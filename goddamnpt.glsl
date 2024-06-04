#define PI 3.14159
#define SAMPLES 5.0
#define MAX_RENDER_DISTANCE 10000.0
#define EPSILON 0.001
#define PATH_LENGTH 12

uint baseHash( uvec2 p ) {
    p = 1103515245U*((p >> 1U)^(p.yx));
    uint h32 = 1103515245U*((p.x)^(p.y>>3U));
    return h32^(h32 >> 16);
}

float hash1( inout float seed ) {
    uint n = baseHash(floatBitsToUint(vec2(seed+=.1,seed+=.1)));
    return float(n)/float(0xffffffffU);
}

vec2 hash2( inout float seed ) {
    uint n = baseHash(floatBitsToUint(vec2(seed+=.1,seed+=.1)));
    uvec2 rz = uvec2(n, n*48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU))/float(0x7fffffff);
}

mat3 GetTangentSpace(vec3 normal)
{
    vec3 helper = vec3(1, 0, 0);
    if (abs(normal.x) > 0.99) helper = vec3(0, 0, 1);
    vec3 tangent = normalize(cross(normal, helper));
    vec3 binormal = normalize(cross(normal, tangent));
    return mat3(tangent, binormal, normal);
}

vec3 SampleHemisphere(float seed, vec3 normal, float alpha) {
	float cosTheta = pow(hash1(seed), 1.0 / (alpha + 1.0));
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
    float phi = 2.0 * PI * hash1(seed);
    vec3 tangentSpaceDir = vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
    
    return GetTangentSpace(normal) * tangentSpaceDir;
}

float CosineSamplingPDF(float NdotL) {
    return NdotL / PI;
}

vec3 FresnelSchlick(float cosTheta, vec3 F0) {
    return clamp(F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0), 0.0, 1.0);
}

float DistributionGGX(vec3 N, vec3 H, float roughness) {
    float a = roughness * roughness;
    float a2 = a * a;
    float NdotH = max(dot(N, H), 0.0);
    float NdotH2 = NdotH * NdotH;

    float nom = a2;
    float denom = max(NdotH2 * (a - 1.0) + 1.0, 0.0000001);
    denom = PI * denom * denom;

    return nom / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness) {
    float r = (roughness + 1.0);
    float k = (r * r) / 8.0;

    float nom = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return nom / denom;
}

float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness) {
    float NdotV = abs(dot(N, V));
    float NdotL = abs(dot(N, L));
    float ggx2 = GeometrySchlickGGX(NdotV, roughness);
    float ggx1 = GeometrySchlickGGX(NdotL, roughness);

    return ggx1 * ggx2;
}

vec3 ImportanceSampleGGX(float x, float y, vec3 N, vec3 V, float roughness) {
    float a = roughness * roughness;

    float phi = 2.0 * PI * x;
    float cosTheta = sqrt((1.0 - y) / (1.0 + (a * a - 1.0) * y));
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

    // from spherical coordinates to cartesian coordinates
    vec3 H;
    H.x = cos(phi) * sinTheta;
    H.y = sin(phi) * sinTheta;
    H.z = cosTheta;

    // from tangent-space vector to world-space sample vector
    vec3 up = abs(N.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
    vec3 tangent = normalize(cross(up, N));
    vec3 bitangent = cross(N, tangent);

    vec3 halfVec = tangent * H.x + bitangent * H.y + N * H.z;
    halfVec = normalize(halfVec);
    
    return halfVec;

}

float ImportanceSampleGGX_PDF(float NDF, float NdotH, float VdotH) {
    return NDF * NdotH / (4.0 * VdotH);
}

float CalculateFresnel(vec3 I, vec3 N, float ior) {
    float kr;
    float cosi = clamp(-1.0, 1.0, dot(I, N));
    float etai = 1.0;
	float etat = ior;
    if (cosi > 0.0) {
        float temp = etai;
        etai = etat;
        etat = temp;
    }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrt(max(0.0, 1.0 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1.0) {
        kr = 1.0;
    } else {
        float cost = sqrt(max(0.f, 1.0 - sint * sint));
        cosi = abs(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        kr = (Rs * Rs + Rp * Rp) / 2.0;
    }
    return kr;
}

vec3 RefractionBTDF(float D, float G, vec3 F, vec3 V, vec3 L, vec3 N, vec3 H, float etaIn, float etaOut) {  
    float NdotL = abs(dot(N, L));
    float NdotV = abs(dot(N, V));
            
    float VdotH = abs(dot(V, H));
    float LdotH = abs(dot(L, H));
            
    
    float term1 = VdotH * LdotH / (NdotV * NdotL);
    vec3 term2 = etaOut * etaOut * (1.0 - F) * G * D;
    float term3 = (etaIn * VdotH + etaOut * LdotH) * (etaIn * VdotH + etaOut * LdotH) + 0.001f;
    vec3 refractionBrdf = term1 * term2 / term3;
    
    return refractionBrdf;
}

vec3 SpecularBRDF(float D, float G, vec3 F, vec3 V, vec3 L, vec3 N) {        
    float NdotL = abs(dot(N, L));
    float NdotV = abs(dot(N, V));
    vec3 nominator = D * G * F;
    float denominator = 4.0 * NdotV * NdotL + 0.001;
    vec3 specularBrdf = nominator / denominator;
    
    return specularBrdf;
}

vec3 DiffuseBRDF(vec3 albedo) {
    return albedo / PI;
}

vec3 Reinhard(vec3 color) {
    return color / (color + 1.0);
}

#define SUN_DIR normalize(vec3(-0.4,0.3,-0.6))

vec3 getSkyColor(vec3 rd) {
    vec3 col = mix(vec3(1),vec3(1,0.7,0.5), 0.5+0.5*rd.y);
    float sun = clamp(dot(SUN_DIR,rd), 0.0, 1.0);
    float sun2 = pow(sun, 32.0);
    col += vec3(1,0.6,0.1) * (pow(sun2, 12.0) + 2.0*pow(sun2, 16.0));
    return max(col, vec3(0));
}

struct Material {
    vec3 albedo;
    vec4 emissive;
    float metallic;
    float roughness;
    vec3 specular;
    float specTrans;
    float ior;
};


struct Sphere {
    vec3 pos;
    float radius;
    vec3 color;
    uint mat;
};

struct Plane {
    vec3 pos;
    vec3 normal;
    vec3 color;
    uint mat;
};

struct Box {
    vec3 min;
    vec3 max;
    vec3 color;
    uint mat;
};

struct Ray {
    vec3 origin;
    vec3 dir;
};

struct HitData {
    bool hit;
    vec3 pos;
    float dist;
    vec3 normal;
    vec3 albedo;
    Material mat;
};

Material[4] g_materials;
Plane[1] g_planes;
Sphere[3] g_spheres;
Box[1] g_boxes;

HitData createHitData() {
    HitData hit;
    hit.hit = false;
    hit.pos = vec3(0);
    hit.dist = MAX_RENDER_DISTANCE;
    hit.normal = vec3(1.0);
    hit.albedo = vec3(0.0);
    
    return hit;
}

Ray setCamera(vec3 orig, vec3 forward, vec2 uv) {
    const float focal_dist = 1.2;
    const vec3 up_vec = vec3(0, 1, 0);
    float aspect_ratio = iResolution.x / iResolution.y;
    
    float i = (uv.x - 0.5) * focal_dist * aspect_ratio;
    float j = (uv.y - 0.5) * focal_dist;
    
    vec3 right = cross(forward, up_vec);
    vec3 up = cross(right, forward);
    vec3 dir = normalize(forward + (right * i) + (up * j));
    
    return Ray(orig, dir);
}

void hitMaterial(Ray ray, inout HitData hit, float dist, vec3 color, uint matId) {
    hit.hit = true;
    hit.pos = ray.origin + ray.dir * dist;
    hit.dist = dist;
    hit.mat = g_materials[matId];
    hit.albedo = color * hit.mat.albedo;
}

void hitSphere(Ray ray, inout HitData hit, Sphere sphere) {
    vec3 base = sphere.pos - ray.origin;
    float d1 = dot(base, ray.dir);
    float d2Sqrt = d1 * d1 - dot(base, base) + sphere.radius * sphere.radius;
    if (d2Sqrt < 0.0) {
        return;
    }
    float d2 = sqrt(d2Sqrt);
    float t = (d1 - d2) > 0.0 ? d1 - d2 : d1 + d2;
    if (t > EPSILON && t < hit.dist) {
        hitMaterial(ray, hit, t, sphere.color, sphere.mat);
        hit.normal = normalize((hit.pos - sphere.pos) / sphere.radius);
    }
}

void hitPlane(Ray ray, inout HitData hit, Plane plane) {
    float d = dot(plane.normal, ray.dir);
    if (d != 0.0) {
        float t = dot(plane.pos - ray.origin, plane.normal) / d;
        if (t > EPSILON && t < hit.dist) {
            hitMaterial(ray, hit, t, plane.color, plane.mat);
            hit.normal = plane.normal;
            hit.albedo = hit.mat.albedo * texture(iChannel1, hit.pos.xz/ 2.0, 0.0).rgb;
        }
    }
}

float safe_inv(float v) {
    return v == 0.0 ? v : 1.0 / v;
}

vec3 safe_inv3(vec3 v) {
    return vec3(safe_inv(v.x), safe_inv(v.y), safe_inv(v.z));
}

void hitBox(Ray ray, inout HitData hit, Box box) {
    vec3 t1 = safe_inv3(ray.dir) * (box.min - ray.origin);
    vec3 t2 = safe_inv3(ray.dir) * (box.max - ray.origin);
    vec3 tminv = min(t1, t2);
    vec3 tmaxv = max(t1, t2);

    float tmin = max(max(tminv.x, tminv.y), tminv.z);
    float tmax = min(min(tmaxv.x, tmaxv.y), tmaxv.z);

    if (tmin <= 0.0) tmin = tmax; // we are inside of box
    
    if (tmax > max(tmin, 0.0) && tmin < hit.dist) {
        hitMaterial(ray, hit, tmin, box.color, box.mat);
        vec3 norm = -sign(ray.dir) * step(tminv.yzx, tminv.xyz) * step(tminv.zxy, tminv.xyz);
        hit.normal = norm;
    }
}

HitData hitScene(Ray ray) {
    HitData hit = createHitData();
    for (int i = 0; i < g_spheres.length(); ++i) {
        hitSphere(ray, hit, g_spheres[i]);
    }
    for (int i = 0; i < g_planes.length(); ++i) {
        hitPlane(ray, hit, g_planes[i]);
    }
    for (int i = 0; i < g_boxes.length(); ++i) {
        hitBox(ray, hit, g_boxes[i]);
    }
    return hit;
}

vec3 tracePath(Ray ray, bool daytime, inout float seed, inout vec3 radi, inout vec3 irradi) {
    vec3 direct = vec3(0.0);
    vec3 energy = vec3(1.0);
    bool is_direct = true;
    
    for (int i = 0; i < PATH_LENGTH; ++i) {
        HitData hit = hitScene(ray);
        
        if (!hit.hit) {
            vec3 sky = daytime ? getSkyColor(ray.dir) : vec3(0);
            radi += sky;
            irradi += sky * energy;
            return direct;
        }
        
        vec4 e = hit.mat.emissive;
        radi += e.rgb * e.w * e.w;
        
        float roulette = hash1(seed);
        float blender = hash1(seed); //used to blend BSDF and BRDF
        
        if (blender < 1.0 - hit.mat.specTrans) {
            if (is_direct) {
                direct += hit.albedo * energy;
            } else {
                irradi += hit.albedo * energy;
            }
            vec3 reflectionDir;
            
            float diffuseRatio = 0.5 * (1.0 - hit.mat.metallic);
            float specularRatio = 1.0 - diffuseRatio;
            vec3 V = normalize(-ray.dir);
            
            if (roulette < diffuseRatio) { //sample diffuse
                reflectionDir = SampleHemisphere(seed, hit.normal, 1.0);
            } else {
                vec3 halfVec = ImportanceSampleGGX(hash1(seed), hash1(seed), hit.normal, V, hit.mat.roughness);
                reflectionDir = reflect(ray.dir, halfVec);
                reflectionDir = normalize(reflectionDir);
            }
            
            vec3 L = normalize(reflectionDir);
            vec3 H = normalize(V + L);
        
            float NdotL = abs(dot(hit.normal, L));
            float NdotH = abs(dot(hit.normal, H));
            float VdotH = abs(dot(V, H));
            
            float NdotV = abs(dot(hit.normal, V));
            
            vec3 F0 = mix(hit.mat.specular, hit.albedo, hit.mat.metallic);
        
            float NDF = DistributionGGX(hit.normal, H, hit.mat.roughness);
            float G = GeometrySmith(hit.normal, V, L, hit.mat.roughness);
            vec3 F = FresnelSchlick(max(dot(H, V), 0.0), F0);
        
            vec3 kS = F;
            vec3 kD = (1.0 - kS) * (1.0 - hit.mat.metallic);
        
            vec3 specularBrdf = SpecularBRDF(NDF, G, F, V, L, hit.normal);
            float specularPdf = ImportanceSampleGGX_PDF(NDF, NdotH, VdotH);
            vec3 diffuseBrdf = DiffuseBRDF(hit.albedo);
            float diffusePdf = CosineSamplingPDF(NdotL);

            vec3 totalBrdf = (diffuseBrdf * kD + specularBrdf) * NdotL;
            float totalPdf = diffuseRatio * diffusePdf + specularRatio * specularPdf;
                
            ray.origin = hit.pos + hit.normal * EPSILON;
            ray.dir = reflectionDir;
            if (totalPdf > 0.0) {
                energy *= totalBrdf / totalPdf;
            }
        } else {
            bool fromOutside = dot(ray.dir, hit.normal) < 0.0;
            vec3 N = fromOutside ? hit.normal : -hit.normal;
            vec3 bias = N * EPSILON;
            
            float etai = 1.0;
            float etat = hit.mat.ior;
            
            vec3 V = normalize(-ray.dir);
            vec3 H = ImportanceSampleGGX(hash1(seed), hash1(seed), N, V, hit.mat.roughness);
            
            vec3 F0 = hit.mat.specular;
            vec3 F = FresnelSchlick(max(dot(H, V), 0.0), F0);
            
            float kr = CalculateFresnel(ray.dir, hit.normal, hit.mat.ior);
            
            float specularRatio = kr;
            float refractionRatio = 1.0 - kr;
            
            vec3 L; 
            if (roulette <= specularRatio) {
                ray.origin = hit.pos + bias;
                L = normalize(reflect(ray.dir, H));
                ray.dir = L;
            } else {
                float eta = fromOutside ? etai / etat : etat / etai;
                L = normalize(refract(ray.dir, H, eta));
                ray.origin = hit.pos - bias;
                ray.dir = L;
                L = N;
                if (!fromOutside) {
                    //since the BTDF is not reciprocal, we need to invert the direction of our vectors.
                    vec3 temp = L;
                    L = V;
                    V = temp;
                        
                    N = -N;
                    H = -H;
                }
            }
            
            float NdotL = abs(dot(N, L));
            float NdotV = abs(dot(N, V));
            
            float NdotH = abs(dot(N, H));
            float VdotH = abs(dot(V, H));
            float LdotH = abs(dot(L, H));
            
            float NDF = DistributionGGX(N, H, hit.mat.roughness);
            float G = GeometrySmith(N, V, L, hit.mat.roughness);

            vec3 specularBrdf = SpecularBRDF(NDF, G, F, V, L, N);
            float specularPdf = ImportanceSampleGGX_PDF(NDF, NdotH, VdotH);
            
            //refraction
            float etaOut = etat;
            float etaIn = etai;
            
            vec3 refractionBtdf = RefractionBTDF(NDF, G, F, V, L, N, H, etaIn, etaOut);
            float refractionPdf = ImportanceSampleGGX_PDF(NDF, NdotH, VdotH);
            
            //BSDF = BRDF + BTDF
            vec3 totalBrdf = (specularBrdf + refractionBtdf * hit.albedo) * NdotL;
            float totalPdf = specularRatio * specularPdf + refractionRatio * refractionPdf;
            if (totalPdf > 0.0) {
                energy *= totalBrdf / totalPdf;
            }
        }
        is_direct = false;
    }
    return direct;
}

void initScene() {
    Material mate1;
    mate1.albedo = vec3(1.0);
    mate1.emissive = vec4(0);
    mate1.metallic = 0.0;
    mate1.roughness = 0.2;
    mate1.specular = vec3(0.2);
    mate1.specTrans = 1.0;
    mate1.ior = 1.05;
    
    Material mate2;
    mate2.albedo = vec3(1.0);
    mate2.emissive = vec4(0);
    mate2.metallic = 0.6;
    mate2.roughness = 0.35;
    mate2.specular = vec3(0.2);
    mate2.specTrans = 0.0;
    mate2.ior = 1.0;
    
    Material mate3;
    mate3.albedo = vec3(1.0);
    mate3.emissive = vec4(0);
    mate3.metallic = 1.0;
    mate3.roughness = 0.01;
    mate3.specular = vec3(1.0);
    mate3.specTrans = 0.0;
    mate3.ior = 1.0;
    
    Material mate4;
    mate4.albedo = vec3(1.0);
    mate4.emissive = vec4(1, 1, 1, 2);
    mate4.metallic = 0.0;
    mate4.roughness = 1.0;
    mate4.specular = vec3(0.0);
    mate4.specTrans = 0.0;
    mate4.ior = 1.0;
    
    g_materials[0] = mate1;
    g_materials[1] = mate2;
    g_materials[2] = mate3;
    g_materials[3] = mate4;
    
    g_planes[0] = Plane(vec3(0, 0, 0), vec3(0, 1, 0), vec3(1, 0, 0), 1u);
    g_spheres[0] = Sphere(vec3(0, 5.01, 0.0), 5.0, vec3(0.9, 0.1, 0.1), 0u);
    g_spheres[1] = Sphere(vec3(-14.3, 6.01, 8.7), 6.0, vec3(1, 0.1, 0), 2u);
    g_spheres[2] = Sphere(vec3(13.2, 4.01, -6.2), 4.0, vec3(1, 1, 1), 3u);
    
    g_boxes[0] = Box(vec3(-15, 0, -15), vec3(-5, 10, -5), vec3(0.4, 0.5, 0.02), 0u);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 p = (-iResolution.xy + 2.*fragCoord - 1.)/iResolution.y;
    float seed = float(baseHash(floatBitsToUint(p - iTime)))/float(0xffffffffU);
    initScene();
    
    const vec3 scene_center = vec3(0.0, 5.0, 0.0);
    const float camera_dist = 30.0;
    vec2 rot = abs(iMouse.xy) / iResolution.xy - 0.5;
    vec3 ss = vec3(cos(1.5 + 6.0 * rot.x), 1.0 + 2.0 * rot.y, sin(1.5 + 6.0 * rot.x));
    vec3 origin = scene_center + camera_dist * ss;
    vec3 forward = normalize(scene_center - origin);
    
    bool reset = false;
    vec4 data = texelFetch(iChannel0, ivec2(0), 0);
    float daytime = mod(floor(iTime / 5.0), 2.0);
    
    if (daytime != data.w
        || round(rot * iResolution.xy) != round(data.xy)
        || round(data.z) != round(iResolution.x)) {
        reset = true;
    }
    
    if(all(equal(ivec2(fragCoord), ivec2(0)))) {
		fragColor = vec4(rot * iResolution.xy, iResolution.x, daytime);
        return;
    }
    
    vec3 color = vec3(0);
    vec3 radi = vec3(0);
    vec3 irradi = vec3(0);
    for (float i = 0.0; i < SAMPLES; ++i) {
        Ray ray = setCamera(origin, forward, (fragCoord + hash2(seed)) / iResolution.xy);
        
        vec3 r = vec3(0);
        vec3 ir = vec3(0);
        vec3 t = tracePath(ray, daytime == 0.0, seed, r, ir);
        radi += r;
        irradi += ir;
        color += t;
    }
    color /= SAMPLES;
    radi /= SAMPLES;
    irradi /= SAMPLES;
    
    color = (color + irradi) * radi;
    
    if (reset) {
       fragColor = vec4(color, 1.0); 
    } else {
       fragColor = vec4(color, 1.0) + texelFetch(iChannel0, ivec2(fragCoord), 0);
    }
}
