float vdf[8 * 8];

void build_voxel_df() {
    // hardcoded voxel distance field, zero means filled
    vdf = float[](
        0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
        1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0,
        1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0,
        1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0
    );

            
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            if (vdf[i + j * 8] == 0.0) continue;
            float mind = 1000.0;
            for (int x = 0; x < 8; x++)
                for (int y = 0; y < 8; y++)
                    if (vdf[x + y * 8] == 0.0) mind = min(mind, floor(distance(vec2(x, y), vec2(i, j))));
            vdf[i + j * 8] = mind;
        }
    }
}

float trace(vec2 origin, vec2 dir, out int steps) {
    const float unit = 1.0 / 8.0;
    origin *= 8.0;
    vec2 pos = floor(origin);
	vec2 rs = sign(dir);
    vec2 ri = 1.0 / dir;
    
	vec2 dis = (pos-origin + 0.5 + rs*0.5) * ri;
    vec2 mm = vec2(0.0);
    steps = 0;
	for (int i = 0; i < 32; i++) {
        float df = vdf[int(pos.x) + int(pos.y) * 8];
		if (df == 0.0) break;
        
        mm = step(dis.xy, dis.yx);
        dis += mm * rs * ri;
        pos += mm * rs;
        steps += 1;
    }
    vec2 mini = (floor(pos) - origin + 0.5 - 0.5 * rs) * ri;
    return max(mini.x, mini.y) * unit + 0.0001;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    build_voxel_df();
    vec2 invR = 1.0 / iResolution.xy;
    vec2 uv = fragCoord.xy * invR;
    float aspect = iResolution.x / iResolution.y;
    vec2 pixel = vec2(uv.x, uv.y * aspect);
    
    vec2 ray_orig = vec2(15.0 * invR.x, 0.5);
    vec3 col = vec3(0);
    ivec2 pos = ivec2(floor(vec2(uv.x, uv.y) * 8.0));
    if (pos.x >= 0 && pos.x < 8 && pos.y >= 0 && pos.y < 8) {
        float d = max(0.0, 1.0-vdf[pos.x + pos.y * 8]);
        col += (0.5 + 0.5 * cos(iTime + uv.xyx + vec3(0,2,4))) * d;
    }

    col += vec3(1.0, 0.0, 0.0) * (distance(ray_orig, uv) <= 0.009 ? 1.0 : 0.0);

    vec2 ray_dir = normalize(vec2(sin(iTime), cos(iTime)));
    int steps;
    float dist = trace(ray_orig, ray_dir, steps);
    
    // Draw Line: https://www.shadertoy.com/view/4ljfRD 
    vec2 p1 = ray_orig;
    vec2 p2 = ray_orig + ray_dir * dist;
    
    float step_size = dist / float(steps);
    for (float i = 0.0; i < dist; i += step_size) {
        col += vec3(0.0, 1.0, 0.0) * (distance(ray_orig + ray_dir * i, uv) <= 0.008 ? 1.0 : 0.0);
    }
    
    col += vec3(1.0, 0.0, 0.0) * (distance(p2, uv) <= 0.01 ? 1.0 : 0.0);
    
    float thickness = 1.0;
    float t1 = distance(p1, p2);
    float t2 = distance(p1, uv);
    float t3 = 1.0 - floor(1.0 - (thickness * invR.x) + distance(mix(p1, p2, clamp(t2 / t1, 0.0, 1.0)), uv));
    
    col += vec3(1) * max(0.0, t3);

    fragColor = vec4(col,1.0);
}
