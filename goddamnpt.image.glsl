#define GAMME 2.2

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    ivec2 o = ivec2(fragCoord);
    vec4 data = texelFetch(iChannel0, o, 0);
    vec3 col = data.rgb / data.w;
    
    col = pow(col, vec3(1.0 / GAMMA)); // gamma correction
    
    // Output to screen
    fragColor = vec4(col, 1.0);
}
