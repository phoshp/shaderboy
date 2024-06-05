// https://www.pcg-random.org/
uint pcg(uint v)
{
	uint state = v * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    uvec2 seed = uvec2(fragCoord + iTime * 100.0);
    
    fragColor = vec4(vec3(float(pcg(seed.x ^ seed.y)) * (1.0 / float(0xffffffffu))), 1.0);
}
