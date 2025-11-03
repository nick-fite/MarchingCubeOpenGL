#version 400 core

in vec2 uv;
in vec3 normal;

uniform sampler2D tex;

void main(void)
{
	vec4 ambientLight = vec4(.1, .1, .1, 1);
	vec4 lightColor = vec4(1, .9, .5, 1);
	vec3 lightDir = vec3(-1, -1, -2);

	// calculate diffuse lighting and clamp between 0 and 1
	float ndotl = clamp(-dot(normalize(lightDir), normalize(normal)), 0, 1); 

	// add diffuse lighting to ambient lighting and clamp a second time
	vec4 lightValue = clamp(lightColor * ndotl + ambientLight, 0, 1);

	// finally, sample from the texuture and multiply in the light.
	gl_FragColor = texture(tex, uv) * lightValue;
}