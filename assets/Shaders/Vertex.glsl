#version 400 core

// Vertex attribute for position
layout(location = 0) in vec3 in_position;
layout(location = 1) in vec2 in_uv;
layout(location = 2) in vec3 in_normal;

// uniform will contain the world matrix.

uniform mat4 worldMatrix;
uniform mat4 cameraView;

out vec2 uv;
out vec3 normal;

void main(void)
{
	//transform the vector
	vec4 worldPosition = worldMatrix * vec4(in_position, 1);
	vec4 viewPosition = cameraView * worldPosition;

	// output the transformed vector
	gl_Position = viewPosition;
	normal = mat3(worldMatrix) * in_normal;
	uv = in_uv;
}