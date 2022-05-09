#version 460 core

out vec4 FragColor;

uniform vec3 color;
uniform float transparent;

void main()
{
	FragColor=vec4(color, transparent);	
}