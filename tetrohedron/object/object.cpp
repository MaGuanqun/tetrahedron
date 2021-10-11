#include"object.h"

void Object::genBuffer()
{
	glGenVertexArrays(1, &VAO);
	glGenBuffers(3, VBO);
	glGenBuffers(1, &EBO);
}


void Object::genShader()
{
	object_shader_front = new Shader("./shader/object.vs", "./shader/object.fs");
	*object_shader_front = Shader("./shader/object.vs", "./shader/object.fs");
	object_shader_back = new Shader("./shader/object.vs", "./shader/object.fs");
	*object_shader_back = Shader("./shader/object.vs", "./shader/object.fs");
	wireframe_shader = new Shader("./shader/wireframe.vs", "./shader/wireframe.fs");
	*wireframe_shader = Shader("./shader/wireframe.vs", "./shader/wireframe.fs");
}

void Object::setWireframwColor(double* color)
{
	wireframe_color = glm::vec3(color[0], color[1], color[2]);
}