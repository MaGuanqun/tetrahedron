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

void Object::getAABB(AABB& target, AABB& aabb0, AABB& aabb1, AABB& aabb2)
{
	for (int i = 0; i < 3; ++i) {
		target.min[i] = myMin(aabb0.min[i], aabb1.min[i]);
		target.min[i] = myMin(target.min[i], aabb2.min[i]);

		target.max[i] = myMax(aabb0.max[i], aabb1.max[i]);
		target.max[i] = myMax(target.max[i], aabb2.max[i]);
	}	
}

void Object::getAABB(AABB& target, AABB& aabb0, AABB& aabb1)
{
	for (int i = 0; i < 3; ++i) {
		target.min[i] = myMin(aabb0.min[i], aabb1.min[i]);
		target.max[i] = myMax(aabb0.max[i], aabb1.max[i]);
	}
}

void Object::getAABB(AABB& target, AABB& aabb0, AABB& aabb1, double radius)
{
	for (int i = 0; i < 3; ++i) {
		target.min[i] = myMin(aabb0.min[i], aabb1.min[i])-radius;
		target.max[i] = myMax(aabb0.max[i], aabb1.max[i])+radius;
	}
}