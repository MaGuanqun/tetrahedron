#ifndef MOUSE_SELECT_H
#define MOUSE_SELECT_H

#include<vector>
#include"../external/glm/glm.hpp"
#include"../external/glm/gtc/matrix_transform.hpp"
#include"../external/shader.h"
#include"global.h"
class MouseSelect
{
public:
	MouseSelect();
	void draw();
	void setFirstCorner(float x, float y);
	void setLastCorner(float x, float y);
private:
	glm::vec3 camera_pos;
	unsigned int VAO, VBO;
	std::vector<float>position;
	Shader* shader;
	void genBuffer();
	void setBufferData();
	glm::mat4 view;
	glm::mat4 projection;
	float first_corner[3];
	float last_corner[3];
	void transferTo3D(float x, float y, float* point);
	void setPosition();
};






#endif // !MOUSE_SELECT_H


#pragma once
