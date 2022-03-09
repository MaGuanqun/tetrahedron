#pragma once
#include<vector>
#include"global.h"
#include"../external/glm/glm.hpp"
#include"../external/shader.h"
#include"camera.h"

class ObjectChosenIndicator
{
public:
	ObjectChosenIndicator();
	void updatePosition(double* AABB);
	void draw(Shader* shader, Camera* camera);

private:
	unsigned int vertex_num;
	std::vector<std::vector<double>> circle_vertices;

	unsigned int angle_num;



	unsigned int VAO[3], VBO[3];

	void setBuffer();

	//void setArrow(int type, double center);
	void genBuffer();
	std::vector<glm::vec3> circle_color;
};

