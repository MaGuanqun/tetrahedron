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
	void draw(Shader* shader, Camera* camera, unsigned int dimension, float line_width);
	bool pickAxes(Shader* shader, Camera* camera, unsigned int& dimension, int* pos);

private:
	unsigned int vertex_num;
	std::vector<std::vector<double>> circle_vertices;

	unsigned int angle_num;

	unsigned int picking_FBO;
	unsigned int picking_RBO_color;

	unsigned int VAO[3], VBO[3];

	void initalFBO();

	void setBuffer();

	//void setArrow(int type, double center);
	void genBuffer();
	std::vector<glm::vec3> circle_color;
	void writingFBO(Camera* camera, Shader* shader);
	void decideAxe(unsigned int& dimension, unsigned int* FBO, int* pos);
	void readPixel(std::vector<unsigned char>* pixel_value, unsigned int* FBO, int* pos);
};

