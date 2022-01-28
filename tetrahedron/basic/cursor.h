#ifndef CURSOR_H
#define CURSOR_H
#include<vector>
#include"../external/shader.h"
#include"camera.h"
#include"global.h"
#include<array>
#undef DEG_2_RADIANS
#define DEG_2_RADIANS(a) (a/180.0*3.141592653589793)

class Cursor
{
public:
	Cursor();
	void draw(Camera* camera);
	void createVertices(double radius, double camera_center[3]);
	void translate(double u[3], double cursor_pos[3]);
private:
	double R;	
	double center[3];
	int longitude_num;
	int latitude_num;
	int vertice_num;
	std::vector<std::array<double, 3>>ori_vertices_pos;
	std::vector<std::array<double,3>>vertices_pos;
	std::vector<std::array<double,3>>cursor_vertices_pos;
	std::vector<int>indices;
	std::vector<std::array<double, 3>>normal;

	void getNormal();
	void genBuffer();
	void setBufferData();
	unsigned int VAO1, VBO1[2], EBO1;
	unsigned int VAO2, VBO2[2], EBO2;
	unsigned int VAO3, VBO3;
	Shader* shader;
	Light light;
	std::vector<std::array<double, 3>>lines;	
};

#endif // !CURSOR_H

#pragma once
