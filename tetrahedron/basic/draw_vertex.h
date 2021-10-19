#ifndef DRAW_VERTEX_H
#define DRAW_VERTEX_H
#include"../external/shader.h"
#include<vector>
#include<array>
#include"camera.h"
#include"global.h"
class DrawVertex
{
public:
	DrawVertex();
	void initialVertex(int vertex_num);
	void draw(Camera* camera, glm::vec3& color, float transparence = 1.0f);
	void setVertex(std::vector<std::array<double,3>>&v_list, std::vector<bool>& index);
	void setVertex(std::vector<std::array<double,3>>&v_list, std::vector<int>& index);
	void setVertex(std::vector<std::array<double, 3>>& v_list, double tolerance);	
	void setVertex(double* vertex, double tolerance);
	void setVertexAccumulate(std::vector<std::array<double, 3>>& v_list, std::vector<int>& index);
private:
	Shader* shader;
	std::vector<std::array<double, 3>> sphere;
	std::vector<std::array<double, 3>> sphere_normal;
	std::vector<int>basic_index;
	void setBasicSphere();
	void genBuffer();
	unsigned int VAO, VBO[2], EBO;
	void setBuffer();	
	std::vector<int>indices;
	int draw_vertex_num;
	std::vector<std::array<double, 3>> draw_vertex;
	std::vector<std::array<double, 3>> draw_normal;
	int sphere_vertex_num;
	int longitude_num = 10;
	int latitude_num = 10;
	Light light;
	void setSphereRadius(double R);
};


#endif // !DRAW_VERTEX_H





#pragma once
