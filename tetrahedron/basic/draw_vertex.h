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
	
	void setInObjNum(unsigned int obj_num);

	void draw(Camera* camera, glm::vec3& color, float transparence = 1.0f);
	void setVertex(std::vector<std::array<double,3>>&v_list, std::vector<bool>& index);
	void setVertex(std::vector<std::array<double,3>>&v_list, std::vector<int>& index);
	void setVertex(std::vector<std::array<double, 3>>& v_list, double tolerance);	
	void setVertex(double* vertex, double tolerance);
	void setVertexAccumulate(std::vector<std::array<double, 3>>& v_list, std::vector<int>& index);

	void setCollisionVertexData(std::vector<std::array<double, 3>*>& v_list, std::vector<std::vector<unsigned int>>& vertex_index);
	void setCollisionVertexDataAllCell(std::vector<std::array<double, 3>*>& v_list, std::vector<std::vector<unsigned int>>& vertex_index);
	void setCollisionVertexDataInACell(std::vector<std::array<double, 3>*>& v_list, std::vector<std::vector<unsigned int>>& vertex_index);
	
	void setShaderData(Camera* camera);

	void drawAllCollisionVertex(unsigned int obj_index, glm::vec3 color, float transparence);
	void drawAllCellVertex(unsigned int obj_index, glm::vec3 color, float transparence);
	void drawOneCellVertex(unsigned int obj_index, glm::vec3 color, float transparence);
	

private:
	Shader* shader;
	std::vector<std::array<double, 3>> sphere;
	std::vector<std::array<double, 3>> sphere_normal;
	std::vector<int>basic_index;
	void setBasicSphere();
	void genBuffer();
	unsigned int VAO, VBO[2], EBO;

	std::vector<unsigned int> VAO1, VBO1, EBO1;
	std::vector<unsigned int> VAO_all_cell, VBO_all_cell, EBO_all_cell;
	std::vector<unsigned int> VAO_in_a_cell, VBO_in_a_cell, EBO_in_a_cell;

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

	std::vector<std::vector<std::array<double, 3>>> vertex_pos;
	std::vector < std::vector<int>>vertex_sphere_indices;
	std::vector<std::vector<std::array<double, 3>>> normal;

	std::vector<std::vector<std::array<double, 3>>> vertex_pos_in_a_cell;
	std::vector < std::vector<int>>vertex_sphere_indices_in_a_cell;
	std::vector<std::vector<std::array<double, 3>>> normal_in_a_cell;

	std::vector<std::vector<std::array<double, 3>>> vertex_pos_all_cell;
	std::vector < std::vector<int>>vertex_sphere_indices_all_cell;
	std::vector<std::vector<std::array<double, 3>>> normal_all_cell;


	void setBuffer1();
	void genBuffer1();
	void setBuffer1(unsigned int obj_No, unsigned int* VAO1, unsigned int* VBO1, unsigned int* EBO1, std::vector<std::array<double, 3>>* vertex_pos,
		std::vector<int>* vertex_sphere_indices, std::vector<std::array<double, 3>>* normal);

	void setBufferAllCell();
	void setBufferOneCell();

	void setVertexAccumulate(std::vector<std::array<double, 3>*>& v_list, std::vector<std::vector<unsigned int>>& vertex_index,
		std::vector < std::array<double, 3>>* vertex_pos, std::vector<int>* vertex_sphere_indices, std::vector<std::array<double, 3>>* normal);
	void initialVertex(unsigned int obj_num, std::vector<std::vector<unsigned int>>& vertex_index,
		std::vector<std::array<double, 3>>* vertex_pos, std::vector < std::array<double, 3>>* normal, std::vector<int>* vertex_sphere_indices);
	void drawCollisionVertex(unsigned int obj_index, glm::vec3 color, float transparence, std::vector<int>* vertex_sphere_indices,
		unsigned int* VAO1);
};


#endif // !DRAW_VERTEX_H





#pragma once
