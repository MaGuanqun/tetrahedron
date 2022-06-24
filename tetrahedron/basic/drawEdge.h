#pragma once
#include<vector>
#include<array>
#include"../external/shader.h"
#include"camera.h"


class Draw_Edge
{
public:
	Draw_Edge();

	void drawEdge(Camera* camera, Shader* wireframe_shader, std::vector<std::array<double, 3>>& vertex_for_render,
		std::vector<unsigned int>& edge_vertex_index, glm::vec3 color);
	void drawEdge(Camera* camera, Shader* wireframe_shader, std::vector<std::array<double, 3>>& vertex_for_render,
		std::vector<unsigned int>& edge_index, glm::vec3 color, std::vector<unsigned int>& edge_vertex_index);
private:

	unsigned int EE_VAO;
	unsigned int EE_VBO;
	unsigned int EE_EBO;

	Shader* shader;
	void genBuffer();

	void setBuffer(std::vector<std::array<double, 3>>&vertex_for_render,
		std::vector<unsigned int>& edge_vertex_index);
};


