#ifndef DRAW_TRIANGLE_H
#define DRAW_TRIANGLE_H
#include"../external/shader.h"
#include<vector>
#include<array>
#include"camera.h"
#include"global.h"

class DrawTriangle
{
public:
	DrawTriangle();
	void drawTriangle(Camera* camera, Shader* object_shader_front, std::vector<std::array<double, 3>>& vertex_for_render,
		std::vector<std::array<int, 3>>& triangle_vertex_index, std::vector<std::array<double, 3>>& norm_for_render, std::vector<unsigned int>& triangle_index, glm::vec3 color);

private:
	unsigned int VAO;
	unsigned int VBO[2];
	unsigned int EBO;

	void genBuffer();
	void setBuffer(std::vector<int>& triangle_vertex_index, std::vector<std::array<double, 3>>& vertex_for_render,
		std::vector<std::array<double, 3>>& norm_for_render);
};




#endif
