#pragma once
#include"../external/glm/glm.hpp"
#include"../external/shader.h"
#include"camera.h"
#include"global.h"
#include"../object/cloth.h"
#include"../object/collider.h"
#include"../object/tetrohedron.h"

class PickTriangle
{
public:
	struct PickIndex
	{
		std::vector<int>index;
		int cloth_No;
	};

	PickTriangle();
	void pickTriangle(std::vector<Cloth>* cloth, std::vector<Collider>* collider, Camera* camera, std::vector<std::vector<bool>>& hide, int* triangle_index, int* pos);

private:
	unsigned int picking_FBO;
	unsigned int picking_RBO_color;
	Shader* shader;
	Shader* collider_shader;
	int picking_base_ID;
	void readPixel(std::vector<unsigned char>* pixel_value, unsigned int* FBO, int* pos);
	void initalFBO();
	void writingFBO(std::vector<Cloth>* cloth, std::vector<Collider>* collider, Camera* camera, std::vector<std::vector<bool>>& hide, Shader* shader);
	void decideTriangle(int& triangle_index, int total_cloth_triangle_num, unsigned int* FBO, int* pos);
	void decideFinalIndi(std::vector<Cloth>* cloth, int sum_triangle_index, int* triangle_index);
};


