#pragma once
#include"spatial_hashing.h"
#include <unordered_set>

class DrawSpatialHashing
{
public:
	DrawSpatialHashing();
	void setCellData(std::vector<std::vector<unsigned int>>* hash, double cell_length, unsigned int* cell_number,
		double* initial_aabb);
	void drawCell(Camera* camera, Shader* shader);
private:
	void setCellIndexBasic();
	unsigned int cell_index_basic[24];
	unsigned int VAO1;
	unsigned int VBO1;
	unsigned int EBO1;
	void genBuffer();

	std::vector<double> point_value;
	std::vector<unsigned int> edge_index;

	unsigned int cell_num0_cell_num1;

	void reverseHash(unsigned int hash_index, double cell_length, unsigned int* cell_number,
		double* initial_aabb, double* point);
	void obtianHashValue(std::vector<std::vector<unsigned int>>* hash);
	std::unordered_set<unsigned int> hash_index;

};