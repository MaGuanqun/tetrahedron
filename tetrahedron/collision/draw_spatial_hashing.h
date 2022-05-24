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
	void drawCellSelect(Camera* camera, Shader* shader);
	void drawCellSelectOne(Camera* camera, Shader* shader);
	void setCellData(std::vector<unsigned int>& hash_index, double cell_length, unsigned int* cell_number,
		double* initial_aabb);
	void setCellData(unsigned int hash_index, double cell_length, unsigned int* cell_number,
		double* initial_aabb);
private:
	void setCellIndexBasic();
	unsigned int cell_index_basic[24];
	unsigned int VAO1;
	unsigned int VBO1;
	unsigned int EBO1;

	unsigned int VAO2;
	unsigned int VBO2;
	unsigned int EBO2;

	unsigned int VAO3;
	unsigned int VBO3;
	unsigned int EBO3;

	void genBuffer();

	std::vector<double> point_value;
	std::vector<unsigned int> edge_index;

	std::vector<double> point_value_select;
	std::vector<unsigned int> edge_index_select;

	std::vector<double> point_value_select_one;
	std::vector<unsigned int> edge_index_select_one;

	unsigned int cell_num0_cell_num1;

	void reverseHash(unsigned int hash_index, double cell_length, unsigned int* cell_number,
		double* initial_aabb, double* point);
	void obtianHashValue(std::vector<std::vector<unsigned int>>* hash);
	std::unordered_set<unsigned int> hash_index;

};