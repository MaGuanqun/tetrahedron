#pragma once
#include"mesh_struct.h"
#include<algorithm>
#include<random>

class GraphColor
{
public:
	
	void graphColorBending(MeshStruct& mesh_struct);
	void graphColorEdgeLength(MeshStruct& mesh_struct);
	void graphColor(std::vector<std::vector<unsigned int>>& element_element, std::vector<std::vector<unsigned int>>& element_not_connect);



	void testTet(MeshStruct& mesh_struct, std::vector<std::vector<unsigned int>>& element_not_connect, std::vector<std::array<int, 4>>& indices);
	void testEdge(std::vector<std::vector<unsigned int>>&element_not_connect, MeshStruct& mesh_struct, 
		std::vector<unsigned int>& edge_vertices);
	void testBend(std::vector<std::vector<unsigned int>>& element_not_connect, MeshStruct& mesh_struct);


	void graphColorTet(MeshStruct& mesh_struct, int different_color_strategy_num);

private:

	void graphColorLoopNode(std::vector<std::vector<unsigned int>>& element_element, int max_color_number,
		std::vector<std::vector<unsigned int>>& element_not_connect, int* element_order, std::vector<int>&color_score);
	
	void decideGroup(unsigned int max_color, std::vector<std::vector<unsigned int>>& unconnected_vertex_index, int* color, int size);
	void testEdgeGroup(int size, std::vector<std::vector<unsigned int>>& unconnected_vertex_index, MeshStruct& mesh_struct);

	void findMinMaxDegree(std::vector<std::vector<unsigned int>>& element_element, unsigned int& max_degree, unsigned int& min_degree);

	void getMaxPaletteSize(int& size, int* palette, unsigned int max_size, unsigned int max_array_size);

	void orderByDegree(std::vector<std::vector<unsigned int>>& tet_tet, std::vector<int>& element_order);
	void orderByScore(std::vector<int>& tet_score, std::vector<int>& element_order);

};
