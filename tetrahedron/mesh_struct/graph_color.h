#pragma once
#include"mesh_struct.h"

class GraphColor
{
public:
	
	void graphColorBending(MeshStruct& mesh_struct);
	void graphColorEdgeLength(MeshStruct& mesh_struct);
	void graphColor(std::vector<std::vector<unsigned int>>& element_element, std::vector<std::vector<unsigned int>>& element_not_connect);


	void testEdge(MeshStruct& mesh_struct, std::vector<std::array<int, 4>>& indices);

private:
	void findEdgeLengthMinMaxDegree(MeshStruct& mesh_struct, unsigned int& max_degree, unsigned int& min_degree);
	void findBendingMinMaxDegree(MeshStruct& mesh_struct, unsigned int& max_degree, unsigned int& min_degree);
	void decideGroup(unsigned int max_color, std::vector<std::vector<unsigned int>>& unconnected_vertex_index, int* color, int size);
	void testEdgeGroup(int size, std::vector<std::vector<unsigned int>>& unconnected_vertex_index, MeshStruct& mesh_struct);

	void findMinMaxDegree(std::vector<std::vector<unsigned int>>& element_element, unsigned int& max_degree, unsigned int& min_degree);

	void getMaxPaletteSize(int& size, int* palette, unsigned int max_size, unsigned int max_array_size);

};
