#pragma once
#include"mesh_struct.h"

class GraphColor
{
public:
	
	void graphColorBending(MeshStruct& mesh_struct);
	void graphColorEdgeLength(MeshStruct& mesh_struct);

private:
	void findEdgeLengthMinMaxDegree(MeshStruct& mesh_struct, unsigned int& max_degree, unsigned int& min_degree);
	void findBendingMinMaxDegree(MeshStruct& mesh_struct, unsigned int& max_degree, unsigned int& min_degree);
	void decideGroup(unsigned int max_color, std::vector<std::vector<unsigned int>>& unconnected_vertex_index, int* color, int size);
	void testEdgeGroup(int size, std::vector<std::vector<unsigned int>>& unconnected_vertex_index, MeshStruct& mesh_struct);
};
