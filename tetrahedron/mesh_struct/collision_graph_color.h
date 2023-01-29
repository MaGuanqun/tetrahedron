#pragma once
#include<algorithm>
#include<vector>
#include<array>

class CollisionGraphColor
{
public:
	void graphColorLoopNode(std::vector<std::vector<std::vector<unsigned int>>>& element_element,
		std::vector<std::vector<unsigned int>>& element_not_connect, unsigned int* start);

	void initial(std::vector<std::vector<unsigned int>>& unconnected_index);
	int max_color_number = 30;

	void testColor(std::vector<unsigned int*>& edge_vertices, std::vector<std::array<int, 3>*>& indices, std::vector<std::vector<unsigned int>>& element_not_connect,
		std::vector<std::vector<unsigned int>*>& pair, std::vector<std::vector<std::vector<unsigned int>>>& element_element, unsigned int total_obj_num, int* vertex_num_obj);

private:

	std::vector<std::vector<int>>color;

	void graphColorLoopNodePerType(int type, std::vector<std::vector<unsigned int>>& element_element, unsigned int start);
	void decideGroup(std::vector<std::vector<unsigned int>>& unconnected_index, std::vector<std::vector<std::vector<unsigned int>>>& element_element, 
		unsigned int* start);
};

