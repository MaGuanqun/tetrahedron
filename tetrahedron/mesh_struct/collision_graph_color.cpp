#include"collision_graph_color.h"
#include<iostream>


void CollisionGraphColor::initial(std::vector<std::vector<unsigned int>>& unconnected_index)
{
	color.resize(5);
	unconnected_index.resize(5);
	for (int i = 0; i < 5; ++i) {
		color[i].clear();
		unconnected_index[i].clear();
	}

}

void CollisionGraphColor::graphColorLoopNode(std::vector<std::vector<std::vector<unsigned int>>>& element_element,
	std::vector<std::vector<unsigned int>>& element_not_connect, unsigned int* start)
{
	for (int i = 0; i < 5; ++i) {
		color[i].resize(element_element[i].size(), -1);
	}
	for (int i = 0; i <5; ++i) {
		graphColorLoopNodePerType(i, element_element[i], start[i]);
	}

	decideGroup(element_not_connect, element_element, start);
}



void CollisionGraphColor::graphColorLoopNodePerType(int type, std::vector<std::vector<unsigned int>>& element_element, unsigned int start)
{
	unsigned int element_number = element_element.size();
	unsigned int unique_color_number = max_color_number - 1;
	std::vector<bool>is_color_used(max_color_number, false);

	for (int i = start; i < element_number; ++i) {
		std::fill(is_color_used.begin(), is_color_used.end(), false);
		for (auto j = element_element[i].begin(); j < element_element[i].end(); j+=2) {
			if (color[*j][*(j+1)] != -1) {
				is_color_used[color[*j][*(j + 1)]] = true;
			}
		}
		for (int j = 0; j < unique_color_number; ++j) {
			if (!is_color_used[j]) {
				color[type][i] = j;
				goto color_assigned;
			}
		}
		color[type][i] = unique_color_number;
	color_assigned:;
	}
}

void CollisionGraphColor::decideGroup(std::vector<std::vector<unsigned int>>& unconnected_index, std::vector<std::vector<std::vector<unsigned int>>>& element_element,
	unsigned int* start)
{
	if (unconnected_index.size() != max_color_number) {
		unconnected_index.resize(max_color_number);
	}	
	int size = 0;
	for (int i = 0; i < 5; ++i) {
		size += element_element[i].size();
	}
	size = size / max_color_number * 2;
	for (unsigned int i = 0; i < max_color_number; ++i) {
		if (unconnected_index[i].capacity() < size) {
			unconnected_index[i].reserve(size);
		}		
	}
	for (int i = 0; i < 5; ++i) {
		for (int j = start[i]; j < element_element[i].size(); ++j) {
			unconnected_index[color[i][j]].emplace_back(i);
			unconnected_index[color[i][j]].emplace_back(j);
		}
	}
}

void CollisionGraphColor::testColor(std::vector<unsigned int*>&edge_vertices, std::vector<std::array<int,3>*>& indices, std::vector<std::vector<unsigned int>>& element_not_connect,
	std::vector<std::vector<unsigned int>*>& pair, std::vector<std::vector<std::vector<unsigned int>>>& element_element, unsigned int total_obj_num, int* vertex_num_obj)
{
	//graphColor(mesh_struct.tet_tet_index, element_not_connect);
	//graphColorLoopNode(mesh_struct.tet_tet_index,25, element_not_connect);
	//std::cout << "==out put " << std::endl;
	//for (int i = 0; i< element_element[3][11].size(); i += 2) {
	//	std::cout << element_element[3][11][i] << " " << element_element[3][11][i + 1] << std::endl;
	//}
	

	std::vector<std::vector<bool>> is_used(element_element.size());
	for (int i = 0; i < 5; ++i) {
		is_used[i].resize(element_element[i].size(), false);
	}

	for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
		for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); j+=2) {
			if (is_used[*j][*(j+1)]) {
				std::cout << "error tet duplicate between different group " << *j<<" "<<*(j+1) << std::endl;
			}
			is_used[*j][*(j + 1)] = true;
		}
	}

	is_used.resize(total_obj_num);
	for (int i = 0; i < is_used.size(); ++i) {
		is_used[i].resize(vertex_num_obj[i], false);
	}

	unsigned int* pair_;
	int* index;
	unsigned int* edge_vertex;
	for (int i = 0; i < element_not_connect.size() - 1; ++i) {
		for (int j = 0; j < is_used.size(); ++j) {
			std::fill(is_used[j].begin(), is_used[j].end(), false);
		}
		for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); j+=2) {
			switch (*j)
			{
			case 0:
				pair_ = pair[*j]->data() + ((*(j + 1)) << 2);
				if (is_used[pair_[0]][pair_[1]]) {
					std::cout << "error case 0 vertex " << *j<<" "<<*(j+1) << std::endl;
				}
				is_used[pair_[0]][pair_[1]] = true;
				index = indices[pair_[2]][pair_[3]].data();
				for (unsigned int k = 0; k < 3; ++k) {
					if (is_used[pair_[2]][index[k]]) {
						std::cout << "error case 0 triangle " << *j<<" "<<*(j+1) << std::endl;
					}
					is_used[pair_[2]][index[k]] = true;
				}
				break;
			case 1:
				pair_ = pair[*j]->data() + ((*(j + 1)) << 2);
				edge_vertex = edge_vertices[pair_[0]]+pair_[1]*2;
				for (unsigned int k = 0; k < 2; ++k) {
					if (is_used[pair_[0]][edge_vertex[k]]) {
						std::cout << "error case 1 e0 " << *j << " " << *(j + 1) << std::endl;
					}
					is_used[pair_[0]][edge_vertex[k]] = true;
				}
				edge_vertex = edge_vertices[pair_[2]] + pair_[3] * 2;
				for (unsigned int k = 0; k < 2; ++k) {
					if (is_used[pair_[2]][edge_vertex[k]]) {
						std::cout << "error case 1 e1 " << *j << " " << *(j + 1) << std::endl;
					}
					is_used[pair_[2]][edge_vertex[k]] = true;
				}
				break;
			case 2:
				pair_ = pair[*j]->data() + ((*(j + 1)) << 1);
				index = indices[pair_[0]][pair_[1]].data();
				for (unsigned int k = 0; k < 3; ++k) {
					if (is_used[pair_[0]][index[k]]) {
						std::cout << "error case 2 triangle " << *j << " " << *(j + 1) << std::endl;
					}
					is_used[pair_[0]][index[k]] = true;
				}
				break;
			case 3:
				pair_ = pair[*j]->data() + ((*(j + 1)) << 1);
				edge_vertex = edge_vertices[pair_[0]] + pair_[1] * 2;
				for (unsigned int k = 0; k < 2; ++k) {
					if (is_used[pair_[0]][edge_vertex[k]]) {
						std::cout << "error case 3 e0 " <<i<<" "<< * j << " " << *(j + 1)<<" "<<j- element_not_connect[i].begin()<<" "<< edge_vertex[0]<<" "<< edge_vertex[1] << std::endl;
					}
					is_used[pair_[0]][edge_vertex[k]] = true;
				}
				break;
			case 4:
				pair_ = pair[*j]->data() + ((*(j + 1)) << 1);
				if (is_used[pair_[0]][pair_[1]]) {
					std::cout << "error case 4 vertex " << *j << " " << *(j + 1) << std::endl;
				}
				is_used[pair_[0]][pair_[1]] = true;
				break;
			}		
		}

	}

	//pair_ = pair[3]->data() + (7<< 1);
	//edge_vertex = edge_vertices[pair_[0]] + pair_[1] * 2;
	//std::cout << edge_vertex[0] << " " << edge_vertex[1] << std::endl;


	unsigned int num = 0;
	for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
		num += element_not_connect[i].size()/2;
	}

	unsigned int actual_num=0;
	for (int i = 0; i < 5; ++i) {
		if (i < 2) {
			actual_num += pair[i]->size() / 4;
		}
		else {
			actual_num += pair[i]->size() / 2;
		}		
	}

	if (num != actual_num) {
		std::cout << "lost some element "<<num<<" " <<actual_num<< std::endl;
	}

	//std::cout << "group of tet "<< element_not_connect.size() << std::endl;
	//for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
	//	//for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); ++j) {
	//	//	std::cout << *j << " ";
	//	//}
	//	std::cout<< element_not_connect[i].size() << std::endl;
	//}
}
