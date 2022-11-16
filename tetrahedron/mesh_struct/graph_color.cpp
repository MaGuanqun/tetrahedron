#include"graph_color.h"


void GraphColor::findMinMaxDegree(std::vector<std::vector<unsigned int>>& element_element, unsigned int& max_degree, unsigned int& min_degree)
{
	max_degree = 0;
	min_degree = UINT_MAX;
	for (auto i = element_element.begin(); i < element_element.end(); ++i) {
		if (!i->empty()) {
			if (max_degree < i->size()) {
				max_degree = i->size();
			}
			if (min_degree > i->size()) {
				min_degree = i->size();
			}
		}
	}
}


void GraphColor::graphColor(std::vector<std::vector<unsigned int>>& element_element, std::vector<std::vector<unsigned int>>& element_not_connect)
{
	unsigned int max_degree, min_degree;
	findMinMaxDegree(element_element, max_degree, min_degree);

	std::cout << "max degree " << max_degree << std::endl;

	//initial the palette array
	unsigned int max_size = max_degree + 3;
	unsigned int max_array_size = element_element.size() * max_size;
	int* palette = new int[max_array_size];
	memset(palette, 0, 4 * max_array_size);
	unsigned int basic_palette_size = max_degree/min_degree;
	int* palette_unit = new int[basic_palette_size+2]; //here we need to record the actual num & all color number we have used
	unsigned int unit_size = basic_palette_size+2;
	for (unsigned int i = 0; i < basic_palette_size; ++i) {
		palette_unit[i + 2] = i;
	}
	palette_unit[0] = basic_palette_size;
	palette_unit[1] = basic_palette_size;

	for (unsigned int i = 0; i < max_array_size; i += max_size) {
		memcpy(palette + i, palette_unit, unit_size << 2);
	}

	unsigned int vertex_number = element_element.size();
	unsigned int ori_vertex_number = vertex_number;

	int* color = new int[vertex_number];
	std::vector<bool> U(vertex_number, true);
	std::vector<int> record_index;
	record_index.reserve(vertex_number / 2);

	bool not_exist;

	int* current_vertex_paletee_start;
	int current_rand;

	bool first_time = true;


	int indicate;

	int* palette_address;

	while (vertex_number > 0)
	{
		//tentative coloring
		if (first_time) {
			for (unsigned int i = 0; i < ori_vertex_number; ++i) {
				color[i] = palette[max_size * i + (rand() % palette[max_size * i]) + 2];
			}
		}
		else {
			for (unsigned int i = 0; i < ori_vertex_number; ++i) {
				if (U[i]) {
					indicate = -1;
					current_vertex_paletee_start = palette + max_size * i ;
					current_rand = rand() % *current_vertex_paletee_start;
					current_vertex_paletee_start ++;
					for (unsigned int j = 0; j < *current_vertex_paletee_start; ++j) {
						if (*(current_vertex_paletee_start+1+ j) > -1) {
							indicate++;
						}
						if (indicate == current_rand) {
							color[i] = *(current_vertex_paletee_start + 1 + j);
							break;
						}

					}
				}
			}
		}
		record_index.clear();
		//conflict resolution
		for (unsigned int i = 0; i < ori_vertex_number; ++i) {
			if (U[i]) {
				not_exist = true;
				for (auto j = element_element[i].begin(); j < element_element[i].end(); ++j) {
					if (color[i] == color[*j] && i < *j) {
						not_exist = false;
						break;						
					}
				}
				if (not_exist) {
					record_index.emplace_back(i);
					for (auto j = element_element[i].begin(); j < element_element[i].end(); ++j) {
						if (U[*j]) {
							if (palette[max_size * (*j) + color[i] + 2] != -1) {
								if (palette[max_size * (*j) + 1] > color[i]) {
									palette[max_size * (*j)]--;
								}
								palette[max_size * (*j) + color[i] + 2] = -1;
							}
						}
					}
				}
			}
		}
		if (!record_index.empty()) {
			first_time = false;
		}
		for (auto i = record_index.begin(); i < record_index.end(); ++i) {
			U[*i] = false;
		}
		vertex_number -= record_index.size();
		//feed the hungry

		for (unsigned int i = 0; i < ori_vertex_number; ++i) {
			if (U[i]) {
				if (palette[max_size * i] == 0) {
					palette_address = palette + max_size * i;
					(* palette_address)++;

					while (*(palette_address + (*(palette_address + 1)) + 2) == -1) {
						( *(palette_address + 1))++;
					}
					*(palette_address + (*(palette_address + 1)) + 2) = (*(palette_address + 1));
					(*(palette_address+1))++;				
				}
			}
		}
	}

	int max_palette_size;
	getMaxPaletteSize(max_palette_size, palette, max_size, max_array_size);

	delete[] palette;
	delete[] palette_unit;

	decideGroup(max_palette_size, element_not_connect, color, ori_vertex_number);
	delete[] color;
}



void GraphColor::graphColorEdgeLength(MeshStruct& mesh_struct)
{
	std::vector<std::vector<unsigned int>>edge_connect_edge;
	edge_connect_edge = mesh_struct.edge_around_edge;
	for (unsigned int i = 0; i < edge_connect_edge.size(); ++i) {
		for (unsigned int j = 0; j < edge_connect_edge[i].size(); ++j) {
			if (edge_connect_edge[i][j] == i) {
				edge_connect_edge[i][j] = *(edge_connect_edge[i].end() - 1);
				edge_connect_edge[i].pop_back();
				break;
			}
		}
	}
	graphColor(edge_connect_edge, mesh_struct.unconnected_edge_index);
	//testEdgeGroup(ori_edge_number, mesh_struct.unconnected_edge_index, mesh_struct);
}

void GraphColor::graphColorBending(MeshStruct& mesh_struct)
{

	std::vector<std::vector<unsigned int>>vertex_vertex(mesh_struct.vertex_position.size());
	for (unsigned int i = 0; i < vertex_vertex.size(); ++i) {
		vertex_vertex[i] = mesh_struct.vertices[i].around_vertex;
	}
	graphColor(vertex_vertex, mesh_struct.unconnected_vertex_index);
	//testEdgeGroup(ori_edge_number, mesh_struct.unconnected_edge_index, mesh_struct);
}

void GraphColor::decideGroup(unsigned int max_color, std::vector<std::vector<unsigned int>>& unconnected_vertex_index, int* color, int size)
{
	unconnected_vertex_index.resize(max_color);
	for (unsigned int i = 0; i < max_color; ++i) {
		unconnected_vertex_index[i].reserve(size / max_color * 2);
	}
	for (unsigned int i = 0; i < size; ++i) {
		unconnected_vertex_index[color[i]].emplace_back(i);
	}
	for (unsigned int i = 0; i < max_color; ++i) {
		unconnected_vertex_index[i].shrink_to_fit();
	}
}



void GraphColor::testBend(std::vector<std::vector<unsigned int>>& element_not_connect, MeshStruct& mesh_struct)
{
	std::vector<bool> is_used(mesh_struct.vertex_position.size(), false);

	for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
		for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); ++j) {
			if (is_used[*j]) {
				std::cout << "error bend duplicate between different group " << *j << std::endl;
			}
			is_used[*j] = true;
		}
	}
	is_used.resize(mesh_struct.vertex_position.size());
	for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
		std::fill(is_used.begin(), is_used.end(), false);
		for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); ++j) {
			if (is_used[*j]) {
				std::cout << "error bend duplicate vertex in a group " << *j << std::endl;
			}
			is_used[*j] = true;
			for (unsigned int k = 0; k < mesh_struct.vertices[*j].neighbor_vertex.size(); ++k) {
				if (is_used[mesh_struct.vertices[*j].neighbor_vertex[k]]) {
					std::cout << "error bend duplicate vertex in a group " << mesh_struct.vertices[*j].neighbor_vertex[k] << std::endl;
				}
				is_used[mesh_struct.vertices[*j].neighbor_vertex[k]] = true;
			}
		}
	}

	unsigned int num = 0;
	for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
		num += element_not_connect[i].size();
	}
	if (num != mesh_struct.vertex_position.size()) {
		std::cout << "bend lost some element " << std::endl;
	}

	std::cout << "group of bend " << element_not_connect.size() << std::endl;
	//for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
	//	//for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); ++j) {
	//	//	std::cout << *j << " ";
	//	//}
	//	std::cout<< element_not_connect[i].size() << std::endl;
	//}
}



void GraphColor::testEdge(std::vector<std::vector<unsigned int>>&element_not_connect, MeshStruct& mesh_struct, std::vector<unsigned int>& edge_vertices)
{
	std::vector<bool> is_used(mesh_struct.edge_around_edge.size(), false);

	for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
		for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); ++j) {
			if (is_used[*j]) {
				std::cout << "error edge duplicate between different group " << *j << std::endl;
			}
			is_used[*j] = true;
		}
	}
	is_used.resize(mesh_struct.vertex_position.size());
	for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
		std::fill(is_used.begin(), is_used.end(), false);
		for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); ++j) {
			for (unsigned int k = 0; k < 2; ++k) {
				if (is_used[edge_vertices[2*(* j) + k]]) {
					std::cout << "error tet duplicate vertex in a group " << *j << std::endl;
				}
				is_used[edge_vertices[2 * (*j) + k]] = true;
			}
		}
	}

	unsigned int num = 0;
	for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
		num += element_not_connect[i].size();
	}
	if (num != edge_vertices.size()/2) {
		std::cout << "edge lost some element " << std::endl;
	}

	std::cout << "group of edge " << element_not_connect.size() << std::endl;
	//for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
	//	//for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); ++j) {
	//	//	std::cout << *j << " ";
	//	//}
	//	std::cout<< element_not_connect[i].size() << std::endl;
	//}
}




void GraphColor::testTet(MeshStruct& mesh_struct, std::vector<std::array<int, 4>>& indices)
{
	std::vector<std::vector<unsigned int>>element_not_connect;
	graphColor(mesh_struct.tet_tet_index, element_not_connect);

	std::vector<bool> is_used(mesh_struct.tet_tet_index.size(), false);

	for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
		for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); ++j) {
			if (is_used[*j]) {
				std::cout << "error tet duplicate between different group " << *j << std::endl;
			}
			is_used[*j] = true;
		}
	}

	is_used.resize(mesh_struct.vertex_position.size());
	for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
		std::fill(is_used.begin(), is_used.end(), false);
		for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); ++j) {
			for(unsigned int k=0;k<4;++k){
				if (is_used[indices[*j][k]]) {
					std::cout << "error tet duplicate vertex in a group " << *j << std::endl;
				}
				is_used[indices[*j][k]] = true;
			}			
		}
	}

	unsigned int num = 0;
	for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
		num += element_not_connect[i].size();
	}
	if (num != indices.size()) {
		std::cout << "lost some element " << std::endl;
	}

	std::cout << "group of tet "<< element_not_connect.size() << std::endl;
	//for (unsigned int i = 0; i < element_not_connect.size(); ++i) {
	//	//for (auto j = element_not_connect[i].begin(); j < element_not_connect[i].end(); ++j) {
	//	//	std::cout << *j << " ";
	//	//}
	//	std::cout<< element_not_connect[i].size() << std::endl;
	//}
}



void GraphColor::testEdgeGroup(int size, std::vector<std::vector<unsigned int>>& unconnected_vertex_index, MeshStruct& mesh_struct)
{
	std::vector<bool> is_used(size, false);
	for (unsigned int i = 0; i < unconnected_vertex_index.size(); ++i) {
		for (auto j = unconnected_vertex_index[i].begin(); j < unconnected_vertex_index[i].end(); ++j) {
			if (is_used[*j]) {
				std::cout << "error " << *j << std::endl;
			}
			is_used[*j] = true;
		}
	}

	unsigned int num = 0;
	for (unsigned int i = 0; i < unconnected_vertex_index.size(); ++i) {
		num += unconnected_vertex_index[i].size();
	}
	if (num != size) {
		std::cout << "lost some element " << std::endl;
	}



	MeshStruct::Edge* edges = mesh_struct.edges.data();
	MeshStruct::Vertex* vertex = mesh_struct.vertices.data();
	unsigned int* edge_vertex = mesh_struct.edge_vertices.data();

	for (unsigned int i = 0; i < unconnected_vertex_index.size(); ++i) {
		std::fill(is_used.begin(), is_used.end(), false);
		for (auto j = unconnected_vertex_index[i].begin(); j < unconnected_vertex_index[i].end(); ++j) {
			if (is_used[*j]) {
				std::cout << "group repeat " << *j << std::endl;
			}
			is_used[*j] = true;

			/*for (auto k = vertex[edge_vertex[(*j) << 1]].edge.begin(); k < vertex[edge_vertex[(*j) << 1]].edge.end(); ++k) {
				if ((*k) != *j) {
					if (is_used[*k]) {
						std::cout << "group repeat around" << *j << " " << *k << std::endl;
					}
					is_used[*k] = true;
				}
			}
			for (auto k = vertex[edge_vertex[((*j) << 1)+1]].edge.begin(); k < vertex[edge_vertex[((*j) << 1)+1]].edge.end(); ++k) {
				if ((*k) != *j) {
					if (is_used[*k]) {
						std::cout << "group repeat around2" << *j << " " << *k << std::endl;
					}
					is_used[*k] = true;
				}
			}*/
		}



	}


}


void GraphColor::getMaxPaletteSize(int& size, int* palette, unsigned int max_size, unsigned int max_array_size)
{
	palette++;
	size = 0;
	for (unsigned int i = 0; i < max_array_size; i += max_size) {
		if (palette[i] > size) {
			size = palette[i];
		}
	}
}