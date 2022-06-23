#include"graph_color.h"

void GraphColor::findMinMaxDegree(MeshStruct& mesh_struct, unsigned int& max_degree, unsigned int& min_degree)
{
	unsigned int max_degree = 0;
	unsigned int min_degree = UINT_MAX;
	MeshStruct::Vertex* vertex=mesh_struct.vertices.data();
	for (unsigned int i = 0; i < mesh_struct.vertex_position.size(); ++i) {
		if (max_degree < vertex[i].around_vertex.size()) {
			max_degree = vertex[i].around_vertex.size();
		}
		if (min_degree > vertex[i].around_vertex.size()) {
			min_degree = vertex[i].around_vertex.size();
		}
	}
}


void GraphColor::graphColorEdgeLength(MeshStruct& mesh_struct, unsigned int max_degree, unsigned int min_degree)
{
	//initial the palette array
	unsigned int max_size = max_degree + 2;
	unsigned int max_array_size = mesh_struct.edges.size() * (max_degree + 2);
	int* palette = new int[max_array_size];
	memset(palette, -1, 4 * max_array_size);
	unsigned int basic_palette_size = max_degree / min_degree + 1;
	int* palette_unit = new int[basic_palette_size + 1];
	unsigned int unit_size = basic_palette_size + 1;
	for (unsigned int i = 0; i < basic_palette_size; ++i) {
		palette_unit[i + 1] = i;
	}
	palette_unit[0] = basic_palette_size;

	for (unsigned int i = 0; i < max_array_size; i += max_size) {
		memcpy(palette + i, palette_unit, unit_size << 2);
	}

	unsigned int edge_number = mesh_struct.edges.size();
	unsigned int ori_edge_number = edge_number;

	int* color = new int[edge_number];
	std::vector<bool> U(edge_number, true);
	std::vector<int> record_index;
	record_index.reserve(edge_number / 2);

	MeshStruct::Edge* edges = mesh_struct.edges.data();


	MeshStruct::Vertex* vertex = mesh_struct.vertices.data();
	
	unsigned int* edge_vertex = mesh_struct.edge_vertices.data();


	bool not_exist;

	int current_vertex_paletee_start = 0;
	int current_rand;

	bool first_time = true;
	bool need_to_add_paletee_unit;


	int indicate;
	while (edge_number > 0)
	{
		//tentative coloring
		if (first_time) {
			for (unsigned int i = 0; i < ori_edge_number; ++i) {
				color[i] = palette[max_size * i + (rand() % palette[max_size * i])];
			}
		}
		else {
			for (unsigned int i = 0; i < ori_edge_number; ++i) {
				if (U[i]) {
					indicate = -1;
					current_vertex_paletee_start = max_size * i;
					current_rand = rand() % palette[current_vertex_paletee_start];
					for (unsigned int j = 0; j < basic_palette_size; ++j) {
						if (palette[current_vertex_paletee_start + j + 1] > -1) {
							indicate++;
						}
						if (indicate == current_rand) {
							color[i] = palette[current_vertex_paletee_start + j + 1];
							break;
						}

					}
				}
			}
		}
		record_index.clear();
		//conflict resolution
		unsigned int i_2 = 0;
		for (unsigned int i = 0; i < ori_edge_number; ++i) {
			if (U[i]) {
				not_exist = true;
				for (auto j = vertex[edge_vertex[i_2]].edge.begin(); j < vertex[edge_vertex[i_2]].edge.end(); ++j) {
					if (color[i] == color[*j]) {
						not_exist = false;
						break;
					}
				}
				if (not_exist) {
					for (auto j = vertex[edge_vertex[i_2+1]].edge.begin(); j < vertex[edge_vertex[i_2+1]].edge.end(); ++j) {
						if (color[i] == color[*j]) {
							not_exist = false;
							break;
						}
					}
				}
				if (not_exist) {
					record_index.emplace_back(i);
					for (auto j = vertex[edge_vertex[i_2]].edge.begin(); j < vertex[edge_vertex[i_2]].edge.end(); ++j) {
						if (palette[max_size * (*j) + color[i] + 1] != -1) {
							palette[max_size * (*j)]--;
							palette[max_size * (*j) + color[i] + 1] = -1;
						}
					}
					for (auto j = vertex[edge_vertex[i_2 + 1]].edge.begin(); j < vertex[edge_vertex[i_2 + 1]].edge.end(); ++j) {
						if (palette[max_size * (*j) + color[i] + 1] != -1) {
							palette[max_size * (*j)]--;
							palette[max_size * (*j) + color[i] + 1] = -1;
						}
					}
				}
			}
			i_2 += 2;
		}
		if (!record_index.empty()) {
			first_time = false;
		}
		for (auto i = record_index.begin(); i < record_index.end(); ++i) {
			U[*i] = false;
		}
		edge_number -= record_index.size();
		//feed the hungry
		need_to_add_paletee_unit = false;
		for (unsigned int i = 0; i < ori_edge_number; ++i) {
			if (U[i]) {
				if (palette[max_size * i] == 0) {
					palette[max_size * i]++;
					palette[max_size * i + basic_palette_size + 1] = basic_palette_size;
					need_to_add_paletee_unit = true;
				}
			}
		}
		if (need_to_add_paletee_unit) {
			basic_palette_size++;
		}
	}
	delete[] palette;
	delete[] palette_unit;

	decideGroup(basic_palette_size, mesh_struct.unconnected_edge_index, color, ori_edge_number);
	delete[] color;
}


void GraphColor::graphColorBending(MeshStruct& mesh_struct, unsigned int max_degree, unsigned int min_degree)
{

	//initial the palette array
	unsigned int max_size = max_degree + 2;
	unsigned int max_array_size = mesh_struct.vertex_position.size() * (max_degree + 2);
	int* palette = new int[max_array_size];
	memset(palette, -1, 4 * max_array_size);
	unsigned int basic_palette_size = max_degree / min_degree + 1;
	int* palette_unit = new int[basic_palette_size +1];
	unsigned int unit_size = basic_palette_size + 1;
	for (unsigned int i = 0; i < basic_palette_size; ++i) {
		palette_unit[i + 1] = i;
	}
	palette_unit[0] = basic_palette_size;

	for (unsigned int i = 0; i < max_array_size; i += max_size) {
		memcpy(palette + i, palette_unit,  unit_size<<2);
	}

	unsigned int vertex_number = mesh_struct.vertex_position.size();
	unsigned int ori_vertex_number= vertex_number;

	int* color = new int[vertex_number];
	std::vector<bool> U(vertex_number,true);
	std::vector<int> record_index;
	record_index.reserve(vertex_number / 2);

	MeshStruct::Vertex* vertices = mesh_struct.vertices.data();
	
	bool not_exist;

	int current_vertex_paletee_start=0;
	int current_rand;

	bool first_time = true;
	bool need_to_add_paletee_unit;


	int indicate;
	while (vertex_number>0)
	{
		//tentative coloring
		if (first_time) {
			for (unsigned int i = 0; i < ori_vertex_number; ++i) {
				color[i] = palette[max_size * i + (rand() % palette[max_size * i])];
			}
		}
		else {
			for (unsigned int i = 0; i < ori_vertex_number; ++i) {
				if (U[i]) {
					indicate = -1;
					current_vertex_paletee_start = max_size * i;
					current_rand = rand() % palette[current_vertex_paletee_start];
					for (unsigned int j = 0; j < basic_palette_size; ++j) {
						if (palette[current_vertex_paletee_start + j + 1] > -1) {
							indicate++;
						}
						if (indicate ==current_rand) {
							color[i] = palette[current_vertex_paletee_start + j + 1];
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
				not_exist =true;
				for (unsigned int j = 0; j < vertices[i].around_vertex.size(); ++j) {
					if (color[i] == color[vertices[i].around_vertex[j]]) {
						not_exist = false;
						break;
					}
				}
				if (not_exist) {
					record_index.emplace_back(i);
					for (auto j = vertices[i].around_vertex.begin(); j < vertices[i].around_vertex.end(); ++j) {
						if (palette[max_size * (*j) + color[i] + 1] != -1) {
							palette[max_size * (*j)]--;
							palette[max_size * (*j) + color[i] + 1] = -1;
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
		need_to_add_paletee_unit = false;
		for (unsigned int i = 0; i < ori_vertex_number; ++i) {
			if (U[i]) {
				if (palette[max_size * i] == 0) {
					palette[max_size * i]++;
					palette[max_size * i + basic_palette_size + 1] = basic_palette_size;
					need_to_add_paletee_unit = true;
				}			
			}
		}
		if (need_to_add_paletee_unit) {
			basic_palette_size++;
		}	
	}
	delete[] palette;
	delete[] palette_unit;

	decideGroup(basic_palette_size, mesh_struct.unconnected_vertex_index, color, ori_vertex_number);
	delete[] color;
}

void GraphColor::decideGroup(unsigned int max_color, std::vector<std::vector<unsigned int>>& unconnected_vertex_index, int* color, int size)
{
	unconnected_vertex_index.resize(max_color);
	for (unsigned int i = 0; i < max_color; ++i) {
		unconnected_vertex_index.reserve(size / max_color * 2);
	}
	for (unsigned int i = 0; i < size; ++i) {
		unconnected_vertex_index[color[i]].emplace_back(i);
	}
}

