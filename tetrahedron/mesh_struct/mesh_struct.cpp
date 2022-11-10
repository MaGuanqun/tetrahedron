#include"mesh_struct.h"
void MeshStruct::initialNormalSize()
{
	face_normal_for_render.resize(triangle_indices.size(), { 0.0,0.0,0.0 });
	vertex_normal_for_render.resize(vertex_position.size());
	triangle_normal_magnitude_reciprocal.resize(triangle_indices.size());
	ori_face_normal_for_render.resize(triangle_indices.size(), { 0.0,0.0,0.0 });
	cross_for_approx_CCD.resize(triangle_indices.size(), { 0.0,0.0,0.0 });
	face_normal = face_normal_for_render;
	ori_face_normal = ori_face_normal_for_render;
	vertex_normal = vertex_normal_for_render;

	f_face_normal_for_render.resize(triangle_indices.size());
	f_face_normal.resize(triangle_indices.size());
	f_cross_for_approx_CCD.resize(triangle_indices.size());
}


void MeshStruct::setAnchorPosition()
{
	anchor_position.resize(anchor_vertex.size());
	for (int i = 0; i < anchor_vertex.size(); ++i) {
		anchor_position[i] = vertex_position[anchor_vertex[i]];
	}
}

void MeshStruct::resetMassInv()
{
	for (unsigned int i = 0; i < mass_inv.size(); ++i) {
		mass_inv[i] = 1.0 / mass[i];
	}
	for (unsigned int i = 0; i < anchor_vertex.size(); ++i) {
		mass_inv[anchor_vertex[i]] = 0.0;
	}
	std::fill(is_vertex_fixed.begin(), is_vertex_fixed.end(), false);
	for (unsigned int i = 0; i < anchor_vertex.size(); ++i) {
		is_vertex_fixed[anchor_vertex[i]] = true;
	}
}

void MeshStruct::updateAnchorPerThread(int total_thread_num)
{
	arrangeIndex(total_thread_num, anchor_vertex.size(), anchor_index_begin_per_thread.data());
}

void MeshStruct::setVertex()
{
	vertices.resize(vertex_position.size());
	mass.resize(vertex_position.size(), 0.0);
	mass_inv.resize(vertex_position.size(), 0.0);
	int face_num = triangle_indices.size();
	for (unsigned int i = 0; i < face_num; ++i) {
		vertices[triangle_indices[i][0]].face.push_back(i);
		vertices[triangle_indices[i][1]].face.push_back(i);
		vertices[triangle_indices[i][2]].face.push_back(i);
	}
	is_vertex_fixed.resize(vertex_position.size(), false);
}

//SORT_TRIANGLE_AROUND_VERTEX_EDGE
void MeshStruct::sortTriangleAroundVertexEdge(int thread_id)
{
	auto k = vertices.begin() + vertex_index_begin_per_thread[thread_id + 1];
	for (auto i = vertices.begin() + vertex_index_begin_per_thread[thread_id]; i < k; ++i) {
		std::sort(i->face.begin(), i->face.end());
		std::sort(i->edge.begin(), i->edge.end());
	}


}


void MeshStruct::setFace()
{
	int face_num = triangle_indices.size();
	faces.resize(face_num);
	surface_triangle_index_in_order = triangle_indices;
	face_around_face.resize(triangle_indices.size());
	edge_around_face.resize(triangle_indices.size());
	//for (int i = 0; i < face_num; ++i) {
	//	faces[i].vertex[0] = triangle_indices[i][0];
	//	faces[i].vertex[1] = triangle_indices[i][1];
	//	faces[i].vertex[2] = triangle_indices[i][2];
	//}
}

void MeshStruct::setEdgeForSpring()
{
	if (vertex_position.size() == 2) {
		edges.resize(1);
		edge_vertices.resize(2);
		edge_vertices[0] = 0;
		edge_vertices[1] = 1;
		double v[3];
		SUB(v, vertex_position[0], vertex_position[1]);
		edge_length.emplace_back(sqrt(DOT(v,  v)));
		vertices[0].edge.emplace_back(0);
		vertices[1].edge.emplace_back(0);
		addNeighborVertex();
		vertices[0].neighbor_vertex.emplace_back(1);
		vertices[1].neighbor_vertex.emplace_back(0);
		vertices[0].around_vertex.emplace_back(1);
		vertices[1].around_vertex.emplace_back(0);
	}
}


//SORT_TRIANGLE_EDGE_AROUND_TRIANGLE_EDGE
void MeshStruct::setFaceEdgeAroundFace(int thread_id)
{
	std::vector<unsigned int> commen_triangle;
	commen_triangle.reserve(4); 
	unsigned int* edge_index;
	int* vertex_index;

	std::vector<unsigned int>* vertex_connnect_edge;

	for (auto i = face_index_begin_per_thread[thread_id]; i < face_index_begin_per_thread[thread_id + 1]; ++i) {
		//face
		vertex_index = triangle_indices[i].data();
		edge_index = face_edges.data() + 3 * i;
		commen_triangle.emplace_back(i);

		for (unsigned int j = 0; j < 3; ++j) {
			for (auto k = edges[edge_index[j]].face.begin(); k < edges[edge_index[j]].face.end(); ++k) {
				if (*k != i) {
					commen_triangle.emplace_back(*k);
				}
			}
		}
		face_around_face[i].reserve(6);
		for (int j = 0; j < 3; ++j) {
			for (auto k = vertices[vertex_index[j]].face.begin(); k < vertices[vertex_index[j]].face.end(); ++k) {
				if (!isCommonUsed(*k, &commen_triangle)) {
					face_around_face[i].emplace_back(*k);
				}
			}
		}

		face_around_face[i].insert(face_around_face[i].end(), commen_triangle.begin() + 1, commen_triangle.end());

		std::sort(face_around_face[i].begin(), face_around_face[i].end());

		//edge
		edge_around_face[i].reserve(8);
		for (unsigned int j = 0; j < 3; ++j) {
			vertex_connnect_edge = &vertices[vertex_index[j]].edge;
			for (auto k = vertex_connnect_edge->begin(); k < vertex_connnect_edge->end(); ++k) {
				if ((*k) != edge_index[0] && (*k) != edge_index[1] && (*k) != edge_index[2]) {
					edge_around_face[i].emplace_back(*k);
				}
			}
		}
		edge_around_face[i].emplace_back(edge_index[0]);
		edge_around_face[i].emplace_back(edge_index[1]);
		edge_around_face[i].emplace_back(edge_index[2]);
		std::sort(edge_around_face[i].begin(), edge_around_face[i].end());		
	}

	unsigned int* edge_vertex_index;
	for (auto i = edge_index_begin_per_thread[thread_id]; i < edge_index_begin_per_thread[thread_id + 1]; ++i) {
		//face
		edge_vertex_index = edge_vertices.data()+i*2;
		for (auto k = edges[i].face.begin(); k < edges[i].face.end(); ++k) {
				commen_triangle.emplace_back(*k);			
		}		
		face_around_edge[i].reserve(6);
		for (int j = 0; j < 2; ++j) {
			for (auto k = vertices[edge_vertex_index[j]].face.begin(); k < vertices[edge_vertex_index[j]].face.end(); ++k) {
				if (!isCommonUsed(*k, &commen_triangle)) {
					face_around_edge[i].emplace_back(*k);
				}
			}
		}
		face_around_edge[i].insert(face_around_edge[i].end(), commen_triangle.begin(), commen_triangle.end());
		std::sort(face_around_edge[i].begin(), face_around_edge[i].end());

		//edge

		edge_around_edge[i].reserve(8);
		for (unsigned int j = 0; j < 2; ++j) {
			vertex_connnect_edge = &vertices[edge_vertex_index[j]].edge;
			for (auto k = vertex_connnect_edge->begin(); k < vertex_connnect_edge->end(); ++k) {
				if ((*k) !=i) {
					edge_around_edge[i].emplace_back(*k);
				}
			}
		}
		edge_around_edge[i].emplace_back(i);
		std::sort(edge_around_edge[i].begin(), edge_around_edge[i].end());
	}


}



bool MeshStruct::isCommonUsed(unsigned int tet_num, std::vector<unsigned int>* tet_index)
{
	for (auto i = tet_index->begin(); i < tet_index->end(); ++i) {
		if (*i == tet_num) {
			return true;
		}
	}
	return false;
}


void MeshStruct::setEdge()
{
	if (triangle_indices.empty()) {
		return;
	}
	unsigned int face_num = triangle_indices.size();
	edges.reserve(face_num * 3 / 2);
	edge_length.reserve(face_num * 3 / 2);
	face_edges.resize(face_num * 3);
	edge_vertices.reserve(face_num * 3);
	unsigned int edgeNo;
	for (unsigned int i = 0; i < face_num; ++i) {
		if (isEdgeExist(triangle_indices[i][0], triangle_indices[i][1], edgeNo)) {
			edges[edgeNo].face.push_back(i);
			edges[edgeNo].opposite_vertex.push_back(triangle_indices[i][2]);
			face_edges[3 * i] = edgeNo;
			//faces[i].edge.push_back(edgeNo);
		}
		else {
			addEdge(triangle_indices[i][0], triangle_indices[i][1], i, triangle_indices[i][2], 0);
		}
		if (isEdgeExist(triangle_indices[i][1], triangle_indices[i][2], edgeNo)) {
			edges[edgeNo].face.push_back(i);
			edges[edgeNo].opposite_vertex.push_back(triangle_indices[i][0]);
			//faces[i].edge.push_back(edgeNo);
			face_edges[3 * i + 2] = edgeNo;
		}
		else {
			addEdge(triangle_indices[i][1], triangle_indices[i][2], i, triangle_indices[i][0], 2);
		}
		if (isEdgeExist(triangle_indices[i][2], triangle_indices[i][0], edgeNo)) {
			edges[edgeNo].face.push_back(i);
			edges[edgeNo].opposite_vertex.push_back(triangle_indices[i][1]);
			face_edges[3 * i + 1] = edgeNo;
			//faces[i].edge.push_back(edgeNo);
		}
		else {
			addEdge(triangle_indices[i][2], triangle_indices[i][0], i, triangle_indices[i][1], 1);
		}
	}

	//unsigned int index[112] = { 0,1,0,5,0,6,1,2,1,6,2,3,2,6,2,7,2,8,3,4,3,8,4,8,4,9,5,6,5,10,6,7,6,10,6,11,6,12,7,8,7,12,8,9,8,12,8,13,8,14,9,14,10,11,10,15,10,16,11,12,11,16,12,13,12,16,12,17,12,18,13,14,13,18,14,18,14,19,15,16,15,20,16,17,16,20,16,21,16,22,17,18,17,22,18,19,18,22,18,23,18,24,19,24,20,21,21,22,22,23,23,24 };


	//for (unsigned int i = 0; i < 112; i+=2) {
	//	addEdge(index[i], index[i + 1], 3, -1, 1);
	//}

	addNeighborVertex();
	edge_vertices.shrink_to_fit();
	edges.shrink_to_fit();
	edge_length.shrink_to_fit();

	for (int i = 0; i < edges.size(); i++) {
		if (edges[i].face.size() == 1) {
			vertices[edge_vertices[i << 1]].on_border = true;
			vertices[edge_vertices[(i << 1) + 1]].on_border = true;
			edges[i].is_border = true;
		}
	}
	//for (int i = 0; i < edges.size(); ++i) {
	//	if (vertices[edges[i].vertex[0]].on_border || vertices[edges[i].vertex[1]].on_border) {
	//		edges[i].near_border = true;
	//	}
	//}
	std::vector<bool>triangle_is_used(faces.size(), false);
	for (int i = 0; i < vertices.size(); i++) {
		//if (vertices[i].on_border) {
		//	for (int j = 0; j < vertices[i].face.size(); j++) {
		//		faces[vertices[i].face[j]].on_border = true;
		//	}
		//}
		if (!vertices[i].edge.empty()) {
			std::fill(triangle_is_used.begin(), triangle_is_used.end(), false);
			for (int k = 0; k < vertices[i].edge.size(); ++k) {
				if (edge_vertices[(vertices[i].edge[k] << 1)] == i) {
					for (int j = 0; j < vertices[edge_vertices[(vertices[i].edge[k] << 1) + 1]].face.size(); ++j) {
						if (!triangle_is_used[vertices[edge_vertices[(vertices[i].edge[k] << 1) + 1]].face[j]]) {
							triangle_is_used[vertices[edge_vertices[(vertices[i].edge[k] << 1) + 1]].face[j]] = true;
							vertices[i].around_face.push_back(vertices[edge_vertices[(vertices[i].edge[k] << 1) + 1]].face[j]);
						}
					}
				}
				else {
					for (int j = 0; j < vertices[edge_vertices[vertices[i].edge[k] << 1]].face.size(); ++j) {
						if (!triangle_is_used[vertices[edge_vertices[vertices[i].edge[k] << 1]].face[j]]) {
							triangle_is_used[vertices[edge_vertices[vertices[i].edge[k] << 1]].face[j]] = true;
							vertices[i].around_face.push_back(vertices[edge_vertices[vertices[i].edge[k] << 1]].face[j]);
						}
					}
				}
			}
		}
	}

	face_around_edge.resize(edge_length.size());
	edge_around_edge.resize(edge_length.size());
	tet_around_edge.resize(edge_length.size());

}




void MeshStruct::addNeighborVertex()
{
	for (int i = 0; i < vertices.size(); ++i) {
		if (!vertices[i].edge.empty()) {
			vertices[i].neighbor_vertex.reserve(vertices[i].edge.size());
			for (int j = 0; j < vertices[i].edge.size(); ++j)
			{
				if (i == edge_vertices[vertices[i].edge[j] << 1]) {
					vertices[i].neighbor_vertex.push_back(edge_vertices[(vertices[i].edge[j] << 1) + 1]);
				}
				else {
					vertices[i].neighbor_vertex.push_back(edge_vertices[vertices[i].edge[j] << 1]);
				}
			}
		}
	}
}

bool MeshStruct::isEdgeExist(unsigned int v0, unsigned int v1, unsigned int& edge_index)
{
	unsigned int* edge_index_address;
	for (int i = 0; i < vertices[v0].edge.size(); ++i) {
		edge_index_address = edge_vertices.data() + (vertices[v0].edge[i] << 1);
		if ((edge_index_address[0] == v0) && (edge_index_address[1] == v1)) {
			edge_index = vertices[v0].edge[i];
			return true;
		}
		if ((edge_index_address[0] == v1) && (edge_index_address[1] == v0)) {
			edge_index = vertices[v0].edge[i];
			return true;
		}
	}
	return false;
}

void MeshStruct::addEdge(int v0, int v1, int face, int opposite_vertex, int edge_index_indicator)
{
	edges.emplace_back(Edge());
	unsigned int edge_index = edges.size() - 1;
	edge_vertices.emplace_back(v0);
	edge_vertices.emplace_back(v1);
	//edges[edge_index].vertex[0] = v0;
	//edges[edge_index].vertex[1] = v1;
	edges[edge_index].face.emplace_back(face);
	edges[edge_index].opposite_vertex.emplace_back(opposite_vertex);
	double tem[3];
	SUB(tem, vertex_position[v0], vertex_position[v1]);
	edge_length.emplace_back(sqrt(DOT(tem, tem)));
	face_edges[3 * face + edge_index_indicator] = edge_index;
	//faces[face].edge.push_back(edge_index);
	vertices[v0].edge.emplace_back(edge_index);
	vertices[v1].edge.emplace_back(edge_index);
	//vertices[v0].neighbor_vertex.push_back(v1);
	//vertices[v1].neighbor_vertex.push_back(v0);
}


void MeshStruct::addArounVertex()
{
	std::vector<bool>point_is_used(vertices.size(), false);
	for (int i = 0; i < vertices.size(); ++i) {
		if (!vertices[i].face.empty()) {
			std::fill(point_is_used.begin(), point_is_used.end(), false);
			point_is_used[i] = true;
			vertices[i].around_vertex.reserve(3 * vertices[i].neighbor_vertex.size());
			vertices[i].around_vertex.insert(vertices[i].around_vertex.end(),vertices[i].neighbor_vertex.begin(), vertices[i].neighbor_vertex.end());
			for (auto j = vertices[i].neighbor_vertex.begin(); j < vertices[i].neighbor_vertex.end(); ++j) {
				point_is_used[*j] = true;
			}
			for (auto k = vertices[i].neighbor_vertex.begin(); k < vertices[i].neighbor_vertex.end(); ++k) {
				for (auto m = vertices[*k].neighbor_vertex.begin(); m < vertices[*k].neighbor_vertex.end(); ++m) {
					if (!point_is_used[*m]) {
						vertices[i].around_vertex.emplace_back(*m);
						point_is_used[*m] = true;
					}
				}
			}
		}
	}
}




//FACE_NORMAL_RENDER
void MeshStruct::getRenderFaceNormalPerThread(int thread_id)
{
	double e0[3], e2[3];
	for (int j = face_index_begin_per_thread[thread_id]; j < face_index_begin_per_thread[thread_id + 1]; ++j) {
		SUB(e2, vertex_for_render[triangle_indices[j][1]], vertex_for_render[triangle_indices[j][0]]);
		SUB(e0, vertex_for_render[triangle_indices[j][2]], vertex_for_render[triangle_indices[j][0]]);
		CROSS(face_normal_for_render[j].data(), e2, e0);
		memcpy(ori_face_normal_for_render[j].data(), face_normal_for_render[j].data(), 24);
		normalize(face_normal_for_render[j].data());
	}
}

void MeshStruct::getRenderFaceNormal()
{
	double e0[3], e2[3];
	for (int j = 0; j < face_normal.size(); ++j) {
		SUB(e2, vertex_for_render[triangle_indices[j][1]], vertex_for_render[triangle_indices[j][0]]);
		SUB(e0, vertex_for_render[triangle_indices[j][2]], vertex_for_render[triangle_indices[j][0]]);
		CROSS(face_normal_for_render[j].data(), e2, e0);
		memcpy(ori_face_normal_for_render[j].data(), face_normal_for_render[j].data(), 24);
		normalize(face_normal_for_render[j].data());
	}
}


void MeshStruct::getFaceNormal()
{
	double e0[3], e2[3];
	double e3[3], e4[3];
	double* current_face_normal;
	int* triangle_vertex;
	double f_cross[3];
	for (int j = 0; j < face_normal.size(); ++j) {
		triangle_vertex = triangle_indices[j].data();
		current_face_normal = face_normal[j].data();
		SUB(e2, vertex_position[triangle_vertex[1]], vertex_position[triangle_vertex[0]]);
		SUB(e0, vertex_position[triangle_vertex[2]], vertex_position[triangle_vertex[0]]);
		CROSS(current_face_normal, e2, e0);
		memcpy(ori_face_normal[j].data(), current_face_normal, 24);
		triangle_normal_magnitude_reciprocal[j] = 1.0 / sqrt(DOT(current_face_normal, current_face_normal));
		MULTI_(current_face_normal, triangle_normal_magnitude_reciprocal[j]);
		//normalize(current_face_normal);

		SUB(e3, vertex_for_render[triangle_vertex[1]], vertex_for_render[triangle_vertex[0]]);
		SUB(e4, vertex_for_render[triangle_vertex[2]], vertex_for_render[triangle_vertex[0]]);

		current_face_normal = cross_for_approx_CCD[j].data();

		CROSS(current_face_normal, e3, e0);
		CROSS(f_cross, e2, e4);
		SUM_(current_face_normal, f_cross);
	}
}



//FACE_NORMAL
void MeshStruct::getFaceNormalPerThread(int thread_id)
{
	double e0[3], e2[3];
	double e3[3], e4[3];
	double* current_face_normal;
	int* triangle_vertex;
	double f_cross[3];
	for (int j = face_index_begin_per_thread[thread_id]; j < face_index_begin_per_thread[thread_id + 1]; ++j) {
		triangle_vertex = triangle_indices[j].data();
		current_face_normal = face_normal[j].data();
		SUB(e2, vertex_position[triangle_vertex[1]], vertex_position[triangle_vertex[0]]);
		SUB(e0, vertex_position[triangle_vertex[2]], vertex_position[triangle_vertex[0]]);
		CROSS(current_face_normal, e2, e0);
		memcpy(ori_face_normal[j].data(), current_face_normal, 24);
		triangle_normal_magnitude_reciprocal[j] = 1.0 / sqrt(DOT(current_face_normal, current_face_normal));
		MULTI_(current_face_normal, triangle_normal_magnitude_reciprocal[j]);
		//normalize(current_face_normal);

		SUB(e3, vertex_for_render[triangle_vertex[1]], vertex_for_render[triangle_vertex[0]]);
		SUB(e4, vertex_for_render[triangle_vertex[2]], vertex_for_render[triangle_vertex[0]]);

		current_face_normal = cross_for_approx_CCD[j].data();

		CROSS(current_face_normal, e3, e0);
		CROSS(f_cross, e2, e4);
		SUM_(current_face_normal, f_cross);
	}
}