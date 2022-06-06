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

}

void MeshStruct::setFace()
{
	int face_num = triangle_indices.size();
	faces.resize(face_num);
	surface_triangle_index_in_order = triangle_indices;
	//for (int i = 0; i < face_num; ++i) {
	//	faces[i].vertex[0] = triangle_indices[i][0];
	//	faces[i].vertex[1] = triangle_indices[i][1];
	//	faces[i].vertex[2] = triangle_indices[i][2];
	//}
}

void MeshStruct::setEdge()
{
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
		if (isEdgeExist(triangle_indices[i][0], triangle_indices[i][2], edgeNo)) {
			edges[edgeNo].face.push_back(i);
			edges[edgeNo].opposite_vertex.push_back(triangle_indices[i][1]);
			face_edges[3 * i + 1] = edgeNo;
			//faces[i].edge.push_back(edgeNo);
		}
		else {
			addEdge(triangle_indices[i][2], triangle_indices[i][0], i, triangle_indices[i][1], 1);
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
	}

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


//void MeshStruct::addArounVertex()
//{
//	int vertex_index;
//	std::vector<bool>point_is_used;
//	point_is_used.resize(vertices.size(), true);
//	for (int i = 0; i < vertices.size(); ++i) {
//		if (!vertices[i].face.empty()) {
//			std::fill(point_is_used.begin(), point_is_used.end(), false);
//			vertices[i].around_vertex.reserve(3 * vertices[i].neighbor_vertex.size());
//			for (int k = 0; k < vertices[i].neighbor_vertex.size(); ++k) {
//				vertex_index = vertices[i].neighbor_vertex[k];
//				for (int m = 0; m < vertices[vertex_index].neighbor_vertex.size(); ++m) {
//					if (!point_is_used[vertices[vertex_index].neighbor_vertex[m]]) {
//						vertices[i].around_vertex.push_back(vertices[vertex_index].neighbor_vertex[m]);
//					}
//				}
//			}
//		}
//	}
//}




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