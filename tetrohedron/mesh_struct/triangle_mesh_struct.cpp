#include"triangle_mesh_struct.h"
#include"../basic/enum_setting.h"

TriangleMeshStruct::TriangleMeshStruct()
{
	type = TRIANGLE;
}

void TriangleMeshStruct::setVertex()
{
	vertices.resize(vertex_position.size());
	int face_num = triangle_indices.size() / 3;
	for (int i = 0; i < face_num; ++i) {
		vertices[triangle_indices[3 * i]].face.push_back(i);
		vertices[triangle_indices[3 * i + 1]].face.push_back(i);
		vertices[triangle_indices[3 * i + 2]].face.push_back(i);
	}
}

//VERTEX_NORMAL_RENDER
void TriangleMeshStruct::getRenderVertexNormalPerThread(int thread_id)
{
	for (int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id + 1]; ++i) {
		for (int k = 0; k < vertices[i].face.size(); ++k) {
			SUM(vertex_norm_for_render[i],
				vertex_norm_for_render[i], temp_face_norm[vertices[i].face[k]]);
		}
		normalize(vertex_norm_for_render[i].data());
	}
}


//FACE_NORMAL
void TriangleMeshStruct::getFaceNormalPerThread(int thread_id)
{
	double e0[3], e2[3];
	for (int j = face_index_begin_per_thread[thread_id]; j < face_index_begin_per_thread[thread_id+1]; ++j) {
		SUB(e2, vertex_position[triangle_indices[3*j+1]], vertex_position[triangle_indices[3 * j]]);
		SUB(e0, vertex_position[triangle_indices[3* j + 2]], vertex_position[triangle_indices[3 * j]]);
		CROSS(temp_face_norm[j], e2, e0);
		memcpy(face_normal[j].data(), temp_face_norm[j].data(), 24);
		normalize(face_normal[j].data());
	}
}

void TriangleMeshStruct::getRenderNormal()
{
	thread->assignTask(this, FACE_NORMAL_RENDER);
	memset(vertex_norm_for_render[0].data(), 0, 24 * vertex_position.size());
	thread->assignTask(this, VERTEX_NORMAL_RENDER);
}

//FACE_NORMAL_RENDER
void TriangleMeshStruct::getRenderFaceNormalPerThread(int thread_id)
{
	double e0[3], e2[3];
	for (int j = face_index_begin_per_thread[thread_id]; j < face_index_begin_per_thread[thread_id + 1]; ++j) {
		SUB(e2, vertex_for_render[triangle_indices[3 * j + 1]], vertex_for_render[triangle_indices[3 * j]]);
		SUB(e0, vertex_for_render[triangle_indices[3 * j + 2]], vertex_for_render[triangle_indices[3 * j]]);
		CROSS(temp_face_norm[j].data(), e2, e0);
		memcpy(face_norm_for_render[j].data(), temp_face_norm[j].data(), 24);
		normalize(face_norm_for_render[j].data());
	}
}

void TriangleMeshStruct::getNormal()
{
	thread->assignTask(this, FACE_NORMAL);
	memset(vertex_normal[0].data(), 0, 24 * vertex_position.size());
	thread->assignTask(this, VERTEX_NORMAL);
}

//VERTEX_NORMAL
void TriangleMeshStruct::getVertexNormalPerThread(int thread_id)
{
	for (int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id+1]; ++i) {
		for (int k = 0; k < vertices[i].face.size(); ++k) {
			SUM(vertex_normal[i],
				vertex_normal[i], temp_face_norm[vertices[i].face[k]]);
		}
		normalize(vertex_normal[i].data());
	}
}


void TriangleMeshStruct::initialInfo()
{
	vertex_normal= vertex_norm_for_render;
	face_normal= face_norm_for_render;
}


void TriangleMeshStruct::setThreadIndex(int total_thread_num_) 
{
	int total_thread_num = total_thread_num_ + 1;

	vertex_index_begin_per_thread.resize(total_thread_num, 0);
	anchor_index_begin_per_thread.resize(total_thread_num, 0);
	face_index_begin_per_thread.resize(total_thread_num, 0);

	arrangeIndex(total_thread_num_, vertices.size(), vertex_index_begin_per_thread);
	arrangeIndex(total_thread_num_, anchor_vertex.size(), anchor_index_begin_per_thread);
	arrangeIndex(total_thread_num_, triangle_indices.size() / 3, face_index_begin_per_thread);

	if (!edges.empty()) {
		edge_index_begin_per_thread.resize(total_thread_num, 0);
		arrangeIndex(total_thread_num_, edges.size(), edge_index_begin_per_thread);
	}	
}

void TriangleMeshStruct::setFace()
{
	int face_num = triangle_indices.size() / 3;
	faces.resize(face_num);
	for (int i = 0; i < face_num; ++i) {
		faces[i].vertex[0] = triangle_indices[3 * i];
		faces[i].vertex[1] = triangle_indices[3 * i + 1];
		faces[i].vertex[2] = triangle_indices[3 * i + 2];
	}
}

void TriangleMeshStruct::setEdge()
{
	int face_num = triangle_indices.size() / 3;
	edges.reserve(face_num * 3 / 2);
	int edgeNo;
	for (int i = 0; i < face_num; ++i) {
		if (isEdgeExist(faces[i].vertex[0], faces[i].vertex[1], edgeNo)) {
			edges[edgeNo].face.push_back(i);
			edges[edgeNo].opposite_vertex.push_back(faces[i].vertex[2]);
			faces[i].edge.push_back(edgeNo);
		}
		else {
			addEdge(faces[i].vertex[0], faces[i].vertex[1], i, faces[i].vertex[2]);
		}
		if (isEdgeExist(faces[i].vertex[0], faces[i].vertex[2], edgeNo)) {
			edges[edgeNo].face.push_back(i);
			edges[edgeNo].opposite_vertex.push_back(faces[i].vertex[1]);
			faces[i].edge.push_back(edgeNo);
		}
		else {
			addEdge(faces[i].vertex[2], faces[i].vertex[0], i, faces[i].vertex[1]);
		}
		if (isEdgeExist(faces[i].vertex[1], faces[i].vertex[2], edgeNo)) {
			edges[edgeNo].face.push_back(i);
			edges[edgeNo].opposite_vertex.push_back(faces[i].vertex[0]);
			faces[i].edge.push_back(edgeNo);
		}
		else {
			addEdge(faces[i].vertex[1], faces[i].vertex[2], i, faces[i].vertex[0]);
		}
	}
	for (int i = 0; i < edges.size(); i++) {
		if (edges[i].face.size() == 1) {
			vertices[edges[i].vertex[0]].on_border = true;
			vertices[edges[i].vertex[1]].on_border = true;
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
		std::fill(triangle_is_used.begin(), triangle_is_used.end(), false);
		for (int k = 0; k < vertices[i].edge.size(); ++k) {
			if (edges[vertices[i].edge[k]].vertex[0] == i) {
				for (int j = 0; j < vertices[edges[vertices[i].edge[k]].vertex[1]].face.size(); ++j) {
					if (!triangle_is_used[vertices[edges[vertices[i].edge[k]].vertex[1]].face[j]]) {
						triangle_is_used[vertices[edges[vertices[i].edge[k]].vertex[1]].face[j]] = true;
						vertices[i].around_face.push_back(vertices[edges[vertices[i].edge[k]].vertex[1]].face[j]);
					}
				}
			}
			else {
				for (int j = 0; j < vertices[edges[vertices[i].edge[k]].vertex[0]].face.size(); ++j) {
					if (!triangle_is_used[vertices[edges[vertices[i].edge[k]].vertex[0]].face[j]]) {
						triangle_is_used[vertices[edges[vertices[i].edge[k]].vertex[0]].face[j]] = true;
						vertices[i].around_face.push_back(vertices[edges[vertices[i].edge[k]].vertex[0]].face[j]);
					}
				}
			}
		}
	}
}

bool TriangleMeshStruct::isEdgeExist(int v0, int v1, int& edge_index)
{

	for (int i = 0; i < vertices[v0].edge.size(); ++i) {
		if (edges[vertices[v0].edge[i]].isSame(v0, v1)) {
			edge_index = vertices[v0].edge[i];
			return true;
		}
	}
	return false;
}

void TriangleMeshStruct::addEdge(int v0, int v1, int face, int opposite_vertex)
{
	edges.push_back(Edge());
	int edge_index = edges.size() - 1;
	edges[edge_index].vertex[0] = v0;
	edges[edge_index].vertex[1] = v1;
	edges[edge_index].face.push_back(face);
	edges[edge_index].opposite_vertex.push_back(opposite_vertex);
	double tem[3];
	SUB(tem, vertex_position[v0], vertex_position[v1]);
	edges[edge_index].length = sqrt(DOT(tem, tem));
	faces[face].edge.push_back(edge_index);
	vertices[v0].edge.push_back(edge_index);
	vertices[v1].edge.push_back(edge_index);
	vertices[v0].neighbor_vertex.push_back(v1);
	vertices[v1].neighbor_vertex.push_back(v0);
}
