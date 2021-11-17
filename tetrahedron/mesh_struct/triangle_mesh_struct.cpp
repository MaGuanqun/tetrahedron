#include"triangle_mesh_struct.h"
#include"../basic/enum_setting.h"

TriangleMeshStruct::TriangleMeshStruct()
{
	type = TRIANGLE;
}

void TriangleMeshStruct::setVertex()
{
	vertices.resize(vertex_position.size());
	mass.resize(vertex_position.size());
	int face_num = triangle_indices.size();
	for (int i = 0; i < face_num; ++i) {
		vertices[triangle_indices[i][0]].face.push_back(i);
		vertices[triangle_indices[i][1]].face.push_back(i);
		vertices[triangle_indices[i][2]].face.push_back(i);
	}
	
}

//VERTEX_NORMAL_RENDER
void TriangleMeshStruct::getRenderVertexNormalPerThread(int thread_id)
{
	std::vector<int>* face_vertex;
	double* current_vertex_normal;
	memset(vertex_normal_for_render[vertex_index_begin_per_thread[thread_id]].data(), 
		0, 24 * (vertex_index_begin_per_thread[thread_id+1]- vertex_index_begin_per_thread[thread_id]));
	for (int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id + 1]; ++i) {
		face_vertex = &vertices[i].face;
		current_vertex_normal = vertex_normal_for_render[i].data();
		for (int k = 0; k < face_vertex->size(); ++k) {
			SUM(current_vertex_normal,
				current_vertex_normal, face_normal_for_render[(*face_vertex)[k]]);
		}
		normalize(current_vertex_normal);		
	}
	
}

//void TriangleMeshStruct::getVertexNormal()
//{
//	std::vector<int>* face_vertex;
//	double* current_vertex_normal;
//	for (int i = 0; i < vertices.size(); ++i) {
//		face_vertex = &vertices[i].face;
//		current_vertex_normal = vertex_normal[i].data();
//		for (int k = 0; k < face_vertex->size(); ++k) {
//			SUM(current_vertex_normal,
//				current_vertex_normal, face_normal[(*face_vertex)[k]]);
//		}
//		normalize(current_vertex_normal);
//	}
//}


//VERTEX_NORMAL
void TriangleMeshStruct::getVertexNormalPerThread(int thread_id)
{
	std::vector<int>* face_vertex;
	double* current_vertex_normal;
	memset(vertex_normal[vertex_index_begin_per_thread[thread_id]].data(), 0, 24 * (vertex_index_begin_per_thread[thread_id + 1]- vertex_index_begin_per_thread[thread_id]));
	for (int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id + 1]; ++i) {
		face_vertex = &vertices[i].face;
		current_vertex_normal = vertex_normal[i].data();
		for (int k = 0; k < face_vertex->size(); ++k) {
			SUM(current_vertex_normal,
				current_vertex_normal, face_normal[(*face_vertex)[k]]);
			
		}
		normalize(current_vertex_normal);
	}
}

//FACE_NORMAL_RENDER
void TriangleMeshStruct::getRenderFaceNormalPerThread(int thread_id)
{
	double e0[3], e2[3];
	double* current_face_normal;
	int* triangle_vertex;
	for (int j = face_index_begin_per_thread[thread_id]; j < face_index_begin_per_thread[thread_id + 1]; ++j) {
		triangle_vertex = triangle_indices[j].data();
		current_face_normal = face_normal_for_render[j].data();
		SUB(e2, vertex_for_render[triangle_vertex[1]], vertex_for_render[triangle_vertex[0]]);
		SUB(e0, vertex_for_render[triangle_vertex[2]], vertex_for_render[triangle_vertex[0]]);
		CROSS(current_face_normal, e2, e0);
		memcpy(ori_face_normal_for_render[j].data(), current_face_normal, 24);
		normalize(current_face_normal);
	}
}

//FACE_NORMAL
void TriangleMeshStruct::getFaceNormalPerThread(int thread_id)
{
	double e0[3], e2[3];
	double e3[3], e4[3], cross[3];
	double* current_face_normal;
	int* triangle_vertex;
	for (int j = face_index_begin_per_thread[thread_id]; j < face_index_begin_per_thread[thread_id+1]; ++j) {
		triangle_vertex = triangle_indices[j].data();
		current_face_normal = face_normal[j].data();
		SUB(e2, vertex_position[triangle_vertex[1]], vertex_position[triangle_vertex[0]]);
		SUB(e0, vertex_position[triangle_vertex[2]], vertex_position[triangle_vertex[0]]);
		CROSS(current_face_normal, e2, e0);
		memcpy(ori_face_normal[j].data(), current_face_normal, 24);
		triangle_normal_magnitude_reciprocal[j] = 1.0 / sqrt(DOT(current_face_normal, current_face_normal));
		normalize(current_face_normal);

		SUB(e3, vertex_for_render[triangle_vertex[1]], vertex_for_render[triangle_vertex[0]]);
		SUB(e4, vertex_for_render[triangle_vertex[2]], vertex_for_render[triangle_vertex[0]]);
		current_face_normal=cross_for_approx_CCD[j].data();
		CROSS(current_face_normal, e3, e0);
		CROSS(cross, e2, e4);
		SUM(current_face_normal, cross, current_face_normal);		
	}
}

void TriangleMeshStruct::getRenderNormal()
{
	thread->assignTask(this, FACE_NORMAL_RENDER);
	thread->assignTask(this, VERTEX_NORMAL_RENDER);
}



void TriangleMeshStruct::getNormal()
{
	thread->assignTask(this, FACE_NORMAL);	
	thread->assignTask(this, VERTEX_NORMAL);
}



void TriangleMeshStruct::initialInfo()
{
	vertex_normal= vertex_normal_for_render;
	face_normal= face_normal_for_render;
}


void TriangleMeshStruct::setThreadIndex(int total_thread_num_) 
{
	int total_thread_num = total_thread_num_ + 1;

	vertex_index_begin_per_thread.resize(total_thread_num, 0);
	anchor_index_begin_per_thread.resize(total_thread_num, 0);
	face_index_begin_per_thread.resize(total_thread_num, 0);

	arrangeIndex(total_thread_num_, vertices.size(), vertex_index_begin_per_thread);
	arrangeIndex(total_thread_num_, anchor_vertex.size(), anchor_index_begin_per_thread);
	arrangeIndex(total_thread_num_, triangle_indices.size(), face_index_begin_per_thread);

	if (!edges.empty()) {
		edge_index_begin_per_thread.resize(total_thread_num, 0);
		arrangeIndex(total_thread_num_, edges.size(), edge_index_begin_per_thread);
	}	
}

void TriangleMeshStruct::setFace()
{
	int face_num = triangle_indices.size();
	faces.resize(face_num);
	for (int i = 0; i < face_num; ++i) {
		faces[i].vertex[0] = triangle_indices[i][0];
		faces[i].vertex[1] = triangle_indices[i][1];
		faces[i].vertex[2] = triangle_indices[i][2];
	}
}

void TriangleMeshStruct::setEdge()
{
	int face_num = triangle_indices.size();
	edges.reserve(face_num * 3 / 2);
	int edgeNo;
	for (int i = 0; i < face_num; ++i) {
		if (isEdgeExist(triangle_indices[i][0], triangle_indices[i][1], edgeNo)) {
			edges[edgeNo].face.push_back(i);
			edges[edgeNo].opposite_vertex.push_back(triangle_indices[i][2]);
			faces[i].edge.push_back(edgeNo);
		}
		else {
			addEdge(triangle_indices[i][0], triangle_indices[i][1], i, triangle_indices[i][2]);
		}
		if (isEdgeExist(triangle_indices[i][0], triangle_indices[i][2], edgeNo)) {
			edges[edgeNo].face.push_back(i);
			edges[edgeNo].opposite_vertex.push_back(triangle_indices[i][1]);
			faces[i].edge.push_back(edgeNo);
		}
		else {
			addEdge(triangle_indices[i][2], triangle_indices[i][0], i, triangle_indices[i][1]);
		}
		if (isEdgeExist(triangle_indices[i][1], triangle_indices[i][2], edgeNo)) {
			edges[edgeNo].face.push_back(i);
			edges[edgeNo].opposite_vertex.push_back(triangle_indices[i][0]);
			faces[i].edge.push_back(edgeNo);
		}
		else {
			addEdge(triangle_indices[i][1], triangle_indices[i][2], i, triangle_indices[i][0]);
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
