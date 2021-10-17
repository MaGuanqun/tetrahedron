#include"tetrohedron_mesh_struct.h"


void TetrohedronMeshStruct::findSurface()
{
	triangle_indices.reserve(indices.size());
	std::map<TetrohedronFace, int> face_in_tet;
	int* index;
	for (int i = 0; i < indices.size(); ++i) {
		index = indices[i].data();
		buildMap(face_in_tet, index[0], index[2], index[1]);
		buildMap(face_in_tet, index[0], index[3], index[2]);
		buildMap(face_in_tet, index[0], index[1], index[3]);
		buildMap(face_in_tet, index[1], index[2], index[3]);
	}
	for (auto i = face_in_tet.begin(); i != face_in_tet.end(); ++i) {
		if (i->second == 1) {
			triangle_indices.push_back(i->first.index);
		}
	}
	triangle_indices.shrink_to_fit();
}

void TetrohedronMeshStruct::buildMap(std::map<TetrohedronFace, int>& face_in_tet, int v0, int v1, int v2)
{
	TetrohedronFace A(v0, v1, v2);
	auto ret = face_in_tet.insert({ A,1 });
	if (!ret.second) {
		++ret.first->second;
	}
}

void TetrohedronMeshStruct::setThreadIndex(int total_thread_num_)
{
	int total_thread_num = total_thread_num_ + 1;

	vertex_index_begin_per_thread.resize(total_thread_num, 0);
	anchor_index_begin_per_thread.resize(total_thread_num, 0);
	face_index_begin_per_thread.resize(total_thread_num, 0);
	tetrohedron_index_begin_per_thread.resize(total_thread_num, 0);

	arrangeIndex(total_thread_num_, vertex_position.size(), vertex_index_begin_per_thread);
	arrangeIndex(total_thread_num_, anchor_vertex.size(), anchor_index_begin_per_thread);
	arrangeIndex(total_thread_num_, triangle_indices.size(), face_index_begin_per_thread);
	arrangeIndex(total_thread_num_, indices.size(), tetrohedron_index_begin_per_thread);
}

void TetrohedronMeshStruct::getRenderNormal()
{
	thread->assignTask(this, FACE_NORMAL_RENDER);
	memset(vertex_norm_for_render[0].data(), 0, 24 * vertex_position.size());
	thread->assignTask(this, VERTEX_NORMAL_RENDER);
}

//FACE_NORMAL_RENDER
void TetrohedronMeshStruct::getRenderFaceNormalPerThread(int thread_id)
{
	double e0[3], e2[3];
	for (int j = face_index_begin_per_thread[thread_id]; j < face_index_begin_per_thread[thread_id + 1]; ++j) {
		SUB(e2, vertex_for_render[triangle_indices[j][1]], vertex_for_render[triangle_indices[j][0]]);
		SUB(e0, vertex_for_render[triangle_indices[j][2]], vertex_for_render[triangle_indices[j][0]]);
		CROSS(temp_face_norm[j].data(), e2, e0);
		memcpy(face_norm_for_render[j].data(), temp_face_norm[j].data(), 24);
		normalize(face_norm_for_render[j].data());
	}
}

//VERTEX_NORMAL_RENDER
void TetrohedronMeshStruct::getRenderVertexNormalPerThread(int thread_id)
{
	for (int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id + 1]; ++i) {
		for (int k = 0; k < vertices[i].face.size(); ++k) {
			SUM(vertex_norm_for_render[i],
				vertex_norm_for_render[i], temp_face_norm[vertices[i].face[k]]);
		}
		normalize(vertex_norm_for_render[i].data());
	}
}

void TetrohedronMeshStruct::setVertex()
{
	vertices.resize(vertex_position.size());
	for (int i = 0; i < triangle_indices.size(); ++i) {
		vertices[triangle_indices[i][0]].face.push_back(i);
		vertices[triangle_indices[i][1]].face.push_back(i);
		vertices[triangle_indices[i][2]].face.push_back(i);
	}
}