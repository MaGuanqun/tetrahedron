#include"tetrahedron_mesh_struct.h"


void TetrahedronMeshStruct::findSurface()
{
	triangle_indices.reserve(indices.size());
	std::map<TetrahedronFace, int> face_in_tet;
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

	vertex_on_surface.resize(vertex_position.size(),false);
	for (int i = 0; i < triangle_indices.size(); ++i) {
		for (int j = 0; j < 3; ++j) {
			vertex_on_surface[triangle_indices[i][j]] = true;
		}		
	}
}

void TetrahedronMeshStruct::buildMap(std::map<TetrahedronFace, int>& face_in_tet, int v0, int v1, int v2)
{
	TetrahedronFace A(v0, v1, v2);
	auto ret = face_in_tet.insert({ A,1 });
	if (!ret.second) {
		++ret.first->second;
	}
}

//SET_VOLUME
void TetrahedronMeshStruct::setVolume(int thread_No)
{
	for (int i = tetrahedron_index_begin_per_thread[thread_No]; i < tetrahedron_index_begin_per_thread[thread_No + 1]; ++i) {
		volume[i] = getTetrahedronVolume(vertex_position[indices[i][0]].data(), vertex_position[indices[i][1]].data(), vertex_position[indices[i][2]].data(), vertex_position[indices[i][3]].data());
	}
}

double TetrahedronMeshStruct::getTetrahedronVolume(double* v1, double* v2, double* v3, double* v4)
{
	double a[3], b[3], c[3];
	SUB(a, v2, v1);
	SUB(b, v3, v1);
	CROSS(c, a, b);
	SUB(a, v4, v1);
	return abs(DOT(c, a)) / 6.0;
}

double TetrahedronMeshStruct::setVolumeMass(double density)
{
	volume.resize(indices.size());
	mass.resize(vertex_position.size(),0.0);
	thread->assignTask(this, SET_VOLUME);
	double total_mass = 0.0;
	double tetrahedron_mass;
	for (int i = 0; i < indices.size(); ++i) {
		tetrahedron_mass = volume[i] * density * 0.25;
		mass[indices[i][0]] += tetrahedron_mass;
		mass[indices[i][1]] += tetrahedron_mass;
		mass[indices[i][2]] += tetrahedron_mass;
		mass[indices[i][3]] += tetrahedron_mass;
	}
	for (int i = 0; i < vertex_position.size(); ++i) {
		total_mass += mass[i];
	}
	return total_mass;
}

void TetrahedronMeshStruct::setThreadIndex(int total_thread_num_)
{
	int total_thread_num = total_thread_num_ + 1;

	vertex_index_begin_per_thread.resize(total_thread_num, 0);
	anchor_index_begin_per_thread.resize(total_thread_num, 0);
	face_index_begin_per_thread.resize(total_thread_num, 0);
	tetrahedron_index_begin_per_thread.resize(total_thread_num, 0);

	arrangeIndex(total_thread_num_, vertex_position.size(), vertex_index_begin_per_thread);
	arrangeIndex(total_thread_num_, anchor_vertex.size(), anchor_index_begin_per_thread);
	arrangeIndex(total_thread_num_, triangle_indices.size(), face_index_begin_per_thread);
	arrangeIndex(total_thread_num_, indices.size(), tetrahedron_index_begin_per_thread);
}

void TetrahedronMeshStruct::getRenderNormal()
{
	thread->assignTask(this, FACE_NORMAL_RENDER);
	memset(vertex_normal_for_render[0].data(), 0, 24 * vertex_position.size());
	thread->assignTask(this, VERTEX_NORMAL_RENDER);
}

//FACE_NORMAL_RENDER
void TetrahedronMeshStruct::getRenderFaceNormalPerThread(int thread_id)
{
	double e0[3], e2[3];
	for (int j = face_index_begin_per_thread[thread_id]; j < face_index_begin_per_thread[thread_id + 1]; ++j) {
		SUB(e2, vertex_for_render[triangle_indices[j][1]], vertex_for_render[triangle_indices[j][0]]);
		SUB(e0, vertex_for_render[triangle_indices[j][2]], vertex_for_render[triangle_indices[j][0]]);
		CROSS(face_normal_for_render[j].data(), e2, e0);
		normalize(face_normal_for_render[j].data());
	}
}

//VERTEX_NORMAL_RENDER
void TetrahedronMeshStruct::getRenderVertexNormalPerThread(int thread_id)
{
	for (int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id + 1]; ++i) {
		for (int k = 0; k < vertices[i].face.size(); ++k) {
			SUM(vertex_normal_for_render[i],
				vertex_normal_for_render[i], face_normal_for_render[vertices[i].face[k]]);
		}
		normalize(vertex_normal_for_render[i].data());
	}
}

void TetrahedronMeshStruct::setVertex()
{
	vertices.resize(vertex_position.size());
	for (int i = 0; i < triangle_indices.size(); ++i) {
		vertices[triangle_indices[i][0]].face.push_back(i);
		vertices[triangle_indices[i][1]].face.push_back(i);
		vertices[triangle_indices[i][2]].face.push_back(i);
	}
}