#include"tetrahedron_mesh_struct.h"


void TetrahedronMeshStruct::setTetEdges()
{
	std::vector<std::vector<unsigned int>> edge_vertex(vertex_position.size());
	
	for (unsigned int i = 0; i < indices.size(); ++i) {
		addTetEdges(indices[i][0], indices[i][1], edge_vertex);
		addTetEdges(indices[i][0], indices[i][2], edge_vertex);
		addTetEdges(indices[i][0], indices[i][3], edge_vertex);
		addTetEdges(indices[i][1], indices[i][2], edge_vertex);
		addTetEdges(indices[i][1], indices[i][3], edge_vertex);
		addTetEdges(indices[i][2], indices[i][3], edge_vertex);
	}
	tet_edge_vertices.reserve(20 * vertex_position.size());
	for (unsigned int i = 0; i < vertex_position.size(); ++i) {
		for (unsigned int j = 0; j < edge_vertex[i].size(); ++j) {
			tet_edge_vertices.emplace_back(i);
			tet_edge_vertices.emplace_back(edge_vertex[i][j]);
		}		
	}
	tet_edge_vertices.shrink_to_fit();

	tet_edge_index_begin_per_thread.resize(thread->thread_num + 1);
	arrangeIndex(thread->thread_num, tet_edge_vertices.size()>>1, tet_edge_index_begin_per_thread.data());
	tet_rest_edge_length.resize(tet_edge_vertices.size() >> 1);
	
	double temp[3];
	for (unsigned int i = 0; i < tet_edge_vertices.size(); i += 2) {
		SUB(temp, vertex_position[tet_edge_vertices[i]], vertex_position[tet_edge_vertices[i + 1]]);
		tet_rest_edge_length[i >> 1] = sqrt(DOT(temp, temp));
	}	

	//for (unsigned int i = 0; i < tet_edge_vertices.size(); i += 2) {
	//	std::cout << tet_edge_vertices[i] << " " << tet_edge_vertices[i + 1] << std::endl;
	//}


}

void TetrahedronMeshStruct::addTetEdges(unsigned int p0, unsigned int p1, std::vector<std::vector<unsigned int>>& edge_vertex)
{
	unsigned int v0, v1;
	if (p0 < p1) {
		v0 = p0; v1 = p1;
	}
	else {
		v0 = p1; v1 = p0;
	}
	if (edge_vertex[v0].empty()) {
		edge_vertex[v0].reserve(10);
		edge_vertex[v0].emplace_back(v1);
	}
	else {
		for (unsigned int i = 0; i < edge_vertex[v0].size(); ++i) {
			if (v1 == edge_vertex[v0][i]) {
				return;
			}
		}
		edge_vertex[v0].emplace_back(v1);
	}
}

void TetrahedronMeshStruct::findSurface()
{
	triangle_indices.reserve(indices.size());
	tet_index_of_surface_face.reserve(indices.size());
	std::map<TetrahedronFace, int> face_in_tet;
	std::vector<int> face_tet_index;
	face_tet_index.reserve(indices.size());
	int* index;
	for (int i = 0; i < indices.size(); ++i) {
		index = indices[i].data();
		buildMap(face_in_tet, index[0], index[2], index[1], face_tet_index, i);
		buildMap(face_in_tet, index[0], index[3], index[2], face_tet_index, i);
		buildMap(face_in_tet, index[0], index[1], index[3], face_tet_index, i);
		buildMap(face_in_tet, index[1], index[2], index[3], face_tet_index, i);
	}
	int j = 0;
	for (auto i = face_in_tet.begin(); i != face_in_tet.end(); ++i) {
		if (i->second == 1) {
			triangle_indices.push_back(i->first.index);
			tet_index_of_surface_face.push_back(face_tet_index[j]);
		}
		j++;
	}

	triangle_indices.shrink_to_fit();
	tet_index_of_surface_face.shrink_to_fit();

	vertex_on_surface.resize(vertex_position.size(), false);
	for (int i = 0; i < triangle_indices.size(); ++i) {
		for (int j = 0; j < 3; ++j) {
			vertex_on_surface[triangle_indices[i][j]] = true;
		}
	}
	vertex_index_on_sureface.reserve(vertex_position.size());
	vertex_surface_index.resize(vertex_position.size(), -1);

	//std::cout << "vertex_surface_index " << vertex_surface_index.size() << std::endl;

	for (unsigned int i = 0; i < vertex_on_surface.size(); ++i) {
		if (vertex_on_surface[i]) {
			vertex_index_on_sureface.push_back(i);
			vertex_surface_index[i] = vertex_index_on_sureface.size() - 1;
		}
	}
	vertex_index_on_sureface.shrink_to_fit();
}


//void TetrahedronMeshStruct::setVertexIndexOnSurfaceEdgeTriangle()
//{
	//edge_vertex_index_on_surface.resize(edges.size());
	//for (int i = 0; i < edge_vertex_index_on_surface.size(); ++i) {
	//	edge_vertex_index_on_surface[i][0] = vertex_surface_index[edges[i].vertex[0]];
	//	edge_vertex_index_on_surface[i][1] = vertex_surface_index[edges[i].vertex[1]];
	//}
	//for (int i = 0; i < surface_triangle_index_in_order.size(); ++i) {
	//	surface_triangle_index_in_order[i][0] = vertex_surface_index[surface_triangle_index_in_order[i][0]];
	//	surface_triangle_index_in_order[i][1] = vertex_surface_index[surface_triangle_index_in_order[i][1]];
	//	surface_triangle_index_in_order[i][2] = vertex_surface_index[surface_triangle_index_in_order[i][2]];
	//}
//}


void TetrahedronMeshStruct::recordTetIndexForSurfaceIndex()
{
	for (int i = 0; i < indices.size(); ++i) {
		for (int j = 0; j < 4; ++j) {
			if (vertex_on_surface[indices[i][j]]) {
				if (vertices[indices[i][j]].tetrahedron.empty()) {
					vertices[indices[i][j]].tetrahedron.reserve(5);
				}
				vertices[indices[i][j]].tetrahedron.push_back(i);
			}
		}
	}
}

void TetrahedronMeshStruct::buildMap(std::map<TetrahedronFace, int>& face_in_tet, int v0, int v1, int v2, std::vector<int>& face_tet_index, int tet_index)
{
	TetrahedronFace A(v0, v1, v2);
	auto ret = face_in_tet.insert({ A,1 });
	if (!ret.second) {
		++ret.first->second;
	}
	else {
		face_tet_index.push_back(tet_index);
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
	return fabs(DOT(c, a)) / 6.0;
}

double TetrahedronMeshStruct::setMass(double density)
{
	volume.resize(indices.size());
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
		mass_inv[i] = 1.0 / mass[i];
	}
	initial_mass_inv = mass_inv;
	return total_mass;
}

void TetrahedronMeshStruct::setThreadIndex(int total_thread_num_)
{
	unsigned int total_thread_num = total_thread_num_ + 1;

	vertex_index_on_surface_begin_per_thread.resize(total_thread_num, 0);
	anchor_index_begin_per_thread.resize(total_thread_num, 0);
	face_index_begin_per_thread.resize(total_thread_num, 0);
	tetrahedron_index_begin_per_thread.resize(total_thread_num, 0);
	vertex_index_begin_per_thread.resize(total_thread_num, 0);


	arrangeIndex(total_thread_num_, vertex_position.size(), vertex_index_begin_per_thread.data());
	arrangeIndex(total_thread_num_, vertex_index_on_sureface.size(), vertex_index_on_surface_begin_per_thread.data());
	arrangeIndex(total_thread_num_, anchor_vertex.size(), anchor_index_begin_per_thread.data());
	arrangeIndex(total_thread_num_, triangle_indices.size(), face_index_begin_per_thread.data());
	arrangeIndex(total_thread_num_, indices.size(), tetrahedron_index_begin_per_thread.data());

	edge_index_begin_per_thread.resize(total_thread_num, 0);
	arrangeIndex(total_thread_num_, edges.size(), edge_index_begin_per_thread.data());

}

void TetrahedronMeshStruct::getRenderNormal()
{
	thread->assignTask(this, FACE_NORMAL_RENDER);
	thread->assignTask(this, VERTEX_NORMAL_RENDER);
}


void TetrahedronMeshStruct::getNormal()
{
	thread->assignTask(this, FACE_NORMAL);
	thread->assignTask(this, VERTEX_NORMAL);
}






void TetrahedronMeshStruct::prepareForDeformationGradient()
{
	//PT_position.resize(indices.size());
	////PPT_determinant.resize(indices.size());
	//PT.resize(indices.size());
	//PT_PPT_inv.resize(indices.size());

	//P_inv.resize(indices.size());
	A.resize(indices.size());
	Matrix<double, 3, 3> p;
	//Matrix<double, 3, 4> A;
	//Matrix3d ppt;
	for (int i = 0; i < indices.size(); ++i)
	{
		p = constructMatrixP(i).inverse();
		//A = constructDeformGradientA(P_inv[i]);
		A[i] = constructDeformGradientA(p);

		//std::cout << (ppt.inverse() * ppt) << std::endl;

		//
		//PT_times_PPT_inv[i] = PT[i] * ppt.inverse();
		//inverse3X3(ppt.data(), PPT_inv[i].data());

	}
}



Matrix<double, 3, 4> TetrahedronMeshStruct::constructDeformGradientA(Matrix3d& p)
{
	Matrix<double, 3, 4> A;
	for (unsigned int i = 0; i < 3; ++i) {
		A.data()[i] = -p.col(i).sum();
		A.data()[3 + i] = p.data()[3 * i];
		A.data()[6 + i] = p.data()[3 * i + 1];
		A.data()[9 + i] = p.data()[3 * i + 2];
	}
	//for (unsigned int i= 0; i < 3; ++i) {
	//	std::cout << A.row(i).sum();
	//}

	return A;
}



Matrix<double, 3, 3> TetrahedronMeshStruct::constructMatrixP(int tetra_index)
{
	Matrix<double, 3, 3> p;
	for (int i = 1; i < 4; ++i)
	{
		SUB((p.data() + 3 * (i-1)), vertex_position[indices[tetra_index][i]].data(), vertex_position[indices[tetra_index][0]].data());
	}
	return p;
	//for (int i = 0; i < 4; ++i)
	//{
	//	memcpy(p.data() + 3 * i, vertex_position[indices[tetra_index][i]].data(), 24);
	//}
	

	//Vector3d center = 0.25 * p.rowwise().sum();
	//return p.colwise() - center;
}

Matrix<double, 3, 4> TetrahedronMeshStruct::constructMatrixP_pos(int tetra_index)
{
	Matrix<double, 3, 4> p;
	for (int i = 0; i < 4; ++i)
	{
		memcpy(p.data() + 3 * i, vertex_position[indices[tetra_index][i]].data(), 24);
	}
	Vector3d center = 0.25 * p.rowwise().sum();
	return p.colwise() - center;
}

//VERTEX_NORMAL
void TetrahedronMeshStruct::getVertexNormalPerThread(int thread_id)
{
	std::vector<unsigned int>* face_vertex;
	double* current_vertex_normal;
	double dot;
	for (int i = vertex_index_on_surface_begin_per_thread[thread_id]; i < vertex_index_on_surface_begin_per_thread[thread_id + 1]; ++i) {
		face_vertex = &vertices[vertex_index_on_sureface[i]].face;
		current_vertex_normal = vertex_normal[vertex_index_on_sureface[i]].data();
		memcpy(current_vertex_normal, face_normal[(*face_vertex)[0]].data(), 24);
		for (int k = 1; k < face_vertex->size(); ++k) {
			SUM_(current_vertex_normal,
				face_normal[(*face_vertex)[k]]);

		}
		normalize(current_vertex_normal);
	}
}

//VERTEX_NORMAL_RENDER
void TetrahedronMeshStruct::getRenderVertexNormalPerThread(int thread_id)
{
	std::vector<unsigned int>* face_vertex;
	double* current_vertex_normal;
	double dot;
	for (int i = vertex_index_on_surface_begin_per_thread[thread_id]; i < vertex_index_on_surface_begin_per_thread[thread_id + 1]; ++i) {
		face_vertex = &vertices[vertex_index_on_sureface[i]].face;
		current_vertex_normal = vertex_normal_for_render[vertex_index_on_sureface[i]].data();
		memcpy(current_vertex_normal, ori_face_normal_for_render[(*face_vertex)[0]].data(), 24);
		for (int k = 1; k < face_vertex->size(); ++k) {
			SUM_(current_vertex_normal,
				ori_face_normal_for_render[(*face_vertex)[k]]);
		}
		dot = DOT(current_vertex_normal, current_vertex_normal);
		if (dot < 1e-20) {
			memcpy(current_vertex_normal, ori_face_normal_for_render[(*face_vertex)[0]].data(), 24);
		}
		else {
			dot = sqrt(dot);
			DEV_(current_vertex_normal, dot);
		}
	}

}