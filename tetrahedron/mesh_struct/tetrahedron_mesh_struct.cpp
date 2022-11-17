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

	unfixed_edge_index_begin_per_thread.resize(thread->thread_num + 1);
	arrangeIndex(thread->thread_num, tet_edge_vertices.size()>>1, unfixed_edge_index_begin_per_thread.data());
	tet_edge_rest_length.resize(tet_edge_vertices.size() >> 1);
	
	double temp[3];
	for (unsigned int i = 0; i < tet_edge_vertices.size(); i += 2) {
		SUB(temp, vertex_position[tet_edge_vertices[i]], vertex_position[tet_edge_vertices[i + 1]]);
		tet_edge_rest_length[i >> 1] = sqrt(DOT(temp, temp));
	}	
	unfixed_rest_edge_length = tet_edge_rest_length;

	only_one_vertex_fixed_edge_index_begin_per_thread.resize(thread->thread_num + 1);
	unfixed_vertex_index_begin_per_thread.resize(thread->thread_num + 1);
	arrangeIndex(thread->thread_num, 0, only_one_vertex_fixed_edge_index_begin_per_thread.data());
	arrangeIndex(thread->thread_num, vertex_position.size(), unfixed_vertex_index_begin_per_thread.data());

	tet_edge_index_begin_per_thread.resize(thread->thread_num + 1);
	arrangeIndex(thread->thread_num, tet_edge_vertices.size() >> 1, tet_edge_index_begin_per_thread.data());

	//for (unsigned int i = 0; i < tet_edge_vertices.size(); i += 2) {
	//	std::cout << tet_edge_vertices[i] << " " << tet_edge_vertices[i + 1] << std::endl;
	//}


}


void TetrahedronMeshStruct::initialUnfixedIndex()
{
	setTetEdges();
	unfixed_point_index.resize(vertex_position.size());
	for (unsigned int i = 0; i < unfixed_point_index.size(); ++i) {
		unfixed_point_index[i] = i;
	}
	unfixed_edge_vertex_index =tet_edge_vertices;
	real_index_to_unfixed_index.resize(vertex_position.size());
	for (unsigned int i = 0; i < unfixed_point_index.size(); ++i) {
		real_index_to_unfixed_index[unfixed_point_index[i]] = i;
	}
	for (unsigned int i = 0; i < anchor_vertex.size(); ++i) {
		real_index_to_unfixed_index[anchor_vertex[i]] = vertex_position.size();
	}
}

void TetrahedronMeshStruct::updateUnfixedPointData()
{
	unfixed_point_index.clear();
	std::vector<bool> is_vertex_fixed(vertex_position.size(),false);
	for (unsigned int i = 0; i < anchor_vertex.size(); ++i) {
		is_vertex_fixed[anchor_vertex[i]] = true;
	}
	for (unsigned int i = 0; i < vertex_position.size(); ++i) {
		if (!is_vertex_fixed[i]) {
			unfixed_point_index.emplace_back(i);
		}
	}

	
	for (unsigned int i = 0; i < unfixed_point_index.size(); ++i) {
		real_index_to_unfixed_index[unfixed_point_index[i]] = i;
	}
	for (unsigned int i = 0; i < anchor_vertex.size(); ++i) {
		real_index_to_unfixed_index[anchor_vertex[i]] = vertex_position.size();
	}

	unfixed_rest_edge_length.clear();
	fixed_one_vertex_rest_edge_length.clear();

	//edge_vertex_index;
	unfixed_edge_vertex_index.clear();
	only_one_vertex_fix_edge.clear();
	only_one_vertex_fix_edge.reserve(anchor_position.size());
	for (unsigned int i = 0; i < tet_edge_vertices.size(); i += 2) {
		if (is_vertex_fixed[tet_edge_vertices[i]]) {
			if (!is_vertex_fixed[tet_edge_vertices[i+1]]) {
				only_one_vertex_fix_edge.emplace_back(real_index_to_unfixed_index[tet_edge_vertices[i + 1]]);
				only_one_vertex_fix_edge.emplace_back(tet_edge_vertices[i]);
				fixed_one_vertex_rest_edge_length.emplace_back(tet_edge_rest_length[i >> 1]);
			}
		}
		else {
			if (is_vertex_fixed[tet_edge_vertices[i + 1]]) {
				only_one_vertex_fix_edge.emplace_back(real_index_to_unfixed_index[tet_edge_vertices[i]]);
				only_one_vertex_fix_edge.emplace_back(tet_edge_vertices[i + 1]);
				fixed_one_vertex_rest_edge_length.emplace_back(tet_edge_rest_length[i >> 1]);
			}
			else {
				unfixed_edge_vertex_index.emplace_back(real_index_to_unfixed_index[tet_edge_vertices[i]]);
				unfixed_edge_vertex_index.emplace_back(real_index_to_unfixed_index[tet_edge_vertices[i + 1]]);
				unfixed_rest_edge_length.emplace_back(tet_edge_rest_length[i >> 1]);
			}
		}
	}



	arrangeIndex(thread->thread_num, unfixed_edge_vertex_index.size() >> 1, unfixed_edge_index_begin_per_thread.data());
	arrangeIndex(thread->thread_num, only_one_vertex_fix_edge.size() >> 1, only_one_vertex_fixed_edge_index_begin_per_thread.data());
	arrangeIndex(thread->thread_num, unfixed_point_index.size(), unfixed_vertex_index_begin_per_thread.data());
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


//SORT_TRIANGLE_AROUND_VERTEX_EDGE
void TetrahedronMeshStruct::sortTetAroundVertexEdge(int thread_id)
{
	auto k = vertex_index_begin_per_thread[thread_id + 1];
	for (auto i = vertex_index_begin_per_thread[thread_id]; i < k; ++i) {
		std::sort(vertex_tet_index[i].begin(), vertex_tet_index[i].end());

	}


}



void TetrahedronMeshStruct::testTetAroundFaceEdge()
{
	std::vector<unsigned int> commen_triangle;
	commen_triangle.reserve(4);
	int* vertex_index;

	std::vector<unsigned int>* vertex_connnect_edge;

	for (auto i = 0; i < triangle_indices.size(); ++i) {
		//face
		vertex_index = triangle_indices[i].data();
		tet_around_face[i].reserve(6);
		for (int j = 0; j < 3; ++j) {
			for (auto k = vertex_tet_index[vertex_index[j]].begin(); k < vertex_tet_index[vertex_index[j]].end(); ++k) {
					tet_around_face[i].emplace_back(*k);				
			}
		}
		std::sort(tet_around_face[i].begin(), tet_around_face[i].end());
		tet_around_face[i].erase(std::unique(tet_around_face[i].begin(), tet_around_face[i].end()), tet_around_face[i].end());
	}


	unsigned int* edge_vertex_index;
	for (auto i = 0; i < edge_length.size(); ++i) {
		edge_vertex_index = edge_vertices.data() + (i << 1);
		tet_around_edge[i].reserve(6);
		for (int j = 0; j < 2; ++j) {
			for (auto k = vertex_tet_index[edge_vertex_index[j]].begin(); k < vertex_tet_index[edge_vertex_index[j]].end(); ++k) {		
					tet_around_edge[i].emplace_back(*k);				
			}
		}
		std::sort(tet_around_edge[i].begin(), tet_around_edge[i].end());
		tet_around_edge[i].erase(std::unique(tet_around_edge[i].begin(), tet_around_edge[i].end()), tet_around_edge[i].end());
	}

	for (unsigned int i = 0; i < triangle_indices.size(); ++i) {
		if (tet_around_face[i] != this->tet_around_face[i]) {
			std::cout << "tet_a_face error " << i << std::endl;
		}
	}
	for (auto i = 0; i < edge_length.size(); ++i) {
		if (tet_around_edge[i] != this->tet_around_edge[i]) {
			std::cout << "tet_a_edge error " << i << std::endl;
		}
	}
	std::cout << "test right" << std::endl;
}


// SORT_TRIANGLE_EDGE_AROUND_TRIANGLE
void TetrahedronMeshStruct::setTetAroundFace(int thread_id)
{
	std::vector<unsigned int> commen_triangle;

	commen_triangle.reserve(4);
	int* vertex_index;

	std::vector<unsigned int>* vertex_connnect_edge;

	for (auto i = face_index_begin_per_thread[thread_id]; i < face_index_begin_per_thread[thread_id + 1]; ++i) {
		//face
		vertex_index = triangle_indices[i].data();
		tet_around_face[i].reserve(6);
		for (int j = 0; j < 3; ++j) {
			tet_around_face[i].insert(tet_around_face[i].end(), vertex_tet_index[vertex_index[j]].begin(),vertex_tet_index[vertex_index[j]].end());
		}
		std::sort(tet_around_face[i].begin(), tet_around_face[i].end());
		tet_around_face[i].erase(std::unique(tet_around_face[i].begin(), tet_around_face[i].end()), tet_around_face[i].end());
	}


	unsigned int* edge_vertex_index;

	for (auto i = edge_index_begin_per_thread[thread_id]; i < edge_index_begin_per_thread[thread_id + 1]; ++i) {
		edge_vertex_index = edge_vertices.data() + (i << 1);
		tet_around_edge[i].reserve(6);
		for (int j = 0; j < 2; ++j) {
			tet_around_edge[i].insert(tet_around_edge[i].end(), vertex_tet_index[edge_vertex_index[j]].begin(), vertex_tet_index[edge_vertex_index[j]].end());
		}
		std::sort(tet_around_edge[i].begin(), tet_around_edge[i].end());
		tet_around_edge[i].erase(std::unique(tet_around_edge[i].begin(), tet_around_edge[i].end()), tet_around_edge[i].end());
	}

}



void TetrahedronMeshStruct::recordTetIndexForVertex()
{
	vertex_tet_index.resize(vertex_position.size());
	for (int i = 0; i < indices.size(); ++i) {
		for (int j = 0; j < 4; ++j) {
			if (vertex_tet_index[indices[i][j]].empty()) {
				vertex_tet_index[indices[i][j]].reserve(8);
			}
			vertex_tet_index[indices[i][j]].emplace_back(i);
		}		
	}
}

void TetrahedronMeshStruct::recordTetIndexForTet()
{
	std::vector<bool> is_used(indices.size(),false);
	tet_tet_index.resize(indices.size());
	for (unsigned int i = 0; i < indices.size(); ++i) {
		std::fill(is_used.begin(), is_used.end(), false);
		is_used[i] = true;
		if (tet_tet_index[i].empty()) {
			tet_tet_index[i].reserve(8);
		}
		for (unsigned int j = 0; j < 4; ++j) {
			for (unsigned int k = 0; k < vertex_tet_index[indices[i][j]].size(); ++k) {
				if (!is_used[vertex_tet_index[indices[i][j]][k]]) {
					tet_tet_index[i].emplace_back(vertex_tet_index[indices[i][j]][k]);
					is_used[vertex_tet_index[indices[i][j]][k]] = true;
				}
			}
		}
	}
}


void TetrahedronMeshStruct::updateTetNeighborInfo()
{
	tet_unfixed_vertex_num.resize(indices.size());
	memset(tet_unfixed_vertex_num.data(), 0, 4 * tet_unfixed_vertex_num.size());
	unsigned int j = 0;
	for (auto itr = indices.begin(); itr < indices.end(); ++itr) {
		for (unsigned int i = 0; i < 4; ++i) {
			if (mass_inv[itr->data()[i]] != 0.0) {
				tet_unfixed_vertex_num[j]++;
			}
		}
		j++;
	}
	updateUnfixedTetVertexIndexInfo();
	updateTetNeighborTetVertexIndex();

//	std::cout << "anchor vertex ";
//for (unsigned int i = 0; i < anchor_vertex.size(); ++i) {
//	std::cout << anchor_vertex[i] << " ";
//}
//
//	std::cout << std::endl;
//	std::cout << "tet neighbor " << std::endl;
//
//	for (unsigned int i = 0; i < indices.size(); ++i) {
//		std::cout << i << std::endl;
//		unsigned int* index = tet_neighbor_tet_vertex_order[i].data();
//		for (unsigned int j = 0; j < tet_tet_index[i].size(); ++j) {
//			int num = *index;
//			index++;
//			for (int i = 0; i < 2*num; ++i) {
//				std::cout << *(index + i) << " ";
//			}
//			std::cout << std::endl;
//			std::cout << indices[tet_tet_index[i][j]][0] << " " << indices[tet_tet_index[i][j]][1] << " " << indices[tet_tet_index[i][j]][2] << " " << indices[tet_tet_index[i][j]][3] << " ";
//			for (unsigned int k = 0; k < 4; ++k) {
//				if (mass_inv[indices[i][k]] != 0) {
//					std::cout << indices[i][k] << " ";
//				}
//			}
//			std::cout<<std::endl;
//			index += 2* num;
//		}
//		std::cout << std::endl;
//	}




	//std::cout << "tet indices " << std::endl;
	//for (unsigned int i = 0; i < indices.size(); ++i) {
	//	std::cout << i << " " << indices[i][0] << " " << indices[i][1] << " " << indices[i][2] << " " <<
	//		indices[i][3] << " " << std::endl;
	//	std::cout <<"= "<< tet_unfixed_vertex_num[i]<<" "<< unfixied_indices[i][0] << " " << unfixied_indices[i][1] << " " << unfixied_indices[i][2] << " " <<
	//		unfixied_indices[i][3] << " " << std::endl;
	//}
}


void TetrahedronMeshStruct::updateUnfixedTetVertexIndexInfo()
{
	unfixied_indices.resize(indices.size(), {-1,-1,-1,-1});
	unfixied_actual_indices.resize(indices.size(), {-1,-1,-1,-1});
	int* vertex_start;
	int* actual_vertex_start;
	for (unsigned int i = 0; i < indices.size(); ++i) {
		vertex_start = unfixied_indices[i].data();
		actual_vertex_start = unfixied_actual_indices[i].data();
		for (unsigned int j = 0; j < 4; ++j) {
			if (mass_inv[indices[i][j]] != 0) {
				*(vertex_start ++)= j;
				*(actual_vertex_start++)= indices[i][j];
			}
		}
	}
}


void TetrahedronMeshStruct::updateTetNeighborTetVertexIndex()
{
	tet_neighbor_tet_vertex_order.resize(indices.size());
	

	thread->assignTask(this, TET_NEIGHBOR_TET_VERTEX_INDEX);
}


//TET_NEIGHBOR_TET_VERTEX_INDEX
void TetrahedronMeshStruct::updateTetNeighborTetVertexIndex(int thread_id)
{
	for (unsigned int i = tetrahedron_index_begin_per_thread[thread_id]; i < tetrahedron_index_begin_per_thread[thread_id+1]; ++i) {
		tet_neighbor_tet_vertex_order[i].clear();
		tet_neighbor_tet_vertex_order[i].reserve(3 * tet_tet_index[i].size());
		for (auto neighbor_tet = tet_tet_index[i].begin(); neighbor_tet < tet_tet_index[i].end(); ++neighbor_tet) {
			findCommonVertexInOrder(indices[i].data(), indices[*neighbor_tet].data(), &tet_neighbor_tet_vertex_order[i],mass_inv.data());
		}
	}
}

void TetrahedronMeshStruct::findCommonVertexInOrder(int* tet_0_index, int* tet_1_index, std::vector<unsigned int>* vertex_in_order, double* inv_mass)
{
	vertex_in_order->emplace_back(0);
	std::vector<unsigned int> vertex;
	vertex.reserve(4);
	unsigned int num = 0;
	for (unsigned int i = 0; i < 4; ++i) {
		if (inv_mass[tet_0_index[i]] != 0.0) {
			for (unsigned int j = 0; j < 4; ++j) {
				if (tet_0_index[i] == tet_1_index[j]) {
					vertex_in_order->emplace_back(j);
					vertex.emplace_back(unfixedVertexIndexInATet(tet_0_index, i, inv_mass));
					num++;
					break;
				}
			}
		}
	}
	*(vertex_in_order->end() - num - 1) = num;
	vertex_in_order->insert(vertex_in_order->end(),vertex.begin(), vertex.end());
}

unsigned int TetrahedronMeshStruct::unfixedVertexIndexInATet(int* tet, int index, double* inv_mass)
{
	unsigned int unfixed_index = index;
	for (unsigned int i = 0; i < index; ++i) {
		if (inv_mass[tet[i]] == 0.0) {
			unfixed_index--;
		}
	}
	return unfixed_index;
}

void TetrahedronMeshStruct::recordTetIndexForSurfaceIndex()
{
	for (int i = 0; i < indices.size(); ++i) {
		for (int j = 0; j < 4; ++j) {
			if (vertex_on_surface[indices[i][j]]) {
				if (vertices[indices[i][j]].tetrahedron.empty()) {
					vertices[indices[i][j]].tetrahedron.reserve(5);
				}
				vertices[indices[i][j]].tetrahedron.emplace_back(i);
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
		if (volume[i] < 0) {
			std::cout << "error " << std::endl;
		}
	}
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
	//initial_mass_inv = mass_inv;
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

//VERTEX_NORMAL_FROM_RENDER
void TetrahedronMeshStruct::getVertexNormalFromRenderPerThread(int thread_id)
{
	std::vector<unsigned int>* face_vertex;
	double* current_vertex_normal;
	double dot;
	for (int i = vertex_index_on_surface_begin_per_thread[thread_id]; i < vertex_index_on_surface_begin_per_thread[thread_id + 1]; ++i) {
		face_vertex = &vertices[vertex_index_on_sureface[i]].face;
		current_vertex_normal = vertex_normal[vertex_index_on_sureface[i]].data();
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


void TetrahedronMeshStruct::sortTriangleAroundElement()
{
	tet_around_face.resize(triangle_indices.size());
	tet_around_edge.resize(edge_length.size());


	thread->assignTask(this, SORT_TRIANGLE_EDGE_AROUND_TRIANGLE_EDGE);
	thread->assignTask(this, SORT_TRIANGLE_AROUND_VERTEX_EDGE);


	std::cout <<"tet around "<< tet_around_face[0].size() << " " << tet_around_edge[0].size() << std::endl;

}


void TetrahedronMeshStruct::recordTriangleIndexOfATet()
{
	triangle_index_of_a_tet.resize(indices.size());
	std::vector<bool>triangle_used(triangle_indices.size(), false);
	std::vector<unsigned int>* vertex_face;
	for (unsigned int i = 0; i < triangle_index_of_a_tet.size(); ++i) {
		std::fill(triangle_used.begin(), triangle_used.end(), false);
		for (unsigned int j = 0; j < 4; ++j) {
			if (vertex_on_surface[indices[i][j]]) {
				vertex_face = &vertices[indices[i][j]].face;
				for (auto k = vertex_face->begin(); k < vertex_face->end(); ++k) {
					if (!triangle_used[*k]) {
						triangle_index_of_a_tet[i].emplace_back(*k);
						triangle_used[*k] = true;
					}
				}
			}
		}
	}
}


void TetrahedronMeshStruct::recordPrimitiveIndexOfATet()
{
	recordTriangleIndexOfATet();
	recordEdgeIndexOfATet();
	//recordEdgeIndexInATet();
}

void TetrahedronMeshStruct::recordEdgeIndexOfATet()
{
	edge_index_of_a_tet.resize(indices.size());
	std::vector<bool>edge_used(edge_vertices.size()>>1, false);
	std::vector<unsigned int>* vertex_edge;
	for (unsigned int i = 0; i < edge_index_of_a_tet.size(); ++i) {
		std::fill(edge_used.begin(), edge_used.end(), false);
		for (unsigned int j = 0; j < 4; ++j) {
			if (vertex_on_surface[indices[i][j]]) {
				vertex_edge = &vertices[indices[i][j]].edge;
				for (auto k = vertex_edge->begin(); k < vertex_edge->end(); ++k) {
					if (!edge_used[*k]) {
						edge_index_of_a_tet[i].emplace_back(*k);
						edge_used[*k] = true;
					}
				}
			}
		}
	}
}

//void TetrahedronMeshStruct::recordEdgeIndexInATet()
//{
//	edge_index_in_a_tet.resize(indices.size());
//	for (unsigned int i = 0; i < edge_index_in_a_tet.size(); ++i) {
//
//		for (auto j = edge_index_of_a_tet[i].begin(); j < edge_index_of_a_tet[i].end(); ++j) {
//			if (checkEdgeInATet(edge_vertices.data() + (*j) * 2, indices[i].data())) {
//				edge_index_in_a_tet[i].emplace_back(*j);
//			}
//		}
//	}
//}

bool TetrahedronMeshStruct::checkEdgeInATet(unsigned int* edge_vertices, int* tet_indices)
{
	for (unsigned int i = 0; i < 4; ++i) {
		if (edge_vertices[0] == tet_indices[i]) {
			for (unsigned int j = 0; j < 4; ++j) {
				if (edge_vertices[1] == tet_indices[j]) {
					return true;
				}
			}
		}
	}
	return false;
}