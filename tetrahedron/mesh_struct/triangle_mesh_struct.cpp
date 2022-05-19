#include"triangle_mesh_struct.h"
#include"../basic/enum_setting.h"

TriangleMeshStruct::TriangleMeshStruct()
{
	type = TRIANGLE;
}



double TriangleMeshStruct::setMass(double density)
{
	double mass_ = 0.0;
	double m;
	if (!faces.empty()) {
		for (int i = 0; i < faces.size(); ++i) {
			m= faces[i].area * density /3.0;
			mass[triangle_indices[i][0]] += m;
			mass[triangle_indices[i][1]] += m;
			mass[triangle_indices[i][2]] += m;
			mass_ += faces[i].area * density;
		}
		//std::cout << mesh_struct.faces[0].area + mesh_struct.faces[1].area << std::endl;
	}
	else {
		for (int i = 0; i < vertices.size(); ++i) {
			mass[i] = 1.25;
			mass_ += mass[i];
		}
	}
	for (int i = 0; i < vertices.size(); ++i) {
		mass_inv[i] = 1.0 / mass[i];
	}
	initial_mass_inv = mass_inv;
	return mass_;
}


//VERTEX_NORMAL
void TriangleMeshStruct::getVertexNormalPerThread(int thread_id)
{
	std::vector<unsigned int>* face_vertex;
	double* current_vertex_normal;
	double dot;
	for (int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id + 1]; ++i) {
		face_vertex = &vertices[i].face;
		if (!face_vertex->empty()) {
			current_vertex_normal = vertex_normal[i].data();
			memcpy(current_vertex_normal, face_normal[(*face_vertex)[0]].data(), 24);
			for (int k = 1; k < face_vertex->size(); ++k) {
				SUM_(current_vertex_normal,
					face_normal[(*face_vertex)[k]]);

			}
			normalize(current_vertex_normal);
		}
	}
}

//VERTEX_NORMAL_RENDER
void TriangleMeshStruct::getRenderVertexNormalPerThread(int thread_id)
{
	std::vector<unsigned int>* face_vertex;
	double* current_vertex_normal;
	double dot;
	for (int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id + 1]; ++i) {
		face_vertex = &vertices[i].face;
		if (!face_vertex->empty()) {
			current_vertex_normal = vertex_normal_for_render[i].data();
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

}

void TriangleMeshStruct::initialUnfixedIndex()
{
	int total_thread_num_ = thread->thread_num;
	unfixed_point_index.resize(vertex_position.size());
	for (unsigned int i = 0; i < unfixed_point_index.size(); ++i) {
		unfixed_point_index[i] = i;
	}
	unfixed_edge_vertex_index = edge_vertices;
	unfixed_rest_edge_length = edge_length;

	only_one_vertex_fixed_edge_index_begin_per_thread.resize(total_thread_num_+1, 0);
	unfixed_edge_index_begin_per_thread.resize(total_thread_num_+1, 0);
	arrangeIndex(total_thread_num_, edges.size(), unfixed_edge_index_begin_per_thread.data());
	arrangeIndex(total_thread_num_, 0, only_one_vertex_fixed_edge_index_begin_per_thread.data());
	unfixed_vertex_index_begin_per_thread.resize(total_thread_num_+1, 0);
	arrangeIndex(total_thread_num_, vertex_position.size(), unfixed_vertex_index_begin_per_thread.data());
	updateUnfixedPointData();

}


void TriangleMeshStruct::updateUnfixedPointData()
{
	unfixed_point_index.clear();
	std::vector<bool> is_vertex_fixed(vertex_position.size(), false);
	for (unsigned int i = 0; i < anchor_vertex.size(); ++i) {
		is_vertex_fixed[anchor_vertex[i]] = true;
	}
	for (unsigned int i = 0; i < vertex_position.size(); ++i) {
		if (!is_vertex_fixed[i]) {
			unfixed_point_index.emplace_back(i);
		}
	}

	std::vector<unsigned int> real_index_to_unfixed_index(vertex_position.size());
	for (unsigned int i = 0; i < unfixed_point_index.size(); ++i) {
		real_index_to_unfixed_index[unfixed_point_index[i]] = i;
	}

	unfixed_rest_edge_length.clear();
	fixed_one_vertex_rest_edge_length.clear();

	//edge_vertex_index;
	unfixed_edge_vertex_index.clear();
	only_one_vertex_fix_edge.clear();
	only_one_vertex_fix_edge.reserve(anchor_position.size());
	for (unsigned int i = 0; i < edge_vertices.size(); i += 2) {
		if (is_vertex_fixed[edge_vertices[i]]) {
			if (!is_vertex_fixed[edge_vertices[i + 1]]) {
				only_one_vertex_fix_edge.emplace_back(real_index_to_unfixed_index[edge_vertices[i + 1]]);
				only_one_vertex_fix_edge.emplace_back(edge_vertices[i]);
				fixed_one_vertex_rest_edge_length.emplace_back(edge_length[i >> 1]);
			}
		}
		else {
			if (is_vertex_fixed[edge_vertices[i + 1]]) {
				only_one_vertex_fix_edge.emplace_back(real_index_to_unfixed_index[edge_vertices[i]]);
				only_one_vertex_fix_edge.emplace_back(edge_vertices[i + 1]);
				fixed_one_vertex_rest_edge_length.emplace_back(edge_length[i >> 1]);
			}
			else {
				unfixed_edge_vertex_index.emplace_back(real_index_to_unfixed_index[edge_vertices[i]]);
				unfixed_edge_vertex_index.emplace_back(real_index_to_unfixed_index[edge_vertices[i + 1]]);
				unfixed_rest_edge_length.emplace_back(edge_length[i >> 1]);
			}
		}
	}
	arrangeIndex(thread->thread_num, only_one_vertex_fix_edge.size() >> 1, only_one_vertex_fixed_edge_index_begin_per_thread.data());
	arrangeIndex(thread->thread_num, unfixed_edge_vertex_index.size() >> 1, unfixed_edge_index_begin_per_thread.data());
	arrangeIndex(thread->thread_num, unfixed_point_index.size(), unfixed_vertex_index_begin_per_thread.data());

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




////FACE_NORMAL_RENDER
//void TriangleMeshStruct::getRenderFaceNormalPerThread(int thread_id)
//{
//	double e0[3], e2[3];
//	double* current_face_normal;
//	int* triangle_vertex;
//	floating* f_current_face_normal;
//	floating f_e2[3];
//	floating f_e0[3];
//	for (int j = face_index_begin_per_thread[thread_id]; j < face_index_begin_per_thread[thread_id + 1]; ++j) {
//		triangle_vertex = triangle_indices[j].data();
//		current_face_normal = face_normal_for_render[j].data();
//		f_current_face_normal = f_face_normal_for_render[j].data();
//		SUB(e2, vertex_for_render[triangle_vertex[1]], vertex_for_render[triangle_vertex[0]]);
//		SUB(e0, vertex_for_render[triangle_vertex[2]], vertex_for_render[triangle_vertex[0]]);
//		for (int i = 0; i < 3; ++i) {
//			f_e2[i].v = e2[i]; 	f_e2[i].sigma = 0.0;
//			f_e0[i].v = e0[i]; 	f_e0[i].sigma = 0.0;
//		}
//		CROSS(f_current_face_normal, f_e2, f_e0);
//		for (int i = 0; i < 3; ++i) {
//			current_face_normal[i] = f_current_face_normal[i].v;
//		}
//		memcpy(ori_face_normal_for_render[j].data(), current_face_normal, 24);
//		normalize(current_face_normal);
//	}
//}

////FACE_NORMAL
//void TriangleMeshStruct::getFaceNormalPerThread(int thread_id)
//{
//	double e0[3], e2[3];
//	double e3[3], e4[3];
//	double* current_face_normal;
//	floating* f_current_face_normal;
//	int* triangle_vertex;//
//	floating f_e0[3], f_e2[3];
//	floating f_e3[3], f_e4[3], f_cross[3];//
//	for (int j = face_index_begin_per_thread[thread_id]; j < face_index_begin_per_thread[thread_id+1]; ++j) {
//		triangle_vertex = triangle_indices[j].data();
//		current_face_normal = face_normal[j].data();
//		f_current_face_normal = f_face_normal[j].data();
//		SUB(e2, vertex_position[triangle_vertex[1]], vertex_position[triangle_vertex[0]]);
//		SUB(e0, vertex_position[triangle_vertex[2]], vertex_position[triangle_vertex[0]]);
//
//		for (int i = 0; i < 3; ++i) {
//			f_e2[i].v = e2[i]; 	f_e2[i].sigma = 0.0;
//			f_e0[i].v = e0[i]; 	f_e0[i].sigma = 0.0;
//		}
//
//		CROSS(f_current_face_normal, f_e2, f_e0);
//		for (int i = 0; i < 3; ++i) {
//			current_face_normal[i] = f_current_face_normal[i].v;
//		}
//
//		memcpy(ori_face_normal[j].data(), current_face_normal, 24);
//		triangle_normal_magnitude_reciprocal[j] = 1.0 / sqrt(DOT(current_face_normal, current_face_normal));
//		normalize(current_face_normal);
//
//		SUB(e3, vertex_for_render[triangle_vertex[1]], vertex_for_render[triangle_vertex[0]]);
//		SUB(e4, vertex_for_render[triangle_vertex[2]], vertex_for_render[triangle_vertex[0]]);
//
//		for (int i = 0; i < 3; ++i) {
//			f_e3[i].v = e3[i]; 	f_e3[i].sigma = 0.0;
//			f_e4[i].v = e4[i]; 	f_e4[i].sigma = 0.0;
//		}
//
//
//		f_current_face_normal = f_cross_for_approx_CCD[j].data();
//		current_face_normal=cross_for_approx_CCD[j].data();
//
//		CROSS(f_current_face_normal, f_e3, f_e0);
//		CROSS(f_cross, f_e2, f_e4);
//		SUM_(f_current_face_normal, f_cross);
//		for (int i = 0; i < 3; ++i) {
//			current_face_normal[i] = f_current_face_normal[i].v;
//		}
//	}
//}



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
	vertex_normal = vertex_normal_for_render;
	face_normal = face_normal_for_render;
}


void TriangleMeshStruct::setThreadIndex(int total_thread_num_)
{
	int total_thread_num = total_thread_num_ + 1;

	vertex_index_begin_per_thread.resize(total_thread_num, 0);
	anchor_index_begin_per_thread.resize(total_thread_num, 0);
	face_index_begin_per_thread.resize(total_thread_num, 0);

	arrangeIndex(total_thread_num_, vertices.size(), vertex_index_begin_per_thread.data());
	arrangeIndex(total_thread_num_, anchor_vertex.size(), anchor_index_begin_per_thread.data());
	arrangeIndex(total_thread_num_, triangle_indices.size(), face_index_begin_per_thread.data());

	if (!edges.empty()) {
		edge_index_begin_per_thread.resize(total_thread_num, 0);
		arrangeIndex(total_thread_num_, edges.size(), edge_index_begin_per_thread.data());

	}
}

