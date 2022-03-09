#include"triangle_mesh_struct.h"
#include"../basic/enum_setting.h"

TriangleMeshStruct::TriangleMeshStruct()
{
	type = TRIANGLE;
}


//VERTEX_NORMAL
void TriangleMeshStruct::getVertexNormalPerThread(int thread_id)
{
	std::vector<int>* face_vertex;
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
	std::vector<int>* face_vertex;
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

