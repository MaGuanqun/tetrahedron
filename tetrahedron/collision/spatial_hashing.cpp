#include"spatial_hashing.h"
//#include"../basic/write_txt.h"


void SpatialHashing::initialHashCellLength(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, double& cell_length, double* tolerance_ratio)
{
	double max_length = 0;
	double ave_length = 0;
	int edge_num = 0;
	double* edge_length;
	unsigned int edge_size;
	for (int i = 0; i < cloth->size(); ++i) {
		edge_size = cloth->data()[i].mesh_struct.edge_length.size();
		edge_length = cloth->data()[i].mesh_struct.edge_length.data();
		for (int j = 0; j < edge_size; ++j) {
			if (max_length < edge_length[j]) {
				max_length = edge_length[j];
			}
			ave_length += edge_length[j];
		}
		edge_num += edge_size;
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		edge_size = tetrahedron->data()[i].mesh_struct.edge_length.size();
		edge_length = tetrahedron->data()[i].mesh_struct.edge_length.data();
		for (int j = 0; j < edge_size; ++j) {
			if (max_length < edge_length[j]) {
				max_length = edge_length[j];
			}
			ave_length += edge_length[j];
		}
		edge_num += edge_size;
	}

	ave_length /= (double)edge_num;
	cell_length = max_length + 2.0 * tolerance_ratio[AABB_SELF_PT] * ave_length;
	//cell_length = 1.0 * max_length +2.0 * tolerance_ratio[SELF_POINT_TRIANGLE] * ave_length;
	std::cout << "AABB tolerance__ " << tolerance_ratio[AABB_SELF_PT] << std::endl;
	std::cout << "ave_length " << ave_length << " max length " << max_length << " ratio " << (double)max_length / (double)ave_length << std::endl;
	triangle_pair_number_thread.resize(thread_num);


}


void SpatialHashing::setSpatialHashingInitialCount()
{
	//getSceneAABB();
	for (unsigned int i = 0; i < 3; ++i) {
		cell_number[i] = (unsigned int)floor((scene_aabb[i + 3] - scene_aabb[i]) / cell_length) + 1;
	}

	cell_num0_cell_num1 = cell_number[0] * cell_number[1];
//	std::cout << "spatial hashing size " << hash_table_size << " " << cell_number[0] << " " << cell_number[1] << " " << cell_number[2] << std::endl;

}


void SpatialHashing::initialHashCell(unsigned int total_triangle_num, unsigned int max_index_number_in_one_cell,
	unsigned int max_index_number_in_one_cell_collider, unsigned int estimate_coeff_for_vt_pair_num, unsigned int estimate_coeff_for_vt_collider_pair_num,
	unsigned int estimate_coeff_for_tv_pair_num, unsigned int estimate_coeff_for_ee_pair_num)
{
	//indicator = new std::vector<unsigned int>[thread_num];
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	indicator[i].resize(1);
	//}

	unsigned int temp_hash_count = total_triangle_num / 4;
	if (temp_hash_count < 50) {
		temp_hash_count = 50;
	}
	hash_cell_count = findNeareastPrimeNumber(temp_hash_count);// 90191;// 24877;// 500009; //;// // ;// 32999;//  // 11987;//84947
	//hash_cell_count = 90191;
	std::cout << "total triangle " << total_triangle_num << " hash cell count " << hash_cell_count << std::endl;
	this->max_index_number_in_one_cell = max_index_number_in_one_cell;// = 801;
	this->max_index_number_in_one_cell_collider = max_index_number_in_one_cell_collider;// = 401;
	this->max_index_number_in_one_cell_edge = max_index_number_in_one_cell;

	if (searchPairByCell) {
		this->max_index_number_in_one_cell_vertex = max_index_number_in_one_cell / 2 + 1;// 	
		this->max_index_number_in_one_cell_collider_vertex = max_index_number_in_one_cell_collider / 2 + 1;//
		this->max_index_number_in_one_cell_collider_edge = max_index_number_in_one_cell_collider;//
	}
	else {
		hash_count_start_per_thread.resize(thread_num + 1);
		arrangeIndex(thread_num, hash_cell_count, hash_count_start_per_thread.data());
	}

	spatial_hashing_cell_triangle = new unsigned int* [thread_num];
	spatial_hashing_cell_collider_triangle = new unsigned int* [thread_num];
	spatial_hashing_cell_edge = new unsigned int* [thread_num];
	for (unsigned int i = 0; i < thread_num; ++i) {
		spatial_hashing_cell_triangle[i] = new unsigned int[hash_cell_count * max_index_number_in_one_cell];
		memset(spatial_hashing_cell_triangle[i], 0, 4 * max_index_number_in_one_cell * hash_cell_count);
		spatial_hashing_cell_edge[i] = new unsigned int[hash_cell_count * max_index_number_in_one_cell_edge];
		memset(spatial_hashing_cell_edge[i], 0, 4 * max_index_number_in_one_cell_edge * hash_cell_count);
		if (has_collider) {
			spatial_hashing_cell_collider_triangle[i] = new unsigned int[hash_cell_count * max_index_number_in_one_cell_collider];
			memset(spatial_hashing_cell_collider_triangle[i], 0, 4 * hash_cell_count * max_index_number_in_one_cell_collider);
		}
	}

	if (!record_pair_by_element) {
		vertex_triangle_pair = new unsigned int* [thread_num];
		vertex_obj_triangle_collider_pair = new unsigned int* [thread_num];
		//vertex_collider_triangle_obj_pair = new unsigned int* [thread_num];
		edge_edge_pair = new unsigned int* [thread_num];
		//edge_edge_pair_collider = new unsigned int* [thread_num];
		for (unsigned int i = 0; i < thread_num; ++i) {
			edge_edge_pair[i] = new unsigned int[2 * estimate_coeff_for_ee_pair_num * total_triangle_num];// 
			memset(edge_edge_pair[i], 0, 4 * 2 * estimate_coeff_for_ee_pair_num * total_triangle_num);
			vertex_triangle_pair[i] = new unsigned int[estimate_coeff_for_ee_pair_num * total_triangle_num / 2];// 
			memset(vertex_triangle_pair[i], 0, 4 * (estimate_coeff_for_ee_pair_num * total_triangle_num / 2));
			if (has_collider) {
				//edge_edge_pair_collider[i] = new unsigned int[estimate_coeff_for_pair_num * total_triangle_num];
				//memset(edge_edge_pair_collider[i], 0, 4 * estimate_coeff_for_pair_num * total_triangle_num);
				vertex_obj_triangle_collider_pair[i] = new unsigned int[estimate_coeff_for_ee_pair_num * total_triangle_num / 2];
				memset(vertex_obj_triangle_collider_pair[i], 0, 4 * (estimate_coeff_for_ee_pair_num * total_triangle_num / 2));
				//vertex_collider_triangle_obj_pair[i] = new unsigned int[estimate_coeff_for_pair_num * total_triangle_num / 2];
				//memset(vertex_collider_triangle_obj_pair[i], 0, 4 * (estimate_coeff_for_pair_num * total_triangle_num / 2));			
			}
			else {
				//edge_edge_pair_collider[i] = new unsigned int[1];
				vertex_obj_triangle_collider_pair[i] = new unsigned int[1];
				//vertex_collider_triangle_obj_pair[i] = new unsigned int[1];
			}
		}
	}
	else {
		vertex_triangle_pair_by_vertex = new unsigned int* [collider_begin_obj_index];
		triangle_vertex_pair_by_triangle = new unsigned int* [collider_begin_obj_index];
		edge_edge_pair_by_edge = new unsigned int* [collider_begin_obj_index];
		vertex_obj_triangle_collider_pair_by_vertex = new unsigned int* [collider_begin_obj_index];

		vertex_triangle_pair_num_record = new unsigned int* [collider_begin_obj_index];
		triangle_vertex_pair_num_record = new unsigned int* [collider_begin_obj_index];
		edge_edge_pair_num_record = new unsigned int* [collider_begin_obj_index];
		vertex_obj_triangle_collider_num_record = new unsigned int* [collider_begin_obj_index];


		total_length_every_element_for_vertex_triange = estimate_coeff_for_vt_pair_num;
		total_length_every_element_for_vertex_triange_collider = estimate_coeff_for_vt_collider_pair_num;
		total_length_every_element_for_triangle_vertex = estimate_coeff_for_tv_pair_num;
		total_length_every_element_for_edge_edge = estimate_coeff_for_ee_pair_num;

		unsigned int vertex_num,triangle_num,edge_num;

		for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
			edge_edge_pair_by_edge[i] = new unsigned int[total_length_every_element_for_edge_edge * mesh_struct[i]->edge_length.size()];// 
			memset(edge_edge_pair_by_edge[i], 0, 4 * total_length_every_element_for_edge_edge * mesh_struct[i]->edge_length.size());

			if (i < cloth->size()) {
				vertex_triangle_pair_by_vertex[i] = new unsigned int[total_length_every_element_for_vertex_triange * mesh_struct[i]->vertex_position.size()];// 
				memset(vertex_triangle_pair_by_vertex[i], 0, 4 * (total_length_every_element_for_vertex_triange * mesh_struct[i]->vertex_position.size()));
				vertex_triangle_pair_num_record[i] = new unsigned int[mesh_struct[i]->vertex_position.size()];// 
				memset(vertex_triangle_pair_num_record[i], 0, 4 * mesh_struct[i]->vertex_position.size());
			}
			else {
				vertex_triangle_pair_by_vertex[i] = new unsigned int[total_length_every_element_for_vertex_triange * tetrahedron->data()[i-cloth->size()].mesh_struct.vertex_index_on_sureface.size()];// 
				memset(vertex_triangle_pair_by_vertex[i], 0, 4 * (total_length_every_element_for_vertex_triange * tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_index_on_sureface.size()));
				vertex_triangle_pair_num_record[i] = new unsigned int[tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_index_on_sureface.size()];// 
				memset(vertex_triangle_pair_num_record[i], 0, 4 * tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_index_on_sureface.size());			
			}
			triangle_vertex_pair_by_triangle[i] = new unsigned int[total_length_every_element_for_triangle_vertex * mesh_struct[i]->triangle_indices.size()];// 
			memset(triangle_vertex_pair_by_triangle[i], 0, 4 * (total_length_every_element_for_triangle_vertex * mesh_struct[i]->triangle_indices.size()));

			edge_edge_pair_num_record[i] = new unsigned int[mesh_struct[i]->edge_length.size()];// 
			memset(edge_edge_pair_num_record[i], 0, 4 * mesh_struct[i]->edge_length.size());
			triangle_vertex_pair_num_record[i] = new unsigned int[mesh_struct[i]->triangle_indices.size()];// 
			memset(triangle_vertex_pair_num_record[i], 0, 4 * mesh_struct[i]->triangle_indices.size());

			if (has_collider) {
				vertex_obj_triangle_collider_pair_by_vertex[i] = new unsigned int[total_length_every_element_for_vertex_triange_collider * mesh_struct[i]->vertex_position.size()];
				memset(vertex_obj_triangle_collider_pair_by_vertex[i], 0, 4 * (total_length_every_element_for_vertex_triange_collider * mesh_struct[i]->vertex_position.size()));

				vertex_obj_triangle_collider_num_record[i] = new unsigned int[mesh_struct[i]->vertex_position.size()];// 
				memset(vertex_obj_triangle_collider_num_record[i], 0, 4 * mesh_struct[i]->vertex_position.size());
			}
			else {
				//edge_edge_pair_collider[i] = new unsigned int[1];
				vertex_obj_triangle_collider_pair_by_vertex[i] = new unsigned int[1];
				vertex_obj_triangle_collider_num_record[i] = new unsigned int[1];
				//vertex_collider_triangle_obj_pair[i] = new unsigned int[1];
			}
		}
	}


	if (searchPairByCell) {		
		spatial_hashing_cell_vertex = new unsigned int* [thread_num];
		spatial_hashing_cell_collider_edge = new unsigned int* [thread_num];
		spatial_hashing_cell_collider_vertex = new unsigned int* [thread_num];
		for (unsigned int i = 0; i < thread_num; ++i) {		
			spatial_hashing_cell_vertex[i] = new unsigned int[hash_cell_count * max_index_number_in_one_cell_vertex];
			memset(spatial_hashing_cell_vertex[i], 0, 4 * max_index_number_in_one_cell_vertex * hash_cell_count);
			if (has_collider) {
				spatial_hashing_cell_collider_edge[i] = new unsigned int[hash_cell_count * max_index_number_in_one_cell_collider_edge];
				spatial_hashing_cell_collider_vertex[i] = new unsigned int[hash_cell_count * max_index_number_in_one_cell_collider_vertex];
				memset(spatial_hashing_cell_collider_edge[i], 0, 4 * hash_cell_count * max_index_number_in_one_cell_collider_edge);
				memset(spatial_hashing_cell_collider_vertex[i], 0, 4 * max_index_number_in_one_cell_collider_vertex);
			}
		}
	}

	spatial_hashing_cell_triangle_size = new unsigned int* [thread_num];
	spatial_hashing_cell_collider_triangle_size = new unsigned int* [thread_num];
	spatial_hashing_cell_edge_size = new unsigned int* [thread_num];
	for (unsigned int i = 0; i < thread_num; ++i) {
		spatial_hashing_cell_triangle_size[i] = new unsigned int[hash_cell_count];
		spatial_hashing_cell_collider_triangle_size[i] = new unsigned int[hash_cell_count];
		memset(spatial_hashing_cell_triangle_size[i], 0, hash_cell_count << 2);
		memset(spatial_hashing_cell_collider_triangle_size[i], 0, hash_cell_count << 2);

		spatial_hashing_cell_edge_size[i] = new unsigned int[hash_cell_count];
		memset(spatial_hashing_cell_edge_size[i], 0, hash_cell_count << 2);

	}

	if (searchPairByCell) {
	
		spatial_hashing_cell_vertex_size = new unsigned int* [thread_num];
		spatial_hashing_cell_collider_edge_size = new unsigned int* [thread_num];
		spatial_hashing_cell_collider_vertex_size = new unsigned int* [thread_num];

		for (unsigned int i = 0; i < thread_num; ++i) {

			spatial_hashing_cell_collider_edge_size[i] = new unsigned int[hash_cell_count];
			memset(spatial_hashing_cell_collider_edge_size[i], 0, hash_cell_count << 2);

			spatial_hashing_cell_vertex_size[i] = new unsigned int[hash_cell_count];
			spatial_hashing_cell_collider_vertex_size[i] = new unsigned int[hash_cell_count];
			memset(spatial_hashing_cell_vertex_size[i], 0, hash_cell_count << 2);
			memset(spatial_hashing_cell_collider_vertex_size[i], 0, hash_cell_count << 2);
		}
	}

	cell_begin_per_thread.resize(thread_num + 1);
	arrangeIndex(thread_num, hash_cell_count, cell_begin_per_thread.data());

	P1 = 73856093; //2147483647
	P2 = 19349663; //500000003
	P3 = 83492791; //900001961


	//spatial_hashing_actual_hash_value_for_test = new std::vector<unsigned int>*[thread_num];
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	spatial_hashing_actual_hash_value_for_test[i] = new std::vector<unsigned int>[hash_cell_count];
	//	for (unsigned int j = 0; j < hash_cell_count; ++j) {
	//		spatial_hashing_actual_hash_value_for_test[i][j].reserve(max_index_number_in_one_cell);

	//	}
	//}

	//pair_vector_size_in_a_thread = new unsigned int[thread_num];
	//pair_count_per_thread = new unsigned int[thread_num];
	//memset(pair_count_per_thread, 0, 4 * thread_num);
	//pair_count_with_collider_per_thread = new unsigned int[thread_num];
	//memset(pair_count_with_collider_per_thread, 0, 4 * thread_num);


	non_empty_cell_index_vertex_triangle = new unsigned int* [thread_num];
	non_empty_cell_index_edge = new unsigned int* [thread_num];
	for (unsigned int i = 0; i < thread_num; ++i) {
		non_empty_cell_index_vertex_triangle[i] = new unsigned int[hash_cell_count + 1];
		non_empty_cell_index_edge[i] = new unsigned int[hash_cell_count + 1];
	}

	non_empty_cell_index_begin_per_thread_vertex_triangle = new unsigned int[thread_num + 1];
	non_empty_cell_index_begin_per_thread_edge = new unsigned int[thread_num + 1];
	non_empty_cell_index_ave_begin_per_thread_vertex_triangle = new unsigned int[thread_num + 1];
	non_empty_cell_index_ave_begin_per_thread_edge = new unsigned int[thread_num + 1];
	memset(non_empty_cell_index_ave_begin_per_thread_vertex_triangle, 0, 4 * (thread_num + 1));
	memset(non_empty_cell_index_ave_begin_per_thread_edge, 0, 4 * (thread_num + 1));
	hash_cell_triangle_count.reserve(hash_cell_count);

	hash_cell_pair_num_prefix_vertex_triangle = new unsigned int[hash_cell_count + 2];
	hash_cell_pair_num_prefix_edge = new unsigned int[hash_cell_count + 2];
	//hash_cell_triangle_num.reserve(hash_cell_count);


	prefix_sum_record_per_thread_start_vertex_triangle = new unsigned int[thread_num];
	prefix_sum_record_per_thread_start_edge = new unsigned int[thread_num];
	memset(prefix_sum_record_per_thread_start_vertex_triangle, 0, thread_num << 2);
	memset(prefix_sum_record_per_thread_start_edge, 0, thread_num << 2);

	//max_pair_num_in_one_loop_per_thread = 1000000;
	global_cell_start_vertex_triangle = new unsigned int[thread_num << 1];
	global_cell_start_edge = new unsigned int[thread_num << 1];
	//max_num_to_loop_to_find_pair = 400;
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	global_cell_start[i] = new unsigned int[max_num_to_loop_to_find_pair + 2];
	//	global_cell_start[i][0] = 0;
	//}

	primitive_index_record = new unsigned int* [thread_num];
	for (unsigned int i = 0; i < thread_num; ++i) {
		primitive_index_record[i] = new unsigned int[1600];
	}

}



void SpatialHashing::setInObject(std::vector<Cloth>* cloth, std::vector<Collider>* collider,
	std::vector<Tetrahedron>* tetrahedron, Thread* thread, double* tolerance_ratio,
	unsigned int max_cell_count,
	unsigned int max_index_number_in_one_cell,
	unsigned int max_index_number_in_one_cell_collider, unsigned int estimate_coeff_for_vt_pair_num, unsigned int estimate_coeff_for_vt_collider_pair_num,
	unsigned int estimate_coeff_for_tv_pair_num, unsigned int estimate_coeff_for_ee_pair_num, bool record_pair_by_element)
{
	this->record_pair_by_element = record_pair_by_element;
	has_collider = !collider->empty();

	this->for_construct_patch = for_construct_patch;
	tetrahedron_begin_obj_index = cloth->size();
	this->cloth = cloth;
	this->collider = collider;
	this->tetrahedron = tetrahedron;
	this->thread = thread;
	thread_num = thread->thread_num;
	scene_aabb_thread.resize(thread_num);
	total_obj_num = cloth->size() + tetrahedron->size() + collider->size();

	collider_begin_obj_index = cloth->size() + tetrahedron->size();

	this->max_cell_count = max_cell_count;
	reorganzieDataOfObjects();

	initialHashCellLength(cloth, tetrahedron, cell_length, tolerance_ratio);


	int total_triangle_num = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		total_triangle_num += cloth->data()[i].mesh_struct.triangle_indices.size();
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		total_triangle_num += tetrahedron->data()[i].mesh_struct.triangle_indices.size();
	}

	for (int i = 0; i < collider->size(); ++i) {
		total_triangle_num += collider->data()[i].mesh_struct.triangle_indices.size();
	}

	initialHashCell(total_triangle_num, max_index_number_in_one_cell, max_index_number_in_one_cell_collider, estimate_coeff_for_vt_pair_num, estimate_coeff_for_vt_collider_pair_num,
		estimate_coeff_for_tv_pair_num, estimate_coeff_for_ee_pair_num);

	initialTriangleHash();


	ori_hash_value.resize(thread_num);
}


void SpatialHashing::initialOriHashValue()
{
	int total_triangle_num = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		total_triangle_num += cloth->data()[i].mesh_struct.triangle_indices.size();
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		total_triangle_num += tetrahedron->data()[i].mesh_struct.triangle_indices.size();
	}

	for (int i = 0; i < collider->size(); ++i) {
		total_triangle_num += collider->data()[i].mesh_struct.triangle_indices.size();
	}
	for (unsigned int i = 0; i < thread_num; ++i) {
		ori_hash_value.reserve(total_obj_num);
	}


}


void SpatialHashing::reorganzieDataOfObjects()
{
	obj_tri_aabb.resize(collider_begin_obj_index);
	obj_vertex_aabb.resize(collider_begin_obj_index);
	obj_edge_aabb.resize(collider_begin_obj_index);
	obj_triangle_index_begin_per_thread.resize(collider_begin_obj_index);
	obj_vertex_index_begin_per_thread.resize(collider_begin_obj_index);
	obj_edge_index_begin_per_thread.resize(collider_begin_obj_index);
	triangle_vertex_index.resize(collider_begin_obj_index);
	edge_vertex_index.resize(collider_begin_obj_index);

	vertices.resize(collider_begin_obj_index);
	edges.resize(collider_begin_obj_index);
	mesh_struct.resize(collider_begin_obj_index);
	
	representative_edge_num.resize(collider_begin_obj_index);
	face_edge.resize(collider_begin_obj_index);


	tetrahedron_vertex_index_on_surface.resize(tetrahedron->size());

	for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
		if (i < cloth->size()) {
			obj_tri_aabb[i] = cloth->data()[i].triangle_AABB.data();
			obj_triangle_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.face_index_begin_per_thread.data();

			obj_vertex_aabb[i] = cloth->data()[i].vertex_AABB.data();
			obj_vertex_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();

			obj_edge_aabb[i] = cloth->data()[i].edge_AABB.data();
			obj_edge_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.edge_index_begin_per_thread.data();

			triangle_vertex_index[i] = cloth->data()[i].mesh_struct.triangle_indices.data();
			edge_vertex_index[i] = cloth->data()[i].mesh_struct.edge_vertices.data();

			vertices[i] = cloth->data()[i].mesh_struct.vertices.data();
			edges[i] = cloth->data()[i].mesh_struct.edges.data();

			representative_edge_num[i] = cloth->data()[i].representative_edge_num.data();
			face_edge[i] = cloth->data()[i].mesh_struct.face_edges.data();

			mesh_struct[i] = &cloth->data()[i].mesh_struct;
		}
		else {
			obj_tri_aabb[i] = tetrahedron->data()[i - cloth->size()].triangle_AABB.data();
			obj_triangle_index_begin_per_thread[i] = tetrahedron->data()[i - cloth->size()].mesh_struct.face_index_begin_per_thread.data();

			obj_vertex_aabb[i] = tetrahedron->data()[i - cloth->size()].vertex_AABB.data();
			obj_vertex_index_begin_per_thread[i] = tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_index_on_surface_begin_per_thread.data();

			obj_edge_aabb[i] = tetrahedron->data()[i - cloth->size()].edge_AABB.data();
			obj_edge_index_begin_per_thread[i] = tetrahedron->data()[i - cloth->size()].mesh_struct.edge_index_begin_per_thread.data();

			//std::cout << obj_triangle_index_begin_per_thread[i][0] << " " << obj_triangle_index_begin_per_thread[i][thread_num] << std::endl;

			triangle_vertex_index[i] = tetrahedron->data()[i - cloth->size()].mesh_struct.triangle_indices.data();
			edge_vertex_index[i] = tetrahedron->data()[i - cloth->size()].mesh_struct.edge_vertices.data();

			vertices[i] = tetrahedron->data()[i - cloth->size()].mesh_struct.vertices.data();
			edges[i] = tetrahedron->data()[i - cloth->size()].mesh_struct.edges.data();

			representative_edge_num[i] = tetrahedron->data()[i - cloth->size()].representative_edge_num.data();
			face_edge[i] = tetrahedron->data()[i - cloth->size()].mesh_struct.face_edges.data();
			mesh_struct[i] = &tetrahedron->data()[i - cloth->size()].mesh_struct;
		}
	}





	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		tetrahedron_vertex_index_on_surface[i] = tetrahedron->data()[i].mesh_struct.vertex_index_on_sureface.data();
	}


	if (has_collider) {
		collider_tri_aabb.resize(collider->size());
		collider_vertex_aabb.resize(collider->size());
		collider_edge_aabb.resize(collider->size());
		triangle_vertex_index_collider.resize(collider->size());
		edge_vertex_index_collider.resize(collider->size());
		collider_vertex_index_begin_per_thread.resize(collider->size());
		collider_edge_index_begin_per_thread.resize(collider->size());
		representative_edge_num_collider.resize(collider->size());
		face_edge_collider.resize(collider->size());
		for (unsigned int i = 0; i < collider->size(); ++i) {
			collider_tri_aabb[i] = collider->data()[i].triangle_AABB.data();
			collider_vertex_aabb[i] = collider->data()[i].vertex_AABB.data();
			collider_edge_aabb[i] = collider->data()[i].edge_AABB.data();
			triangle_vertex_index_collider[i] = collider->data()[i].mesh_struct.triangle_indices.data();
			edge_vertex_index_collider[i] = collider->data()[i].mesh_struct.edge_vertices.data();
			collider_vertex_index_begin_per_thread[i] = collider->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
			collider_edge_index_begin_per_thread[i] = collider->data()[i].mesh_struct.edge_index_begin_per_thread.data();
			representative_edge_num_collider[i] = collider->data()[i].representative_edge_num.data();
			face_edge_collider[i] = collider->data()[i].mesh_struct.face_edges.data();
		}
	}
}


void SpatialHashing::deleteArray()
{
	//delete[] spatial_hashing_value;
	////delete[] spatial_hashing_obj_index;
	//delete[] spatial_hashing_triangle_index;
	//delete[] spatial_hashing_value_collider;
	////delete[] spatial_hashing_obj_index_collider;
	//delete[] spatial_hashing_triangle_index_collider;

	//delete[] radix_sort;
	//delete[] radix_sort_collider;


	//delete[] vertex_triangle_pair;
}


void SpatialHashing::initialTriangleHash()
{
	max_triangle_cell_size = 28;
	max_vertex_cell_size = max_triangle_cell_size>>1;
	obj_edge_hash = new unsigned int* [cloth->size() + tetrahedron->size()];
	for (int i = 0; i < cloth->size(); ++i) {
		obj_edge_hash[i] = new unsigned int[max_triangle_cell_size * (*cloth)[i].mesh_struct.edges.size()];
		memset(obj_edge_hash[i], 0, 4 * max_triangle_cell_size * (*cloth)[i].mesh_struct.edges.size());
	}
	int obj_no;
	for (int i = 0; i < tetrahedron->size(); ++i) {
		obj_no = i + tetrahedron_begin_obj_index;
		obj_edge_hash[obj_no] = new unsigned int[max_triangle_cell_size * (*tetrahedron)[i].mesh_struct.edges.size()];
		memset(obj_edge_hash[obj_no], 0, 4 * max_triangle_cell_size * (*tetrahedron)[i].mesh_struct.edges.size());
	}

	

	obj_vertex_hash = new unsigned int* [cloth->size() + tetrahedron->size()];
	for (int i = 0; i < cloth->size(); ++i) {
		obj_vertex_hash[i] = new unsigned int[max_vertex_cell_size * (*cloth)[i].mesh_struct.vertex_position.size()];
		memset(obj_vertex_hash[i], 0, 4 * max_vertex_cell_size * (*cloth)[i].mesh_struct.vertex_position.size());
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		obj_no = i + tetrahedron_begin_obj_index;
		obj_vertex_hash[obj_no] = new unsigned int[max_vertex_cell_size * (*tetrahedron)[i].mesh_struct.vertex_position.size()];
		memset(obj_vertex_hash[obj_no], 0, 4 * max_vertex_cell_size * (*tetrahedron)[i].mesh_struct.vertex_position.size());
	}

	obj_is_used = new bool** [thread_num];
	unsigned int obj_num = std::max(collider_begin_obj_index, (unsigned int)collider->size());
	unsigned int index_size;

	for (int i = 0; i < thread_num; ++i) {
		obj_is_used[i] = new bool* [obj_num];
		for (int j = 0; j < cloth->size(); ++j) {
			index_size = std::max((*cloth)[j].mesh_struct.triangle_indices.size(), (*cloth)[j].mesh_struct.edges.size());
			if (j < collider->size()) {
				if (index_size < collider->data()[j].mesh_struct.triangle_indices.size()) {
					index_size = collider->data()[j].mesh_struct.triangle_indices.size();
				}
				if (index_size < collider->data()[j].mesh_struct.edges.size()) {
					index_size = collider->data()[j].mesh_struct.edges.size();
				}
			}
			obj_is_used[i][j] = new bool[index_size];
			memset(obj_is_used[i][j], 0, index_size);
		}
		for (int j = 0; j < tetrahedron->size(); ++j) {
			index_size = std::max((*tetrahedron)[j].mesh_struct.triangle_indices.size(), (*tetrahedron)[j].mesh_struct.edges.size());
			if (j + tetrahedron_begin_obj_index < collider->size()) {
				if (index_size < collider->data()[j + tetrahedron_begin_obj_index].mesh_struct.triangle_indices.size()) {
					index_size = collider->data()[j + tetrahedron_begin_obj_index].mesh_struct.triangle_indices.size();
				}
				if (index_size < collider->data()[j + tetrahedron_begin_obj_index].mesh_struct.edges.size()) {
					index_size = collider->data()[j + tetrahedron_begin_obj_index].mesh_struct.edges.size();
				}
			}
			obj_is_used[i][j + tetrahedron_begin_obj_index] = new bool[index_size];
			memset(obj_is_used[i][j + tetrahedron_begin_obj_index], 0, index_size);
		}
		if (collider->size()>collider_begin_obj_index) {
			for (unsigned int j = collider_begin_obj_index; j < collider->size(); ++j) {
				obj_is_used[i][j] = new bool[std::max((*collider)[j].mesh_struct.triangle_indices.size(), (*collider)[j].mesh_struct.edges.size())];
				memset(obj_is_used[i][j], 0, std::max((*collider)[j].mesh_struct.triangle_indices.size(), (*collider)[j].mesh_struct.edges.size()));
			}
		}
	}

	if (has_collider) {
		collider_vertex_hash = new unsigned int* [collider->size()];
		for (int i = 0; i < collider->size(); ++i) {
			collider_vertex_hash[i] = new unsigned int[max_vertex_cell_size * (*collider)[i].mesh_struct.vertex_position.size()];
			memset(collider_vertex_hash[i], 0, 4*max_vertex_cell_size * (*collider)[i].mesh_struct.vertex_position.size());
		}
	}
}



void SpatialHashing::buildSpatialHashingForOri(double* scene_aabb)
{
	memcpy(this->scene_aabb, scene_aabb, 48);
	setSpatialHashingInitialCount();
	thread->assignTask(this, ORI_TRIANGLE_HASHING);
}


void SpatialHashing::buildSpatialHashing(double* scene_aabb)
{
	memcpy(this->scene_aabb, scene_aabb, 48);
	time_t t = clock();
	time_t t1 = clock();

	//std::cout<<"spatial hashing aabb "<< scene_aabb[0] << " " << scene_aabb[1] << " " << scene_aabb[2] << " " << scene_aabb[3] << " " << scene_aabb[4] << " " << scene_aabb[5] << std::endl;

	//t = clock();
	thread->assignTask(this, TRIANGLE_HASHING_SMALLER_HASH_TABLE);

	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	triangleHashingSmallerHashTable(i);
	//}
	if (searchPairByCell) {
		recordNonEmptyCell();
			//use this for loop by cell
		thread->assignTask(this, FIND_ALL_PAIRS_HASH_TABLE);
	}
	else {
		thread->assignTask(this, COMBINE_HASH_TABLE);
		thread->assignTask(this, FIND_ALL_PAIRS_HASH_TABLE_BY_ELEMENT);

		if (record_pair_by_element) {
			findAllTriangleVertexPairByTriangle();
		}
	}
		//for (unsigned int i = 0; i < thread_num; ++i) {
		//	findAllEdgeEdgePairs(i);
		//}
		

		//use this for loop by element
		//thread->assignTask(this, FIND_ALL_TRIANGLE_PAIRS_HASH_TABLE_ELEMENTWISE);
	//}
//t1 = clock();
	//std::cout << "find all triangle pairs multi thread " << t1 - t << std::endl;

}


// FIND_ALL_TRIANGLE_PAIRS_HASH_TABLE_ELEMENTWISE, use this for loop by element
void SpatialHashing::findAllPairsHashTableElementwise(int thread_No)
{
	unsigned int* primitive_pair_ = edge_edge_pair[thread_No] + 1;
	unsigned int* primitive_pair_collider_ = edge_edge_pair[thread_No] + 1;// = edge_edge_pair_collider[thread_No] + 1;
	unsigned int* spatial_hashing_cell_= spatial_hashing_cell_edge[0];
	unsigned int* spatial_hashing_cell_size_= spatial_hashing_cell_edge_size[0];
	unsigned int* spatial_hashing_cell_collider_ = spatial_hashing_cell_collider_edge[0];
	unsigned int* spatial_hashing_cell_size_collider_ = spatial_hashing_cell_collider_edge_size[0];

	unsigned int max_index_number_in_one_cell_= max_index_number_in_one_cell_edge;
	unsigned int max_index_number_in_one_cell_collider_ = max_index_number_in_one_cell_collider_edge;

	for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {		
		for (unsigned int j = obj_edge_index_begin_per_thread[i][thread_No]; j < obj_edge_index_begin_per_thread[i][thread_No + 1]; ++j) {
			searchPrimitive(obj_edge_aabb[i][j].data(), obj_edge_hash[i] + j * max_triangle_cell_size + 1, obj_edge_hash[i][j * max_triangle_cell_size],
				i, j, primitive_pair_, primitive_pair_collider_, thread_No, spatial_hashing_cell_, spatial_hashing_cell_size_, spatial_hashing_cell_collider_,
				spatial_hashing_cell_size_collider_,true, &vertices[i][edge_vertex_index[i][j<<1]].edge, &vertices[i][edge_vertex_index[i][(j << 1)+1]].edge,
				max_index_number_in_one_cell_, max_index_number_in_one_cell_collider_);
		}
	}
	edge_edge_pair[thread_No][0] = primitive_pair_ - edge_edge_pair[thread_No] - 1;
	//edge_edge_pair_collider[thread_No][0] = primitive_pair_collider_ - edge_edge_pair_collider[thread_No] - 1;

	primitive_pair_ = vertex_triangle_pair[thread_No] + 1;
	primitive_pair_collider_ = vertex_obj_triangle_collider_pair[thread_No] + 1;
	spatial_hashing_cell_ = spatial_hashing_cell_triangle[0];
	spatial_hashing_cell_size_ = spatial_hashing_cell_triangle_size[0];

	spatial_hashing_cell_collider_ = spatial_hashing_cell_collider_triangle[0];
	spatial_hashing_cell_size_collider_ = spatial_hashing_cell_collider_triangle_size[0];


	max_index_number_in_one_cell_ = max_index_number_in_one_cell;
	max_index_number_in_one_cell_collider_ = max_index_number_in_one_cell_collider;

	unsigned int* vertex_index_on_surface;
	unsigned int cloth_size = cloth->size();
	unsigned int vertex_index;
	for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
		if (i < cloth_size) {
			for (unsigned int j = obj_vertex_index_begin_per_thread[i][thread_No]; j < obj_vertex_index_begin_per_thread[i][thread_No + 1]; ++j) {
				searchPrimitive(obj_vertex_aabb[i][j].data(), obj_vertex_hash[i] + j * max_vertex_cell_size + 1, obj_vertex_hash[i][j * max_vertex_cell_size],
					i, j, primitive_pair_, primitive_pair_collider_, thread_No, spatial_hashing_cell_, spatial_hashing_cell_size_, spatial_hashing_cell_collider_,
					spatial_hashing_cell_size_collider_, false, &vertices[i][j].face, &vertices[i][j].face,
					max_index_number_in_one_cell_, max_index_number_in_one_cell_collider_);
			}
		}
		else {
			vertex_index_on_surface = tetrahedron_vertex_index_on_surface[i - cloth_size];
			for (unsigned int j = obj_vertex_index_begin_per_thread[i][thread_No]; j < obj_vertex_index_begin_per_thread[i][thread_No + 1]; ++j) {
				vertex_index = vertex_index_on_surface[j];
				searchPrimitive(obj_vertex_aabb[i][vertex_index].data(), obj_vertex_hash[i] + vertex_index * max_vertex_cell_size + 1, obj_vertex_hash[i][vertex_index * max_vertex_cell_size],
					i, vertex_index, primitive_pair_, primitive_pair_collider_, thread_No, spatial_hashing_cell_, spatial_hashing_cell_size_, spatial_hashing_cell_collider_,
					spatial_hashing_cell_size_collider_, false, &vertices[i][vertex_index].face, &vertices[i][vertex_index].face,
					max_index_number_in_one_cell_, max_index_number_in_one_cell_collider_);
			}

		}

		
	}
	vertex_triangle_pair[thread_No][0]= primitive_pair_ - vertex_triangle_pair[thread_No] - 1;
	vertex_obj_triangle_collider_pair[thread_No][0]= primitive_pair_collider_ - vertex_obj_triangle_collider_pair[thread_No] - 1;

	//if (has_collider) {
	//	primitive_pair_collider_ = vertex_collider_triangle_obj_pair[thread_No] + 1;
	//	spatial_hashing_cell_ = spatial_hashing_cell_triangle[0];
	//	spatial_hashing_cell_size_ = spatial_hashing_cell_triangle_size[0];
	//	for (unsigned int i = 0; i < collider->size(); ++i) {
	//		for (unsigned int j = collider_vertex_index_begin_per_thread[i][thread_No]; j < collider_vertex_index_begin_per_thread[i][thread_No + 1]; ++j) {
	//			searchPrimitiveOnCollider(collider_vertex_aabb[i][j].data(), collider_vertex_hash[i] + j * max_vertex_cell_size + 1, collider_vertex_hash[i][j * max_vertex_cell_size],
	//				i, j, primitive_pair_collider_, thread_No, spatial_hashing_cell_, spatial_hashing_cell_size_, max_index_number_in_one_cell_);
	//		}
	//	}
	//}
}



void SpatialHashing::findAllEdgeEdgePairs(int thread_No)
{
	//triangle_pair_number_thread[thread_No] = 0;
	unsigned int last_cell_index = global_cell_start_edge[(thread_No << 1) + 1];
	//unsigned int last_cell_index = non_empty_cell_index_begin_per_thread[thread_No + 1];
	unsigned int start = global_cell_start_edge[thread_No << 1];
	//unsigned int start= non_empty_cell_index_begin_per_thread[thread_No];

	unsigned int* cell_edge_index;
	unsigned int* cell_collider_edge_index;

	unsigned int* primitive_pair_;
	//unsigned int* edge_edge_collider_pair_;

	primitive_pair_ = edge_edge_pair[thread_No] + 1;
	//edge_edge_collider_pair_ = edge_edge_pair_collider[thread_No] + 1;

	unsigned int max_index_number_in_one_cell_ = max_index_number_in_one_cell_edge;
	unsigned int max_index_number_in_one_cell_collider_ = max_index_number_in_one_cell_collider_edge;

	double* aabb_1;
	bool has_collider_ = has_collider;
	unsigned int obj_index_0, primitive_index_0;// , obj_index_1, triangle_index_1;

	unsigned int hash_cell_element_size;
	unsigned int hash_cell_vertex_size;
	unsigned int hash_cell_element_size_collider;
	unsigned int hash_cell_vertex_size_collider;

	for (unsigned int cell_index = start; cell_index < last_cell_index; ++cell_index) {
		//cell_index = non_empty_cell_index_[kk];
		hash_cell_element_size = spatial_hashing_cell_edge_size[0][cell_index];
		cell_edge_index = spatial_hashing_cell_edge[0]+cell_index* max_index_number_in_one_cell_;
		if (hash_cell_element_size > 2) {
			for (unsigned int i = 0; i < hash_cell_element_size; i += 2) {
				obj_index_0 = cell_edge_index[i + 1]; primitive_index_0 = cell_edge_index[i];
				aabb_1 = obj_edge_aabb[obj_index_0][primitive_index_0].data();
				for (unsigned int j = i + 2; j < hash_cell_element_size; j += 2) {
					if (AABB::AABB_intersection(aabb_1, obj_edge_aabb[cell_edge_index[j+1]][cell_edge_index[j]].data())) {
						if (obj_index_0 == cell_edge_index[j + 1]) {
							if (!edgeEdgeconnected(edge_vertex_index[obj_index_0] + (primitive_index_0 << 1),
								edge_vertex_index[obj_index_0] + (cell_edge_index[j] << 1))) {
								memcpy(primitive_pair_, cell_edge_index + i, 8);
								primitive_pair_ += 2;
								memcpy(primitive_pair_, cell_edge_index + j, 8);
								primitive_pair_ += 2;
							}
						}
						else {
							memcpy(primitive_pair_, cell_edge_index + i, 8);
							primitive_pair_ += 2;
							memcpy(primitive_pair_, cell_edge_index + j, 8);
							primitive_pair_ += 2;
						}
					}
				}
			}
		}
		//if (has_collider_) {
		//	hash_cell_element_size_collider = spatial_hashing_cell_collider_edge_size[0][cell_index];
		//	cell_collider_edge_index = spatial_hashing_cell_collider_edge[0]+cell_index* max_index_number_in_one_cell_collider_;
		//	if (hash_cell_element_size > 0) {
		//		for (unsigned int i = 0; i < hash_cell_element_size; i += 2) {
		//			obj_index_0 = cell_edge_index[i + 1]; primitive_index_0 = cell_edge_index[i];
		//			aabb_1 = obj_edge_aabb[obj_index_0][primitive_index_0].data();
		//			for (unsigned int j = 0; j < hash_cell_element_size_collider; j += 2) {
		//				if (AABB::AABB_intersection(aabb_1, 
		//					collider_edge_aabb[cell_collider_edge_index[j+1]][cell_collider_edge_index[j]].data())) {
		//					memcpy(edge_edge_collider_pair_, cell_edge_index + i, 8);
		//					edge_edge_collider_pair_ += 2;
		//					memcpy(edge_edge_collider_pair_, cell_collider_edge_index + j, 8);
		//					edge_edge_collider_pair_ += 2;
		//				}
		//			}
		//		}
		//	}
		//}
	}

	edge_edge_pair[thread_No][0] = primitive_pair_ - edge_edge_pair[thread_No] - 1;
	//edge_edge_pair_collider[thread_No][0] = edge_edge_collider_pair_ - edge_edge_pair_collider[thread_No] - 1;
}


void SpatialHashing::findAllVertexTrianglePairsByPrimitive(int thread_No)
{
	unsigned int* primitive_pair_;

	//obj vertex obj triangle
	primitive_pair_ = vertex_triangle_pair[thread_No] + 1;
	for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
		findAllVertexTrianglePairsByPrimitiveSingleObj(thread_No, i, primitive_pair_, obj_vertex_index_begin_per_thread[i][thread_No],
			obj_vertex_index_begin_per_thread[i][thread_No + 1], obj_vertex_aabb[i], max_index_number_in_one_cell, spatial_hashing_cell_triangle[0],
			spatial_hashing_cell_triangle_size[0], obj_tri_aabb,true,i>=cloth->size());
	}
	vertex_triangle_pair[thread_No][0] = primitive_pair_ - vertex_triangle_pair[thread_No] - 1;

	if (has_collider) {
		//obj vertex collider triangle
		primitive_pair_ = vertex_obj_triangle_collider_pair[thread_No] + 1;
		for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
			findAllVertexTrianglePairsByPrimitiveSingleObj(thread_No, i, primitive_pair_, obj_vertex_index_begin_per_thread[i][thread_No],
				obj_vertex_index_begin_per_thread[i][thread_No + 1], obj_vertex_aabb[i], max_index_number_in_one_cell_collider, spatial_hashing_cell_collider_triangle[0],
				spatial_hashing_cell_collider_triangle_size[0], collider_tri_aabb, false, i >= cloth->size());
		}
		vertex_obj_triangle_collider_pair[thread_No][0] = primitive_pair_ - vertex_obj_triangle_collider_pair[thread_No] - 1;

		//collider vertex obj triangle
		//primitive_pair_ = vertex_collider_triangle_obj_pair[thread_No] + 1;
		//for (unsigned int i = 0; i < collider->size(); ++i) {
		//	findAllVertexTrianglePairsByPrimitiveSingleObj(thread_No, i, primitive_pair_, collider_vertex_index_begin_per_thread[i][thread_No],
		//		collider_vertex_index_begin_per_thread[i][thread_No + 1], collider_vertex_aabb[i], max_index_number_in_one_cell, spatial_hashing_cell_triangle[0],
		//		spatial_hashing_cell_triangle_size[0], obj_tri_aabb, false,false);
		//}
		//vertex_collider_triangle_obj_pair[thread_No][0] = primitive_pair_ - vertex_collider_triangle_obj_pair[thread_No] - 1;
	}
}



void SpatialHashing::initialPairByElement()
{
	for (unsigned int i = 0; i < collider_begin_obj_index;++i){
		memset(triangle_vertex_pair_num_record[i], 0, 4 * mesh_struct[i]->triangle_indices.size());
	}
}



void SpatialHashing::findAllTriangleVertexPairByTriangle()
{
	initialPairByElement();
	for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
		findAllTriangleVertexPairByTriangleSingleObj(i, vertex_triangle_pair_by_vertex[i], vertex_triangle_pair_num_record[i], triangle_vertex_pair_by_triangle,
			triangle_vertex_pair_num_record, total_length_every_element_for_vertex_triange, total_length_every_element_for_triangle_vertex);
	}

}



void SpatialHashing::findAllVertexTrianglePairsByPrimitiveByVertex(int thread_No)
{
	unsigned int* primitive_pair_;
	unsigned int* primitive_pair_num_record;

	//obj vertex obj triangle
	for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
		primitive_pair_ = vertex_triangle_pair_by_vertex[i];
		primitive_pair_num_record = vertex_triangle_pair_num_record[i];
		findAllVertexTrianglePairsByPrimitiveSingleObj_ByVertex(thread_No, i, primitive_pair_, primitive_pair_num_record, obj_vertex_index_begin_per_thread[i][thread_No],
			obj_vertex_index_begin_per_thread[i][thread_No + 1], obj_vertex_aabb[i], max_index_number_in_one_cell, spatial_hashing_cell_triangle[0],
			spatial_hashing_cell_triangle_size[0], obj_tri_aabb, true, i >= cloth->size());
	}
	if (has_collider) {
		for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
			primitive_pair_ = vertex_obj_triangle_collider_pair_by_vertex[i];
			primitive_pair_num_record = vertex_obj_triangle_collider_num_record[i];
			findAllVertexTrianglePairsByPrimitiveSingleObj_ByVertex(thread_No, i, primitive_pair_, primitive_pair_num_record, obj_vertex_index_begin_per_thread[i][thread_No],
				obj_vertex_index_begin_per_thread[i][thread_No + 1], obj_vertex_aabb[i], max_index_number_in_one_cell_collider, spatial_hashing_cell_collider_triangle[0],
				spatial_hashing_cell_collider_triangle_size[0], collider_tri_aabb, false, i >= cloth->size());
		}
	}
}



void SpatialHashing::findAllEdgeEdgePairsByPrimitive(int thread_No)
{
	unsigned int* primitive_pair_;
	//obj edge obj edge
	primitive_pair_ = edge_edge_pair[thread_No] + 1;
	for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
		findAllEdgeEdgePairsByPrimitiveSingleObj(thread_No, i, primitive_pair_,
			obj_edge_aabb, max_index_number_in_one_cell_edge, spatial_hashing_cell_edge[0],
			spatial_hashing_cell_edge_size[0], true);
	}
	edge_edge_pair[thread_No][0] = primitive_pair_ - edge_edge_pair[thread_No] - 1;

	//obj edge collider edge
	if (has_collider) {
		//primitive_pair_ = edge_edge_pair_collider[thread_No] + 1;
		//for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
		//	findAllEdgeEdgePairsByPrimitiveSingleObj(thread_No, i, primitive_pair_, representative_edge_num_collider,
		//		face_edge_collider, collider_edge_aabb, max_index_number_in_one_cell_collider, spatial_hashing_cell_collider_triangle[0],
		//		spatial_hashing_cell_collider_triangle_size[0], false);
		//}
		//edge_edge_pair_collider[thread_No][0] = primitive_pair_ - edge_edge_pair_collider[thread_No] - 1;
	}
}

void SpatialHashing::findAllEdgeEdgePairsByPrimitiveByEdge(int thread_No)
{
	unsigned int* primitive_pair_;
	unsigned int* primitive_pair_num_record_;
	for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
		primitive_pair_ = edge_edge_pair_by_edge[i];
		primitive_pair_num_record_ = edge_edge_pair_num_record[i];
		findAllEdgeEdgePairsByPrimitiveSingleObjByEdge(thread_No, i, primitive_pair_, primitive_pair_num_record_,
			obj_edge_aabb, max_index_number_in_one_cell_edge, spatial_hashing_cell_edge[0],
			spatial_hashing_cell_edge_size[0], true);
	}
}


void SpatialHashing::findAllVertexTrianglePairs(int thread_No)
{
	//triangle_pair_number_thread[thread_No] = 0;
	unsigned int last_cell_index = global_cell_start_vertex_triangle[(thread_No << 1) + 1];
	//unsigned int last_cell_index = non_empty_cell_index_begin_per_thread[thread_No + 1];
	unsigned int start = global_cell_start_vertex_triangle[thread_No << 1];
	//unsigned int start= non_empty_cell_index_begin_per_thread[thread_No];

	unsigned int* cell_triangle_index;
	unsigned int* cell_collider_triangle_index;

	unsigned int* cell_vertex_index;
	unsigned int* cell_collider_vertex_index;

	unsigned int* primitive_pair_;
	//unsigned int* vertex_collider_triangle_pair_;
	unsigned int* vertex_triangle_collider_pair_;

	primitive_pair_ = vertex_triangle_pair[thread_No] + 1;
	//vertex_collider_triangle_pair_ = vertex_collider_triangle_obj_pair[thread_No] + 1;
	vertex_triangle_collider_pair_ = vertex_obj_triangle_collider_pair[thread_No] + 1;
	double* aabb_1;

	unsigned int obj_index_0, primitive_index_0;// , obj_index_1, triangle_index_1;
	//unsigned int cell_index;
	bool has_collider_ = has_collider;

	//unsigned int a = 0;
	unsigned int hash_cell_element_size;
	unsigned int hash_cell_vertex_size;
	unsigned int hash_cell_element_size_collider;
	unsigned int hash_cell_vertex_size_collider;

	//vertex-triangle
	
	unsigned int max_index_number_in_one_cell_triangle_ = max_index_number_in_one_cell;
	unsigned int max_index_number_in_one_cell_collider_triangle_ = max_index_number_in_one_cell_collider;

	unsigned int max_index_number_in_one_cell_vertex_ = max_index_number_in_one_cell_vertex;
	unsigned int max_index_number_in_one_cell_collider_vertex_ = max_index_number_in_one_cell_collider_vertex;


	for (unsigned int cell_index = start; cell_index < last_cell_index; ++cell_index) {
		//cell_index = non_empty_cell_index_[kk];
		hash_cell_element_size = spatial_hashing_cell_triangle_size[0][cell_index];
		cell_triangle_index = spatial_hashing_cell_triangle[0]+cell_index* max_index_number_in_one_cell_triangle_;

		hash_cell_vertex_size = spatial_hashing_cell_vertex_size[0][cell_index];
		cell_vertex_index = spatial_hashing_cell_vertex[0]+cell_index* max_index_number_in_one_cell_vertex_;

		for (unsigned int i = 0; i < hash_cell_vertex_size; i += 2) {
			obj_index_0 = cell_vertex_index[i + 1]; primitive_index_0 = cell_vertex_index[i];
			aabb_1 = obj_vertex_aabb[obj_index_0][primitive_index_0].data();
			for (unsigned int j = 0; j < hash_cell_element_size; j += 2) {
				if (AABB::AABB_intersection(aabb_1, obj_tri_aabb[cell_triangle_index[j + 1]][cell_triangle_index[j]].data())) {
					if (obj_index_0 == cell_triangle_index[j + 1]) {
						//if (primitive_index_0 == 491) {
						//	std::cout << cell_triangle_index[j] << std::endl;
						//}

						if (!vertexInTriangle(triangle_vertex_index[cell_triangle_index[j + 1]][cell_triangle_index[j]].data(),
							primitive_index_0)) {
							memcpy(primitive_pair_, cell_vertex_index + i, 8);
							primitive_pair_ += 2;
							memcpy(primitive_pair_, cell_triangle_index + j, 8);
							primitive_pair_ += 2;
						}
					}
					else {
						memcpy(primitive_pair_, cell_vertex_index + i, 8);
						primitive_pair_ += 2;
						memcpy(primitive_pair_, cell_triangle_index + j, 8);
						primitive_pair_ += 2;
					}
				}
			}
		}
		if (has_collider_) {
			//std::cout << "--" << std::endl;
			hash_cell_element_size_collider = spatial_hashing_cell_collider_triangle_size[0][cell_index];
			cell_collider_triangle_index = spatial_hashing_cell_collider_triangle[0]+cell_index* max_index_number_in_one_cell_collider_triangle_;

			hash_cell_vertex_size_collider = spatial_hashing_cell_collider_vertex_size[0][cell_index];
			cell_collider_vertex_index = spatial_hashing_cell_collider_vertex[0] +cell_index* max_index_number_in_one_cell_collider_vertex_;

			//if (hash_cell_vertex_size_collider) {
			//	for (unsigned int j = 0; j < hash_cell_vertex_size_collider; j += 2) {
			//		obj_index_0 = cell_collider_vertex_index[j + 1]; primitive_index_0 = cell_collider_vertex_index[j];
			//		aabb_1 = collider_vertex_aabb[obj_index_0][primitive_index_0].data();
			//		for (unsigned int i = 0; i < hash_cell_element_size; i += 2) {
			//			if (AABB::AABB_intersection(aabb_1, obj_tri_aabb[cell_triangle_index[i + 1]][cell_triangle_index[i]].data())) {
			//				memcpy(vertex_collider_triangle_pair_, cell_collider_vertex_index + j, 8);
			//				vertex_collider_triangle_pair_ += 2;
			//				memcpy(vertex_collider_triangle_pair_, cell_triangle_index + i, 8);
			//				vertex_collider_triangle_pair_ += 2;
			//			}
			//		}
			//	}
			//}

			if (hash_cell_element_size_collider) {
				for (unsigned int i = 0; i < hash_cell_vertex_size; i += 2) {
					obj_index_0 = cell_vertex_index[i + 1]; primitive_index_0 = cell_vertex_index[i];
					aabb_1 = obj_vertex_aabb[obj_index_0][primitive_index_0].data();
					for (unsigned int j = 0; j < hash_cell_element_size_collider; j += 2) {
						if (AABB::AABB_intersection(aabb_1, collider_tri_aabb[cell_collider_triangle_index[j + 1]][cell_collider_triangle_index[j]].data())) {
							memcpy(vertex_triangle_collider_pair_, cell_vertex_index + i, 8);
							vertex_triangle_collider_pair_ += 2;
							memcpy(vertex_triangle_collider_pair_, cell_collider_triangle_index + j, 8);
							vertex_triangle_collider_pair_ += 2;

						}
					}
				}
			}

		}
	}

	// triangle_pair_number_thread[thread_No] = a;
	vertex_triangle_pair[thread_No][0] = primitive_pair_ - vertex_triangle_pair[thread_No] - 1;
	vertex_obj_triangle_collider_pair[thread_No][0] = vertex_triangle_collider_pair_ - vertex_obj_triangle_collider_pair[thread_No] - 1;
	//vertex_collider_triangle_obj_pair[thread_No][0] = vertex_collider_triangle_pair_ - vertex_collider_triangle_obj_pair[thread_No] - 1;

}


//FIND_ALL_PAIRS_HASH_TABLE_BY_ELEMENT
void SpatialHashing::findAllPairsByPrimitive(int thread_No)
{
	if (record_pair_by_element) {
		findAllVertexTrianglePairsByPrimitiveByVertex(thread_No);
		findAllEdgeEdgePairsByPrimitiveByEdge(thread_No);
	}
	else {
		findAllVertexTrianglePairsByPrimitive(thread_No);
		findAllEdgeEdgePairsByPrimitive(thread_No);
	}
}






//FIND_ALL_PAIRS_HASH_TABLE,  use this for loop by cell
void SpatialHashing::findAllPairsHashTable(int thread_No)
{
	findAllVertexTrianglePairs(thread_No);
	findAllEdgeEdgePairs(thread_No);
}




void SpatialHashing::searchPrimitiveOnCollider(double* aabb, unsigned int* hash_index, unsigned int hash_cell_size, 
	unsigned int input_obj_No, unsigned int vertex_index,
	unsigned int*& triangle_pair_with_collider_, unsigned int thread_No,
	unsigned int* spatial_hashing_cell_, unsigned int* spatial_hash_size_, unsigned int max_index_number_in_one_cell_collider)
{
	unsigned int* primitive_index_record_ = primitive_index_record[thread_No] + 1;

	bool** obj_is_used_ = obj_is_used[thread_No];
	unsigned int hash_cell_index_end;
	unsigned int obj_No; unsigned int current_triangle_index;
	unsigned int obj_num = cloth->size() + tetrahedron->size();

	for (unsigned int i = 0; i < hash_cell_size; ++i) {
		for (unsigned int j = 0; j < spatial_hash_size_[hash_index[i]]; j += 2) {
			obj_No = spatial_hashing_cell_[hash_index[i]* max_index_number_in_one_cell_collider + j + 1];
			current_triangle_index = spatial_hashing_cell_[hash_index[i]* max_index_number_in_one_cell_collider+j];
			if (!obj_is_used_[obj_No][current_triangle_index]) {
				obj_is_used_[obj_No][current_triangle_index] = true;
				memcpy(primitive_index_record_, spatial_hashing_cell_+hash_index[i] * max_index_number_in_one_cell_collider + j, 8);
				primitive_index_record_ += 2;
				if (AABB::AABB_intersection(aabb, obj_tri_aabb[obj_No][current_triangle_index].data())) {
					*(triangle_pair_with_collider_++) = vertex_index;
					*(triangle_pair_with_collider_++) = input_obj_No;
					*(triangle_pair_with_collider_++) = current_triangle_index;
					*(triangle_pair_with_collider_++) = obj_No;
				}
			}
		}
	}
	primitive_index_record[thread_No][0] = primitive_index_record_ - primitive_index_record[thread_No] - 1;
	primitive_index_record_ = primitive_index_record[thread_No];
	for (unsigned int i = 0; i < primitive_index_record_[0]; i += 2) {
		obj_is_used_[primitive_index_record_[i + 2]][primitive_index_record_[i + 1]] = false;
	}
}

void SpatialHashing::findAllEdgeEdgePairsByPrimitiveSingleObj(int thread_No, int obj_No, unsigned int*& primitive_pair_,
	std::vector<std::array<double, 6>*>& obj_edge_aabb_,
	unsigned int max_index_number_in_one_cell_edge_, unsigned int* spatial_hashing_cell_edge, unsigned int* spatial_hashing_cell_edge_size,
	bool is_self)
{
	double hash_cell_length = cell_length;
	unsigned int hash_cell_count_ = hash_cell_count;
	uint64_t p1 = P1;
	uint64_t p2 = P2;
	uint64_t p3 = P3;
	unsigned int edge_start = obj_edge_index_begin_per_thread[obj_No][thread_No];
	unsigned int edge_end = obj_edge_index_begin_per_thread[obj_No][thread_No + 1];
	std::array<double, 6>* edge_aabb = obj_edge_aabb[obj_No];

	double* scene_aabb_ = scene_aabb;
	std::vector<unsigned int> spatial_hash_index;
	spatial_hash_index.reserve(16);

	bool** obj_is_used_ = obj_is_used[thread_No];

	unsigned int* cell_edge_index;
	unsigned int edge_size;

	std::vector<unsigned int>* neighbor_edge;

	std::vector<unsigned int> primitive_index_record_;
	primitive_index_record_.reserve(40);

	unsigned int* edge_vertex_index_;
	//unsigned int* face_edge_;
	edge_vertex_index_ = edge_vertex_index[obj_No];
	//face_edge_ = this->face_edge[obj_No];


	std::vector<unsigned int>* vertex_neighbor_triangle;

	for (unsigned int i = edge_start; i < edge_end; ++i) {
		primitive_index_record_.clear();
		triangleHashValue(edge_aabb[i].data(), &spatial_hash_index, scene_aabb_, hash_cell_length, hash_cell_count_, p1, p2, p3);
		if (is_self) {
			neighbor_edge = &vertices[obj_No][edge_vertex_index_[i << 1]].edge;
			for (unsigned int j = 0; j < neighbor_edge->size(); ++j) {
				obj_is_used_[obj_No][neighbor_edge->data()[j]] = true;
				primitive_index_record_.emplace_back(obj_No);
				primitive_index_record_.emplace_back(neighbor_edge->data()[j]);
			}
			neighbor_edge = &vertices[obj_No][edge_vertex_index_[(i << 1) + 1]].edge;
			for (unsigned int j = 0; j < neighbor_edge->size(); ++j) {
				obj_is_used_[obj_No][neighbor_edge->data()[j]] = true;
				primitive_index_record_.emplace_back(obj_No);
				primitive_index_record_.emplace_back(neighbor_edge->data()[j]);
			}
			for (unsigned int j = 0; j < spatial_hash_index.size(); ++j) {
				cell_edge_index = spatial_hashing_cell_edge + spatial_hash_index[j] * max_index_number_in_one_cell_edge_;
				edge_size = spatial_hashing_cell_edge_size[spatial_hash_index[j]];
				for (unsigned int k = 0; k < edge_size; k += 2) {
					if (!obj_is_used_[cell_edge_index[k + 1]][cell_edge_index[k]]) {
						obj_is_used_[cell_edge_index[k + 1]][cell_edge_index[k]] = true;
						primitive_index_record_.emplace_back(cell_edge_index[k + 1]);
						primitive_index_record_.emplace_back(cell_edge_index[k]);
						if (cell_edge_index[k + 1] == obj_No) {
							//	for (unsigned int m = 0; m < representative_edge_num[cell_edge_index[k + 1]][cell_edge_index[k]]; ++m) {
							if (i < cell_edge_index[k]) {
								if (AABB::AABB_intersection(edge_aabb[i].data(),
									edge_aabb[cell_edge_index[k]].data())) {
									*(primitive_pair_++) = i;
									*(primitive_pair_++) = obj_No;
									*(primitive_pair_++) = cell_edge_index[k];
									*(primitive_pair_++) = obj_No;
								}
							}
						}
						//}
						else {
							//	for (unsigned int m = 0; m < representative_edge_num[cell_edge_index[k + 1]][cell_edge_index[k]]; ++m) {
							if (AABB::AABB_intersection(edge_aabb[i].data(),
								obj_edge_aabb_[cell_edge_index[k + 1]][cell_edge_index[k]].data())) {
								*(primitive_pair_++) = i;
								*(primitive_pair_++) = obj_No;
								*(primitive_pair_++) = cell_edge_index[k];
								*(primitive_pair_++) = cell_edge_index[k + 1];
							}
							//	}
						}
					}
				}
			}
		}
		else {
			for (unsigned int j = 0; j < spatial_hash_index.size(); ++j) {
				cell_edge_index = spatial_hashing_cell_edge + spatial_hash_index[j] * max_index_number_in_one_cell_edge_;
				edge_size = spatial_hashing_cell_edge_size[spatial_hash_index[j]];
				for (unsigned int k = 0; k < edge_size; k += 2) {
					if (!obj_is_used_[cell_edge_index[k + 1]][cell_edge_index[k]]) {
						obj_is_used_[cell_edge_index[k + 1]][cell_edge_index[k]] = true;
						primitive_index_record_.emplace_back(cell_edge_index[k + 1]);
						primitive_index_record_.emplace_back(cell_edge_index[k]);
						//for (unsigned int m = 0; m < representative_edge_num[cell_edge_index[k + 1]][cell_edge_index[k]]; ++m) {
						if (AABB::AABB_intersection(edge_aabb[i].data(),
							obj_edge_aabb_[cell_edge_index[k + 1]][cell_edge_index[k]].data())) {
							*(primitive_pair_++) = i;
							*(primitive_pair_++) = obj_No;
							*(primitive_pair_++) = cell_edge_index[k];
							*(primitive_pair_++) = cell_edge_index[k + 1];
						}
						//}
					}
				}
			}
		}
		for (unsigned int j = 0; j < primitive_index_record_.size(); j += 2) {
			obj_is_used_[primitive_index_record_[j]][primitive_index_record_[j + 1]] = false;
		}
	}
}



void SpatialHashing::findAllEdgeEdgePairsByPrimitiveSingleObjByEdge(int thread_No, int obj_No, unsigned int* primitive_pair, unsigned int* primitive_pair_num_record_,
	std::vector<std::array<double, 6>*>& obj_edge_aabb_,
	unsigned int max_index_number_in_one_cell_edge_, unsigned int* spatial_hashing_cell_edge, unsigned int* spatial_hashing_cell_edge_size,
	bool is_self)
{
	double hash_cell_length = cell_length;
	unsigned int hash_cell_count_ = hash_cell_count;
	uint64_t p1 = P1;
	uint64_t p2 = P2;
	uint64_t p3 = P3;
	unsigned int edge_start = obj_edge_index_begin_per_thread[obj_No][thread_No];
	unsigned int edge_end = obj_edge_index_begin_per_thread[obj_No][thread_No+1];
	std::array<double, 6>* edge_aabb = obj_edge_aabb[obj_No];

	double* scene_aabb_ = scene_aabb;
	std::vector<unsigned int> spatial_hash_index;
	spatial_hash_index.reserve(16);

	bool** obj_is_used_ = obj_is_used[thread_No];

	unsigned int* cell_edge_index;
	unsigned int edge_size;

	std::vector<unsigned int>* neighbor_edge;

	std::vector<unsigned int> primitive_index_record_;
	primitive_index_record_.reserve(40);

	unsigned int* edge_vertex_index_;
	//unsigned int* face_edge_;
	edge_vertex_index_ = edge_vertex_index[obj_No];
	//face_edge_ = this->face_edge[obj_No];


	std::vector<unsigned int>* vertex_neighbor_triangle;

	unsigned int* primitive_pair_;
	unsigned int* primitive_pair_initial_;

	for (unsigned int i = edge_start; i < edge_end; ++i) {
		primitive_index_record_.clear();
		triangleHashValue(edge_aabb[i].data(), &spatial_hash_index, scene_aabb_, hash_cell_length, hash_cell_count_, p1, p2, p3);

		primitive_pair_ = primitive_pair + i * total_length_every_element_for_edge_edge;
		primitive_pair_initial_ = primitive_pair_;
		neighbor_edge = &vertices[obj_No][edge_vertex_index_[i<<1]].edge;
		for (unsigned int j = 0; j < neighbor_edge->size(); ++j) {
			obj_is_used_[obj_No][neighbor_edge->data()[j]] = true;
			primitive_index_record_.emplace_back(obj_No);
			primitive_index_record_.emplace_back(neighbor_edge->data()[j]);
		}
		neighbor_edge = &vertices[obj_No][edge_vertex_index_[(i << 1)+1]].edge;
		for (unsigned int j = 0; j < neighbor_edge->size(); ++j) {
			obj_is_used_[obj_No][neighbor_edge->data()[j]] = true;
			primitive_index_record_.emplace_back(obj_No);
			primitive_index_record_.emplace_back(neighbor_edge->data()[j]);
		}
		for (unsigned int j = 0; j < spatial_hash_index.size(); ++j) {
			cell_edge_index = spatial_hashing_cell_edge + spatial_hash_index[j] * max_index_number_in_one_cell_edge_;
			edge_size = spatial_hashing_cell_edge_size[spatial_hash_index[j]];
			for (unsigned int k = 0; k < edge_size; k += 2) {
				if (!obj_is_used_[cell_edge_index[k + 1]][cell_edge_index[k]]) {
					obj_is_used_[cell_edge_index[k + 1]][cell_edge_index[k]] = true;
					primitive_index_record_.emplace_back(cell_edge_index[k + 1]);
					primitive_index_record_.emplace_back(cell_edge_index[k]);			
					if (AABB::AABB_intersection(edge_aabb[i].data(),
						obj_edge_aabb_[cell_edge_index[k + 1]][cell_edge_index[k]].data())) {
						*(primitive_pair_++) = cell_edge_index[k + 1];
						*(primitive_pair_++) = cell_edge_index[k];
					}								
				}
			}
		}
		
		for (unsigned int j = 0; j < primitive_index_record_.size(); j += 2) {
			obj_is_used_[primitive_index_record_[j]][primitive_index_record_[j + 1]] = false;
		}
		primitive_pair_num_record_[i] = primitive_pair_ - primitive_pair_initial_;
	}
}




void SpatialHashing::findAllVertexTrianglePairsByPrimitiveSingleObj(int thread_No, int obj_No, unsigned int*& primitive_pair_,
	unsigned int vertex_start, unsigned int vertex_end,  std::array<double, 6>* vertex_aabb, 
	unsigned int max_index_number_in_one_cell_triangle_, unsigned int* spatial_hashing_cell_triangle, unsigned int* spatial_hashing_cell_triangle_size,
	std::vector<std::array<double, 6>*>& obj_tri_aabb, bool is_self, bool is_tet)
{
	//unsigned int end = obj_vertex_index_begin_per_thread[obj_No][thread_No + 1];
	//unsigned int start = obj_vertex_index_begin_per_thread[obj_No][thread_No];
	//std::array<double, 6>* vertex_aabb = obj_vertex_aabb[obj_No];

	double hash_cell_length = cell_length;
	unsigned int hash_cell_count_ = hash_cell_count;
	uint64_t p1 = P1;
	uint64_t p2 = P2;
	uint64_t p3 = P3;

	double* scene_aabb_ = scene_aabb;
	std::vector<unsigned int> spatial_hash_index;
	spatial_hash_index.reserve(16);

	bool** obj_is_used_ = obj_is_used[thread_No];

	//unsigned int max_index_number_in_one_cell_triangle_ = max_index_number_in_one_cell;

	unsigned int* cell_triangle_index;
	unsigned int triangle_size;

	std::vector<unsigned int>* neighbor_vertex;

	std::vector<unsigned int> primitive_index_record_;
	primitive_index_record_.reserve(40);
	unsigned int i;


	for (unsigned int ind = vertex_start; ind < vertex_end; ++ind) {
		primitive_index_record_.clear();
		if (is_tet) {
			i = tetrahedron_vertex_index_on_surface[obj_No - cloth->size()][ind];
		}
		else {
			i = ind;
		}
		triangleHashValue(vertex_aabb[i].data(), &spatial_hash_index, scene_aabb_, hash_cell_length, hash_cell_count_, p1, p2, p3);
		if (is_self) {
			neighbor_vertex = &vertices[obj_No][i].face;

			for (unsigned int j = 0; j < neighbor_vertex->size(); ++j) {
				obj_is_used_[obj_No][neighbor_vertex->data()[j]] = true;
				primitive_index_record_.emplace_back(obj_No);
				primitive_index_record_.emplace_back(neighbor_vertex->data()[j]);
			}
		}
		for (unsigned int j = 0; j < spatial_hash_index.size(); ++j) {
			cell_triangle_index = spatial_hashing_cell_triangle + spatial_hash_index[j] * max_index_number_in_one_cell_triangle_;
			triangle_size= spatial_hashing_cell_triangle_size[spatial_hash_index[j]];
			for (unsigned int k = 0; k < triangle_size; k+=2) {
				if (!obj_is_used_[cell_triangle_index[k + 1]][cell_triangle_index[k]]) {
					obj_is_used_[cell_triangle_index[k + 1]][cell_triangle_index[k]] = true;
					primitive_index_record_.emplace_back(cell_triangle_index[k + 1]);
					primitive_index_record_.emplace_back(cell_triangle_index[k]);					
					if (AABB::AABB_intersection(vertex_aabb[i].data(), obj_tri_aabb[cell_triangle_index[k + 1]][cell_triangle_index[k]].data())) {
						*(primitive_pair_++) = i;
						*(primitive_pair_++) = obj_No;
						*(primitive_pair_++) = cell_triangle_index[k];
						*(primitive_pair_++) = cell_triangle_index[k + 1];						
					}
				}				
			}
		}
		for (unsigned int j = 0; j < primitive_index_record_.size(); j += 2) {
			obj_is_used_[primitive_index_record_[j]][primitive_index_record_[j + 1]] = false;
		}
	}
}


void SpatialHashing::findAllTriangleVertexPairByTriangleSingleObj(int obj_No, unsigned int* vt_pair_initial, unsigned int* vt_pair_num, unsigned int** tv_pair,
	unsigned int** tv_pair_num, unsigned int total_length_every_element_vt, unsigned int total_length_every_element_tv)
{
	unsigned int vertex_end,i;
	unsigned int total_num;

	bool is_tet = obj_No >= cloth->size();
	unsigned int* temp_address;

	unsigned int* vt_pair;
	if (obj_No < cloth->size()) {
		vertex_end = mesh_struct[obj_No]->vertex_position.size();
	}
	else {
		vertex_end = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_index_on_sureface.size();
	}

	for (unsigned int ind = 0; ind < vertex_end; ++ind) {
		if (is_tet) {
			i = tetrahedron_vertex_index_on_surface[obj_No - cloth->size()][ind];
		}
		else {
			i = ind;
		}
		total_num = vt_pair_num[ind];
		vt_pair = vt_pair_initial + ind * total_length_every_element_vt;
		for (unsigned int j = 0; j < total_num; j += 2) {
			temp_address = tv_pair[vt_pair[j]] + vt_pair[j + 1] * total_length_every_element_tv + tv_pair_num[vt_pair[j]][vt_pair[j + 1]];
			*(temp_address++) = obj_No;
			*temp_address = i;
			tv_pair_num[vt_pair[j]][vt_pair[j + 1]]++;
		}
	}
}



void SpatialHashing::findAllVertexTrianglePairsByPrimitiveSingleObj_ByVertex(int thread_No, int obj_No, unsigned int* primitive_pair, 
	unsigned int* primitive_pair_num_record,
	unsigned int vertex_start, unsigned int vertex_end, std::array<double, 6>* vertex_aabb,
	unsigned int max_index_number_in_one_cell_triangle_, unsigned int* spatial_hashing_cell_triangle, unsigned int* spatial_hashing_cell_triangle_size,
	std::vector<std::array<double, 6>*>& obj_tri_aabb, bool is_self, bool is_tet)
{

	double hash_cell_length = cell_length;
	unsigned int hash_cell_count_ = hash_cell_count;
	uint64_t p1 = P1;
	uint64_t p2 = P2;
	uint64_t p3 = P3;

	double* scene_aabb_ = scene_aabb;
	std::vector<unsigned int> spatial_hash_index;
	spatial_hash_index.reserve(16);

	bool** obj_is_used_ = obj_is_used[thread_No];

	unsigned int* cell_triangle_index;
	unsigned int triangle_size;

	std::vector<unsigned int>* neighbor_vertex;

	std::vector<unsigned int> primitive_index_record_;
	primitive_index_record_.reserve(40);
	unsigned int i;


	unsigned int* primitive_pair_;
	unsigned int* primitive_pair_initial_;

	for (unsigned int ind = vertex_start; ind < vertex_end; ++ind) {
		primitive_index_record_.clear();
		if (is_tet) {
			i = tetrahedron_vertex_index_on_surface[obj_No - cloth->size()][ind];
		}
		else {
			i = ind;
		}
		triangleHashValue(vertex_aabb[i].data(), &spatial_hash_index, scene_aabb_, hash_cell_length, hash_cell_count_, p1, p2, p3);
		if (is_self) {
			primitive_pair_ = primitive_pair + ind * total_length_every_element_for_vertex_triange;
			neighbor_vertex = &vertices[obj_No][i].face;
			for (unsigned int j = 0; j < neighbor_vertex->size(); ++j) {
				obj_is_used_[obj_No][neighbor_vertex->data()[j]] = true;
				primitive_index_record_.emplace_back(obj_No);
				primitive_index_record_.emplace_back(neighbor_vertex->data()[j]);
			}			
		}
		else {
			primitive_pair_ = primitive_pair + ind * total_length_every_element_for_vertex_triange_collider;
		}
		primitive_pair_initial_ = primitive_pair_;
		for (unsigned int j = 0; j < spatial_hash_index.size(); ++j) {
			cell_triangle_index = spatial_hashing_cell_triangle + spatial_hash_index[j] * max_index_number_in_one_cell_triangle_;
			triangle_size = spatial_hashing_cell_triangle_size[spatial_hash_index[j]];
			for (unsigned int k = 0; k < triangle_size; k += 2) {
				if (!obj_is_used_[cell_triangle_index[k + 1]][cell_triangle_index[k]]) {
					obj_is_used_[cell_triangle_index[k + 1]][cell_triangle_index[k]] = true;
					primitive_index_record_.emplace_back(cell_triangle_index[k + 1]);
					primitive_index_record_.emplace_back(cell_triangle_index[k]);
					if (AABB::AABB_intersection(vertex_aabb[i].data(), obj_tri_aabb[cell_triangle_index[k + 1]][cell_triangle_index[k]].data())) {
						*(primitive_pair_++) = cell_triangle_index[k + 1];
						*(primitive_pair_++) = cell_triangle_index[k];
					}
				}
			}
		}
		for (unsigned int j = 0; j < primitive_index_record_.size(); j += 2) {
			obj_is_used_[primitive_index_record_[j]][primitive_index_record_[j + 1]] = false;
		}
		primitive_pair_num_record[ind] = primitive_pair_ - primitive_pair_initial_;
	}
}

void SpatialHashing::searchPrimitive(double* aabb, unsigned int* hash_index, unsigned int hash_cell_size, unsigned int input_obj_No, unsigned int triangle_index,
	unsigned int*& triangle_pair_, unsigned int*& triangle_pair_with_collider_, unsigned int thread_No, 
	unsigned int* spatial_hashing_cell_, unsigned int* spatial_hash_size,	
	unsigned int* spatial_hashing_cell_collider_, unsigned int* spatial_hash_size_collider, bool is_edge,
	std::vector<unsigned int>* neighbor_primitive_0, std::vector<unsigned int>* neighbor_primitive_1,
	unsigned int max_index_number_in_one_cell, unsigned int max_index_number_in_one_cell_collider)
{
	unsigned int* primitive_index_record_ = primitive_index_record[thread_No] + 1;
	bool** obj_is_used_ = obj_is_used[thread_No];
	unsigned int hash_cell_index_end;
	unsigned int obj_No; unsigned int current_triangle_index;
	unsigned int obj_num = cloth->size() + tetrahedron->size();


	for (unsigned int i = 0; i < neighbor_primitive_0->size(); ++i) {
		obj_is_used_[input_obj_No][neighbor_primitive_0->data()[i]] = true;
		*(primitive_index_record_++) = neighbor_primitive_0->data()[i];
		*(primitive_index_record_++) = input_obj_No;
	}		
	if (is_edge) {
		for (unsigned int i = 0; i < neighbor_primitive_1->size(); ++i) {
			obj_is_used_[input_obj_No][neighbor_primitive_1->data()[i]] = true;
			*(primitive_index_record_++) = neighbor_primitive_1->data()[i];
			*(primitive_index_record_++) = input_obj_No;
		}
	}
	if (is_edge) {
		for (unsigned int i = 0; i < hash_cell_size; ++i) {
			for (unsigned int j = 0; j < spatial_hash_size[hash_index[i]]; j += 2) {
				obj_No = spatial_hashing_cell_[hash_index[i]* max_index_number_in_one_cell + j + 1];
				current_triangle_index = spatial_hashing_cell_[hash_index[i]* max_index_number_in_one_cell +j];
				if (!obj_is_used_[obj_No][current_triangle_index]) {
					obj_is_used_[obj_No][current_triangle_index] = true;
					memcpy(primitive_index_record_, spatial_hashing_cell_+hash_index[i]* max_index_number_in_one_cell + j, 8);
					primitive_index_record_ += 2;
					if (obj_No > input_obj_No ||
						(obj_No == input_obj_No && current_triangle_index > triangle_index)) {
						if (AABB::AABB_intersection(aabb, obj_edge_aabb[obj_No][current_triangle_index].data())) {
							*(triangle_pair_++) = triangle_index;
							*(triangle_pair_++) = input_obj_No;
							*(triangle_pair_++) = current_triangle_index;
							*(triangle_pair_++) = obj_No;
						}
					}
				}
			}
		}
	}
	else {
		for (unsigned int i = 0; i < hash_cell_size; ++i) {
			for (unsigned int j = 0; j < spatial_hash_size[hash_index[i]]; j += 2) {
				obj_No = spatial_hashing_cell_[hash_index[i]* max_index_number_in_one_cell +j + 1];
				current_triangle_index = spatial_hashing_cell_[hash_index[i]* max_index_number_in_one_cell+j];
				if (!obj_is_used_[obj_No][current_triangle_index]) {
					obj_is_used_[obj_No][current_triangle_index] = true;
					memcpy(primitive_index_record_, spatial_hashing_cell_+hash_index[i]* max_index_number_in_one_cell + j, 8);
					primitive_index_record_ += 2;				
					if (AABB::AABB_intersection(aabb, obj_tri_aabb[obj_No][current_triangle_index].data())) {
						*(triangle_pair_++) = triangle_index;
						*(triangle_pair_++) = input_obj_No;
						*(triangle_pair_++) = current_triangle_index;
						*(triangle_pair_++) = obj_No;
					}					
				}
			}
		}
	}
	

	primitive_index_record[thread_No][0] = primitive_index_record_ - primitive_index_record[thread_No] - 1;
	primitive_index_record_ = primitive_index_record[thread_No];
	for (unsigned int i = 0; i < primitive_index_record_[0]; i += 2) {
		obj_is_used_[primitive_index_record_[i + 2]][primitive_index_record_[i + 1]] = false;
	}

	if (has_collider) {
		
		if (is_edge) {
			//*(primitive_index_record_++) = 0;
			//for (unsigned int i = 0; i < hash_cell_size; ++i) {
			//	for (unsigned int j = 0; j < spatial_hash_size_collider[hash_index[i]]; j += 2) {
			//		obj_No = spatial_hashing_cell_collider_[hash_index[i]* max_index_number_in_one_cell_collider + j + 1];
			//		current_triangle_index = spatial_hashing_cell_collider_[hash_index[i]* max_index_number_in_one_cell_collider + j];
			//		if (!obj_is_used_[obj_No][current_triangle_index]) {
			//			obj_is_used_[obj_No][current_triangle_index] = true;
			//			memcpy(primitive_index_record_, spatial_hashing_cell_collider_+hash_index[i]* max_index_number_in_one_cell_collider + j, 8);
			//			primitive_index_record_ += 2;
			//			if (AABB::AABB_intersection(aabb, collider_edge_aabb[obj_No][current_triangle_index].data())) {
			//				*(triangle_pair_with_collider_++) = triangle_index;
			//				*(triangle_pair_with_collider_++) = input_obj_No;
			//				*(triangle_pair_with_collider_++) = current_triangle_index;
			//				*(triangle_pair_with_collider_++) = obj_No;
			//			}
			//		}
			//	}
			//}
		}
		else {
			*(primitive_index_record_++) = 0;
			for (unsigned int i = 0; i < hash_cell_size; ++i) {
				for (unsigned int j = 0; j < spatial_hash_size_collider[hash_index[i]]; j += 2) {
					obj_No = spatial_hashing_cell_collider_[hash_index[i]* max_index_number_in_one_cell_collider + j + 1];
					current_triangle_index = spatial_hashing_cell_collider_[hash_index[i]* max_index_number_in_one_cell_collider + j];
					if (!obj_is_used_[obj_No][current_triangle_index]) {
						obj_is_used_[obj_No][current_triangle_index] = true;
						memcpy(primitive_index_record_, spatial_hashing_cell_collider_ + hash_index[i]* max_index_number_in_one_cell_collider + j, 8);
						primitive_index_record_ += 2;
						if (AABB::AABB_intersection(aabb, collider_tri_aabb[obj_No][current_triangle_index].data())) {
							*(triangle_pair_with_collider_++) = triangle_index;
							*(triangle_pair_with_collider_++) = input_obj_No;
							*(triangle_pair_with_collider_++) = current_triangle_index;
							*(triangle_pair_with_collider_++) = obj_No;
						}
					}
				}
			}
		}

		primitive_index_record[thread_No][0] = primitive_index_record_ - primitive_index_record[thread_No] - 1;
		primitive_index_record_ = primitive_index_record[thread_No];
		for (unsigned int i = 0; i < primitive_index_record_[0]; i += 2) {
			obj_is_used_[primitive_index_record_[i + 2]][primitive_index_record_[i + 1]] = false;
		}
	}
}







//TRIANGLE_HASHING_SMALLER_HASH_TABLE
void SpatialHashing::triangleHashingSmallerHashTable(int thread_No)
{
	double scene_aabb_[6];
	memcpy(scene_aabb_, scene_aabb, 48);
	double hash_cell_length = cell_length;
	unsigned int hash_cell_count_ = hash_cell_count;

	uint64_t p1 = P1;
	uint64_t p2 = P2;
	uint64_t p3 = P3;
	bool has_collider_ = has_collider;
	unsigned int collider_begin_obj_index_ = collider_begin_obj_index;
	unsigned int collider_size = collider->size();
	unsigned int primitive_begin;
	unsigned int primitive_end;
	std::array<double, 6>* aabb;
	int vector_size;
	unsigned int largest_count_in_hash_value_list;
	unsigned int* spatial_hashing_cell_ = spatial_hashing_cell_triangle[thread_No];
	unsigned int* spatial_hashing_cell_size_ = spatial_hashing_cell_triangle_size[thread_No];
	memset(spatial_hashing_cell_size_, 0, hash_cell_count_ << 2);	

	unsigned int* spatial_hashing_cell_collider_;
	unsigned int* spatial_hashing_cell_collider_size_;

	unsigned int cloth_size = cloth->size();

	unsigned int max_index_number_in_one_cell_= max_index_number_in_one_cell;
	unsigned int max_index_number_in_one_cell_collider_= max_index_number_in_one_cell_collider;


	unsigned int primitive_index[2];//first triangle_index, second obj index
	for (unsigned int i = 0; i < collider_begin_obj_index_; ++i) {
		aabb = obj_tri_aabb[i];
		primitive_begin = obj_triangle_index_begin_per_thread[i][thread_No];
		primitive_end = obj_triangle_index_begin_per_thread[i][thread_No + 1];
		primitive_index[1] = i;
		for (unsigned int j = primitive_begin; j < primitive_end; ++j) {
			primitive_index[0] = j;
			triangleHashValueWithoutRecord(aabb[j].data(), spatial_hashing_cell_, primitive_index, scene_aabb_, hash_cell_length,
				hash_cell_count_, p1, p2, p3, spatial_hashing_cell_size_,max_index_number_in_one_cell_);
		}
	}
	if (has_collider_) {
		spatial_hashing_cell_collider_ = spatial_hashing_cell_collider_triangle[thread_No];
		spatial_hashing_cell_collider_size_ = spatial_hashing_cell_collider_triangle_size[thread_No];
		memset(spatial_hashing_cell_collider_size_, 0, hash_cell_count_ << 2);
		for (unsigned int i = 0; i < collider_size; ++i) {
			aabb = (*collider)[i].triangle_AABB.data();
			primitive_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread[thread_No];
			primitive_end = (*collider)[i].mesh_struct.face_index_begin_per_thread[thread_No + 1];
			primitive_index[1] = i;
			for (unsigned int j = primitive_begin; j < primitive_end; ++j) {
				primitive_index[0] = j;
				triangleHashValueWithoutRecord(aabb[j].data(), spatial_hashing_cell_collider_, primitive_index, scene_aabb_,
					hash_cell_length, hash_cell_count_, p1, p2, p3, spatial_hashing_cell_collider_size_, max_index_number_in_one_cell_collider_);
			}
		}
	}


		memset(spatial_hashing_cell_edge_size[thread_No], 0, hash_cell_count_ << 2);
	
		//edge 
		spatial_hashing_cell_ = spatial_hashing_cell_edge[thread_No];
		spatial_hashing_cell_size_ = spatial_hashing_cell_edge_size[thread_No];
		max_index_number_in_one_cell_ = max_index_number_in_one_cell_edge;
		max_index_number_in_one_cell_collider_ = max_index_number_in_one_cell_collider_edge;

		//memset(spatial_hashing_cell_size_, 0, hash_cell_count_ << 2);
		for (unsigned int i = 0; i < collider_begin_obj_index_; ++i) {
			aabb = obj_edge_aabb[i];
			primitive_begin = obj_edge_index_begin_per_thread[i][thread_No];
			primitive_end = obj_edge_index_begin_per_thread[i][thread_No + 1];
			primitive_index[1] = i;
			for (unsigned int j = primitive_begin; j < primitive_end; ++j) {
				primitive_index[0] = j;

				//use this for loop by every element
			/*	triangleHashValueWithRecord(aabb[j].data(), spatial_hashing_cell_, primitive_index, scene_aabb_, hash_cell_length,
					hash_cell_count_, p1, p2, p3, spatial_hashing_cell_size_, obj_edge_hash[i] + j * max_triangle_cell_size, max_index_number_in_one_cell_);*/
					//use this for loop by every cell
				triangleHashValueWithoutRecord(aabb[j].data(), spatial_hashing_cell_, primitive_index, scene_aabb_, hash_cell_length,
					hash_cell_count_, p1, p2, p3, spatial_hashing_cell_size_, max_index_number_in_one_cell_);
			}
		}
		if (searchPairByCell) {
		if (has_collider_) {
			spatial_hashing_cell_collider_ = spatial_hashing_cell_collider_edge[thread_No];
			spatial_hashing_cell_collider_size_ = spatial_hashing_cell_collider_edge_size[thread_No];
			memset(spatial_hashing_cell_collider_size_, 0, hash_cell_count_ << 2);
			for (unsigned int i = 0; i < collider_size; ++i) {
				aabb = (*collider)[i].edge_AABB.data();
				primitive_begin = (*collider)[i].mesh_struct.edge_index_begin_per_thread[thread_No];
				primitive_end = (*collider)[i].mesh_struct.edge_index_begin_per_thread[thread_No + 1];
				primitive_index[1] = i;
				for (unsigned int j = primitive_begin; j < primitive_end; ++j) {
					primitive_index[0] = j;
					triangleHashValueWithoutRecord(aabb[j].data(), spatial_hashing_cell_collider_, primitive_index, scene_aabb_,
						hash_cell_length, hash_cell_count_, p1, p2, p3, spatial_hashing_cell_collider_size_, max_index_number_in_one_cell_collider_);
				}
			}
		}


		memset(spatial_hashing_cell_vertex_size[thread_No], 0, hash_cell_count_ << 2);
		////vertex
		spatial_hashing_cell_ = spatial_hashing_cell_vertex[thread_No];
		spatial_hashing_cell_size_ = spatial_hashing_cell_vertex_size[thread_No];
		//memset(spatial_hashing_cell_size_, 0, hash_cell_count_ << 2);
		unsigned int* vertex_index_on_surface;
		max_index_number_in_one_cell_ = max_index_number_in_one_cell_vertex;
		max_index_number_in_one_cell_collider_ = max_index_number_in_one_cell_collider_vertex;
		for (unsigned int i = 0; i < collider_begin_obj_index_; ++i) {
			aabb = obj_vertex_aabb[i];
			primitive_begin = obj_vertex_index_begin_per_thread[i][thread_No];
			primitive_end = obj_vertex_index_begin_per_thread[i][thread_No + 1];
			primitive_index[1] = i;
			if (i < cloth_size) {
				for (unsigned int j = primitive_begin; j < primitive_end; ++j) {
					primitive_index[0] = j;
					//use this for loop by every element
					//triangleHashValueWithRecord(aabb[j].data(), spatial_hashing_cell_, primitive_index, scene_aabb_, hash_cell_length,
					//	hash_cell_count_, p1, p2, p3, spatial_hashing_cell_size_, obj_vertex_hash[i] + j * max_vertex_cell_size, max_index_number_in_one_cell_);
					//use this for loop by every cell
					triangleHashValueWithoutRecord(aabb[j].data(), spatial_hashing_cell_, primitive_index, scene_aabb_, hash_cell_length,
						hash_cell_count_, p1, p2, p3, spatial_hashing_cell_size_, max_index_number_in_one_cell_);
				}
			}
			else {
				vertex_index_on_surface = tetrahedron_vertex_index_on_surface[i - cloth_size];
				for (unsigned int j = primitive_begin; j < primitive_end; ++j) {
					primitive_index[0] = vertex_index_on_surface[j];
					//use this for loop by every element
					/*triangleHashValueWithRecord(aabb[vertex_index_on_surface[j]].data(), spatial_hashing_cell_, primitive_index, scene_aabb_, hash_cell_length,
						hash_cell_count_, p1, p2, p3, spatial_hashing_cell_size_, obj_vertex_hash[i] + vertex_index_on_surface[j] * max_vertex_cell_size, max_index_number_in_one_cell_);*/
						//use this for loop by every cell
					triangleHashValueWithoutRecord(aabb[vertex_index_on_surface[j]].data(), spatial_hashing_cell_, primitive_index, scene_aabb_, hash_cell_length,
						hash_cell_count_, p1, p2, p3, spatial_hashing_cell_size_, max_index_number_in_one_cell_);
				}
			}
		}
		if (has_collider_) {
			spatial_hashing_cell_collider_ = spatial_hashing_cell_collider_vertex[thread_No];
			spatial_hashing_cell_collider_size_ = spatial_hashing_cell_collider_vertex_size[thread_No];
			memset(spatial_hashing_cell_collider_size_, 0, hash_cell_count_ << 2);
			for (unsigned int i = 0; i < collider_size; ++i) {
				aabb = (*collider)[i].vertex_AABB.data();
				primitive_begin = (*collider)[i].mesh_struct.vertex_index_begin_per_thread[thread_No];
				primitive_end = (*collider)[i].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
				primitive_index[1] = i;
				for (unsigned int j = primitive_begin; j < primitive_end; ++j) {
					primitive_index[0] = j;
					triangleHashValueWithRecord(aabb[j].data(), spatial_hashing_cell_collider_, primitive_index, scene_aabb_,
						hash_cell_length, hash_cell_count_, p1, p2, p3, spatial_hashing_cell_collider_size_, collider_vertex_hash[i] + j * max_vertex_cell_size,
						max_index_number_in_one_cell_collider_);
				}
			}
		}
	}
}

void SpatialHashing::countHashIndexPerThread(int thread_No, unsigned int actual_count, std::vector<unsigned int>* hash_index_count,
	unsigned int* spatial_hashing_value_)
{
	hash_index_count->clear();
	unsigned int hash_index = spatial_hashing_value_[0];
	unsigned int* end = spatial_hashing_value_ + actual_count;

	while (true) {
		hash_index_count->emplace_back(0);
		while (*spatial_hashing_value_ == hash_index)
		{
			hash_index_count->back()++;
			spatial_hashing_value_++;
			if (end == spatial_hashing_value_) {
				return;
			}
		}
		hash_index = *spatial_hashing_value_;
	}
}


//ORI_TRIANGLE_HASHING
void SpatialHashing::oriTriangleHashing(int thread_No)
{
	double scene_aabb_[6];
	memcpy(scene_aabb_, scene_aabb, 48);
	double hash_cell_length = cell_length;
	unsigned int collider_begin_obj_index_ = collider_begin_obj_index;
	unsigned int collider_size = collider->size();
	unsigned int triangle_begin;
	unsigned int triangle_end;
	std::array<double, 6>* aabb;
	std::vector<unsigned int>* spatial_hashing_index_ = &ori_hash_value[thread_No];

	spatial_hashing_index_->clear();


	for (unsigned int i = 0; i < collider_begin_obj_index_; ++i) {
		aabb = obj_tri_aabb[i];
		triangle_begin = obj_triangle_index_begin_per_thread[i][thread_No];
		triangle_end = obj_triangle_index_begin_per_thread[i][thread_No + 1];
		for (unsigned int j = triangle_begin; j < triangle_end; ++j) {
			initialTriangleHashValue(aabb[j].data(), spatial_hashing_index_, scene_aabb_, hash_cell_length);
		}
	}
	if (has_collider) {
		for (unsigned int i = 0; i < collider_size; ++i) {
			aabb = (*collider)[i].triangle_AABB.data();
			triangle_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread[thread_No];
			triangle_end = (*collider)[i].mesh_struct.face_index_begin_per_thread[thread_No + 1];
			for (unsigned int j = triangle_begin; j < triangle_end; ++j) {
				initialTriangleHashValue(aabb[j].data(), spatial_hashing_index_, scene_aabb_,	hash_cell_length);
			}
		}
	}
}



void SpatialHashing::initialTriangleHashValue(double* aabb,
	std::vector<unsigned int>* spatial_hashing_cell, double* scene_aabb, double cell_length)
{
	std::uint64_t spatial_index[6];
	unsigned int hash_index;
	for (unsigned int j = 0; j < 3; ++j) {
		spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j]) / cell_length);
	}
	for (unsigned int j = 3; j < 6; ++j) {
		spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j - 3]) / cell_length) + 1;
	}
	unsigned int index_y_multi;
	unsigned int index_z_multi;

	for (unsigned int index_y = spatial_index[1]; index_y < spatial_index[4]; ++index_y) {
		index_y_multi = index_y * cell_number[0];
		for (unsigned int index_z = spatial_index[2]; index_z < spatial_index[5]; ++index_z) {
			index_z_multi = index_z * cell_num0_cell_num1;
			for (unsigned int index_x = spatial_index[0]; index_x < spatial_index[3]; ++index_x) {
				spatial_hashing_cell->emplace_back(index_x + index_y_multi + index_z_multi);
			}
		}
	}
}


void SpatialHashing::triangleHashValueWithRecord(double* aabb,
	unsigned int* spatial_hashing_cell, unsigned int* triangle_index, double* scene_aabb, double cell_length,
	unsigned int hash_cell_count, uint64_t P1, uint64_t P2, uint64_t P3, unsigned int* spatial_hashing_cell_triangle_size,
	unsigned int* obj_triangle_hash, unsigned int max_index_number_in_one_cell)
{
	std::uint64_t spatial_index[6];
	unsigned int hash_index;
	for (unsigned int j = 0; j < 3; ++j) {
		spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j]) / cell_length);
	}
	for (unsigned int j = 3; j < 6; ++j) {
		spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j - 3]) / cell_length) + 1;
	}

	unsigned int size = 0;
	for (std::uint64_t index_y = spatial_index[1]; index_y < spatial_index[4]; ++index_y) {
		//index_y_multi = P2 * index_y;
		for (std::uint64_t index_z = spatial_index[2]; index_z < spatial_index[5]; ++index_z) {
			//index_z_multi = index_z * P3;
			for (std::uint64_t index_x = spatial_index[0]; index_x < spatial_index[3]; ++index_x) {
				hash_index = ((index_x * P1) ^ (P2 * index_y) ^ (index_z * P3)) % hash_cell_count;
				//cell_size_index = spatial_hashing_cell[hash_index];
				if (spatial_hashing_cell_triangle_size[hash_index]) {
					if (memcmp(spatial_hashing_cell+hash_index* max_index_number_in_one_cell + spatial_hashing_cell_triangle_size[hash_index] - 2, triangle_index, 8) != 0) {
						memcpy(spatial_hashing_cell+hash_index* max_index_number_in_one_cell + spatial_hashing_cell_triangle_size[hash_index], triangle_index, 8);
						spatial_hashing_cell_triangle_size[hash_index] += 2;

						size++;
						obj_triangle_hash[size] = hash_index;

					}
				}
				else {
					memcpy(spatial_hashing_cell+hash_index* max_index_number_in_one_cell + spatial_hashing_cell_triangle_size[hash_index], triangle_index, 8);
					spatial_hashing_cell_triangle_size[hash_index] += 2;

					size++;
					obj_triangle_hash[size] = hash_index;
				}
			}
		}

	}
	obj_triangle_hash[0] = size;
}

void SpatialHashing::triangleHashValue(double* aabb,
	std::vector<unsigned int>* spatial_hashing_index, double* scene_aabb, double cell_length,
	unsigned int hash_cell_count, uint64_t P1, uint64_t P2, uint64_t P3)
{
	std::uint64_t spatial_index[6];
	//unsigned int hash_index;
	for (unsigned int j = 0; j < 3; ++j) {
		spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j]) / cell_length);
	}
	for (unsigned int j = 3; j < 6; ++j) {
		spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j - 3]) / cell_length) + 1;
	}
	spatial_hashing_index->clear();
	unsigned int size = 0;
	for (std::uint64_t index_y = spatial_index[1]; index_y < spatial_index[4]; ++index_y) {
		//index_y_multi = P2 * index_y;
		for (std::uint64_t index_z = spatial_index[2]; index_z < spatial_index[5]; ++index_z) {
			//index_z_multi = index_z * P3;
			for (std::uint64_t index_x = spatial_index[0]; index_x < spatial_index[3]; ++index_x) {
				//hash_index = ((index_x * P1) ^ (P2 * index_y) ^ (index_z * P3)) % hash_cell_count;
				//cell_size_index = spatial_hashing_cell[hash_index];
				spatial_hashing_index->emplace_back(((index_x * P1) ^ (P2 * index_y) ^ (index_z * P3)) % hash_cell_count);
			}
		}
	}
}



void SpatialHashing::triangleHashValueWithoutRecord(double* aabb,
	unsigned int* spatial_hashing_cell, unsigned int* triangle_index, double* scene_aabb, double cell_length,
	unsigned int hash_cell_count, uint64_t P1, uint64_t P2, uint64_t P3, unsigned int* spatial_hashing_cell_triangle_size,
	unsigned int max_index_number_in_one_cell)
{
	std::uint64_t spatial_index[6];
	//std::uint64_t index_z_multi;
	//std::uint64_t index_y_multi;

	//unsigned int* cell_size_index;

	unsigned int hash_index;
	for (unsigned int j = 0; j < 3; ++j) {
		spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j]) / cell_length);
	}
	for (unsigned int j = 3; j < 6; ++j) {
		spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j - 3]) / cell_length)+ 1;
	}

	//if (test_type == 0) {
	//	std::cout << (spatial_index[4] - spatial_index[1]) * (spatial_index[5] - spatial_index[2]) * (spatial_index[3] - spatial_index[0]) << std::endl;
	//}
	for (std::uint64_t index_y = spatial_index[1]; index_y < spatial_index[4]; ++index_y) {
		//index_y_multi = P2 * index_y;
		for (std::uint64_t index_z = spatial_index[2]; index_z < spatial_index[5]; ++index_z) {
			//index_z_multi = index_z * P3;
			for (std::uint64_t index_x = spatial_index[0]; index_x < spatial_index[3]; ++index_x) {
				hash_index = ((index_x * P1) ^ (P2 * index_y) ^ (index_z * P3)) % hash_cell_count;
				//cell_size_index = spatial_hashing_cell[hash_index];
				if (spatial_hashing_cell_triangle_size[hash_index]) {
					if (memcmp(spatial_hashing_cell + hash_index* max_index_number_in_one_cell + spatial_hashing_cell_triangle_size[hash_index] - 2, triangle_index, 8) != 0) {
						memcpy(spatial_hashing_cell + hash_index * max_index_number_in_one_cell + spatial_hashing_cell_triangle_size[hash_index], triangle_index, 8);
						spatial_hashing_cell_triangle_size[hash_index] += 2;
						//cell_size_index[*cell_size_index] = triangle_index;
						//(*cell_size_index)++;
						//cell_size_index[*cell_size_index] = obj_index;
					}
				}
				else {
					memcpy(spatial_hashing_cell + hash_index* max_index_number_in_one_cell + spatial_hashing_cell_triangle_size[hash_index], triangle_index, 8);
					spatial_hashing_cell_triangle_size[hash_index] += 2;
					//(*cell_size_index)++;
					//cell_size_index[*cell_size_index] = triangle_index;
					//(*cell_size_index)++;
					//cell_size_index[*cell_size_index] = obj_index;
				}
			}
		}
	}
}



void SpatialHashing::getSceneAABB()
{

	thread->assignTask(this, SCENE_AABB);
	memcpy(scene_aabb, scene_aabb_thread[0].data(), 48);
	for (int i = 1; i < thread_num; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (scene_aabb[j] > scene_aabb_thread[i][j]) {
				scene_aabb[j] = scene_aabb_thread[i][j];
			}
		}
		for (int j = 3; j < 6; ++j) {
			if (scene_aabb[j] < scene_aabb_thread[i][j]) {
				scene_aabb[j] = scene_aabb_thread[i][j];
			}
		}
	}
	for (int i = 0; i < 3; ++i) {
		scene_aabb[i + 3] += 0.1;
		scene_aabb[i] -= 0.1;
	}
}


//SCENE_AABB
void SpatialHashing::getSceneAABB(int thread_No)
{
	double* scene_aabb_ = scene_aabb_thread[thread_No].data();
	memset(scene_aabb_ + 3, 0xFE, 24); //set double to -5.31401e+303
	memset(scene_aabb_, 0x7F, 24); //set double to 1.38242e+306
	std::array<double, 6>* obj_aabb;
	int vertex_index_begin;
	int vertex_index_end;
	double* vertex;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		obj_aabb = (*cloth)[i].vertex_AABB.data();
		vertex_index_begin = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_index_end = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
		for (unsigned int j = vertex_index_begin; j < vertex_index_end; ++j) {
			vertex = obj_aabb[j].data();
			for (unsigned int k = 0; k < 3; ++k) {
				if (scene_aabb_[k] > vertex[k]) {
					scene_aabb_[k] = vertex[k];
				}
			}
			for (unsigned int k = 3; k < 6; ++k) {
				if (scene_aabb_[k] < vertex[k]) {
					scene_aabb_[k] = vertex[k];
				}
			}
		}
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		obj_aabb = (*collider)[i].vertex_AABB.data();
		vertex_index_begin = (*collider)[i].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_index_end = (*collider)[i].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
		for (unsigned int j = vertex_index_begin; j < vertex_index_end; ++j) {
			vertex = obj_aabb[j].data();
			for (unsigned int k = 0; k < 3; ++k) {
				if (scene_aabb_[k] > vertex[k]) {
					scene_aabb_[k] = vertex[k];
				}
			}
			for (unsigned int k = 3; k < 6; ++k) {
				if (scene_aabb_[k] < vertex[k]) {
					scene_aabb_[k] = vertex[k];
				}
			}
		}
	}

	for (int i = 0; i < tetrahedron->size(); ++i) {
		obj_aabb = (*tetrahedron)[i].vertex_AABB.data();
		vertex_index_begin = (*tetrahedron)[i].mesh_struct.vertex_index_on_surface_begin_per_thread[thread_No];
		vertex_index_end = (*tetrahedron)[i].mesh_struct.vertex_index_on_surface_begin_per_thread[thread_No + 1];
		for (int j = vertex_index_begin; j < vertex_index_end; ++j) {
			vertex = obj_aabb[j].data();
			for (unsigned int k = 0; k < 3; ++k) {
				if (scene_aabb_[k] > vertex[k]) {
					scene_aabb_[k] = vertex[k];
				}
			}
			for (unsigned int k = 3; k < 6; ++k) {
				if (scene_aabb_[k] < vertex[k]) {
					scene_aabb_[k] = vertex[k];
				}
			}
		}
	}
}


//find the nearset prime number that is not less than the input
unsigned int SpatialHashing::findNeareastPrimeNumber(unsigned int index)
{
	if (index < 4) {
		return index;
	}
	if (index < 6) {
		return 5;
	}
	if (index % 2 == 0) {
		index++;
	}
	unsigned int coeff;
	bool is_prime;
	while (true)
	{
		if (index % 10 == 5) {
			index += 2;
			continue;
		}
		is_prime = true;
		coeff = sqrt(index);
		for (unsigned int i = 3; i <= coeff; ++i) {
			if (index % i == 0) {
				is_prime = false;
				break;
			}
		}
		if (is_prime) {
			return index;
		}
		else {
			index += 2;
		}

	}
}



void SpatialHashing::recordNonEmptyCell()
{
	thread->assignTask(this, RECORD_NONEMPTY_CELL);

	for (unsigned int i = 1; i < thread_num; ++i) {
		memcpy(non_empty_cell_index_vertex_triangle[0] + 1 + non_empty_cell_index_vertex_triangle[0][0],
			non_empty_cell_index_vertex_triangle[i] + 1,
			4 * non_empty_cell_index_vertex_triangle[i][0]);
		non_empty_cell_index_vertex_triangle[0][0] += non_empty_cell_index_vertex_triangle[i][0];

		memcpy(non_empty_cell_index_edge[0] + 1 + non_empty_cell_index_edge[0][0],
			non_empty_cell_index_edge[i] + 1,
			4 * non_empty_cell_index_edge[i][0]);
		non_empty_cell_index_edge[0][0] += non_empty_cell_index_edge[i][0];
	}

	//std::cout << "non empty cell " << non_empty_cell_index[0][0] << std::endl;
	
	arrangeIndex(thread_num, non_empty_cell_index_vertex_triangle[0][0], non_empty_cell_index_begin_per_thread_vertex_triangle);
	arrangeIndex(thread_num, non_empty_cell_index_edge[0][0], non_empty_cell_index_begin_per_thread_edge);


	hash_cell_pair_num_prefix_vertex_triangle[0] = non_empty_cell_index_vertex_triangle[0][0];
	hash_cell_pair_num_prefix_edge[0] = non_empty_cell_index_edge[0][0];
	hash_cell_pair_num_prefix_vertex_triangle[1] = 0;
	hash_cell_pair_num_prefix_edge[1] = 0;

	thread->assignTask(this, OBTAIN_PAIR_COUNT);


	for (unsigned int i = 1; i < thread_num; ++i) {
		prefix_sum_record_per_thread_start_vertex_triangle[i] = prefix_sum_record_per_thread_start_vertex_triangle[i - 1] + hash_cell_pair_num_prefix_vertex_triangle[non_empty_cell_index_begin_per_thread_vertex_triangle[i] + 1];
		prefix_sum_record_per_thread_start_edge[i] = prefix_sum_record_per_thread_start_edge[i - 1] + hash_cell_pair_num_prefix_edge[non_empty_cell_index_begin_per_thread_edge[i] + 1];
	}

	thread->assignTask(this, SET_HASH_CELL_PAIR_NUM_PREFIX_SUM_TOGETHER);


	if (thread_num > 1) {
		thread->assignTask(this, SET_PAIR_AVE);
	}
	else {
		for (unsigned int i = 0; i < 2; ++i) {
			setPairAveInThread(i);
		}
	}

}



//RECORD_NONEMPTY_CELL
void SpatialHashing::recordNonEmptyCell(int thread_No)
{
	unsigned int cell_collider_triangle_index_size;
	unsigned int last_cell_index = cell_begin_per_thread[thread_No + 1];
	unsigned int start_index = cell_begin_per_thread[thread_No];
	//non_empty_cell_index_and_pair_vector_size[thread_No][1] = 0;
	unsigned int* non_empty_cell_index_triangle_ = non_empty_cell_index_vertex_triangle[thread_No] + 1;
	unsigned int* non_empty_cell_index_edge_ = non_empty_cell_index_edge[thread_No] + 1;
	bool contain_obj_tri;
	bool contain_obj_vertex;
	unsigned int thread_num_ = thread_num;
	bool has_collider_ = has_collider;

	for (unsigned int cell_index = start_index;
		cell_index < last_cell_index; ++cell_index) {
		contain_obj_vertex = false;
		for (unsigned int t = 0; t < thread_num_; ++t) {
			if (spatial_hashing_cell_vertex_size[t][cell_index]) {
				contain_obj_vertex = true;
				break;
			}
		}
		if (contain_obj_vertex) {
			*(non_empty_cell_index_triangle_++) = cell_index;
			continue;
		}
		else {
			if (has_collider_) {
				contain_obj_tri = false;
				for (unsigned int t = 0; t < thread_num_; ++t) {
					if (spatial_hashing_cell_triangle_size[t][cell_index]) {
						contain_obj_tri = true;						
						break;
					}
				}
				if (contain_obj_tri) {
					for (unsigned int t = 0; t < thread_num_; ++t) {
						if (spatial_hashing_cell_collider_vertex_size[t][cell_index]) {
							*(non_empty_cell_index_triangle_++) = cell_index;
							break;
						}
					}
				}
			}

		}
	}

	bool is_not_empty;
	unsigned int cell_triangle_index_size;
	for (unsigned int cell_index = start_index;
		cell_index < last_cell_index; ++cell_index) {
		cell_triangle_index_size = 0;
		is_not_empty = false;
		for (unsigned int t = 0; t < thread_num_; ++t) {
			if (spatial_hashing_cell_edge_size[t][cell_index]) {
				cell_triangle_index_size += spatial_hashing_cell_edge_size[t][cell_index];
			}
			if (cell_triangle_index_size > 2) {
				is_not_empty = true;
				break;
			}
		}
		if (is_not_empty) {
			*(non_empty_cell_index_edge_++) = cell_index;
			continue;
		}
		else {
			if (cell_triangle_index_size == 0) {
				continue;
			}
		}
		if (has_collider_) {
			for (unsigned int t = 0; t < thread_num_; ++t) {
				if (spatial_hashing_cell_collider_edge_size[t][cell_index]) {
					*(non_empty_cell_index_edge_++) = cell_index;
					break;
				}
			}
		}
	}
	non_empty_cell_index_edge[thread_No][0] = non_empty_cell_index_edge_ - non_empty_cell_index_edge[thread_No] - 1;
	non_empty_cell_index_vertex_triangle[thread_No][0] = non_empty_cell_index_triangle_ - non_empty_cell_index_vertex_triangle[thread_No] - 1;
}


//COMBINE_HASH_TABLE
void SpatialHashing::combineHashTable(int thread_No)
{
	unsigned int max_index_number_in_one_cell_ = max_index_number_in_one_cell;
	unsigned int max_index_number_in_one_cell_collider_ = max_index_number_in_one_cell_collider;
	unsigned int thread_num_ = thread_num;
	
	unsigned int ori_first_thread_num;

	unsigned int* spatial_hashing_cell_ = spatial_hashing_cell_triangle[0];
	unsigned int* spatial_hashing_cell_size_ = spatial_hashing_cell_triangle_size[0];

	unsigned int* spatial_hashing_cell_collider_ = spatial_hashing_cell_collider_triangle[0];
	unsigned int* spatial_hashing_cell_collider_size_ = spatial_hashing_cell_collider_triangle_size[0];

	bool has_collider_ = has_collider;

	unsigned int* spatial_hashing_cell_address;

	unsigned int start = hash_count_start_per_thread[thread_No];
	unsigned int last_cell_index = hash_count_start_per_thread[thread_No + 1];

	for (unsigned int cell_index = start; cell_index < last_cell_index; ++cell_index) {
		ori_first_thread_num = spatial_hashing_cell_size_[cell_index];
		spatial_hashing_cell_address = spatial_hashing_cell_ + max_index_number_in_one_cell_ * cell_index + ori_first_thread_num;
		for (unsigned int t = 1; t < thread_num_; ++t) {
			if (spatial_hashing_cell_triangle_size[t][cell_index]) {
				memcpy(spatial_hashing_cell_address, spatial_hashing_cell_triangle[t] + max_index_number_in_one_cell_ * cell_index, spatial_hashing_cell_triangle_size[t][cell_index] << 2);
				spatial_hashing_cell_address += spatial_hashing_cell_triangle_size[t][cell_index];
			}
		}
		spatial_hashing_cell_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_ - max_index_number_in_one_cell_ * cell_index;

		if (has_collider_) {
			ori_first_thread_num = spatial_hashing_cell_collider_size_[cell_index];
			spatial_hashing_cell_address = spatial_hashing_cell_collider_ + max_index_number_in_one_cell_collider_ * cell_index + ori_first_thread_num;
			for (unsigned int t = 1; t < thread_num_; ++t) {
				if (spatial_hashing_cell_collider_triangle_size[t][cell_index]) {
					memcpy(spatial_hashing_cell_address, spatial_hashing_cell_collider_triangle[t] + max_index_number_in_one_cell_collider_ * cell_index, spatial_hashing_cell_collider_triangle_size[t][cell_index] << 2);
					spatial_hashing_cell_address += spatial_hashing_cell_collider_triangle_size[t][cell_index];
				}
			}
			spatial_hashing_cell_collider_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_collider_ - max_index_number_in_one_cell_collider_ * cell_index;
		}
	}


	//edge
	max_index_number_in_one_cell_ = max_index_number_in_one_cell_edge;
	spatial_hashing_cell_ = spatial_hashing_cell_edge[0];
	spatial_hashing_cell_size_ = spatial_hashing_cell_edge_size[0];

	for (unsigned int cell_index = start; cell_index < last_cell_index; ++cell_index) {
		ori_first_thread_num = spatial_hashing_cell_size_[cell_index];
		spatial_hashing_cell_address = spatial_hashing_cell_ + max_index_number_in_one_cell_ * cell_index + ori_first_thread_num;
		for (unsigned int t = 1; t < thread_num_; ++t) {
			if (spatial_hashing_cell_edge_size[t][cell_index]) {
				memcpy(spatial_hashing_cell_address, spatial_hashing_cell_edge[t] + max_index_number_in_one_cell_ * cell_index, spatial_hashing_cell_edge_size[t][cell_index] << 2);
				spatial_hashing_cell_address += spatial_hashing_cell_edge_size[t][cell_index];
			}
		}
		spatial_hashing_cell_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_ - max_index_number_in_one_cell_ * cell_index;

	}


}



//OBTAIN_PAIR_COUNT
void SpatialHashing::obtainPairCount(int thread_No)
{
	unsigned int max_index_number_in_one_cell_ = max_index_number_in_one_cell;
	unsigned int max_index_number_in_one_cell_collider_ = max_index_number_in_one_cell_collider;

	unsigned int max_index_number_in_one_cell_vertex_ = max_index_number_in_one_cell_vertex;
	unsigned int max_index_number_in_one_cell_collider_vertex_= max_index_number_in_one_cell_collider_vertex;

	unsigned int thread_num_ = thread_num;
	unsigned int* hash_cell_pair_num_ = hash_cell_pair_num_prefix_vertex_triangle + 2;
	unsigned int last_cell_index = non_empty_cell_index_begin_per_thread_vertex_triangle[thread_No + 1];
	unsigned int start = non_empty_cell_index_begin_per_thread_vertex_triangle[thread_No];
	unsigned int* non_empty_cell_index_ = non_empty_cell_index_vertex_triangle[0] + 1;
	//unsigned int pair_num_with_collider = 0;
	unsigned int ori_first_thread_num;

	unsigned int* spatial_hashing_cell_ = spatial_hashing_cell_triangle[0];
	unsigned int* spatial_hashing_cell_size_ = spatial_hashing_cell_triangle_size[0];

	unsigned int* spatial_hashing_cell_vertex_ = spatial_hashing_cell_vertex[0];
	unsigned int* spatial_hashing_cell_vertex_size_ = spatial_hashing_cell_vertex_size[0];

	unsigned int* spatial_hashing_cell_collider_ = spatial_hashing_cell_collider_triangle[0];
	unsigned int* spatial_hashing_cell_collider_size_ = spatial_hashing_cell_collider_triangle_size[0];

	unsigned int* spatial_hashing_cell_collider_vertex_ = spatial_hashing_cell_collider_vertex[0];
	unsigned int* spatial_hashing_cell_collider_vertex_size_ = spatial_hashing_cell_collider_vertex_size[0];


	bool has_collider_ = has_collider;
	unsigned int cell_index;





	unsigned int* spatial_hashing_cell_address;
	unsigned int pair_num_cell;
	for (unsigned int i = start; i < last_cell_index; ++i) {
		cell_index = non_empty_cell_index_[i];
		pair_num_cell = 0;
		ori_first_thread_num = spatial_hashing_cell_size_[cell_index];
		spatial_hashing_cell_address = spatial_hashing_cell_ + max_index_number_in_one_cell_ * cell_index + ori_first_thread_num;
		for (unsigned int t = 1; t < thread_num_; ++t) {
			if (spatial_hashing_cell_triangle_size[t][cell_index]) {
				memcpy(spatial_hashing_cell_address, spatial_hashing_cell_triangle[t] + max_index_number_in_one_cell_ * cell_index, spatial_hashing_cell_triangle_size[t][cell_index] << 2);
				spatial_hashing_cell_address += spatial_hashing_cell_triangle_size[t][cell_index];
			}
		}
		spatial_hashing_cell_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_- max_index_number_in_one_cell_* cell_index;

		ori_first_thread_num = spatial_hashing_cell_vertex_size_[cell_index];
		spatial_hashing_cell_address = spatial_hashing_cell_vertex_+ max_index_number_in_one_cell_vertex_ *cell_index + ori_first_thread_num;

		for (unsigned int t = 1; t < thread_num_; ++t) {
			if (spatial_hashing_cell_vertex_size[t][cell_index]) {
				memcpy(spatial_hashing_cell_address, spatial_hashing_cell_vertex[t]+ max_index_number_in_one_cell_vertex_ *cell_index, spatial_hashing_cell_vertex_size[t][cell_index] << 2);
				spatial_hashing_cell_address += spatial_hashing_cell_vertex_size[t][cell_index];
			}
		}
		spatial_hashing_cell_vertex_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_vertex_- max_index_number_in_one_cell_vertex_ * cell_index;

		pair_num_cell = (spatial_hashing_cell_size_[cell_index] * spatial_hashing_cell_vertex_size_[cell_index]) >> 2;
		
		if (has_collider_) {
			ori_first_thread_num = spatial_hashing_cell_collider_size_[cell_index];
			spatial_hashing_cell_address = spatial_hashing_cell_collider_+ max_index_number_in_one_cell_collider_ * cell_index + ori_first_thread_num;
			for (unsigned int t = 1; t < thread_num_; ++t) {
				if (spatial_hashing_cell_collider_triangle_size[t][cell_index]) {
					memcpy(spatial_hashing_cell_address, spatial_hashing_cell_collider_triangle[t]+ max_index_number_in_one_cell_collider_ * cell_index, spatial_hashing_cell_collider_triangle_size[t][cell_index] << 2);
					spatial_hashing_cell_address += spatial_hashing_cell_collider_triangle_size[t][cell_index];
				}
			}
			spatial_hashing_cell_collider_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_collider_- max_index_number_in_one_cell_collider_*cell_index;

			ori_first_thread_num = spatial_hashing_cell_collider_vertex_size_[cell_index];
			spatial_hashing_cell_address = spatial_hashing_cell_collider_vertex_ + max_index_number_in_one_cell_collider_vertex_ * cell_index + ori_first_thread_num;
			for (unsigned int t = 1; t < thread_num_; ++t) {
				if (spatial_hashing_cell_collider_vertex_size[t][cell_index]) {
					memcpy(spatial_hashing_cell_address, spatial_hashing_cell_collider_vertex[t]+ max_index_number_in_one_cell_collider_vertex_*cell_index, spatial_hashing_cell_collider_vertex_size[t][cell_index] << 2);
					spatial_hashing_cell_address += spatial_hashing_cell_collider_vertex_size[t][cell_index];
				}
			}
			spatial_hashing_cell_collider_vertex_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_collider_vertex_- max_index_number_in_one_cell_collider_vertex_*cell_index;

			
			pair_num_cell += ((spatial_hashing_cell_size_[cell_index] * spatial_hashing_cell_collider_vertex_size_[cell_index]
				+ spatial_hashing_cell_vertex_size_[cell_index] * spatial_hashing_cell_collider_size_[cell_index]) >> 2);
		}

		if (i == start) {
			hash_cell_pair_num_[i] = pair_num_cell;
		}
		else {
			hash_cell_pair_num_[i] = hash_cell_pair_num_[i - 1] + pair_num_cell;
		}
	}

	hash_cell_pair_num_ = hash_cell_pair_num_prefix_edge + 2;
	last_cell_index = non_empty_cell_index_begin_per_thread_edge[thread_No + 1];
	start = non_empty_cell_index_begin_per_thread_edge[thread_No];
	non_empty_cell_index_ = non_empty_cell_index_edge[0] + 1;
	
	spatial_hashing_cell_ = spatial_hashing_cell_edge[0];
	spatial_hashing_cell_size_ = spatial_hashing_cell_edge_size[0];
	spatial_hashing_cell_collider_ = spatial_hashing_cell_collider_edge[0];
	spatial_hashing_cell_collider_size_ = spatial_hashing_cell_collider_edge_size[0];

	max_index_number_in_one_cell_ = max_index_number_in_one_cell_edge;
	max_index_number_in_one_cell_collider_ = max_index_number_in_one_cell_collider_edge;



	for (unsigned int i = start; i < last_cell_index; ++i) {
		cell_index = non_empty_cell_index_[i];
		pair_num_cell = 0;
		ori_first_thread_num = spatial_hashing_cell_size_[cell_index];
		spatial_hashing_cell_address = spatial_hashing_cell_+ max_index_number_in_one_cell_ *cell_index + ori_first_thread_num;
		for (unsigned int t = 1; t < thread_num_; ++t) {
			if (spatial_hashing_cell_edge_size[t][cell_index]) {
				memcpy(spatial_hashing_cell_address, spatial_hashing_cell_edge[t]+ max_index_number_in_one_cell_ * cell_index, spatial_hashing_cell_edge_size[t][cell_index] << 2);
				spatial_hashing_cell_address += spatial_hashing_cell_edge_size[t][cell_index];
			}
		}
		spatial_hashing_cell_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_ - max_index_number_in_one_cell_*cell_index;
		if (spatial_hashing_cell_size_[cell_index] > 2) {
			pair_num_cell = (spatial_hashing_cell_size_[cell_index] * (spatial_hashing_cell_size_[cell_index] - 2)) >> 3;
		}
		if (has_collider_) {
			ori_first_thread_num = spatial_hashing_cell_collider_size_[cell_index];
			spatial_hashing_cell_address = spatial_hashing_cell_collider_+ max_index_number_in_one_cell_collider_*cell_index + ori_first_thread_num;
			for (unsigned int t = 1; t < thread_num_; ++t) {
				if (spatial_hashing_cell_collider_edge_size[t][cell_index]) {
					memcpy(spatial_hashing_cell_address, spatial_hashing_cell_collider_edge[t]+ max_index_number_in_one_cell_collider_ *cell_index, spatial_hashing_cell_collider_edge_size[t][cell_index] << 2);
					spatial_hashing_cell_address += spatial_hashing_cell_collider_edge_size[t][cell_index];
				}
			}
			spatial_hashing_cell_collider_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_collider_- max_index_number_in_one_cell_collider_*cell_index;
			if (spatial_hashing_cell_collider_size_[cell_index]) {
				pair_num_cell += (spatial_hashing_cell_size_[cell_index] * spatial_hashing_cell_collider_size_[cell_index]) >> 2;
			}
		}

		if (i == start) {
			hash_cell_pair_num_[i] = pair_num_cell;
		}
		else {
			hash_cell_pair_num_[i] = hash_cell_pair_num_[i - 1] + pair_num_cell;
		}
	}
}


//SET_HASH_CELL_PAIR_NUM_PREFIX_SUM_TOGETHER
void SpatialHashing::setHashCellPairNumPrefixSumTogether(int thread_No)
{
	if (thread_No == 0) {
		return;
	}
	unsigned int start = non_empty_cell_index_begin_per_thread_vertex_triangle[thread_No] + 2;
	unsigned int end = non_empty_cell_index_begin_per_thread_vertex_triangle[thread_No + 1] + 2;
	unsigned int add_last_thread_pair_num = prefix_sum_record_per_thread_start_vertex_triangle[thread_No];
	for (unsigned int i = start; i < end; ++i) {
		hash_cell_pair_num_prefix_vertex_triangle[i] += add_last_thread_pair_num;
	}

	start = non_empty_cell_index_begin_per_thread_edge[thread_No] + 2;
	end = non_empty_cell_index_begin_per_thread_edge[thread_No + 1] + 2;
	add_last_thread_pair_num = prefix_sum_record_per_thread_start_edge[thread_No];
	for (unsigned int i = start; i < end; ++i) {
		hash_cell_pair_num_prefix_edge[i] += add_last_thread_pair_num;
	}
}


//SET_PAIR_AVE
void SpatialHashing::setPairAveInThread(int thread_No)
{
	if (thread_No == 0) {
		setPairAveInThreadEdgeEdge();
	}
	else if (thread_No == 1) {
		setPairAveInThreadVertexTriangle();
	}
}


void SpatialHashing::setPairAveInThreadEdgeEdge()
{
	memset(non_empty_cell_index_ave_begin_per_thread_edge, 0, 4 * (thread_num + 1));
	std::vector<unsigned int>prefix_sum_pair_per_thread(thread_num + 1, 0);
	for (unsigned int i = 0; i < thread_num; ++i) {
		prefix_sum_pair_per_thread[i + 1] = hash_cell_pair_num_prefix_edge[1 + non_empty_cell_index_begin_per_thread_edge[i + 1]];//  prefix_sum_pair_per_thread[i] + pair_count_per_thread[i] + pair_count_with_collider_per_thread[i];
	}
	std::vector<unsigned int> ave_pair_num_start_per_thread(thread_num + 1);
	arrangeIndex(thread_num, prefix_sum_pair_per_thread[thread_num], ave_pair_num_start_per_thread.data());
	unsigned int record_prefix_sum_thread = 0;
	unsigned int start_index = 0;
	unsigned int pair_count = 0;
	unsigned int total_cell_num = non_empty_cell_index_edge[0][0];
	//unsigned int pre_pair_count = 0;
	unsigned int*  hash_cell_pair_num_prefix_start = hash_cell_pair_num_prefix_edge + 2;

	for (unsigned int t = 1; t <= thread_num; ++t) {
		for (unsigned int i = record_prefix_sum_thread; i < thread_num; ++i) {
			if (prefix_sum_pair_per_thread[i + 1] > ave_pair_num_start_per_thread[t]) {
				record_prefix_sum_thread = i;
				break;
			}
		}
		if (start_index < non_empty_cell_index_begin_per_thread_edge[record_prefix_sum_thread]) {
			start_index = non_empty_cell_index_begin_per_thread_edge[record_prefix_sum_thread];
			pair_count = prefix_sum_pair_per_thread[record_prefix_sum_thread];
		}
		for (unsigned int i = start_index; i < total_cell_num; ++i) {
			if (hash_cell_pair_num_prefix_start[i] >= ave_pair_num_start_per_thread[t]) {
				pair_count = hash_cell_pair_num_prefix_start[i];
				non_empty_cell_index_ave_begin_per_thread_edge[t] = i + 1;
				start_index = i + 1;
				break;
			}
		}
	}

	for (unsigned int i = 0; i < thread_num; ++i) {
		global_cell_start_edge[i << 1] = non_empty_cell_index_edge[0][non_empty_cell_index_ave_begin_per_thread_edge[i] + 1];
		global_cell_start_edge[(i << 1) + 1] = non_empty_cell_index_edge[0][non_empty_cell_index_ave_begin_per_thread_edge[i + 1]] + 1;
	}
}


void SpatialHashing::setPairAveInThreadVertexTriangle()
{
	memset(non_empty_cell_index_ave_begin_per_thread_vertex_triangle, 0, 4 * (thread_num + 1));
	std::vector<unsigned int>prefix_sum_pair_per_thread(thread_num + 1, 0);
	for (unsigned int i = 0; i < thread_num; ++i) {
		prefix_sum_pair_per_thread[i + 1] = hash_cell_pair_num_prefix_vertex_triangle[1 + non_empty_cell_index_begin_per_thread_vertex_triangle[i + 1]];//  prefix_sum_pair_per_thread[i] + pair_count_per_thread[i] + pair_count_with_collider_per_thread[i];
	}
	std::vector<unsigned int> ave_pair_num_start_per_thread(thread_num + 1);
	arrangeIndex(thread_num, prefix_sum_pair_per_thread[thread_num], ave_pair_num_start_per_thread.data());
	unsigned int record_prefix_sum_thread = 0;
	unsigned int start_index = 0;
	unsigned int pair_count = 0;
	unsigned int total_cell_num = non_empty_cell_index_vertex_triangle[0][0];
	//unsigned int pre_pair_count = 0;
	unsigned int* hash_cell_pair_num_prefix_start = hash_cell_pair_num_prefix_vertex_triangle + 2;

	for (unsigned int t = 1; t <= thread_num; ++t) {
		for (unsigned int i = record_prefix_sum_thread; i < thread_num; ++i) {
			if (prefix_sum_pair_per_thread[i + 1] > ave_pair_num_start_per_thread[t]) {
				record_prefix_sum_thread = i;
				break;
			}
		}
		if (start_index < non_empty_cell_index_begin_per_thread_vertex_triangle[record_prefix_sum_thread]) {
			start_index = non_empty_cell_index_begin_per_thread_vertex_triangle[record_prefix_sum_thread];
			pair_count = prefix_sum_pair_per_thread[record_prefix_sum_thread];
		}
		for (unsigned int i = start_index; i < total_cell_num; ++i) {
			if (hash_cell_pair_num_prefix_start[i] >= ave_pair_num_start_per_thread[t]) {
				pair_count = hash_cell_pair_num_prefix_start[i];
				non_empty_cell_index_ave_begin_per_thread_vertex_triangle[t] = i + 1;
				start_index = i + 1;
				break;
			}
		}
	}


	for (unsigned int i = 0; i < thread_num; ++i) {
		global_cell_start_vertex_triangle[i << 1] = non_empty_cell_index_vertex_triangle[0][non_empty_cell_index_ave_begin_per_thread_vertex_triangle[i] + 1];
		global_cell_start_vertex_triangle[(i << 1) + 1] = non_empty_cell_index_vertex_triangle[0][non_empty_cell_index_ave_begin_per_thread_vertex_triangle[i + 1]] + 1;
	}
}


void SpatialHashing::selectCell(double* start_pos, double* dir, std::vector<unsigned int>& select_hash_index, std::vector<unsigned int>&	initial_hash_index)
{
	double start_position[3];
	double end[3];
	memcpy(start_position, start_pos, 24);
	select_hash_index.clear();
	initial_hash_index.clear();
	if (intersectSpatialHasingCube(scene_aabb, start_position, dir, end)) {
		collectAllVoxel(select_hash_index, initial_hash_index, start_position, end, scene_aabb);
	}
}

void SpatialHashing::findAllElementsInOneCell(unsigned int select_hash_index, unsigned int ori_hash_index,
	std::vector<std::vector<unsigned int>>& vertex_index, std::vector<std::vector<unsigned int>>& triangle_index,
	std::vector<std::vector<unsigned int>>& edge_index)
{
	double aabb[6];
	unsigned int index[3];
	index[2] = ori_hash_index / cell_num0_cell_num1;
	index[1] = (ori_hash_index % cell_num0_cell_num1) / cell_number[0];
	index[0] = (ori_hash_index % cell_num0_cell_num1) % cell_number[0];
	for (unsigned int i = 0; i < 3; ++i) {
		aabb[i] = cell_length * index[i] + scene_aabb[i];
		aabb[i + 3] = aabb[i] + cell_length;
	}
	unsigned int hash_location;
	for (unsigned int i = 0; i < spatial_hashing_cell_triangle_size[0][select_hash_index]; i += 2) {
		hash_location = select_hash_index * max_index_number_in_one_cell + i;
		if (AABB::AABB_intersection(aabb, obj_tri_aabb[spatial_hashing_cell_triangle[0][hash_location + 1]][spatial_hashing_cell_triangle[0][hash_location]].data())) {
			triangle_index[spatial_hashing_cell_triangle[0][hash_location + 1]].emplace_back(spatial_hashing_cell_triangle[0][hash_location]);
		}
	}
	for (unsigned int i = 0; i < spatial_hashing_cell_vertex_size[0][select_hash_index]; i += 2) {
		hash_location = select_hash_index * max_index_number_in_one_cell_vertex + i;
		if (AABB::AABB_intersection(aabb, obj_vertex_aabb[spatial_hashing_cell_vertex[0][hash_location + 1]][spatial_hashing_cell_vertex[0][hash_location]].data())) {
			vertex_index[spatial_hashing_cell_vertex[0][hash_location + 1]].emplace_back(spatial_hashing_cell_vertex[0][hash_location]);
			//for (unsigned int i = 0; i < 6; ++i) {
			//	std::cout << aabb[i] << " ";
			//}
			//std::cout << std::endl;
			//for (unsigned int i = 0; i < 6; ++i) {
			//	std::cout << obj_vertex_aabb[spatial_hashing_cell_vertex[0][hash_location + 1]][spatial_hashing_cell_vertex[0][hash_location]][i] << " ";
			//}
			//std::cout << std::endl;
		}
	}
	for (unsigned int i = 0; i < spatial_hashing_cell_edge_size[0][select_hash_index]; i += 2) {
		hash_location = select_hash_index * max_index_number_in_one_cell_edge + i;
		if (AABB::AABB_intersection(aabb, obj_edge_aabb[spatial_hashing_cell_edge[0][hash_location + 1]][spatial_hashing_cell_edge[0][hash_location]].data())) {
			edge_index[spatial_hashing_cell_edge[0][hash_location + 1]].emplace_back(spatial_hashing_cell_edge[0][hash_location]);
		}
	}
}


bool SpatialHashing::intersectSpatialHasingCube(double* AABB, double* start, double* direction, double* end)
{
	double t_xmin = (AABB[0] -start[0]) / direction[0];
	double t_xmax = (AABB[3] - start[0]) / direction[0];
	double t_ymin = (AABB[1] -start[1]) / direction[1];
	double t_ymax = (AABB[4] - start[1]) / direction[1];
	double t_zmin = (AABB[2] -start[2]) / direction[2];
	double t_zmax = (AABB[5] - start[2]) / direction[2];


	double t;
	if (direction[0] < 0)
	{
		t = t_xmax;
		t_xmax = t_xmin;
		t_xmin = t;
	}
	if (direction[1] < 0)
	{
		t = t_ymax;
		t_ymax = t_ymin;
		t_ymin = t;
	}
	if (direction[2] < 0)
	{
		t = t_zmax;
		t_zmax = t_zmin;
		t_zmin = t;
	}

	double t_enter = myMax(t_xmin, t_ymin);
	if (t_enter < t_zmin) {
		t_enter = t_zmin;
	}
	double t_exit = myMin(t_xmax, t_ymax);
	if (t_exit > t_zmax) {
		t_exit = t_zmax;
	}
	if (t_enter <= t_exit && t_exit >= 0) {
		end[0] = start[0] + direction[0] * t_exit;
		end[1] = start[1] + direction[1] * t_exit;
		end[2] = start[2] + direction[2] * t_exit;
		if (t_enter > 0) {//the start point(camera pos) is in the spatial cube
			start[0] += direction[0] * t_enter;
			start[1] += direction[1] * t_enter;
			start[2] += direction[2] * t_enter;
		}
		return true;
	}
	else {
		return false;
	}
}


void SpatialHashing::collectAllVoxel(std::vector<unsigned int>& select_hash_index, std::vector<unsigned int>&
	initial_hash_index, double start__[3], double end__[3], double* scene_aabb)
//collect all voxels that possible be selected by the mouse
{
	double start[3];
	double end[3];
	int d[3];
	double tDelta[3];
	double tMax[3];
	int voxel[3];
	unsigned int hash_value;
	unsigned int ori_hash_value;

	double start_[3];
	double end_[3];
	SUB(start_, start__, scene_aabb);
	SUB(end_, end__, scene_aabb);

	//int around[8][3] = { 1,1,1,1,1,-1,
	//	1,-1,1,1,-1,-1,
	//	-1,1,1,-1,1,-1,
	//	-1,-1,1,-1,-1,-1 };
	start[0] = start_[0] / cell_length;
	start[1] = start_[1] / cell_length;
	start[2] = start_[2] / cell_length;
	end[0] = end_[0] / cell_length;
	end[1] = end_[1] / cell_length;
	end[2] = end_[2] / cell_length;
	d[0] = SIGN(end[0] - start[0]);
	d[1] = SIGN(end[1] - start[1]);
	d[2] = SIGN(end[2] - start[2]);
	if (d[0] != 0) {
		tDelta[0] = (double)d[0] / (end[0] - start[0]);
	}
	else tDelta[0] = 1e7;
	if (d[0] > 0) {
		tMax[0] = tDelta[0] * frac1(start[0]);
	}
	else tMax[0] = tDelta[0] * frac0(start[0]);
	voxel[0] = (int)floor(start[0]);
	if (d[1] != 0) {
		tDelta[1] = d[1] / (end[1] - start[1]);
	}
	else tDelta[1] = 1e7;
	if (d[1] > 0) {
		tMax[1] = tDelta[1] * frac1(start[1]);
	}
	else tMax[1] = tDelta[1] * frac0(start[1]);
	voxel[1] = (int)floor(start[1]);
	if (d[2] != 0) {
		tDelta[2] = d[2] / (end[2] - start[2]);
	}
	else tDelta[2] = 1e7;
	if (d[2] > 0) {
		tMax[2] = tDelta[2] * frac1(start[2]);
	}
	else tMax[2] = tDelta[2] * frac0(start[2]);
	voxel[2] = (int)floor(start[2]);					
	hash_value = ((voxel[0] * P1) ^ (P2 * voxel[1]) ^ (voxel[2] * P3)) % hash_cell_count;
	ori_hash_value = voxel[1] * cell_number[0] + cell_num0_cell_num1 * voxel[2] + voxel[0];
	select_hash_index.emplace_back(hash_value);
	initial_hash_index.emplace_back(ori_hash_value);
	while (true) {
		if (tMax[0] > 1 && tMax[1] > 1 && tMax[2] > 1) {
			break;
		}
		if (tMax[0] < tMax[1]) {
			if (tMax[0] < tMax[2]) {
				voxel[0] += d[0];
				tMax[0] += tDelta[0];
			}
			else {
				voxel[2] += d[2];
				tMax[2] += tDelta[2];
			}
		}
		else {
			if (tMax[1] < tMax[2]) {
				voxel[1] += d[1];
				tMax[1] += tDelta[1];
			}
			else {
				voxel[2] += d[2];
				tMax[2] += tDelta[2];
			}
		}
		if (voxel[0] > -1 && voxel[0] < cell_number[0] &&
			voxel[1] > -1 && voxel[1] < cell_number[1] &&
			voxel[2] > -1 && voxel[2] < cell_number[2]) {
			hash_value = ((voxel[0] * P1) ^ (P2 * voxel[1]) ^ (voxel[2] * P3)) % hash_cell_count;
			ori_hash_value = voxel[1] * cell_number[0] + cell_num0_cell_num1 * voxel[2] + voxel[0];
			select_hash_index.emplace_back(hash_value);
			initial_hash_index.emplace_back(ori_hash_value);
		}
	}
}

inline double SpatialHashing::frac0(double x)
{
	return x - floor(x);
}
inline double SpatialHashing::frac1(double x)
{
	return 1.0 - x + floor(x);
}