#include"spatial_hashing.h"
//#include"../basic/write_txt.h"


void SpatialHashing::initialHashCellLength(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, double& cell_length, double* tolerance_ratio)
{
	double max_length = 0;
	double ave_length = 0;
	int edge_num = 0;
	std::vector<MeshStruct::Edge>* edge;
	for (int i = 0; i < cloth->size(); ++i) {
		edge = &cloth->data()[i].mesh_struct.edges;
		for (int j = 0; j < edge->size(); ++j) {
			if (max_length < (*edge)[j].length) {
				max_length = (*edge)[j].length;
			}
			ave_length += (*edge)[j].length;
		}
		edge_num += edge->size();
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		edge = &tetrahedron->data()[i].mesh_struct.edges;
		for (int j = 0; j < edge->size(); ++j) {
			if (max_length < (*edge)[j].length) {
				max_length = (*edge)[j].length;
			}
			ave_length += (*edge)[j].length;
		}
		edge_num += edge->size();
	}

	ave_length /= (double)edge_num;
	cell_length =1.0* max_length + 2.0 * tolerance_ratio[SELF_POINT_TRIANGLE] * ave_length;
	//cell_length = 1.0 * max_length +2.0 * tolerance_ratio[SELF_POINT_TRIANGLE] * ave_length;
	std::cout << "tolerance__ " << tolerance_ratio[SELF_POINT_TRIANGLE] << std::endl;
	std::cout << "ave_length " << ave_length << " max length " << max_length << " ratio " << (double)max_length / (double)ave_length << std::endl;
	triangle_pair_number_thread.resize(thread_num);
}



void SpatialHashing::initialHashCell(unsigned int total_triangle_num, unsigned int max_index_number_in_one_cell,
	unsigned int max_index_number_in_one_cell_collider)
{
	//indicator = new std::vector<unsigned int>[thread_num];
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	indicator[i].resize(1);
	//}
	hash_cell_count = findNeareastPrimeNumber(0.5 * total_triangle_num);// 90191;// 24877;// 500009; //;// // ;// 32999;//  // 11987;//84947
	//hash_cell_count = 90191;
	std::cout << "total triangle " << total_triangle_num << " hash cell count " << hash_cell_count << std::endl;
	this->max_index_number_in_one_cell = max_index_number_in_one_cell;// = 801;
	this->max_index_number_in_one_cell_collider = max_index_number_in_one_cell_collider;// = 401;
	vertex_triangle_pair = new unsigned int* [thread_num];
	vertex_obj_triangle_collider_pair = new unsigned int* [thread_num];
	vertex_collider_triangle_obj_pair = new unsigned int* [thread_num];
	edge_edge_pair = new unsigned int* [thread_num];
	edge_edge_pair_collider = new unsigned int* [thread_num];
	spatial_hashing_cell_triangle = new unsigned int** [thread_num];
	spatial_hashing_cell_edge = new unsigned int** [thread_num];
	spatial_hashing_cell_vertex = new unsigned int** [thread_num];
	spatial_hashing_cell_collider_triangle = new unsigned int** [thread_num];
	spatial_hashing_cell_collider_edge = new unsigned int** [thread_num];
	spatial_hashing_cell_collider_vertex = new unsigned int** [thread_num];

	for (unsigned int i = 0; i < thread_num; ++i) {
		edge_edge_pair[i] = new unsigned int[max_index_number_in_one_cell * total_triangle_num];// 
		memset(edge_edge_pair[i], 0, 4 * max_index_number_in_one_cell * total_triangle_num);

		vertex_triangle_pair[i] = new unsigned int[max_index_number_in_one_cell * total_triangle_num / 2];// 
		memset(vertex_triangle_pair[i], 0, 4 * (max_index_number_in_one_cell * total_triangle_num / 2));

		spatial_hashing_cell_triangle[i] = new unsigned int* [hash_cell_count];
		spatial_hashing_cell_edge[i] = new unsigned int* [hash_cell_count];
		spatial_hashing_cell_vertex[i] = new unsigned int* [hash_cell_count];

		for (unsigned int j = 0; j < hash_cell_count; ++j) {
			spatial_hashing_cell_triangle[i][j] = new unsigned int[max_index_number_in_one_cell];
			memset(spatial_hashing_cell_triangle[i][j], 0, 4 * max_index_number_in_one_cell);

			spatial_hashing_cell_edge[i][j] = new unsigned int[ max_index_number_in_one_cell];
			memset(spatial_hashing_cell_edge[i][j], 0, 4 * max_index_number_in_one_cell);

			spatial_hashing_cell_vertex[i][j] = new unsigned int[max_index_number_in_one_cell/2+1];
			memset(spatial_hashing_cell_vertex[i][j], 0, (max_index_number_in_one_cell / 2 + 1) * 4);
		}

		if (has_collider) {
			edge_edge_pair_collider[i] = new unsigned int[max_index_number_in_one_cell * total_triangle_num];
			memset(edge_edge_pair_collider[i], 0, 4 * max_index_number_in_one_cell * total_triangle_num);

			vertex_obj_triangle_collider_pair[i] = new unsigned int[max_index_number_in_one_cell * total_triangle_num / 2];
			memset(vertex_obj_triangle_collider_pair[i], 0, 4 * (max_index_number_in_one_cell * total_triangle_num / 2));

			vertex_collider_triangle_obj_pair[i] = new unsigned int[max_index_number_in_one_cell * total_triangle_num / 2];
			memset(vertex_collider_triangle_obj_pair[i], 0, 4 * (max_index_number_in_one_cell * total_triangle_num / 2));


			spatial_hashing_cell_collider_triangle[i] = new unsigned int* [hash_cell_count];
			spatial_hashing_cell_collider_edge[i] = new unsigned int* [hash_cell_count];
			spatial_hashing_cell_collider_vertex[i] = new unsigned int* [hash_cell_count];
			for (unsigned int j = 0; j < hash_cell_count; ++j) {
				spatial_hashing_cell_collider_triangle[i][j] = new unsigned int[max_index_number_in_one_cell_collider];
				memset(spatial_hashing_cell_collider_triangle[i][j], 0, 4 * max_index_number_in_one_cell_collider);

				spatial_hashing_cell_collider_edge[i][j] = new unsigned int[max_index_number_in_one_cell_collider];
				memset(spatial_hashing_cell_collider_edge[i][j], 0, 4 * max_index_number_in_one_cell_collider);

				spatial_hashing_cell_collider_vertex[i][j] = new unsigned int[max_index_number_in_one_cell_collider/2+1];
				memset(spatial_hashing_cell_collider_vertex[i][j], 0, 4 * (max_index_number_in_one_cell_collider/2+1));

			}
		}
		else {
			edge_edge_pair_collider[i] = new unsigned int[1];
			vertex_obj_triangle_collider_pair[i] = new unsigned int[1];
			vertex_collider_triangle_obj_pair[i] = new unsigned int[1];
		}
	}

	spatial_hashing_cell_triangle_size = new unsigned int* [thread_num];
	spatial_hashing_cell_edge_size = new unsigned int* [thread_num];
	spatial_hashing_cell_vertex_size = new unsigned int* [thread_num];
	spatial_hashing_cell_collider_triangle_size = new unsigned int* [thread_num];
	spatial_hashing_cell_collider_edge_size = new unsigned int* [thread_num];
	spatial_hashing_cell_collider_vertex_size = new unsigned int* [thread_num];

	for (unsigned int i = 0; i < thread_num; ++i) {
		spatial_hashing_cell_triangle_size[i] = new unsigned int[hash_cell_count];
		spatial_hashing_cell_collider_triangle_size[i] = new unsigned int[hash_cell_count];
		memset(spatial_hashing_cell_triangle_size[i], 0, hash_cell_count << 2);
		memset(spatial_hashing_cell_collider_triangle_size[i], 0, hash_cell_count << 2);


		spatial_hashing_cell_edge_size[i] = new unsigned int[hash_cell_count];
		spatial_hashing_cell_collider_edge_size[i] = new unsigned int[hash_cell_count];
		memset(spatial_hashing_cell_edge_size[i], 0, hash_cell_count << 2);
		memset(spatial_hashing_cell_collider_edge_size[i], 0, hash_cell_count << 2);

		spatial_hashing_cell_vertex_size[i] = new unsigned int[hash_cell_count];
		spatial_hashing_cell_collider_vertex_size[i] = new unsigned int[hash_cell_count];
		memset(spatial_hashing_cell_vertex_size[i], 0, hash_cell_count << 2);
		memset(spatial_hashing_cell_collider_vertex_size[i], 0, hash_cell_count << 2);
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



}

void SpatialHashing::findListOfSpatialHashingCell()
{
	//std::vector<unsigned int*> index_per_thread(thread_num);
	//std::vector<unsigned int> start_index(thread_num, 0);
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	index_per_thread[i] = spatial_hashing_value[i];
	//}
	//unsigned int current_hash_index = UINT_MAX;
	//bool end_loop;
	//hash_value_for_test.clear();
	//while (true)
	//{
	//	end_loop = (index_per_thread[0] == actual_hash_value_end_index_ref[0]);
	//	for (unsigned int i = 1; i < thread_num; ++i) {
	//		end_loop = end_loop && (index_per_thread[i] == actual_hash_value_end_index_ref[i]);
	//	}
	//	if (end_loop) {
	//		break;
	//	}
	//	current_hash_index = UINT_MAX;
	//	for (unsigned int i = 0; i < thread_num; ++i) {
	//		if (index_per_thread[i] != actual_hash_value_end_index_ref[i]) {
	//			if (current_hash_index > *(index_per_thread[i])) {
	//				current_hash_index = *(index_per_thread[i]);
	//			}
	//		}
	//	}
	//	hash_value_for_test.emplace_back(current_hash_index);
	//	for (unsigned int i = 0; i < thread_num; ++i) {
	//		if (index_per_thread[i] != actual_hash_value_end_index_ref[i]) {
	//			if (*(index_per_thread[i]) == current_hash_index)
	//			{
	//				index_per_thread[i] += hash_index_count_per_thread[i][start_index[i]];
	//				start_index[i]++;
	//			}
	//		}
	//	}
	//}
}


void SpatialHashing::setInObject(std::vector<Cloth>* cloth, std::vector<Collider>* collider,
	std::vector<Tetrahedron>* tetrahedron, Thread* thread, double* tolerance_ratio,
	unsigned int max_cell_count, bool for_construct_patch,
	unsigned int max_index_number_in_one_cell,
	unsigned int max_index_number_in_one_cell_collider)
{
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
	triangle_begin_per_obj_per_thread.resize(thread_num);
	triangle_begin_per_obj_per_thread_collider.resize(thread_num);
	collider_begin_obj_index = cloth->size() + tetrahedron->size();
	for (unsigned int i = 0; i < thread_num; ++i) {
		triangle_begin_per_obj_per_thread[i].resize(collider_begin_obj_index);
		triangle_begin_per_obj_per_thread_collider[i].resize(collider->size());
	}
	for (unsigned int j = 0; j < thread_num; ++j) {
		unsigned int total_triangle_num = 0;
		unsigned int total_triangle_num_collider = 0;
		for (unsigned int i = 0; i < cloth->size(); ++i) {
			triangle_begin_per_obj_per_thread[j][i] = total_triangle_num;
			total_triangle_num += cloth->data()[i].mesh_struct.face_index_begin_per_thread[j + 1] - cloth->data()[i].mesh_struct.face_index_begin_per_thread[j];
		}
		for (int i = 0; i < tetrahedron->size(); ++i) {
			triangle_begin_per_obj_per_thread[j][i + tetrahedron_begin_obj_index] = total_triangle_num;
			total_triangle_num += tetrahedron->data()[i].mesh_struct.face_index_begin_per_thread[j + 1] - tetrahedron->data()[i].mesh_struct.face_index_begin_per_thread[j];
		}
		for (int i = 0; i < collider->size(); ++i) {
			triangle_begin_per_obj_per_thread_collider[j][i] = total_triangle_num_collider;
			total_triangle_num_collider += collider->data()[i].mesh_struct.face_index_begin_per_thread[j + 1] - collider->data()[i].mesh_struct.face_index_begin_per_thread[j];
		}
	}
	for (unsigned int j = 0; j < thread_num; ++j) {
		for (unsigned int i = 1; i < collider_begin_obj_index; ++i) {
			triangle_begin_per_obj_per_thread[j][i] += triangle_begin_per_obj_per_thread[j][i - 1];
		}
		for (unsigned int i = 1; i < collider->size(); ++i) {
			triangle_begin_per_obj_per_thread_collider[j][i] += triangle_begin_per_obj_per_thread_collider[j][i - 1];
		}
		for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
			if (i < cloth->size()) {
				triangle_begin_per_obj_per_thread[j][i] -= cloth->data()[i].mesh_struct.face_index_begin_per_thread[j];
			}
			else {
				triangle_begin_per_obj_per_thread[j][i] -= tetrahedron->data()[i - cloth->size()].mesh_struct.face_index_begin_per_thread[j];
			}
			//std::cout <<"could be negative "<< triangle_begin_per_obj_per_thread[j][i] << std::endl;
		}
		for (unsigned int i = 0; i < collider->size(); ++i) {
			triangle_begin_per_obj_per_thread_collider[j][i] -= collider->data()[i].mesh_struct.face_index_begin_per_thread[j];
		}
	}

	this->max_cell_count = max_cell_count;

	initialHashCellLength(cloth, tetrahedron, cell_length, tolerance_ratio);


	if (!for_construct_patch) {
		//actual_hash_count_start_per_thread = new unsigned int[thread_num + 1];
		//total_hash_count_start_per_thread = new unsigned int[thread_num + 1];
		//prefix_sum_thread_start = new unsigned int[thread_num];
		//memset(prefix_sum_thread_start, 0, 4 * thread_num);

		//initialParallePrefix();

		//prefix_sum.resize(thread_num);
		//prefix_sum_collider.resize(thread_num);
		//hash_index_count_per_thread.resize(thread_num);
		//hash_index_count_per_thread_collider.resize(thread_num);
		//for (unsigned int i = 0; i < thread_num; ++i) {
		//	prefix_sum[i].reserve(100 * 100 * 100 / thread_num);
		//	if (has_collider) {
		//		prefix_sum_collider[i].reserve(100 * 100 * 100 / thread_num);
		//		hash_index_count_per_thread_collider[i].reserve(100 * 100 * 100 / thread_num);
		//	}
		//	hash_index_count_per_thread[i].reserve(100 * 100 * 100 / thread_num);
		//}

		//radix_sort = new RadixSort[thread_num];
		//radix_sort_collider = new RadixSort[thread_num];
		//for (unsigned int i = 0; i < thread_num; ++i) {
		//	radix_sort[i].initial(thread, true);
		//	radix_sort_collider[i].initial(thread, true);
		//}
		//testRadixSort();
		//testPrefixSumTime();
	}

	//total_hash_size = max_cell_count * total_triangle_num;
	triangle_number_per_thread.resize(thread_num, 0);
	triangle_number_per_thread_collider.resize(thread_num, 0);
	for (unsigned int i = 0; i < thread_num; ++i) {
		for (int j = 0; j < cloth->size(); ++j) {
			triangle_number_per_thread[i] += cloth->data()[j].mesh_struct.face_index_begin_per_thread[i + 1]
				- cloth->data()[j].mesh_struct.face_index_begin_per_thread[i];
		}
		for (int j = 0; j < tetrahedron->size(); ++j) {
			triangle_number_per_thread[i] += tetrahedron->data()[j].mesh_struct.face_index_begin_per_thread[i + 1]
				- tetrahedron->data()[j].mesh_struct.face_index_begin_per_thread[i];
		}
		for (int j = 0; j < collider->size(); ++j) {
			triangle_number_per_thread_collider[i] += collider->data()[j].mesh_struct.face_index_begin_per_thread[i + 1]
				- collider->data()[j].mesh_struct.face_index_begin_per_thread[i];
		}
	}
	total_hash_size_per_thread.resize(thread_num);
	total_hash_size_per_thread_collider.resize(thread_num);
	for (int i = 0; i < thread_num; ++i) {
		total_hash_size_per_thread[i] = triangle_number_per_thread[i] * max_cell_count;
		total_hash_size_per_thread_collider[i] = triangle_number_per_thread_collider[i] * max_cell_count;
	}
	//spatial_hashing_triangle_value = new unsigned int* [thread_num];
	//spatial_hashing_vertex_value = new unsigned int* [thread_num];
	//spatial_hashing_edge_value = new unsigned int* [thread_num];
	////spatial_hashing_obj_index = new unsigned int* [thread_num];
	//spatial_hashing_triangle_index = new unsigned int* [thread_num];
	//spatial_hashing_vertex_index = new unsigned int* [thread_num];
	//spatial_hashing_edge_index = new unsigned int* [thread_num];

	//spatial_hashing_triangle_collider = new unsigned int* [thread_num];
	//spatial_hashing_vertex_collider = new unsigned int* [thread_num];
	//spatial_hashing_edge_collider = new unsigned int* [thread_num];
	////spatial_hashing_obj_index_collider = new unsigned int* [thread_num];
	//spatial_hashing_triangle_index_collider = new unsigned int* [thread_num];
	//spatial_hashing_vertex_index_collider = new unsigned int* [thread_num];
	//spatial_hashing_edge_index_collider = new unsigned int* [thread_num];

	//for (int i = 0; i < thread_num; ++i) {
	//	spatial_hashing_value[i] = new unsigned int[total_hash_size_per_thread[i]];
	//	//spatial_hashing_obj_index[i] = new unsigned int[total_hash_size_per_thread[i]];
	//	spatial_hashing_triangle_index[i] = new unsigned int[2 * total_hash_size_per_thread[i]];
	//	memset(spatial_hashing_value[i], -1, 4 * total_hash_size_per_thread[i]);
	//	memset(spatial_hashing_triangle_index[i], -1, 8 * total_hash_size_per_thread[i]);
	//	if (has_collider && triangle_number_per_thread_collider[i] != 0) {
	//		spatial_hashing_value_collider[i] = new unsigned int[total_hash_size_per_thread_collider[i]];
	//		//spatial_hashing_obj_index_collider[i] = new unsigned int[total_hash_size_per_thread_collider[i]];
	//		spatial_hashing_triangle_index_collider[i] = new unsigned int[2 * total_hash_size_per_thread_collider[i]];
	//		memset(spatial_hashing_value_collider[i], -1, 4 * total_hash_size_per_thread_collider[i]);
	//	}
	//}

	if (!for_construct_patch) {
		//for (unsigned int i = 0; i < thread_num; ++i) {
		//	radix_sort[i].initialArray(total_hash_size_per_thread[i]);
		//}
		//if (has_collider) {
		//	for (unsigned int i = 0; i < thread_num; ++i) {
		//		if (triangle_number_per_thread_collider[i] != 0) {
		//			radix_sort_collider[i].initialArray(total_hash_size_per_thread_collider[i]);
		//		}
		//	}
		//}

		//obj_is_used0 = new bool** [thread_num];
		//obj_is_used1 = new bool** [thread_num];
		//for (int i = 0; i < thread_num; ++i) {
		//	obj_is_used0[i] = new bool* [collider_begin_obj_index];
		//	obj_is_used1[i] = new bool* [collider_begin_obj_index];
		//	for (int j = 0; j < cloth->size(); ++j) {
		//		obj_is_used0[i][j] = new bool[(*cloth)[j].mesh_struct.triangle_indices.size()];
		//		obj_is_used1[i][j] = new bool[(*cloth)[j].mesh_struct.triangle_indices.size()];
		//		memset(obj_is_used0[i][j], 0, (*cloth)[j].mesh_struct.triangle_indices.size());
		//		memset(obj_is_used1[i][j], 0, (*cloth)[j].mesh_struct.triangle_indices.size());
		//	}
		//	for (int j = 0; j < tetrahedron->size(); ++j) {
		//		obj_is_used0[i][j + tetrahedron_begin_obj_index] = new bool[(*tetrahedron)[j].mesh_struct.triangle_indices.size()];
		//		obj_is_used1[i][j + tetrahedron_begin_obj_index] = new bool[(*tetrahedron)[j].mesh_struct.triangle_indices.size()];
		//		memset(obj_is_used0[i][j + tetrahedron_begin_obj_index], 0, (*tetrahedron)[j].mesh_struct.triangle_indices.size());
		//		memset(obj_is_used1[i][j + tetrahedron_begin_obj_index], 0, (*tetrahedron)[j].mesh_struct.triangle_indices.size());
		//	}
		//	/*if (for_construct_patch) {
		//		for (int j = 0; j < collider->size(); ++j) {
		//			obj_is_used[i][j + collider_begin_obj_index] = new bool[(*collider)[j].mesh_struct.triangle_indices.size()];
		//			memset(obj_is_used[i][j + collider_begin_obj_index], 0, (*collider)[j].mesh_struct.triangle_indices.size());
		//		}
		//	}*/
		//}
		//collider_is_used0 = new bool** [thread_num];
		//if (has_collider) {
		//	for (int i = 0; i < thread_num; ++i) {
		//		collider_is_used0[i] = new bool* [collider->size()];
		//		for (int j = 0; j < collider->size(); ++j) {
		//			collider_is_used0[i][j] = new bool[(*collider)[j].mesh_struct.triangle_indices.size()];
		//			memset(collider_is_used0[i][j], 0, (*collider)[j].mesh_struct.triangle_indices.size());
		//		}
		//	}
		//}
	}
	//
	//actual_hash_value_end_index_ref.resize(thread_num);
	//actual_hash_value_end_index_ref_collider.resize(thread_num);
	//actual_exist_cell_begin_per_thread.resize(thread_num + 1);
	//actual_hash_value_count_per_thread.resize(thread_num, 0);
	//actual_hash_value_count_per_thread_collider.resize(thread_num, 0);



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

	//hash_value_for_test.reserve(total_triangle_num);

	initialHashCell(total_triangle_num, max_index_number_in_one_cell, max_index_number_in_one_cell_collider);

	//vertex_tet_pair.resize(thread_num);
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	vertex_tet_pair[i].reserve(total_triangle_num);
	//}

	//RadixSort radix_sort;
	//radix_sort.initial(thread, true);
	//radix_sort.initialArray(5);
	//unsigned int largerst_count;
	//unsigned int test_[5] = { 3,1,1,7,300 };
	//unsigned int test_2[10] = { 333,1,212,0,11,4,654,7,300,2};
	//radix_sort.radixSort(8748, test_,test_2, largerst_count);
	//for (unsigned int i = 0; i < 5; ++i) {
	//	std::cout << test_[i] << " ";
	//}
	//std::cout << std::endl;
	//for (unsigned int i = 0; i < 10; ++i) {
	//	std::cout << test_2[i] << " ";
	//}
	initialTriangleHash();
	reorganzieDataOfObjects();
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
		for (unsigned int i = 0; i < collider->size(); ++i) {
			collider_tri_aabb[i] = collider->data()[i].triangle_AABB.data();
			collider_vertex_aabb[i] = collider->data()[i].triangle_AABB.data();
			collider_edge_aabb[i] = collider->data()[i].triangle_AABB.data();
			triangle_vertex_index_collider[i] = collider->data()[i].mesh_struct.triangle_indices.data();
			edge_vertex_index_collider[i] = collider->data()[i].mesh_struct.edge_vertices.data();
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
	//obj_triangle_hash.resize(cloth->size() + tetrahedron->size());
	//for (int i = 0; i < cloth->size(); ++i) {
	//	obj_triangle_hash[i].resize((*cloth)[i].mesh_struct.triangle_indices.size());
	//	for (int j = 0; j < obj_triangle_hash[i].size(); ++j) {
	//		obj_triangle_hash[i][j].reserve(32);
	//	}
	//}
	//int obj_no;
	//for (int i = 0; i < tetrahedron->size(); ++i) {
	//	obj_no = i + tetrahedron_begin_obj_index;
	//	obj_triangle_hash[obj_no].resize((*tetrahedron)[i].mesh_struct.triangle_indices.size());
	//	for (int j = 0; j < obj_triangle_hash[obj_no].size(); ++j) {
	//		obj_triangle_hash[obj_no][j].reserve(32);
	//	}
	//}

	//obj_is_used = new bool** [thread_num];
	//collider_is_used = new bool** [thread_num];
	//for (int i = 0; i < thread_num; ++i) {
	//	obj_is_used[i] = new bool* [collider_begin_obj_index];
	//	for (int j = 0; j < cloth->size(); ++j) {
	//		obj_is_used[i][j] = new bool[(*cloth)[j].mesh_struct.triangle_indices.size()];
	//		memset(obj_is_used[i][j], 0, (*cloth)[j].mesh_struct.triangle_indices.size());
	//	}
	//	for (int j = 0; j < tetrahedron->size(); ++j) {
	//		obj_is_used[i][j + tetrahedron_begin_obj_index] = new bool[(*tetrahedron)[j].mesh_struct.triangle_indices.size()];
	//		memset(obj_is_used[i][j + tetrahedron_begin_obj_index], 0, (*tetrahedron)[j].mesh_struct.triangle_indices.size());
	//	}
	//	if (!collider->empty()) {
	//		collider_is_used[i] = new bool* [collider->size()];
	//		for (int j = 0; j < collider->size(); ++j) {
	//			collider_is_used[i][j] = new bool[(*collider)[j].mesh_struct.triangle_indices.size()];
	//			memset(collider_is_used[i][j], 0, (*collider)[j].mesh_struct.triangle_indices.size());
	//		}
	//	}
	//}
}


//void SpatialHashing::initialParallePrefix()
//{
//	unsigned int stage = 32;
//	d_for_prefix_sum = new unsigned int[stage];
//	for (int i = 0; i < stage; ++i) {
//		d_for_prefix_sum[i] = pow(2, i);
//	}
//}




void SpatialHashing::testRadixSort()
{
	unsigned int vec_size = 4e5;
	std::vector<unsigned int> value;
	std::vector<unsigned int> triangle_index;
	std::vector<unsigned int> hash_cloth_No;
	value.reserve(vec_size);
	triangle_index.reserve(vec_size);
	hash_cloth_No.reserve(vec_size);
	time_t t1 = clock();
	for (unsigned int j = 0; j < 100; ++j) {
		int size = 100 * 100 * 100;
		value.clear();
		triangle_index.clear();
		hash_cloth_No.clear();
		for (unsigned int i = 0; i < vec_size; ++i) {
			value.push_back(rand() % size);
			triangle_index.push_back(i);
			hash_cloth_No.push_back(0);
		}
		for (unsigned int i = 0; i < vec_size; ++i) {
			if (i % 2 == 0) {
				value[i] += 100000;
			}
		}
	}
	t1 = clock() - t1;
	//std::cout << "index " << t1 << std::endl;

	//radix_sort.initialArray(vec_size);
	//time_t t = clock();
	//unsigned int count;
	//for (unsigned int j = 0; j < 100; ++j) {
	//	unsigned int size = 100 * 100 * 100;
	//	value.clear();
	//	triangle_index.clear();
	//	hash_cloth_No.clear();
	//	for (unsigned int i = 0; i < vec_size; ++i) {
	//		value.push_back(rand() % size);
	//		triangle_index.push_back(i);
	//		hash_cloth_No.push_back(0);
	//	}
	//	for (unsigned int i = 0; i < vec_size; ++i) {
	//		if (i % 2 == 0) {
	//			value[i] += 100000;
	//		}
	//	}
	//	//for (int i = 0; i < vec_size; ++i) {
	//	//	std::cout << cell[i][1] << " ";
	//	//}
	//	//std::cout << std::endl;
	//	radix_sort.radixSort(size, value.data(), triangle_index.data(), hash_cloth_No.data(), count);
	//	//for (unsigned int i = 0; i < value.size() - 1; ++i) {
	//	//	if (value[i] > value[i + 1]) {
	//	//		std::cout << "error"<<value[i]<<" "<<value[i+1]<< std::endl;
	//	//	}
	//	//}
	//}
	//std::cout << "time " << clock() - t << " " << clock() - t - t1 << std::endl;
	//radix_sort.deleteArray();
	//for (int i = 0; i < vec_size; ++i) {
	//	std::cout << cell[i][1] << " ";
	//}
	//std::cout << std::endl;
	//}

}

//void SpatialHashing::setSpatialHashingPatch(double* scene_aabb)
//{
	//memcpy(this->scene_aabb, scene_aabb,48);
	//for (unsigned int i = 0; i < 3; ++i) {
	//	cell_number[i] = (unsigned int)floor((scene_aabb[i + 3] - scene_aabb[i]) / cell_length) + 1;
	//	hash_max_index[i] = cell_number[i] - 1;
	//}
	//cell_num0_cell_num1 = cell_number[0] * cell_number[1];
	//hash_table_size = cell_number[0] * cell_number[1] * cell_number[2];
	//memset(spatial_hashing_value, -1, 4 * actual_hash_value_count);
//}

void SpatialHashing::setSpatialHashing()
{
	//getSceneAABB();
	for (unsigned int i = 0; i < 3; ++i) {
		cell_number[i] = (unsigned int)floor((scene_aabb[i + 3] - scene_aabb[i]) / cell_length) + 1;
		hash_max_index[i] = cell_number[i] - 1;
	}

	cell_num0_cell_num1 = cell_number[0] * cell_number[1];
	hash_table_size = cell_number[0] * cell_number[1] * cell_number[2];
	std::cout << "spatial hashing size " << hash_table_size << " " << cell_number[0] << " " << cell_number[1] << " " << cell_number[2] << std::endl;

}

void SpatialHashing::buildSpatialHashing(double* scene_aabb)
{
	memcpy(this->scene_aabb, scene_aabb, 48);
	setSpatialHashing();

	time_t t = clock();
	time_t t1 = clock();
	t = clock();
	//for (unsigned int i = 0; i < 100; ++i) {
		//thread->assignTask(this, TRIANGLE_HASHING);
		//setPrifixSum();
	//}

	//t1 = clock();
	//std::cout << " prefix sum with triangle_hashing " <<t1 - t << std::endl;
	//record_time1.push_back(t1 - t);
	//t = clock();
//	for (unsigned int i = 0; i < 100; ++i) {
		//thread->assignTask(this, SH_FIND_ALL_TRIANGLE_PAIRS);	
		//for (unsigned int j = 0; j < thread_num; ++j) {
		//	findAllTrianglePairs(j);
		//}
//	}
	//t1 = clock();
	//std::cout << " find all triangle pairs " << t1 - t << std::endl;
	//record_time2.push_back(t1 - t);
	//std::cout <<"exist cell size "<< prefix_sum[0].size() - 1 << std::endl;
	//findTheMaxNumberOfTriangleInCell();	
	//if (record_time0.size() == 5) {
	//	time_t t_0=0, t_1=0, t_2=0;
	//	for (unsigned int i = 0; i < record_time0.size(); ++i) {
	//		t_0 += record_time0[i];
	//		t_1 += record_time1[i];
	//		t_2 += record_time2[i];
	//	}
	//	std::cout << "ave 5: triangle_hashing " << (double)t_0 / 5000.0 << std::endl;
	//	std::cout << "ave 5: prefix sum " << (double)(t_1-t_0) / 5000.0 << std::endl;
	//	std::cout << "ave 5: find all triangle pairs " << (double)t_0 / 500.0 << std::endl;
	//}
	//findListOfSpatialHashingCell();
	//for (unsigned int j = 0; j < thread_num; ++j) {
	//	findAllTrianglePairsTest(j);
	//}
	//testHashCollisionSP();
	//t = clock();
	//for (unsigned int i = 0; i < 100; ++i) {
	//	clearHashCell();
	//}
	//t1 = clock();
	//std::cout << "clear triangle hashing " << t1 - t << std::endl;

	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	triangleHashingSmallerHashTable(i);
	//}


	t = clock();
	for (unsigned int i = 0; i < 100; ++i) {
		thread->assignTask(this, TRIANGLE_HASHING_SMALLER_HASH_TABLE);
	}
	t1 = clock();
	std::cout << "triangle hashing " << t1 - t << std::endl;
	t = clock();
	for (unsigned int i = 0; i < 100; ++i) {
		for (unsigned int j = 0; j < thread_num; ++j) {
			triangleHashingSmallerHashTable(j);
		}
	}		
	t1 = clock();
	std::cout << "triangle hashing single thread " << t1 - t << std::endl;

	t = clock();
	for (unsigned int i = 0; i < 100; ++i) {
		thread->assignTask(this, TRIANGLE_HASHING_SMALLER_HASH_TABLE);
		recordNonEmptyCell();
	}
	t1 = clock();
	std::cout << "record nonemoty cell " << t1 - t << std::endl;

	t = clock();
	for (unsigned int i = 0; i < 10; ++i) {
		for (unsigned int j = 0; j < thread_num; ++j) {
			findAllPairsHashTable(j);
		}
	}
	t1 = clock();
	std::cout << "find all triangle pairs single thread " << " : " << t1 - t << std::endl;


	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	std::cout << pair_count_per_thread[i] << std::endl;
	//}
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	std::cout << triangle_count_per_thread[i] << std::endl;
	//}

	//	t = clock();
	//for (unsigned int i = 0; i < 100; ++i) {
	//	recordNonEmptyCell();
	//	}
	//	std::cout << "record none empty cell multi thread " << clock() - t << std::endl;

	//t = clock();
	//for (unsigned int i = 0; i < 100; ++i) {
	//	thread->assignTask(this, FIND_ALL_TRIANGLE_PAIRS_HASH_TABLE_COMPARE);
	//}
	//t1 = clock();
	//std::cout << "compare find all triangle pairs multi thread " << t1 - t << std::endl;


	t = clock();
	for (unsigned int i = 0; i < 10; ++i) {
		thread->assignTask(this, FIND_ALL_PAIRS_HASH_TABLE);
	}
	t1 = clock();
	std::cout << "find all triangle pairs multi thread " << t1 - t << std::endl;





	//t = clock();	
	//thread->assignTask(this, TRIANGLE_HASHING_RECORD_REAL_HASH_VALUE);
	//testHashCollision();
	//findTheNumberOfTriangleInCellSP();

	//findMaxTriangleNumInATCell();

	//t = clock();
	//for (unsigned int i = 0; i < 100; ++i) {
	//	testSpatialHashingCellSet0();
	//}
	//t1 = clock();
	//std::cout << "test cell intial for one thread " << t1 - t << std::endl;
	//for (unsigned int i = 0; i < 100; ++i) {
		//for (unsigned int j = 0; j < thread_num; ++j) {
		//	findAllTrianglePairsHashTableTest(j);
		//}
	//}
	//unsigned int triangle_pair_number = 0;
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	std::cout << triangle_pair_number_thread[i] << std::endl;
	//}
	//std::cout << "loop triangle_pair_number " << triangle_pair_number << std::endl;

	//std::vector<unsigned int>* spatial_hashing_cell_ = spatial_hashing_cell[0];
	//t = clock();
	//unsigned int kk = 0;
	//for (unsigned int j = 0; j < 1000; ++j) {
	//	for (unsigned int i = 0; i < hash_cell_count; ++i) {
	//		spatial_hashing_cell_[i].clear();
	//	}
	//	kk++;
	//}
	//std::cout << "clear time " << clock() - t<<" "<<kk << std::endl;

	//t = clock();
	//for (unsigned int kk = 0; kk < 100; ++kk) {
	//	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
	//		thread->assignTask(&tetrahedron->data()[i], TETRAHEDRON_AABB);
	//	}
	//	thread->assignTask(this, TET_HASHING_SMALLER_HASH_TABLE);
	//	thread->assignTask(this, FIND_VERTEX_TET_PAIRS_HASH_TABLE);
	//}
	//std::cout << "time consumed for vertex tet intersection " << clock() - t << " " << std::endl;

	//for (unsigned int kk = 0; kk < 100; ++kk) {
	//	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
	//		for (unsigned int j = 0; j < thread_num; ++j) {
	//			tetrahedron->data()[i].getTetAABBPerThread(j);
	//		}
	//	}
	//	for (unsigned int i = 0; i < thread_num; ++i) {
	//		tetHashingSmallerHashTable(i);
	//	}
	//	for (unsigned int i = 0; i < thread_num; ++i) {
	//		findVertexTetPairHashTable(i);
	//	}
	//}

}


void SpatialHashing::findMaxTriangleNumInATCell()
{
	unsigned int max = 0;
	for (unsigned int i = 0; i < hash_cell_count; ++i) {
		if (max < spatial_hashing_cell_triangle_size[0][i]) {
			max = spatial_hashing_cell_triangle_size[0][i];
		}
	}
	std::cout << "max triangle num " << max << std::endl;
}

//void SpatialHashing::buildSpatialHashingPatchTriangle(double* scene_aabb)
//{
//	//setSpatialHashingPatch(scene_aabb);
//
//}

//void SpatialHashing::findPatch(std::vector<std::vector<std::vector<unsigned int>>>* triangle_patch_)
//{
//	std::vector<std::vector<std::vector<unsigned int>>>triangle_patch;
//	triangle_patch.resize(total_obj_num);
//	std::vector<unsigned int>is_hash_cell_used(hash_table_size,0);
//	for (unsigned int i = 0; i < total_hash_size; ++i) {
//		if (spatial_hashing_value[i] < UINT_MAX) {
//			is_hash_cell_used[spatial_hashing_value[i]]++;
//		}
//	}
//
//	unsigned int actual_exist_hash_cell_number = 0;
//	for (unsigned int i = 0; i < hash_table_size; ++i) {
//		if (is_hash_cell_used[i]!=0) {
//			is_hash_cell_used[i] = actual_exist_hash_cell_number;
//			actual_exist_hash_cell_number++;
//		}
//	}
//	for (unsigned int i = 0; i < total_obj_num; ++i) {
//		triangle_patch[i].resize(actual_exist_hash_cell_number);
//		for (unsigned int j = 0; j < actual_exist_hash_cell_number; ++j) {
//			triangle_patch[i][j].reserve(20);
//		}
//	}
//	unsigned int triangle_index;
//	unsigned int index_of_hash_array;
//
//	unsigned int choosen_hash_value;
//	unsigned int record_size;
//	for (unsigned int i = 0; i < cloth->size(); ++i) {
//		for (unsigned int j = 0; j < cloth->data()[i].mesh_struct.triangle_indices.size(); ++j) {
//			triangle_index = triangle_begin_per_obj[i] + j;
//			index_of_hash_array = max_cell_count * triangle_index;
//			choosen_hash_value = spatial_hashing_value[index_of_hash_array];
//			record_size = triangle_patch[i][is_hash_cell_used[choosen_hash_value]].size();
//			for (unsigned int k = 1; k < max_cell_count; ++k) {
//				index_of_hash_array++;
//				if (spatial_hashing_value[index_of_hash_array] < UINT_MAX) {
//					if (record_size > triangle_patch[i][is_hash_cell_used[spatial_hashing_value[index_of_hash_array]]].size()) {
//						choosen_hash_value = spatial_hashing_value[index_of_hash_array];
//						record_size = triangle_patch[i][is_hash_cell_used[choosen_hash_value]].size();
//						
//					}
//				}				
//			}
//			triangle_patch[i][is_hash_cell_used[choosen_hash_value]].push_back(j);
//		}
//	}
//
//	int obj_No;
//	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
//		obj_No = tetrahedron_begin_obj_index + i;
//		for (unsigned int j = 0; j < tetrahedron->data()[i].mesh_struct.triangle_indices.size(); ++j) {
//			triangle_index = triangle_begin_per_obj[obj_No] + j;
//			index_of_hash_array = max_cell_count * triangle_index;
//			choosen_hash_value = spatial_hashing_value[index_of_hash_array];
//			record_size = triangle_patch[obj_No][is_hash_cell_used[choosen_hash_value]].size();
//			for (unsigned int k = 1; k < max_cell_count; ++k) {
//				index_of_hash_array++;
//				if (spatial_hashing_value[index_of_hash_array] < UINT_MAX) {
//					if (record_size > triangle_patch[obj_No][is_hash_cell_used[spatial_hashing_value[index_of_hash_array]]].size()) {
//						choosen_hash_value = spatial_hashing_value[index_of_hash_array];
//						record_size = triangle_patch[obj_No][is_hash_cell_used[choosen_hash_value]].size();						
//					}
//				}
//			}
//			triangle_patch[obj_No][is_hash_cell_used[choosen_hash_value]].push_back(j);
//		}
//	}
//
//	for (unsigned int i = 0; i < collider->size(); ++i) {
//		obj_No = collider_begin_obj_index + i;
//		for (unsigned int j = 0; j < collider->data()[i].mesh_struct.triangle_indices.size(); ++j) {
//			triangle_index = triangle_begin_per_obj[obj_No] + j;
//			index_of_hash_array = max_cell_count * triangle_index;
//			choosen_hash_value = spatial_hashing_value[index_of_hash_array];
//			record_size = triangle_patch[obj_No][is_hash_cell_used[choosen_hash_value]].size();
//			for (unsigned int k = 1; k < max_cell_count; ++k) {
//				index_of_hash_array++;
//				if (spatial_hashing_value[index_of_hash_array] < UINT_MAX) {
//					if (record_size > triangle_patch[obj_No][is_hash_cell_used[spatial_hashing_value[index_of_hash_array]]].size()) {
//						choosen_hash_value = spatial_hashing_value[index_of_hash_array];
//						record_size = triangle_patch[obj_No][is_hash_cell_used[choosen_hash_value]].size();
//					}
//				}
//			}
//			triangle_patch[obj_No][is_hash_cell_used[choosen_hash_value]].push_back(j);
//		}
//	}
//
//	triangle_patch_->resize(total_obj_num);
//	for (unsigned int i = 0; i < total_obj_num; ++i) {
//		for (unsigned int j = 0; j < triangle_patch[i].size(); ++j) {
//			if (!triangle_patch[i][j].empty()) {
//				triangle_patch_->data()[i].push_back(triangle_patch[i][j]);
//			}
//		}
//	}
//
//	deleteArray();
//}



//void SpatialHashing::prefixSumParallel()
//{
//	unsigned int stage = log2(hash_table_size);
//	//std::cout << "stage " << stage << std::endl;
//	for (unsigned int i = 0; i < stage; ++i) {
//		thread->assignTask(this, PREFIX_SUM_UP, i);
//	}
//	for (int i = stage - 1; i >= 0; --i) {
//		thread->assignTask(this, PREFIX_SUM_DOWN, i);
//	}
//}

////PREFIX_SUM_UP
//void SpatialHashing::prefixSumParallelUp(int thread_No, unsigned int stage)
//{
//
//	unsigned int d_1 = d_for_prefix_sum[stage + 1];
//	unsigned int d = d_for_prefix_sum[stage];
//
//	unsigned int size = (hash_table_size / d_1) / thread_num;
//
//	if ((hash_table_size / d_1) % thread_num != 0) {
//		size++;
//	}
//	unsigned int start = d_1 * thread_No * size + d_1 - 1;
//	unsigned int end = d_1 * (thread_No+1) * size + d_1 - 1;
//
//	if (end > hash_table_size) {//
//		end = hash_table_size;
//	}
//	//std::cout << start << " " << end << " " << d_1<<" "<< hash_table_size << std::endl;
//	for (unsigned int i = start; i < end; i += d_1) {
//		//if (d_1 == 1024) {
//		//	std::cout << i + d_1 - 1 << std::endl;
//		//}
//		*(prefix_sum_1_address + i ) += *(prefix_sum_1_address + i - d);
//	}	
//}

////PREFIX_SUM_DOWN
//void SpatialHashing::prefixSumParallelDown(int thread_No, unsigned int stage)
//{
//	unsigned int d_1 = d_for_prefix_sum[stage + 1];
//	unsigned int d = d_for_prefix_sum[stage];
//	unsigned int size = (hash_table_size / d_1) / thread_num;
//	if ((hash_table_size / d_1) % thread_num != 0) {
//		size++;
//	}
//	unsigned int start = d_1 * thread_No * size + d_1 - 1;
//	unsigned int end = d_1 * (thread_No+1) * size + d_1 - 1;
//
//	if (end + d > hash_table_size) {
//		end = hash_table_size - d;
//	}
//
//	for (unsigned int i = start; i < end; i += d_1) {
//		*(prefix_sum_1_address + i + d) += *(prefix_sum_1_address + i);
//	}
//	
//}

//void SpatialHashing::testPrefixSumTime()
//{
	//delete[] spatial_hashing_value;
	//unsigned int hash_size_dimension = 120;
	//hash_table_size = hash_size_dimension * hash_size_dimension * hash_size_dimension;
	//unsigned int triangle_number = 200000;
	//unsigned int max_per_triangle = 27;
	//unsigned int max_count_hash_triangle = max_per_triangle * triangle_number;
	//radix_sort.initialArray(max_count_hash_triangle);
	//spatial_hashing_value = new unsigned int[max_count_hash_triangle];
	//memset(spatial_hashing_value, -1, 4 * max_count_hash_triangle);
	//unsigned int size[3];
	//unsigned int start[3];

	//std::vector<unsigned int> triangle_index(max_count_hash_triangle, 0);
	//std::vector<unsigned int> hash_cloth_No(max_count_hash_triangle, 0);
	//for (unsigned int i = 0; i < triangle_number; ++i) {
	//	size[0] = 1 + rand() % 3;
	//	size[1] = 1 + rand() % 3;
	//	size[2] = 1 + rand() % 3;
	//	start[0] = 3 + rand() % (hash_size_dimension - 20);
	//	start[1] = 3 + rand() % (hash_size_dimension - 20);
	//	start[2] = 3 + rand() % (hash_size_dimension - 20);
	//	unsigned int note = 0;
	//	for (unsigned int j = start[0]; j < start[0] + size[0]; ++j) {
	//		for (unsigned int k = start[1]; k < start[1] + size[1]; ++k) {
	//			for (unsigned int l = start[2]; l < start[2] + size[2]; ++l) {
	//				spatial_hashing_value[max_per_triangle * i + note] = j + hash_size_dimension * k + hash_size_dimension * hash_size_dimension * l;
	//				note++;
	//			}
	//		}
	//	}
	//}
	//radix_sort.radixSort(hash_table_size, spatial_hashing_value, triangle_index.data(), hash_cloth_No.data(),
	//	largest_count_in_hash_value_list);
	//radix_sort.deleteArray();

	//actual_hash_value_count = max_count_hash_triangle - largest_count_in_hash_value_list;
	//for (unsigned int i = actual_hash_value_count; i < max_count_hash_triangle; ++i) {
	//	if (spatial_hashing_value[i] == UINT_MAX) {
	//		actual_hash_value_count = i;
	//		break;
	//	}
	//}
	//std::cout << "actual_hash_value_count " << actual_hash_value_count << std::endl;

	//prefix_sum.resize(hash_table_size + 1);
	//prefix_sum_1_address = prefix_sum.data() + 1;
	//arrangeIndex(thread_num, hash_table_size, total_hash_count_start_per_thread);

	////

	//time_t t1 = clock();
	//
	//for (int j = 0; j < 1000; ++j) {		
	//	arrangeIndex(thread_num, actual_hash_value_count, actual_hash_count_start_per_thread);
	//	thread->assignTask(this, MEMSET_PREFIX);
	//	//thread->assignTask(this, PREPARE_FOR_ACTUAL_HASH_VALUE_COUNT_THREAD);		
	//	prepareForActualHashValueCount();
	//	thread->assignTask(this, ADD_COUNT_FOR_prefix_sum);
	//}
	//std::cout << "time parallel sum count " << clock() - t1 << std::endl;
	//t1 = clock();
	//for (int j = 0; j < 1000; ++j) {
	//	arrangeIndex(thread_num, actual_hash_value_count, actual_hash_count_start_per_thread);
	//	thread->assignTask(this, MEMSET_PREFIX);
	//	//thread->assignTask(this, PREPARE_FOR_ACTUAL_HASH_VALUE_COUNT_THREAD);
	//	prepareForActualHashValueCount();
	//	thread->assignTask(this, ADD_COUNT_FOR_prefix_sum);
	//	thread->assignTask(this, PREFIX_SUM_THREAD_1);
	//	for (int i = 1; i < thread_num; ++i) {
	//		prefix_sum_thread_start[i] = *(prefix_sum_1_address + total_hash_count_start_per_thread[i] - 1) + prefix_sum_thread_start[i - 1];
	//	}
	//	thread->assignTask(this, PREFIX_SUM_THREAD_2);
	//	//prefixSumParallel();
	//}
	//std::cout << "time parallel sum count " << clock() - t1 << std::endl;
	//testPrifixSum1();

	//delete[] spatial_hashing_value;

//}


// SH_FIND_ALL_TRIANGLE_PAIRS
void SpatialHashing::findAllTrianglePairs(int thread_No)
{
	//unsigned int last_cell_index = actual_exist_cell_begin_per_thread[thread_No + 1];
	//std::vector<unsigned int> cell_triangle_index;
	////std::vector<unsigned int> cell_obj_index;
	//std::vector<unsigned int> cell_collider_triangle_index;
	////std::vector<unsigned int> cell_collider_index;
	//cell_triangle_index.reserve(800);
	////cell_obj_index.reserve(30);
	//cell_collider_triangle_index.reserve(800);
	////cell_collider_index.reserve(30);

	//unsigned int previous_tri_num;
	//unsigned int current_thread_triangle_num;

	//unsigned int* triangle_pair_;
	//unsigned int* triangle_pair_with_collider_;

	//triangle_pair_ = triangle_pair[thread_No];
	//triangle_pair_with_collider_ = triangle_pair_with_collider[thread_No];

	//double* aabb_1;

	//unsigned int obj_index_0, obj_index_1, triangle_index_0, triangle_index_1;
	////triangle_pair_number_thread[thread_No] = 0;
	//for (unsigned int cell_index = actual_exist_cell_begin_per_thread[thread_No];
	//	cell_index < last_cell_index; ++cell_index) {
	//	cell_triangle_index.clear();
	//	//cell_obj_index.clear();
	//	cell_collider_triangle_index.clear();
	//	//cell_collider_index.clear();
	//	for (unsigned int t = 0; t < thread_num; ++t) {
	//		current_thread_triangle_num = prefix_sum[t][cell_index + 1] - prefix_sum[t][cell_index];
	//		if (current_thread_triangle_num > 0) {
	//			previous_tri_num = cell_triangle_index.size();
	//			cell_triangle_index.resize(previous_tri_num + current_thread_triangle_num);
	//			//cell_obj_index.resize(previous_tri_num + current_thread_triangle_num);
	//			memcpy(cell_triangle_index.data() + previous_tri_num, spatial_hashing_triangle_index[t] + prefix_sum[t][cell_index], 4 * current_thread_triangle_num);
	//			//memcpy(cell_obj_index.data() + previous_tri_num, spatial_hashing_obj_index[t] + prefix_sum[t][cell_index], 4 * current_thread_triangle_num);
	//		}
	//	}
	//	if (has_collider) {
	//		for (unsigned int t = 0; t < thread_num; ++t) {
	//			current_thread_triangle_num = prefix_sum_collider[t][(cell_index << 1) + 1] - prefix_sum_collider[t][cell_index << 1];
	//			if (current_thread_triangle_num > 0) {
	//				previous_tri_num = cell_collider_triangle_index.size();
	//				cell_collider_triangle_index.resize(previous_tri_num + current_thread_triangle_num);
	//				//cell_collider_index.resize(previous_tri_num + current_thread_triangle_num);
	//				memcpy(cell_collider_triangle_index.data() + previous_tri_num,
	//					spatial_hashing_triangle_index_collider[t] + prefix_sum_collider[t][cell_index << 1], 4 * current_thread_triangle_num);
	//				//memcpy(cell_collider_index.data() + previous_tri_num,
	//				//	spatial_hashing_obj_index_collider[t] + prefix_sum_collider[t][cell_index << 1], 4 * current_thread_triangle_num);
	//			}
	//		}
	//	}
	//	//there exist collision

	//	if (cell_triangle_index.size() > 2 || (!cell_collider_triangle_index.empty())) {
	//		//if (thread_No == 4) {
	//		//	std::cout << cell_triangle_index.size() << std::endl;
	//		//	for (unsigned int k = 0; k < cell_triangle_index.size(); ++k) {
	//		//		std::cout << cell_triangle_index[k] << " ";
	//		//	}
	//		//	std::cout << std::endl;
	//		//}

	//		for (unsigned int i = 0; i < cell_triangle_index.size(); i+=2) {
	//			obj_index_0 = cell_triangle_index[i + 1];  triangle_index_0 = cell_triangle_index[i];
	//			aabb_1 = obj_tri_aabb[obj_index_0][triangle_index_0].data();
	//			for (unsigned int j = i + 2; j < cell_triangle_index.size(); j+=2) {
	//				obj_index_1 = cell_triangle_index[j+1];
	//				triangle_index_1 = cell_triangle_index[j];
	//				//triangle_pair_number_thread[thread_No]++;
	//				if (AABB::AABB_intersection(aabb_1, obj_tri_aabb[obj_index_1][triangle_index_1].data())) {

	//					triangle_pair_->emplace_back(triangle_index_0);
	//					triangle_pair_->emplace_back(obj_index_0);
	//					triangle_pair_->emplace_back(triangle_index_1);
	//					triangle_pair_->emplace_back(obj_index_1);

	//				}
	//			}
	//			for (unsigned int j = 0; j < cell_collider_triangle_index.size(); j+=2) {
	//				obj_index_1 = cell_collider_triangle_index[j + 1];
	//				triangle_index_1 = cell_collider_triangle_index[j];
	//				if (AABB::AABB_intersection(aabb_1, collider_tri_aabb[obj_index_1][triangle_index_1].data())) {
	//					triangle_pair_with_collider_->emplace_back(triangle_index_0);
	//					triangle_pair_with_collider_->emplace_back(obj_index_0);
	//					triangle_pair_with_collider_->emplace_back(triangle_index_1);
	//					triangle_pair_with_collider_->emplace_back(obj_index_1);
	//				}
	//			}
	//		}
	//	}
	//}

}

void SpatialHashing::findAllTrianglePairsTest(int thread_No)
{
	//unsigned int last_cell_index = actual_exist_cell_begin_per_thread[thread_No + 1];

	//std::vector<unsigned int> cell_triangle_index;
	////std::vector<unsigned int> cell_obj_index;
	//std::vector<unsigned int> cell_collider_triangle_index;
	////std::vector<unsigned int> cell_collider_index;

	//cell_triangle_index.reserve(60);
	////cell_obj_index.reserve(30);
	//cell_collider_triangle_index.reserve(60);
	////cell_collider_index.reserve(30);

	//unsigned int previous_tri_num;
	//unsigned int current_thread_triangle_num;

	//std::vector<unsigned int>* triangle_pair_;
	//std::vector<unsigned int>* triangle_pair_with_collider_;

	//triangle_pair_ = &triangle_pair[thread_No];
	//triangle_pair_with_collider_ = &triangle_pair_with_collider[thread_No];

	//triangle_pair_->clear();
	//triangle_pair_with_collider_->clear();

	//double* aabb_1;

	//unsigned int obj_index_0, obj_index_1, triangle_index_0, triangle_index_1;

	//triangle_pair_number_thread[thread_No] = 0;

	//for (unsigned int cell_index = actual_exist_cell_begin_per_thread[thread_No];
	//	cell_index < last_cell_index; ++cell_index) {
	//	cell_triangle_index.clear();
	//	//cell_obj_index.clear();
	//	cell_collider_triangle_index.clear();
	//	//cell_collider_index.clear();

	//	for (unsigned int t = 0; t < thread_num; ++t) {
	//		current_thread_triangle_num = prefix_sum[t][cell_index + 1] - prefix_sum[t][cell_index];
	//		if (current_thread_triangle_num > 0) {
	//			previous_tri_num = cell_triangle_index.size();
	//			cell_triangle_index.resize(previous_tri_num + current_thread_triangle_num);
	//			//cell_obj_index.resize(previous_tri_num + current_thread_triangle_num);
	//			memcpy(cell_triangle_index.data() + previous_tri_num, spatial_hashing_triangle_index[t] + prefix_sum[t][cell_index], 4 * current_thread_triangle_num);
	//			//memcpy(cell_obj_index.data() + previous_tri_num, spatial_hashing_obj_index[t] + prefix_sum[t][cell_index], 4 * current_thread_triangle_num);
	//		}
	//	}
	//	if (has_collider) {
	//		for (unsigned int t = 0; t < thread_num; ++t) {
	//			current_thread_triangle_num = prefix_sum_collider[t][(cell_index << 1) + 1] - prefix_sum_collider[t][cell_index << 1];
	//			if (current_thread_triangle_num > 0) {
	//				previous_tri_num = cell_collider_triangle_index.size();
	//				cell_collider_triangle_index.resize(previous_tri_num + current_thread_triangle_num);
	//				//cell_collider_index.resize(previous_tri_num + current_thread_triangle_num);
	//				memcpy(cell_collider_triangle_index.data() + previous_tri_num,
	//					spatial_hashing_triangle_index_collider[t] + prefix_sum_collider[t][cell_index << 1], 4 * current_thread_triangle_num);
	//				//memcpy(cell_collider_index.data() + previous_tri_num,
	//				//	spatial_hashing_obj_index_collider[t] + prefix_sum_collider[t][cell_index << 1], 4 * current_thread_triangle_num);
	//			}
	//		}
	//	}
	//	//there exist collision

	//	if (cell_triangle_index.size() > 2 || (!cell_collider_triangle_index.empty())) {
	//		//if (thread_No == 4) {
	//		//	std::cout << cell_triangle_index.size() << std::endl;
	//		//	for (unsigned int k = 0; k < cell_triangle_index.size(); ++k) {
	//		//		std::cout << cell_triangle_index[k] << " ";
	//		//	}
	//		//	std::cout << std::endl;
	//		//}

	//		for (unsigned int i = 0; i < cell_triangle_index.size(); i += 2) {
	//			obj_index_0 = cell_triangle_index[i + 1];  triangle_index_0 = cell_triangle_index[i];
	//			aabb_1 = obj_tri_aabb[obj_index_0][triangle_index_0].data();
	//			for (unsigned int j = i + 2; j < cell_triangle_index.size(); j += 2) {
	//				obj_index_1 = cell_triangle_index[j + 1];
	//				triangle_index_1 = cell_triangle_index[j];
	//				triangle_pair_number_thread[thread_No]++;
	//				if (AABB::AABB_intersection(aabb_1, obj_tri_aabb[obj_index_1][triangle_index_1].data())) {

	//					triangle_pair_->emplace_back(triangle_index_0);
	//					triangle_pair_->emplace_back(obj_index_0);
	//					triangle_pair_->emplace_back(triangle_index_1);
	//					triangle_pair_->emplace_back(obj_index_1);

	//				}
	//			}
	//			for (unsigned int j = 0; j < cell_collider_triangle_index.size(); j += 2) {
	//				obj_index_1 = cell_collider_triangle_index[j + 1];
	//				triangle_index_1 = cell_collider_triangle_index[j];
	//				if (AABB::AABB_intersection(aabb_1, collider_tri_aabb[obj_index_1][triangle_index_1].data())) {
	//					triangle_pair_with_collider_->emplace_back(triangle_index_0);
	//					triangle_pair_with_collider_->emplace_back(obj_index_0);
	//					triangle_pair_with_collider_->emplace_back(triangle_index_1);
	//					triangle_pair_with_collider_->emplace_back(obj_index_1);
	//				}
	//			}
	//		}
	//	}
	//}

}


// FIND_VERTEX_TET_PAIRS_HASH_TABLE
void SpatialHashing::findVertexTetPairHashTable(int thread_No)
{
	//unsigned int vertex_end;
	//int* vertex_index_on_surface;
	//std::array<double, 3>* pos;
	//unsigned int obj_index;
	//vertex_tet_pair[thread_No].clear();

	//for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
	//	obj_index = tetrahedron_begin_obj_index + i;
	//	vertex_index_on_surface = tetrahedron->data()[i].mesh_struct.vertex_index_on_sureface.data();
	//	vertex_end = tetrahedron->data()[i].mesh_struct.vertex_index_on_surface_begin_per_thread[thread_No + 1];
	//	pos = tetrahedron->data()[i].mesh_struct.vertex_position.data();
	//	for (unsigned int j = tetrahedron->data()[i].mesh_struct.vertex_index_on_surface_begin_per_thread[thread_No];
	//		j < vertex_end; ++j) {
	//		vertexTetIntersect(thread_No, pos[vertex_index_on_surface[j]].data(), vertex_index_on_surface[j], obj_index);
	//	}
	//}
}





void SpatialHashing::vertexTetIntersect(int thread_No, double* pos, unsigned int vertex_index, unsigned int obj_index)
{
	//std::vector<unsigned int>* vertex_tet_pair_ = &vertex_tet_pair[thread_No];
	//unsigned int  spatial_index[3];
	//for (unsigned int j = 0; j < 3; ++j) {
	//	spatial_index[j] = (unsigned int)floor((pos[j] - scene_aabb[j]) / cell_length);
	//}
	//unsigned int hash_index = ((spatial_index[0] * P1) ^ (P2 * spatial_index[1]) ^ (spatial_index[2] * P3)) % hash_cell_count;
	//unsigned int compare_tet_obj_index; unsigned int compare_tet_index;
	//int* indices;
	//std::array<double, 3>* tet_pos;
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	for (unsigned int j = 0; j < spatial_hashing_cell[i][hash_index].size(); j+=2) {
	//		compare_tet_obj_index = spatial_hashing_cell[i][hash_index][j + 1];
	//		compare_tet_index= spatial_hashing_cell[i][hash_index][j];
	//		indices = tetrahedron->data()[compare_tet_obj_index].mesh_struct.indices[compare_tet_index].data();
	//		tet_pos = tetrahedron->data()[compare_tet_obj_index].mesh_struct.vertex_position.data();
	//		if (compare_tet_obj_index != obj_index) {
	//			if (testIntersect(pos, tet_pos[indices[0]].data(), tet_pos[indices[1]].data(), tet_pos[indices[2]].data(),
	//				tet_pos[indices[3]].data())) {
	//				vertex_tet_pair_->emplace_back(vertex_index);
	//				vertex_tet_pair_->emplace_back(obj_index);
	//				vertex_tet_pair_->emplace_back(compare_tet_index);
	//				vertex_tet_pair_->emplace_back(compare_tet_obj_index);
	//			}
	//		}
	//		else {
	//			if (!vertexFromTet(vertex_index, indices))
	//			{
	//				if (testIntersect(pos, tet_pos[indices[0]].data(), tet_pos[indices[1]].data(), tet_pos[indices[2]].data(),
	//					tet_pos[indices[3]].data())) {
	//					vertex_tet_pair_->emplace_back(vertex_index);
	//					vertex_tet_pair_->emplace_back(obj_index);
	//					vertex_tet_pair_->emplace_back(compare_tet_index);
	//					vertex_tet_pair_->emplace_back(compare_tet_obj_index);
	//				}
	//			}
	//		}
	//	}
	//}
}


bool SpatialHashing::testIntersect(double* pos, double* p0, double* p1, double* p2, double* p3)
{
	double A_inverse[9];
	double A[9];
	SUB(A, p1, p0);
	SUB((A + 3), p2, p0);
	SUB((A + 6), p3, p0);
	inverse3X3(A, A_inverse);

	double p_x0[3];
	SUB(p_x0, pos, p0);
	double beta[3];
	for (unsigned int i = 0; i < 3; ++i) {
		beta[i] = A_inverse[i] * p_x0[0] + A_inverse[i + 3] * p_x0[1] + A_inverse[i + 6] * p_x0[2];
	}
	if (beta[0] >= 0 && beta[1] >= 0 && beta[2] >= 0 && (1.0 - beta[0] - beta[1] - beta[2]) >= 0) {
		return true;
	}
	return false;
}

bool SpatialHashing::vertexFromTet(int vertex_index, int* tet_indices)
{
	for (unsigned int i = 0; i < 4; ++i) {
		if (vertex_index == tet_indices[i]) {
			return true;
		}
	}
	return false;
}

// FIND_ALL_TRIANGLE_PAIRS_HASH_TABLE
//void SpatialHashing::findAllTrianglePairsHashTable(int thread_No)
//{
//	for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
//		for (unsigned int j = obj_triangle_index_begin_per_thread[i][thread_No]; j < obj_triangle_index_begin_per_thread[i][thread_No + 1]; ++j) {
//			searchTriangle(obj_tri_aabb[i][j].data(), obj_triangle_hash[i][j].data(), obj_triangle_hash[i][j].size(),
//				i, j, &triangle_pair[thread_No], &triangle_pair_with_collider[thread_No], thread_No);
//		}
//	}
//}

void SpatialHashing::findTheNumberOfTriangleInCellSP()
{
	unsigned int num1 = 0;
	unsigned int max_num = 0;
	unsigned int num = 0;
	unsigned int non_empty_count = 0;

	//hash_cell_triangle_count.clear();
	//for (unsigned int i = 0; i < hash_cell_count; ++i) {
	//	num1 = 0;
	//	for (unsigned int j = 0; j < thread_num; ++j) {
	//		num1 += spatial_hashing_cell[j][i].size();
	//	}
	//	if (num1 > max_num) {
	//		max_num = num1;
	//	}
	//	num += num1;
	//	if (num1 > 0) {
	//		hash_cell_triangle_count.emplace_back(num1);
	//		non_empty_count++;
	//	}
	//}
	//std::cout << "max_num " << max_num << " " << (double)num / (double)non_empty_count << std::endl;
	//std::string name = "cell num";
	//WriteTxt::writeTxt(hash_cell_triangle_count, name);
	//std::cout << "non empty cell " << non_empty_cell_index_begin_per_thread[8] << std::endl;

	//hash_cell_triangle_count.clear();
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	num = 0;
	//	for (unsigned int j = non_empty_cell_index_begin_per_thread[i]; j < non_empty_cell_index_begin_per_thread[i + 1]; ++j) {
	//		num1 = 0;
	//		for (unsigned int k = 0; k < thread_num; ++k) {
	//			num1 += (spatial_hashing_cell[k][non_empty_cell_index[0][j]].size()>>1);
	//		}
	//		//num += num1;
	//		hash_cell_triangle_count.emplace_back(num1);
	//		num += (num1 * (num1 - 1)) >> 1;
	//	}	
	//}
	//std::string name = "cell num";
	//WriteTxt::writeTxt(hash_cell_triangle_count, name);
	//std::cout << "=====" << std::endl;
	////for (unsigned int i = 0; i < 8; ++i) {
	////	std::cout << indicator[i][0] << std::endl;
	////}
	//std::cout << "=====" << std::endl;
}






bool SpatialHashing::testA(unsigned int& a) {
	if (a % 2 == 1) {
		return true;
	}
	return false;
}


//void SpatialHashing::findAllTrianglePairsHashTable()
//{
//	for (unsigned int i = 0; i < global_cell_start[0][0]; ++i) {
//		thread->assignTask(this, FIND_ALL_TRIANGLE_PAIRS_HASH_TABLE, i);
//	}	
//}


//void SpatialHashing::findAllTrianglePairsHashTableSingleThread()
//{
//	for (unsigned int i = 0; i < global_cell_start[0][0]; ++i) {
//		for (unsigned int j = 0; j < thread_num; ++j) {
//			findAllTrianglePairsHashTable(j, i);
//		}	
//	}
//}

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
	unsigned int* edge_edge_collider_pair_;

	primitive_pair_ = edge_edge_pair[thread_No] + 1;
	edge_edge_collider_pair_ = edge_edge_pair_collider[thread_No] + 1;

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
		cell_edge_index = spatial_hashing_cell_edge[0][cell_index];
		if (hash_cell_element_size > 2) {
			for (unsigned int i = 0; i < hash_cell_element_size; i += 2) {
				obj_index_0 = cell_edge_index[i + 1]; primitive_index_0 = cell_edge_index[i];
				aabb_1 = obj_edge_aabb[obj_index_0][primitive_index_0].data();
				for (unsigned int j = i + 2; j < hash_cell_element_size; j += 2) {
					if (AABB::AABB_intersection(aabb_1, obj_edge_aabb[cell_edge_index[j+1]][cell_edge_index[j]].data())) {
						if (obj_index_0 == cell_edge_index[j + 1]) {
							if (!edgeEdgeconnected(edge_vertex_index[obj_index_0] + (primitive_index_0 << 1),
								edge_vertex_index[obj_index_0] + (cell_edge_index[j] << 1))) {

								//if (primitive_index_0 < cell_edge_index[j]) {
								//	if (primitive_index_0 == 12106) {
								//		std::cout << cell_edge_index[j] << std::endl;
								//	}
								//}
								//else {
								//	if (cell_edge_index[j] == 12106) {
								//		std::cout << primitive_index_0 << std::endl;
								//	}
								//}

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
		if (has_collider_) {
			hash_cell_element_size_collider = spatial_hashing_cell_collider_edge_size[0][cell_index];
			cell_collider_edge_index = spatial_hashing_cell_collider_edge[0][cell_index];
			if (hash_cell_element_size > 0) {
				for (unsigned int i = 0; i < hash_cell_element_size; i += 2) {
					obj_index_0 = cell_edge_index[i + 1]; primitive_index_0 = cell_edge_index[i];
					aabb_1 = obj_edge_aabb[obj_index_0][primitive_index_0].data();
					for (unsigned int j = 0; j < hash_cell_element_size_collider; j += 2) {
						if (AABB::AABB_intersection(aabb_1, 
							collider_edge_aabb[cell_collider_edge_index[j+1]][cell_collider_edge_index[j]].data())) {
							memcpy(edge_edge_collider_pair_, cell_edge_index + i, 8);
							edge_edge_collider_pair_ += 2;
							memcpy(edge_edge_collider_pair_, cell_collider_edge_index + j, 8);
							edge_edge_collider_pair_ += 2;

						}
					}
				}
			}
		}
	}

	edge_edge_pair[thread_No][0] = primitive_pair_ - edge_edge_pair[thread_No] - 1;
	edge_edge_pair_collider[thread_No][0] = edge_edge_collider_pair_ - edge_edge_pair_collider[thread_No] - 1;
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
	unsigned int* vertex_collider_triangle_pair_;
	unsigned int* vertex_triangle_collider_pair_;

	primitive_pair_ = vertex_triangle_pair[thread_No] + 1;
	vertex_collider_triangle_pair_ = vertex_collider_triangle_obj_pair[thread_No] + 1;
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

	for (unsigned int cell_index = start; cell_index < last_cell_index; ++cell_index) {
		//cell_index = non_empty_cell_index_[kk];
		hash_cell_element_size = spatial_hashing_cell_triangle_size[0][cell_index];
		cell_triangle_index = spatial_hashing_cell_triangle[0][cell_index];

		hash_cell_vertex_size = spatial_hashing_cell_vertex_size[0][cell_index];
		cell_vertex_index = spatial_hashing_cell_vertex[0][cell_index];

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
			hash_cell_element_size_collider = spatial_hashing_cell_collider_triangle_size[0][cell_index];
			cell_collider_triangle_index = spatial_hashing_cell_collider_triangle[0][cell_index];

			hash_cell_vertex_size_collider = spatial_hashing_cell_collider_vertex_size[0][cell_index];
			cell_collider_vertex_index = spatial_hashing_cell_collider_vertex[0][cell_index];

			if (hash_cell_vertex_size_collider) {
				for (unsigned int j = 0; j < hash_cell_vertex_size_collider; j += 2) {
					obj_index_0 = cell_collider_vertex_index[j + 1]; primitive_index_0 = cell_collider_vertex_index[j];
					aabb_1 = collider_vertex_aabb[obj_index_0][primitive_index_0].data();
					for (unsigned int i = 0; i < hash_cell_element_size; i += 2) {
						if (AABB::AABB_intersection(aabb_1, obj_tri_aabb[cell_triangle_index[i + 1]][cell_triangle_index[i]].data())) {
							memcpy(vertex_collider_triangle_pair_, cell_collider_vertex_index + j, 8);
							vertex_collider_triangle_pair_ += 2;
							memcpy(vertex_collider_triangle_pair_, cell_triangle_index + i, 8);
							vertex_collider_triangle_pair_ += 2;
						}
					}
				}
			}

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
	vertex_collider_triangle_obj_pair[thread_No][0] = vertex_collider_triangle_pair_ - vertex_collider_triangle_obj_pair[thread_No] - 1;

}


//FIND_ALL_PAIRS_HASH_TABLE
void SpatialHashing::findAllPairsHashTable(int thread_No)
{
	findAllVertexTrianglePairs(thread_No);
	findAllEdgeEdgePairs(thread_No);
}



//void SpatialHashing::findAllTrianglePairsHashTableSingleThread()
//{
//	//triangle_pair_number_thread[thread_No] = 0;
//	unsigned int last_cell_index = 35724;// global_cell_start[thread_num];
//	//unsigned int last_cell_index = non_empty_cell_index_begin_per_thread[thread_No + 1];
//	unsigned int start = 0;// global_cell_start[0];
//	//unsigned int start= non_empty_cell_index_begin_per_thread[thread_No];
//	
//	unsigned int* cell_triangle_index;
//	unsigned int* cell_collider_triangle_index;
//	unsigned int* non_empty_cell_index_ = non_empty_cell_index[0] + 1;
//	unsigned int* triangle_pair_;
//	unsigned int* triangle_pair_with_collider_;
//
//	triangle_pair_ = triangle_pair[0] + 1;
//	triangle_pair_with_collider_ = triangle_pair_with_collider[0] + 1;
//	double* aabb_1;
//
//	unsigned int obj_index_0, obj_index_1, triangle_index_0, triangle_index_1;
//	//triangle_pair_number_thread[thread_No] = 0;
//	unsigned int cell_index;
//
//	unsigned int* cell_triangle_index_;
//	unsigned int cell_vector_size;
//	unsigned int cell_vector_size_collider;
//
//
//	unsigned int* cell_triangle_index_address;
//
//	unsigned int hash_cell_element_size;
//	unsigned int hash_cell_element_size_collider;
//
//	unsigned int a = 0;
//
//
//	for (unsigned int cell_index = start; cell_index < last_cell_index; ++cell_index) {
//		//cell_index = non_empty_cell_index_[kk];
//		//cell_triangle_index.clear();
//		//cell_collider_triangle_index.clear();
//		//cell_vector_size = 0;
//		//for (unsigned int t = 0; t < thread_num; ++t) {
//		//	if (!spatial_hashing_cell[t][cell_index].empty()) {
//		//		cell_vector_size += spatial_hashing_cell[t][cell_index].size();
//		//		//cell_triangle_index.insert(cell_triangle_index.end(), spatial_hashing_cell[t][cell_index].begin(),
//		//		//	spatial_hashing_cell[t][cell_index].end());
//		//	}
//		//}
//		//cell_triangle_index.resize(cell_vector_size);
//		//cell_triangle_index_address = cell_triangle_index.data();
//		//for (unsigned int t = 0; t < thread_num; ++t) {
//		//	if (!spatial_hashing_cell[t][cell_index].empty()) {
//		//		memcpy(cell_triangle_index_address, spatial_hashing_cell[t][cell_index].data(), 
//		//			4 * spatial_hashing_cell[t][cell_index].size());
//		//		cell_triangle_index_address += spatial_hashing_cell[t][cell_index].size();
//		//	}
//		//}
//		//if (has_collider) {
//		//	cell_vector_size_collider = 0;
//		//	for (unsigned int t = 0; t < thread_num; ++t) {
//		//		if (!spatial_hashing_cell_collider[t][cell_index].empty()) {
//		//			cell_vector_size_collider += spatial_hashing_cell_collider[t][cell_index].size();
//		//	/*		cell_collider_triangle_index.insert(cell_collider_triangle_index.end(), spatial_hashing_cell_collider[t][cell_index].begin(),
//		//				spatial_hashing_cell_collider[t][cell_index].end());*/
//		//		}
//		//	}
//		//	cell_collider_triangle_index.resize(cell_vector_size_collider);
//		//	cell_triangle_index_address = cell_collider_triangle_index.data();
//		//	for (unsigned int t = 0; t < thread_num; ++t) {
//		//		if (!spatial_hashing_cell_collider[t][cell_index].empty()) {
//		//			memcpy(cell_triangle_index_address, spatial_hashing_cell_collider[t][cell_index].data(),
//		//				4 * spatial_hashing_cell_collider[t][cell_index].size());
//		//			cell_triangle_index_address += spatial_hashing_cell_collider[t][cell_index].size();
//		//		}
//		//	}
//		//}		
//
//		hash_cell_element_size = spatial_hashing_cell[0][cell_index][0];
//		cell_triangle_index = spatial_hashing_cell[0][cell_index] + 1;
//		if (hash_cell_element_size > 2) {
//			for (unsigned int i = 0; i < hash_cell_element_size; i += 2) {
//				obj_index_0 = cell_triangle_index[i + 1]; triangle_index_0 = cell_triangle_index[i];
//				aabb_1 = obj_tri_aabb[obj_index_0][triangle_index_0].data();
//				for (unsigned int j = i + 2; j < hash_cell_element_size; j += 2) {
//					obj_index_1 = cell_triangle_index[j + 1];
//					triangle_index_1 = cell_triangle_index[j];
//					//a++;
//					//
//					//if (AABB::AABB_intersection(aabb_1, obj_tri_aabb[obj_index_1][triangle_index_1].data())) {
//					*(triangle_pair_++) = triangle_index_0;
//					*(triangle_pair_++) = obj_index_0;
//					*(triangle_pair_++) = triangle_index_1;
//					*(triangle_pair_++) = obj_index_1;
//					//}
//				}
//			}
//		}
//		if (has_collider) {
//			hash_cell_element_size_collider = spatial_hashing_cell_collider[0][cell_index][0];
//			cell_collider_triangle_index = spatial_hashing_cell_collider[0][cell_index] + 1;
//			if (hash_cell_element_size > 0) {
//				for (unsigned int i = 0; i < hash_cell_element_size; i += 2) {
//					obj_index_0 = cell_triangle_index[i + 1]; triangle_index_0 = cell_triangle_index[i];
//					aabb_1 = obj_tri_aabb[obj_index_0][triangle_index_0].data();
//					for (unsigned int j = 0; j < hash_cell_element_size_collider; j += 2) {
//						obj_index_1 = cell_collider_triangle_index[j + 1];
//						triangle_index_1 = cell_collider_triangle_index[j];
//						if (AABB::AABB_intersection(aabb_1, collider_tri_aabb[obj_index_1][triangle_index_1].data())) {
//							*(triangle_pair_with_collider_++) = triangle_index_0;
//							*(triangle_pair_with_collider_++) = obj_index_0;
//							*(triangle_pair_with_collider_++) = triangle_index_1;
//							*(triangle_pair_with_collider_++) = obj_index_1;
//						}
//					}
//				}
//			}
//		}
//		//for (unsigned int t0 = 0; t0 < thread_num; ++t0) {
//		//	for (unsigned int k0 = 0; k0 < spatial_hashing_cell[t0][cell_index].size(); k0 += 2) {
//		//		obj_index_0 = spatial_hashing_cell[t0][cell_index][k0+1]; 
//		//		triangle_index_0 = spatial_hashing_cell[t0][cell_index][k0];
//		//		aabb_1 = obj_tri_aabb[obj_index_0][triangle_index_0].data();
//		//		for (unsigned int k1 = k0 + 2; k1 < spatial_hashing_cell[t0][cell_index].size(); k1 += 2) {
//		//			obj_index_1 = spatial_hashing_cell[t0][cell_index][k1 + 1];
//		//			triangle_index_1= spatial_hashing_cell[t0][cell_index][k1];
//		//			//indicator[thread_No] ++;
//		//			if (AABB::AABB_intersection(aabb_1, obj_tri_aabb[obj_index_1][triangle_index_1].data())) {
//		//				triangle_pair_->emplace_back(triangle_index_0);
//		//				triangle_pair_->emplace_back(obj_index_0);
//		//				triangle_pair_->emplace_back(triangle_index_1);
//		//				triangle_pair_->emplace_back(obj_index_1);
//		//				//indicator[thread_No] ++;
//		//			}
//		//		}
//		//		for (unsigned int t1 = t0 + 1; t1 < thread_num; ++t1) {
//		//			for (unsigned int k1 = 0; k1 < spatial_hashing_cell[t1][cell_index].size(); k1 += 2) {
//		//				obj_index_1 = spatial_hashing_cell[t1][cell_index][k1 + 1];
//		//				triangle_index_1 = spatial_hashing_cell[t1][cell_index][k1];
//		//				//indicator[thread_No] ++;
//		//				if (AABB::AABB_intersection(aabb_1, obj_tri_aabb[obj_index_1][triangle_index_1].data())) {
//		//					triangle_pair_->emplace_back(triangle_index_0);
//		//					triangle_pair_->emplace_back(obj_index_0);
//		//					triangle_pair_->emplace_back(triangle_index_1);
//		//					triangle_pair_->emplace_back(obj_index_1);
//		//					//indicator[thread_No] ++;
//		//				}
//		//			}
//		//		}
//		//		if (has_collider) {
//		//			for (unsigned int t1 = 0; t1 < thread_num; ++t1) {
//		//				for (unsigned int k1 = 0; k1 < spatial_hashing_cell_collider[t1][cell_index].size(); k1 += 2) {
//		//					obj_index_1 = spatial_hashing_cell_collider[t1][cell_index][k1 + 1];
//		//					triangle_index_1 = spatial_hashing_cell_collider[t1][cell_index][k1];
//		//					if (AABB::AABB_intersection(aabb_1, collider_tri_aabb[obj_index_1][triangle_index_1].data())) {
//		//						triangle_pair_with_collider_->emplace_back(triangle_index_0);
//		//						triangle_pair_with_collider_->emplace_back(obj_index_0);
//		//						triangle_pair_with_collider_->emplace_back(triangle_index_1);
//		//						triangle_pair_with_collider_->emplace_back(obj_index_1);
//		//					}
//		//				}
//		//			}
//		//		}
//		//	}
//		//}
//	}
//	triangle_pair[0][0] = triangle_pair_ - triangle_pair[0] - 1;
//	if (has_collider) {
//		triangle_pair_with_collider[0][0] = triangle_pair_with_collider_ - triangle_pair_with_collider[0] - 1;
//	}	
//	//triangle_pair_number_thread[thread_No] = a;
//}




//void SpatialHashing::searchTriangle(double* aabb,unsigned int* hash_index, unsigned int hash_cell_size,  unsigned int input_obj_No, unsigned int triangle_index,
//	std::vector<unsigned int>* triangle_pair_, std::vector<unsigned int>* triangle_pair_with_collider_, unsigned int thread_No)
//{
//	std::vector<unsigned int>triangle_index_record;
//	triangle_index_record.reserve(60);
//	bool** obj_is_used_ = obj_is_used[thread_No];
//	
//	unsigned int hash_cell_index_end;
//	unsigned int obj_No; unsigned int current_triangle_index;
//	unsigned int obj_num = cloth->size() + tetrahedron->size();
//	obj_is_used_[input_obj_No][triangle_index] = true;
//	triangle_index_record.emplace_back(triangle_index);
//	triangle_index_record.emplace_back(input_obj_No);
//	for (unsigned int i = 0; i < hash_cell_size; ++i) {
//		for (unsigned int k = 0; k < thread_num; ++k) {
//			for (unsigned int j = 0; j < spatial_hashing_cell[k][hash_index[i]].size(); j+=2) {
//				obj_No = spatial_hashing_cell[k][hash_index[i]][j + 1];
//				current_triangle_index = spatial_hashing_cell[k][hash_index[i]][j];
//				if (!obj_is_used_[obj_No][current_triangle_index]) {
//					obj_is_used_[obj_No][current_triangle_index] = true;
//					triangle_index_record.emplace_back(current_triangle_index);
//					triangle_index_record.emplace_back(obj_No);
//					if (obj_No > input_obj_No ||
//						(obj_No == input_obj_No && current_triangle_index > triangle_index)) {					
//						if (AABB::AABB_intersection(aabb, obj_tri_aabb[obj_No][current_triangle_index].data())) {
//							triangle_pair_->emplace_back(triangle_index);
//							triangle_pair_->emplace_back(input_obj_No);
//							triangle_pair_->emplace_back(current_triangle_index);
//							triangle_pair_->emplace_back(obj_No);
//						}						
//					}			
//				}
//			}			
//		}
//	}	
//	for (unsigned int i = 0; i < triangle_index_record.size(); i += 2) {
//		obj_is_used_[triangle_index_record[i + 1]][triangle_index_record[i]] = false;
//	}
//	
//	if (has_collider) {
//		triangle_index_record.clear();
//		bool** collider_is_used_ = collider_is_used[thread_No];
//		for (unsigned int i = 0; i < hash_cell_size; ++i) {
//			for (unsigned int k = 0; k < thread_num; ++k) {
//				for (unsigned int j = 0; j < spatial_hashing_cell_collider[k][hash_index[i]].size(); j += 2) {
//					obj_No = spatial_hashing_cell_collider[k][hash_index[i]][j + 1];
//					current_triangle_index = spatial_hashing_cell_collider[k][hash_index[i]][j];
//					if (!collider_is_used_[obj_No][current_triangle_index]) {
//						collider_is_used_[obj_No][current_triangle_index] = true;
//						triangle_index_record.emplace_back(current_triangle_index);
//						triangle_index_record.emplace_back(obj_No);
//						if (AABB::AABB_intersection(aabb, collider_tri_aabb[obj_No][current_triangle_index].data())) {
//							triangle_pair_with_collider_->emplace_back(triangle_index);
//							triangle_pair_with_collider_->emplace_back(input_obj_No);
//							triangle_pair_with_collider_->emplace_back(current_triangle_index);
//							triangle_pair_with_collider_->emplace_back(obj_No);
//						}
//					}
//				}
//			}
//		}
//		for (unsigned int i = 0; i < triangle_index_record.size(); i += 2) {
//			collider_is_used_[triangle_index_record[i + 1]][triangle_index_record[i]] = false;
//		}
//	}
//}

//void SpatialHashing::findAllTrianglePairsHashTableTest(int thread_No)
//{
//	unsigned int last_cell_index = cell_begin_per_thread[thread_No + 1];
//
//	std::vector<unsigned int> cell_triangle_index;
//	std::vector<unsigned int> cell_collider_triangle_index;
//
//	cell_triangle_index.reserve(200);
//	cell_collider_triangle_index.reserve(200);
//
//	unsigned int previous_tri_num;
//	unsigned int current_thread_triangle_num;
//
//	std::vector<unsigned int>* triangle_pair_;
//	std::vector<unsigned int>* triangle_pair_with_collider_;
//
//	triangle_pair_ = &triangle_pair[thread_No];
//	triangle_pair_with_collider_ = &triangle_pair_with_collider[thread_No];
//
//	triangle_pair_->clear();
//	triangle_pair_with_collider_->clear();
//
//	double* aabb_1;
//
//	unsigned int obj_index_0, obj_index_1, triangle_index_0, triangle_index_1;
//	triangle_pair_number_thread[thread_No] = 0;
//
//	for (unsigned int cell_index = cell_begin_per_thread[thread_No];
//		cell_index < last_cell_index; ++cell_index) {
//
//		cell_triangle_index.clear();
//		cell_collider_triangle_index.clear();
//
//		for (unsigned int t = 0; t < thread_num; ++t) {
//			if (!spatial_hashing_cell[t][cell_index].empty()) {
//				cell_triangle_index.insert(cell_triangle_index.end(), spatial_hashing_cell[t][cell_index].begin(),
//					spatial_hashing_cell[t][cell_index].end());
//			}
//		}
//		if (has_collider) {
//			for (unsigned int t = 0; t < thread_num; ++t) {
//				if (!spatial_hashing_cell_collider[t][cell_index].empty()) {
//					cell_collider_triangle_index.insert(cell_collider_triangle_index.end(), spatial_hashing_cell_collider[t][cell_index].begin(),
//						spatial_hashing_cell_collider[t][cell_index].end());
//				}
//			}
//		}
//		//there exist collision
//
//		if (cell_triangle_index.size() > 2 || (!cell_collider_triangle_index.empty())) {
//			for (unsigned int i = 0; i < cell_triangle_index.size(); i += 2) {
//				obj_index_0 = cell_triangle_index[i + 1]; triangle_index_0 = cell_triangle_index[i];
//				aabb_1 = obj_tri_aabb[obj_index_0][triangle_index_0].data();
//				for (unsigned int j = i + 2; j < cell_triangle_index.size(); j += 2) {
//					obj_index_1 = cell_triangle_index[j + 1];
//					triangle_index_1 = cell_triangle_index[j];
//					triangle_pair_number_thread[thread_No]++;
//					if (AABB::AABB_intersection(aabb_1, obj_tri_aabb[obj_index_1][triangle_index_1].data())) {
//
//						triangle_pair_->emplace_back(triangle_index_0);
//						triangle_pair_->emplace_back(obj_index_0);
//						triangle_pair_->emplace_back(triangle_index_1);
//						triangle_pair_->emplace_back(obj_index_1);
//					}
//				}
//				if (!cell_collider_triangle_index.empty()) {
//					for (unsigned int j = 0; j < cell_collider_triangle_index.size(); j += 2) {
//						obj_index_1 = cell_collider_triangle_index[j + 1];
//						triangle_index_1 = cell_collider_triangle_index[j];
//						if (AABB::AABB_intersection(aabb_1, collider_tri_aabb[obj_index_1][triangle_index_1].data())) {
//							triangle_pair_with_collider_->emplace_back(triangle_index_0);
//							triangle_pair_with_collider_->emplace_back(obj_index_0);
//							triangle_pair_with_collider_->emplace_back(triangle_index_1);
//							triangle_pair_with_collider_->emplace_back(obj_index_1);
//						}
//					}
//				}
//				//if (triangle_index_0 == 32384) {
//				//	for (unsigned int i = 0; i < cell_triangle_index.size(); i += 2) {
//				//		std::cout << cell_triangle_index[i] << " ";
//				//	}
//				//	std::cout << std::endl;
//				//}
//			}
//		}
//	}
//
//}




void SpatialHashing::findTheMaxNumberOfTriangleInCell()
{
	//unsigned int num1 = 0;
	//unsigned int max_num = 0;
	//unsigned int num = 0;
	//for (unsigned int i = 1; i < prefix_sum[0].size(); ++i) {
	//	num1 = 0;
	//	for (unsigned int j = 0; j < thread_num; ++j) {
	//		num += (prefix_sum[j][i] - prefix_sum[j][i - 1]) >> 1;
	//		num1 += (prefix_sum[j][i] - prefix_sum[j][i - 1]) >> 1;
	//	}
	//	if (num1 > max_num) {
	//		max_num = num1;
	//	}
	//}
	//std::cout << "max_num " << max_num << " " << (double)num / (double)(prefix_sum[0].size() - 1) << std::endl;
}

void SpatialHashing::setPrifixSum()
{
	//std::cout << "running... " << std::endl;
	//for (int i = 1; i < spatial_hashing_value.size(); ++i) {
	//	if (spatial_hashing_value[i] - spatial_hashing_value[i - 1] < 0) {
	//		std::cout << "order error " << spatial_hashing_value[i] << " " << spatial_hashing_value[i - 1] << std::endl;
	//	}
	//	if (spatial_hashing_triangle_index[i]< 0 || spatial_hashing_triangle_index[i]>(*cloth)[0].mesh_struct.triangle_indices.size()) {
	//		std::cout << "triangle index error " << spatial_hashing_triangle_index[i]<< std::endl;
	//	}
	//	if (spatial_hashing_cloth_index[i] != 0) {
	//		std::cout << "cloth index error " << spatial_hashing_cloth_index[i] << std::endl;
	//	}
	//

	////initial prefix sum in TRIANGLE_HASHING, clear() & push_back(0)
	//std::vector<unsigned int*> index_per_thread(thread_num);
	//std::vector<unsigned int*> index_per_thread_collider(thread_num);
	//std::vector<unsigned int> start_index(thread_num, 0);
	//std::vector<unsigned int> start_index_collider(thread_num, 0);
	//std::vector<unsigned int> record_collider_index_start(thread_num, 0);
	////std::vector<unsigned int> record_collider_index_end(thread_num, 0);

	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	index_per_thread[i] = spatial_hashing_value[i];
	//}
	//if (has_collider) {
	//	for (unsigned int i = 0; i < thread_num; ++i) {
	//		index_per_thread_collider[i] = spatial_hashing_value_collider[i];
	//	}
	//}
	//unsigned int current_hash_index = UINT_MAX;
	//bool end_loop;
	//while (true)
	//{
	//	end_loop = (index_per_thread[0] == actual_hash_value_end_index_ref[0]);
	//	for (unsigned int i = 1; i < thread_num; ++i) {
	//		end_loop = end_loop && (index_per_thread[i] == actual_hash_value_end_index_ref[i]);
	//	}
	//	if (end_loop) {
	//		break;
	//	}
	//	current_hash_index = UINT_MAX;
	//	for (unsigned int i = 0; i < thread_num; ++i) {
	//		if (index_per_thread[i] != actual_hash_value_end_index_ref[i]) {
	//			if (current_hash_index > *(index_per_thread[i])) {
	//				current_hash_index = *(index_per_thread[i]);
	//			}
	//		}

	//	}
	//	for (unsigned int i = 0; i < thread_num; ++i) {
	//		prefix_sum[i].emplace_back(prefix_sum[i].back());
	//		if (index_per_thread[i] != actual_hash_value_end_index_ref[i]) {
	//			if (*(index_per_thread[i]) == current_hash_index)
	//			{
	//				prefix_sum[i].back() += hash_index_count_per_thread[i][start_index[i]];
	//				index_per_thread[i] += hash_index_count_per_thread[i][start_index[i]];
	//				start_index[i]++;
	//			}
	//		}
	//	}

	//	if (has_collider)
	//	{
	//		for (unsigned int i = 0; i < thread_num; ++i) {
	//			if (start_index_collider[i] < hash_index_count_per_thread_collider[i].size()) {
	//				while (*(index_per_thread_collider[i]) < current_hash_index)
	//				{
	//					index_per_thread_collider[i] += hash_index_count_per_thread_collider[i][start_index_collider[i]];
	//					record_collider_index_start[i] += hash_index_count_per_thread_collider[i][start_index_collider[i]];
	//					start_index_collider[i]++;
	//					if (start_index_collider[i] == hash_index_count_per_thread_collider[i].size()) {
	//						break;
	//					}
	//				}
	//			}
	//			prefix_sum_collider[i].emplace_back(record_collider_index_start[i]);

	//			if (start_index_collider[i] < hash_index_count_per_thread_collider[i].size()) {
	//				if (*(index_per_thread_collider[i]) == current_hash_index)
	//				{
	//					index_per_thread_collider[i] += hash_index_count_per_thread_collider[i][start_index_collider[i]];
	//					record_collider_index_start[i] += hash_index_count_per_thread_collider[i][start_index_collider[i]];
	//					start_index_collider[i]++;
	//				}
	//			}
	//			prefix_sum_collider[i].emplace_back(record_collider_index_start[i]);
	//		}

	//	}

	//}
	//thread->assignTask(this, PREFIX_SUM_MULTI_2);
	//arrangeIndex(thread_num, prefix_sum[0].size() - 1, actual_exist_cell_begin_per_thread.data());

	//for (unsigned int i = 0; i < prefix_sum[7].size(); ++i) {
	//	std::cout << prefix_sum[7][i] << std::endl;
	//}

	//test();
	//if (prefix_sum[0].size() != prefix_sum_collider[0].size() / 2 + 1) {
	//	std::cout << prefix_sum[0].size() << " " << prefix_sum_collider[0].size() << std::endl;
	//	std::cout << "prefix sum not consistent with prefix sum collider " << std::endl;
	//}

	//testPrifixSum1();

}


//PREFIX_SUM_MULTI_2
void SpatialHashing::prefixSumMulti2(int thread_No)
{
	//std::vector<unsigned int>* prefix = &prefix_sum[thread_No];
	//for (unsigned int i = 0; i < prefix->size(); ++i) {
	//	prefix->data()[i] <<= 1;
	//}
	//if (has_collider) {
	//	std::vector<unsigned int>* prefix_c = &prefix_sum_collider[thread_No];
	//	for (unsigned int i = 0; i < prefix_c->size(); ++i) {
	//		prefix_c->data()[i] <<= 1;
	//	}
	//}
}

//single thread
//void SpatialHashing::setPrifixSum()
//{
//	//initial prefix sum in TRIANGLE_HASHING, clear() & push_back(0)
//	std::vector<unsigned int*> index_per_thread(thread_num);
//	for (unsigned int i = 0; i < thread_num; ++i) {
//		index_per_thread[i] = spatial_hashing_value[i];
//	}
//	unsigned int current_hash_index = 0;
//	bool end_loop;
//	while (true)
//	{
//		end_loop = (index_per_thread[0] == actual_hash_value_end_index_ref[0]);
//		current_hash_index = *(index_per_thread[0]);
//		for (unsigned int i = 1; i < thread_num; ++i) {
//			if (current_hash_index > *(index_per_thread[i])) {
//				current_hash_index = *(index_per_thread[i]);
//			}
//			end_loop = end_loop && (index_per_thread[i] == actual_hash_value_end_index_ref[i]);
//		}
//		if (end_loop) {
//			break;
//		}
//		for (unsigned int i = 0; i < thread_num; ++i) {
//			prefix_sum[i].push_back(prefix_sum[i].back());
//			while (*(index_per_thread[i]) == current_hash_index)
//			{
//				prefix_sum[i].back()++;
//				if (index_per_thread[i] == actual_hash_value_end_index_ref[i]) {
//					break;
//				}
//				index_per_thread[i]++;
//			}
//		}
//	}
//	//testPrifixSum1();
//}

//void SpatialHashing::testPrifixSum1()
//{
	//std::vector<int> prefix_sum_;
	//prefix_sum_.resize(hash_table_size + 1);
	//time_t t1 = clock();
	//memset(prefix_sum_.data(), 0, 4 * prefix_sum_.size());
	//for (int j = 0; j < 1000; ++j) {		
	//	for (int i = 0; i < actual_hash_value_count; ++i) {
	//		prefix_sum_[spatial_hashing_value[i]]++;
	//	}
	//}
	//std::cout << "time sequence sum count " << clock() - t1 << std::endl;
	////std::cout << "run this" << std::endl;
	//t1 = clock();
	//for (int j = 0; j < 1000; ++j) {
	//	memset(prefix_sum_.data(), 0, 4 * prefix_sum_.size());
	//	for (int i = 0; i < actual_hash_value_count; ++i) {
	//		prefix_sum_[spatial_hashing_value[i]]++;
	//	}
	//	int s = 0; int t;
	//	for (int i = 0; i < prefix_sum_.size(); ++i) {
	//		t = s + prefix_sum_[i];
	//		prefix_sum_[i] = s;
	//		s = t;
	//	}
	//	prefix_sum_[prefix_sum_.size() - 1] = t;
	//}
	//std::cout << "total time sequence sum count " << clock() - t1 << std::endl;
	//for (int i = 0; i < prefix_sum_.size(); ++i) {
	//	if (prefix_sum_[i] != prefix_sum[i]) {
	//		std::cout << "prefix sum error " << prefix_sum_[i] << " " << prefix_sum[i] << std::endl;
	//	}
	//}
//}

//void SpatialHashing::prepareForActualHashValueCount()
//{
//	unsigned int* value_address;
//	unsigned int k;
//	for (unsigned int i = 1; i < thread_num; ++i) {
//		k = 1;
//		value_address = spatial_hashing_value + actual_hash_count_start_per_thread[i];
//		while (*(value_address - k) == *(value_address)) {
//			k++;
//		}
//		if (k < actual_hash_count_start_per_thread[i]) {
//			actual_hash_count_start_per_thread[i] -= k - 1;
//		}
//		else {
//			actual_hash_count_start_per_thread[i] = 0;
//		}
//	}
//}

// PREFIX_SUM_THREAD_1
//void SpatialHashing::prifixSum2(int thread_No)
//{
//	unsigned int end = total_hash_count_start_per_thread[thread_No + 1];
//	for (unsigned int i = total_hash_count_start_per_thread[thread_No]+1; i < end; ++i) {
//		*(prefix_sum_1_address + i) += *(prefix_sum_1_address + i - 1);
//	}
//}

// PREFIX_SUM_THREAD_2
//void SpatialHashing::prifixSum3(int thread_No)
//{
//	unsigned int start_index = prefix_sum_thread_start[thread_No];
//	unsigned int end = total_hash_count_start_per_thread[thread_No + 1];
//	//*(prefix_sum_1_address + total_hash_count_start_per_thread[thread_No]) =
//	//	start_index + hash_value_count_start_thread[thread_No];
//	for (unsigned int i = total_hash_count_start_per_thread[thread_No]; i < end; ++i) {
//		*(prefix_sum_1_address + i) += start_index;
//	}
//}



//PREPARE_FOR_ACTUAL_HASH_VALUE_COUNT_THREAD
//void SpatialHashing::prepareForActualHashValueCountThread(int thread_No)
//{
//	if (thread_No == 0) {
//		return;
//	}
//	unsigned int* value_address;
//	unsigned int k=1;
//	value_address = spatial_hashing_value + actual_hash_count_start_per_thread[thread_No];
//	while (*(value_address - k) == *(value_address)) {
//		k++;
//	}
//	if (k < actual_hash_count_start_per_thread[thread_No]) {
//		actual_hash_count_start_per_thread[thread_No] -= k - 1;
//	}
//	else {
//		actual_hash_count_start_per_thread[thread_No] = 0;
//	}
//}

// ADD_COUNT_FOR_prefix_sum
//void SpatialHashing::prifixSum1(int thread_No)
//{
//	unsigned int end = actual_hash_count_start_per_thread[thread_No + 1];
//	for (unsigned int i = actual_hash_count_start_per_thread[thread_No]; i < end; ++i) {
//		(*(prefix_sum_1_address + spatial_hashing_value[i]))++;
//	}
//}

// MEMSET_PREFIX
//void SpatialHashing::memsetThread(int thread_No)
//{
//	memset(prefix_sum_1_address + total_hash_count_start_per_thread[thread_No], 0, 4 * (total_hash_count_start_per_thread[thread_No + 1] - total_hash_count_start_per_thread[thread_No]));
//}





void SpatialHashing::setHashTogether()
{
	//unsigned int total_hash_num = spatial_hashing_value_per_thread[0].size();
	//for (unsigned int i = 1; i < thread_num; ++i) {
	//	hash_value_begin[i] = total_hash_num;
	//	total_hash_num+= spatial_hashing_value_per_thread[i].size();
	//	
	//}
	//spatial_hashing_value.resize(total_hash_num);
	//spatial_hashing_obj_index.resize(total_hash_num);
	//spatial_hashing_triangle_index.resize(total_hash_num);
	//thread->assignTask(this, SET_HASH_TOGETHER);
}

//SET_HASH_TOGETHER
void SpatialHashing::setHashTogether(int thread_No)
{
	//memcpy(&spatial_hashing_value[hash_value_begin[thread_No]], spatial_hashing_value_per_thread[thread_No].data(), 4 * spatial_hashing_value_per_thread[thread_No].size());
	//memcpy(&spatial_hashing_obj_index[hash_value_begin[thread_No]], spatial_hashing_obj_index_per_thread[thread_No].data(), 4 * spatial_hashing_value_per_thread[thread_No].size());
	//memcpy(&spatial_hashing_triangle_index[hash_value_begin[thread_No]], spatial_hashing_triangle_index_per_thread[thread_No].data(), 4 * spatial_hashing_value_per_thread[thread_No].size());
}

//PATCH_TRIANGLE_HASHING
//void SpatialHashing::patchTriangleHashing(int thread_No)
//{
//
//}

//TRIANGLE_HASHING
void SpatialHashing::triangleHashing(int thread_No)
{
	//unsigned int* triangle_begin;
	//std::array<double, 6>* aabb;
	//int vector_size;
	//unsigned int triangle_index_total;
	//unsigned int triangle_end;
	//unsigned int largest_count_in_hash_value_list;
	//unsigned int* spatial_hashing_triangle_index_ = spatial_hashing_triangle_index[thread_No];
	//unsigned int* spatial_hashing_value_ = spatial_hashing_value[thread_No];
	////unsigned int* spatial_hashing_obj_index_ = spatial_hashing_obj_index[thread_No];
	//memset(spatial_hashing_value_, -1, 4 * actual_hash_value_count_per_thread[thread_No]);
	//for (unsigned int i = 0; i < cloth->size(); ++i) {
	//	aabb = (*cloth)[i].triangle_AABB.data();
	//	triangle_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
	//	triangle_end = triangle_begin[thread_No + 1];
	//	for (unsigned int j = triangle_begin[thread_No]; j < triangle_end; ++j) {
	//		triangle_index_total = max_cell_count * (triangle_begin_per_obj_per_thread[thread_No][i] + j);
	//		triangleHashValue(aabb[j].data(), spatial_hashing_triangle_index_ + 2 * triangle_index_total,
	//			spatial_hashing_value_ + triangle_index_total, j, i);
	//	}
	//}
	//unsigned int obj_No;
	//for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
	//	obj_No = i + tetrahedron_begin_obj_index;
	//	aabb = (*tetrahedron)[i].triangle_AABB.data();
	//	triangle_begin = (*tetrahedron)[i].mesh_struct.face_index_begin_per_thread.data();
	//	triangle_end = triangle_begin[thread_No + 1];
	//	for (unsigned int j = triangle_begin[thread_No]; j < triangle_end; ++j) {
	//		triangle_index_total = max_cell_count * (triangle_begin_per_obj_per_thread[thread_No][obj_No] + j);
	//		triangleHashValue(aabb[j].data(), spatial_hashing_triangle_index_ + 2 * triangle_index_total,
	//			spatial_hashing_value_ + triangle_index_total, j, obj_No);
	//	}
	//}
	//radix_sort[thread_No].radixSort(hash_table_size, spatial_hashing_value_, spatial_hashing_triangle_index_,
	//	largest_count_in_hash_value_list);
	//testIfRadixSortIsRight(spatial_hashing_value_, total_hash_size_per_thread[thread_No]);
	//actual_hash_value_count_per_thread[thread_No] = total_hash_size_per_thread[thread_No];
	//actual_hash_value_end_index_ref[thread_No] = spatial_hashing_value[thread_No] + total_hash_size_per_thread[thread_No];
	//for (unsigned int i = total_hash_size_per_thread[thread_No] - largest_count_in_hash_value_list;
	//	i < total_hash_size_per_thread[thread_No]; ++i) {
	//	if (spatial_hashing_value_[i] == UINT_MAX) {
	//		actual_hash_value_end_index_ref[thread_No] = spatial_hashing_value[thread_No] + i;
	//		actual_hash_value_count_per_thread[thread_No] = i;
	//		break;
	//	}
	//}
	//prefix_sum[thread_No].clear();
	//prefix_sum[thread_No].push_back(0);
	//countHashIndexPerThread(thread_No, actual_hash_value_count_per_thread[thread_No], &hash_index_count_per_thread[thread_No],
	//	spatial_hashing_value_);
	//if (has_collider) {
	//	unsigned int* spatial_hashing_triangle_index_collider_ = spatial_hashing_triangle_index_collider[thread_No];
	//	unsigned int* spatial_hashing_value_collider_ = spatial_hashing_value_collider[thread_No];
	//	//unsigned int* spatial_hashing_obj_index_collider_ = spatial_hashing_obj_index_collider[thread_No];
	//	memset(spatial_hashing_value_collider_, -1, 4 * actual_hash_value_count_per_thread_collider[thread_No]);
	//	for (unsigned int i = 0; i < collider->size(); ++i) {
	//		aabb = (*collider)[i].triangle_AABB.data();
	//		triangle_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread.data();
	//		triangle_end = triangle_begin[thread_No + 1];
	//		for (unsigned int j = triangle_begin[thread_No]; j < triangle_end; ++j) {
	//			triangle_index_total = max_cell_count * (triangle_begin_per_obj_per_thread_collider[thread_No][i] + j);
	//			triangleHashValue(aabb[j].data(), spatial_hashing_triangle_index_collider_ + 2 * triangle_index_total,
	//				spatial_hashing_value_collider_ + triangle_index_total, j, i);
	//		}
	//	}
	//	radix_sort_collider[thread_No].radixSort(hash_table_size, spatial_hashing_value_collider_, spatial_hashing_triangle_index_collider_,
	//		largest_count_in_hash_value_list);
	//	testIfRadixSortIsRight(spatial_hashing_value_collider_, total_hash_size_per_thread_collider[thread_No]);
	//	actual_hash_value_count_per_thread_collider[thread_No] = total_hash_size_per_thread_collider[thread_No];
	//	actual_hash_value_end_index_ref_collider[thread_No] = spatial_hashing_value_collider[thread_No] + total_hash_size_per_thread_collider[thread_No];
	//	for (unsigned int i = total_hash_size_per_thread_collider[thread_No] - largest_count_in_hash_value_list;
	//		i < total_hash_size_per_thread_collider[thread_No]; ++i) {
	//		if (spatial_hashing_value_collider_[i] == UINT_MAX) {
	//			actual_hash_value_end_index_ref_collider[thread_No] = spatial_hashing_value_collider[thread_No] + i;
	//			actual_hash_value_count_per_thread_collider[thread_No] = i;
	//			break;
	//		}
	//	}
	//	prefix_sum_collider[thread_No].clear();
	//	countHashIndexPerThread(thread_No, actual_hash_value_count_per_thread_collider[thread_No],
	//		&hash_index_count_per_thread_collider[thread_No], spatial_hashing_value_collider_);
	//}

}

//TET_HASHING_SMALLER_HASH_TABLE
void SpatialHashing::tetHashingSmallerHashTable(int thread_No)
{
	//double scene_aabb_[6];
	//memcpy(scene_aabb_, scene_aabb, 48);
	//double hash_cell_length = cell_length;
	//unsigned int hash_cell_count_ = hash_cell_count;
	//uint64_t p1 = P1;
	//uint64_t p2 = P2;
	//uint64_t p3 = P3;
	//unsigned int* triangle_begin;
	//std::array<double, 6>* aabb;
	//int vector_size;
	//unsigned int triangle_end;
	//unsigned int largest_count_in_hash_value_list;
	//unsigned int** spatial_hashing_cell_ = spatial_hashing_cell[thread_No];
	//unsigned int** spatial_hashing_cell_collider_;
	//unsigned int* spatial_hashing_cell_triangle_size_ = spatial_hashing_cell_triangle_size[thread_No];
	//memset(spatial_hashing_cell_triangle_size_, 0, hash_cell_count_ << 2);
	////for (unsigned int i = 0; i < hash_cell_count; ++i) {
	////	spatial_hashing_cell_[i][0] = 0;
	////}
	////std::cout << hash_cell_count << std::endl;
	//unsigned int triangle_index[2];//first triangle_index, second obj index
	//for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
	//	triangle_index[1] = i + tetrahedron_begin_obj_index;
	//	aabb = (*tetrahedron)[i].tet_AABB.data();
	//	triangle_begin = (*tetrahedron)[i].mesh_struct.tetrahedron_index_begin_per_thread.data();
	//	triangle_end = triangle_begin[thread_No + 1];
	//	//std::cout << obj_No << " " << triangle_end << std::endl;
	//	for (unsigned int j = triangle_begin[thread_No]; j < triangle_end; ++j) {
	//		triangle_index[0] = j;
	//		triangleHashValueWithoutRecord(aabb[j].data(), spatial_hashing_cell_, triangle_index, scene_aabb_,
	//			hash_cell_length, hash_cell_count_, p1, p2, p3, spatial_hashing_cell_triangle_size_);
	//	}
	//}
}


void SpatialHashing::testSpatialHashingCellSet0()
{
	//unsigned int* spatial_hashing_cell_ = spatial_hashing_cell_triangle_size[0];
	//for (unsigned int i = 0; i < hash_cell_count; ++i) {
	//	if (spatial_hashing_cell_[i][0]) {
	//		spatial_hashing_cell_[i][0] = 0;
	//	}
	//}
	memset(spatial_hashing_cell_triangle_size[0], 0, hash_cell_count << 2);
}



void SpatialHashing::clearHashCell()
{
	for (unsigned int i = 0; i < thread_num; ++i) {
		memset(spatial_hashing_cell_triangle_size[i], 0, hash_cell_count << 2);
	}
	if (has_collider) {
		for (unsigned int i = 0; i < thread_num; ++i) {
			memset(spatial_hashing_cell_collider_triangle_size[i], 0, hash_cell_count << 2);
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
	unsigned int** spatial_hashing_cell_ = spatial_hashing_cell_triangle[thread_No];
	unsigned int* spatial_hashing_cell_size_ = spatial_hashing_cell_triangle_size[thread_No];
	memset(spatial_hashing_cell_size_, 0, hash_cell_count_ << 2);
	memset(spatial_hashing_cell_edge_size[thread_No], 0, hash_cell_count_ << 2);
	memset(spatial_hashing_cell_vertex_size[thread_No], 0, hash_cell_count_ << 2);

	unsigned int** spatial_hashing_cell_collider_;
	unsigned int* spatial_hashing_cell_collider_size_;

	unsigned int cloth_size = cloth->size();


	unsigned int primitive_index[2];//first triangle_index, second obj index
	for (unsigned int i = 0; i < collider_begin_obj_index_; ++i) {
		aabb = obj_tri_aabb[i];
		primitive_begin = obj_triangle_index_begin_per_thread[i][thread_No];
		primitive_end = obj_triangle_index_begin_per_thread[i][thread_No + 1];
		primitive_index[1] = i;
		for (unsigned int j = primitive_begin; j < primitive_end; ++j) {
			primitive_index[0] = j;
			triangleHashValueWithoutRecord(aabb[j].data(), spatial_hashing_cell_, primitive_index, scene_aabb_, hash_cell_length,
				hash_cell_count_, p1, p2, p3, spatial_hashing_cell_size_);
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
					hash_cell_length, hash_cell_count_, p1, p2, p3, spatial_hashing_cell_collider_size_);
			}
		}
	}

	//edge 
	spatial_hashing_cell_ = spatial_hashing_cell_edge[thread_No];
	spatial_hashing_cell_size_ = spatial_hashing_cell_edge_size[thread_No];
	//memset(spatial_hashing_cell_size_, 0, hash_cell_count_ << 2);


	


	for (unsigned int i = 0; i < collider_begin_obj_index_; ++i) {
		aabb = obj_edge_aabb[i];
		primitive_begin = obj_edge_index_begin_per_thread[i][thread_No];
		primitive_end = obj_edge_index_begin_per_thread[i][thread_No + 1];
		primitive_index[1] = i;
		for (unsigned int j = primitive_begin; j < primitive_end; ++j) {
			primitive_index[0] = j;
			triangleHashValueWithoutRecord(aabb[j].data(), spatial_hashing_cell_, primitive_index, scene_aabb_, hash_cell_length,
				hash_cell_count_, p1, p2, p3, spatial_hashing_cell_size_);
		}
	}
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
					hash_cell_length, hash_cell_count_, p1, p2, p3, spatial_hashing_cell_collider_size_);
			}
		}
	}
	//vertex
	spatial_hashing_cell_ = spatial_hashing_cell_vertex[thread_No];
	spatial_hashing_cell_size_ = spatial_hashing_cell_vertex_size[thread_No];
	//memset(spatial_hashing_cell_size_, 0, hash_cell_count_ << 2);

	unsigned int* vertex_index_on_surface;

	for (unsigned int i = 0; i < collider_begin_obj_index_; ++i) {
		aabb = obj_vertex_aabb[i];
		primitive_begin = obj_vertex_index_begin_per_thread[i][thread_No];
		primitive_end = obj_vertex_index_begin_per_thread[i][thread_No + 1];
		primitive_index[1] = i;

		if (i < cloth_size) {
			for (unsigned int j = primitive_begin; j < primitive_end; ++j) {
				primitive_index[0] = j;
				triangleHashValueWithoutRecord(aabb[j].data(), spatial_hashing_cell_, primitive_index, scene_aabb_, hash_cell_length,
					hash_cell_count_, p1, p2, p3, spatial_hashing_cell_size_);
			}
		}
		else {
			vertex_index_on_surface = tetrahedron_vertex_index_on_surface[i - cloth_size];
			for (unsigned int j = primitive_begin; j < primitive_end; ++j) {
				primitive_index[0] = vertex_index_on_surface[j];
				triangleHashValueWithoutRecord(aabb[vertex_index_on_surface[j]].data(), spatial_hashing_cell_, primitive_index, scene_aabb_, hash_cell_length,
					hash_cell_count_, p1, p2, p3, spatial_hashing_cell_size_);
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
				triangleHashValueWithoutRecord(aabb[j].data(), spatial_hashing_cell_collider_, primitive_index, scene_aabb_,
					hash_cell_length, hash_cell_count_, p1, p2, p3, spatial_hashing_cell_collider_size_);
			}
		}
	}
}


//TRIANGLE_HASHING_RECORD_REAL_HASH_VALUE
void SpatialHashing::recordRealTriangleHashValue(int thread_No)
{
	//unsigned int* triangle_begin;
	//std::array<double, 6>* aabb;
	//int vector_size;
	//unsigned int triangle_end;
	//unsigned int largest_count_in_hash_value_list;
	//std::vector<unsigned int>* spatial_hashing_cell_ = spatial_hashing_actual_hash_value_for_test[thread_No];
	//for (unsigned int i = 0; i < hash_cell_count; ++i) {
	//	spatial_hashing_cell_[i].clear();
	//}
	//for (unsigned int i = 0; i < cloth->size(); ++i) {
	//	aabb = (*cloth)[i].triangle_AABB.data();
	//	triangle_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
	//	triangle_end = triangle_begin[thread_No + 1];
	//	for (unsigned int j = triangle_begin[thread_No]; j < triangle_end; ++j) {
	//		recordRealTriangleHashValue(aabb[j].data(), spatial_hashing_cell_, j, i);
	//	}
	//}
	//unsigned int obj_No;
	//for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
	//	obj_No = i + tetrahedron_begin_obj_index;
	//	aabb = (*tetrahedron)[i].triangle_AABB.data();
	//	triangle_begin = (*tetrahedron)[i].mesh_struct.face_index_begin_per_thread.data();
	//	triangle_end = triangle_begin[thread_No + 1];
	//	for (unsigned int j = triangle_begin[thread_No]; j < triangle_end; ++j) {
	//		recordRealTriangleHashValue(aabb[j].data(), spatial_hashing_cell_, j, obj_No);
	//	}
	//}
}

void SpatialHashing::test()
{

	for (unsigned int i = 0; i < thread_num; ++i) {
		//std::cout << actual_hash_value_count_per_thread[i] << std::endl;
		//for (unsigned int j = 0; j < actual_hash_value_count_per_thread[i]+2; ++j) {
		//	std::cout <<j<<" "<< spatial_hashing_value[i][j] << " " << spatial_hashing_triangle_index[i][j] << " " << std::endl;
		//}
		//std::cout << prefix_sum[i].size() << std::endl;
		//for (unsigned int j = 0; j < prefix_sum[i].size(); ++j) {
		//	std::cout<< prefix_sum[i][j]<<" ";
		//}
		//std::cout << std::endl;
		//std::cout << actual_exist_cell_begin_per_thread[i] << std::endl;
	}

}


void SpatialHashing::testIfRadixSortIsRight(unsigned int* has_value, unsigned int count)
{
	for (unsigned int i = 1; i < count; ++i) {
		if (has_value[i] < has_value[i - 1]) {
			std::cout << "hash value is not in increasing order" << std::endl;
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

void SpatialHashing::obtainTriangleHashingValue(double* aabb, std::vector<unsigned int>* hash_value)
{
	unsigned int spatial_index[6];
	unsigned int index_z_multi;
	unsigned int index_y_multi;
	for (unsigned int j = 0; j < 3; ++j) {
		spatial_index[j] = (unsigned int)floor((aabb[j] - scene_aabb[j]) / cell_length);
	}
	for (unsigned int j = 3; j < 6; ++j) {
		spatial_index[j] = (unsigned int)floor((aabb[j] - scene_aabb[j - 3]) / cell_length) + 1;
	}
	unsigned int size = (spatial_index[3] - spatial_index[0]) * (spatial_index[4] - spatial_index[1]) * (spatial_index[5] - spatial_index[2]);
	//hash_value.clear();
	hash_value->resize(size);
	unsigned int* value = hash_value->data();
	for (unsigned int index_y = spatial_index[1]; index_y < spatial_index[4]; ++index_y) {
		index_y_multi = index_y * cell_number[0];
		for (unsigned int index_z = spatial_index[2]; index_z < spatial_index[5]; ++index_z) {
			index_z_multi = index_z * cell_num0_cell_num1;
			for (unsigned int index_x = spatial_index[0]; index_x < spatial_index[3]; ++index_x) {
				*(value++) = index_x + index_y_multi + index_z_multi;
			}
		}
	}
}

void SpatialHashing::triangleHashingValue(double* aabb, std::vector<unsigned int>* spatial_hashing_triangle_index,
	std::vector<unsigned int>* spatial_hashing_value, unsigned int triangle_index, std::vector<unsigned int>* hash_value)
{
	unsigned int spatial_index[6];
	unsigned int index_z_multi;
	unsigned int index_y_multi;
	for (unsigned int j = 0; j < 3; ++j) {
		spatial_index[j] = (unsigned int)floor((aabb[j] - scene_aabb[j]) / cell_length);
	}
	for (unsigned int j = 3; j < 6; ++j) {
		spatial_index[j] = (unsigned int)floor((aabb[j] - scene_aabb[j - 3]) / cell_length) + 1;
	}
	//int size = (max_index[0] - min_index[0]) * (max_index[1] - min_index[1]) * (max_index[2] - min_index[2]);
	hash_value->clear();
	unsigned int value;
	for (unsigned int index_y = spatial_index[1]; index_y < spatial_index[4]; ++index_y) {
		index_y_multi = index_y * cell_number[0];
		for (unsigned int index_z = spatial_index[2]; index_z < spatial_index[5]; ++index_z) {
			index_z_multi = index_z * cell_num0_cell_num1;
			for (unsigned int index_x = spatial_index[0]; index_x < spatial_index[3]; ++index_x) {
				value = index_x + index_y_multi + index_z_multi;
				hash_value->push_back(value);
			}
		}
	}
	spatial_hashing_value->insert(spatial_hashing_value->end(), hash_value->begin(), hash_value->end());
	spatial_hashing_triangle_index->insert(spatial_hashing_triangle_index->end(), hash_value->size(), triangle_index);
}

void SpatialHashing::triangleHashValue(double* aabb, unsigned int* spatial_hashing_triangle_index,
	unsigned int* spatial_hashing_value, unsigned int triangle_index, unsigned int obj_index)
{
	unsigned int spatial_index[6];
	unsigned int index_z_multi;
	unsigned int index_y_multi;
	for (unsigned int j = 0; j < 3; ++j) {
		spatial_index[j] = (unsigned int)floor((aabb[j] - scene_aabb[j]) / cell_length);
	}
	for (unsigned int j = 3; j < 6; ++j) {
		spatial_index[j] = (unsigned int)floor((aabb[j] - scene_aabb[j - 3]) / cell_length) + 1;
	}
	//if ((spatial_index[4] - spatial_index[1]) * (spatial_index[5] - spatial_index[2]) * (spatial_index[3] - spatial_index[0]) > max_cell_count) {
	//	std::cout << aabb[0] << " " << aabb[1] << " " << aabb[2] << std::endl;
	//	std::cout << scene_aabb[0] << " " << scene_aabb[1] << " " << scene_aabb[2] << std::endl;
	//	std::cout << "cell count is not large enough " << (spatial_index[4] - spatial_index[1]) * (spatial_index[5] - spatial_index[2]) * (spatial_index[3] - spatial_index[0]) << std::endl;
	//}

	for (unsigned int index_y = spatial_index[1]; index_y < spatial_index[4]; ++index_y) {
		index_y_multi = index_y * cell_number[0];
		for (unsigned int index_z = spatial_index[2]; index_z < spatial_index[5]; ++index_z) {
			index_z_multi = index_z * cell_num0_cell_num1;
			for (unsigned int index_x = spatial_index[0]; index_x < spatial_index[3]; ++index_x) {
				*(spatial_hashing_value++) = index_x + index_y_multi + index_z_multi;
				*(spatial_hashing_triangle_index++) = triangle_index;
				*(spatial_hashing_triangle_index++) = obj_index;
				//*(spatial_hashing_obj_index++) = obj_index;
			}
		}
	}

}



void SpatialHashing::triangleHashValueWithRecord(double* aabb,
	std::vector<unsigned int>* spatial_hashing_cell, unsigned int triangle_index, unsigned int obj_index)
{
	//std::uint64_t spatial_index[6];
	//std::uint64_t index_z_multi;
	//std::uint64_t index_y_multi;
	//unsigned int hash_index;
	//for (unsigned int j = 0; j < 3; ++j) {
	//	spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j]) / cell_length);
	//}
	//for (unsigned int j = 3; j < 6; ++j) {
	//	spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j - 3]) / cell_length) + 1;
	//}
	//std::vector<unsigned int>* obj_triangle_hash_ = &obj_triangle_hash[obj_index][triangle_index];
	//obj_triangle_hash_->clear();
	//for (std::uint64_t index_y = spatial_index[1]; index_y < spatial_index[4]; ++index_y) {
	//	index_y_multi = P2 * index_y;
	//	for (std::uint64_t index_z = spatial_index[2]; index_z < spatial_index[5]; ++index_z) {
	//		index_z_multi = index_z * P3;
	//		for (std::uint64_t index_x = spatial_index[0]; index_x < spatial_index[3]; ++index_x) {
	//			hash_index = ((index_x * P1) ^ index_y_multi ^ index_z_multi) % hash_cell_count;
	//			if (spatial_hashing_cell[hash_index].back() != obj_index
	//				|| *(spatial_hashing_cell[hash_index].end() - 2) != triangle_index) {
	//				spatial_hashing_cell[hash_index].emplace_back(triangle_index);
	//				spatial_hashing_cell[hash_index].emplace_back(obj_index);
	//				obj_triangle_hash_->emplace_back(hash_index);
	//			}
	//		}
	//	}
	//}
}

void SpatialHashing::triangleHashValueWithoutRecord(double* aabb,
	unsigned int** spatial_hashing_cell, unsigned int* triangle_index, double* scene_aabb, double cell_length,
	unsigned int hash_cell_count, uint64_t P1, uint64_t P2, uint64_t P3, unsigned int* spatial_hashing_cell_triangle_size)
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

	for (std::uint64_t index_y = spatial_index[1]; index_y < spatial_index[4]; ++index_y) {
		//index_y_multi = P2 * index_y;
		for (std::uint64_t index_z = spatial_index[2]; index_z < spatial_index[5]; ++index_z) {
			//index_z_multi = index_z * P3;
			for (std::uint64_t index_x = spatial_index[0]; index_x < spatial_index[3]; ++index_x) {
				hash_index = ((index_x * P1) ^ (P2 * index_y) ^ (index_z * P3)) % hash_cell_count;
				//cell_size_index = spatial_hashing_cell[hash_index];
				if (spatial_hashing_cell_triangle_size[hash_index]) {
					if (memcmp(spatial_hashing_cell[hash_index] + spatial_hashing_cell_triangle_size[hash_index] - 2, triangle_index, 8) != 0) {
						memcpy(spatial_hashing_cell[hash_index] + spatial_hashing_cell_triangle_size[hash_index], triangle_index, 8);
						spatial_hashing_cell_triangle_size[hash_index] += 2;
						//cell_size_index[*cell_size_index] = triangle_index;
						//(*cell_size_index)++;
						//cell_size_index[*cell_size_index] = obj_index;
					}
				}
				else {
					memcpy(spatial_hashing_cell[hash_index] + spatial_hashing_cell_triangle_size[hash_index], triangle_index, 8);
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

void SpatialHashing::recordRealTriangleHashValue(double* aabb,
	std::vector<unsigned int>* spatial_hashing_cell, unsigned int triangle_index, unsigned int obj_index)
{
	std::uint64_t spatial_index[6];
	std::uint64_t index_z_multi;
	std::uint64_t index_y_multi;

	std::uint64_t index_z_multi1;
	std::uint64_t index_y_multi1;

	unsigned int hash_index;
	for (unsigned int j = 0; j < 3; ++j) {
		spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j]) / cell_length);
	}
	for (unsigned int j = 3; j < 6; ++j) {
		spatial_index[j] = (std::uint64_t)floor((aabb[j] - scene_aabb[j - 3]) / cell_length) + 1;
	}
	for (std::uint64_t index_y = spatial_index[1]; index_y < spatial_index[4]; ++index_y) {
		index_y_multi = P2 * index_y;
		index_y_multi1 = index_y * cell_number[0];
		for (std::uint64_t index_z = spatial_index[2]; index_z < spatial_index[5]; ++index_z) {
			index_z_multi = index_z * P3;
			index_z_multi1 = index_z * cell_num0_cell_num1;
			for (std::uint64_t index_x = spatial_index[0]; index_x < spatial_index[3]; ++index_x) {
				hash_index = ((index_x * P1) ^ index_y_multi ^ index_z_multi) % hash_cell_count;
				spatial_hashing_cell[hash_index].emplace_back(index_x + index_y_multi1 + index_z_multi1);
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

void SpatialHashing::testHashCollision()
{
	unsigned int count = 0;
	unsigned int not_empty_count = 0;
	for (unsigned int i = 0; i < hash_cell_count; ++i) {
		if (hashCellIsNotEmpty(i)) {
			if (hashCollisionHappened(i)) {
				count++;
			}
			not_empty_count++;
		}

	}
	std::cout << "hash collision " << count << " not emput count " << not_empty_count << " " << (double)count / (double)not_empty_count << " total count " << hash_cell_count << std::endl;

	//for (unsigned int j = 0; j < hash_cell_count; ++j) {
	//	count = 0;
	//	for (unsigned int i = 0; i < thread_num; ++i) {
	//		count += spatial_hashing_cell[i][j].size();
	//	}
	//	std::cout << "count " << count << std::endl;
	//}
	//
}


void SpatialHashing::testHashCollisionSP()
{
	//unsigned int count = 0;
	//unsigned int not_empty_count = 0;
	//for (unsigned int j = 0; j < prefix_sum[0].size() - 1; ++j) {
	//	if (hashCellIsNotEmptyTest(j)) {
	//		if (hashCollisionHappenedTest(j)) {
	//			count++;
	//		}
	//		not_empty_count++;
	//	}
	//}
	//std::cout << "hash collision " << count << " not emput count " << not_empty_count << " " << (double)count / (double)not_empty_count << std::endl;
}

bool SpatialHashing::hashCollisionHappened(unsigned int hash_index)
{
	//int record_index = -1;
	//for (unsigned int j = 0; j < thread_num; ++j) {
	//	for (unsigned int k = 0; k < spatial_hashing_actual_hash_value_for_test[j][hash_index].size(); k++) {
	//		if ((int)spatial_hashing_actual_hash_value_for_test[j][hash_index][k] != record_index) {
	//			if (record_index == -1) {
	//				record_index = spatial_hashing_actual_hash_value_for_test[j][hash_index][k];
	//			}
	//			else {
	//				//std::cout << record_index << " " << spatial_hashing_actual_hash_value_for_test[j][hash_index][k] << std::endl;
	//				return true;
	//			}
	//		}
	//	}
	//}
	return false;
}

bool SpatialHashing::hashCollisionHappenedTest(unsigned int hash_index)
{
	//int record_index = -1;
	//for (unsigned int j = 0; j < thread_num; ++j) {
	//	for (unsigned int k = prefix_sum[j][hash_index] + 1; k < prefix_sum[j][hash_index + 1]; k += 2) {
	//		if ((int)spatial_hashing_triangle_index[j][k] != record_index) {
	//			if (record_index == -1) {
	//				record_index = spatial_hashing_triangle_index[j][k];
	//			}
	//			else {
	//				//std::cout << record_index << " " << spatial_hashing_actual_hash_value_for_test[j][hash_index][k] << std::endl;
	//				return true;
	//			}
	//		}
	//	}
	//}
	//return false;
}


bool SpatialHashing::hashCellIsNotEmpty(unsigned int hash_index)
{
	//for (unsigned int j = 0; j < thread_num; ++j) {
	//	if (!spatial_hashing_actual_hash_value_for_test[j][hash_index].empty()) {
	//		return true;
	//	}		
	//}
	return false;
}

bool SpatialHashing::hashCellIsNotEmptyTest(unsigned int hash_index)
{
	//for (unsigned int j = 0; j < thread_num; ++j) {
	//	if ((prefix_sum[j][hash_index + 1] - prefix_sum[j][hash_index]) > 0) {
	//		return true;
	//	}
	//}
	//return false;
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

	//unsigned int pair_count_ = 0;
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	pair_count_ += pair_count_per_thread[i];
	//	std::cout << pair_count_per_thread[i] << " ";
	//}
	//std::cout << std::endl;
	//std::cout << "total pair count " << pair_count_ << std::endl;

	for (unsigned int i = 1; i < thread_num; ++i) {
		prefix_sum_record_per_thread_start_vertex_triangle[i] = prefix_sum_record_per_thread_start_vertex_triangle[i - 1] + hash_cell_pair_num_prefix_vertex_triangle[non_empty_cell_index_begin_per_thread_vertex_triangle[i] + 1];
		prefix_sum_record_per_thread_start_edge[i] = prefix_sum_record_per_thread_start_edge[i - 1] + hash_cell_pair_num_prefix_edge[non_empty_cell_index_begin_per_thread_edge[i] + 1];
	}

	thread->assignTask(this, SET_HASH_CELL_PAIR_NUM_PREFIX_SUM_TOGETHER);

	//std::cout << "prefix sum pair " << hash_cell_pair_num_prefix[non_empty_cell_index[0][0]+1]<<std::endl;

	thread->assignTask(this, SET_PAIR_AVE);
	//std::cout << "--------" << std::endl;
	//setGlobalCellStartPerThreadPerLoop();

	//std::cout << std::endl;
	//std::cout << "cell start " << std::endl;
	//for (unsigned int i = 0; i <= thread_num ; ++i) {
	//	std::cout << non_empty_cell_index_ave_begin_per_thread[i] << std::endl;
	//}
	//for (unsigned int i = 0; i < thread_num*2; ++i) {
	//	std::cout << global_cell_start[i] << " ";
	//}
	//std::cout << std::endl;
}


//void SpatialHashing::setGlobalCellStartPerThreadPerLoop()
//{
//	unsigned int cell_num_in_one_thread;
//	unsigned int max_loop_num;
//	unsigned int* global_cell_start_index;
//	unsigned int previous_loop_end;
//	max_loop_num = (non_empty_cell_index_ave_begin_per_thread[1] - non_empty_cell_index_ave_begin_per_thread[0])/ max_pair_num_in_one_loop_per_thread;
//	if (max_loop_num > max_num_to_loop_to_find_pair) {
//		max_loop_num = max_num_to_loop_to_find_pair;
//	}
//	for (unsigned int i = 0; i < thread_num; ++i) {
//		global_cell_start_index = global_cell_start[i] + 1;
//		cell_num_in_one_thread = non_empty_cell_index_ave_begin_per_thread[i + 1] - non_empty_cell_index_ave_begin_per_thread[i];
//		//max_loop_num = cell_num_in_one_thread / max_pair_num_in_one_loop_per_thread;
//		arrangeIndex(max_loop_num, cell_num_in_one_thread, global_cell_start_index);
//		global_cell_start[i][0] = max_loop_num;
//		if (i > 0) {
//			previous_loop_end = global_cell_start[i - 1][global_cell_start[i - 1][0] + 1];
//			for (unsigned int j = 0; j <= global_cell_start[i][0]; ++j) {
//				global_cell_start_index[j] += previous_loop_end;
//			}
//		}
//	}
//	for (unsigned int i = 0; i <= thread_num; ++i) {
//		std::cout << non_empty_cell_index_ave_begin_per_thread[i] << " ";		
//	}
//	std::cout << std::endl;
//	for (unsigned int i = 0; i <= thread_num; ++i) {
//		for (unsigned int j = 0; j <= global_cell_start[i][0]; ++j) {
//			std::cout << global_cell_start[i][1 + j] << " ";
//		}
//		std::cout << std::endl;
//	}
//}


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

//OBTAIN_PAIR_COUNT
void SpatialHashing::obtainPairCount(int thread_No)
{
	unsigned int thread_num_ = thread_num;
	unsigned int* hash_cell_pair_num_ = hash_cell_pair_num_prefix_vertex_triangle + 2;
	unsigned int last_cell_index = non_empty_cell_index_begin_per_thread_vertex_triangle[thread_No + 1];
	unsigned int start = non_empty_cell_index_begin_per_thread_vertex_triangle[thread_No];
	unsigned int* non_empty_cell_index_ = non_empty_cell_index_vertex_triangle[0] + 1;
	//unsigned int pair_num_with_collider = 0;
	unsigned int ori_first_thread_num;

	unsigned int** spatial_hashing_cell_ = spatial_hashing_cell_triangle[0];
	unsigned int* spatial_hashing_cell_size_ = spatial_hashing_cell_triangle_size[0];

	unsigned int** spatial_hashing_cell_vertex_ = spatial_hashing_cell_vertex[0];
	unsigned int* spatial_hashing_cell_vertex_size_ = spatial_hashing_cell_vertex_size[0];

	unsigned int** spatial_hashing_cell_collider_ = spatial_hashing_cell_collider_triangle[0];
	unsigned int* spatial_hashing_cell_collider_size_ = spatial_hashing_cell_collider_triangle_size[0];

	unsigned int** spatial_hashing_cell_collider_vertex_ = spatial_hashing_cell_collider_vertex[0];
	unsigned int* spatial_hashing_cell_collider_vertex_size_ = spatial_hashing_cell_collider_vertex_size[0];


	bool has_collider_ = has_collider;
	unsigned int cell_index;

	unsigned int* spatial_hashing_cell_address;
	unsigned int pair_num_cell;
	for (unsigned int i = start; i < last_cell_index; ++i) {
		cell_index = non_empty_cell_index_[i];
		pair_num_cell = 0;
		ori_first_thread_num = spatial_hashing_cell_size_[cell_index];
		spatial_hashing_cell_address = spatial_hashing_cell_[cell_index] + ori_first_thread_num;
		for (unsigned int t = 1; t < thread_num_; ++t) {
			if (spatial_hashing_cell_triangle_size[t][cell_index]) {
				memcpy(spatial_hashing_cell_address, spatial_hashing_cell_triangle[t][cell_index], spatial_hashing_cell_triangle_size[t][cell_index] << 2);
				spatial_hashing_cell_address += spatial_hashing_cell_triangle_size[t][cell_index];
			}
		}
		spatial_hashing_cell_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_[cell_index];

		ori_first_thread_num = spatial_hashing_cell_vertex_size_[cell_index];
		spatial_hashing_cell_address = spatial_hashing_cell_vertex_[cell_index] + ori_first_thread_num;

		for (unsigned int t = 1; t < thread_num_; ++t) {
			if (spatial_hashing_cell_vertex_size[t][cell_index]) {
				memcpy(spatial_hashing_cell_address, spatial_hashing_cell_vertex[t][cell_index], spatial_hashing_cell_vertex_size[t][cell_index] << 2);
				spatial_hashing_cell_address += spatial_hashing_cell_vertex_size[t][cell_index];
			}
		}
		spatial_hashing_cell_vertex_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_vertex_[cell_index];

		pair_num_cell = (spatial_hashing_cell_size_[cell_index] * spatial_hashing_cell_vertex_size_[cell_index]) >> 2;
		
		if (has_collider_) {
			ori_first_thread_num = spatial_hashing_cell_collider_size_[cell_index];
			spatial_hashing_cell_address = spatial_hashing_cell_collider_[cell_index] + ori_first_thread_num;
			for (unsigned int t = 1; t < thread_num_; ++t) {
				if (spatial_hashing_cell_collider_triangle_size[t][cell_index]) {
					memcpy(spatial_hashing_cell_address, spatial_hashing_cell_collider_triangle[t][cell_index], spatial_hashing_cell_collider_triangle_size[t][cell_index] << 2);
					spatial_hashing_cell_address += spatial_hashing_cell_collider_triangle_size[t][cell_index];
				}
			}
			spatial_hashing_cell_collider_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_collider_[cell_index];

			ori_first_thread_num = spatial_hashing_cell_collider_vertex_size_[cell_index];
			spatial_hashing_cell_address = spatial_hashing_cell_collider_vertex_[cell_index] + ori_first_thread_num;
			for (unsigned int t = 1; t < thread_num_; ++t) {
				if (spatial_hashing_cell_collider_vertex_size[t][cell_index]) {
					memcpy(spatial_hashing_cell_address, spatial_hashing_cell_collider_vertex[t][cell_index], spatial_hashing_cell_collider_vertex_size[t][cell_index] << 2);
					spatial_hashing_cell_address += spatial_hashing_cell_collider_vertex_size[t][cell_index];
				}
			}
			spatial_hashing_cell_collider_vertex_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_collider_vertex_[cell_index];

			
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
	for (unsigned int i = start; i < last_cell_index; ++i) {
		cell_index = non_empty_cell_index_[i];
		pair_num_cell = 0;
		ori_first_thread_num = spatial_hashing_cell_size_[cell_index];
		spatial_hashing_cell_address = spatial_hashing_cell_[cell_index] + ori_first_thread_num;
		for (unsigned int t = 1; t < thread_num_; ++t) {
			if (spatial_hashing_cell_edge_size[t][cell_index]) {
				memcpy(spatial_hashing_cell_address, spatial_hashing_cell_edge[t][cell_index], spatial_hashing_cell_edge_size[t][cell_index] << 2);
				spatial_hashing_cell_address += spatial_hashing_cell_edge_size[t][cell_index];
			}
		}
		spatial_hashing_cell_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_[cell_index];
		if (spatial_hashing_cell_size_[cell_index] > 2) {
			pair_num_cell = (spatial_hashing_cell_size_[cell_index] * (spatial_hashing_cell_size_[cell_index] - 2)) >> 3;
		}
		if (has_collider_) {
			ori_first_thread_num = spatial_hashing_cell_collider_size_[cell_index];
			spatial_hashing_cell_address = spatial_hashing_cell_collider_[cell_index] + ori_first_thread_num;
			for (unsigned int t = 1; t < thread_num_; ++t) {
				if (spatial_hashing_cell_collider_edge_size[t][cell_index]) {
					memcpy(spatial_hashing_cell_address, spatial_hashing_cell_collider_edge[t][cell_index], spatial_hashing_cell_collider_edge_size[t][cell_index] << 2);
					spatial_hashing_cell_address += spatial_hashing_cell_collider_edge_size[t][cell_index];
				}
			}
			spatial_hashing_cell_collider_size_[cell_index] = spatial_hashing_cell_address - spatial_hashing_cell_collider_[cell_index];
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