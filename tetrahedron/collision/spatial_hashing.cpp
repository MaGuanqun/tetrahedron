#include"spatial_hashing.h"



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
	if (!for_construct_patch) {
		cell_length = 2.0 * max_length + 2.0 * tolerance_ratio[SELF_POINT_TRIANGLE] * ave_length;
	}
	else {
		cell_length = tolerance_ratio[SELF_POINT_TRIANGLE] * max_length;
	}
	//std::cout << "ave_length" << " " << ave_length << std::endl;
}

void SpatialHashing::setInObject(std::vector<Cloth>* cloth, std::vector<Collider>* collider,
	std::vector<Tetrahedron>* tetrahedron, Thread* thread, double* tolerance_ratio,
	unsigned int max_cell_count, bool for_construct_patch)
{
	has_collider = !collider->empty();

	this->for_construct_patch = for_construct_patch;
	tetrahedron_begin_obj_index = cloth->size();
	initialHashCellLength(cloth, tetrahedron, cell_length, tolerance_ratio);
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

	if (!for_construct_patch) {
		//actual_hash_count_start_per_thread = new unsigned int[thread_num + 1];
		//total_hash_count_start_per_thread = new unsigned int[thread_num + 1];
		//prefix_sum_thread_start = new unsigned int[thread_num];
		//memset(prefix_sum_thread_start, 0, 4 * thread_num);

		//initialParallePrefix();

		prefix_sum.resize(thread_num);
		prefix_sum_collider.resize(thread_num);
		hash_index_count_per_thread.resize(thread_num);
		hash_index_count_per_thread_collider.resize(thread_num);
		for (unsigned int i = 0; i < thread_num; ++i) {
			prefix_sum[i].reserve(100 * 100 * 100 / thread_num);
			if (has_collider) {
				prefix_sum_collider[i].reserve(100 * 100 * 100 / thread_num);
				hash_index_count_per_thread_collider[i].reserve(100 * 100 * 100 / thread_num);
			}
			hash_index_count_per_thread[i].reserve(100 * 100 * 100 / thread_num);
		}

		radix_sort = new RadixSort[thread_num];
		radix_sort_collider = new RadixSort[thread_num];
		for (unsigned int i = 0; i < thread_num; ++i) {
			radix_sort[i].initial(thread, true);
			radix_sort_collider[i].initial(thread, true);
		}
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
	spatial_hashing_value = new unsigned int* [thread_num];
	spatial_hashing_obj_index = new unsigned int* [thread_num];
	spatial_hashing_triangle_index = new unsigned int* [thread_num];

	spatial_hashing_value_collider = new unsigned int* [thread_num];
	spatial_hashing_obj_index_collider = new unsigned int* [thread_num];
	spatial_hashing_triangle_index_collider = new unsigned int* [thread_num];

	for (int i = 0; i < thread_num; ++i) {
		spatial_hashing_value[i] = new unsigned int[total_hash_size_per_thread[i]];
		spatial_hashing_obj_index[i] = new unsigned int[total_hash_size_per_thread[i]];
		spatial_hashing_triangle_index[i] = new unsigned int[total_hash_size_per_thread[i]];
		memset(spatial_hashing_value[i], -1, 4 * total_hash_size_per_thread[i]);

		if (has_collider && triangle_number_per_thread_collider[i] != 0) {
			spatial_hashing_value_collider[i] = new unsigned int[total_hash_size_per_thread_collider[i]];
			spatial_hashing_obj_index_collider[i] = new unsigned int[total_hash_size_per_thread_collider[i]];
			spatial_hashing_triangle_index_collider[i] = new unsigned int[total_hash_size_per_thread_collider[i]];
			memset(spatial_hashing_value_collider[i], -1, 4 * total_hash_size_per_thread_collider[i]);
		}
	}

	if (!for_construct_patch) {
		for (unsigned int i = 0; i < thread_num; ++i) {
			radix_sort[i].initialArray(total_hash_size_per_thread[i]);
		}
		if (has_collider) {
			for (unsigned int i = 0; i < thread_num; ++i) {
				if (triangle_number_per_thread_collider[i] != 0) {
					radix_sort_collider[i].initialArray(total_hash_size_per_thread_collider[i]);
				}
			}
		}

		obj_is_used0 = new bool** [thread_num];
		obj_is_used1 = new bool** [thread_num];
		for (int i = 0; i < thread_num; ++i) {
			obj_is_used0[i] = new bool* [collider_begin_obj_index];
			obj_is_used1[i] = new bool* [collider_begin_obj_index];
			for (int j = 0; j < cloth->size(); ++j) {
				obj_is_used0[i][j] = new bool[(*cloth)[j].mesh_struct.triangle_indices.size()];
				obj_is_used1[i][j] = new bool[(*cloth)[j].mesh_struct.triangle_indices.size()];
				memset(obj_is_used0[i][j], 0, (*cloth)[j].mesh_struct.triangle_indices.size());
				memset(obj_is_used1[i][j], 0, (*cloth)[j].mesh_struct.triangle_indices.size());
			}
			for (int j = 0; j < tetrahedron->size(); ++j) {
				obj_is_used0[i][j + tetrahedron_begin_obj_index] = new bool[(*tetrahedron)[j].mesh_struct.triangle_indices.size()];
				obj_is_used1[i][j + tetrahedron_begin_obj_index] = new bool[(*tetrahedron)[j].mesh_struct.triangle_indices.size()];
				memset(obj_is_used0[i][j + tetrahedron_begin_obj_index], 0, (*tetrahedron)[j].mesh_struct.triangle_indices.size());
				memset(obj_is_used1[i][j + tetrahedron_begin_obj_index], 0, (*tetrahedron)[j].mesh_struct.triangle_indices.size());
			}
			/*if (for_construct_patch) {
				for (int j = 0; j < collider->size(); ++j) {
					obj_is_used[i][j + collider_begin_obj_index] = new bool[(*collider)[j].mesh_struct.triangle_indices.size()];
					memset(obj_is_used[i][j + collider_begin_obj_index], 0, (*collider)[j].mesh_struct.triangle_indices.size());
				}
			}*/
		}

		collider_is_used0 = new bool** [thread_num];
		if (has_collider) {
			for (int i = 0; i < thread_num; ++i) {
				collider_is_used0[i] = new bool* [collider->size()];
				for (int j = 0; j < collider->size(); ++j) {
					collider_is_used0[i][j] = new bool[(*collider)[j].mesh_struct.triangle_indices.size()];
					memset(collider_is_used0[i][j], 0, (*collider)[j].mesh_struct.triangle_indices.size());
				}
			}
		}
	}
	//
	actual_hash_value_end_index_ref.resize(thread_num);
	actual_hash_value_end_index_ref_collider.resize(thread_num);
	actual_exist_cell_begin_per_thread.resize(thread_num + 1);
	actual_hash_value_count_per_thread.resize(thread_num, 0);
	actual_hash_value_count_per_thread_collider.resize(thread_num, 0);

	triangle_pair.resize(thread_num);
	triangle_pair_with_collider.resize(thread_num);

	int total_triangle_num = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		total_triangle_num += cloth->data()[i].mesh_struct.triangle_indices.size();
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		total_triangle_num += tetrahedron->data()[i].mesh_struct.triangle_indices.size();
	}
	for (unsigned int i = 0; i < thread_num; ++i) {
		triangle_pair[i].reserve(2 * total_triangle_num);
		triangle_pair_with_collider[i].reserve(2 * total_triangle_num);
	}

	obj_tri_aabb.resize(collider_begin_obj_index);
	for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
		if (i < cloth->size()) {
			obj_tri_aabb[i] = cloth->data()[i].triangle_AABB.data();
		}
		else {
			obj_tri_aabb[i] = tetrahedron->data()[i].triangle_AABB.data();
		}
	}

	if (has_collider) {
		collider_tri_aabb.resize(collider->size());
		for (unsigned int i = 0; i < collider->size(); ++i) {
			collider_tri_aabb[i] = collider->data()[i].triangle_AABB.data();
		}
	}


}



void SpatialHashing::deleteArray()
{
	delete[] spatial_hashing_value;
	delete[] spatial_hashing_obj_index;
	delete[] spatial_hashing_triangle_index;
	delete[] spatial_hashing_value_collider;
	delete[] spatial_hashing_obj_index_collider;
	delete[] spatial_hashing_triangle_index_collider;

	delete[] radix_sort;
	delete[] radix_sort_collider;
}


void SpatialHashing::initialTriangleHash()
{
	//obj_triangle_hash.resize(cloth->size()+ tetrahedron->size());
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
	getSceneAABB();
	for (unsigned int i = 0; i < 3; ++i) {
		cell_number[i] = (unsigned int)floor((scene_aabb[i + 3] - scene_aabb[i]) / cell_length) + 1;
		hash_max_index[i] = cell_number[i] - 1;
	}

	cell_num0_cell_num1 = cell_number[0] * cell_number[1];
	hash_table_size = cell_number[0] * cell_number[1] * cell_number[2];
	std::cout << "spatial hashing size " << hash_table_size << " " << cell_number[0] << " " << cell_number[1] << " " << cell_number[2] << std::endl;

}

void SpatialHashing::buildSpatialHashing()
{
	setSpatialHashing();
	thread->assignTask(this, TRIANGLE_HASHING);
	//if (!for_construct_patch) {
		//memcpy(obj_triangle_hash, spatial_hashing_value, total_hash_size * 4);	

	setPrifixSum();
	thread->assignTask(this, SH_FIND_ALL_TRIANGLE_PAIRS);
	//}
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
	unsigned int last_cell_index = actual_exist_cell_begin_per_thread[thread_No + 1];

	std::vector<unsigned int> cell_triangle_index;
	std::vector<unsigned int> cell_obj_index;
	std::vector<unsigned int> cell_collider_triangle_index;
	std::vector<unsigned int> cell_collider_index;

	cell_triangle_index.reserve(30);
	cell_obj_index.reserve(30);
	cell_collider_triangle_index.reserve(30);
	cell_collider_index.reserve(30);

	unsigned int previous_tri_num;
	unsigned int current_thread_triangle_num;

	bool** obj_is_used0_ = obj_is_used0[thread_No];
	bool** obj_is_used1_ = obj_is_used1[thread_No];

	bool** collider_is_used0_ = collider_is_used0[thread_No];

	for (unsigned int i = 0; i < collider_begin_obj_index; ++i) {
		if (i < cloth->size()) {
			memset(obj_is_used0_[i], 0, cloth->data()[i].mesh_struct.triangle_indices.size());
			memset(obj_is_used1_[i], 0, cloth->data()[i].mesh_struct.triangle_indices.size());
		}
		else {
			memset(obj_is_used0_[i], 0, tetrahedron->data()[i - cloth->size()].mesh_struct.triangle_indices.size());
			memset(obj_is_used1_[i], 0, tetrahedron->data()[i - cloth->size()].mesh_struct.triangle_indices.size());
		}
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		memset(collider_is_used0_[i], 0, collider->data()[i].mesh_struct.triangle_indices.size());
	}

	std::vector<unsigned int>* triangle_pair_;
	std::vector<unsigned int>* triangle_pair_with_collider_;

	triangle_pair_ = &triangle_pair[thread_No];
	triangle_pair_with_collider_ = &triangle_pair_with_collider[thread_No];

	triangle_pair_->clear();
	triangle_pair_with_collider_->clear();

	double* aabb_1;

	unsigned int obj_index_0, obj_index_1, triangle_index_0, triangle_index_1;

	for (unsigned int cell_index = actual_exist_cell_begin_per_thread[thread_No];
		cell_index < last_cell_index; ++cell_index) {
		cell_triangle_index.clear();
		cell_obj_index.clear();
		cell_collider_triangle_index.clear();
		cell_collider_index.clear();

		for (unsigned int t = 0; t < thread_num; ++t) {
			current_thread_triangle_num = prefix_sum[t][cell_index + 1] - prefix_sum[t][cell_index];
			if (current_thread_triangle_num > 0) {
				previous_tri_num = cell_triangle_index.size();
				cell_triangle_index.resize(previous_tri_num + current_thread_triangle_num);
				cell_obj_index.resize(previous_tri_num + current_thread_triangle_num);
				memcpy(cell_triangle_index.data() + previous_tri_num, spatial_hashing_triangle_index[t] + prefix_sum[t][cell_index], 4 * current_thread_triangle_num);
				memcpy(cell_obj_index.data() + previous_tri_num, spatial_hashing_obj_index[t] + prefix_sum[t][cell_index], 4 * current_thread_triangle_num);
			}
		}
		if (has_collider) {
			for (unsigned int t = 0; t < thread_num; ++t) {
				current_thread_triangle_num = prefix_sum_collider[t][(cell_index << 1) + 1] - prefix_sum_collider[t][cell_index << 1];
				if (current_thread_triangle_num > 0) {
					previous_tri_num = cell_collider_triangle_index.size();
					cell_collider_triangle_index.resize(previous_tri_num + current_thread_triangle_num);
					cell_collider_index.resize(previous_tri_num + current_thread_triangle_num);
					memcpy(cell_collider_triangle_index.data() + previous_tri_num,
						spatial_hashing_triangle_index_collider[t] + prefix_sum_collider[t][cell_index << 1], 4 * current_thread_triangle_num);
					memcpy(cell_collider_index.data() + previous_tri_num,
						spatial_hashing_obj_index_collider[t] + prefix_sum_collider[t][cell_index << 1], 4 * current_thread_triangle_num);
				}
			}
		}
		//there exist collision

		if (cell_triangle_index.size() > 1 || (!cell_collider_triangle_index.empty())) {
			//if (thread_No == 4) {
			//	std::cout << cell_triangle_index.size() << std::endl;
			//	for (unsigned int k = 0; k < cell_triangle_index.size(); ++k) {
			//		std::cout << cell_triangle_index[k] << " ";
			//	}
			//	std::cout << std::endl;
			//}

			for (unsigned int i = 0; i < cell_triangle_index.size(); ++i) {
				obj_index_0 = cell_obj_index[i]; triangle_index_0 = cell_triangle_index[i];
				aabb_1 = obj_tri_aabb[obj_index_0][triangle_index_0].data();
				for (unsigned int j = i + 1; j < cell_triangle_index.size(); ++j) {
					obj_index_1 = cell_obj_index[j];
					triangle_index_1 = cell_triangle_index[j];

					//if (obj_index_0 < obj_index_1 || (obj_index_0 == obj_index_1 && triangle_index_0 < triangle_index_1)) {
						//if (!(obj_is_used0_[obj_index_0][triangle_index_0] && obj_is_used1_[obj_index_1][triangle_index_1])) {
							//obj_is_used0_[obj_index_0][triangle_index_0] = true;
							//obj_is_used1_[obj_index_1][triangle_index_1] = true;

							//if (triangle_index_0 == 1) {
							//	std::cout << triangle_index_1 << " " << triangle_index_0 << std::endl;
							//}

					if (AABB::AABB_intersection(aabb_1, obj_tri_aabb[obj_index_1][triangle_index_1].data())) {

						triangle_pair_->emplace_back(triangle_index_0);
						triangle_pair_->emplace_back(obj_index_0);
						triangle_pair_->emplace_back(triangle_index_1);
						triangle_pair_->emplace_back(obj_index_1);

					}
					//}
				//}
				//else {
				//	if (!(obj_is_used0_[obj_index_1][triangle_index_1] && obj_is_used1_[obj_index_0][triangle_index_0])) {
				//		obj_is_used0_[obj_index_1][triangle_index_1] = true;
				//		obj_is_used1_[obj_index_0][triangle_index_0] = true;
				//		if (AABB::AABB_intersection(aabb_1, obj_tri_aabb[obj_index_1][triangle_index_1].data())) {
				//			triangle_pair_->emplace_back(triangle_index_1);
				//			triangle_pair_->emplace_back(obj_index_1);
				//			triangle_pair_->emplace_back(triangle_index_0);
				//			triangle_pair_->emplace_back(obj_index_0);						
				//		}
				//	}
				//}
				}
				for (unsigned int j = 0; j < cell_collider_triangle_index.size(); ++j) {
					obj_index_1 = cell_collider_index[j];
					triangle_index_1 = cell_collider_triangle_index[j];
					//if (!(obj_is_used0_[obj_index_0][triangle_index_0] && collider_is_used0_[obj_index_1][triangle_index_1])) {
					//	obj_is_used0_[obj_index_0][triangle_index_0] = true;
					//	collider_is_used0_[obj_index_1][triangle_index_1] = true;
					if (AABB::AABB_intersection(aabb_1, collider_tri_aabb[obj_index_1][triangle_index_1].data())) {
						triangle_pair_with_collider_->emplace_back(triangle_index_0);
						triangle_pair_with_collider_->emplace_back(obj_index_0);
						triangle_pair_with_collider_->emplace_back(triangle_index_1);
						triangle_pair_with_collider_->emplace_back(obj_index_1);
					}
					//}
				}
			}
		}
	}

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

	//initial prefix sum in TRIANGLE_HASHING, clear() & push_back(0)
	std::vector<unsigned int*> index_per_thread(thread_num);
	std::vector<unsigned int*> index_per_thread_collider(thread_num);
	std::vector<unsigned int> start_index(thread_num, 0);
	std::vector<unsigned int> start_index_collider(thread_num, 0);
	std::vector<unsigned int> record_collider_index_start(thread_num, 0);
	//std::vector<unsigned int> record_collider_index_end(thread_num, 0);

	for (unsigned int i = 0; i < thread_num; ++i) {
		index_per_thread[i] = spatial_hashing_value[i];
	}
	if (has_collider) {
		for (unsigned int i = 0; i < thread_num; ++i) {
			index_per_thread_collider[i] = spatial_hashing_value_collider[i];
		}
	}
	unsigned int current_hash_index = UINT_MAX;
	bool end_loop;
	while (true)
	{
		end_loop = (index_per_thread[0] == actual_hash_value_end_index_ref[0]);
		for (unsigned int i = 1; i < thread_num; ++i) {
			end_loop = end_loop && (index_per_thread[i] == actual_hash_value_end_index_ref[i]);
		}
		if (end_loop) {
			break;
		}
		current_hash_index = UINT_MAX;
		for (unsigned int i = 0; i < thread_num; ++i) {
			if (index_per_thread[i] != actual_hash_value_end_index_ref[i]) {
				if (current_hash_index > *(index_per_thread[i])) {
					current_hash_index = *(index_per_thread[i]);
				}
			}

		}
		for (unsigned int i = 0; i < thread_num; ++i) {
			prefix_sum[i].emplace_back(prefix_sum[i].back());
			if (index_per_thread[i] != actual_hash_value_end_index_ref[i]) {
				if (*(index_per_thread[i]) == current_hash_index)
				{
					prefix_sum[i].back() += hash_index_count_per_thread[i][start_index[i]];
					index_per_thread[i] += hash_index_count_per_thread[i][start_index[i]];
					start_index[i]++;
				}
			}
		}

		if (has_collider)
		{
			for (unsigned int i = 0; i < thread_num; ++i) {
				if (start_index_collider[i] < hash_index_count_per_thread_collider[i].size()) {
					while (*(index_per_thread_collider[i]) < current_hash_index)
					{
						index_per_thread_collider[i] += hash_index_count_per_thread_collider[i][start_index_collider[i]];
						record_collider_index_start[i] += hash_index_count_per_thread_collider[i][start_index_collider[i]];
						start_index_collider[i]++;
						if (start_index_collider[i] == hash_index_count_per_thread_collider[i].size()) {
							break;
						}
					}
				}
				prefix_sum_collider[i].emplace_back(record_collider_index_start[i]);

				if (start_index_collider[i] < hash_index_count_per_thread_collider[i].size()) {
					if (*(index_per_thread_collider[i]) == current_hash_index)
					{
						index_per_thread_collider[i] += hash_index_count_per_thread_collider[i][start_index_collider[i]];
						record_collider_index_start[i] += hash_index_count_per_thread_collider[i][start_index_collider[i]];
						start_index_collider[i]++;
					}
				}
				prefix_sum_collider[i].emplace_back(record_collider_index_start[i]);
			}

		}

	}




	arrangeIndex(thread_num, prefix_sum[0].size() - 1, actual_exist_cell_begin_per_thread.data());
	//test();
	//if (prefix_sum[0].size() != prefix_sum_collider[0].size() / 2 + 1) {
	//	std::cout << prefix_sum[0].size() << " " << prefix_sum_collider[0].size() << std::endl;
	//	std::cout << "prefix sum not consistent with prefix sum collider " << std::endl;
	//}

	//testPrifixSum1();
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



//void SpatialHashing::searchTriangle(double* aabb, unsigned int input_obj_No, unsigned int triangle_index, 
//	std::vector<unsigned int>* obj_neighbor_index, bool is_collider, unsigned int thread_No)
//{
//	std::vector<unsigned int>cloth_index_record;	std::vector<unsigned int>triangle_index_record;
//	cloth_index_record.reserve(20);
//	triangle_index_record.reserve(20);
//	bool** obj_is_used_ = obj_is_used[thread_No];
//	unsigned int hash_cell_index_end;
//	unsigned int obj_No; unsigned int current_triangle_index;
//	unsigned int obj_num = cloth->size() + tetrahedron->size();
//	for (unsigned int i = 0; i < obj_num; ++i) {
//		obj_neighbor_index[i].clear();
//		//obj_neighbor_index[i].reserve(10);
//	}
//	if (is_collider) {
//		std::vector<unsigned int>hash_value;
//		obtainTriangleHashingValue(aabb, &hash_value);
//		for (unsigned int i = 0; i < hash_value.size(); ++i) {
//			hash_cell_index_end = prefix_sum[hash_value[i] + 1];
//			for (unsigned int j = prefix_sum[hash_value[i]]; j < hash_cell_index_end; ++j) {
//				obj_No = spatial_hashing_obj_index[j];
//				current_triangle_index = spatial_hashing_triangle_index[j];
//				if (!obj_is_used_[obj_No][current_triangle_index]) {
//					obj_is_used_[obj_No][current_triangle_index] = true;					
//					cloth_index_record.push_back(obj_No);
//					triangle_index_record.push_back(current_triangle_index);
//					if (obj_No < tetrahedron_begin_obj_index) {
//						if (AABB::AABB_intersection(aabb, (*cloth)[obj_No].triangle_AABB[current_triangle_index].data())) {
//							obj_neighbor_index[obj_No].push_back(current_triangle_index);
//							//std::cout << "around " << obj_No << " " << current_triangle_index << std::endl;
//						}
//					}
//					else {
//						if (AABB::AABB_intersection(aabb, (*tetrahedron)[obj_No - tetrahedron_begin_obj_index].triangle_AABB[current_triangle_index].data())) {
//							obj_neighbor_index[obj_No].push_back(current_triangle_index);
//						}
//					}
//				}
//			}
//			//if (input_obj_No == 0 && triangle_index == 0) {
//			//	std::cout <<"prefix "<< prefix_sum[hash_value[i]] << " " << hash_cell_index_end << std::endl;
//			//}
//		}
//		//if (thread_No == 0) {
//		//	if (!triangle_index_record.empty()) {
//		//		std::cout << "in func " << obj_neighbor_index[0].size() << std::endl;
//		//	}
//		//}
//	}
//	else {
//		obj_is_used_[input_obj_No][triangle_index] = true;
//		cloth_index_record.push_back(input_obj_No);
//		triangle_index_record.push_back(triangle_index);
//
//		unsigned int* hash_value;
//		hash_value = obj_triangle_hash + max_cell_count*(triangle_begin_per_obj[input_obj_No]+triangle_index);
//		for (unsigned int i = 0; i < max_cell_count; ++i) {			
//			if (*hash_value == UINT_MAX) {
//				break;
//			}
//			hash_cell_index_end = prefix_sum[*hash_value + 1];
//
//			for (unsigned int j = prefix_sum[*hash_value]; j < hash_cell_index_end; ++j) {
//				
//				obj_No = spatial_hashing_obj_index[j];
//				//std::cout << "prefix_sum "<< *hash_value<<" "<< triangle_index<<" "<<i << " " << j << " " << obj_No << std::endl;
//				current_triangle_index = spatial_hashing_triangle_index[j];
//
//				if (!obj_is_used_[obj_No][current_triangle_index]) {
//					obj_is_used_[obj_No][current_triangle_index] = true;
//					cloth_index_record.push_back(obj_No);
//					triangle_index_record.push_back(current_triangle_index);
//
//					if (obj_No < tetrahedron_begin_obj_index) {
//						if (AABB::AABB_intersection(aabb, (*cloth)[obj_No].triangle_AABB[current_triangle_index].data())) {
//							obj_neighbor_index[obj_No].push_back(current_triangle_index);
//						}
//					}
//					else {
//						if (AABB::AABB_intersection(aabb, 
//							(*tetrahedron)[obj_No - tetrahedron_begin_obj_index].triangle_AABB[current_triangle_index].data())) {
//							obj_neighbor_index[obj_No].push_back(current_triangle_index);
//						}
//					}
//				}
//			}
//			hash_value++;		
//		}
//	}
//	
//	
//	for (unsigned int i = 0; i < triangle_index_record.size(); ++i) {
//		obj_is_used_[cloth_index_record[i]][triangle_index_record[i]] = false;
//	}
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
	unsigned int* triangle_begin;
	std::array<double, 6>* aabb;
	int vector_size;
	unsigned int triangle_index_total;
	unsigned int triangle_end;
	unsigned int largest_count_in_hash_value_list;

	unsigned int* spatial_hashing_triangle_index_ = spatial_hashing_triangle_index[thread_No];
	unsigned int* spatial_hashing_value_ = spatial_hashing_value[thread_No];
	unsigned int* spatial_hashing_obj_index_ = spatial_hashing_obj_index[thread_No];


	memset(spatial_hashing_value_, -1, 4 * actual_hash_value_count_per_thread[thread_No]);

	for (unsigned int i = 0; i < cloth->size(); ++i) {
		aabb = (*cloth)[i].triangle_AABB.data();
		triangle_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
		triangle_end = triangle_begin[thread_No + 1];
		for (unsigned int j = triangle_begin[thread_No]; j < triangle_end; ++j) {
			triangle_index_total = max_cell_count * (triangle_begin_per_obj_per_thread[thread_No][i] + j);
			triangleHashValue(aabb[j].data(), spatial_hashing_triangle_index_ + triangle_index_total,
				spatial_hashing_value_ + triangle_index_total, j, spatial_hashing_obj_index_ + triangle_index_total, i);
		}
	}

	unsigned int obj_No;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		obj_No = i + tetrahedron_begin_obj_index;
		aabb = (*tetrahedron)[i].triangle_AABB.data();
		triangle_begin = (*tetrahedron)[i].mesh_struct.face_index_begin_per_thread.data();
		triangle_end = triangle_begin[thread_No + 1];
		for (unsigned int j = triangle_begin[thread_No]; j < triangle_end; ++j) {
			triangle_index_total = max_cell_count * (triangle_begin_per_obj_per_thread[thread_No][obj_No] + j);
			triangleHashValue(aabb[j].data(), spatial_hashing_triangle_index_ + triangle_index_total,
				spatial_hashing_value_ + triangle_index_total, j, spatial_hashing_obj_index_ + triangle_index_total, obj_No);
		}
	}

	radix_sort[thread_No].radixSort(hash_table_size, spatial_hashing_value_, spatial_hashing_triangle_index_, spatial_hashing_obj_index_,
		largest_count_in_hash_value_list);



	testIfRadixSortIsRight(spatial_hashing_value_, total_hash_size_per_thread[thread_No]);
	actual_hash_value_count_per_thread[thread_No] = total_hash_size_per_thread[thread_No];
	actual_hash_value_end_index_ref[thread_No] = spatial_hashing_value[thread_No] + total_hash_size_per_thread[thread_No];
	for (unsigned int i = total_hash_size_per_thread[thread_No] - largest_count_in_hash_value_list;
		i < total_hash_size_per_thread[thread_No]; ++i) {
		if (spatial_hashing_value_[i] == UINT_MAX) {
			actual_hash_value_end_index_ref[thread_No] = spatial_hashing_value[thread_No] + i;
			actual_hash_value_count_per_thread[thread_No] = i;
			break;
		}
	}
	prefix_sum[thread_No].clear();
	prefix_sum[thread_No].push_back(0);
	countHashIndexPerThread(thread_No, actual_hash_value_count_per_thread[thread_No], &hash_index_count_per_thread[thread_No],
		spatial_hashing_value_);

	if (has_collider) {
		unsigned int* spatial_hashing_triangle_index_collider_ = spatial_hashing_triangle_index_collider[thread_No];
		unsigned int* spatial_hashing_value_collider_ = spatial_hashing_value_collider[thread_No];
		unsigned int* spatial_hashing_obj_index_collider_ = spatial_hashing_obj_index_collider[thread_No];
		memset(spatial_hashing_value_collider_, -1, 4 * actual_hash_value_count_per_thread_collider[thread_No]);
		for (unsigned int i = 0; i < collider->size(); ++i) {
			aabb = (*collider)[i].triangle_AABB.data();
			triangle_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread.data();
			triangle_end = triangle_begin[thread_No + 1];
			for (unsigned int j = triangle_begin[thread_No]; j < triangle_end; ++j) {
				triangle_index_total = max_cell_count * (triangle_begin_per_obj_per_thread_collider[thread_No][i] + j);
				triangleHashValue(aabb[j].data(), spatial_hashing_triangle_index_collider_ + triangle_index_total,
					spatial_hashing_value_collider_ + triangle_index_total, j, spatial_hashing_obj_index_collider_ + triangle_index_total, i);
			}
		}
		radix_sort_collider[thread_No].radixSort(hash_table_size, spatial_hashing_value_collider_, spatial_hashing_triangle_index_collider_, spatial_hashing_obj_index_collider_,
			largest_count_in_hash_value_list);
		testIfRadixSortIsRight(spatial_hashing_value_collider_, total_hash_size_per_thread_collider[thread_No]);

		actual_hash_value_count_per_thread_collider[thread_No] = total_hash_size_per_thread_collider[thread_No];
		actual_hash_value_end_index_ref_collider[thread_No] = spatial_hashing_value_collider[thread_No] + total_hash_size_per_thread_collider[thread_No];
		for (unsigned int i = total_hash_size_per_thread_collider[thread_No] - largest_count_in_hash_value_list;
			i < total_hash_size_per_thread_collider[thread_No]; ++i) {
			if (spatial_hashing_value_collider_[i] == UINT_MAX) {
				actual_hash_value_end_index_ref_collider[thread_No] = spatial_hashing_value_collider[thread_No] + i;
				actual_hash_value_count_per_thread_collider[thread_No] = i;
				break;
			}
		}
		prefix_sum_collider[thread_No].clear();
		countHashIndexPerThread(thread_No, actual_hash_value_count_per_thread_collider[thread_No],
			&hash_index_count_per_thread_collider[thread_No], spatial_hashing_value_collider_);

	}

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
		std::cout << actual_exist_cell_begin_per_thread[i] << std::endl;
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
	unsigned int* spatial_hashing_value, unsigned int triangle_index, unsigned int* spatial_hashing_obj_index, unsigned int obj_index)
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
				*(spatial_hashing_obj_index++) = obj_index;
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