#include"spatial_hashing.h"



void SpatialHashing::initialHashCellLength(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, double& cell_length, double* tolerance_ratio)
{
	double max_length = 0;
	double ave_length = 0;
	int edge_num = 0;
	std::vector<MeshStruct::Edge>* edge;
	for (int i = 0; i < cloth->size(); ++i) {
		edge = &cloth->data()[i].mesh_struct.edges;
		for(int j=0;j< edge->size();++j){
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
	cell_length = max_length + 2.0 *tolerance_ratio[SELF_POINT_TRIANGLE]* ave_length;
	//std::cout << "ave_length" << " " << ave_length << std::endl;
}

void SpatialHashing::setInObject(std::vector<Cloth>* cloth, std::vector<Collider>* collider, 
	std::vector<Tetrahedron>* tetrahedron, Thread* thread, double* tolerance_ratio)
{
	tetrahedron_begin_obj_index = cloth->size();
	initialHashCellLength(cloth, tetrahedron, cell_length, tolerance_ratio);
	this->cloth = cloth;
	this->collider = collider;
	this->tetrahedron = tetrahedron;
	this->thread = thread;
	thread_num = thread->thread_num;
	scene_aabb_max_thread.resize(thread_num);
	scene_aabb_min_thread.resize(thread_num);

	//spatial_hashing_value_per_thread.resize(thread_num);
	//spatial_hashing_obj_index_per_thread.resize(thread_num);
	//spatial_hashing_triangle_index_per_thread.resize(thread_num);
	//spatial_hashing_triangle_per_thread.resize(thread_num);

	triangle_begin_per_obj.resize(cloth->size()+ tetrahedron->size() +1);

	int total_triangle_num = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		triangle_begin_per_obj[i] = total_triangle_num;
		total_triangle_num += cloth->data()[i].mesh_struct.triangle_indices.size();
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		triangle_begin_per_obj[i+tetrahedron_begin_obj_index] = total_triangle_num;
		total_triangle_num += tetrahedron->data()[i].mesh_struct.triangle_indices.size();
	}
	triangle_begin_per_obj[triangle_begin_per_obj.size()-1] = total_triangle_num;

	max_cell_count = 27;

	total_hash_size = max_cell_count * total_triangle_num;
	spatial_hashing_value = new unsigned int[total_hash_size];
	spatial_hashing_obj_index = new unsigned int[total_hash_size];
	spatial_hashing_triangle_index= new unsigned int[total_hash_size];

	memset(spatial_hashing_value, -1, 4 * total_hash_size);
	//std::cout << "spatial_hashing_value"<< spatial_hashing_value[0] << std::endl;
	//for (int i = 0; i < thread_num; ++i) {
	//	spatial_hashing_triangle_per_thread[i].reserve(64 * total_cloth_num / thread_num);
	//}
	obj_triangle_hash = new unsigned int[total_hash_size];
	//initialTriangleHash();
	hash_value_begin.resize(thread_num,0);
	radix_sort.initial(thread);
	//testRadixSort();
	radix_sort.initialArray(total_hash_size);
	obj_is_used = new bool** [thread_num];
	for (int i = 0; i < thread_num; ++i) {
		obj_is_used[i] = new bool* [cloth->size()+tetrahedron->size()];
		for (int j = 0; j < cloth->size(); ++j) {
			obj_is_used[i][j] = new bool[(*cloth)[j].mesh_struct.triangle_indices.size()];
			memset(obj_is_used[i][j], 0, (*cloth)[j].mesh_struct.triangle_indices.size());
		}
		for (int j = 0; j < tetrahedron->size(); ++j) {
			obj_is_used[i][j+ tetrahedron_begin_obj_index] = new bool[(*tetrahedron)[j].mesh_struct.triangle_indices.size()];
			memset(obj_is_used[i][j+ tetrahedron_begin_obj_index], 0, (*tetrahedron)[j].mesh_struct.triangle_indices.size());
		}
	}
	//
	actual_hash_value_count = 1;

	actual_hash_count_start_per_thread = new unsigned int[thread_num + 1];
	total_hash_count_start_per_thread = new unsigned int[thread_num + 1];
	total_hash_count_start_per_thread_move_1 = new unsigned int[thread_num + 1];
	prefix_sum_thread_start = new unsigned int[thread_num];
	hash_value_count_start_thread = new unsigned int[thread_num];
	memset(prefix_sum_thread_start,0,4*thread_num);
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
	
	radix_sort.initialArray(vec_size);
	time_t t = clock();
	unsigned int count;
	for (unsigned int j = 0; j < 100; ++j) {
		unsigned int size = 100 * 100 * 100;
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
		//for (int i = 0; i < vec_size; ++i) {
		//	std::cout << cell[i][1] << " ";
		//}
		//std::cout << std::endl;
		radix_sort.radixSort(size, value.data(), triangle_index.data(), hash_cloth_No.data(), value.size(), count);
		//for (unsigned int i = 0; i < value.size() - 1; ++i) {
		//	if (value[i] > value[i + 1]) {
		//		std::cout << "error"<<value[i]<<" "<<value[i+1]<< std::endl;
		//	}
		//}
	}
	std::cout << "time " << clock() - t << " " << clock() - t - t1 << std::endl;
	radix_sort.deleteArray();
	//for (int i = 0; i < vec_size; ++i) {
	//	std::cout << cell[i][1] << " ";
	//}
	//std::cout << std::endl;
	//}
	
}

void SpatialHashing::setSpatialHashing()
{
	getSceneAABB();
	for (unsigned int i = 0; i < 3; ++i) {
		cell_number[i] = (unsigned int)floor((scene_aabb.max[i] - scene_aabb.min[i]) / cell_length) + 1;
		hash_max_index[i] = cell_number[i] - 1;
	}
	cell_num0_cell_num1 = cell_number[0] * cell_number[1];
	hash_table_size = cell_number[0] * cell_number[1] * cell_number[2];

	memset(spatial_hashing_value, -1, 4 * actual_hash_value_count);
}

void SpatialHashing::buildSpatialHashing()
{
	setSpatialHashing();
	thread->assignTask(this, TRIANGLE_HASHING);
	//setHashTogether();

	memcpy(obj_triangle_hash, spatial_hashing_value, total_hash_size * 4);

	radix_sort.radixSort(hash_table_size, spatial_hashing_value, spatial_hashing_triangle_index, spatial_hashing_obj_index,
		total_hash_size, largest_count_in_hash_value_list);
	setPrifixSum();
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
	//}
	actual_hash_value_count = total_hash_size - largest_count_in_hash_value_list;
	for (unsigned int i = actual_hash_value_count; i < total_hash_size; ++i) {
		if (spatial_hashing_value[i] == UINT_MAX) {
			actual_hash_value_count = i;
			break;
		}
	}

	

	prifix_sum.resize(cell_number[0] * cell_number[1] * cell_number[2] + 1);
	memset(prifix_sum.data(), 0, 4 * prifix_sum.size());

	arrangeIndex(thread_num, actual_hash_value_count, actual_hash_count_start_per_thread);
	thread->assignTask(this, PREPARE_FOR_ACTUAL_HASH_VALUE_COUNT_THREAD);

	thread->assignTask(this, ADD_COUNT_FOR_PRIFIX_SUM);

	//unsigned int* value_address;
	//unsigned int k;
	//for (unsigned int i = 1; i < thread_num; ++i) {
	//	k = 1;
	//	value_address = spatial_hashing_value + actual_hash_count_start_per_thread[i];
	//	while (*(value_address - k) == *(value_address)) {
	//		k++;
	//	}
	//	if (k < actual_hash_count_start_per_thread[i]) {
	//		actual_hash_count_start_per_thread[i] -= k - 1;
	//	}
	//	else {
	//		actual_hash_count_start_per_thread[i] = 0;
	//	}
	//}

	arrangeIndex(thread_num, prifix_sum.size(), total_hash_count_start_per_thread);
	for (int i = 0; i < thread_num; ++i) {
		total_hash_count_start_per_thread_move_1[i] = total_hash_count_start_per_thread[i] + 1;
	}
	total_hash_count_start_per_thread_move_1[thread_num] = total_hash_count_start_per_thread[thread_num];

	thread->assignTask(this, PREFIX_SUM_THREAD_1);
	hash_value_count_start_thread[0] = prifix_sum[0];
	for (int i = 1; i < thread_num; ++i) {
		hash_value_count_start_thread[i] = prifix_sum[total_hash_count_start_per_thread[i]];
		prefix_sum_thread_start[i] = prifix_sum[total_hash_count_start_per_thread[i] - 1] + prefix_sum_thread_start[i - 1];
	}

	thread->assignTask(this, PREFIX_SUM_THREAD_2);

	//for (int i = 0; i < actual_hash_value_count; ++i) {
	//	prifix_sum[spatial_hashing_value[i]]++;
	//}
	//int s = 0; int t;
	//for (int i = 0; i < prifix_sum.size(); ++i) {
	//	t = s + prifix_sum[i];
	//	prifix_sum[i] = s;
	//	s = t;
	//}
	//prifix_sum[prifix_sum.size() - 1] = t;
}

// PREFIX_SUM_THREAD_1
void SpatialHashing::prifixSum2(int thread_No)
{
	unsigned int end = total_hash_count_start_per_thread[thread_No + 1];
	for (unsigned int i = total_hash_count_start_per_thread[thread_No]+1; i < end; ++i) {
		prifix_sum[i] += prifix_sum[i - 1];
	}
}

// PREFIX_SUM_THREAD_2
void SpatialHashing::prifixSum3(int thread_No)
{
	unsigned int start_index = prefix_sum_thread_start[thread_No];
	unsigned int end = total_hash_count_start_per_thread_move_1[thread_No + 1];
	prifix_sum[total_hash_count_start_per_thread_move_1[thread_No]] =
		start_index + hash_value_count_start_thread[thread_No];
	for (unsigned int i = total_hash_count_start_per_thread_move_1[thread_No]+1; i < end; ++i) {
		prifix_sum[i] = start_index + prifix_sum[i - 1];
	}
}



//PREPARE_FOR_ACTUAL_HASH_VALUE_COUNT_THREAD
void SpatialHashing::prepareForActualHashValueCountThread(int thread_No)
{
	if (thread_No == 0) {
		return;
	}
	unsigned int* value_address;
	unsigned int k=1;
	value_address = spatial_hashing_value + actual_hash_count_start_per_thread[thread_No];
	while (*(value_address - k) == *(value_address)) {
		k++;
	}
	if (k < actual_hash_count_start_per_thread[thread_No]) {
		actual_hash_count_start_per_thread[thread_No] -= k - 1;
	}
	else {
		actual_hash_count_start_per_thread[thread_No] = 0;
	}
}

// ADD_COUNT_FOR_PRIFIX_SUM
void SpatialHashing::prifixSum1(int thread_No)
{
	unsigned int end = actual_hash_count_start_per_thread[thread_No + 1];
	for (unsigned int i = actual_hash_count_start_per_thread[thread_No]; i < end; ++i) {
		prifix_sum[spatial_hashing_value[i]]++;
	}
}




void SpatialHashing::searchTriangle(AABB& aabb, unsigned int input_obj_No, unsigned int triangle_index, 
	std::vector<int>* obj_neighbor_index, bool is_collider, unsigned int thread_No)
{
	std::vector<unsigned int>cloth_index_record;	std::vector<unsigned int>triangle_index_record;
	cloth_index_record.reserve(20);
	triangle_index_record.reserve(20);
	bool** obj_is_used_ = obj_is_used[thread_No];
	unsigned int hash_cell_index_end;
	unsigned int obj_No; unsigned int current_triangle_index;
	unsigned int obj_num = cloth->size() + tetrahedron->size();
	for (unsigned int i = 0; i < obj_num; ++i) {
		obj_neighbor_index[i].clear();
		//obj_neighbor_index[i].reserve(10);
	}
	if (is_collider) {
		std::vector<unsigned int>hash_value;
		obtainTriangleHashingValue(aabb, &hash_value);
		for (unsigned int i = 0; i < hash_value.size(); ++i) {
			hash_cell_index_end = prifix_sum[hash_value[i] + 1];
			for (unsigned int j = prifix_sum[hash_value[i]]; j < hash_cell_index_end; ++j) {
				obj_No = spatial_hashing_obj_index[j];
				current_triangle_index = spatial_hashing_triangle_index[j];
				if (!obj_is_used_[obj_No][current_triangle_index]) {
					obj_is_used_[obj_No][current_triangle_index] = true;					
					cloth_index_record.push_back(obj_No);
					triangle_index_record.push_back(current_triangle_index);
					if (obj_No < tetrahedron_begin_obj_index) {
						if (aabb.AABB_intersection((*cloth)[obj_No].triangle_AABB[current_triangle_index])) {
							obj_neighbor_index[obj_No].push_back(current_triangle_index);
							//std::cout << "around " << obj_No << " " << current_triangle_index << std::endl;
						}
					}
					else {
						if (aabb.AABB_intersection((*tetrahedron)[obj_No - tetrahedron_begin_obj_index].triangle_AABB[current_triangle_index])) {
							obj_neighbor_index[obj_No].push_back(current_triangle_index);
						}
					}
				}
			}
			//if (input_obj_No == 0 && triangle_index == 0) {
			//	std::cout <<"prefix "<< prifix_sum[hash_value[i]] << " " << hash_cell_index_end << std::endl;
			//}
		}
		//if (thread_No == 0) {
		//	if (!triangle_index_record.empty()) {
		//		std::cout << "in func " << obj_neighbor_index[0].size() << std::endl;
		//	}
		//}
	}
	else {
		obj_is_used_[input_obj_No][triangle_index] = true;
		cloth_index_record.push_back(input_obj_No);
		triangle_index_record.push_back(triangle_index);

		unsigned int* hash_value;
		hash_value = obj_triangle_hash + max_cell_count*(triangle_begin_per_obj[input_obj_No]+triangle_index);
		for (unsigned int i = 0; i < max_cell_count; ++i) {			
			if (*hash_value == UINT_MAX) {
				break;
			}
			hash_cell_index_end = prifix_sum[*hash_value + 1];

			for (unsigned int j = prifix_sum[*hash_value]; j < hash_cell_index_end; ++j) {
				
				obj_No = spatial_hashing_obj_index[j];
				//std::cout << "prefix_sum "<< *hash_value<<" "<< triangle_index<<" "<<i << " " << j << " " << obj_No << std::endl;
				current_triangle_index = spatial_hashing_triangle_index[j];

				if (!obj_is_used_[obj_No][current_triangle_index]) {
					obj_is_used_[obj_No][current_triangle_index] = true;
					cloth_index_record.push_back(obj_No);
					triangle_index_record.push_back(current_triangle_index);

					if (obj_No < tetrahedron_begin_obj_index) {
						if (aabb.AABB_intersection((*cloth)[obj_No].triangle_AABB[current_triangle_index])) {
							obj_neighbor_index[obj_No].push_back(current_triangle_index);
						}
					}
					else {
						if (aabb.AABB_intersection((*tetrahedron)[obj_No - tetrahedron_begin_obj_index].triangle_AABB[current_triangle_index])) {
							obj_neighbor_index[obj_No].push_back(current_triangle_index);
						}
					}
				}
			}
			hash_value++;		
		}
	}
	
	
	for (unsigned int i = 0; i < triangle_index_record.size(); ++i) {
		obj_is_used_[cloth_index_record[i]][triangle_index_record[i]] = false;
	}
}

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


//TRIANGLE_HASHING
void SpatialHashing::triangleHashing(int thread_No)
{
	unsigned int* triangle_begin;
	AABB* aabb;
	//std::vector<unsigned int>* spatial_hashing_value_;
	//std::vector<unsigned int>* spatial_hashing_obj_;
	//std::vector<unsigned int>* spatial_hashing_triangle_;
	//std::vector<unsigned int>* hash_value;
	//spatial_hashing_obj_ = &spatial_hashing_obj_index_per_thread[thread_No];
	//spatial_hashing_triangle_ = &spatial_hashing_triangle_index_per_thread[thread_No];
	//spatial_hashing_value_ = &spatial_hashing_value_per_thread[thread_No];
	//spatial_hashing_obj_->clear();
	//spatial_hashing_triangle_->clear();
	//spatial_hashing_value_->clear();
	int vector_size;
	unsigned int triangle_index_total;
	unsigned int triangle_end;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		aabb = (*cloth)[i].triangle_AABB.data();		
		triangle_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
		//hash_value = obj_triangle_hash[i].data();
		//vector_size = spatial_hashing_triangle_->size();
		triangle_end = triangle_begin[thread_No + 1];
		for (unsigned int j = triangle_begin[thread_No]; j < triangle_end; ++j) {
			triangle_index_total = max_cell_count * (triangle_begin_per_obj[i] + j);
			triangleHashValue(aabb[j], spatial_hashing_triangle_index+ triangle_index_total,
				spatial_hashing_value+ triangle_index_total, j, spatial_hashing_obj_index + triangle_index_total,i);
		}
	}

	unsigned int obj_No;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		obj_No = i + tetrahedron_begin_obj_index;
		aabb = (*tetrahedron)[i].triangle_AABB.data();
		triangle_begin = (*tetrahedron)[i].mesh_struct.face_index_begin_per_thread.data();
		//hash_value = obj_triangle_hash[i+tetrahedron_begin_obj_index].data();
		//vector_size = spatial_hashing_triangle_->size();
		triangle_end = triangle_begin[thread_No + 1];
		for (unsigned int j = triangle_begin[thread_No]; j < triangle_end; ++j) {
			triangle_index_total = max_cell_count * (triangle_begin_per_obj[obj_No] + j);
			triangleHashValue(aabb[j], spatial_hashing_triangle_index + triangle_index_total,
				spatial_hashing_value + triangle_index_total, j, spatial_hashing_obj_index + triangle_index_total, obj_No);
		}
	}
}

void SpatialHashing::obtainTriangleHashingValue(AABB& aabb, std::vector<unsigned int>* hash_value)
{
	unsigned int min_index[3];
	unsigned int max_index[3];
	unsigned int index_z_multi;
	unsigned int index_y_multi;
	for (unsigned int j = 0; j < 3; ++j) {
		min_index[j] = (unsigned int)floor((aabb.min[j] - scene_aabb.min[j]) / cell_length);
		max_index[j] = (unsigned int)floor((aabb.max[j] - scene_aabb.min[j]) / cell_length) + 1;
	}
	unsigned int size = (max_index[0] - min_index[0]) * (max_index[1] - min_index[1]) * (max_index[2] - min_index[2]);
	//hash_value.clear();
	hash_value->resize(size);
	unsigned int* value = hash_value->data();
	for (unsigned int index_y = min_index[1]; index_y < max_index[1]; ++index_y) {
		index_y_multi = index_y * cell_number[0];
		for (unsigned int index_z = min_index[2]; index_z < max_index[2]; ++index_z) {
			index_z_multi = index_z * cell_num0_cell_num1;
			for (unsigned int index_x = min_index[0]; index_x < max_index[0]; ++index_x) {
				*(value++) = index_x + index_y_multi + index_z_multi;
			}
		}
	}
}

void SpatialHashing::triangleHashingValue(AABB& aabb, std::vector<unsigned int>* spatial_hashing_triangle_index, 
	std::vector<unsigned int>* spatial_hashing_value, unsigned int triangle_index, std::vector<unsigned int>* hash_value)
{
	unsigned int min_index[3];
	unsigned int max_index[3];
	unsigned int index_z_multi;
	unsigned int index_y_multi;
	for (unsigned int j = 0; j < 3; ++j) {
		min_index[j] = (unsigned int)floor((aabb.min[j]- scene_aabb.min[j]) / cell_length);
		max_index[j] = (unsigned int)floor((aabb.max[j] - scene_aabb.min[j]) / cell_length) + 1;
	}
	//int size = (max_index[0] - min_index[0]) * (max_index[1] - min_index[1]) * (max_index[2] - min_index[2]);
	hash_value->clear();
	unsigned int value;
	for (unsigned int index_y = min_index[1]; index_y < max_index[1]; ++index_y) {
		index_y_multi = index_y * cell_number[0];
		for (unsigned int index_z = min_index[2]; index_z < max_index[2]; ++index_z) {
			index_z_multi = index_z * cell_num0_cell_num1;
			for (unsigned int index_x = min_index[0]; index_x < max_index[0]; ++index_x) {
				value = index_x + index_y_multi + index_z_multi;
				hash_value->push_back(value);
			}
		}
	}	
	spatial_hashing_value->insert(spatial_hashing_value->end(), hash_value->begin(),hash_value->end());
	spatial_hashing_triangle_index->insert(spatial_hashing_triangle_index->end(), hash_value->size(), triangle_index);
}

void SpatialHashing::triangleHashValue(AABB& aabb, unsigned int* spatial_hashing_triangle_index,
	unsigned int* spatial_hashing_value, unsigned int triangle_index, unsigned int* spatial_hashing_obj_index, unsigned int obj_index)
{
	unsigned int min_index[3];
	unsigned int max_index[3];
	unsigned int index_z_multi;
	unsigned int index_y_multi;
	for (unsigned int j = 0; j < 3; ++j) {
		min_index[j] = (unsigned int)floor((aabb.min[j] - scene_aabb.min[j]) / cell_length);
		max_index[j] = (unsigned int)floor((aabb.max[j] - scene_aabb.min[j]) / cell_length) + 1;
	}
	for (unsigned int index_y = min_index[1]; index_y < max_index[1]; ++index_y) {
		index_y_multi = index_y * cell_number[0];
		for (unsigned int index_z = min_index[2]; index_z < max_index[2]; ++index_z) {
			index_z_multi = index_z * cell_num0_cell_num1;
			for (unsigned int index_x = min_index[0]; index_x < max_index[0]; ++index_x) {
				*(spatial_hashing_value++)= index_x + index_y_multi + index_z_multi;
				*(spatial_hashing_triangle_index++) = triangle_index;
				*(spatial_hashing_obj_index++) = obj_index;
			}
		}
	}

}

void SpatialHashing::getSceneAABB()
{
	thread->assignTask(this, SCENE_AABB);
	double* max = scene_aabb.max;
	double* min = scene_aabb.min;
	memcpy(max, scene_aabb_max_thread[0].data(), 24); //set double to -5.31401e+303
	memcpy(min, scene_aabb_min_thread[0].data(), 24); //set double to 1.38242e+306
	for (int i = 1; i < thread_num; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (max[j] < scene_aabb_max_thread[i][j]) {
				max[j] = scene_aabb_max_thread[i][j];
			}
			if (min[j] > scene_aabb_min_thread[i][j]) {
				min[j] = scene_aabb_min_thread[i][j];
			}
		}
	}
	for (int i = 0; i < 3; ++i) {
		max[i] += 0.1;
		min[i] -= 0.1;
	}
}


//SCENE_AABB
void SpatialHashing::getSceneAABB(int thread_No)
{
	double* max = scene_aabb_max_thread[thread_No].data();
	double* min = scene_aabb_min_thread[thread_No].data();
	memset(max, 0xFE, 24); //set double to -5.31401e+303
	memset(min, 0x7F, 24); //set double to 1.38242e+306
	AABB* obj_aabb;
	int vertex_index_begin;
	int vertex_index_end;
	double* min_vertex;
	double* max_vertex;
	for (int i = 0; i < cloth->size(); ++i) {
		obj_aabb = (*cloth)[i].vertex_AABB.data();
		vertex_index_begin = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_index_end = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
		for (int j = vertex_index_begin; j < vertex_index_end; ++j) {
			min_vertex = obj_aabb[j].min;
			max_vertex = obj_aabb[j].max;
			for (int k = 0; k < 3; ++k) {
				if (max[k] < max_vertex[k]) {
					max[k] = max_vertex[k];
				}
				if (min[k] > min_vertex[k]) {
					min[k] = min_vertex[k];
				}				
			}
		}
	}
	for (int i = 0; i < collider->size(); ++i) {
		obj_aabb = (*collider)[i].vertex_AABB.data();
		vertex_index_begin = (*collider)[i].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_index_end = (*collider)[i].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
		for (int j = vertex_index_begin; j < vertex_index_end; ++j) {
			min_vertex = obj_aabb[j].min;
			max_vertex = obj_aabb[j].max;
			for (int k = 0; k < 3; ++k) {
				if (max[k] < max_vertex[k]) {
					max[k] = max_vertex[k];
				}
				if (min[k] > min_vertex[k]) {
					min[k] = min_vertex[k];
				}
			}
		}
	}

	for (int i = 0; i < tetrahedron->size(); ++i) {
		obj_aabb = (*tetrahedron)[i].vertex_AABB.data();
		vertex_index_begin = (*tetrahedron)[i].mesh_struct.vertex_index_on_surface_begin_per_thread[thread_No];
		vertex_index_end = (*tetrahedron)[i].mesh_struct.vertex_index_on_surface_begin_per_thread[thread_No + 1];
		for (int j = vertex_index_begin; j < vertex_index_end; ++j) {
			min_vertex = obj_aabb[j].min;
			max_vertex = obj_aabb[j].max;
			for (int k = 0; k < 3; ++k) {
				if (max[k] < max_vertex[k]) {
					max[k] = max_vertex[k];
				}
				if (min[k] > min_vertex[k]) {
					min[k] = min_vertex[k];
				}
			}
		}
	}
}