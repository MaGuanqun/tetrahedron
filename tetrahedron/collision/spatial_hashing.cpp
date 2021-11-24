#include"spatial_hashing.h"



void SpatialHashing::initialHashCellLength(std::vector<Cloth>* cloth, double& cell_length)
{
	double ave_length = 0;
	int edge_num = 0;
	std::vector<MeshStruct::Edge>* edge;
	for (int i = 0; i < cloth->size(); ++i) {
		edge = &(*cloth)[i].mesh_struct.edges;
		for(int j=0;j< edge->size();++j){
			ave_length += (*edge)[j].length;
		}
		edge_num += edge->size();
	}
	ave_length /= (double)edge_num;
	cell_length = ave_length;
}

void SpatialHashing::setInObject(std::vector<Cloth>* cloth, std::vector<Collider>* collider, 
	std::vector<Tetrahedron>* tetrahedron, Thread* thread)
{
	initialHashCellLength(cloth, cell_length);
	this->cloth = cloth;
	this->collider = collider;
	this->tetrahedron = tetrahedron;
	this->thread = thread;
	thread_num = thread->thread_num;
	scene_aabb_max_thread.resize(thread_num);
	scene_aabb_min_thread.resize(thread_num);
	spatial_hashing_value_per_thread.resize(thread_num);
	spatial_hashing_cloth_index_per_thread.resize(thread_num);
	spatial_hashing_triangle_index_per_thread.resize(thread_num);
	//spatial_hashing_triangle_per_thread.resize(thread_num);

	int total_cloth_num = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		total_cloth_num += cloth->size();
	}
	int estimate_hash_num= 64 * total_cloth_num / thread_num;
	spatial_hashing_value.reserve(estimate_hash_num);
	spatial_hashing_cloth_index.reserve(estimate_hash_num);
	spatial_hashing_triangle_index.reserve(estimate_hash_num);
	for (int i = 0; i < thread_num; ++i) {
		spatial_hashing_value_per_thread[i].reserve(estimate_hash_num);
		spatial_hashing_cloth_index_per_thread[i].reserve(estimate_hash_num);
		spatial_hashing_triangle_index_per_thread[i].reserve(estimate_hash_num);
	}

	//for (int i = 0; i < thread_num; ++i) {
	//	spatial_hashing_triangle_per_thread[i].reserve(64 * total_cloth_num / thread_num);
	//}

	initialTriangleHash();
	hash_value_begin.resize(thread_num,0);
	radix_sort.initial(thread);

	cloth_is_used = new bool** [thread_num];
	for (int i = 0; i < thread_num; ++i) {
		cloth_is_used[i] = new bool* [cloth->size()];
		for (int j = 0; j < cloth->size(); ++j) {
			cloth_is_used[i][j] = new bool[(*cloth)[j].mesh_struct.triangle_indices.size()];
			memset(cloth_is_used[i][j], 0, (*cloth)[j].mesh_struct.triangle_indices.size());
		}		
	}
	//testRadixSort();
	
}


void SpatialHashing::initialTriangleHash()
{
	cloth_triangle_hash.resize(cloth->size());
	for (int i = 0; i < cloth->size(); ++i) {
		cloth_triangle_hash[i].resize((*cloth)[i].mesh_struct.triangle_indices.size());
		for (int j = 0; j < cloth_triangle_hash[i].size(); ++j) {
			cloth_triangle_hash[i][j].reserve(32);
		}
	}
}

void SpatialHashing::testRadixSort()
{
	int vec_size = 4e5;
	std::vector<int> value;
	std::vector<int> triangle_index;
	std::vector<int> hash_cloth_No;
	value.reserve(vec_size);
	triangle_index.reserve(vec_size);
	hash_cloth_No.reserve(vec_size);
	time_t t1 = clock();
	for (int j = 0; j < 100; ++j) {
		int size = 100 * 100 * 100;
		value.clear();
		triangle_index.clear();
		hash_cloth_No.clear();
		for (int i = 0; i < vec_size; ++i) {
			value.push_back(rand() % size);
			triangle_index.push_back(i);
			hash_cloth_No.push_back(0);
		}
		for (int i = 0; i < vec_size; ++i) {
			if (i % 2 == 0) {
				value[i] += 100000;
			}
		}
	}
	t1 = clock() - t1;
	std::cout << "index " << t1 << std::endl;
	time_t t = clock();
	for (int j = 0; j < 100; ++j) {
		int size = 100 * 100 * 100;
		value.clear();
		triangle_index.clear();
		hash_cloth_No.clear();
		for (int i = 0; i < vec_size; ++i) {
			value.push_back(rand() % size);
			triangle_index.push_back(i);
			hash_cloth_No.push_back(0);
		}
		for (int i = 0; i < vec_size; ++i) {
			if (i % 2 == 0) {
				value[i] += 100000;
			}
		}
		//for (int i = 0; i < vec_size; ++i) {
		//	std::cout << cell[i][1] << " ";
		//}
		//std::cout << std::endl;
		radix_sort.radixSort(size, &value, &triangle_index, &hash_cloth_No);
	}
	//for (int i = 0; i < vec_size; ++i) {
	//	std::cout << cell[i][1] << " ";
	//}
	//std::cout << std::endl;
	//}
	std::cout <<"time "<< clock() - t<<" "<< clock() - t - t1 << std::endl;
}

void SpatialHashing::setSpatialHashing()
{
	getSceneAABB();
	for (int i = 0; i < 3; ++i) {
		cell_number[i] = (int)floor((scene_aabb.max[i] - scene_aabb.min[i]) / cell_length) + 1;
		hash_max_index[i] = cell_number[i] - 1;
	}
	cell_num0_cell_num1 = cell_number[0] * cell_number[1];
}

void SpatialHashing::buildSpatialHashing()
{
	setSpatialHashing();
	thread->assignTask(this, TRIANGLE_HASHING);
	setHashTogether();
	radix_sort.radixSort(cell_number[0] * cell_number[1] * cell_number[2], &spatial_hashing_value, &spatial_hashing_triangle_index, &spatial_hashing_cloth_index);
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

	prifix_sum.resize(cell_number[0] * cell_number[1] * cell_number[2] + 1);
	memset(prifix_sum.data(), 0, 4 * prifix_sum.size());
	for (int i = 0; i < spatial_hashing_value.size(); ++i) {
		prifix_sum[spatial_hashing_value[i]]++;
	}
	int s = 0; int t;
	for (int i = 0; i < prifix_sum.size(); ++i) {
		t = s + prifix_sum[i];
		prifix_sum[i] = s;
		s = t;
	}
	prifix_sum[prifix_sum.size() - 1] = t;
}

void SpatialHashing::searchTriangle(AABB& aabb,int obj_No, int triangle_index, std::vector<int>* cloth_neighbor_index, bool is_collider, int thread_No)
{
	std::vector<int>cloth_index_record;	std::vector<int>triangle_index_record;
	cloth_index_record.reserve(20);
	triangle_index_record.reserve(20);
	bool** cloth_is_used_ = cloth_is_used[thread_No];
	int hash_cell_index_end;
	int cloth_No; int current_triangle_index;
	for (int i = 0; i < cloth->size(); ++i) {
		cloth_neighbor_index[i].clear();
		cloth_neighbor_index[i].reserve(10);
	}
	if (is_collider) {
		std::vector<int>hash_value;
		obtainTriangleHashingValue(aabb, &hash_value);
		for (int i = 0; i < hash_value.size(); ++i) {
			hash_cell_index_end = prifix_sum[hash_value[i]+1];
			for (int j = prifix_sum[hash_value[i]]; j < hash_cell_index_end; ++j) {
				cloth_No = spatial_hashing_cloth_index[j];
				current_triangle_index= spatial_hashing_triangle_index[j];
				if (!cloth_is_used_[cloth_No][current_triangle_index]) {
					cloth_is_used_[cloth_No][current_triangle_index] = true;
					cloth_index_record.push_back(cloth_No);
					triangle_index_record.push_back(current_triangle_index);
					if (aabb.AABB_intersection((*cloth)[cloth_No].triangle_AABB[current_triangle_index])) {
						cloth_neighbor_index[cloth_No].push_back(current_triangle_index);
					}				
				}				
			}
		}
	}
	else {
		std::vector<int>* hash_value;
		hash_value = &cloth_triangle_hash[obj_No][triangle_index];
		for (int i = 0; i < hash_value->size(); ++i) {
			hash_cell_index_end = prifix_sum[(*hash_value)[i]+1];
			for (int j = prifix_sum[(*hash_value)[i]]; j < hash_cell_index_end; ++j) {
				cloth_No = spatial_hashing_cloth_index[j];
				current_triangle_index = spatial_hashing_triangle_index[j];
				if (!cloth_is_used_[cloth_No][current_triangle_index]) {
					cloth_is_used_[cloth_No][current_triangle_index] = true;
					cloth_index_record.push_back(cloth_No);
					triangle_index_record.push_back(current_triangle_index);
					if (aabb.AABB_intersection((*cloth)[cloth_No].triangle_AABB[current_triangle_index])) {
						cloth_neighbor_index[cloth_No].push_back(current_triangle_index);
					}
				}
			}
		}
	}
	for (int i = 0; i < triangle_index_record.size(); ++i) {
		cloth_is_used_[cloth_index_record[i]][triangle_index_record[i]] = false;
	}
}

void SpatialHashing::setHashTogether()
{
	int total_hash_num = spatial_hashing_value_per_thread[0].size();
	for (int i = 1; i < thread_num; ++i) {
		hash_value_begin[i] = total_hash_num;
		total_hash_num+= spatial_hashing_value_per_thread[i].size();
		
	}
	spatial_hashing_value.resize(total_hash_num);
	spatial_hashing_cloth_index.resize(total_hash_num);
	spatial_hashing_triangle_index.resize(total_hash_num);
	thread->assignTask(this, SET_HASH_TOGETHER);
}

//SET_HASH_TOGETHER
void SpatialHashing::setHashTogether(int thread_No)
{
	memcpy(&spatial_hashing_value[hash_value_begin[thread_No]], spatial_hashing_value_per_thread[thread_No].data(), 4 * spatial_hashing_value_per_thread[thread_No].size());
	memcpy(&spatial_hashing_cloth_index[hash_value_begin[thread_No]], spatial_hashing_cloth_index_per_thread[thread_No].data(), 4 * spatial_hashing_value_per_thread[thread_No].size());
	memcpy(&spatial_hashing_triangle_index[hash_value_begin[thread_No]], spatial_hashing_triangle_index_per_thread[thread_No].data(), 4 * spatial_hashing_value_per_thread[thread_No].size());
}


//TRIANGLE_HASHING
void SpatialHashing::triangleHashing(int thread_No)
{
	int* triangle_begin;
	AABB* aabb;
	std::vector<int>* spatial_hashing_value_;
	std::vector<int>* spatial_hashing_cloth_;
	std::vector<int>* spatial_hashing_triangle_;
	std::vector<int>* hash_value;
	spatial_hashing_cloth_ = &spatial_hashing_cloth_index_per_thread[thread_No];
	spatial_hashing_triangle_ = &spatial_hashing_triangle_index_per_thread[thread_No];
	spatial_hashing_value_ = &spatial_hashing_value_per_thread[thread_No];
	spatial_hashing_cloth_->clear();
	spatial_hashing_triangle_->clear();
	spatial_hashing_value_->clear();
	int vector_size;
	for (int i = 0; i < cloth->size(); ++i) {
		aabb = (*cloth)[i].triangle_AABB.data();		
		triangle_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
		hash_value = cloth_triangle_hash[i].data();
		vector_size = spatial_hashing_triangle_->size();
		for (int j = triangle_begin[thread_No]; j < triangle_begin[thread_No + 1]; ++j) {
			clothTriangleHashingValue(aabb[j], spatial_hashing_triangle_,spatial_hashing_value_, j,&hash_value[j]);
		}
		spatial_hashing_cloth_->insert(spatial_hashing_cloth_->end(), spatial_hashing_triangle_->size() - vector_size, i);
	}
}

void SpatialHashing::obtainTriangleHashingValue(AABB& aabb, std::vector<int>* hash_value)
{
	int min_index[3];
	int max_index[3];
	int index_z_multi;
	int index_y_multi;
	for (int j = 0; j < 3; ++j) {
		min_index[j] = (int)floor((aabb.min[j] - scene_aabb.min[j]) / cell_length);
		max_index[j] = (int)floor((aabb.max[j] - scene_aabb.min[j]) / cell_length) + 1;
	}
	int size = (max_index[0] - min_index[0]) * (max_index[1] - min_index[1]) * (max_index[2] - min_index[2]);
	//hash_value.clear();
	hash_value->reserve(size);
	for (int index_y = min_index[1]; index_y < max_index[1]; ++index_y) {
		index_y_multi = index_y * cell_number[0];
		for (int index_z = min_index[2]; index_z < max_index[2]; ++index_z) {
			index_z_multi = index_z * cell_num0_cell_num1;
			for (int index_x = min_index[0]; index_x < max_index[0]; ++index_x) {
				hash_value->push_back(index_x + index_y_multi + index_z_multi);
			}
		}
	}
}

void SpatialHashing::clothTriangleHashingValue(AABB& aabb, std::vector<int>* spatial_hashing_triangle_index, std::vector<int>* spatial_hashing_value, int triangle_index, std::vector<int>* hash_value)
{
	int min_index[3];
	int max_index[3];
	int index_z_multi;
	int index_y_multi;	
	for (int j = 0; j < 3; ++j) {
		min_index[j] = (int)floor((aabb.min[j]- scene_aabb.min[j]) / cell_length);
		max_index[j] = (int)floor((aabb.max[j] - scene_aabb.min[j]) / cell_length) + 1;
	}
	int size = (max_index[0] - min_index[0]) * (max_index[1] - min_index[1]) * (max_index[2] - min_index[2]);
	hash_value->clear();
	int value;
	for (int index_y = min_index[1]; index_y < max_index[1]; ++index_y) {
		index_y_multi = index_y * cell_number[0];
		for (int index_z = min_index[2]; index_z < max_index[2]; ++index_z) {
			index_z_multi = index_z * cell_num0_cell_num1;
			for (int index_x = min_index[0]; index_x < max_index[0]; ++index_x) {
				value = index_x + index_y_multi + index_z_multi;
				hash_value->push_back(value);
			}
		}
	}	
	spatial_hashing_value->insert(spatial_hashing_value->end(), hash_value->begin(),hash_value->end());
	spatial_hashing_triangle_index->insert(spatial_hashing_triangle_index->end(), hash_value->size(), triangle_index);
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
}