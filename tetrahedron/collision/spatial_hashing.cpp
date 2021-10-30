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

	spatial_hashing_cloth_triangle.resize(thread_num);
	spatial_hashing_collider_triangle.resize(thread_num);
	for (int i = 0; i < thread_num; ++i) {
		spatial_hashing_cloth_triangle[i].resize(cloth->size());
		spatial_hashing_collider_triangle[i].resize(collider->size());
		for (int j = 0; j < cloth->size(); ++j) {
			spatial_hashing_cloth_triangle[i][j].reserve(8 * (*cloth)[j].mesh_struct.triangle_indices.size() / thread_num);
		}
		for (int j = 0; j < collider->size(); ++j) {
			spatial_hashing_collider_triangle[i][j].reserve(8 * (*collider)[j].mesh_struct.triangle_indices.size() / thread_num);
		}
	}
	radix_sort.initial(thread);
}


void SpatialHashing::setSpatialHashing()
{
	getSceneAABB();
	for (int i = 0; i < 3; ++i) {
		cell_number[i] = (int)floor((scene_aabb.max[i] - scene_aabb.min[i]) / cell_length) + 1;
		hash_max_index[i] = cell_number[i] - 1;
	}
	cell_num0_cell_num1 = cell_number[0] * cell_number[1];
	initialSpatialHashingList();
}

void SpatialHashing::culling()
{
	setSpatialHashing();
	thread->assignTask(this, TRIANGLE_HASHING);

}






void SpatialHashing::initialSpatialHashingList()
{
	for (int i = 0; i < thread_num; ++i) {
		for (int j = 0; j < cloth->size(); ++j) {
			spatial_hashing_cloth_triangle[i][j].clear();
		}
		for (int j = 0; j < collider->size(); ++j) {
			spatial_hashing_collider_triangle[i][j].clear();
		}
	}
}

//TRIANGLE_HASHING
void SpatialHashing::triangleHashing(int thread_No)
{
	int* triangle_begin;
	AABB* aabb;
	std::vector<std::array<int, 2>>* spatial_hashing_triangle;
	for (int i = 0; i < cloth->size(); ++i) {
		aabb = (*cloth)[i].triangle_AABB.data();
		spatial_hashing_triangle = &spatial_hashing_cloth_triangle[thread_No][i];
		triangle_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
		for (int j = triangle_begin[thread_No]; j < triangle_begin[thread_No + 1]; ++j) {
			triangleHashingValue(aabb[j], spatial_hashing_triangle, j);
		}
	}
	for (int i = 0; i < collider->size(); ++i) {
		aabb = (*collider)[i].triangle_AABB.data();
		spatial_hashing_triangle = &spatial_hashing_collider_triangle[thread_No][i];
		triangle_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread.data();
		for (int j = triangle_begin[thread_No]; j < triangle_begin[thread_No + 1]; ++j) {
			triangleHashingValue(aabb[j], spatial_hashing_triangle, j);
		}
	}
}


void SpatialHashing::triangleHashingValue(AABB& aabb, std::vector<std::array<int, 2>>* spatial_hashing_triangle, int triangle_index)
{
	int min_index[3];
	int max_index[3];
	int index_z_multi;
	int index_y_multi;	

	for (int j = 0; j < 3; ++j) {
		min_index[j] = (int)floor((aabb.min[j]- scene_aabb.min[j]) / cell_length);
		max_index[j] = (int)floor((aabb.max[j] - scene_aabb.min[j]) / cell_length);
	}
	for (int index_y = min_index[1]; index_y <= max_index[1]; ++index_y) {
		index_y_multi = index_y * cell_number[0];
		for (int index_z = min_index[2]; index_z <= max_index[2]; ++index_z) {
			index_z_multi = index_z * cell_num0_cell_num1;
			for (int index_x = min_index[0]; index_x <= max_index[0]; ++index_x) {
				spatial_hashing_triangle->push_back({ triangle_index, index_x + index_y_multi + index_z_multi });
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
			max[j] = myMax(max[j], scene_aabb_max_thread[i][j]);
			min[j] = myMin(min[j], scene_aabb_min_thread[i][j]);
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
	int vertex_num;
	int vertex_index_begin;
	int vertex_index_end;
	double* min_vertex;
	double* max_vertex;
	for (int i = 0; i < cloth->size(); ++i) {
		vertex_num = (*cloth)[i].mesh_struct.vertex_position.size();
		obj_aabb = (*cloth)[i].vertex_AABB.data();
		vertex_index_begin = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_index_end = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
		for (int j = vertex_index_begin; j < vertex_index_end; ++j) {
			min_vertex = obj_aabb[j].min;
			max_vertex = obj_aabb[j].max;
			for (int k = 0; k < 3; ++k) {
				max[k] = myMax(max[k], max_vertex[k]);
				min[k] = myMin(min[k], min_vertex[k]);
			}
		}
	}
	for (int i = 0; i < collider->size(); ++i) {
		vertex_num = (*collider)[i].mesh_struct.vertex_position.size();
		obj_aabb = (*collider)[i].vertex_AABB.data();
		vertex_index_begin = (*collider)[i].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_index_end = (*collider)[i].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
		for (int j = vertex_index_begin; j < vertex_index_end; ++j) {
			min_vertex = obj_aabb[j].min;
			max_vertex = obj_aabb[j].max;
			for (int k = 0; k < 3; ++k) {
				max[k] = myMax(max[k], max_vertex[k]);
				min[k] = myMin(min[k], min_vertex[k]);
			}
		}
	}
}