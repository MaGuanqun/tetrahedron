//Here we list the <element,spatial hashing index> pairs, avoid creating cells.
//use radix sorting to sort pairs by index -> for searching
#pragma once
#include"../thread.h"
#include"../object/cloth.h"
#include"../object/collider.h"
#include"../object/tetrahedron.h"
#include"parallel_radix_sort.h"

class SpatialHashing
{
public:
	void setInObject(std::vector<Cloth>* cloth, std::vector<Collider>* collider,std::vector<Tetrahedron>* tetrahedron, 
		Thread* thread, double* tolerance_ratio);
	void getSceneAABB(int thread_No);
	void triangleHashing(int thread_No);
	void buildSpatialHashing();
	void setHashTogether(int thread_No);
	void searchTriangle(AABB& aabb, unsigned int input_obj_No, unsigned int triangle_index,
		std::vector<int>* obj_neighbor_index, bool is_collider, unsigned int thread_No);
private:
	Thread* thread;
	std::vector<Cloth>* cloth;
	std::vector<Collider>* collider;
	std::vector<Tetrahedron>* tetrahedron;

	unsigned int tetrahedron_begin_obj_index;

	void initialHashCellLength(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, double& cell_length, double* tolerance_ratio);
	double cell_length;//this is the size of spatial hashing cell, = average edge length
	double thread_num;

	AABB scene_aabb;
	std::vector<std::array<double, 3>>scene_aabb_max_thread;
	std::vector<std::array<double, 3>>scene_aabb_min_thread;
	void getSceneAABB();
	unsigned int cell_number[3];
	unsigned int cell_num0_cell_num1;
	unsigned int hash_max_index[3];
	void setSpatialHashing();
	//std::vector<std::vector<std::array<int, 3>>> spatial_hashing_triangle_per_thread;
	//std::vector<std::array<int, 3>> spatial_hashing_triangle;

	//std::vector<std::vector<unsigned int>> spatial_hashing_value_per_thread;
	unsigned int* spatial_hashing_value;
	//std::vector<std::vector<unsigned int>>spatial_hashing_obj_index_per_thread;
	unsigned int* spatial_hashing_obj_index;
	//std::vector<std::vector<unsigned int>>spatial_hashing_triangle_index_per_thread;
	unsigned int* spatial_hashing_triangle_index;
	
	unsigned int max_cell_count;
	
	std::vector<unsigned int> triangle_begin_per_obj;
	//std::vector<unsigned int> triangle_begin_per_tetrahedron;

	void triangleHashingValue(AABB& aabb, std::vector<unsigned int>* spatial_hashing_triangle_index,
		std::vector<unsigned int>* spatial_hashing_value, unsigned int triangle_index, std::vector<unsigned int>* hash_value);

	void triangleHashValue(AABB& aabb, unsigned int* spatial_hashing_triangle_index,
		unsigned int* spatial_hashing_value, unsigned int triangle_index, unsigned int* spatial_hashing_obj_index, unsigned int obj_index);

	void obtainTriangleHashingValue(AABB& aabb, std::vector<unsigned int>* hash_value);
	RadixSort radix_sort;
	void testRadixSort();

	unsigned int total_hash_size;

	unsigned int hash_table_size;

	unsigned int* obj_triangle_hash;
	//std::vector<std::vector<std::vector<unsigned int>>> obj_triangle_hash;
	void initialTriangleHash();
	void setHashTogether();
	std::vector<int>hash_value_begin;
	std::vector<int> prifix_sum;
	void setPrifixSum();
	bool*** obj_is_used;

	unsigned int largest_count_in_hash_value_list;
	unsigned int actual_hash_value_count;//total_hash_size-largest_count_in_hash_value_list
};

