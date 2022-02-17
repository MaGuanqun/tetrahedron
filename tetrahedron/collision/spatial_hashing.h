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
		Thread* thread, double* tolerance_ratio, unsigned int max_cell_count, bool for_construct_patch);
	void getSceneAABB(int thread_No);
	void triangleHashing(int thread_No);
	void buildSpatialHashing();
	void setHashTogether(int thread_No);
	//void searchTriangle(double* aabb, unsigned int input_obj_No, unsigned int triangle_index,
	//	std::vector<unsigned int>* obj_neighbor_index, bool is_collider, unsigned int thread_No);
	//void prifixSum1(int thread_No);
	//void prifixSum2(int thread_No);
	//void prifixSum3(int thread_No);
	//void prepareForActualHashValueCountThread(int thread_No);
	//void memsetThread(int thread_No);

	//void prefixSumParallelUp(int thread_No, unsigned int stage);
	//void prefixSumParallelDown(int thread_No, unsigned int stage);

	//void findPatch(std::vector<std::vector<std::vector<unsigned int>>>* triangle_patch_);	

	//void buildSpatialHashingPatchTriangle(double* scene_aabb);
	std::vector<unsigned int> actual_exist_cell_begin_per_thread;

	void findAllTrianglePairs(int thread_No);

	std::vector<std::vector<unsigned int>> triangle_pair; //except collider, inner vector store tri_1 index, obj_1_index, tri_2_index, obj_2_index
	std::vector<std::vector<unsigned int>> triangle_pair_with_collider; //inner vector store tri_1 index, obj_1_index, tri_2_index, obj_2_index

private:
	void deleteArray();
	Thread* thread;
	std::vector<Cloth>* cloth;
	std::vector<Collider>* collider;
	std::vector<Tetrahedron>* tetrahedron;

	unsigned int tetrahedron_begin_obj_index;
	unsigned int collider_begin_obj_index;

	void initialHashCellLength(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, double& cell_length, double* tolerance_ratio);
	double cell_length;//this is the size of spatial hashing cell, = average edge length
	unsigned int thread_num;

	double scene_aabb[6];
	std::vector<std::array<double, 6>>scene_aabb_thread;
	void getSceneAABB();
	unsigned int cell_number[3];
	unsigned int cell_num0_cell_num1;
	unsigned int hash_max_index[3];
	void setSpatialHashing();
	//void setSpatialHashingPatch(double* scene_aabb);
	//std::vector<std::vector<std::array<int, 3>>> spatial_hashing_triangle_per_thread;
	//std::vector<std::array<int, 3>> spatial_hashing_triangle;

	//std::vector<std::vector<unsigned int>> spatial_hashing_value_per_thread;
	unsigned int** spatial_hashing_value;
	//std::vector<std::vector<unsigned int>>spatial_hashing_obj_index_per_thread;
	unsigned int** spatial_hashing_obj_index;
	//std::vector<dstd::vector<unsigned int>>spatial_hashing_triangle_index_per_thread;
	unsigned int** spatial_hashing_triangle_index;

	unsigned int** spatial_hashing_value_collider;
	//std::vector<std::vector<unsigned int>>spatial_hashing_obj_index_per_thread;
	unsigned int** spatial_hashing_obj_index_collider;
	//std::vector<dstd::vector<unsigned int>>spatial_hashing_triangle_index_per_thread;
	unsigned int** spatial_hashing_triangle_index_collider;

	 
	unsigned int max_cell_count;
	
	std::vector<unsigned int> triangle_number_per_thread;
	std::vector<unsigned int> triangle_number_per_thread_collider;

	std::vector<std::vector<int>> triangle_begin_per_obj_per_thread;
	std::vector<std::vector<int>> triangle_begin_per_obj_per_thread_collider;
	//for instance, for obj 2 first index 7 at thread 3, if there are 5 indices before ,the start index is (5-7)=-2, so that 
	//we can directly get start index by 7+(-2)=5;
	std::vector<unsigned int> total_hash_size_per_thread;
	std::vector<unsigned int> total_hash_size_per_thread_collider;



	//std::vector<unsigned int> triangle_begin_per_tetrahedron;

	void triangleHashingValue(double* aabb, std::vector<unsigned int>* spatial_hashing_triangle_index,
		std::vector<unsigned int>* spatial_hashing_value, unsigned int triangle_index, std::vector<unsigned int>* hash_value);

	void triangleHashValue(double* aabb, unsigned int* spatial_hashing_triangle_index,
		unsigned int* spatial_hashing_value, unsigned int triangle_index, unsigned int* spatial_hashing_obj_index, unsigned int obj_index);

	void obtainTriangleHashingValue(double* aabb, std::vector<unsigned int>* hash_value);
	RadixSort* radix_sort;
	RadixSort* radix_sort_collider;
	void testRadixSort();

	//unsigned int total_hash_size;

	unsigned int hash_table_size;
	//std::vector<std::vector<std::vector<unsigned int>>> obj_triangle_hash;
	void initialTriangleHash();
	void setHashTogether();
	//std::vector<int>hash_value_begin;
	std::vector<std::vector<unsigned int>> prifix_sum;
	std::vector<std::vector<unsigned int>> prifix_sum_collider; //record start & end respectively
	void setPrifixSum();
	bool*** obj_is_used0;
	bool*** obj_is_used1;
	bool*** collider_is_used0;

	std::vector<unsigned int*> actual_hash_value_end_index_ref;//total_hash_size-largest_count_in_hash_value_list-1
	std::vector<unsigned int*> actual_hash_value_end_index_ref_collider;//total_hash_size-largest_count_in_hash_value_list-1
	std::vector<unsigned int> actual_hash_value_count_per_thread;//total_hash_size-largest_count_in_hash_value_list-1
	std::vector<unsigned int> actual_hash_value_count_per_thread_collider;//total_hash_size-largest_count_in_hash_value_list-1

	std::vector<std::vector<unsigned int>> hash_index_count_per_thread;
	std::vector<std::vector<unsigned int>> hash_index_count_per_thread_collider;

	//unsigned int* actual_hash_count_start_per_thread;
	//unsigned int* total_hash_count_start_per_thread;
	//unsigned int* total_hash_count_start_per_thread_move_1;//all index add 1

	//unsigned int* prefix_sum_thread_start;//record the count before the start index

	//unsigned int* hash_value_count_start_thread;//we record the count of every start of total_hash_count_start_per_thread of prefix_sum
	//void prepareForActualHashValueCount();
	//void testPrefixSumTime();
	//void testPrifixSum1();

	//void prefixSumParallel();
	//void initialParallePrefix();

	unsigned int* d_for_prefix_sum;//2^d

	//unsigned int* prefix_sum_1_address;

	int total_obj_num;
	bool for_construct_patch;

	void countHashIndexPerThread(int thread_No, unsigned int actual_count, std::vector<unsigned int>* hash_index_count,
		unsigned int* spatial_hashing_value_);

	

	bool has_collider;


	//reorganize the data in every object
	std::vector<std::array<double, 6>*> obj_tri_aabb;
	std::vector<std::array<double, 6>*> collider_tri_aabb;


	

};

