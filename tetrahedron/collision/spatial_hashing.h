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
	void setInObject(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
		Thread* thread, double* tolerance_ratio, unsigned int max_cell_count, bool for_construct_patch,
		unsigned int max_index_number_in_one_cell,
		unsigned int max_index_number_in_one_cell_collider);
	void getSceneAABB(int thread_No);
	void triangleHashing(int thread_No);
	void buildSpatialHashing(double* scene_aabb);
	void setHashTogether(int thread_No);

	void triangleHashingSmallerHashTable(int thread_No);


	//void prifixSum1(int thread_No);
	//void prifixSum2(int thread_No);
	//void prifixSum3(int thread_No);
	//void prepareForActualHashValueCountThread(int thread_No);
	//void memsetThread(int thread_No);

	//void prefixSumParallelUp(int thread_No, unsigned int stage);
	//void prefixSumParallelDown(int thread_No, unsigned int stage);

	//void findPatch(std::vector<std::vector<std::vector<unsigned int>>>* triangle_patch_);	

	//void buildSpatialHashingPatchTriangle(double* scene_aabb);
	//std::vector<unsigned int> actual_exist_cell_begin_per_thread;
	std::vector<unsigned int> cell_begin_per_thread;

	void findAllTrianglePairs(int thread_No);
	void findAllPairsHashTable(int thread_No);

	unsigned int** vertex_triangle_pair; //except collider, inner vector store vertex_1 index, obj_1_index, tri_2_index, obj_2_index
	unsigned int** vertex_obj_triangle_collider_pair; //inner vector store vertex_1 index, obj_1_index, tri_2_index, obj_2_index
	unsigned int** vertex_collider_triangle_obj_pair; //inner vector store vertex_1 index, obj_1_index, tri_2_index, obj_2_index

	unsigned int** edge_edge_pair;
	unsigned int** edge_edge_pair_collider; //first obj, second collider


	void test();


	////std::vector<std::vector<unsigned int>> spatial_hashing_value_per_thread;
	//unsigned int** spatial_hashing_triangle_value;
	//unsigned int** spatial_hashing_vertex_value;
	//unsigned int** spatial_hashing_edge_value;
	////std::vector<std::vector<unsigned int>>spatial_hashing_obj_index_per_thread;
	////unsigned int** spatial_hashing_obj_index;
	////std::vector<dstd::vector<unsigned int>>spatial_hashing_triangle_index_per_thread;
	//unsigned int** spatial_hashing_triangle_index; //in every vector, we store triangle index, obj_index respectively
	//unsigned int** spatial_hashing_vertex_index; //in every vector, we store triangle index, obj_index respectively
	//unsigned int** spatial_hashing_edge_index; //in every vector, we store triangle index, obj_index respectively

	//unsigned int** spatial_hashing_triangle_collider;
	//unsigned int** spatial_hashing_vertex_collider;
	//unsigned int** spatial_hashing_edge_collider;
	////std::vector<std::vector<unsigned int>>spatial_hashing_obj_index_per_thread;
	////unsigned int** spatial_hashing_obj_index_collider;
	////std::vector<dstd::vector<unsigned int>>spatial_hashing_triangle_index_per_thread;
	//unsigned int** spatial_hashing_triangle_index_collider;//in every vector, we store triangle index, obj_index respectively
	//unsigned int** spatial_hashing_vertex_index_collider;//in every vector, we store triangle index, obj_index respectively
	//unsigned int** spatial_hashing_edge_index_collider;//in every vector, we store triangle index, obj_index respectively

	//std::vector<std::vector<unsigned int>> prefix_sum;
	//std::vector<std::vector<unsigned int>> prefix_sum_collider; //record start & end respectively

	std::vector<unsigned int>hash_value_for_test;
	double scene_aabb[6];
	unsigned int cell_number[3];
	unsigned int cell_num0_cell_num1;
	double cell_length;//this is the size of spatial hashing cell, = average edge length




	unsigned int*** spatial_hashing_cell_triangle;//in every vector, we store triangle index,obj_index respectively
	unsigned int*** spatial_hashing_cell_edge;//in every vector, we store triangle index,obj_index respectively
	unsigned int*** spatial_hashing_cell_vertex;//in every vector, we store triangle index,obj_index respectively
	//unsigned int*** spatial_hashing_actual_hash_value_for_test;//in every vector, we store triangle index,obj_index respectively
	unsigned int*** spatial_hashing_cell_collider_triangle;//in every vector, we store triangle index, obj_index respectively
	unsigned int*** spatial_hashing_cell_collider_edge;//in every vector, we store triangle index, obj_index respectively
	unsigned int*** spatial_hashing_cell_collider_vertex;//in every vector, we store triangle index, obj_index respectively
	unsigned int hash_cell_count;


	//unsigned int*** spatial_hashing_cell;
	//unsigned int*** spatial_hashing_cell_collider;


	void recordRealTriangleHashValue(int thread_No);
	void prefixSumMulti2(int thread_No);

	void recordNonEmptyCell(int thread_No);

	void findVertexTetPairHashTable(int thread_No);
	void tetHashingSmallerHashTable(int thread_No);

	std::vector<std::vector<unsigned int>> vertex_tet_pair; //except collider, inner vector store vertex_1 index, obj_1_index, tet_2_index, obj_2_index
	void obtainPairCount(int thread_No);
	void setHashCellPairNumPrefixSumTogether(int thread_No);

	//void findAllTrianglePairsHashTableCompare(int thread_No);

private:

	void deleteArray();
	Thread* thread;
	std::vector<Cloth>* cloth;
	std::vector<Collider>* collider;
	std::vector<Tetrahedron>* tetrahedron;

	unsigned int tetrahedron_begin_obj_index;
	unsigned int collider_begin_obj_index;

	void initialHashCellLength(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, double& cell_length, double* tolerance_ratio);

	unsigned int thread_num;


	std::vector<std::array<double, 6>>scene_aabb_thread;
	void getSceneAABB();

	unsigned int hash_max_index[3];
	void setSpatialHashing();
	//void setSpatialHashingPatch(double* scene_aabb);
	//std::vector<std::vector<std::array<int, 3>>> spatial_hashing_triangle_per_thread;
	//std::vector<std::array<int, 3>> spatial_hashing_triangle;




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
		unsigned int* spatial_hashing_value, unsigned int triangle_index, unsigned int obj_index);

	void obtainTriangleHashingValue(double* aabb, std::vector<unsigned int>* hash_value);
	//RadixSort* radix_sort;
	//RadixSort* radix_sort_collider;
	void testRadixSort();

	//unsigned int total_hash_size;

	unsigned int hash_table_size;

	void initialTriangleHash();
	void setHashTogether();
	//std::vector<int>hash_value_begin;

	void setPrifixSum();
	//bool*** obj_is_used;
	////bool*** obj_is_used1;
	//bool*** collider_is_used;

	//std::vector<unsigned int*> actual_hash_value_end_index_ref;//total_hash_size-largest_count_in_hash_value_list-1
	//std::vector<unsigned int*> actual_hash_value_end_index_ref_collider;//total_hash_size-largest_count_in_hash_value_list-1
	//std::vector<unsigned int> actual_hash_value_count_per_thread;//total_hash_size-largest_count_in_hash_value_list-1
	//std::vector<unsigned int> actual_hash_value_count_per_thread_collider;//total_hash_size-largest_count_in_hash_value_list-1

	//std::vector<std::vector<unsigned int>> hash_index_count_per_thread;
	//std::vector<std::vector<unsigned int>> hash_index_count_per_thread_collider;

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

	std::vector<std::array<double, 6>*> obj_vertex_aabb;
	std::vector<std::array<double, 6>*> collider_vertex_aabb;

	std::vector<std::array<double, 6>*> obj_edge_aabb;
	std::vector<std::array<double, 6>*> collider_edge_aabb;

	std::vector<unsigned int*> obj_triangle_index_begin_per_thread;
	std::vector<unsigned int*> obj_vertex_index_begin_per_thread;
	std::vector<unsigned int*> obj_edge_index_begin_per_thread;

	std::vector<unsigned int*> tetrahedron_vertex_index_on_surface;


	std::vector<std::array<int, 3>*> triangle_vertex_index;
	std::vector<std::array<int, 3>*> triangle_vertex_index_collider;

	std::vector<unsigned int*> edge_vertex_index;
	std::vector<unsigned int*> edge_vertex_index_collider;




	void testIfRadixSortIsRight(unsigned int* has_value, unsigned int count);

	std::vector<time_t> record_time0;
	std::vector<time_t> record_time1;
	std::vector<time_t> record_time2;

	void findTheMaxNumberOfTriangleInCell();

	void findListOfSpatialHashingCell();
	void initialHashCell(unsigned int total_triangle_num, unsigned int max_index_number_in_one_cell,
		unsigned int max_index_number_in_one_cell_collider);


	unsigned int max_index_number_in_one_cell;
	unsigned int max_index_number_in_one_cell_collider;

	void triangleHashValueWithRecord(double* aabb,
		std::vector<unsigned int>* spatial_hashing_cell, unsigned int triangle_index, unsigned int obj_index);

	void triangleHashValueWithoutRecord(double* aabb,
		unsigned int** spatial_hashing_cell, unsigned int* triangle_index, double* scene_aabb, double cell_length, unsigned int hash_cell_count,
		uint64_t P1, uint64_t P2, uint64_t P3, unsigned int* spatial_hashing_cell_triangle_size);

	std::uint64_t P1, P2, P3;
	void testHashCollision();
	bool hashCollisionHappened(unsigned int hash_index);

	void recordRealTriangleHashValue(double* aabb,
		std::vector<unsigned int>* spatial_hashing_cell, unsigned int triangle_index, unsigned int obj_index);

	bool hashCellIsNotEmpty(unsigned int hash_index);

	std::vector<unsigned int> triangle_pair_number_thread;
	void findTheNumberOfTriangleInCellSP();

	void findAllTrianglePairsHashTableTest(int thread_No);

	void findAllTrianglePairsTest(int thread_No);

	void testHashCollisionSP();

	bool hashCellIsNotEmptyTest(unsigned int hash_index);
	bool hashCollisionHappenedTest(unsigned int hash_index);

	unsigned int findNeareastPrimeNumber(unsigned int index);

	//size_t** time_stamp_cell;
	//size_t** time_stamp_cell_collider;
	void recordNonEmptyCell();
	unsigned int** non_empty_cell_index_vertex_triangle; //the first element store the actual element number 
	unsigned int** non_empty_cell_index_edge;

	void vertexTetIntersect(int thread_No, double* pos, unsigned int vertex_index, unsigned int obj_index);


	bool vertexFromTet(int vertex_index, int* tet_indices);
	bool testIntersect(double* pos, double* p0, double* p1, double* p2, double* p3);

	unsigned int* non_empty_cell_index_begin_per_thread_vertex_triangle;//
	unsigned int* non_empty_cell_index_begin_per_thread_edge;//

	//std::vector<std::vector<std::vector<unsigned int>>> obj_triangle_hash;
	void searchTriangle(double* aabb, unsigned int* hash_index, unsigned int hash_cell_size, unsigned int input_obj_No, unsigned int triangle_index,
		std::vector<unsigned int>* triangle_pair_, std::vector<unsigned int>* triangle_pair_with_collider_, unsigned int thread_No);

	void reorganzieDataOfObjects();

	std::vector<unsigned int> hash_cell_triangle_count;

	//unsigned int* pair_count_per_thread;//actually record the size of pair vector which = 4*pair count
	//unsigned int* pair_count_with_collider_per_thread;//actually record the size of pair vector which = 4*pair count
	//std::vector<unsigned int>* indicator;

	//unsigned int* prefix_sum_pair_count
	unsigned int* hash_cell_pair_num_prefix_vertex_triangle;
	unsigned int* hash_cell_pair_num_prefix_edge;
	//std::vector<unsigned int> hash_cell_triangle_num;

	void setPairAveInThread();

	unsigned int* non_empty_cell_index_ave_begin_per_thread_vertex_triangle;
	unsigned int* non_empty_cell_index_ave_begin_per_thread_edge;
	unsigned int* global_cell_start_vertex_triangle;//the size of size of outter is 2*thread_num
	unsigned int* global_cell_start_edge;//the size of size of outter is 2*thread_num

	//unsigned int max_num_to_loop_to_find_pair;//the size of inner array in global_cell_start

//	unsigned int max_pair_num_in_one_loop_per_thread;

	//unsigned int triangle_count_per_thread[8];
//	unsigned int triangle_count_with_collider_per_thread[8];

	//void findAllTrianglePairsHashTableSingleThread();


	void findMaxTriangleNumInATCell();



	unsigned int* prefix_sum_record_per_thread_start_vertex_triangle;
	unsigned int* prefix_sum_record_per_thread_start_edge;

	bool testA(unsigned int& a);

	//void setGlobalCellStartPerThreadPerLoop();

	//void findAllTrianglePairsHashTable();

	//void findAllTrianglePairsHashTableSingleThread();

	void testSpatialHashingCellSet0();

	unsigned int** spatial_hashing_cell_triangle_size;
	unsigned int** spatial_hashing_cell_edge_size;
	unsigned int** spatial_hashing_cell_vertex_size;
	unsigned int** spatial_hashing_cell_collider_triangle_size;
	unsigned int** spatial_hashing_cell_collider_edge_size;
	unsigned int** spatial_hashing_cell_collider_vertex_size;

	void clearHashCell();


	void findAllVertexTrianglePairs(int thread_No);
	void findAllEdgeEdgePairs(int thread_No);
};

