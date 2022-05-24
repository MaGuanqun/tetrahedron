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
		Thread* thread, double* tolerance_ratio, unsigned int max_cell_count, 
		unsigned int max_index_number_in_one_cell,
		unsigned int max_index_number_in_one_cell_collider, unsigned int estimate_coeff_for_pair_num);
	void getSceneAABB(int thread_No);

	void buildSpatialHashing(double* scene_aabb);


	void triangleHashingSmallerHashTable(int thread_No);

	std::vector<unsigned int> cell_begin_per_thread;




	void findAllPairsHashTable(int thread_No);

	unsigned int** vertex_triangle_pair; //except collider, inner vector store vertex_1 index, obj_1_index, tri_2_index, obj_2_index
	unsigned int** vertex_obj_triangle_collider_pair; //inner vector store vertex_1 index, obj_1_index, tri_2_index, obj_2_index
	unsigned int** vertex_collider_triangle_obj_pair; //inner vector store vertex_1 index, obj_1_index, tri_2_index, obj_2_index

	unsigned int** edge_edge_pair;
	unsigned int** edge_edge_pair_collider; //first obj, second collider



	std::vector<unsigned int>hash_value_for_test;
	double scene_aabb[6];

	double cell_length;//this is the size of spatial hashing cell, = average edge length




	unsigned int** spatial_hashing_cell_triangle;//in every vector, we store triangle index,obj_index respectively
	unsigned int** spatial_hashing_cell_edge;//in every vector, we store triangle index,obj_index respectively
	unsigned int** spatial_hashing_cell_vertex;//in every vector, we store triangle index,obj_index respectively
	//unsigned int*** spatial_hashing_actual_hash_value_for_test;//in every vector, we store triangle index,obj_index respectively
	unsigned int** spatial_hashing_cell_collider_triangle;//in every vector, we store triangle index, obj_index respectively
	unsigned int** spatial_hashing_cell_collider_edge;//in every vector, we store triangle index, obj_index respectively
	unsigned int** spatial_hashing_cell_collider_vertex;//in every vector, we store triangle index, obj_index respectively
	unsigned int hash_cell_count;

	void recordNonEmptyCell(int thread_No);

	std::vector<std::vector<unsigned int>> vertex_tet_pair; //except collider, inner vector store vertex_1 index, obj_1_index, tet_2_index, obj_2_index
	void obtainPairCount(int thread_No);
	void setHashCellPairNumPrefixSumTogether(int thread_No);

	//void findAllTrianglePairsHashTableCompare(int thread_No);

	void setPairAveInThread(int thread_No);

	void findAllPairsHashTableElementwise(int thread_No);


	std::vector<std::vector<unsigned int>> ori_hash_value;// the hash value computed by length * width * height
	void initialOriHashValue();

	unsigned int cell_number[3];
	unsigned int cell_num0_cell_num1;

	void oriTriangleHashing(int thread_No);

	void buildSpatialHashingForOri(double* scene_aabb);

	void selectCell(double* start_pos, double* dir, std::vector<unsigned int>& select_hash_index
		, std::vector<unsigned int>& initial_hash_index);

	void 	findAllElementsInOneCell(unsigned int select_hash_index, unsigned int ori_hash_index,
		std::vector<std::vector<unsigned int>>& vertex_index, std::vector<std::vector<unsigned int>>& triangle_index,
		std::vector<std::vector<unsigned int>>& edge_index);

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
	unsigned int max_cell_count;

	void initialTriangleHash();

	//std::vector<int>hash_value_begin;

	bool*** obj_is_used;


	unsigned int* d_for_prefix_sum;//2^d

	//unsigned int* prefix_sum_1_address;

	int total_obj_num;
	bool for_construct_patch;

	void countHashIndexPerThread(int thread_No, unsigned int actual_count, std::vector<unsigned int>* hash_index_count,
		unsigned int* spatial_hashing_value_);



	bool has_collider;


	std::vector<std::array<double, 6>*> obj_tri_aabb;
	std::vector<std::array<double, 6>*> collider_tri_aabb;

	std::vector<std::array<double, 6>*> obj_vertex_aabb;
	std::vector<std::array<double, 6>*> collider_vertex_aabb;

	std::vector<std::array<double, 6>*> obj_edge_aabb;
	std::vector<std::array<double, 6>*> collider_edge_aabb;

	std::vector<unsigned int*> obj_triangle_index_begin_per_thread;
	std::vector<unsigned int*> obj_vertex_index_begin_per_thread;
	std::vector<unsigned int*> obj_edge_index_begin_per_thread;


	std::vector<unsigned int*> collider_vertex_index_begin_per_thread;

	std::vector<unsigned int*> tetrahedron_vertex_index_on_surface;


	std::vector<std::array<int, 3>*> triangle_vertex_index;
	std::vector<std::array<int, 3>*> triangle_vertex_index_collider;

	std::vector<unsigned int*> edge_vertex_index;
	std::vector<unsigned int*> edge_vertex_index_collider;



	void initialHashCell(unsigned int total_triangle_num, unsigned int max_index_number_in_one_cell,
		unsigned int max_index_number_in_one_cell_collider, unsigned int estimate_coeff_for_pair_num);


	unsigned int max_index_number_in_one_cell;
	unsigned int max_index_number_in_one_cell_edge;
	unsigned int max_index_number_in_one_cell_vertex;
	unsigned int max_index_number_in_one_cell_collider;
	unsigned int max_index_number_in_one_cell_collider_vertex;
	unsigned int max_index_number_in_one_cell_collider_edge;

	void triangleHashValueWithRecord(double* aabb,
		unsigned int* spatial_hashing_cell, unsigned int* triangle_index, double* scene_aabb, double cell_length,
		unsigned int hash_cell_count, uint64_t P1, uint64_t P2, uint64_t P3, unsigned int* spatial_hashing_cell_triangle_size,
		unsigned int* obj_triangle_hash, unsigned int max_index_number_in_one_cell);

	void triangleHashValueWithoutRecord(double* aabb,
		unsigned int* spatial_hashing_cell, unsigned int* triangle_index, double* scene_aabb, double cell_length, unsigned int hash_cell_count,
		uint64_t P1, uint64_t P2, uint64_t P3, unsigned int* spatial_hashing_cell_triangle_size, unsigned int max_index_number_in_one_cell,
		int test_type);

	std::uint64_t P1, P2, P3;



	std::vector<unsigned int> triangle_pair_number_thread;




	unsigned int findNeareastPrimeNumber(unsigned int index);

	void recordNonEmptyCell();
	unsigned int** non_empty_cell_index_vertex_triangle; //the first element store the actual element number 
	unsigned int** non_empty_cell_index_edge;


	unsigned int* non_empty_cell_index_begin_per_thread_vertex_triangle;//
	unsigned int* non_empty_cell_index_begin_per_thread_edge;//

	unsigned int** obj_edge_hash;
	unsigned int** obj_vertex_hash;
	unsigned int** collider_vertex_hash;



	void searchPrimitive(double* aabb, unsigned int* hash_index, unsigned int hash_cell_size, unsigned int input_obj_No, unsigned int triangle_index,
		unsigned int*& triangle_pair_, unsigned int*& triangle_pair_with_collider_, unsigned int thread_No,
		unsigned int* spatial_hashing_cell_, unsigned int* spatial_hash_size,
		unsigned int* spatial_hashing_cell_collider_, unsigned int* spatial_hash_size_collider, bool is_edge,
		std::vector<unsigned int>* neighbor_primitive_0, std::vector<unsigned int>* neighbor_primitive_1,
		unsigned int max_index_number_in_one_cell, unsigned int max_index_number_in_one_cell_collider);

	void searchPrimitiveOnCollider(double* aabb, unsigned int* hash_index, unsigned int hash_cell_size,
		unsigned int input_obj_No, unsigned int vertex_index,
		unsigned int*& triangle_pair_with_collider_, unsigned int thread_No,
		unsigned int* spatial_hashing_cell_, unsigned int* spatial_hash_size_, unsigned int max_index_number_in_one_cell_collider);

	void reorganzieDataOfObjects();

	std::vector<unsigned int> hash_cell_triangle_count;



	//unsigned int* prefix_sum_pair_count
	unsigned int* hash_cell_pair_num_prefix_vertex_triangle;
	unsigned int* hash_cell_pair_num_prefix_edge;
	//std::vector<unsigned int> hash_cell_triangle_num;

	void setPairAveInThreadVertexTriangle();
	void setPairAveInThreadEdgeEdge();

	unsigned int* non_empty_cell_index_ave_begin_per_thread_vertex_triangle;
	unsigned int* non_empty_cell_index_ave_begin_per_thread_edge;
	unsigned int* global_cell_start_vertex_triangle;//the size of size of outter is 2*thread_num
	unsigned int* global_cell_start_edge;//the size of size of outter is 2*thread_num



	unsigned int* prefix_sum_record_per_thread_start_vertex_triangle;
	unsigned int* prefix_sum_record_per_thread_start_edge;


	unsigned int** spatial_hashing_cell_triangle_size;
	unsigned int** spatial_hashing_cell_edge_size;
	unsigned int** spatial_hashing_cell_vertex_size;
	unsigned int** spatial_hashing_cell_collider_triangle_size;
	unsigned int** spatial_hashing_cell_collider_edge_size;
	unsigned int** spatial_hashing_cell_collider_vertex_size;

	void findAllVertexTrianglePairs(int thread_No);
	void findAllEdgeEdgePairs(int thread_No);

	unsigned int** primitive_index_record;
	unsigned int max_triangle_cell_size;
	unsigned int max_vertex_cell_size;


	std::vector<MeshStruct::Vertex*> vertices;
	std::vector<MeshStruct::Edge*> edges;
	
	void initialTriangleHashValue(double* aabb,
		std::vector<unsigned int>* spatial_hashing_cell, double* scene_aabb, double cell_length);
	void setSpatialHashingInitialCount();

	void 	collectAllVoxel(std::vector<unsigned int>& select_hash_index, std::vector<unsigned int>&
		initial_hash_index, double start__[3], double end__[3],
		double* scene_aabb);
	bool intersectSpatialHasingCube(double* AABB, double* start, double* direction, double* end);

	inline double frac0(double x);
	inline double frac1(double x);
};

