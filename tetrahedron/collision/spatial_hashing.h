//Here we list the <element,spatial hashing index> pairs, avoid creating cells.
//use radix sorting to sort pairs by index -> for searching
#pragma once
#include"../thread.h"
#include"../object/cloth.h"
#include"../object/collider.h"
#include"../object/tetrahedron.h"
#include"parallel_radix_sort.h"
#include<atomic>
#include<unordered_set>

class SpatialHashing
{
public:
	void setInObject(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
		Thread* thread, double* tolerance_ratio, unsigned int max_cell_count, 
		unsigned int max_index_number_in_one_cell,
		unsigned int max_index_number_in_one_cell_collider, unsigned int estimate_coeff_for_vt_pair_num, unsigned int estimate_coeff_for_vt_collider_pair_num,
		unsigned int estimate_coeff_for_tv_pair_num, unsigned int estimate_coeff_for_ee_pair_num, bool record_pair_by_element,
		unsigned int estimate_coeff_for_tv_collider_pair_num, unsigned int estimate_coeff_for_ee_collider_pair_num);
	void getSceneAABB(int thread_No);

	void buildSpatialHashing(double* scene_aabb);


	void triangleHashingSmallerHashTable(int thread_No);

	std::vector<unsigned int> cell_begin_per_thread;


	void findAllPairsByPrimitive(int thread_No);

	void findAllPairsHashTable(int thread_No);

	unsigned int** vertex_triangle_pair; //except collider, inner vector store vertex_1 index, obj_1_index, tri_2_index, obj_2_index
	unsigned int** vertex_obj_triangle_collider_pair; //inner vector store vertex_1 index, obj_1_index, tri_2_index, obj_2_index
	unsigned int** vertex_collider_triangle_obj_pair; //inner vector store vertex_1 index, obj_1_index, tri_2_index, obj_2_index
	unsigned int** edge_edge_pair; 
	unsigned int** edge_edge_pair_collider; //

	unsigned int** vertex_triangle_pair_by_vertex;//store pair by every vertex. (obj_index,triangle_index)
	unsigned int** vertex_triangle_pair_num_record;//record the number of triangle pairs for every vertex. For fast initialize, we recoed it in this variable
	//unsigned int** triangle_vertex_pair_by_triangle; //store pair by triangle. for every triangle, store vertex: (obj_index,vertex_index)
	//unsigned int** triangle_vertex_pair_num_record;// record the number of vertex pairs for every triangle. For fast initialize, we recoed it in this variable
	unsigned int** edge_edge_pair_by_edge; //for coordinate descent, we should store both (e1,e2) & (e2,e1)
	unsigned int** edge_edge_pair_num_record;//record the number of triangle pairs for every vertex. For fast initialize, we recoed it in this variable
	unsigned int** vertex_obj_triangle_collider_pair_by_vertex;// store pair by every vertex. (obj_index, triangle_index)
	unsigned int** vertex_obj_triangle_collider_num_record;//record the number of triangle pairs for every vertex. For fast initialize, we record it in this variable

	unsigned int** triangle_obj_vertex_collider_pair_by_triangle;// store pair by every triangle. (obj_index, triangle_index)
	unsigned int** triangle_obj_vertex_collider_num_record;//record the number of vertex pairs for every triangle. For fast initialize, we record it in this variable
	unsigned int** edge_obj_edge_collider_pair_by_edge; //for coordinate descent, we should store both (e1,e2) & (e2,e1)
	unsigned int** edge_obj_edge_collider_num_record;//record the number of triangle pairs for every vertex. For fast initialize, we recoed it in this variable


	bool** is_used_vertex_triangle_pair_by_vertex;
	bool** is_used_edge_edge_pair_by_edge;
	bool** is_used_vertex_obj_triangle_collider_pair_by_vertex;
	bool** is_used_triangle_obj_vertex_collider_pair_by_triangle;
	bool** is_used_edge_obj_edge_collider_pair_by_edge;





	std::vector<unsigned int>hash_value_for_test;
	double scene_aabb[6];

	double cell_length;//this is the size of spatial hashing cell, = average edge length




	unsigned int* spatial_hashing_cell_triangle;//in every vector, we store triangle index,obj_index respectively
	unsigned int* spatial_hashing_cell_edge;//in every vector, we store triangle index,obj_index respectively
	unsigned int* spatial_hashing_cell_vertex;//in every vector, we store triangle index,obj_index respectively
	//unsigned int*** spatial_hashing_actual_hash_value_for_test;//in every vector, we store triangle index,obj_index respectively
	unsigned int* spatial_hashing_cell_collider_triangle;//in every vector, we store triangle index, obj_index respectively
	unsigned int* spatial_hashing_cell_collider_edge;//in every vector, we store triangle index, obj_index respectively
	unsigned int* spatial_hashing_cell_collider_vertex;//in every vector, we store triangle index, obj_index respectively
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

	//void combineHashTable(int thread_No);

	void testColliderPair();

private:

	std::unordered_set<unsigned int>* unduplicate_cell_index_for_element;


	unsigned int total_length_every_element_for_vertex_triange;
	unsigned int total_length_every_element_for_vertex_triange_collider;
	unsigned int total_length_every_element_for_triangle_vertex;
	unsigned int total_length_every_element_for_edge_edge;

	unsigned int total_length_every_element_for_triange_vertex_collider;
	unsigned int total_length_every_element_for_edge_edge_collider;


	unsigned int total_length_every_element_for_vertex_triange_exist;
	unsigned int total_length_every_element_for_vertex_triange_collider_exist;
	unsigned int total_length_every_element_for_triangle_vertex_exist;
	unsigned int total_length_every_element_for_edge_edge_exist;
	unsigned int total_length_every_element_for_triange_vertex_collider_exist;
	unsigned int total_length_every_element_for_edge_edge_collider_exist;

	bool record_pair_by_element;//if true, record pairs by vertex, triangle, or edge

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
	std::vector<unsigned int*> obj_vertex_index_begin_per_thread;// for tet, size is surface vertex size
	std::vector<unsigned int*> obj_edge_index_begin_per_thread;


	std::vector<unsigned int*> collider_vertex_index_begin_per_thread;
	std::vector<unsigned int*> collider_edge_index_begin_per_thread;

	std::vector<unsigned int*> tetrahedron_vertex_index_on_surface;


	std::vector<std::array<int, 3>*> triangle_vertex_index;
	std::vector<std::array<int, 3>*> triangle_vertex_index_collider;

	std::vector<unsigned int*> edge_vertex_index;
	std::vector<unsigned int*> edge_vertex_index_collider;



	void initialHashCell(unsigned int total_triangle_num, unsigned int max_index_number_in_one_cell,
		unsigned int max_index_number_in_one_cell_collider, unsigned int estimate_coeff_for_vt_pair_num, unsigned int estimate_coeff_for_vt_collider_pair_num,
		unsigned int estimate_coeff_for_tv_pair_num, unsigned int estimate_coeff_for_ee_pair_num,
		unsigned int estimate_coeff_for_tv_collider_pair_num, unsigned int estimate_coeff_for_ee_collider_pair_num);


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
		uint64_t P1, uint64_t P2, uint64_t P3, std::atomic_uint* spatial_hashing_cell_triangle_size, unsigned int max_index_number_in_one_cell,
		std::unordered_set<unsigned int>* record_hash_index);



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
		unsigned int* spatial_hashing_cell_, std::atomic_uint* spatial_hash_size,
		unsigned int* spatial_hashing_cell_collider_, std::atomic_uint* spatial_hash_size_collider, bool is_edge,
		std::vector<unsigned int>* neighbor_primitive_0, std::vector<unsigned int>* neighbor_primitive_1,
		unsigned int max_index_number_in_one_cell, unsigned int max_index_number_in_one_cell_collider);

	void searchPrimitiveOnCollider(double* aabb, unsigned int* hash_index, unsigned int hash_cell_size,
		unsigned int input_obj_No, unsigned int vertex_index,
		unsigned int*& triangle_pair_with_collider_, unsigned int thread_No,
		unsigned int* spatial_hashing_cell_, unsigned int* spatial_hash_size_, unsigned int max_index_number_in_one_cell_collider);

	void reorganzieDataOfObjects();





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


	std::atomic_uint* spatial_hashing_cell_triangle_size;
	std::atomic_uint* spatial_hashing_cell_edge_size;
	std::atomic_uint* spatial_hashing_cell_vertex_size;
	std::atomic_uint* spatial_hashing_cell_collider_triangle_size;
	std::atomic_uint* spatial_hashing_cell_collider_edge_size;
	std::atomic_uint* spatial_hashing_cell_collider_vertex_size;

	void findAllVertexTrianglePairs(int thread_No);

	void findAllEdgeEdgePairs(int thread_No);

	std::vector<std::vector<unsigned int>> primitive_index_record;
	unsigned int max_triangle_cell_size;
	unsigned int max_vertex_cell_size;


	std::vector<MeshStruct::Vertex*> vertices;
	std::vector<MeshStruct::Edge*> edges;

	std::vector<MeshStruct*>mesh_struct;

	
	void initialTriangleHashValue(double* aabb,
		std::vector<unsigned int>* spatial_hashing_cell, double* scene_aabb, double cell_length);
	void setSpatialHashingInitialCount();

	void 	collectAllVoxel(std::vector<unsigned int>& select_hash_index, std::vector<unsigned int>&
		initial_hash_index, double start__[3], double end__[3],
		double* scene_aabb);
	bool intersectSpatialHasingCube(double* AABB, double* start, double* direction, double* end);

	inline double frac0(double x);
	inline double frac1(double x);

	bool searchPairByCell = false;//this means we need to do spatial hashing for vertex, triangle and edge, and find all pairs in a cell, may have duplicate pairs

	void findAllVertexTrianglePairsByPrimitive(int thread_No);
	void findAllEdgeEdgePairsByPrimitive(int thread_No);

	void findAllVertexTrianglePairsByPrimitiveByVertex(int thread_No);
	void findAllEdgeEdgePairsByPrimitiveByEdge(int thread_No);
	void findAllTriangleVertexPairsByPrimitiveByTriangle(int thread_No);

	void findAllVertexTrianglePairsByPrimitiveSingleObj(int thread_No, int obj_No, unsigned int* &primitive_pair_,
		unsigned int vertex_start, unsigned int vertex_end, std::array<double, 6>* vertex_aabb, 
		unsigned int max_index_number_in_one_cell_triangle_, unsigned int* spatial_hashing_cell_triangle, std::atomic_uint* spatial_hashing_cell_triangle_size,
		std::vector<std::array<double, 6>*>& obj_tri_aabb, bool is_self, bool is_tet);

	void findAllVertexTrianglePairsByPrimitiveSingleObj_ByVertex(int thread_No, int obj_No, unsigned int* primitive_pair,
		unsigned int* primitive_pair_num_record,
		unsigned int vertex_start, unsigned int vertex_end, std::array<double, 6>* vertex_aabb,
		unsigned int max_index_number_in_one_cell_triangle_, unsigned int* spatial_hashing_cell_triangle, std::atomic_uint* spatial_hashing_cell_triangle_size,
		std::vector<std::array<double, 6>*>& obj_tri_aabb, bool is_self, bool is_tet, bool* exist_flag); //this is for vertex_triangle_pair_by_vertex etc

	void findAllTriangleVertexColliderPairsByPrimitiveSingleObj_ByTriangle(int thread_No, int obj_No, unsigned int* primitive_pair,
		unsigned int* primitive_pair_num_record,
		std::array<double, 6>* triangle_aabb,
		unsigned int max_index_number_in_one_cell_vertex_, unsigned int* spatial_hashing_cell_vertex, std::atomic_uint* spatial_hashing_cell_vertex_size,
		std::vector<std::array<double, 6>*>& obj_vertex_aabb, bool* exist_flag);

	void triangleHashValue(double* aabb,
		std::vector<unsigned int>* spatial_hashing_index, double* scene_aabb, double cell_length,
		unsigned int hash_cell_count, uint64_t P1, uint64_t P2, uint64_t P3);

	std::vector<unsigned int>hash_count_start_per_thread;
	void findAllEdgeEdgePairsByPrimitiveSingleObjByEdge(int thread_No, int obj_No, unsigned int* primitive_pair, unsigned int* primitive_pair_num_record_,
		std::vector<std::array<double, 6>*>& obj_edge_aabb_,
		unsigned int max_index_number_in_one_cell_edge_, unsigned int* spatial_hashing_cell_edge, std::atomic_uint* spatial_hashing_cell_edge_size,
		bool is_self, unsigned int total_length_every_element_for_edge_edge, bool* exist_flag);
	void findAllEdgeEdgePairsByPrimitiveSingleObj(int thread_No, int obj_No, unsigned int*& primitive_pair_,
		std::vector<std::array<double, 6>*>& obj_edge_aabb_,
		unsigned int max_index_number_in_one_cell_edge_, unsigned int* spatial_hashing_cell_edge, std::atomic_uint* spatial_hashing_cell_edge_size,
		bool is_self);
	std::vector<unsigned int*>representative_edge_num;
	std::vector<unsigned int*>representative_edge_num_collider;
	std::vector<unsigned int*>face_edge;
	std::vector<unsigned int*>face_edge_collider;

	void findAllTriangleVertexPairByTriangleSingleObj(int obj_No, unsigned int* vt_pair_initial, unsigned int* vt_pair_num, unsigned int** tv_pair,
		unsigned int** tv_pair_num, unsigned int total_length_every_element_vt, unsigned int total_length_every_element_tv);

	//void findAllTriangleVertexPairByTriangle();

	//initial triangle_vertex_pair_num_record
	//void initialPairByElement();

	void initialExistFlag(bool* flag, unsigned int num);

};

