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
	void searchTriangle(AABB& aabb, int input_obj_No, int triangle_index, std::vector<int>* obj_neighbor_index, bool is_collider, int thread_No);
private:
	Thread* thread;
	std::vector<Cloth>* cloth;
	std::vector<Collider>* collider;
	std::vector<Tetrahedron>* tetrahedron;

	int tetrahedron_begin_obj_index;

	void initialHashCellLength(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, double& cell_length, double* tolerance_ratio);
	double cell_length;//this is the size of spatial hashing cell, = average edge length
	double thread_num;

	AABB scene_aabb;
	std::vector<std::array<double, 3>>scene_aabb_max_thread;
	std::vector<std::array<double, 3>>scene_aabb_min_thread;
	void getSceneAABB();
	int cell_number[3];
	int cell_num0_cell_num1;
	int hash_max_index[3];
	void setSpatialHashing();
	//std::vector<std::vector<std::array<int, 3>>> spatial_hashing_triangle_per_thread;
	//std::vector<std::array<int, 3>> spatial_hashing_triangle;

	std::vector<std::vector<int>> spatial_hashing_value_per_thread;
	std::vector<int> spatial_hashing_value;
	std::vector<std::vector<int>>spatial_hashing_obj_index_per_thread;
	std::vector<int> spatial_hashing_obj_index;
	std::vector<std::vector<int>>spatial_hashing_triangle_index_per_thread;
	std::vector<int> spatial_hashing_triangle_index;


	void triangleHashingValue(AABB& aabb, std::vector<int>* spatial_hashing_triangle_index, std::vector<int>* spatial_hashing_value, int triangle_index, std::vector<int>* hash_value);
	void obtainTriangleHashingValue(AABB& aabb, std::vector<int>* hash_value);
	RadixSort radix_sort;
	void testRadixSort();


	std::vector<std::vector<std::vector<int>>> obj_triangle_hash;
	void initialTriangleHash();
	void setHashTogether();
	std::vector<int>hash_value_begin;
	std::vector<int> prifix_sum;
	void setPrifixSum();
	bool*** obj_is_used;
};

