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
	void setInObject(std::vector<Cloth>* cloth, std::vector<Collider>* collider,std::vector<Tetrahedron>* tetrahedron, Thread* thread);
	void getSceneAABB(int thread_No);
	void triangleHashing(int thread_No);
	void culling();
private:
	Thread* thread;
	std::vector<Cloth>* cloth;
	std::vector<Collider>* collider;
	std::vector<Tetrahedron>* tetrahedron;
	void initialHashCellLength(std::vector<Cloth>* cloth, double& cell_length);
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
	void initialSpatialHashingList();
	std::vector<std::vector<std::vector<std::array<int, 2>>>> spatial_hashing_cloth_triangle;
	std::vector<std::vector<std::vector<std::array<int, 2>>>> spatial_hashing_collider_triangle;
	void triangleHashingValue(AABB& aabb, std::vector<std::array<int, 2>>* spatial_hashing_triangle, int triangle_index);
	RadixSort radix_sort;
};

