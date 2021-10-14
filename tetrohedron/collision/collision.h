#pragma once
#include"BVH.h"
#include"../object/cloth.h"
#include"../object/collider.h"
#include"../object/tetrohedron.h"
#include"../thread.h"

class Collision
{
public:
	void initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrohedron>* tetrohedron, Thread* thread);
	
	void buildBVH();
	void searchTriangle(AABB& aabb, int compare_index, std::vector<std::vector<int>>* cloth_neighbor_index, std::vector<std::vector<int>>* collider_neighbor_index);
private:
	std::vector<BVH> cloth_BVH;
	std::vector<BVH> collider_BVH;
	std::vector<BVH> tetrohedron_BVH;

	std::vector<Cloth>* cloth;
	std::vector<Collider>* collider;
	std::vector<Tetrohedron>* tetrohedron;

	Thread* thread;

	void initialBVH(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrohedron>* tetrohedron, Thread* thread);
	
};
