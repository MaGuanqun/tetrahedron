#pragma once
#include"BVH.h"
#include"../object/cloth.h"
#include"../object/collider.h"
#include"../object/tetrahedron.h"
#include"../thread.h"

class Collision
{
public:
	void initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron, Thread* thread);

	void findAllTrianglePairs(int thread_No);
	void findAllNeighborPairs();
	void findPrimitivesAround(int thread_No);
	void collisionDetection(int thread_No);

	void test();
private:
	std::vector<BVH> cloth_BVH;
	std::vector<BVH> collider_BVH;
	std::vector<BVH> tetrahedron_BVH;

	std::vector<Cloth>* cloth;
	std::vector<Collider>* collider;
	std::vector<Tetrahedron>* tetrahedron;

	Thread* thread;
	void buildBVH();
	void initialBVH(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron, Thread* thread);
	void searchTriangle(AABB& aabb, int compare_index, int cloth_No, std::vector<std::vector<int>>* cloth_neighbor_index, std::vector<std::vector<int>>* collider_neighbor_index);
	void findTriangleAroundVertex(int thread_No);
	void findEdgeAroundEdge(int thread_No);
	inline bool vertexInTriangle(int* face_index, int vertex_index);
	inline bool edgeEdgeconnected(int* edge1, int* edge2);
	void pointTriangleCollisionDetection(int thread_No, std::vector<std::vector<std::vector<int>>>& vertex_neighbor_cloth_traingle, std::vector<int>& vertex_index_begin);
	bool checkPointTriangleCollision(double* initial_position, double* current_position, int* triangle_vertex_index, std::vector<std::array<double, 3>>& initial_triangle_position, std::vector<std::array<double, 3>>& current_triangle_position);
	void getAABB();
};
