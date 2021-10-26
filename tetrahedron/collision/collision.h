#pragma once
#include"BVH.h"
#include"../object/cloth.h"
#include"../object/collider.h"
#include"../object/tetrahedron.h"
#include"../thread.h"
#include"predictive_contact.h"

class Collision
{
public:
	void initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron, Thread* thread);

	void findAllTrianglePairs(int thread_No);
	void globalCollision();
	void findPrimitivesAround(int thread_No);
	void collisionDetection(int thread_No);
	void sumTargetPositionPerThread(int thread_id);
	void collisionReDetection(int thread_No);
	void resumTargetPositionPerThread(int thread_id);
	void updateCollisionPosition();
	void test();
	struct TargetPosition
	{
		std::vector<std::vector<std::array<double, 3>>>b_sum; // add on the b vector in projectDynamic.cpp
		std::vector<std::vector<double>>stiffness;//add on the system matrix
		std::vector<double>collision_energy;
		bool** need_update;//to indicate if that item is nonzeor(used)
		void initialSet(int cloth_num)
		{
			b_sum.resize(cloth_num);
			stiffness.resize(cloth_num);
			collision_energy.resize(cloth_num, 0.0);
			need_update = new bool* [cloth_num];
		}

		void initialSet2(int cloth_No, int num)
		{
			b_sum[cloth_No].resize(num, std::array{ 0.0,0.0,0.0 });
			stiffness[cloth_No].resize(num, 0.0);
			need_update[cloth_No] = new bool[num];
		}

		void initial()
		{
			for (int i = 0; i < b_sum.size(); ++i) {
				memset(b_sum[i][0].data(), 0, 24 * b_sum[i].size());
				memset(&stiffness[i][0], 0, 8 * stiffness[i].size());
				memset(need_update[i], 0, b_sum[i].size());
			}
			memset(&collision_energy[0], 0, 8 * collision_energy.size());
		}

		void partialInitial()
		{
			for (int i = 0; i < b_sum.size(); ++i) {
				memset(b_sum[i][0].data(), 0, 24 * b_sum[i].size());
			}
			memset(&collision_energy[0], 0, 8 * collision_energy.size());
		}
	};
	TargetPosition cloth_target_pos;

private:
	int thread_num;
	std::vector<BVH> cloth_BVH;
	std::vector<BVH> collider_BVH;
	std::vector<BVH> tetrahedron_BVH;

	std::vector<Cloth>* cloth;
	std::vector<Collider>* collider;
	std::vector<Tetrahedron>* tetrahedron;

	Thread* thread;
	PredictiveContact predictive_contact;
	void buildBVH();
	void initialBVH(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron, Thread* thread);
	void searchTriangle(AABB& aabb, int compare_index, int cloth_No, std::vector<std::vector<int>>* cloth_neighbor_index, std::vector<std::vector<int>>* collider_neighbor_index);
	void findTriangleAroundVertex(int thread_No);
	void findEdgeAroundEdge(int thread_No);
	inline bool vertexInTriangle(int* face_index, int vertex_index);
	inline bool edgeEdgeconnected(int* edge1, int* edge2);
	void pointSelfTriangleCollisionDetection(int thread_No, int vertex_index, int cloth_No,
		std::vector<std::vector<int>>* vertex_neighbor_triangle, std::vector<std::vector<int>>* collide_vertex_triangle, MeshStruct* vertex_mesh, double radius0, TargetPosition* target_pos);
	bool checkPointTriangleCollision(MeshStruct* vertex_mesh, MeshStruct* triangle_mesh,
		double radius, int vertex_index, int triangle_index, int vertex_cloth_No, int triangle_cloth_No, TargetPosition* target_position, bool new_collision_registration,
		std::vector<double>& vertex_collision_stiffness, std::vector<double>& triangle_collision_stiffness);
	void getAABB();
	std::vector<TargetPosition> cloth_target_pos_per_thread;
	void initialTargetPos(std::vector<Cloth>* cloth,std::vector<Tetrahedron>* tetrahedron, Thread* thread);
	void addTargetPosToSystem(double* b_sum, double& energy, double* current_pos, double* target_pos, double stiffness);
	void pointColliderTriangleCollisionDetection(int thread_No, int vertex_index, int cloth_No,
		std::vector<std::vector<int>>* vertex_neighbor_triangle, std::vector<std::vector<int>>* collide_vertex_triangle, MeshStruct* vertex_mesh, double radius0, TargetPosition* target_pos);
	bool checkPointColliderTriangleCollision(MeshStruct* vertex_mesh, MeshStruct* triangle_mesh,
		double radius, int vertex_index, int triangle_index, int vertex_cloth_No, int triangle_obj_No, TargetPosition* target_position, bool new_collision_registration, std::vector<double>& vertex_collision_stiffness);
	void edgeSelfEdgeCollisionDetection(int thread_No, int edge_index, int cloth_No,
		std::vector<std::vector<int>>* edge_neighbor_edge, std::vector<std::vector<int>>* collide_edge_edge, TriangleMeshStruct* edge_mesh, double radius0, TargetPosition* target_pos);
	bool checkEdgeEdgeCollision(TriangleMeshStruct* edge_mesh, TriangleMeshStruct* compare_mesh,
		double radius, int edge_index, int compare_edge_index, int edge_cloth_No, int compare_cloth_No, TargetPosition* target_position, bool new_collision_registration,
		std::vector<double>& collision_stiffness, std::vector<double>& compare_collision_stiffness);
	void sumTargetPosition();
	void pointSelfTriangleCollisionReDetection(int thread_No, int vertex_index, int cloth_No,std::vector<std::vector<int>>* collide_vertex_triangle,
		MeshStruct* vertex_mesh, double radius0, std::vector<double>* collision_stiffness, TargetPosition* target_postion_);
	void pointColliderTriangleCollisionReDetection(int thread_No, int vertex_index, int cloth_No, std::vector<std::vector<int>>* collide_vertex_triangle,
		MeshStruct* vertex_mesh, double radius0, std::vector<double>* collision_stiffness, TargetPosition* target_postion_);
	void edgeSelfEdgeCollisionReDetection(int thread_No, int edge_index, int cloth_No, std::vector<std::vector<int>>* collide_edge_edge, TriangleMeshStruct* edge_mesh,
		double radius0, std::vector<double>& collision_stiffness, TargetPosition* target_postion_);
	void resumTargetPosition();
};
