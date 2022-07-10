#pragma once
#include"BVH.h"
#include"../object/cloth.h"
#include"../object/collider.h"
#include"../object/tetrahedron.h"
#include"../thread.h"
#include"predictive_contact.h"
#include"spatial_hashing.h"
#include"ApproxCCD.h"
#include"../external/Eigen/Dense"
#include"collision_constraint.h"
#include"CCD.h"
//#include"drawCulling.h"
#include"primitive_distance.h"
#include"DCD.h"
#include"../basic/floor.h"



//#include"mesh_patch.h"

using namespace Eigen;

class Collision
{
public:
	double collision_time;

	void initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron, Thread* thread, Floor* floor,
		double* tolerance_ratio, unsigned int use_method);
	void initialDHatTolerance(double ave_edge_length);
	void findAllTrianglePairs(int thread_No);
	void globalCollision();
	void findPrimitivesAround(int thread_No);
	void collisionDetection(int thread_No);
	void sumTargetPositionPerThread(int thread_id);
	void collisionReDetection(int thread_No);
	void resumTargetPositionPerThread(int thread_id);
	void updateCollisionPosition();
	void collisionTime(int thread_No);
	void collisionConstraint(int thread_No);
	void re_collisionConstraint(int thread_No);
	void solveCollisionConstraint();
	void test();
	void collisionCulling();
	void globalCollisionTime();
	//void findAllPatchPairs(int thread_No);

	//std::vector<std::vector<double>> target_position_and_stiffness; //thus, the actual number is 4*target_position_index[0]
	//std::vector<std::vector<unsigned int>>target_position_index; //the first element store the number of primitives. thus, the actual number from [1] is 2*[0]
	std::vector<std::vector<unsigned int>>point_triangle_target_pos_index; 
	std::vector<std::vector<unsigned int>>point_triangle_collider_target_pos_index;
	//std::vector<std::vector<unsigned int>>point_collider_triangle_target_pos_index;
	std::vector<std::vector<unsigned int>>edge_edge_target_pos_index;
	//std::vector<std::vector<unsigned int>>edge_edge_collider_target_pos_index;


	unsigned int ave_pair_num[5];//vertex_triangle_pair,edge_edge_pair,vertex_obj_triangle_collider_pair,vertex_collider_triangle_obj_pair,edge_edge_pair_collider.


	struct TargetPosition
	{
		std::vector<std::vector<std::array<double, 3>>>b_sum; // add on the b vector in projectDynamic.cpp
		std::vector<std::vector<double>>stiffness;//add on the system matrix
		double collision_energy;
		bool** need_update;//to indicate if that item is nonzero(used)
		void initialSet(int cloth_num)
		{
			b_sum.resize(cloth_num);
			stiffness.resize(cloth_num);
			collision_energy = 0.0;
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
			collision_energy = 0.0;
		}

		void partialInitial()
		{
			for (int i = 0; i < b_sum.size(); ++i) {
				memset(b_sum[i][0].data(), 0, 24 * b_sum[i].size());
			}
			collision_energy = 0.0;
		}
	};
	TargetPosition obj_target_pos;

	size_t* time_stamp;
	//void drawMeshPatch(Camera* camera);

	void getSceneAABB();

	void findPointTriangleEdgeEdgePair(int thread_No);

	void getAABB();

	void findAllVertexVertexEdgePairs(int thread_No);

	void collisionTimeCompare(int thread_No);

	void collisionEnergy(int thread_No);
	void collisionEnergy();

	void solveCollisionConstraintDCD();
	void reSolveCollisionConstraintDCD();

	//DrawCulling* draw_culling;

	void getCollisionPair();
	void getCollisionPair(int thread_No);


	unsigned int collisionConstraintNumber(unsigned int* point_triangle_collider_constraint, unsigned int* point_triangle_constraint, unsigned int* edge_edge_constraint);

	void XPBDsolveCollisionConstraint();

	void re_XPBDsolveCollisionConstraint();


	//void setParameter(std::vector<double>* lambda, double* floor_lambda, std::vector<unsigned int>* collision_lambda_index_start, double damp_stiffness,double dt);



	SpatialHashing spatial_hashing;
	void buildSpatialHashingForOri();

	double* energy;

	double* friction_coe;

private:

	

	//double* floor_lambda;

	unsigned int** vertex_edge_pair;
	unsigned int** vertex_obj_edge_collider_pair;
	unsigned int** vertex_collider_edge_obj_pair;

	unsigned int** vertex_vertex_pair;
	unsigned int** vertex_vertex_pair_collider;


	unsigned int max_index_number_in_one_cell;
	unsigned int max_index_number_in_one_cell_collider;


	unsigned int estimate_coeff_for_pair_num;

	unsigned int tetrahedron_begin_obj_index;
	unsigned int total_obj_num;//except collider
	unsigned int total_obj_with_collider;// include collider

	unsigned int thread_num;
	std::vector<BVH> obj_BVH; //cloth + tetrahedron
	std::vector<BVH> collider_BVH;

	std::vector<Cloth>* cloth;
	std::vector<Collider>* collider;
	std::vector<Tetrahedron>* tetrahedron;

	Thread* thread;
	PredictiveContact predictive_contact;


	bool use_BVH;

	std::vector<unsigned int> vertex_triangle_pair_index_start_per_thread;//thread_No, index, respectively
	std::vector<unsigned int> vertex_obj_triangle_collider_pair_index_start_per_thread;//thread_No, index, respectively
	std::vector<unsigned int> vertex_collider_triangle_obj_pair_index_start_per_thread;//thread_No, index, respectively
	std::vector<unsigned int> edge_edge_pair_index_start_per_thread;//thread_No, index, respectively
	std::vector<unsigned int> edge_edge_pair_collider_index_start_per_thread;//thread_No, index, respectively


	std::vector<unsigned int> vertex_edge_pair_index_start_per_thread;//thread_No, index, respectively
	std::vector<unsigned int> vertex_obj_edge_collider_pair_index_start_per_thread;//thread_No, index, respectively
	std::vector<unsigned int> vertex_collider_edge_obj_pair_index_start_per_thread;//thread_No, index, respectively
	std::vector<unsigned int> vertex_vertex_pair_index_start_per_thread;//thread_No, index, respectively
	std::vector<unsigned int> vertex_vertex_pair_collider_index_start_per_thread;//thread_No, index, respectively

	//std::vector<unsigned int> target_position_element_start_per_thread;
	std::vector<unsigned int> point_triangle_target_position_element_start_per_thread;
	std::vector<unsigned int> point_collider_triangle_target_position_element_start_per_thread;
	std::vector<unsigned int> edge_edge_target_position_element_start_per_thread;
	std::vector<unsigned int> edge_edge_collider_target_position_element_start_per_thread;

	std::vector<unsigned int> point_triangle_collider_target_position_element_start_per_thread;
	//std::vector<unsigned int> target_position_start_per_thread;


	void setPairIndexEveryThread(unsigned int** pair, std::vector<unsigned int>& pair_index_start_per_thread,
		unsigned int& ave_pair_num);
	void setPairIndexEveryThread();
	void setVertexVertexEdgePairIndexEveryThread();

	std::vector<double>collision_time_thread;

	ApproxCCD approx_CCD;
	CollisionConstraint collision_constraint;

	double d_hat;
	double tolerance_2;// distance to report a collision

	double epsilon;


	void buildBVH();
	void initialBVH(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron, Thread* thread);
	void searchTriangle(double* aabb, unsigned int compare_index, unsigned int obj_No, std::vector<std::vector<unsigned int>>* obj_neighbor_index,
		std::vector<std::vector<unsigned int>>* collider_neighbor_index, bool is_collider);
	void findObjTriangleAroundVertex(int thread_No);
	void findColliderTriangleAroundVertex(int thread_No);
	void findVertexAroundColliderTriangle(int thread_No);
	void findEdgeAroundEdge(int thread_No);
	//inline bool vertexInTriangle(int* face_index, int vertex_index);
	//inline bool edgeEdgeconnected(int* edge1, int* edge2);
	void pointSelfTriangleCollisionDetection(int thread_No, int vertex_index, int cloth_No,
		std::vector<int>* vertex_neighbor_triangle, std::vector<int>* collide_vertex_triangle, MeshStruct* vertex_mesh,
		double radius0, TargetPosition* target_pos, double vertex_collision_stiffness);
	bool checkPointTriangleCollision(MeshStruct* vertex_mesh, MeshStruct* triangle_mesh,
		double radius, int vertex_index, int triangle_index, int vertex_cloth_No, int triangle_cloth_No, TargetPosition* target_position, bool new_collision_registration,
		double vertex_collision_stiffness, double triangle_collision_stiffness);
	
	std::vector<TargetPosition> obj_target_pos_per_thread;
	void initialTargetPos(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, Thread* thread);
	void addTargetPosToSystem(double* b_sum, double& energy, double* current_pos, double* target_pos, double stiffness);
	void pointColliderTriangleCollisionDetection(int thread_No, int vertex_index, int cloth_No,
		std::vector<int>* vertex_neighbor_triangle, std::vector<int>* collide_vertex_triangle, MeshStruct* vertex_mesh, double radius0, TargetPosition* target_pos,
		double vertex_collision_stiffness);
	bool checkPointColliderTriangleCollision(MeshStruct* vertex_mesh, MeshStruct* triangle_mesh,
		double radius, int vertex_index, int triangle_index, int vertex_cloth_No, int triangle_obj_No,
		TargetPosition* target_position, bool new_collision_registration, double vertex_collision_stiffness);
	void edgeSelfEdgeCollisionDetection(int thread_No, int edge_index, int cloth_No,
		std::vector<int>* edge_neighbor_edge, std::vector<int>* collide_edge_edge, MeshStruct* edge_mesh, double radius0,
		TargetPosition* target_pos, double stiffness_0);
	bool checkEdgeEdgeCollision(MeshStruct* edge_mesh, MeshStruct* compare_mesh,
		double radius, int edge_index, int compare_edge_index, int edge_cloth_No, int compare_cloth_No, TargetPosition* target_position, bool new_collision_registration,
		double collision_stiffness, double compare_collision_stiffness);
	void sumTargetPosition();
	void pointSelfTriangleCollisionReDetection(int thread_No, int vertex_index, int cloth_No, std::vector<int>* collide_vertex_triangle,
		MeshStruct* vertex_mesh, double radius0, double collision_stiffness, TargetPosition* target_postion_);
	void pointColliderTriangleCollisionReDetection(int thread_No, int vertex_index, int cloth_No, std::vector<std::vector<int>>* collide_vertex_triangle,
		MeshStruct* vertex_mesh, double radius0, std::vector<double>* collision_stiffness, TargetPosition* target_postion_);
	void edgeSelfEdgeCollisionReDetection(int thread_No, int edge_index, int cloth_No, std::vector<std::vector<int>>* collide_edge_edge, TriangleMeshStruct* edge_mesh,
		double radius0, std::vector<double>& collision_stiffness, TargetPosition* target_postion_);
	void resumTargetPosition();
	void initialSpatialHashing(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
		Thread* thread, double* tolerance_ratio);
	void initialNeighborPrimitive();
	void colliderTriangleVertexCollisionDetection(int thread_No, int triangle_index, int collider_No,
		std::vector<std::vector<int>>* triangle_neighbor_vertex, std::vector<std::vector<int>>* collide_triangle_vertex, MeshStruct* triangle_mesh, double radius0,
		TargetPosition* target_pos);
	void colliderTriangleVertexCollisionReDetection(int thread_No, int triangle_index, int collider_No,
		std::vector<int>* collide_triangle_vertex, MeshStruct* triangle_mesh, double radius0, TargetPosition* target_pos);
	void testCollision();


	void pointSelfTriangleCollisionTime(double* collision_time, std::vector<int>* vertex_neighbor_triangle, double* initial_vertex_pos, double* current_vertex_pos, int vertex_index);
	void edgeEdgeCollisionTime(double* collision_time, std::vector<int>* edge_neighbor_edge, double* initial_edge_vertex_0, double* initial_edge_vertex_1, double* current_edge_vertex_0, double* current_edge_vertex_1);
	void pointColliderTriangleCollisionTime(double* collision_time, int* triangle_vertex_index, std::vector<int>* triangle_neighbor_vertex,
		std::array<double, 3>* initial_position, std::array<double, 3>* current_position,
		double* initial_ori_face_normal, double* current_ori_face_normal, double* cross_for_CCD,
		floating* f_initial_normal, floating* f_current_normal, floating* f_cross_for_CCD, int triangle_index);
	void pointSelfTriangleClose(std::vector<int>* vertex_neighbor_triangle, double* initial_vertex_pos, double* current_vertex_pos,
		int vertex_index, int cloth_No, double mass, TargetPosition* target_position);
	void pointColliderTriangleClose(int* triangle_vertex_index, std::vector<int>* triangle_neighbor_vertex,
		std::vector<double*>& current_position, double* current_face_normal, TargetPosition* target_position);
	void addTargetPosToSystemTotal(double* b_sum, double& energy, double* current_pos, double* target_pos, double stiffness, double& sum_stiffness, bool& update);
	void edgeEdgeClose(std::vector<int>* edge_neighbor_edge, double* initial_edge_vertex_0, double* initial_edge_vertex_1, double* current_edge_vertex_0, double* current_edge_vertex_1,
		int cloth_No, int edge_vertex_index_0, int edge_vertex_index_1, double* mass, TargetPosition* target_position);
	void testCulling();
	void testIfBuildCollisionConstraint();
	double tolerance;
	double* tolerance_ratio;
	double tolerance_radius[4];
	double eta;//for setting gap in ccd
	void testTwoVectorsAreSame(std::vector<std::vector<int>>& vec1, std::vector<std::vector<int>>& vec2, unsigned int obj_index,
		unsigned int triangle_index);

	void findPatchOfObjects();

	//MeshPatch mesh_patch;

	void getAABBWithoutTolerance();

	//void testBVHUpdate();
	void searchPatch(double* aabb, unsigned int compare_index, unsigned int obj_No, bool& intersect);

	double scene_aabb[6];


	void initialCollidePairInfo();

	//unsigned int** point_triangle_pair; //except collider, inner vector store vertex_1 index, obj_1_index, tri_2_index, obj_2_index
	//unsigned int** point_collider_triangle_obj_pair;  //, inner vector store vertex_1 index, obj_1_index, tri_2_index, obj_2_index
	//unsigned int** point_obj_triangle_collider_pair;
	//unsigned int** edge_edge_pair;  //except collider, inner vector store edge_1 index, obj_1_index, edge_2_index, obj_2_index
	//unsigned int** edge_edge_collider_pair;

	//reorganize the data of different objects
	std::vector<std::array<double, 6>*> obj_tri_aabb;
	std::vector<std::array<double, 6>*> vertex_aabb;
	std::vector<std::array<double, 6>*> edge_aabb;

	std::vector<std::array<double, 6>*> obj_tri_aabb_collider;
	std::vector<std::array<double, 6>*> vertex_aabb_collider;
	std::vector<std::array<double, 6>*> edge_aabb_collider;


	std::vector<std::array<double, 3>*> vertex_for_render;
	std::vector<std::array<double, 3>*> vertex_position;
	std::vector<std::array<double, 3>*> triangle_normal_render;
	std::vector<std::array<double, 3>*> triangle_normal;

	std::vector<std::array<double, 3>*> triangle_normal_not_normalized;
	std::vector<std::array<double, 3>*> triangle_normal_render_not_normalized;
	std::vector<std::array<double, 3>*> cross_for_approx_CCD;

	std::vector<std::array<int, 3>*> triangle_indices;
	std::vector<std::array<int, 3>*> triangle_indices_collider;


	std::vector<std::array<double, 3>*> vertex_for_render_collider;
	std::vector<std::array<double, 3>*> vertex_position_collider;
	std::vector<std::array<double, 3>*> triangle_normal_render_collider;
	std::vector<std::array<double, 3>*> triangle_normal_collider;

	std::vector<std::array<double, 3>*> triangle_normal_not_normalized_collider;
	std::vector<std::array<double, 3>*> triangle_normal_render_not_normalized_collider;
	std::vector<std::array<double, 3>*> cross_for_approx_CCD_collider;


	std::vector<unsigned int*> representative_vertex_num;
	std::vector<unsigned int*> representative_edge_num;

	std::vector<unsigned int*> representative_vertex_num_collider;
	std::vector<unsigned int*> representative_edge_num_collider;

	std::vector<std::array<int, 3>*> triangle_index_in_order;
	std::vector<std::array<int, 3>*> triangle_index_in_order_collider;
	std::vector<MeshStruct::Face*> faces;
	std::vector<MeshStruct::Edge*> edges;

	std::vector<unsigned int*> face_edges;
	std::vector<unsigned int*> collider_face_edges;

	std::vector<unsigned int*> edge_vertices;
	std::vector<unsigned int*> collider_edge_vertices;


	std::vector<double*>triangle_normal_magnitude_reciprocal;
	std::vector<double*>collider_triangle_normal_magnitude_reciprocal;


	std::vector<unsigned int*> vertex_index_start_per_thread;


	std::vector<double*>mass;
	std::vector<double*>mass_inv;

	void reorganzieDataOfObjects();
	bool has_collider;

	void testIfSPRight();
	void findInBVH(unsigned int obj_0, unsigned int tri_index_0, unsigned int obj_1, unsigned int tri_index1, bool with_collider);
	void findInSP(std::vector<std::vector<std::vector<unsigned int>>>* tri_obj, unsigned int obj_index);
	void findInSPCollider(std::vector<std::vector<std::vector<unsigned int>>>* tri_obj, unsigned int obj_index);
	void testRepeatability();

	struct TriangleElementPair
	{
		unsigned int index[4];

		bool operator<(const TriangleElementPair& t1) const
		{
			if (index[1] < t1.index[1])
				return true;
			else if (index[1] == t1.index[1]) {
				if (index[0] < t1.index[0])
					return true;
				else if (index[0] == t1.index[0]) {
					if (index[3] < t1.index[3]) {
						return true;
					}
					else if (index[3] == t1.index[3]) {
						if (index[2] < t1.index[2]) {
							return true;
						}
					}
					return false;
				}
			}
			return false;
		}

		bool operator==(const TriangleElementPair& rhs)
		{
			for (unsigned int i = 0; i < 4; ++i) {
				if (index[i] != rhs.index[i]) {
					return false;
				}
			}
			return true;
		}
	};

	std::vector<time_t> record_time;
	std::vector<time_t> record_time1;

	std::vector<unsigned int> edge_edge_count;
	std::vector<unsigned int> vertex_triangle_count;

	void totalCount();

	void findVertexTriangleInBVH(unsigned int obj_0, unsigned int vertex_index_0, unsigned int obj_1, unsigned int tri_index1);
	void findEdgeEdgeInBVH(unsigned int obj_0, unsigned int edge_index_0, unsigned int obj_1, unsigned int edge_index_1);

	void vertexTriangleCollisionTime(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index, unsigned int end_pair_index,
		double& collision_time);
	void edgeEdgeCollisionTime(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, double& collision_time);
	void vertexColliderTriangleCollisionTime(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, double& collision_time);
	void vertexTriangleColliderCollisionTime(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, double& collision_time);
	void edgeEdgeColliderCollisionTime(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, double& collision_time);
	double conservative_rescaling;

	void initialPair();
	void findAllVertexEdgePairs(int thread_No);
	void findAllVertexVertexPairs(int thread_No);

	void findVertexEdgePair(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, CollisionIndicateType type);

	void findVertexVertexPair(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, CollisionIndicateType type);

	void vertexVertexCollisionTime(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, double& collision_time, CollisionIndicateType type);
	void vertexEdgeCollisionTime(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, double& collision_time, CollisionIndicateType type);

	void vertexTriangleCollisionTimeCompare(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, double& collision_time);

	std::vector<std::vector<Eigen::Vector3d>> vertex_for_render_eigen;
	std::vector<std::vector<Eigen::Vector3d>> vertex_position_eigen;

	void edgeEdgeCollisionTimeCompare(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, double& collision_time);

	void vertexVertexCollisionTimeCompare(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, double& collision_time, CollisionIndicateType type);

	void vertexEdgeCollisionTimeCompare(int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, double& collision_time, CollisionIndicateType type);

	void updateEigenPosition();
	void pointTriangleResponse(int thread_No, TargetPosition* target_pos);
	void pointTriangleColliderResponse(int thread_No, TargetPosition* target_pos);
	void pointColliderTriangleResponse(int thread_No, TargetPosition* target_pos);
	void edgeEdgeResponse(int thread_No, TargetPosition* target_pos);
	void edgeEdgeColliderResponse(int thread_No, TargetPosition* target_pos);

	void pointTriangleResponse(unsigned int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, TargetPosition* target_pos);
	void pointTriangleColliderResponse(unsigned int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, TargetPosition* target_pos);
	void pointColliderTriangleResponse(unsigned int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, TargetPosition* target_pos);
	void edgeEdgeResponse(unsigned int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, TargetPosition* target_pos);
	void edgeEdgeColliderResponse(unsigned int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, TargetPosition* target_pos);

	void checkTargetPosSize(int thread_No);


	std::vector<unsigned int> test_point_triangle_record_true_number;
	std::vector<unsigned int> test_edge_edge_record_true_number;

	

	void setTargetPositionEven();

	void setPairIndexEveryThread(std::vector<unsigned int>* pair, std::vector<unsigned int>& pair_index_start_per_thread);

	void computeEnergy(unsigned int pair_thread_No, unsigned int index_start, unsigned int index_end, double& energy);

	double collision_stiffness[4];//we assume the initial stiffness of every objects are same

	DCD dcd;

	void re_pointColliderTriangleResponse(unsigned int pair_thread_No, unsigned int index_start, unsigned int index_end, TargetPosition* target_pos);

	void re_pointTriangleResponse(unsigned int pair_thread_No, unsigned int index_start, unsigned int index_end, TargetPosition* target_pos);
	void setIndexEveryThread(std::vector<unsigned int>* pair, std::vector<unsigned int>& pair_index_start_per_thread);
	void re_pointTriangleColliderResponse(unsigned int pair_thread_No, unsigned int index_start, unsigned int index_end, TargetPosition* target_pos);

	void construct_b_sum(double* b_sum, double* target_pos, double stiffness);

	void re_pointTriangleResponse(int thread_No, TargetPosition* target_pos);
	void re_pointTriangleColliderResponse(int thread_No, TargetPosition* target_pos);
	void re_pointColliderTriangleResponse(int thread_No, TargetPosition* target_pos);

	void testPairEven();

	void re_edgeEdgeResponse(unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, TargetPosition* target_pos);
	void re_edgeEdgeColliderResponse(unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index, TargetPosition* target_pos);
	void re_edgeEdgeResponse(int thread_No, TargetPosition* target_pos);
	void re_edgeEdgeColliderResponse(int thread_No, TargetPosition* target_pos);




	void testNearestPoint();

	void pointTriangleCollisionPair(unsigned int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index);
	void pointTriangleCollisionPair(int thread_No);
	void edgeEdgeCollisionPair(int thread_No);
	void edgeEdgeCollisionPair(unsigned int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index);

	void pointTriangleColliderPair(int thread_No);
	void pointTriangleColliderPair(unsigned int thread_No, unsigned int pair_thread_No, unsigned int start_pair_index,
		unsigned int end_pair_index);

	void solveXPBDpointTriangleResponse(unsigned int thread_No, unsigned int pair_thread_No, unsigned int index_start, unsigned int index_end);
	void solveXPBDedgeEdgeResponse(unsigned int thread_No, unsigned int pair_thread_No, unsigned int index_start, unsigned int index_end);
	void solveXPBDpointTriangleColliderResponse(unsigned int thread_No, unsigned int pair_thread_No, unsigned int index_start, unsigned int index_end);


	void re_solveXPBDedgeEdgeResponse(unsigned int pair_thread_No, unsigned int index_start, unsigned int index_end);
	void re_solveXPBDpointTriangleResponse(unsigned int pair_thread_No, unsigned int index_start, unsigned int index_end);
	void re_solveXPBDpointTriangleColliderResponse(unsigned int pair_thread_No, unsigned int index_start, unsigned int index_end);
	//std::vector<double>* lambda;
	//std::vector<unsigned int>* collision_lambda_index_start;
	//double damp_stiffness;
	//double dt;
	void solveXPBDpointTriangleResponse(int thread_No);
	void solveXPBDpointTriangleColliderResponse(int thread_No);
	void solveXPBDedgeEdgeResponse(int thread_No);
	void re_solveXPBDedgeEdgeResponse(int thread_No);
	void re_solveXPBDpointTriangleResponse(int thread_No);
	void re_solveXPBDpointTriangleColliderResponse(int thread_No);
	Floor* floor;

	std::vector<double*> damp_collision;

	void XPBDfloorCollisionResponse();

	std::vector<std::vector<unsigned int>> floor_collision_vertex;
	void floorCollisionVertex(int thread_No);
	void re_FloorCollisionVertex(int thread_No);

	unsigned int use_method;
	bool CCD_compare=false;
};
