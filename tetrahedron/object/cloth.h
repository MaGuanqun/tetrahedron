#pragma once

#include"triangle_object.h"


class Cloth:public TriangleObject
{
public:
	void setSceneShader(Light& light, Camera* camera, float& far_plane, Shader* object_shader_front);
	void loadMesh(OriMesh& ori_mesh, double density, Thread* thread);
	void draw(Camera* camera, Shader* object_shader_front);
	void setArea();
	double tolerance;
	
	std::vector<int>update_stiffness_iteration_number;
	void recordInitialMesh(SingleClothInfo& single_cloth_info_ref);
	SingleClothInfo single_cloth_info_ref;

	std::vector<double>length_stiffness;
	double bend_stiffness;
	double collision_stiffness[4];
	//std::vector<std::vector<double>>collision_stiffness_time_step_starts;
	//std::vector<std::vector<double>>collision_stiffness_time_step_starts_indicator;
	double collision_stiffness_update_indicator;
	double position_stiffness;
	double collision_stiffness_initial[4];

	void initial();
	void initialMouseChosenVertex();

	//void getTriangleAABBPerThread(int thread_No);
	void getEdgeTriangleAABBPerThread(int thread_No);
	void getVertexAABBPerThread(int thread_No);
	void obtainAABB();
	



	void findAllNeighborVertex(int face_index, double cursor_pos[3], double average_edge_length);

	void setTolerance(double* tolerance_ratio, double ave_edge_length);


	std::vector<std::vector<std::vector<int>>>triangle_neighbor_collider_triangle;
	void initialNeighborPrimitiveRecording(int cloth_num, int tetrahedron_num, int collider_num, bool use_BVH);

	std::vector<std::vector<std::vector<int>>>vertex_neighbor_obj_triangle;//except collider
	std::vector<std::vector<std::vector<int>>>vertex_neighbor_collider_triangle;
	std::vector<std::vector<std::vector<int>>>edge_neighbor_obj_edge;//except collider
	std::vector<std::vector<std::vector<int>>>collide_vertex_obj_triangle;//except collider
	std::vector<std::vector<std::vector<int>>>collide_vertex_collider_triangle;
	std::vector<std::vector<std::vector<int>>>collide_edge_obj_edge;//except collider
	std::vector<int> vertex_from_rep_triangle_index;
	

private:
	void setMeshStruct(double density, OriMesh& ori_mesh);

	void setMass(double density);
	void setAnchor();
	
	void setRepresentativePrimitve();
	void findNeighborVertex(int vertex_index, int recursion_deepth, std::vector<bool>& is_vertex_used);
	void initialHashAABB();

};

