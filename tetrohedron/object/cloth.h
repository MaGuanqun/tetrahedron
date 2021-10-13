#pragma once

#include"triangle_object.h"


class Cloth:public TriangleObject
{
public:
	void setSceneShader(Light& light, Camera* camera, float& far_plane);
	void loadMesh(OriMesh& ori_mesh, double density, Thread* thread);
	void draw(Camera* camera);
	void setAnchor(std::vector<int>& anchor_vertex);
	void setArea();
	std::vector<double> PC_radius;
	std::vector<std::vector<int>>hash_index_for_edge;
	std::vector<std::vector<int>>hash_index_for_vertex;
	std::vector<int>update_stiffness_iteration_number;
	void recordInitialMesh(SingleClothInfo& single_cloth_info_ref);
	SingleClothInfo single_cloth_info_ref;

	std::vector<double>length_stiffness;
	double bend_stiffness;
	std::vector<std::array<double, 4>>collision_stiffness;
	std::vector<std::array<double, 4>>collision_stiffness_time_step_starts;
	std::vector<std::array<double, 4>>collision_stiffness_time_step_starts_indicator;
	double collision_stiffness_update_indicator;
	double position_stiffness;

	void initial();
	void initialMouseChosenVertex();

	void getTriangleAABBPerThread(int thread_No);
	void getEdgeAABBPerThread(int thread_No);
	void getVertexAABBPerThread(int thread_No);
	void obtainAABB();
	double tolerance;

	std::vector<int> representative_vertex_num;
	std::vector<int> representative_edge_num;

	
private:
	void setMeshStruct(double density, OriMesh& ori_mesh);
	std::vector<bool>is_vertex_used;
	double PC_radius_coe;
	void setMass(double density);
	void setAnchor();
	void setOrder(bool* in_this_triangle, int count, int* index);
	void setOrderEdge(bool* in_this_triangle, int count, int* index);
	void setRepresentativePrimitve();
	void setRepresentativeVertex(std::vector<MeshStruct::Face>& face, std::vector<MeshStruct::Vertex>& vertex);
	void setRepresentativeEdge(std::vector<MeshStruct::Face>& face, std::vector<MeshStruct::Edge>& edge);
};

