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
	std::vector<std::vector<AABB>> aabb;
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
private:
	void setMeshStruct(double density, OriMesh& ori_mesh);
	std::vector<bool>is_vertex_used;
	double PC_radius_coe;
	void setMass(double density);
};

