#pragma once
#include<cstring>
#include<string>
#include"basic/global.h"

struct MeshMaterial {
	std::string material_name;
	float Kd[3] = { 0.5,0.5,0.5 };
	float Ka[3] = { 0.1,0.1,0.1 };
	float Ks[3] = { 0.5,0.5,0.5 };
	float Tf[3] = { 1.0,1.0,1.0 };
	float Ni = 1.0;
	float Ns = 32.0;
	int illum = 4;
	double density = 50;
	MeshMaterial& operator=(MeshMaterial const& mesh_material)
	{
		Ni = mesh_material.Ni;
		Ns = mesh_material.Ns;
		illum = mesh_material.illum;
		density = mesh_material.density;
		memcpy(Kd, mesh_material.Kd, 24);
		memcpy(Ka, mesh_material.Ka, 24);
		memcpy(Ks, mesh_material.Ks, 24);
		memcpy(Tf, mesh_material.Tf, 24);
		material_name = mesh_material.material_name;
		return *this;
	}
};


struct OriMesh {
	int type;
	std::vector<std::array<double, 3>> vertices;
	std::vector<std::array<double, 2>> texture;
	std::vector<int> indices;
	MeshMaterial front_material;
	MeshMaterial back_material;
	OriMesh& operator=(OriMesh const& mesh)
	{
		type = mesh.type;
		vertices = mesh.vertices;
		if (!texture.empty()) {
			texture = mesh.texture;
		}
		indices=mesh.indices;
		front_material=mesh.front_material;
		back_material=mesh.back_material;
		return *this;
	}
};


struct SingleClothInfo {
	double density;
	double length_stiffness;			// stiffness of length constraint
	double position_stiffness;			// stiffness of position constraint
	double bending_stiffness;			// stiffness of bending constraint
	double collision_stiffness[4];			// stiffness of collision constraint //=0 body point triangle, =1 point-triangle =2 edge-edge =3 point-point
	double friction_stiffness_tangent;
	double friction_stiffness_normal;
	double virtual_length_stiffness;
	SingleClothInfo() {};
	SingleClothInfo(double density, double length_stiffness, double position_stiffness,
		double bending_stiffness, double* collision_stiffness, double friction_stiffness_tangent,
		double friction_stiffness_normal, double virtual_length_stiffness) {
		this->density = density;
		this->length_stiffness = length_stiffness;
		this->position_stiffness = position_stiffness;
		this->bending_stiffness = bending_stiffness;
		memcpy(this->collision_stiffness, collision_stiffness, 32);
		this->friction_stiffness_tangent = friction_stiffness_tangent;
		this->friction_stiffness_normal = friction_stiffness_normal;
		this->virtual_length_stiffness = virtual_length_stiffness;
	};
	SingleClothInfo& operator=(SingleClothInfo const& single_cloth_info)
	{
		this->density = single_cloth_info.density;
		this->length_stiffness = single_cloth_info.length_stiffness;
		this->position_stiffness = single_cloth_info.position_stiffness;
		this->bending_stiffness = single_cloth_info.bending_stiffness;
		memcpy(this->collision_stiffness, single_cloth_info.collision_stiffness, 32);
		this->friction_stiffness_tangent = single_cloth_info.friction_stiffness_tangent;
		this->friction_stiffness_normal = single_cloth_info.friction_stiffness_normal;
		this->virtual_length_stiffness = single_cloth_info.virtual_length_stiffness;
		return *this;
	}
};

struct SingleTetrahedronInfo {
	double density;
	double position_stiffness;			// stiffness of position constraint
	double ARAP_stiffness;		
	double collision_stiffness[4];			// stiffness of collision constraint //=0 body point triangle, =1 point-triangle =2 edge-edge =3 point-point
	SingleTetrahedronInfo() {};
	SingleTetrahedronInfo(double density, double position_stiffness,
		double ARAP_stiffness, double* collision_stiffness) {
		this->density = density;
		this->position_stiffness = position_stiffness;
		this->ARAP_stiffness = ARAP_stiffness;
		memcpy(this->collision_stiffness, collision_stiffness, 32);
	};
	SingleTetrahedronInfo& operator=(SingleTetrahedronInfo const& single_cloth_info)
	{
		this->density = single_cloth_info.density;
		this->position_stiffness = single_cloth_info.position_stiffness;
		this->ARAP_stiffness = single_cloth_info.ARAP_stiffness;
		memcpy(this->collision_stiffness, single_cloth_info.collision_stiffness, 32);
		return *this;
	}
};

struct UpdateClothStiffness
{
	bool update_length;
	double length_stiffness;
	bool update_bend;
	double bend_stiffness;

	bool  update_collision[4];
	double collision_stiffness[4];

	UpdateClothStiffness() {
		update_length = false;
		update_bend = false;
		memset(update_collision, 0, 4);
	}
};
