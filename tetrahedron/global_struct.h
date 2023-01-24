#pragma once
#include<cstring>
#include<string>
#include"basic/global.h"
#include"basic/enum_setting.h"

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
		memcpy(Kd, mesh_material.Kd, 12);
		memcpy(Ka, mesh_material.Ka, 12);
		memcpy(Ks, mesh_material.Ks, 12);
		memcpy(Tf, mesh_material.Tf, 12);
		material_name = mesh_material.material_name;
		return *this;
	}
	void setMaterial(float* kd, float* ka, float* ks)
	{
		memcpy(this->Kd, kd, 12);
		memcpy(this->Ka, ka, 12);
		memcpy(this->Ks, ks, 12);
	}
};


struct  StoreHessianWithOrderInConstraint
{
	std::vector<unsigned int> order_in_constraint;//store type(0 vt,1ee,2 tv_c, 3 ee_c,4 vt_c ), constraint index, start index in that constraint hessian
	std::array<double, 9>hessian;
	bool is_update;
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
		indices = mesh.indices;
		front_material = mesh.front_material;
		back_material = mesh.back_material;
		return *this;
	}
};


struct SingleClothInfo {
	double density;
	double length_stiffness[2];			// stiffness of length constraint
	double position_stiffness;			// stiffness of position constraint
	double bending_stiffness[2];			// stiffness of bending constraint
	double collision_stiffness[8];			// stiffness of collision constraint //=0 body point triangle, =1 point-triangle =2 edge-edge =3 point-point
	double friction_stiffness_tangent;
	double friction_stiffness_normal;
	double virtual_length_stiffness;
	SingleClothInfo() {};
	SingleClothInfo(double density, double length_stiffness, double position_stiffness,
		double bending_stiffness, double* collision_stiffness, double friction_stiffness_tangent,
		double friction_stiffness_normal, double virtual_length_stiffness, double damp_length, double damp_bending) {
		this->density = density;
		this->length_stiffness[0] = length_stiffness;
		this->length_stiffness[1] = damp_length;
		this->position_stiffness = position_stiffness;
		this->bending_stiffness[0] = bending_stiffness;
		this->bending_stiffness[1] = damp_bending;
		memcpy(this->collision_stiffness, collision_stiffness, 32);
		memcpy(this->collision_stiffness+4, collision_stiffness+ DAMP_BODY_POINT_TRIANGLE, 32);
		this->friction_stiffness_tangent = friction_stiffness_tangent;
		this->friction_stiffness_normal = friction_stiffness_normal;
		this->virtual_length_stiffness = virtual_length_stiffness;
	};
	SingleClothInfo& operator=(SingleClothInfo const& single_cloth_info)
	{
		this->density = single_cloth_info.density;
		memcpy(this->length_stiffness, single_cloth_info.length_stiffness,16);
		this->position_stiffness = single_cloth_info.position_stiffness;
		memcpy(this->bending_stiffness, single_cloth_info.bending_stiffness,16);
		memcpy(this->collision_stiffness, single_cloth_info.collision_stiffness, 64);
		this->friction_stiffness_tangent = single_cloth_info.friction_stiffness_tangent;
		this->friction_stiffness_normal = single_cloth_info.friction_stiffness_normal;
		this->virtual_length_stiffness = single_cloth_info.virtual_length_stiffness;
		return *this;
	}
};

struct SingleTetrahedronInfo {
	double density;
	double edge_length_stiffness;
	double position_stiffness;			// stiffness of position constraint
	double ARAP_stiffness[2];
	double volume_preserve_stiffness[2];
	double collision_stiffness[8];			// stiffness of collision constraint //=0 body point triangle, =1 point-triangle =2 edge-edge =3 point-point
	double sigma_limit[2];//max min sigma for volume preserve
	double youngs_modulus;
	double poisson_ratio;
	SingleTetrahedronInfo() {};
	SingleTetrahedronInfo(double density, double position_stiffness,
		double ARAP_stiffness, double volume_preserve_stiffness, double* collision_stiffness, double* sigma_limit,
		double youngs_modulus, double poisson_ratio, double edge_length_stiffness,
		double damp_ARAP, double damp_volume_preserve) {
		this->density = density;
		this->position_stiffness = position_stiffness;
		this->volume_preserve_stiffness[0] = volume_preserve_stiffness;
		this->volume_preserve_stiffness[1] = damp_volume_preserve;
		this->ARAP_stiffness[0] = ARAP_stiffness;
		this->ARAP_stiffness[1] = damp_ARAP;
		memcpy(this->collision_stiffness, collision_stiffness, 32);
		memcpy(this->collision_stiffness + 4, collision_stiffness + DAMP_BODY_POINT_TRIANGLE, 32);
		memcpy(this->sigma_limit, sigma_limit, 16);
		this->youngs_modulus = youngs_modulus;
		this->poisson_ratio = poisson_ratio;
		this->edge_length_stiffness = edge_length_stiffness;
	};
	SingleTetrahedronInfo& operator=(SingleTetrahedronInfo const& single_cloth_info)
	{
		this->density = single_cloth_info.density;
		this->position_stiffness = single_cloth_info.position_stiffness;
		memcpy(this->volume_preserve_stiffness, single_cloth_info.volume_preserve_stiffness, 16);
		memcpy(this->ARAP_stiffness, single_cloth_info.ARAP_stiffness,16);
		this->youngs_modulus = single_cloth_info.youngs_modulus;
		this->poisson_ratio = single_cloth_info.poisson_ratio;
		this->edge_length_stiffness = single_cloth_info.edge_length_stiffness;
		memcpy(this->collision_stiffness, single_cloth_info.collision_stiffness, 64);
		memcpy(this->sigma_limit, single_cloth_info.sigma_limit, 16);
		return *this;
	}
};

struct UpdateObjStiffness
{
	bool update_length;
	double length_stiffness[2];
	bool update_bend;
	double bend_stiffness[2];
	bool update_ARAP;
	double ARAP_stiffness[2];

	double youngs_modulus;
	double poisson_ratio;

	bool  update_collision[4];
	double collision_stiffness[8];
	bool update_tet_edge_length;
	double tet_edge_length_stiffness;
	UpdateObjStiffness() {
		update_length = false;
		update_bend = false;
		update_ARAP = false;
		update_tet_edge_length = false;
		memset(update_collision, 0, 4);
	}
};
