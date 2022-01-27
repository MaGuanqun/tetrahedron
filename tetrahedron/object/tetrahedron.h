#pragma once
#include"object.h"
#include"../mesh_struct/tetrahedron_mesh_struct.h"
class Tetrahedron:public Object
{
private:
	void setMeshStruct(double density, OriMesh& ori_mesh);
	void initialHashAABB();

public:
	MeshMaterial material;
	TetrahedronMeshStruct mesh_struct;
	void drawShadow(Camera* camera, Shader* shader);
	void drawWireframe(Camera* camera, Shader* wireframe_shader);
	void setBuffer();
	void simpDraw(Camera* camera, Shader* shader);

	int tetrahedron_num;
	void loadMesh(OriMesh& ori_mesh, double density, Thread* thread);
	void draw(Camera* camera, Shader* object_shader_front);
	void setSceneShader(Light& light, Camera* camera, float& far_plane, Shader* object_shader_front);
	
	double density;
	double ARAP_stiffness;
	double volume_preserve_stiffness;
	double position_stiffness;
	double collision_stiffness[4];
	double sigma_limit[2];//min max sigma for volume preserve
	void recordInitialMesh(SingleTetrahedronInfo& single_tetrahedron_info_ref);
	void initial();
	void reset();
	SingleTetrahedronInfo single_tetrahedron_info_ref;
	void findAllNeighborVertex(int face_index, double cursor_pos[3], double average_edge_length);
	void findNeighborVertex(int vertex_index, int recursion_deepth, std::vector<bool>& is_vertex_used);
	void findInnerVertex(std::vector<bool>& is_vertex_used);
	void setTolerance(double* tolerance_ratio, double ave_edge_length);
	double tolerance;
};

