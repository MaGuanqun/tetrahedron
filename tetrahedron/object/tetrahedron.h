#pragma once
#include"object.h"
#include"../mesh_struct/tetrahedron_mesh_struct.h"
class Tetrahedron:public Object
{
private:
	void setMeshStruct(double density, OriMesh& ori_mesh);


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
	SingleTetrahedronInfo single_tetrahedron_info_ref;
};

