#pragma once
#include"object.h"
#include"../mesh_struct/tetrohedron_mesh_struct.h"
class Tetrohedron:public Object
{
private:
	void genShader();
	void setMeshStruct(double density, OriMesh& ori_mesh);
public:
	MeshMaterial material;
	TetrohedronMeshStruct mesh_struct;
	void drawShadow(Camera* camera, Shader* shader);
	void drawWireframe(Camera* camera);
	void setBuffer();
	void simpDraw(Camera* camera, Shader* shader);

	int tetrohedron_num;
	void loadMesh(OriMesh& ori_mesh, double density, Thread* thread);
	void draw(Camera* camera);
	void setSceneShader(Light& light, Camera* camera, float& far_plane);
};

