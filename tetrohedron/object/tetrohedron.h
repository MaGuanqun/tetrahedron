#pragma once
#include"object.h"
#include"../mesh_struct/tetrohedron_mesh_struct.h"
class Tetrohedron:public Object
{
private:
	void setMaterial(OriMesh& ori_mesh);
public:
	MeshMaterial material;
	TetrohedronMeshStruct mesh_struct;
	void drawShadow(Camera* camera, Shader* shader);
	void drawWireframe(Camera* camera);
	void setBuffer();
	void simpDraw(Camera* camera, Shader* shader);
};

