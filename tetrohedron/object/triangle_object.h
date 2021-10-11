#pragma once
#include"object.h"
#include"../mesh_struct/triangle_mesh_struct.h"

class TriangleObject:public Object
{
protected:
	struct Material
	{
		MeshMaterial front_material;
		MeshMaterial back_material;
	};
	void setMaterial(OriMesh& ori_mesh);
public:
	TriangleMeshStruct mesh_struct;	
	void drawShadow(Camera* camera, Shader* shader);
	void drawWireframe(Camera* camera);
	void setBuffer();
	void simpDraw(Camera* camera, Shader* shader);
	Material material;
	void reset();
private:

};

