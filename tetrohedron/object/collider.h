#pragma once
#include"triangle_object.h"

class Collider:public TriangleObject
{
public:
	void draw(Camera* camera);
	void setSceneShader(Light& light, Camera* camera, float& far_plane);
	void loadMesh(OriMesh& ori_mesh, Thread* thread);
	std::vector<AABB> aabb;

private:
	void setMeshStruct(OriMesh& ori_mesh);
};

