#pragma once
#include"triangle_object.h"

class Collider:public TriangleObject
{
public:
	double tolerance;
	void draw(Camera* camera);
	void setSceneShader(Light& light, Camera* camera, float& far_plane);
	void loadMesh(OriMesh& ori_mesh, Thread* thread);
	void getTriangleAABBPerThread(int thread_No);
	void obtainAABB();
	void setTolerance(double* tolerance_ratio, double ave_edge_length);

private:
	void setMeshStruct(OriMesh& ori_mesh);
};

