#pragma once
#include"triangle_object.h"

class Collider:public TriangleObject
{
public:
	double tolerance;
	void draw(Camera* camera, Shader* object_shader_front);
	void setSceneShader(Light& light, Camera* camera, float& far_plane, Shader* object_shader_front);
	void loadMesh(OriMesh& ori_mesh, Thread* thread);
	void getTriangleAABBPerThread(int thread_No);
	void getVertexAABBPerThread(int thread_No);
	void obtainAABB();
	void setTolerance(double* tolerance_ratio, double ave_edge_length);

	std::vector<std::vector<std::vector<int>>>triangle_neighbor_cloth_vertex;
	std::vector<std::vector<std::vector<int>>>collider_triangle_cloth_vertex;
	void initialNeighborPrimitiveRecording(int cloth_num, int tetrahedron_num);

private:
	void setMeshStruct(OriMesh& ori_mesh);
};

