#pragma once
#include"object.h"
#include"../mesh_struct/triangle_mesh_struct.h"

class TriangleObject :public Object
{
protected:
	struct Material
	{
		MeshMaterial front_material;
		MeshMaterial back_material;
	};
	void setMaterial(OriMesh& ori_mesh);
	void setRepresentativePrimitve();
public:
	TriangleMeshStruct mesh_struct;
	void drawShadow(Camera* camera, Shader* shader);
	void drawWireframe(Camera* camera, Shader* wireframe_shader);
	void setBuffer();
	void simpDraw(Camera* camera, Shader* shader);
	Material material;
	void reset();
	std::vector<int> vertex_from_rep_triangle_index;
	void getEdgeTriangleAABBPerThread(int thread_No);
	void getVertexAABBPerThread(int thread_No, bool has_tolerance);

private:

};

