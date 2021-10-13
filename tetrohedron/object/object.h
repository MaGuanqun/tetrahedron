#pragma once
#include"../mesh_struct/mesh_struct.h"
#include"../basic/camera.h"
#include"../external\shader.h"
#include"../global_struct.h"
#include"../thread.h"
#include<array>
class Object
{
protected:
	unsigned int VAO;
	unsigned int VBO[3];
	unsigned int EBO;
	Shader* object_shader_back;
	Shader* object_shader_front;
	Shader* wireframe_shader;
	int total_thread_num;
	void genShader();
	void genBuffer();
	virtual void setMaterial(OriMesh& ori_mesh) =0;

	std::vector<std::array<double, 3>> ori_vertices;
	glm::vec3 wireframe_color;

	void getAABB(AABB& target, AABB& aabb0, AABB& aabb1, AABB& aabb2);
	void getAABB(AABB& target, AABB& aabb0, AABB& aabb1, double radius);
	void getAABB(AABB& target, AABB& aabb0, AABB& aabb1);
	Thread* thread;

private:
public:
	double mass;
	std::vector<double>coe_neighbor_vertex_force;
	std::vector<int>neighbor_vertex;
	void setWireframwColor(double* color);
	double density;
	std::vector<AABB>triangle_AABB;
	std::vector<AABB>edge_AABB;
	std::vector<AABB>vertex_AABB;
};
