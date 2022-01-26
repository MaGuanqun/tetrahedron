#pragma once
#include"../mesh_struct/mesh_struct.h"
#include"../basic/camera.h"
#include"../external\shader.h"
#include"../global_struct.h"
#include"../thread.h"
#include<array>
#include"../basic/aabb.h"

class Object
{
protected:
	unsigned int VAO;
	unsigned int VBO[3];
	unsigned int EBO;

	int total_thread_num;
	void genBuffer();

	glm::vec3 wireframe_color;

	void getAABB(AABB& target, AABB& aabb0, AABB& aabb1, AABB& aabb2);
	void getAABB(AABB& target, AABB& aabb0, AABB& aabb1, double radius);
	void getAABB(AABB& target, AABB& aabb0, AABB& aabb1);
	Thread* thread;
	void setOrder(bool* in_this_triangle, int count, int* index);
	void setOrderEdge(bool* in_this_triangle, int count, int* index);
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
	std::vector<std::array<double, 3>> ori_vertices;

	std::vector<double> PC_radius;
	std::vector<std::vector<int>>hash_index_for_edge;
	std::vector<std::vector<int>>hash_index_for_vertex;

	std::vector<int> representative_vertex_num;
	std::vector<int> representative_edge_num;
	std::vector<int> vertex_from_rep_triangle_index;
	std::vector<int> edge_from_rep_triangle_index;

};
