#pragma once
#include"spatial_hashing.h"

class MeshPatch
{
public:
	~MeshPatch();
	void initialPatch(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
		Thread* thread);//, std::vector<std::vector<std::vector<unsigned int>>>* triangle_patch,std::vector<std::vector<std::vector<unsigned int>>>* patch_vertex

	void draw(Camera* camera);
	void setBuffer(unsigned int obj_index, unsigned int tetrahedron_start_index);
	void findVertex(int thread_No);
	void obtainAABB(int thread_No);

	std::vector<std::vector<std::array<double,6>>> patch_AABB;
	bool** patch_is_intersect;
	std::vector<std::vector<unsigned int>> patch_index_start_per_thread;

	std::vector<std::vector<unsigned int>> triangle_index_noted_as_true;
	std::vector<std::vector<unsigned int>> triangle_index_noted_begin_per_thread;

	void findAllNotedTrueTriangle();
	void decideTriangleIndexSize(int thread_No);
	void findNotedTrueTriangleIndex(int thread_No);
	void findAllNotedTrueTriangleSingleThread();
	void test();
private:
	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;
	std::vector<std::vector<std::vector<unsigned int>>> triangle_patch;
	std::vector<std::vector<std::vector<unsigned int>>> patch_vertex; //for tetrahedron, index in toal not surface

	void setInObjectData();

	unsigned int VAO;
	unsigned int VBO[2];
	void genBuffer();

	Shader* shader;
	
	unsigned int draw_element_num;
	void test(unsigned int obj_index, unsigned int tetrahedron_start_index);
	void findVertex();
	

	void findVertex(std::vector<bool>& is_vertex_used,std::array<int, 3>* triangle_indices,
		std::vector<unsigned int>* patch, unsigned int start, unsigned int end, std::vector<unsigned int>* patch_vertex);
	void findTetVertex(std::vector<bool>& is_vertex_used, std::array<int, 3>* triangle_indices,
		std::vector<unsigned int>* patch, unsigned int start, unsigned int end, std::vector<unsigned int>* patch_vertex,
		std::vector<int>& surface_vertex_index);
	unsigned int total_obj_num; //include collider
	unsigned int tetrahedron_end_index;
	void getAABB(double* target, double* aabb0);

	void obtainAABB(std::array<double,6>* aabb, unsigned int start, unsigned int end, std::vector<unsigned int>* vertex_patch, std::array<double, 6>* vertex_aabb);

	std::vector<std::vector<unsigned int>> triangle_size_start_per_thread;//triangle size noted as true;
	unsigned int total_thread_num;

	//std::vector<std::vector<std::vector<unsigned int>>> triangle_patch_noted_true;
};


