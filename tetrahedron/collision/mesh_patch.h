#pragma once
#include"spatial_hashing.h"

class MeshPatch
{
public:

	void initialPatch(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
		Thread* thread, std::vector<std::vector<std::vector<unsigned int>>>* triangle_patch);

	void draw(Camera* camera);
	void setBuffer(unsigned int obj_index, unsigned int tetrahedron_start_index);
private:
	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;
	std::vector<std::vector<std::vector<unsigned int>>>* triangle_patch;
	void setInObjectData();

	unsigned int VAO;
	unsigned int VBO[2];
	void genBuffer();

	Shader* shader;
	
	unsigned int draw_element_num;
	void test();
};


