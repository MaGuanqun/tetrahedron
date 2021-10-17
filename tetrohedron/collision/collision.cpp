#include"collision.h"

void Collision::initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, 
	std::vector<Tetrohedron>* tetrohedron, Thread* thread)
{
	this->cloth = cloth;
	this->collider = collider;
	this->tetrohedron = tetrohedron;
	this->thread = thread;
	initialBVH(cloth, collider, tetrohedron, thread);
}

void Collision::initialBVH(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrohedron>* tetrohedron, Thread* thread)
{
	cloth_BVH.resize(cloth->size());
	collider_BVH.resize(collider->size());
	tetrohedron_BVH.resize(tetrohedron->size());
	for (int i = 0; i < cloth->size(); ++i) {
		cloth_BVH[i].init((*cloth)[i].mesh_struct.faces.size(), (*cloth)[i].mesh_struct.face_index_begin_per_thread, thread);
	}
	for (int i = 0; i < collider->size(); ++i) {
		collider_BVH[i].init((*collider)[i].mesh_struct.faces.size(), (*collider)[i].mesh_struct.face_index_begin_per_thread, thread);
	}
	//for (int i = 0; i < tetrohedron->size(); ++i) {
	//	tetrohedron_BVH[i].init((*tetrohedron)[i].mesh_struct.triangle_indices.size()/3, (*tetrohedron)[i].mesh_struct.face_index_begin_per_thread, thread);
	//}
}



void Collision::buildBVH()
{
	for (int i = 0; i < cloth->size(); ++i) {
		(*cloth)[i].obtainAABB();
		cloth_BVH[i].buildBVH(&(*cloth)[i].triangle_AABB);
	}
	for (int i = 0; i < collider->size(); ++i) {
		(*collider)[i].obtainAABB();
		collider_BVH[i].buildBVH(&(*collider)[i].triangle_AABB);
	}
}



void Collision::searchTriangle(AABB& aabb, int compare_index, std::vector<std::vector<int>>* cloth_neighbor_index, std::vector<std::vector<int>>* collider_neighbor_index)
{
	for (int i = 0; i < cloth->size(); ++i) {
		cloth_BVH[i].search(aabb, compare_index, &((*cloth_neighbor_index)[i]),1,0, (*cloth)[i].triangle_AABB.size());
	}
	for (int i = 0; i < collider->size(); ++i) {
		collider_BVH[i].search(aabb, compare_index, &((*collider_neighbor_index)[i]), 1, 0, (*collider)[i].triangle_AABB.size());
	}

}

//void Collision::findAllTrianglePairs()
//{
//
//}