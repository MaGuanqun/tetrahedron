#include"collision.h"

void Collision::initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, 
	std::vector<Tetrahedron>* tetrahedron, Thread* thread)
{
	this->cloth = cloth;
	this->collider = collider;
	this->tetrahedron = tetrahedron;
	this->thread = thread;
	initialBVH(cloth, collider, tetrahedron, thread);
}

void Collision::initialBVH(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron, Thread* thread)
{
	cloth_BVH.resize(cloth->size());
	collider_BVH.resize(collider->size());
	tetrahedron_BVH.resize(tetrahedron->size());
	for (int i = 0; i < cloth->size(); ++i) {
		cloth_BVH[i].init((*cloth)[i].mesh_struct.faces.size(), (*cloth)[i].mesh_struct.face_index_begin_per_thread, thread);
	}
	for (int i = 0; i < collider->size(); ++i) {
		collider_BVH[i].init((*collider)[i].mesh_struct.faces.size(), (*collider)[i].mesh_struct.face_index_begin_per_thread, thread);
	}
	//for (int i = 0; i < tetrahedron->size(); ++i) {
	//	tetrahedron_BVH[i].init((*tetrahedron)[i].mesh_struct.triangle_indices.size()/3, (*tetrahedron)[i].mesh_struct.face_index_begin_per_thread, thread);
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


void Collision::findAllTrianglePairs()
{
	buildBVH();
	thread->assignTask(this, FIND_TRIANGLE_PAIRS);
	thread->assignTask(this, FIND_PRIMITIVE_AROUND);

}


void Collision::test()
{
	findAllTrianglePairs();

	std::cout << "vertex " << std::endl;
	for (int i = 0; i < (*cloth)[0].vertex_neighbor_cloth_traingle[0][0].size(); ++i) {
		std::cout << (*cloth)[0].vertex_neighbor_cloth_traingle[0][0][i] << " ";
	}
	std::cout << std::endl;
	std::cout << "triangle " << std::endl;
	for (int i = 0; i < (*cloth)[0].triangle_neighbor_cloth_traingle[0][0].size(); ++i) {
		std::cout << (*cloth)[0].triangle_neighbor_cloth_traingle[0][0][i] << " ";
	}
	std::cout << std::endl;
	std::cout << "edge " << std::endl;
	for (int i = 0; i < (*cloth)[0].edge_neighbor_cloth_edge[0][0].size(); ++i) {
		std::cout << (*cloth)[0].edge_neighbor_cloth_edge[0][0][i] << " ";
	}
	std::cout << std::endl;
}

void Collision::searchTriangle(AABB& aabb, int compare_index, int cloth_No, std::vector<std::vector<int>>* cloth_neighbor_index, std::vector<std::vector<int>>* collider_neighbor_index)
{
	for (int i = cloth_No; i < cloth->size(); ++i) {
		(*cloth_neighbor_index)[i].reserve(5);
		(*cloth_neighbor_index)[i].clear();
		cloth_BVH[i].search(aabb, compare_index,i==cloth_No, &((*cloth_neighbor_index)[i]),1,0, (*cloth)[i].triangle_AABB.size());
	}
	for (int i = 0; i < collider->size(); ++i) {
		(*collider_neighbor_index)[i].reserve(5);
		(*collider_neighbor_index)[i].clear();
		collider_BVH[i].search(aabb, compare_index, false, &((*collider_neighbor_index)[i]), 1, 0, (*collider)[i].triangle_AABB.size());
	}
}


//FIND_TRIANGLE_PAIRS
void Collision::findAllTrianglePairs(int thread_No)
{
	int* thread_begin;
	for (int i = 0; i < cloth->size(); ++i) {
		thread_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
		for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
			searchTriangle((*cloth)[i].triangle_AABB[j], j, i, &(*cloth)[i].triangle_neighbor_cloth_traingle[j],
				&(*cloth)[i].triangle_neighbor_collider_triangle[j]);
		}
	}
}

//FIND_PRIMITIVE_AROUND
void Collision::findPrimitivesAround(int thread_No)
{
	findTriangleAroundVertex(thread_No);
	findEdgeAroundEdge(thread_No);
}



void Collision::findTriangleAroundVertex(int thread_No)
{
	int* thread_begin;
	std::vector<std::vector<std::vector<int>>>* triangle_neighbor_traingle;
	std::vector<std::vector<std::vector<int>>>* vertex_neighbor_traingle; //triangle index near vertex
	std::vector<int>* vertex_from_rep_triangle_index; // record of vertex's representative triangle index

	std::vector<std::array<int, 3>>* face_indices;//the record of every triangle's index
	std::vector<int>* vertex_neighbor_this_cloth_traingle; //vector to record the triangle index near vertex which triangle and vertex in same cloth
	std::vector<int>* triangle_neighbor_this_cloth_traingle;//vector to record the triangle index near triangle in same cloth

	std::vector<AABB>* vertex_aabb;
	std::vector<AABB>* compared_triangle_aabb;

	for (int i = 0; i < cloth->size(); ++i) {
		thread_begin = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread.data();
		triangle_neighbor_traingle = &(*cloth)[i].triangle_neighbor_cloth_traingle;
		vertex_neighbor_traingle = &(*cloth)[i].vertex_neighbor_cloth_traingle;
		vertex_from_rep_triangle_index = &(*cloth)[i].vertex_from_rep_triangle_index;

		vertex_aabb = &(*cloth)[i].vertex_AABB;

		for (int k = 0; k < cloth->size(); ++k) {
			compared_triangle_aabb= &(*cloth)[k].triangle_AABB;

			if (i != k) {
				for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
					vertex_neighbor_this_cloth_traingle = &(*vertex_neighbor_traingle)[j][k];
					triangle_neighbor_this_cloth_traingle = &(*triangle_neighbor_traingle)[(*vertex_from_rep_triangle_index)[j]][k];//vertex's represent triangle index = (*vertex_from_rep_triangle_index)[j];
					vertex_neighbor_this_cloth_traingle->clear();
					vertex_neighbor_this_cloth_traingle->reserve(triangle_neighbor_this_cloth_traingle->size());
					for (int m = 0; m < triangle_neighbor_this_cloth_traingle->size(); ++m) {						
						if ((*vertex_aabb)[j].AABB_intersection((*compared_triangle_aabb)[(*triangle_neighbor_this_cloth_traingle)[m]])) {
							(*vertex_neighbor_this_cloth_traingle).push_back((*triangle_neighbor_this_cloth_traingle)[m]);
						}
					}
				}
			}
			else {
				face_indices = &(*cloth)[k].mesh_struct.triangle_indices;
				for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
					vertex_neighbor_this_cloth_traingle = &(*vertex_neighbor_traingle)[j][k];
					triangle_neighbor_this_cloth_traingle = &(*triangle_neighbor_traingle)[(*vertex_from_rep_triangle_index)[j]][k];//vertex's represent triangle index = (*vertex_from_rep_triangle_index)[j];
					vertex_neighbor_this_cloth_traingle->clear();
					vertex_neighbor_this_cloth_traingle->reserve(triangle_neighbor_this_cloth_traingle->size());
					for (int m = 0; m < triangle_neighbor_this_cloth_traingle->size(); ++m) {
						if (!vertexInTriangle((*face_indices)[(*triangle_neighbor_this_cloth_traingle)[m]].data(), j)) {
							if ((*vertex_aabb)[j].AABB_intersection((*compared_triangle_aabb)[(*triangle_neighbor_this_cloth_traingle)[m]])) {
								(*vertex_neighbor_this_cloth_traingle).push_back((*triangle_neighbor_this_cloth_traingle)[m]);
							}
						}
					}
				}
			}
		}
		vertex_neighbor_traingle = &(*cloth)[i].vertex_neighbor_collider_triangle;
		triangle_neighbor_traingle = &(*cloth)[i].triangle_neighbor_collider_triangle;
		for (int k = 0; k < collider->size(); ++k) {
			compared_triangle_aabb = &(*collider)[k].triangle_AABB;
			for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
				vertex_neighbor_this_cloth_traingle = &(*vertex_neighbor_traingle)[j][k];
				triangle_neighbor_this_cloth_traingle = &(*triangle_neighbor_traingle)[(*vertex_from_rep_triangle_index)[j]][k];//vertex's represent triangle index = (*vertex_from_rep_triangle_index)[j];
				vertex_neighbor_this_cloth_traingle->clear();
				vertex_neighbor_this_cloth_traingle->reserve(triangle_neighbor_this_cloth_traingle->size());
				for (int m = 0; m < triangle_neighbor_this_cloth_traingle->size(); ++m) {
					if ((*vertex_aabb)[j].AABB_intersection((*compared_triangle_aabb)[(*triangle_neighbor_this_cloth_traingle)[m]])) {
						(*vertex_neighbor_this_cloth_traingle).push_back((*triangle_neighbor_this_cloth_traingle)[m]);
					}
				}
			}
		}
	}
}


void Collision::findEdgeAroundEdge(int thread_No)
{
	std::vector<std::array<int, 3>>* face_indices;//the record of every triangle's index
	int* thread_begin;
	std::vector<std::vector<std::vector<int>>>* triangle_neighbor_traingle;
	std::vector<std::vector<std::vector<int>>>* edge_neighbor_edge; //edge index near edge
	std::vector<int>* edge_from_rep_triangle_index; // record of edge's representative triangle index

	std::vector<int>* edge_neighbor_one_cloth_edge; //vector to record the edge index near edge which triangle and vertex in same cloth
	std::vector<int>* triangle_neighbor_one_cloth_traingle;//vector to record the triangle index near triangle in same cloth

	std::vector<int>* representative_edge_num;

	std::vector<MeshStruct::Face>* face;
	int face_index;

	std::vector<MeshStruct::Edge>* edge;

	std::vector<AABB>* edge_aabb;
	std::vector<AABB>* compared_edge_aabb;
	int compare_edge_index;

	for (int i = 0; i < cloth->size(); ++i) {
		thread_begin = (*cloth)[i].mesh_struct.edge_index_begin_per_thread.data();
		triangle_neighbor_traingle = &(*cloth)[i].triangle_neighbor_cloth_traingle;
		edge_neighbor_edge = &(*cloth)[i].edge_neighbor_cloth_edge;
		edge_from_rep_triangle_index = &(*cloth)[i].edge_from_rep_triangle_index;
		edge = &(*cloth)[i].mesh_struct.edges;
		edge_aabb = &(*cloth)[i].edge_AABB;

		for (int k = 0; k < cloth->size(); ++k) {
			compared_edge_aabb= &(*cloth)[k].edge_AABB;
			representative_edge_num = &(*cloth)[k].representative_edge_num;
			face = &(*cloth)[k].mesh_struct.faces;
			if (i != k) {
				for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
					edge_neighbor_one_cloth_edge = &(*edge_neighbor_edge)[j][k];
					triangle_neighbor_one_cloth_traingle = &(*triangle_neighbor_traingle)[(*edge_from_rep_triangle_index)[j]][k];
					//vertex's represent triangle index = (*vertex_from_rep_triangle_index)[j];
					edge_neighbor_one_cloth_edge->clear();
					edge_neighbor_one_cloth_edge->reserve(triangle_neighbor_one_cloth_traingle->size());					
					for (int m = 0; m < triangle_neighbor_one_cloth_traingle->size(); ++m) {
						face_index = (*triangle_neighbor_one_cloth_traingle)[m];
						for (int n = 0; n < (*representative_edge_num)[face_index]; ++n) {
							if ((*edge_aabb)[j].AABB_intersection((*compared_edge_aabb)[(*face)[face_index].edge[n]])) {
								(*edge_neighbor_one_cloth_edge).push_back((*face)[face_index].edge[n]);
							}
						}						
					}
				}
			}
			else {
				for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
					edge_neighbor_one_cloth_edge = &(*edge_neighbor_edge)[j][k];
					triangle_neighbor_one_cloth_traingle = &(*triangle_neighbor_traingle)[(*edge_from_rep_triangle_index)[j]][k];
					//vertex's represent triangle index = (*vertex_from_rep_triangle_index)[j];
					edge_neighbor_one_cloth_edge->clear();
					edge_neighbor_one_cloth_edge->reserve(triangle_neighbor_one_cloth_traingle->size());
					for (int m = 0; m < triangle_neighbor_one_cloth_traingle->size(); ++m) {
						face_index = (*triangle_neighbor_one_cloth_traingle)[m];
						for (int n = 0; n < (*representative_edge_num)[face_index]; ++n) {
							compare_edge_index = (*face)[face_index].edge[n];
							if (!edgeEdgeconnected((*edge)[j].vertex, (*edge)[compare_edge_index].vertex)) {
								if ((*edge_aabb)[j].AABB_intersection((*compared_edge_aabb)[compare_edge_index])) {
									(*edge_neighbor_one_cloth_edge).push_back(compare_edge_index);
								}								
							}							
						}
					}
				}
			}
		}
	}

}


bool Collision::edgeEdgeconnected(int* edge1, int* edge2)
{
	if (edge1[0] == edge2[0])
		return true;
	if (edge1[1] == edge2[1])
		return true;
	if (edge1[0] == edge2[1])
		return true;
	if (edge1[1] == edge2[0])
		return true;
	return false;
}

bool Collision::vertexInTriangle(int* face_index, int vertex_index)
{
	if (face_index[0] == vertex_index)
		return true;
	if (face_index[1] == vertex_index)
		return true;
	if (face_index[2] == vertex_index)
		return true;
	return false;
}
