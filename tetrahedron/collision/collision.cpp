#include"collision.h"

void Collision::initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, 
	std::vector<Tetrahedron>* tetrahedron, Thread* thread)
{
	this->cloth = cloth;
	this->collider = collider;
	this->tetrahedron = tetrahedron;
	this->thread = thread;
	initialBVH(cloth, collider, tetrahedron, thread);
	initialTargetPos(cloth, tetrahedron, thread);
	initialSpatialHashing(cloth, collider, tetrahedron, thread);
	std::cout << "6" << std::endl;
}

void Collision::initialTargetPos(std::vector<Cloth>* cloth,	std::vector<Tetrahedron>* tetrahedron, Thread* thread)
{
	thread_num = thread->thread_num;
	cloth_target_pos_per_thread.resize(thread_num);
	for (int i = 0; i < cloth_target_pos_per_thread.size(); ++i) {
		cloth_target_pos_per_thread[i].initialSet(cloth->size());
		for (int j = 0; j < cloth->size(); ++j) {
			cloth_target_pos_per_thread[i].initialSet2(j, (*cloth)[j].ori_vertices.size());
		}
	}
	cloth_target_pos.initialSet(cloth->size());
	for (int j = 0; j < cloth->size(); ++j) {
		cloth_target_pos.initialSet2(j, (*cloth)[j].ori_vertices.size());
	}
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

void Collision::initialSpatialHashing(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron, Thread* thread)
{
	spatial_hashing.setInObject(cloth, collider, tetrahedron, thread);
}


void Collision::getAABB()
{
	for (int i = 0; i < cloth->size(); ++i) {
		(*cloth)[i].obtainAABB();
	}
	for (int i = 0; i < collider->size(); ++i) {
		(*collider)[i].obtainAABB();
	}
}

void Collision::buildBVH()
{
	for (int i = 0; i < cloth->size(); ++i) {
		cloth_BVH[i].buildBVH(&(*cloth)[i].triangle_AABB);
	}
	for (int i = 0; i < collider->size(); ++i) {
		collider_BVH[i].buildBVH(&(*collider)[i].triangle_AABB);
	}
}


void Collision::globalCollision()
{
	getAABB();
	buildBVH();
	for (int i = 0; i < cloth->size(); ++i) {
		(*cloth)[i].mesh_struct.getNormal();
	}
	thread->assignTask(this, FIND_TRIANGLE_PAIRS);
	thread->assignTask(this, FIND_PRIMITIVE_AROUND);
	thread->assignTask(this, GLOBAL_COLLISION_DETECTION);
	sumTargetPosition();
}

void Collision::updateCollisionPosition()
{
	for (int i = 0; i < cloth->size(); ++i) {
		(*cloth)[i].mesh_struct.getNormal();
	}
	thread->assignTask(this, RE_DETECTION);
	resumTargetPosition();
}



void Collision::resumTargetPosition()
{
	cloth_target_pos.partialInitial();
	thread->assignTask(this, RESUM_TARGET_POSITION);
	for (int i = 0; i < thread_num; ++i) {
		for (int j = 0; j < cloth->size(); ++j) {
			cloth_target_pos.collision_energy[j] += cloth_target_pos_per_thread[i].collision_energy[j];
		}
	}
}

void Collision::sumTargetPosition()
{
	cloth_target_pos.initial();
	for (int i = 0; i < thread_num; ++i) {
		for (int j = 0; j < cloth->size(); ++j) {
			cloth_target_pos.collision_energy[j] += cloth_target_pos_per_thread[i].collision_energy[j];
		}
	}
	thread->assignTask(this, SUM_TARGET_POSITION);
}

//RESUM_TARGET_POSITION
void Collision::resumTargetPositionPerThread(int thread_id)
{
	std::vector<int>* index_begin;
	bool* need_update;
	std::vector<std::array<double, 3>>* b_sum;
	std::vector<std::array<double, 3>>* b_sum_per_thread;
	for (int j = 0; j < cloth->size(); ++j) {
		b_sum = &cloth_target_pos.b_sum[j];
		index_begin = &(*cloth)[j].mesh_struct.vertex_index_begin_per_thread;
		for (int i = 0; i < thread_num; ++i) {			
			need_update = cloth_target_pos_per_thread[i].need_update[j];
			b_sum_per_thread = &cloth_target_pos_per_thread[i].b_sum[j];
			for (int k = (*index_begin)[thread_id]; k < (*index_begin)[thread_id + 1]; ++k) {
				if (need_update[k]) {
					SUM((*b_sum)[k], (*b_sum)[k], (*b_sum_per_thread)[k]);
				}
			}
		}
	}
}
//SUM_TARGET_POSITION
void Collision::sumTargetPositionPerThread(int thread_id)
{
	std::vector<int>* index_begin;
	bool* need_update;
	std::vector<std::array<double, 3>>* b_sum;
	std::vector<std::array<double, 3>>* b_sum_per_thread;
	bool* global_need_update;
	double* stiffness;
	double* global_stiffness;

	for (int j = 0; j < cloth->size(); ++j) {
		global_need_update = cloth_target_pos.need_update[j];
		b_sum = &cloth_target_pos.b_sum[j];
		index_begin = &(*cloth)[j].mesh_struct.vertex_index_begin_per_thread;
		global_stiffness = cloth_target_pos.stiffness[j].data();
		for (int i = 0; i < thread_num; ++i) {
			need_update = cloth_target_pos_per_thread[i].need_update[j];
			b_sum_per_thread = &cloth_target_pos_per_thread[i].b_sum[j];
			stiffness = cloth_target_pos_per_thread[i].stiffness[j].data();
			for (int k = (*index_begin)[thread_id]; k < (*index_begin)[thread_id + 1]; ++k) {
				if (need_update[k]) {
					global_need_update[k] = true;
					SUM((*b_sum)[k], (*b_sum)[k],(*b_sum_per_thread)[k]);
					global_stiffness[k] += stiffness[k];
				}
			}
		}
	}
}


//GLOBAL_COLLISION_DETECTION
void Collision::collisionDetection(int thread_No)
{
	TargetPosition* target_pos = &cloth_target_pos_per_thread[thread_No];
	target_pos->initial();
	std::vector<int>* index_begin;
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		index_begin = &(*cloth)[cloth_No].mesh_struct.vertex_index_begin_per_thread;
		for (int i =(*index_begin)[thread_No]; i < (*index_begin)[thread_No + 1]; ++i) {
			pointSelfTriangleCollisionDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].vertex_neighbor_cloth_traingle[i],
				&(*cloth)[cloth_No].collide_vertex_cloth_traingle[i], &(*cloth)[cloth_No].mesh_struct, (*cloth)[cloth_No].PC_radius[SELF_POINT_TRIANGLE], target_pos);
			pointColliderTriangleCollisionDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].vertex_neighbor_collider_triangle[i],
				&(*cloth)[cloth_No].collide_vertex_collider_triangle[i], &(*cloth)[cloth_No].mesh_struct, (*cloth)[cloth_No].PC_radius[BODY_POINT_TRIANGLE], target_pos);
		}
		index_begin = &(*cloth)[cloth_No].mesh_struct.edge_index_begin_per_thread;
		for (int i = (*index_begin)[thread_No]; i < (*index_begin)[thread_No + 1]; ++i) {
			edgeSelfEdgeCollisionDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].edge_neighbor_cloth_edge[i], &(*cloth)[cloth_No].collide_edge_cloth_edge[i], 
				&(*cloth)[cloth_No].mesh_struct, (*cloth)[cloth_No].PC_radius[SELF_EDGE_EDGE], target_pos);
		}
	}
}

//RE_DETECTION
void Collision::collisionReDetection(int thread_No)
{
	TargetPosition* target_pos = &cloth_target_pos_per_thread[thread_No];
	target_pos->partialInitial();
	std::vector<int>* index_begin;
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		index_begin = &(*cloth)[cloth_No].mesh_struct.vertex_index_begin_per_thread;
		for (int i = (*index_begin)[thread_No]; i < (*index_begin)[thread_No + 1]; ++i) {
			pointSelfTriangleCollisionReDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].collide_vertex_cloth_traingle[i], &(*cloth)[cloth_No].mesh_struct,
				(*cloth)[cloth_No].PC_radius[SELF_POINT_TRIANGLE], &(*cloth)[cloth_No].collision_stiffness[SELF_POINT_TRIANGLE], target_pos);
			pointColliderTriangleCollisionReDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].collide_vertex_collider_triangle[i], &(*cloth)[cloth_No].mesh_struct,
				(*cloth)[cloth_No].PC_radius[BODY_POINT_TRIANGLE], &(*cloth)[cloth_No].collision_stiffness[BODY_POINT_TRIANGLE], target_pos);
		}
		index_begin = &(*cloth)[cloth_No].mesh_struct.edge_index_begin_per_thread;
		for (int i = (*index_begin)[thread_No]; i < (*index_begin)[thread_No + 1]; ++i) {
			edgeSelfEdgeCollisionReDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].collide_edge_cloth_edge[i], &(*cloth)[cloth_No].mesh_struct, 
				(*cloth)[cloth_No].PC_radius[SELF_EDGE_EDGE], (*cloth)[cloth_No].collision_stiffness[SELF_EDGE_EDGE], target_pos);
		}
	}
}



void Collision::pointColliderTriangleCollisionDetection(int thread_No, int vertex_index, int cloth_No,
	std::vector<std::vector<int>>* vertex_neighbor_triangle, std::vector<std::vector<int>>* collide_vertex_triangle, MeshStruct* vertex_mesh, double radius0,
	TargetPosition* target_pos)
{
	MeshStruct* triangle_mesh;
	std::vector<int>* neighbor_triangle;
	std::vector<int>* collide_triangle;
	double radius1;
	for (int i = 0; i < collider->size(); ++i) {
		triangle_mesh = &(*collider)[i].mesh_struct;
		neighbor_triangle = &(*vertex_neighbor_triangle)[i];
		collide_triangle = &(*collide_vertex_triangle)[i];
		collide_triangle->clear();
		collide_triangle->reserve(neighbor_triangle->size());
		radius1 = (*collider)[i].tolerance;
		for (int k = 0; k < neighbor_triangle->size(); ++k) {
			if (checkPointColliderTriangleCollision(vertex_mesh, triangle_mesh, radius0 + radius1, vertex_index, (*neighbor_triangle)[k], cloth_No,
				i, target_pos,true, (*cloth)[cloth_No].collision_stiffness[BODY_POINT_TRIANGLE])) {
				collide_triangle->push_back((*neighbor_triangle)[k]);
			}
		}
	}
}



void Collision::pointSelfTriangleCollisionDetection(int thread_No,int vertex_index, int cloth_No,
	std::vector<std::vector<int>>* vertex_neighbor_triangle, std::vector<std::vector<int>>* collide_vertex_triangle, MeshStruct* vertex_mesh, double radius0, TargetPosition* target_pos)
{
	MeshStruct* triangle_mesh;
	std::vector<int>* neighbor_triangle;
	std::vector<int>* collide_triangle;
	double radius1;
	for (int i = 0; i < cloth->size(); ++i) {
		triangle_mesh = &(*cloth)[i].mesh_struct;
		neighbor_triangle = &(*vertex_neighbor_triangle)[i];
		collide_triangle = &(*collide_vertex_triangle)[i];
		collide_triangle->clear();
		collide_triangle->reserve(neighbor_triangle->size());
		radius1= (*cloth)[i].PC_radius[SELF_POINT_TRIANGLE];
		for (int k = 0; k < neighbor_triangle->size(); ++k) {
			if (checkPointTriangleCollision(vertex_mesh, triangle_mesh, radius0+radius1,vertex_index, (*neighbor_triangle)[k],cloth_No,
				i, target_pos,true, (*cloth)[cloth_No].collision_stiffness[SELF_POINT_TRIANGLE], (*cloth)[i].collision_stiffness[SELF_POINT_TRIANGLE])) {
				collide_triangle->push_back((*neighbor_triangle)[k]);
			}
		}
	}	
}

void Collision::pointColliderTriangleCollisionReDetection(int thread_No, int vertex_index, int cloth_No, std::vector<std::vector<int>>* collide_vertex_triangle, MeshStruct* vertex_mesh,
	double radius0, std::vector<double>* collision_stiffness, TargetPosition* target_postion_)
{
	MeshStruct* triangle_mesh;
	std::vector<int>* collide_triangle;	
	double radius1;
	for (int i = 0; i < collider->size(); ++i) {
		triangle_mesh = &(*collider)[i].mesh_struct;
		collide_triangle = &(*collide_vertex_triangle)[i];
		radius1 = (*collider)[i].tolerance;
		for (int k = 0; k < collide_triangle->size(); ++k) {
			if (!checkPointColliderTriangleCollision(vertex_mesh, triangle_mesh, radius0 + radius1, vertex_index, (*collide_triangle)[k], cloth_No,
				i, target_postion_, false, (*collision_stiffness))) {
				addTargetPosToSystem(target_postion_->b_sum[cloth_No][vertex_index].data(),
					target_postion_->collision_energy[cloth_No], vertex_mesh->vertex_position[vertex_index].data(), vertex_mesh->vertex_position[vertex_index].data(), (*collision_stiffness)[vertex_index]);
			}
		}
	}
}

void Collision::pointSelfTriangleCollisionReDetection(int thread_No, int vertex_index, int cloth_No,std::vector<std::vector<int>>* collide_vertex_triangle, MeshStruct* vertex_mesh,
	double radius0, std::vector<double>* collision_stiffness, TargetPosition* target_postion_)
{
	MeshStruct* triangle_mesh;
	std::vector<int>* collide_triangle;
	double radius1;
	int* triangle_vertex_index;
	int triangle_vertex;
	std::vector<double>* triangle_collision_stiffness;
	for (int i = 0; i < cloth->size(); ++i) {
		triangle_mesh = &(*cloth)[i].mesh_struct;
		collide_triangle = &(*collide_vertex_triangle)[i];
		radius1 = (*cloth)[i].PC_radius[SELF_POINT_TRIANGLE];
		triangle_collision_stiffness = &(*cloth)[i].collision_stiffness[SELF_POINT_TRIANGLE];
		for (int k = 0; k < collide_triangle->size(); ++k) {
			if (!checkPointTriangleCollision(vertex_mesh, triangle_mesh, radius0 + radius1, vertex_index, (*collide_triangle)[k], cloth_No,
				i, target_postion_, false, (*collision_stiffness), (*triangle_collision_stiffness))) {
				addTargetPosToSystem(target_postion_->b_sum[cloth_No][vertex_index].data(),
					target_postion_->collision_energy[cloth_No], vertex_mesh->vertex_position[vertex_index].data(), vertex_mesh->vertex_position[vertex_index].data(), (*collision_stiffness)[vertex_index]);
				
				triangle_vertex_index = triangle_mesh->triangle_indices[(*collide_triangle)[k]].data();
				for (int j = 0; j < 3; ++j) {
					triangle_vertex = triangle_vertex_index[j];
					addTargetPosToSystem(target_postion_->b_sum[i][triangle_vertex].data(),
						target_postion_->collision_energy[i], triangle_mesh->vertex_position[triangle_vertex].data(), triangle_mesh->vertex_position[triangle_vertex].data(), (*triangle_collision_stiffness)[triangle_vertex]);
				}
			}
		}
	}
}

void Collision::edgeSelfEdgeCollisionReDetection(int thread_No, int edge_index, int cloth_No, std::vector<std::vector<int>>* collide_edge_edge, TriangleMeshStruct* edge_mesh, 
	double radius0, std::vector<double>& collision_stiffness, TargetPosition* target_postion_)
{
	TriangleMeshStruct* compare_edge_mesh;
	std::vector<int>* collide_edge;
	double radius1;
	std::vector<double>* compare_collision_stiffness;
	int* edge_vertex = edge_mesh->edges[edge_index].vertex;
	int* compare_edge_index;
	for (int i = 0; i < cloth->size(); ++i) {
		compare_edge_mesh = &(*cloth)[i].mesh_struct;
		collide_edge = &(*collide_edge_edge)[i];
		radius1 = (*cloth)[i].PC_radius[SELF_EDGE_EDGE];
		compare_collision_stiffness= &(*cloth)[i].collision_stiffness[SELF_EDGE_EDGE];
		for (int k = 0; k < collide_edge->size(); ++k) {
			compare_edge_index = compare_edge_mesh->edges[(*collide_edge)[k]].vertex;
			if (!checkEdgeEdgeCollision(edge_mesh, compare_edge_mesh, radius0 + radius1, edge_index, (*collide_edge)[k], cloth_No,
				i, target_postion_, false, collision_stiffness, (*compare_collision_stiffness))) {
				for (int j = 0; j < 2; ++j) {
					addTargetPosToSystem(target_postion_->b_sum[cloth_No][edge_vertex[j]].data(),target_postion_->collision_energy[cloth_No], 
						edge_mesh->vertex_position[edge_vertex[j]].data(), edge_mesh->vertex_position[edge_vertex[j]].data(), collision_stiffness[edge_vertex[j]]);
					addTargetPosToSystem(target_postion_->b_sum[i][compare_edge_index[j]].data(),target_postion_->collision_energy[i],
						compare_edge_mesh->vertex_position[compare_edge_index[j]].data(), compare_edge_mesh->vertex_position[compare_edge_index[j]].data(), (*compare_collision_stiffness)[compare_edge_index[j]]);
				}
			}
		}
	}
}

void Collision::edgeSelfEdgeCollisionDetection(int thread_No, int edge_index, int cloth_No,
	std::vector<std::vector<int>>* edge_neighbor_edge, std::vector<std::vector<int>>* collide_edge_edge, TriangleMeshStruct* edge_mesh, double radius0, TargetPosition* target_pos)
{
	TriangleMeshStruct* compare_edge_mesh;
	std::vector<int>* neighbor_edge;
	std::vector<int>* collide_edge;
	double radius1;
	for (int i = 0; i < cloth->size(); ++i) {
		compare_edge_mesh = &(*cloth)[i].mesh_struct;
		neighbor_edge = &(*edge_neighbor_edge)[i];
		collide_edge = &(*collide_edge_edge)[i];
		collide_edge->clear();
		collide_edge->reserve(neighbor_edge->size());
		radius1 = (*cloth)[i].PC_radius[SELF_EDGE_EDGE];
		for (int k = 0; k < neighbor_edge->size(); ++k) {
			if (checkEdgeEdgeCollision(edge_mesh, compare_edge_mesh, radius0 + radius1, edge_index, (*neighbor_edge)[k], cloth_No,
				i, target_pos, true,(*cloth)[cloth_No].collision_stiffness[SELF_EDGE_EDGE], (*cloth)[i].collision_stiffness[SELF_EDGE_EDGE])) {
				collide_edge->push_back((*neighbor_edge)[k]);
			}
		}
	}
}


bool Collision::checkPointTriangleCollision(MeshStruct* vertex_mesh, MeshStruct* triangle_mesh,
	double radius, int vertex_index, int triangle_index, int vertex_cloth_No, int triangle_cloth_No, TargetPosition* target_position, bool new_collision_registration,
	std::vector<double>& vertex_collision_stiffness, std::vector<double>& triangle_collision_stiffness)
{	
	std::vector<double*>initial_triangle_pos(3); std::vector<double*>current_triangle_pos(3);
	double triangle_mass[3];
	int* triangle_vertex_index = triangle_mesh->triangle_indices[triangle_index].data();
	for (int i = 0; i < 3; ++i) {
		initial_triangle_pos[i] = triangle_mesh->vertex_for_render[triangle_vertex_index[i]].data();
		current_triangle_pos[i] = triangle_mesh->vertex_position[triangle_vertex_index[i]].data();
		triangle_mass[i] = triangle_mesh->mass[triangle_vertex_index[i]];
	}
	double vertex_target_pos[3];
	std::vector<std::array<double, 3>> triangle_target_pos;
	if (predictive_contact.pointTriangleCollision(vertex_mesh->vertex_for_render[vertex_index].data(), vertex_mesh->vertex_position[vertex_index].data(), initial_triangle_pos, current_triangle_pos,
		triangle_mesh->face_norm_for_render[triangle_index].data(), triangle_mesh->face_norm[triangle_index].data(), vertex_target_pos, triangle_target_pos, radius, triangle_mesh->triangle_normal_magnitude_reciprocal[triangle_index],
		vertex_mesh->mass[vertex_index], triangle_mass)) {		
		addTargetPosToSystem(target_position->b_sum[vertex_cloth_No][vertex_index].data(),
			target_position->collision_energy[vertex_cloth_No], vertex_mesh->vertex_position[vertex_index].data(), vertex_target_pos, vertex_collision_stiffness[vertex_index]);
		int triangle_vertex;
		for (int i = 0; i < 3; ++i) {
			triangle_vertex = triangle_vertex_index[i];
			addTargetPosToSystem(target_position->b_sum[triangle_cloth_No][triangle_vertex].data(),
				target_position->collision_energy[triangle_cloth_No], triangle_mesh->vertex_position[triangle_vertex].data(), triangle_target_pos[i].data(), triangle_collision_stiffness[triangle_vertex]);
		}
		if (new_collision_registration) {
			target_position->stiffness[vertex_cloth_No][vertex_index] += vertex_collision_stiffness[vertex_index];
			target_position->need_update[vertex_cloth_No][vertex_index] = true;
			for (int i = 0; i < 3; ++i) {
				target_position->stiffness[triangle_cloth_No][triangle_vertex_index[i]] += triangle_collision_stiffness[triangle_vertex];
				target_position->need_update[triangle_cloth_No][triangle_vertex_index[i]] = true;
			}
		}

		return true;
	}
	return false;
}

bool Collision::checkEdgeEdgeCollision(TriangleMeshStruct* edge_mesh, TriangleMeshStruct* compare_mesh,
	double radius, int edge_index, int compare_edge_index, int edge_cloth_No, int compare_cloth_No, TargetPosition* target_position, bool new_collision_registration,
	std::vector<double>& collision_stiffness, std::vector<double>& compare_collision_stiffness)
{
	double mass[4];
	int* edge_vertex_index = edge_mesh->edges[edge_index].vertex;
	int* compare_edge_vertex_index = compare_mesh->edges[compare_edge_index].vertex;
	for (int i = 0; i < 2; ++i) {		
		mass[i] = edge_mesh->mass[edge_vertex_index[i]];
		mass[2+i] = edge_mesh->mass[compare_edge_vertex_index[i]];
	}
	std::vector<std::array<double, 3>> target_pos;
	std::vector<std::array<double, 3>> compare_target_pos;
	if (predictive_contact.edgeEdgeCollision(target_pos, compare_target_pos,radius, edge_mesh->vertex_position[edge_vertex_index[0]].data(), edge_mesh->vertex_position[edge_vertex_index[1]].data(),
		edge_mesh->vertex_for_render[edge_vertex_index[0]].data(), edge_mesh->vertex_for_render[edge_vertex_index[1]].data(), compare_mesh->vertex_position[compare_edge_vertex_index[0]].data(),
		compare_mesh->vertex_position[compare_edge_vertex_index[1]].data(), compare_mesh->vertex_for_render[compare_edge_vertex_index[0]].data(), compare_mesh->vertex_for_render[compare_edge_vertex_index[1]].data(),
		mass)) {
		for (int i = 0; i < 2; ++i) {
			addTargetPosToSystem(target_position->b_sum[edge_cloth_No][edge_vertex_index[i]].data(),
				target_position->collision_energy[edge_cloth_No], edge_mesh->vertex_position[edge_vertex_index[i]].data(), target_pos[i].data(), collision_stiffness[edge_vertex_index[i]]);
			addTargetPosToSystem(target_position->b_sum[compare_cloth_No][compare_edge_vertex_index[i]].data(),
				target_position->collision_energy[compare_cloth_No], compare_mesh->vertex_position[compare_edge_vertex_index[i]].data(), compare_target_pos[i].data(), compare_collision_stiffness[compare_edge_vertex_index[i]]);
		}	
		if (new_collision_registration) {
			for (int i = 0; i < 2; ++i) {
				target_position->stiffness[edge_cloth_No][edge_vertex_index[i]] += collision_stiffness[edge_vertex_index[i]];
				target_position->need_update[edge_cloth_No][edge_vertex_index[i]] = true;
				target_position->stiffness[compare_cloth_No][compare_edge_vertex_index[i]] += compare_collision_stiffness[compare_edge_vertex_index[i]];
				target_position->need_update[edge_cloth_No][compare_edge_vertex_index[i]] = true;
			}
		}
		return true;
	}
	return false;
}

bool Collision::checkPointColliderTriangleCollision(MeshStruct* vertex_mesh, MeshStruct* triangle_mesh,
	double radius, int vertex_index, int triangle_index, int vertex_cloth_No, int triangle_obj_No, TargetPosition* target_position, bool new_collision_registration,
	std::vector<double>& vertex_collision_stiffness)
{
	std::vector<double*>initial_triangle_pos(3); std::vector<double*>current_triangle_pos(3);
	int* triangle_vertex_index = triangle_mesh->triangle_indices[triangle_index].data();
	for (int i = 0; i < 3; ++i) {
		initial_triangle_pos[i] = triangle_mesh->vertex_for_render[triangle_vertex_index[i]].data();
		current_triangle_pos[i] = triangle_mesh->vertex_position[triangle_vertex_index[i]].data();
	}
	double vertex_target_pos[3];
	if (predictive_contact.pointColliderTriangleCollision(vertex_mesh->vertex_for_render[vertex_index].data(), vertex_mesh->vertex_position[vertex_index].data(), initial_triangle_pos, current_triangle_pos,
		triangle_mesh->face_norm_for_render[triangle_index].data(), triangle_mesh->face_norm[triangle_index].data(), vertex_target_pos, radius)) {
		if (new_collision_registration) {
			target_position->stiffness[vertex_cloth_No][vertex_index] += vertex_collision_stiffness[vertex_index];
			target_position->need_update[vertex_cloth_No][vertex_index] = true;
		}
		addTargetPosToSystem(target_position->b_sum[vertex_cloth_No][vertex_index].data(),
			target_position->collision_energy[vertex_cloth_No], vertex_mesh->vertex_position[vertex_index].data(), vertex_target_pos, vertex_collision_stiffness[vertex_index]);
		return true;
	}
	return false;
}

void Collision::addTargetPosToSystem(double* b_sum, double& energy, double* current_pos, double* target_pos, double stiffness)
{
	double temp[3];
	SUB(temp, current_pos, target_pos);
	energy += 0.5 * stiffness * DOT(temp, temp);
	b_sum[0] += target_pos[0] * stiffness;
	b_sum[1] += target_pos[2] * stiffness;
	b_sum[2] += target_pos[1] * stiffness;
}
void Collision::test()
{
	//findAllNeighborPairs();

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

		for (int k = i; k < cloth->size(); ++k) {
			compared_edge_aabb= &(*cloth)[k].edge_AABB;
			representative_edge_num = &(*cloth)[k].representative_edge_num;
			face = &(*cloth)[k].mesh_struct.faces;
			if (i != k) {
				for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
					edge_neighbor_one_cloth_edge = &(*edge_neighbor_edge)[j][k];
					triangle_neighbor_one_cloth_traingle = &(*triangle_neighbor_traingle)[(*edge_from_rep_triangle_index)[j]][k];
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
					edge_neighbor_one_cloth_edge->clear();
					edge_neighbor_one_cloth_edge->reserve(triangle_neighbor_one_cloth_traingle->size());
					for (int m = 0; m < triangle_neighbor_one_cloth_traingle->size(); ++m) {
						face_index = (*triangle_neighbor_one_cloth_traingle)[m];
						for (int n = 0; n < (*representative_edge_num)[face_index]; ++n) {
							compare_edge_index = (*face)[face_index].edge[n];
							if (j < compare_edge_index) {
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
