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
	use_BVH = false;
	initialNeighborPrimitive();
	collision_time_thread.resize(thread->thread_num);
	//collision_constraint.testPT();
	//approx_CCD.test();
	//CCD::test();
	//std::cout <<"floor coordinate "<< (*collider)[0].ori_vertices[0][1] << std::endl;
}

void Collision::initialDHatTolerance(double ave_edge_length)
{
	d_hat =1e-2 * ave_edge_length;
	eta = 0.01;
	tolerance = 1e-3 * d_hat;
	tolerance_2 = tolerance* tolerance;
	//std::cout << "d_hat_2 " << d_hat_2 << std::endl;
	
}

void Collision::initialNeighborPrimitive()
{
	for (int i = 0; i < cloth->size(); ++i) {
		(*cloth)[i].initialNeighborPrimitiveRecording(cloth->size(), tetrahedron->size(), collider->size(),use_BVH);
	}
	if (!use_BVH) {
		for (int i = 0; i < collider->size(); ++i) {
			(*collider)[i].initialNeighborPrimitiveRecording(cloth->size(), tetrahedron->size());
		}		
	}
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

	tet_target_pos_per_thread.resize(thread_num);
	for (int i = 0; i < tet_target_pos_per_thread.size(); ++i) {
		tet_target_pos_per_thread[i].initialSet(tetrahedron->size());
		for (int j = 0; j < tetrahedron->size(); ++j) {
			tet_target_pos_per_thread[i].initialSet2(j, (*tetrahedron)[j].ori_vertices.size());
		}
	}
	tet_target_pos.initialSet(tetrahedron->size());
	for (int j = 0; j < tetrahedron->size(); ++j) {
		tet_target_pos.initialSet2(j, (*tetrahedron)[j].ori_vertices.size());
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
	//time_t t1 = clock();
	if (use_BVH) {
		//for (int i = 0; i < 1000; ++i) {
			buildBVH();
		//}
	}
	else {		
		//for (int i = 0; i < 1000; ++i) {
			spatial_hashing.buildSpatialHashing();
		//}				
	}
	////std::cout << "build " << clock() - t1 << std::endl;
	//t1 = clock();
	//for (int i = 0; i < 1000; ++i) {
		thread->assignTask(this, FIND_TRIANGLE_PAIRS);
	//}
	////std::cout << "find triangle pair " << clock() - t1 << std::endl;

	//testCollision();

	//t1 = clock();
	//for (int i = 0; i < 1000; ++i) {
		thread->assignTask(this, FIND_PRIMITIVE_AROUND);

	//}
	////std::cout << "find around primitive " << clock() - t1 << std::endl;
	////std::cout << "search " << clock() - t1 << std::endl;
	thread->assignTask(this, GLOBAL_COLLISION_DETECTION);
	sumTargetPosition();
}

void Collision::testCulling()
{
	for (int i = 0; i < (*collider)[0].triangle_AABB.size(); ++i) {
		for (int j = 0; j < (*collider)[0].triangle_neighbor_cloth_vertex[i][0].size(); ++j) {
			if ((*collider)[0].triangle_neighbor_cloth_vertex[i][0][j] == 13) {
				std::cout << i << std::endl;
				}
			} 		
	}
}

void Collision::collisionCulling()
{
	getAABB();
	if (use_BVH) {
		buildBVH();
	}
	else {
		spatial_hashing.buildSpatialHashing();			
	}
	//for (int i = 0; i < collider->size(); ++i) {
	//	thread->assignTask(&(*collider)[i].mesh_struct, FACE_NORMAL);
	//}
	//for (int i = 0; i < cloth->size(); ++i) {
	//	thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	//}
	//for (int i = 0; i < tetrahedron->size(); ++i) {
	//	thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL);
	//}
	thread->assignTask(this, FIND_TRIANGLE_PAIRS);
	thread->assignTask(this, FIND_PRIMITIVE_AROUND);
	//testCulling();
}


void Collision::globalCollisionTime()
{
	for (int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL);
	}

	thread->assignTask(this, GLOBAL_COLLISION_TIME);
	//std::cout << "finish" << std::endl;
	collision_time = collision_time_thread[0];
	for (int i = 1; i < thread_num; ++i) {
		if (collision_time > collision_time_thread[i]) {
			collision_time = collision_time_thread[i];
		}
	}
	collision_time *= 0.9;

	if (collision_time > 1.0) {
		collision_time = 1.0;
	}

	//std::cout <<"collision time "<< collision_time << std::endl;
}


void Collision::solveCollisionConstraint()
{
	for (int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
	thread->assignTask(this, COLLISION_CONSTRAINT);
	sumTargetPosition();

	testIfBuildCollisionConstraint();
}


void Collision::testIfBuildCollisionConstraint()
{
	for (int i = 0; i < cloth_target_pos.b_sum[0].size(); ++i) {
		if (cloth_target_pos.need_update[0][i]) {
			//std::cout << "build collision constraint "<<i << std::endl;
		}
	}

}

//COLLISION_CONSTRAINT
void Collision::collisionConstraint(int thread_No)
{
	TargetPosition* target_pos = &cloth_target_pos_per_thread[thread_No];
	target_pos->initial();
	TriangleMeshStruct* mesh_struct;
	std::vector<std::vector<int>>* neighbor_primitve;
	int index_begin;
	int index_end;
	double* mass;
	std::array<double, 3>* current_pos;
	std::array<double, 3>* initial_pos;
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		mesh_struct = &(*cloth)[cloth_No].mesh_struct;
		index_begin = mesh_struct->vertex_index_begin_per_thread[thread_No];
		index_end = mesh_struct->vertex_index_begin_per_thread[thread_No + 1];
		neighbor_primitve = (*cloth)[cloth_No].vertex_neighbor_cloth_traingle.data();
		mass = mesh_struct->mass.data();
		current_pos = mesh_struct->vertex_position.data();
		initial_pos = mesh_struct->vertex_for_render.data();
		for (int i = index_begin; i < index_end; ++i) {
			pointSelfTriangleClose(neighbor_primitve[i].data(), initial_pos[i].data(), current_pos[i].data(), i, cloth_No, mass[i], target_pos);
		}
	}
	std::array<int, 3>* triangle_vertex_index;
	std::vector<double*> triangle_pos(3);
	std::array<double, 3>* triangle_normal;
	for (int collider_No = 0; collider_No < collider->size(); ++collider_No) {
		mesh_struct = &(*collider)[collider_No].mesh_struct;
		index_begin = mesh_struct->face_index_begin_per_thread[thread_No];
		index_end = mesh_struct->face_index_begin_per_thread[thread_No + 1];
		current_pos = mesh_struct->vertex_position.data();
		neighbor_primitve = (*collider)[collider_No].triangle_neighbor_cloth_vertex.data();
		triangle_vertex_index = mesh_struct->triangle_indices.data();
		triangle_normal = mesh_struct->face_normal.data();
		for (int i = index_begin; i < index_end; ++i) {
			for (int j = 0; j < 3; ++j) {
				triangle_pos[j] = current_pos[triangle_vertex_index[i][j]].data();
			}
			pointColliderTriangleClose(triangle_vertex_index[i].data(), neighbor_primitve[i].data(), triangle_pos, triangle_normal[i].data(), target_pos);
		}
	}
	int* edge_vertex_index;
	double mass_[4];
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		mesh_struct = &(*cloth)[cloth_No].mesh_struct;
		index_begin = mesh_struct->edge_index_begin_per_thread[thread_No];
		index_end = mesh_struct->edge_index_begin_per_thread[thread_No + 1];
		neighbor_primitve = (*cloth)[cloth_No].edge_neighbor_cloth_edge.data();
		initial_pos = mesh_struct->vertex_for_render.data();
		current_pos = mesh_struct->vertex_position.data();
		for (int i = index_begin; i < index_end; ++i) {
			edge_vertex_index = mesh_struct->edges[i].vertex;
			mass_[0] = mesh_struct->mass[edge_vertex_index[0]];
			mass_[1] = mesh_struct->mass[edge_vertex_index[1]];
			edgeEdgeClose(neighbor_primitve[i].data(), initial_pos[edge_vertex_index[0]].data(), initial_pos[edge_vertex_index[1]].data(),
				current_pos[edge_vertex_index[0]].data(), current_pos[edge_vertex_index[1]].data(), cloth_No, edge_vertex_index[0],
				edge_vertex_index[1], mass_, target_pos);
		}
	}

	target_pos== &tet_target_pos_per_thread[thread_No];
	target_pos->initial();
}




//GLOBAL_COLLISION_TIME
void Collision::collisionTime(int thread_No)
{
	int thread_test=2;
	double* collision_time = &collision_time_thread[thread_No];
	(*collision_time) = 2.0;
	TriangleMeshStruct* mesh_struct;
	std::vector<std::vector<int>>* neighbor_primitve;
	int index_begin;
	int index_end;
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		mesh_struct = &(*cloth)[cloth_No].mesh_struct;
		index_begin = mesh_struct->vertex_index_begin_per_thread[thread_No];
		index_end= mesh_struct->vertex_index_begin_per_thread[thread_No+1];
		neighbor_primitve = (*cloth)[cloth_No].vertex_neighbor_cloth_traingle.data();
		for (int i = index_begin; i < index_end; ++i) {
			pointSelfTriangleCollisionTime(collision_time, neighbor_primitve[i].data(),mesh_struct->vertex_for_render[i].data(), mesh_struct->vertex_position[i].data(),i);
		}
	}
	std::array<double,3>* initial_pos;
	std::array<double,3>* current_pos;
	for (int collider_No = 0; collider_No < collider->size(); ++collider_No) {
		mesh_struct = &(*collider)[collider_No].mesh_struct;
		index_begin = mesh_struct->face_index_begin_per_thread[thread_No];
		index_end = mesh_struct->face_index_begin_per_thread[thread_No+1];		
		initial_pos = mesh_struct->vertex_for_render.data();
		current_pos = mesh_struct->vertex_position.data();
		neighbor_primitve = (*collider)[collider_No].triangle_neighbor_cloth_vertex.data();
		for (int i = index_begin; i < index_end; ++i) {			
			pointColliderTriangleCollisionTime(collision_time, mesh_struct->triangle_indices[i].data(), neighbor_primitve[i].data(), initial_pos,current_pos, 
				mesh_struct->ori_face_normal_for_render[i].data(), mesh_struct->ori_face_normal[i].data(), mesh_struct->cross_for_approx_CCD[i].data(),
				mesh_struct->f_face_normal_for_render[i].data(),mesh_struct->f_face_normal[i].data(),mesh_struct->f_cross_for_approx_CCD[i].data(),i);
			
		}

	}	
	
	int* edge_vertex_index;
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		mesh_struct = &(*cloth)[cloth_No].mesh_struct;
		index_begin = mesh_struct->edge_index_begin_per_thread[thread_No];
		index_end = mesh_struct->edge_index_begin_per_thread[thread_No+1];		
		neighbor_primitve = (*cloth)[cloth_No].edge_neighbor_cloth_edge.data();
		initial_pos= mesh_struct->vertex_for_render.data();
		current_pos = mesh_struct->vertex_position.data();
		for (int i = index_begin; i < index_end; ++i) {
			edge_vertex_index = mesh_struct->edges[i].vertex;
			edgeEdgeCollisionTime(collision_time, neighbor_primitve[i].data(), initial_pos[edge_vertex_index[0]].data(), initial_pos[edge_vertex_index[1]].data(), 
				current_pos[edge_vertex_index[0]].data(), current_pos[edge_vertex_index[1]].data());
		}
	}
}


void Collision::testCollision()
{
	int k = 0;
	for (int i = 0; i < (*cloth)[0].triangle_neighbor_cloth_triangle.size(); ++i) {
		//(*cloth)[0].triangle_neighbor_cloth_triangle[i][0].push_back(i);
		k += (*cloth)[0].triangle_neighbor_cloth_triangle[i][0].size();
	}
	//std::cout << k + (*cloth)[0].triangle_neighbor_cloth_triangle.size()<<std::endl;

	//for (int i = 0; i < (*cloth)[0].mesh_struct.vertices.size(); ++i) {
	//int i = 0;

	//for (int j = 0; j < (*cloth)[0].vertex_neighbor_collider_triangle[i][0].size(); ++j) {
	//	//std::cout << (*cloth)[0].vertex_AABB[i].min[0] << " " << (*cloth)[0].vertex_AABB[i].min[1] << " " << (*cloth)[0].vertex_AABB[i].min[2] << " "
	//		<< (*cloth)[0].vertex_AABB[i].max[0] << " " << (*cloth)[0].vertex_AABB[i].max[1] << " " << (*cloth)[0].vertex_AABB[i].max[2] << std::endl;
	//	//std::cout << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].min[0] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].min[1] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].min[2] << " "
	//		<< (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].max[0] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].max[1] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].max[2] << std::endl;
	//	//std::cout << i << " " << (*cloth)[0].vertex_neighbor_collider_triangle[i][0][j] << std::endl;
	//	/*	if (!(*cloth)[0].vertex_AABB[i].AABB_intersection((*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]])) {
	//			//std::cout << i << " " << (*cloth)[0].vertex_neighbor_collider_triangle[i][0][j] << std::endl;
	//		}*/
	//}

	//if (!(*cloth)[0].vertex_neighbor_collider_triangle[i][0].empty()) {
	//	system("pause");
	//}
	//}
	//for (int i = 0; i < (*collider)[0].mesh_struct.triangle_indices.size(); ++i) {
	//	for (int j = 0; j < (*collider)[0].triangle_neighbor_cloth_vertex[i][0].size(); ++j) {
	//		if ((*collider)[0].triangle_neighbor_cloth_vertex[i][0][j] == 0) {
	//			if (!(*collider)[0].triangle_AABB[i].AABB_intersection((*cloth)[0].vertex_AABB[(*collider)[0].triangle_neighbor_cloth_vertex[i][0][j]])) {
	//				//std::cout << (*collider)[0].triangle_neighbor_cloth_vertex[i][0][j] << " " << i << std::endl;
	//			}
	//			system("pause");
	//		}
	//	}
	//}
	//if (!(*cloth)[0].vertex_neighbor_collider_triangle[i][0].empty()) {
	//	
	//}
	//for (int i = 0; i < (*cloth)[0].mesh_struct.triangle_indices.size(); ++i) {
	//	for (int j = 0; j < (*cloth)[0].triangle_neighbor_collider_triangle[i][0].size(); ++j) {
	//		if (!(*cloth)[0].triangle_AABB[i].AABB_intersection((*collider)[0].triangle_AABB[(*cloth)[0].triangle_neighbor_collider_triangle[i][0][j]])) {
	//			//std::cout << i << " " << (*cloth)[0].triangle_neighbor_collider_triangle[i][0][j] << std::endl;
	//		}
	//		
	//	}
	//}
	
}

void Collision::updateCollisionPosition()
{
	for (int i = 0; i < cloth->size(); ++i) {
		(*cloth)[i].mesh_struct.getNormal();
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		(*tetrahedron)[i].mesh_struct.getNormal();
	}
	thread->assignTask(this, RE_DETECTION);
	resumTargetPosition();
}



void Collision::resumTargetPosition()
{
	cloth_target_pos.partialInitial();
	tet_target_pos.partialInitial();
	thread->assignTask(this, RESUM_TARGET_POSITION);
	for (int i = 0; i < thread_num; ++i) {
		for (int j = 0; j < cloth->size(); ++j) {
			cloth_target_pos.collision_energy[j] += cloth_target_pos_per_thread[i].collision_energy[j];
		}
		for (int j = 0; j < tetrahedron->size(); ++j) {
			tet_target_pos.collision_energy[j] += tet_target_pos_per_thread[i].collision_energy[j];
		}
	}
}

void Collision::sumTargetPosition()
{
	cloth_target_pos.initial();
	tet_target_pos.initial();
	for (int i = 0; i < thread_num; ++i) {
		for (int j = 0; j < cloth->size(); ++j) {
			cloth_target_pos.collision_energy[j] += cloth_target_pos_per_thread[i].collision_energy[j];
		}
		for (int j = 0; j < tetrahedron->size(); ++j) {
			tet_target_pos.collision_energy[j] += tet_target_pos_per_thread[i].collision_energy[j];
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

	for (int j = 0; j < tetrahedron->size(); ++j) {
		b_sum = &tet_target_pos.b_sum[j];
		index_begin = &(*tetrahedron)[j].mesh_struct.vertex_index_begin_per_thread;
		for (int i = 0; i < thread_num; ++i) {
			need_update = tet_target_pos_per_thread[i].need_update[j];
			b_sum_per_thread = &tet_target_pos_per_thread[i].b_sum[j];
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
					SUM_((*b_sum)[k],(*b_sum_per_thread)[k]);
					global_stiffness[k] += stiffness[k];
				}
			}
		}
	}
	for (int j = 0; j < tetrahedron->size(); ++j) {
		global_need_update = tet_target_pos.need_update[j];
		b_sum = &tet_target_pos.b_sum[j];
		index_begin = &(*tetrahedron)[j].mesh_struct.vertex_index_begin_per_thread;
		global_stiffness = tet_target_pos.stiffness[j].data();
		for (int i = 0; i < thread_num; ++i) {
			need_update = tet_target_pos_per_thread[i].need_update[j];
			b_sum_per_thread = &tet_target_pos_per_thread[i].b_sum[j];
			stiffness = tet_target_pos_per_thread[i].stiffness[j].data();
			for (int k = (*index_begin)[thread_id]; k < (*index_begin)[thread_id + 1]; ++k) {
				if (need_update[k]) {
					global_need_update[k] = true;
					SUM_((*b_sum)[k], (*b_sum_per_thread)[k]);
					global_stiffness[k] += stiffness[k];
				}
			}
		}
	}
	////std::cout << "collision vertex ";
	//for (int i = 0; i < (*cloth)[0].ori_vertices.size(); ++i) {
	//	if (cloth_target_pos.need_update[0][i]) {
	//		//std::cout << i << " ";
	//	}
	//}
	////std::cout << std::endl;
}


//GLOBAL_COLLISION_DETECTION
void Collision::collisionDetection(int thread_No)
{
	TargetPosition* target_pos = &cloth_target_pos_per_thread[thread_No];
	target_pos->initial();
	int* index_begin;
	if (use_BVH) {
		for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
			index_begin = (*cloth)[cloth_No].mesh_struct.vertex_index_begin_per_thread.data();
			for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
				pointSelfTriangleCollisionDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].vertex_neighbor_cloth_traingle[i],
					&(*cloth)[cloth_No].collide_vertex_cloth_triangle[i], &(*cloth)[cloth_No].mesh_struct, (*cloth)[cloth_No].PC_radius[SELF_POINT_TRIANGLE], target_pos);
				pointColliderTriangleCollisionDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].vertex_neighbor_collider_triangle[i],
					&(*cloth)[cloth_No].collide_vertex_collider_triangle[i], &(*cloth)[cloth_No].mesh_struct, (*cloth)[cloth_No].PC_radius[BODY_POINT_TRIANGLE], target_pos);
			}
		}
	}
	else {
		for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
			index_begin = (*cloth)[cloth_No].mesh_struct.vertex_index_begin_per_thread.data();
			for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
				pointSelfTriangleCollisionDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].vertex_neighbor_cloth_traingle[i],
					&(*cloth)[cloth_No].collide_vertex_cloth_triangle[i], &(*cloth)[cloth_No].mesh_struct, (*cloth)[cloth_No].PC_radius[SELF_POINT_TRIANGLE], target_pos);
			}
		}
		for (int collider_No = 0; collider_No < collider->size(); ++collider_No) {
			index_begin = (*collider)[collider_No].mesh_struct.face_index_begin_per_thread.data();
			for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
				colliderTriangleVertexCollisionDetection(thread_No, i, collider_No, &(*collider)[collider_No].triangle_neighbor_cloth_vertex[i],
					&(*collider)[collider_No].collider_triangle_cloth_vertex[i], &(*collider)[collider_No].mesh_struct, (*collider)[collider_No].tolerance, target_pos);
			}
		}
	}
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		index_begin = (*cloth)[cloth_No].mesh_struct.edge_index_begin_per_thread.data();
		for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
			edgeSelfEdgeCollisionDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].edge_neighbor_cloth_edge[i], &(*cloth)[cloth_No].collide_edge_cloth_edge[i], 
				&(*cloth)[cloth_No].mesh_struct, (*cloth)[cloth_No].PC_radius[SELF_EDGE_EDGE], target_pos);
		}
	}

	target_pos = &tet_target_pos_per_thread[thread_No];
	target_pos->initial();
}

void Collision::colliderTriangleVertexCollisionDetection(int thread_No, int triangle_index, int collider_No,
	std::vector<std::vector<int>>* triangle_neighbor_vertex, std::vector<std::vector<int>>* collide_triangle_vertex, MeshStruct* triangle_mesh, double radius0,
	TargetPosition* target_pos)
{
	MeshStruct* vertex_mesh;
	std::vector<int>* neighbor_vertex;
	std::vector<int>* collide_vertex;
	double radius1;
	for (int i = 0; i < cloth->size(); ++i) {
		vertex_mesh = &(*cloth)[i].mesh_struct;
		neighbor_vertex = &(*triangle_neighbor_vertex)[i];
		collide_vertex = &(*collide_triangle_vertex)[i];
		collide_vertex->clear();
		collide_vertex->reserve(neighbor_vertex->size());
		radius1 = (*cloth)[i].PC_radius[BODY_POINT_TRIANGLE];
		for (int k = 0; k < neighbor_vertex->size(); ++k) {
			if (checkPointColliderTriangleCollision(vertex_mesh, triangle_mesh, radius0 + radius1, (*neighbor_vertex)[k], triangle_index,  i,
				collider_No, target_pos, true, (*cloth)[i].collision_stiffness[BODY_POINT_TRIANGLE])) {
				collide_vertex->push_back((*neighbor_vertex)[k]);
			}
		}
	}
}

//RE_DETECTION
void Collision::collisionReDetection(int thread_No)
{
	TargetPosition* target_pos = &cloth_target_pos_per_thread[thread_No];
	target_pos->partialInitial();
	int* index_begin;
	if (use_BVH) {
		for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
			index_begin = (*cloth)[cloth_No].mesh_struct.vertex_index_begin_per_thread.data();
			for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
				pointSelfTriangleCollisionReDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].collide_vertex_cloth_triangle[i], &(*cloth)[cloth_No].mesh_struct,
					(*cloth)[cloth_No].PC_radius[SELF_POINT_TRIANGLE], &(*cloth)[cloth_No].collision_stiffness[SELF_POINT_TRIANGLE], target_pos);
			}
			for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
				pointColliderTriangleCollisionReDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].collide_vertex_collider_triangle[i], &(*cloth)[cloth_No].mesh_struct,
					(*cloth)[cloth_No].PC_radius[BODY_POINT_TRIANGLE], &(*cloth)[cloth_No].collision_stiffness[BODY_POINT_TRIANGLE], target_pos);
			}
		}
	}
	else {
		for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
			index_begin = (*cloth)[cloth_No].mesh_struct.vertex_index_begin_per_thread.data();
			for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
				pointSelfTriangleCollisionReDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].collide_vertex_cloth_triangle[i], &(*cloth)[cloth_No].mesh_struct,
					(*cloth)[cloth_No].PC_radius[SELF_POINT_TRIANGLE], &(*cloth)[cloth_No].collision_stiffness[SELF_POINT_TRIANGLE], target_pos);
			}
		}
		for (int collider_No = 0; collider_No < collider->size(); ++collider_No) {
			index_begin = (*collider)[collider_No].mesh_struct.face_index_begin_per_thread.data();
			for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
				colliderTriangleVertexCollisionReDetection(thread_No, i, collider_No,(*collider)[collider_No].collider_triangle_cloth_vertex[i].data(), 
					&(*collider)[collider_No].mesh_struct, (*collider)[collider_No].tolerance, target_pos);
			}
		}
	}
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		index_begin = (*cloth)[cloth_No].mesh_struct.edge_index_begin_per_thread.data();
		for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
			edgeSelfEdgeCollisionReDetection(thread_No, i, cloth_No, &(*cloth)[cloth_No].collide_edge_cloth_edge[i], &(*cloth)[cloth_No].mesh_struct,
				(*cloth)[cloth_No].PC_radius[SELF_EDGE_EDGE], (*cloth)[cloth_No].collision_stiffness[SELF_EDGE_EDGE], target_pos);
		}
	}

	target_pos = &tet_target_pos_per_thread[thread_No];
	target_pos->partialInitial();
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



void Collision::pointSelfTriangleClose(std::vector<int>* vertex_neighbor_triangle, double* initial_vertex_pos, double* current_vertex_pos,
	int vertex_index, int cloth_No, double mass, TargetPosition* target_position)
{
	TriangleMeshStruct* triangle_mesh;
	std::vector<int>* neighbor_triangle;
	int* triangle_vertex_index;
	std::array<double, 3>* current_position;
	std::array<double, 3>* initial_position;
	std::array<int, 3>* triangle_indices;
	std::array<double, 3>* initial_face_normal;
	int triangle_index;

	std::vector<double*>triangle_initial_pos(3);
	std::vector<double*>triangle_current_pos(3);

	double target_pos[3];
	std::vector<std::array<double, 3>> triangle_target_pos(3);
	double stiffness;
	double triangle_mass[3];
	double record_stiffness;
	double* vertex_b_sum = target_position->b_sum[cloth_No][vertex_index].data();
	double* vertex_stiffness = &target_position->stiffness[cloth_No][vertex_index];
	bool* vetex_need_update = &target_position->need_update[cloth_No][vertex_index];

	for (int i = 0; i < cloth->size(); ++i) {
		triangle_mesh = &(*cloth)[i].mesh_struct;
		neighbor_triangle = &(vertex_neighbor_triangle[i]);
		current_position = triangle_mesh->vertex_position.data();
		initial_position = triangle_mesh->vertex_for_render.data();
		triangle_indices = triangle_mesh->triangle_indices.data();
		initial_face_normal = triangle_mesh->face_normal_for_render.data();
		record_stiffness = (*cloth)[i].collision_stiffness_initial[SELF_POINT_TRIANGLE];
		for (int k = 0; k < neighbor_triangle->size(); ++k) {
			triangle_index = (*neighbor_triangle)[k];
			triangle_vertex_index = triangle_indices[triangle_index].data();

			for (int j = 0; j < 3; ++j) {
				triangle_initial_pos[j] = initial_position[triangle_vertex_index[j]].data();
				triangle_current_pos[j] = current_position[triangle_vertex_index[j]].data();
				triangle_mass[j] = triangle_mesh->mass[triangle_vertex_index[j]];
			}
			stiffness = record_stiffness;
			if (collision_constraint.pointSelfTriangle(initial_vertex_pos, current_vertex_pos, triangle_initial_pos, triangle_current_pos,
				initial_face_normal[triangle_index].data(), target_pos, triangle_target_pos.data(), d_hat, stiffness,mass, triangle_mass)) {
				addTargetPosToSystemTotal(vertex_b_sum, target_position->collision_energy[cloth_No], initial_vertex_pos,
					target_pos, stiffness, *vertex_stiffness, *vetex_need_update);
				int triangle_vertex;
				for (int j = 0; j < 3; ++j) {
					triangle_vertex = triangle_vertex_index[j];
					addTargetPosToSystemTotal(target_position->b_sum[i][triangle_vertex].data(), target_position->collision_energy[i], initial_position[triangle_vertex].data(),
						triangle_target_pos[j].data(), stiffness,
						target_position->stiffness[i][triangle_vertex], target_position->need_update[i][triangle_vertex]);
					//if (vertex_index == 1) {
					//	//std::cout << triangle_initial_pos[j][0] << " " << triangle_initial_pos[j][1] << " " << triangle_initial_pos[j][2] << std::endl;
					//	//std::cout << triangle_current_pos[j][0] << " " << triangle_current_pos[j][1] << " " << triangle_current_pos[j][2] << std::endl;
					//	//std::cout << triangle_target_pos[j][0] << " " << triangle_target_pos[j][1] << " " << triangle_target_pos[j][2] << std::endl;
					//	//std::cout << initial_face_normal[triangle_index][0] << " " << initial_face_normal[triangle_index][1] << " " << initial_face_normal[triangle_index][2] << std::endl;
					//}
					//
					////std::cout << target_position->collision_energy[i] << std::endl;
				}
				
			}
		}
	}
}

void Collision::edgeEdgeClose(std::vector<int>* edge_neighbor_edge, double* initial_edge_vertex_0, double* initial_edge_vertex_1, double* current_edge_vertex_0, double* current_edge_vertex_1,
	int cloth_No, int edge_vertex_index_0, int edge_vertex_index_1, double* mass, TargetPosition* target_position)
{
	TriangleMeshStruct* compare_mesh;
	std::vector<int>* neighbor_edge;
	double current_collision_time;
	std::array<double, 3>* current_position;
	std::array<double, 3>* initial_position;
	std::vector<std::array<double, 3>> target_pos(2);
	std::vector<std::array<double, 3>> compare_target_pos(2);
	int compare_edge_index_0;
	int compare_edge_index_1;
	double stiffness;
	double record_stiffness;
	for (int i = 0; i < cloth->size(); ++i) {
		compare_mesh = &(*cloth)[i].mesh_struct;
		neighbor_edge = &(edge_neighbor_edge[i]);
		current_position = compare_mesh->vertex_position.data();
		initial_position = compare_mesh->vertex_for_render.data();
		record_stiffness = (*cloth)[i].collision_stiffness_initial[SELF_EDGE_EDGE];
		for (int k = 0; k < neighbor_edge->size(); ++k) {
			compare_edge_index_0= compare_mesh->edges[(*neighbor_edge)[k]].vertex[0];
			compare_edge_index_1= compare_mesh->edges[(*neighbor_edge)[k]].vertex[1];
			mass[2] = compare_mesh->mass[compare_edge_index_0];
			mass[3] = compare_mesh->mass[compare_edge_index_1];
			stiffness = record_stiffness;
			if (collision_constraint.edgeEdgeCollision(target_pos, compare_target_pos,current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0, 
				initial_edge_vertex_1, current_position[compare_edge_index_0].data(), current_position[compare_edge_index_1].data(),
				initial_position[compare_edge_index_0].data(), initial_position[compare_edge_index_1].data(),mass,d_hat, stiffness)) {
				addTargetPosToSystemTotal(target_position->b_sum[cloth_No][edge_vertex_index_0].data(),
					target_position->collision_energy[cloth_No], initial_edge_vertex_0,
					target_pos[0].data(), stiffness, target_position->stiffness[cloth_No][edge_vertex_index_0], target_position->need_update[cloth_No][edge_vertex_index_0]);
				addTargetPosToSystemTotal(target_position->b_sum[cloth_No][edge_vertex_index_1].data(),
					target_position->collision_energy[cloth_No], initial_edge_vertex_1,
					target_pos[1].data(), stiffness, target_position->stiffness[cloth_No][edge_vertex_index_1], target_position->need_update[cloth_No][edge_vertex_index_1]);
				addTargetPosToSystemTotal(target_position->b_sum[i][compare_edge_index_0].data(),
					target_position->collision_energy[i], initial_position[compare_edge_index_0].data(),
					compare_target_pos[0].data(), stiffness, target_position->stiffness[i][compare_edge_index_0], target_position->need_update[i][compare_edge_index_0]);
				addTargetPosToSystemTotal(target_position->b_sum[i][compare_edge_index_1].data(),
					target_position->collision_energy[i], initial_position[compare_edge_index_1].data(),
					compare_target_pos[1].data(), stiffness, target_position->stiffness[i][compare_edge_index_1], target_position->need_update[i][compare_edge_index_1]);
			}
		}
	}
}



void Collision::pointColliderTriangleClose(int* triangle_vertex_index, std::vector<int>* triangle_neighbor_vertex,
	std::vector<double*>& current_position, double* current_face_normal, TargetPosition* target_position)
{
	MeshStruct* vertex_mesh;
	std::vector<int>* neighbor_vertex;
	double current_collision_time;
	double target_pos[3];
	double stiffness;
	double record_stiffness;
	std::array<double,3>* vertex_b_sum;
	double* vertex_stiffness;
	bool* vertex_need_update;
	int vertex_index;
	double* energy;
	for (int i = 0; i < cloth->size(); ++i) {
		vertex_mesh = &(*cloth)[i].mesh_struct;
		neighbor_vertex = &(triangle_neighbor_vertex[i]);
		vertex_b_sum = target_position->b_sum[i].data();
		vertex_stiffness = target_position->stiffness[i].data();
		vertex_need_update = target_position->need_update[i];
		record_stiffness = (*cloth)[i].collision_stiffness_initial[BODY_POINT_TRIANGLE];
		energy = &target_position->collision_energy[i];
		for (int k = 0; k < neighbor_vertex->size(); ++k) {
			vertex_index = (*neighbor_vertex)[k];
			stiffness = record_stiffness;
			if (collision_constraint.pointColliderTriangle(vertex_mesh->vertex_for_render[vertex_index].data(), vertex_mesh->vertex_position[vertex_index].data(),
				current_position, current_face_normal, target_pos, d_hat, stiffness)) {
				addTargetPosToSystemTotal(vertex_b_sum[vertex_index].data(), 
					*energy, vertex_mesh->vertex_for_render[vertex_index].data(),
					target_pos, stiffness,
					vertex_stiffness[vertex_index], vertex_need_update[vertex_index]);
				////std::cout << vertex_index << " " << stiffness << std::endl;
			}
		}
	}
}



void Collision::pointSelfTriangleCollisionTime(double* collision_time, std::vector<int>* vertex_neighbor_triangle, double* initial_vertex_pos, double* current_vertex_pos, int vertex_index)
{
	MeshStruct* triangle_mesh;
	std::vector<int>* neighbor_triangle;
	double current_collision_time;
	int* triangle_vertex_index;
	std::array<double, 3>* current_position;
	std::array<double, 3>* initial_position;
	std::array<int, 3>* triangle_indices;
	std::array<double, 3>* current_ori_face_normal;
	std::array<double, 3>* initial_ori_face_normal;
	std::array<double, 3>* cross_for_CCD;
	std::array<floating, 3>* f_current_face_normal;
	std::array<floating, 3>* f_initial_face_normal;
	std::array<floating, 3>* f_cross_for_CCD;

	int triangle_index;
	for (int i = 0; i < cloth->size(); ++i) {
		triangle_mesh = &(*cloth)[i].mesh_struct;
		neighbor_triangle = &(vertex_neighbor_triangle[i]);
		current_position = triangle_mesh->vertex_position.data();
		initial_position = triangle_mesh->vertex_for_render.data();
		triangle_indices = triangle_mesh->triangle_indices.data();
		current_ori_face_normal = triangle_mesh->ori_face_normal.data();
		initial_ori_face_normal = triangle_mesh->ori_face_normal_for_render.data();
		cross_for_CCD = triangle_mesh->cross_for_approx_CCD.data();

		f_current_face_normal=triangle_mesh->f_face_normal.data();
		f_initial_face_normal =triangle_mesh->f_face_normal_for_render.data();
		f_cross_for_CCD= triangle_mesh->f_cross_for_approx_CCD.data();

		for (int k = 0; k < neighbor_triangle->size(); ++k) {
			triangle_index = (*neighbor_triangle)[k];
			triangle_vertex_index = triangle_indices[triangle_index].data();

			current_collision_time = CCD::pointTriangleCcd(initial_vertex_pos, initial_position[triangle_vertex_index[0]].data(), initial_position[triangle_vertex_index[1]].data(), initial_position[triangle_vertex_index[2]].data(),
				current_vertex_pos, current_position[triangle_vertex_index[0]].data(), current_position[triangle_vertex_index[1]].data(), current_position[triangle_vertex_index[2]].data(), eta, tolerance);
			if ((*collision_time) > current_collision_time) {
				(*collision_time) = current_collision_time;
			}
			//if (approx_CCD.pointTriangleCollisionTime(current_collision_time, initial_vertex_pos, current_vertex_pos, initial_position[triangle_vertex_index[0]].data(), current_position[triangle_vertex_index[0]].data(),
			//	initial_position[triangle_vertex_index[1]].data(), current_position[triangle_vertex_index[1]].data(), initial_position[triangle_vertex_index[2]].data(), current_position[triangle_vertex_index[2]].data(),
			//	initial_ori_face_normal[triangle_index].data(),
			//	current_ori_face_normal[triangle_index].data(), cross_for_CCD[triangle_index].data(), tolerance_2,
			//	f_initial_face_normal[triangle_index].data(), f_current_face_normal[triangle_index].data(),f_cross_for_CCD[triangle_index].data(),vertex_index)) {//
			//	//if (current_collision_time < 1e-4) {
			//	//	std::cout << current_collision_time << " " << vertex_index << " " << triangle_vertex_index[0] << " " << triangle_vertex_index[1] << " " << triangle_vertex_index[2] << std::endl;
			//	//}
			//	if ((*collision_time) > current_collision_time) {
			//		(*collision_time) = current_collision_time;
			//	}
			//}
		}
	}
}

void Collision::pointColliderTriangleCollisionTime(double* collision_time, int* triangle_vertex_index, std::vector<int>* triangle_neighbor_vertex,
	std::array<double, 3>* initial_position, std::array<double, 3>* current_position,
	double* initial_ori_face_normal, double* current_ori_face_normal, double* cross_for_CCD,
	floating* f_initial_normal, floating* f_current_normal, floating* f_cross_for_CCD, int triangle_index)
{
	MeshStruct* vertex_mesh;
	std::vector<int>* neighbor_vertex;
	double current_collision_time;
	for (int i = 0; i < cloth->size(); ++i) {
		vertex_mesh = &(*cloth)[i].mesh_struct;
		neighbor_vertex = &(triangle_neighbor_vertex[i]);	
		//if (*time_stamp == 15 && triangle_index == 9104)
		//{
		//	std::cout << "k11" << std::endl;
		//}
		for (int k = 0; k < neighbor_vertex->size(); ++k) {
			//if (*time_stamp == 15 && triangle_index == 9104 && k==3)
			//{
			//	std::cout << vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]][0]<<" "<< vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]][1]<<" "
			//		<< vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]][2] << std::endl;
			//	std::cout << vertex_mesh->vertex_position[(*neighbor_vertex)[k]][0] << " " << vertex_mesh->vertex_position[(*neighbor_vertex)[k]][1] << " "
			//		<< vertex_mesh->vertex_position[(*neighbor_vertex)[k]][2] << std::endl;
			//	std::cout << initial_position[triangle_vertex_index[0]][0] << " " << initial_position[triangle_vertex_index[0]][1] << " " <<
			//		initial_position[triangle_vertex_index[0]][2] << std::endl;
			//	std::cout << initial_position[triangle_vertex_index[1]][0] << " " << initial_position[triangle_vertex_index[1]][1] << " " <<
			//		initial_position[triangle_vertex_index[1]][2] << std::endl;
			//	std::cout << initial_position[triangle_vertex_index[2]][0] << " " << initial_position[triangle_vertex_index[2]][1] << " " <<
			//		initial_position[triangle_vertex_index[2]][2] << std::endl;
			//	std::cout << current_position[triangle_vertex_index[0]][0] << " " << current_position[triangle_vertex_index[0]][1] << " " <<
			//		current_position[triangle_vertex_index[0]][2] << std::endl;
			//	std::cout << current_position[triangle_vertex_index[1]][0] << " " << current_position[triangle_vertex_index[1]][1] << " " <<
			//		current_position[triangle_vertex_index[1]][2] << std::endl;
			//	std::cout << current_position[triangle_vertex_index[2]][0] << " " << current_position[triangle_vertex_index[2]][1] << " " <<
			//		current_position[triangle_vertex_index[2]][2] << std::endl;
			//}
			current_collision_time = CCD::pointTriangleCcd(vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]].data(), initial_position[triangle_vertex_index[0]].data(), initial_position[triangle_vertex_index[1]].data(), initial_position[triangle_vertex_index[2]].data(),
				vertex_mesh->vertex_position[(*neighbor_vertex)[k]].data(), current_position[triangle_vertex_index[0]].data(), current_position[triangle_vertex_index[1]].data(), current_position[triangle_vertex_index[2]].data(), eta, tolerance);
			if ((*collision_time) > current_collision_time) {
				(*collision_time) = current_collision_time;
			}
			//if (*time_stamp == 15 && triangle_index == 9104)
			//{
			//	std::cout << k << " " << (*neighbor_vertex)[k] << std::endl;
			//}
			/*if (approx_CCD.pointTriangleCollisionTime(current_collision_time, vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]].data(), vertex_mesh->vertex_position[(*neighbor_vertex)[k]].data(), initial_position[triangle_vertex_index[0]].data(), current_position[triangle_vertex_index[0]].data(),
				initial_position[triangle_vertex_index[1]].data(), current_position[triangle_vertex_index[1]].data(),
				initial_position[triangle_vertex_index[2]].data(), current_position[triangle_vertex_index[2]].data(),
				initial_ori_face_normal, current_ori_face_normal, cross_for_CCD, tolerance_2,
				f_initial_normal, f_current_normal, f_cross_for_CCD, (*neighbor_vertex)[k])) {
				if ((*collision_time) > current_collision_time) {					
					(*collision_time) = current_collision_time;
				}
			}*/
		}
		//if (*time_stamp == 15 && triangle_index == 9104)
		//{
		//	std::cout << "k22" << std::endl;
		//}
	}
}

void Collision::edgeEdgeCollisionTime(double* collision_time, std::vector<int>* edge_neighbor_edge, double* initial_edge_vertex_0, double* initial_edge_vertex_1, double* current_edge_vertex_0, double* current_edge_vertex_1)
{
	TriangleMeshStruct* compare_mesh;
	std::vector<int>* neighbor_edge;
	double current_collision_time;
	std::array<double, 3>* current_position;
	std::array<double, 3>* initial_position;
	int* edge_vertex_index;
	for (int i = 0; i < cloth->size(); ++i) {
		compare_mesh = &(*cloth)[i].mesh_struct;
		neighbor_edge = &(edge_neighbor_edge[i]);
		current_position = compare_mesh->vertex_position.data();
		initial_position = compare_mesh->vertex_for_render.data();
		for (int k = 0; k < neighbor_edge->size(); ++k) {
			edge_vertex_index = compare_mesh->edges[(*neighbor_edge)[k]].vertex;

			current_collision_time = CCD::edgeEdgeCcd(initial_edge_vertex_0, initial_edge_vertex_1, initial_position[edge_vertex_index[0]].data(), initial_position[edge_vertex_index[1]].data(),
				current_edge_vertex_0, current_edge_vertex_1, current_position[edge_vertex_index[0]].data(), current_position[edge_vertex_index[1]].data(), eta, tolerance);
			if ((*collision_time) > current_collision_time) {
				(*collision_time) = current_collision_time;
			}

			/*if (approx_CCD.edgeEdgeCollisionTime(current_collision_time, current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0, initial_edge_vertex_1, current_position[edge_vertex_index[0]].data(), current_position[edge_vertex_index[1]].data(),
				initial_position[edge_vertex_index[0]].data(), initial_position[edge_vertex_index[1]].data(), tolerance_2)) {
				if ((*collision_time) > current_collision_time) {
					(*collision_time) = current_collision_time;
				}
			}*/
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
		radius1 = (*cloth)[i].PC_radius[SELF_POINT_TRIANGLE]+ radius0;
		triangle_collision_stiffness = &(*cloth)[i].collision_stiffness[SELF_POINT_TRIANGLE];
		for (int k = 0; k < collide_triangle->size(); ++k) {
			if (!checkPointTriangleCollision(vertex_mesh, triangle_mesh, radius1, vertex_index, (*collide_triangle)[k], cloth_No,
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

void Collision::colliderTriangleVertexCollisionReDetection(int thread_No, int triangle_index, int collider_No,
	std::vector<int>* collide_triangle_vertex, MeshStruct* triangle_mesh, double radius0, TargetPosition* target_pos)
{
	MeshStruct* vertex_mesh;
	std::vector<int>* collide_vertex;
	double radius1;
	std::vector<double>* collision_stiffness;
	for (int i = 0; i < cloth->size(); ++i) {
		vertex_mesh = &(*cloth)[i].mesh_struct;
		collide_vertex = &collide_triangle_vertex[i];
		radius1 = (*cloth)[i].PC_radius[BODY_POINT_TRIANGLE]+radius0;
		collision_stiffness = &(*cloth)[i].collision_stiffness[BODY_POINT_TRIANGLE];
		for (int k = 0; k < collide_vertex->size(); ++k) {
			if (!checkPointColliderTriangleCollision(vertex_mesh, triangle_mesh, radius1, (*collide_vertex)[k], triangle_index, i,
				collider_No, target_pos, false, (*collision_stiffness))) {
				addTargetPosToSystem(target_pos->b_sum[i][(*collide_vertex)[k]].data(),
					target_pos->collision_energy[i], vertex_mesh->vertex_position[(*collide_vertex)[k]].data(), vertex_mesh->vertex_position[(*collide_vertex)[k]].data(), (*collision_stiffness)[(*collide_vertex)[k]]);
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
		triangle_mesh->face_normal_for_render[triangle_index].data(), triangle_mesh->face_normal[triangle_index].data(), vertex_target_pos, triangle_target_pos, radius, triangle_mesh->triangle_normal_magnitude_reciprocal[triangle_index],
		vertex_mesh->mass[vertex_index], triangle_mass)) {	
		//if (vertex_index == 0) {
		//	//std::cout << "++" << std::endl;
		//	//std::cout << vertex_index << " " << triangle_index << " " << vertex_mesh->vertex_for_render[vertex_index][0] << " " << vertex_mesh->vertex_for_render[vertex_index][1] << " " << vertex_mesh->vertex_for_render[vertex_index][2] << " " << std::endl;
		//	for (int i = 0; i < triangle_target_pos.size(); ++i) {
		//		//std::cout << vertex_index << " " << triangle_index << " " << initial_triangle_pos[i][0] << " " << initial_triangle_pos[i][1] << " " << initial_triangle_pos[i][2] << " " << std::endl;
		//	}
		//	//std::cout << vertex_mesh->vertex_position[vertex_index][0] << " " << vertex_mesh->vertex_position[vertex_index][1] << " " << vertex_mesh->vertex_position[vertex_index][2] << " " << std::endl;
		//	for (int i = 0; i < triangle_target_pos.size(); ++i) {
		//		//std::cout << current_triangle_pos[i][0] << " " << current_triangle_pos[i][1] << " " << current_triangle_pos[i][2] << " " << std::endl;
		//	}
		//	//std::cout << vertex_index << " " << triangle_index <<" "<< vertex_target_pos[0] << " " << vertex_target_pos[1] << " " << vertex_target_pos[2] << " " << std::endl;
		//	for (int i = 0; i < triangle_target_pos.size(); ++i) {
		//		//std::cout << vertex_index<<" "<< triangle_index << " " << triangle_target_pos[i][0] << " " << triangle_target_pos[i][1] << " " << triangle_target_pos[i][2] << " " << std::endl;
		//	}
		//	//std::cout << "++" << std::endl;
		//}

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
		triangle_mesh->face_normal_for_render[triangle_index].data(), triangle_mesh->face_normal[triangle_index].data(), vertex_target_pos, radius)) {
		if (new_collision_registration) {
			target_position->stiffness[vertex_cloth_No][vertex_index] += vertex_collision_stiffness[vertex_index];
			target_position->need_update[vertex_cloth_No][vertex_index] = true;
			//if (vertex_index == 0) {
			//	//std::cout << vertex_target_pos[0] << " " << vertex_target_pos[1] << " " << vertex_target_pos[2] << " "
			//		<< vertex_mesh->vertex_position[vertex_index][0] << " " << vertex_mesh->vertex_position[vertex_index][1] << " " << vertex_mesh->vertex_position[vertex_index][2] << std::endl;
			//	//std::cout << current_triangle_pos[0][1] << std::endl;
			//}
		}
		else {
			/*if (vertex_index == 0) {
				//std::cout << vertex_target_pos[0] << " " << vertex_target_pos[1] << " " << vertex_target_pos[2] << " "
					<< vertex_mesh->vertex_position[vertex_index][0] << " " << vertex_mesh->vertex_position[vertex_index][1] << " " << vertex_mesh->vertex_position[vertex_index][2] << std::endl;
			}*/
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
	b_sum[1] += target_pos[1] * stiffness;
	b_sum[2] += target_pos[2] * stiffness;
}

void Collision::addTargetPosToSystemTotal(double* b_sum, double& energy, double* current_pos, double* target_pos, double stiffness, double& sum_stiffness, bool& update)
{
	double temp[3];
	SUB(temp, current_pos, target_pos);
	energy += 0.5 * stiffness * DOT(temp, temp);
	b_sum[0] += target_pos[0] * stiffness;
	b_sum[1] += target_pos[1] * stiffness;
	b_sum[2] += target_pos[2] * stiffness;
	update = true;
	sum_stiffness += stiffness;
}

void Collision::test()
{
	//findAllNeighborPairs();

	//std::cout << "vertex " << std::endl;
	for (int i = 0; i < (*cloth)[0].vertex_neighbor_cloth_traingle[0][0].size(); ++i) {
		//std::cout << (*cloth)[0].vertex_neighbor_cloth_traingle[0][0][i] << " ";
	}
	//std::cout << std::endl;
	//std::cout << "triangle " << std::endl;
	for (int i = 0; i < (*cloth)[0].triangle_neighbor_cloth_triangle[0][0].size(); ++i) {
		//std::cout << (*cloth)[0].triangle_neighbor_cloth_triangle[0][0][i] << " ";
	}
	//std::cout << std::endl;
	//std::cout << "edge " << std::endl;
	for (int i = 0; i < (*cloth)[0].edge_neighbor_cloth_edge[0][0].size(); ++i) {
		//std::cout << (*cloth)[0].edge_neighbor_cloth_edge[0][0][i] << " ";
	}
	//std::cout << std::endl;
}



void Collision::searchTriangle(AABB& aabb, int compare_index, int cloth_No, std::vector<std::vector<int>>* cloth_neighbor_index, std::vector<std::vector<int>>* collider_neighbor_index)
{
	for (int i = cloth_No; i < cloth->size(); ++i) {
		(*cloth_neighbor_index)[i].clear();
		cloth_BVH[i].search(aabb, compare_index,i==cloth_No, &((*cloth_neighbor_index)[i]),1,0, (*cloth)[i].triangle_AABB.size());
	}
	for (int i = 0; i < collider->size(); ++i) {
		(*collider_neighbor_index)[i].clear();
		collider_BVH[i].search(aabb, compare_index, false, &((*collider_neighbor_index)[i]), 1, 0, (*collider)[i].triangle_AABB.size());
	}
}


//FIND_TRIANGLE_PAIRS
void Collision::findAllTrianglePairs(int thread_No)
{
	int* thread_begin;
	if (use_BVH) {
		for (int i = 0; i < cloth->size(); ++i) {
			thread_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
			for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
				searchTriangle((*cloth)[i].triangle_AABB[j], j, i, &(*cloth)[i].triangle_neighbor_cloth_triangle[j],
					&(*cloth)[i].triangle_neighbor_collider_triangle[j]);
			}
		}
	}
	else {
		for (int i = 0; i < cloth->size(); ++i) { 
			thread_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
			for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
				spatial_hashing.searchTriangle((*cloth)[i].triangle_AABB[j],i, j,(*cloth)[i].triangle_neighbor_cloth_triangle[j].data(),
					false,thread_No);
				//if (j == 0) {
				//	//std::cout << j << " ";
				//	for (int k = 0; k < (*cloth)[i].triangle_neighbor_cloth_triangle[j][0].size(); ++k) {
				//		//std::cout << (*cloth)[i].triangle_neighbor_cloth_triangle[j][0][k] << " ";
				//	}
				//	//std::cout << std::endl;					
				//}
			}
		}
		for (int i = 0; i < collider->size(); ++i) {
			thread_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread.data();
			for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
				spatial_hashing.searchTriangle((*collider)[i].triangle_AABB[j], i, j, (*collider)[i].triangle_neighbor_cloth_triangle[j].data(),
					true, thread_No);
/*				for (int k = 0; k < (*collider)[i].triangle_neighbor_cloth_triangle[j][0].size(); ++k) {
					//std::cout << (*collider)[i].triangle_neighbor_cloth_triangle[j][0][k] << " ";
				}
				//std::cout << std::endl;	*/								
			}
		}
	}
}

//FIND_PRIMITIVE_AROUND
void Collision::findPrimitivesAround(int thread_No)
{
	findClothTriangleAroundVertex(thread_No);
	if (use_BVH) {
		findColliderTriangleAroundVertex(thread_No);
	}
	else {
		findVertexAroundColliderTriangle(thread_No);
	}
	findEdgeAroundEdge(thread_No);
}

void Collision::findVertexAroundColliderTriangle(int thread_No)
{
	int* thread_begin;
	std::vector<std::vector<std::vector<int>>>* triangle_neighbor_triangle;
	std::vector<std::vector<std::vector<int>>>* triangle_neighbor_vertex; //triangle index near vertex
	int* rep_vertex_num; // record of vertex's representative triangle index

	std::vector<int>* triangle_neighbor_this_cloth_vertex; //vector to record the triangle index near vertex which triangle and vertex in same cloth
	std::vector<int>* triangle_neighbor_this_cloth_triangle;//vector to record the triangle index near triangle in same cloth

	std::vector<AABB>* triangle_aabb;
	std::vector<AABB>* compared_vertex_aabb;
	int cloth_triangle_index;

	MeshStruct::Face* faces;

	for (int i = 0; i < collider->size(); ++i) {
		thread_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread.data();
		triangle_aabb = &(*collider)[i].triangle_AABB;
		triangle_neighbor_vertex = &(*collider)[i].triangle_neighbor_cloth_vertex;
		triangle_neighbor_triangle = &(*collider)[i].triangle_neighbor_cloth_triangle;
		for (int k = 0; k < cloth->size(); ++k) {
			compared_vertex_aabb = &(*cloth)[k].vertex_AABB;
			rep_vertex_num = (*cloth)[k].representative_vertex_num.data();
			faces = (*cloth)[k].mesh_struct.faces.data();
			for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
				triangle_neighbor_this_cloth_vertex = &(*triangle_neighbor_vertex)[j][k];
				triangle_neighbor_this_cloth_triangle = &(*triangle_neighbor_triangle)[j][k];
				triangle_neighbor_this_cloth_vertex->clear();				
				triangle_neighbor_this_cloth_vertex->reserve(triangle_neighbor_this_cloth_triangle->size());
				for (int m = 0; m < triangle_neighbor_this_cloth_triangle->size(); ++m) {
					cloth_triangle_index = (*triangle_neighbor_this_cloth_triangle)[m];
					for (int n = 0; n < rep_vertex_num[cloth_triangle_index]; ++n) {
						if ((*triangle_aabb)[j].AABB_intersection((*compared_vertex_aabb)[faces[cloth_triangle_index].vertex[n]])) {
							triangle_neighbor_this_cloth_vertex->push_back(faces[cloth_triangle_index].vertex[n]);
						}
					}
				}
			}
		}
	}
}

void Collision::findClothTriangleAroundVertex(int thread_No)
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
		triangle_neighbor_traingle = &(*cloth)[i].triangle_neighbor_cloth_triangle;
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
							vertex_neighbor_this_cloth_traingle->push_back((*triangle_neighbor_this_cloth_traingle)[m]);
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
								vertex_neighbor_this_cloth_traingle->push_back((*triangle_neighbor_this_cloth_traingle)[m]);
							}
						}
					}
				}
			}
		}
	}
}


void Collision::findColliderTriangleAroundVertex(int thread_No)
{
	int* thread_begin;
	std::vector<std::vector<std::vector<int>>>* triangle_neighbor_traingle;
	std::vector<std::vector<std::vector<int>>>* vertex_neighbor_traingle; //triangle index near vertex
	std::vector<int>* vertex_from_rep_triangle_index; // record of vertex's representative triangle index

	std::vector<int>* vertex_neighbor_this_cloth_traingle; //vector to record the triangle index near vertex which triangle and vertex in same cloth
	std::vector<int>* triangle_neighbor_this_cloth_traingle;//vector to record the triangle index near triangle in same cloth

	std::vector<AABB>* vertex_aabb;
	std::vector<AABB>* compared_triangle_aabb;

	for (int i = 0; i < cloth->size(); ++i) {
		thread_begin = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread.data();
		vertex_from_rep_triangle_index = &(*cloth)[i].vertex_from_rep_triangle_index;
		vertex_aabb = &(*cloth)[i].vertex_AABB;
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
		triangle_neighbor_traingle = &(*cloth)[i].triangle_neighbor_cloth_triangle;
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
