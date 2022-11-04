#include"XPBD_IPC.h"

XPBD_IPC::XPBD_IPC()
{
	gravity_ = 9.8;
	sub_step_num = 1;
	iteration_number = 300;

	damping_coe = 0.0;

	perform_collision = false;
	max_iteration_number = 10;
	outer_max_iteration_number = 20;
	XPBD_constraint.epsilon_for_bending = 1e-10;

	velocity_damp = 0.995;
	energy_converge_ratio = 2e-3;

	min_inner_iteration = 10;
	min_outer_iteration = 1;



}
void XPBD_IPC::initial()
{
	resetExternalForce();
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
	}
}

void XPBD_IPC::reset()
{
	resetExternalForce();
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
	}
}

void XPBD_IPC::updateItrInfo(int* iteration_num)
{
	iteration_num[LOCAL_GLOBAL] = iteration_number;
	iteration_num[OUTER] = outer_itr_num;
	//outer_iteration_number = iteration_num[OUTER];
	//sub_step_num = iteration_num[OUTER];
	sub_time_step = time_step / (double)sub_step_num;
}


void XPBD_IPC::setForXPBD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
	Thread* thread, double* tolerance_ratio)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	sub_time_step = time_step / (double)sub_step_num;
	this->thread = thread;
	total_thread_num = thread->thread_num;

	total_obj_num = cloth->size() + tetrahedron->size();

	this->floor = floor;

	reorganzieDataOfObjects();
	initialVariable();
	initialClothBending();
	setConstraintIndex();
	//energy_per_thread.resize(thread->thread_num,0.0);
	if (perform_collision) {
		//collision.energy = energy_per_thread.data();
		collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio, XPBD_,true);
		collision.setCollisionFreeVertex(&record_gloabl_CCD_vertex_position);
		//collision.setParameter(&lambda_collision,lambda.data()+ constraint_index_start[3], collision_constraint_index_start.data(), damping_coe, sub_time_step);
	}

	setConvergeCondition();

}

void XPBD_IPC::setConvergeCondition()
{
	converge_condition_ratio = 1e-2;
	double edge_length = calEdgeLength();
	max_move_standard = converge_condition_ratio * edge_length;
	converge_condition_ratio = 1e-4;
	max_move_standard_inner_itr = converge_condition_ratio * edge_length;
}

void XPBD_IPC::initialClothBending()
{
	lbo_weight.resize(cloth->size());
	vertex_lbo.resize(cloth->size());
	rest_mean_curvature_norm.resize(cloth->size());
	//rest_Aq.resize(cloth->size());
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		XPBD_constraint.initial_LBO_EdgeCotWeight(cloth->data()[i].mesh_struct, lbo_weight[i], vertex_lbo[i], rest_mean_curvature_norm[i]);
	}
}

void XPBD_IPC::setConstraintIndex()
{
	constraint_index_start.resize(5); //bending, edge_length, ARAP, floor collision
	constraint_index_start[0] = 0;
	unsigned int constraint_number = 0;
	if (use_bending_based_on_vertex) {
		for (unsigned int i = 0; i < cloth->size(); ++i) {
			constraint_number += mesh_struct[i]->vertex_position.size();
		}
	}
	else {
		for (unsigned int i = 0; i < cloth->size(); ++i) {
			constraint_number += mesh_struct[i]->edge_vertices.size() >> 1;
		}
	}
	constraint_index_start[1] = constraint_number;
	constraint_number = 0;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		constraint_number += mesh_struct[i]->edges.size();
	}
	constraint_index_start[2] = constraint_number + constraint_index_start[1];
	constraint_number = 0;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		constraint_number += tetrahedron->data()[i].mesh_struct.indices.size();
	}
	constraint_index_start[3] = constraint_number + constraint_index_start[2];
	constraint_number = 0;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		constraint_number += mesh_struct[i]->vertex_position.size();
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		constraint_number += tetrahedron->data()[i].mesh_struct.vertex_index_on_sureface.size();
	}
	constraint_index_start[4] = constraint_number + constraint_index_start[3];
	lambda.resize(constraint_index_start[4]);

	constraint_number = 0;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		constraint_number += tetrahedron->data()[i].mesh_struct.vertex_position.size();
	}

	lambda_collision.reserve(constraint_index_start[1] + constraint_number);
	collision_constraint_index_start.resize(3);
	for (unsigned int i = 0; i < collision_constraint_index_start.size(); ++i) {
		collision_constraint_index_start[i].resize(total_thread_num + 1, 0);
	}
}


void XPBD_IPC::initialVariable()
{
	f_ext.resize(total_obj_num);
	velocity.resize(total_obj_num);
	gravity[0] = 0;
	gravity[1] = -gravity_;
	gravity[2] = 0;

	sn.resize(total_obj_num);


	residual.resize(total_obj_num);

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		f_ext[i].resize(mesh_struct[i]->vertex_position.size());
		velocity[i].resize(mesh_struct[i]->vertex_position.size());
		sn[i].resize(mesh_struct[i]->vertex_position.size());
		residual[i].resize(mesh_struct[i]->vertex_position.size());
	}
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
		memset(f_ext[i][0].data(), 0, 24 * f_ext[i].size());
		memset(sn[i][0].data(), 0, 24 * sn[i].size());
	}

}

void XPBD_IPC::reorganzieDataOfObjects()
{
	vertex_position.resize(total_obj_num);
	initial_vertex_position.resize(total_obj_num);
	mesh_struct.resize(total_obj_num);
	vertex_index_begin_per_thread.resize(total_obj_num);
	record_vertex_position.resize(total_obj_num);
	record_gloabl_CCD_vertex_position.resize(total_obj_num);

	//record_outer_vertex_position.resize(total_obj_num);
	//unfixed_vertex.resize(total_obj_num);
	triangle_indices.resize(total_obj_num);
	edge_vertices.resize(total_obj_num);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position.data();
		initial_vertex_position[i] = cloth->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i] = &cloth->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		record_gloabl_CCD_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		triangle_indices[i] = cloth->data()[i].mesh_struct.triangle_indices.data();
		//record_outer_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		//unfixed_vertex[i] = &cloth->data()[i].mesh_struct.unfixed_point_index;
		edge_vertices[i] = cloth->data()[i].mesh_struct.edge_vertices.data();
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		initial_vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		record_gloabl_CCD_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		triangle_indices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.triangle_indices.data();
		//record_outer_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		//unfixed_vertex[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.unfixed_point_index;
		edge_vertices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_vertices.data();
	}
	if (!collider->empty()) {
		triangle_indices_collider.resize(collider->size());
		collider_mesh_struct.resize(collider->size());
		collider_edge_vertices.resize(collider->size());
		vertex_position_collider.resize(collider->size());

		for (unsigned int i = 0; i < collider->size(); ++i) {
			collider_mesh_struct[i] = &collider->data()[i].mesh_struct;
			triangle_indices_collider[i] = collider->data()[i].mesh_struct.triangle_indices.data();
			collider_edge_vertices[i] = collider->data()[i].mesh_struct.edge_vertices.data();
			vertex_position_collider[i] = collider->data()[i].mesh_struct.vertex_position.data();
		}
	}
	

}

void XPBD_IPC::saveScene(double* force_direction, int obj_No, bool have_force)
{
	std::string file_name;
	save_scene.save_scene_XPBD(*time_stamp, *time_indicate_for_simu, mesh_struct, &velocity,
		collider_mesh_struct, file_name);
	if (have_force) {
		if (obj_No < cloth->size()) {
			save_scene.save_force(file_name, force_direction, cloth->data()[obj_No].coe_neighbor_vertex_force,
				cloth->data()[obj_No].neighbor_vertex, obj_No);
		}
		else {
			save_scene.save_force(file_name, force_direction, tetrahedron->data()[obj_No - cloth->size()].coe_neighbor_vertex_force,
				tetrahedron->data()[obj_No - cloth->size()].neighbor_vertex, obj_No);
		}
	}

}

void XPBD_IPC::readScene(const char* file_name, double* force_direction, int& obj_No)
{
	std::vector<double>force_coe;
	std::vector<int>neighbor_vertex_index;
	save_scene.read_scene_XPBD(file_name, time_stamp, time_indicate_for_simu, mesh_struct,
		&velocity, collider_mesh_struct, *has_force, force_direction, force_coe, neighbor_vertex_index, obj_No);
	for (unsigned int i = 0; i < mesh_struct.size(); ++i) {
		memcpy(mesh_struct[i]->vertex_for_render[0].data(), mesh_struct[i]->vertex_position[0].data(), 24 * mesh_struct[i]->vertex_position.size());
	}
	for (unsigned int i = 0; i < collider_mesh_struct.size(); ++i) {
		memcpy(collider_mesh_struct[i]->vertex_for_render[0].data(), collider_mesh_struct[i]->vertex_position[0].data(), 24 * collider_mesh_struct[i]->vertex_position.size());
	}
	if (*has_force) {
		if (obj_No < cloth->size()) {
			cloth->data()[obj_No].coe_neighbor_vertex_force = force_coe;
			cloth->data()[obj_No].neighbor_vertex = neighbor_vertex_index;
		}
		else {
			tetrahedron->data()[obj_No - cloth->size()].coe_neighbor_vertex_force = force_coe;
			tetrahedron->data()[obj_No - cloth->size()].neighbor_vertex = neighbor_vertex_index;
		}
	}
	updateRenderNormal();
	updateNormal();
	updateRenderVertexNormal();
}

double XPBD_IPC::calEdgeLength()
{
	double max_length = 0;
	double ave_length = 0;
	int edge_num = 0;
	double* edge_length;
	unsigned int edge_size;
	for (int i = 0; i < total_obj_num; ++i) {
		edge_size = mesh_struct[i]->edge_length.size();
		edge_length = mesh_struct[i]->edge_length.data();
		for (int j = 0; j < edge_size; ++j) {
			if (max_length < edge_length[j]) {
				max_length = edge_length[j];
			}
			ave_length += edge_length[j];
		}
		edge_num += edge_size;
	}
	return ave_length / (double)edge_num;
}

void XPBD_IPC::updateCollisionFreePosition()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(record_vertex_position[i][0].data(), vertex_position[i][0].data(), 24 * record_vertex_position[i].size());
		memcpy(record_gloabl_CCD_vertex_position[i][0].data(), vertex_position[i][0].data(), 24 * record_gloabl_CCD_vertex_position[i].size());
	}
}

void XPBD_IPC::initialCollisionConstriantNum()
{
	lambda_collision.resize(collision.collisionConstraintNumber(collision_constraint_index_start[0].data(), collision_constraint_index_start[1].data(), collision_constraint_index_start[2].data()));
}


void XPBD_IPC::XPBD_IPC_Position_Solve()
{
	updateCollisionFreePosition();
	thread->assignTask(this, SET_POS_PREDICT_);
	updateSn();
	firstNewtonCD();
	iteration_number = 0;
	if (perform_collision) {
		collision.collisionCulling();
	}
	outer_itr_num = 0;

	while (!convergeCondition(outer_itr_num)) {

		if (perform_collision) {
			collision.globalCollisionTime();
			thread->assignTask(this, COLLISION_FREE_POSITION_);
			updateCollisionFreePosition();
			collision.findClosePair();
			collision.saveCollisionPairVolume();
		}		
		nearly_not_move = false;

		while (!innerConvergeCondition(inner_iteration_number))
		{
			nearly_not_move = true;
			newtonCDTetWithCollision();
			inner_iteration_number++;
		}
		outer_itr_num++;
		//std::cout << inner_iteration_number << std::endl;
		iteration_number += inner_iteration_number;

	}


	collision.globalCollisionTime();
	thread->assignTask(this, COLLISION_FREE_POSITION_);

	thread->assignTask(this, XPBD_IPC_VELOCITY);
	updatePosition();
	updateRenderNormal();

	updateRenderVertexNormal();
}



void XPBD_IPC::XPBD_IPC_Block_Solve()
{
	std::cout << "====" << std::endl;

	record_energy.clear();
	updateCollisionFreePosition();
	thread->assignTask(this, SET_POS_PREDICT_);
	updateSn();
	iteration_number = 0;
	if (perform_collision) {
		collision.collisionCulling();
	}
	outer_itr_num = 0;
	computeCurrentEnergy();

	memset(lambda.data(), 0, 8 * lambda.size());

	while (!convergeCondition(outer_itr_num)) {
		if (perform_collision) {
			collision.globalCollisionTime();
			thread->assignTask(this, COLLISION_FREE_POSITION_);
			updateCollisionFreePosition();
			collision.findClosePair();
			//collision.saveCollisionPairVolume();
			//if (outer_itr_num == 0) {
			//	firstOnlyInertialCollision();
			//}			
		}
		inner_iteration_number = 0;
		nearly_not_move = false;
		previous_energy = energy;
		while (!innerConvergeCondition(inner_iteration_number))
		{
			previous_energy = energy;
			nearly_not_move = true;
			newtonCDTetBlock();
			computeCurrentEnergy();
			inner_iteration_number++;

			std::cout << "finish one itr " << inner_iteration_number<<" "<< energy << std::endl;

		}
		outer_itr_num++;
		//std::cout << inner_iteration_number << std::endl;
		iteration_number += inner_iteration_number;
	}

	if (perform_collision) {
		collision.globalCollisionTime();
		thread->assignTask(this, COLLISION_FREE_POSITION_);
	}
	thread->assignTask(this, XPBD_IPC_VELOCITY);
	updatePosition();
	updateRenderNormal();

	updateRenderVertexNormal();



	for (int i = 0; i < record_energy.size(); ++i) {
		std::cout << record_energy[i] << std::endl;
	}
}

void XPBD_IPC::XPBD_IPCSolve()
{
	updateCollisionFreePosition();
	thread->assignTask(this, SET_POS_PREDICT_);

	updateSn();
	firstNewtonCD();
	iteration_number = 0;
	if (perform_collision) {
		collision.collisionCulling();
	}
	outer_itr_num = 0;
	computeCurrentEnergy();

	memset(lambda.data(), 0, 8 * lambda.size());
	previous_energy = energy;
	//while (!convergeCondition(outer_itr_num)) {

	//	if (perform_collision) {
	//		collision.globalCollisionTime();
	//		thread->assignTask(this, COLLISION_FREE_POSITION_);
	//		updateCollisionFreePosition();
	//		collision.findClosePair();
	//		collision.saveCollisionPairVolume();
	//		firstOnlyInertialCollision();
	//	}
	//	inner_iteration_number = 0;
	//	nearly_not_move = false;	
		
		while (!innerConvergeCondition(inner_iteration_number))
		{
			previous_energy = energy;
			nearly_not_move = true;
			newtonCDTetWithCollision();
			inner_iteration_number++;
		}		
		outer_itr_num++;
		//std::cout << inner_iteration_number << std::endl;
		iteration_number += inner_iteration_number;
	//}

	if (perform_collision) {
		collision.globalCollisionTime();
		thread->assignTask(this, COLLISION_FREE_POSITION_);
	}
	thread->assignTask(this, XPBD_IPC_VELOCITY);
	updatePosition();
	updateRenderNormal();
	
	updateRenderVertexNormal();
}

bool XPBD_IPC::innerConvergeCondition(unsigned int iteration_num)
{
	if (iteration_num < min_inner_iteration) {//max_iteration_number
		return false;
	}

	if (abs(energy - previous_energy) / previous_energy < energy_converge_ratio) {
		return true;
	}


	if (iteration_num > max_iteration_number) {
		return true;
	}

	return false;

	//if (!nearly_not_move) {
	//	return false;
	//}

	//return true;
}



bool XPBD_IPC::convergeCondition(unsigned int iteration_num)
{
	if (iteration_num < min_outer_iteration) {//max_iteration_number
		return false;
	}

	//return true;

	if (iteration_num > outer_max_iteration_number - 1) {
		return true;
	}

	//if (abs(energy - previous_energy) / previous_energy < energy_converge_ratio) {
	//	return true;
	//}
	//return false;

	unsigned int num;
	std::array<double, 3>* current_pos;
	std::array<double, 3>* previous_pos;
	double* mass_;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		previous_pos = record_gloabl_CCD_vertex_position[i].data();
		current_pos = vertex_position[i];		
		num = mesh_struct[i]->vertex_position.size();
		mass_ = mesh_struct[i]->mass_inv.data();
		for (unsigned int j = 0; j < num; ++j) {
			if (mass_[j] != 0.0) {
				for (unsigned int k = 0; k < 3; ++k) {
					if (abs(previous_pos[j][k] - current_pos[j][k]) > max_move_standard) {
						return false;
					}
				}
			}		
		}
	}
	//std::cout << iteration_num << std::endl;
	//std::cout << abs(energy - previous_energy) / previous_energy << std::endl;
	return true;

}

void XPBD_IPC::updatePosition()
{

	unsigned int vertex_num;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_num = mesh_struct[i]->vertex_position.size();
		memcpy(initial_vertex_position[i][0].data(), vertex_position[i][0].data(), 24 * vertex_num);
	}

	for (unsigned int i = 0; i < collider->size(); ++i) {
		memcpy(collider->data()[i].mesh_struct.vertex_for_render[0].data(), collider->data()[i].mesh_struct.vertex_position[0].data(), 24 * collider->data()[i].mesh_struct.vertex_position.size());
	}
}

void XPBD_IPC::computeCurrentEnergy()
{
	energy = 0.0;
	energy += 0.5 * computeInertialEnergy();
	energy += computeCurrentARAPEnergy();
	record_energy.emplace_back(energy);
}

void XPBD_IPC::firstNewtonCD()
{
	previous_energy = energy;
	newtonCDTet();
}



void XPBD_IPC::newtonCDTetWithCCD()
{
	unsigned int size;
	std::array<int, 4>* indices;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	double* mass_inv;
	double stiffness;
	Matrix<double, 3, 4>* A;
	std::array<double, 3>* sn_;
	double* mass;
	double colliision_stiffness;
	int* vertex_index;
	std::vector<bool>* is_surface_vertex;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		mesh_struct_ = mesh_struct[i + cloth->size()];
		size = tetrahedron->data()[i].mesh_struct.vertex_position.size();
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i + cloth->size()];
		initial_vertex_pos = initial_vertex_position[i + cloth->size()];
		stiffness = tetrahedron->data()[i].ARAP_stiffness;
		A = tetrahedron->data()[i].mesh_struct.A.data();
		mass_inv = mesh_struct_->mass_inv.data();
		mass = mesh_struct_->mass.data();
		sn_ = sn[i + cloth->size()].data();
		vertex_index = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
		is_surface_vertex = &tetrahedron->data()[i].mesh_struct.vertex_on_surface;
		colliision_stiffness = tetrahedron->data()[i].collision_stiffness[0];

		for (unsigned int j = 0; j < size; ++j) {
			if (mass_inv[j] != 0.0) {
				solveNewtonCDTetWithCCD(vertex_pos, initial_vertex_pos, stiffness, sub_time_step, A,
					mesh_struct_->vertex_tet_index[j], indices, mass, volume, j, sn_,
					colliision_stiffness, i + cloth->size(), (*is_surface_vertex)[j], vertex_index[j]);
			}
		}
	}
}


void XPBD_IPC::firstOnlyInertialCollision()
{
	unsigned int size;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* last_step_vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	double* mass_inv;
	std::array<double, 3>* sn_;
	double* mass;
	double colliision_stiffness;
	int* vertex_index;
	std::vector<bool>* is_surface_vertex;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		mesh_struct_ = mesh_struct[i + cloth->size()];
		size = tetrahedron->data()[i].mesh_struct.vertex_position.size();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i + cloth->size()];
		last_step_vertex_pos = initial_vertex_position[i + cloth->size()];
		initial_vertex_pos = record_vertex_position[i + cloth->size()].data();
		mass_inv = mesh_struct_->mass_inv.data();
		mass = mesh_struct_->mass.data();
		sn_ = sn[i + cloth->size()].data();
		vertex_index = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
		is_surface_vertex = &tetrahedron->data()[i].mesh_struct.vertex_on_surface;
		colliision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
		//std::cout << "=== " << std::endl;
		for (unsigned int j = 0; j < size; ++j) {
			if (mass_inv[j] != 0.0) {
				solveInertialCollision(vertex_pos, initial_vertex_pos[j].data(),
					last_step_vertex_pos[j].data(), sub_time_step, mass, volume, j, sn_,
					colliision_stiffness, i + cloth->size(), (*is_surface_vertex)[j], vertex_index[j]);
			}
		}
	}
}


void XPBD_IPC::newtonCDTetWithCollision()
{
	unsigned int size;
	std::array<int, 4>* indices;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* last_step_vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	double* mass_inv;
	double stiffness;
	Matrix<double, 3, 4>* A;
	std::array<double, 3>* sn_;
	double* mass;
	double colliision_stiffness;
	int* vertex_index;
	std::vector<bool>* is_surface_vertex;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		mesh_struct_ = mesh_struct[i + cloth->size()];
		size = tetrahedron->data()[i].mesh_struct.vertex_position.size();
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i + cloth->size()];
		last_step_vertex_pos = initial_vertex_position[i + cloth->size()];
		initial_vertex_pos = record_vertex_position[i + cloth->size()].data();
		stiffness = tetrahedron->data()[i].ARAP_stiffness;
		A = tetrahedron->data()[i].mesh_struct.A.data();
		mass_inv = mesh_struct_->mass_inv.data();
		mass = mesh_struct_->mass.data();
		sn_ = sn[i + cloth->size()].data();
		vertex_index = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
		is_surface_vertex = &tetrahedron->data()[i].mesh_struct.vertex_on_surface;
		colliision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
		//std::cout << "=== " << std::endl;
		for (unsigned int j = 0; j < size; ++j) {
			if (mass_inv[j] != 0.0) {
					solveNewtonCDTetWithCollision(vertex_pos, initial_vertex_pos[j].data(),
					last_step_vertex_pos[j].data(),	stiffness, sub_time_step, A,
						mesh_struct_->vertex_tet_index[j], indices, mass, volume, j, sn_,
						colliision_stiffness, i + cloth->size(), (*is_surface_vertex)[j],vertex_index[j]);
			}
		}
	}
}


void XPBD_IPC::newtonCDTetBlock()
{
	unsigned int size;
	std::array<int, 4>* indices;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;
	double* mass_inv;
	double stiffness;
	Matrix<double, 3, 4>* A;
	std::array<double, 3>* sn_;
	double* mass;
	unsigned int* unfixed_vertex_num;
	std::vector<unsigned int>* tet_neightbor_tet;
	std::vector<unsigned int>* tet_neightbor_tet_common_vertex;
	std::array<int,4>* unfixed_vertex_index;

	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		mesh_struct_ = mesh_struct[i + cloth->size()];
		size = tetrahedron->data()[i].mesh_struct.indices.size();
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i + cloth->size()];
		stiffness =0.5* tetrahedron->data()[i].ARAP_stiffness;
		A = tetrahedron->data()[i].mesh_struct.A.data();
		mass_inv = mesh_struct_->mass_inv.data();
		mass = mesh_struct_->mass.data();
		sn_ = sn[i + cloth->size()].data();
		unfixed_vertex_num = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();
		tet_neightbor_tet = tetrahedron->data()[i].mesh_struct.tet_tet_index.data();
		tet_neightbor_tet_common_vertex= tetrahedron->data()[i].mesh_struct.tet_neighbor_tet_vertex_order.data();
		unfixed_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_indices.data();
		for (unsigned int j = 0; j < size; ++j) {
			solveNewtonCD_tetBlock(vertex_pos, stiffness, sub_time_step, mass, A, tet_neightbor_tet[j],
				indices, volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
				unfixed_vertex_index[j].data(), unfixed_vertex_num[j]);
		}
	}
}



void XPBD_IPC::newtonCDTet()
{
	unsigned int size;
	std::array<int, 4>* indices;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;
	double* mass_inv;
	double stiffness;
	Matrix<double, 3, 4>* A;
	std::array<double, 3>* sn_;
	double* mass;
	double* lambda_= lambda.data() + constraint_index_start[2];
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		mesh_struct_ = mesh_struct[i + cloth->size()];
		size = tetrahedron->data()[i].mesh_struct.vertex_position.size();
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i + cloth->size()];
		stiffness = tetrahedron->data()[i].ARAP_stiffness;
		A = tetrahedron->data()[i].mesh_struct.A.data();
		mass_inv = mesh_struct_->mass_inv.data();
		mass = mesh_struct_->mass.data();
		sn_ = sn[i + cloth->size()].data();
		for (unsigned int j = 0; j < size; ++j) {
			if (mass_inv[j] != 0.0) {
				solveNewtonCD_tet(vertex_pos, stiffness, sub_time_step, A,
					mesh_struct_->vertex_tet_index[j], indices, mass, volume, j, sn_,lambda_);
			}
		}
		lambda_ += tetrahedron->data()[i].mesh_struct.indices.size();
	}

}


bool XPBD_IPC::getFloorHessian(double& Hessian, double& grad, double* vertex_position, double floor_value,
	double* last_step_position, unsigned int dimension, double collision_stiffness, bool direction, double tolerance)
{
	if(direction){
		if (vertex_position[dimension] > floor_value+ tolerance) {
			return false;
		}
	}
	else {
		if (vertex_position[dimension] < floor_value- tolerance) {
			return false;
		}
	}
	double d_current = vertex_position[dimension] - floor_value;
	double d_ori = abs(last_step_position[dimension] - floor_value);// 
	if (d_ori < tolerance) {
		d_ori = tolerance;
	}
	double ln_;

	if (d_current < 0) {
		d_current = -d_current;
		ln_ = log(d_current / d_ori);
		grad =-collision_stiffness * (d_ori - d_current) * (2.0 * ln_ - d_ori / d_current + 1.0);
		Hessian =- collision_stiffness *(-2.0 * ln_ + (d_ori - d_current) * (d_ori + 3.0 * d_current) / (d_current * d_current));
	}
	else {
		ln_ = log(d_current / d_ori);
		grad = collision_stiffness * (d_ori - d_current) * (2.0 * ln_ - d_ori / d_current + 1.0);
		Hessian = collision_stiffness * (-2.0 * ln_ + (d_ori - d_current) * (d_ori + 3.0 * d_current) / (d_current * d_current));
	}
	return true;
}





void XPBD_IPC::getCollisionHessian(Matrix3d& Hessian, Vector3d& grad, std::array<double, 3>* vertex_position, 
	double* last_step_vertex_position,
	double collision_stiffness, unsigned int obj_No,
	unsigned int vertex_index,	unsigned int vertex_index_on_surface)
{
	Matrix3d a = Hessian;
	getVTCollisionHessain(Hessian, grad, vertex_position[vertex_index].data(), collision_stiffness,
		collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * vertex_index_on_surface, collision.vertex_triangle_pair_num_record[obj_No][vertex_index_on_surface],
		collision.VT_volume[obj_No].data() + collision.VT_start_index[obj_No][vertex_index_on_surface],obj_No,vertex_index);
	int* triangle_;
	std::vector<unsigned int>* triangle = &mesh_struct[obj_No]->vertices[vertex_index].face;
	for (unsigned int i = 0; i < triangle->size(); ++i) {
		triangle_ = triangle_indices[obj_No][triangle->data()[i]].data();
		getTVCollisionHessain(Hessian, grad, vertex_position[triangle_[0]].data(), vertex_position[triangle_[1]].data(),
			vertex_position[triangle_[2]].data(),
			findVertexNo(vertex_index, triangle_, 3),
			collision_stiffness, collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * triangle->data()[i],
			collision.triangle_vertex_pair_num_record[obj_No][triangle->data()[i]],
			collision.TV_volume[obj_No].data() + collision.TV_start_index[obj_No][triangle->data()[i]]);
	}
	std::vector<unsigned int>* edge = &mesh_struct[obj_No]->vertices[vertex_index].edge;
	unsigned int* edge_;
	for (unsigned int i = 0; i < edge->size(); ++i) {
		edge_ = edge_vertices[obj_No] + (edge->data()[i] << 1);
		getEECollisionHessian(Hessian, grad, vertex_position[edge_[0]].data(), vertex_position[edge_[1]].data(),
			collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * edge->data()[i],
			collision.edge_edge_pair_num_record[obj_No][edge->data()[i]], 
			collision.EE_volume[obj_No].data() + collision.EE_start_index[obj_No][edge->data()[i]], 
			collision_stiffness, obj_No, edge->data()[i],
			findVertexNo(vertex_index, edge_, 2));
	}
	if (!collider->empty()) {
		getVT_ColiderCollisionHessain(Hessian, grad, vertex_position[vertex_index].data(), collision_stiffness,
			collision.vertex_obj_triangle_collider_pair_by_vertex[obj_No] + collision.close_vt_collider_pair_num * vertex_index_on_surface,
			collision.vertex_obj_triangle_collider_num_record[obj_No][vertex_index_on_surface],
			collision.VT_collider_volume[obj_No].data() + collision.VT_collider_start_index[obj_No][vertex_index_on_surface]);
	}

	//floor
	if (floor->exist) {
		double Hessian_, grad_;
		if (getFloorHessian(Hessian_, grad_, vertex_position[vertex_index].data(), floor->value,
			last_step_vertex_position, floor->dimension, collision_stiffness, floor->normal_direction, collision.d_hat)) {
			grad.data()[floor->dimension] -= grad_;
			Hessian.data()[floor->dimension << 2] += Hessian_;
		}	
	}

}


void XPBD_IPC::solveNewtonCDTetWithCCD(std::array<double, 3>* vertex_position, std::array<double, 3>* initial_vertex_position,
	double ARAP_stiffness, double dt,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* tet_vertex_indices, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn, double collision_stiffness, unsigned int obj_No,
	bool vertex_on_surface, unsigned int vertex_index_on_surface)
{
	if (vertex_on_surface) {
		collision.collisionTimeSingleVertex(obj_No, vertex_index, vertex_index_on_surface,  initial_vertex_position,
			vertex_position);

	}
	
}

void XPBD_IPC::solveInertialCollision(std::array<double, 3>* vertex_position,
	double* record_vertex_position_,
	double* last_step_vertex_position,
	double dt, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn, double collision_stiffness, unsigned int obj_No,
	bool vertex_on_surface, unsigned int vertex_index_on_surface)
{
	Matrix3d Hessian;
	Vector3d grad;
	Hessian.setZero();
	grad.setZero();
	//if (vertex_on_surface) {
	//	getCollisionHessian(Hessian, grad, vertex_position, last_step_vertex_position, collision_stiffness, obj_No, vertex_index, vertex_index_on_surface);
	//}
	double mass_dt_2 = mass[vertex_index] / (dt * dt);

	Hessian.data()[0] += mass_dt_2;
	Hessian.data()[4] += mass_dt_2;
	Hessian.data()[8] += mass_dt_2;
	grad.data()[0] -= mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	grad.data()[1] -= mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	grad.data()[2] -= mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <Matrix3d> linear(Hessian);
	Vector3d result = linear.solve(grad);

	SUM_(vertex_position[vertex_index], result);
	if (vertex_on_surface) {
		collision.collisionFreeOneVertex(obj_No, vertex_index, vertex_index_on_surface,
			record_vertex_position_, vertex_position[vertex_index].data(),
			record_vertex_position[obj_No].data(), vertex_position, this->vertex_position.data());
	}
	memcpy(record_vertex_position_, vertex_position[vertex_index].data(), 24);
}




void XPBD_IPC::solveNewtonCDTetWithCollision(std::array<double, 3>* vertex_position, 
	double* record_vertex_position_,
	double* last_step_vertex_position,
	double ARAP_stiffness, double dt,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* tet_vertex_indices, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn, double collision_stiffness, unsigned int obj_No,
	bool vertex_on_surface, unsigned int vertex_index_on_surface)
{
	Matrix3d Hessian;
	Vector3d grad;
	Hessian.setZero();
	grad.setZero();

	getARAPHessian(Hessian, grad, vertex_position, ARAP_stiffness, A, tet_indices, tet_vertex_indices, volume, vertex_index, obj_No);

	if (perform_collision) {
		if (vertex_on_surface) {
			getCollisionHessian(Hessian, grad, vertex_position, last_step_vertex_position, collision_stiffness, obj_No, vertex_index, vertex_index_on_surface);
		}
	}


	double mass_dt_2 = mass[vertex_index] / (dt * dt);

	Hessian.data()[0] += mass_dt_2;
	Hessian.data()[4] += mass_dt_2;
	Hessian.data()[8] += mass_dt_2;
	grad.data()[0] -= mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	grad.data()[1] -= mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	grad.data()[2] -= mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <Matrix3d> linear(Hessian);
	Vector3d result = linear.solve(grad);

	SUM_(vertex_position[vertex_index], result);
	
//	if (obj_No == 0) {
//		if (inner_iteration_number == 1) {
//			if (abs(result[1]) < 1e-4) {
//				std::cout << "index " << vertex_index << std::endl;
//			}
//		}
//		
///*		if (vertex_index==11 || vertex_index == 8) {
//			std::cout <<vertex_index<<" "<< sn[vertex_index][1] << " " << vertex_position[vertex_index][1] << std::endl;
//			std::cout << result[1] << " " << vertex_index<<" mdt2 "<<mass_dt_2 << std::endl;
//		}	*/	
//	}

	if (nearly_not_move) {
		for (unsigned int i = 0; i < 3; ++i) {
			if (abs(result[i]) > max_move_standard_inner_itr) {
				nearly_not_move = false;
				break;
			}
		}
	}	

	if (perform_collision) {
		if (vertex_on_surface) {
			collision.collisionFreeOneVertex(obj_No, vertex_index, vertex_index_on_surface,
				record_vertex_position_, vertex_position[vertex_index].data(),
				record_vertex_position[obj_No].data(), vertex_position, this->vertex_position.data());
		}
		memcpy(record_vertex_position_, vertex_position[vertex_index].data(), 24);
	}



}





void XPBD_IPC::getEECollisionHessian(Matrix3d& Hessian, Vector3d& grad, double* pos0, double* pos1, unsigned int* EE, unsigned int num,
	double* ori_volume, double stiffness, unsigned int obj_index, unsigned int edge_index, unsigned int vertex_no)
{
	Matrix3d Hessian_single;
	Vector3d grad_single;
	unsigned int* edge_vertex;
	double volume;
	for (unsigned int i = 0; i < num; i += 2) {
		volume = (ori_volume[i >> 1] > collision.volume_boundary ? ori_volume[i >> 1] : collision.volume_boundary);
		edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
		if (obj_index < EE[i] || (obj_index == EE[i] && edge_index < EE[i + 1])) {
			if (second_order_constraint.getCollisionPairHessian(pos0, pos1, vertex_position[EE[i]][edge_vertex[0]].data(),
				vertex_position[EE[i]][edge_vertex[1]].data(),
				volume, Hessian_single, grad_single, vertex_no)) {
				Hessian_single *= stiffness;
				grad_single *= stiffness;
				grad -= grad_single;
				Hessian += Hessian_single;
			}
		}
		else {
			if (second_order_constraint.getCollisionPairHessian(vertex_position[EE[i]][edge_vertex[0]].data(),
				vertex_position[EE[i]][edge_vertex[1]].data(), pos0, pos1,
				volume, Hessian_single, grad_single, vertex_no + 2)) {
				Hessian_single *= stiffness;
				grad_single *= stiffness;
				grad -= grad_single;
				Hessian += Hessian_single;
			}
		}
	}
}

void XPBD_IPC::getTVCollisionHessain(Matrix3d& Hessian, Vector3d& grad,
	double* pos_0, double* pos_1, double* pos_2,
	unsigned int vertex_no, double stiffness, unsigned int* TV, unsigned int num, double* ori_volume)
{
	Matrix3d Hessian_single;
	Vector3d grad_single;
	double volume;
	for (unsigned int i = 0; i < num; i += 2) {
		volume = (ori_volume[i >> 1] > collision.volume_boundary ? ori_volume[i >> 1] : collision.volume_boundary);
		if (second_order_constraint.getCollisionPairHessian(pos_0, pos_1,  pos_2,
			vertex_position[TV[i]][TV[i+1]].data(),
			volume, Hessian_single, grad_single, vertex_no)) {
			Hessian_single *= stiffness;
			grad_single *= stiffness;
			grad -= grad_single;
			Hessian += Hessian_single;
		}
	}
}



void XPBD_IPC::getVT_ColiderCollisionHessain(Matrix3d& Hessian, Vector3d& grad, double* vertex_position_, double stiffness,
	unsigned int* VT, unsigned int num, double* ori_volume)
{
	Matrix3d Hessian_single;
	Vector3d grad_single;
	int* triangle_vertex;
	double volume;
	for (unsigned int i = 0; i < num; i += 2) {
		triangle_vertex = triangle_indices_collider[VT[i]][VT[i + 1]].data();
		volume = (ori_volume[i >> 1] > collision.volume_boundary ? ori_volume[i >> 1] : collision.volume_boundary);
		if (second_order_constraint.getCollisionPairHessian(vertex_position_collider[VT[i]][triangle_vertex[0]].data(),
			vertex_position_collider[VT[i]][triangle_vertex[1]].data(), vertex_position_collider[VT[i]][triangle_vertex[2]].data(), vertex_position_,
			volume, Hessian_single, grad_single, 3)) {
			Hessian_single *= stiffness;
			grad_single *= stiffness;
			grad -= grad_single;
			Hessian += Hessian_single;
		}
	}
}



void XPBD_IPC::getVTCollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
	unsigned int* VT, unsigned int num, unsigned int obj_No,unsigned int vertex_order_in_matrix,
	unsigned int vertex_index_on_surface, unsigned int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num )
{
	collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * vertex_index_on_surface;
	collision.vertex_triangle_pair_num_record[obj_No][vertex_index_on_surface];
	int* triangle_vertex;
	std::vector<int>triangle_vertex_order_in_tet, triangle_vertex_order_in_pair;
	triangle_vertex_order_in_tet.reserve(4);
	triangle_vertex_order_in_pair.reserve(4);
	triangle_vertex_order_in_tet.emplace_back(vertex_order_in_matrix);
	triangle_vertex_order_in_pair.emplace_back(0);
	for (unsigned int i = 0; i < num; i += 2) {
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		triangle_vertex_order_in_tet.resize(1);
		triangle_vertex_order_in_pair.resize(1);
		checkPairIndexInTriangle(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, &triangle_vertex_order_in_tet,
			&triangle_vertex_order_in_pair);
		second_order_constraint.getCollisionHessianVertex(vertex_position_, vertex_position[VT[i]][triangle_vertex[0]].data(),
			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(),
			Hessian, grad, triangle_vertex_order_in_tet, triangle_vertex_order_in_pair);
	}

}






void XPBD_IPC::checkPairIndexInTriangle(int unfixed_tet_vertex_num, unsigned int* tet_unfixed_vertex_indices, int* triangle_indices, std::vector<int>* triangle_vertex_order_in_system,
	std::vector<int>* triangle_vertex_order_in_pair)
{
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < unfixed_tet_vertex_num; ++j) {
			if (tet_unfixed_vertex_indices[j] == triangle_indices[i]) {
				triangle_vertex_order_in_system->emplace_back(j);
				triangle_vertex_order_in_pair->emplace_back(i);
				break;
			}
		}
	}

}


void XPBD_IPC::getVTCollisionHessain(Matrix3d& Hessian, Vector3d& grad, double* vertex_position_, double stiffness, 
	unsigned int* VT, unsigned int num, double* ori_volume, unsigned int obj_No, unsigned int vertex_index)
{
	Matrix3d Hessian_single;
	Vector3d grad_single;
	int* triangle_vertex;
	double volume;
	for (unsigned int i = 0; i < num; i += 2) {
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		volume =( ori_volume[i >> 1]>collision.volume_boundary?ori_volume[i >> 1] :collision.volume_boundary);
		if (second_order_constraint.getCollisionPairHessian( vertex_position[VT[i]][triangle_vertex[0]].data(),
			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(), vertex_position_,
			volume, Hessian_single, grad_single, 3)) {
			Hessian_single *= stiffness;
			grad_single *= stiffness;
			grad -= grad_single;
			Hessian += Hessian_single;
		}
	}


}


void XPBD_IPC::getARAPHessian(Matrix3d& Hessian, Vector3d& grad, std::array<double, 3>* vertex_position, double stiffness,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, 
	double* volume, unsigned int vertex_index, unsigned int obj_No)
{
	unsigned int tet_index, vertex_no;
	Matrix3d Hessian_single;
	Vector3d grad_single;
	double C;
	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		tet_index = tet_indices[i];
		vertex_no = findVertexNo(vertex_index, indices[tet_index].data(),4);

		if (second_order_constraint.getARAPGradHessianNewton(vertex_position[indices[tet_index][0]].data(), vertex_position[indices[tet_index][1]].data(),
			vertex_position[indices[tet_index][2]].data(), vertex_position[indices[tet_index][3]].data(),
			A[tet_index], Hessian_single, grad_single, C, vertex_no)) {
			Hessian += (0.5 * stiffness * volume[tet_index]) * Hessian_single;
			grad -= (0.5 * stiffness * volume[tet_index]) * grad_single;
		}
		//if (obj_No == 0 && vertex_index == 11) {
		//	std::cout << C << " " << grad_single.transpose() << std::endl;
		//	std::cout << Hessian_single << std::endl;
		//}

	}
}


void XPBD_IPC::solveNewtonCD_tetBlock(std::array<double, 3>* vertex_position, double stiffness, double dt,
	double* mass,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices, std::array<int, 4>* indices,
	double* volume, unsigned int tet_index, std::array<double, 3>* sn, unsigned int* common_vertex_in_order,
	int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num)
{
	if (unfixed_vertex_num == 0) {
		return;
	}

	MatrixXd Hessian;
	VectorXd grad;
	grad.resize(3 * unfixed_vertex_num);
	Hessian.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);

	second_order_constraint.solveCD_ARAP_block(Hessian, grad, vertex_position, stiffness, A, neighbor_tet_indices,
		indices, volume, tet_index, common_vertex_in_order, tet_vertex_index,
		unfixed_tet_vertex_index, unfixed_vertex_num);

	double mass_dt_2;
	int vertex_index;
	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_vertex_index[unfixed_tet_vertex_index[i]];
		mass_dt_2 = mass[vertex_index] / (dt * dt);
		Hessian(3 * i, 3 * i) += mass_dt_2;
		Hessian(3 * i + 1, 3 * i + 1) += mass_dt_2;
		Hessian(3 * i + 2, 3 * i + 2) += mass_dt_2;
		grad.data()[3 * i] += mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
		grad.data()[3 * i + 1] += mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
		grad.data()[3 * i + 2] += mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);
	}


	LDLT <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(grad);

	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_vertex_index[unfixed_tet_vertex_index[i]];
		SUB_(vertex_position[vertex_index], (result.data() + 3 * i));

	}

}



void XPBD_IPC::solveNewtonCD_tet(std::array<double, 3>* vertex_position, double stiffness, double dt,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn, double* lambda)
{
	Matrix3d Hessian;
	Vector3d grad;
	Hessian.setZero();
	grad.setZero();

	getARAPHessian(Hessian, grad, vertex_position, stiffness, A, tet_indices, indices, volume, vertex_index,0);


	double mass_dt_2 = mass[vertex_index] / (dt * dt);

	Hessian.data()[0] += mass_dt_2;
	Hessian(1, 1) += mass_dt_2;
	Hessian(2, 2) += mass_dt_2;
	grad.data()[0] -= mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	grad.data()[1] -= mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	grad.data()[2] -= mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <Matrix3d> linear(Hessian);
	Vector3d result = linear.solve(grad);
	SUM_(vertex_position[vertex_index], result);

}


//XPBD_IPC_VELOCITY
void XPBD_IPC::computeVelocity(int thread_No)
{
	double damp_coe = velocity_damp;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	unsigned int vertex_end = 0;
	std::array<double, 3>* velocity_;
	double delta_t = sub_time_step / damp_coe;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		initial_vertex_pos = initial_vertex_position[i];
		vertex_end = vertex_index_begin_per_thread[i][thread_No + 1];
		velocity_ = velocity[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < vertex_end; ++j) {
			for (unsigned int k = 0; k < 3; ++k) {
				velocity_[j][k] = (vertex_pos[j][k] - initial_vertex_pos[j][k]) / delta_t;
			}
		}
	}
}


void XPBD_IPC::updateSn()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(sn[i][0].data(), vertex_position[i][0].data(), sn[i].size() * 24);
	}
}


double XPBD_IPC::computeInertialEnergy()
{
	double energy = 0.0;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* sn_;
	unsigned int vertex_end;
	double* mass;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		sn_ = sn[i].data();
		vertex_end = vertex_index_begin_per_thread[i][total_thread_num];
		mass = mesh_struct[i]->mass.data();
		for (unsigned int j = 0; j < vertex_end; ++j) {
			energy += mass[j] * (EDGE_LENGTH(vertex_pos[j], sn_[j]));
		}
	}
	return energy / (sub_time_step * sub_time_step);
}
//SET_POS_PREDICT_
void XPBD_IPC::setPosPredict(int thread_No)
{
	std::array<double, 3>* vertex_pos;
	unsigned int vertex_end = 0;

	double delta_t = time_step;
	double delta_t_2 = delta_t * delta_t;
	double* mass_inv;
	std::array<double, 3>* f_ext_;
	std::array<double, 3>* velocity_;
	double gravity__[3];
	memcpy(gravity__, gravity, 24);
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		vertex_end = vertex_index_begin_per_thread[i][thread_No + 1];
		mass_inv = mesh_struct[i]->mass_inv.data();
		f_ext_ = f_ext[i].data();
		velocity_ = velocity[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < vertex_end; ++j) {
			if (mass_inv[j] != 0) {
				for (unsigned int k = 0; k < 3; ++k) {
					vertex_pos[j][k] += delta_t * velocity_[j][k] + delta_t_2 * (mass_inv[j] * f_ext_[j][k] + gravity__[k]);
				}
			}
		}
	}
}


void XPBD_IPC::resetExternalForce()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(f_ext[i][0].data(), 0, 24 * f_ext[i].size());
	}
}

void XPBD_IPC::initialDHatTolerance(double ave_edge_length)
{
	if (perform_collision) {
		collision.initialDHatTolerance(ave_edge_length);
	}
}

void XPBD_IPC::updateTetrahedronAnchorVertices()
{
	double* mass_inv;
	int* anchor_vertex;
	unsigned int anchor_vertex_size;
	for (unsigned int i = cloth->size(); i < cloth->size() + tetrahedron->size(); ++i) {
		mass_inv = mesh_struct[i]->mass_inv.data();
		anchor_vertex_size = mesh_struct[i]->anchor_vertex.size();
		anchor_vertex = mesh_struct[i]->anchor_vertex.data();
		mesh_struct[i]->resetMassInv();

	}
}



void XPBD_IPC::addExternalForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int obj_No)
{
	if (!coe.empty()) {
		for (unsigned int i = 0; i < coe.size(); ++i) {
			for (unsigned int j = 0; j < 3; ++j) {
				f_ext[obj_No][neighbor_vertex[i]][j] += coe[i] * neighbor_vertex_force_direction[j];
			}
		}
	}
}


void XPBD_IPC::updateRenderVertexNormal()
{
	for (unsigned int i = 0; i < collider->size(); ++i) {
		thread->assignTask(&(*collider)[i].mesh_struct, VERTEX_NORMAL_FROM_RENDER);
	}
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, VERTEX_NORMAL_FROM_RENDER);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, VERTEX_NORMAL_FROM_RENDER);
	}
}


void XPBD_IPC::updateRenderNormal()
{
	for (unsigned int i = 0; i < collider->size(); ++i) {
		thread->assignTask(&(*collider)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
}

void XPBD_IPC::updateNormal()
{
	for (unsigned int i = 0; i < collider->size(); ++i) {
		thread->assignTask(&(*collider)[i].mesh_struct, FACE_NORMAL);
	}
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL);
	}
}

double XPBD_IPC::computeCurrentARAPEnergy()
{
	double energy = 0.0;
	unsigned int size;
	std::array<int, 4>* indices;
	std::array<double, 3>* vertex_pos;
	Matrix<double, 3, 4>* A;
	double* volume;
	double stiffness;
	double* mass_inv_;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		size = tetrahedron->data()[i].mesh_struct.indices.size();
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		vertex_pos = vertex_position[i + cloth->size()];
		A = tetrahedron->data()[i].mesh_struct.A.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		stiffness = tetrahedron->data()[i].ARAP_stiffness;
		mass_inv_ = tetrahedron->data()[i].mesh_struct.mass_inv.data();
		for (unsigned int j = 0; j < size; ++j) {
			if (mass_inv_[indices[j][0]] != 0.0 || mass_inv_[indices[j][1]] != 0.0 || mass_inv_[indices[j][2]] != 0.0 || mass_inv_[indices[j][3]] != 0.0) {
				energy += compute_energy.computeARAPEnergy(vertex_pos[indices[j][0]].data(), vertex_pos[indices[j][1]].data(),
					vertex_pos[indices[j][2]].data(), vertex_pos[indices[j][3]].data(), A[j], volume[j], stiffness);
			}
		}

	}
	return energy;
}

//COLLISION_FREE_POSITION_
void XPBD_IPC::computeCollisionFreePosition(int thread_No)
{
	unsigned int index_end;

	double collision_time = collision.collision_time;
	std::array<double, 3>* q_pre;
	std::array<double, 3>* q_end;

	for (unsigned int i = 0; i <total_obj_num; ++i) {
		index_end = vertex_index_begin_per_thread[i][thread_No + 1];
		q_end = vertex_position[i];
		q_pre = record_gloabl_CCD_vertex_position[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			q_end[j][0] = q_pre[j][0] + collision_time * (q_end[j][0] - q_pre[j][0]);
			q_end[j][1] = q_pre[j][1] + collision_time * (q_end[j][1] - q_pre[j][1]);
			q_end[j][2] = q_pre[j][2] + collision_time * (q_end[j][2] - q_pre[j][2]);
		}
	}
}

