#include"XPBD_IPC.h"
#include<algorithm>
#include"XPBD/FEM_relate.h"

XPBD_IPC::XPBD_IPC()
{
	gravity_ = 9.8;
	sub_step_num = 1;
	iteration_number = 300;

	damping_coe = 0.0;

	perform_collision = true;

	XPBD_constraint.epsilon_for_bending = 1e-10;

	velocity_damp = 0.995;
	energy_converge_ratio = 1e-3;

	min_inner_iteration = 3;
	min_outer_iteration = 2;

	max_move_standard = 1e-3;
	max_move_standard_inner_itr = 1e-4;
	max_iteration_number = 50;
	outer_max_iteration_number = 40;
	energy_converge_standard = 1e-11;

	second_order_constraint.solve_exact_ARAP_hessian = false;
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
	has_collider = !collider->empty();
	this->floor = floor;

	reorganzieDataOfObjects();
	initialVariable();
	initialClothBending();
	setConstraintIndex();
	//energy_per_thread.resize(thread->thread_num,0.0);
	if (perform_collision) {
		//collision.energy = energy_per_thread.data();
		collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio, XPBD_IPC_,true);
		collision.setCollisionFreeVertex(&record_collision_free_vertex_position_address, &record_vertex_position);
		//collision.setParameter(&lambda_collision,lambda.data()+ constraint_index_start[3], collision_constraint_index_start.data(), damping_coe, sub_time_step);
	}

	initalARAPHessianStorages();
	setColorNum();

}


void XPBD_IPC::setColorNum()
{
	int color_num = 1;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		if (color_num < tetrahedron->data()[i].mesh_struct.tet_color_group[0].size()) {
			color_num = tetrahedron->data()[i].mesh_struct.tet_color_group[0].size();
		}
	}
	max_tet_color_num = color_num;
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
	vertex_position_render.resize(total_obj_num);
	mesh_struct.resize(total_obj_num);
	vertex_index_begin_per_thread.resize(total_obj_num);
	record_vertex_position.resize(total_obj_num);
	record_vertex_position_address.resize(total_obj_num);
	record_collision_free_vertex_position_address.resize(total_obj_num);
	record_collision_free_vertex_position.resize(total_obj_num);

	//record_outer_vertex_position.resize(total_obj_num);
	//unfixed_vertex.resize(total_obj_num);
	triangle_indices.resize(total_obj_num);
	tet_indices.resize(total_obj_num);
	tet_A.resize(total_obj_num);
	tet_volume.resize(total_obj_num);
	edge_vertices.resize(total_obj_num);
	is_vertex_fixed.resize(total_obj_num);

	vertex_index_surface.resize(total_obj_num);
	mass.resize(total_obj_num);


	triangle_around_triangle.resize(total_obj_num);
	edge_around_triangle.resize(total_obj_num);
	vertices.resize(total_obj_num);
	tet_around_vertex.resize(total_obj_num);
	tet_around_triangle.resize(total_obj_num);
	vertex_surface_to_global.resize(total_obj_num);

	triangle_around_edge.resize(total_obj_num);
	edge_around_edge.resize(total_obj_num);
	tet_around_edge.resize(total_obj_num);


	rest_edge_length.resize(total_obj_num);



	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position.data();
		vertex_position_render[i] = cloth->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i] = &cloth->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		record_vertex_position_address[i]= record_vertex_position[i].data();
		record_collision_free_vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position;
		record_collision_free_vertex_position_address[i] = record_collision_free_vertex_position[i].data();
		triangle_indices[i] = cloth->data()[i].mesh_struct.triangle_indices.data();
		//record_outer_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		//unfixed_vertex[i] = &cloth->data()[i].mesh_struct.unfixed_point_index;
		edge_vertices[i] = cloth->data()[i].mesh_struct.edge_vertices.data();

		is_vertex_fixed[i] = &cloth->data()[i].mesh_struct.is_vertex_fixed;
		vertex_index_surface[i] = cloth->data()[i].mesh_struct.vertex_surface_index.data();
		mass[i] = cloth->data()[i].mesh_struct.mass.data();

		triangle_around_triangle[i] = cloth->data()[i].mesh_struct.face_around_face.data();
		edge_around_triangle[i] = cloth->data()[i].mesh_struct.edge_around_face.data();
		vertices[i] = cloth->data()[i].mesh_struct.vertices.data();

		triangle_around_edge[i] = cloth->data()[i].mesh_struct.face_around_edge.data();
		edge_around_edge[i] = cloth->data()[i].mesh_struct.edge_around_edge.data();

		rest_edge_length[i] = cloth->data()[i].mesh_struct.edge_length.data();

	}

	unfix_tet_index.resize(tetrahedron->size());
	unfixed_tet_vertex_num.resize(tetrahedron->size());

	tet_color_groups.resize(tetrahedron->size());

	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		vertex_position_render[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		record_vertex_position_address[i + cloth->size()] = record_vertex_position[i + cloth->size()].data();
		record_collision_free_vertex_position[i + cloth->size()]=tetrahedron->data()[i].mesh_struct.vertex_position;
		record_collision_free_vertex_position_address[i + cloth->size()] = record_collision_free_vertex_position[i + cloth->size()].data();
		triangle_indices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.triangle_indices.data();
		tet_indices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.indices.data();
		tet_A[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.A.data();
		tet_volume[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.volume.data();
		//record_outer_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		//unfixed_vertex[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.unfixed_point_index;
		edge_vertices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_vertices.data();

		is_vertex_fixed[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.is_vertex_fixed;
		vertex_index_surface[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
		mass[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.mass.data();

		triangle_around_triangle[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.face_around_face.data();
		edge_around_triangle[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_around_face.data();
		vertices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertices.data();
		tet_around_vertex[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_tet_index.data();
		tet_around_triangle[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tet_around_face.data();
		vertex_surface_to_global[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_on_sureface.data();

		triangle_around_edge[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.face_around_edge.data();
		edge_around_edge[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_around_edge.data();
		tet_around_edge[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tet_around_edge.data();
		rest_edge_length[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_length.data();

		
		unfix_tet_index[i] = tetrahedron->data()[i].mesh_struct.unfixied_indices.data();
		unfixed_tet_vertex_num[i] = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();

		tet_color_groups[i] =&tetrahedron->data()[i].mesh_struct.tet_color_group;
	}
	if (!collider->empty()) {
		triangle_indices_collider.resize(collider->size());
		collider_mesh_struct.resize(collider->size());
		collider_edge_vertices.resize(collider->size());
		vertex_position_collider.resize(collider->size());
		rest_edge_length_collider.resize(collider->size());
		


		for (unsigned int i = 0; i < collider->size(); ++i) {
			collider_mesh_struct[i] = &collider->data()[i].mesh_struct;
			triangle_indices_collider[i] = collider->data()[i].mesh_struct.triangle_indices.data();
			collider_edge_vertices[i] = collider->data()[i].mesh_struct.edge_vertices.data();
			vertex_position_collider[i] = collider->data()[i].mesh_struct.vertex_position.data();
			rest_edge_length_collider[i] = collider->data()[i].mesh_struct.edge_length.data();
		}
	}
	


	int total_vertex_num = 0;
	vertex_num_on_surface_prefix_sum.resize(total_obj_num+1);
	vertex_num_on_surface_prefix_sum[0] = 0;
	for (int i = 1; i <= total_obj_num; ++i) {
		if (i <= cloth->size()) {
			vertex_num_on_surface_prefix_sum[i] = vertex_num_on_surface_prefix_sum[i - 1] + vertex_index_begin_per_thread[i - 1][total_thread_num];
		}
		else {
			vertex_num_on_surface_prefix_sum[i] = vertex_num_on_surface_prefix_sum[i - 1] + tetrahedron->data()[i].mesh_struct.vertex_index_on_sureface.size();
		}
	}




}


void XPBD_IPC::initalARAPHessianStorages()
{
	prefix_sum_of_every_tet_index.resize(tetrahedron->size()+1);
	prefix_sum_of_every_tet_index[0] = 0;
	for (int i = 0; i < tetrahedron->size(); ++i) {
		prefix_sum_of_every_tet_index[i + 1] = prefix_sum_of_every_tet_index[i] + tetrahedron->data()[i].mesh_struct.indices.size();
	}
	store_tet_arap_hessian.resize(16 * prefix_sum_of_every_tet_index[tetrahedron->size()],0.0);
	store_tet_arap_grad.resize(12 * prefix_sum_of_every_tet_index[tetrahedron->size()],0.0);

	//auto t0 = std::chrono::system_clock::now();

	//for (int i = 0; i < total_thread_num; ++i) {
	//	tetHessian(i);
	//}

	//auto t1 = std::chrono::system_clock::now();

	thread->assignTask(this, UPDATE_TET_HESSIAN);

	//auto t2 = std::chrono::system_clock::now();

	//double time_0 = double(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den * 1000.0;
	//double time_1 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den * 1000.0;

	//std::cout << "time " << time_0 << " " << time_1 << std::endl;

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

void XPBD_IPC::recordInitialPosition()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(vertex_position_render[i][0].data(), vertex_position[i][0].data(), 24 * record_vertex_position[i].size());
	}
}

void XPBD_IPC::updateCollisionFreePosition()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(record_vertex_position[i][0].data(), vertex_position[i][0].data(), 24 * record_vertex_position[i].size());
		memcpy(record_collision_free_vertex_position[i][0].data(), vertex_position[i][0].data(), 24 * record_vertex_position[i].size());
	}
}


void XPBD_IPC::tempSavePos()
{
	temp_save_pos.resize(total_obj_num);
	for (int i = 0; i < total_obj_num; ++i) {
		temp_save_pos[i].resize(record_vertex_position[i].size());
		memcpy(temp_save_pos[i][0].data(), vertex_position[i][0].data(), 24 * record_vertex_position[i].size());
	}
}

void XPBD_IPC::tempRestorePos()
{
	for (int i = 0; i < total_obj_num; ++i) {
		memcpy(vertex_position[i][0].data(), temp_save_pos[i][0].data(),  24 * record_vertex_position[i].size());
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



void XPBD_IPC::XPBD_IPC_Block_Solve_Multithread()
{
	//recordInitialPosition();
	updateCollisionFreePosition();
	thread->assignTask(this, SET_POS_PREDICT_);
	updateSn();
	iteration_number = 0;
	if (perform_collision) {
		collision.collisionCulling();
	}
	outer_itr_num = 0;
	displacement_satisfied = false;
	//vertex_trace.clear();
	//vertex_trace.push_back(vertex_position[0][6]);
	while (!convergeCondition(outer_itr_num)) {
		if (perform_collision) {
			//collision.collisionCulling();
			collision.globalCollisionTime();
			thread->assignTask(this, COLLISION_FREE_POSITION_);
			updateCollisionFreePosition();
			collision.findClosePair();			
		}
		inner_iteration_number = 0;
		nearly_not_move = false;
		previous_energy = energy;
		while (!innerConvergeCondition(inner_iteration_number))
		{
			previous_energy = energy;
			nearly_not_move = true;
			solveNewtonCD_tetBlock();
			//newtonCDBlock();

			if (inner_iteration_number >= min_inner_iteration - 1) {
				computeCurrentEnergy();
			}
			//std::cout << "finish one itr " << inner_iteration_number << " " << energy << std::endl;
			inner_iteration_number++;
			//vertex_trace.push_back(vertex_position[0][6]);
		}

		//std::cout << "outer " << outer_itr_num << std::endl;

		outer_itr_num++;
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

	std::cout <<"time stamp "<< * time_stamp << std::endl;
}



void XPBD_IPC::XPBD_IPC_Block_Solve()
{
	//std::cout << "====" << std::endl;

	record_energy.clear();
	recordInitialPosition();
	updateCollisionFreePosition();
	thread->assignTask(this, SET_POS_PREDICT_);
	updateSn();
	iteration_number = 0;
	if (perform_collision) {
		collision.collisionCulling();
	}
	outer_itr_num = 0;
	displacement_satisfied = false;
	while (!convergeCondition(outer_itr_num)) {
		if (perform_collision) {
			collision.collisionCulling();
			collision.globalCollisionTime();
			thread->assignTask(this, COLLISION_FREE_POSITION_);
			updateCollisionFreePosition();
			collision.findClosePair();

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
			newtonCDBlock();

			if (inner_iteration_number >= min_inner_iteration - 1) {
				computeCurrentEnergy();
			}

			inner_iteration_number++;
			std::cout << "finish one itr " << inner_iteration_number<<" "<< energy << std::endl;

		}

		outer_itr_num++;
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



	//for (int i = 0; i < record_energy.size(); ++i) {
	//	std::cout << record_energy[i] << std::endl;
	//}
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

	if (abs(energy - previous_energy) < energy_converge_standard) {
		return true;
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


bool XPBD_IPC::checkMaxDisplacement()
{
	unsigned int num;
	std::array<double, 3>* current_pos;
	std::array<double, 3>* previous_pos;
	double* mass_;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		previous_pos = record_collision_free_vertex_position[i].data();
		current_pos = vertex_position[i];
		num = mesh_struct[i]->vertex_position.size();
		mass_ = mesh_struct[i]->mass_inv.data();
		for (unsigned int j = 0; j < num; ++j) {
			if (mass_[j] != 0.0) {
				for (unsigned int k = 0; k < 3; ++k) {
					if (abs(previous_pos[j][k] - current_pos[j][k]) > max_move_standard) {
						//std::cout << outer_itr_num<<" "<< abs(previous_pos[j][k] - current_pos[j][k]) << std::endl;

						return false;
					}
				}
			}
		}
	}
	return true;
}

//UPDATE_TET_HESSIAN 
void 	XPBD_IPC::tetHessian(int thread_No)
{
	unsigned int tet_around_tet_color_group_end;
	unsigned int tet_around_a_group_start;

	MatrixXd Hessian;

	Hessian.resize(4, 4);
	double stiffness = 0.0;
	Matrix<double, 3, 4>* A;
	double* volume;

	unsigned int prefix_sum_start;

	int k = 0;
	for (int i = 0; i < tetrahedron->size(); ++i) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		prefix_sum_start = prefix_sum_of_every_tet_index[i];

		A = tetrahedron->data()[i].mesh_struct.A.data();

		volume = tet_volume[i + cloth->size()];
		tet_around_a_group_start = tetrahedron->data()[i].mesh_struct.tetrahedron_index_begin_per_thread[thread_No];
		tet_around_tet_color_group_end = tetrahedron->data()[i].mesh_struct.tetrahedron_index_begin_per_thread[thread_No + 1];

		for (auto j = tet_around_a_group_start;
			j < tet_around_tet_color_group_end; ++j) {
			second_order_constraint.setARAPHessian(Hessian, stiffness, A[j],// ,
				volume[j]); //,
			memcpy(store_tet_arap_hessian.data() + 16 * (prefix_sum_start + j), Hessian.data(), 128);//
		}
	}

}

//UPDATE_TET_GRAD_SHARED
void XPBD_IPC::tetGradForColor(int thread_No, unsigned int color_No)
{

	unsigned int tet_around_tet_color_group_end;
	unsigned int* tet_around_a_group;

	MatrixXd grad;
	grad.resize(3, 4);
	std::array<int, 4>* tet_indices_;

	std::array<double, 3>* vertex_pos;
	double stiffness = 0.0;
	Matrix<double, 3, 4>* A;
	double* volume;

	unsigned int prefix_sum_start;

	int k=0;

	int color_group_index;

	for (int i = 0; i < tetrahedron->size(); ++i) {

		color_group_index = inner_iteration_number % tet_color_groups[i]->size();

		stiffness =0.5* tetrahedron->data()[0].ARAP_stiffness;
		prefix_sum_start = prefix_sum_of_every_tet_index[i];

		vertex_pos = vertex_position[i + cloth->size()];
		A = tetrahedron->data()[i].mesh_struct.A.data();

		volume = tet_volume[i + cloth->size()];
		tet_indices_ = tet_indices[i + cloth->size()];

		tet_around_a_group = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group[color_No].data();
		tet_around_tet_color_group_end = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group_start_per_thread[color_No][thread_No + 1];

		for (auto j = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group_start_per_thread[color_No][thread_No]; 	
			j < tet_around_tet_color_group_end; ++j) {
			second_order_constraint.setARAPGrad(grad, vertex_pos, stiffness, A[tet_around_a_group[j]],// ,
				volume[tet_around_a_group[j]], tet_indices_[tet_around_a_group[j]].data()); //,
			memcpy(store_tet_arap_grad.data() + 12 * (prefix_sum_start + tet_around_a_group[j]), grad.data(), 96);//
		}
	}
}



void XPBD_IPC::computeTetHessianInAdvance(int thread_No, int color_No)
{
	unsigned int tet_around_tet_color_group_end;
	unsigned int* tet_around_a_group;
	MatrixXd Hessian; 
	MatrixXd grad;
	Hessian.resize(4, 4);
	grad.resize(3, 4);


	std::array<double, 3>* vertex_pos;
	double stiffness=0.0;
	if (!tetrahedron->empty()) {
		stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		stiffness = cloth->data()[0].collision_stiffness[0];
	}
	Matrix<double, 3, 4>* A;
	double* volume;

	std::array<int, 4>*tet_indices_;
	unsigned int prefix_sum_start;



	for (int i = 0; i < tetrahedron->size(); ++i) {
		if (color_No >= unconnected_tet_index[i]->size()) {
			continue;
		}
		tet_around_tet_color_group_end = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group_start_per_thread[color_No][thread_No + 1];
		tet_around_a_group = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group[color_No].data();

		vertex_pos = vertex_position[i + cloth->size()];
		A = tetrahedron->data()[i].mesh_struct.A.data();

		volume = tet_volume[i + cloth->size()];
		tet_indices_ = tet_indices[i + cloth->size()];

		prefix_sum_start = prefix_sum_of_every_tet_index[i];

		//std::cout << "color " << color_No << " " << tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group_start_per_thread[color_No][thread_No] 
		//	<< " " << tet_around_tet_color_group_end<<" "<< prefix_sum_start << std::endl;

		for (int j = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group_start_per_thread[color_No][thread_No];
			j < tet_around_tet_color_group_end; ++j) {
			if (tet_around_tet_color_group_end > tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group[color_No].size()) {
				std::cout << "error to large size " << std::endl;
			}

			if (tet_around_a_group[j]>= tetrahedron->data()[i].mesh_struct.indices.size()) {
				std::cout << "erro  tet index " << std::endl;
			}

			second_order_constraint.setARAPGrad(grad, vertex_pos, stiffness, A[tet_around_a_group[j]],
				volume[tet_around_a_group[j]], tet_indices_[tet_around_a_group[j]].data());
			//memcpy(store_tet_arap_hessian.data() + 16 * (prefix_sum_start + tet_around_a_group[j]), Hessian.data(),128);
			if (12 * (prefix_sum_start + tet_around_a_group[j]) >= store_tet_arap_grad.size()) {
				std::cout<<"error "<<std::endl;
			}
			//memcpy(store_tet_arap_grad.data() + 12 * (prefix_sum_start + tet_around_a_group[j]), grad.data(), 96);
		}


		//for (int j = 0; j < tetrahedron->data()[i].mesh_struct.indices.size(); ++j) {
		//	second_order_constraint.setARAPHessianGrad(Hessian, grad, vertex_pos, stiffness, A[j],
		//		volume[j], tet_indices_[j].data());
		//	memcpy(store_tet_arap_hessian.data() + 16 * (prefix_sum_start + j), Hessian.data(), 128);
		//	memcpy(store_tet_arap_grad.data() + 12 * (prefix_sum_start + j), grad.data(), 96);
		//}
	}
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

	return checkMaxDisplacement();

}

void XPBD_IPC::updatePosition()
{

	unsigned int vertex_num;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_num = mesh_struct[i]->vertex_position.size();
		memcpy(vertex_position_render[i][0].data(), vertex_position[i][0].data(), 24 * vertex_num);
	}

	for (unsigned int i = 0; i < collider->size(); ++i) {
		memcpy(collider->data()[i].mesh_struct.vertex_for_render[0].data(), collider->data()[i].mesh_struct.vertex_position[0].data(), 24 * collider->data()[i].mesh_struct.vertex_position.size());
	}
}

void XPBD_IPC::computeCurrentEnergy()
{
	energy = 1e-15;
	double inertial_energy = 0.5 * computeInertialEnergy();
	energy += inertial_energy;
	double ARAP_energy = computeCurrentARAPEnergy();
	energy += ARAP_energy;
	double barrier_energy = 0.0;
	if (has_collider) {
		barrier_energy = computeBarrierEnergy();
		energy += barrier_energy;
	}
	record_energy.emplace_back(energy);

	if (energy > 1e9) {
		std::cout <<inner_iteration_number<<" "<< outer_itr_num<<"  inertial " << inertial_energy << " ARAP " << ARAP_energy << " barrier " << barrier_energy << std::endl;
	}
}

void XPBD_IPC::firstNewtonCD()
{
	previous_energy = energy;
	newtonCDTet();
}


//void XPBD_IPC::firstOnlyInertialCollision()
//{
//	unsigned int size;
//	MeshStruct* mesh_struct_;
//	double* volume;
//	std::array<double, 3>* vertex_pos;
//	std::array<double, 3>* last_step_vertex_pos;
//	std::array<double, 3>* initial_vertex_pos;
//	double* mass_inv;
//	std::array<double, 3>* sn_;
//	double* mass;
//	double colliision_stiffness;
//	int* vertex_index;
//	std::vector<bool>* is_surface_vertex;
//	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
//		mesh_struct_ = mesh_struct[i + cloth->size()];
//		size = tetrahedron->data()[i].mesh_struct.vertex_position.size();
//		volume = tetrahedron->data()[i].mesh_struct.volume.data();
//		vertex_pos = vertex_position[i + cloth->size()];
//		last_step_vertex_pos = initial_vertex_position[i + cloth->size()];
//		initial_vertex_pos = record_vertex_position[i + cloth->size()].data();
//		mass_inv = mesh_struct_->mass_inv.data();
//		mass = mesh_struct_->mass.data();
//		sn_ = sn[i + cloth->size()].data();
//		vertex_index = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
//		is_surface_vertex = &tetrahedron->data()[i].mesh_struct.vertex_on_surface;
//		colliision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
//		//std::cout << "=== " << std::endl;
//		for (unsigned int j = 0; j < size; ++j) {
//			if (mass_inv[j] != 0.0) {
//				solveInertialCollision(vertex_pos, initial_vertex_pos[j].data(),
//					last_step_vertex_pos[j].data(), sub_time_step, mass, volume, j, sn_,
//					colliision_stiffness, i + cloth->size(), (*is_surface_vertex)[j], vertex_index[j]);
//			}
//		}
//	}
//}


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
		last_step_vertex_pos = record_collision_free_vertex_position[i + cloth->size()].data();
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



double XPBD_IPC::computeBarrierEnergy()
{
	double energy = 0.0;

	energy+=computeVTCollisionEnergy(collision.vertex_triangle_pair_by_vertex, collision.vertex_triangle_pair_num_record,
		vertex_position.data(), vertex_position.data(), collision.close_vt_pair_num, false, triangle_indices.data());


	energy += computeEECollisionEnergy(collision.edge_edge_pair_by_edge, collision.edge_edge_pair_num_record, vertex_position.data(),
		vertex_position.data(), collision.close_ee_pair_num, edge_vertices.data(), edge_vertices.data());

	if (has_collider) {
		energy += computeVTCollisionEnergy(collision.vertex_obj_triangle_collider_pair_by_vertex, collision.vertex_obj_triangle_collider_num_record,
			vertex_position_collider.data(), vertex_position.data(), collision.close_vt_collider_pair_num, false, triangle_indices_collider.data());


		energy += computeVTCollisionEnergy(collision.triangle_vertex_collider_pair_by_triangle, collision.triangle_vertex_collider_pair_num_record,
			vertex_position.data(), vertex_position_collider.data(), collision.close_tv_collider_pair_num, true, triangle_indices.data());

		energy += computeEECollisionEnergy(collision.edge_edge_collider_pair_by_edge, collision.edge_edge_collider_pair_num_record, vertex_position.data(),
			vertex_position_collider.data(), collision.close_ee_collider_pair_num, edge_vertices.data(), collider_edge_vertices.data());
	}

	return energy;
}


double XPBD_IPC::computeEECollisionEnergy(unsigned int**edge_edge_pair_by_vertex_, unsigned int** edge_edge_pair_num_record_,
	std::array<double, 3>** edge_0_position, std::array<double, 3>** edge_1_position, unsigned int close_pair_num,
	unsigned int** edge_0_vertex, unsigned int** edge_1_vertex)
{
	double energy = 0.0;
	double	collision_stiffness;
	if (!tetrahedron->empty()) {
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	int size;
	unsigned int* edge_edge_pair_by_vertex;
	unsigned int* edge_edge_pair_num_record;


	unsigned int* pair_index;
	unsigned int pair_num;

	unsigned int* edge_vertex_index;
	unsigned int* edge_vertex_1_index;
	
		for (int i = 0; i < total_obj_num; ++i) {
			edge_edge_pair_by_vertex = edge_edge_pair_by_vertex_[i];
			edge_edge_pair_num_record = edge_edge_pair_num_record_[i];
			size = mesh_struct[i]->triangle_indices.size();

			for (int j = 0; j < size; ++j) {
				pair_index = edge_edge_pair_by_vertex + close_pair_num * j;
				pair_num = edge_edge_pair_num_record[j];
				edge_vertex_index = edge_0_vertex[i] + (j << 1);
				for (int k = 0; k < pair_num; k += 2) {
					edge_vertex_1_index = edge_1_vertex[pair_index[k]] + (pair_index[k + 1] << 1);
					energy+=compute_energy.computeBarrierEnergy(edge_0_position[i][edge_vertex_index[0]].data(),
						edge_0_position[i][edge_vertex_index[1]].data(),
						edge_1_position[pair_index[k]][edge_vertex_1_index[0]].data(),
						edge_1_position[pair_index[k]][edge_vertex_1_index[1]].data(), collision_stiffness, collision.d_hat_2, false);
				}
			}
		}
		return energy;
}



double XPBD_IPC::computeVTCollisionEnergy(unsigned int** vertex_triangle_pair_by_vertex_, unsigned int** vertex_triangle_pair_num_record_,
	std::array<double,3>** triangle_position, std::array<double,3>** vertex_position, unsigned int close_pair_num, bool is_TV,
	std::array<int, 3>** triangle_vertex)
{
	double energy=0.0;
	double	collision_stiffness;
	if (!tetrahedron->empty()) {
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	int size;
	unsigned int* vertex_triangle_pair_by_vertex;
	unsigned int* vertex_triangle_pair_num_record;


	unsigned int* pair_index;
	unsigned int pair_num;

	int* triangle_vertex_index;
	if (!is_TV) {
		for (int i = 0; i < total_obj_num; ++i) {
			vertex_triangle_pair_by_vertex = vertex_triangle_pair_by_vertex_[i];
			vertex_triangle_pair_num_record = vertex_triangle_pair_num_record_[i];
			if (i < cloth->size()) {
				size = mesh_struct[i]->vertex_position.size();
				for (int j = 0; j < size; ++j) {
					pair_index = vertex_triangle_pair_by_vertex + close_pair_num * j;
					pair_num = vertex_triangle_pair_num_record[j];
					for (int k = 0; k < pair_num; k += 2) {
						triangle_vertex_index = triangle_vertex[pair_index[k]][pair_index[k + 1]].data();
						energy+= compute_energy.computeBarrierEnergy(vertex_position[i][j].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[0]].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[1]].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[2]].data(), collision_stiffness, collision.d_hat_2, true);
					}
				}
			}
			else {
				vertex_triangle_pair_by_vertex = vertex_triangle_pair_by_vertex_[i];
				vertex_triangle_pair_num_record = vertex_triangle_pair_num_record_[i];
				int j;
				size = tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_index_on_sureface.size();
				for (int m = 0; m < size; ++m) {
					j = vertex_surface_to_global[i][m];
					pair_index = vertex_triangle_pair_by_vertex + close_pair_num * m;
					pair_num = vertex_triangle_pair_num_record[m];
					for (int k = 0; k < pair_num; k += 2) {
						triangle_vertex_index = triangle_vertex[pair_index[k]][pair_index[k + 1]].data();
						energy += compute_energy.computeBarrierEnergy(vertex_position[i][j].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[0]].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[1]].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[2]].data(), collision_stiffness, collision.d_hat_2, true);
					}
				}
			}
		}
	}
	else {
		for (int i = 0; i < total_obj_num; ++i) {
			vertex_triangle_pair_by_vertex = vertex_triangle_pair_by_vertex_[i];
			vertex_triangle_pair_num_record = vertex_triangle_pair_num_record_[i];
			size = mesh_struct[i]->triangle_indices.size();

			for (int j = 0; j < size; ++j) {
				pair_index = vertex_triangle_pair_by_vertex + close_pair_num * j;
				pair_num = vertex_triangle_pair_num_record[j];
				triangle_vertex_index = triangle_vertex[i][j].data();
				for (int k = 0; k < pair_num; k += 2) {
					energy += compute_energy.computeBarrierEnergy(vertex_position[pair_index[k]][pair_index[k+1]].data(),
						triangle_position[i][triangle_vertex_index[0]].data(),
						triangle_position[i][triangle_vertex_index[1]].data(),
						triangle_position[i][triangle_vertex_index[2]].data(), collision_stiffness, collision.d_hat_2, true);
				}
			}

		}
	}
	return energy;
}




void XPBD_IPC::newtonVTCollisionBlock()
{
	int size;
	int pair_num;
	unsigned int* pair_index;
	unsigned int* vertex_triangle_pair_by_vertex;
	unsigned int* vertex_triangle_pair_num_record;
	double stiffness = 0;
	double collision_stiffness =0.0;

	MeshStruct::Vertex* vertex;

	if (!tetrahedron->empty()) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}	
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	std::vector<unsigned int>tet_around;

	std::vector<unsigned int>* tet_around_vertex_;
	std::vector<unsigned int>* tet_around_triangle_;

	for (int i = 0; i <total_obj_num ; ++i) {
		vertex_triangle_pair_by_vertex = collision.vertex_triangle_pair_by_vertex[i];
		vertex_triangle_pair_num_record = collision.vertex_triangle_pair_num_record[i];
		vertex = mesh_struct[i]->vertices.data();
		if (i < cloth->size()) {
			size = mesh_struct[i]->vertex_position.size();
			tet_around_vertex_ = &tet_around;
			for (int j = 0; j < size; ++j) {
				pair_index = vertex_triangle_pair_by_vertex + collision.close_vt_pair_num * j;
				pair_num = vertex_triangle_pair_num_record[j];
				for (int k = 0; k < pair_num; k += 2) {
					if(pair_index[k]<cloth->size()){
						tet_around_triangle_ = &tet_around;;
					}
					else {
						tet_around_triangle_ = &tet_around_triangle[pair_index[k]][pair_index[k + 1]];
					}
					solveVT_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &vertex[j].face, &triangle_around_triangle[pair_index[k]][pair_index[k + 1]],
						&vertex[j].edge, &edge_around_triangle[pair_index[k]][pair_index[k + 1]], tet_around_vertex_, tet_around_triangle_, collision.d_hat_2,
						false,false);
				}				
			}
		}
		else {
			int j;
			size = tetrahedron->data()[i-cloth->size()].mesh_struct.vertex_index_on_sureface.size();
			for (int m = 0; m < size; ++m) {
				j = vertex_surface_to_global[i][m];
				tet_around_vertex_ = &tet_around_vertex[i][j];
				pair_index = vertex_triangle_pair_by_vertex + collision.close_vt_pair_num * m;
				pair_num = vertex_triangle_pair_num_record[m];
				for (int k = 0; k < pair_num; k += 2) {
					if (pair_index[k] < cloth->size()) {
						tet_around_triangle_ = &tet_around;
					}
					else {
						tet_around_triangle_ = &tet_around_triangle[pair_index[k]][pair_index[k + 1]];
					}
					solveVT_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &vertex[j].face, &triangle_around_triangle[pair_index[k]][pair_index[k + 1]],
						&vertex[j].edge, &edge_around_triangle[pair_index[k]][pair_index[k + 1]], tet_around_vertex_, tet_around_triangle_, collision.d_hat_2,
						false,false);
				}
			}		
		}		
	}
}


void XPBD_IPC::newtonTVColliderCollisionBlock()
{
	int size;
	int pair_num;
	unsigned int* pair_index;
	unsigned int* vertex_triangle_pair_by_triangle;
	unsigned int* vertex_triangle_pair_num_record;
	double stiffness = 0;
	double collision_stiffness = 0.0;

	MeshStruct::Vertex* vertex;

	if (!tetrahedron->empty()) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	std::vector<unsigned int>tet_around;

	std::vector<unsigned int>* tet_around_v_;
	std::vector<unsigned int>* tet_around_t_;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_triangle_pair_by_triangle = collision.triangle_vertex_collider_pair_by_triangle[i];
		vertex_triangle_pair_num_record = collision.triangle_vertex_collider_pair_num_record[i];
		size = mesh_struct[i]->triangle_indices.size();
		if (i < cloth->size()) {
			tet_around_t_ = &tet_around;
			tet_around_v_ = &tet_around;
			for (unsigned int j = 0; j < size; ++j) {
				pair_index = vertex_triangle_pair_by_triangle + collision.close_tv_collider_pair_num * j;
				pair_num = vertex_triangle_pair_num_record[j];
				for (unsigned int k = 0; k < pair_num; k += 2) {
					solveVT_collisionBlock(pair_index[k], pair_index[k + 1], i, j,
						stiffness, sub_time_step, collision_stiffness, &tet_around , &triangle_around_triangle[i][j],
						&tet_around, &edge_around_triangle[i][j], tet_around_v_, tet_around_t_, collision.d_hat_2,
						true, false);
				}
			}
		}
		else {
			tet_around_v_ = &tet_around;
			for (unsigned int j = 0; j < size; ++j) {
				pair_index = vertex_triangle_pair_by_triangle + collision.close_tv_collider_pair_num * j;
				pair_num = vertex_triangle_pair_num_record[j];
				tet_around_t_ = &tet_around_triangle[i][j];
				for (unsigned int k = 0; k < pair_num; k += 2) {
					solveVT_collisionBlock(pair_index[k], pair_index[k + 1], i, j,
						stiffness, sub_time_step, collision_stiffness, &tet_around, &triangle_around_triangle[i][j],
						&tet_around, &edge_around_triangle[i][j], tet_around_v_, tet_around_t_, collision.d_hat_2,
						true, false);
				}
			}
			
		}
	}
}


void XPBD_IPC::newtonVTColliderCollisionBlock()
{
	int size;
	int pair_num;
	unsigned int* pair_index;
	unsigned int* vertex_triangle_pair_by_vertex;
	unsigned int* vertex_triangle_pair_num_record;
	double stiffness = 0;
	double collision_stiffness = 0.0;

	MeshStruct::Vertex* vertex;

	if (!tetrahedron->empty()) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	std::vector<unsigned int>tet_around;

	std::vector<unsigned int>* tet_around_vertex_;
	std::vector<unsigned int>* tet_around_triangle_;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_triangle_pair_by_vertex = collision.vertex_obj_triangle_collider_pair_by_vertex[i];
		vertex_triangle_pair_num_record = collision.vertex_obj_triangle_collider_num_record[i];
		vertex = mesh_struct[i]->vertices.data();
		if (i < cloth->size()) {
			size = mesh_struct[i]->vertex_position.size();
			tet_around_vertex_ = &tet_around;
			tet_around_triangle_ = &tet_around;
			for (unsigned int j = 0; j < size; ++j) {
				pair_index = vertex_triangle_pair_by_vertex + collision.close_vt_collider_pair_num * j;
				pair_num = vertex_triangle_pair_num_record[j];
				for (unsigned int k = 0; k < pair_num; k += 2) {
					solveVT_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &vertex[j].face, &triangle_around_triangle[pair_index[k]][pair_index[k + 1]],
						&vertex[j].edge, &edge_around_triangle[pair_index[k]][pair_index[k + 1]], tet_around_vertex_, tet_around_triangle_, collision.d_hat_2,
						false,true);
				}
			}
		}
		else {
			int j;
			tet_around_triangle_ = &tet_around;
			size = tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_index_on_sureface.size();
			for (unsigned int m = 0; m < size; ++m) {
				j = vertex_surface_to_global[i][m];
				tet_around_vertex_ = &tet_around_vertex[i][j];
				pair_index = vertex_triangle_pair_by_vertex + collision.close_vt_collider_pair_num * m;
				pair_num = vertex_triangle_pair_num_record[m];
				for (unsigned int k = 0; k < pair_num; k += 2) {

					solveVT_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &vertex[j].face, &tet_around,
						&vertex[j].edge, &tet_around, tet_around_vertex_, tet_around_triangle_,
						collision.d_hat_2,false ,true);
				}
			}
		}
	}
}



void XPBD_IPC::newtonEEColliderCollisionBlock()
{
	int size;
	int pair_num;
	unsigned int* pair_index;
	unsigned int* edge_edge_pair_by_edge;
	unsigned int* edge_edge_pair_num_record;
	double stiffness = 0;
	double collision_stiffness = 0.0;



	if (!tetrahedron->empty()) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	std::vector<unsigned int>* tet_around_edge_0;

	std::vector<unsigned int>* tet_around_edge_1;

	std::vector<unsigned int>* triangle_around_edge_;
	std::vector<unsigned int>* edge_around_edge_;

	std::vector<unsigned int>tet_around;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		edge_edge_pair_by_edge = collision.edge_edge_collider_pair_by_edge[i];
		edge_edge_pair_num_record = collision.edge_edge_collider_pair_num_record[i];
		triangle_around_edge_ = triangle_around_edge[i];
		edge_around_edge_ = edge_around_edge[i];
		size = mesh_struct[i]->edge_length.size();
		if (i < cloth->size()) {
			tet_around_edge_0 = &tet_around;
			tet_around_edge_1 = &tet_around;
			for (unsigned int j = 0; j < size; ++j) {
				pair_index = edge_edge_pair_by_edge + collision.close_ee_collider_pair_num * j;
				pair_num = edge_edge_pair_num_record[j];
				for (int k = 0; k < pair_num; k += 2) {
					solveEE_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &triangle_around_edge_[j], &tet_around,
						&edge_around_edge_[j], &tet_around,
						tet_around_edge_0, tet_around_edge_1, collision.d_hat_2, false, true);
				}
			}
		}
		else {
			tet_around_edge_1 = &tet_around;
			for (unsigned int j = 0; j < size; ++j) {

				tet_around_edge_0 = &tet_around_edge[i][j];
				pair_index = edge_edge_pair_by_edge + collision.close_ee_collider_pair_num * j;
				pair_num = edge_edge_pair_num_record[j];
				for (int k = 0; k < pair_num; k += 2) {
					solveEE_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &triangle_around_edge_[j], &tet_around,
						&edge_around_edge_[j], &tet_around,
						tet_around_edge_0, tet_around_edge_1, collision.d_hat_2, false, true);
				}
			}
		}
	}

}


void XPBD_IPC::newtonEECollisionBlock()
{
	int size;
	int pair_num;
	unsigned int* pair_index;
	unsigned int* edge_edge_pair_by_edge;
	unsigned int* edge_edge_pair_num_record;
	double stiffness = 0;
	double collision_stiffness = 0.0;



	if (!tetrahedron->empty()) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	std::vector<unsigned int>* tet_around_edge_0;

	std::vector<unsigned int>* tet_around_edge_1;

	std::vector<unsigned int>* triangle_around_edge_;
	std::vector<unsigned int>* edge_around_edge_;

	std::vector<unsigned int>tet_around;

	for (int i = 0; i < total_obj_num; ++i) {
		edge_edge_pair_by_edge = collision.edge_edge_pair_by_edge[i];
		edge_edge_pair_num_record = collision.edge_edge_pair_num_record[i];
		triangle_around_edge_ = triangle_around_edge[i];
		edge_around_edge_ = edge_around_edge[i];
		size = mesh_struct[i]->edge_length.size();
		if (i < cloth->size()) {
			tet_around_edge_0 = &tet_around;
			for (int j = 0; j < size; ++j) {
				pair_index = edge_edge_pair_by_edge + collision.close_ee_pair_num * j;
				pair_num = edge_edge_pair_num_record[j];
				for (int k = 0; k < pair_num; k += 2) {
					if (i > pair_index[k]) {
						continue;
					}
					else if (i == pair_index[k]) {
						if (j > pair_index[k + 1]) {
							continue;
						}
					}
					if (pair_index[k] < cloth->size()) {
						tet_around_edge_1 = &tet_around;
					}
					else {
						tet_around_edge_1 = &tet_around_edge[pair_index[k]][pair_index[k + 1]];
					}
					solveEE_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &triangle_around_edge_[j], &triangle_around_edge[pair_index[k]][pair_index[k + 1]],
						&edge_around_edge_[j], &edge_around_edge[pair_index[k]][pair_index[k + 1]],
						tet_around_edge_0, tet_around_edge_1, collision.d_hat_2,false,false);
				}
			}
		}
		else {
			for (int j = 0; j < size; ++j) {
				tet_around_edge_0 = &tet_around_edge[i][j];
				pair_index = edge_edge_pair_by_edge + collision.close_ee_pair_num * j;
				pair_num = edge_edge_pair_num_record[j];
				for (int k = 0; k < pair_num; k += 2) {
					if (i > pair_index[k]) {
						continue;
					}
					else if (i == pair_index[k]) {
						if (j > pair_index[k + 1]) {
							continue;
						}
					}
					if (pair_index[k] < cloth->size()) {
						tet_around_edge_1 = &tet_around;
					}
					else {
						tet_around_edge_1 = &tet_around_edge[pair_index[k]][pair_index[k + 1]];
					}
					solveEE_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &triangle_around_edge_[j], &triangle_around_edge[pair_index[k]][pair_index[k + 1]],
						&edge_around_edge_[j], &edge_around_edge[pair_index[k]][pair_index[k + 1]],
						tet_around_edge_0, tet_around_edge_1, collision.d_hat_2, false, false);
				}
			}
		}
	}
}






void XPBD_IPC::newtonCDBlock()
{
	//newtonEECollisionBlock();
	//newtonVTCollisionBlock();
	//if (has_collider) {
	//	newtonEEColliderCollisionBlock();
	//	newtonVTColliderCollisionBlock();
	//	newtonTVColliderCollisionBlock();
	//}
	newtonCDTetBlock();
}




void XPBD_IPC::solveNewtonCD_tetBlock()
{
	//for (int i = 0; i < max_tet_color_num; ++i) {
	//	for (int j = 0; j < thread->thread_num; ++j) {
	//		tetHessianForColor(j, i);
	//	}
	//}
	//std::cout << store_tet_arap_hessian.size() << std::endl;
	//for (int i = 0; i < store_tet_arap_hessian.size(); ++i) {
	//	if (store_tet_arap_hessian[i] != 0.0) {
	//		system("pause");
	//		//std::cout << "error hessian not equal to zero " << std::endl;
	//	}
	//}

	//for (int i = 0; i < store_tet_arap_grad.size(); ++i) {
	//	if (store_tet_arap_grad[i] != 0.0) {
	//		system("pause");
	//	}
	//}
	//for (int i = 0; i < max_tet_color_num; ++i) {		//
	//	newtonCDTetBlockTest(i);
	//}
	//newtonCDTetBlock();

	std::vector<double>time;
	for (int i = 0; i < max_tet_color_num; ++i) {		//
	
		//update shared tet, collision hessian 
		thread->assignTask(this, UPDATE_TET_GRAD_SHARED, i);	
		//for (int j = 0; j < total_thread_num; ++j) {
		//	tetGradForColor(j, i);
		//}
		//tetHessianForColor(0, i);
		//testComputeTetHessian(0, i);
		if (perform_collision) {
			collision.computeHessian(i);
		}		
	
		//auto t0 = std::chrono::system_clock::now();

		thread->assignTask(this, SOLVE_TET_BLOCK, i);

		//for (unsigned int j = 0; j < total_thread_num; ++j) {
		//	newtonCDTetBlockAGroup(j, i);
		//}	

	//	//local ccd
		if (perform_collision) {
			collision.collisionTimeColor(i);
		}
		//auto t1 = std::chrono::system_clock::now();
	
		//double time_0 = double(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den * 1000.0;
		////double time_1 = double(std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den * 1000.0;
		//time.push_back(time_0);
	}
	
	//for (int i = 0; i < time.size(); ++i) {
	//	std::cout << time[i] << " ";
	//}
	//std::cout << std::endl;

	//system("pause");
	//newtonCDTetBlockAGroupTest();

}



//SOLVE_TET_BLOCK
void XPBD_IPC::newtonCDTetBlockAGroup(int thread_No, int color)
{

	unsigned int* tet_group;
	int start, end;

	std::array<int, 4>* indices;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;

	double stiffness;
	double collision_stiffness;;
	Matrix<double, 3, 4>* A;
	std::array<double, 3>* sn_;
	double* mass;
	unsigned int* unfixed_vertex_num;
	std::vector<unsigned int>* tet_neightbor_tet;
	std::vector<unsigned int>* tet_neightbor_tet_common_vertex;
	std::array<int, 4>* unfixed_vertex_index;
	std::array<int, 4>* unfixed_actual_vertex_index;


	std::vector<unsigned int>* triangle_of_a_tet;
	std::vector<unsigned int>* edge_of_a_tet;
	unsigned int obj_No;
	int* vertex_index_on_surface;
	int j;

	double* hessian_record;
	double* grad_record;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		if (color >= tetrahedron->data()[i].mesh_struct.unconnected_tet_index.size()) {
			continue;
		}
		collision_stiffness = tetrahedron->data()[i].collision_stiffness[0];

		tet_group = mesh_struct[i + cloth->size()]->unconnected_tet_index[color].data();
		end = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread[color][thread_No + 1];
		start = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread[color][thread_No];

		mesh_struct_ = mesh_struct[i + cloth->size()];
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i + cloth->size()];
		stiffness = 0.5 * tetrahedron->data()[i].ARAP_stiffness;

		A = tetrahedron->data()[i].mesh_struct.A.data();
		mass = mesh_struct_->mass.data();
		sn_ = sn[i + cloth->size()].data();
		unfixed_vertex_num = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();
		tet_neightbor_tet = tetrahedron->data()[i].mesh_struct.tet_tet_index.data();
		tet_neightbor_tet_common_vertex = tetrahedron->data()[i].mesh_struct.tet_neighbor_tet_vertex_order.data();
		unfixed_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_indices.data();
		unfixed_actual_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_actual_indices.data();
		triangle_of_a_tet = tetrahedron->data()[i].mesh_struct.triangle_index_of_a_tet.data();
		edge_of_a_tet = tetrahedron->data()[i].mesh_struct.edge_index_of_a_tet.data();
		obj_No = i + cloth->size();
		vertex_index_on_surface = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();

		hessian_record = store_tet_arap_hessian.data() + prefix_sum_of_every_tet_index[i] * 16;
		grad_record = store_tet_arap_grad.data() + prefix_sum_of_every_tet_index[i] * 12;


		//std::cout << thread_No<<" "<< color<<" " << start << " " << end << std::endl;

		for (int k = start; k < end; ++k) {
			j = tet_group[k];
			solveTetBlock(vertex_pos, stiffness, sub_time_step, mass, A, tet_neightbor_tet[j],
				volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
				unfixed_vertex_index[j].data(), unfixed_vertex_num[j], &(triangle_of_a_tet[j]), &(edge_of_a_tet[j]), collision_stiffness,
				obj_No, unfixed_actual_vertex_index[j].data(), vertex_index_on_surface, hessian_record, grad_record);
			//
			//
		}
	}
}


//void XPBD_IPC::newtonCDTetBlockAGroupTest(int  thread_No, int color)
//{
//	unsigned int* tet_group;
//	int start, end;
//
//	std::array<int, 4>* indices;
//	MeshStruct* mesh_struct_;
//	double* volume;
//	std::array<double, 3>* vertex_pos;
//
//	double stiffness;
//	double collision_stiffness;;
//	Matrix<double, 3, 4>* A;
//	std::array<double, 3>* sn_;
//	double* mass;
//	unsigned int* unfixed_vertex_num;
//	std::vector<unsigned int>* tet_neightbor_tet;
//	std::vector<unsigned int>* tet_neightbor_tet_common_vertex;
//	std::array<int, 4>* unfixed_vertex_index;
//	std::array<int, 4>* unfixed_actual_vertex_index;
//
//
//	std::vector<unsigned int>* triangle_of_a_tet;
//	std::vector<unsigned int>* edge_of_a_tet;
//	unsigned int obj_No;
//	int* vertex_index_on_surface;
//
//	std::array<double, 3>* ori_vertex_pos;
//
//	int end_index;
//
//	for (int i = 0; i < tetrahedron->size(); ++i) {
//
//		collision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
//
//		mesh_struct_ = mesh_struct[i + cloth->size()];
//		indices = tetrahedron->data()[i].mesh_struct.indices.data();
//		volume = tetrahedron->data()[i].mesh_struct.volume.data();
//		vertex_pos = vertex_position[i + cloth->size()];
//		stiffness = 0.5 * tetrahedron->data()[i].ARAP_stiffness;
//
//		A = tetrahedron->data()[i].mesh_struct.A.data();
//		mass = mesh_struct_->mass.data();
//		sn_ = sn[i + cloth->size()].data();
//		unfixed_vertex_num = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();
//		tet_neightbor_tet = tetrahedron->data()[i].mesh_struct.tet_tet_index.data();
//		tet_neightbor_tet_common_vertex = tetrahedron->data()[i].mesh_struct.tet_neighbor_tet_vertex_order.data();
//		unfixed_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_indices.data();
//		unfixed_actual_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_actual_indices.data();
//		triangle_of_a_tet = tetrahedron->data()[i].mesh_struct.triangle_index_of_a_tet.data();
//		edge_of_a_tet = tetrahedron->data()[i].mesh_struct.edge_index_of_a_tet.data();
//		obj_No = i + cloth->size();
//		vertex_index_on_surface = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
//		ori_vertex_pos = initial_vertex_position[i + cloth->size()];
//		int j;
//		end_index = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread[color][thread_No + 1];
//		for (int m = 0; m < end_index; ++m) {
//			j = tetrahedron->data()[i].mesh_struct.unconnected_tet_index[color][m];
//			solveTetBlock(vertex_pos, stiffness, sub_time_step, mass, A, tet_neightbor_tet[j],
//				volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
//				unfixed_vertex_index[j].data(), unfixed_vertex_num[j], &(triangle_of_a_tet[j]), &(edge_of_a_tet[j]), collision_stiffness,
//				obj_No, unfixed_actual_vertex_index[j].data(), vertex_index_on_surface, store_tet_arap_hessian.data() + prefix_sum_of_every_tet_index[i] * 16, 
//				store_tet_arap_grad.data() + prefix_sum_of_every_tet_index[i] * 12);//
//			//
//		}
//
//
//
//
//		//for (unsigned int j = 0; j < tetrahedron->data()[i].mesh_struct.indices.size(); ++j) {
//		//	solveNewtonCD_tetBlock(vertex_pos, stiffness, sub_time_step, mass, A, tet_neightbor_tet[j],
//		//		indices, volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
//		//		unfixed_vertex_index[j].data(), unfixed_vertex_num[j], &(triangle_of_a_tet[j]), &(edge_of_a_tet[j]), collision_stiffness,
//		//		obj_No, unfixed_actual_vertex_index[j].data(), vertex_index_on_surface, ori_vertex_pos);
//		//}
//	}
//}

void XPBD_IPC::newtonCDTetBlockTest(int color_No)
{
	unsigned int size;
	std::array<int, 4>* indices;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* ori_vertex_pos;
	double* mass_inv;
	double stiffness;
	double collision_stiffness;;
	Matrix<double, 3, 4>* A;
	std::array<double, 3>* sn_;
	double* mass;
	unsigned int* unfixed_vertex_num;
	std::vector<unsigned int>* tet_neightbor_tet;
	std::vector<unsigned int>* tet_neightbor_tet_common_vertex;
	std::array<int, 4>* unfixed_vertex_index;
	std::array<int, 4>* unfixed_actual_vertex_index;


	std::vector<unsigned int>* triangle_of_a_tet;
	std::vector<unsigned int>* edge_of_a_tet;
	unsigned int obj_No;
	int* vertex_index_on_surface;


	int j;

	unsigned int* color;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		collision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
		mesh_struct_ = mesh_struct[i + cloth->size()];
		size = tetrahedron->data()[i].mesh_struct.indices.size();
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i + cloth->size()];
		stiffness = 0.5 * tetrahedron->data()[i].ARAP_stiffness;

		A = tetrahedron->data()[i].mesh_struct.A.data();
		mass_inv = mesh_struct_->mass_inv.data();
		mass = mesh_struct_->mass.data();
		sn_ = sn[i + cloth->size()].data();
		unfixed_vertex_num = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();
		tet_neightbor_tet = tetrahedron->data()[i].mesh_struct.tet_tet_index.data();
		tet_neightbor_tet_common_vertex = tetrahedron->data()[i].mesh_struct.tet_neighbor_tet_vertex_order.data();
		unfixed_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_indices.data();
		unfixed_actual_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_actual_indices.data();
		triangle_of_a_tet = tetrahedron->data()[i].mesh_struct.triangle_index_of_a_tet.data();
		edge_of_a_tet = tetrahedron->data()[i].mesh_struct.edge_index_of_a_tet.data();
		obj_No = i + cloth->size();
		vertex_index_on_surface = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
		ori_vertex_pos = record_vertex_position[i + cloth->size()].data();

		color = tetrahedron->data()[i].mesh_struct.unconnected_tet_index[color_No].data();

		for (unsigned int k = 0; k < tetrahedron->data()[i].mesh_struct.unconnected_tet_index[color_No].size(); ++k) {
			j = color[k];
			solveNewtonCD_tetBlock(vertex_pos, stiffness, sub_time_step, mass, A, tet_neightbor_tet[j],
				indices, volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
				unfixed_vertex_index[j].data(), unfixed_vertex_num[j], &(triangle_of_a_tet[j]), &(edge_of_a_tet[j]), collision_stiffness,
				obj_No, unfixed_actual_vertex_index[j].data(), vertex_index_on_surface, ori_vertex_pos,
				store_tet_arap_hessian.data()+16*prefix_sum_of_every_tet_index[i]);
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
	std::array<double, 3>* ori_vertex_pos;
	double* mass_inv;
	double stiffness;
	double collision_stiffness;;
	Matrix<double, 3, 4>* A;
	std::array<double, 3>* sn_;
	double* mass;
	unsigned int* unfixed_vertex_num;
	std::vector<unsigned int>* tet_neightbor_tet;
	std::vector<unsigned int>* tet_neightbor_tet_common_vertex;
	std::array<int,4>* unfixed_vertex_index;
	std::array<int,4>* unfixed_actual_vertex_index;


	std::vector<unsigned int>* triangle_of_a_tet;
	std::vector<unsigned int>* edge_of_a_tet;
	unsigned int obj_No;
	int* vertex_index_on_surface;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		collision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
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
		unfixed_actual_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_actual_indices.data();
		triangle_of_a_tet = tetrahedron->data()[i].mesh_struct.triangle_index_of_a_tet.data();
		edge_of_a_tet = tetrahedron->data()[i].mesh_struct.edge_index_of_a_tet.data();
		obj_No = i + cloth->size();
		vertex_index_on_surface = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
		ori_vertex_pos = record_vertex_position[i + cloth->size()].data();
		for (unsigned int j = 0; j < size; ++j) {
			solveNewtonCD_tetBlock(vertex_pos, stiffness, sub_time_step, mass, A, tet_neightbor_tet[j],
				indices, volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
				unfixed_vertex_index[j].data(), unfixed_vertex_num[j],&(triangle_of_a_tet[j]),&(edge_of_a_tet[j]), collision_stiffness,
				obj_No, unfixed_actual_vertex_index[j].data(), vertex_index_on_surface, ori_vertex_pos,
				store_tet_arap_hessian.data()+16* prefix_sum_of_every_tet_index[i]);
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


double XPBD_IPC::getCollisionTime(std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet,
	int* tet_actual_unfixed_vertex_indices,
	int unfixed_tet_vertex_num,
	int** vertex_index_on_surface,
	std::array<double, 3>** current_vertex_position,
	std::array<double, 3>** initial_vertex_position)
{
	double collision_time = 1.0;
	int* triangle_;		unsigned int* edge_;
	int surface_index;

	for (int i = 0; i < unfixed_tet_vertex_num; i+=2) {
		surface_index = vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i + 1]];
		if (surface_index == -1) {
			continue;
		}
		collision.VTCollisionTimeOneVertex(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data(),
			current_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data(), collision_time,
			collision.vertex_triangle_pair_num_record[tet_actual_unfixed_vertex_indices[i]][surface_index],
			collision.vertex_triangle_pair_by_vertex[tet_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * surface_index,
			initial_vertex_position, vertex_position.data(), triangle_indices.data());

		if (has_collider) {
			collision.VTCollisionTimeOneVertex(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data(),
				current_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data(), collision_time,
				collision.vertex_obj_triangle_collider_num_record[tet_actual_unfixed_vertex_indices[i]][surface_index],
				collision.vertex_obj_triangle_collider_pair_by_vertex[tet_actual_unfixed_vertex_indices[i]] + 
				collision.close_vt_collider_pair_num * surface_index,
				vertex_position_collider.data(), vertex_position_collider.data(), triangle_indices_collider.data());
		}
		if (floor->exist) {
			collision.floorCollisionTime(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data(),
				current_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data(), floor->dimension,
				floor->normal_direction, floor->value, collision_time, collision.tolerance);
		}
	}

	for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); i+=2) {
		triangle_ = triangle_indices[*i][*(i+1)].data();
		collision.TVCollisionTimeOneTriangle(initial_vertex_position[*i][triangle_[0]].data(), initial_vertex_position[*i][triangle_[1]].data(),
			initial_vertex_position[*i][triangle_[2]].data(),
			current_vertex_position[*i][triangle_[0]].data(), current_vertex_position[*i][triangle_[1]].data(),
			current_vertex_position[*i][triangle_[2]].data(), collision_time,
			collision.triangle_vertex_pair_num_record[*i][*(i+1)],
			collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i+1)),
			initial_vertex_position, vertex_position.data(),3);
		if (has_collider) {
			collision.TVCollisionTimeOneTriangle(initial_vertex_position[*i][triangle_[0]].data(), initial_vertex_position[*i][triangle_[1]].data(),
				initial_vertex_position[*i][triangle_[2]].data(),
				current_vertex_position[*i][triangle_[0]].data(), current_vertex_position[*i][triangle_[1]].data(),
				current_vertex_position[*i][triangle_[2]].data(), collision_time,
				collision.triangle_vertex_collider_pair_num_record[*i][*(i + 1)],
				collision.triangle_vertex_collider_pair_by_triangle[*i] + collision.close_tv_collider_pair_num * (*(i + 1)),
				vertex_position_collider.data(), vertex_position_collider.data(),2);
		}
	}

	for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); i+=2) {
		edge_ = edge_vertices[*i] + ((*(i+1)) << 1);
		collision.EECollisionTimeOneEdgeAll(initial_vertex_position[*i][*edge_].data(),
			initial_vertex_position[*i][*(edge_ + 1)].data(),
			current_vertex_position[*i][*edge_].data(),
			current_vertex_position[*i][*(edge_ + 1)].data(), collision_time,
			collision.edge_edge_pair_num_record[*i][*(i+1)],
			collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i+1)),
			initial_vertex_position, vertex_position.data(), edge_vertices.data());
		if (has_collider) {
			collision.EECollisionTimeOneEdgeAll(initial_vertex_position[*i][*edge_].data(),
				initial_vertex_position[*i][*(edge_ + 1)].data(),
				current_vertex_position[*i][*edge_].data(),
				current_vertex_position[*i][*(edge_ + 1)].data(), collision_time,
				collision.edge_edge_collider_pair_num_record[*i][*(i + 1)],
				collision.edge_edge_collider_pair_by_edge[*i] + collision.close_ee_collider_pair_num * (*(i + 1)),
				vertex_position_collider.data(), vertex_position_collider.data(), collider_edge_vertices.data());
		}
	}
	return collision_time;
}


double XPBD_IPC::getCollisionTime(std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet,
	unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	int unfixed_tet_vertex_num,
	int* vertex_index_on_surface, std::array<double, 3>* current_vertex_position,
	std::array<double, 3>* initial_vertex_position)
{
	double collision_time = 1.0;
	int* triangle_;		unsigned int* edge_;

		for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
			if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
				continue;
			}
			collision.VTCollisionTimeOneVertex(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]].data(),
				current_vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), collision_time,
				collision.vertex_triangle_pair_num_record[obj_No][vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]]],
				collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]],
				this->record_collision_free_vertex_position_address.data(), vertex_position.data(), triangle_indices.data());
			if (has_collider) {
				collision.VTCollisionTimeOneVertex(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]].data(),
					current_vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), collision_time,
					collision.vertex_obj_triangle_collider_num_record[obj_No][vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]]],
					collision.vertex_obj_triangle_collider_pair_by_vertex[obj_No] + collision.close_vt_collider_pair_num * vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]],
					vertex_position_collider.data(), vertex_position_collider.data(), triangle_indices_collider.data());
			}
			if (floor->exist) {
				collision.floorCollisionTime(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]].data(),
					current_vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), floor->dimension,
					floor->normal_direction, floor->value, collision_time, collision.tolerance);
			}
		}

		for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
			triangle_ = triangle_indices[obj_No][*i].data();
			collision.TVCollisionTimeOneTriangle(initial_vertex_position[triangle_[0]].data(), initial_vertex_position[triangle_[1]].data(),
				initial_vertex_position[triangle_[2]].data(),
				current_vertex_position[triangle_[0]].data(), current_vertex_position[triangle_[1]].data(),
				current_vertex_position[triangle_[2]].data(), collision_time,
				collision.triangle_vertex_pair_num_record[obj_No][*i],
				collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * (*i),
				this->record_collision_free_vertex_position_address.data(), vertex_position.data(),3);
			if (has_collider) {
				collision.TVCollisionTimeOneTriangle(initial_vertex_position[triangle_[0]].data(), initial_vertex_position[triangle_[1]].data(),
					initial_vertex_position[triangle_[2]].data(),
					current_vertex_position[triangle_[0]].data(), current_vertex_position[triangle_[1]].data(),
					current_vertex_position[triangle_[2]].data(), collision_time,
					collision.triangle_vertex_collider_pair_num_record[obj_No][*i],
					collision.triangle_vertex_collider_pair_by_triangle[obj_No] + collision.close_tv_collider_pair_num * (*i),
					vertex_position_collider.data(), vertex_position_collider.data(),2);
			}
		}

		for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
			edge_ = edge_vertices[obj_No] + ((*i) << 1);
			collision.EECollisionTimeOneEdgeAll(initial_vertex_position[*edge_].data(),
				initial_vertex_position[*(edge_ + 1)].data(),
				current_vertex_position[*edge_].data(),
				current_vertex_position[*(edge_ + 1)].data(), collision_time,
				collision.edge_edge_pair_num_record[obj_No][*i],
				collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * (*i),
				this->record_collision_free_vertex_position_address.data(), vertex_position.data(), edge_vertices.data());
			if (has_collider) {
				collision.EECollisionTimeOneEdgeAll(initial_vertex_position[*edge_].data(),
					initial_vertex_position[*(edge_ + 1)].data(),
					current_vertex_position[*edge_].data(),
					current_vertex_position[*(edge_ + 1)].data(), collision_time,
					collision.edge_edge_collider_pair_num_record[obj_No][*i],
					collision.edge_edge_collider_pair_by_edge[obj_No] + collision.close_ee_collider_pair_num * (*i),
					vertex_position_collider.data(), vertex_position_collider.data(), collider_edge_vertices.data());
			}
		}

		return collision_time;
}


void XPBD_IPC::getCollisionBlockTetHessian(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* tets, double stiffness, int* pair_actual_unfixed_vertex_indices, //pair_actual_unfixed_vertex_indices first obj, second primitive index
	int unfixed_vertex_num)//unfixed_vertex_num size double )
{
	for (auto i = tets->begin(); i < tets->end(); i += 2) {
		getARAPCollisionHessianForPair(Hessian, grad, stiffness, *i, *(i + 1), tet_indices[*i][*(i + 1)].data(),
			pair_actual_unfixed_vertex_indices, unfixed_vertex_num, vertex_position[*i], tet_A[*i], tet_volume[*i]);
	}
}

void XPBD_IPC::getCollisionBlockCollisionHessianTest(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangles,
	std::vector<unsigned int>* edges,
	double collision_stiffness, int* pair_actual_unfixed_vertex_indices, //pair_actual_unfixed_vertex_indices first obj, second primitive index
	int unfixed_vertex_num, double d_hat_2, int** vertex_index_on_surface, unsigned int vertex_obj_No, unsigned int tri_obj_No)//unfixed_vertex_num size double
{
	if (has_collider) {
		for (int i = 0; i < unfixed_vertex_num; i += 2) {
			getVTCollisionHessainForPairTest(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]].data(),
				collision_stiffness, collision.vertex_triangle_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_triangle_pair_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]],
				i >> 1, vertex_obj_No, tri_obj_No, pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				collision.vertex_obj_triangle_collider_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_collider_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_obj_triangle_collider_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (int i = 0; i < unfixed_vertex_num; i += 2) {
			getVTCollisionHessainForPairTest(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]].data(),
				collision_stiffness, collision.vertex_triangle_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_triangle_pair_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]],
				i >> 1, vertex_obj_No, tri_obj_No, pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				emp, 0);
		}
	}
	int* triangle_;
	if (has_collider) {
		for (auto i = triangles->begin(); i < triangles->end(); i += 2) {
			triangle_ = triangle_indices[*i][*(i + 1)].data();
			getTVCollisionHessainForPairTest(Hessian, grad, vertex_position[*i][triangle_[0]].data(), vertex_position[*i][triangle_[1]].data(),
				vertex_position[*i][triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i + 1)),
				collision.triangle_vertex_pair_num_record[*i][*(i + 1)], *i,
				vertex_obj_No, tri_obj_No, triangle_, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
				d_hat_2,
				collision.triangle_vertex_collider_pair_by_triangle[*i] + collision.close_tv_collider_pair_num * (*(i + 1)),
				collision.triangle_vertex_collider_pair_num_record[*i][*(i + 1)]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = triangles->begin(); i < triangles->end(); i += 2) {
			triangle_ = triangle_indices[*i][*(i + 1)].data();
			getTVCollisionHessainForPairTest(Hessian, grad, vertex_position[*i][triangle_[0]].data(), vertex_position[*i][triangle_[1]].data(),
				vertex_position[*i][triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i + 1)),
				collision.triangle_vertex_pair_num_record[*i][*(i + 1)], *i,
				vertex_obj_No, tri_obj_No, triangle_, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
				d_hat_2,
				emp, 0);
		}
	}
	unsigned int* edge_;
	if (has_collider) {
		for (auto i = edges->begin(); i < edges->end(); i += 2) {
			edge_ = edge_vertices[*i] + ((*(i + 1)) << 1);
			getEECollisionHessainForPairTest(Hessian, grad, *i, edge_,
				collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i + 1)),
				collision.edge_edge_pair_num_record[*i][*(i + 1)],
				pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				vertex_obj_No, tri_obj_No,
				i - edges->begin(), edges, vertex_position[*i][*edge_].data(),
				vertex_position[*i][*(edge_ + 1)].data(), collision_stiffness, collision.edge_edge_collider_pair_by_edge[*i] + collision.close_ee_collider_pair_num * (*(i + 1)),
				collision.edge_edge_collider_pair_num_record[*i][*(i + 1)], rest_edge_length[*i][*(i + 1)]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = edges->begin(); i < edges->end(); i += 2) {
			edge_ = edge_vertices[*i] + ((*(i + 1)) << 1);
			getEECollisionHessainForPairTest(Hessian, grad, *i, edge_,
				collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i + 1)),
				collision.edge_edge_pair_num_record[*i][*(i + 1)],
				pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				vertex_obj_No, tri_obj_No,
				i - edges->begin(), edges, vertex_position[*i][*edge_].data(),
				vertex_position[*i][*(edge_ + 1)].data(), collision_stiffness, emp, 0, rest_edge_length[*i][*(i + 1)]);
		}
	}

	if (floor->exist) {
		for (int i = 0; i < unfixed_vertex_num; i += 2) {
			getFloorHessianForTet(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]].data(),
				floor->value,
				floor->dimension, collision_stiffness, floor->normal_direction, collision.d_hat_2, i >> 1, unfixed_vertex_num >> 1);
		}
	}
}


void XPBD_IPC::getCollisionBlockCollisionHessian(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangles,
	std::vector<unsigned int>* edges,
	double collision_stiffness, int* pair_actual_unfixed_vertex_indices, //pair_actual_unfixed_vertex_indices first obj, second primitive index
	int unfixed_vertex_num, double d_hat_2, int** vertex_index_on_surface, unsigned int vertex_obj_No, unsigned int tri_obj_No)//unfixed_vertex_num size double
{
	if (has_collider) {
		for (int i = 0; i < unfixed_vertex_num; i += 2) {
			getVTCollisionHessainForPair(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]].data(),
				collision_stiffness, collision.vertex_triangle_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_triangle_pair_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]],
				i >> 1, vertex_obj_No, tri_obj_No, pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				collision.vertex_obj_triangle_collider_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_collider_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_obj_triangle_collider_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (int i = 0; i < unfixed_vertex_num; i += 2) {
			getVTCollisionHessainForPair(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]].data(),
				collision_stiffness, collision.vertex_triangle_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_triangle_pair_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]],
				i >> 1, vertex_obj_No, tri_obj_No, pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				emp, 0);
		}
	}
	int* triangle_;
	if (has_collider) {
		for (auto i = triangles->begin(); i < triangles->end(); i += 2) {
			triangle_ = triangle_indices[*i][*(i + 1)].data();
			getTVCollisionHessainForPair(Hessian, grad, vertex_position[*i][triangle_[0]].data(), vertex_position[*i][triangle_[1]].data(),
				vertex_position[*i][triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i + 1)),
				collision.triangle_vertex_pair_num_record[*i][*(i + 1)], *i,
				vertex_obj_No, tri_obj_No, triangle_, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
				d_hat_2,
				collision.triangle_vertex_collider_pair_by_triangle[*i] + collision.close_tv_collider_pair_num * (*(i + 1)),
				collision.triangle_vertex_collider_pair_num_record[*i][*(i + 1)]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = triangles->begin(); i < triangles->end(); i += 2) {
			triangle_ = triangle_indices[*i][*(i + 1)].data();
			getTVCollisionHessainForPair(Hessian, grad, vertex_position[*i][triangle_[0]].data(), vertex_position[*i][triangle_[1]].data(),
				vertex_position[*i][triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i + 1)),
				collision.triangle_vertex_pair_num_record[*i][*(i + 1)], *i,
				vertex_obj_No, tri_obj_No, triangle_, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
				d_hat_2,
				emp,0);
		}
	}
	unsigned int* edge_;
	if (has_collider) {
		for (auto i = edges->begin(); i < edges->end(); i += 2) {
			edge_ = edge_vertices[*i] + ((*(i + 1)) << 1);
			getEECollisionHessainForPair(Hessian, grad, *i, edge_,
				collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i + 1)),
				collision.edge_edge_pair_num_record[*i][*(i + 1)],
				pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				vertex_obj_No, tri_obj_No,
				i - edges->begin(), edges, vertex_position[*i][*edge_].data(),
				vertex_position[*i][*(edge_ + 1)].data(), collision_stiffness, collision.edge_edge_collider_pair_by_edge[*i] + collision.close_ee_collider_pair_num * (*(i + 1)),
				collision.edge_edge_collider_pair_num_record[*i][*(i + 1)], rest_edge_length[*i][*(i+1)]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = edges->begin(); i < edges->end(); i += 2) {
			edge_ = edge_vertices[*i] + ((*(i + 1)) << 1);
			getEECollisionHessainForPair(Hessian, grad, *i, edge_,
				collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i + 1)),
				collision.edge_edge_pair_num_record[*i][*(i + 1)],
				pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				vertex_obj_No, tri_obj_No,
				i - edges->begin(), edges, vertex_position[*i][*edge_].data(),
				vertex_position[*i][*(edge_ + 1)].data(), collision_stiffness,emp,	0, rest_edge_length[*i][*(i + 1)]);
		}
	}

	if (floor->exist) {
		for (int i = 0; i < unfixed_vertex_num; i+=2) {
			getFloorHessianForTet(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i+1]].data(), 
				floor->value,
				floor->dimension, collision_stiffness, floor->normal_direction, collision.d_hat_2, i>>1, unfixed_vertex_num>>1);
		}
	}
}



void XPBD_IPC::getCollisionHessianFromRecord(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet,
	double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	int unfixed_tet_vertex_num, 
	int* vertex_index_on_surface, unsigned int* vt_prefix_sum, unsigned int* ee_prefix_sum, unsigned int* tv_collider_prefix_sum, unsigned int* ee_collider_prefix_sum, 
	std::array<double, 3>* vertex_position)
{
	int surface_index;
	if (has_collider) {
		for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
			if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
				continue;
			}
			surface_index = vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]];
			getVTCollisionHessainForTetFromRecord(Hessian, grad, collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * surface_index,
				collision.vertex_triangle_pair_num_record[obj_No][surface_index],
				i, obj_No, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num,
				collision.vertex_obj_triangle_collider_pair_by_vertex[obj_No] + collision.close_vt_collider_pair_num * surface_index,
				collision.vertex_obj_triangle_collider_num_record[obj_No][surface_index],
				collision.vt_hessian_record.data() + 144 * vt_prefix_sum[surface_index],
				collision.vt_grad_record.data() + 12 * vt_prefix_sum[surface_index],
				collision.vt_hessian_record_index.data() + 5 * vt_prefix_sum[surface_index],
				collision.vt_colldier_hessian_record.data() + (collision.vertex_num_on_surface_prefix_sum[obj_No] + surface_index) * 9,
				collision.vt_colldier_grad_record.data() + (collision.vertex_num_on_surface_prefix_sum[obj_No] + surface_index) * 3);
		}
	}
	else {
		unsigned int emp[1] = { 0 }; 
		double temp_double[1] = { 0.0 };
		for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
			if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
				continue;
			}
			surface_index = vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]];
			getVTCollisionHessainForTetFromRecord(Hessian, grad, collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * surface_index,
				collision.vertex_triangle_pair_num_record[obj_No][surface_index],
				i, obj_No, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num,
				emp, 0,
				collision.vt_hessian_record.data() + 144 * vt_prefix_sum[surface_index],
				collision.vt_grad_record.data() + 12 * vt_prefix_sum[surface_index],
				collision.vt_hessian_record_index.data() + 5 * vt_prefix_sum[surface_index],
				temp_double, temp_double);
		}
	}

	int* triangle_;
	if (has_collider) {
		for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
			triangle_ = triangle_indices[obj_No][*i].data();
			getTVCollisionHessainForTetFromRecord(Hessian, grad,
				collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * (*i),
				collision.triangle_vertex_pair_num_record[obj_No][*i], obj_No, triangle_, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num,
				collision.triangle_vertex_collider_pair_by_triangle[obj_No] + collision.close_tv_collider_pair_num * (*i),
				collision.triangle_vertex_collider_pair_num_record[obj_No][*i],
				collision.tv_colldier_hessian_record_index.data() + 4 * tv_collider_prefix_sum[*i],
				collision.tv_colldier_hessian_record.data() + 81 * tv_collider_prefix_sum[*i],
				collision.tv_colldier_grad_record.data() + 9 * tv_collider_prefix_sum[*i]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		double temp_double[1] = { 0.0 }; int temp_int[1] = { 0 };
		for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
			triangle_ = triangle_indices[obj_No][*i].data();
			getTVCollisionHessainForTetFromRecord(Hessian, grad,
				collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * (*i),
				collision.triangle_vertex_pair_num_record[obj_No][*i], obj_No, triangle_, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num,
				emp,0,
				temp_int,temp_double, temp_double);
		}
	}


	unsigned int* edge_;
	if (has_collider) {
		for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
			edge_ = edge_vertices[obj_No] + ((*i) << 1);
			getEECollisionHessainForTetFromRecord(Hessian, grad, obj_No, edge_,
				collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * (*i),
				collision.edge_edge_pair_num_record[obj_No][*i],
				tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, 
				i - edge_of_a_tet->begin(), edge_of_a_tet, collision.edge_edge_collider_pair_by_edge[obj_No] + collision.close_ee_collider_pair_num * (*i),
				collision.edge_edge_collider_pair_num_record[obj_No][*i],
				collision.ee_hessian_record.data()+144*ee_prefix_sum[*i], collision.ee_grad_record.data()+12*ee_prefix_sum[*i],
				collision.ee_hessian_record_index.data()+5*ee_prefix_sum[*i],
				collision.ee_collider_hessian_record_index.data()+3*ee_collider_prefix_sum[*i],
				collision.ee_collider_hessian_record.data()+36*ee_collider_prefix_sum[*i],
				collision.ee_collider_grad_record.data()+6*ee_collider_prefix_sum[*i]	);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		double temp_double[1] = { 0.0 }; int temp_int[1] = { 0 };
		for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
			edge_ = edge_vertices[obj_No] + ((*i) << 1);

			getEECollisionHessainForTetFromRecord(Hessian, grad, obj_No, edge_,
				collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * (*i),
				collision.edge_edge_pair_num_record[obj_No][*i],
				tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num,
				i - edge_of_a_tet->begin(), edge_of_a_tet, emp,0,
				collision.ee_hessian_record.data() + 144 * ee_prefix_sum[*i], collision.ee_grad_record.data() + 12 * ee_prefix_sum[*i],
				collision.ee_hessian_record_index.data() + 5 * ee_prefix_sum[*i],
				temp_int, temp_double,temp_double);
		}
	}

	//if (Hessian.norm() != 0.0) {
	//	std::cout << "actual " << std::endl;
	//	std::cout << grad.transpose() << std::endl;
	//	std::cout << Hessian << std::endl;
	//}


	if (floor->exist) {
		for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
			if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
				continue;
			}
			grad[3 * i + floor->dimension] += collision.floor_grad_record[vertex_num_on_surface_prefix_sum[obj_No]+ vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]]];
			Hessian.data()[(3 * i + floor->dimension) * (3 * unfixed_tet_vertex_num + 1)] += collision.floor_hessian_record[vertex_num_on_surface_prefix_sum[obj_No] + vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]]];
		}
	}
}




void XPBD_IPC::getCollisionHessian(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet,
	double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	int unfixed_tet_vertex_num, double d_hat_2,
	int* vertex_index_on_surface, std::array<double, 3>* vertex_position)
{
	if (has_collider) {
		for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
			if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
				continue;
			}
			getVTCollisionHessainForTet(Hessian, grad, vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), collision_stiffness,
				collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]],
				collision.vertex_triangle_pair_num_record[obj_No][vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]]],
				i, obj_No, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				collision.vertex_obj_triangle_collider_pair_by_vertex[obj_No] + collision.close_vt_collider_pair_num * vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]],
				collision.vertex_obj_triangle_collider_num_record[obj_No][vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]]]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
			if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
				continue;
			}
			getVTCollisionHessainForTet(Hessian, grad, vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), collision_stiffness,
				collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]],
				collision.vertex_triangle_pair_num_record[obj_No][vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]]],
				i, obj_No, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				emp, 0);
		}
	}
	
	int* triangle_;
	if (has_collider) {
		for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
			triangle_ = triangle_indices[obj_No][*i].data();
			getTVCollisionHessainForTet(Hessian, grad, vertex_position[triangle_[0]].data(), vertex_position[triangle_[1]].data(),
				vertex_position[triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * (*i),
				collision.triangle_vertex_pair_num_record[obj_No][*i], obj_No,
				triangle_, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				collision.triangle_vertex_collider_pair_by_triangle[obj_No] + collision.close_tv_collider_pair_num * (*i),
				collision.triangle_vertex_collider_pair_num_record[obj_No][*i]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
			triangle_ = triangle_indices[obj_No][*i].data();
			getTVCollisionHessainForTet(Hessian, grad, vertex_position[triangle_[0]].data(), vertex_position[triangle_[1]].data(),
				vertex_position[triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * (*i),
				collision.triangle_vertex_pair_num_record[obj_No][*i], obj_No,
				triangle_, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				emp, 0);
		}
	}




	unsigned int* edge_;
	if (has_collider) {
		for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
			edge_ = edge_vertices[obj_No] + ((*i) << 1);
			getEECollisionHessainForTet(Hessian, grad, obj_No, edge_,
				collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * (*i),
				collision.edge_edge_pair_num_record[obj_No][*i],
				tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				i - edge_of_a_tet->begin(), edge_of_a_tet, vertex_position[*edge_].data(),
				vertex_position[*(edge_ + 1)].data(), collision_stiffness, collision.edge_edge_collider_pair_by_edge[obj_No] + collision.close_ee_collider_pair_num * (*i), 
				collision.edge_edge_collider_pair_num_record[obj_No][*i],rest_edge_length[obj_No][*i]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
			edge_ = edge_vertices[obj_No] + ((*i) << 1);
			getEECollisionHessainForTet(Hessian, grad, obj_No, edge_,
				collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * (*i),
				collision.edge_edge_pair_num_record[obj_No][*i],
				tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				i - edge_of_a_tet->begin(), edge_of_a_tet, vertex_position[*edge_].data(),
				vertex_position[*(edge_ + 1)].data(), collision_stiffness, emp,
				0,rest_edge_length[obj_No][*i]);
		}
	}

	//if (Hessian.norm() != 0.0) {
	//	std::cout << "compare correct " << std::endl;
	//	std::cout << grad.transpose() << std::endl;
	//	std::cout << Hessian << std::endl;
	//}


	if (floor->exist) {
		for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
			if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
				continue;
			}
			getFloorHessianForTet(Hessian, grad, vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), floor->value,
				floor->dimension, collision_stiffness, floor->normal_direction, collision.d_hat_2,i, unfixed_tet_vertex_num);
		}		
	}

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
		if (second_order_constraint.getCollisionPairHessian(pos_0, pos_1, pos_2,
			vertex_position[TV[i]][TV[i + 1]].data(),
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


void XPBD_IPC::getEECollisionHessainForPairTest(MatrixXd& Hessian, VectorXd& grad, unsigned int ee_obj_No, unsigned int* edge_vertex_index,
	unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2,
	unsigned int obj_No_0, unsigned int obj_No_1,
	int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, double* ea0, double* ea1, double stiffness,
	unsigned int* EE_collider, int num_collider, double edge_length_0)//edge_order_in_tet should x2
{
	unsigned int* edge_vertex;
	int vertex_order_in_tet[4];
	memset(vertex_order_in_tet, 0xff, 8);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex_index, ee_obj_No, vertex_order_in_tet);
	for (int i = 0; i < num; i += 2) {
		double barrier_=0.0;
		edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
		memset(vertex_order_in_tet + 2, 0xff, 8);
		if (EE[i] == obj_No_0 || EE[i] == obj_No_1) {
			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, EE[i + 1], EE[i])) {
				continue;
			}
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex, EE[i], vertex_order_in_tet + 2);
		}
		second_order_constraint.computeEEBarrierGradientHessianTest(ea0, ea1, vertex_position[EE[i]][*edge_vertex].data(),
			vertex_position[EE[i]][edge_vertex[1]].data(), Hessian, grad,
			vertex_order_in_tet, stiffness, d_hat_2, edge_length_0, rest_edge_length[EE[i]][EE[i + 1]],
			barrier_);

		if (barrier_ != 0.0) {
			std::cout << "EE " << barrier_ << std::endl;
		}

	}

	if (has_collider) {
		memset(vertex_order_in_tet + 2, 0xff, 8);
		for (int i = 0; i < num_collider; i += 2) {
			edge_vertex = collider_edge_vertices[EE_collider[i]] + (EE_collider[i + 1] << 1);
			second_order_constraint.computeEEBarrierGradientHessian(ea0, ea1, vertex_position_collider[EE_collider[i]][*edge_vertex].data(),
				vertex_position_collider[EE_collider[i]][edge_vertex[1]].data(), Hessian, grad,
				vertex_order_in_tet, stiffness, d_hat_2, edge_length_0,
				rest_edge_length_collider[EE_collider[i]][EE_collider[i + 1]]);
		}
	}
}





void XPBD_IPC::getEECollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, unsigned int ee_obj_No, unsigned int* edge_vertex_index,
	unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2,
	unsigned int obj_No_0, unsigned int obj_No_1,
	int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, double* ea0, double* ea1, double stiffness,
	unsigned int* EE_collider, int num_collider, double edge_length_0)//edge_order_in_tet should x2
{
	unsigned int* edge_vertex;
	int vertex_order_in_tet[4];
	memset(vertex_order_in_tet, 0xff, 8);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex_index,ee_obj_No, vertex_order_in_tet);
	for (int i = 0; i < num; i += 2) {
		edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
		memset(vertex_order_in_tet + 2, 0xff, 8);
		if (EE[i] == obj_No_0 || EE[i]==obj_No_1) {
			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, EE[i + 1],EE[i])) {
				continue;
			}
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex, EE[i], vertex_order_in_tet + 2);
		}
		second_order_constraint.computeEEBarrierGradientHessian(ea0, ea1, vertex_position[EE[i]][*edge_vertex].data(),
			vertex_position[EE[i]][edge_vertex[1]].data(), Hessian, grad,
			vertex_order_in_tet, stiffness, d_hat_2, edge_length_0,rest_edge_length[EE[i]][EE[i+1]]);

	}

	if (has_collider) {
		memset(vertex_order_in_tet + 2, 0xff, 8);
		for (int i = 0; i < num_collider; i += 2) {
			edge_vertex = collider_edge_vertices[EE_collider[i]] + (EE_collider[i + 1] << 1);
			second_order_constraint.computeEEBarrierGradientHessian(ea0, ea1, vertex_position_collider[EE_collider[i]][*edge_vertex].data(),
				vertex_position_collider[EE_collider[i]][edge_vertex[1]].data(), Hessian, grad,
				vertex_order_in_tet, stiffness, d_hat_2, edge_length_0,
				rest_edge_length_collider[EE_collider[i]][EE_collider[i+1]]);
		}
	}
}

void XPBD_IPC::getARAPCollisionHessianForPair(MatrixXd& Hessian, VectorXd& grad, double stiffness, int tet_obj, int tet_index, int* tet_vertex, int* tet_unfixed_vertex_indices,
	int unfixed_tet_vertex_num, std::array<double,3>* vertex_position, Matrix<double, 3, 4>* A, double* volume)
{
	int triangle_vertex_order_in_tet[4];
	memset(triangle_vertex_order_in_tet, 0xff, 16);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, tet_vertex, tet_obj, triangle_vertex_order_in_tet, 4);
	second_order_constraint.solveTetCertainVertices(vertex_position, stiffness, A[tet_index], triangle_vertex_order_in_tet, tet_vertex, Hessian,
		volume[tet_index], grad);

}


void XPBD_IPC::getTVCollisionHessainForPairTest(MatrixXd& Hessian, VectorXd& grad, double* t0, double* t1, double* t2, double stiffness,
	unsigned int* TV, int num, unsigned int tri_obj, unsigned int obj_No_0, unsigned int obj_No_1, int* triangle_indices,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* TV_collider, int collider_num)
{
	int triangle_vertex_order_in_tet[4];
	memset(triangle_vertex_order_in_tet, 0xff, 16);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_indices, tri_obj, triangle_vertex_order_in_tet + 1, 3);

	for (int i = 0; i < num; i += 2) {
		if (TV[i] == obj_No_0 || TV[i] == obj_No_1) {
			if (vertexInPair(unfixed_tet_vertex_num, TV[i + 1], tet_unfixed_vertex_indices, TV[i])) {
				continue;
			}
		}

		double barrier = 0.0;

		second_order_constraint.computeVTBarrierGradientHessianTest(Hessian, grad, vertex_position[TV[i]][TV[i + 1]].data(),
			t0, t1, t2,
			d_hat_2, triangle_vertex_order_in_tet, stiffness, barrier);
		if (barrier != 0.0) {
			std::cout << "TV " << barrier << std::endl;
		}
	}

	if (has_collider) {
		for (int i = 0; i < collider_num; i += 2) {
			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_collider[TV_collider[i]][TV_collider[i + 1]].data(),
				t0, t1, t2,
				d_hat_2, triangle_vertex_order_in_tet, stiffness);
		}
	}

}




void XPBD_IPC::getTVCollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, double* t0, double* t1, double* t2, double stiffness,
	unsigned int* TV, int num, unsigned int tri_obj, unsigned int obj_No_0, unsigned int obj_No_1, int* triangle_indices,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* TV_collider, int collider_num)
{
	int triangle_vertex_order_in_tet[4];
	memset(triangle_vertex_order_in_tet, 0xff, 16);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_indices, tri_obj, triangle_vertex_order_in_tet+1,3);

	for (int i = 0; i < num; i += 3) {
		if (TV[i] == obj_No_0 || TV[i] == obj_No_1) {
			if (vertexInPair(unfixed_tet_vertex_num, TV[i + 1], tet_unfixed_vertex_indices,TV[i])) {
				continue;
			}
		}
		second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position[TV[i]][TV[i + 1]].data(),
			t0, t1, t2,
			d_hat_2, triangle_vertex_order_in_tet, stiffness);
	}

	if (has_collider) {
		for (int i = 0; i < collider_num; i += 2) {
			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_collider[TV_collider[i]][TV_collider[i + 1]].data(),
				t0, t1, t2,
				d_hat_2, triangle_vertex_order_in_tet, stiffness);
		}
	}

}

void XPBD_IPC::getVTCollisionHessainForPairTest(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int vertex_obj_No, unsigned int tri_obj_No,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* VT_collider, int num_collider)
{
	int* triangle_vertex;
	int triangle_vertex_order_in_tet[4];
	triangle_vertex_order_in_tet[0] = vertex_order_in_matrix;
	
	for (int i = 0; i < num; i += 2) {
		double barrier = 0.0;
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		if (VT[i] == tri_obj_No || VT[i] == vertex_obj_No) {
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, VT[i], triangle_vertex_order_in_tet + 1, 3);
		}
		second_order_constraint.computeVTBarrierGradientHessianTest(Hessian, grad, vertex_position_, vertex_position[VT[i]][triangle_vertex[0]].data(),
			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(),
			d_hat_2, triangle_vertex_order_in_tet, stiffness, barrier);

		if (barrier != 0.0) {
			std::cout << barrier << std::endl;
		}

	}

	if (has_collider) {
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		for (int i = 0; i < num_collider; i += 2) {
			triangle_vertex = triangle_indices_collider[VT_collider[i]][VT_collider[i + 1]].data();
			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position_collider[VT_collider[i]][triangle_vertex[0]].data(),
				vertex_position_collider[VT_collider[i]][triangle_vertex[1]].data(), vertex_position_collider[VT_collider[i]][triangle_vertex[2]].data(),
				d_hat_2, triangle_vertex_order_in_tet, stiffness);
		}
	}

}



void XPBD_IPC::getVTCollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int vertex_obj_No, unsigned int tri_obj_No,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* VT_collider, int num_collider)
{
	int* triangle_vertex;
	int triangle_vertex_order_in_tet[4];
	triangle_vertex_order_in_tet[0] = vertex_order_in_matrix;
	for (int i = 0; i < num; i += 2) {
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		if (VT[i] == tri_obj_No || VT[i] == vertex_obj_No) {
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, VT[i], triangle_vertex_order_in_tet+1,3);
		}
		second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position[VT[i]][triangle_vertex[0]].data(),
			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(),
			d_hat_2, triangle_vertex_order_in_tet, stiffness);
	}

	if (has_collider) {
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		for (int i = 0; i < num_collider; i += 2) {
			triangle_vertex = triangle_indices_collider[VT_collider[i]][VT_collider[i + 1]].data();
			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position_collider[VT_collider[i]][triangle_vertex[0]].data(),
				vertex_position_collider[VT_collider[i]][triangle_vertex[1]].data(), vertex_position_collider[VT_collider[i]][triangle_vertex[2]].data(),
				d_hat_2, triangle_vertex_order_in_tet, stiffness);
		}
	}
}


void XPBD_IPC::getTVCollisionHessainForTetFromRecord(MatrixXd& Hessian, VectorXd& grad,
	unsigned int* TV, int num, unsigned int obj_No, int* triangle_indices,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, unsigned int* TV_collider, int collider_num,
	int* tv_collider_hessian_record_index, double* tv_collider_hessian_record, double* tv_collider_grad_record)
{
	int triangle_vertex_order_in_tet[4];
	memset(triangle_vertex_order_in_tet, 0xff, 16);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_indices,
		triangle_vertex_order_in_tet + 1, 3);
	int pair_index;
	int* vt_record_index;

	for (int i = 0; i < num; i += 3) {
		if (TV[i] == obj_No) {
			if (vertexInTet(unfixed_tet_vertex_num, TV[i + 1], tet_unfixed_vertex_indices)) {
				continue;
			}
		}

		if (TV[i] < cloth->size()) {
			pair_index = collision.vertex_triangle_pair_num_record_prefix_sum[TV[i]][TV[i + 1]] + TV[i + 2];
		}
		else {
			pair_index = collision.vertex_triangle_pair_num_record_prefix_sum[TV[i]][vertex_index_surface[TV[i]][TV[i + 1]]]+TV[i + 2];
		}
		vt_record_index = collision.vt_hessian_record_index.data() + 5 * pair_index;
		if (*vt_record_index != 0) {
			setTetHessianFromBarrierHessian(Hessian, grad.data(), collision.vt_hessian_record.data() + 144 * pair_index, 
				collision.vt_grad_record.data() + 12 * pair_index, triangle_vertex_order_in_tet,
				vt_record_index + 1, *vt_record_index);
		}
	}
	if (has_collider) {
		for (int i = 0; i < collider_num; i += 2) {
			if ((*tv_collider_hessian_record_index) != 0) {
				setTetHessianFromBarrierHessian(Hessian, grad.data(), tv_collider_hessian_record, tv_collider_grad_record, triangle_vertex_order_in_tet,
					tv_collider_hessian_record_index + 1, *tv_collider_hessian_record_index);
			}
			tv_collider_hessian_record_index += 4;
			tv_collider_hessian_record += 9;
			tv_collider_grad_record += 81;
		}
	}

}


void XPBD_IPC::getVTCollisionHessainForTetFromRecord(MatrixXd& Hessian, VectorXd& grad,
	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int obj_No, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num,
	unsigned int* VT_collider, int num_collider, double* vt_hessian_record, double* vt_grad_record, int* vt_hessian_record_index,
	double* vt_collider_hessian_record, double* vt_collider_grad_record)
{
	int* triangle_vertex;
	int triangle_vertex_order_in_tet[4];
	triangle_vertex_order_in_tet[0] = vertex_order_in_matrix;
	for (unsigned int i = 0; i < num; i += 2) {
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		if (VT[i] == obj_No) {
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, triangle_vertex_order_in_tet + 1, 3);
		}
		if ((*vt_hessian_record_index) != 0) {
			setTetHessianFromBarrierHessian(Hessian, grad.data(), vt_hessian_record, vt_grad_record, triangle_vertex_order_in_tet,
				vt_hessian_record_index + 1, *vt_hessian_record_index);
		}
		vt_hessian_record_index += 5;
		vt_grad_record += 12;
		vt_hessian_record += 144;
	}

	if (has_collider) {
		int vertex_pos = 3 * vertex_order_in_matrix * (Hessian.cols() + 1);
		for (unsigned int i = 0; i < 3; ++i) {
			grad.data()[3* vertex_order_in_matrix + i] += vt_collider_grad_record[i];
			Hessian.data()[vertex_pos + i] += vt_collider_hessian_record[i];
			Hessian.data()[vertex_pos + Hessian.cols() + i] += vt_collider_hessian_record[i+3];
			Hessian.data()[vertex_pos + Hessian.cols() + Hessian.cols()+i] += vt_collider_hessian_record[i+6];
		}
	}
}



void XPBD_IPC::setTetHessianFromBarrierHessian(MatrixXd& Hessian_system, double* grad_system, double* Hessian_, double* grad_,
	int* triangle_vertex_order_in_system, int* vertex_in_pair, int vertex_in_use)
{
	double* sys_start;
	double* record_start;

	int size_sys = Hessian_system.cols();
	int size_record = vertex_in_use * 3;

	for (int i = 0; i < vertex_in_use; ++i) {
		if (triangle_vertex_order_in_system[vertex_in_pair[i]] != -1) {
			for (int j = 0; j < vertex_in_use; ++j) {
				if (triangle_vertex_order_in_system[vertex_in_pair[j]] != -1) {
					sys_start = &Hessian_system(3 * triangle_vertex_order_in_system[vertex_in_pair[i]], 3 * triangle_vertex_order_in_system[vertex_in_pair[j]]);
					record_start = Hessian_ + 9 * j * vertex_in_use + 3 * i;

					*sys_start += *record_start;
					*(sys_start + 1) += *(record_start + 1);
					*(sys_start + 2) += *(record_start + 2);

					*(sys_start + size_sys) += *(record_start + size_record);
					*(sys_start + size_sys + 1) += *(record_start + size_record + 1);
					*(sys_start + size_sys + 2) += *(record_start + size_record + 2);

					*(sys_start + (size_sys<<1)) += *(record_start + (size_record<<1));
					*(sys_start + (size_sys << 1) + 1) += *(record_start + (size_record << 1) + 1);
					*(sys_start + (size_sys << 1) + 2) += *(record_start + (size_record << 1) + 2);
				}

			}
			sys_start = grad_system + 3 * triangle_vertex_order_in_system[vertex_in_pair[i]];
			record_start = grad_ + 3 * i;

			*sys_start += *record_start;
			*(sys_start + 1) += *(record_start + 1);
			*(sys_start + 2) += *(record_start + 2);
		}
	}
}


void XPBD_IPC::getVTCollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int obj_No,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* VT_collider, int num_collider)
{
	int* triangle_vertex;
	int triangle_vertex_order_in_tet[4];
	triangle_vertex_order_in_tet[0] = vertex_order_in_matrix;
	for (unsigned int i = 0; i < num; i += 2) {
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		memset(triangle_vertex_order_in_tet +1, 0xff, 12);
		if (VT[i] == obj_No) {
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, triangle_vertex_order_in_tet+1,3);
		}
		second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position[VT[i]][triangle_vertex[0]].data(),
			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(),
			d_hat_2, triangle_vertex_order_in_tet, stiffness);
	}
	if (has_collider) {
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		for (unsigned int i = 0; i < num_collider; i += 2) {
			triangle_vertex = triangle_indices_collider[VT_collider[i]][VT_collider[i + 1]].data();
			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position_collider[VT_collider[i]][triangle_vertex[0]].data(),
				vertex_position_collider[VT_collider[i]][triangle_vertex[1]].data(), vertex_position_collider[VT_collider[i]][triangle_vertex[2]].data(),
				d_hat_2, triangle_vertex_order_in_tet, stiffness);
		}
	}
}


//void XPBD_IPC::getVTCollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
//	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, 
//	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2)
//{
//	int* triangle_vertex;
//	int triangle_vertex_order_in_tet[4];
//	triangle_vertex_order_in_tet[0] = vertex_order_in_matrix;
//	for (unsigned int i = 0; i < num; i += 2) {
//		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
//		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
//		checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, VT[i], triangle_vertex_order_in_tet);		
//		second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position[VT[i]][triangle_vertex[0]].data(),
//			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(),
//			d_hat_2, triangle_vertex_order_in_tet, stiffness);
//	}
//}




void XPBD_IPC::getFloorHessianForTet(MatrixXd& Hessian, VectorXd& grad, double* vertex_position, double floor_value,
	unsigned int dimension, double collision_stiffness, bool direction, double d_hat_2, unsigned int vertex_order_in_matrix, unsigned int unfixed_vertex_num)
{
	double distance = (vertex_position[dimension] - floor_value) * (vertex_position[dimension] - floor_value);

	if (distance > d_hat_2) {
		return;
	}
	
	double h, g;
	barrierGradHessian(distance, d_hat_2, g, h);

	//if (vertex_order_in_matrix == 1) {
		//std::cout <<"vertex_order_in_matrix "<< vertex_order_in_matrix<<" "<< g << " " << h << std::endl;
	//}

	double grad_d, h_d;
	grad_d = 2 * (vertex_position[dimension] - floor_value);
	h_d = 2;
	grad[3 * vertex_order_in_matrix + dimension] += collision_stiffness * g * grad_d;
	Hessian.data()[(3 * vertex_order_in_matrix + dimension) * (3 * unfixed_vertex_num + 1)] += collision_stiffness * (h * grad_d * grad_d + g * h_d);

}




void XPBD_IPC::getTVCollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, double* t0, double* t1, double* t2, double stiffness,
	unsigned int* TV, int num, unsigned int obj_No, int* triangle_indices,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* TV_collider, int collider_num)
{
	int triangle_vertex_order_in_tet[4];
	memset(triangle_vertex_order_in_tet, 0xff, 16);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_indices, triangle_vertex_order_in_tet+1,3);

	for (int i = 0; i < num; i += 3) {
		if (TV[i] == obj_No) {
			if (vertexInTet(unfixed_tet_vertex_num, TV[i + 1], tet_unfixed_vertex_indices)) {
				continue;
			}
		}
		second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position[TV[i]][TV[i + 1]].data(),
			t0, t1, t2,
			d_hat_2, triangle_vertex_order_in_tet, stiffness);
	}
	
	if (has_collider) {
		for (int i = 0; i < collider_num; i += 2) {
			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_collider[TV_collider[i]][TV_collider[i + 1]].data(),
				t0, t1, t2,
				d_hat_2, triangle_vertex_order_in_tet, stiffness);
		}
	}

}



//void XPBD_IPC::checkIfRight(int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, unsigned int obj_No, int* triangle_vertex_order_in_tet)
//{
//
//}



void XPBD_IPC ::getEECollisionHessainForTetFromRecord(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_No, unsigned int* edge_vertex_index,
	unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, 
	int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, 
	unsigned int* EE_collider, int num_collider, 
	double* ee_hessian_record, double* ee_grad_record, int* ee_hessian_record_index,
	int* ee_collider_hessian_record_index, double* ee_collider_hessian_record, double* ee_collider_grad_record)
{
	unsigned int* edge_vertex;
	int vertex_order_in_tet[4];
	memset(vertex_order_in_tet, 0xff, 8);

	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex_index, vertex_order_in_tet);

	for (int i = 0; i < num; i += 2) {
		edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
		memset(vertex_order_in_tet + 2, 0xff, 8);
		if (EE[i] == obj_No) {
			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, EE[i + 1])) {
				continue;
			}
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex, vertex_order_in_tet + 2);
		}

		if ((*ee_hessian_record_index) != 0) {
			setTetHessianFromBarrierHessian(Hessian, grad.data(), ee_hessian_record, ee_grad_record, vertex_order_in_tet,
				ee_hessian_record_index + 1, *ee_hessian_record_index);
		}
		ee_hessian_record_index += 5;
		ee_grad_record += 12;
		ee_hessian_record += 144;
	}

	if (has_collider) {
		memset(vertex_order_in_tet + 2, 0xff, 8);
		for (int i = 0; i < num_collider; i += 2) {
			if ((*ee_collider_hessian_record_index) != 0) {
				setTetHessianFromBarrierHessian(Hessian, grad.data(), ee_collider_hessian_record, ee_collider_grad_record, vertex_order_in_tet,
					ee_collider_hessian_record_index + 1, *ee_collider_hessian_record_index);
			}
			ee_collider_hessian_record_index += 3;
			ee_collider_grad_record += 6;
			ee_collider_hessian_record += 36;
		}

	}

}


void XPBD_IPC::getEECollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_No, unsigned int* edge_vertex_index,
	unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2,
	int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, double*ea0, double* ea1, double stiffness,
	unsigned int* EE_collider, int num_collider, double edge_length_0)
{
	unsigned int* edge_vertex;
	int vertex_order_in_tet[4];
	memset(vertex_order_in_tet, 0xff, 8);	

	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex_index, vertex_order_in_tet);

	for (int i = 0; i < num; i += 2) {
		edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
		memset(vertex_order_in_tet + 2, 0xff, 8);
		if (EE[i] == obj_No) {
			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, EE[i + 1])) {
				continue;
			}
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex, vertex_order_in_tet+2);
		}
		second_order_constraint.computeEEBarrierGradientHessian(ea0, ea1,  vertex_position[EE[i]][*edge_vertex].data(), 
			vertex_position[EE[i]][edge_vertex[1]].data(), Hessian, grad,
			vertex_order_in_tet, stiffness, d_hat_2, edge_length_0,rest_edge_length[EE[i]][EE[i+1]]);

	}

	if (has_collider) {
		
		memset(vertex_order_in_tet + 2, 0xff, 8);
		for (int i = 0; i < num_collider; i += 2) {
			std::cout << EE_collider[i] << " " << EE_collider[i + 1] << std::endl;
			//std::cout << "test " << std::endl;
			//std::cout << grad.transpose() << std::endl;
			edge_vertex = collider_edge_vertices[EE_collider[i]] + (EE_collider[i + 1] << 1);
			second_order_constraint.computeEEBarrierGradientHessian(ea0, ea1, vertex_position_collider[EE_collider[i]][*edge_vertex].data(),
				vertex_position_collider[EE_collider[i]][edge_vertex[1]].data(), Hessian, grad,
				vertex_order_in_tet, stiffness, d_hat_2, edge_length_0,rest_edge_length_collider[EE[i]][EE[i+1]]);
			//std::cout << grad.transpose() << std::endl;

		}
	}
}


bool XPBD_IPC::edgeInSameTetDuplicate(int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet,
	unsigned int compare_edge_index)
{
	for (auto i = edge_of_a_tet->begin() + edge_order_in_tet + 1; i < edge_of_a_tet->end(); ++i) {
		if (*i == compare_edge_index) {
			return true;
		}
	}
	return false;
}


bool XPBD_IPC::edgeInSameTetDuplicate(int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet,
	unsigned int compare_edge_index, unsigned int compare_edge_obj)
{
	for (auto i = edge_of_a_tet->begin() + edge_order_in_tet + 2; i < edge_of_a_tet->end(); i+=2) {
		if (*i == compare_edge_obj &&  *(i+1) == compare_edge_index) {
			return true;
		}
	}
	return false;
}

bool XPBD_IPC::vertexInTet(int unfixed_tet_vertex_num, int vertex_No, int* tet_unfixed_vertex_indices)
{
	for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
		if (vertex_No == tet_unfixed_vertex_indices[i]) {
			return true;
		}
	}
	return false;
}

bool XPBD_IPC::vertexInPair(int unfixed_tet_vertex_num, int vertex_No, int* tet_unfixed_vertex_indices, int obj_No)
{
	for (int i = 0; i < unfixed_tet_vertex_num; i+=2) {
		if (vertex_No == tet_unfixed_vertex_indices[i+1] && obj_No == tet_unfixed_vertex_indices[i]) {
			return true;
		}
	}
	return false;
}


void XPBD_IPC::checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, int* element_indices, int obj_No,
	int* triangle_vertex_order_in_system, int size_num)
{
	for (int i = 0; i < size_num; ++i) {
		for (int j = 0; j < unfixed_tet_vertex_num; j+=2) {
			if (tet_unfixed_vertex_indices[j+1] == element_indices[i] && tet_unfixed_vertex_indices[j]==obj_No) {
				triangle_vertex_order_in_system[i] = (j>>1);
				break;
			}
		}
	}

}


void XPBD_IPC::checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, int* element_indices,
	int* triangle_vertex_order_in_system, int size_num)
{
	for (int i = 0; i < size_num; ++i) {
		for (int j = 0; j < unfixed_tet_vertex_num; ++j) {
			if (tet_unfixed_vertex_indices[j] == element_indices[i]) {
				triangle_vertex_order_in_system[i] = j;
				break;
			}
		}
	}

}

void XPBD_IPC::checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, unsigned int* element_indices,
	int* triangle_vertex_order_in_system)
{
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < unfixed_tet_vertex_num; ++j) {
			if (tet_unfixed_vertex_indices[j] == element_indices[i]) {
				triangle_vertex_order_in_system[i] = j;
				break;
			}
		}
	}

}

void XPBD_IPC::checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, unsigned int* element_indices, int obj_No,
	int* triangle_vertex_order_in_system)
{
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < unfixed_tet_vertex_num; j+=2) {
			if (tet_unfixed_vertex_indices[j+1] == element_indices[i] && tet_unfixed_vertex_indices[j]==obj_No) {
				triangle_vertex_order_in_system[i] = (j>>1);
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


void XPBD_IPC::solveEE_collisionBlock(unsigned int obj_No_0, unsigned int primitive_0_index, unsigned int obj_No_1, unsigned int primitive_1_index,
	double stiffness, double dt, double collision_stiffne, std::vector<unsigned int>* triangle_around_0, std::vector<unsigned int>* triangle_around_1,
	std::vector<unsigned int>* edge_around_0, std::vector<unsigned int>* edge_around_1,
	std::vector<unsigned int>* tet_around_0, std::vector<unsigned int>* tet_around_1,
	double d_hat_2, bool edge_0_collider, bool edge_1_collider)
{
	MatrixXd Hessian;
	VectorXd grad;
	int unfixed_pair_vertex_index[8];
	memset(unfixed_pair_vertex_index, 0xff, 32);
	int unfixed_num = 0;

	if (!edge_0_collider) {
		unsigned int* edge_vertex_0 = edge_vertices[obj_No_0] + (primitive_0_index << 1);
		if (!(*is_vertex_fixed[obj_No_0])[*edge_vertex_0]) {
			unfixed_pair_vertex_index[unfixed_num] = obj_No_0;
			unfixed_pair_vertex_index[unfixed_num + 1] = *edge_vertex_0;
			unfixed_num += 2;
		}
		if (!(*is_vertex_fixed[obj_No_0])[edge_vertex_0[1]]) {
			unfixed_pair_vertex_index[unfixed_num] = obj_No_0;
			unfixed_pair_vertex_index[unfixed_num + 1] = edge_vertex_0[1];
			unfixed_num += 2;
		}
	}

	if (!edge_1_collider) {
		unsigned int* edge_vertex_1 = edge_vertices[obj_No_1] + (primitive_1_index << 1);
		if (!(*is_vertex_fixed[obj_No_1])[*edge_vertex_1]) {
			unfixed_pair_vertex_index[unfixed_num] = obj_No_1;
			unfixed_pair_vertex_index[unfixed_num + 1] = *edge_vertex_1;
			unfixed_num += 2;
		}
		if (!(*is_vertex_fixed[obj_No_1])[edge_vertex_1[1]]) {
			unfixed_pair_vertex_index[unfixed_num] = obj_No_1;
			unfixed_pair_vertex_index[unfixed_num + 1] = edge_vertex_1[1];
			unfixed_num += 2;
		}
	}

	if (unfixed_num == 0) {
		return;
	}

	grad.resize(3 * (unfixed_num >> 1));
	grad.setZero();
	Hessian.resize(3 * (unfixed_num >> 1), 3 * (unfixed_num >> 1));
	Hessian.setZero();

	std::vector<unsigned int> around_triangle, around_edge, around_tet;

	getCollisionPairHessian(Hessian, grad, obj_No_0, obj_No_1, collision_stiffne, triangle_around_0, triangle_around_1,
		edge_around_0, edge_around_1,
		tet_around_0, tet_around_1,
		d_hat_2, unfixed_pair_vertex_index, unfixed_num, &around_triangle, &around_edge, &around_tet, stiffness,
		edge_0_collider, edge_1_collider);

	double mass_dt_2;
	int obj_index_, vertex_index_;

	int i;

	for (int j = 0; j < unfixed_num; j += 2) {
		i = (j >> 1);
		vertex_index_ = unfixed_pair_vertex_index[j + 1];
		obj_index_ = unfixed_pair_vertex_index[j];
		mass_dt_2 = mass[obj_index_][vertex_index_] / (dt * dt);
		Hessian(3 * i, 3 * i) += mass_dt_2;
		Hessian(3 * i + 1, 3 * i + 1) += mass_dt_2;
		Hessian(3 * i + 2, 3 * i + 2) += mass_dt_2;
		grad.data()[3 * i] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][0] - sn[obj_index_][vertex_index_][0]);
		grad.data()[3 * i + 1] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][1] - sn[obj_index_][vertex_index_][1]);
		grad.data()[3 * i + 2] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][2] - sn[obj_index_][vertex_index_][2]);
	}

	//HouseholderQR <MatrixXd> linear(Hessian);
	LLT <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(grad);

	double tt = result.norm();


	//if (edge_1_collider == false && primitive_0_index == 3533 && primitive_1_index == 6302) {
	//	//std::cout<<ed
	//	//double e0[3]; double e1[3];
	//	//std::vector<unsigned int> around_triangle_, around_edge_, around_tet_;
	//	//MatrixXd hessen_record = Hessian; hessen_record.setZero();
	//	//VectorXd grad_record = grad; grad_record.setZero();
	//	//getCollisionPairHessianTest(hessen_record, grad_record, obj_No_0, obj_No_1, collision_stiffne, triangle_around_0, triangle_around_1,
	//	//	edge_around_0, edge_around_1,
	//	//	tet_around_0, tet_around_1,
	//	//	d_hat_2, unfixed_pair_vertex_index, unfixed_num, &around_triangle_, &around_edge_, &around_tet_, stiffness,
	//	//	edge_0_collider, edge_1_collider);
	//}

	//if (tt > 1e-2) {
	//	std::cout <<"displacement ee"<< edge_1_collider<<" "<< primitive_0_index << " " << primitive_1_index << " " << tt << std::endl;
	//}
	

	for (int i = 0; i < unfixed_num; i += 2) {
		vertex_index_ = unfixed_pair_vertex_index[i + 1];
		obj_index_ = unfixed_pair_vertex_index[i];
		SUB_(vertex_position[obj_index_][vertex_index_], (result.data() + 3 * (i >> 1)));
	}

	double t = getCollisionTime(&around_triangle, &around_edge, unfixed_pair_vertex_index, unfixed_num, vertex_index_surface.data(),
		vertex_position.data(), record_vertex_position_address.data());

	if (t < 1.0) {
		t *= 0.9;
		double* p_c; double* p_i;
		for (int i = 0; i < unfixed_num; i += 2) {
			p_i = record_vertex_position_address[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data();
			p_c = vertex_position[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data();
			COLLISION_POS(p_c, t, p_i, p_c);
			memcpy(p_i, p_c, 24);
		}
	}
	else {
		for (int i = 0; i < unfixed_num; i += 2) {
			memcpy(record_vertex_position_address[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data(),
				vertex_position[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data(), 24);
		}
	}

}


void XPBD_IPC::solveVT_collisionBlock(unsigned int vertex_obj_no, unsigned int vertex_index, unsigned int triangle_obj_No, unsigned int triangle_index,
	double stiffness, double dt, double collision_stiffne, std::vector<unsigned int>* triangle_around_vertex, std::vector<unsigned int>* triangle_around_triangle,
	std::vector<unsigned int>* edge_around_vertex, std::vector<unsigned int>* edge_around_triangle, 
	std::vector<unsigned int>* tet_around_vertex, std::vector<unsigned int>*tet_around_triangle,
	double d_hat_2, bool vertex_collider, bool triangle_collider)
{
	MatrixXd Hessian;
	VectorXd grad;
	int unfixed_pair_vertex_index[8];
	memset(unfixed_pair_vertex_index, 0xff, 32);
	int unfixed_num = 0;
	int* tri_indices = triangle_indices[triangle_obj_No][triangle_index].data();
	if (!vertex_collider) {
		if (!(*is_vertex_fixed[vertex_obj_no])[vertex_index]) {
			unfixed_pair_vertex_index[unfixed_num] = vertex_obj_no;
			unfixed_pair_vertex_index[unfixed_num + 1] = vertex_index;
			unfixed_num += 2;
		}
	}
	if (!triangle_collider) {
		for (unsigned int i = 0; i < 3; ++i) {
			if (!(*is_vertex_fixed[triangle_obj_No])[tri_indices[i]]) {
				unfixed_pair_vertex_index[unfixed_num] = triangle_obj_No;
				unfixed_pair_vertex_index[unfixed_num + 1] = tri_indices[i];
				unfixed_num += 2;
			}
		}
	}
	if (unfixed_num == 0) {
		return;
	}

	grad.resize(3*(unfixed_num >> 1));
	grad.setZero();
	Hessian.resize(3 * (unfixed_num >> 1), 3 * (unfixed_num >> 1));
	Hessian.setZero();
	
	std::vector<unsigned int> around_triangle, around_edge, around_tet;

	getCollisionPairHessian(Hessian, grad, vertex_obj_no, triangle_obj_No, collision_stiffne, triangle_around_vertex, triangle_around_triangle,
		edge_around_vertex, edge_around_triangle,
		tet_around_vertex, tet_around_triangle,
		d_hat_2,unfixed_pair_vertex_index, unfixed_num,&around_triangle,&around_edge,&around_tet, stiffness,vertex_collider, triangle_collider);


	double mass_dt_2;
	int obj_index_,vertex_index_;

	int i;
	for (int j = 0; j < unfixed_num; j+=2) {
		i =( j >> 1);
		vertex_index_ = unfixed_pair_vertex_index[j+1];
		obj_index_= unfixed_pair_vertex_index[j];
		mass_dt_2 = mass[obj_index_][vertex_index_] / (dt * dt);
		Hessian(3 * i, 3 * i) += mass_dt_2;
		Hessian(3 * i + 1, 3 * i + 1) += mass_dt_2;
		Hessian(3 * i + 2, 3 * i + 2) += mass_dt_2;
		grad.data()[3 * i] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][0] - sn[obj_index_][vertex_index_][0]);
		grad.data()[3 * i + 1] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][1] - sn[obj_index_][vertex_index_][1]);
		grad.data()[3 * i + 2] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][2] - sn[obj_index_][vertex_index_][2]);
	}

	//HouseholderQR <MatrixXd> linear(Hessian);
	LLT <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(grad);
	//std::cout << Hessian << std::endl;
	//std::cout << "++" << std::endl;
	//std::cout << grad.transpose() << std::endl;

	double tt = result.norm();

	//if (tt > 1e-2) {
	//	std::cout << "displacement vt" << triangle_collider <<" "<< vertex_collider << " " << vertex_index << " " << triangle_index << " " << tt << std::endl;
	//}

	for (int i = 0; i < unfixed_num; i+=2) {
		vertex_index_ = unfixed_pair_vertex_index[i + 1];
		obj_index_ = unfixed_pair_vertex_index[i];
		SUB_(vertex_position[obj_index_][vertex_index_], (result.data() + 3 * (i>>1)));
	}

	double t = getCollisionTime(&around_triangle, &around_edge, unfixed_pair_vertex_index, unfixed_num, vertex_index_surface.data(),
		vertex_position.data(), record_vertex_position_address.data());

	if (t < 1.0) {
		t *= 0.9;
		double* p_c; double* p_i;
		for (int i = 0; i < unfixed_num; i += 2) {
			p_i = record_vertex_position_address[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data();
			p_c = vertex_position[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data();
			COLLISION_POS(p_c, t, p_i, p_c);
			memcpy(p_i, p_c, 24);
		}
	}
	else {
		for (int i = 0; i < unfixed_num; i += 2) {
			memcpy(record_vertex_position_address[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i+1]].data(),
				vertex_position[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i+1]].data(), 24);
		}
	}
}




void XPBD_IPC::getCollisionPairHessianTest(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_0, unsigned int obj_1,
	double collision_stiffness, std::vector<unsigned int>* triangle_around_0, std::vector<unsigned int>* triangle_around_1,
	std::vector<unsigned int>* edge_around_0, std::vector<unsigned int>* edge_around_1,
	std::vector<unsigned int>* tet_around_0, std::vector<unsigned int>* tet_around_1,
	double d_hat_2, int* unfixed_pair_vertex_index, int unfixed_num, std::vector<unsigned int>* around_triangle, std::vector<unsigned int>* around_edge,
	std::vector<unsigned int>* around_tet, double arap_stiffness, bool obj_0_collider, bool obj_1_collider)
{
	if (obj_0_collider) {
		for (auto i = triangle_around_1->begin(); i < triangle_around_1->end(); ++i) {
			around_triangle->emplace_back(obj_1);
			around_triangle->emplace_back(*i);
		}
		for (auto i = edge_around_1->begin(); i < edge_around_1->end(); ++i) {
			around_edge->emplace_back(obj_1);
			around_edge->emplace_back(*i);
		}
		for (auto i = tet_around_1->begin(); i < tet_around_1->end(); ++i) {
			around_tet->emplace_back(obj_1);
			around_tet->emplace_back(*i);
		}
	}
	else if (obj_1_collider) {
		for (auto i = triangle_around_0->begin(); i < triangle_around_0->end(); ++i) {
			around_triangle->emplace_back(obj_0);
			around_triangle->emplace_back(*i);
		}
		for (auto i = edge_around_0->begin(); i < edge_around_0->end(); ++i) {
			around_edge->emplace_back(obj_0);
			around_edge->emplace_back(*i);
		}
		for (auto i = tet_around_0->begin(); i < tet_around_0->end(); ++i) {
			around_tet->emplace_back(obj_0);
			around_tet->emplace_back(*i);
		}
	}
	else {
		comparePrimitiveAroundPrimitveTogether(triangle_around_0, triangle_around_1, obj_0, obj_1,
			around_triangle);

		comparePrimitiveAroundPrimitveTogether(edge_around_0, edge_around_1, obj_0, obj_1,
			around_edge);

		comparePrimitiveAroundPrimitveTogether(tet_around_0, tet_around_1, obj_0, obj_1,
			around_tet);
	}
	getCollisionBlockCollisionHessianTest(Hessian, grad, around_triangle, around_edge, collision_stiffness,
		unfixed_pair_vertex_index, unfixed_num, d_hat_2, vertex_index_surface.data(), obj_0, obj_1);

	getCollisionBlockTetHessian(Hessian, grad, around_tet, arap_stiffness, unfixed_pair_vertex_index, unfixed_num);
}


void XPBD_IPC::getCollisionPairHessian(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_0, unsigned int obj_1,
	double collision_stiffness, std::vector<unsigned int>* triangle_around_0, std::vector<unsigned int>* triangle_around_1,
	std::vector<unsigned int>* edge_around_0, std::vector<unsigned int>* edge_around_1, 
	std::vector<unsigned int>* tet_around_0, std::vector<unsigned int>* tet_around_1,
	double d_hat_2, int* unfixed_pair_vertex_index, int unfixed_num, std::vector<unsigned int>* around_triangle, std::vector<unsigned int>* around_edge,
	std::vector<unsigned int>* around_tet, double arap_stiffness, bool obj_0_collider, bool obj_1_collider)
{
	if (obj_0_collider) {
		for (auto i = triangle_around_1->begin(); i < triangle_around_1->end(); ++i) {
			around_triangle->emplace_back(obj_1);
			around_triangle->emplace_back(*i);
		}
		for (auto i = edge_around_1->begin(); i < edge_around_1->end(); ++i) {
			around_edge->emplace_back(obj_1);
			around_edge->emplace_back(*i);
		}
		for (auto i = tet_around_1->begin(); i < tet_around_1->end(); ++i) {
			around_tet->emplace_back(obj_1);
			around_tet->emplace_back(*i);
		}
	}
	else if(obj_1_collider) {
		for (auto i = triangle_around_0->begin(); i < triangle_around_0->end(); ++i) {
			around_triangle->emplace_back(obj_0);
			around_triangle->emplace_back(*i);
		}
		for (auto i = edge_around_0->begin(); i < edge_around_0->end(); ++i) {
			around_edge->emplace_back(obj_0);
			around_edge->emplace_back(*i);
		}
		for (auto i = tet_around_0->begin(); i < tet_around_0->end(); ++i) {
			around_tet->emplace_back(obj_0);
			around_tet->emplace_back(*i);
		}
	}
	else {
		comparePrimitiveAroundPrimitveTogether(triangle_around_0, triangle_around_1, obj_0, obj_1,
			around_triangle);

		comparePrimitiveAroundPrimitveTogether(edge_around_0, edge_around_1, obj_0, obj_1,
			around_edge);

		comparePrimitiveAroundPrimitveTogether(tet_around_0, tet_around_1, obj_0, obj_1,
			around_tet);
	}
	getCollisionBlockCollisionHessian(Hessian, grad, around_triangle, around_edge, collision_stiffness,
		unfixed_pair_vertex_index, unfixed_num, d_hat_2, vertex_index_surface.data(), obj_0, obj_1);

	//FEM::SPDprojection(Hessian);

	getCollisionBlockTetHessian(Hessian, grad, around_tet, arap_stiffness, unfixed_pair_vertex_index, unfixed_num);
}





void XPBD_IPC::comparePrimitiveAroundPrimitveTogether(std::vector<unsigned int>* primitive_around_1, std::vector<unsigned int>* primitive_around_2,
	unsigned int obj_1, unsigned int obj_2, std::vector<unsigned int>* primitive_together)
{
	primitive_together->reserve((primitive_around_1->size() + primitive_around_2->size()) << 1);
	if (obj_1 != obj_2) {
		for (auto i = primitive_around_1->begin(); i < primitive_around_1->end(); ++i) {
			primitive_together->emplace_back(obj_1);
			primitive_together->emplace_back(*i);
		}
		for (auto i = primitive_around_2->begin(); i < primitive_around_2->end(); ++i) {
			primitive_together->emplace_back(obj_2);
			primitive_together->emplace_back(*i);
		}
	}
	else {
		std::set_union(primitive_around_1->begin(), primitive_around_1->end(),
			primitive_around_2->begin(), primitive_around_2->end(), std::back_inserter(*primitive_together));
		primitive_together->resize(primitive_together->size() << 1);
		for (int i = primitive_together->size() - 1; i > -1; i -= 2) {
			primitive_together->data()[i] = primitive_together->data()[i >> 1];
			primitive_together->data()[i - 1] = obj_1;
		}
	}
}


void XPBD_IPC::solveTetBlock(std::array<double, 3>* vertex_position, double stiffness, double dt,
	double* mass,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices, 
	double* volume, unsigned int tet_index, std::array<double, 3>* sn, unsigned int* common_vertex_in_order,
	int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num, std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet, double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	int* vertex_index_on_surface, double* hessian_record, double* grad_record)//
{
	if (unfixed_vertex_num == 0) {
		return;
	}

	MatrixXd Hessian;
	VectorXd grad;
	grad.resize(3 * unfixed_vertex_num);
	Hessian.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);
	Hessian.setZero();
	grad.setZero();


	MatrixXd Hessian_collision_test;
	Hessian_collision_test.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);
	Hessian_collision_test.setZero();
	VectorXd grad_collision_test;
	grad_collision_test.resize(3 * unfixed_vertex_num);
	grad_collision_test.setZero();


	second_order_constraint.solveCD_ARAP_block(Hessian, grad, vertex_position, stiffness, A[tet_index],
		volume[tet_index], tet_vertex_index,
		unfixed_tet_vertex_index, unfixed_vertex_num, hessian_record+16*tet_index);


	for (auto i = neighbor_tet_indices.begin(); i < neighbor_tet_indices.end(); ++i) {
		second_order_constraint.solveHessianForNeighborTet(common_vertex_in_order,
			Hessian, grad, hessian_record + 16 * (*i), grad_record + 12 * (*i));
	}

	getCollisionHessianFromRecord(Hessian, grad, triangle_of_a_tet, edge_of_a_tet, collision_stiffness, obj_No, tet_actual_unfixed_vertex_indices,
		unfixed_vertex_num,
		vertex_index_on_surface, collision.vertex_triangle_pair_num_record_prefix_sum[obj_No], collision.edge_edge_pair_num_record_prefix_sum[obj_No],
		collision.triangle_vertex_collider_num_record_prefix_sum[obj_No], collision.edge_edge_collider_num_record_prefix_sum[obj_No], vertex_position);

	//if (Hessian_collision_test.norm() != 0.0) {
	//	std::cout << "show" << std::endl;
	//	std::cout << grad_collision_test.transpose() << std::endl;
	//	std::cout << Hessian_collision_test << std::endl;
	//}


	//if (perform_collision)
	//{
	//	MatrixXd Hessian_collision;
	//	VectorXd grad_collision;
	//	Hessian_collision.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);
	//	Hessian_collision.setZero();
	//	grad_collision.resize(3 * unfixed_vertex_num);
	//	grad_collision.setZero();
	//	getCollisionHessian(Hessian_collision, grad_collision, triangle_of_a_tet, edge_of_a_tet, collision_stiffness, obj_No, tet_actual_unfixed_vertex_indices,
	//		unfixed_vertex_num,
	//		collision.d_hat_2, vertex_index_on_surface, vertex_position);
	//	//if (Hessian_collision.norm() != 0.0) {
	//	//	std::cout << "compare correct " << std::endl;
	//	//	std::cout << grad_collision.transpose() << std::endl;
	//	//	std::cout << Hessian_collision << std::endl;
	//	//}
	//	//FEM::SPDprojection(Hessian_collision);
	//	Hessian += Hessian_collision;
	//	grad += grad_collision;
	//}


	double mass_dt_2;
	int vertex_index;
	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_actual_unfixed_vertex_indices[i];
		mass_dt_2 = mass[vertex_index] / (dt * dt);
		Hessian(3 * i, 3 * i) += mass_dt_2;
		Hessian(3 * i + 1, 3 * i + 1) += mass_dt_2;
		Hessian(3 * i + 2, 3 * i + 2) += mass_dt_2;
		grad.data()[3 * i] += mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
		grad.data()[3 * i + 1] += mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
		grad.data()[3 * i + 2] += mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);
	}

	//ColPivHouseholderQR<MatrixXd> linear(Hessian);


	LLT <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(grad);


	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_actual_unfixed_vertex_indices[i];
		SUB_(vertex_position[vertex_index], (result.data() + 3 * i));
	}


	

	////unsigned int* common_vertex_in_order_ = common_vertex_in_order;


	////if (result.norm() > 1e-1) {
	////	std::cout << result.norm() << std::endl;
	////	second_order_constraint.solveCD_ARAP_blockTest(Hessian, grad, vertex_position, stiffness, A[tet_index],
	////		volume[tet_index], tet_vertex_index,
	////		unfixed_tet_vertex_index, unfixed_vertex_num);

	////	for (auto i = neighbor_tet_indices.begin(); i < neighbor_tet_indices.end(); ++i) {
	////		std::cout << "i " << i- neighbor_tet_indices.begin() << std::endl;
	////		second_order_constraint.solveCertainHessianForNeighborTetTest(vertex_position, stiffness, A[*i], common_vertex_in_order_,
	////			tet_indices[obj_No][*i].data(),
	////			Hessian, volume[*i], grad);
	////	}
	////}


}


void XPBD_IPC::solveNewtonCD_tetBlock(std::array<double, 3>* vertex_position, double stiffness, double dt,
	double* mass,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices, std::array<int, 4>* indices,
	double* volume, unsigned int tet_index, std::array<double, 3>* sn, unsigned int* common_vertex_in_order,
	int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num, std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet, double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	int* vertex_index_on_surface, std::array<double, 3>* record_ori_pos, double* hessian_record)
{
	if (unfixed_vertex_num == 0) {
		return;
	}

	MatrixXd Hessian;
	VectorXd grad;
	grad.resize(3 * unfixed_vertex_num);
	Hessian.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);

	second_order_constraint.solveCD_ARAP_block(Hessian, grad, vertex_position, stiffness, A[tet_index],
		volume[tet_index],  tet_vertex_index,
		unfixed_tet_vertex_index, unfixed_vertex_num, hessian_record+16* tet_index);

	for (auto i = neighbor_tet_indices.begin(); i < neighbor_tet_indices.end(); ++i) {
		second_order_constraint.solveCertainHessianForNeighborTet(vertex_position, stiffness, A[*i], common_vertex_in_order, indices[*i].data(),
			Hessian, volume[*i], grad);
	}

	//MatrixXd test = Hessian;

	if (perform_collision)
	{
		MatrixXd Hessian_collision;
		Hessian_collision.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);
		Hessian_collision.setZero();
		getCollisionHessian(Hessian_collision, grad, triangle_of_a_tet, edge_of_a_tet, collision_stiffness, obj_No, tet_actual_unfixed_vertex_indices,
			unfixed_vertex_num,
			collision.d_hat_2, vertex_index_on_surface, vertex_position);
		//FEM::SPDprojection(Hessian_collision);
		Hessian += Hessian_collision;
	}



	//if (tet_index == 11607) {
	//	std::cout << Hessian - test << std::endl;
	//}
	double mass_dt_2;
	int vertex_index;
	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_actual_unfixed_vertex_indices[i];
		mass_dt_2 = mass[vertex_index] / (dt * dt);
		Hessian(3 * i, 3 * i) += mass_dt_2;
		Hessian(3 * i + 1, 3 * i + 1) += mass_dt_2;
		Hessian(3 * i + 2, 3 * i + 2) += mass_dt_2;
		grad.data()[3 * i] += mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
		grad.data()[3 * i + 1] += mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
		grad.data()[3 * i + 2] += mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);
	}

	//ColPivHouseholderQR<MatrixXd> linear(Hessian);
	LLT <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(grad);


	//if (tet_index == 11607) {
	//	JacobiSVD<MatrixXd > k(Hessian);
	//	std::cout << "svd "<< k.singularValues().transpose() << std::endl;
	//
	//}


	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_actual_unfixed_vertex_indices[i];
		SUB_(vertex_position[vertex_index], (result.data() + 3 * i));
	}

	if (perform_collision) {

		double t = getCollisionTime(triangle_of_a_tet, edge_of_a_tet, obj_No, tet_actual_unfixed_vertex_indices,
			unfixed_vertex_num, vertex_index_on_surface, vertex_position, record_ori_pos);

		if (t < 1.0) {
			t *= 0.9;
			for (int i = 0; i < unfixed_vertex_num; ++i) {
				vertex_index = tet_actual_unfixed_vertex_indices[i];
				COLLISION_POS(vertex_position[vertex_index], t, record_ori_pos[vertex_index], vertex_position[vertex_index]);
				memcpy(record_ori_pos[vertex_index].data(), vertex_position[vertex_index].data(), 24);
			}
		}
		else {
			for (int i = 0; i < unfixed_vertex_num; ++i) {
				memcpy(record_ori_pos[tet_actual_unfixed_vertex_indices[i]].data(), vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), 24);
			}
		}
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
		initial_vertex_pos = vertex_position_render[i];
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
	double tet_energy;
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
				tet_energy = compute_energy.computeARAPEnergy(vertex_pos[indices[j][0]].data(), vertex_pos[indices[j][1]].data(),
					vertex_pos[indices[j][2]].data(), vertex_pos[indices[j][3]].data(), A[j], volume[j], stiffness);
				energy += tet_energy;				
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
		q_pre = record_collision_free_vertex_position[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			q_end[j][0] = q_pre[j][0] + collision_time * (q_end[j][0] - q_pre[j][0]);
			q_end[j][1] = q_pre[j][1] + collision_time * (q_end[j][1] - q_pre[j][1]);
			q_end[j][2] = q_pre[j][2] + collision_time * (q_end[j][2] - q_pre[j][2]);
		}
	}
}

