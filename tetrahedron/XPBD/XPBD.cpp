#include"XPBD.h"


XPBD::XPBD()
{
	gravity_ = 9.8;
	sub_step_num =18;
	iteration_number =100;

	damping_coe = 0.0;

	perform_collision = true;
	max_iteration_number = 3;
	outer_max_iteration_number = 100;
	XPBD_constraint.epsilon_for_bending = 1e-10;

	velocity_damp = 0.995;
	//energy_converge_ratio = 5e-3;
	
}


void XPBD::initial()
{
	resetExternalForce();
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
	}
}

void XPBD::reset()
{
	resetExternalForce();
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
	}
}

void XPBD::updateItrInfo(int* iteration_num)
{
	iteration_num[LOCAL_GLOBAL] = iteration_number;
	//outer_iteration_number = iteration_num[OUTER];
	//sub_step_num = iteration_num[OUTER];
	sub_time_step = time_step / (double)sub_step_num;
}


void XPBD::setForXPBD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
	Thread* thread, double* tolerance_ratio)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	sub_time_step = time_step / (double)sub_step_num;
	this->thread = thread;
	total_thread_num =thread->thread_num;

	total_obj_num = cloth->size() + tetrahedron->size();
	reorganzieDataOfObjects();
	initialVariable();
	initialClothBending();
	setConstraintIndex();
	//energy_per_thread.resize(thread->thread_num,0.0);
	if (perform_collision) {
		//collision.energy = energy_per_thread.data();
		collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio, XPBD_);
		//collision.setParameter(&lambda_collision,lambda.data()+ constraint_index_start[3], collision_constraint_index_start.data(), damping_coe, sub_time_step);
	}

	setConvergeCondition();

}


void XPBD::setConvergeCondition()
{
	converge_condition_ratio = 1e-3;
	double edge_length = calEdgeLength();
	max_move_standard = converge_condition_ratio * edge_length;
	outer_max_move_standard = 5.0 * converge_condition_ratio * edge_length;
}


void XPBD::initialClothBending()
{
	lbo_weight.resize(cloth->size());
	vertex_lbo.resize(cloth->size());
	rest_mean_curvature_norm.resize(cloth->size());
	//rest_Aq.resize(cloth->size());
	for (unsigned int i = 0; i < cloth->size(); ++i){
		XPBD_constraint.initial_LBO_EdgeCotWeight(cloth->data()[i].mesh_struct, lbo_weight[i], vertex_lbo[i], rest_mean_curvature_norm[i]);
	}

}

void XPBD::setConstraintIndex()
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
			constraint_number += mesh_struct[i]->edge_vertices.size()>>1;
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

void XPBD::initialVariable()
{
	f_ext.resize(total_obj_num);
	velocity.resize(total_obj_num);
	gravity[0] = 0;
	gravity[1] = -gravity_;
	gravity[2] = 0;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		f_ext[i].resize(mesh_struct[i]->vertex_position.size());
		velocity[i].resize(mesh_struct[i]->vertex_position.size());
	}
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
		memset(f_ext[i][0].data(), 0, 24 * f_ext[i].size());
	}
}

void XPBD::reorganzieDataOfObjects()
{
	vertex_position.resize(total_obj_num);
	initial_vertex_position.resize(total_obj_num);
	mesh_struct.resize(total_obj_num);
	vertex_index_begin_per_thread.resize(total_obj_num);
	record_vertex_position.resize(total_obj_num);
	//record_outer_vertex_position.resize(total_obj_num);
	unfixed_vertex.resize(total_obj_num);
	
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_position[i]= cloth->data()[i].mesh_struct.vertex_position.data();
		initial_vertex_position[i]= cloth->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i]= &cloth->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i]= cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		//record_outer_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		unfixed_vertex[i] = &cloth->data()[i].mesh_struct.unfixed_point_index;
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		initial_vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		//record_outer_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		unfixed_vertex[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.unfixed_point_index;
	}
	recordVertexPosition();

	collider_mesh_struct.resize(collider->size());
	for (unsigned int i = 0; i < collider->size(); ++i) {
		collider_mesh_struct[i] = &collider->data()[i].mesh_struct;
	}

}


void XPBD::saveScene()
{
	save_scene.save_scene_XPBD(*time_stamp, *time_indicate_for_simu, mesh_struct, &velocity, collider_mesh_struct);
}

void XPBD::readScene(const char* file_name)
{
	save_scene.read_scene_XPBD(file_name, time_stamp, time_indicate_for_simu, mesh_struct, &velocity, collider_mesh_struct);
	for (unsigned int i = 0; i < mesh_struct.size(); ++i) {
		memcpy(mesh_struct[i]->vertex_for_render[0].data(), mesh_struct[i]->vertex_position[0].data(), 24 * mesh_struct[i]->vertex_position.size());
	}
	for (unsigned int i = 0; i < collider_mesh_struct.size(); ++i) {
		memcpy(collider_mesh_struct[i]->vertex_for_render[0].data(), collider_mesh_struct[i]->vertex_position[0].data(), 24 * collider_mesh_struct[i]->vertex_position.size());
	}
	updateRenderNormal();
	updateNormal();
	updateRenderVertexNormal();
}



double XPBD::calEdgeLength()
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


void XPBD::recordVertexPosition()
{
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		memcpy(record_vertex_position[i][0].data(), cloth->data()[i].mesh_struct.vertex_position[0].data(), 24 * record_vertex_position[i].size());
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		memcpy(record_vertex_position[i + cloth->size()][0].data(), tetrahedron->data()[i].mesh_struct.vertex_position[0].data(), 
			24 * record_vertex_position[i + cloth->size()].size());
	}
}


//void XPBD::recordOuterVertexPosition()
//{
//	for (unsigned int i = 0; i < cloth->size(); ++i) {
//		memcpy(record_outer_vertex_position[i][0].data(), cloth->data()[i].mesh_struct.vertex_position[0].data(), 24 * record_outer_vertex_position[i].size());
//	}
//	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
//		memcpy(record_outer_vertex_position[i + cloth->size()][0].data(), tetrahedron->data()[i].mesh_struct.vertex_position[0].data(),
//			24 * record_outer_vertex_position[i + cloth->size()].size());
//	}
//}


void XPBD::initialCollisionConstriantNum()
{
	lambda_collision.resize(collision.collisionConstraintNumber(collision_constraint_index_start[0].data(), collision_constraint_index_start[1].data(), collision_constraint_index_start[2].data()));
}


void XPBD::PBD_IPCSolve()
{
	thread->assignTask(this, SET_POS_PREDICT);

	for (unsigned int i = 0; i < iteration_number; ++i) {



	}

}


void XPBD::solveByXPBD()
{
	//if (sub_step_num == 1) {
	if (perform_collision) {
		thread->assignTask(this, SET_POS_PREDICT);
		//time_t t1 = clock();
		//for (unsigned int j = 0; j < 10; ++j) {
			collision.collisionCulling();
		//}
	//	time_t t2 = clock() - t1;
	//	std::cout << "t2 " << t2 << std::endl;
	//	unsigned int vt_pair_num = 0;
	//	unsigned int ee_pair_num = 0;
	//	unsigned int vt_c_pair_num=0;
	//	for (unsigned int j = 0; j < total_thread_num; ++j) {
	//		vt_pair_num += collision.spatial_hashing.vertex_triangle_pair[j][0];
	//		ee_pair_num += collision.spatial_hashing.edge_edge_pair[j][0];
	//		vt_c_pair_num += collision.spatial_hashing.vertex_obj_triangle_collider_pair[j][0];
	//		//spatial_hashing.vertex_triangle_pair[j][0] = 0;
	//		//spatial_hashing.edge_edge_pair[j][0] = 0;
	////		spatial_hashing.vertex_obj_triangle_collider_pair[j][0] = 0;
	//	}
	//	std::cout << "vt " << vt_pair_num << " ee " << ee_pair_num << " vt_c " << vt_c_pair_num << std::endl;
	}
	//}
	iteration_number = 0;
	//time_t t3 = clock();
	for (unsigned int sub_step = 0; sub_step < sub_step_num; ++sub_step) {
		memset(lambda.data(), 0, 8 * lambda.size());
		inner_iteration_number = 0;
		if (sub_step_num >1) {
			thread->assignTask(this, SET_POS_PREDICT_SUB_TIME_STEP);
			if (control_parameter[START_TEST]) {
				//if (!collider->empty()) {
				//	move_model->moveSphere(*time_indicate_for_simu, collider->data()[0].mesh_struct.vertex_for_render, collider->data()[0].mesh_struct.vertex_position, sub_step_num);
				//}
				//move_model->moveSkirt(*time_indicate_for_simu, mesh_struct, false, sub_step_num);
				//rorate band capsule
				if (!collider->empty()) {
					move_model->sceneRotateCapsule(*time_indicate_for_simu, collider->data()[0].mesh_struct.vertex_for_render, collider->data()[0].mesh_struct.vertex_position, mesh_struct[0], false, sub_step_num);
				}
			}
		}	
		while (!convergeCondition(inner_iteration_number)) {
			//recordVertexPosition();
			if (perform_collision) {
				updateNormal();
			}
			solveConstraint((inner_iteration_number==0 ) && sub_step% *sub_step_per_detection ==0);//sub_step % prediction_sub_step_size//|| inner_iteration_number== (max_iteration_number/2+1)
			inner_iteration_number++;
		}
		iteration_number += inner_iteration_number;
		thread->assignTask(this, XPBD_VELOCITY);
		updatePosition();
		updateRenderNormal();
	}

	//second_order_constraint.test(*(mesh_struct[0]), lbo_weight[0], vertex_lbo[0], rest_mean_curvature_norm[0]);

	//time_t t4 = clock() - t3;
	//std::cout << "t4 " << t4 << std::endl;
	updateRenderVertexNormal();
}


void XPBD::solveByPBD()
{
	thread->assignTask(this, SET_POS_PREDICT);

	for (unsigned int sub_step = 0; sub_step < sub_step_num; ++sub_step) {
		iteration_number = 0;
		if (sub_step_num > 1) {
			thread->assignTask(this, SET_POS_PREDICT_SUB_TIME_STEP);
		}
		if (perform_collision) {
			collision.collisionCulling();
			collision.getCollisionPair();
		}
		while (iteration_number < max_iteration_number) {
			recordVertexPosition();
			if (perform_collision) {
				updateNormal();
			}
			solveConstraint(iteration_number==0);
			iteration_number++;
		}
		thread->assignTask(this, XPBD_VELOCITY);
		updatePosition();
		updateRenderNormal();
	}
	updateRenderVertexNormal();
}

void XPBD::PBDsolve()
{
	if (use_PBD) {
		solveByPBD();
	}
	else {
		solveByXPBD();
	}
}


//bool XPBD::outerConvergeCondition(unsigned int iteration_num)
//{
//	if (iteration_num < 1) {
//		return false;
//	}
//	if (iteration_num > outer_max_iteration_number) {
//		return true;
//	}
//	unsigned int* unfixed_vertex_index;
//	std::array<double, 3>* current_pos;
//	std::array<double, 3>* previous_pos;
//
//	for (unsigned int i = 0; i < total_obj_num; ++i) {
//		unfixed_vertex_index = unfixed_vertex[i]->data();
//		previous_pos = record_outer_vertex_position[i].data();
//		current_pos = vertex_position[i];
//
//		for (unsigned int j = 0; j < unfixed_vertex[i]->size(); ++j) {
//			for (unsigned int k = 0; k < 3; ++k) {
//				if (abs(previous_pos[unfixed_vertex_index[j]][k] - current_pos[unfixed_vertex_index[j]][k]) > max_move_standard) {
//					return false;
//				}
//			}
//		}
//	}
//
//	return true;
//}

bool XPBD::convergeCondition(unsigned int iteration_num)
{
	if (iteration_num < max_iteration_number) {
		return false;
	}

	return true;
	
	//if (iteration_num > max_iteration_number-1) {
	//	return true;
	//}
	//unsigned int* unfixed_vertex_index;
	//std::array<double, 3>* current_pos;
	//std::array<double, 3>* previous_pos;
	//if (abs(energy - previous_energy) / previous_energy < energy_converge_ratio) {
	//	return true;
	//}
	//for (unsigned int i = 0; i < total_obj_num; ++i) {
	//	unfixed_vertex_index = unfixed_vertex[i]->data();
	//	previous_pos = record_vertex_position[i].data();
	//	current_pos = vertex_position[i];
	//	
	//	for (unsigned int j = 0; j < unfixed_vertex[i]->size(); ++j) {
	//		for(unsigned int k=0;k<3;++k){
	//			if (abs(previous_pos[unfixed_vertex_index[j]][k] - current_pos[unfixed_vertex_index[j]][k])> max_move_standard) {
	//				return false;
	//			}
	//		}
	//	}
	//}

	//return false;

}

void XPBD::updatePosition()
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

void XPBD::solveConstraint(bool need_detection)
{
	//previous_energy = energy;
	//energy = 0.0;
	solveBendingConstraint();
	solveEdgeLengthConstraint();
	solveTetStrainConstraint();
	if (need_detection) {
		collision.XPBDsolveCollisionConstraint();
	}
	else {
		collision.re_XPBDsolveCollisionConstraint();
	}
}


void XPBD::updateNormal()
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


void XPBD::updateRenderVertexNormal()
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


void XPBD::updateRenderNormal()
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


void XPBD::solveEdgeLengthConstraint()
{
	unsigned int size;
	MeshStruct* mesh_struct_;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	double stiffness;
	unsigned int* edge_vertex_index;
	double* mass_inv;
	//double delta_t = sub_time_step;
	double* lambda_ = lambda.data()+ constraint_index_start[1];
	//double damping_stiffness = damping_coe;
	//for (unsigned int i = 0; i < mesh_struct[0]->vertex_position.size(); ++i) {
	//	std::cout << vertex_position[0][i].data()[0]<<" "<< vertex_position[0][i][0] << std::endl;
	//}
	double damp_stiffness;

	//double energy_;

	for (unsigned int i = 0; i < cloth->size(); ++i) {
		mesh_struct_ = mesh_struct[i];
		size = mesh_struct_->edge_length.size();
		vertex_pos = vertex_position[i];
		initial_vertex_pos = initial_vertex_position[i];

		//std::cout << vertex_position[i] << " " << vertex_pos << std::endl;
		damp_stiffness = cloth->data()[i].damp_length_stiffness;
		stiffness = cloth->data()[i].length_stiffness;
		edge_vertex_index = mesh_struct_->edge_vertices.data();
		mass_inv = mesh_struct_->mass_inv.data();
		//for (unsigned int j = 0; j < size; ++j) {			
		for (auto k = mesh_struct_->unconnected_edge_index.begin(); k < mesh_struct_->unconnected_edge_index.end(); ++k) {
			for (auto j = k->begin(); j < k->end(); ++j) {
				//std::cout << edge_vertex_index[j << 1] << " " << edge_vertex_index[(j << 1) + 1] << std::endl;
				//std::cout << *lambda_ << " " << vertex_position[i][edge_vertex_index[j << 1]].data()[0] << " " << vertex_position[i][edge_vertex_index[(j << 1) + 1]].data()[0] << std::endl;
				XPBD_constraint.solveEdgeLengthConstraint(vertex_pos[edge_vertex_index[(*j) << 1]].data(),
					vertex_pos[edge_vertex_index[((*j) << 1) + 1]].data(), mesh_struct_->edge_length[*j], stiffness, sub_time_step, mass_inv[edge_vertex_index[(*j) << 1]],
					mass_inv[edge_vertex_index[((*j) << 1) + 1]], *lambda_, damp_stiffness, initial_vertex_pos[edge_vertex_index[(*j) << 1]].data(),
					initial_vertex_pos[edge_vertex_index[((*j) << 1) + 1]].data());
				//std::cout << *lambda_ << " " << vertex_position[i][edge_vertex_index[j << 1]].data()[0] << " " << vertex_position[i][edge_vertex_index[(j << 1) + 1]].data()[0] << std::endl;
				//std::cout << vertex_pos[edge_vertex_index[j << 1]].data()[0]<<" "<< vertex_pos[edge_vertex_index[(j << 1) + 1]].data()[0]<<" "<< edge_vertex_index[j << 1] << " " << edge_vertex_index[(j << 1) + 1] << std::endl;		
				lambda_++;

				//energy += energy_;
			}
		}
	}
}



void XPBD::solveTetStrainConstraint()
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
	double* lambda_ = lambda.data() + constraint_index_start[2];
	double* sigma_limit;
	double youngs_modulus, poisson_ratio;
	std::array<double, 3>* original_vertex_pos;
	double damp_stiffness;

	double iteration_num_inverse = 1.0 /(double)(max_iteration_number + 1);
	double energy_;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		mesh_struct_ = mesh_struct[i+cloth->size()];
		size =tetrahedron->data()[i].mesh_struct.indices.size();
		indices =tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i+cloth->size()];
		initial_vertex_pos = initial_vertex_position[i + cloth->size()];
		stiffness = tetrahedron->data()[i].ARAP_stiffness;
		A = tetrahedron->data()[i].mesh_struct.A.data();
		sigma_limit = (*tetrahedron)[i].sigma_limit;
		mass_inv = mesh_struct_->mass_inv.data();
		youngs_modulus = tetrahedron->data()[i].youngs_modulus;
		poisson_ratio = tetrahedron->data()[i].poisson_ratio;
		damp_stiffness = tetrahedron->data()[i].damp_ARAP_stiffness;
		if (use_PBD) {
			for (unsigned int j = 0; j < size; ++j) {
				XPBD_constraint.PBDsolveARAPConstraint(vertex_pos, initial_vertex_pos, stiffness, sub_time_step, A[j], indices[j].data(), mass_inv,
					volume[j], iteration_num_inverse);
				lambda_++;
			}
		}
		else {
			for (unsigned int j = 0; j < size; ++j) {
				XPBD_constraint.solveARAPConstraint(vertex_pos, initial_vertex_pos, stiffness, sub_time_step, A[j], indices[j].data(), mass_inv,
					*lambda_, damp_stiffness, sigma_limit[0], sigma_limit[1], volume[j], energy_);

				lambda_++;
				//energy += energy_;
			}
		}
		//original_vertex_pos = tetrahedron->data()[i].ori_vertices.data();
		//for (unsigned int j = 0; j < size; ++j) {
		//	XPBD_constraint.solveARAPConstraint2(original_vertex_pos, vertex_pos, initial_vertex_pos, stiffness, sub_time_step, A[j], indices[j].data(), mass_inv,
		//		*lambda_, damping_coe, sigma_limit[0], sigma_limit[1], volume[j]);
		//	lambda_++;
		//}

		//for (unsigned int j = 0; j < size; ++j) {
		//	XPBD_constraint.solveTetStrainConstraint(vertex_pos, initial_vertex_pos, stiffness, sub_time_step, A[j], indices[j].data(), mass_inv,
		//		*lambda_, damping_coe, volume[j], youngs_modulus, poisson_ratio);
		//	lambda_ ++;
		//}
	}
}


void XPBD::solveBendingConstraint()
{
	unsigned int size;
	double* mass_inv;
	double* rest_mean_curvature_norm_;
	double* lbo_weight_;
	VectorXd* vertex_lbo_;
	MeshStruct* mesh_struct_;
	std::array<double, 3>*vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	//Vector3d* rest_Aq_;
	double stiffness;
	double damp_stiffness;
	double* lambda_ = lambda.data();
	//double delta_t = sub_time_step;
	//double energy_;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		stiffness = cloth->data()[i].bend_stiffness;
		if (stiffness > 0.0) {
			damp_stiffness = cloth->data()[i].damp_bend_stiffness;
			mesh_struct_ = mesh_struct[i];
			mass_inv = mesh_struct_->mass_inv.data();
			size = mesh_struct_->vertex_position.size();
			rest_mean_curvature_norm_ = rest_mean_curvature_norm[i].data();
			//rest_Aq_ = rest_Aq[i].data();
			lbo_weight_ = lbo_weight[i].data();
			vertex_lbo_ = vertex_lbo[i].data();
			vertex_pos = vertex_position[i];
			initial_vertex_pos = initial_vertex_position[i];

			//for (unsigned int j = 0; j < size; ++j) {
			for(auto k= mesh_struct_->unconnected_vertex_index.begin();k<mesh_struct_->unconnected_vertex_index.end(); ++k) {
				for (auto j = k->begin(); j < k->end(); ++j) {
					XPBD_constraint.solveBendingConstraint(vertex_pos[*j].data(), mass_inv[*j], vertex_pos, mesh_struct_->vertices[*j].neighbor_vertex.data(), mesh_struct_->vertices[*j].neighbor_vertex.size(),
						rest_mean_curvature_norm_[*j], lbo_weight_[*j], vertex_lbo_[*j], stiffness, sub_time_step, mass_inv, *lambda_);//, damp_stiffness,				initial_vertex_pos[*j].data(), initial_vertex_pos
					lambda_++;
					//energy += energy_;
				}			
			}
		}
		else {
			lambda_+= mesh_struct[i]->vertex_position.size();
		}
	}
}

//XPBD_VELOCITY
void XPBD::computeVelocity(int thread_No)
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


//SET_POS_PREDICT
void XPBD::setPosPredict(int thread_No)
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


//SET_POS_PREDICT_SUB_TIME_STEP_FOR_CULLING
//SET_POS_PREDICT_SUB_TIME_STEP
void XPBD::setPosPredictSubTimeStep(int thread_No, bool predictLargerStep)
{
	std::array<double, 3>*vertex_pos;
	std::array<double, 3>*vertex_pos_initial;
	unsigned int vertex_end=0;
	double* mass_inv;
	std::array<double, 3>* f_ext_;
	std::array<double, 3>* velocity_;
	double gravity__[3];
	memcpy(gravity__, gravity, 24);
	double delta_t = sub_time_step;
	double delta_t_2 = delta_t * delta_t;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		vertex_pos_initial = initial_vertex_position[i];
		vertex_end = vertex_index_begin_per_thread[i][thread_No + 1];
		mass_inv = mesh_struct[i]->mass_inv.data();
		f_ext_ = f_ext[i].data();
		velocity_ = velocity[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < vertex_end; ++j) {
			if (mass_inv[j] != 0) {
				for (unsigned int k = 0; k < 3; ++k) {
					vertex_pos[j][k] = vertex_pos_initial[j][k] + delta_t * velocity_[j][k] + delta_t_2 * (mass_inv[j] * f_ext_[j][k] + gravity__[k]);
				}
			}	
		}
	}
}






void XPBD::resetExternalForce()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(f_ext[i][0].data(), 0, 24 * f_ext[i].size());
	}
}

void XPBD::initialDHatTolerance(double ave_edge_length)
{
	if (perform_collision) {
		collision.initialDHatTolerance(ave_edge_length);
	}
}

void XPBD::updateTetrahedronAnchorVertices()
{
	double* mass_inv;
	int* anchor_vertex;
	unsigned int anchor_vertex_size;
	for (unsigned int i = cloth->size(); i < cloth->size()+tetrahedron->size(); ++i) {
		mass_inv = mesh_struct[i]->mass_inv.data();
		anchor_vertex_size = mesh_struct[i]->anchor_vertex.size();
		anchor_vertex = mesh_struct[i]->anchor_vertex.data();
		mesh_struct[i]->resetMassInv();
	}
}

void XPBD::addExternalForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int obj_No)
{
	if (!coe.empty()) {
		for (unsigned int i = 0; i < coe.size(); ++i) {
			for (unsigned int j = 0; j < 3; ++j) {
				f_ext[obj_No][neighbor_vertex[i]][j] += coe[i] * neighbor_vertex_force_direction[j];
			}
		}
	}
}