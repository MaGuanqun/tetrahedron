#include"XPBD_large_system.h"
#include"basic/write_txt.h"
#include"XPBD/FEM_relate.h"

#include"testHessian.h"

SecondOrderLargeSystem::SecondOrderLargeSystem()
{
	gravity_ = 9.8;

	iteration_number = 10;

	time_step = 1.0 / 30.0;
	perform_collision = false;
	time_step_square = time_step * time_step;
	conv_rate = time_step * 1e-5;

	max_itr_num = 100;
	velocity_damp = 0.99;
	beta = 0.25;
	gamma = 0.5;

	energy_conv_rate = 1e-3;
	//TEST_HESSIAN::testARAPHessianMulti();
}



void SecondOrderLargeSystem::setForNewtonMethod(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
	Thread* thread, double* tolerance_ratio)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	this->thread = thread;
	total_obj_num = cloth->size() + tetrahedron->size();
	total_thread_num = thread->thread_num;
	reorganzieDataOfObjects();
	energy_per_thread.resize(total_thread_num);
	off_diagonal_hessian_nnz_index_begin_per_thread.resize(total_thread_num + 1);
	unfixed_gradC_hessian_index_begin_per_thread.resize(total_thread_num + 1);
	fixed_gradC_hessian_index_begin_per_thread.resize(total_thread_num + 1);
	unfixed_constraint_start_per_thread_in_system.resize(total_thread_num + 1);
	fixed_vertex_constraint_start_per_thread_in_system.resize(total_thread_num + 1);
	hessian_coeff_diagonal.resize(total_thread_num);
	initialHessianNnz();
	computeGravity();
	initialLambda();
	setHessian();

	if (perform_collision) {
		collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio, NEWTON_);
	}
}





void SecondOrderLargeSystem::initialHessianNnz()
{
	unsigned int edge_num = 0;
	off_diagonal_hessian_nnz_index_begin_per_thread[0] = 9 * vertex_begin_per_obj[total_obj_num];

	for (unsigned int i = 0; i < total_thread_num; ++i) {
		edge_num = 0;
		for (unsigned int j = 0; j < cloth->size(); ++j) {
			edge_num += edge_index_begin_per_thread_for_mass_spring[j][i + 1] - edge_index_begin_per_thread_for_mass_spring[j][i];
		}
		off_diagonal_hessian_nnz_index_begin_per_thread[i + 1] = off_diagonal_hessian_nnz_index_begin_per_thread[i] + 18 * edge_num;
	}

	unfixed_gradC_hessian_index_begin_per_thread[0] = off_diagonal_hessian_nnz_index_begin_per_thread[total_thread_num];
	unfixed_constraint_start_per_thread_in_system[0] = 3 * vertex_begin_per_obj[total_obj_num];
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		edge_num = 0;
		for (unsigned int j = 0; j < cloth->size(); ++j) {
			edge_num += edge_index_begin_per_thread_for_mass_spring[j][i + 1] - edge_index_begin_per_thread_for_mass_spring[j][i];
		}
		unfixed_gradC_hessian_index_begin_per_thread[i + 1] = unfixed_gradC_hessian_index_begin_per_thread[i] + 13 * edge_num;
		unfixed_constraint_start_per_thread_in_system[i + 1] = unfixed_constraint_start_per_thread_in_system[i] + edge_num;
	}

	fixed_gradC_hessian_index_begin_per_thread[0] = unfixed_gradC_hessian_index_begin_per_thread[total_thread_num];
	fixed_vertex_constraint_start_per_thread_in_system[0] = unfixed_constraint_start_per_thread_in_system[total_thread_num];
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		edge_num = 0;
		for (unsigned int j = 0; j < cloth->size(); ++j) {
			edge_num += only_one_vertex_fixed_edge_index_begin_per_thread[j][i + 1] - only_one_vertex_fixed_edge_index_begin_per_thread[j][i];
		}
		fixed_gradC_hessian_index_begin_per_thread[i + 1] = fixed_gradC_hessian_index_begin_per_thread[i] + 7 * edge_num;
		fixed_vertex_constraint_start_per_thread_in_system[i + 1] = fixed_vertex_constraint_start_per_thread_in_system[i] + edge_num;	
	}

	ARAP_constrant_in_system.resize(tetrahedron->size()+1);

	ARAP_constrant_in_system[0] = fixed_vertex_constraint_start_per_thread_in_system[total_thread_num];


	std::cout << "ARAP_constrant_in_system[0] " << ARAP_constrant_in_system[0] << std::endl;

	
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		hessian_coeff_diagonal[i].resize(9 * vertex_begin_per_obj[total_obj_num]);
	}

	unsigned int total_tet_num=0;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		total_tet_num += tetrahedron->data()[i].mesh_struct.indices.size();
	}

	tet_hessian_index[0] = fixed_gradC_hessian_index_begin_per_thread[total_thread_num];
	tet_hessian_index[1] = fixed_gradC_hessian_index_begin_per_thread[total_thread_num] + 133 * total_tet_num;

	hessian_nnz.reserve(tet_hessian_index[1]);
	hessian_nnz.resize(tet_hessian_index[0]);






}


void SecondOrderLargeSystem::setSizeOfSys()
{
	sys_total_size = ARAP_constrant_in_system[tetrahedron->size()];
	sys_total_size_index = sys_total_size - 1;
	Hessian.resize(sys_total_size, sys_total_size);
	Sn.resize(3 * vertex_begin_per_obj[total_obj_num]);
	f_ext.resize(3 * vertex_begin_per_obj[total_obj_num]);
	velocity.resize(3 * vertex_begin_per_obj[total_obj_num]);
	velocity.setZero();

	b_thread.resize(total_thread_num);
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		b_thread[i].resize(sys_total_size);
	}

	b_test.resize(sys_total_size);
	Matrix_test.resize(sys_total_size, sys_total_size);

	residual.resize(total_obj_num);
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		residual[i].resize(total_vertex_num[i],Vector3d::Zero());
	}
}

void SecondOrderLargeSystem::computeGravity()
{
	gravity.resize(3 * vertex_begin_per_obj[total_obj_num]);
	//double gravity_accerlation[3] = { 0,0.0,gravity_};
	//double gravity_accerlation[3] = { gravity_, 0,0.0};
	double gravity_accerlation[3] = { 0.0, -gravity_, 0.0 };
	double* mass_;
	unsigned int vertex_index_start;

	unsigned int* unfixed_index_to_normal_index;

	//set cloth
	unsigned int size;
	for (unsigned int j = 0; j < total_obj_num; ++j) {
		mass_ = mass[j];
		vertex_index_start = vertex_begin_per_obj[j];
		size = vertex_begin_per_obj[j + 1] - vertex_index_start;
		unfixed_index_to_normal_index = unfixed_vertex[j]->data();
		for (unsigned int k = 0; k < size; ++k) {
			gravity[3 * (vertex_index_start + k)] = gravity_accerlation[0] * mass_[unfixed_index_to_normal_index[k]];
			gravity[3 * (vertex_index_start + k) + 1] = gravity_accerlation[1] * mass_[unfixed_index_to_normal_index[k]];
			gravity[3 * (vertex_index_start + k) + 2] = gravity_accerlation[2] * mass_[unfixed_index_to_normal_index[k]];
		}
	}
	//Sn = gravity;
}

void SecondOrderLargeSystem::initial()
{
	velocity.setZero();

}


//void SecondOrderLargeSystem::recordEdgeHessian()
//{
//	hessian_coeff_off_diagonal.resize(total_thread_num);
//	hessian_coeff_diagonal.resize(total_thread_num);
//	unsigned int edge_num = 0;
//	unsigned int vertex_num = 0;
//	for (unsigned int i = 0; i < total_thread_num; ++i) {
//		edge_num = 0;
//		for (unsigned int j = 0; j < cloth->size(); ++j) {
//			edge_num += cloth->data()[j].mesh_struct.edge_index_begin_per_thread[i + 1] - cloth->data()[j].mesh_struct.edge_index_begin_per_thread[i];
//		}
//		for (unsigned int j = 0; j < tetrahedron->size(); ++j) {
//			edge_num += tetrahedron->data()[j].mesh_struct.tet_edge_index_begin_per_thread[i + 1] - tetrahedron->data()[j].mesh_struct.tet_edge_index_begin_per_thread[i];
//		}
//		hessian_coeff_off_diagonal[i].resize(edge_num);
//		hessian_coeff_diagonal[i].resize(vertex_begin_per_obj[total_obj_num],0.0);
//	}
//}

////INITIAL_HESSIAN_COEFF
//void SecondOrderLargeSystem::initial_hessian_coeff(int thread_No)
//{
//	unsigned int vertex_begin;
//	unsigned int vertex_end;
//	unsigned int* edge_vertex;
//	std::array<double, 3>* vertex_pos;
//
//	HessianCoeff* hessian_coeff_off_diagonal_ = hessian_coeff_off_diagonal[thread_No].data();
//
//	unsigned int j = 0;
//	for (unsigned int obj_No = 0; obj_No < total_obj_num; ++obj_No) {
//		vertex_pos = vertex_position[obj_No];
//		edge_vertex = edge_vertices_mass_spring[obj_No];
//		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
//		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
//		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
//			hessian_coeff_off_diagonal_[j].row = edge_vertex[i] + vertex_begin_per_obj[obj_No];
//			hessian_coeff_off_diagonal_[j].col = edge_vertex[i+1] + vertex_begin_per_obj[obj_No];
//			j ++;
//		}
//	}
//}






void SecondOrderLargeSystem::setHessian()
{
	initialCoeffDiagonal();
	thread->assignTask(this, SET_MASS_SPRING);
	setARAP();
	thread->assignTask(this, SET_HESSIAN_DIAGONAL);
	setSizeOfSys();
	Hessian.setFromTriplets(hessian_nnz.begin(), hessian_nnz.end());

	Hessian_coeff_address.reserve(hessian_nnz.size());
	Hessian_coeff_address.resize(fixed_gradC_hessian_index_begin_per_thread[total_thread_num]);


	//Hessian_coeff_address.resize(Hessian.nonZeros());
	//if (Hessian_coeff_address.size() != hessian_nnz.size()) {
	//	std::cout << "error, nonzeros of hessian does not equal to hessian_nnz size " << std::endl;
	//}


	thread->assignTask(this, GET_COEFF_ADDRESS);
	getARAPCoeffAddress();
	//std::cout << "k2" << std::endl;
	global_llt.analyzePattern(Hessian);
	//std::cout << "k3" << std::endl;
	//setK();
	//std::cout << Hessian.block(111, 111, Hessian.rows() - 111, Hessian.rows() - 111) << std::endl;
	//std::cout << Hessian.rows() << " " << Hessian.cols() << std::endl;
}

//void SecondOrderLargeSystem::setK()
//{
//	thread->assignTask(this, UPDATE_HESSIAN_FIXED_STRUCTURE);
//	thread->assignTask(this, UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE_INITIAL_STIFFNESS);
//
//	memcpy(store_ori_value.data(), Hessian.valuePtr(), 8 * Hessian.nonZeros());
//
//	//std::cout << ori_stiffness_matrix_record << std::endl;
//}



bool SecondOrderLargeSystem::edgeLengthStiffnessHasChanged()
{
	for (unsigned int i = 0; i < previous_frame_edge_length_stiffness.size(); ++i) {
		if (previous_frame_edge_length_stiffness[i] != *edge_length_stiffness[i]) {
			return false;
		}
	}
	return true;
}

void SecondOrderLargeSystem::updateHessianFixedStructure()
{
	memset(Hessian.valuePtr(), 0, 8*Hessian.nonZeros());

	for (unsigned int i = 0; i < total_thread_num; ++i) {
		b_thread[i].setZero();
	}
	thread->assignTask(this, UPDATE_HESSIAN_FIXED_STRUCTURE);

	updateARAPHessianFixedStructure();

	//std::cout << Hessian.block(111, 111, Hessian.rows() - 111, Hessian.rows() - 111) << std::endl;
	//std::cout << Hessian.rows() << " " << Hessian.cols() << std::endl;

	thread->assignTask(this, UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE);



	thread->assignTask(this, UPDATE_INTERNAL_FORCE);
	thread->assignTask(this, SUM_B);

	//checkIfSampleRight();

	//std::cout << (b_test - b_thread[0])[135] << std::endl;
	//std::cout << temp_record_0 - temp_record_1 << std::endl;
	//std::cout << "=========" << std::endl;

	//

	//thread->assignTask(this, UPDATE_DAMP);
	//thread->assignTask(this, UPDATE_ANCHOR_POINT_HESSIAN);
}



void SecondOrderLargeSystem::addExternalForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int obj_No)
{
	if (!coe.empty()) {
		unsigned int vertex_start = vertex_begin_per_obj[obj_No];
		unsigned int* real_index_to_unfixed_index_ = real_index_to_unfixed_index[obj_No];
		unsigned int size = total_vertex_num[obj_No];
		unsigned int index;
		for (unsigned int i = 0; i < coe.size(); ++i) {
			if (real_index_to_unfixed_index_[neighbor_vertex[i]] < size) {
				index = 3 * (vertex_start + real_index_to_unfixed_index_[neighbor_vertex[i]]);
				for (unsigned int j = 0; j < 3; ++j) {
					f_ext[index + j] += coe[i] * neighbor_vertex_force_direction[j];
				}
			}
		}
	}
}






void SecondOrderLargeSystem::resetExternalForce()
{
	f_ext = gravity;
}


void  SecondOrderLargeSystem::solveNewtonMethod()
{
	solveNewtonMethod_();
}


void SecondOrderLargeSystem::resetLambda()
{
	for (unsigned int i = 0; i < lambda_unfixed.size(); ++i) {
		memset(lambda_unfixed[i].data(), 0, 8 * lambda_unfixed[i].size());
	}
	for (unsigned int i = 0; i < lambda_fixed_one_vertex.size(); ++i) {
		memset(lambda_fixed_one_vertex[i].data(), 0, 8 * lambda_fixed_one_vertex[i].size());
	}

	for (unsigned int i = 0; i < tet_lambda.size(); ++i) {
		memset(tet_lambda[i].data(), 0, 8 * tet_lambda[i].size());
	}

	for (unsigned int i = 0; i < lambda_test.size(); ++i) {
		memset(lambda_test[i].data(), 0, 8 * lambda_test[i].size());
		//std::fill(lambda_test[i].begin(), lambda_test[i].end(), 1.0);
	}
}


void SecondOrderLargeSystem::solveNewtonMethod_()
{
	//storeInitialPosition();
	thread->assignTask(this, SET_S_N);
	updatePositionFromSn();


	resetLambda();

	iteration_number = 0;

	//std::cout << "====" << std::endl;
	computeEnergy();
	//std::cout << "energy " << total_energy << std::endl;
	computeResidual();
	//std::cout <<"== "<< total_residual << std::endl;
	previous_residual = total_residual;
	previous_energy = total_energy;
	while (convergenceCondition())
	{
		updateTest();
		//updateHessianFixedStructure();

		//global_llt.factorize(Hessian);
		//delta_x = global_llt.solve(b_thread[0]);

		displacement_coe = 1.0;
		thread->assignTask(this, UPDATE_POSITION_NEWTON);
		updateARAPLambda();
		computeEnergy();
		//std::cout << "energy0 " << total_energy << std::endl;
		computeResidual();

		if (iteration_number != 0) {
			if (abs(total_energy) >abs(previous_energy)) {
				displacement_coe *= 0.5;
				change_direction = false;
				//std::cout << "activate " << std::endl;
				while (abs(total_energy) > abs(previous_energy)) 
				{
					thread->assignTask(this, UPDATE_POSITION_NEWTON_FROM_ORI);
					updateARAPLambdaFromOri();
					computeEnergy();
					computeResidual();
					if (abs(displacement_coe) < 1e-4) {
						//std::cout << "too many times " << std::endl;
						//std::cout <<"determinent "<< global_llt.determinant()<<" "<< global_llt.logAbsDeterminant() << std::endl;
						//std::cout << (Matrix_test * delta_x - b_test).norm()<<" "<<b_test.norm() << std::endl;
						//total_residual
						//std::cout << "displacement_coe too small " << displacement_coe << std::endl;
						//if (!change_direction) {
						//	//displacement_coe = 1.0;
						//	change_direction = true;
						//}
						//else {
							break;
						//}
					}
					displacement_coe *= 0.5;
				}
			}
		}
		
		//if (total_residual > previous_residual) {
		//	std::cout << total_residual << " " << previous_residual << std::endl;
		//}

		//std::cout << total_energy << std::endl;

		iteration_number++;
		previous_energy = total_energy;
		previous_residual = total_residual;
		//std::cout << "residual " << residual << std::endl;
	}
	thread->assignTask(this, VELOCITY_NEWTON);
	updateNormal();
	updateRenderPosition();
	//thread->assignTask(this, VELOCITY_NEWTON_2);
	//log_store_residual.resize(store_residual.size());
	//if (iteration_number > 10000) {
		//std::cout << iteration_number << std::endl;
		//system("pause");
		//std::string txt_file_name = "residual result " + std::to_string(iteration_number);
		//WriteTxt::writeTxt(store_residual, txt_file_name);
		//for (unsigned int i = 0; i < store_residual.size(); ++i) {
		//	log_store_residual[i] = log10(store_residual[i]);
		//}
		//txt_file_name = "residual result log" + std::to_string(iteration_number);
		//WriteTxt::writeTxt(log_store_residual, txt_file_name);
	//}


}

void SecondOrderLargeSystem::updateNormal()
{
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	}
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, VERTEX_NORMAL);
	}
}


void SecondOrderLargeSystem::initialDHatTolerance(double ave_edge_length)
{
	if (perform_collision) {
		collision.initialDHatTolerance(ave_edge_length);
	}
}


bool SecondOrderLargeSystem::convergenceCondition()
{
	if (iteration_number < 2) {
		return true;
	}

	if (iteration_number > max_itr_num) {
		return false;
	}

	//if (abs((total_energy- previous_energy)/ previous_energy)< energy_conv_rate) {
	//	return false;
	//}

	//return true;

	double a = 0.0;
	double b = 0.0;
	for (unsigned int i = 0; i < delta_x.size(); i++) {
		if (a < abs(delta_x.data()[i])) {
			a = abs(delta_x.data()[i]);
			//b = delta_x.data()[i];
		}
	}

	////std::cout << b << std::endl;

	if (a > conv_rate) {
		return true;
	}
	return false;
}

void SecondOrderLargeSystem::updateRenderPosition()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(render_position[i][0].data(), vertex_position[i][0].data(), 24 * total_vertex_num[i]);
	}
}




//UPDATE_INTERNAL_FORCE
void SecondOrderLargeSystem::updateInternalForce(int thread_No)
{

	unsigned int vertex_begin;
	unsigned int vertex_end;
	unsigned int* edge_vertex;
	std::array<double, 3>* vertex_pos;

	double* rest_length_;

	unsigned int* unfixed_vertex_index;

	unsigned int index_start;

	double* f_ = b_thread[thread_No].data();

	double* h = b_thread[0].data() + unfixed_constraint_start_per_thread_in_system[thread_No];
	double* h_fixed_one_vertex = b_thread[0].data() + fixed_vertex_constraint_start_per_thread_in_system[thread_No];

	double alpha;
	double* lambda;

	std::array<double, 3>* ori_vertex_pos;

	double damp_ = rayleigh_damp_stiffness[0] * time_step;

	for (unsigned int obj_No = 0; obj_No < cloth->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No]->data();
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		rest_length_ = unfixed_rest_length[obj_No]->data();

		unfixed_vertex_index = unfixed_vertex[obj_No]->data();
		index_start = vertex_begin_per_obj[obj_No];
		alpha = 1.0 / (*edge_length_stiffness[obj_No] * time_step_square);
		lambda = lambda_unfixed[obj_No].data();
		ori_vertex_pos = render_position[obj_No];
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			updateInternalForce(vertex_pos[unfixed_vertex_index[edge_vertex[i]]].data(), vertex_pos[unfixed_vertex_index[edge_vertex[i + 1]]].data(),
				f_ + 3 * (edge_vertex[i] + index_start), f_ + 3 * (edge_vertex[i + 1] + index_start),
				rest_length_[i >> 1], h, lambda[i>>1], alpha, ori_vertex_pos[unfixed_vertex_index[edge_vertex[i]]].data(), ori_vertex_pos[unfixed_vertex_index[edge_vertex[i + 1]]].data(), 
				damp_);
			h++;
		}
		//only one vertex fixed in the edge
		edge_vertex = only_one_vertex_fix_edge_vertices[obj_No]->data();
		vertex_begin = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No] << 1;
		vertex_end = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No + 1] << 1;
		rest_length_ = fixed_one_vertices_rest_length[obj_No]->data();
		lambda = lambda_fixed_one_vertex[obj_No].data();
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			updateInternalForceOnlyOneEdgeFixed(vertex_pos[unfixed_vertex_index[edge_vertex[i]]].data(), vertex_pos[edge_vertex[i + 1]].data(),
				f_ + 3 * (edge_vertex[i]+ index_start),
				alpha, rest_length_[i >> 1], h_fixed_one_vertex, lambda[i >> 1],
				ori_vertex_pos[unfixed_vertex_index[edge_vertex[i]]].data(), ori_vertex_pos[edge_vertex[i + 1]].data(), damp_);
			h_fixed_one_vertex++;
		}
	}
}


//void SecondOrderLargeSystem::updateEnergy



void SecondOrderLargeSystem::computeMassSpringEnergy(int thread_No)
{
	unsigned int vertex_begin;
	unsigned int vertex_end;
	unsigned int* edge_vertex;
	std::array<double, 3>* vertex_pos;

	double* rest_length_;
	double edge_length_stiffness_;

	unsigned int* unfixed_vertex_index;
	unsigned int index_start;
	double energy = 0;
	double* lambda;
	for (unsigned int obj_No = 0; obj_No < cloth->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No]->data();
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		rest_length_ = unfixed_rest_length[obj_No]->data();
		unfixed_vertex_index = unfixed_vertex[obj_No]->data();
		lambda = lambda_unfixed[obj_No].data();
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			energy += compute_energy.computeMassSpringConstraint(vertex_pos[unfixed_vertex_index[edge_vertex[i]]].data(), vertex_pos[unfixed_vertex_index[edge_vertex[i + 1]]].data(),
				rest_length_[i >> 1], *edge_length_stiffness[obj_No], lambda[i >> 1],time_step_square);
		}
		//only one vertex fixed in the edge
		edge_vertex = only_one_vertex_fix_edge_vertices[obj_No]->data();
		vertex_begin = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No] << 1;
		vertex_end = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No + 1] << 1;
		rest_length_ = fixed_one_vertices_rest_length[obj_No]->data();
		lambda = lambda_fixed_one_vertex[obj_No].data();
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			energy += compute_energy.computeMassSpringConstraint(vertex_pos[unfixed_vertex_index[edge_vertex[i]]].data(), vertex_pos[edge_vertex[i + 1]].data(),
				rest_length_[i >> 1], *edge_length_stiffness[obj_No], lambda[i >> 1], time_step_square);
		}
	}
	energy_per_thread[thread_No] += energy;
}



//UPDATE_HESSIAN_FIXED_STRUCTURE
void SecondOrderLargeSystem::updateHessianFixedStructure(int thread_No)
{
	unsigned int vertex_begin;
	unsigned int vertex_end;
	unsigned int* edge_vertex;
	std::array<double, 3>* vertex_pos;

	memset(hessian_coeff_diagonal[thread_No].data(), 0, 8 * hessian_coeff_diagonal[thread_No].size());
	double* hessian_coeff_diag = hessian_coeff_diagonal[thread_No].data();

	double** hessian_nnz_ = Hessian_coeff_address.data() + off_diagonal_hessian_nnz_index_begin_per_thread[thread_No];
	unsigned int edge_no = 0;
	unsigned int edge_no_fixed_one_vertex = 0;

	double time_step_square_ = time_step_square;
	unsigned int vertex_begin_in_obj;
	unsigned int* unfixed_index_to_normal_index;
	double alpha;

	double** address_for_grad_c = Hessian_coeff_address.data() + unfixed_gradC_hessian_index_begin_per_thread[thread_No];
	double** address_for_fixed_vertex_grad_C = Hessian_coeff_address.data() + fixed_gradC_hessian_index_begin_per_thread[thread_No];

	double* lambda_;
	std::array<double, 3>* ori_vertex_pos;
	double damp_stiff = rayleigh_damp_stiffness[0] * time_step;
	for (unsigned int obj_No = 0; obj_No < cloth->size(); ++obj_No) {
		alpha = 1.0 / (*edge_length_stiffness[obj_No] * time_step_square_);
		vertex_begin_in_obj = vertex_begin_per_obj[obj_No];
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No]->data();
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		lambda_ = lambda_unfixed[obj_No].data();
		unfixed_index_to_normal_index = unfixed_vertex[obj_No]->data();

		ori_vertex_pos = render_position[obj_No];
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			computeHessianFixedStructure(vertex_pos[unfixed_index_to_normal_index[edge_vertex[i]]].data(),
				vertex_pos[unfixed_index_to_normal_index[edge_vertex[i + 1]]].data(),
				hessian_coeff_diag + 9 * (edge_vertex[i] + vertex_begin_in_obj), hessian_coeff_diag + 9 * (edge_vertex[i + 1] + vertex_begin_in_obj),
				hessian_nnz_ + 18 * edge_no, alpha, address_for_grad_c +13*edge_no, lambda_[i>>1],
				ori_vertex_pos[unfixed_index_to_normal_index[edge_vertex[i]]].data(),
				ori_vertex_pos[unfixed_index_to_normal_index[edge_vertex[i + 1]]].data(), damp_stiff);
			edge_no++;
		}
		//only one vertex fixed in the edge
		vertex_begin = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No] << 1;
		vertex_end = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No + 1] << 1;
		edge_vertex = only_one_vertex_fix_edge_vertices[obj_No]->data();
		lambda_ = lambda_fixed_one_vertex[obj_No].data();
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			computeHessianOnlyOneVertexFixedEdge(vertex_pos[unfixed_index_to_normal_index[edge_vertex[i]]].data(),
				vertex_pos[edge_vertex[i + 1]].data(),
				hessian_coeff_diag + 9 * (edge_vertex[i] + vertex_begin_in_obj),
				alpha, address_for_fixed_vertex_grad_C+7* edge_no_fixed_one_vertex,
				lambda_[i>>1], ori_vertex_pos[unfixed_index_to_normal_index[edge_vertex[i]]].data(),
				ori_vertex_pos[edge_vertex[i + 1]].data(), damp_stiff);
			edge_no_fixed_one_vertex++;
		}
	}
}

void SecondOrderLargeSystem::initialLambda()
{
	lambda_unfixed.resize(cloth->size());
	lambda_fixed_one_vertex.resize(cloth->size());
	tet_lambda.resize(tetrahedron->size());
	lambda_test.resize(tetrahedron->size());

	unsigned int constraint_number = 0;
	unsigned int constraint_number1 = 0;
	for (unsigned int i = 0; i < cloth->size(); ++i) {

		lambda_unfixed[i].resize(cloth->data()[i].mesh_struct.unfixed_rest_edge_length.size());
		lambda_fixed_one_vertex[i].resize( cloth->data()[i].mesh_struct.fixed_one_vertex_rest_edge_length.size());
	}

	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		tet_lambda[i].resize(tetrahedron->data()[i].mesh_struct.indices.size());
		lambda_test[i].resize(tetrahedron->data()[i].mesh_struct.indices.size());
	}
	resetLambda();
}






//GET_COEFF_ADDRESS
void SecondOrderLargeSystem::hessianCoeffAddress(int thread_No)
{
	unsigned int vertex_begin;
	unsigned int vertex_end;
	unsigned int* edge_vertex;



	double** address = Hessian_coeff_address.data() + off_diagonal_hessian_nnz_index_begin_per_thread[thread_No];
	double** address_for_grad_c = Hessian_coeff_address.data() + unfixed_gradC_hessian_index_begin_per_thread[thread_No];
	double** address_for_fixed_vertex_grad_C = Hessian_coeff_address.data() + fixed_gradC_hessian_index_begin_per_thread[thread_No];

	unsigned int unfixed_constraint_start = unfixed_constraint_start_per_thread_in_system[thread_No];
	unsigned int fixed_vertex_constraint_start = fixed_vertex_constraint_start_per_thread_in_system[thread_No];

	unsigned int edge_no = 0;
	unsigned int fixed_edge_no = 0;
	unsigned int vertex_index_begin_in_system;

	//off diagonal mass spring
	for (unsigned int obj_No = 0; obj_No < cloth->size(); ++obj_No) {
		edge_vertex = edge_vertices_mass_spring[obj_No]->data();
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		vertex_index_begin_in_system = vertex_begin_per_obj[obj_No];
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			getHessianCoeffAddress(address + 18 * edge_no, 3 * (vertex_index_begin_in_system + edge_vertex[i]),
				3 * (vertex_index_begin_in_system + edge_vertex[i + 1]), address_for_grad_c+13*edge_no, unfixed_constraint_start);
			edge_no++;
			unfixed_constraint_start++;
		}

		vertex_begin = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No] << 1;
		vertex_end = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No + 1] << 1;
		edge_vertex = only_one_vertex_fix_edge_vertices[obj_No]->data();
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			getHessianCoeffAddressFixedOneAddress(3 * (vertex_index_begin_in_system + edge_vertex[i]), address_for_fixed_vertex_grad_C + 7 * fixed_edge_no, fixed_vertex_constraint_start);
			fixed_edge_no++;
			fixed_vertex_constraint_start++;
		}
	}
	//diagonal
	unsigned int index_end;
	unsigned int hessian_index_start;
	unsigned int global_matrix_row_start;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		for (unsigned int j = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			hessian_index_start = 9 * j;
			global_matrix_row_start = 3 * j;
			for (unsigned int k = 0; k < 3; ++k) {
				Hessian_coeff_address[hessian_index_start] = &Hessian.coeffRef(global_matrix_row_start, global_matrix_row_start + k);
				Hessian_coeff_address[hessian_index_start + 1] = Hessian_coeff_address[hessian_index_start] + 1;
				Hessian_coeff_address[hessian_index_start + 2] = Hessian_coeff_address[hessian_index_start] + 2;
				hessian_index_start += 3;
			}
		}
	}

}
void SecondOrderLargeSystem::initialCoeffDiagonal()
{
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		memset(hessian_coeff_diagonal[i].data(), 0, 8 * hessian_coeff_diagonal[i].size());
	}
}



void SecondOrderLargeSystem::getARAPCoeffAddress()
{
	unsigned int size;
	int vertex_position_in_system[4];
	bool is_unfixed[4];
	std::array<int, 4>* indices;
	double* inv_mass_;
	unsigned int vertex_index_start;
	unsigned int* real_index_to_system_index;
	unsigned int constraint_index_in_system_ = ARAP_constrant_in_system[0];

	unsigned int record_size = 0;
	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
		size = tet_mesh_struct[obj_No]->indices.size();
		inv_mass_ = inv_mass[cloth->size() + obj_No];
		indices = tet_indices[obj_No];
		vertex_index_start = vertex_begin_per_obj[cloth->size() + obj_No];
		real_index_to_system_index = real_index_to_unfixed_index[cloth->size() + obj_No];

		for (unsigned int i = 0; i < size; ++i) {
			for (unsigned int j = 0; j < 4; ++j) {
				is_unfixed[j] = (inv_mass_[indices[i][j]] != 0.0);
				vertex_position_in_system[j] = 3 * (vertex_index_start + real_index_to_system_index[indices[i][j]]);				
			}

			if (is_unfixed[0] || is_unfixed[1] || is_unfixed[2] || is_unfixed[3]) {

				record_size = Hessian_coeff_address.size();

				getARAPHessianCoeffAddress(vertex_position_in_system, &Hessian_coeff_address, constraint_index_in_system_, is_unfixed);

				constraint_index_in_system_++;
			}
		}
	}

	//std::cout << Hessian_coeff_address.size() << std::endl;;
}


void SecondOrderLargeSystem::computeARAPEnergy()
{
	unsigned int tet_end;

	std::array<double, 3>* vertex_pos;
	std::array<int, 4>* tet_index;

	Matrix<double, 3, 4>* A;
	double* volume;
	double stiffness;
	double energy = 0.0;
	double* mass_inv_;
	double* lambda_;
	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No+cloth->size()];
		tet_end = tet_index_begin_per_thread[obj_No][total_thread_num];
		tet_index = tet_indices[obj_No];
		A = tet_mesh_struct[obj_No]->A.data();
		volume = tet_mesh_struct[obj_No]->volume.data();
		stiffness = tetrahedron->data()[obj_No].ARAP_stiffness;
		mass_inv_ = inv_mass[obj_No + cloth->size()];
		lambda_ = lambda_test[obj_No].data();
		for (unsigned int i =0; i < tet_end; i ++) {
			if (mass_inv_[tet_index[i][0]] != 0.0 || mass_inv_[tet_index[i][1]] != 0.0 || mass_inv_[tet_index[i][2]] != 0.0 || mass_inv_[tet_index[i][3]] != 0.0) {
				//energy += compute_energy.computeARAPConstraint(vertex_pos[tet_index[i][0]].data(), vertex_pos[tet_index[i][1]].data(), vertex_pos[tet_index[i][2]].data(),
				//	vertex_pos[tet_index[i][3]].data(), A[i], volume[i], stiffness, *lambda_,time_step_square);
				//lambda_++;
				energy += compute_energy.computeARAPEnergy(vertex_pos[tet_index[i][0]].data(), vertex_pos[tet_index[i][1]].data(), vertex_pos[tet_index[i][2]].data(),
					vertex_pos[tet_index[i][3]].data(), A[i], volume[i], stiffness);
			}

		}
	}
	total_energy += energy;
}

//also compute internal force
void SecondOrderLargeSystem::updateARAPHessianFixedStructure()
{

	std::array<double, 3>* vertex_pos;
	double* hessian_coeff_diag = hessian_coeff_diagonal[0].data();
	double alpha;
	unsigned int size;
	double* lambda;

	double* inv_mass_;

	bool is_unfixed[4];

	std::array<int, 4>* indices;
	unsigned int vertex_index_start;

	int vertex_position_in_system[4];
	unsigned int* real_index_to_system_index;

	double use_for_temp[9];
	double* diagonal_coeff[4];
	double* g_tet[4];
	double* volume;
	Matrix<double, 3, 4>* A;

	double** address = Hessian_coeff_address.data() +3* ARAP_constrant_in_system[0];

	double* g_ = b_thread[0].data();
	double* h = b_thread[0].data() + ARAP_constrant_in_system[0];


	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No + cloth->size()];
		alpha = 1.0 / (tetrahedron->data()[obj_No].ARAP_stiffness * time_step_square);
		size = tet_mesh_struct[obj_No]->indices.size();
		lambda = tet_lambda[obj_No].data();
		//lambda = lambda_test[obj_No].data();
		inv_mass_ = inv_mass[cloth->size() + obj_No];
		indices = tet_indices[obj_No];
		vertex_index_start = vertex_begin_per_obj[cloth->size() + obj_No];
		real_index_to_system_index = real_index_to_unfixed_index[cloth->size() + obj_No];
		volume = tet_mesh_struct[obj_No]->volume.data();
		A = tet_mesh_struct[obj_No]->A.data();

		for (unsigned int i = 0; i < size; ++i) {
			for (unsigned int j = 0; j < 4; ++j) {
				is_unfixed[j] = (inv_mass_[indices[i][j]] != 0.0);
				vertex_position_in_system[j] = 3 * (vertex_index_start + real_index_to_system_index[indices[i][j]]);
				if (is_unfixed[j]) {
					diagonal_coeff[j] = hessian_coeff_diag + 3 * vertex_position_in_system[j];
					g_tet[j] = g_+vertex_position_in_system[j];
				}
				else {
					diagonal_coeff[j] = use_for_temp;
					g_tet[j] = use_for_temp;
				}
			}
			if (is_unfixed[0] || is_unfixed[1] || is_unfixed[2] || is_unfixed[3]) {
				//if (i == 25) {
				//	std::cout <<"read h start "<<  h - b_thread[0].data() << std::endl;
				//}
				computeARAPHessianFixedStructure(vertex_pos[indices[i][0]].data(), vertex_pos[indices[i][1]].data(), vertex_pos[indices[i][2]].data(), vertex_pos[indices[i][3]].data(),
					diagonal_coeff[0], diagonal_coeff[1], diagonal_coeff[2], diagonal_coeff[3], address, alpha / volume[i], A[i], is_unfixed,
					lambda, h, g_tet[0], g_tet[1], g_tet[2], g_tet[3],i);
				h++;
				lambda++;

				//if (i == 25) {
				//	std::cout << "real " << indices[i][0] << " " << indices[i][1] << " " << indices[i][2] << " " << indices[i][3] << std::endl;
				//}
			}
		}
	}


}


void SecondOrderLargeSystem::checkIfSampleRight()
{
	std::array<double, 3>* vertex_pos;
	std::vector<Triplet<double>>* hessian_nnz_ = &hessian_nnz;
	double alpha;
	unsigned int size;
	double* lambda;
	double* inv_mass_;

	bool is_unfixed[4];
	std::array<int, 4>* indices;
	unsigned int vertex_index_start;

	int vertex_position_in_system[4];
	unsigned int* real_index_to_system_index;

	unsigned int constraint_index_in_system_ = ARAP_constrant_in_system[0];
	double* volume;
	Matrix<double, 3, 4>* A;
	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No + cloth->size()];
		alpha = 1.0 / (tetrahedron->data()[obj_No].ARAP_stiffness * time_step_square);
		size = tet_mesh_struct[obj_No]->indices.size();
		lambda = tet_lambda[obj_No].data();
		inv_mass_ = inv_mass[cloth->size() + obj_No];
		indices = tet_indices[obj_No];
		vertex_index_start = vertex_begin_per_obj[cloth->size() + obj_No];
		real_index_to_system_index = real_index_to_unfixed_index[cloth->size() + obj_No];
		volume = tet_mesh_struct[obj_No]->volume.data();
		A = tet_mesh_struct[obj_No]->A.data();

		for (unsigned int i = 0; i < size; ++i) {
			for (unsigned int j = 0; j < 4; ++j) {
				is_unfixed[j] = (inv_mass_[indices[i][j]] != 0.0);
				vertex_position_in_system[j] = 3 * (vertex_index_start + real_index_to_system_index[indices[i][j]]);
			}
			if (is_unfixed[0] || is_unfixed[1] || is_unfixed[2] || is_unfixed[3]) {
				checkIfSampleRight(vertex_pos[indices[i][0]].data(), vertex_pos[indices[i][1]].data(), vertex_pos[indices[i][2]].data(), vertex_pos[indices[i][3]].data(),
					vertex_position_in_system,
					constraint_index_in_system_, alpha / volume[i], A[i], is_unfixed, lambda);
				constraint_index_in_system_++;
				lambda++;
			}
		}
	}
}

void SecondOrderLargeSystem::updateTest()
{
	b_test.setZero();
	test_nnz.clear();

	setARAP_ForTest();
	Matrix_test.setFromTriplets(test_nnz.begin(), test_nnz.end());


	

	global_llt.compute(Matrix_test);

	
	//global_llt1.compute(Matrix_test);

	//std::cout <<"determinent "<< global_llt1.determinant()<<" "<< global_llt1.logAbsDeterminant() << std::endl;

	//std::cout << Matrix_test << std::endl;

	delta_x = global_llt.solve(b_test);

	//std::cout << b_test << std::endl;
	//std::cout << "====" << std::endl;

	//std::cout << "error " << (Matrix_test * delta_x - b_test).norm() << std::endl;


}

void SecondOrderLargeSystem::setARAP_ForTest()
{

	std::array<double, 3>* vertex_pos;
	double alpha;
	unsigned int size;
	double* lambda;
	double* inv_mass_;

	bool is_unfixed[4];
	std::array<int, 4>* indices;
	unsigned int vertex_index_start;

	int vertex_position_in_system[4];
	unsigned int* real_index_to_system_index;

	double use_for_temp[9];
	double* diagonal_coeff[4];

	unsigned int constraint_index_in_system_ = ARAP_constrant_in_system[0];

	//std::cout << "ARAP_constrant_in_system[0] " << ARAP_constrant_in_system[0] << std::endl;

	double* volume;
	Matrix<double, 3, 4>* A;

	std::array<double, 3>* ori_vertex_pos;
	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No + cloth->size()];
		alpha = 1.0 / (tetrahedron->data()[obj_No].ARAP_stiffness * time_step_square);
		size = tet_mesh_struct[obj_No]->indices.size();
		lambda = lambda_test[obj_No].data();
		inv_mass_ = inv_mass[cloth->size() + obj_No];
		indices = tet_indices[obj_No];
		vertex_index_start = vertex_begin_per_obj[cloth->size() + obj_No];
		real_index_to_system_index = real_index_to_unfixed_index[cloth->size() + obj_No];
		volume = tet_mesh_struct[obj_No]->volume.data();
		A = tet_mesh_struct[obj_No]->A.data();
		ori_vertex_pos = render_position[obj_No + cloth->size()];

		for (unsigned int i = 0; i < size; ++i) {
			for (unsigned int j = 0; j < 4; ++j) {
				is_unfixed[j] = (inv_mass_[indices[i][j]] != 0.0);
				vertex_position_in_system[j] = 3 * (vertex_index_start + real_index_to_system_index[indices[i][j]]);
			}
			if (is_unfixed[0] || is_unfixed[1] || is_unfixed[2] || is_unfixed[3]) {
				setARAPHessianForTest(vertex_pos[indices[i][0]].data(), vertex_pos[indices[i][1]].data(), vertex_pos[indices[i][2]].data(),
					vertex_pos[indices[i][3]].data(),
					&test_nnz, vertex_position_in_system,
					constraint_index_in_system_, alpha / volume[i], A[i], is_unfixed, lambda,b_test.data()+ vertex_position_in_system[0],
					b_test.data() + vertex_position_in_system[1], b_test.data() + vertex_position_in_system[2], b_test.data() + vertex_position_in_system[3],
					b_test.data() + constraint_index_in_system_,i, ori_vertex_pos[indices[i][0]].data(), ori_vertex_pos[indices[i][1]].data(), ori_vertex_pos[indices[i][2]].data(),
					ori_vertex_pos[indices[i][3]].data(), rayleigh_damp_stiffness[0]*time_step);

				//if (i == 25) {
					//std::cout << *lambda << std::endl;
					//std::cout << "test h start " << constraint_index_in_system_ << std::endl;
					//std::cout << "test " << constraint_index_in_system_ - 1 << " " << indices[i][0] << " " << indices[i][1] << " " << indices[i][2] << " " << indices[i][3] << std::endl;
				//}
				constraint_index_in_system_++;
				lambda++;

		

			}
		}

		double mass_;
		unsigned int actual_index;
		unsigned int* unfixed_index_to_normal_index = unfixed_vertex[cloth->size()+obj_No]->data();
		unsigned int index_end = vertex_begin_per_obj[cloth->size() + obj_No] + unfixed_vertex[cloth->size() + obj_No]->size();



		for (unsigned int j = vertex_begin_per_obj[cloth->size() + obj_No]; j < index_end; ++j) {
			actual_index = unfixed_index_to_normal_index[j - vertex_begin_per_obj[cloth->size() + obj_No]];
			mass_ = mass[cloth->size() + obj_No][actual_index];
			test_nnz.emplace_back(Triplet<double>(3*j, 3*j, mass_));
			test_nnz.emplace_back(Triplet<double>(3 * j + 1, 3 * j + 1, mass_));
			test_nnz.emplace_back(Triplet<double>(3 * j + 2, 3 * j + 2, mass_));

			b_test.data()[3 * j] -= mass_ * (vertex_pos[actual_index][0] - Sn.data()[3 * j]);
			b_test.data()[3 * j+1] -= mass_ * (vertex_pos[actual_index][1] - Sn.data()[3 * j+1]);
			b_test.data()[3 * j+2] -= mass_ * (vertex_pos[actual_index][2] - Sn.data()[3 * j+2]);
		}


		ARAP_constrant_in_system[obj_No + 1] = constraint_index_in_system_;
		lambda_test[obj_No].resize(constraint_index_in_system_ - ARAP_constrant_in_system[obj_No]);
	}


}


void SecondOrderLargeSystem::setARAP()
{
	std::array<double, 3>* vertex_pos;
	double* hessian_coeff_diag = hessian_coeff_diagonal[0].data();
	std::vector<Triplet<double>>* hessian_nnz_ = &hessian_nnz;
	double alpha;
	unsigned int size;
	double* lambda;	
	double* inv_mass_;

	bool is_unfixed[4];
	std::array<int, 4>* indices;
	unsigned int vertex_index_start;

	int vertex_position_in_system[4];
	unsigned int* real_index_to_system_index;

	double use_for_temp[9];
	double* diagonal_coeff[4];

	unsigned int constraint_index_in_system_= ARAP_constrant_in_system[0];
	double* volume;
	Matrix<double, 3, 4>* A;
	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No+cloth->size()];
		alpha= 1.0 /( tetrahedron->data()[obj_No].ARAP_stiffness * time_step_square);
		size = tet_mesh_struct[obj_No]->indices.size();
		lambda = tet_lambda[obj_No].data();
		inv_mass_ = inv_mass[cloth->size() + obj_No];
		indices = tet_indices[obj_No];
		vertex_index_start = vertex_begin_per_obj[cloth->size()+ obj_No];
		real_index_to_system_index = real_index_to_unfixed_index[cloth->size() + obj_No];
		volume = tet_mesh_struct[obj_No]->volume.data();
		A = tet_mesh_struct[obj_No]->A.data();
		
		for (unsigned int i = 0; i < size; ++i) {
			for (unsigned int j = 0; j < 4; ++j) {
				is_unfixed[j] = (inv_mass_[indices[i][j]] != 0.0);
				vertex_position_in_system[j] = 3 * (vertex_index_start + real_index_to_system_index[indices[i][j]]);
				if (is_unfixed[j]) {
					diagonal_coeff[j] = hessian_coeff_diag + 3 * vertex_position_in_system[j];
				}
				else {
					diagonal_coeff[j] = use_for_temp;
				}
			}			
			if (is_unfixed[0] || is_unfixed[1] || is_unfixed[2] || is_unfixed[3]) {
				computeARAPHessian(vertex_pos[indices[i][0]].data(), vertex_pos[indices[i][1]].data(), vertex_pos[indices[i][2]].data(), vertex_pos[indices[i][3]].data(),
					diagonal_coeff[0], diagonal_coeff[1], diagonal_coeff[2], diagonal_coeff[3], hessian_nnz_, vertex_position_in_system,
					constraint_index_in_system_, alpha / volume[i], A[i], is_unfixed, lambda);
				constraint_index_in_system_++;
				lambda++;
			}
		}
		ARAP_constrant_in_system[obj_No + 1] = constraint_index_in_system_;
		tet_lambda[obj_No].resize(constraint_index_in_system_ - ARAP_constrant_in_system[obj_No]);
	}
	
}



// SET_MASS_SPRING
void SecondOrderLargeSystem::massSpring(int thread_No)
{
	unsigned int edge_begin;
	unsigned int edge_end;
	unsigned int* edge_vertex;
	std::array<double, 3>* vertex_pos;

	
	double* hessian_coeff_diag = hessian_coeff_diagonal[thread_No].data();

	Triplet<double>* hessian_nnz_ = hessian_nnz.data() + off_diagonal_hessian_nnz_index_begin_per_thread[thread_No];
	Triplet<double>* hessian_nnz_grad_c = hessian_nnz.data() + unfixed_gradC_hessian_index_begin_per_thread[thread_No];
	Triplet<double>* hessian_nnz_grad_c_fixed = hessian_nnz.data() + fixed_gradC_hessian_index_begin_per_thread[thread_No];

	unsigned int edge_no = 0;
	unsigned int fixed_one_edge = 0;
	double time_step_square_ = time_step_square;

	double* rest_length_;

	double alpha;
	unsigned int vertex_index_begin_in_system;

	unsigned int* unfixed_index_to_normal_index;

	double* lambda;
	unsigned int sys_start_for_constraint;

	unsigned int unfixed_constraint_start = unfixed_constraint_start_per_thread_in_system[thread_No];
	unsigned int fixed_vertex_constraint_start = fixed_vertex_constraint_start_per_thread_in_system[thread_No];


	for (unsigned int obj_No = 0; obj_No < cloth->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No]->data();
		edge_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		edge_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		rest_length_ = unfixed_rest_length[obj_No]->data();
		alpha = 1.0 /(*edge_length_stiffness[obj_No]* time_step_square_);
		vertex_index_begin_in_system = vertex_begin_per_obj[obj_No];
		unfixed_index_to_normal_index = unfixed_vertex[obj_No]->data();
		lambda = lambda_unfixed[obj_No].data();
		for (unsigned int i = edge_begin; i < edge_end; i += 2) {
			computeHessian(vertex_pos[unfixed_index_to_normal_index[edge_vertex[i]]].data(), vertex_pos[unfixed_index_to_normal_index[edge_vertex[i + 1]]].data(),
				hessian_coeff_diag + 9 * (edge_vertex[i] + vertex_index_begin_in_system), hessian_coeff_diag + 9 * (edge_vertex[i + 1] + vertex_index_begin_in_system),
				hessian_nnz_ + 18 * edge_no, alpha, rest_length_[i >> 1], 3 * (vertex_index_begin_in_system + edge_vertex[i]),
				3 * (vertex_index_begin_in_system + edge_vertex[i + 1]), hessian_nnz_grad_c + 13 * edge_no, lambda[i>>1], 
				unfixed_constraint_start);
			unfixed_constraint_start++;
			edge_no++;
		}
		//only one vertex fixed in the edge
		edge_begin = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No] << 1;
		edge_end = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No + 1] << 1;
		edge_vertex = only_one_vertex_fix_edge_vertices[obj_No]->data();
		rest_length_ = fixed_one_vertices_rest_length[obj_No]->data();
		lambda = lambda_fixed_one_vertex[obj_No].data();
		for (unsigned int i = edge_begin; i < edge_end; i += 2) {
			computeHessianFixedOneVertex(vertex_pos[unfixed_index_to_normal_index[edge_vertex[i]]].data(),
				vertex_pos[edge_vertex[i + 1]].data(),
				hessian_coeff_diag + 9 * (edge_vertex[i] + vertex_index_begin_in_system),
				alpha, rest_length_[i >> 1],
				3 * (edge_vertex[i] + vertex_index_begin_in_system), hessian_nnz_grad_c_fixed + 7 * fixed_one_edge, lambda[i>>1],
				fixed_vertex_constraint_start);
				fixed_one_edge++;
				fixed_vertex_constraint_start++;
		}
	}
}


void SecondOrderLargeSystem::updatePositionFromSn()
{
	unsigned int index_end;
	unsigned int index_start;

	double* vertex_pos;

	unsigned int vertex_start;
	unsigned int* unfixed_index_to_normal_index;
	unsigned int j;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_start = vertex_begin_per_obj[i];
		vertex_pos = vertex_position[i][0].data();		
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		index_end = unfixed_vertex[i]->size();
		for (unsigned int l = 0; l < index_end; ++l) {
			memcpy(vertex_pos + 3 * unfixed_index_to_normal_index[l], Sn.data() + 3 * (l + vertex_start), 24);
		}
	}
}

//SET_S_N
void SecondOrderLargeSystem::setSn(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	double time_step_ = time_step;
	double time_step_square_ = time_step_square;
	double* vertex_pos;
	double* mass_;

	unsigned int vertex_start;
	unsigned int* unfixed_index_to_normal_index;

	unsigned int j;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		mass_ = mass[i];
		vertex_start = vertex_begin_per_obj[i];
		vertex_pos = vertex_position[i][0].data();
		index_end = vertex_start + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_start + unfixed_vertex_begin_per_thread[i][thread_No];

		unfixed_index_to_normal_index = unfixed_vertex[i]->data();

		for (unsigned int l = index_start; l < index_end; ++l) {
			j = 3 * l;
			for (unsigned int k = 0; k < 3; ++k) {
				Sn.data()[j + k] =
					(vertex_pos[3 * unfixed_index_to_normal_index[l - vertex_start] + k]
						+ time_step_ * velocity.data()[j + k]) + time_step_square_ * f_ext.data()[j + k] / mass_[unfixed_index_to_normal_index[l - vertex_start]];
			}
		}
	}
}


void SecondOrderLargeSystem::computeEnergy()
{
	thread->assignTask(this, NEWTON_METHOD_ENERGY);
	total_energy = 0;
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		total_energy += energy_per_thread[i];
	}
	computeARAPEnergy();

}

//NEWTON_METHOD_ENERGY
void SecondOrderLargeSystem::computeEnergy(int thread_No)
{
	energy_per_thread[thread_No] = 0;
	computeMassSpringEnergy(thread_No);
	computeInertial(thread_No);

}

void SecondOrderLargeSystem::computeResidual()
{
	double* mass_;
	Vector3d* residual_;
	unsigned int vertex_size;


	unsigned int index_end;
	unsigned int index_start;
	std::array<double,3>* vertex_pos;

	unsigned int vertex_start;
	unsigned int* unfixed_index_to_normal_index;

	unsigned int j;
	unsigned int start;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		mass_ = mass[i];
		vertex_pos = vertex_position[i];
		residual_ = residual[i].data();
		vertex_size = unfixed_vertex_begin_per_thread[i][total_thread_num];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		for (unsigned int k = 0; k < vertex_size; ++k) {
			start = 3 * (vertex_begin_per_obj[i] + k);
			j = unfixed_index_to_normal_index[k];
			residual_[j][0] = mass_[j] * (vertex_pos[j][0] - Sn[start]) / (time_step_square);
			residual_[j][1] = mass_[j] * (vertex_pos[j][1] - Sn[start+1]) / (time_step_square);
			residual_[j][2] = mass_[j] * (vertex_pos[j][2] - Sn[start+2]) / (time_step_square);
		}
	}

	computeEdgeLenthResidual();
	computeARAPResidual();
	double residual_norm = 0.0;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		residual_ = residual[i].data();
		for (unsigned int j = 0; j < residual[i].size(); ++j) {
			residual_norm += residual_[j].squaredNorm();
		}
	}

	total_residual = sqrt(residual_norm);

	//std::cout << "residual norm "<<residual_norm << std::endl;

}


void SecondOrderLargeSystem::computeInertial(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	double* vertex_pos;
	double* mass_;

	unsigned int vertex_start;
	unsigned int* unfixed_index_to_normal_index;

	unsigned int j;

	double energy = 0;
	unsigned int start;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		mass_ = mass[i];
		vertex_start = vertex_begin_per_obj[i];
		vertex_pos = vertex_position[i][0].data();
		index_end = vertex_start + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_start + unfixed_vertex_begin_per_thread[i][thread_No];

		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		energy_per_thread[thread_No] += compute_energy.computeInertial(time_step,index_start, index_end, vertex_pos, mass_, Sn, vertex_start, unfixed_index_to_normal_index);
	}
}





void SecondOrderLargeSystem::computeEdgeLenthResidual()
{

	unsigned int size;
	double stiffness;
	unsigned int* edge_vertex_index;
	Vector3d force_0;
	Vector3d force_1;

	std::array<double, 3>* v_p;
	Vector3d* residual_;
	double* edge_length;
	double* inv_mass_;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		size = cloth->data()[i].mesh_struct.edge_length.size();
		edge_vertex_index = cloth->data()[i].mesh_struct.edge_vertices.data();
		stiffness = cloth->data()[i].length_stiffness;
		v_p = vertex_position[i];
		residual_ = residual[i].data();
		edge_length = cloth->data()[i].mesh_struct.edge_length.data();
		inv_mass_ = inv_mass[i];
		for (unsigned int j = 0; j < size; ++j) {
			if (inv_mass_[edge_vertex_index[(j) << 1]] != 0 || inv_mass_[edge_vertex_index[((j) << 1) + 1]] != 0) {
				second_order_constraint.computeEdgeLengthForce(v_p[edge_vertex_index[(j) << 1]].data(),
					v_p[edge_vertex_index[((j) << 1) + 1]].data(), stiffness, force_0.data(), force_1.data(),
					edge_length[j]);
				if (inv_mass_[edge_vertex_index[(j) << 1]] != 0) {
					residual_[edge_vertex_index[(j) << 1]] += force_0;
				}
				if (inv_mass_[edge_vertex_index[((j) << 1) + 1]] != 0) {
					residual_[edge_vertex_index[((j) << 1) + 1]] += force_1;
				}
			}
		}
	}
}



void SecondOrderLargeSystem::computeARAPResidual()
{
	unsigned int size;
	double stiffness;
	unsigned int* edge_vertex_index;

	std::array<double, 3>* vertex_pos;
	Vector3d* residual_;
	std::array<int, 4>* indices;
	double* volume;

	Matrix<double, 3, 4>* A;

	Matrix<double, 3, 4> force;
	double* inv_mass_;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		size = tetrahedron->data()[i].mesh_struct.indices.size();
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		vertex_pos = vertex_position[i + cloth->size()];
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		stiffness = tetrahedron->data()[i].ARAP_stiffness;
		A = tetrahedron->data()[i].mesh_struct.A.data();
		residual_ = residual[i + cloth->size()].data();
		inv_mass_ = inv_mass[cloth->size() + i];
		for (unsigned int j = 0; j < size; ++j) {
			if (inv_mass_[indices[j][0]] != 0.0 || inv_mass_[indices[j][1]] != 0.0 || inv_mass_[indices[j][2]] != 0.0 || inv_mass_[indices[j][3]] != 0.0) {
				second_order_constraint.computeARAPForce(vertex_pos[indices[j][0]].data(), vertex_pos[indices[j][1]].data(),
					vertex_pos[indices[j][2]].data(), vertex_pos[indices[j][3]].data(), stiffness, A[j], volume[j], force);
				if (inv_mass_[indices[j][0]] != 0.0) {
					residual_[indices[j][0]] += force.col(0);
				}
				if (inv_mass_[indices[j][1]] != 0.0) {
					residual_[indices[j][1]] += force.col(1);
				}
				if (inv_mass_[indices[j][2]] != 0.0) {
					residual_[indices[j][2]] += force.col(2);
				}
				if (inv_mass_[indices[j][3]] != 0.0) {
					residual_[indices[j][3]] += force.col(3);
				}
			}		
		}
	}

}

//
void SecondOrderLargeSystem::updateARAPLambda()
{
	unsigned int size;
	double* lambda;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		//size = tet_lambda[i].size();
		//lambda = tet_lambda[i].data();
		size = lambda_test[i].size();
		lambda = lambda_test[i].data();
		for (unsigned int j = 0; j < size; ++j) {
			lambda[j] += delta_x[ARAP_constrant_in_system[i] + j];
		}
	}

}

//UPDATE_POSITION_NEWTON
void SecondOrderLargeSystem::updatePosition(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	double* vertex_pos;
	unsigned int index;

	unsigned int* unfixed_index_to_normal_index;

	unsigned int vertex_start;
	unsigned int initial_start= unfixed_constraint_start_per_thread_in_system[thread_No];
	unsigned int initial_start_fixed_vertex= fixed_vertex_constraint_start_per_thread_in_system[thread_No];

	double* lambda;
	unsigned int edge_No = 0;
	unsigned int edge_No_fixed_one_edge = 0;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i][0].data();
		index = unfixed_vertex_begin_per_thread[i][thread_No];
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		for (unsigned int j = index_start; j < index_end; ++j) {
			vertex_start = 3 * unfixed_index_to_normal_index[index];
			vertex_pos[vertex_start] += delta_x[3 * j];
			vertex_pos[vertex_start + 1] += delta_x[3 * j + 1];
			vertex_pos[vertex_start + 2] += delta_x[3 * j + 2];
			index++;
		}
	}
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		index_start = edge_index_begin_per_thread_for_mass_spring[i][thread_No];
		index_end = edge_index_begin_per_thread_for_mass_spring[i][thread_No+1];
		lambda = lambda_unfixed[i].data();

		for (unsigned int j = index_start; j < index_end; ++j) {
			lambda[j] += delta_x[initial_start + edge_No];
			edge_No++;
		}

		index_start = only_one_vertex_fixed_edge_index_begin_per_thread[i][thread_No];
		index_end = only_one_vertex_fixed_edge_index_begin_per_thread[i][thread_No + 1];
		lambda = lambda_fixed_one_vertex[i].data();
		for (unsigned int j = index_start; j < index_end; ++j) {
			lambda[j] += delta_x[initial_start_fixed_vertex + edge_No_fixed_one_edge];
			edge_No_fixed_one_edge++;
		}

	}

	
}

//
void SecondOrderLargeSystem::updateARAPLambdaFromOri()
{
	unsigned int size;
	double* lambda;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		size = lambda_test[i].size();
		lambda = lambda_test[i].data();
		for (unsigned int j = 0; j < size; ++j) {
			lambda[j] -= displacement_coe*delta_x[ARAP_constrant_in_system[i] + j];
		}
	}

}

//UPDATE_POSITION_NEWTON_FROM_ORI
void SecondOrderLargeSystem::updatePositionFromOri(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	double* vertex_pos;

	unsigned int index;

	unsigned int* unfixed_index_to_normal_index;
	unsigned int vertex_start;

	double coe = displacement_coe;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i][0].data();
		index = unfixed_vertex_begin_per_thread[i][thread_No];
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		for (unsigned int j = index_start; j < index_end; ++j) {
			vertex_start = 3 * unfixed_index_to_normal_index[index];
			for (unsigned int k = 0; k < 3; ++k) {
				vertex_pos[vertex_start + k] -= coe * delta_x[3 * j + k];
			}
			index++;
		}
	}
}




//VELOCITY_NEWTON
void SecondOrderLargeSystem::updateVelocity(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	double* vertex_pos;
	unsigned int index;

	unsigned int* unfixed_index_to_normal_index;
	unsigned int vertex_start;

	double time_step_ = time_step;
	double* vertex_render_pos;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i][0].data();
		vertex_render_pos = render_position[i][0].data();
		index = unfixed_vertex_begin_per_thread[i][thread_No];
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		for (unsigned int j = index_start; j < index_end; ++j) {
			vertex_start = 3 * unfixed_index_to_normal_index[index];
			velocity.data()[3 * j] = velocity_damp * (vertex_pos[vertex_start] - vertex_render_pos[vertex_start]) / time_step_;
			velocity.data()[3 * j+1] = velocity_damp * (vertex_pos[vertex_start+1] - vertex_render_pos[vertex_start+1]) / time_step_;
			velocity.data()[3 * j+2] = velocity_damp * (vertex_pos[vertex_start+2] - vertex_render_pos[vertex_start+2]) / time_step_;
			index++;
		}
	}
}


//SUM_B
void SecondOrderLargeSystem::sumB(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	unsigned int total_thread_num_ = total_thread_num;
	double* vertex_pos;
	double* mass_;

	unsigned int index;
	unsigned int* unfixed_index_to_normal_vertex;
	unsigned int j;



	for (unsigned int i = 0; i < total_obj_num; ++i) {
		mass_ = mass[i];
		vertex_pos = vertex_position[i][0].data();
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No];
		unfixed_index_to_normal_vertex = unfixed_vertex[i]->data();

		for (unsigned int l = index_start; l < index_end; ++l) {
			index = unfixed_index_to_normal_vertex[l - vertex_begin_per_obj[i]];
			j = 3 * l;
			for (unsigned int m = 0; m < 3; ++m) {
				for (unsigned int k = 1; k < total_thread_num_; ++k) {
					b_thread[0][j + m] += b_thread[k][j + m];
				}
				b_thread[0][j + m] += mass_[index] * (Sn.data()[j + m] - vertex_pos[3 * index + m]);
			}
		}
	}
	
}

//UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE_INITIAL_STIFFNESS
void SecondOrderLargeSystem::setHessianDiagonalFixedStructureInitialStiffness(int thread_No)
{
	unsigned int index_end;
	unsigned int vertex_start;

	double* hessian_coeff_diagonal_0 = hessian_coeff_diagonal[0].data();
	unsigned int total_thread_num_ = total_thread_num;

	unsigned int* unfixed_index_to_normal_index;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		for (unsigned int j = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			vertex_start = 9 * j;
			for (unsigned int m = 0; m < 9; ++m) {
				for (unsigned int k = 1; k < total_thread_num_; ++k) {
					hessian_coeff_diagonal_0[vertex_start + m] += hessian_coeff_diagonal[k][vertex_start + m];
				}
			}
			for (unsigned int k = 0; k < 9; ++k) {
				*Hessian_coeff_address[vertex_start + k] = hessian_coeff_diagonal_0[vertex_start + k];
			}
		}
	}
}

//UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE
void SecondOrderLargeSystem::setHessianDiagonalFixedStructure(int thread_No)
{
	unsigned int index_end;
	unsigned int vertex_start;

	double* hessian_coeff_diagonal_0 = hessian_coeff_diagonal[0].data();
	unsigned int total_thread_num_ = total_thread_num;

	unsigned int* unfixed_index_to_normal_index;

	double mass_;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		for (unsigned int j = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			vertex_start = 9 * j;
			for (unsigned int m = 0; m < 9; ++m) {
				for (unsigned int k = 1; k < total_thread_num_; ++k) {
					hessian_coeff_diagonal_0[vertex_start + m] += hessian_coeff_diagonal[k][vertex_start + m];
				}
			}
			mass_ =mass[i][unfixed_index_to_normal_index[j - vertex_begin_per_obj[i]]];

			hessian_coeff_diagonal_0[vertex_start] += mass_;
			hessian_coeff_diagonal_0[vertex_start + 4] += mass_;
			hessian_coeff_diagonal_0[vertex_start + 8] += mass_;

			for (unsigned int k = 0; k < 9; ++k) {
				*Hessian_coeff_address[vertex_start + k] = hessian_coeff_diagonal_0[vertex_start + k];
			}

		}
	}
}


//SET_HESSIAN_DIAGONAL
void SecondOrderLargeSystem::setHessianDiagonal(int thread_No)
{
	unsigned int index_end;
	unsigned int vertex_start;
	unsigned int hessian_index_start;

	double* hessian_coeff_diagonal_0 = hessian_coeff_diagonal[0].data();
	unsigned int* unfixed_index_to_normal_index;
	unsigned int global_matrix_row_start;
	double mass_;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		for (unsigned int j = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			vertex_start = 9 * j;
			for (unsigned int m = 0; m < 9; ++m) {
				for (unsigned int k = 1; k < total_thread_num; ++k) {
					hessian_coeff_diagonal_0[vertex_start + m] += hessian_coeff_diagonal[k][vertex_start + m];
				}
			}
			mass_ = mass[i][unfixed_index_to_normal_index[j - vertex_begin_per_obj[i]]];
			hessian_coeff_diagonal_0[vertex_start] += mass_;
			hessian_coeff_diagonal_0[vertex_start + 4] += mass_;
			hessian_coeff_diagonal_0[vertex_start + 8] += mass_;
			hessian_index_start = 9 * j;
			global_matrix_row_start = 3 * j;
			for (unsigned int col = 0; col < 3; ++col) {
				for (unsigned int row = 0; row < 3; ++row) {
					hessian_nnz[hessian_index_start++] = Triplet<double>(global_matrix_row_start + row, global_matrix_row_start + col, hessian_coeff_diagonal_0[vertex_start++]);
				}
			}
		}
	}
}


void SecondOrderLargeSystem::getHessianCoeffAddressFixedOneAddress(unsigned int start_index_in_system_0, double** grad_C_address,
	unsigned int constraint_start_index_0)
{
	grad_C_address[0] = &Hessian.coeffRef(start_index_in_system_0, constraint_start_index_0);
	grad_C_address[1] = grad_C_address[0] + 1;
	grad_C_address[2] = grad_C_address[1] + 1;

	grad_C_address[3] = &Hessian.coeffRef(constraint_start_index_0, start_index_in_system_0);
	grad_C_address[4] = &Hessian.coeffRef(constraint_start_index_0, start_index_in_system_0 + 1);
	grad_C_address[5] = &Hessian.coeffRef(constraint_start_index_0, start_index_in_system_0 + 2);

	grad_C_address[6] = &Hessian.coeffRef(constraint_start_index_0, constraint_start_index_0);
}


void SecondOrderLargeSystem::getARAPHessianCoeffAddress(int* start_index_in_system, std::vector<double*>* address,
	unsigned int constraint_start_index, bool* is_unfixed)
{
	for (unsigned int i = 1; i < 4; ++i) {
		if (is_unfixed[i]) {
			for (unsigned int j = 0; j < i; ++j) {
				if (is_unfixed[j]) {
					address->emplace_back(&Hessian.coeffRef(start_index_in_system[j], start_index_in_system[i]));
					address->emplace_back(address->back() + 1);
					address->emplace_back(address->back() + 1);

					address->emplace_back(&Hessian.coeffRef(start_index_in_system[j], start_index_in_system[i] + 1));
					address->emplace_back(address->back() + 1);
					address->emplace_back(address->back() + 1);

					address->emplace_back(&Hessian.coeffRef(start_index_in_system[j], start_index_in_system[i] + 2));
					address->emplace_back(address->back() + 1);
					address->emplace_back(address->back() + 1);
				}
			}
		}
	}

	for (unsigned int i = 0; i < 3; ++i) {
		if (is_unfixed[i]) {
			for (unsigned int j = i + 1; j < 4; ++j) {
				if (is_unfixed[j]) {
					address->emplace_back(&Hessian.coeffRef(start_index_in_system[j], start_index_in_system[i]));
					address->emplace_back(address->back() + 1);
					address->emplace_back(address->back() + 1);

					address->emplace_back(&Hessian.coeffRef(start_index_in_system[j], start_index_in_system[i] + 1));
					address->emplace_back(address->back() + 1);
					address->emplace_back(address->back() + 1);

					address->emplace_back(&Hessian.coeffRef(start_index_in_system[j], start_index_in_system[i] + 2));
					address->emplace_back(address->back() + 1);
					address->emplace_back(address->back() + 1);
				}
			}
		}
	}

	//record grad_c
	for (unsigned int i = 0; i < 4; ++i) {
		if (is_unfixed[i]) {
			address->emplace_back(&Hessian.coeffRef(start_index_in_system[i], constraint_start_index));
			address->emplace_back(address->back() + 1);
			address->emplace_back(address->back() + 1);
		}
	}

	for (unsigned int i = 0; i < 4; ++i) {
		if (is_unfixed[i]) {
			address->emplace_back(&Hessian.coeffRef(constraint_start_index, start_index_in_system[i]));
			address->emplace_back(&Hessian.coeffRef(constraint_start_index, start_index_in_system[i] + 1));
			address->emplace_back(&Hessian.coeffRef(constraint_start_index, start_index_in_system[i] + 2));
		}
	}

	address->emplace_back(&Hessian.coeffRef(constraint_start_index, constraint_start_index));
}






void SecondOrderLargeSystem::getHessianCoeffAddress(double** address, unsigned int start_index_in_system_0, unsigned int start_index_in_system_1, double** grad_C_address,
	unsigned int constraint_start_index_0)
{
	unsigned int k = 0;
	for (unsigned int j = start_index_in_system_1; j < start_index_in_system_1 + 3; ++j) {
		address[k] = &Hessian.coeffRef(start_index_in_system_0, j);
		address[k + 1] = address[k] + 1;
		address[k + 2] = address[k] + 2;
		k += 3;
	}
	for (unsigned int j = start_index_in_system_0; j < start_index_in_system_0 + 3; ++j) {
		address[k] = &Hessian.coeffRef(start_index_in_system_1, j);
		address[k + 1] = address[k] + 1;
		address[k + 2] = address[k] + 2;
		k += 3;
	}

	grad_C_address[0] = &Hessian.coeffRef(start_index_in_system_0, constraint_start_index_0);
	grad_C_address[1] = grad_C_address[0]+1;
	grad_C_address[2] = grad_C_address[1]+1;

	grad_C_address[3]= &Hessian.coeffRef(start_index_in_system_1, constraint_start_index_0);
	grad_C_address[4] = grad_C_address[3] + 1;
	grad_C_address[5] = grad_C_address[4] + 1;

	grad_C_address[6] = &Hessian.coeffRef(constraint_start_index_0,start_index_in_system_0);
	grad_C_address[7] = &Hessian.coeffRef(constraint_start_index_0, start_index_in_system_0+1);
	grad_C_address[8] = &Hessian.coeffRef(constraint_start_index_0, start_index_in_system_0 + 2);

	grad_C_address[9] = &Hessian.coeffRef(constraint_start_index_0, start_index_in_system_1);
	grad_C_address[10] = &Hessian.coeffRef(constraint_start_index_0, start_index_in_system_1+1);
	grad_C_address[11] = &Hessian.coeffRef(constraint_start_index_0, start_index_in_system_1+2);

	grad_C_address[12] = &Hessian.coeffRef(constraint_start_index_0, constraint_start_index_0);
}



void SecondOrderLargeSystem::updateInternalForce(double* vertex_position_0, double* vertex_position_1, double* force_0,
	double* force_1,  double rest_length, double* h, double lambda, double alpha, double* ori_0, double* ori_1, double beta)
{
	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double coe = sqrt(DOT(Ax, Ax));
	double C = coe - rest_length;
	coe = 1.0 / coe;
	MULTI_(Ax, coe);

	Vector3d ori_vec;
	for (unsigned int i = 0; i < 3; ++i) {
		ori_vec[i] = vertex_position_0[i] - ori_0[i] - vertex_position_1[i] + ori_1[i];
	}

	*h = C + lambda * alpha+alpha*beta*DOT(Ax,ori_vec);

	MULTI_(Ax, lambda);
	SUM_(force_0, Ax);
	SUB_(force_1, Ax);




}

void SecondOrderLargeSystem::updateInternalForceOnlyOneEdgeFixed(double* vertex_position_0, double* vertex_position_1, double* force_0,
	double alpha, double rest_length, double* h, double lambda, double* ori_0, double* ori_1, double beta)
{

	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double coe = sqrt(DOT(Ax, Ax));
	double C = coe - rest_length;
	coe = 1.0 / coe;
	MULTI_(Ax, coe);

	Vector3d ori_vec;
	for (unsigned int i = 0; i < 3; ++i) {
		ori_vec[i] = vertex_position_0[i] - ori_0[i] - vertex_position_1[i] + ori_1[i];
	}
	*h = C + lambda * alpha + alpha * beta * DOT(Ax, ori_vec);

	MULTI_(Ax, lambda);
	SUM_(force_0, Ax);
	


}

void SecondOrderLargeSystem::computeHessianFixedOneVertex(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
	double alpha, double rest_length, unsigned int start_index_in_system_0,
	Triplet<double>* hessian_nnz_grad_C, double lambda, unsigned int constraint_start_in_sys)
{
	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length = sqrt(DOT(Ax, Ax));
	DEV_(Ax, length);

	double coe = lambda / length;

	//matrix store in row_major
	double matrix[9] = { coe * (1.0 - Ax[0] * Ax[0]), 	-coe * Ax[0] * Ax[1], 	-coe * (Ax[2] * Ax[0]),
		-coe * Ax[1] * Ax[0], coe * (1.0 - Ax[1] * Ax[1]), 	-coe * Ax[2] * Ax[1],
		-coe * Ax[2] * Ax[0], -coe * Ax[2] * Ax[1], coe * (1.0 - Ax[2] * Ax[2]) };

	for (unsigned int i = 0; i < 9; ++i) {
		diagonal_coeff_0[i] -= matrix[i];
	}

	for (unsigned int i = 0; i < 3; ++i) {
		hessian_nnz_grad_C[i] = Triplet<double>(start_index_in_system_0 + i, constraint_start_in_sys, -Ax[i]);
		hessian_nnz_grad_C[i + 3] = Triplet<double>(constraint_start_in_sys, start_index_in_system_0 + i, Ax[i]);
	}
	hessian_nnz_grad_C[6] = Triplet<double>(constraint_start_in_sys, constraint_start_in_sys, alpha);
}

void SecondOrderLargeSystem::computeHessianOnlyOneVertexFixedEdge(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
	double alpha, double** grad_C_address, double lambda, double* ori_0, double* ori_1, double beta)
{
	Vector3d Ax;
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length = sqrt(DOT(Ax, Ax));
	DEV_(Ax, length);

	double coe = 1.0 / length;

	Matrix3d matrix;
	matrix << coe * (1.0 - Ax[0] * Ax[0]), 	-coe * Ax[0] * Ax[1], 	-coe * (Ax[2] * Ax[0]),
		-coe * Ax[1] * Ax[0], coe * (1.0 - Ax[1] * Ax[1]), 	-coe * Ax[2] * Ax[1],
		-coe * Ax[2] * Ax[0], -coe * Ax[2] * Ax[1], coe * (1.0 - Ax[2] * Ax[2]);



	Vector3d ori_vec;
	for (unsigned int i = 0; i < 3; ++i) {
		ori_vec[i] = vertex_position_0[i] - ori_0[i] - vertex_position_1[i] + ori_1[i];
	}
	Vector3d damp_vec = (1.0 + alpha * beta) * Ax + alpha * beta * matrix * ori_vec;

	matrix *= lambda;
	for (unsigned int i = 0; i < 9; ++i) {
		diagonal_coeff_0[i] -= matrix.data()[i];// +damp_derivative[i];
	}

	for (unsigned int i = 0; i < 3; ++i) {
		*(grad_C_address[i]) = -Ax[i];
		*(grad_C_address[i + 3]) = -damp_vec[i];
	}
	*(grad_C_address[6]) = -alpha;
}




void SecondOrderLargeSystem::computeHessianFixedStructure(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
	double* diagonal_coeff_1, double** hessian_coeff_address, double alpha, double** grad_C_address,
	 double lambda, double* ori_0, double* ori_1, double beta)
{
	Vector3d Ax;
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length = sqrt(DOT(Ax, Ax));
	DEV_(Ax, length);

	double coe = 1.0 / length;

	//matrix store in row_major
	Matrix3d matrix;
	matrix << coe * (1.0 - Ax[0] * Ax[0]), -coe * Ax[0] * Ax[1], -coe * (Ax[2] * Ax[0]),
		-coe * Ax[1] * Ax[0], coe * (1.0 - Ax[1] * Ax[1]), 	-coe * Ax[2] * Ax[1],
		-coe * Ax[2] * Ax[0], -coe * Ax[2] * Ax[1], coe * (1.0 - Ax[2] * Ax[2]);

	Vector3d ori_vec;
	for (unsigned int i = 0; i < 3; ++i) {
		ori_vec[i] = vertex_position_0[i] - ori_0[i] - vertex_position_1[i] + ori_1[i];
	}
	Vector3d damp_vec = (1.0 + alpha * beta) * Ax + alpha * beta * matrix * ori_vec;

	matrix *= lambda;

	for (unsigned int i = 0; i < 9; ++i) {
		*hessian_coeff_address[i] = matrix.data()[i];// +damp_derivative[i]);
		*hessian_coeff_address[i + 9] = matrix.data()[i];// +damp_derivative[i]);
	}
	for (unsigned int i = 0; i < 9; ++i) {
		diagonal_coeff_0[i] -= matrix.data()[i];// +damp_derivative[i];
		diagonal_coeff_1[i] -= matrix.data()[i];// +damp_derivative[i];
	}



	for (unsigned int i = 0; i < 3; ++i) {
		*grad_C_address[i] = -Ax[i];
		*grad_C_address[i+3] = Ax[i];
		*grad_C_address[i+6] = -damp_vec[i];
		*grad_C_address[i+9] = damp_vec[i];
	}
	*grad_C_address[12] = -alpha;
}

void SecondOrderLargeSystem::computeARAPHessianFixedStructure(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3, 
	double* diagonal_coeff_0, double* diagonal_coeff_1, double* diagonal_coeff_2, double* diagonal_coeff_3, double** &address,
	double alpha, Matrix<double, 3, 4>& A, bool* is_unfixed, double* lambda, double* h, double * g_0, double* g_1, double* g_2, double* g_3, unsigned int index)
{
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d U, V, rotation;
	Matrix<double, 12, 12> Hessian;
	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);

	FEM::extractRotation(deformation_gradient, eigen_value, U, V, rotation);
	double C = (deformation_gradient - rotation).norm();
	Matrix<double, 12, 1> grad;

	if (C < 1e-8) {
		Hessian.setZero();
		grad.setZero();
		*h = -C - alpha * (*lambda);
	}
	else {
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (1.0 / C) * (deformation_gradient - rotation) * A;//

		memcpy(grad.data(), grad_C_transpose.data(), 96);

		Matrix<double, 9, 9> dPdF;
		FEM::getdPdF(U, V, eigen_value, dPdF);
		FEM::backpropagateElementHessian(Hessian, dPdF, A);
		Hessian *= (0.5 / C);
		Hessian -= ((1.0 / C) * grad) * grad.transpose();
		Hessian *= -*lambda;

		//Hessian.setZero();
		*h = C + alpha * (*lambda);
	}

		//if (index == 25) {
		//	std::cout <<"real "<< C-temp_record_0 << " " << alpha-temp_record_2 << " " << *lambda-temp_record_1 << std::endl;
		//	//temp_record_1 = *h;
		//}

	//record hessian
	for (unsigned int i = 1; i < 4; ++i) {
		if (is_unfixed[i]) {
			for (unsigned int j = 0; j < i; ++j) {
				if (is_unfixed[j]) {
					**address += Hessian(3 * j, 3 * i);
					address++;
					**address += Hessian(3 * j+1, 3 * i);
					address++;
					**address += Hessian(3 * j+2, 3 * i);
					address++;
					**address += Hessian(3 * j, 3 * i+1);
					address++;
					**address += Hessian(3 * j+1, 3 * i+1);
					address++;
					**address += Hessian(3 * j+2, 3 * i+1);
					address++;
					**address += Hessian(3 * j, 3 * i + 2);
					address++;
					**address += Hessian(3 * j+1, 3 * i + 2);
					address++;
					**address += Hessian(3 * j+2, 3 * i + 2);
					address++;
				}
			}
		}
	}

	for (unsigned int i = 0; i < 3; ++i) {
		if (is_unfixed[i]) {
			for (unsigned int j = i + 1; j < 4; ++j) {
				if (is_unfixed[j]) {
					**address += Hessian(3 * j, 3 * i);
					address++;
					**address += Hessian(3 * j + 1, 3 * i);
					address++;
					**address += Hessian(3 * j + 2, 3 * i);
					address++;
					**address += Hessian(3 * j, 3 * i + 1);
					address++;
					**address += Hessian(3 * j + 1, 3 * i + 1);
					address++;
					**address += Hessian(3 * j + 2, 3 * i + 1);
					address++;
					**address += Hessian(3 * j, 3 * i + 2);
					address++;
					**address += Hessian(3 * j + 1, 3 * i + 2);
					address++;
					**address += Hessian(3 * j + 2, 3 * i + 2);
					address++;
				}
			}
		}
	}
	//record grad_c
	for (unsigned int i = 0; i < 4; ++i) {
		if (is_unfixed[i]) {
			**address = -grad.data()[3 * i];
			address++;
			**address = -grad.data()[3 * i+1];
			address++;
			**address = -grad.data()[3 * i+2];
			address++;
		}
	}

	for (unsigned int i = 0; i < 4; ++i) {
		if (is_unfixed[i]) {
			**address = -grad.data()[3 * i];
			address++;
			**address = -grad.data()[3 * i+1];
			address++;
			**address = -grad.data()[3 * i+2];
			address++;
		}
	}
	**address = -alpha;
	address++;

	if (is_unfixed[0]) {
		for (unsigned int i = 0; i < 3; ++i) {//column
			for (unsigned int j = 0; j < 3; ++j) {//row
				diagonal_coeff_0[3 * i + j] += Hessian(j, i);
			}
		}
	}

	if (is_unfixed[1]) {
		for (unsigned int i = 0; i < 3; ++i) {//column
			for (unsigned int j = 0; j < 3; ++j) {//row
				diagonal_coeff_1[3 * i + j] += Hessian(3 + j, 3 + i);
			}
		}
	}

	if (is_unfixed[2]) {
		for (unsigned int i = 0; i < 3; ++i) {//column
			for (unsigned int j = 0; j < 3; ++j) {//row
				diagonal_coeff_2[3 * i + j] += Hessian(6 + j, 6 + i);
			}
		}
	}

	if (is_unfixed[3]) {
		for (unsigned int i = 0; i < 3; ++i) {//column
			for (unsigned int j = 0; j < 3; ++j) {//row
				diagonal_coeff_3[3 * i + j] += Hessian(9 + j, 9 + i);
			}
		}
	}

	for (unsigned int i = 0; i < 3; ++i) {
		g_0[i] += *lambda * grad.data()[i];
		g_1[i] += *lambda * grad.data()[i+3];
		g_2[i] += *lambda * grad.data()[i+6];
		g_3[i] += *lambda * grad.data()[i+9];
	}


}

void SecondOrderLargeSystem::checkIfSampleRight(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	int* vertex_index, unsigned int constraint_start_in_sys,
	double alpha, Matrix<double, 3, 4>& A, bool* is_unfixed, double* lambda)
{
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d U, V, rotation;
	Matrix<double, 12, 12> Hessian_;
	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);

	FEM::extractRotation(deformation_gradient, eigen_value, U, V, rotation);
	double C = (deformation_gradient - rotation).norm();
	Matrix<double, 12, 1> grad;
	double h;

	if (C < 1e-8) {
		Hessian_.setZero();
		grad.setZero();
		h = 0;
	}
	else {
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (1.0 / C) * (deformation_gradient - rotation) * A;

		memcpy(grad.data(), grad_C_transpose.data(), 96);

		//Matrix<double, 9, 9> dPdF;
		//FEM::getdPdF(U, V, eigen_value, dPdF);
		//FEM::backpropagateElementHessian_(Hessian_, dPdF, A);
		//Hessian_ *= (0.5 / C);
		//Hessian_ -= ((1.0 / C) * grad) * grad.transpose();
		//Hessian_ *= -*lambda;

		Hessian_.setZero();
		h=  -C - alpha * (*lambda);
	}

	//record grad_c
	for (unsigned int i = 0; i < 4; ++i) {
		if (is_unfixed[i]) {
			if (Hessian.coeff(vertex_index[i], constraint_start_in_sys) != -grad.data()[3 * i]) {
				std::cout << "error " << Hessian.coeff(vertex_index[i], constraint_start_in_sys) + grad.data()[3 * i] << std::endl;
			}
			if (Hessian.coeff(vertex_index[i]+1, constraint_start_in_sys) != -grad.data()[3 * i+1]) {
				std::cout << "error " << Hessian.coeff(vertex_index[i]+1, constraint_start_in_sys) + grad.data()[3 * i+1] << std::endl;
			}
			if (Hessian.coeff(vertex_index[i] + 2, constraint_start_in_sys) != -grad.data()[3 * i + 2]) {
				std::cout << "error " << Hessian.coeff(vertex_index[i] + 2, constraint_start_in_sys) + grad.data()[3 * i + 2] << std::endl;
			}
		}
	}

	for (unsigned int i = 0; i < 4; ++i) {
		if (is_unfixed[i]) {
			if (Hessian.coeff(constraint_start_in_sys, vertex_index[i]) != grad.data()[3 * i]) {
				std::cout << "error " << Hessian.coeff(constraint_start_in_sys,vertex_index[i]) - grad.data()[3 * i] << std::endl;
			}
			if (Hessian.coeff(constraint_start_in_sys, vertex_index[i]+1) != grad.data()[3 * i+1]) {
				std::cout << "error " << Hessian.coeff(constraint_start_in_sys, vertex_index[i]+1) - grad.data()[3 * i+1] << std::endl;
			}
			if (Hessian.coeff(constraint_start_in_sys, vertex_index[i] + 2) != grad.data()[3 * i + 2]) {
				std::cout << "error " << Hessian.coeff(constraint_start_in_sys, vertex_index[i] + 2) - grad.data()[3 * i + 2] << std::endl;
			}
		}
	}

	if (Hessian.coeff(constraint_start_in_sys, constraint_start_in_sys) != alpha) {
		std::cout << "error " << Hessian.coeff(constraint_start_in_sys, constraint_start_in_sys) - alpha<< std::endl;
	}

	if (b_thread[0][constraint_start_in_sys] != h) {
		std::cout << "error " << std::endl;
	}
}

void SecondOrderLargeSystem::setARAPHessianForTest(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	std::vector<Triplet<double>>* hessian_nnz,
	int* vertex_index, unsigned int constraint_start_in_sys,
	double alpha, Matrix<double, 3, 4>& A, bool* is_unfixed, double* lambda, double* g_0, double* g_1,double* g_2, double* g_3, double* h, unsigned int index,
	double* ori_0, double* ori_1, double* ori_2, double* ori_3, double beta)
{
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d S, rotation;
	Matrix<double, 12, 12> Hessian;
	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);
	FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);
	double C = (deformation_gradient - rotation).norm();
	Matrix<double, 12, 1> grad;

	//if (iteration_number == 0) {
	//}

	if (C < 1e-8) {
		Hessian.setZero();
		grad.setZero();
		*h = 0.0;
		//*lambda = -C / alpha;
		hessian_nnz->emplace_back(Triplet<double>(constraint_start_in_sys, constraint_start_in_sys, -alpha));
		return;
	}
	else {
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (1.0 / C) * (deformation_gradient - rotation) * A;// 
		memcpy(grad.data(), grad_C_transpose.data(), 96);

		Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
		FEM::getHessian(Hessian, S, rotation, Dm, A);
		Hessian *= (0.5 / C);
		Hessian -= ((1.0 / C) * grad) * grad.transpose();

		Matrix<double, 12, 1> x_dis;
		for (unsigned int i = 0; i < 3; ++i) {
			x_dis[i] = vertex_position_0[i] - ori_0[i];
			x_dis[i+3] = vertex_position_1[i] - ori_1[i];
			x_dis[i+6] = vertex_position_2[i] - ori_2[i];
			x_dis[i+9] = vertex_position_3[i] - ori_3[i];
		}

		Matrix<double, 12, 1> grad_damp = (1 + alpha * beta) * grad + alpha * beta * (Hessian * x_dis);


		Hessian *= -*lambda;
		*h =  C + alpha * (*lambda) + alpha * beta * grad.dot(x_dis);

		//std::cout << *h << std::endl;

		//if (index == 25) {
			//std::cout << C << " " << alpha << " " << *lambda << std::endl;
			//temp_record_0 = C;
			//temp_record_1 = *lambda;
			//temp_record_2 = alpha;
		//}



	////record hessian
		for (unsigned int i = 0; i < 4; ++i) {
			if (is_unfixed[i]) {
				for (unsigned int j = 0; j < 4; ++j) {
					if (is_unfixed[j]) {
						hessian_nnz->emplace_back(Triplet<double>(vertex_index[j], vertex_index[i], Hessian(3 * j, 3 * i)));
						hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 1, vertex_index[i], Hessian(3 * j + 1, 3 * i)));
						hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 2, vertex_index[i], Hessian(3 * j + 2, 3 * i)));
						hessian_nnz->emplace_back(Triplet<double>(vertex_index[j], vertex_index[i] + 1, Hessian(3 * j, 3 * i + 1)));
						hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 1, vertex_index[i] + 1, Hessian(3 * j + 1, 3 * i + 1)));
						hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 2, vertex_index[i] + 1, Hessian(3 * j + 2, 3 * i + 1)));
						hessian_nnz->emplace_back(Triplet<double>(vertex_index[j], vertex_index[i] + 2, Hessian(3 * j, 3 * i + 2)));
						hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 1, vertex_index[i] + 2, Hessian(3 * j + 1, 3 * i + 2)));
						hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 2, vertex_index[i] + 2, Hessian(3 * j + 2, 3 * i + 2)));
					}
				}
			}
		}

		//record grad_c
		for (unsigned int i = 0; i < 4; ++i) {
			if (is_unfixed[i]) {
				hessian_nnz->emplace_back(Triplet<double>(vertex_index[i], constraint_start_in_sys, -grad.data()[3 * i]));
				hessian_nnz->emplace_back(Triplet<double>(vertex_index[i] + 1, constraint_start_in_sys, -grad.data()[3 * i + 1]));
				hessian_nnz->emplace_back(Triplet<double>(vertex_index[i] + 2, constraint_start_in_sys, -grad.data()[3 * i + 2]));
			}
		}

		for (unsigned int i = 0; i < 4; ++i) {
			if (is_unfixed[i]) {
				hessian_nnz->emplace_back(Triplet<double>(constraint_start_in_sys, vertex_index[i], -grad_damp.data()[3 * i]));
				hessian_nnz->emplace_back(Triplet<double>(constraint_start_in_sys, vertex_index[i] + 1,-grad_damp.data()[3 * i + 1]));
				hessian_nnz->emplace_back(Triplet<double>(constraint_start_in_sys, vertex_index[i] + 2, -grad_damp.data()[3 * i + 2]));
			}
		}

		hessian_nnz->emplace_back(Triplet<double>(constraint_start_in_sys, constraint_start_in_sys, -alpha));


		for (unsigned int i = 0; i < 3; ++i) {
			g_0[i] += *lambda * grad.data()[i];
			g_1[i] += *lambda * grad.data()[i + 3];
			g_2[i] += *lambda * grad.data()[i + 6];
			g_3[i] += *lambda * grad.data()[i + 9];
		}
	}
}



//vertex_index already times 3
void SecondOrderLargeSystem::computeARAPHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	double* diagonal_coeff_0, double* diagonal_coeff_1, double* diagonal_coeff_2, double* diagonal_coeff_3, std::vector<Triplet<double>>* hessian_nnz,
	int* vertex_index, unsigned int constraint_start_in_sys,
	double alpha, Matrix<double, 3, 4>& A, bool* is_unfixed, double* lambda)

{
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d U, V, rotation;
	Matrix<double, 12, 12> Hessian;
	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);

	FEM::extractRotation(deformation_gradient, eigen_value, U, V, rotation);
	double C = (deformation_gradient - rotation).norm();
	Matrix<double, 12, 1> grad;

	if (C < 1e-8) {
		Hessian.setZero();
		grad.setZero();
	}
	else {
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (1.0 / C) * (deformation_gradient - rotation) * A;

		memcpy(grad.data(), grad_C_transpose.data(), 96);

		Matrix<double, 9, 9> dPdF;

		FEM::getdPdF(U, V, eigen_value, dPdF);


		FEM::backpropagateElementHessian(Hessian, dPdF, A);
		Hessian *= (0.5 / C);
		Hessian -= ((1.0 / C) * grad) * grad.transpose();
		Hessian *= -*lambda;


	}
	//record hessian
	for (unsigned int i = 1; i < 4; ++i) {
		if (is_unfixed[i]) {
			for (unsigned int j = 0; j < i; ++j) {
				if (is_unfixed[j]) {
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j], vertex_index[i], Hessian(3 * j, 3 * i)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 1, vertex_index[i], Hessian(3 * j + 1, 3 * i)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 2, vertex_index[i], Hessian(3 * j + 2, 3 * i)));

					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j], vertex_index[i] + 1, Hessian(3 * j, 3 * i + 1)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 1, vertex_index[i] + 1, Hessian(3 * j + 1, 3 * i + 1)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 2, vertex_index[i] + 1, Hessian(3 * j + 2, 3 * i + 1)));

					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j], vertex_index[i] + 2, Hessian(3 * j, 3 * i + 2)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 1, vertex_index[i] + 2, Hessian(3 * j + 1, 3 * i + 2)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 2, vertex_index[i] + 2, Hessian(3 * j + 2, 3 * i + 2)));
				}
			}
		}
	}

	for (unsigned int i = 0; i < 3; ++i) {
		if (is_unfixed[i]) {
			for (unsigned int j = i+1; j < 4; ++j) {
				if (is_unfixed[j]) {
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j], vertex_index[i], Hessian(3 * j, 3 * i)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 1, vertex_index[i], Hessian(3 * j + 1, 3 * i)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 2, vertex_index[i], Hessian(3 * j + 2, 3 * i)));

					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j], vertex_index[i] + 1, Hessian(3 * j, 3 * i + 1)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 1, vertex_index[i] + 1, Hessian(3 * j + 1, 3 * i + 1)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 2, vertex_index[i] + 1, Hessian(3 * j + 2, 3 * i + 1)));

					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j], vertex_index[i] + 2, Hessian(3 * j, 3 * i + 2)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 1, vertex_index[i] + 2, Hessian(3 * j + 1, 3 * i + 2)));
					hessian_nnz->emplace_back( Triplet<double>(vertex_index[j] + 2, vertex_index[i] + 2, Hessian(3 * j + 2, 3 * i + 2)));
				}
			}
		}
	}

	//record grad_c
	for (unsigned int i = 0; i < 4; ++i) {
		if (is_unfixed[i]) {
			hessian_nnz->emplace_back( Triplet<double>(vertex_index[i], constraint_start_in_sys, -grad.data()[3 * i]));
			hessian_nnz->emplace_back( Triplet<double>(vertex_index[i] + 1, constraint_start_in_sys, -grad.data()[3 * i + 1]));
			hessian_nnz->emplace_back( Triplet<double>(vertex_index[i] + 2, constraint_start_in_sys, -grad.data()[3 * i + 2]));
		}
	}
	
	for (unsigned int i = 0; i < 4; ++i) {
		if (is_unfixed[i]) {
			hessian_nnz->emplace_back( Triplet<double>(constraint_start_in_sys, vertex_index[i], grad.data()[3 * i]));
			hessian_nnz->emplace_back( Triplet<double>(constraint_start_in_sys, vertex_index[i] + 1, grad.data()[3 * i + 1]));
			hessian_nnz->emplace_back( Triplet<double>(constraint_start_in_sys, vertex_index[i] + 2, grad.data()[3 * i + 2]));
		}
	}

	hessian_nnz->emplace_back(Triplet<double>(constraint_start_in_sys, constraint_start_in_sys, alpha));

	if (is_unfixed[0]) {
		for (unsigned int i = 0; i < 3; ++i) {//column
			for (unsigned int j = 0; j < 3; ++j) {//row
				diagonal_coeff_0[3 * i + j] += Hessian(j, i);
			}
		}
	}

	if (is_unfixed[1]) {
		for (unsigned int i = 0; i < 3; ++i) {//column
			for (unsigned int j = 0; j < 3; ++j) {//row
				diagonal_coeff_1[3 * i + j] += Hessian(3 + j, 3 + i);
			}
		}
	}

	if (is_unfixed[2]) {
		for (unsigned int i = 0; i < 3; ++i) {//column
			for (unsigned int j = 0; j < 3; ++j) {//row
				diagonal_coeff_2[3 * i + j] += Hessian(6 + j, 6 + i);
			}
		}
	}

	if (is_unfixed[3]) {
		for (unsigned int i = 0; i < 3; ++i) {//column
			for (unsigned int j = 0; j < 3; ++j) {//row
				diagonal_coeff_3[3 * i + j] += Hessian(9 + j, 9 + i);
			}
		}
	}

}


void SecondOrderLargeSystem::computeHessian(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
	double* diagonal_coeff_1, Triplet<double>* hessian_nnz, double alpha, double rest_length, unsigned int start_index_in_system_0,
	unsigned int start_index_in_system_1, Triplet<double>* hessian_nnz_grad_C, double lambda, unsigned int constraint_start_in_sys)
{
	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length = sqrt(DOT(Ax, Ax));
	DEV_(Ax, length);

	double coe = lambda / length;

	//matrix store in row_major
	double matrix[9] = {coe * (1.0- Ax[0] * Ax[0]), 	-coe*Ax[0]*Ax[1], 	-coe * ( Ax[2] * Ax[0]),
		-coe* Ax[1] * Ax[0], coe * (1.0- Ax[1] * Ax[1]), 	-coe * Ax[2] * Ax[1],
		-coe * Ax[2] * Ax[0], -coe *  Ax[2] * Ax[1], coe * (1.0- Ax[2] * Ax[2]) };

	hessian_nnz[0] = Triplet<double>(start_index_in_system_0, start_index_in_system_1, matrix[0]);
	hessian_nnz[1] = Triplet<double>(start_index_in_system_0 + 1, start_index_in_system_1, matrix[1]);
	hessian_nnz[2] = Triplet<double>(start_index_in_system_0 + 2, start_index_in_system_1, matrix[2]);
	hessian_nnz[3] = Triplet<double>(start_index_in_system_0, start_index_in_system_1 + 1, matrix[3]);
	hessian_nnz[4] = Triplet<double>(start_index_in_system_0 + 1, start_index_in_system_1 + 1, matrix[4]);
	hessian_nnz[5] = Triplet<double>(start_index_in_system_0 + 2, start_index_in_system_1 + 1, matrix[5]);
	hessian_nnz[6] = Triplet<double>(start_index_in_system_0, start_index_in_system_1 + 2, matrix[6]);
	hessian_nnz[7] = Triplet<double>(start_index_in_system_0 + 1, start_index_in_system_1 + 2, matrix[7]);
	hessian_nnz[8] = Triplet<double>(start_index_in_system_0 + 2, start_index_in_system_1 + 2, matrix[8]);

	hessian_nnz[9] = Triplet<double>(start_index_in_system_1, start_index_in_system_0, matrix[0]);
	hessian_nnz[10] = Triplet<double>(start_index_in_system_1 + 1, start_index_in_system_0, matrix[1]);
	hessian_nnz[11] = Triplet<double>(start_index_in_system_1 + 2, start_index_in_system_0, matrix[2]);
	hessian_nnz[12] = Triplet<double>(start_index_in_system_1, start_index_in_system_0 + 1, matrix[3]);
	hessian_nnz[13] = Triplet<double>(start_index_in_system_1 + 1, start_index_in_system_0 + 1, matrix[4]);
	hessian_nnz[14] = Triplet<double>(start_index_in_system_1 + 2, start_index_in_system_0 + 1, matrix[5]);
	hessian_nnz[15] = Triplet<double>(start_index_in_system_1, start_index_in_system_0 + 2, matrix[6]);
	hessian_nnz[16] = Triplet<double>(start_index_in_system_1 + 1, start_index_in_system_0 + 2, matrix[7]);
	hessian_nnz[17] = Triplet<double>(start_index_in_system_1 + 2, start_index_in_system_0 + 2, matrix[8]);

	for (unsigned int i = 0; i < 9; ++i) {
		diagonal_coeff_0[i] -= matrix[i];
		diagonal_coeff_1[i] -= matrix[i];
	}
	for (unsigned int i = 0; i < 3; ++i) {
		hessian_nnz_grad_C[i]= Triplet<double>(start_index_in_system_0 +i, constraint_start_in_sys, -Ax[i]);
		hessian_nnz_grad_C[i+3]= Triplet<double>(start_index_in_system_1 +i, constraint_start_in_sys, Ax[i]);
		hessian_nnz_grad_C[i +6] = Triplet<double>(constraint_start_in_sys, start_index_in_system_0 + i, Ax[i]);
		hessian_nnz_grad_C[i +9] = Triplet<double>(constraint_start_in_sys, start_index_in_system_1 + i, -Ax[i]);
	}
	hessian_nnz_grad_C[12]= Triplet<double>(constraint_start_in_sys, constraint_start_in_sys,alpha);

}

void SecondOrderLargeSystem::updateIndexBeginPerObj()
{
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_begin_per_obj[i + 1] = vertex_begin_per_obj[i] + cloth->data()[i].mesh_struct.unfixed_point_index.size();
		//edge_begin_per_obj[i + 1] = edge_begin_per_obj[i] + (cloth->data()[i].mesh_struct.unfixed_edge_vertex_index.size() >> 1);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_begin_per_obj[i + 1 + cloth->size()] = vertex_begin_per_obj[i + cloth->size()] + tetrahedron->data()[i].mesh_struct.unfixed_point_index.size();
		//edge_begin_per_obj[i + 1 + cloth->size()] = edge_begin_per_obj[i + cloth->size()] + (tetrahedron->data()[i].mesh_struct.unfixed_edge_vertex_index.size() >> 1);
	}
	initialHessianNnz();
	computeGravity();
	setHessian();
}

void SecondOrderLargeSystem::reorganzieDataOfObjects()
{
	vertex_position.resize(total_obj_num);
	render_position.resize(total_obj_num);
	edge_index_begin_per_thread_for_mass_spring.resize(total_obj_num);
	only_one_vertex_fixed_edge_index_begin_per_thread.resize(total_obj_num);
	vertex_begin_per_obj.resize(total_obj_num + 1, 0);
	//edge_begin_per_obj.resize(total_obj_num+1,0);
	edge_vertices_mass_spring.resize(total_obj_num);
	only_one_vertex_fix_edge_vertices.resize(total_obj_num);
	unfixed_rest_length.resize(total_obj_num);
	fixed_one_vertices_rest_length.resize(total_obj_num);
	edge_length_stiffness.resize(total_obj_num);
	//vertex_index_begin_per_thread.resize(total_obj_num);
	mass.resize(total_obj_num);
	inv_mass.resize(total_obj_num);

	unfixed_vertex_begin_per_thread.resize(total_obj_num);
	unfixed_vertex.resize(total_obj_num);
	total_vertex_num.resize(total_obj_num);
	real_index_to_unfixed_index.resize(total_obj_num);

	anchor_stiffness.resize(total_obj_num);

	previous_frame_edge_length_stiffness.resize(total_obj_num);

	//vertex_pos_max_step.resize(total_obj_num);
	//for (unsigned int i = 0; i < cloth->size(); ++i) {
	//	vertex_pos_max_step[i]= cloth->data()[i].mesh_struct.vertex_position
	//}

	tet_indices.resize(tetrahedron->size());
	tet_mesh_struct.resize(tetrahedron->size());

	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position.data();
		render_position[i] = cloth->data()[i].mesh_struct.vertex_for_render.data();
		edge_index_begin_per_thread_for_mass_spring[i] = cloth->data()[i].mesh_struct.unfixed_edge_index_begin_per_thread.data();
		only_one_vertex_fixed_edge_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.only_one_vertex_fixed_edge_index_begin_per_thread.data();
		edge_vertices_mass_spring[i] = &cloth->data()[i].mesh_struct.unfixed_edge_vertex_index;
		only_one_vertex_fix_edge_vertices[i] = &cloth->data()[i].mesh_struct.only_one_vertex_fix_edge;
		unfixed_rest_length[i] = &cloth->data()[i].mesh_struct.unfixed_rest_edge_length;
		fixed_one_vertices_rest_length[i] = &cloth->data()[i].mesh_struct.fixed_one_vertex_rest_edge_length;

		edge_length_stiffness[i] = &cloth->data()[i].length_stiffness;
		anchor_stiffness[i] = cloth->data()[i].position_stiffness;
		//vertex_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		mass[i] = cloth->data()[i].mesh_struct.mass.data();
		inv_mass[i] = cloth->data()[i].mesh_struct.mass_inv.data();

		unfixed_vertex_begin_per_thread[i] = cloth->data()[i].mesh_struct.unfixed_vertex_index_begin_per_thread.data();


		unfixed_vertex[i] = &cloth->data()[i].mesh_struct.unfixed_point_index;
		vertex_begin_per_obj[i + 1] = vertex_begin_per_obj[i] + cloth->data()[i].mesh_struct.unfixed_point_index.size();
		//edge_begin_per_obj[i+1] = edge_begin_per_obj[i]+ (cloth->data()[i].mesh_struct.edge_vertices.size()>>1);
		total_vertex_num[i] = cloth->data()[i].mesh_struct.vertex_position.size();
		real_index_to_unfixed_index[i] = cloth->data()[i].mesh_struct.real_index_to_unfixed_index.data();

		previous_frame_edge_length_stiffness[i] = cloth->data()[i].length_stiffness;
	}

	tet_index_begin_per_thread.resize(tetrahedron->size());
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		render_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.data();
		edge_length_stiffness[i + cloth->size()] = &tetrahedron->data()[i].edge_length_stiffness;
		anchor_stiffness[i + cloth->size()] = tetrahedron->data()[i].position_stiffness;
		edge_index_begin_per_thread_for_mass_spring[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.unfixed_edge_index_begin_per_thread.data();
		only_one_vertex_fixed_edge_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.only_one_vertex_fixed_edge_index_begin_per_thread.data();
		edge_vertices_mass_spring[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.unfixed_edge_vertex_index;
		only_one_vertex_fix_edge_vertices[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.only_one_vertex_fix_edge;
		unfixed_rest_length[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.unfixed_rest_edge_length;
		fixed_one_vertices_rest_length[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.fixed_one_vertex_rest_edge_length;
		//vertex_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		inv_mass[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.mass_inv.data();
		mass[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.mass.data();
		unfixed_vertex_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.unfixed_vertex_index_begin_per_thread.data();

		unfixed_vertex[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.unfixed_point_index;
		vertex_begin_per_obj[i + 1 + cloth->size()] = vertex_begin_per_obj[i + cloth->size()] + tetrahedron->data()[i].mesh_struct.unfixed_point_index.size();
		//edge_begin_per_obj[i + 1 + cloth->size()] = edge_begin_per_obj[i + cloth->size()] + (tetrahedron->data()[i].mesh_struct.tet_edge_vertices.size()>>1);
		total_vertex_num[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.size();
		real_index_to_unfixed_index[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.real_index_to_unfixed_index.data();

		previous_frame_edge_length_stiffness[i] = tetrahedron->data()[i].edge_length_stiffness;


		tet_indices[i] = tetrahedron->data()[i].mesh_struct.indices.data();
		tet_mesh_struct[i] = &tetrahedron->data()[i].mesh_struct;

		tet_index_begin_per_thread[i] = tetrahedron->data()[i].mesh_struct.tetrahedron_index_begin_per_thread.data();
	}

	//std::cout << edge_vertices_mass_spring[0]->size() << " " << only_one_vertex_fix_edge_vertices[0]->size() << " " <<
	//	unfixed_rest_length[0]->size() << " " << fixed_one_vertices_rest_length[0]->size() << " " << unfixed_vertex[0]->size() << std::endl;

	//for (unsigned int i = 0; i < unfixed_vertex[0]->size(); ++i) {
	//	std::cout << unfixed_vertex[0]->data()[i] << std::endl;
	//}



	//for (unsigned int i = 0; i < total_thread_num+1; ++i) {
	//	std::cout << edge_index_begin_per_thread_for_mass_spring[0][i] << " " << only_one_vertex_fixed_edge_index_begin_per_thread[0][i] <<
	//		" " << unfixed_vertex_begin_per_thread[0][i] << std::endl;
	//}
}



