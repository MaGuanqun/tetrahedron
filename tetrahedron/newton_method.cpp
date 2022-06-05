#include"newton_method.h"
#include"./basic/write_txt.h"
#include"test_assemble_matrix_newton.h"

NewtonMethod::NewtonMethod()
{
	gravity_ = 9.8;
	total_thread_num = std::thread::hardware_concurrency();
	iteration_number = 1000;

	time_step = 1.0 / 100.0;
	perform_collision = false;
	time_step_square = time_step * time_step;
	conv_rate = time_step * 1e-5;

	max_itr_num = 5e4;
	//damp_coe = 0.99;
	beta = 0.25;
	gamma = 0.5;
}



void NewtonMethod::setForNewtonMethod(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
	Thread* thread, double* tolerance_ratio)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	this->thread = thread;
	total_obj_num = cloth->size() + tetrahedron->size();

	reorganzieDataOfObjects();
	energy_per_thread.resize(total_thread_num);
	off_diagonal_hessian_nnz_index_begin_per_thread.resize(total_thread_num + 1);
	hessian_coeff_diagonal.resize(total_thread_num);
	initialHessianNnz();
	//recordEdgeHessian();
	computeGravity();

	setHessian();

	if (perform_collision) {
		collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio);
	}
}





void NewtonMethod::initialHessianNnz()
{
	unsigned int edge_num = 0;
	off_diagonal_hessian_nnz_index_begin_per_thread[0] = 9 * vertex_begin_per_obj[total_obj_num];

	for (unsigned int i = 0; i < total_thread_num; ++i) {
		edge_num = 0;
		for (unsigned int j = 0; j < total_obj_num; ++j) {
			edge_num += edge_index_begin_per_thread_for_mass_spring[j][i + 1] - edge_index_begin_per_thread_for_mass_spring[j][i];
		}
		off_diagonal_hessian_nnz_index_begin_per_thread[i + 1] = off_diagonal_hessian_nnz_index_begin_per_thread[i] + 18 * edge_num;
	}
	hessian_nnz.resize(off_diagonal_hessian_nnz_index_begin_per_thread[total_thread_num]);

	for (unsigned int i = 0; i < total_thread_num; ++i) {
		hessian_coeff_diagonal[i].resize(9 * vertex_begin_per_obj[total_obj_num]);
	}

	Hessian.resize(3 * vertex_begin_per_obj[total_obj_num], 3 * vertex_begin_per_obj[total_obj_num]);
	b.resize(3 * vertex_begin_per_obj[total_obj_num]);
	Sn.resize(3 * vertex_begin_per_obj[total_obj_num]);
	f_ext.resize(3 * vertex_begin_per_obj[total_obj_num]);
	velocity.resize(3 * vertex_begin_per_obj[total_obj_num]);
	velocity.setZero();
	acceleration.resize(3 * vertex_begin_per_obj[total_obj_num]);
	acceleration.setZero();

	pos_dis.resize(3 * vertex_begin_per_obj[total_obj_num]);
	pos_dis.setZero();

	position_of_beginning.resize(3 * vertex_begin_per_obj[total_obj_num]);

	b_thread.resize(total_thread_num);
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		b_thread[i].resize(3 * vertex_begin_per_obj[total_obj_num]);
	}

}


void NewtonMethod::computeGravity()
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

void NewtonMethod::initial()
{
	velocity.setZero();
	acceleration.setZero();

}


//void NewtonMethod::recordEdgeHessian()
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
//void NewtonMethod::initial_hessian_coeff(int thread_No)
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


void NewtonMethod::setHessian()
{
	thread->assignTask(this, SET_MASS_SPRING);
	thread->assignTask(this, SET_HESSIAN_DIAGONAL);
	Hessian.setFromTriplets(hessian_nnz.begin(), hessian_nnz.end());

	Hessian_coeff_address.resize(Hessian.nonZeros());
	if (Hessian_coeff_address.size() != hessian_nnz.size()) {
		std::cout << "error, nonzeros of hessian does not equal to hessian_nnz size " << std::endl;
	}
	thread->assignTask(this, GET_COEFF_ADDRESS);

	global_llt.analyzePattern(Hessian);
	setK();
}

void NewtonMethod::setK()
{
	thread->assignTask(this, UPDATE_HESSIAN_FIXED_STRUCTURE);
	thread->assignTask(this, UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE_INITIAL_STIFFNESS);
	store_ori_value.resize(Hessian.nonZeros());
	memcpy(store_ori_value.data(), Hessian.valuePtr(), 8 * Hessian.nonZeros());

	if (is_newmark) {
		ori_stiffness_matrix_record = rayleigh_damp_stiffness[1] * Hessian;
	}
	else {
		ori_stiffness_matrix_record = (rayleigh_damp_stiffness[1] / time_step) * Hessian;
	}

	previous_frame_rayleigh_damp_stiffness_beta = rayleigh_damp_stiffness[1];
	//std::cout << ori_stiffness_matrix_record << std::endl;
}


void NewtonMethod::updateK()
{
	if (!edgeLengthStiffnessHasChanged()) {
		setK();
		for (unsigned int i = 0; i < previous_frame_edge_length_stiffness.size(); ++i) {
			previous_frame_edge_length_stiffness[i] = *edge_length_stiffness[i];
		}
	}


	if (previous_frame_rayleigh_damp_stiffness_beta != rayleigh_damp_stiffness[1]) {
		memcpy(ori_stiffness_matrix_record.valuePtr(), store_ori_value.data(), 8 * Hessian.nonZeros());
		if (is_newmark) {
			ori_stiffness_matrix_record *= rayleigh_damp_stiffness[1];
		}
		else {
			ori_stiffness_matrix_record *= (rayleigh_damp_stiffness[1] / time_step);
		}
		previous_frame_rayleigh_damp_stiffness_beta = rayleigh_damp_stiffness[1];
	}
}


bool NewtonMethod::edgeLengthStiffnessHasChanged()
{
	for (unsigned int i = 0; i < previous_frame_edge_length_stiffness.size(); ++i) {
		if (previous_frame_edge_length_stiffness[i] != *edge_length_stiffness[i]) {
			return false;
		}
	}
	return true;
}

void NewtonMethod::updateHessianFixedStructure()
{
	thread->assignTask(this, UPDATE_HESSIAN_FIXED_STRUCTURE);
	thread->assignTask(this, UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE);

	if (is_newmark) {
		double coe = gamma / (beta * time_step);
		for (unsigned int i = 0; i < Hessian.nonZeros(); ++i) {
			Hessian.valuePtr()[i] += coe * ori_stiffness_matrix_record.valuePtr()[i];
		}
	}
	else {
		for (unsigned int i = 0; i < Hessian.nonZeros(); ++i) {
			Hessian.valuePtr()[i] += ori_stiffness_matrix_record.valuePtr()[i];
		}
	}

	thread->assignTask(this, UPDATE_INTERNAL_FORCE);
	thread->assignTask(this, SUM_B);
	if (!is_newmark) {
		b += ori_stiffness_matrix_record * pos_dis;
	}

	//std::cout << Hessian << std::endl;
	//std::cout << b << std::endl;



	//thread->assignTask(this, UPDATE_DAMP);
	//thread->assignTask(this, UPDATE_ANCHOR_POINT_HESSIAN);
}



void NewtonMethod::addExternalForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int obj_No)
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






void NewtonMethod::resetExternalForce()
{
	f_ext = gravity;
}


void  NewtonMethod::solveNewtonMethod()
{
	updateK();
	if (is_newmark) {
		newmarkBetaMethod();
	}
	else {
		solveNewtonMethod_();
	}
}


void NewtonMethod::newmarkBetaMethod()
{
	storeInitialPosition();
	iteration_number = 0;
	store_residual.clear();
	while (convergenceCondition())
	{
		setBn();
		updateHessianFixedStructure();
		global_llt.factorize(Hessian);
		delta_x = global_llt.solve(b);
		displacement_coe = 1.0;
		thread->assignTask(this, UPDATE_POSITION_NEWTON);
		residual = sqrt(b.dot(b));
		store_residual.emplace_back(residual);
		iteration_number++;
	}

	thread->assignTask(this, UPDATEVELOCITY_ACCELERATION_NEWMARK);
}


void NewtonMethod::solveNewtonMethod_()
{
	storeInitialPosition();
	thread->assignTask(this, SET_S_N);
	//if (iteration_number > 1000) {
	//	system("pause");
	//}

	iteration_number = 0;

	//std::cout << "====" << std::endl;
//	computeEnergy();
	//std::cout << "energy " << total_energy << std::endl;
	previous_energy = total_energy;
	store_residual.clear();
	while (convergenceCondition())
	{
		updateHessianFixedStructure();
		//compareTwoMatrix(Hessian, vertex_position, edge_vertices_mass_spring, only_one_vertex_fix_edge_vertices,
		//	edge_length_stiffness[0][0], unfixed_rest_length, fixed_one_vertices_rest_length, unfixed_vertex, mass, time_step);
		//compareTwoForce(Sn, b, vertex_position, edge_vertices_mass_spring, only_one_vertex_fix_edge_vertices,
		//	edge_length_stiffness[0][0], unfixed_rest_length, fixed_one_vertices_rest_length, unfixed_vertex, mass, time_step);
		global_llt.factorize(Hessian);
		delta_x = global_llt.solve(b);

		//std::cout << Hessian << std::endl;
		//if (*time_stamp == 92 && iteration_number<5) {
		//	std::cout << "====" << std::endl;
		//	std::cout << delta_x.segment(0,10) << std::endl;
		//std::cout << Hessian << std::endl;
		//}

		displacement_coe = 1.0;
		thread->assignTask(this, UPDATE_POSITION_NEWTON);
		//thread->assignTask(this, VELOCITY_NEWTON);
		//computeEnergy();
		//std::cout << "energy0 " << total_energy << std::endl;

		//if (total_energy > previous_energy) {
		//	updateRenderPosition();
		//	displacement_coe *= 0.5;
		//	change_direction = false;
		//	while (total_energy > previous_energy)
		//	{				
		//		thread->assignTask(this, UPDATE_POSITION_NEWTON_FROM_ORI);
		//		computeEnergy();
		//		if (abs(displacement_coe) < 1e-4) {
		//			//std::cout << "displacement_coe too small " << displacement_coe << std::endl;
		//			//if (!change_direction) {
		//			//	displacement_coe = -1.0;
		//			//	change_direction = true;
		//			//}
		//			//else {
		//				break;
		//			//}
		//		}
		//		displacement_coe *= 0.5;
		//	}
		//}

		residual = sqrt(b.dot(b));
		store_residual.emplace_back(residual);
		iteration_number++;
		previous_energy = total_energy;
		//std::cout << "residual " << residual << std::endl;
	}
	thread->assignTask(this, VELOCITY_NEWTON);
	//thread->assignTask(this, VELOCITY_NEWTON_2);
	log_store_residual.resize(store_residual.size());
	if (iteration_number > 10000) {
		//std::cout << iteration_number << std::endl;
		//system("pause");
		//std::string txt_file_name = "residual result " + std::to_string(iteration_number);
		//WriteTxt::writeTxt(store_residual, txt_file_name);
		//for (unsigned int i = 0; i < store_residual.size(); ++i) {
		//	log_store_residual[i] = log10(store_residual[i]);
		//}
		//txt_file_name = "residual result log" + std::to_string(iteration_number);
		//WriteTxt::writeTxt(log_store_residual, txt_file_name);
	}


}




void NewtonMethod::initialDHatTolerance(double ave_edge_length)
{
	if (perform_collision) {
		collision.initialDHatTolerance(ave_edge_length);
	}
}


bool NewtonMethod::convergenceCondition()
{
	if (iteration_number < 2) {
		return true;
	}

	if (iteration_number > max_itr_num) {
		return false;
	}

	double a = 0.0;
	double b = 0.0;
	for (unsigned int i = 0; i < delta_x.size(); i++) {
		if (a < abs(delta_x.data()[i])) {
			a = abs(delta_x.data()[i]);
			//b = delta_x.data()[i];
		}
	}

	//std::cout << b << std::endl;

	if (a > conv_rate) {
		return true;
	}
	return false;
}

void NewtonMethod::updateRenderPosition()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(render_position[i][0].data(), vertex_position[i][0].data(), 24 * total_vertex_num[i]);
	}
}




//UPDATE_INTERNAL_FORCE
void NewtonMethod::updateInternalForce(int thread_No)
{
	//memset(hessian_coeff_diagonal[thread_No].data(), 0, 24 * vertex_begin_per_obj[total_obj_num]);
	//double* hessian_coeff_diag = hessian_coeff_diagonal[thread_No].data();

	unsigned int vertex_begin;
	unsigned int vertex_end;
	unsigned int* edge_vertex;
	std::array<double, 3>* vertex_pos;

	double* rest_length_;
	double edge_length_stiffness_;

	unsigned int* unfixed_vertex_index;

	unsigned int index_start;
	double damp_stiff = *damp_stiffness / time_step;

	double* f_ = b_thread[thread_No].data();

	for (unsigned int obj_No = 0; obj_No < total_obj_num; ++obj_No) {
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No]->data();
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		rest_length_ = unfixed_rest_length[obj_No]->data();
		edge_length_stiffness_ = *edge_length_stiffness[obj_No];
		unfixed_vertex_index = unfixed_vertex[obj_No]->data();
		index_start = vertex_begin_per_obj[obj_No];


		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			updateInternalForce(vertex_pos[unfixed_vertex_index[edge_vertex[i]]].data(), vertex_pos[unfixed_vertex_index[edge_vertex[i + 1]]].data(),
				f_ + 3 * (edge_vertex[i] + index_start), f_ + 3 * (edge_vertex[i + 1] + index_start),
				edge_length_stiffness_, rest_length_[i >> 1], position_of_beginning.data() + 3 * (edge_vertex[i] + index_start),
				position_of_beginning.data() + 3 * (edge_vertex[i + 1] + index_start), damp_stiff);
		}
		//only one vertex fixed in the edge
		edge_vertex = only_one_vertex_fix_edge_vertices[obj_No]->data();
		vertex_begin = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No] << 1;
		vertex_end = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No + 1] << 1;
		rest_length_ = fixed_one_vertices_rest_length[obj_No]->data();
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			updateInternalForceOnlyOneEdgeFixed(vertex_pos[unfixed_vertex_index[edge_vertex[i]]].data(), vertex_pos[edge_vertex[i + 1]].data(),
				f_ + 3 * edge_vertex[i],
				edge_length_stiffness_, rest_length_[i >> 1], position_of_beginning.data() + 3 * (edge_vertex[i] + index_start),
				vertex_pos[edge_vertex[i + 1]].data(),
				damp_stiff);
		}
	}
}


//void NewtonMethod::updateEnergy



void NewtonMethod::computeMassSpringEnergy(int thread_No)
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
	for (unsigned int obj_No = 0; obj_No < total_obj_num; ++obj_No) {
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No]->data();
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		rest_length_ = unfixed_rest_length[obj_No]->data();
		edge_length_stiffness_ = *edge_length_stiffness[obj_No];
		unfixed_vertex_index = unfixed_vertex[obj_No]->data();
		index_start = vertex_begin_per_obj[obj_No];
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			energy += computeMassSpringEnergy(vertex_pos[unfixed_vertex_index[edge_vertex[i]]].data(), vertex_pos[unfixed_vertex_index[edge_vertex[i + 1]]].data(),
				rest_length_[i >> 1], edge_length_stiffness_);
		}
	}
	energy_per_thread[thread_No] += energy;
}


double NewtonMethod::computeMassSpringEnergy(double* position_0, double* position_1, double rest_length, double stiffness)
{
	double current_length;
	current_length = sqrt((position_0[0] - position_1[0]) * (position_0[0] - position_1[0])
		+ (position_0[1] - position_1[1]) * (position_0[1] - position_1[1])
		+ (position_0[2] - position_1[2]) * (position_0[2] - position_1[2]));
	current_length -= rest_length;
	return 0.5 * stiffness * current_length * current_length;
}


//UPDATE_HESSIAN_FIXED_STRUCTURE
void NewtonMethod::updateHessianFixedStructure(int thread_No)
{
	unsigned int vertex_begin;
	unsigned int vertex_end;
	unsigned int* edge_vertex;
	std::array<double, 3>* vertex_pos;

	double* f_ = b_thread[thread_No].data();
	memset(f_, 0, 8 * b_thread[thread_No].size());

	memset(hessian_coeff_diagonal[thread_No].data(), 0, 8 * hessian_coeff_diagonal[thread_No].size());
	double* hessian_coeff_diag = hessian_coeff_diagonal[thread_No].data();

	double** hessian_nnz_ = Hessian_coeff_address.data() + off_diagonal_hessian_nnz_index_begin_per_thread[thread_No];
	unsigned int edge_no = 0;
	double* rest_length_;
	double edge_length_stiffness_;

	double time_step_square_ = time_step_square;
	unsigned int vertex_begin_in_obj;
	unsigned int* unfixed_index_to_normal_index;
	double damp_stiff = *damp_stiffness * time_step;

	for (unsigned int obj_No = 0; obj_No < total_obj_num; ++obj_No) {
		vertex_begin_in_obj = vertex_begin_per_obj[obj_No];
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No]->data();
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		rest_length_ = unfixed_rest_length[obj_No]->data();
		edge_length_stiffness_ = *edge_length_stiffness[obj_No];
		unfixed_index_to_normal_index = unfixed_vertex[obj_No]->data();
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			computeHessianFixedStructure(vertex_pos[unfixed_index_to_normal_index[edge_vertex[i]]].data(),
				vertex_pos[unfixed_index_to_normal_index[edge_vertex[i + 1]]].data(),
				hessian_coeff_diag + 9 * (edge_vertex[i] + vertex_begin_in_obj), hessian_coeff_diag + 9 * (edge_vertex[i + 1] + vertex_begin_in_obj),
				hessian_nnz_ + 18 * edge_no, edge_length_stiffness_, rest_length_[i >> 1], time_step_square_, damp_stiff,
				f_ + 3 * (edge_vertex[i] + vertex_begin_in_obj), f_ + 3 * (edge_vertex[i + 1] + vertex_begin_in_obj));
			edge_no++;
		}
		//only one vertex fixed in the edge
		vertex_begin = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No] << 1;
		vertex_end = only_one_vertex_fixed_edge_index_begin_per_thread[obj_No][thread_No + 1] << 1;
		edge_vertex = only_one_vertex_fix_edge_vertices[obj_No]->data();
		rest_length_ = fixed_one_vertices_rest_length[obj_No]->data();
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			computeHessianOnlyOneVertexFixedEdge(vertex_pos[unfixed_index_to_normal_index[edge_vertex[i]]].data(),
				vertex_pos[edge_vertex[i + 1]].data(),
				hessian_coeff_diag + 9 * (edge_vertex[i] + vertex_begin_in_obj),
				edge_length_stiffness_, rest_length_[i >> 1], time_step_square_,
				damp_stiff, f_ + 3 * (edge_vertex[i] + vertex_begin_in_obj));
		}
	}


}
// SET_MASS_SPRING
void NewtonMethod::massSpring(int thread_No)
{
	unsigned int edge_begin;
	unsigned int edge_end;
	unsigned int* edge_vertex;
	std::array<double, 3>* vertex_pos;

	memset(hessian_coeff_diagonal[thread_No].data(), 0, 8 * hessian_coeff_diagonal[thread_No].size());
	double* hessian_coeff_diag = hessian_coeff_diagonal[thread_No].data();

	Triplet<double>* hessian_nnz_ = hessian_nnz.data() + off_diagonal_hessian_nnz_index_begin_per_thread[thread_No];
	unsigned int edge_no = 0;

	double time_step_square_ = time_step_square;

	double* rest_length_;

	double edge_length_stiffness_;
	unsigned int vertex_index_begin_in_system;

	unsigned int* unfixed_index_to_normal_index;

	for (unsigned int obj_No = 0; obj_No < total_obj_num; ++obj_No) {
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No]->data();
		edge_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		edge_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		rest_length_ = unfixed_rest_length[obj_No]->data();
		edge_length_stiffness_ = *edge_length_stiffness[obj_No];
		vertex_index_begin_in_system = vertex_begin_per_obj[obj_No];
		unfixed_index_to_normal_index = unfixed_vertex[obj_No]->data();
		for (unsigned int i = edge_begin; i < edge_end; i += 2) {
			computeHessian(vertex_pos[unfixed_index_to_normal_index[edge_vertex[i]]].data(), vertex_pos[unfixed_index_to_normal_index[edge_vertex[i + 1]]].data(),
				hessian_coeff_diag + 9 * (edge_vertex[i] + vertex_index_begin_in_system), hessian_coeff_diag + 9 * (edge_vertex[i + 1] + vertex_index_begin_in_system),
				hessian_nnz_ + 18 * edge_no, edge_length_stiffness_, rest_length_[i >> 1], 3 * (vertex_index_begin_in_system + edge_vertex[i]),
				3 * (vertex_index_begin_in_system + edge_vertex[i + 1]), time_step_square_);
			edge_no++;
		}
	}
}

//SET_S_N
void NewtonMethod::setSn(int thread_No)
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

void NewtonMethod::setBn()
{
	thread->assignTask(this, SET_B_N);
	Sn += ori_stiffness_matrix_record * b;
}

//SET_B_N //for newmark
void NewtonMethod::setBn(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	double time_step_ = time_step;
	double time_step_square_ = time_step_square;
	double* vertex_pos;
	unsigned int vertex_start;
	unsigned int* unfixed_index_to_normal_index;
	unsigned int j;


	double coe0 = gamma / (beta * time_step_);
	double coe1 = 1.0 - gamma / beta;
	double coe2 = (1.0 - gamma / (2.0 * beta)) * time_step_;
	double* mass_;

	double alpha = rayleigh_damp_stiffness[0] * time_step_square_;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_start = vertex_begin_per_obj[i];
		vertex_pos = vertex_position[i][0].data();
		index_end = vertex_start + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_start + unfixed_vertex_begin_per_thread[i][thread_No];

		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		mass_ = mass[i];
		for (unsigned int l = index_start; l < index_end; ++l) {
			j = 3 * l;

			for (unsigned int k = 0; k < 3; ++k) {
				b.data()[j + k] =
					coe0 * (vertex_pos[3 * unfixed_index_to_normal_index[l - vertex_start] + k] - position_of_beginning[j + k])
					+ coe1 * velocity.data()[j + k] + coe2 * acceleration.data()[j + k];
				Sn.data()[j + k] = mass_[l - vertex_start] * alpha * b.data()[j + k];
			}
		}
	}
	//use b to store Bn temporarily. Sn store alpha M Bn
}


void NewtonMethod::computeEnergy()
{
	thread->assignTask(this, NEWTON_METHOD_ENERGY);
	total_energy = 0;
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		total_energy += energy_per_thread[i];
	}
}

//NEWTON_METHOD_ENERGY
void NewtonMethod::computeEnergy(int thread_No)
{
	energy_per_thread[thread_No] = 0;
	computeMassSpringEnergy(thread_No);
	computeInertial(thread_No);
}


void NewtonMethod::computeInertial(int thread_No)
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

		for (unsigned int l = index_start; l < index_end; ++l) {
			j = 3 * l;
			start = 3 * unfixed_index_to_normal_index[l - vertex_start];
			energy += mass_[unfixed_index_to_normal_index[l - vertex_start]] *
				((vertex_pos[start] - Sn.data()[j]) * (vertex_pos[start] - Sn.data()[j]) +
					(vertex_pos[start + 1] - Sn.data()[j + 1]) * (vertex_pos[start + 1] - Sn.data()[j + 1]) +
					(vertex_pos[start + 2] - Sn.data()[j + 2]) * (vertex_pos[start + 2] - Sn.data()[j + 2]));
		}
	}
	energy /= (2.0 * time_step_square);
	energy_per_thread[thread_No] += energy;
}


//UPDATEVELOCITY_ACCELERATION_NEWMARK
void NewtonMethod::updateVelocityAccelerationNewMark(int thread_No)
{
	double alpha_1 = 1.0 / (beta * time_step_square);
	double alpha_2 = 1.0 / (beta * time_step);
	double alpha_3 = 1.0 / (2.0 * beta) - 1.0;
	double alpha_4 = gamma / (beta * time_step);
	double alpha_5 = 1.0 - gamma / beta;
	double alpha_6 = (1.0 - gamma / (2.0 * beta)) * time_step;

	unsigned int index_end;
	unsigned int index_start;
	double* vertex_pos;
	unsigned int index;

	unsigned int* unfixed_index_to_normal_index;
	unsigned int vertex_start;

	double time_step_ = time_step;

	double ori_velocity;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i][0].data();
		index = unfixed_vertex_begin_per_thread[i][thread_No];
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		for (unsigned int j = index_start; j < index_end; ++j) {
			vertex_start = 3 * unfixed_index_to_normal_index[index];
			for (unsigned int k = 0; k < 3; ++k) {
				ori_velocity = velocity.data()[3 * j + k];
				velocity.data()[3 * j + k] = alpha_4 * (vertex_pos[vertex_start + k] - position_of_beginning.data()[3 * j + k])
					+ alpha_5 * ori_velocity + alpha_6 * acceleration.data()[3 * j + k];
				acceleration.data()[3 * j + k] = alpha_1 * (vertex_pos[vertex_start + k] - position_of_beginning.data()[3 * j + k])
					- alpha_2 * ori_velocity - alpha_3 * acceleration.data()[3 * j + k];
			}
			index++;
		}
	}
}


//UPDATE_POSITION_NEWTON
void NewtonMethod::updatePosition(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	double* vertex_pos;
	double* ori_vertex_pos;
	unsigned int index;

	unsigned int* unfixed_index_to_normal_index;

	unsigned int vertex_start;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i][0].data();
		ori_vertex_pos = render_position[i][0].data();
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
}



//UPDATE_POSITION_NEWTON_FROM_ORI
void NewtonMethod::updatePositionFromOri(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	double* vertex_pos;
	double* ori_vertex_pos;
	unsigned int index;

	unsigned int* unfixed_index_to_normal_index;

	unsigned int vertex_start;

	double coe = displacement_coe;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i][0].data();
		ori_vertex_pos = render_position[i][0].data();
		index = unfixed_vertex_begin_per_thread[i][thread_No];
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		for (unsigned int j = index_start; j < index_end; ++j) {
			vertex_start = 3 * unfixed_index_to_normal_index[index];
			for (unsigned int k = 0; k < 3; ++k) {
				vertex_pos[vertex_start + k] = ori_vertex_pos[vertex_start + k] + coe * delta_x[3 * j + k];
			}
			index++;
		}
	}
}




//VELOCITY_NEWTON
void NewtonMethod::updateVelocity(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	double* vertex_pos;
	unsigned int index;

	unsigned int* unfixed_index_to_normal_index;
	unsigned int vertex_start;

	double time_step_ = time_step;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i][0].data();
		index = unfixed_vertex_begin_per_thread[i][thread_No];
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		for (unsigned int j = index_start; j < index_end; ++j) {
			vertex_start = 3 * unfixed_index_to_normal_index[index];
			velocity.data()[3 * j] = (vertex_pos[vertex_start] - position_of_beginning.data()[3 * j]) / time_step_;
			velocity.data()[3 * j + 1] = (vertex_pos[vertex_start + 1] - position_of_beginning.data()[3 * j + 1]) / time_step_;
			velocity.data()[3 * j + 2] = (vertex_pos[vertex_start + 2] - position_of_beginning.data()[3 * j + 2]) / time_step_;
			index++;
		}
	}
}


//VELOCITY_NEWTON_2
void NewtonMethod::updateVelocity2(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	double* vertex_pos;
	unsigned int index;

	unsigned int* unfixed_index_to_normal_index;
	unsigned int vertex_start;

	double time_step_ = time_step;
	double damp_coe_ = 0.98;// *damp_stiffness;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i][0].data();
		index = unfixed_vertex_begin_per_thread[i][thread_No];
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		for (unsigned int j = index_start; j < index_end; ++j) {
			vertex_start = 3 * unfixed_index_to_normal_index[index];
			velocity.data()[3 * j] = damp_coe_ * (vertex_pos[vertex_start] - position_of_beginning.data()[3 * j]) / time_step_;
			velocity.data()[3 * j + 1] = damp_coe_ * (vertex_pos[vertex_start + 1] - position_of_beginning.data()[3 * j + 1]) / time_step_;
			velocity.data()[3 * j + 2] = damp_coe_ * (vertex_pos[vertex_start + 2] - position_of_beginning.data()[3 * j + 2]) / time_step_;
			index++;
		}
	}
}


void NewtonMethod::storeInitialPosition()
{
	unsigned int* unfixed_vertex_to_normal;
	unsigned int vertex_size;
	std::array<double, 3>* vertex_pos;
	unsigned int vertex_start;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		unfixed_vertex_to_normal = unfixed_vertex[i]->data();
		vertex_size = unfixed_vertex[i]->size();
		vertex_pos = vertex_position[i];
		vertex_start = vertex_begin_per_obj[i];
		for (unsigned int j = 0; j < vertex_size; ++j) {
			memcpy(position_of_beginning.data() + 3 * (vertex_start + j),
				vertex_pos[unfixed_vertex_to_normal[j]].data(), 24);
		}
	}
}

void NewtonMethod::computeResidual()
{

}



//SUM_B
void NewtonMethod::sumB(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	unsigned int total_thread_num_ = total_thread_num;
	double time_step_square_ = time_step_square;

	double* vertex_pos;
	double* mass_;

	unsigned int index;
	unsigned int* unfixed_index_to_normal_vertex;
	unsigned int j;
	double coe = rayleigh_damp_stiffness[0] * time_step;
	double coe1 = 1.0 / beta;
	double coe2 = time_step / beta;
	double coe3 = (1 - 2.0 * beta) / (2.0 * beta) * time_step_square_;



	for (unsigned int i = 0; i < total_obj_num; ++i) {
		mass_ = mass[i];
		vertex_pos = vertex_position[i][0].data();
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No];
		unfixed_index_to_normal_vertex = unfixed_vertex[i]->data();
		if (is_newmark) {
			for (unsigned int l = index_start; l < index_end; ++l) {
				index = unfixed_index_to_normal_vertex[l - vertex_begin_per_obj[i]];
				j = 3 * l;
				for (unsigned int m = 0; m < 3; ++m) {
					for (unsigned int k = 1; k < total_thread_num_; ++k) {
						b_thread[0][j + m] += b_thread[k][j + m];
					}
					b.data()[j + m] = (b_thread[0][j + m] + f_ext.data()[j + m]) * time_step_square_ +
						mass_[index] * (coe1 * (position_of_beginning.data()[j + m] - vertex_pos[3 * index + m]) + coe2 * velocity.data()[j + m] + coe3 * acceleration.data()[j + m])
						- Sn.data()[j + m];
				}
			}
		}
		else {
			for (unsigned int l = index_start; l < index_end; ++l) {
				index = unfixed_index_to_normal_vertex[l - vertex_begin_per_obj[i]];
				j = 3 * l;
				for (unsigned int m = 0; m < 3; ++m) {
					for (unsigned int k = 1; k < total_thread_num_; ++k) {
						b_thread[0][j + m] += b_thread[k][j + m];
					}
					b.data()[j + m] = b_thread[0][j + m] * time_step_square_ +
						mass_[index] * (Sn.data()[j + m] - vertex_pos[3 * index + m])
						+ coe * mass_[index] * (position_of_beginning[j + m] - vertex_pos[3 * index + m]);
					pos_dis.data()[j + m] = position_of_beginning[j + m] - vertex_pos[3 * index + m];
				}
			}
		}
	}
}


//GET_COEFF_ADDRESS
void NewtonMethod::hessianCoeffAddress(int thread_No)
{
	unsigned int vertex_begin;
	unsigned int vertex_end;
	unsigned int* edge_vertex;
	double** address = Hessian_coeff_address.data() + off_diagonal_hessian_nnz_index_begin_per_thread[thread_No];
	unsigned int edge_no = 0;
	unsigned int vertex_index_begin_in_system;

	//off diagonal
	for (unsigned int obj_No = 0; obj_No < total_obj_num; ++obj_No) {
		edge_vertex = edge_vertices_mass_spring[obj_No]->data();
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		vertex_index_begin_in_system = vertex_begin_per_obj[obj_No];
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			getHessianCoeffAddress(address + 18 * edge_no, 3 * (vertex_index_begin_in_system + edge_vertex[i]),
				3 * (vertex_index_begin_in_system + edge_vertex[i + 1]));
			edge_no++;
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


//UPDATE_DAMP
void NewtonMethod::updateHessianForDamp(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;

	double* vertex_pos;

	unsigned int index;
	unsigned int* unfixed_index_to_normal_vertex;
	unsigned int hessian_index_start;

	double damp_stiff = *damp_stiffness * time_step;

	unsigned int force_start_index;
	double temp[3];


	for (unsigned int i = 0; i < total_obj_num; ++i) {
		index_end = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No + 1];
		index_start = vertex_begin_per_obj[i] + unfixed_vertex_begin_per_thread[i][thread_No];

		unfixed_index_to_normal_vertex = unfixed_vertex[i]->data();

		vertex_pos = vertex_position[i][0].data();

		for (unsigned int l = index_start; l < index_end; ++l) {
			index = 3 * unfixed_index_to_normal_vertex[l - vertex_begin_per_obj[i]];

			hessian_index_start = 9 * l;
			*Hessian_coeff_address[hessian_index_start] += damp_stiff;
			*Hessian_coeff_address[hessian_index_start + 4] += damp_stiff;
			*Hessian_coeff_address[hessian_index_start + 8] += damp_stiff;

			force_start_index = 3 * l;

			temp[0] = damp_stiff * (vertex_pos[index] - position_of_beginning[force_start_index]);
			temp[1] = damp_stiff * (vertex_pos[index + 1] - position_of_beginning[force_start_index + 1]);
			temp[2] = damp_stiff * (vertex_pos[index + 2] - position_of_beginning[force_start_index + 2]);

			b[force_start_index] -= temp[0];
			b[force_start_index + 1] -= temp[1];
			b[force_start_index + 2] -= temp[2];
		}

	}
}

//UPDATE_ANCHOR_POINT_HESSIAN
void NewtonMethod::updateHessianForFixPoint(int thread_No)
{
	unsigned int index_end;
	unsigned int hessian_index_start;
	int* anchor_vertex_;
	std::array<double, 3>* anchor_pos_;
	std::array<double, 3>* vertex_pos;
	double anchor_stiffness_;

	unsigned int force_start_index;

	double temp[3];

	unsigned int vertex_index;
	double time_step_square_ = time_step_square;


	for (unsigned int i = 0; i < total_obj_num; ++i) {
		anchor_vertex_ = anchor_vertex[i]->data();
		index_end = anchor_vertex_begin_per_thread[i][thread_No + 1];
		anchor_stiffness_ = anchor_stiffness[i];
		anchor_pos_ = anchor_position[i]->data();
		vertex_pos = vertex_position[i];
		for (unsigned int j = anchor_vertex_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			//std::cout << "has anchor " << std::endl;

			hessian_index_start = 9 * (vertex_begin_per_obj[i] + anchor_vertex_[j]);
			*Hessian_coeff_address[hessian_index_start] += anchor_stiffness_;
			*Hessian_coeff_address[hessian_index_start + 4] += anchor_stiffness_;
			*Hessian_coeff_address[hessian_index_start + 8] += anchor_stiffness_;

			force_start_index = 3 * (vertex_begin_per_obj[i] + anchor_vertex_[j]);
			temp[0] = time_step_square_ * anchor_stiffness_ * (vertex_pos[anchor_vertex_[j]][0] - anchor_pos_[j][0]);
			temp[1] = time_step_square_ * anchor_stiffness_ * (vertex_pos[anchor_vertex_[j]][1] - anchor_pos_[j][1]);
			temp[2] = time_step_square_ * anchor_stiffness_ * (vertex_pos[anchor_vertex_[j]][2] - anchor_pos_[j][2]);

			b[force_start_index] -= temp[0];
			b[force_start_index + 1] -= temp[1];
			b[force_start_index + 2] -= temp[2];
		}

	}
}


//UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE_INITIAL_STIFFNESS
void NewtonMethod::setHessianDiagonalFixedStructureInitialStiffness(int thread_No)
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
void NewtonMethod::setHessianDiagonalFixedStructure(int thread_No)
{
	unsigned int index_end;
	unsigned int vertex_start;

	double* hessian_coeff_diagonal_0 = hessian_coeff_diagonal[0].data();
	unsigned int total_thread_num_ = total_thread_num;

	unsigned int* unfixed_index_to_normal_index;

	double mass_;

	double stiff_coe;
	if (is_newmark) {
		stiff_coe = 1.0 / beta + rayleigh_damp_stiffness[0] * time_step * gamma / beta;
	}
	else {
		stiff_coe = 1.0 + rayleigh_damp_stiffness[0] * time_step;
	}

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
			mass_ = stiff_coe * mass[i][unfixed_index_to_normal_index[j - vertex_begin_per_obj[i]]];

			hessian_coeff_diagonal_0[vertex_start] += mass_;
			hessian_coeff_diagonal_0[vertex_start + 4] += mass_;
			hessian_coeff_diagonal_0[vertex_start + 8] += mass_;

			for (unsigned int k = 0; k < 9; ++k) {
				*Hessian_coeff_address[vertex_start + k] = hessian_coeff_diagonal_0[vertex_start + k];
				//*Hessian_coeff_address[hessian_index_start + 1] = hessian_coeff_diagonal_0[vertex_start + 1];
				//*Hessian_coeff_address[hessian_index_start + 2] = hessian_coeff_diagonal_0[vertex_start + 2];
				//*Hessian_coeff_address[hessian_index_start + 3] = hessian_coeff_diagonal_0[vertex_start + 1];
				//*Hessian_coeff_address[hessian_index_start + 4] = hessian_coeff_diagonal_0[vertex_start + 3];
				//*Hessian_coeff_address[hessian_index_start + 5] = hessian_coeff_diagonal_0[vertex_start + 4];
				//*Hessian_coeff_address[hessian_index_start + 6] = hessian_coeff_diagonal_0[vertex_start + 2];
				//*Hessian_coeff_address[hessian_index_start + 7] = hessian_coeff_diagonal_0[vertex_start + 4];
				//*Hessian_coeff_address[hessian_index_start + 8] = hessian_coeff_diagonal_0[vertex_start + 5];
			}

		}
	}
}


//SET_HESSIAN_DIAGONAL
void NewtonMethod::setHessianDiagonal(int thread_No)
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

			/*hessian_nnz[hessian_index_start + 3] = Triplet<double>(global_matrix_row_start, global_matrix_row_start + 1, hessian_coeff_diagonal_0[vertex_start + 1]);
			hessian_nnz[hessian_index_start + 4] = Triplet<double>(global_matrix_row_start + 1, global_matrix_row_start + 1, hessian_coeff_diagonal_0[vertex_start + 3]);
			hessian_nnz[hessian_index_start + 5] = Triplet<double>(global_matrix_row_start + 2, global_matrix_row_start + 1, hessian_coeff_diagonal_0[vertex_start + 4]);
			hessian_nnz[hessian_index_start + 6] = Triplet<double>(global_matrix_row_start, global_matrix_row_start + 2, hessian_coeff_diagonal_0[vertex_start + 2]);
			hessian_nnz[hessian_index_start + 7] = Triplet<double>(global_matrix_row_start + 1, global_matrix_row_start + 2, hessian_coeff_diagonal_0[vertex_start + 4]);
			hessian_nnz[hessian_index_start + 8] = Triplet<double>(global_matrix_row_start + 2, global_matrix_row_start + 2, hessian_coeff_diagonal_0[vertex_start + 5]);*/
		}
	}
}

void NewtonMethod::getHessianCoeffAddress(double** address, unsigned int start_index_in_system_0, unsigned int start_index_in_system_1)
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
}



void NewtonMethod::updateInternalForce(double* vertex_position_0, double* vertex_position_1, double* force_0,
	double* force_1, double stiffness, double rest_length, double* ori_position_0, double* ori_position_1, double damp_stiffness)
{
	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length_ = DOT(Ax, Ax);
	double coe = stiffness - stiffness * rest_length / sqrt(length_);// +coe3 / length_;

	double Ax_[3];
	SUB(Ax_, ori_position_0, ori_position_1);
	SUB(Ax_, Ax, Ax_);
	double x0AAx;
	//x0AAx = damp_stiffness * DOT(Ax, Ax_) / length_;
	//coe += x0AAx;	
	MULTI_(Ax, coe);
	SUM_(force_1, Ax);
	SUB_(force_0, Ax);

}

void NewtonMethod::updateInternalForceOnlyOneEdgeFixed(double* vertex_position_0, double* vertex_position_1, double* force_0,
	double stiffness, double rest_length, double* ori_position_0, double* ori_position_1, double damp_stiffness)
{
	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length_ = DOT(Ax, Ax);

	//double coe3 = stiffness * damp_stiffness * damp_stiffness * (DOT(Ax, velocity_0));
	double coe = stiffness - stiffness * rest_length / sqrt(length_);// +coe3 / length_;

	double Ax_[3];
	SUB(Ax_, ori_position_0, ori_position_1);
	SUB(Ax_, Ax, Ax_);

	//double x0AAx = damp_stiffness * DOT(Ax, Ax_) / length_;
	//coe += x0AAx;
	force_0[0] -= Ax[0] * coe;
	force_0[1] -= Ax[1] * coe;
	force_0[2] -= Ax[2] * coe;


}

void NewtonMethod::computeHessianOnlyOneVertexFixedEdge(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
	double stiffness, double rest_length, double time_step_square, double damp_stiffness,
	double* damp_force_0)
{
	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length_ = DOT(Ax, Ax);
	double length = sqrt(length_);

	//double coe3 = stiffness * damp_stiffness* damp_stiffness * (DOT(Ax, velocity_0));
	double coe = stiffness - stiffness * rest_length / length;// +coe3 / length_;
	double coe2 = stiffness * rest_length / (length_ * length);// +coe3 / (length_ * length_);




	//if (length < 1e-8) {
	//	std::cout << "too small "<<length << std::endl;
	//}

	//double Ax_[3];
	//SUB(Ax_, ori_position_0, ori_position_1);



	double matrix[9] = { (coe + coe2 * Ax[0] * Ax[0]),  (coe2 * Ax[1] * Ax[0]),  (coe2 * Ax[2] * Ax[0]),
		(coe2 * Ax[1] * Ax[0]),	 (coe + coe2 * Ax[1] * Ax[1]),  (coe2 * Ax[2] * Ax[1]),
		(coe2 * Ax[2] * Ax[0]),  (coe2 * Ax[2] * Ax[1]), 	(coe + coe2 * Ax[2] * Ax[2]) };
	//double x0AAx = DOT(Ax, Ax_) / length_;
	//double coe_grad_damp = damp_stiffness * (1.0 - x0AAx);
	//x0AAx /= length_;
	//double temp_vec[3];
	//temp_vec[0] = damp_stiffness * (Ax_[0] / length_ - 2.0 * x0AAx * Ax[0]);
	//temp_vec[1] = damp_stiffness * (Ax_[1] / length_ - 2.0 * x0AAx * Ax[1]);
	//temp_vec[2] = damp_stiffness * (Ax_[2] / length_ - 2.0 * x0AAx * Ax[2]);
	//double damp_derivative[9] = {
	// coe_grad_damp - Ax[0] * temp_vec[0], -Ax[1] * temp_vec[0],-Ax[2] * temp_vec[0],
	// -Ax[0] * temp_vec[1], coe_grad_damp - Ax[1] * temp_vec[1],-Ax[2] * temp_vec[1],
	//  -Ax[0] * temp_vec[2],  -Ax[1] * temp_vec[2], coe_grad_damp - Ax[2] * temp_vec[2]
	//};

	//std::cout << "temp vec " << temp_vec[0] << " " << temp_vec[1] << " " << temp_vec[2] << std::endl;
	//std::cout << "Ax " << Ax[0] << " " << Ax[1] << " " << Ax[2] << std::endl;
	//for (unsigned int i = 0; i < 3; ++i) {
	//	std::cout << damp_derivative[i] << " " << damp_derivative[i + 3] << " " << damp_derivative[i + 6] << std::endl;
	//}
	for (unsigned int i = 0; i < 9; ++i) {
		diagonal_coeff_0[i] += time_step_square * matrix[i];// +damp_derivative[0];
	}

}




void NewtonMethod::computeHessianFixedStructure(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
	double* diagonal_coeff_1, double** hessian_coeff_address, double stiffness, double rest_length, double time_step_square,
	double damp_stiffness, double* damp_force_0, double* damp_force_1)
{
	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length_ = DOT(Ax, Ax);
	double length = sqrt(length_);

	//if (length < 1e-8) {
	//	std::cout << "too small "<<length << std::endl;
	//}



	//double coe3 = stiffness * damp_stiffness * damp_stiffness * (DOT(Ax, velocity_0) - DOT(Ax, velocity_1));

	double coe = stiffness - stiffness * rest_length / length;// +coe3 / length_;
	double coe2 = stiffness * rest_length / (length_ * length);// +coe3 / (length_ * length_);

	//matrix store in row_major
	double matrix[9] = { (coe + coe2 * Ax[0] * Ax[0]),  (coe2 * Ax[1] * Ax[0]),  (coe2 * Ax[2] * Ax[0]),
		 (coe2 * Ax[1] * Ax[0]),	 (coe + coe2 * Ax[1] * Ax[1]),  (coe2 * Ax[2] * Ax[1]),
		 (coe2 * Ax[2] * Ax[0]),  (coe2 * Ax[2] * Ax[1]), 	(coe + coe2 * Ax[2] * Ax[2]) };


	/*	double Ax_[3];
	SUB(Ax_, ori_position_0, ori_position_1);
	double x0AAx =  DOT(Ax, Ax_) / length_;
	double coe_grad_damp = damp_stiffness * (1.0 - x0AAx);
	x0AAx /= length_;
	double temp_vec[3];
	temp_vec[0] = damp_stiffness * (Ax_[0] / length_ - 2.0 * x0AAx * Ax[0]);
	temp_vec[1] = damp_stiffness * (Ax_[1] / length_ - 2.0 * x0AAx * Ax[1]);
	temp_vec[2] = damp_stiffness * (Ax_[2] / length_ - 2.0 * x0AAx * Ax[2]);
	double damp_derivative[9] = {
	 coe_grad_damp - Ax[0] * temp_vec[0], -Ax[1] * temp_vec[0],-Ax[2] * temp_vec[0],
	 -Ax[0] * temp_vec[1], coe_grad_damp - Ax[1] * temp_vec[1],-Ax[2] * temp_vec[1],
	  -Ax[0] * temp_vec[2],  -Ax[1] * temp_vec[2], coe_grad_damp - Ax[2] * temp_vec[2]
	};*/



	for (unsigned int i = 0; i < 9; ++i) {
		matrix[i] *= time_step_square;
	}

	for (unsigned int i = 0; i < 9; ++i) {
		*hessian_coeff_address[i] = -(matrix[i]);// +damp_derivative[i]);
		*hessian_coeff_address[i + 9] = -(matrix[i]);// +damp_derivative[i]);
	}
	//*hessian_coeff_address[1] = -matrix[1];
	//*hessian_coeff_address[2] = -matrix[2];
	//*hessian_coeff_address[3] = -matrix[1];
	//*hessian_coeff_address[4] = -matrix[3];
	//*hessian_coeff_address[5] = -matrix[4];
	//*hessian_coeff_address[6] = -matrix[2];
	//*hessian_coeff_address[7] = -matrix[4];
	//*hessian_coeff_address[8] = -matrix[5];

	//*hessian_coeff_address[9] = -matrix[0];
	//*hessian_coeff_address[10] = -matrix[1];
	//*hessian_coeff_address[11] = -matrix[2];
	//*hessian_coeff_address[12] = -matrix[1];
	//*hessian_coeff_address[13] = -matrix[3];
	//*hessian_coeff_address[14] = -matrix[4];
	//*hessian_coeff_address[15] = -matrix[2];
	//*hessian_coeff_address[16] = -matrix[4];
	//*hessian_coeff_address[17] = -matrix[5];

	for (unsigned int i = 0; i < 9; ++i) {
		diagonal_coeff_0[i] += matrix[i];// +damp_derivative[i];
		diagonal_coeff_1[i] += matrix[i];// +damp_derivative[i];
	}

}

void NewtonMethod::computeHessian(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
	double* diagonal_coeff_1, Triplet<double>* hessian_nnz, double stiffness, double rest_length, unsigned int start_index_in_system_0,
	unsigned int start_index_in_system_1, double time_step_square)
{
	double beta = rayleigh_damp_stiffness[1];
	double stiff_coe = time_step_square + beta * time_step;
	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length_ = DOT(Ax, Ax);
	double length = sqrt(length_);
	length_ *= length;
	double coe = stiffness - stiffness * rest_length / length;
	double coe2 = stiffness * rest_length / length_;
	//matrix store in row_major
	double matrix[9] = { stiff_coe * (coe + coe2 * Ax[0] * Ax[0]), 	stiff_coe * (coe2 * Ax[1] * Ax[0]), 	stiff_coe * (coe2 * Ax[2] * Ax[0]),
		stiff_coe * (coe2 * Ax[1] * Ax[0]), stiff_coe * (coe + coe2 * Ax[1] * Ax[1]), 	stiff_coe * (coe2 * Ax[2] * Ax[1]),
		stiff_coe * (coe2 * Ax[2] * Ax[0]),stiff_coe * (coe2 * Ax[2] * Ax[1]), 	stiff_coe * (coe + coe2 * Ax[2] * Ax[2]) };

	hessian_nnz[0] = Triplet<double>(start_index_in_system_0, start_index_in_system_1, -matrix[0]);
	hessian_nnz[1] = Triplet<double>(start_index_in_system_0 + 1, start_index_in_system_1, -matrix[1]);
	hessian_nnz[2] = Triplet<double>(start_index_in_system_0 + 2, start_index_in_system_1, -matrix[2]);
	hessian_nnz[3] = Triplet<double>(start_index_in_system_0, start_index_in_system_1 + 1, -matrix[3]);
	hessian_nnz[4] = Triplet<double>(start_index_in_system_0 + 1, start_index_in_system_1 + 1, -matrix[4]);
	hessian_nnz[5] = Triplet<double>(start_index_in_system_0 + 2, start_index_in_system_1 + 1, -matrix[5]);
	hessian_nnz[6] = Triplet<double>(start_index_in_system_0, start_index_in_system_1 + 2, -matrix[6]);
	hessian_nnz[7] = Triplet<double>(start_index_in_system_0 + 1, start_index_in_system_1 + 2, -matrix[7]);
	hessian_nnz[8] = Triplet<double>(start_index_in_system_0 + 2, start_index_in_system_1 + 2, -matrix[8]);
	//hessian_nnz[3] = Triplet<double>(start_index_in_system_0, start_index_in_system_1+1, -matrix[1]);
	//hessian_nnz[4] = Triplet<double>(start_index_in_system_0+1, start_index_in_system_1+1, -matrix[3]);
	//hessian_nnz[5] = Triplet<double>(start_index_in_system_0+2, start_index_in_system_1+1, -matrix[4]);
	//hessian_nnz[6] = Triplet<double>(start_index_in_system_0, start_index_in_system_1+2, -matrix[2]);
	//hessian_nnz[7] = Triplet<double>(start_index_in_system_0+1, start_index_in_system_1+2, -matrix[4]);
	//hessian_nnz[8] = Triplet<double>(start_index_in_system_0+2, start_index_in_system_1+2, -matrix[5]);

	hessian_nnz[9] = Triplet<double>(start_index_in_system_1, start_index_in_system_0, -matrix[0]);
	hessian_nnz[10] = Triplet<double>(start_index_in_system_1 + 1, start_index_in_system_0, -matrix[1]);
	hessian_nnz[11] = Triplet<double>(start_index_in_system_1 + 2, start_index_in_system_0, -matrix[2]);
	hessian_nnz[12] = Triplet<double>(start_index_in_system_1, start_index_in_system_0 + 1, -matrix[3]);
	hessian_nnz[13] = Triplet<double>(start_index_in_system_1 + 1, start_index_in_system_0 + 1, -matrix[4]);
	hessian_nnz[14] = Triplet<double>(start_index_in_system_1 + 2, start_index_in_system_0 + 1, -matrix[5]);
	hessian_nnz[15] = Triplet<double>(start_index_in_system_1, start_index_in_system_0 + 2, -matrix[6]);
	hessian_nnz[16] = Triplet<double>(start_index_in_system_1 + 1, start_index_in_system_0 + 2, -matrix[7]);
	hessian_nnz[17] = Triplet<double>(start_index_in_system_1 + 2, start_index_in_system_0 + 2, -matrix[8]);

	for (unsigned int i = 0; i < 9; ++i) {
		diagonal_coeff_0[i] += matrix[i];
		diagonal_coeff_1[i] += matrix[i];
	}

}

void NewtonMethod::updateIndexBeginPerObj()
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

void NewtonMethod::reorganzieDataOfObjects()
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
	anchor_vertex_begin_per_thread.resize(total_obj_num);
	unfixed_vertex_begin_per_thread.resize(total_obj_num);
	unfixed_vertex.resize(total_obj_num);
	total_vertex_num.resize(total_obj_num);
	real_index_to_unfixed_index.resize(total_obj_num);


	anchor_vertex.resize(total_obj_num);
	anchor_stiffness.resize(total_obj_num);
	anchor_position.resize(total_obj_num);


	previous_frame_edge_length_stiffness.resize(total_obj_num);

	//vertex_pos_max_step.resize(total_obj_num);
	//for (unsigned int i = 0; i < cloth->size(); ++i) {
	//	vertex_pos_max_step[i]= cloth->data()[i].mesh_struct.vertex_position
	//}


	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position.data();
		render_position[i] = cloth->data()[i].mesh_struct.vertex_for_render.data();
		edge_index_begin_per_thread_for_mass_spring[i] = cloth->data()[i].mesh_struct.unfixed_edge_index_begin_per_thread.data();
		only_one_vertex_fixed_edge_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.only_one_vertex_fixed_edge_index_begin_per_thread.data();
		edge_vertices_mass_spring[i] = &cloth->data()[i].mesh_struct.unfixed_edge_vertex_index;
		only_one_vertex_fix_edge_vertices[i] = &cloth->data()[i].mesh_struct.only_one_vertex_fix_edge;
		unfixed_rest_length[i] = &cloth->data()[i].mesh_struct.unfixed_rest_edge_length;
		fixed_one_vertices_rest_length[i] = &cloth->data()[i].mesh_struct.fixed_one_vertex_rest_edge_length;

		edge_length_stiffness[i] = cloth->data()[i].length_stiffness.data();
		anchor_stiffness[i] = cloth->data()[i].position_stiffness;
		//vertex_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		mass[i] = cloth->data()[i].mesh_struct.mass.data();
		anchor_vertex_begin_per_thread[i] = cloth->data()[i].mesh_struct.anchor_index_begin_per_thread.data();
		unfixed_vertex_begin_per_thread[i] = cloth->data()[i].mesh_struct.unfixed_vertex_index_begin_per_thread.data();

		anchor_vertex[i] = &cloth->data()[i].mesh_struct.anchor_vertex;
		anchor_position[i] = &cloth->data()[i].mesh_struct.anchor_position;
		unfixed_vertex[i] = &cloth->data()[i].mesh_struct.unfixed_point_index;
		vertex_begin_per_obj[i + 1] = vertex_begin_per_obj[i] + cloth->data()[i].mesh_struct.unfixed_point_index.size();
		//edge_begin_per_obj[i+1] = edge_begin_per_obj[i]+ (cloth->data()[i].mesh_struct.edge_vertices.size()>>1);
		total_vertex_num[i] = cloth->data()[i].mesh_struct.vertex_position.size();
		real_index_to_unfixed_index[i] = cloth->data()[i].mesh_struct.real_index_to_unfixed_index.data();

		previous_frame_edge_length_stiffness[i] = cloth->data()[i].length_stiffness[0];
	}
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
		mass[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.mass.data();
		anchor_vertex_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.anchor_index_begin_per_thread.data();
		unfixed_vertex_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.unfixed_vertex_index_begin_per_thread.data();
		anchor_vertex[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.anchor_vertex;
		anchor_position[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.anchor_position;
		unfixed_vertex[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.unfixed_point_index;
		vertex_begin_per_obj[i + 1 + cloth->size()] = vertex_begin_per_obj[i + cloth->size()] + tetrahedron->data()[i].mesh_struct.unfixed_point_index.size();
		//edge_begin_per_obj[i + 1 + cloth->size()] = edge_begin_per_obj[i + cloth->size()] + (tetrahedron->data()[i].mesh_struct.tet_edge_vertices.size()>>1);
		total_vertex_num[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.size();
		real_index_to_unfixed_index[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.real_index_to_unfixed_index.data();

		previous_frame_edge_length_stiffness[i] = tetrahedron->data()[i].edge_length_stiffness;
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



