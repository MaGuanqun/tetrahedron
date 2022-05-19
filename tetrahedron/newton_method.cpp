#include"newton_method.h"

NewtonMethod::NewtonMethod()
{
	gravity_ = 9.8;
	total_thread_num = std::thread::hardware_concurrency();
	iteration_number = 1000;

	time_step = 1.0 / 100.0;
	perform_collision = false;
	time_step_square = time_step * time_step;
	conv_rate = time_step * 1e-5;

	max_itr_num = 100;
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
	off_diagonal_hessian_nnz_index_begin_per_thread.resize(total_thread_num+1);

	unsigned int edge_num = 0;
	off_diagonal_hessian_nnz_index_begin_per_thread[0] = 9 * vertex_begin_per_obj[total_obj_num];

	for(unsigned int i=0;i<total_thread_num;++i){
		edge_num = 0;
		for (unsigned int j = 0; j <total_obj_num; ++j) {
			edge_num += edge_index_begin_per_thread_for_mass_spring[j][i + 1] - edge_index_begin_per_thread_for_mass_spring[j][i];
		}
		off_diagonal_hessian_nnz_index_begin_per_thread[i + 1] = off_diagonal_hessian_nnz_index_begin_per_thread[i] + 18 * edge_num;
	}
	hessian_nnz.resize(off_diagonal_hessian_nnz_index_begin_per_thread[total_thread_num]);	


	hessian_coeff_diagonal.resize(total_thread_num);
	 for (unsigned int i = 0; i < total_thread_num; ++i) {
		 hessian_coeff_diagonal[i].resize(6*vertex_begin_per_obj[total_obj_num]);
	}

	 Hessian.resize(3 * vertex_begin_per_obj[total_obj_num], 3 * vertex_begin_per_obj[total_obj_num]);
	 b.resize(3 * vertex_begin_per_obj[total_obj_num]);
	 Sn.resize(3 * vertex_begin_per_obj[total_obj_num]);
	 velocity.resize(3 * vertex_begin_per_obj[total_obj_num]);
	 velocity.setZero();
	 position_of_beginning.resize(3 * vertex_begin_per_obj[total_obj_num]);
}


void NewtonMethod::computeGravity()
{
	gravity.resize(3 * vertex_begin_per_obj[total_obj_num]);
	//double gravity_accerlation[3] = { 0,0.0,gravity_};
	//double gravity_accerlation[3] = { gravity_, 0,0.0};
	double gravity_accerlation[3] = { 0.0, -gravity_, 0.0 };
	double* mass_;
	unsigned int vertex_index_start;

	//set cloth
	unsigned int size;
	for (unsigned int j = 0; j < total_obj_num; ++j) {
		mass_ = mass[j];
		vertex_index_start = vertex_begin_per_obj[j];
		size= vertex_begin_per_obj[j+1]- vertex_index_start;
		for (unsigned int k = 0; k < size; ++k) {
			gravity[3 * (vertex_index_start + k)] = gravity_accerlation[0] * mass_[k];
			gravity[3 * (vertex_index_start + k) + 1] = gravity_accerlation[1] * mass_[k];
			gravity[3 * (vertex_index_start + k) + 2] = gravity_accerlation[2] * mass_[k];
		}		
	}
	Sn = gravity;
}

void NewtonMethod::initial()
{
	velocity.setZero();

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
}



void NewtonMethod::updateHessianFixedStructure()
{
	thread->assignTask(this, UPDATE_HESSIAN_FIXED_STRUCTURE);
	thread->assignTask(this, UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE);
	thread->assignTask(this, UPDATE_INTERNAL_FORCE);
	thread->assignTask(this, SUM_B);
	thread->assignTask(this, UPDATE_ANCHOR_POINT_HESSIAN);
}


//UPDATE_INTERNAL_FORCE
void NewtonMethod::updateInternalForce(int thread_No)
{
	memset(hessian_coeff_diagonal[thread_No].data(), 0, 24 * vertex_begin_per_obj[total_obj_num]);
	double* hessian_coeff_diag = hessian_coeff_diagonal[thread_No].data();

	unsigned int vertex_begin;
	unsigned int vertex_end;
	unsigned int* edge_vertex;
	std::array<double, 3>* vertex_pos;
	double* rest_length_;
	double edge_length_stiffness_;

	for (unsigned int obj_No = 0; obj_No < total_obj_num; ++obj_No) {
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No];
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		rest_length_ = rest_length[obj_No];
		edge_length_stiffness_ = edge_length_stiffness[obj_No];

		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			updateInternalForce(vertex_pos[edge_vertex[i]].data(), vertex_pos[edge_vertex[i + 1]].data(),
				hessian_coeff_diag + 3 * edge_vertex[i], hessian_coeff_diag + 3 * edge_vertex[i + 1],
				edge_length_stiffness_, rest_length_[i >> 1]);
		}
	}
}



//UPDATE_HESSIAN_FIXED_STRUCTURE
void NewtonMethod::updateHessianFixedStructure(int thread_No)
{
	unsigned int vertex_begin;
	unsigned int vertex_end;
	unsigned int* edge_vertex;
	std::array<double, 3>* vertex_pos;

	memset(hessian_coeff_diagonal[thread_No].data(), 0, 8 * hessian_coeff_diagonal[thread_No].size());
	double* hessian_coeff_diag = hessian_coeff_diagonal[thread_No].data();

	double** hessian_nnz_ = Hessian_coeff_address.data() + off_diagonal_hessian_nnz_index_begin_per_thread[thread_No];
	unsigned int edge_no = 0;
	double* rest_length_;
	double edge_length_stiffness_;

	double time_step_square_ = time_step_square;
	unsigned int vertex_begin_in_obj;

	for (unsigned int obj_No = 0; obj_No < total_obj_num; ++obj_No) {
		vertex_begin_in_obj = vertex_begin_per_obj[obj_No];
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No];
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		rest_length_ = rest_length[obj_No];
		edge_length_stiffness_ = edge_length_stiffness[obj_No];
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			computeHessianFixedStructure(vertex_pos[edge_vertex[i]].data(), vertex_pos[edge_vertex[i + 1]].data(),
				hessian_coeff_diag + 6 * (edge_vertex[i] + vertex_begin_in_obj), hessian_coeff_diag + 6 * (edge_vertex[i + 1]+ vertex_begin_in_obj),
				hessian_nnz_ + 18 * edge_no, edge_length_stiffness_, rest_length_[i >> 1],  time_step_square_);
			edge_no++;
		}
	}

}


void NewtonMethod::setExternalForce()
{
	Sn = gravity;
}



void NewtonMethod::solveNewtonMethod()
{
	setExternalForce();
	storeInitialPosition();
	thread->assignTask(this, SET_S_N);
	iteration_number = 0;
	while (convergenceCondition())
	{
		updateRenderPosition();
		updateHessianFixedStructure();
		global_llt.factorize(Hessian);
		delta_x = global_llt.solve(b);
		thread->assignTask(this, UPDATE_POSITION_NEWTON);
		iteration_number++;
	}

	thread->assignTask(this, VELOCITY_NEWTON);	
}



void NewtonMethod::initialDHatTolerance(double ave_edge_length)
{
	if (perform_collision) {
		collision.initialDHatTolerance(ave_edge_length);
	}
}


//void NewtonMethod::addExternalClothForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int cloth_No)
//{
//	if (!coe.empty()) {
//		unsigned int vertex_index_start = vertex_begin_per_obj[cloth_No];
//		for (unsigned int i = 0; i < coe.size(); ++i) {
//			for (unsigned int j = 0; j < 3; ++j) {
//				f_ext.data()[neighbor_vertex[i] + vertex_index_start] += coe[i] * neighbor_vertex_force_direction[j];
//			}
//		}
//	}
//}

bool NewtonMethod::convergenceCondition()
{
	if (iteration_number < 1) {
		return true;
	}

	if (iteration_number > max_itr_num) {
		return false;
	}

	double a = 0.0;
	for (unsigned int i = 0; i < delta_x.size(); i++) {
		if (a < abs(delta_x.data()[i])) {
			a = abs(delta_x.data()[i]);
		}
	}

	//std::cout << conv_rate<<" "<< a << std::endl;

	if (a > conv_rate) {
		return true;
	}
	return false;
}

void NewtonMethod::updateRenderPosition()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(render_position[i][0].data(), vertex_position[i][0].data(), 24 * (vertex_begin_per_obj[i + 1] - vertex_begin_per_obj[i]));
	}
}


// SET_MASS_SPRING
void NewtonMethod::massSpring(int thread_No)
{
	unsigned int vertex_begin;
	unsigned int vertex_end;
	unsigned int* edge_vertex;
	std::array<double, 3>* vertex_pos ;

	memset(hessian_coeff_diagonal[thread_No].data(), 0, 8 * hessian_coeff_diagonal[thread_No].size());
	double* hessian_coeff_diag = hessian_coeff_diagonal[thread_No].data();

	Triplet<double>* hessian_nnz_ = hessian_nnz.data() + off_diagonal_hessian_nnz_index_begin_per_thread[thread_No];
	unsigned int edge_no=0;

	double time_step_square_ = time_step_square;

	double* rest_length_;

	double edge_length_stiffness_;
	unsigned int vertex_index_begin_in_system;

	for (unsigned int obj_No = 0; obj_No < total_obj_num; ++obj_No) {
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No];
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		rest_length_ = rest_length[obj_No];
		edge_length_stiffness_ = edge_length_stiffness[obj_No];
		vertex_index_begin_in_system = vertex_begin_per_obj[obj_No];
		for (unsigned int i = vertex_begin; i < vertex_end; i += 2) {
			computeHessian(vertex_pos[edge_vertex[i]].data(), vertex_pos[edge_vertex[i + 1]].data(),
				hessian_coeff_diag + 6 * (edge_vertex[i]+ vertex_index_begin_in_system), hessian_coeff_diag + 6 * (edge_vertex[i + 1]+ vertex_index_begin_in_system),
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

	unsigned int index;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		mass_ = mass[i];
		vertex_pos = vertex_position[i][0].data();
		index_end = (vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No + 1]) * 3;
		index_start = 3 * (vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No]);

		index = 3*vertex_index_begin_per_thread[i][thread_No];

		for (unsigned int j = index_start; j < index_end; ++j) {
			Sn.data()[j] = mass_[index/3] * (vertex_pos[index] + time_step_ * velocity.data()[j]) + time_step_square_* Sn.data()[j];
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
	unsigned int index;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i][0].data();
		index = 3 * vertex_index_begin_per_thread[i][thread_No];
		index_end = (vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No + 1]) * 3;
		index_start = 3 * (vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No]);
		for (unsigned int j = index_start; j < index_end; ++j) {
			vertex_pos[index] += delta_x[j];
			index++;
		}
	}
}


void NewtonMethod::storeInitialPosition()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(position_of_beginning.data() + 3 * vertex_begin_per_obj[i], vertex_position[i][0].data(), 24 * (vertex_begin_per_obj[i + 1] - vertex_begin_per_obj[i]));
	}
}


//VELOCITY_NEWTON
void NewtonMethod::updateVelocity(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;
	double* vertex_pos;
	unsigned int index;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i][0].data();
		index = 3 * vertex_index_begin_per_thread[i][thread_No];
		index_end = (vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No + 1]) * 3;
		index_start = 3 * (vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No]);
		for (unsigned int j = index_start; j < index_end; ++j) {
			velocity.data()[j] = vertex_pos[index] - position_of_beginning.data()[j];
			index++;
		}
	}
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

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		mass_ = mass[i];
		vertex_pos = vertex_position[i][0].data();

		index = 3 * vertex_index_begin_per_thread[i][thread_No];

		index_end = (vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No + 1])*3;
		index_start = 3 * (vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No]);
		memcpy(b.data() + index_start, hessian_coeff_diagonal[0].data() + index_start, 8 * (index_end - index_start));
		for (unsigned int j = index_start; j < index_end; ++j) {
			for (unsigned int k = 1; k < total_thread_num_; ++k) {
				b.data()[j] += hessian_coeff_diagonal[k][j];
			}			
			b.data()[j] = time_step_square_ * b.data()[j] + Sn.data()[j] - mass_[index / 3] * vertex_pos[index];
			index++;
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
		edge_vertex = edge_vertices_mass_spring[obj_No];
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
		index_end = vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No + 1];
		for (unsigned int j = vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			hessian_index_start = 9 * j;
			global_matrix_row_start = 3 * j;
			for (unsigned int k = 0; k < 3; ++k) {
				Hessian_coeff_address[hessian_index_start] = &Hessian.coeffRef(global_matrix_row_start, global_matrix_row_start+k);
				Hessian_coeff_address[hessian_index_start + 1] = Hessian_coeff_address[hessian_index_start] + 1;
				Hessian_coeff_address[hessian_index_start + 2] = Hessian_coeff_address[hessian_index_start] + 2;
				hessian_index_start += 3;
			}
		}
	}

}

//UPDATE_ANCHOR_POINT_HESSIAN
void NewtonMethod::updateHessianForFixPoint(int thread_No)
{
	unsigned int index_end;
	unsigned int hessian_index_start;
	int* anchor_vertex_;
	std::array<double,3>* anchor_pos_;
	std::array<double,3>* vertex_pos;
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

			force_start_index = 3* (vertex_begin_per_obj[i] + anchor_vertex_[j]);

			temp[0] = time_step_square_ * anchor_stiffness_ * (vertex_pos[anchor_vertex_[j]][0] - anchor_pos_[j][0]);
			temp[1] = time_step_square_ * anchor_stiffness_ * (vertex_pos[anchor_vertex_[j]][1] - anchor_pos_[j][1]);
			temp[2] = time_step_square_ * anchor_stiffness_ * (vertex_pos[anchor_vertex_[j]][2] - anchor_pos_[j][2]);

			b[force_start_index] -= temp[0];
			b[force_start_index + 1] -= temp[1];
			b[force_start_index + 2] -= temp[2];
		}

	}
}

//UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE
void NewtonMethod::setHessianDiagonalFixedStructure(int thread_No)
{
	unsigned int index_end;
	unsigned int vertex_start;
	unsigned int hessian_index_start;

	double* hessian_coeff_diagonal_0 = hessian_coeff_diagonal[0].data();
	unsigned int total_thread_num_ = total_thread_num;

	double mass_;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		index_end = vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No + 1];
		for (unsigned int j = vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			vertex_start = 6 * j;
			for (unsigned int m = 0; m < 6; ++m) {
				for (unsigned int k = 1; k < total_thread_num_; ++k) {
					hessian_coeff_diagonal_0[vertex_start + m] += hessian_coeff_diagonal[k][vertex_start + m];
				}
			}
			mass_ = mass[i][j - vertex_begin_per_obj[i]];
			hessian_coeff_diagonal_0[vertex_start] += mass_;
			hessian_coeff_diagonal_0[vertex_start + 3] += mass_;
			hessian_coeff_diagonal_0[vertex_start + 5] += mass_;
			hessian_index_start = 9 * j;

			*Hessian_coeff_address[hessian_index_start] = hessian_coeff_diagonal_0[vertex_start];
			*Hessian_coeff_address[hessian_index_start + 1] = hessian_coeff_diagonal_0[vertex_start + 1];
			*Hessian_coeff_address[hessian_index_start + 2] = hessian_coeff_diagonal_0[vertex_start + 2];
			*Hessian_coeff_address[hessian_index_start + 3] = hessian_coeff_diagonal_0[vertex_start + 1];
			*Hessian_coeff_address[hessian_index_start + 4] = hessian_coeff_diagonal_0[vertex_start + 3];
			*Hessian_coeff_address[hessian_index_start + 5] = hessian_coeff_diagonal_0[vertex_start + 4];
			*Hessian_coeff_address[hessian_index_start + 6] = hessian_coeff_diagonal_0[vertex_start + 2];
			*Hessian_coeff_address[hessian_index_start + 7] = hessian_coeff_diagonal_0[vertex_start + 4];
			*Hessian_coeff_address[hessian_index_start + 8] = hessian_coeff_diagonal_0[vertex_start + 5];


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

	unsigned int global_matrix_row_start;
	double mass_;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		index_end =vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No + 1];
		for (unsigned int j = vertex_begin_per_obj[i] + vertex_index_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			vertex_start = 6 * j;
			for (unsigned int m = 0; m < 6; ++m) {
				for (unsigned int k = 1; k < total_thread_num; ++k) {
					hessian_coeff_diagonal_0[vertex_start + m] += hessian_coeff_diagonal[k][vertex_start + m];
				}
			}
			mass_ = mass[i][j - vertex_begin_per_obj[i]];
			hessian_coeff_diagonal_0[vertex_start] += mass_;
			hessian_coeff_diagonal_0[vertex_start + 3] += mass_;
			hessian_coeff_diagonal_0[vertex_start + 5] += mass_;
			hessian_index_start = 9 * j;
			global_matrix_row_start = 3 * j;
			hessian_nnz[hessian_index_start] = Triplet<double>(global_matrix_row_start, global_matrix_row_start, hessian_coeff_diagonal_0[vertex_start]);
			hessian_nnz[hessian_index_start + 1] = Triplet<double>(global_matrix_row_start + 1, global_matrix_row_start, hessian_coeff_diagonal_0[vertex_start + 1]);
			hessian_nnz[hessian_index_start + 2] = Triplet<double>(global_matrix_row_start + 2, global_matrix_row_start, hessian_coeff_diagonal_0[vertex_start + 2]);
			hessian_nnz[hessian_index_start + 3] = Triplet<double>(global_matrix_row_start, global_matrix_row_start+1, hessian_coeff_diagonal_0[vertex_start + 1]);
			hessian_nnz[hessian_index_start + 4] = Triplet<double>(global_matrix_row_start+1, global_matrix_row_start+1, hessian_coeff_diagonal_0[vertex_start + 3]);
			hessian_nnz[hessian_index_start + 5] = Triplet<double>(global_matrix_row_start+2, global_matrix_row_start+1, hessian_coeff_diagonal_0[vertex_start + 4]);
			hessian_nnz[hessian_index_start + 6] = Triplet<double>(global_matrix_row_start, global_matrix_row_start+2, hessian_coeff_diagonal_0[vertex_start + 2]);
			hessian_nnz[hessian_index_start + 7] = Triplet<double>(global_matrix_row_start+1, global_matrix_row_start+2, hessian_coeff_diagonal_0[vertex_start + 4]);
			hessian_nnz[hessian_index_start + 8] = Triplet<double>(global_matrix_row_start+2, global_matrix_row_start+2, hessian_coeff_diagonal_0[vertex_start + 5]);
		}		
	}
}

void NewtonMethod::getHessianCoeffAddress(double** address, unsigned int start_index_in_system_0, unsigned int start_index_in_system_1)
{
	unsigned int k = 0;
	for (unsigned int j = start_index_in_system_1; j < start_index_in_system_1 + 3;++j) {
		address[k] = &Hessian.coeffRef(start_index_in_system_0, j);
		address[k +1] = address[k] +1;
		address[k +2] = address[k] + 2;
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
	double* force_1, double stiffness, double rest_length)
{
	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length = sqrt(DOT(Ax, Ax));
	double coe = stiffness - stiffness * rest_length / length;
	MULTI_(Ax, coe);
	SUM_(force_1, Ax);
	SUB_(force_0, Ax);

}


void NewtonMethod::computeHessianFixedStructure(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
	double* diagonal_coeff_1, double** hessian_coeff_address, double stiffness, double rest_length, double time_step_square)
{
	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length_ = DOT(Ax, Ax);
	double length = sqrt(length_);
	double coe = stiffness - stiffness*rest_length / length;
	double coe2 = stiffness*rest_length / (length_*length);
	//matrix store in row_major
	double matrix[6] = { time_step_square * (coe + coe2 * Ax[0] * Ax[0]), time_step_square * (coe2 * Ax[1] * Ax[0]), time_step_square * (coe2 * Ax[2] * Ax[0]),
		time_step_square * (coe + coe2 * Ax[1] * Ax[1]), time_step_square * (coe2 * Ax[2] * Ax[1]), time_step_square * (coe + coe2 * Ax[2] * Ax[2]) };

	*hessian_coeff_address[0] = -matrix[0];
	*hessian_coeff_address[1] = -matrix[1];
	*hessian_coeff_address[2] = -matrix[2];
	*hessian_coeff_address[3] = -matrix[1];
	*hessian_coeff_address[4] = -matrix[3];
	*hessian_coeff_address[5] = -matrix[4];
	*hessian_coeff_address[6] = -matrix[2];
	*hessian_coeff_address[7] = -matrix[4];
	*hessian_coeff_address[8] = -matrix[5];

	*hessian_coeff_address[9] = -matrix[0];
	*hessian_coeff_address[10] = -matrix[1];
	*hessian_coeff_address[11] = -matrix[2];
	*hessian_coeff_address[12] = -matrix[1];
	*hessian_coeff_address[13] = -matrix[3];
	*hessian_coeff_address[14] = -matrix[4];
	*hessian_coeff_address[15] = -matrix[2];
	*hessian_coeff_address[16] = -matrix[4];
	*hessian_coeff_address[17] = -matrix[5];

	for (unsigned int i = 0; i < 6; ++i) {
		diagonal_coeff_0[i] += matrix[i];
		diagonal_coeff_1[i] += matrix[i];
	}

}

void NewtonMethod::computeHessian(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0, 
	double* diagonal_coeff_1, Triplet<double>* hessian_nnz, double stiffness, double rest_length, unsigned int start_index_in_system_0,
	unsigned int start_index_in_system_1, double time_step_square)
{
	double Ax[3];
	SUB(Ax, vertex_position_0, vertex_position_1);
	double length_ = DOT(Ax, Ax);
	double length = sqrt(length_);
	length_ *= length;
	double coe = stiffness - stiffness*rest_length / length;
	double coe2 = stiffness*rest_length / length_;
	//matrix store in row_major
	double matrix[6] = { time_step_square * (coe + coe2 * Ax[0] * Ax[0]), 
		time_step_square *(coe2 * Ax[1] * Ax[0]), 
		time_step_square *(coe2 * Ax[2] * Ax[0]),
		time_step_square * (coe + coe2 * Ax[1] * Ax[1]), 
		time_step_square *(coe2 * Ax[2] * Ax[1]), 
		time_step_square *(coe + coe2 * Ax[2] * Ax[2]) };

	hessian_nnz[0] = Triplet<double>(start_index_in_system_0, start_index_in_system_1, -matrix[0]);
	hessian_nnz[1] = Triplet<double>(start_index_in_system_0+1, start_index_in_system_1, -matrix[1]);
	hessian_nnz[2] = Triplet<double>(start_index_in_system_0+2, start_index_in_system_1, -matrix[2]);
	hessian_nnz[3] = Triplet<double>(start_index_in_system_0, start_index_in_system_1+1, -matrix[1]);
	hessian_nnz[4] = Triplet<double>(start_index_in_system_0+1, start_index_in_system_1+1, -matrix[3]);
	hessian_nnz[5] = Triplet<double>(start_index_in_system_0+2, start_index_in_system_1+1, -matrix[4]);
	hessian_nnz[6] = Triplet<double>(start_index_in_system_0, start_index_in_system_1+2, -matrix[2]);
	hessian_nnz[7] = Triplet<double>(start_index_in_system_0+1, start_index_in_system_1+2, -matrix[4]);
	hessian_nnz[8] = Triplet<double>(start_index_in_system_0+2, start_index_in_system_1+2, -matrix[5]);

	hessian_nnz[9] = Triplet<double>(start_index_in_system_1, start_index_in_system_0, -matrix[0]);
	hessian_nnz[10] = Triplet<double>(start_index_in_system_1 + 1, start_index_in_system_0, -matrix[1]);
	hessian_nnz[11] = Triplet<double>(start_index_in_system_1 + 2, start_index_in_system_0, -matrix[2]);
	hessian_nnz[12] = Triplet<double>(start_index_in_system_1, start_index_in_system_0 + 1, -matrix[1]);
	hessian_nnz[13] = Triplet<double>(start_index_in_system_1 + 1, start_index_in_system_0 + 1, -matrix[3]);
	hessian_nnz[14] = Triplet<double>(start_index_in_system_1 + 2, start_index_in_system_0 + 1, -matrix[4]);
	hessian_nnz[15] = Triplet<double>(start_index_in_system_1, start_index_in_system_0 + 2, -matrix[2]);
	hessian_nnz[16] = Triplet<double>(start_index_in_system_1 + 1, start_index_in_system_0 + 2, -matrix[4]);
	hessian_nnz[17] = Triplet<double>(start_index_in_system_1 + 2, start_index_in_system_0 + 2, -matrix[5]);

	for (unsigned int i = 0; i<6; ++i) {
		diagonal_coeff_0[i] += matrix[i];
		diagonal_coeff_1[i] += matrix[i];
	}

}

void NewtonMethod::reorganzieDataOfObjects()
{
	vertex_position.resize(total_obj_num);
	render_position.resize(total_obj_num);
	edge_index_begin_per_thread_for_mass_spring.resize(total_obj_num);
	vertex_begin_per_obj.resize(total_obj_num+1,0);
	edge_begin_per_obj.resize(total_obj_num+1,0);
	edge_vertices_mass_spring.resize(total_obj_num);
	rest_length.resize(total_obj_num);
	edge_length_stiffness.resize(total_obj_num);
	vertex_index_begin_per_thread.resize(total_obj_num);
	mass.resize(total_obj_num);
	anchor_vertex_begin_per_thread.resize(total_obj_num);
	anchor_vertex.resize(total_obj_num);
	anchor_stiffness.resize(total_obj_num);
	anchor_position.resize(total_obj_num);

	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position.data();
		render_position[i] = cloth->data()[i].mesh_struct.vertex_for_render.data();
		edge_index_begin_per_thread_for_mass_spring[i] = cloth->data()[i].mesh_struct.edge_index_begin_per_thread.data();
		edge_vertices_mass_spring[i] = cloth->data()[i].mesh_struct.edge_vertices.data();
		rest_length[i] = cloth->data()[i].mesh_struct.edge_length.data();
		edge_length_stiffness[i] = cloth->data()[i].length_stiffness[0];
		anchor_stiffness[i] = cloth->data()[i].position_stiffness;
		vertex_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		mass[i] = cloth->data()[i].mesh_struct.mass.data();
		anchor_vertex_begin_per_thread[i] = cloth->data()[i].mesh_struct.anchor_index_begin_per_thread.data();
		anchor_vertex[i] = &cloth->data()[i].mesh_struct.anchor_vertex;
		anchor_position[i] = &cloth->data()[i].mesh_struct.anchor_position;
		vertex_begin_per_obj[i+1] = vertex_begin_per_obj[i]+ cloth->data()[i].mesh_struct.vertex_position.size();
		edge_begin_per_obj[i+1] = edge_begin_per_obj[i]+ (cloth->data()[i].mesh_struct.edge_vertices.size()>>1);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		render_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.data();
		edge_length_stiffness[i + cloth->size()] = tetrahedron->data()[i].edge_length_stiffness;		
		anchor_stiffness[i + cloth->size()] = tetrahedron->data()[i].position_stiffness;
		edge_index_begin_per_thread_for_mass_spring[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tet_edge_index_begin_per_thread.data();
		edge_vertices_mass_spring[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tet_edge_vertices.data();
		rest_length[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tet_rest_edge_length.data();
		vertex_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		mass[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.mass.data();
		anchor_vertex_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.anchor_index_begin_per_thread.data();
		anchor_vertex[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.anchor_vertex;
		anchor_position[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.anchor_position;
		vertex_begin_per_obj[i + 1 + cloth->size()] = vertex_begin_per_obj[i + cloth->size()] + tetrahedron->data()[i].mesh_struct.vertex_position.size();
		edge_begin_per_obj[i + 1 + cloth->size()] = edge_begin_per_obj[i + cloth->size()] + (tetrahedron->data()[i].mesh_struct.tet_edge_vertices.size()>>1);
	}
}



