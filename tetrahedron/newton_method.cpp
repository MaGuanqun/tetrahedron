#include"newton_method.h"
#include"./basic/write_txt.h"
#include"test_assemble_matrix_newton.h"
#include"XPBD/FEM_relate.h"


NewtonMethod::NewtonMethod()
{
	gravity_ = 9.8;

	iteration_number = 1000;

	time_step = 1.0 / 30.0;
	perform_collision = false;
	time_step_square = time_step * time_step;
	conv_rate = time_step * 1e-2;

	max_itr_num = 100;
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
	total_thread_num =thread->thread_num;
	reorganzieDataOfObjects();
	energy_per_thread.resize(total_thread_num);
	off_diagonal_hessian_nnz_index_begin_per_thread.resize(total_thread_num + 1);
	hessian_coeff_diagonal.resize(total_thread_num);
	initialHessianNnz();
	//recordEdgeHessian();
	computeGravity();

	setHessian();

	if (perform_collision) {
		collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio, NEWTON_);
	}
}





void NewtonMethod::initialHessianNnz()
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

	tet_hessian_index[0] = off_diagonal_hessian_nnz_index_begin_per_thread[total_thread_num];


	unsigned int total_tet_num = 0;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		total_tet_num += tetrahedron->data()[i].mesh_struct.indices.size();
	}

	tet_hessian_index[1] = off_diagonal_hessian_nnz_index_begin_per_thread[total_thread_num] + 108 * total_tet_num;
	hessian_nnz.reserve(tet_hessian_index[1]);
	hessian_nnz.resize(tet_hessian_index[0]);


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


	velocity_total.resize(total_obj_num);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		velocity_total[i].resize(mesh_struct[i]->vertex_position.size());
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		velocity_total[i].resize(mesh_struct[i + cloth->size()]->vertex_position.size());
	}

	//b_for_test.resize(3 * vertex_begin_per_obj[total_obj_num]);
	//Matrix_test.resize(3 * vertex_begin_per_obj[total_obj_num], 3 * vertex_begin_per_obj[total_obj_num]);
	//only_constraint.resize(3 * vertex_begin_per_obj[total_obj_num], 3 * vertex_begin_per_obj[total_obj_num]);
}





void NewtonMethod::updateVelocity()
{
	unsigned int vertex_index_start;
	unsigned int* unfixed_index_to_normal_index;
	unsigned int size;
	std::array<double, 3>* velocity_;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_index_start = vertex_begin_per_obj[i];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		size = unfixed_vertex[i]->size();
		velocity_ = velocity_total[i].data();
		for (unsigned int k = 0; k < size; ++k) {
			memcpy(velocity.data() + 3 * (vertex_index_start + k), velocity_[unfixed_index_to_normal_index[k]].data(), 24);
		}
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_index_start = vertex_begin_per_obj[i + cloth->size()];
		unfixed_index_to_normal_index = unfixed_vertex[i + cloth->size()]->data();
		size = unfixed_vertex[i + cloth->size()]->size();
		velocity_ = velocity_total[i + cloth->size()].data();
		for (unsigned int k = 0; k < size; ++k) {
			memcpy(velocity.data() + 3 * (vertex_index_start + k), velocity_[unfixed_index_to_normal_index[k]].data(), 24);
		}
	}
}



void NewtonMethod::updateTotalVelocity()
{
	unsigned int vertex_index_start;
	unsigned int* unfixed_index_to_normal_index;
	unsigned int size;
	std::array<double,3>* velocity_;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		memset(velocity_total[i][0].data(), 0, 24 * velocity_total[i].size());
		vertex_index_start = vertex_begin_per_obj[i];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		size = unfixed_vertex[i]->size();
		velocity_ = velocity_total[i].data();
		for (unsigned int k = 0; k < size; ++k) {
			memcpy(velocity_[unfixed_index_to_normal_index[k]].data(), velocity.data() + 3 * (vertex_index_start + k), 24);
		}
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		memset(velocity_total[i+cloth->size()][0].data(), 0, 24 * velocity_total[i + cloth->size()].size());
		vertex_index_start = vertex_begin_per_obj[i + cloth->size()];
		unfixed_index_to_normal_index = unfixed_vertex[i + cloth->size()]->data();
		size = unfixed_vertex[i + cloth->size()]->size();
		velocity_ = velocity_total[i+cloth->size()].data();
		for (unsigned int k = 0; k < size; ++k) {
			memcpy(velocity_[unfixed_index_to_normal_index[k]].data(), velocity.data() + 3 * (vertex_index_start + k), 24);
		}
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

void NewtonMethod::initialCoeffDiagonal()
{
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		memset(hessian_coeff_diagonal[i].data(), 0, 8 * hessian_coeff_diagonal[i].size());
	}
}



void NewtonMethod::setHessian()
{
	initialCoeffDiagonal();
	thread->assignTask(this, SET_MASS_SPRING);
	setARAP();
	thread->assignTask(this, SET_HESSIAN_DIAGONAL);
	Hessian.setFromTriplets(hessian_nnz.begin(), hessian_nnz.end());

	Hessian_coeff_address.reserve(hessian_nnz.size());
	Hessian_coeff_address.resize(tet_hessian_index[0]);
	//if (Hessian_coeff_address.size() != hessian_nnz.size()) {
	//	std::cout << "error, nonzeros of hessian does not equal to hessian_nnz size " << std::endl;
	//}
	thread->assignTask(this, GET_COEFF_ADDRESS);
	getARAPCoeffAddress();
	global_llt.analyzePattern(Hessian);
	setK();
}



void NewtonMethod::setARAP()
{
	std::array<double, 3>* vertex_pos;
	double* hessian_coeff_diag = hessian_coeff_diagonal[0].data();
	std::vector<Triplet<double>>* hessian_nnz_ = &hessian_nnz;
	double stiffness;
	unsigned int size;
	double* inv_mass_;

	bool is_unfixed[4];
	std::array<int, 4>* indices;
	unsigned int vertex_index_start;

	int vertex_position_in_system[4];
	unsigned int* real_index_to_system_index;

	double use_for_temp[9];
	double* diagonal_coeff[4];

	double* volume;
	Matrix<double, 3, 4>* A;
	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No + cloth->size()];
		stiffness = tetrahedron->data()[obj_No].ARAP_stiffness;
		size = tet_mesh_struct[obj_No]->indices.size();
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
				}
				else {
					diagonal_coeff[j] = use_for_temp;
				}
			}
			if (is_unfixed[0] || is_unfixed[1] || is_unfixed[2] || is_unfixed[3]) {
				computeARAPHessian(vertex_pos[indices[i][0]].data(), vertex_pos[indices[i][1]].data(), vertex_pos[indices[i][2]].data(), vertex_pos[indices[i][3]].data(),
					diagonal_coeff[0], diagonal_coeff[1], diagonal_coeff[2], diagonal_coeff[3], hessian_nnz_, vertex_position_in_system,
					stiffness* volume[i], A[i], is_unfixed);
			}
		}
	}


}


void NewtonMethod::computeARAPHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	double* diagonal_coeff_0, double* diagonal_coeff_1, double* diagonal_coeff_2, double* diagonal_coeff_3, std::vector<Triplet<double>>* hessian_nnz,
	int* vertex_index, 	double stiffness, Matrix<double, 3, 4>& A, bool* is_unfixed)
{
	Vector3d eigen_value;
	Matrix3d deformation_gradient;
	Matrix3d S, rotation;
	Matrix<double, 12, 12> Hessian;
	Matrix<double, 12, 1> grad;

	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);

	FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);
	

	Matrix<double, 3, 4> grad_C_transpose;
	grad_C_transpose =(2.0* stiffness) * (deformation_gradient - rotation) * A;
	memcpy(grad.data(), grad_C_transpose.data(), 96);

	Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
	FEM::getHessian(Hessian, S, rotation, Dm, A);

	//Matrix4d result = 2.0 * A.transpose() * A;
	//Hessian.setZero();
	//for (unsigned int i = 0; i < 4; ++i) {
	//	for (unsigned int j = 0; j < 4; ++j) {
	//		Hessian(3 * i, 3 * j) = result(i, j);
	//		Hessian(3 * i+1, 3 * j+1) = result(i, j);
	//		Hessian(3 * i+2, 3 * j+2) = result(i, j);
	//	}
	//}

	Hessian *= stiffness;


	for (unsigned int i = 1; i < 4; ++i) {
		if (is_unfixed[i]) {
			for (int j = 0; j < i; ++j) {
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

	for (unsigned int i = 0; i < 3; ++i) {
		if (is_unfixed[i]) {
			for (unsigned int j = i + 1; j < 4; ++j) {
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

void NewtonMethod::setK()
{
	memset(Hessian.valuePtr(), 0, 8 * Hessian.nonZeros());
	thread->assignTask(this, UPDATE_HESSIAN_FIXED_STRUCTURE);
	updateARAPHessianFixedStructure();
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

	memset(Hessian.valuePtr(), 0, 8 * Hessian.nonZeros());
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		b_thread[i].setZero();
	}


	thread->assignTask(this, UPDATE_HESSIAN_FIXED_STRUCTURE);

	updateARAPHessianFixedStructure();



	thread->assignTask(this, UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE);

	//std::cout << only_constraint << std::endl;
	//std::cout << "/////" << std::endl;
	//std::cout << Hessian << std::endl;
	//std::cout << "error " << (only_constraint - Hessian).squaredNorm() << std::endl;
	//for (unsigned int i = 0; i < Hessian.nonZeros(); ++i) {
	//	if (Hessian.valuePtr()[i] != only_constraint.valuePtr()[i]) {
	//		std::cout << i << " " << Hessian.innerIndexPtr()[i] << " " << Hessian.valuePtr()[i] << " " << only_constraint.valuePtr()[i] << std::endl;
	//	}
	//}


	///temp comment damp
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
	////////////

	thread->assignTask(this, UPDATE_INTERNAL_FORCE);
	thread->assignTask(this, SUM_B);

	//temp comment damp
	if (!is_newmark) {
		b += ori_stiffness_matrix_record * pos_dis;
	}
	////////////

	//std::cout << Hessian << std::endl;
	//std::cout << "/////" << std::endl;
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





//void NewtonMethod::construct_b_Test()
//{
//	double* mass_;
//	unsigned int vertex_position_in_system;
//	unsigned int read_index;
//	std::array<double, 3>* vetex_pos;
//	unsigned int* unfixed;
//	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
//		mass_ = mass[obj_No];
//		unfixed = tetrahedron->data()[obj_No].mesh_struct.unfixed_point_index.data();
//		vetex_pos = vertex_position[obj_No + cloth->size()];
//		for (unsigned int i = 0; i < tetrahedron->data()[obj_No].mesh_struct.unfixed_point_index.size(); ++i) {
//			vertex_position_in_system = 3 * (vertex_begin_per_obj[cloth->size() + obj_No] + i);
//			read_index = unfixed[i];
//			b_for_test.data()[vertex_position_in_system] -= mass_[read_index] / time_step_square * (vetex_pos[read_index][0] - Sn.data()[vertex_position_in_system]);
//			b_for_test.data()[vertex_position_in_system+1] -= mass_[read_index] / time_step_square * (vetex_pos[read_index][1] - Sn.data()[vertex_position_in_system+1]);
//			b_for_test.data()[vertex_position_in_system+2] -= mass_[read_index] / time_step_square * (vetex_pos[read_index][2] - Sn.data()[vertex_position_in_system+2]);
//		}
//	}
//}


//void NewtonMethod::testSetHessian()
//{
//	std::array<double, 3>* vertex_pos;
//	double stiffness;
//	unsigned int size;
//	double* inv_mass_;
//
//	bool is_unfixed[4];
//
//	std::array<int, 4>* indices;
//	unsigned int vertex_index_start;
//
//	int vertex_position_in_system[4];
//	unsigned int* real_index_to_system_index;
//
//	double use_for_temp[9];
//	double* diagonal_coeff[4];
//	double* g_tet[4];
//	double* volume;
//	Matrix<double, 3, 4>* A;
//	double* g_;
//	test_nnz.clear();
//	double* mass_;
//	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
//		vertex_pos = vertex_position[obj_No + cloth->size()];
//		stiffness = tetrahedron->data()[obj_No].ARAP_stiffness;
//		size = tet_mesh_struct[obj_No]->indices.size();
//		inv_mass_ = inv_mass[cloth->size() + obj_No];
//		indices = tet_indices[obj_No];
//		vertex_index_start = vertex_begin_per_obj[cloth->size() + obj_No];
//		real_index_to_system_index = real_index_to_unfixed_index[cloth->size() + obj_No];
//		volume = tet_mesh_struct[obj_No]->volume.data();
//		A = tet_mesh_struct[obj_No]->A.data();
//		g_ = b_for_test.data();
//		for (unsigned int i = 0; i < size; ++i) {
//			for (unsigned int j = 0; j < 4; ++j) {
//				is_unfixed[j] = (inv_mass_[indices[i][j]] != 0.0);
//				vertex_position_in_system[j] = 3 * (vertex_index_start + real_index_to_system_index[indices[i][j]]);
//				if (is_unfixed[j]) {
//					g_tet[j] = g_ + vertex_position_in_system[j];
//				}
//				else {
//					g_tet[j] = use_for_temp;
//				}
//			}
//			if (is_unfixed[0] || is_unfixed[1] || is_unfixed[2] || is_unfixed[3]) {
//				setARAPHessianForTest(vertex_pos[indices[i][0]].data(), vertex_pos[indices[i][1]].data(), vertex_pos[indices[i][2]].data(),
//					vertex_pos[indices[i][3]].data(), vertex_position_in_system,
//					&test_nnz, stiffness * volume[i], A[i], is_unfixed,
//					g_tet[0], g_tet[1], g_tet[2], g_tet[3]);
//			}
//		}
//		only_constraint.setFromTriplets(test_nnz.begin(), test_nnz.end());
//		only_constraint *= time_step_square;
//		mass_ = mass[obj_No];
//		unsigned int* unfixed = tetrahedron->data()[obj_No].mesh_struct.unfixed_point_index.data();
//		for (unsigned int i = 0; i < tetrahedron->data()[obj_No].mesh_struct.unfixed_point_index.size(); ++i) {
//			test_nnz.emplace_back(Triplet<double>(3 * (vertex_index_start + i), 3 * (vertex_index_start + i), mass_[unfixed[i]]/time_step_square));
//			test_nnz.emplace_back(Triplet<double>(3 * (vertex_index_start + i)+1, 3 * (vertex_index_start + i)+1, mass_[unfixed[i]] / time_step_square));
//			test_nnz.emplace_back(Triplet<double>(3 * (vertex_index_start + i)+2, 3 * (vertex_index_start + i)+2, mass_[unfixed[i]] / time_step_square));
//		}
//	}
//}






void NewtonMethod::updateARAPHessianFixedStructure()
{
	std::array<double, 3>* vertex_pos;
	double* hessian_coeff_diag = hessian_coeff_diagonal[0].data();
	double stiffness;
	unsigned int size;
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

	double** address = Hessian_coeff_address.data() + tet_hessian_index[0];

	double* g_ = b_thread[0].data();


	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No + cloth->size()];
		stiffness = tetrahedron->data()[obj_No].ARAP_stiffness* time_step_square;
		size = tet_mesh_struct[obj_No]->indices.size();
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
					g_tet[j] = g_ + vertex_position_in_system[j];
				}
				else {
					diagonal_coeff[j] = use_for_temp;
					g_tet[j] = use_for_temp;
				}
			}
			if (is_unfixed[0] || is_unfixed[1] || is_unfixed[2] || is_unfixed[3]) {

				computeARAPHessianFixedStructure(vertex_pos[indices[i][0]].data(), vertex_pos[indices[i][1]].data(), vertex_pos[indices[i][2]].data(), vertex_pos[indices[i][3]].data(),
					diagonal_coeff[0], diagonal_coeff[1], diagonal_coeff[2], diagonal_coeff[3], address, stiffness* volume[i], A[i], is_unfixed,
					g_tet[0], g_tet[1], g_tet[2], g_tet[3]);
			}
		}
	}
}


//void NewtonMethod::setARAPHessianForTest(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3, int* vertex_index,
//	std::vector<Triplet<double>>* hessian_nnz, double stiffness, Matrix<double, 3, 4>& A, bool* is_unfixed, double* g_0, double* g_1, double* g_2, double* g_3)
//{
//
//	Vector3d eigen_value;
//	Vector3d position;
//	Matrix3d deformation_gradient;
//	Matrix3d U, V, rotation;
//	Matrix<double, 12, 12> Hessian;
//	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);
//	FEM::extractRotation(deformation_gradient, eigen_value, U, V, rotation);
//	Matrix<double, 3, 4> grad_C_transpose;
//	grad_C_transpose = (2.0 * stiffness) * (deformation_gradient - rotation) * A;
//
//	Matrix4d result = 2.0 * A.transpose() * A;
//	Hessian.setZero();
//	for (unsigned int i = 0; i < 4; ++i) {
//		for (unsigned int j = 0; j < 4; ++j) {
//			Hessian(3 * i, 3 * j) = result(i, j);
//			Hessian(3 * i + 1, 3 * j + 1) = result(i, j);
//			Hessian(3 * i + 2, 3 * j + 2) = result(i, j);
//		}
//	}
//	Hessian *= stiffness;
//
//	for (unsigned int i = 0; i < 4; ++i) {
//		if (is_unfixed[i]) {
//			for (int j = 0; j < 4; ++j) {
//				if (is_unfixed[j]) {
//					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j], vertex_index[i], Hessian(3 * j, 3 * i)));
//					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 1, vertex_index[i], Hessian(3 * j + 1, 3 * i)));
//					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 2, vertex_index[i], Hessian(3 * j + 2, 3 * i)));
//
//					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j], vertex_index[i] + 1, Hessian(3 * j, 3 * i + 1)));
//					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 1, vertex_index[i] + 1, Hessian(3 * j + 1, 3 * i + 1)));
//					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 2, vertex_index[i] + 1, Hessian(3 * j + 2, 3 * i + 1)));
//
//					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j], vertex_index[i] + 2, Hessian(3 * j, 3 * i + 2)));
//					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 1, vertex_index[i] + 2, Hessian(3 * j + 1, 3 * i + 2)));
//					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 2, vertex_index[i] + 2, Hessian(3 * j + 2, 3 * i + 2)));
//				}
//			}
//		}
//	}
//
//
//	for (unsigned int i = 0; i < 3; ++i) {
//		g_0[i] -= grad_C_transpose.data()[i];
//		g_1[i] -= grad_C_transpose.data()[i + 3];
//		g_2[i] -= grad_C_transpose.data()[i + 6];
//		g_3[i] -= grad_C_transpose.data()[i + 9];
//	}
//
//}


void NewtonMethod::computeARAPHessianFixedStructure(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	double* diagonal_coeff_0, double* diagonal_coeff_1, double* diagonal_coeff_2, double* diagonal_coeff_3, double**& address,
	double stiffness, Matrix<double, 3, 4>& A, bool* is_unfixed, double* g_0, double* g_1, double* g_2, double* g_3)
{
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d S, rotation;
	Matrix<double, 12, 12> Hessian;
	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);
	FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);
	Matrix<double, 3, 4> grad_C_transpose;
	grad_C_transpose =(2.0*stiffness) * (deformation_gradient - rotation) * A;

	Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
	FEM::getHessian(Hessian, S, rotation, Dm, A);
	Hessian *=stiffness;

	//record hessian
	for (unsigned int i = 1; i < 4; ++i) {
		if (is_unfixed[i]) {
			for (int j = 0; j < i; ++j) {
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
		g_0[i] -=  grad_C_transpose.data()[i];
		g_1[i] -=  grad_C_transpose.data()[i+3];
		g_2[i] -=  grad_C_transpose.data()[i+6];
		g_3[i] -=  grad_C_transpose.data()[i+9];
	}

}


void NewtonMethod::saveScene(double* force_direction, int obj_No, bool have_force)
{
	updateTotalVelocity();
	std::string file_name;
	save_scene.save_scene_XPBD(*time_stamp, *time_indicate_for_simu, mesh_struct, 
		&velocity_total, collider_mesh_struct, file_name);
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

void NewtonMethod::readScene(const char* file_name, double* force_direction, int& obj_No)
{
	std::vector<double>force_coe;
	std::vector<int>neighbor_vertex_index;
	save_scene.read_scene_XPBD(file_name, time_stamp, time_indicate_for_simu, mesh_struct, 
		&velocity_total, collider_mesh_struct,*has_force,force_direction, force_coe, neighbor_vertex_index, obj_No);
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
	updateVelocity();
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


//void NewtonMethod::updateTest()
//{
//	b_for_test.setZero();
//	testSetHessian();
//
//	Matrix_test.setFromTriplets(test_nnz.begin(), test_nnz.end());
//	construct_b_Test();
//	Matrix_test *= time_step_square;
//	b_for_test *= time_step_square;
//
//	global_llt.compute(Matrix_test);
//
//	delta_x = global_llt.solve(b_for_test);
//	std::cout << (Matrix_test - Hessian).squaredNorm() << " " << (b_for_test - b).squaredNorm() << std::endl;
//}

void NewtonMethod::solveNewtonMethod_()
{
	storeInitialPosition();
	thread->assignTask(this, SET_S_N);
	
	//if (iteration_number > 1000) {
	//	system("pause");
	//}

	iteration_number = 0;
	std::cout << "===" << std::endl;
	computeEnergy();
	//std::cout << "energy " << total_energy << std::endl;
	previous_energy = total_energy;
	store_residual.clear();
	while (convergenceCondition())
	{
		//updateTest();
		updateHessianFixedStructure();
		////compareTwoMatrix(Hessian, vertex_position, edge_vertices_mass_spring, only_one_vertex_fix_edge_vertices,
		////	edge_length_stiffness[0][0], unfixed_rest_length, fixed_one_vertices_rest_length, unfixed_vertex, mass, time_step);
		////compareTwoForce(Sn, b, vertex_position, edge_vertices_mass_spring, only_one_vertex_fix_edge_vertices,
		////	edge_length_stiffness[0][0], unfixed_rest_length, fixed_one_vertices_rest_length, unfixed_vertex, mass, time_step);
		global_llt.factorize(Hessian);
		delta_x = global_llt.solve(b);



		displacement_coe = 1.0;
		thread->assignTask(this, UPDATE_POSITION_NEWTON);
		//thread->assignTask(this, VELOCITY_NEWTON);
		computeEnergy();
		//std::cout << total_energy << std::endl;

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
	std::cout << total_energy << std::endl;
	if (iteration_number <50) {
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

	for (unsigned int obj_No = 0; obj_No < cloth->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No];
		edge_vertex = edge_vertices_mass_spring[obj_No]->data();
		vertex_begin = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No] << 1;
		vertex_end = edge_index_begin_per_thread_for_mass_spring[obj_No][thread_No + 1] << 1;
		rest_length_ = unfixed_rest_length[obj_No]->data();
		edge_length_stiffness_ = *edge_length_stiffness[obj_No]*time_step_square;
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
	for (unsigned int obj_No = 0; obj_No < cloth->size(); ++obj_No) {
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

	for (unsigned int obj_No = 0; obj_No < cloth->size(); ++obj_No) {
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

	double* hessian_coeff_diag = hessian_coeff_diagonal[thread_No].data();

	Triplet<double>* hessian_nnz_ = hessian_nnz.data() + off_diagonal_hessian_nnz_index_begin_per_thread[thread_No];
	unsigned int edge_no = 0;

	double time_step_square_ = time_step_square;

	double* rest_length_;

	double edge_length_stiffness_;
	unsigned int vertex_index_begin_in_system;

	unsigned int* unfixed_index_to_normal_index;

	for (unsigned int obj_No = 0; obj_No < cloth->size(); ++obj_No) {
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
					vertex_pos[3 * unfixed_index_to_normal_index[l - vertex_start] + k]
						+= time_step_ * velocity.data()[j + k] + time_step_square_ * f_ext.data()[j + k] / mass_[unfixed_index_to_normal_index[l - vertex_start]];
					Sn.data()[j + k] = vertex_pos[3 * unfixed_index_to_normal_index[l - vertex_start] + k];
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
	//computeInertialEnergy();
	computeARAPEnergy();
}

//NEWTON_METHOD_ENERGY
void NewtonMethod::computeEnergy(int thread_No)
{
	energy_per_thread[thread_No] = 0;
	computeMassSpringEnergy(thread_No);
	computeInertial(thread_No);
}



void NewtonMethod::computeARAPEnergy()
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
		vertex_pos = vertex_position[obj_No + cloth->size()];
		tet_end = tet_mesh_struct[obj_No]->indices.size();
		tet_index = tet_indices[obj_No];
		A = tet_mesh_struct[obj_No]->A.data();
		volume = tet_mesh_struct[obj_No]->volume.data();
		stiffness = tetrahedron->data()[obj_No].ARAP_stiffness;
		mass_inv_ = inv_mass[obj_No + cloth->size()];
		for (unsigned int i = 0; i < tet_end; i++) {
			if (mass_inv_[tet_index[i][0]] != 0.0 || mass_inv_[tet_index[i][1]] != 0.0 || mass_inv_[tet_index[i][2]] != 0.0 || mass_inv_[tet_index[i][3]] != 0.0) {
				energy += compute_energy.computeARAPEnergy(vertex_pos[tet_index[i][0]].data(), vertex_pos[tet_index[i][1]].data(), vertex_pos[tet_index[i][2]].data(),
					vertex_pos[tet_index[i][3]].data(), A[i], volume[i], stiffness);
			}

		}
	}
	total_energy += energy;
}



double NewtonMethod::computeInertialEnergy()
{
	double energy = 0.0;
	std::array<double, 3>* vertex_pos;
	unsigned int vertex_index_start, size;
	double* mass;
	unsigned int* unfixed_index_to_normal_index;
	double* sn_single;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		vertex_index_start = vertex_begin_per_obj[i];
		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		size = unfixed_vertex[i]->size();
		mass = mesh_struct[i]->mass.data();
		for (unsigned int k = 0; k < size; ++k) {
			sn_single = Sn.data() + 3 * (vertex_index_start + k);
			energy += mass[unfixed_index_to_normal_index[k]] * (EDGE_LENGTH(vertex_pos[unfixed_index_to_normal_index[k]], sn_single));
		}
	}
	return energy / (2.0 * time_step * time_step);
}
void NewtonMethod::computeInertial(int thread_No)
{

	std::array<double, 3>* vertex_pos;
	double* mass_;

	unsigned int vertex_start;
	unsigned int* unfixed_index_to_normal_index;

	unsigned int j;

	double energy = 0;
	unsigned int end;
	double* sn_single;
	
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		mass_ = mesh_struct[i]->mass.data();
		vertex_start = vertex_begin_per_obj[i];
		vertex_pos = vertex_position[i];

		unfixed_index_to_normal_index = unfixed_vertex[i]->data();
		end = unfixed_vertex_begin_per_thread[i][thread_No + 1];
		for (unsigned int k = unfixed_vertex_begin_per_thread[i][thread_No]; k < end; ++k) {
			sn_single = Sn.data() + 3 * (vertex_start + k);
			energy += mass_[unfixed_index_to_normal_index[k]] *
				(EDGE_LENGTH(vertex_pos[unfixed_index_to_normal_index[k]], sn_single));
		}
	}
	energy /= (2.0 * time_step*time_step);
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
					b.data()[j + m] = b_thread[0][j + m] + f_ext.data()[j + m] * time_step_square_ +
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
					b.data()[j + m] = b_thread[0][j + m] +
						mass_[index] * (Sn.data()[j + m] - vertex_pos[3 * index + m])
						+ coe * mass_[index] * (position_of_beginning[j + m] - vertex_pos[3 * index + m]);
					pos_dis.data()[j + m] = position_of_beginning[j + m] - vertex_pos[3 * index + m];
				}
			}
		}
	}
}

void NewtonMethod::getARAPCoeffAddress()
{
	unsigned int size;
	int vertex_position_in_system[4];
	bool is_unfixed[4];
	std::array<int, 4>* indices;
	double* inv_mass_;
	unsigned int vertex_index_start;
	unsigned int* real_index_to_system_index;


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

				getARAPHessianCoeffAddress(vertex_position_in_system, &Hessian_coeff_address, is_unfixed);

			}
		}
	}
}

void NewtonMethod::getARAPHessianCoeffAddress(int* start_index_in_system, std::vector<double*>* address, bool* is_unfixed)
{
	for (unsigned int i = 1; i < 4; ++i) {
		if (is_unfixed[i]) {
			for (int j = 0; j < i; ++j) {
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
	for (unsigned int obj_No = 0; obj_No < cloth->size(); ++obj_No) {
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


	mesh_struct.resize(total_obj_num);
	collider_mesh_struct.resize(collider->size());


	anchor_vertex.resize(total_obj_num);
	anchor_stiffness.resize(total_obj_num);
	anchor_position.resize(total_obj_num);


	previous_frame_edge_length_stiffness.resize(total_obj_num);

	//vertex_pos_max_step.resize(total_obj_num);
	//for (unsigned int i = 0; i < cloth->size(); ++i) {
	//	vertex_pos_max_step[i]= cloth->data()[i].mesh_struct.vertex_position
	//}

	tet_indices.resize(tetrahedron->size());
	tet_mesh_struct.resize(tetrahedron->size());
	inv_mass.resize(total_obj_num);


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
		anchor_vertex_begin_per_thread[i] = cloth->data()[i].mesh_struct.anchor_index_begin_per_thread.data();
		unfixed_vertex_begin_per_thread[i] = cloth->data()[i].mesh_struct.unfixed_vertex_index_begin_per_thread.data();

		anchor_vertex[i] = &cloth->data()[i].mesh_struct.anchor_vertex;
		anchor_position[i] = &cloth->data()[i].mesh_struct.anchor_position;
		unfixed_vertex[i] = &cloth->data()[i].mesh_struct.unfixed_point_index;
		vertex_begin_per_obj[i + 1] = vertex_begin_per_obj[i] + cloth->data()[i].mesh_struct.unfixed_point_index.size();
		//edge_begin_per_obj[i+1] = edge_begin_per_obj[i]+ (cloth->data()[i].mesh_struct.edge_vertices.size()>>1);
		total_vertex_num[i] = cloth->data()[i].mesh_struct.vertex_position.size();
		real_index_to_unfixed_index[i] = cloth->data()[i].mesh_struct.real_index_to_unfixed_index.data();

		previous_frame_edge_length_stiffness[i] = cloth->data()[i].length_stiffness;

		inv_mass[i] = cloth->data()[i].mesh_struct.mass_inv.data();

		mesh_struct[i] = &cloth->data()[i].mesh_struct;



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


		inv_mass[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.mass_inv.data();
		tet_indices[i] = tetrahedron->data()[i].mesh_struct.indices.data();
		tet_mesh_struct[i] = &tetrahedron->data()[i].mesh_struct;
		mesh_struct[i+cloth->size()] = &tetrahedron->data()[i].mesh_struct;
	}


	for (unsigned int i = 0; i < collider->size(); ++i) {
		collider_mesh_struct[i]= &collider->data()[i].mesh_struct;
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



