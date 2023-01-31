#include"IPC.h"
#include"./basic/write_txt.h"
#include"XPBD/FEM_relate.h"
#include"tet_inversion.h"
#include"NeoHookean.h"

IPC::IPC()
{
	gravity_ = 9.8;
	sub_step_num = 1;
	iteration_number = 300;

	damping_coe = 0.0;

	perform_collision = true;

	velocity_damp = 0.995;
	energy_converge_ratio = 1e-3;

	min_inner_iteration = 4;
	min_outer_iteration = 2;

	max_move_standard = 1e-3;

	max_iteration_number = 100;
	outer_max_iteration_number = 40;
	energy_converge_standard = 1e-10;

	second_order_constraint.solve_exact_ARAP_hessian = false;
	time_step = 0.01;

}

void IPC::setForIPC(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
	Thread* thread, double* tolerance_ratio)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	this->thread = thread;
	total_obj_num = cloth->size() + tetrahedron->size();
	total_thread_num = thread->thread_num;

	has_collider = !collider->empty();

	reorganzieDataOfObjects();
	energy_per_thread.resize(total_thread_num);

	initialVariable();
	setInertial();

	//recordEdgeHessian();
	//computeGravity();
	//setHessian();
	this->floor = floor;

	collision.common_hessian = &common_hessian;
	if (perform_collision) {
		collision.inner_iteration_number = &inner_iteration_number;
		collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio, IPC_, false);
		collision.setCollisionFreeVertex(&record_collision_free_vertex_position_address, &record_vertex_position);
		collision.inner_itr_num_standard = &inner_itr_num_standard;
		collision.min_collision_time = &min_collision_time;
	}
}

bool IPC::convergeCondition(unsigned int iteration_num, bool for_warm_start, double max_move_standard, unsigned int outer_max_iteration_number)
{

	if (iteration_num < min_outer_iteration) {//max_iteration_number
		std::cout << iteration_num << " " << max_displacement << std::endl;
		return false;
	}
	
	if (iteration_num > outer_max_iteration_number - 1) {
		return true;
	}
	if (for_warm_start) {
		std::cout << "warm start " << iteration_num << " " << max_displacement << std::endl;
	}
	else {
		std::cout << iteration_num << " " << max_displacement << std::endl;
	}


	if (max_displacement > max_move_standard) {

		if (iteration_num > 100) {
			if (dis_record.size() > 3) {
				dis_record.clear();
			}
			for (int i = 0; i < dis_record.size(); ++i) {
				if (abs(dis_record[i] - max_displacement) < 1e-8) {
					return true;
				}
			}
			dis_record.emplace_back(max_displacement);
		}
		return false;
	}
	return true;
}



void IPC::initialDHatTolerance(double ave_edge_length)
{
	if (perform_collision) {
		collision.initialDHatTolerance(ave_edge_length);
	}
}

void IPC::addExternalForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int obj_No)
{
	if (!coe.empty()) {
		for (unsigned int i = 0; i < coe.size(); ++i) {
			for (unsigned int j = 0; j < 3; ++j) {
				f_ext[obj_No][neighbor_vertex[i]][j] += coe[i] * neighbor_vertex_force_direction[j];
			}
		}
	}
}


void IPC::updateTetrahedronAnchorVertices()
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

void IPC::IPC_solve()
{
	updateCollisionFreePosition();
	setPosPredict();
	updateSn();
	iteration_number = 0;
	if (perform_collision) {
		collision.initialPairRecordInfo();
	}
	outer_itr_num = 0;
	displacement_satisfied = false;
	inner_itr_num_standard = min_inner_iteration;
	inner_iteration_number = 0;
	if (perform_collision) {
		collision.collisionCulling();
		if (add_pair_by_distance) {
			collision.globalCollisionTime();
		}
		else {
			collision.collisionTimeWithPair(false);
		}
		
		inversionTest();
		computeCollisionFreePosition();
		if (add_pair_by_distance) {
			collision.addPairByDistance();
		}
		updateCollisionFreePosition();
	}
	double ori_energy = computeCurrentEnergy();

	while (!convergeCondition(outer_itr_num, false, max_move_standard, outer_max_iteration_number))
	{
		if (outer_itr_num > 0) {
			updateCollisionFreePosition();
		}
		//solve
		solveSystem();
		if (perform_collision) {
			collision.collisionCulling();
			if (add_pair_by_distance) {
				collision.globalCollisionTime();
			}
			else {
				collision.collisionTimeWithPair(false);
			}
			inversionTest();
			computeCollisionFreePosition();
			if (add_pair_by_distance) {
				collision.addPairByDistance();
			}
			//line search
			lineSearch(ori_energy);
		}
		outer_itr_num++;
	}
	iteration_number += inner_iteration_number;
	computeVelocity();
	updatePosition();
	updateRenderNormal();

	updateRenderVertexNormal();
}



double IPC::computeInertialEnergy()
{
	double energy = 0.0;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* sn_;
	unsigned int vertex_end;
	double* mass;
	unsigned int start;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		sn_ = sn[i].data();
		vertex_end = vertex_index_begin_per_thread[i][total_thread_num];
		mass = mesh_struct[i]->mass.data();
		for (unsigned int j = 0; j < vertex_end; ++j) {
			energy += mass[j] * (EDGE_LENGTH(vertex_pos[j], sn_[j]));
		}
	}
	return  0.5 * energy / (time_step * time_step);
}


//ARAP_ENERGY
double IPC::computeElasticEnergy()
{
	double energy = 0.0;
	unsigned int end, start;
	std::array<int, 4>* indices;
	std::array<double, 3>* vertex_pos;
	Matrix<double, 3, 4>* A;
	double* volume;
	//double stiffness;
	double* mass_inv_;
	double tet_energy;
	double mu; double lambda;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		end = tet_index_begin_per_thread[i + cloth->size()][total_thread_num];
		start = 0;
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		vertex_pos = vertex_position[i + cloth->size()];
		A = tetrahedron->data()[i].mesh_struct.A.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		//stiffness = tetrahedron->data()[i].ARAP_stiffness;
		mass_inv_ = tetrahedron->data()[i].mesh_struct.mass_inv.data();
		LameCoeff(tetrahedron->data()[i].youngs_modulus, tetrahedron->data()[i].poisson_ratio,
			mu, lambda);
		for (unsigned int j = start; j < end; ++j) {
			if (mass_inv_[indices[j][0]] != 0.0 || mass_inv_[indices[j][1]] != 0.0 || mass_inv_[indices[j][2]] != 0.0 || mass_inv_[indices[j][3]] != 0.0) {
				tet_energy = NeoHookean::energy(vertex_pos[indices[j][0]].data(), vertex_pos[indices[j][1]].data(),
					vertex_pos[indices[j][2]].data(), vertex_pos[indices[j][3]].data(), A[j], volume[j], mu, lambda);
				energy += tet_energy;
				//if (!NeoHookean::testDeterminent(vertex_pos[indices[j][0]].data(), vertex_pos[indices[j][1]].data(),
				//	vertex_pos[indices[j][2]].data(), vertex_pos[indices[j][3]].data(), A[j])) {
				//	std::cout << "error tet index " << j << std::endl;
				//	system("pause");
				//}
			}
		}
	}
	return energy;
}


double IPC::computeCurrentEnergy()
{
	double energy = 0;
	double inertial_energy = computeInertialEnergy();
	energy += inertial_energy;
	double elastic_energy = computeElasticEnergy();
	energy += elastic_energy;
	double barrier_energy = 0.0;
	if (perform_collision) {
		barrier_energy = computeBarrierEnergy();
		energy += barrier_energy;
	}

	std::cout <<"energy " << energy << " " << inertial_energy << " " << elastic_energy << " " << barrier_energy << std::endl;

	return energy;
}

double IPC::computeBarrierEnergy()
{
	double energy = 0.0;
	double stiffness = 0.0;
	if (!tetrahedron->empty()) {
		stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		stiffness = cloth->data()[0].collision_stiffness[0];
	}
	//vt
	energy += computeVTEnergy(&collision.record_vt_pair_sum_all_thread, triangle_indices.data(), vertex_position.data(), vertex_position.data(),
		stiffness, collision.record_vt_pair_d_hat.data(),
		0, collision.record_vt_pair_sum_all_thread.size());
	//ee
	energy += computeEEEnergy(&collision.record_ee_pair_sum_all_thread, edge_vertices.data(), edge_vertices.data(),
		vertex_position.data(), vertex_position.data(), stiffness, collision.record_ee_pair_d_hat.data(),
		0,
		collision.record_ee_pair_sum_all_thread.size());
	if (has_collider) {
		//vt_c
		energy += computeVTEnergy(&collision.record_vt_collider_pair_sum_all_thread, triangle_indices_collider.data(), vertex_position.data(), vertex_position_collider.data(),
			stiffness, collision.record_vt_collider_pair_d_hat.data(),
			0,
			collision.record_vt_collider_pair_sum_all_thread.size());
		//tv_c
		energy += computeVTEnergy(&collision.record_tv_collider_pair_sum_all_thread, triangle_indices.data(), vertex_position_collider.data(), vertex_position.data(),
			stiffness, collision.record_tv_collider_pair_d_hat.data(),
			0,
			collision.record_tv_collider_pair_sum_all_thread.size());
		//ee c
		energy += computeEEEnergy(&collision.record_ee_collider_pair_sum_all_thread, edge_vertices.data(), collider_edge_vertices.data(),
			vertex_position.data(), vertex_position_collider.data(), stiffness, collision.record_ee_collider_pair_d_hat.data(),
			0,
			collision.record_ee_collider_pair_sum_all_thread.size());
	}
	if (floor->exist) {
		energy += computeFloorEnergy(0, stiffness, &collision.record_vertex_collide_with_floor_sum_all_thread, 0,
			collision.record_vertex_collide_with_floor_sum_all_thread.size());
	}
	return energy;
}
double IPC::computeFloorEnergy(int type, double collision_stiffness, std::vector<unsigned int>* record_vertex_collide_with_floor, int start, int end)
{
	auto start_ = record_vertex_collide_with_floor->begin() + start;
	auto end_ = record_vertex_collide_with_floor->begin() + end;

	double energy = 0.0;
	if (type == 0) {
		for (auto i = start_; i < end_; i += 2) {
			energy += compute_energy.computeFloorBarrierEnergy(vertex_position[*i][*(i + 1)][floor->dimension], collision.record_vertex_collide_with_floor_d_hat[*i][*(i + 1)],
				collision_stiffness, floor->value);
		}
	}
	else {
		for (auto i = start_; i < end_; i += 2) {
			if (collision.vertex_belong_to_color_group[*i][*(i + 1)]) {
				energy += compute_energy.computeFloorBarrierEnergy(vertex_position[*i][*(i + 1)][floor->dimension], collision.record_vertex_collide_with_floor_d_hat[*i][*(i + 1)],
					collision_stiffness, floor->value);
			}
		}
	}

	return energy;

}


double IPC::computeVTEnergy(std::vector<unsigned int>* record_vt_pair, std::array<int, 3>** triangle_indices,
	std::array<double, 3>** v_current_pos, std::array<double, 3>** t_current_pos, double collision_stiffness, double* d_hat,
	unsigned int start, unsigned int end)
{
	double energy = 0.0;
	int* indices;
	auto pair_end = record_vt_pair->begin() + end;
	auto pair_start = record_vt_pair->begin() + start;
	double* d_hat_ = d_hat + (start >> 2);
	for (auto i = pair_start; i < pair_end; i += 4) {
		indices = triangle_indices[*(i + 2)][*(i + 3)].data();
		energy += compute_energy.computeBarrierEnergy(v_current_pos[*i][*(i + 1)].data(),
			t_current_pos[*(i + 2)][indices[0]].data(),
			t_current_pos[*(i + 2)][indices[1]].data(),
			t_current_pos[*(i + 2)][indices[2]].data(), collision_stiffness, *d_hat_, true);
		d_hat_++;
	}
	return energy;
}

double IPC::computeEEEnergy(std::vector<unsigned int>* record_pair, unsigned int** edge_v_0, unsigned int** edge_v_1,
	std::array<double, 3>** e0_current_pos, std::array<double, 3>** e1_current_pos, double collision_stiffness, double* d_hat,
	unsigned int start, unsigned int end)
{
	double energy = 0.0;
	unsigned int* edge_0_vertex;
	unsigned int* edge_1_vertex;

	auto pair_end = record_pair->begin() + end;
	auto pair_start = record_pair->begin() + start;
	double* d_hat_ = d_hat + (start >> 2);

	for (auto i = pair_start; i < pair_end; i += 4) {
		edge_0_vertex = edge_v_0[*i] + ((*(i + 1)) << 1);
		edge_1_vertex = edge_v_1[*(i + 2)] + ((*(i + 3)) << 1);
		energy += compute_energy.computeBarrierEnergy(e0_current_pos[*i][*edge_0_vertex].data(),
			e0_current_pos[*i][*(edge_0_vertex + 1)].data(), e1_current_pos[*(i + 2)][*edge_1_vertex].data(),
			e1_current_pos[*(i + 2)][*(edge_1_vertex + 1)].data(), collision_stiffness, *d_hat_, false);
		d_hat_++;
	}
	return energy;
}


void IPC::lineSearch(double& ori_energy)
{
	double current_energy = 0.0;
	if (perform_collision) {
		current_energy = computeCurrentEnergy();
		if (current_energy > energy_converge_standard) {
			double record_collision_time = collision.collision_time;
			if (current_energy > ori_energy) {
				if (record_collision_time < 1e-6) {
					std::cout << "global ccd original record_collision_time is too small " << ori_energy << " " << current_energy << std::endl;
				}
			}
			collision.collision_time = 0.5;
			int itr_num = 0;
			while (current_energy > ori_energy)
			{
				record_collision_time *= 0.5;
				computeCollisionFreePosition();
				if (add_pair_by_distance) {
					collision.addPairByDistance();
				}
				itr_num++;
				if (itr_num > 6) {
					//std::cout << "global record_collision_time is too small " << ori_energy << " " << current_energy << std::endl;
					break;
				}
				current_energy = computeCurrentEnergy();
			}
		}
		ori_energy = current_energy;
	}
}

void IPC::solveSystem()
{
	b.setZero();
	hessian_nnz.resize(record_size_of_inertial_triplet);

	setNeoHookean();
	if (perform_collision) {
		setCollisionHessian();
	}
	setInertialGrad(b.data());
	Hessian.setFromTriplets(hessian_nnz.begin(), hessian_nnz.end());
	global_llt.compute(Hessian);
	//std::cout << "result " << std::endl;
	//std::cout << delta_x.transpose() << std::endl;
	//std::cout << Hessian << std::endl;

	delta_x = global_llt.solve(b);
	max_displacement = delta_x.lpNorm<Eigen::Infinity>();
	updatePos();
}

void IPC::setInertialGrad(double* grad)
{
	int num;
	int prefix;
	double* mass_;
	double m_t;

	double* mass_inv;
	for (int i = 0; i < total_obj_num; ++i) {
		prefix = vertex_index_prefix_sum_obj[i];
		num = mesh_struct[i]->vertex_position.size();
		mass_ = mass[i];
		mass_inv = mesh_struct[i]->mass_inv.data();
		for (int j = 0; j < num; ++j) {
			if (mass_inv[j]!=0.0) {
				m_t = mass_[j] / (time_step * time_step);
				grad[3 * (prefix + j)] += m_t * (vertex_position[i][j][0] - sn[i][j][0]);
				grad[3 * (prefix + j) + 1] += m_t * (vertex_position[i][j][1] - sn[i][j][1]);
				grad[3 * (prefix + j) + 2] += m_t * (vertex_position[i][j][2] - sn[i][j][2]);
			}			
		}
	}
}

void IPC::updateRenderVertexNormal()
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



void IPC::updatePosition()
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


void IPC::updateRenderNormal()
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


void IPC::updatePos()
{
	std::array<double, 3>* vertex_pos;
	unsigned int vertex_end = 0;
	unsigned int prefix_;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		prefix_ = vertex_index_prefix_sum_obj[i];
		vertex_end = vertex_index_begin_per_thread[i][total_thread_num];
		for (unsigned int j = 0; j < vertex_end; ++j) {
			for (unsigned int k = 0; k < 3; ++k) {
				vertex_pos[j][k] -= delta_x[3 * (prefix_ + j) + k];
			}
		}
	}
}

void IPC::computeVelocity()
{
	double damp_coe = velocity_damp;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	unsigned int vertex_end = 0;
	std::array<double, 3>* velocity_;
	double delta_t = time_step / damp_coe;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		initial_vertex_pos = vertex_position_render[i];
		vertex_end = vertex_index_begin_per_thread[i][total_thread_num];
		velocity_ = velocity[i].data();
		for (unsigned int j = 0; j < vertex_end; ++j) {
			for (unsigned int k = 0; k < 3; ++k) {
				velocity_[j][k] = (vertex_pos[j][k] - initial_vertex_pos[j][k]) / delta_t;
			}
		}
	}
}


void IPC::computeCollisionFreePosition()
{
	double collision_time = collision.collision_time;

	if (collision_time == 1.0) {
		return;
	}
	unsigned int index_end;
	std::array<double, 3>* q_pre;
	std::array<double, 3>* q_end;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		index_end = vertex_index_begin_per_thread[i][total_thread_num];
		q_end = vertex_position[i];
		q_pre = record_vertex_position[i].data();

		for (unsigned int j =0; j < index_end; ++j) {
			q_end[j][0] = q_pre[j][0] + collision_time * (q_end[j][0] - q_pre[j][0]);
			q_end[j][1] = q_pre[j][1] + collision_time * (q_end[j][1] - q_pre[j][1]);
			q_end[j][2] = q_pre[j][2] + collision_time * (q_end[j][2] - q_pre[j][2]);
		}
	}
}

void IPC::updateSn()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(sn[i][0].data(), vertex_position[i][0].data(), mesh_struct[i]->vertex_position.size() * 24);
	}
}

void IPC::inversionTest()
{
	unsigned int start, end;
	unsigned int* tet_around_a_group;

	MatrixXd grad;
	grad.resize(3, 4);
	std::array<int, 4>* tet_indices_;

	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* ori_vertex_pos;
	double stiffness = 0.0;
	Matrix<double, 3, 4>* A;
	double* volume;

	int prefix_sum_start;

	int k = 0;

	//int color_group_index;
	unsigned int* tet_in_a_group;

	//solve all colors except the last color
	double* grad_address;

	double collision_time = 1.0;
	double inversion_time;
	int* tet_vertex;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		vertex_pos = vertex_position[i + cloth->size()];
		ori_vertex_pos = record_vertex_position[i + cloth->size()].data();

		tet_indices_ = tet_indices[i + cloth->size()];
		start = 0;
		end = tet_index_begin_per_thread[i + cloth->size()][total_thread_num];
		for (auto j = start; j < end; ++j) {
			tet_vertex = tet_indices_[j].data();
			if (inversionTest::TetInversionTest(ori_vertex_pos[tet_vertex[0]].data(), ori_vertex_pos[tet_vertex[1]].data(),
				ori_vertex_pos[tet_vertex[2]].data(), ori_vertex_pos[tet_vertex[3]].data(), vertex_pos[tet_vertex[0]].data(),
				vertex_pos[tet_vertex[1]].data(), vertex_pos[tet_vertex[2]].data(), vertex_pos[tet_vertex[3]].data(), &inversion_time)) {
				if (inversion_time < collision_time) {
					collision_time = inversion_time;
				}
			}
		}
	}

	if (collision_time < 1.0) {
		collision_time *= 0.9;
	}
	if (collision_time < 0.0) {
		std::cout << "error for inversion test time inversionTest()" << std::endl;
	}
	if (collision_time == 0.0) {
		std::cout << "inversion test time inversionTest() equals 0" << std::endl;
	}

	if (collision.collision_time > collision_time) {
		collision.collision_time = collision_time;
	}
}


void IPC::setPosPredict()
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
		vertex_end = mesh_struct[i]->vertex_position.size();
		mass_inv = mesh_struct[i]->mass_inv.data();
		f_ext_ = f_ext[i].data();
		velocity_ = velocity[i].data();
		for (unsigned int j = 0; j < vertex_end; ++j) {
			if (mass_inv[j] != 0) {
				for (unsigned int k = 0; k < 3; ++k) {
					vertex_pos[j][k] += delta_t * velocity_[j][k] + delta_t_2 * (mass_inv[j] * f_ext_[j][k] + gravity__[k]);
				}
			}
		}
	}
}



void IPC::setCollisionHessian()
{
	//vt
	computeVTHessian(0, collision.record_vt_pair_sum_all_thread.size(), collision.record_vt_pair_sum_all_thread.data(),
		collision.record_vt_pair_d_hat.data(), 0, b.data(), hessian_nnz);
	//ee
	computeEEHessian(0, collision.record_ee_pair_sum_all_thread.size(), collision.record_ee_pair_sum_all_thread.data(),
		collision.record_ee_pair_d_hat.data(), b.data(), hessian_nnz,false);

	if (has_collider) {
		//tv_c
		computeVTHessian(0, collision.record_tv_collider_pair_sum_all_thread.size(), collision.record_tv_collider_pair_sum_all_thread.data(),
			collision.record_tv_collider_pair_d_hat.data(), 2, b.data(), hessian_nnz);
		//ee_c
		computeEEHessian(0, collision.record_ee_collider_pair_sum_all_thread.size(), collision.record_ee_collider_pair_sum_all_thread.data(),
			collision.record_ee_collider_pair_d_hat.data(), b.data(), hessian_nnz, true);
		//vt_c
		computeVTHessian(0, collision.record_vt_collider_pair_sum_all_thread.size(), collision.record_vt_collider_pair_sum_all_thread.data(),
			collision.record_vt_collider_pair_d_hat.data(), 1, b.data(), hessian_nnz);
	}
	if (floor->exist) {
		computeFloorHessian();
	}
}

void IPC::computeFloorHessian()
{
	double stiffness;
	double floor_value = floor->value;
	int floor_dimension = floor->dimension;
	unsigned int global_index;
	double hessian, grad;
	if (!tetrahedron->empty()) {
		stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		stiffness = cloth->data()[0].collision_stiffness[0];
	}
	for (auto i = collision.record_vertex_collide_with_floor_sum_all_thread.begin(); i < collision.record_vertex_collide_with_floor_sum_all_thread.end(); i+=2) {
		if (collision.computeFloorHessian(collision.d_hat_2, stiffness, floor_value, &hessian, &grad,
			vertex_position[*i][*(i + 1)][floor_dimension], true)) {
			global_index = vertex_index_prefix_sum_obj[*i] + *(i + 1);
			hessian_nnz.emplace_back(Triplet<double>(3 * global_index + floor->dimension,
				3 * global_index + floor->dimension, hessian));
			std::cout << b[3 * global_index + floor_dimension] << std::endl;
			b[3 * global_index + floor_dimension] += grad;
			std::cout << b[3 * global_index + floor_dimension] << std::endl;
		}
	}
}

//0 ee, 1 ee collider
void IPC::computeEEHessian(int start, int end, unsigned int* pair, double* d_hat, double* grad_,
	std::vector<Triplet<double>>& hessian, bool is_collider)
{
	double stiffness = 0.0;

	unsigned int vertex_index_in_sum[4];
	if (!tetrahedron->empty()) {
		stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		stiffness = cloth->data()[0].collision_stiffness[0];
	}

	int obj_1, obj_2;
	unsigned int* edge_0_vertex;
	unsigned int* edge_1_vertex;

	MatrixXd Hessian; VectorXd grad;
	int hessian_record_index[5];
	memset(hessian_record_index, 0, 20);
	char* record_exist;
	if (is_collider) {
		bool not_collider[4] = { true,true,false, false };
		for (int i = start; i < end; i += 4) {
			edge_0_vertex = edge_vertices[pair[i]] + (pair[i + 1] << 1);
			obj_2 = pair[i + 2];
			edge_1_vertex = collider_edge_vertices[obj_2] + (pair[i + 3] << 1);
			vertex_index_in_sum[0] = vertex_index_prefix_sum_obj[pair[i]] + edge_0_vertex[0];
			vertex_index_in_sum[1] = vertex_index_prefix_sum_obj[pair[i]] + edge_0_vertex[1];
			if (second_order_constraint.computeBarrierEEGradientHessian(vertex_position[pair[i]][edge_0_vertex[0]].data(),
				vertex_position[pair[i]][edge_0_vertex[1]].data(),
				vertex_position_collider[obj_2][edge_1_vertex[0]].data(), vertex_position_collider[obj_2][edge_1_vertex[1]].data(), Hessian, grad,
				hessian_record_index, stiffness, d_hat[i >> 2], rest_edge_length[pair[i]][pair[i + 1]],
				rest_edge_length_collider[obj_2][pair[i + 3]], true, is_collider)) {
				setHessian(hessian_record_index, vertex_index_in_sum, Hessian, grad.data(), hessian, grad_, not_collider);
			}
		}
	}
	else {
		bool not_collider[4] = { true,true,true, true };
		for (int i = start; i < end; i += 4) {
			edge_0_vertex = edge_vertices[pair[i]] + (pair[i + 1] << 1);
			obj_2 = pair[i + 2];
			edge_1_vertex = edge_vertices[obj_2] + (pair[i + 3] << 1);
			vertex_index_in_sum[0] = vertex_index_prefix_sum_obj[pair[i]] + edge_0_vertex[0];
			vertex_index_in_sum[1] = vertex_index_prefix_sum_obj[pair[i]] + edge_0_vertex[1];
			vertex_index_in_sum[2] = vertex_index_prefix_sum_obj[obj_2] + edge_1_vertex[0];
			vertex_index_in_sum[3] = vertex_index_prefix_sum_obj[obj_2] + edge_1_vertex[1];
			if (second_order_constraint.computeBarrierEEGradientHessian(vertex_position[pair[i]][edge_0_vertex[0]].data(),
				vertex_position[pair[i]][edge_0_vertex[1]].data(),
				vertex_position[obj_2][edge_1_vertex[0]].data(), vertex_position[obj_2][edge_1_vertex[1]].data(), Hessian, grad,
				hessian_record_index, stiffness, d_hat[i >> 2], rest_edge_length[pair[i]][pair[i + 1]],
				rest_edge_length[obj_2][pair[i + 3]], true, is_collider)) {				
				setHessian(hessian_record_index, vertex_index_in_sum, Hessian, grad.data(), hessian, grad_, not_collider);
			}
		}
	}
}

void IPC::computeVTHessian(int start, int end, unsigned int* pair, double* d_hat, int type, double* grad_, 
	std::vector<Triplet<double>>& hessian)
{
	double stiffness = 0.0;
	unsigned int vertex_index_in_sum[4];
	if (!tetrahedron->empty()) {
		stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		stiffness = cloth->data()[0].collision_stiffness[0];
	}
	int obj_1, obj_2, vertex_index;
	int* triangle_vertex;
	MatrixXd Hessian; VectorXd grad;
	int hessian_record_index[5] = { 0,-1,-1,-1,-1 };
	switch (type)
	{
	case 0: {
		bool not_collider[4] = { true,true, true,true };
		for (int i = start; i < end; i += 4) {
			vertex_index = pair[i + 1];
			obj_2 = pair[i + 2];
			triangle_vertex = triangle_indices[obj_2][pair[i + 3]].data();
			vertex_index_in_sum[0] = vertex_index_prefix_sum_obj[pair[i]] + vertex_index;
			vertex_index_in_sum[1] = vertex_index_prefix_sum_obj[obj_2] + triangle_vertex[0];
			vertex_index_in_sum[2] = vertex_index_prefix_sum_obj[obj_2] + triangle_vertex[1];
			vertex_index_in_sum[3] = vertex_index_prefix_sum_obj[obj_2] + triangle_vertex[2];
			if (second_order_constraint.computeBarrierVTGradientHessian(Hessian, grad, vertex_position[pair[i]][vertex_index].data(), vertex_position[obj_2][triangle_vertex[0]].data(),
				vertex_position[obj_2][triangle_vertex[1]].data(), vertex_position[obj_2][triangle_vertex[2]].data(), d_hat[i >> 2],
				hessian_record_index, stiffness, true, type)) {
				setHessian(hessian_record_index, vertex_index_in_sum, Hessian, grad.data(), hessian, grad_, not_collider);
			}
		}
	}
		break;

	case 1:
	{
		bool not_collider[4] = { true,false,false, false };
		for (int i = start; i < end; i += 4) {
			vertex_index = pair[i + 1];
			obj_2 = pair[i + 2];
			triangle_vertex = triangle_indices_collider[obj_2][pair[i + 3]].data();
			vertex_index_in_sum[0] = vertex_index_prefix_sum_obj[pair[i]] + vertex_index;
			if (second_order_constraint.computeBarrierVTGradientHessian(Hessian, grad, vertex_position[pair[i]][vertex_index].data(), vertex_position_collider[obj_2][triangle_vertex[0]].data(),
				vertex_position_collider[obj_2][triangle_vertex[1]].data(), vertex_position_collider[obj_2][triangle_vertex[2]].data(), d_hat[i >> 2],
				hessian_record_index, stiffness, true, type)) {
				//std::cout << "===" << std::endl;
				//std::cout << grad.transpose() << std::endl;
				//std::cout << Hessian << std::endl;
				setHessian(hessian_record_index, vertex_index_in_sum, Hessian, grad.data(), hessian, grad_, not_collider);

			}
		}
	}
		break;

	case 2:
	{
		bool not_collider[4] = { false,true, true,true };
		for (int i = start; i < end; i += 4) {
			vertex_index = pair[i + 1];
			obj_2 = pair[i + 2];
			triangle_vertex = triangle_indices[obj_2][pair[i + 3]].data();
			vertex_index_in_sum[1] = vertex_index_prefix_sum_obj[obj_2] + triangle_vertex[0];
			vertex_index_in_sum[2] = vertex_index_prefix_sum_obj[obj_2] + triangle_vertex[1];
			vertex_index_in_sum[3] = vertex_index_prefix_sum_obj[obj_2] + triangle_vertex[2];
			if (second_order_constraint.computeBarrierVTGradientHessian(Hessian, grad, vertex_position_collider[pair[i]][vertex_index].data(), vertex_position[obj_2][triangle_vertex[0]].data(),
				vertex_position[obj_2][triangle_vertex[1]].data(), vertex_position[obj_2][triangle_vertex[2]].data(), d_hat[i >> 2],
				hessian_record_index, stiffness, true, type)) {
				setHessian(hessian_record_index, vertex_index_in_sum, Hessian, grad.data(), hessian, grad_, not_collider);

			}
		}
	}
		break;
	}
}

void IPC::updateItrInfo(int* iteration_num)
{
	iteration_num[LOCAL_GLOBAL] = iteration_number;
	iteration_num[OUTER] = outer_itr_num;
	//outer_iteration_number = iteration_num[OUTER];
	//sub_step_num = iteration_num[OUTER];
}


void IPC::setHessian(int* record_index, unsigned int* vertex_index_total, MatrixXd& Hessian, double* grad,
	std::vector<Triplet<double>>& hessian, double* common_grad, bool* not_collider)
{
	int size = record_index[0];
	record_index += 1;
	double* result;

	double* hessian_locate;
	double* grad_locate;

	double* grad_in_global;

	for (int i = 0; i < size; ++i) {
		if (not_collider[record_index[i]]) {
			for (int j = 0; j < size; ++j) {
				if (not_collider[record_index[j]]) {
					for (int m = 0; m < 3; ++m) {
						for (int n = 0; n < 3; ++n) {
							hessian.emplace_back(Triplet<double>(3 * vertex_index_total[record_index[i]] + m,
								3 * vertex_index_total[record_index[j]] + n, Hessian(3 * i + m, 3 * j + n)));
						}
					}
				}
			}			
			grad_locate = grad + 3 * i;
			grad_in_global = common_grad + 3 * vertex_index_total[record_index[i]];
			grad_in_global[0] += *grad_locate;
			grad_in_global[1] += *(grad_locate + 1);
			grad_in_global[2] += *(grad_locate + 2);
		}		
	}
}


void IPC::setNeoHookean()
{
	std::array<double, 3>* vertex_pos;
	unsigned int size;
	double* inv_mass_;

	bool is_unfixed[4];
	std::array<int, 4>* indices;
	unsigned int vertex_index_start;

	int vertex_position_in_system[4];


	double use_for_temp[9];

	double* volume;
	Matrix<double, 3, 4>* A;
	//std::vector<bool>* this_is_vertex_fixed;



	double mu; double lambda;

	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No + cloth->size()];
		size = tetrahedron->data()[obj_No].mesh_struct.indices.size();
		indices = tet_indices[obj_No];
		vertex_index_start = vertex_index_prefix_sum_obj[cloth->size() + obj_No];
		volume = tet_volume[cloth->size() + obj_No];
		A = tet_A[cloth->size() + obj_No];
		inv_mass_ = mesh_struct[obj_No + cloth->size()]->mass_inv.data();
		//this_is_vertex_fixed = is_vertex_fixed[cloth->size() + obj_No];
		LameCoeff(tetrahedron->data()[obj_No].youngs_modulus, tetrahedron->data()[obj_No].poisson_ratio,
			mu, lambda);
		for (unsigned int i = 0; i < size; ++i) {
			for (unsigned int j = 0; j < 4; ++j) {

				is_unfixed[j] = (inv_mass_[indices[i][j]]!=0);
				vertex_position_in_system[j] = 3 * (vertex_index_start + indices[i][j]);
			}
			if (is_unfixed[0] || is_unfixed[1] || is_unfixed[2] || is_unfixed[3]) {
				computeARAPHessian(vertex_pos[indices[i][0]].data(), vertex_pos[indices[i][1]].data(), vertex_pos[indices[i][2]].data(), vertex_pos[indices[i][3]].data(),
					&hessian_nnz, vertex_position_in_system,
					mu, lambda, A[i], is_unfixed, volume[i]);
			}
		}
	}
}

void IPC::setInertial()
{
	std::array<double, 3>* vertex_pos;

	std::vector<Triplet<double>>* hessian_nnz_ = &hessian_nnz;
	double stiffness;
	unsigned int size;
	double* inv_mass_;

	bool is_unfixed[4];
	std::array<int, 4>* indices;
	unsigned int vertex_index_start;

	int vertex_position_in_system[4];


	double use_for_temp[9];

	double* volume;
	Matrix<double, 3, 4>* A;
	std::vector<bool>* this_is_vertex_fixed;

	hessian_nnz.reserve(18*vertex_index_prefix_sum_obj[total_obj_num]);

	double* mass_;
	double mass_dt_2;
	int index;
	for (int i = 0; i < total_obj_num; ++i) {
		size = mesh_struct[i]->vertex_position.size();
		mass_ = mass[i];
		vertex_index_start = vertex_index_prefix_sum_obj[i];
		for (int j = 0; j < size; ++j) {			
			mass_dt_2 = mass_[j] / (time_step * time_step);
			index = 3 * (vertex_index_start + j);
			hessian_nnz.emplace_back(Triplet<double>(index, index, mass_dt_2));
			hessian_nnz.emplace_back(Triplet<double>(index + 1, index + 1, mass_dt_2));
			hessian_nnz.emplace_back(Triplet<double>(index + 2, index + 2, mass_dt_2));
		}
	}
	record_size_of_inertial_triplet = hessian_nnz.size();
}


void IPC::computeARAPHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	std::vector<Triplet<double>>* hessian_nnz,
	int* vertex_index, double mu, double lambda, Matrix<double, 3, 4>& A, bool* is_unfixed, double volume)
{
	MatrixXd Hessian; VectorXd grad(12);
	Hessian.resize(12, 12);

	NeoHookean::gradientHessianMPM(vertex_position_0, vertex_position_1,
		vertex_position_2, vertex_position_3,
		A, volume, mu, lambda, grad, Hessian);

	for (int i = 0; i < 4; ++i) {
		if (is_unfixed[i]) {
			for (int j = 0; j < 4; ++j) {
				if (is_unfixed[j]) {
					for (int m = 0; m < 3; ++m) {
						for (int n = 0; n < 3; ++n) {
							hessian_nnz->emplace_back(Triplet<double>(vertex_index[j]+n, vertex_index[i]+m, Hessian(j*3+n, i*3+m)));
						}
					}				
				}
			}
			b[vertex_index[i]] += grad[3 * i];
			b[vertex_index[i] + 1] += grad[3 * i + 1];
			b[vertex_index[i]+ 2 ] += grad[3 * i + 2];
		}
	}
}



void IPC::updateCollisionFreePosition()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(record_vertex_position[i][0].data(), vertex_position[i][0].data(), 24 * record_vertex_position[i].size());
		memcpy(record_collision_free_vertex_position[i][0].data(), vertex_position[i][0].data(), 24 * record_vertex_position[i].size());
	}
}

void IPC::initialVariable()
{
	gravity[0] = 0;
	gravity[1] = -gravity_;
	gravity[2] = 0;

	Hessian.resize(3 * vertex_index_prefix_sum_obj[total_obj_num], 3 * vertex_index_prefix_sum_obj[total_obj_num]);
	b.resize(3 * vertex_index_prefix_sum_obj[total_obj_num]);

	f_ext.resize(total_obj_num);
	velocity.resize(total_obj_num);
	sn.resize(total_obj_num);
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		f_ext[i].resize(mesh_struct[i]->vertex_position.size(), {0.0,0.0,0.0});
		velocity[i].resize(mesh_struct[i]->vertex_position.size(), { 0.0,0.0,0.0 });
		sn[i].resize(mesh_struct[i]->vertex_position.size(), { 0.0,0.0,0.0 });
	}
	acceleration.resize(3 * vertex_index_prefix_sum_obj[total_obj_num]);
	acceleration.setZero();

	grad_max_store.resize(total_thread_num);
	max_dis_record.resize(total_thread_num);
}

void IPC::reorganzieDataOfObjects()
{
	energy_per_thread.resize(total_thread_num);
	vertex_position.resize(total_obj_num);
	vertex_position_render.resize(total_obj_num);
	mesh_struct.resize(total_obj_num);
	vertex_index_begin_per_thread.resize(total_obj_num);
	tet_index_begin_per_thread.resize(total_obj_num);
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
	vertex_index_of_a_tet_color_group.resize(total_obj_num);

	vertex_index_of_a_tet_color_per_thread_start_group.resize(total_obj_num);

	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position.data();
		vertex_position_render[i] = cloth->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i] = &cloth->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		record_vertex_position_address[i] = record_vertex_position[i].data();
		record_collision_free_vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position;
		record_collision_free_vertex_position_address[i] = record_collision_free_vertex_position[i].data();
		triangle_indices[i] = cloth->data()[i].mesh_struct.triangle_indices.data();
		//record_outer_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		//unfixed_vertex[i] = &cloth->data()[i].mesh_struct.unfixed_point_index;
		edge_vertices[i] = cloth->data()[i].mesh_struct.edge_vertices.data();

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
	//tet_color_groups_label.resize(tetrahedron->size());

	tet_around_tet.resize(tetrahedron->size());

	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		vertex_position_render[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		record_vertex_position_address[i + cloth->size()] = record_vertex_position[i + cloth->size()].data();
		record_collision_free_vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position;
		record_collision_free_vertex_position_address[i + cloth->size()] = record_collision_free_vertex_position[i + cloth->size()].data();
		triangle_indices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.triangle_indices.data();
		tet_indices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.indices.data();
		tet_A[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.A.data();
		tet_volume[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.volume.data();
		//record_outer_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		//unfixed_vertex[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.unfixed_point_index;
		edge_vertices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_vertices.data();

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

		tet_color_groups[i] = &tetrahedron->data()[i].mesh_struct.tet_color_group;
		//tet_color_groups_label[i] =tetrahedron->data()[i].mesh_struct.tet_in_collision.data();

		tet_around_tet[i] = tetrahedron->data()[i].mesh_struct.tet_tet_index.data();
		vertex_index_of_a_tet_color_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_of_a_tet_color_group.data();

		tet_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tetrahedron_index_begin_per_thread.data();

		vertex_index_of_a_tet_color_per_thread_start_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_of_a_tet_color_per_thread_start_group.data();
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
	vertex_num_on_surface_prefix_sum.resize(total_obj_num + 1);
	vertex_index_prefix_sum_obj.resize(total_obj_num + 1, 0);

	vertex_num_on_surface_prefix_sum[0] = 0;
	for (int i = 1; i <= total_obj_num; ++i) {
		if (i <= cloth->size()) {
			vertex_num_on_surface_prefix_sum[i] = vertex_num_on_surface_prefix_sum[i - 1] + vertex_index_begin_per_thread[i - 1][total_thread_num];
		}
		else {
			vertex_num_on_surface_prefix_sum[i] = vertex_num_on_surface_prefix_sum[i - 1] + tetrahedron->data()[i - 1 - cloth->size()].mesh_struct.vertex_index_on_sureface.size();
		}

		vertex_index_prefix_sum_obj[i] = vertex_index_prefix_sum_obj[i - 1] + mesh_struct[i - 1]->vertex_for_render.size();
	}


	global_vertex_index_start_per_thread.resize(total_thread_num + 1, 0);
	arrangeIndex(total_thread_num, vertex_index_prefix_sum_obj[total_obj_num], global_vertex_index_start_per_thread.data());


	int vertex_num;
	record_vertex_position_every_thread.resize(total_obj_num);
	record_vertex_position_num_every_thread.resize(total_obj_num);

	record_vertex_by_thread.resize(total_thread_num);
	record_vertex_update_num_by_thread.resize(total_thread_num);

	for (int i = 0; i < total_obj_num; ++i) {
		vertex_num = mesh_struct[i]->vertex_position.size();
		record_vertex_position_every_thread[i].resize(total_thread_num);
		record_vertex_position_num_every_thread[i].resize(total_thread_num);
		for (int j = 0; j < total_thread_num; ++j) {
			record_vertex_position_every_thread[i][j].resize(vertex_num);
			record_vertex_position_num_every_thread[i][j].resize(vertex_num);
		}
	}


	for (int i = 0; i < total_thread_num; ++i) {
		record_vertex_by_thread[i].resize(total_obj_num);
		record_vertex_update_num_by_thread[i].resize(total_obj_num);
		for (int j = 0; j < total_obj_num; ++j) {
			record_vertex_by_thread[i][j] = record_vertex_position_every_thread[j][i].data();
			record_vertex_update_num_by_thread[i][j] = record_vertex_position_num_every_thread[j][i].data();
		}

	}

	collision.indicate_if_involved_in_last_color = record_vertex_position_num_every_thread.data();
	//collision_compare.indicate_if_involved_in_last_color = record_vertex_position_num_every_thread.data();
	std::cout << "xpbd ipc reorganzieDataOfObjects " << vertex_num_on_surface_prefix_sum[1] << std::endl;
}

void IPC::saveScene(double* force_direction, int obj_No, bool have_force)
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

void IPC::readScene(const char* file_name, double* force_direction, int& obj_No)
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


void IPC::updateNormal()
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


void IPC::initial()
{
	resetExternalForce();
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
	}
}

void IPC::resetExternalForce()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(f_ext[i][0].data(), 0, 24 * f_ext[i].size());
	}
}

void IPC::reset()
{
	resetExternalForce();
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
	}
}