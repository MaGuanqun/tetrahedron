#include"project_dynamic.h"
#include"basic/EigenMatrixIO.h"

ProjectDynamic::ProjectDynamic()
{
	gravity_ = 9.8;
	total_thread_num = std::thread::hardware_concurrency();
	temEnergy.resize(total_thread_num);
	outer_itr_conv_rate = 1e-3;// 7.5e-2; 
	local_global_conv_rate = 5e-2;
	sub_step_num = 1;

	use_dierct_solve_for_coarest_mesh = true;
	super_jacobi_step_size = 3;
	max_it = 30;
	max_jacobi_itr_num = 20;
	displacement_norm_thread.resize(total_thread_num);


	iteration_method.setConvergenceRate(1e-7, 20);
	max_inner_iteration_num = 5;

}

void ProjectDynamic::setForPD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Thread* thread)
{
	sub_time_step = time_step / (double)sub_step_num;
	this->thread = thread;
	setForClothPD(cloth);
	setForTetrahedronPD(tetrahedron);
	setIndexPerThread();
	collision.initial(cloth, collider, tetrahedron, thread);
	total_collider_num = collider->size();
	this->collider = collider;

	iteration_method.setBasicInfo(total_cloth_num, cloth_sys_size, thread, cloth_per_thread_begin);
	iteration_method.initialGlobalDiagonalInv(&cloth_global_mat_diagonal_ref_address);
	initialJacobi();

	ave_iteration.resize(cloth->size());
	for (int i = 0; i < ave_iteration.size(); ++i) {
		ave_iteration[i].resize(3);
	}
	iteration_method.test();
	
}


void ProjectDynamic::initialDHatTolerance(double ave_edge_length)
{
	collision.initialDHatTolerance(ave_edge_length);
	int element_count = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		element_count += (*cloth)[i].mesh_struct.vertex_for_render.size();
	}
	displacement_bound = 2e-3 * ave_edge_length;
	displacement_bound *= displacement_bound;
	displacement_bound *= (double)element_count;
}

void ProjectDynamic::setForClothPD(std::vector<Cloth>* cloth)
{
	total_cloth_num = cloth->size();
	cloth_sys_size.resize(total_cloth_num);
	for (int i = 0; i < total_cloth_num; ++i) {
		cloth_sys_size[i] = (*cloth)[i].mesh_struct.vertices.size();
	}
	this->cloth = cloth;
	std::vector<std::vector<double>> edge_cot_weight;
	computeEdgeCotWeight(edge_cot_weight);
	computeLBOWeight(edge_cot_weight);
	computeClothGravity();
	computeVertexLBO(edge_cot_weight);
	computeGlobalStepMatrix();
	initialClothPDvariable();
	restBendingMeanCurvature();
}

void ProjectDynamic::setForTetrahedronPD(std::vector<Tetrahedron>* tetrahedron)
{
	total_tetrahedron_num = tetrahedron->size();
	tetrahedron_sys_size.resize(total_tetrahedron_num);
	for (int i = 0; i < total_tetrahedron_num; ++i) {
		tetrahedron_sys_size[i] = (*tetrahedron)[i].mesh_struct.vertex_position.size();
	}

}

void ProjectDynamic::restBendingMeanCurvature()
{
	rest_mean_curvature_norm.resize(total_cloth_num);
	int neighbor, ve;
	double vertex_curvature[3];
	for (int k = 0; k < total_cloth_num; ++k) {
		restBendingMeanCurvatureSingleCloth((*cloth)[k].mesh_struct, rest_mean_curvature_norm[k], vertex_around_vertex_for_bending[k],
			vertex_lbo[k]);
	}
}

void ProjectDynamic::restBendingMeanCurvatureSingleCloth(TriangleMeshStruct& mesh_struct, std::vector<double>& rest_mean_curvature_norm,
	std::vector<std::vector<int>>& vertex_around_vertex_for_bending, std::vector<VectorXd>& vertex_lbo)
{
	rest_mean_curvature_norm.resize(mesh_struct.vertices.size());
	VectorXd q;
	int size;
	double vertex_curvature[3];
	for (int i = 0; i < mesh_struct.vertices.size(); ++i) {
		size = vertex_around_vertex_for_bending[i].size();
		q.resize(size);
		for (int j = 0; j < 3; ++j) {
			for (int h = 0; h < size; h++) {
				q[h] = mesh_struct.vertex_position[vertex_around_vertex_for_bending[i][h]][j];
			}
			vertex_curvature[j] = dotProductX_(vertex_lbo[i], q);
		}
		rest_mean_curvature_norm[i] = sqrt(DOT(vertex_curvature, vertex_curvature));
	}
}


void ProjectDynamic::initialTetrahedronPDvariable()
{

}


void ProjectDynamic::initialClothPDvariable()
{
	p_edge_length.resize(total_cloth_num);
	p_bending.resize(total_cloth_num);
	cloth_u.resize(total_cloth_num);
	cloth_v.resize(total_cloth_num);

	for (int j = 0; j < total_cloth_num; ++j) {
		p_edge_length[j].resize((*cloth)[j].mesh_struct.edges.size());
		p_bending[j].resize(cloth_sys_size[j]);
		cloth_u[j].resize(3);
		cloth_v[j].resize(3);
		for (int i = 0; i < 3; ++i) {
			cloth_u[j][i].resize(cloth_sys_size[j]);	cloth_u[j][i].setZero();
			cloth_v[j][i].resize(cloth_sys_size[j]); cloth_v[j][i].setZero();
			for (int k = 0; k < cloth_sys_size[j]; ++k) {
				cloth_u[j][i].data()[k] = (*cloth)[j].mesh_struct.vertex_position[k][i];
			}
		}
		for (int i = 0; i < cloth_sys_size[j]; ++i) {
			p_bending[j][i].resize(3);
		}
	}

	cloth_u_previous_itr = cloth_u;
	cloth_u_ = cloth_u;
	cloth_v_ = cloth_v;
	cloth_b = cloth_v;
	cloth_u_prediction = cloth_u;
	cloth_acceleration = cloth_v;

}
void ProjectDynamic::setIndexPerThread()
{
	cloth_per_thread_begin.resize(total_thread_num + 1);
	if (total_cloth_num > 0) {
		arrangeIndex(total_thread_num, total_cloth_num, cloth_per_thread_begin);
	}
	cloth_dimension_per_thread.resize(total_thread_num);
	int total_num = total_cloth_num * 3;
	for (int i = 0; i < total_num; ++i) {
		cloth_dimension_per_thread[i % total_thread_num].push_back(i);
	}
	tetrahedron_per_thread_begin.resize(total_thread_num + 1);
	if (total_tetrahedron_num > 0) {
		arrangeIndex(total_thread_num, total_tetrahedron_num, tetrahedron_per_thread_begin);
	}
	tetrahedron_dimension_per_thread.resize(total_thread_num);
	total_num = total_tetrahedron_num * 3;
	for (int i = 0; i < total_num; ++i) {
		tetrahedron_dimension_per_thread[i % total_thread_num].push_back(i);
	}
}

void ProjectDynamic::computeVertexLBO(std::vector<std::vector<double>>& edge_cot_weight)
{
	std::vector<std::vector<std::vector<int>>> edge_around_vertex_for_bending;
	edge_around_vertex_for_bending.resize(total_cloth_num);
	vertex_around_vertex_for_bending.resize(total_cloth_num);
	for (int k = 0; k < total_cloth_num; ++k) {
		edge_around_vertex_for_bending[k].resize(cloth_sys_size[k]);
		vertex_around_vertex_for_bending[k].resize(cloth_sys_size[k]);
		setAroundVertexPrimitive((*cloth)[k].mesh_struct, edge_around_vertex_for_bending[k], vertex_around_vertex_for_bending[k]);
	}
	vertex_lbo.resize(total_cloth_num);
	for (int k = 0; k < total_cloth_num; ++k) {
		computeVertexLBOSingleCloth((*cloth)[k].mesh_struct, vertex_lbo[k], edge_cot_weight[k], edge_around_vertex_for_bending[k], lbo_weight[k]);
	}
}

void ProjectDynamic::computeGlobalStepMatrix()
{
	cloth_global_mat.resize(total_cloth_num);
	cloth_global_mat_diagonal_ref.resize(total_cloth_num);
	cloth_global_mat_diagonal_ref_address.resize(total_cloth_num);
	cloth_llt = new SimplicialLLT<SparseMatrix<double>>[total_cloth_num];
	collision_free_cloth_llt = new SimplicialLLT<SparseMatrix<double>>[total_cloth_num];
	for (int i = 0; i < total_cloth_num; ++i) {
		computeGlobalStepMatrixSingleCloth(&cloth_global_mat[i], cloth_global_mat_diagonal_ref[i], cloth_global_mat_diagonal_ref_address[i],
			&cloth_llt[i], (*cloth)[i].mesh_struct, vertex_around_vertex_for_bending[i], lbo_weight[i], vertex_lbo[i], cloth_sys_size[i],
			(*cloth)[i].bend_stiffness, (*cloth)[i].length_stiffness, (*cloth)[i].position_stiffness, &collision_free_cloth_llt[i]);
	}
	initial_cloth_global_mat = cloth_global_mat;
}


void ProjectDynamic::computeOffDiagonal()
{
	iteration_method.offDiagonalSize();	
	for (int i = 0; i < total_cloth_num; ++i) {
		setOffDiagonal(i, vertex_around_vertex_for_bending[i], lbo_weight[i], vertex_lbo[i], cloth_sys_size[i],
			(*cloth)[i].bend_stiffness, (*cloth)[i].length_stiffness, (*cloth)[i].mesh_struct);
	}
}


void ProjectDynamic::setOffDiagonal(int cloth_No, std::vector<std::vector<int>>& vertex_around_vertex_for_bending,
	std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo, int sys_size, double bending_stiffness, std::vector<double>& length_stiffness,
	TriangleMeshStruct& mesh_struct)
{
	std::vector<Triplet<double>>global_mat_nnz;
	int estimate_total_constraint = 4 * mesh_struct.edges.size() + 66 * sys_size; //to estimate the size of global_mat_nnz
	global_mat_nnz.reserve(estimate_total_constraint);
	//bending
	int lbo_length; MatrixXd ATA;
	if (!mesh_struct.faces.empty()) {
		for (int i = 0; i < sys_size; ++i) {
			lbo_length = vertex_lbo[i].size();
			ATA = (bending_stiffness * lbo_weight[i] * vertex_lbo[i]) * vertex_lbo[i].transpose();
			for (int l = 0; l < lbo_length; ++l) {
				for (int k = l + 1; k < lbo_length; ++k) {
					global_mat_nnz.push_back(Triplet<double>(vertex_around_vertex_for_bending[i][l], vertex_around_vertex_for_bending[i][k], ATA.data()[lbo_length * l + k]));
					global_mat_nnz.push_back(Triplet<double>(vertex_around_vertex_for_bending[i][k], vertex_around_vertex_for_bending[i][l], ATA.data()[lbo_length * k + l]));
				}
			}
		}
	}
	//edge length
	int id0, id1;
	for (int i = 0; i < mesh_struct.edges.size(); ++i) {
		id0 = mesh_struct.edges[i].vertex[0];
		id1 = mesh_struct.edges[i].vertex[1];
		global_mat_nnz.push_back(Triplet<double>(id0, id1, -length_stiffness[i]));
		global_mat_nnz.push_back(Triplet<double>(id1, id0, -length_stiffness[i]));
	}
	iteration_method.setOffDiagonal(cloth_No, global_mat_nnz);
}


void ProjectDynamic::computeGlobalStepMatrixSingleCloth(SparseMatrix<double, RowMajor>* global_mat, std::vector<double>& global_mat_diagonal_ref,
	std::vector<double*>& global_collision_mat_diagonal_ref, SimplicialLLT<SparseMatrix<double>>* global_llt, TriangleMeshStruct& mesh_struct,
	std::vector<std::vector<int>>& vertex_around_vertex_for_bending,
	std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo, int sys_size, double bending_stiffness, std::vector<double>& length_stiffness,
	double position_stiffness, SimplicialLLT<SparseMatrix<double>>* collision_free_global_llt)
{
	std::vector<Triplet<double>>global_mat_nnz;
	global_mat->resize(sys_size, sys_size);
	int estimate_total_constraint = 4 * mesh_struct.edges.size() + 66 * sys_size; //to estimate the size of global_mat_nnz
	global_mat_nnz.reserve(estimate_total_constraint);
	//bending
	int lbo_length; MatrixXd ATA;
	if (!mesh_struct.faces.empty()) {
		for (int i = 0; i < sys_size; ++i) {
			lbo_length = vertex_lbo[i].size();
			ATA = (bending_stiffness * lbo_weight[i] * vertex_lbo[i]) * vertex_lbo[i].transpose();
			for (int l = 0; l < lbo_length; ++l) {
				for (int k = l + 1; k < lbo_length; ++k) {
					global_mat_nnz.push_back(Triplet<double>(vertex_around_vertex_for_bending[i][l], vertex_around_vertex_for_bending[i][k], ATA.data()[lbo_length * l + k]));
					global_mat_nnz.push_back(Triplet<double>(vertex_around_vertex_for_bending[i][k], vertex_around_vertex_for_bending[i][l], ATA.data()[lbo_length * k + l]));
				}
				global_mat_nnz.push_back(Triplet<double>(vertex_around_vertex_for_bending[i][l], vertex_around_vertex_for_bending[i][l], ATA.data()[lbo_length * l + l]));
			}
		}
	}
	//edge length
	int id0, id1;
	for (int i = 0; i < mesh_struct.edges.size(); ++i) {
		id0 = mesh_struct.edges[i].vertex[0];
		id1 = mesh_struct.edges[i].vertex[1];

		global_mat_nnz.push_back(Triplet<double>(id0, id0, length_stiffness[i]));
		global_mat_nnz.push_back(Triplet<double>(id1, id1, length_stiffness[i]));
		global_mat_nnz.push_back(Triplet<double>(id0, id1, -length_stiffness[i]));
		global_mat_nnz.push_back(Triplet<double>(id1, id0, -length_stiffness[i]));
	}

	//position
	for (int i = 0; i < mesh_struct.anchor_vertex.size(); ++i) {
		global_mat_nnz.push_back(Triplet<double>(mesh_struct.anchor_vertex[i], mesh_struct.anchor_vertex[i], position_stiffness));
	}
	//mass
	for (int i = 0; i < sys_size; ++i) {
		global_mat_nnz.push_back(Triplet<double>(i, i, mesh_struct.mass[i] / (sub_time_step * sub_time_step)));
	}
	global_mat->setFromTriplets(global_mat_nnz.begin(), global_mat_nnz.end());
	global_collision_mat_diagonal_ref.resize(sys_size);
	global_mat_diagonal_ref.resize(sys_size);
	for (int i = 0; i < sys_size; ++i) {
		global_collision_mat_diagonal_ref[i] = &global_mat->coeffRef(i, i);
	}
	for (int i = 0; i < sys_size; ++i) {
		global_mat_diagonal_ref[i] = *(global_collision_mat_diagonal_ref[i]);
	}
	global_llt->analyzePattern((*global_mat));
	global_llt->factorize((*global_mat));
	collision_free_global_llt->compute((*global_mat));
}

void ProjectDynamic::computeVertexLBOSingleCloth(TriangleMeshStruct& mesh_struct, std::vector<VectorXd>& vertex_lbo, std::vector<double>& edge_cot_weight,
	std::vector<std::vector<int>>& edge_around_vertex_for_bending, std::vector<double>& lbo_weight)
{
	double total;
	double edge_weight;
	vertex_lbo.resize(mesh_struct.vertices.size());
	for (int i = 0; i < mesh_struct.vertices.size(); ++i) {
		vertex_lbo[i].resize(edge_around_vertex_for_bending[i].size() + 1);
		vertex_lbo[i].setZero();
		total = 0.0;
		for (int j = 0; j < edge_around_vertex_for_bending[i].size(); ++j) {
			edge_weight = edge_cot_weight[edge_around_vertex_for_bending[i][j]] / lbo_weight[i];
			total += edge_weight;
			vertex_lbo[i].data()[j + 1] = edge_weight;
		}
		vertex_lbo[i].data()[0] = -total;
	}

}

void ProjectDynamic::setAroundVertexPrimitive(TriangleMeshStruct& mesh_struct, std::vector<std::vector<int>>& edge_around_vertex_for_bending,
	std::vector<std::vector<int>>& vertex_around_vertex_for_bending)
{
	std::vector<TriangleMeshStruct::Vertex>* vertex;
	std::vector<int>* around_edge_index;
	for (int i = 0; i < mesh_struct.vertices.size(); ++i) {
		around_edge_index = &mesh_struct.vertices[i].edge;
		edge_around_vertex_for_bending[i].reserve(around_edge_index->size());
		vertex_around_vertex_for_bending[i].reserve(around_edge_index->size() + 1);
		vertex_around_vertex_for_bending[i].push_back(i);
		for (int j = 0; j < around_edge_index->size(); ++j) {
			if (mesh_struct.edges[(*around_edge_index)[j]].opposite_vertex.size() == 2) {
				edge_around_vertex_for_bending[i].push_back((*around_edge_index)[j]);
				if (i == mesh_struct.edges[(*around_edge_index)[j]].vertex[0]) {
					vertex_around_vertex_for_bending[i].push_back(mesh_struct.edges[(*around_edge_index)[j]].vertex[1]);
				}
				else {
					vertex_around_vertex_for_bending[i].push_back(mesh_struct.edges[(*around_edge_index)[j]].vertex[0]);
				}
			}
		}
	}
}

void ProjectDynamic::computeClothGravity()
{
	cloth_f_ext.resize(total_cloth_num);
	cloth_gravity.resize(total_cloth_num);
	double gravity_accerlation[3] = { 0, -gravity_,0 };
	std::vector<double>* mass_;
	for (int j = 0; j < total_cloth_num; ++j) {
		cloth_f_ext[j].resize(3);
		cloth_gravity[j].resize(3);
		mass_ = &(*cloth)[j].mesh_struct.mass;
		for (int i = 0; i < 3; ++i) {
			if (cloth_sys_size[j] > 1) {
				cloth_gravity[j][i].resize(cloth_sys_size[j]);
				for (int k = 0; k < cloth_sys_size[j]; ++k) {
					cloth_gravity[j][i][k] = gravity_accerlation[i] * (*mass_)[k];
				}
			}
			else {
				cloth_gravity[j][i] = VectorXd::Zero(cloth_sys_size[j]);
			}
			cloth_f_ext[j][i] = cloth_gravity[j][i];
		}
	}
}



void ProjectDynamic::computeLBOWeightSingleCloth(std::vector<double>& edge_cot_weight, std::vector<double>& lbo_weight,
	TriangleMeshStruct& mesh_struct, VectorXd& mass_inv, double density, VectorXd& mass)
{
	double m;
	for (int i = 0; i < mesh_struct.faces.size(); ++i) {
		m = mesh_struct.faces[i].area / 3.0;
		lbo_weight[mesh_struct.triangle_indices[i][0]] += m;
		lbo_weight[mesh_struct.triangle_indices[i][1]] += m;
		lbo_weight[mesh_struct.triangle_indices[i][2]] += m;
	}
	double mass_;
	if (!mesh_struct.faces.empty()) {
		for (int i = 0; i < mass_inv.size(); ++i) {
			mass_ = density * lbo_weight[i];
			mass_inv.data()[i] = 1.0 / mass_;
			mass.data()[i] = mass_;
			mesh_struct.mass[i] = mass_;
		}
	}
	else {
		for (int i = 0; i < mass_inv.size(); ++i) {
			mass_ = 1.25;
			mass_inv.data()[i] = 1.0 / mass_;
			mass.data()[i] = mass_;
			mesh_struct.mass[i] = mass_;
		}
	}
}

void ProjectDynamic::computeLBOWeight(std::vector<std::vector<double>>& edge_cot_weight)
{

	lbo_weight.resize(total_cloth_num);
	cloth_mass_inv.resize(total_cloth_num);
	cloth_mass.resize(total_cloth_num);
	for (int j = 0; j < total_cloth_num; ++j) {
		lbo_weight[j].resize(cloth_sys_size[j], 0.0);
		cloth_mass_inv[j].resize(cloth_sys_size[j]);
		cloth_mass[j].resize(cloth_sys_size[j]);
		computeLBOWeightSingleCloth(edge_cot_weight[j], lbo_weight[j], (*cloth)[j].mesh_struct, cloth_mass_inv[j], (*cloth)[j].density, cloth_mass[j]);
	}
}

void ProjectDynamic::computeEdgeCotWeight(std::vector<std::vector<double>>& edge_cot_weight)
{
	edge_cot_weight.resize(total_cloth_num);
	double cotan0, cotan1;
	double len10, len20, len13, len23;
	double x10[3], x20[3];
	double x13[3], x23[3];
	double theta0, theta1;
	for (int j = 0; j < total_cloth_num; ++j) {
		computeEdgeCotWeightSingleCloth(edge_cot_weight[j], (*cloth)[j].mesh_struct);
	}
}

void ProjectDynamic::computeEdgeCotWeightSingleCloth(std::vector<double>& edge_cot_weight, TriangleMeshStruct& mesh_struct)
{
	int edge_num;
	edge_num = mesh_struct.edges.size();
	edge_cot_weight.resize(edge_num);
	double cotan0, cotan1;
	double len10, len20, len13, len23;
	double x10[3], x20[3];
	double x13[3], x23[3];
	double theta0, theta1;
	for (int i = 0; i < edge_num; ++i) {
		if (mesh_struct.edges[i].opposite_vertex.size() > 1) {
			cotan0 = 0;
			cotan1 = 0;
			SUB(x10, mesh_struct.vertex_position[mesh_struct.edges[i].vertex[0]], mesh_struct.vertex_position[mesh_struct.edges[i].opposite_vertex[0]]);
			SUB(x20, mesh_struct.vertex_position[mesh_struct.edges[i].vertex[1]], mesh_struct.vertex_position[mesh_struct.edges[i].opposite_vertex[0]]);
			len10 = sqrt(DOT(x10, x10));
			len20 = sqrt(DOT(x20, x20));
			theta0 = acos(DOT(x10, x20) / (len10 * len20));
			cotan0 = 1.0 / tan(theta0);

			SUB(x13, mesh_struct.vertex_position[mesh_struct.edges[i].vertex[0]], mesh_struct.vertex_position[mesh_struct.edges[i].opposite_vertex[1]]);
			SUB(x23, mesh_struct.vertex_position[mesh_struct.edges[i].vertex[1]], mesh_struct.vertex_position[mesh_struct.edges[i].opposite_vertex[1]]);
			len13 = sqrt(DOT(x13, x13));
			len23 = sqrt(DOT(x23, x23));
			theta1 = acos(DOT(x13, x23) / (len13 * len23));
			cotan1 = 1.0 / tan(theta1);
			edge_cot_weight[i] = -0.5 * (cotan0 + cotan1);
		}
	}
}

void ProjectDynamic::reset()
{
	for (int j = 0; j < total_cloth_num; ++j) {
		for (int i = 0; i < cloth_sys_size[j]; ++i) {
			for (int k = 0; k < 3; ++k) {
				cloth_u[j][k].data()[i] = (*cloth)[j].ori_vertices[i][k];
			}
		}
		for (int i = 0; i < 3; ++i) {
			cloth_u_[j][i] = cloth_u[j][i];
			cloth_v[j][i].setZero();
			cloth_v_[j][i].setZero();
			cloth_acceleration[j][i].setZero();
		}
	}
	cloth_f_ext = cloth_gravity;

	for (int j = 0; j < total_tetrahedron_num; ++j) {
		for (int i = 0; i < tetrahedron_sys_size[j]; ++i) {
			for (int k = 0; k < 3; ++k) {
				tetrahedron_u[j][k].data()[i] = (*tetrahedron)[j].ori_vertices[i][k];
			}
		}
		for (int i = 0; i < 3; ++i) {
			tetrahedron_u_[j][i] = tetrahedron_u[j][i];
			tetrahedron_v[j][i].setZero();
			tetrahedron_v_[j][i].setZero();
			tetrahedron_acceleration[j][i].setZero();
		}
	}
	tetrahedron_f_ext = tetrahedron_gravity;
}

void ProjectDynamic::initial()
{
	cloth_global_mat = initial_cloth_global_mat;
	for (int j = 0; j < total_cloth_num; ++j) {
		for (int i = 0; i < cloth_sys_size[j]; ++i) {
			cloth_global_mat_diagonal_ref_address[j][i] = &cloth_global_mat[j].coeffRef(i, i);
			cloth_global_mat_diagonal_ref[j][i] = *cloth_global_mat_diagonal_ref_address[j][i];
		}
	}


	tetrahedron_global_mat = initial_tetrahedron_global_mat;

	reset();
}

void ProjectDynamic::updateModelPosition()
{
	std::array<double, 3>* vertex_position;
	for (int j = 0; j < total_cloth_num; ++j) {
		vertex_position = (*cloth)[j].mesh_struct.vertex_position.data();
		for (int i = 0; i < cloth_sys_size[j]; ++i) {
			vertex_position[i][0] = cloth_u[j][0].data()[i];
			vertex_position[i][1] = cloth_u[j][1].data()[i];
			vertex_position[i][2] = cloth_u[j][2].data()[i];
		}
	}
	for (int j = 0; j < total_tetrahedron_num; ++j) {
		vertex_position = (*tetrahedron)[j].mesh_struct.vertex_position.data();
		for (int i = 0; i < tetrahedron_sys_size[j]; ++i) {
			vertex_position[i][0] = tetrahedron_u[j][0].data()[i];
			vertex_position[i][1] = tetrahedron_u[j][1].data()[i];
			vertex_position[i][2] = tetrahedron_u[j][2].data()[i];
		}
	}
}


void ProjectDynamic::updateRenderPositionIPC()
{
	TriangleMeshStruct* mesh_struct;
	for (int j = 0; j < total_cloth_num; ++j) {
		mesh_struct = &(*cloth)[j].mesh_struct;
		thread->assignTask(mesh_struct, FACE_NORMAL_RENDER);
		thread->assignTask(mesh_struct, VERTEX_NORMAL_RENDER);
	}
	for (int j = 0; j < total_tetrahedron_num; ++j) {
		(*tetrahedron)[j].mesh_struct.vertex_for_render = (*tetrahedron)[j].mesh_struct.vertex_position;
		(*tetrahedron)[j].mesh_struct.getRenderNormal();
	}
	for (int j = 0; j < total_collider_num; ++j) {
		mesh_struct = &(*collider)[j].mesh_struct;
		thread->assignTask(mesh_struct, FACE_NORMAL_RENDER);
		thread->assignTask(mesh_struct, VERTEX_NORMAL_RENDER);
		mesh_struct->ori_face_normal_for_render = mesh_struct->ori_face_normal;
	}
}

void ProjectDynamic::updateRenderPosition()
{
	TriangleMeshStruct* mesh_struct;
	for (int j = 0; j < total_cloth_num; ++j) {
		mesh_struct = &(*cloth)[j].mesh_struct;
		mesh_struct->vertex_for_render = mesh_struct->vertex_position;
		mesh_struct->face_normal_for_render = mesh_struct->face_normal;
		thread->assignTask(mesh_struct, VERTEX_NORMAL_RENDER);
	}
	for (int j = 0; j < total_tetrahedron_num; ++j) {
		(*tetrahedron)[j].mesh_struct.vertex_for_render = (*tetrahedron)[j].mesh_struct.vertex_position;
		(*tetrahedron)[j].mesh_struct.getRenderNormal();
	}

	for (int j = 0; j < total_collider_num; ++j) {
		mesh_struct = &(*collider)[j].mesh_struct;
		mesh_struct->vertex_for_render = mesh_struct->vertex_position;
		mesh_struct->face_normal_for_render = mesh_struct->face_normal;
		thread->assignTask(mesh_struct, VERTEX_NORMAL_RENDER);
	}
}

void ProjectDynamic::PDsetPosPredict()
{
	PDClothPredict();
	PDTetrahedronPredict();
	updateModelPosition();
}



void ProjectDynamic::PDTetrahedronPredict()
{

}

void ProjectDynamic::PDClothPredict()
{
	for (int j = 0; j < total_cloth_num; ++j) {
		for (int i = 0; i < 3; ++i) {
			cloth_u[j][i] = cloth_u_[j][i] + sub_time_step * cloth_v_[j][i] + (0.25 * sub_time_step * sub_time_step) * (cloth_acceleration[j][i]);//sub_time_step * sub_time_step * mass_inv[j].cwiseProduct(f_ext[j][i]);
			cloth_u_prediction[j][i] = cloth_u_[j][i] + sub_time_step * cloth_v_[j][i] + (sub_time_step * sub_time_step) * cloth_mass_inv[j].cwiseProduct(cloth_f_ext[j][i]);
			//cloth_u[j][i] = cloth_u_prediction[j][i];
		}
	}
}


void ProjectDynamic::firstPDForIPC()
{
	//face_normal_render
	PDsetPosPredict();
	thread->assignTask(this, LOCAL_PROJECTION);

	current_constraint_energy = temEnergy[0];
	for (int i = 1; i < total_thread_num; ++i) {
		current_constraint_energy += temEnergy[i];
	}

	thread->assignTask(this, SOLVE_SYSYTEM_WITHOUT_COLLISION);
	updateModelPosition();
	collision.collisionCulling();

	current_PD_energy = temEnergy[0];
	for (int i = 1; i < total_thread_num; ++i) {
		current_PD_energy += temEnergy[i];
	}
	current_PD_energy += current_constraint_energy;
	current_collision_energy = 1e-15;
	previous_PD_energy = 1e-15;
	////std::cout << "++++" << std::endl;
	//for (int i = 0; i < cloth_sys_size[0]; ++i) {
	//	//std::cout << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "
	//		<< cloth_u[0][2][i] << std::endl;
	//}
	//std::cout << "first local-global without collision" << std::endl;
	for (int i = 0; i < cloth_sys_size[0]; ++i) {
		//std::cout << "    " << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "<< cloth_u[0][2][i] << std::endl;
	}
}


void ProjectDynamic::PD_IPC_solve(bool& record_matrix)
{
	firstPDForIPC();
	outer_iteration_num = 1;
	reset_ave_iteration_record();
	while (!IPC_PDConvergeCondition()) {
		////std::cout << "==" << std::endl;
		collision.globalCollisionTime();
		thread->assignTask(this, COLLISION_FREE_POSITION);//in document, we use q_n+1, however, here, we use vertices_for_render & cloth_u to store this collision free position.
		cloth_u_previous_itr = cloth_u;
		collision.solveCollisionConstraint();
		PDupdateSystemMatrix();
		//std::cout << "==iteration number " << outer_iteration_num << std::endl;
		//std::cout << "collision free position " << std::endl;
		//for (int i = 0; i < cloth_sys_size[0]; ++i) {
		//	std::cout<<"    " << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "<< cloth_u[0][2][i] << std::endl;
		//}
		local_global_iteration_num = 0;

		if (record_matrix) {
			saveSparseMatrix(cloth_global_mat[0],"global.dat");
			saveSparseMatrix(iteration_method.off_diagonal[0],"off_diagonal.dat");

			for (int i = 0; i < 3; ++i) {
				saveMatrix(i, cloth_u[0][i], "u");
			}			
		}

		while (!innerIterationConvergeCondition()) {
			thread->assignTask(this, LOCAL_PROJECTION_WITHOUT_ENERGY);
			current_constraint_energy = temEnergy[0];
			for (int i = 1; i < total_thread_num; ++i) {
				current_constraint_energy += temEnergy[i];
			}
			thread->assignTask(this, SOLVE_SYSYTEM_WITHOUT_ENERGY);
			local_global_iteration_num++;
			computeInnerEnergyIPCPD();

			if (record_matrix) {
				for (int i = 0; i < 3; ++i) {
					saveMatrix(i, cloth_b[0][i], "b");
				}
				record_matrix = false;
			}
		}
		//std::cout << outer_iteration_num << std::endl;
		updateModelPosition();
		outer_iteration_num++;
		computeEnergyIPCPD();
		displacement_ratio_dif = previous_displacement_norm - displacement_norm;
		previous_displacement_norm = displacement_norm;
		//system("pause");
		//std::cout << "pd position " << std::endl;
		//for (int i = 0; i < cloth_sys_size[0]; ++i) {
		//	std::cout << "    " << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "	<< cloth_u[0][2][i] << std::endl;
		//}
		//std::cout << "+++++" << std::endl;
	}
	collision.globalCollisionTime();
	thread->assignTask(this, COLLISION_FREE_POSITION);
	//std::cout << "final collision free position " << std::endl;
	//for (int i = 0; i < cloth_sys_size[0]; ++i) {
	//	std::cout << "    " << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "	<< cloth_u[0][2][i] << std::endl;
	//}
	thread->assignTask(this, UPDATE_UV);
	updateRenderPositionIPC();
	//std::cout << cloth_v[0][1] << std::endl;
	//std::cout << "========" << std::endl;
}


void ProjectDynamic::computeInnerEnergyIPCPD()
{
	previous_itr_PD_energy = current_PD_energy;
	current_collision_energy = 1e-15;
	for (int k = 0; k < total_cloth_num; ++k) {
		current_collision_energy += collision.cloth_target_pos.collision_energy[k];
	}
	current_PD_energy = temEnergy[0];
	for (int i = 1; i < total_thread_num; ++i) {
		current_PD_energy += temEnergy[i];
	}
	current_PD_energy += current_constraint_energy + current_constraint_energy;
}

void ProjectDynamic::computeEnergyIPCPD()
{
	PD_energy_dif = previous_PD_energy - current_PD_energy;
	collision_energy_dif = previous_collision_energy - current_collision_energy;
	previous_collision_energy = current_collision_energy;
	previous_PD_energy = current_PD_energy;
	current_collision_energy = 1e-15;
	for (int k = 0; k < total_cloth_num; ++k) {
		current_collision_energy += collision.cloth_target_pos.collision_energy[k];
	}
	current_PD_energy = temEnergy[0];
	for (int i = 1; i < total_thread_num; ++i) {
		current_PD_energy += temEnergy[i];
	}
	current_PD_energy += current_constraint_energy + current_constraint_energy;

}

void ProjectDynamic::PDsolve()
{
	PDsetPosPredict();
	int itr_num = 0;
	outer_iteration_num = 0;
	initialEnergy();
	local_global_iteration_num = 0;
	for (int i = 0; i < collider->size(); ++i) {
		thread->assignTask(&(*collider)[i].mesh_struct, FACE_NORMAL);
	}
	for (int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	}
	////std::cout << "==============================================" << std::endl;
	while (!PDConvergeCondition()) {
		collision.globalCollision();
		PDupdateSystemMatrix();
		initialEnergyOuterInteration();
		local_global_itr_in_single_outer = 0;
		while (!PDLocalGlobalConvergeCondition()) {
			initialEnergyLocalGlobal();
			if (local_global_itr_in_single_outer > 0) {
				//time_t t = clock();
				collision.updateCollisionPosition();
				////std::cout << clock() - t << std::endl;
			}
			time_t t = clock();
			//for (int i = 0; i < 1000; ++i) {
				//thread->assignTask(this, LOCAL_EDGE_LENGTH_PROJECTION);
			localProjection();
			for (int i = 0; i < total_thread_num; ++i) {
				current_constraint_energy += temEnergy[i];
			}
			//}
			////std::cout << "local " << clock() - t << std::endl;
			//t = clock();
			//for (int i = 0; i < 1000; ++i) {
				////std::cout << "===" << std::endl;
				//for (int i = 0; i < cloth_sys_size[0]; ++i) {
				//	//std::cout << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "
				//		<< cloth_u[0][2][i] << std::endl;
				//}
			thread->assignTask(this, SOLVE_SYSYTEM);//solve b	
			//for (int i = 0; i < cloth_sys_size[0]; ++i) {
			//	//std::cout << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "
			//		<< cloth_u[0][2][i] << std::endl;
			//}
		//}
		////std::cout << "global " << clock() - t << std::endl;


			for (int k = 0; k < total_cloth_num; ++k) {
				current_collision_energy += collision.cloth_target_pos.collision_energy[k];
			}
			current_constraint_energy += current_collision_energy;
			for (int i = 0; i < total_thread_num; ++i) {
				current_PD_energy += temEnergy[i];
			}
			current_PD_energy += current_constraint_energy;
			updateModelPosition();
			for (int i = 0; i < cloth->size(); ++i) {
				thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
			}
			////std::cout << "result"<< local_global_itr_in_single_outer << std::endl;
			//for (int j = 0; j < total_cloth_num; ++j) {
			//	for (int i = 0; i < cloth_sys_size[j]; ++i) {
			//		//std::cout << cloth_u[j][0].data()[i] << " " << cloth_u[j][1].data()[i] << " " << cloth_u[j][2].data()[i] << std::endl;
			//	}
			//}
			////std::cout << "=======" << std::endl;
			local_global_itr_in_single_outer++;
			local_global_iteration_num++;
			//std::cout << "==iteration number " << local_global_iteration_num << std::endl;
			//std::cout << "pd position " << std::endl;
			for (int i = 0; i < cloth_sys_size[0]; ++i) {
				//std::cout << "    " << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "	<< cloth_u[0][2][i] << std::endl;
			}
		}
		//std::cout << "==local===global=========" << std::endl;
		outer_iteration_num++;

	}
	thread->assignTask(this, UPDATE_UV);
	updateRenderPosition();
}


//UPDATE_UV
void ProjectDynamic::updateUVPerThread(int thread_id)
{
	updateClothUV(thread_id);
	updateTetrahedronUV(thread_id);
}

void ProjectDynamic::updateClothUV(int thread_id)
{
	int object_No;
	int dimension;
	for (int i = 0; i < cloth_dimension_per_thread[thread_id].size(); ++i) {
		object_No = cloth_dimension_per_thread[thread_id][i] / 3;
		dimension = cloth_dimension_per_thread[thread_id][i] % 3;
		cloth_v[object_No][dimension] = (cloth_u[object_No][dimension] - cloth_u_[object_No][dimension]) / sub_time_step;
		cloth_v[object_No][dimension] *= 0.98;
		cloth_acceleration[object_No][dimension] = (cloth_v[object_No][dimension] - cloth_v_[object_No][dimension]) / sub_time_step;
		cloth_u_[object_No][dimension] = cloth_u[object_No][dimension];
		cloth_v_[object_No][dimension] = cloth_v[object_No][dimension];
	}
}


void ProjectDynamic::updateTetrahedronUV(int thread_id)
{
	int object_No;
	int dimension;
	for (int i = 0; i < tetrahedron_dimension_per_thread[thread_id].size(); ++i) {
		object_No = tetrahedron_dimension_per_thread[thread_id][i] / 3;
		dimension = tetrahedron_dimension_per_thread[thread_id][i] % 3;
		tetrahedron_v[object_No][dimension] = (tetrahedron_u[object_No][dimension] - tetrahedron_u_[object_No][dimension]) / sub_time_step;
		tetrahedron_v[object_No][dimension] *= 0.98;
		tetrahedron_acceleration[object_No][dimension] = (tetrahedron_v[object_No][dimension] - tetrahedron_v_[object_No][dimension]) / sub_time_step;
		tetrahedron_u_[object_No][dimension] = tetrahedron_u[object_No][dimension];
		tetrahedron_v_[object_No][dimension] = tetrahedron_v[object_No][dimension];
	}
}




void ProjectDynamic::initialEnergy()
{
	previous_itr_PD_energy = 1e-15;
	previous_itr_constraint_energy = 1e-15;
	previous_itr_collision_energy = 1e-15;
	current_PD_energy = 1e-15;
	current_collision_energy = 1e-15;
	current_constraint_energy = 1e-15;
}

void ProjectDynamic::initialEnergyOuterInteration()
{
	previous_itr_PD_energy = current_PD_energy;
	previous_itr_collision_energy = current_collision_energy;
	previous_itr_constraint_energy = current_constraint_energy;
	previous_PD_energy = 1e-15;
	previous_collision_energy = 1e-15;
	previous_constraint_energy = 1e-15;
}

void ProjectDynamic::initialEnergyLocalGlobal()
{
	previous_PD_energy = current_PD_energy;
	previous_collision_energy = current_collision_energy;
	previous_constraint_energy = current_constraint_energy;
	current_PD_energy = 1e-15;
	current_collision_energy = 1e-15;
	current_constraint_energy = 1e-15;
}

void ProjectDynamic::PDupdateSystemMatrix()
{
	//updateMatrix();
	thread->assignTask(this, UPDATE_MATRIX);
	
	switch (itr_solver_method)
	{
	case DIRECT_SOLVE:
		thread->assignTask(this, MATRIX_DECOMPOSITION);
		break;
	case JACOBI:
		thread->assignTask(&iteration_method, UPDATE_JACOBI_R);
		break;
	case SUPER_JACOBI:
		thread->assignTask(&iteration_method, UPDATE_JACOBI_R);
		break;
	case CHEBYSHEV_SUPER_JACOBI:
		thread->assignTask(&iteration_method, UPDATE_JACOBI_R);
		iteration_method.estimateSuperJacobiEigenValue(cloth_u);
		break;
	case GAUSS_SEIDEL_CHEBYSHEV:
		iteration_method.estimateGaussSeidelEigenValue(cloth_u, cloth_global_mat);
		break;
	case CHEBYSHEV_JACOBI:
		thread->assignTask(&iteration_method, UPDATE_JACOBI_R);
		iteration_method.estimateJacobiEigenValue(cloth_u);
		break;
	case PCG:
		iteration_method.updateGlobalDiagonalInv();
		break;
	case WEIGHTED_JACOBI:
		thread->assignTask(&iteration_method, UPDATE_JACOBI_R);
		break;
	}

}

void ProjectDynamic::updateMatrix()
{
	thread->assignTask(this, UPDATE_MATRIX);
	//updateCollisionMatrix();
}

void ProjectDynamic::updateCollisionMatrix()
{
	std::vector<int>* anchor_index;
	for (int i = 0; i < total_cloth_num; ++i) {
		anchor_index = &(*cloth)[i].mesh_struct.anchor_vertex;
		for (int j = 0; j < anchor_index->size(); ++j) {
			*(cloth_global_mat_diagonal_ref_address[i][(*anchor_index)[j]]) += (*cloth)[i].position_stiffness;
		}
	}
}


//COLLISION_FREE_POSITION
void ProjectDynamic::computeCollisionFreePosition(int thread_No)
{
	int index_end;
	MeshStruct* mesh_struct;
	double collision_time = collision.collision_time;
	std::array<double, 3>* q_pre;
	std::array<double, 3>* q_end;
	double* u_x; double* u_y; double* u_z;
	for (int i = 0; i < cloth->size(); ++i) {
		mesh_struct = &(*cloth)[i].mesh_struct;
		index_end = mesh_struct->vertex_index_begin_per_thread[thread_No + 1];
		q_end = mesh_struct->vertex_position.data();
		q_pre = mesh_struct->vertex_for_render.data();
		u_x = cloth_u[i][0].data();
		u_y = cloth_u[i][1].data();
		u_z = cloth_u[i][2].data();
		for (int j = mesh_struct->vertex_index_begin_per_thread[thread_No]; j < index_end; ++j) {
			q_pre[j][0] += collision_time * (q_end[j][0] - q_pre[j][0]);
			q_pre[j][1] += collision_time * (q_end[j][1] - q_pre[j][1]);
			q_pre[j][2] += collision_time * (q_end[j][2] - q_pre[j][2]);
			u_x[j] = q_pre[j][0];
			u_y[j] = q_pre[j][1];
			u_z[j] = q_pre[j][2];
		}
	}
}


//COMPUTE_DISPLACEMENT
void ProjectDynamic::computeDisplacement(int thread_No)
{
	int index_end;
	double* u_x; double* u_y; double* u_z;
	double* u_pre_x; double* u_pre_y; double* u_pre_z;
	double x, y, z;
	double* displacement_norm_ = &displacement_norm_thread[thread_No];
	*displacement_norm_ = 0;
	double displace_current;
	for (int i = 0; i < cloth->size(); ++i) {
		index_end = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
		u_x = cloth_u[i][0].data();
		u_y = cloth_u[i][1].data();
		u_z = cloth_u[i][2].data();
		u_pre_x = cloth_u_previous_itr[i][0].data();
		u_pre_y = cloth_u_previous_itr[i][1].data();
		u_pre_z = cloth_u_previous_itr[i][2].data();
		for (int j = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread[thread_No]; j < index_end; ++j) {
			x = u_x[j] - u_pre_x[j];
			y = u_y[j] - u_pre_y[j];
			z = u_z[j] - u_pre_z[j];
			displace_current = x * x + y * y + z * z;
			//if (*displacement_norm_ < displace_current) {
				*displacement_norm_ += displace_current;
			//}
		}
	}
}






//UPDATE_MATRIX
void ProjectDynamic::updateMatrixPerThread(int thread_No)
{
	setClothMatrix(thread_No);
}

void ProjectDynamic::setClothMatrix(int thread_No)
{
	int collision_num = 0;
	std::vector<int>* vertex_index_begin_per_thread;
	bool* need_update;
	double* stiffness;
	double* diagonal_ref;
	double** diagonal_ref_address;
	for (int j = 0; j < total_cloth_num; ++j) {
		need_update = collision.cloth_target_pos.need_update[j];
		stiffness = collision.cloth_target_pos.stiffness[j].data();
		vertex_index_begin_per_thread = &(*cloth)[j].mesh_struct.vertex_index_begin_per_thread;
		diagonal_ref = cloth_global_mat_diagonal_ref[j].data();
		diagonal_ref_address = cloth_global_mat_diagonal_ref_address[j].data();
		for (int i = (*vertex_index_begin_per_thread)[thread_No]; i < (*vertex_index_begin_per_thread)[thread_No + 1]; ++i) {
			*(diagonal_ref_address[i]) = diagonal_ref[i];
			////std::cout << i << " " << *(diagonal_ref_address[i]) << std::endl;
			if (need_update[i]) {
				*(diagonal_ref_address[i]) += stiffness[i];
				////std::cout<<"update "<<i << " " << *(diagonal_ref_address[i]) << std::endl;
			}
		}
	}


}

//MATRIX_DECOMPOSITION
void ProjectDynamic::matrixDecomposition(int thread_id)
{

	for (int i = cloth_per_thread_begin[thread_id]; i < cloth_per_thread_begin[thread_id + 1]; ++i) {
		cloth_llt[i].factorize(cloth_global_mat[i]);
	}
	for (int i = tetrahedron_per_thread_begin[thread_id]; i < tetrahedron_per_thread_begin[thread_id + 1]; ++i) {
		tetrahedron_llt[i].factorize(cloth_global_mat[i]);
	}
}




bool ProjectDynamic::innerIterationConvergeCondition()
{
	return local_global_iteration_num > max_inner_iteration_num;

	//if (local_global_iteration_num > 0) {
	//	bool energy_changing = fabs(current_PD_energy - previous_PD_energy) / previous_PD_energy < local_global_conv_rate || current_PD_energy < 5e-15;
	//	if (energy_changing) {
	//		return true;
	//	}
	//}
	//return false;
}

bool ProjectDynamic::IPC_PDConvergeCondition()
{
	if (outer_iteration_num > 2) {
		if (outer_iteration_num < max_it) {
			//bool system_energy = fabs(current_PD_energy - previous_PD_energy) / previous_PD_energy < outer_itr_conv_rate || current_PD_energy < 5e-15;
			//bool collision_energy = fabs(previous_collision_energy - current_collision_energy) / previous_collision_energy < local_global_conv_rate;
			//bool energy_changing = fabs(PD_energy_dif + (previous_PD_energy - current_PD_energy))/ current_PD_energy < outer_itr_conv_rate;
			//bool collision_energy_changing= fabs(collision_energy_dif + (previous_collision_energy - current_collision_energy)) / current_collision_energy < 0.1;
			//this is actually the changing ratio between itr-2 and itr;
			//std::cout <<current_PD_energy << " " << current_collision_energy << std::endl;
			//std::cout << fabs(current_PD_energy - previous_PD_energy) / previous_PD_energy << std::endl;
			//if (outer_iteration_num > 990) {
				//std::cout <<current_collision_energy <<" " << fabs(collision_energy_dif + (previous_collision_energy - current_collision_energy)) << std::endl;
				//std::cout << current_PD_energy <<" " << fabs(PD_energy_dif + (previous_PD_energy - current_PD_energy)) << std::endl;
				//std::cout << fabs(current_PD_energy - previous_PD_energy) / previous_PD_energy <<" " << fabs(PD_energy_dif + (previous_PD_energy - current_PD_energy)) / current_PD_energy << std::endl;
			//}
			//if (system_energy || energy_changing){//|| current_collision_energy<1e-6) {//&&(collision_energy|| collision_energy_changing)
				//std::cout << outer_iteration_num << std::endl;
			thread->assignTask(this, COMPUTE_DISPLACEMENT);
			displacement_norm = displacement_norm_thread[0];
			for (int i = 1; i < total_thread_num; ++i) {
				//if (displacement_norm < displacement_norm_thread[i]) {
					displacement_norm += displacement_norm_thread[i];
				//}
			}
	

			bool ratio_changing = fabs(displacement_ratio_dif + (previous_displacement_norm - displacement_norm)) / displacement_norm < 1e-4;
			//if (outer_iteration_num > 990) {
				//std::cout << "displacement ratio " << displacement_norm / displacement_bound <<" "<< fabs(displacement_ratio_dif + (previous_displacement_norm - displacement_norm)) / displacement_norm << std::endl;
			//}
			if (displacement_norm / displacement_bound < 1.0 || ratio_changing) {//  
				return true;
			}
			//}
		}
		else {
			////std::cout << "larger than 1000 " << std::endl;
			return true;
		}
	}
	return false;
}

bool ProjectDynamic::PDConvergeCondition()
{
	bool system_energy = fabs(current_PD_energy - previous_itr_PD_energy) / previous_itr_PD_energy < outer_itr_conv_rate || current_PD_energy < 5e-15;
	bool collision_energy = fabs(previous_itr_collision_energy - current_collision_energy) / previous_itr_collision_energy < local_global_conv_rate || current_collision_energy < 5e-15;

	bool energy_satisfied = system_energy && collision_energy;
	bool need_to_stop = local_global_iteration_num > max_it - 3 || fabs(current_PD_energy - previous_itr_PD_energy) / previous_itr_PD_energy < 5e-6;

	bool standard = (energy_satisfied || need_to_stop) && outer_iteration_num > 1;
	if (standard) {
		return true;
	}
	else {
		////std::cout << fabs(current_PD_energy - previous_itr_PD_energy) / previous_itr_PD_energy << " " << current_PD_energy << std::endl;
		return false;
	}

}

bool ProjectDynamic::PDLocalGlobalConvergeCondition()
{
	bool system_energy = fabs(current_PD_energy - previous_PD_energy) / previous_PD_energy < local_global_conv_rate;
	bool constraint_energy = fabs(current_constraint_energy - previous_constraint_energy) / previous_constraint_energy < local_global_conv_rate || current_constraint_energy < 1e-14;
	bool collision_energy = fabs(previous_collision_energy - current_collision_energy) / previous_collision_energy < local_global_conv_rate || current_collision_energy < 1e-14;
	bool need_to_stop = fabs(current_PD_energy - previous_PD_energy) / previous_PD_energy < 5e-6 || local_global_iteration_num > max_it - 3;
	bool energy_satisfied = system_energy && constraint_energy && collision_energy;
	bool standard = (energy_satisfied || need_to_stop) && local_global_itr_in_single_outer > 0;

	//if (local_global_itr_in_single_outer > 300) {
	//	//std::cout<<"energy " << current_collision_energy<<" "<< abs(previous_collision_energy - current_collision_energy) / previous_collision_energy << std::endl;
	//	//std::cout << current_constraint_energy<<" "<< abs(current_constraint_energy - previous_constraint_energy) / previous_constraint_energy << std::endl;
	//	//std::cout << current_PD_energy <<" "<< abs(current_PD_energy - previous_PD_energy) / previous_PD_energy << std::endl;
	//}

	if (standard) {
		return true;
	}
	else {
		return false;
	}
}




void ProjectDynamic::localProjection()
{
	//time_t t = clock();
	//for (int i = 0; i < 1000; ++i) {
	//	thread->assignTask(this, TEST_LOCAL_PROJECTION);
	//}
	////std::cout <<"time "<< clock() - t << std::endl;;

	thread->assignTask(this, LOCAL_PROJECTION);
}


//TEST_LOCAL_PROJECTION
void ProjectDynamic::testLocalProjectionPerThread(int thread_id)
{
	temEnergy[thread_id] = 0.0;
	localBendingProjectionPerThread(thread_id, true);
}


//LOCAL_PROJECTION
//LOCAL_PROJECTION_WITHOUT_ENERGY
void ProjectDynamic::localProjectionPerThread(int thread_id, bool with_energy)
{
	//edge length
	localEdgeLengthProjectionPerThread(thread_id, with_energy);
	//bending
	localBendingProjectionPerThread(thread_id, with_energy);
	//anchor
	localPositionProjectionPerThread(thread_id, with_energy);
}



void ProjectDynamic::solveClothSystemPerThead(int thread_id, bool with_collision, bool compute_energy)
{
	int cloth_No;
	int dimension;
	temEnergy[thread_id] = 0.0;
	for (int i = 0; i < cloth_dimension_per_thread[thread_id].size(); ++i) {
		cloth_No = cloth_dimension_per_thread[thread_id][i] / 3;
		dimension = cloth_dimension_per_thread[thread_id][i] % 3;
		solveClothSystemPerThead(cloth_b[cloth_No][dimension], cloth_u[cloth_No][dimension], (*cloth)[cloth_No].mesh_struct, (*cloth)[cloth_No].length_stiffness,
			p_edge_length[cloth_No], cloth_No, dimension, vertex_around_vertex_for_bending[cloth_No], vertex_lbo[cloth_No], p_bending[cloth_No],
			cloth_u_prediction[cloth_No][dimension], thread_id, collision.cloth_target_pos.b_sum[cloth_No], collision.cloth_target_pos.need_update[cloth_No],
			with_collision, compute_energy);

	}
}


//LOCAL_EDGE_LENGTH_PROJECTION
//void ProjectDynamic::localEdgeLengthProjectionPerThread(int thread_id)
//{
//	Vector3d q0, q1, q01;
//	double curLen;
//	TriangleMeshStruct* mesh_struct;
//	double ori_edge_length;
//	int edge_index;
//	VectorXd* u;
//	int* neighbor_vertex;
//	for (int j = 0; j < total_cloth_num; ++j) {
//		mesh_struct = &(*cloth)[j].mesh_struct;
//		u = cloth_u[j].data();
//		for (int i = mesh_struct->vertex_index_begin_per_thread[thread_id]; i < mesh_struct->vertex_index_begin_per_thread[thread_id + 1]; ++i) {
//			q0.data()[0] = u[0].data()[i];
//			q0.data()[1] = u[1].data()[i];
//			q0.data()[2] = u[2].data()[i];
//			memset(p_edge_length[j][i].data(), 0, 24);
//			neighbor_vertex = mesh_struct->vertices[i].neighbor_vertex.data();
//			for (int k = 0; k < mesh_struct->vertices[i].neighbor_vertex.size(); ++k) {
//				q1.data()[0] = u[0].data()[neighbor_vertex[k]];
//				q1.data()[1] = u[1].data()[neighbor_vertex[k]];
//				q1.data()[2] = u[2].data()[neighbor_vertex[k]];
//				edge_index = mesh_struct->vertices[i].edge[k];
//				q01 = q0 - q1;		
//				curLen = sqrt(dotProduct(q01.data(), q01.data()));				
//				ori_edge_length = mesh_struct->edges[edge_index].length;
//				temEnergy[thread_id] += 0.25 * (*cloth)[j].length_stiffness[edge_index] * (curLen - ori_edge_length) * (curLen - ori_edge_length);
//				p_edge_length[j][i] +=  q01* ((*cloth)[j].length_stiffness[edge_index] / curLen * ori_edge_length);
//			}
//		}
//	}
//}
//LOCAL_EDGE_LENGTH_PROJECTION
void ProjectDynamic::localEdgeLengthProjectionPerThread(int thread_id, bool with_energy)
{
	if (with_energy) {
		temEnergy[thread_id] = 0.0;
		Vector3d q0, q1, q01;
		double curLen;
		VectorXd* u;
		MeshStruct::Edge* edges;
		int* edge_index_begin_per_thread;
		double* length_stiffness;
		for (int j = 0; j < total_cloth_num; ++j) {
			u = cloth_u[j].data();
			edges = (*cloth)[j].mesh_struct.edges.data();
			edge_index_begin_per_thread = (*cloth)[j].mesh_struct.edge_index_begin_per_thread.data();
			length_stiffness = (*cloth)[j].length_stiffness.data();
			for (int i = edge_index_begin_per_thread[thread_id]; i < edge_index_begin_per_thread[thread_id + 1]; ++i) {
				q0.data()[0] = u[0].data()[edges[i].vertex[0]];
				q1.data()[0] = u[0].data()[edges[i].vertex[1]];
				q0.data()[1] = u[1].data()[edges[i].vertex[0]];
				q1.data()[1] = u[1].data()[edges[i].vertex[1]];
				q0.data()[2] = u[2].data()[edges[i].vertex[0]];
				q1.data()[2] = u[2].data()[edges[i].vertex[1]];
				q01 = q0 - q1;
				curLen = sqrt(dotProduct(q01.data(), q01.data()));
				temEnergy[thread_id] += 0.5 * length_stiffness[i] * (curLen - edges[i].length) * (curLen - edges[i].length);
				p_edge_length[j][i] = q01 * (length_stiffness[i] / curLen * edges[i].length);
			}
		}
	}
	else {
		Vector3d q0, q1, q01;
		double curLen;
		VectorXd* u;
		MeshStruct::Edge* edges;
		int* edge_index_begin_per_thread;
		double* length_stiffness;
		for (int j = 0; j < total_cloth_num; ++j) {
			u = cloth_u[j].data();
			edges = (*cloth)[j].mesh_struct.edges.data();
			edge_index_begin_per_thread = (*cloth)[j].mesh_struct.edge_index_begin_per_thread.data();
			length_stiffness = (*cloth)[j].length_stiffness.data();
			for (int i = edge_index_begin_per_thread[thread_id]; i < edge_index_begin_per_thread[thread_id + 1]; ++i) {
				q0.data()[0] = u[0].data()[edges[i].vertex[0]];
				q1.data()[0] = u[0].data()[edges[i].vertex[1]];
				q0.data()[1] = u[1].data()[edges[i].vertex[0]];
				q1.data()[1] = u[1].data()[edges[i].vertex[1]];
				q0.data()[2] = u[2].data()[edges[i].vertex[0]];
				q1.data()[2] = u[2].data()[edges[i].vertex[1]];
				q01 = q0 - q1;
				curLen = sqrt(dotProduct(q01.data(), q01.data()));
				p_edge_length[j][i] = q01 * (length_stiffness[i] / curLen * edges[i].length);
			}
		}
	}
}


void ProjectDynamic::localBendingProjectionPerThread(int thread_id, bool with_energy)
{
	int size; std::vector<VectorXd>q(3);
	Vector3d aq; double aqnorm; //double a_rest[3];
	VectorXd* u;
	VectorXd* lbo;
	std::vector<std::vector<int>>* vertex_around_vertex_for_bending_;
	double bend_stiffness;
	int* vertex_index_begin_per_thread;
	Vector3d p;
	VectorXd* p_bend;
	if (with_energy) {
		for (int k = 0; k < total_cloth_num; ++k) {
			if (!(*cloth)[k].mesh_struct.faces.empty()) {
				bend_stiffness = (*cloth)[k].bend_stiffness;
				u = cloth_u[k].data();
				vertex_around_vertex_for_bending_ = &vertex_around_vertex_for_bending[k];
				lbo = vertex_lbo[k].data();
				vertex_index_begin_per_thread = (*cloth)[k].mesh_struct.vertex_index_begin_per_thread.data();
				for (int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id + 1]; ++i) {
					p_bend = p_bending[k][i].data();
					size = (*vertex_around_vertex_for_bending_)[i].size();
					for (int j = 0; j < 3; ++j) {
						q[j].resize(size);
					}
					for (int h = 0; h < size; h++) {
						q[0][h] = u[0].data()[(*vertex_around_vertex_for_bending_)[i][h]];
						q[1][h] = u[1].data()[(*vertex_around_vertex_for_bending_)[i][h]];
						q[2][h] = u[2].data()[(*vertex_around_vertex_for_bending_)[i][h]];
					}
					aq.data()[0] = lbo[i].dot(q[0]);
					aq.data()[1] = lbo[i].dot(q[1]);
					aq.data()[2] = lbo[i].dot(q[2]);
					aqnorm = sqrt(dotProduct(aq.data(), aq.data()));
					if (aqnorm < 1e-10)
						p = aq;
					else {
						p = aq * rest_mean_curvature_norm[k][i] / aqnorm;
					}
					temEnergy[thread_id] += 0.5 * bend_stiffness * norm2(aq - p);
					for (int j = 0; j < 3; ++j) {
						p_bend[j] = lbo[i] * (p.data()[j] * bend_stiffness * lbo_weight[k][i]);
					}
				}
			}
		}
	}
	else {
		for (int k = 0; k < total_cloth_num; ++k) {
			if (!(*cloth)[k].mesh_struct.faces.empty()) {
				bend_stiffness = (*cloth)[k].bend_stiffness;
				u = cloth_u[k].data();
				vertex_around_vertex_for_bending_ = &vertex_around_vertex_for_bending[k];
				lbo = vertex_lbo[k].data();
				vertex_index_begin_per_thread = (*cloth)[k].mesh_struct.vertex_index_begin_per_thread.data();
				for (int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id + 1]; ++i) {
					p_bend = p_bending[k][i].data();
					size = (*vertex_around_vertex_for_bending_)[i].size();
					for (int j = 0; j < 3; ++j) {
						q[j].resize(size);
					}
					for (int h = 0; h < size; h++) {
						q[0][h] = u[0].data()[(*vertex_around_vertex_for_bending_)[i][h]];
						q[1][h] = u[1].data()[(*vertex_around_vertex_for_bending_)[i][h]];
						q[2][h] = u[2].data()[(*vertex_around_vertex_for_bending_)[i][h]];
					}
					aq.data()[0] = lbo[i].dot(q[0]);
					aq.data()[1] = lbo[i].dot(q[1]);
					aq.data()[2] = lbo[i].dot(q[2]);
					aqnorm = sqrt(dotProduct(aq.data(), aq.data()));
					if (aqnorm < 1e-10)
						p = aq;
					else {
						p = aq * rest_mean_curvature_norm[k][i] / aqnorm;
					}
					for (int j = 0; j < 3; ++j) {
						p_bend[j] = lbo[i] * (p.data()[j] * bend_stiffness * lbo_weight[k][i]);
					}
				}
			}
		}
	}
}
void ProjectDynamic::localPositionProjectionPerThread(int thread_id, bool with_energy)
{
	double delta_q[3];
	TriangleMeshStruct* mesh_struct;
	if (with_energy) {
		for (int j = 0; j < total_cloth_num; ++j) {
			if (!(*cloth)[j].mesh_struct.anchor_vertex.empty()) {
				mesh_struct = &(*cloth)[j].mesh_struct;
				for (int i = mesh_struct->anchor_index_begin_per_thread[thread_id]; i < mesh_struct->anchor_index_begin_per_thread[thread_id + 1]; ++i)
				{
					delta_q[0] = cloth_u[j][0].data()[mesh_struct->anchor_vertex[i]] - mesh_struct->anchor_position[i][0];
					delta_q[1] = cloth_u[j][1].data()[mesh_struct->anchor_vertex[i]] - mesh_struct->anchor_position[i][1];
					delta_q[2] = cloth_u[j][2].data()[mesh_struct->anchor_vertex[i]] - mesh_struct->anchor_position[i][2];
					temEnergy[thread_id] += 0.5 * (*cloth)[j].position_stiffness * dotProduct(delta_q, delta_q);
				}
			}
		}
	}
	else {
		for (int j = 0; j < total_cloth_num; ++j) {
			if (!(*cloth)[j].mesh_struct.anchor_vertex.empty()) {
				mesh_struct = &(*cloth)[j].mesh_struct;
				for (int i = mesh_struct->anchor_index_begin_per_thread[thread_id]; i < mesh_struct->anchor_index_begin_per_thread[thread_id + 1]; ++i)
				{
					delta_q[0] = cloth_u[j][0].data()[mesh_struct->anchor_vertex[i]] - mesh_struct->anchor_position[i][0];
					delta_q[1] = cloth_u[j][1].data()[mesh_struct->anchor_vertex[i]] - mesh_struct->anchor_position[i][1];
					delta_q[2] = cloth_u[j][2].data()[mesh_struct->anchor_vertex[i]] - mesh_struct->anchor_position[i][2];
				}
			}
		}
	}
}



//SOLVE_SYSYTEM
//SOLVE_SYSYTEM_WITHOUT_COLLISION
//SOLVE_SYSYTEM_WITHOUT_ENERGY
void ProjectDynamic::solveSystemPerThead(int thread_id, bool with_collision, bool compute_energy)
{
	temEnergy[thread_id] = 0;
	solveClothSystemPerThead(thread_id, with_collision, compute_energy);
}


void ProjectDynamic::solveClothSystemPerThead(VectorXd& b, VectorXd& u, TriangleMeshStruct& mesh_struct, std::vector<double>& length_stiffness,
	std::vector<Vector3d>& p_edge_length, int cloth_No, int dimension, std::vector<std::vector<int>>& vertex_around_vertex_for_bending,
	std::vector<VectorXd>& vertex_lbo, std::vector<std::vector<VectorXd>>& p_bending, VectorXd& u_prediction, int thread_id,
	std::vector<std::array<double, 3>>& collision_b_sum, bool* collision_b_need_update, bool with_collision, bool compute_energy)
{
	b.setZero();
	//collision
	if (with_collision) {
		for (int i = 0; i < cloth_sys_size[cloth_No]; ++i) {
			if (collision_b_need_update[i]) {
				b.data()[i] += collision_b_sum[i][dimension];
			}
		}
	}
	//edge length
	for (int i = 0; i < mesh_struct.edges.size(); ++i) {
		b.data()[mesh_struct.edges[i].vertex[0]] += p_edge_length[i].data()[dimension];
		b.data()[mesh_struct.edges[i].vertex[1]] -= p_edge_length[i].data()[dimension];
	}
	//for (int i = 0; i < cloth_sys_size[cloth_No]; ++i) {
	//	b.data()[i] += p_edge_length[i].data()[dimension];
	//}
	//bending	
	if (!mesh_struct.faces.empty()) {
		for (int i = 0; i < cloth_sys_size[cloth_No]; ++i) {
			for (int j = 0; j < vertex_lbo[i].size(); ++j) {
				b.data()[vertex_around_vertex_for_bending[i][j]] += p_bending[i][dimension].data()[j];
			}
		}
	}
	//position
	for (int i = 0; i < (*cloth)[cloth_No].mesh_struct.anchor_vertex.size(); ++i) {
		b.data()[mesh_struct.anchor_vertex[i]] += (*cloth)[cloth_No].position_stiffness * mesh_struct.anchor_position[i][dimension];
	}
	b += (1.0 / (sub_time_step * sub_time_step)) * (cloth_mass[cloth_No].cwiseProduct(u_prediction));
	if (with_collision) {
		int itr_num;
		switch (itr_solver_method)
		{
		case DIRECT_SOLVE:
			u = cloth_llt[cloth_No].solve(b);
			//std::cout << "direct solve " << std::endl;
			break;
		case JACOBI: {			
			iteration_method.solveByJacobi(u, b, cloth_global_mat[cloth_No], cloth_No, itr_num);
			//std::cout << "Jacobi "<< itr_num << std::endl;
		}
			break;
		case SUPER_JACOBI: {
			iteration_method.solveBySuperJacobi(u, b, cloth_global_mat[cloth_No], cloth_No, itr_num);
			//std::cout << "super jacobi " << itr_num << std::endl;
		}
			break;
		case CHEBYSHEV_SUPER_JACOBI: {
			iteration_method.solveByChebyshevSemiIterativeSuperJacobi(u, b, cloth_global_mat[cloth_No], cloth_No, itr_num);
			//std::cout << "chebyshev super jacobi "<< itr_num << std::endl;
		}
			break;
		case GAUSS_SEIDEL: {
			iteration_method.solveByGaussSeidel(u, b, cloth_global_mat[cloth_No], cloth_No, itr_num);
			//std::cout << "gauss_seidel "<< itr_num << std::endl;
		}
			break;
		case GAUSS_SEIDEL_CHEBYSHEV: {
			iteration_method.solveByChebyshevGaussSeidel(u, b, cloth_global_mat[cloth_No], cloth_No, itr_num,0.6);		
			//std::cout << "gauss_seidel_chebysev "<< itr_num << std::endl;
		}
			break;
		case CHEBYSHEV_JACOBI: {
			iteration_method.solveByChebyshevSemiIterativeJacobi(u, b, cloth_global_mat[cloth_No], cloth_No, itr_num);
			//std::cout <<"chebyshev jacobi "<< itr_num << std::endl;
		}
			break;
		case PCG: {
			iteration_method.solveByPCG(u, b, cloth_global_mat[cloth_No], cloth_No, itr_num);
			//std::cout <<"PCG "<< itr_num << std::endl;
		}
			break;
		case WEIGHTED_JACOBI: {
			iteration_method.solveByWeightedJacobi(u, b, cloth_global_mat[cloth_No], cloth_No, itr_num,2.0/3.0);
		}
			break;
		}
		if (compute_energy) {
			VectorXd p0, p1;
			p0 = u - u_prediction;
			p1 = p0.cwiseProduct(cloth_mass[cloth_No]);
			temEnergy[thread_id] += 0.5 / (sub_time_step * sub_time_step) * p0.dot(p1);
		}
		ave_iteration[cloth_No][dimension] += itr_num;
	}
	else {
		u = collision_free_cloth_llt[cloth_No].solve(b);
	}
}




void ProjectDynamic::addExternalClothForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int cloth_No)
{
	if (!coe.empty()) {
		for (int i = 0; i < coe.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				cloth_f_ext[cloth_No][j].data()[neighbor_vertex[i]] += coe[i] * neighbor_vertex_force_direction[j];
			}
		}
	}
}

void ProjectDynamic::resetExternalForce()
{
	cloth_f_ext = cloth_gravity;
	tetrahedron_f_ext = tetrahedron_gravity;
}

void ProjectDynamic::mainProcess()
{
	PDsetPosPredict();
	for (int i = 0; i < total_cloth_num; ++i) {
		(*cloth)[i].mesh_struct.getNormal();
	}


}


void ProjectDynamic::initialJacobi()
{
	computeOffDiagonal();
	iteration_method.initialJacobi();
}


void ProjectDynamic::updateIterateSolverParameter(double conv_rate)
{
	iteration_method.updateConvergenceRate(conv_rate);
}

void ProjectDynamic::saveMatrix(int dimension, VectorXd& vec, std::string name)
{
	std::string file_name_matrix = name+std::to_string(dimension)+".dat";
	EigenMatrixIO::write_binary(file_name_matrix.c_str(), vec);
}

void ProjectDynamic::saveSparseMatrix(SparseMatrix<double,RowMajor>& matrix, std::string file_name_matrix)
{
	EigenMatrixIO::write_sp_binary(file_name_matrix.c_str(), matrix);
}

void ProjectDynamic::reset_ave_iteration_record()
{
	for (int i = 0; i < ave_iteration.size(); ++i) {
		memset(ave_iteration[i].data(), 0, 24);
	}
}

void ProjectDynamic::update_ave_iteration_record(std::vector<std::vector<double>>& ave_itr)
{
	for (int i = 0; i < ave_iteration.size(); ++i) {
		for (int j = 0; j < 3; ++j) {
			ave_iteration[i][j] /= (max_inner_iteration_num * (outer_iteration_num - 1));
		}
	
	}
	ave_itr = ave_iteration;
}