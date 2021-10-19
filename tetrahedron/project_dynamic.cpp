#include"project_dynamic.h"

ProjectDynamic::ProjectDynamic()
{
	gravity_ = 9.8;
	total_thread_num = std::thread::hardware_concurrency();
	temEnergy.resize(total_thread_num);
	outer_itr_conv_rate = 1e-2;// 7.5e-2; 
	local_global_conv_rate = 2e-2;
	sub_step_num = 1;

	use_dierct_solve_for_coarest_mesh = true;	
	super_jacobi_step_size = 3;
	max_it = 1000;
	max_jacobi_itr_num = 20;
}

void ProjectDynamic::setForPD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Thread* thread)
{
	sub_time_step = time_step / (double)sub_step_num;
	this->thread = thread;
	setForClothPD(cloth);
	setForTetrahedronPD(tetrahedron);	
	setIndexPerThread();
	collision.initial(cloth, collider, tetrahedron, thread);
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

void ProjectDynamic::initialClothPDvariable()
{
	p_edge_length.resize(total_cloth_num);
	p_bending.resize(total_cloth_num);
	cloth_u.resize(total_cloth_num);
	cloth_u_.resize(total_cloth_num);
	cloth_v.resize(total_cloth_num);
	cloth_v_.resize(total_cloth_num);
	cloth_b.resize(total_cloth_num);
	cloth_u_prediction.resize(total_cloth_num);
	cloth_acceleration.resize(total_cloth_num);
	for (int j = 0; j < total_cloth_num; ++j) {
		p_edge_length[j].resize((*cloth)[j].mesh_struct.edges.size());
		p_bending[j].resize(cloth_sys_size[j]);
		cloth_u[j].resize(3);
		cloth_u_[j].resize(3);
		cloth_v[j].resize(3);
		cloth_v_[j].resize(3);
		cloth_b[j].resize(3);
		cloth_u_prediction[j].resize(3);
		cloth_acceleration[j].resize(3);
		for (int i = 0; i < 3; ++i) {
			cloth_u[j][i].resize(cloth_sys_size[j]);	cloth_u[j][i].setZero();
			cloth_u_[j][i].resize(cloth_sys_size[j]);	cloth_u_[j][i].setZero();
			cloth_v[j][i].resize(cloth_sys_size[j]); cloth_v[j][i].setZero();
			cloth_v_[j][i].resize(cloth_sys_size[j]); cloth_v_[j][i].setZero();
			cloth_b[j][i].resize(cloth_sys_size[j]); cloth_b[j][i].setZero();
			cloth_u_prediction[j][i].resize(cloth_sys_size[j]); cloth_u_prediction[j][i].setZero();
			cloth_acceleration[j][i].resize(cloth_sys_size[j]); cloth_acceleration[j][i].setZero();
			for (int k = 0; k < cloth_sys_size[j]; ++k) {
				cloth_u[j][i].data()[k] = (*cloth)[j].mesh_struct.vertex_position[k][i];			
			}
			cloth_u_[j][i] = cloth_u[j][i];
		}
	}
}
void ProjectDynamic::setIndexPerThread()
{
	cloth_per_thread_begin.resize(total_thread_num+1);
	if (total_cloth_num > 0) {
		(*cloth)[0].mesh_struct.arrangeIndex(total_thread_num, total_cloth_num, cloth_per_thread_begin);
	}
	cloth_dimension_per_thread.resize(total_thread_num);
	int total_num = total_cloth_num * 3;
	for (int i = 0; i < total_num; ++i) {
		cloth_dimension_per_thread[i % total_thread_num].push_back(i);
	}
	tetrahedron_per_thread_begin.resize(total_thread_num + 1);
	if (total_tetrahedron_num > 0) {
		(*tetrahedron)[0].mesh_struct.arrangeIndex(total_thread_num, total_tetrahedron_num, tetrahedron_per_thread_begin);
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
	for (int i = 0; i < total_cloth_num; ++i) {
		computeGlobalStepMatrixSingleCloth(&cloth_global_mat[i], cloth_global_mat_diagonal_ref[i], cloth_global_mat_diagonal_ref_address[i],
			&cloth_llt[i], (*cloth)[i].mesh_struct, vertex_around_vertex_for_bending[i], lbo_weight[i], vertex_lbo[i], cloth_sys_size[i],
			(*cloth)[i].bend_stiffness, (*cloth)[i].length_stiffness, (*cloth)[i].position_stiffness);
	}
	initial_cloth_global_mat = cloth_global_mat;
}
void ProjectDynamic::computeGlobalStepMatrixSingleCloth(SparseMatrix<double>* global_mat, std::vector<double>& global_mat_diagonal_ref,
	std::vector<double*>& global_collision_mat_diagonal_ref, SimplicialLLT<SparseMatrix<double>>* global_llt, TriangleMeshStruct& mesh_struct,
	std::vector<std::vector<int>>& vertex_around_vertex_for_bending,
	std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo, int sys_size, double bending_stiffness, std::vector<double>& length_stiffness,
	double position_stiffness)
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
	//for (int i = 0; i < mesh_struct.anchor_vertex.size(); ++i) {
	//	global_mat_nnz.push_back(Triplet<double>(mesh_struct.anchor_vertex[i], mesh_struct.anchor_vertex[i], position_stiffness));
	//}
	//mass
	for (int i = 0; i < sys_size; ++i) {
		global_mat_nnz.push_back(Triplet<double>(i, i, mesh_struct.vertices[i].mass / (sub_time_step * sub_time_step)));
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

void ProjectDynamic::setAroundVertexPrimitive(TriangleMeshStruct& mesh_struct,std::vector<std::vector<int>>& edge_around_vertex_for_bending,
	std::vector<std::vector<int>>& vertex_around_vertex_for_bending)
{
	std::vector<TriangleMeshStruct::Vertex>* vertex;
	std::vector<int>* around_edge_index;
	for (int i = 0; i < mesh_struct.vertices.size(); ++i) {
		around_edge_index = &mesh_struct.vertices[i].edge;
		edge_around_vertex_for_bending[i].reserve(around_edge_index->size());
		vertex_around_vertex_for_bending[i].reserve(around_edge_index->size()+1);
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
	std::vector<TriangleMeshStruct::Vertex>* vertex;
	for (int j = 0; j < total_cloth_num; ++j) {
		cloth_f_ext[j].resize(3);
		cloth_gravity[j].resize(3);
		vertex = &(*cloth)[j].mesh_struct.vertices;
		for (int i = 0; i < 3; ++i){
			if (cloth_sys_size[j] > 1) {
				cloth_gravity[j][i].resize(cloth_sys_size[j]);
				for (int k = 0; k < cloth_sys_size[j]; ++k) {
					cloth_gravity[j][i][k] = gravity_accerlation[i] * (*vertex)[k].mass;
				}
			}
			else {
				cloth_gravity[j][i] = VectorXd::Zero(cloth_sys_size[j]);
			}
			cloth_f_ext[j][i] = cloth_gravity[j][i];
		}
	}
}

void ProjectDynamic::setForTetrahedronPD(std::vector<Tetrahedron>* tetrahedron)
{
	total_tetrahedron_num = tetrahedron->size();
	tetrahedron_sys_size.resize(total_tetrahedron_num);
	for (int i = 0; i < total_tetrahedron_num; ++i) {
		tetrahedron_sys_size[i] = (*tetrahedron)[i].mesh_struct.vertex_position.size();
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
			mesh_struct.vertices[i].mass = mass_;
		}
	}
	else {
		for (int i = 0; i < mass_inv.size(); ++i) {
			mass_ = 1.25;
			mass_inv.data()[i] = 1.0 / mass_;
			mass.data()[i] = mass_;
			mesh_struct.vertices[i].mass = mass_;
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
	for (int j = 0; j < total_cloth_num; ++j) {
		for (int i = 0; i < cloth_sys_size[j]; ++i) {
			(*cloth)[j].mesh_struct.vertex_position[i][0] = cloth_u[j][0].data()[i];
			(*cloth)[j].mesh_struct.vertex_position[i][1] = cloth_u[j][1].data()[i];
			(*cloth)[j].mesh_struct.vertex_position[i][2] = cloth_u[j][2].data()[i];
		}
	}
	for (int j = 0; j < total_tetrahedron_num; ++j) {
		for (int i = 0; i < tetrahedron_sys_size[j]; ++i) {
			(*tetrahedron)[j].mesh_struct.vertex_position[i][0] = tetrahedron_u[j][0].data()[i];
			(*tetrahedron)[j].mesh_struct.vertex_position[i][1] = tetrahedron_u[j][1].data()[i];
			(*tetrahedron)[j].mesh_struct.vertex_position[i][2] = tetrahedron_u[j][2].data()[i];
		}
	}
}

void ProjectDynamic::updateRenderPosition()
{
	for (int j = 0; j < total_cloth_num; ++j) {
		(*cloth)[j].mesh_struct.vertex_for_render = (*cloth)[j].mesh_struct.vertex_position;
		(*cloth)[j].mesh_struct.getRenderNormal();

	}
	for (int j = 0; j < total_tetrahedron_num; ++j) {
		(*tetrahedron)[j].mesh_struct.vertex_for_render = (*tetrahedron)[j].mesh_struct.vertex_position;
		(*tetrahedron)[j].mesh_struct.getRenderNormal();
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
			cloth_u_prediction[j][i] = cloth_u_[j][i] + sub_time_step * cloth_v_[j][i] + sub_time_step * sub_time_step * cloth_mass_inv[j].cwiseProduct(cloth_f_ext[j][i]);
		}
	}
}


void ProjectDynamic::PDsolve()
{
	PDsetPosPredict();
	int itr_num = 0;
	outer_iteration_num = 0;
	initialEnergy();
	current_PD_energy = 1e-15;
	current_collision_energy = 1e-15;
	current_constraint_energy = 1e-15;
	PDupdateSystemMatrix();
	local_global_iteration_num = 0;
	while (!PDConvergeCondition())	{
		initialEnergyOuterInteration();
		local_global_itr_in_single_outer = 0;
		//collision.findAllTrianglePairs();
		while (!PDLocalGlobalConvergeCondition()){
			initialEnergyLocalGlobal();
			localProjection();
			current_constraint_energy += current_collision_energy;
			thread->assignTask(this, SOLVE_SYSYTEM);//solve b	
			for (int i = 0; i < total_thread_num; ++i) {
				current_PD_energy += temEnergy[i];
			}
			current_PD_energy += current_constraint_energy;
			updateModelPosition();			
			local_global_itr_in_single_outer++;
		}
		outer_iteration_num++;
		local_global_iteration_num+= local_global_itr_in_single_outer;
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
	updateMatrix();
	thread->assignTask(this, MATRIX_DECOMPOSITION);
}

void ProjectDynamic::updateMatrix()
{
	thread->assignTask(this, UPDATE_MATRIX);
	updateCollisionMatrix();
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

//UPDATE_MATRIX
void ProjectDynamic::updateMatrixPerThread(int thread_No)
{
	setClothMatrix(thread_No);
}

void ProjectDynamic::setClothMatrix(int thread_No)
{
	int collision_num = 0;
	std::vector<int>* vertex_index_begin_per_thread;
	for (int j = 0; j < total_cloth_num; ++j) {
		vertex_index_begin_per_thread = &(*cloth)[j].mesh_struct.vertex_index_begin_per_thread;
		for (int i = (*vertex_index_begin_per_thread)[thread_No]; i < (*vertex_index_begin_per_thread)[thread_No+1]; ++i) {
			*(cloth_global_mat_diagonal_ref_address[j][i]) = cloth_global_mat_diagonal_ref[j][i];
			//if ((*vertex_collision_sum).need_update[j][i]) {
			//	*(global_collision_mat_diagonal_ref[j][i]) += (*vertex_collision_sum).stiffness[j][i];
			//	*(global_mat_diagonal_diagonal_ref[j][i]) += (*vertex_collision_sum).stiffness[j][i];
			//	diagonal_approx_global_matrix[j][i] += (*vertex_collision_sum).stiffness[j][i];
			//}
		}
	}


}

//MATRIX_DECOMPOSITION
void ProjectDynamic::matrixDecomposition(int thread_id)
{	

	for (int i = cloth_per_thread_begin[thread_id]; i < cloth_per_thread_begin[thread_id+1]; ++i) {
		cloth_llt[i].factorize(cloth_global_mat[i]);
	}
	for (int i = tetrahedron_per_thread_begin[thread_id]; i < tetrahedron_per_thread_begin[thread_id + 1]; ++i) {
		tetrahedron_llt[i].factorize(cloth_global_mat[i]);
	}
}

bool ProjectDynamic::PDConvergeCondition()
{
	bool system_energy = abs(current_PD_energy - previous_itr_PD_energy) / previous_itr_PD_energy < outer_itr_conv_rate || current_PD_energy < 1.1e-15;
	bool collision_energy = abs(previous_itr_collision_energy - current_collision_energy) / previous_itr_collision_energy < local_global_conv_rate || current_collision_energy < 1.1e-15;

	bool energy_satisfied = system_energy && collision_energy;
	bool need_to_stop = local_global_iteration_num > max_it - 3 || abs(current_PD_energy - previous_itr_PD_energy) / previous_itr_PD_energy < 5e-6;

	bool standard = (energy_satisfied || need_to_stop) && outer_iteration_num > 0;
	if (standard) {
		return true;
	}
	else {
		//std::cout << abs(current_PD_energy - previous_itr_PD_energy) / previous_itr_PD_energy << " " << current_PD_energy << std::endl;
		return false;
	}

}

bool ProjectDynamic::PDLocalGlobalConvergeCondition()
{
	bool system_energy = abs(current_PD_energy - previous_PD_energy) / previous_PD_energy < local_global_conv_rate || current_PD_energy < 1.1e-15;
	bool constraint_energy = abs(current_constraint_energy - previous_constraint_energy) / previous_constraint_energy < local_global_conv_rate || current_constraint_energy < 1.1e-15;
	bool collision_energy = abs(previous_collision_energy - current_collision_energy) / previous_collision_energy < local_global_conv_rate || current_collision_energy < 1.1e-15;
	bool need_to_stop = abs(current_PD_energy - previous_PD_energy) / previous_PD_energy < 5e-6 || local_global_iteration_num > max_it - 3;
	bool energy_satisfied = system_energy && constraint_energy && collision_energy;
	bool standard = (energy_satisfied || need_to_stop) && local_global_itr_in_single_outer > 0;
	
	if (standard) {
		return true;
	}
	else {
		return false;
	}
}


void ProjectDynamic::localProjection()
{
	thread->assignTask(this, LOCAL_PROJECTION);
	for (int i = 0; i < total_thread_num; ++i) {
		current_constraint_energy += temEnergy[i];
	}
}


//LOCAL_PROJECTION
void ProjectDynamic::localProjectionPerThread(int thread_id)
{
	//edgeLengthConstraint
	temEnergy[thread_id] = 0.0;
	localEdgeLengthProjectionPerThread(thread_id);
	//bending
	localBendingProjectionPerThread(thread_id);
	//anchor
	localPositionProjectionPerThread(thread_id);
}

void ProjectDynamic::localEdgeLengthProjectionPerThread(int thread_id)
{
	Vector3d q0, q1, q01;
	double curLen;
	TriangleMeshStruct* mesh_struct;
	for (int j = 0; j < total_cloth_num; ++j) {
		mesh_struct = &(*cloth)[j].mesh_struct;
		for (int i = mesh_struct->edge_index_begin_per_thread[thread_id]; i < mesh_struct->edge_index_begin_per_thread[thread_id+1]; ++i) {
			q0.data()[0] = cloth_u[j][0].data()[mesh_struct->edges[i].vertex[0]];
			q1.data()[0] = cloth_u[j][0].data()[mesh_struct->edges[i].vertex[1]];
			q0.data()[1] = cloth_u[j][1].data()[mesh_struct->edges[i].vertex[0]];
			q1.data()[1] = cloth_u[j][1].data()[mesh_struct->edges[i].vertex[1]];
			q0.data()[2] = cloth_u[j][2].data()[mesh_struct->edges[i].vertex[0]];
			q1.data()[2] = cloth_u[j][2].data()[mesh_struct->edges[i].vertex[1]];
			q01 = q0 - q1;
			curLen = sqrt(dotProduct(q01.data(), q01.data()));
			temEnergy[thread_id] += 0.5 * (*cloth)[j].length_stiffness[i] * (curLen - mesh_struct->edges[i].length) * (curLen - mesh_struct->edges[i].length);
			p_edge_length[j][i] = q01 / curLen * mesh_struct->edges[i].length;
		}
	}
}

void ProjectDynamic::localBendingProjectionPerThread(int thread_id)
{
	int size; std::vector<VectorXd>q(3);
	Vector3d aq; double aqnorm; //double a_rest[3];
	std::vector<std::vector<int>>* vertex_around_vertex_for_bending_;
	for (int k = 0; k < total_cloth_num; ++k) {
		if (!(*cloth)[k].mesh_struct.faces.empty()) {
			vertex_around_vertex_for_bending_ = &vertex_around_vertex_for_bending[k];
			for (int i = (*cloth)[k].mesh_struct.vertex_index_begin_per_thread[thread_id]; i < (*cloth)[k].mesh_struct.vertex_index_begin_per_thread[thread_id+1]; ++i) {
				size = (*vertex_around_vertex_for_bending_)[i].size();
				for (int j = 0; j < 3; ++j) {
					q[j].resize(size);				
				}
				for (int h = 0; h < size; h++) {
					q[0][h] = cloth_u[k][0].data()[(*vertex_around_vertex_for_bending_)[i][h]];
					q[1][h] = cloth_u[k][1].data()[(*vertex_around_vertex_for_bending_)[i][h]];
					q[2][h] = cloth_u[k][2].data()[(*vertex_around_vertex_for_bending_)[i][h]];
				}
				aq.data()[0] = vertex_lbo[k][i].dot(q[0]);
				aq.data()[1] = vertex_lbo[k][i].dot(q[1]);
				aq.data()[2] = vertex_lbo[k][i].dot(q[2]);
				aqnorm = sqrt(dotProduct(aq.data(), aq.data()));
				if (aqnorm < 1e-10)
					p_bending[k][i] = aq;
				else {
					//memcpy(a_rest, &(rest_mean_curvature[k].data()[3 * i]), 24);
					p_bending[k][i] = aq * rest_mean_curvature_norm[k][i] / aqnorm;
				}
				temEnergy[thread_id] += 0.5 * (*cloth)[k].bend_stiffness * norm2(aq - p_bending[k][i]);
			}
		}
	}
}
void ProjectDynamic::localPositionProjectionPerThread(int thread_id)
{
	double delta_q[3];
	TriangleMeshStruct* mesh_struct;
	for (int j = 0; j < total_cloth_num; ++j) {
		if (!(*cloth)[j].mesh_struct.anchor_vertex.empty()) {
			mesh_struct = &(*cloth)[j].mesh_struct;
			for (int i = mesh_struct->anchor_index_begin_per_thread[thread_id]; i < mesh_struct->anchor_index_begin_per_thread[thread_id+1]; ++i)
			{
				delta_q[0] = cloth_u[j][0].data()[mesh_struct->anchor_vertex[i]] - mesh_struct->anchor_position[i][0];
				delta_q[1] = cloth_u[j][1].data()[mesh_struct->anchor_vertex[i]] - mesh_struct->anchor_position[i][1];
				delta_q[2] = cloth_u[j][2].data()[mesh_struct->anchor_vertex[i]] - mesh_struct->anchor_position[i][2];
				temEnergy[thread_id] += 0.5 * (*cloth)[j].position_stiffness * dotProduct(delta_q, delta_q);
			}
		}
	}
}



//SOLVE_SYSYTEM
void ProjectDynamic::solveSystemPerThead(int thread_id)
{
	solveClothSystemPerThead(thread_id);
}


void ProjectDynamic::solveClothSystemPerThead(VectorXd& b, VectorXd& u, TriangleMeshStruct& mesh_struct, std::vector<double>& length_stiffness,
	std::vector<Vector3d>& p_edge_length,int cloth_No, int dimension, std::vector<std::vector<int>>& vertex_around_vertex_for_bending,
	std::vector<VectorXd>& vertex_lbo, std::vector<Vector3d>& p_bending,std::vector<double>& lbo_weight, VectorXd& u_prediction, int thread_id)
{
	b.setZero();
	//collision
	//for (int i = 0; i < sys_size[cloth_No]; ++i) {
	//	if ((*vertex_collision_sum).need_update[cloth_No][i]) {
	//		b[cloth_No][dimension].data()[i] += (*vertex_collision_sum).b_sum[cloth_No][i][dimension];
	//	}
	//}
	
	//edge length
	for (int i = 0; i < mesh_struct.edges.size(); ++i) {
		b.data()[mesh_struct.edges[i].vertex[0]] += length_stiffness[i] * p_edge_length[i].data()[dimension];
		b.data()[mesh_struct.edges[i].vertex[1]] -= length_stiffness[i] * p_edge_length[i].data()[dimension];
	}
	//bending	
	VectorXd p0;
	if (!mesh_struct.faces.empty()) {
		for (int i = 0; i < cloth_sys_size[cloth_No]; ++i) {		
			p0 = vertex_lbo[i] * (p_bending[i].data()[dimension] * ((*cloth)[cloth_No].bend_stiffness * lbo_weight[i]));
			for (int j = 0; j < vertex_lbo[i].size(); ++j) {
				b.data()[vertex_around_vertex_for_bending[i][j]] += p0.data()[j];
			}
		}
	}
	//position
	for (int i = 0; i < (*cloth)[cloth_No].mesh_struct.anchor_vertex.size(); ++i) {
		b.data()[mesh_struct.anchor_vertex[i]] += (*cloth)[cloth_No].position_stiffness * mesh_struct.anchor_position[i][dimension];
	}
	b += (1.0 / (sub_time_step * sub_time_step)) * (cloth_mass[cloth_No].cwiseProduct(u_prediction));
	u= cloth_llt[cloth_No].solve(b);
	p0 = u - u_prediction;
	VectorXd p1;
	p1 = p0.cwiseProduct(cloth_mass[cloth_No]);
	temEnergy[thread_id] += 0.5 / (sub_time_step * sub_time_step) * p0.dot(p1);
}

void ProjectDynamic::solveClothSystemPerThead(int thread_id)
{
	int cloth_No;
	int dimension;
	temEnergy[thread_id] = 0.0;
	for (int i = 0; i < cloth_dimension_per_thread[thread_id].size(); ++i) {
		cloth_No = cloth_dimension_per_thread[thread_id][i] / 3;
		dimension = cloth_dimension_per_thread[thread_id][i] % 3;
		solveClothSystemPerThead(cloth_b[cloth_No][dimension], cloth_u[cloth_No][dimension], (*cloth)[cloth_No].mesh_struct, (*cloth)[cloth_No].length_stiffness,
			p_edge_length[cloth_No], cloth_No, dimension, vertex_around_vertex_for_bending[cloth_No], vertex_lbo[cloth_No], p_bending[cloth_No],
			lbo_weight[cloth_No], cloth_u_prediction[cloth_No][dimension], thread_id);
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