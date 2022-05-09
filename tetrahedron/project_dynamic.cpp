#include"project_dynamic.h"
#include"basic/EigenMatrixIO.h"

ProjectDynamic::ProjectDynamic()
{
	gravity_ = 9.8;
	total_thread_num = std::thread::hardware_concurrency();
	temEnergy.resize(total_thread_num);
	outer_itr_conv_rate = 1.5e-3;// 7.5e-2; 
	local_global_conv_rate = 2e-3;
	sub_step_num = 1;

	use_dierct_solve_for_coarest_mesh = true;
	super_jacobi_step_size = 3;
	max_it = 50;
	max_jacobi_itr_num = 50;
	displacement_norm_thread.resize(total_thread_num);


	iteration_method.setConvergenceRate(1e-8, 100);
	max_inner_iteration_num = 7;

}

void ProjectDynamic::setForPD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
	Thread* thread,double* tolerance_ratio)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	sub_time_step = time_step / (double)sub_step_num;
	this->thread = thread;
	tetrahedron_begin_obj_index = cloth->size();
	//collision.draw_culling = draw_culling_;
	collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio);
	total_collider_num = collider->size();
	this->collider = collider;
	
	setSystemIndexInfo();
	initialPDvariable();
	iteration_method.setBasicInfo(sys_size, thread, &global_mat);
	setOverRelaxationCoefficient();
	setForClothPD(cloth);
	setForTetrahedronPD();	
	computeGlobalStepMatrix();
	computeGravity();
	setIndexPerThread();
	
	
	//iteration_method.test();
	//iteration_method.testRelativeError();
	
}



void ProjectDynamic::setOverRelaxationCoefficient()
{
	double jacobi_or = 0.5;
	double A_jacobi_or = jacobi_or;
	double A_jaocbi_3_or = jacobi_or;
	double gauss_seidel_or = 0.8;
	double jacobi_chebyshev_or = jacobi_or;
	double A_jacobi_chebyshev_or = jacobi_or;
	double A_jaocbi_3_chebyshev_or = jacobi_or;
	double weight_for_chebyshev_jacobi = 0.8;
	double weight_for_chebyshev_A_jacobi = 0.8;
	iteration_method.setOverRelaxationCoefficient(jacobi_or, A_jacobi_or, A_jaocbi_3_or, gauss_seidel_or, jacobi_chebyshev_or,
		A_jacobi_chebyshev_or, A_jaocbi_3_chebyshev_or, weight_for_chebyshev_jacobi, weight_for_chebyshev_A_jacobi);

}


void ProjectDynamic::updateSystemPos()
{
	unsigned int vertex_index_start;
	std::array<double, 3>* vertex_position;
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
		vertex_index_start = vertex_begin_per_cloth[j];
		vertex_position = (*cloth)[j].mesh_struct.vertex_position.data();
		for (unsigned int i = 0; i < cloth_sys_size[j]; ++i) {
			u[0].data()[i + vertex_index_start] = vertex_position[i][0];
			u[1].data()[i + vertex_index_start] = vertex_position[i][1];
			u[2].data()[i + vertex_index_start] = vertex_position[i][2];
		}
	}
	for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
		vertex_position = (*tetrahedron)[j].mesh_struct.vertex_position.data();
		vertex_index_start = vertex_begin_per_tetrahedron[j];
		for (unsigned int i = 0; i < tetrahedron_sys_size[j]; ++i) {
			u[0].data()[i + vertex_index_start] = vertex_position[i][0];
			u[1].data()[i + vertex_index_start] = vertex_position[i][1];
			u[2].data()[i + vertex_index_start] = vertex_position[i][2];
		}
	}
	u_previous_itr = u;
	u_ = u;
	u_prediction = u;
	for (unsigned int i = 0; i < v.size(); ++i) {
		v[i].setZero();
		v_[i].setZero();
		acceleration[i].setZero();
	}
	TriangleMeshStruct* mesh_struct;
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
		mesh_struct = &(*cloth)[j].mesh_struct;
		thread->assignTask(mesh_struct, FACE_NORMAL_RENDER);
		thread->assignTask(mesh_struct, VERTEX_NORMAL_RENDER);
	}
	TetrahedronMeshStruct* mesh_struct_;
	for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
		mesh_struct_ = &(*tetrahedron)[j].mesh_struct;
		thread->assignTask(mesh_struct_, FACE_NORMAL_RENDER);
		thread->assignTask(mesh_struct_, VERTEX_NORMAL_RENDER);
	}
	for (unsigned int j = 0; j < total_collider_num; ++j) {
		mesh_struct = &(*collider)[j].mesh_struct;
		thread->assignTask(mesh_struct, FACE_NORMAL_RENDER);
		thread->assignTask(mesh_struct, VERTEX_NORMAL_RENDER);
	}
}



void ProjectDynamic::setSystemIndexInfo()
{
	vertex_begin_per_cloth.resize(cloth->size() + 1);
	vertex_begin_per_tetrahedron.resize(tetrahedron->size() + 1);
	vertex_begin_per_cloth[0] = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		vertex_begin_per_cloth[i + 1] = vertex_begin_per_cloth[i] + (*cloth)[i].mesh_struct.vertices.size();
	}
	vertex_begin_per_tetrahedron[0] = vertex_begin_per_cloth[vertex_begin_per_cloth.size() - 1];
	for (int i = 0; i < tetrahedron->size(); ++i){
		vertex_begin_per_tetrahedron[i + 1] = vertex_begin_per_tetrahedron[i] + (*tetrahedron)[i].mesh_struct.vertex_position.size();
	}
	sys_size = vertex_begin_per_tetrahedron[vertex_begin_per_tetrahedron.size()-1];
}

void ProjectDynamic::initialDHatTolerance(double ave_edge_length)
{
	collision.initialDHatTolerance(ave_edge_length);
	int element_count = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		element_count += (*cloth)[i].mesh_struct.vertex_for_render.size();
	}
	displacement_bound = 1e-3 * ave_edge_length;
	displacement_bound *= displacement_bound;
	displacement_bound *= (double)element_count;
}

void ProjectDynamic::setForClothPD(std::vector<Cloth>* cloth)
{

	std::vector<std::vector<double>> edge_cot_weight;
	computeEdgeCotWeight(edge_cot_weight);
	computeLBOWeight(edge_cot_weight);
	computeVertexLBO(edge_cot_weight);
	restBendingMeanCurvature();
	
}

void ProjectDynamic::setForTetrahedronPD()
{
	for (unsigned int i = 0; i < total_tetrahedron_num; ++i) {
		computeTetMass((*tetrahedron)[i].mesh_struct.mass.data(), mass_inv, mass, vertex_begin_per_tetrahedron[i], (*tetrahedron)[i].mesh_struct.vertex_for_render.size());
	}
}


void ProjectDynamic::computeTetMass(double* mass_, VectorXd& mass_inv, VectorXd& mass, int vertex_index_start, int vertex_num)
{
	for (int i = 0; i < vertex_num; ++i) {
		mass.data()[i + vertex_index_start] = mass_[i];
		mass_inv.data()[i + vertex_index_start] = 1.0 / mass_[i];
	}
}

void ProjectDynamic::restBendingMeanCurvature()
{
	rest_mean_curvature_norm.resize(total_cloth_num);
	int neighbor, ve;
	double vertex_curvature[3];
	for (unsigned int k = 0; k < total_cloth_num; ++k) {
		restBendingMeanCurvatureSingleCloth((*cloth)[k].mesh_struct, rest_mean_curvature_norm[k], vertex_lbo[k]);
	}
}

void ProjectDynamic::restBendingMeanCurvatureSingleCloth(TriangleMeshStruct& mesh_struct, std::vector<double>& rest_mean_curvature_norm,
	std::vector<VectorXd>& vertex_lbo)
{
	rest_mean_curvature_norm.resize(mesh_struct.vertices.size());
	VectorXd q;
	unsigned int size;
	double vertex_curvature[3];
	for (unsigned int i = 0; i < rest_mean_curvature_norm.size(); ++i) {
		size = vertex_lbo[i].size();
		q.resize(size);
		for (unsigned int j = 0; j < 3; ++j) {
			q[0] = mesh_struct.vertex_position[i][j];
			for (unsigned int h = 1; h < size; h++) {
				q[h] = mesh_struct.vertex_position[mesh_struct.vertices[i].neighbor_vertex[h-1]][j];
			}
			vertex_curvature[j] = dotProductX_(vertex_lbo[i], q);
		}
		rest_mean_curvature_norm[i] = sqrt(DOT(vertex_curvature, vertex_curvature));
	}
}





void ProjectDynamic::initialPDvariable()
{
	total_cloth_num = cloth->size();
	cloth_sys_size.resize(total_cloth_num);
	for (unsigned int i = 0; i < total_cloth_num; ++i) {
		cloth_sys_size[i] = (*cloth)[i].mesh_struct.vertices.size();
	}
	total_tetrahedron_num = tetrahedron->size();
	tetrahedron_sys_size.resize(total_tetrahedron_num);
	for (unsigned int i = 0; i < total_tetrahedron_num; ++i) {
		tetrahedron_sys_size[i] = (*tetrahedron)[i].mesh_struct.vertex_position.size();
	}
	u.resize(3);
	v.resize(3);
	for (unsigned int i = 0; i < 3; ++i) {
		u[i].resize(sys_size);	u[i].setZero();
		v[i].resize(sys_size); v[i].setZero();
	}
	initialClothPDvariable();
	initialTetrahedronPDvariable();

	u_previous_itr = u;
	u_ = u;
	v_ = v;
	b = v;
	u_prediction = u;
	acceleration = v;
}

void ProjectDynamic::initialTetrahedronPDvariable()
{
	tet_local_A = getARAPmatrix();
	unsigned int vertex_start;
	p_ARAP_volume_preserve.resize(total_tetrahedron_num);
	std::array<double, 3>* position;
	for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
		p_ARAP_volume_preserve[j].resize((*tetrahedron)[j].mesh_struct.indices.size());
		vertex_start = vertex_begin_per_tetrahedron[j];
		position = (*tetrahedron)[j].mesh_struct.vertex_position.data();
		for (unsigned int i = 0; i < 3; ++i) {
			for (unsigned int k = 0; k < tetrahedron_sys_size[j]; ++k) {
				u[i].data()[vertex_start + k] = position[k][i];
			}
		}
	}
}

void ProjectDynamic::initialClothPDvariable()
{
	unsigned int vertex_start;
	p_edge_length.resize(total_cloth_num);
	p_bending.resize(total_cloth_num);	
	std::array<double, 3>* position;
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
		p_edge_length[j].resize((*cloth)[j].mesh_struct.edges.size());
		p_bending[j].resize(cloth_sys_size[j]);
		vertex_start = vertex_begin_per_cloth[j];
		position= (*cloth)[j].mesh_struct.vertex_position.data();
		for (unsigned int i = 0; i < 3; ++i) {
			for (unsigned int k = 0; k < cloth_sys_size[j]; ++k) {
				u[i].data()[vertex_start + k] = position[k][i];
			}
		}
		for (unsigned int i = 0; i < cloth_sys_size[j]; ++i) {
			p_bending[j][i].resize(3);
		}
	}
}
void ProjectDynamic::setIndexPerThread()
{
	cloth_per_thread_begin.resize(total_thread_num + 1);
	if (total_cloth_num > 0) {
		arrangeIndex(total_thread_num, total_cloth_num, cloth_per_thread_begin.data());
	}
	cloth_dimension_per_thread.resize(total_thread_num);
	unsigned int total_num = total_cloth_num * 3;
	for (unsigned int i = 0; i < total_num; ++i) {
		cloth_dimension_per_thread[i % total_thread_num].push_back(i);
	}
	tetrahedron_per_thread_begin.resize(total_thread_num + 1);
	if (total_tetrahedron_num > 0) {
		arrangeIndex(total_thread_num, total_tetrahedron_num, tetrahedron_per_thread_begin.data());
	}
	tetrahedron_dimension_per_thread.resize(total_thread_num);
	total_num = total_tetrahedron_num * 3;
	for (unsigned int i = 0; i < total_num; ++i) {
		tetrahedron_dimension_per_thread[i % total_thread_num].push_back(i);
	}

	system_vertex_index_per_thread.resize(total_thread_num + 1);
	arrangeIndex(total_thread_num, sys_size, system_vertex_index_per_thread.data());

	dimension_per_thread.resize(total_thread_num + 1);
	arrangeIndex(total_thread_num, 3, dimension_per_thread.data());
}

void ProjectDynamic::computeVertexLBO(std::vector<std::vector<double>>& edge_cot_weight)
{
	vertex_lbo.resize(total_cloth_num);
	for (unsigned int k = 0; k < total_cloth_num; ++k) {
		computeVertexLBOSingleCloth((*cloth)[k].mesh_struct, vertex_lbo[k], edge_cot_weight[k], lbo_weight[k]);
	}
}




void ProjectDynamic::computeGlobalStepMatrix()
{
	global_mat.resize(sys_size, sys_size);
	std::vector<Triplet<double>>global_mat_nnz;
	unsigned int estimate_total_constraint = 40 * sys_size; //to estimate the size of global_mat_nnz
	global_mat_nnz.reserve(estimate_total_constraint);
	//set cloth
	for (unsigned int i = 0; i < total_cloth_num; ++i) {
		computeGlobalStepMatrixSingleCloth((*cloth)[i].mesh_struct, global_mat_nnz, lbo_weight[i], vertex_lbo[i], cloth_sys_size[i],
			(*cloth)[i].bend_stiffness, (*cloth)[i].length_stiffness, (*cloth)[i].position_stiffness, vertex_begin_per_cloth[i]);
	}
	//set tetrahedron
//	Matrix4d arap_matrix = tet_local_A.transpose() * tet_local_A;
	for (unsigned int i = 0; i < total_tetrahedron_num; ++i) {
		copmuteGlobalStepMatrixSingleTetrahedron((*tetrahedron)[i].mesh_struct, global_mat_nnz, tetrahedron_sys_size[i],
			(*tetrahedron)[i].ARAP_stiffness, (*tetrahedron)[i].volume_preserve_stiffness, (*tetrahedron)[i].position_stiffness,
			vertex_begin_per_tetrahedron[i]);
	}

	global_mat.setFromTriplets(global_mat_nnz.begin(), global_mat_nnz.end());

	
	//chech A
	SparseMatrix<double, RowMajor> A;
	A = global_mat;
	
	for (int i = 0; i < tetrahedron_sys_size[0]; ++i) {
		//std::cout << tetrahedron->data()[0].mesh_struct.mass[i] << std::endl;
	//	A.coeffRef(i, i) -= tetrahedron->data()[0].mesh_struct.mass[i] / (sub_time_step * sub_time_step);
	//	std::cout << tetrahedron->data()[0].mesh_struct.mass[i] / (sub_time_step * sub_time_step) << std::endl;
	}
	//std::cout.setf(std::ios::right);
	//std::cout.width(5);
	//std::cout << A << std::endl;
	//for (unsigned int i = 0; i < tetrahedron_sys_size[0]; ++i)
	//{
	//	std::cout << A.row(i).sum() << std::endl;
	//}


	global_mat_diagonal_ref_address.resize(sys_size);
	global_mat_diagonal_ref.resize(sys_size);
	for (unsigned int i = 0; i < sys_size; ++i) {
		global_mat_diagonal_ref_address[i] = &global_mat.coeffRef(i, i);
	}
	for (unsigned int i = 0; i < sys_size; ++i) {
		global_mat_diagonal_ref[i] = *(global_mat_diagonal_ref_address[i]);
	}
	ori_global_mat_diagonal_ref = global_mat_diagonal_ref;
	global_llt.analyzePattern(global_mat);
	global_llt.factorize(global_mat);
	collision_free_llt.analyzePattern(global_mat);
	collision_free_llt.factorize(global_mat);
	initial_global_mat = global_mat;


	std::vector<std::array<int, 2>> global_mat_coeff_index;
	std::vector<double> global_mat_coeff;
	global_mat_coeff_index.resize(global_mat_nnz.size());
	global_mat_coeff.resize(global_mat_nnz.size());
	for (int i = 0; i < global_mat_nnz.size(); ++i) {
		global_mat_coeff_index[i].data()[0] = global_mat_nnz[i].row();
		global_mat_coeff_index[i].data()[1] = global_mat_nnz[i].col();
		global_mat_coeff[i] = global_mat_nnz[i].value();
	}


	iteration_method.setOffDiagonal();
	iteration_method.initialGlobalDiagonalInv(&global_mat_diagonal_ref_address);
	iteration_method.initialJacobi();

	iteration_method.createAJacobiOperator(global_mat_coeff_index, global_mat_coeff);


	
}

void ProjectDynamic::copmuteGlobalStepMatrixSingleTetrahedron(TetrahedronMeshStruct& mesh_struct, std::vector<Triplet<double>>& global_mat_nnz, int sys_size, double& ARAP_stiffness,
	double& volume_preserve_stiffness, double position_stiffness, int vertex_index_start)
{
	//ARAP + volume preserve: they have the same matrix.
	//Matrix4d m_ARAP = (ARAP_stiffness+volume_preserve_stiffness) * m_for_ARAP;//
	int index[4];

	//std::cout << "tet " << mesh_struct.indices.size() << std::endl;

	Matrix4d AT_A;


	for (unsigned int i = 0; i < mesh_struct.indices.size(); ++i) {
		memcpy(index, mesh_struct.indices[i].data(), 16);
		for (unsigned int j = 0; j < 4; ++j){
			index[j] += vertex_index_start;
		}

		AT_A = ((ARAP_stiffness + volume_preserve_stiffness) * mesh_struct.volume[i]) * mesh_struct.A[i].transpose() * mesh_struct.A[i];

		for (unsigned int j = 0; j < 4; ++j) {
			for (unsigned int k = j + 1; k < 4; ++k) {
				global_mat_nnz.push_back(Triplet<double>(index[j], index[k], AT_A.data()[4 * j + k]));//mesh_struct.volume[i] *
				global_mat_nnz.push_back(Triplet<double>(index[k], index[j], AT_A.data()[4 * k + j]));// mesh_struct.volume[i] *
			}
			global_mat_nnz.push_back(Triplet<double>(index[j], index[j], AT_A.data()[4 * j + j]));//mesh_struct.volume[i] *

		}
		//std::cout << mesh_struct.volume[i] << std::endl;
	}
	//std::cout << "tet " << mesh_struct.indices.size() << std::endl;
	
	//position
	//for (int i = 0; i < mesh_struct.anchor_vertex.size(); ++i) {
	//	global_mat_nnz.push_back(Triplet<double>(mesh_struct.anchor_vertex[i] + vertex_index_start, mesh_struct.anchor_vertex[i] + vertex_index_start, position_stiffness));
	//}
	//mass
	for (int i = 0; i < sys_size; ++i) {
		global_mat_nnz.push_back(Triplet<double>(i + vertex_index_start, i + vertex_index_start, mesh_struct.mass[i] / (sub_time_step * sub_time_step)));
		//std::cout << mesh_struct.mass[i] / (sub_time_step * sub_time_step) << std::endl;
	}
}


void ProjectDynamic::updateTetrahedronAnchorVertices()
{
	for (unsigned int i = 0; i < total_tetrahedron_num; ++i)	{
		updateTetrahedronAnchorVertices(i, (*tetrahedron)[i].mesh_struct, vertex_begin_per_tetrahedron[i], (*tetrahedron)[i].position_stiffness);
	}
	global_llt.factorize(global_mat);
	collision_free_llt.factorize(global_mat);

	for (unsigned int i = 0; i < total_tetrahedron_num; ++i) {
		iteration_method.updateDiagonalWithAnchorVertices((*tetrahedron)[i].mesh_struct.anchor_vertex.size(), (*tetrahedron)[i].mesh_struct.anchor_vertex.data(),
			tetrahedron_sys_size[i], vertex_begin_per_tetrahedron[i], (*tetrahedron)[i].position_stiffness);
	}
	iteration_method.updateDiagonalWithAnchorVerticesTotal();
}


void ProjectDynamic::updateTetrahedronAnchorVertices(int tetrahedron_index,TetrahedronMeshStruct& mesh_struct, int vertex_index_start,
	double position_stiffness)
{
	unsigned int anchor_vertex_size = mesh_struct.anchor_vertex.size();
	int* anchor_vertex = mesh_struct.anchor_vertex.data();
	unsigned int system_size = tetrahedron_sys_size[tetrahedron_index];
	for (unsigned int i = 0; i < system_size; ++i){
		*(global_mat_diagonal_ref_address[i + vertex_index_start]) = ori_global_mat_diagonal_ref[i + vertex_index_start];
	}

	for (unsigned int i = 0; i < anchor_vertex_size; ++i){
		*(global_mat_diagonal_ref_address[anchor_vertex[i] + vertex_index_start]) += position_stiffness;
	}
	for (unsigned int i = 0; i < system_size; ++i) {
		global_mat_diagonal_ref[i + vertex_index_start] = *(global_mat_diagonal_ref_address[i + vertex_index_start]);
	}	



}

void ProjectDynamic::computeGlobalStepMatrixSingleCloth(TriangleMeshStruct& mesh_struct, std::vector<Triplet<double>>& global_mat_nnz,
	std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo, int sys_size, double bending_stiffness, std::vector<double>& length_stiffness,
	double position_stiffness, int vertex_index_start)
{	
	//bending
	unsigned int lbo_length; MatrixXd ATA;
	std::vector<int> index;
	if (!mesh_struct.faces.empty()) {
		for (unsigned int i = 0; i < sys_size; ++i) {
			lbo_length = vertex_lbo[i].size();
			index.clear();
			index.resize(lbo_length);
			index[0] = i;
			memcpy(&(index[1]), mesh_struct.vertices[i].neighbor_vertex.data(), 4 * (lbo_length - 1));		
			for (unsigned int j = 0; j < index.size(); ++j)
			{
				index[j] += vertex_index_start;
			}
			ATA = (bending_stiffness * lbo_weight[i] * vertex_lbo[i]) * vertex_lbo[i].transpose();
			for (unsigned int l = 0; l < lbo_length; ++l) {
				for (unsigned int k = l + 1; k < lbo_length; ++k) {
					global_mat_nnz.push_back(Triplet<double>(index[l], index[k], ATA.data()[lbo_length * l + k]));
					global_mat_nnz.push_back(Triplet<double>(index[k], index[l], ATA.data()[lbo_length * k + l]));
				}
				global_mat_nnz.push_back(Triplet<double>(index[l], index[l], ATA.data()[lbo_length * l + l]));
			}
		}
	}
	//edge length
	int id0, id1;
	for (int i = 0; i < mesh_struct.edges.size(); ++i) {
		id0 = mesh_struct.edge_vertices[i << 1] + vertex_index_start;
		id1 = mesh_struct.edge_vertices[(i << 1) + 1] + vertex_index_start;
		global_mat_nnz.push_back(Triplet<double>(id0, id0, length_stiffness[i]));
		global_mat_nnz.push_back(Triplet<double>(id1, id1, length_stiffness[i]));
		global_mat_nnz.push_back(Triplet<double>(id0, id1, -length_stiffness[i]));
		global_mat_nnz.push_back(Triplet<double>(id1, id0, -length_stiffness[i]));
	}

	//position
	for (int i = 0; i < mesh_struct.anchor_vertex.size(); ++i) {
		global_mat_nnz.push_back(Triplet<double>(mesh_struct.anchor_vertex[i] + vertex_index_start, mesh_struct.anchor_vertex[i] + vertex_index_start, position_stiffness));
	}
	//mass
	for (int i = 0; i < sys_size; ++i) {
		global_mat_nnz.push_back(Triplet<double>(i+ vertex_index_start, i+ vertex_index_start, mesh_struct.mass[i] / (sub_time_step * sub_time_step)));
	}
}

void ProjectDynamic::computeVertexLBOSingleCloth(TriangleMeshStruct& mesh_struct, std::vector<VectorXd>& vertex_lbo, std::vector<double>& edge_cot_weight,
	std::vector<double>& lbo_weight)
{
	double total;
	double edge_weight;
	vertex_lbo.resize(mesh_struct.vertices.size());
	for (int i = 0; i < mesh_struct.vertices.size(); ++i) {
		vertex_lbo[i].resize(mesh_struct.vertices[i].edge.size() + 1);
		vertex_lbo[i].setZero();
		total = 0.0;
		for (int j = 0; j < mesh_struct.vertices[i].edge.size(); ++j) {
			edge_weight = edge_cot_weight[mesh_struct.vertices[i].edge[j]] / lbo_weight[i];
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
	std::vector<unsigned int>* around_edge_index;
	for (int i = 0; i < mesh_struct.vertices.size(); ++i) {
		around_edge_index = &mesh_struct.vertices[i].edge;
		edge_around_vertex_for_bending[i].reserve(around_edge_index->size());
		vertex_around_vertex_for_bending[i].reserve(around_edge_index->size() + 1);
		vertex_around_vertex_for_bending[i].push_back(i);
		for (int j = 0; j < around_edge_index->size(); ++j) {
			if (mesh_struct.edges[(*around_edge_index)[j]].opposite_vertex.size() == 2) {
				edge_around_vertex_for_bending[i].push_back((*around_edge_index)[j]);
				if (i == mesh_struct.edge_vertices[(*around_edge_index)[j] << 1]) {
					vertex_around_vertex_for_bending[i].push_back(mesh_struct.edge_vertices[((*around_edge_index)[j] << 1) + 1]);
				}
				else {
					vertex_around_vertex_for_bending[i].push_back(mesh_struct.edge_vertices[(*around_edge_index)[j] << 1]);
				}
			}
		}
	}
}

void ProjectDynamic::computeGravity()
{
	f_ext.resize(3);
	total_gravity.resize(3);
	for (int i = 0; i < 3; ++i) {
		total_gravity[i].resize(sys_size);
		total_gravity[i].setZero();
	}
	//double gravity_accerlation[3] = { 0,0.0,gravity_};
	//double gravity_accerlation[3] = { gravity_, 0,0.0};
	double gravity_accerlation[3] = {0.0, -gravity_, 0.0};
	std::vector<double>* mass_;
	unsigned int vertex_index_start;

	//set cloth
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
		mass_ = &(*cloth)[j].mesh_struct.mass;
		vertex_index_start = vertex_begin_per_cloth[j];
		for (unsigned int i = 0; i < 3; ++i) {
			if (cloth_sys_size[j] > 1) {			
				for (unsigned int k = 0; k < cloth_sys_size[j]; ++k) {
					total_gravity[i][vertex_index_start + k] = gravity_accerlation[i] * (*mass_)[k];
				}
			}
		}
	}

	//set tetrahedron
	for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
		mass_ = &(*tetrahedron)[j].mesh_struct.mass;
		vertex_index_start = vertex_begin_per_tetrahedron[j];
		for (unsigned int i = 0; i < 3; ++i) {
			if (tetrahedron_sys_size[j] > 1) {
				for (unsigned int k = 0; k < tetrahedron_sys_size[j]; ++k) {
					total_gravity[i][vertex_index_start + k] = gravity_accerlation[i] * (*mass_)[k];
				}
			}
		}
	}

	f_ext = total_gravity;
}



void ProjectDynamic::computeLBOWeightSingleCloth(std::vector<double>& edge_cot_weight, std::vector<double>& lbo_weight,
	TriangleMeshStruct& mesh_struct, VectorXd& mass_inv, double density, VectorXd& mass, int vertex_index_start)
{
	double m;
	for (int i = 0; i < mesh_struct.faces.size(); ++i) {
		m = mesh_struct.faces[i].area / 3.0;
		lbo_weight[mesh_struct.triangle_indices[i][0]] += m;
		lbo_weight[mesh_struct.triangle_indices[i][1]] += m;
		lbo_weight[mesh_struct.triangle_indices[i][2]] += m;
	}
	double mass_;
	memcpy(mass.data() + vertex_index_start, mesh_struct.mass.data(), 8 * mesh_struct.vertices.size());
	memcpy(mass_inv.data() + vertex_index_start, mesh_struct.mass_inv.data(), 8 * mesh_struct.vertices.size());
	//if (!mesh_struct.faces.empty()) 
	//	for (int i = 0; i < mesh_struct.vertices.size(); ++i) {
	//		mass_ = density * lbo_weight[i];
	//		mass_inv.data()[i + vertex_index_start] = 1.0 / mass_;
	//		mass.data()[i + vertex_index_start] = mass_;
	//	}
	//}
	//else {
	//	for (int i = 0; i < mesh_struct.vertices.size(); ++i) {
	//		mass_ = 1.25;
	//		mass_inv.data()[i + vertex_index_start] = 1.0 / mass_;
	//		mass.data()[i + vertex_index_start] = mass_;
	//	}
	//}
}

void ProjectDynamic::computeLBOWeight(std::vector<std::vector<double>>& edge_cot_weight)
{

	lbo_weight.resize(total_cloth_num);
	mass_inv.resize(sys_size);
	mass.resize(sys_size);
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
		lbo_weight[j].resize(cloth_sys_size[j], 0.0);
		computeLBOWeightSingleCloth(edge_cot_weight[j], lbo_weight[j], (*cloth)[j].mesh_struct, mass_inv, (*cloth)[j].density, mass, vertex_begin_per_cloth[j]);
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
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
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
	unsigned int edge_vertex_0;
	unsigned int edge_vertex_1;
	unsigned int opposite_0;
	unsigned int opposite_1;
	for (int i = 0; i < edge_num; ++i) {
		if (mesh_struct.edges[i].opposite_vertex.size() > 1) {
			edge_vertex_0 = mesh_struct.edge_vertices[i << 1];
			edge_vertex_1 = mesh_struct.edge_vertices[(i << 1) + 1];
			opposite_0 = mesh_struct.edges[i].opposite_vertex[0];
			opposite_1 = mesh_struct.edges[i].opposite_vertex[1];
			cotan0 = 0;
			cotan1 = 0;
			SUB(x10, mesh_struct.vertex_position[edge_vertex_0], mesh_struct.vertex_position[opposite_0]);
			SUB(x20, mesh_struct.vertex_position[edge_vertex_1], mesh_struct.vertex_position[opposite_0]);
			len10 = sqrt(DOT(x10, x10));
			len20 = sqrt(DOT(x20, x20));
			theta0 = acos(DOT(x10, x20) / (len10 * len20));
			cotan0 = 1.0 / tan(theta0);

			SUB(x13, mesh_struct.vertex_position[edge_vertex_0], mesh_struct.vertex_position[opposite_1]);
			SUB(x23, mesh_struct.vertex_position[edge_vertex_1], mesh_struct.vertex_position[opposite_1]);
			len13 = sqrt(DOT(x13, x13));
			len23 = sqrt(DOT(x23, x23));
			theta1 = acos(DOT(x13, x23) / (len13 * len23));
			cotan1 = 1.0 / tan(theta1);
			edge_cot_weight[i] = -0.5 * (cotan0 + cotan1);
		}
		else
		{
			edge_cot_weight[i] = 0.0;
		}
	}
}

void ProjectDynamic::reset()
{
	int vertex_index_start;
	//cloth
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
		vertex_index_start = vertex_begin_per_cloth[j];
		for (unsigned int i = 0; i < cloth_sys_size[j]; ++i) {
			for (unsigned int k = 0; k < 3; ++k) {
				u[k].data()[vertex_index_start + i] = (*cloth)[j].ori_vertices[i][k];
			}
		}		
	}
	//tetrahedron
	for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
		vertex_index_start = vertex_begin_per_tetrahedron[j];
		for (unsigned int i = 0; i < tetrahedron_sys_size[j]; ++i) {
			for (unsigned int k = 0; k < 3; ++k) {
				u[k].data()[vertex_index_start + i] = (*tetrahedron)[j].ori_vertices[i][k];
			}
		}
	}
	f_ext = total_gravity;
	u_ = u;
	for (unsigned int i = 0; i < 3; ++i) {
		v[i].setZero();
		v_[i].setZero();
		acceleration[i].setZero();
	}
}

void ProjectDynamic::initial()
{
	global_mat = initial_global_mat;
	global_llt.factorize(global_mat);
	collision_free_llt.factorize(global_mat);
	for (unsigned int i = 0; i < sys_size; ++i) {
		global_mat_diagonal_ref_address[i] = &global_mat.coeffRef(i, i);
	}
	global_mat_diagonal_ref = ori_global_mat_diagonal_ref;
	
	iteration_method.initialRecordDiagonal_Operator();
	reset();
}

void ProjectDynamic::updateModelPosition()
{
	unsigned int vertex_index_start;
	std::array<double, 3>* vertex_position;
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
		vertex_index_start = vertex_begin_per_cloth[j];
		vertex_position = (*cloth)[j].mesh_struct.vertex_position.data();
		for (unsigned int i = 0; i < cloth_sys_size[j]; ++i) {
			vertex_position[i][0] = u[0].data()[i + vertex_index_start];
			vertex_position[i][1] = u[1].data()[i + vertex_index_start];
			vertex_position[i][2] = u[2].data()[i + vertex_index_start];
		}
	}
	for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
		vertex_position = (*tetrahedron)[j].mesh_struct.vertex_position.data();
		vertex_index_start = vertex_begin_per_tetrahedron[j];
		for (unsigned int i = 0; i < tetrahedron_sys_size[j]; ++i) {
			vertex_position[i][0] = u[0].data()[i + vertex_index_start];
			vertex_position[i][1] = u[1].data()[i + vertex_index_start];
			vertex_position[i][2] = u[2].data()[i + vertex_index_start];
		}
	}
}


void ProjectDynamic::updateRenderPositionIPC()
{
	TriangleMeshStruct* mesh_struct;
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
		mesh_struct = &(*cloth)[j].mesh_struct;
		thread->assignTask(mesh_struct, FACE_NORMAL_RENDER);
		thread->assignTask(mesh_struct, VERTEX_NORMAL_RENDER);
	}
	for (unsigned int j = 0; j < total_collider_num; ++j) {
		mesh_struct = &(*collider)[j].mesh_struct;
		thread->assignTask(mesh_struct, FACE_NORMAL_RENDER);
		thread->assignTask(mesh_struct, VERTEX_NORMAL_RENDER);
		mesh_struct->ori_face_normal_for_render = mesh_struct->ori_face_normal;
	}
	TetrahedronMeshStruct* mesh_struct_;
	for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
		mesh_struct_ = &(*tetrahedron)[j].mesh_struct;
		thread->assignTask(mesh_struct_, FACE_NORMAL_RENDER);
		thread->assignTask(mesh_struct_, VERTEX_NORMAL_RENDER);
	}
	
}

void ProjectDynamic::updateRenderPosition()
{
	TriangleMeshStruct* mesh_struct;
	TetrahedronMeshStruct* mesh_struct_;
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
		mesh_struct = &(*cloth)[j].mesh_struct;
		mesh_struct->vertex_for_render = mesh_struct->vertex_position;
		mesh_struct->face_normal_for_render = mesh_struct->face_normal;
		thread->assignTask(mesh_struct, VERTEX_NORMAL);
	}
	for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
		mesh_struct_ = &(*tetrahedron)[j].mesh_struct;
		mesh_struct_->vertex_for_render = mesh_struct_->vertex_position;
		mesh_struct_->face_normal_for_render = mesh_struct_->face_normal;
		thread->assignTask(mesh_struct_, VERTEX_NORMAL);
	}
	for (unsigned int j = 0; j < total_collider_num; ++j) {
		mesh_struct = &(*collider)[j].mesh_struct;
		mesh_struct->vertex_for_render = mesh_struct->vertex_position;
		mesh_struct->face_normal_for_render = mesh_struct->face_normal;
		thread->assignTask(mesh_struct, VERTEX_NORMAL);
	}
}

void ProjectDynamic::PDsetPosPredict()
{
	PDPredict();
	updateModelPosition();
}

void ProjectDynamic::PDPredict()
{
	
	for (unsigned int i = 0; i < 3; ++i) {
		u[i] = u_[i] +sub_time_step * v_[i] + (0.25 * sub_time_step * sub_time_step) * (acceleration[i]);//sub_time_step * sub_time_step * mass_inv[j].cwiseProduct(f_ext[j][i]);
		u_prediction[i] = u_[i] +sub_time_step * v_[i] + (sub_time_step * sub_time_step) * mass_inv.cwiseProduct(f_ext[i]);
		//u[i] = u_prediction[i];
	}
	
}


void ProjectDynamic::firstPDForIPC(bool& record_matrix)
{
	//face_normal_render
	PDsetPosPredict();
	
	previous_PD_energy = 1e-15;
	current_PD_energy= 1e-15;

	local_global_iteration_num = 0;

	iteration_method.setOperatorCollisionFree();

	while (!innerIterationConvergeCondition()) {
		thread->assignTask(this, LOCAL_PROJECTION);
		current_constraint_energy = temEnergy[0];
		for (unsigned int i = 1; i < total_thread_num; ++i) {
			current_constraint_energy += temEnergy[i];
		}
		if (record_matrix) {
			//	std::string matrix_name = "./save_matrix/global_" + std::to_string(*time_stamp) + ".dat";
			//	saveSparseMatrix(cloth_global_mat[0], matrix_name);
			//
			////std::string off_diagonal_name = "off_diagonal_" + std::to_string(*time_stamp) + ".dat";
			////saveSparseMatrix(iteration_method.off_diagonal[0], off_diagonal_name);
			//std::string u_name = "./save_matrix/u_" + std::to_string(*time_stamp) + "_";
			//for (int i = 0; i < 3; ++i) {
			//	saveMatrix(i, cloth_u[0][i], u_name);
			//}
		}
		
		thread->assignTask(this, CONSTRUCT_B_WITHOUT_COLLISION);
		thread->assignTask(this, SOLVE_WITHOUT_COLLISION);
		thread->assignTask(this, COMPUTE_ENERGY);

		computeInnerCollisionFreeEnergy();
		
		if (record_matrix) {
			//std::string b_name = "./save_matrix/b_" + std::to_string(*time_stamp) + "_";
			//for (int i = 0; i < 3; ++i) {
			//	saveMatrix(i, cloth_b[0][i], b_name);
			//}
			//record_matrix = false;
		}

		local_global_iteration_num++;

		
	}

	updateModelPosition();

	//collision.collisionCulling();

	current_PD_energy = temEnergy[0];
	for (unsigned int i = 1; i < total_thread_num; ++i) {
		current_PD_energy += temEnergy[i];
	}
	current_PD_energy += current_constraint_energy;
	current_collision_energy = 1e-15;
	////std::cout << "++++" << std::endl;
	//for (int i = 0; i < cloth_sys_size[0]; ++i) {
	//	//std::cout << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "
	//		<< cloth_u[0][2][i] << std::endl;
	//}
	//std::cout << "first local-global without collision" << std::endl;
	for (unsigned int i = 0; i < cloth_sys_size[0]; ++i) {
		//std::cout << "    " << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "<< cloth_u[0][2][i] << std::endl;
	}
}


void ProjectDynamic::PD_IPC_solve(bool& record_matrix)
{
	collision.collisionCulling();
	firstPDForIPC(record_matrix);
//	outer_iteration_num = 1;
//	reset_ave_iteration_record();
//	while (!IPC_PDConvergeCondition()) {
//		////std::cout << "==" << std::endl;
		collision.globalCollisionTime();
		thread->assignTask(this, COLLISION_FREE_POSITION);//in document, we use q_n+1, however, here, we use vertices_for_render & cloth_u to store this collision free position.
//		u_previous_itr = u;
		collision.solveCollisionConstraint();
//		PDupdateSystemMatrix();
//		//std::cout << "==iteration number " << outer_iteration_num << std::endl;
//		//std::cout << "collision free position " << std::endl;
//		//for (int i = 0; i < cloth_sys_size[0]; ++i) {
//		//	std::cout<<"    " << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "<< cloth_u[0][2][i] << std::endl;
//		//}
//		local_global_iteration_num = 0;
//
////		if (record_matrix && outer_iteration_num==1) {
////			std::string matrix_name = "./save_matrix/global_" + std::to_string(*time_stamp) + ".dat";
////			//saveSparseMatrix(cloth_global_mat[0], matrix_name);
////			//std::string off_diagonal_name = "off_diagonal_" + std::to_string(*time_stamp) + ".dat";
////			//saveSparseMatrix(iteration_method.off_diagonal[0], off_diagonal_name);
////			std::string u_name = "./save_matrix/u_" + std::to_string(*time_stamp) + "_";
/////*			for (int i = 0; i < 3; ++i) {
////				saveMatrix(i, cloth_u[0][i], u_name);
////			}	*/		
////		}
//
//		while (!innerIterationConvergeCondition()) {
		collision.collisionEnergy();
//			thread->assignTask(this, LOCAL_PROJECTION_WITHOUT_ENERGY);
//			current_constraint_energy = temEnergy[0];
//			for (unsigned int i = 1; i < total_thread_num; ++i) {
//				current_constraint_energy += temEnergy[i];
//			}
//			thread->assignTask(this, CONSTRUCT_B);
//			thread->assignTask(this, SOLVE_WITH_COLLISION);
//			solveClothSystem2(false);
//			local_global_iteration_num++;
//			computeInnerEnergyIPCPD();
//
//			//if (record_matrix && outer_iteration_num == 1 && local_global_iteration_num==1) {
//			//	std::string b_name = "./save_matrix/b_" + std::to_string(*time_stamp) + "_";
//			//	for (int i = 0; i < 3; ++i) {
//			//		//saveMatrix(i, cloth_b[0][i], b_name);
//			//	}
//			//	//record_matrix = false;
//			//}
//		}
//		//std::cout << outer_iteration_num << std::endl;
//		updateModelPosition();
//		outer_iteration_num++;
//		//computeEnergyIPCPD();
//		displacement_ratio_dif = previous_displacement_norm - displacement_norm;
//		previous_displacement_norm = displacement_norm;
//		//system("pause");
//		//std::cout << "pd position " << std::endl;
//		//for (int i = 0; i < cloth_sys_size[0]; ++i) {
//		//	std::cout << "    " << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "	<< cloth_u[0][2][i] << std::endl;
//		//}
//		//std::cout << "+++++" << std::endl;
//	}
//	//collision.globalCollisionTime();
//	thread->assignTask(this, COLLISION_FREE_POSITION);
//	//std::cout << "final collision free position " << std::endl;
//	//for (int i = 0; i < cloth_sys_size[0]; ++i) {
//	//	std::cout << "    " << cloth_u[0][0][i] << " " << cloth_u[0][1][i] << " "	<< cloth_u[0][2][i] << std::endl;
//	//}
//	thread->assignTask(this, UPDATE_UV);
//	updateRenderPositionIPC();
//	//std::cout << cloth_v[0][1] << std::endl;
//	//std::cout << "========" << std::endl;

		std::cout << tetrahedron->data()[0].mesh_struct.vertex_for_render[0][0] << " " <<
			tetrahedron->data()[0].mesh_struct.vertex_for_render[0][1] << " " <<
			tetrahedron->data()[0].mesh_struct.vertex_for_render[0][2] << std::endl;
		std::cout << tetrahedron->data()[0].mesh_struct.vertex_position[0][0] << " " <<
			tetrahedron->data()[0].mesh_struct.vertex_position[0][1] << " " <<
			tetrahedron->data()[0].mesh_struct.vertex_position[0][2] << std::endl;
}


void ProjectDynamic::computeInnerCollisionFreeEnergy()
{
	previous_PD_energy = current_PD_energy;
	current_PD_energy = temEnergy[0];
	for (unsigned int i = 1; i < total_thread_num; ++i) {
		current_PD_energy += temEnergy[i];
	}
	current_PD_energy += current_constraint_energy;
}


void ProjectDynamic::computeInnerEnergyIPCPD()
{
	previous_itr_PD_energy = current_PD_energy;
	current_collision_energy = 1e-15;
	for (unsigned int k = 0; k < total_cloth_num; ++k) {
		current_collision_energy += collision.obj_target_pos.collision_energy;
	}
	current_PD_energy = temEnergy[0];
	for (unsigned int i = 1; i < total_thread_num; ++i) {
		current_PD_energy += temEnergy[i];
	}
	current_PD_energy += current_constraint_energy + current_collision_energy;
}

void ProjectDynamic::computeEnergyIPCPD()
{
	PD_energy_dif = previous_PD_energy - current_PD_energy;
	collision_energy_dif = previous_collision_energy - current_collision_energy;
	previous_collision_energy = current_collision_energy;
	previous_PD_energy = current_PD_energy;
	current_collision_energy = 1e-15;
	for (unsigned int k = 0; k < total_cloth_num; ++k) {
		current_collision_energy += collision.obj_target_pos.collision_energy;
	}
	current_PD_energy = temEnergy[0];
	for (unsigned int i = 1; i < total_thread_num; ++i) {
		current_PD_energy += temEnergy[i];
	}
	current_PD_energy += current_constraint_energy + current_collision_energy;

}

void ProjectDynamic::PDsolve()
{
	PDsetPosPredict();
	//collision.collisionCulling();

	int itr_num = 0;
	outer_iteration_num = 0;
	initialEnergy();
	local_global_iteration_num = 0;

	//for (int i = 0; i < cloth->size(); ++i) {
	//	thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL_RENDER);
	//}
	//for (int i = 0; i < tetrahedron->size(); ++i) {
	//	thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL_RENDER);
	//}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		thread->assignTask(&(*collider)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
//	std::cout << "===" << std::endl;
	while (!PDConvergeCondition()) {	
		//collision.globalCollision()
		initialEnergyOuterInteration();
		local_global_itr_in_single_outer = 0;
		//std::cout << "outer===" << std::endl;
		while (!PDLocalGlobalConvergeCondition()) {
			//std::cout << "inner===" << std::endl;
			initialEnergyLocalGlobal();
			if (local_global_itr_in_single_outer == 0) {				
				collision.collisionCulling();
				collision.solveCollisionConstraintDCD();
				PDupdateSystemMatrix();
			}
			else {
				collision.reSolveCollisionConstraintDCD();
			}
			//time_t t = clock();
			localProjection();
			for (unsigned int i = 0; i < total_thread_num; ++i) {
				current_constraint_energy += temEnergy[i];
			}
			//thread->assignTask(this, CONSTRUCT_B_WITHOUT_COLLISION);
			//thread->assignTask(this, SOLVE_WITHOUT_COLLISION);
			thread->assignTask(this, CONSTRUCT_B);
			//time_t t = clock();
			//for (unsigned int i = 0; i < 10; ++i) {
				thread->assignTask(this, SOLVE_WITH_COLLISION);
				solveClothSystem2(true);
			//}
			//time_t t1 = clock();
			
			current_collision_energy = 1e-15;
			for (unsigned int k = 0; k < total_cloth_num; ++k) {
				current_collision_energy += collision.obj_target_pos.collision_energy;
			}
			current_constraint_energy += current_collision_energy;
			for (unsigned int i = 0; i < total_thread_num; ++i) {
				current_PD_energy += temEnergy[i];
			}
			current_PD_energy += current_constraint_energy;
			updateModelPosition();			
			local_global_itr_in_single_outer++;
			local_global_iteration_num++;
		}
		outer_iteration_num++;

	}
	//updateRenderPosition();
	thread->assignTask(this, UPDATE_UV);

	//std::cout << v[1] << std::endl;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL);
	}
	updateRenderPosition();

	//std::cout << "==========================" << std::endl;
}


//UPDATE_UV
void ProjectDynamic::updateUVPerThread(int thread_id)
{
	VectorXd* this_v;
	VectorXd* this_u;
	VectorXd* this_u_;
	VectorXd* this_v_;
	VectorXd* this_acceleration;
	for (unsigned int j = 0; j < 3; ++j)
	{
		this_v = &v[j]; this_u = &u[j]; this_u_ = &u_[j]; this_v_ = &v_[j];
		this_acceleration = &acceleration[j];
		for (unsigned int i = system_vertex_index_per_thread[thread_id]; i < system_vertex_index_per_thread[thread_id + 1]; ++i) {
			this_v->data()[i] = (this_u->data()[i] - this_u_->data()[i]) / sub_time_step;
			this_v->data()[i] *= 0.98;
			this_acceleration->data()[i] = (this_v->data()[i] - this_v_->data()[i]) / sub_time_step;
			//this_u_->data()[]
		}
		memcpy(this_u_->data() + system_vertex_index_per_thread[thread_id], this_u->data() + system_vertex_index_per_thread[thread_id],
			8 * (system_vertex_index_per_thread[thread_id + 1] - system_vertex_index_per_thread[thread_id]));
		memcpy(this_v_->data() + system_vertex_index_per_thread[thread_id], this_v->data() + system_vertex_index_per_thread[thread_id],
			8 * (system_vertex_index_per_thread[thread_id + 1] - system_vertex_index_per_thread[thread_id]));
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


void ProjectDynamic::testJacobiMatrix()
{
	std::cout << "reference " << std::endl;
	//std::cout << global_mat << std::endl;

	VectorXd a(sys_size);
	for (unsigned int i = 0; i < sys_size; ++i) {
		a[i] = *(global_mat_diagonal_ref_address[i]);
	}
	//std::cout << a << std::endl;
	for (unsigned int i = 0; i < sys_size; ++i) {
		a[i] = 1.0 / a[i];
	}
	SparseMatrix<double, RowMajor> test = global_mat;
	for (unsigned int i = 0; i < sys_size; ++i) {
		test.coeffRef(i, i) = 0;
	}
	std::cout << -1.0*a.asDiagonal() * test << std::endl;
	std::cout << a.asDiagonal()*b[0] << std::endl;
	std::cout << u[0] << std::endl;
	std::cout << "=============" << std::endl;
	


}


void ProjectDynamic::testJacobi()
{
	std::cout << "test for " << std::endl;
	std::cout << global_mat << std::endl;
	VectorXd a(sys_size);
	for (unsigned int i = 0; i < sys_size; ++i) {
		a[i] = *(global_mat_diagonal_ref_address[i]);
	}
	std::cout << a << std::endl;
	for (unsigned int i = 0; i < sys_size; ++i) {
		a[i] = 1.0 / a[i];
	}
	SparseMatrix<double, RowMajor> test = global_mat;
	for (unsigned int i = 0; i < sys_size; ++i) {
		test.coeffRef(i, i) = 0;
	}

	SparseMatrix<double, RowMajor> test_R_jacobi = -1.0 * a.asDiagonal() * test;
	Matrix4d test1 = Matrix4d(test_R_jacobi);
	std::cout << test1 << std::endl;
	JacobiSVD<Matrix4d> svd;
	svd.compute(test1);
	//	rotation_2 = svd.matrixU();

	Vector4d eigen_value = svd.singularValues();
	std::cout << eigen_value << std::endl;
}



void ProjectDynamic::PDupdateSystemMatrix()
{	
	switch (itr_solver_method)
	{
	case DIRECT_SOLVE: {
		//time_t t = clock();
		thread->assignTask(this, UPDATE_MATRIX);
		global_llt.factorize(global_mat);
		//time_t t1 = clock();
		//std::cout << "update matrix &factorize " << t1 - t << std::endl;
	}
		break;
	case JACOBI:
		thread->assignTask(this, UPDATE_MATRIX);
		iteration_method.updateJacobi();
		//testJacobiMatrix();
		break;
	case A_JACOBI:
		thread->assignTask(this, UPDATE_DIAGONAL);
		thread->assignTask(&iteration_method, UPDATE_JACOBI_OPERATOR);
		thread->assignTask(&iteration_method, UPDATE_2_A_JACOBI_ITR_MATRIX);
		//iteration_method.updateJacobi();
		break;
	case CHEBYSHEV_A_JACOBI:
		thread->assignTask(this, UPDATE_DIAGONAL);
		thread->assignTask(&iteration_method, UPDATE_JACOBI_OPERATOR);
		thread->assignTask(&iteration_method, UPDATE_2_A_JACOBI_ITR_MATRIX);
		//iteration_method.updateJacobi();
		iteration_method.estimateSuperJacobiEigenValue(u.data(), 2);
		break;
	case GAUSS_SEIDEL_CHEBYSHEV:
		thread->assignTask(this, UPDATE_MATRIX);
		iteration_method.estimateGaussSeidelEigenValue(u, global_mat);
		break;
	case CHEBYSHEV_JACOBI:
		thread->assignTask(this, UPDATE_MATRIX);
		iteration_method.updateJacobi();
		iteration_method.estimateJacobiEigenValue(u);
		break;
	case PCG:
		thread->assignTask(this, UPDATE_MATRIX);
		iteration_method.updateGlobalDiagonalInv();
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
	int vertex_index_start;
	for (unsigned int i = 0; i < total_cloth_num; ++i) {
		anchor_index = &(*cloth)[i].mesh_struct.anchor_vertex;
		vertex_index_start = vertex_begin_per_cloth[i];
		for (int j = 0; j < anchor_index->size(); ++j) {
			*(global_mat_diagonal_ref_address[(*anchor_index)[j]+ vertex_index_start]) += (*cloth)[i].position_stiffness;
		}
	}
	for (unsigned int i = 0; i < total_tetrahedron_num; ++i) {
		anchor_index = &(*tetrahedron)[i].mesh_struct.anchor_vertex;
		vertex_index_start = vertex_begin_per_tetrahedron[i];
		for (int j = 0; j < anchor_index->size(); ++j) {
			*(global_mat_diagonal_ref_address[(*anchor_index)[j] + vertex_index_start]) += (*tetrahedron)[i].position_stiffness;
		}
	}
}


//COLLISION_FREE_POSITION
void ProjectDynamic::computeCollisionFreePosition(int thread_No)
{
	unsigned int index_end;
	
	double collision_time = collision.collision_time;
	std::array<double, 3>* q_pre;
	std::array<double, 3>* q_end;
	double* u_x; double* u_y; double* u_z;
	u_x = u[0].data();
	u_y = u[1].data();
	u_z = u[2].data();
	unsigned int vertex_index_start;
	//cloth
	MeshStruct* mesh_struct;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		mesh_struct = &(*cloth)[i].mesh_struct;
		index_end = mesh_struct->vertex_index_begin_per_thread[thread_No + 1];
		q_end = mesh_struct->vertex_position.data();
		q_pre = mesh_struct->vertex_for_render.data();
		vertex_index_start = vertex_begin_per_cloth[i];
		for (unsigned int j = mesh_struct->vertex_index_begin_per_thread[thread_No]; j < index_end; ++j) {
			q_pre[j][0] += collision_time * (q_end[j][0] - q_pre[j][0]);
			q_pre[j][1] += collision_time * (q_end[j][1] - q_pre[j][1]);
			q_pre[j][2] += collision_time * (q_end[j][2] - q_pre[j][2]);
			u_x[j + vertex_index_start] = q_pre[j][0];
			u_y[j + vertex_index_start] = q_pre[j][1];
			u_z[j + vertex_index_start] = q_pre[j][2];
		}
	}

	//tetrahedron
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		mesh_struct = &(*tetrahedron)[i].mesh_struct;
		index_end = mesh_struct->vertex_index_begin_per_thread[thread_No + 1];
		q_end = mesh_struct->vertex_position.data();
		q_pre = mesh_struct->vertex_for_render.data();
		vertex_index_start = vertex_begin_per_tetrahedron[i];
		for (unsigned int j = mesh_struct->vertex_index_begin_per_thread[thread_No]; j < index_end; ++j) {
			q_pre[j][0] += collision_time * (q_end[j][0] - q_pre[j][0]);
			q_pre[j][1] += collision_time * (q_end[j][1] - q_pre[j][1]);
			q_pre[j][2] += collision_time * (q_end[j][2] - q_pre[j][2]);
			u_x[j + vertex_index_start] = q_pre[j][0];
			u_y[j + vertex_index_start] = q_pre[j][1];
			u_z[j + vertex_index_start] = q_pre[j][2];
		}
	}
}


//COMPUTE_DISPLACEMENT
void ProjectDynamic::computeDisplacement(int thread_No)
{
	unsigned int index_end;
	double* u_x; double* u_y; double* u_z;
	double* u_pre_x; double* u_pre_y; double* u_pre_z;
	double x, y, z;
	double* displacement_norm_ = &displacement_norm_thread[thread_No];
	*displacement_norm_ = 0;
	double displace_current;
	u_x = u[0].data();
	u_y = u[1].data();
	u_z = u[2].data();
	u_pre_x = u_previous_itr[0].data();
	u_pre_y = u_previous_itr[1].data();
	u_pre_z = u_previous_itr[2].data();
	for (unsigned int j = system_vertex_index_per_thread[thread_No]; j < system_vertex_index_per_thread[thread_No + 1]; ++j) {
		x = u_x[j] - u_pre_x[j];
		y = u_y[j] - u_pre_y[j];
		z = u_z[j] - u_pre_z[j];
		displace_current = x * x + y * y + z * z;
		//if (*displacement_norm_ < displace_current) {
		*displacement_norm_ += displace_current;
	}
}


//UPDATE_DIAGONAL
void ProjectDynamic::updateDiagonalPerThread(int thread_No)
{
	std::vector<int>* vertex_index_begin_per_thread;
	

	//update_cloth
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
		iteration_method.updateClothDiagonalPerThread((*cloth)[j].mesh_struct.vertex_index_begin_per_thread[thread_No],
			(*cloth)[j].mesh_struct.vertex_index_begin_per_thread[thread_No + 1], vertex_begin_per_cloth[j], collision.obj_target_pos.stiffness[j].data(),
			collision.obj_target_pos.need_update[j]);
	}
	//update tetrahedron
	unsigned int tet_index;
	for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
		tet_index = tetrahedron_begin_obj_index + j;
		iteration_method.updateTetDiagonalPerThread((*tetrahedron)[j].mesh_struct.vertex_index_on_surface_begin_per_thread[thread_No],
			(*tetrahedron)[j].mesh_struct.vertex_index_on_surface_begin_per_thread[thread_No + 1], vertex_begin_per_tetrahedron[j],
			collision.obj_target_pos.stiffness[tet_index].data(),collision.obj_target_pos.need_update[tet_index], (*tetrahedron)[j].mesh_struct.vertex_index_on_sureface.data());
	}
}



//UPDATE_MATRIX
void ProjectDynamic::updateMatrixPerThread(int thread_No)
{
	int collision_num = 0;
	unsigned int* vertex_index_begin_per_thread;
	bool* need_update;
	double* stiffness;
	double* diagonal_ref;
	double** diagonal_ref_address;
	unsigned int vertex_index_start;
	diagonal_ref = global_mat_diagonal_ref.data();
	diagonal_ref_address = global_mat_diagonal_ref_address.data();
	unsigned int vertex_end;
	//update_cloth
	for (unsigned int j = 0; j < total_cloth_num; ++j) {
		vertex_index_start = vertex_begin_per_cloth[j];
		need_update = collision.obj_target_pos.need_update[j];
		stiffness = collision.obj_target_pos.stiffness[j].data();
		vertex_index_begin_per_thread = (*cloth)[j].mesh_struct.vertex_index_begin_per_thread.data();		
		vertex_end = vertex_index_begin_per_thread[thread_No + 1];
		for (unsigned int i = vertex_index_begin_per_thread[thread_No]; i < vertex_end; ++i) {
			*(diagonal_ref_address[i + vertex_index_start]) = diagonal_ref[i + vertex_index_start];
			////std::cout << i << " " << *(diagonal_ref_address[i]) << std::endl;
			if (need_update[i]) {
				*(diagonal_ref_address[i + vertex_index_start]) += stiffness[i];
				////std::cout<<"update "<<i << " " << *(diagonal_ref_address[i]) << std::endl;
			}
		}
	}
	//update tetrahedron
	unsigned int* vertex_index_on_surface;
	unsigned int vertex_index;
	for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
		vertex_index_start = vertex_begin_per_tetrahedron[j];
		need_update = collision.obj_target_pos.need_update[j+tetrahedron_begin_obj_index];
		stiffness = collision.obj_target_pos.stiffness[j+ tetrahedron_begin_obj_index].data();
		vertex_index_begin_per_thread = (*tetrahedron)[j].mesh_struct.vertex_index_on_surface_begin_per_thread.data();
		vertex_index_on_surface = (*tetrahedron)[j].mesh_struct.vertex_index_on_sureface.data();
		vertex_end = vertex_index_begin_per_thread[thread_No + 1];
		for (unsigned int i = vertex_index_begin_per_thread[thread_No]; i < vertex_end; ++i) {
			vertex_index = vertex_index_on_surface[i];
			*(diagonal_ref_address[vertex_index + vertex_index_start]) = diagonal_ref[vertex_index + vertex_index_start];
			////std::cout << i << " " << *(diagonal_ref_address[i]) << std::endl;
			if (need_update[vertex_index]) {
				*(diagonal_ref_address[vertex_index + vertex_index_start]) += stiffness[vertex_index];
				////std::cout<<"update "<<i << " " << *(diagonal_ref_address[i]) << std::endl;
			}
		}
	}
}




bool ProjectDynamic::innerIterationConvergeCondition()
{
	//return local_global_iteration_num > max_inner_iteration_num;

	if (local_global_iteration_num > 10) {
		bool energy_changing = fabs(current_PD_energy - previous_PD_energy) / previous_PD_energy < local_global_conv_rate || current_PD_energy < 5e-15;
		if (energy_changing) {
			return true;
		}
	}
	return false;
}

bool ProjectDynamic::IPC_PDConvergeCondition()
{
	if (outer_iteration_num > 3) {
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
			for (unsigned int i = 1; i < total_thread_num; ++i) {
				//if (displacement_norm < displacement_norm_thread[i]) {
					displacement_norm += displacement_norm_thread[i];
				//}
			}
	

			bool ratio_changing = fabs(displacement_ratio_dif + (previous_displacement_norm - displacement_norm)) / displacement_norm < 1e-7;
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
	//std::cout << outer_iteration_num << std::endl;
	bool standard = (energy_satisfied || need_to_stop) && outer_iteration_num > 5;
	//if (outer_iteration_num > 0) {//
	if (standard) {//
		//std::cout << (current_PD_energy - previous_itr_PD_energy) / previous_itr_PD_energy << std::endl;
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
	bool standard = (energy_satisfied || need_to_stop) && local_global_itr_in_single_outer > 5;

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
	temEnergy[thread_id] = 0.0;
	//edge length
	localEdgeLengthProjectionPerThread(thread_id, with_energy);
	//bending
	localBendingProjectionPerThread(thread_id, with_energy);
	//anchor
	localPositionProjectionPerThread(thread_id, with_energy);
	//ARAP
	localARAPProjectionPerThread(thread_id, with_energy);
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
// 

void ProjectDynamic::localARAPProjectionPerThread(int thread_id, bool with_energy)
{
	unsigned int vertex_index_start;
	unsigned int index_end;
	unsigned int index_begin;
	Matrix<double, 3, 4> q;
	TetrahedronMeshStruct* mesh_struct;
	int* vertex_index;
	std::array<double,3>* vertex_pos;
	//Vector3d center;
	int temp_vertex_index;
	Matrix3d transform;

	Matrix3d P_inv;
	//double* det;

	//Matrix<double, 3, 3>* PT;

	//Matrix<double, 4, 3>* PT_pos;

	Matrix3d rotation_2;	
	Matrix<double,4,3>* ATAp_ARAP_volume_preserve;
	double ARAP_stiffness;
	//double volume_preserve_stiffness;

	double sigma_min;
	double sigma_max;
	Vector3d sigma;
	//Matrix3d transform_for_volume_perserve;

	//Vector3d singular_value;

	double determinant;
	//bool keep_transform;
	Matrix<double, 3, 3> q_e;
	//Matrix<double, 3, 3> p_e;

	//Matrix<double, 4, 3> temp_p;

	//Matrix<double, 3, 4> temp_P_;

	Matrix3d deformation_gradient;

	Vector3d eigen_value;
	JacobiSVD<Matrix3d> svd;
	Matrix<double, 3, 4> p;

	if (with_energy) {
		for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
			vertex_index_start = vertex_begin_per_tetrahedron[j];
			mesh_struct = &(*tetrahedron)[j].mesh_struct;
			index_end = mesh_struct->tetrahedron_index_begin_per_thread[thread_id + 1];
			index_begin = mesh_struct->tetrahedron_index_begin_per_thread[thread_id];
			ATAp_ARAP_volume_preserve = p_ARAP_volume_preserve[j].data();
			ARAP_stiffness = (*tetrahedron)[j].ARAP_stiffness;
			//volume_preserve_stiffness = (*tetrahedron)[j].volume_preserve_stiffness;
			sigma_min = (*tetrahedron)[j].sigma_limit[0];
			sigma_max = (*tetrahedron)[j].sigma_limit[1];
			
			//P_inv = (*tetrahedron)[j].mesh_struct.P_inv.data();


			for (unsigned int i = index_begin; i < index_end; ++i) {
				vertex_index = mesh_struct->indices[i].data();
				//set matrix q
				for (unsigned int k = 0; k < 4; ++k) {
					temp_vertex_index = vertex_index[k] + vertex_index_start;
					for (int l = 0; l < 3; ++l) {
						q.data()[3 * k + l] = u[l].data()[temp_vertex_index];
					}					
				}
				for (unsigned int j = 1; j < 4; ++j) 
				{
					q_e.col(j-1) = q.col(j) - q.col(0);	
				}			


				memcpy(P_inv.data(), mesh_struct->A[i].data() + 3, 72);

				deformation_gradient = q_e * P_inv.transpose();
				//transform = (q * PT_times_PPT_inv[i]);
				//transform = PT[i].transpose()* q_e.transpose();
	
				svd.compute(deformation_gradient, ComputeFullU | ComputeFullV);
			//	rotation_2 = svd.matrixU();

				eigen_value = svd.singularValues();
				determinant = eigen_value[0] * eigen_value[1] * eigen_value[2];

				if (determinant < 0) {
					eigen_value[2] = -eigen_value[2];
/*					rotation_2.data()[6] *= -1.0;
					rotation_2.data()[7] *= -1.0;
					rotation_2.data()[8] *= -1.0;			*/	
				}

				for (unsigned int j = 0; j < 3; ++j) {
					if (eigen_value[j]>0) {
						if (eigen_value[j] < sigma_min) {
							eigen_value[j] = sigma_min;
						}
						else if (eigen_value[j] > sigma_max) {
							eigen_value[j] = sigma_max;
						}
					}
					else {
						if (eigen_value[j] > -sigma_min) {
							eigen_value[j] = -sigma_min;
						}
						else if (eigen_value[j] < -sigma_max) {
							eigen_value[j] = -sigma_max;
						}
					}
				}
				
				for (unsigned int j = 0; j < 3; ++j) {
					for (unsigned int k = 0; k < 3; ++k) {
						rotation_2.data()[3*j+k]= eigen_value[j] * svd.matrixU().data()[3*j+k];
					}
				}
				transform = rotation_2 * svd.matrixV().transpose();


				//singular_value = svd.singularValues();
				//keep_transform = false;

				//double k = (svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().transpose()- transform).norm();
				//if (k > 1e-12) {
				//	std::cout << k << std::endl;
				//}

				//if (!getDigonalForVolumePreserve(singular_value, volume_preserve_max, volume_preserve_min, sigma)) {
				//	transform_for_volume_perserve = Matrix3d::Identity();
				//	
				//}
				//else {
				//	for (unsigned int k = 0; k < 3; ++k) {
				//		for (unsigned int kk = 0; kk < 3; ++kk) {
				//			transform_for_volume_perserve.data()[3 * k + kk] = sigma.data()[k] * svd.matrixV().data()[3 * k + kk];
				//		}
				//	}
				//	transform_for_volume_perserve = (rotation_2 * transform_for_volume_perserve.transpose()).eval();
				//}
				////if (i == 0) {
				////	std::cout << (transform - svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().transpose()).squaredNorm() << std::endl;
				////}
				//transform = rotation_2* svd.matrixV().transpose();

				//if ((transform - Matrix3d::Identity()).squaredNorm() > 1e-6) {
				//	std::cout << transform << std::endl;
				//}

				//for (int j = 0; j < 9; ++j) {
				//	if (abs(transform.data()[j]) < 1e-12) {
				//		transform.data()[j] = 0.0;
				//	}
				//}
				//if (i == 0) {
				//	std::cout << "transform " << std::endl;					
				//	std::cout << transform << std::endl;
				//}
				//p = transform * (PT_pos[i].transpose());
				//temp_p = getARAPmatrix() * p.transpose();

				//std::cout << temp_p.colwise().sum() << std::endl;

				//p_volume_preserve = transform_for_volume_perserve * (PT_pos[i].transpose());
				ATAp_ARAP_volume_preserve[i] =
					((ARAP_stiffness * mesh_struct->volume[i]) * transform * mesh_struct->A[i]).transpose();//+(volume_preserve_stiffness * mesh_struct->volume[i])* transform_for_volume_perserve
				//.transpose()//+ volume_preserve_stiffness * p_volume_preserve
																			   //tet_local_A*tet_local_A=tet_local_A, we omit it here   tet_local_A is symmetric, need not to use transpose
				//mesh_struct->volume[i]*
				temEnergy[thread_id] += 0.5 * ((ARAP_stiffness * mesh_struct->volume[i]) * (transform- deformation_gradient).squaredNorm());//+ volume_preserve_stiffness * (p_volume_preserve - q).squaredNorm()
				//mesh_struct->volume[i] *
			}
		}
	}
	else {

		for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
			vertex_index_start = vertex_begin_per_tetrahedron[j];
			mesh_struct = &(*tetrahedron)[j].mesh_struct;
			index_end = mesh_struct->tetrahedron_index_begin_per_thread[thread_id + 1];
			index_begin = mesh_struct->tetrahedron_index_begin_per_thread[thread_id];
			ATAp_ARAP_volume_preserve = p_ARAP_volume_preserve[j].data();
			ARAP_stiffness = (*tetrahedron)[j].ARAP_stiffness;
			//volume_preserve_stiffness = (*tetrahedron)[j].volume_preserve_stiffness;
			sigma_min = (*tetrahedron)[j].sigma_limit[0];
			sigma_max = (*tetrahedron)[j].sigma_limit[1];

			//P_inv = (*tetrahedron)[j].mesh_struct.P_inv.data();

			for (unsigned int i = index_begin; i < index_end; ++i) {
				vertex_index = mesh_struct->indices[i].data();
				//set matrix q
				for (int k = 0; k < 4; ++k) {
					temp_vertex_index = vertex_index[k] + vertex_index_start;
					for (int l = 0; l < 3; ++l) {
						q.data()[3 * k + l] = u[l].data()[temp_vertex_index];
					}
				}
				//center = 0.25 * q.rowwise().sum();
				//q = (q.colwise() - center).eval();
				//for (unsigned int j = 0; j < 3; ++j)
				//{
				//	q_e.col(j) = q.col(j) - q.col(3);
				//}

				for (unsigned int j = 1; j < 4; ++j)
				{
					q_e.col(j-1) = q.col(j) - q.col(0);
					//p_e.col(j) = temp_P_.col(j) - temp_P_.col(3);
				}

				//transform = q_e * p_e.inverse();
				memcpy(P_inv.data(), mesh_struct->A[i].data() + 3, 72);
				deformation_gradient = q_e * P_inv.transpose();

				//transform = PT[i].transpose() * q_e.transpose();
				//transform = (q * PT_times_PPT_inv[i]);				//

				svd.compute(deformation_gradient, ComputeFullU | ComputeFullV);
				//rotation_2 = svd.matrixU();
				determinant = deformation_gradient.determinant();
				if (determinant < 0) {
					eigen_value[2] = -eigen_value[2];
				}

				for (unsigned int j = 0; j < 3; ++j) {
					if (eigen_value[j] > 0) {
						if (eigen_value[j] < sigma_min) {
							eigen_value[j] = sigma_min;
						}
						else if (eigen_value[j] > sigma_max) {
							eigen_value[j] = sigma_max;
						}
					}
					else {
						if (eigen_value[j] > -sigma_min) {
							eigen_value[j] = -sigma_min;
						}
						else if (eigen_value[j] < -sigma_max) {
							eigen_value[j] = -sigma_max;
						}
					}
				}

				for (unsigned int j = 0; j < 3; ++j) {
					for (unsigned int k = 0; k < 3; ++k) {
						rotation_2.data()[3 * j + k] = eigen_value[j] * svd.matrixU().data()[3 * j + k];
					}
				}
				transform = rotation_2 * svd.matrixV().transpose();

				//singular_value = svd.singularValues();
				//keep_transform = false;
				//if (!getDigonalForVolumePreserve(singular_value, volume_preserve_max, volume_preserve_min, sigma)) {
				//	transform_for_volume_perserve = Matrix3d::Identity();
				//}
				//else {
				//	for (unsigned int k = 0; k < 3; ++k) {
				//		for (unsigned int kk = 0; kk < 3; ++kk) {
				//			transform_for_volume_perserve.data()[3 * k + kk] = sigma.data()[k] * svd.matrixV().data()[3 * k + kk];
				//		}
				//	}
				//	transform_for_volume_perserve = (rotation_2 * transform_for_volume_perserve.transpose()).eval();
				//}
				//transform = rotation_2* svd.matrixV().transpose();
			//	p = transform * (PT_pos[i].transpose());
			//	temp_p = getARAPmatrix() * p.transpose();
				//p_volume_preserve = transform_for_volume_perserve * (PT_pos[i].transpose());
				ATAp_ARAP_volume_preserve[i] = ((ARAP_stiffness * mesh_struct->volume[i]) * transform * mesh_struct->A[i]).transpose();//+ (volume_preserve_stiffness * mesh_struct->volume[i]) * transform_for_volume_perserve
				//tet_local_A*tet_local_A=tet_local_A, we omit it here   tet_local_A is symmetric, need not to use transpose
			}
		}
	}
}

bool ProjectDynamic::getDigonalForVolumePreserve(Vector3d& svd_eigen, double max, double min, Vector3d& sigma)
{
	double sigma_ = svd_eigen[0] * svd_eigen[1] * svd_eigen[2];
	double target_sigma;
	if (sigma_ > max) {
		target_sigma = max;
	}
	else if (sigma_ < min) {
		target_sigma = min;
	}
	else {
		sigma = svd_eigen;
		return false;
	}
	Vector3d D;
	D.setZero();
	Vector3d delta_CD;
	double sigma_multi;
	sigma_multi = sigma_;
	sigma = svd_eigen + D;
	Vector3d pre;
	pre.setOnes();
	int itr_num = 0;
	while ((sigma_multi <min || sigma_multi > max) && itr_num < 50) {
		pre = D;
		delta_CD.data()[0] = sigma.data()[1] * sigma.data()[2];
		delta_CD.data()[1] = sigma.data()[0] * sigma.data()[2];
		delta_CD.data()[2] = sigma.data()[0] * sigma.data()[1];		
		D = ((delta_CD.dot(D) - (sigma_multi-target_sigma)) / delta_CD.squaredNorm()) * delta_CD;		
		sigma = svd_eigen + D;
		sigma_multi = sigma.data()[0] * sigma.data()[1] * sigma.data()[2];
		itr_num++;
	}
	//std::cout << itr_num << std::endl;

	return true;

}


//LOCAL_EDGE_LENGTH_PROJECTION
void ProjectDynamic::localEdgeLengthProjectionPerThread(int thread_id, bool with_energy)
{
	unsigned int vertex_index_start;
	if (with_energy) {
		Vector3d q0, q1, q01;
		double curLen;
		double* edges_length;
		unsigned int* edge_vertices;
		unsigned int* edge_index_begin_per_thread;
		double* length_stiffness;
		for (unsigned int j = 0; j < total_cloth_num; ++j) {
			vertex_index_start = vertex_begin_per_cloth[j];
			edges_length = (*cloth)[j].mesh_struct.edge_length.data();
			edge_index_begin_per_thread = (*cloth)[j].mesh_struct.edge_index_begin_per_thread.data();
			length_stiffness = (*cloth)[j].length_stiffness.data();
			edge_vertices = (*cloth)[j].mesh_struct.edge_vertices.data();
			for (unsigned int i = edge_index_begin_per_thread[thread_id]; i < edge_index_begin_per_thread[thread_id + 1]; ++i) {
				q0.data()[0] = u[0].data()[edge_vertices[i << 1] + vertex_index_start];
				q1.data()[0] = u[0].data()[edge_vertices[(i << 1) + 1] + vertex_index_start];
				q0.data()[1] = u[1].data()[edge_vertices[i << 1] + vertex_index_start];
				q1.data()[1] = u[1].data()[edge_vertices[(i << 1) + 1] + vertex_index_start];
				q0.data()[2] = u[2].data()[edge_vertices[i << 1] + vertex_index_start];
				q1.data()[2] = u[2].data()[edge_vertices[(i << 1) + 1] + vertex_index_start];
				q01 = q0 - q1;
				curLen = sqrt(dotProduct(q01.data(), q01.data()));
				temEnergy[thread_id] += 0.5 * length_stiffness[i] * (curLen - edges_length[i]) * (curLen - edges_length[i]);
				p_edge_length[j][i] = q01 * (length_stiffness[i] / curLen * edges_length[i]);
			}
		}
	}
	else {
		Vector3d q0, q1, q01;
		double curLen;
		double* edges_length;
		unsigned int* edge_index_begin_per_thread;
		double* length_stiffness;
		unsigned int* edge_vertices;
		for (unsigned int j = 0; j < total_cloth_num; ++j) {
			vertex_index_start = vertex_begin_per_cloth[j];
			edges_length = (*cloth)[j].mesh_struct.edge_length.data();
			edge_index_begin_per_thread = (*cloth)[j].mesh_struct.edge_index_begin_per_thread.data();
			length_stiffness = (*cloth)[j].length_stiffness.data();
			edge_vertices = (*cloth)[j].mesh_struct.edge_vertices.data();
			for (unsigned int i = edge_index_begin_per_thread[thread_id]; i < edge_index_begin_per_thread[thread_id + 1]; ++i) {
				q0.data()[0] = u[0].data()[edge_vertices[i << 1] + vertex_index_start];
				q1.data()[0] = u[0].data()[edge_vertices[(i << 1) + 1] + vertex_index_start];
				q0.data()[1] = u[1].data()[edge_vertices[i << 1] + vertex_index_start];
				q1.data()[1] = u[1].data()[edge_vertices[(i << 1) + 1] + vertex_index_start];
				q0.data()[2] = u[2].data()[edge_vertices[i << 1] + vertex_index_start];
				q1.data()[2] = u[2].data()[edge_vertices[(i << 1) + 1] + vertex_index_start];
				q01 = q0 - q1;
				curLen = sqrt(dotProduct(q01.data(), q01.data()));
				p_edge_length[j][i] = q01 * (length_stiffness[i] / curLen * edges_length[i]);
			}
		}
	}
}


void ProjectDynamic::localBendingProjectionPerThread(int thread_id, bool with_energy)
{
	unsigned int size; std::vector<VectorXd>q(3);
	Vector3d aq; double aqnorm; //double a_rest[3];
	VectorXd* lbo;
	std::vector<unsigned int>* vertex_around_vertex_for_bending_;
	double bend_stiffness;
	unsigned int* vertex_index_begin_per_thread;
	Vector3d p;
	VectorXd* p_bend;
	TriangleMeshStruct* mesh_struct;
	int vertex_index_start;
	if (with_energy) {
		for (unsigned int k = 0; k < total_cloth_num; ++k) {
			vertex_index_start = vertex_begin_per_cloth[k];
			if (!(*cloth)[k].mesh_struct.faces.empty()) {
				bend_stiffness = (*cloth)[k].bend_stiffness;
				mesh_struct = &(*cloth)[k].mesh_struct;
				lbo = vertex_lbo[k].data();
				vertex_index_begin_per_thread = (*cloth)[k].mesh_struct.vertex_index_begin_per_thread.data();
				for (unsigned int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id + 1]; ++i) {
					vertex_around_vertex_for_bending_ = &mesh_struct->vertices[i].neighbor_vertex;
					p_bend = p_bending[k][i].data();
					size = vertex_around_vertex_for_bending_->size()+1;
					for (unsigned int j = 0; j < 3; ++j) {
						q[j].resize(size);
					}
					q[0][0] = u[0].data()[i + vertex_index_start];
					q[1][0] = u[1].data()[i + vertex_index_start];
					q[2][0] = u[2].data()[i + vertex_index_start];
					for (unsigned int h = 1; h < size; h++) {
						q[0][h] = u[0].data()[vertex_around_vertex_for_bending_->data()[h-1] + vertex_index_start];
						q[1][h] = u[1].data()[vertex_around_vertex_for_bending_->data()[h-1] + vertex_index_start];
						q[2][h] = u[2].data()[vertex_around_vertex_for_bending_->data()[h-1] + vertex_index_start];
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
		for (unsigned int k = 0; k < total_cloth_num; ++k) {
			if (!(*cloth)[k].mesh_struct.faces.empty()) {
				bend_stiffness = (*cloth)[k].bend_stiffness;
				mesh_struct = &(*cloth)[k].mesh_struct;
				lbo = vertex_lbo[k].data();
				vertex_index_begin_per_thread = (*cloth)[k].mesh_struct.vertex_index_begin_per_thread.data();
				vertex_index_start = vertex_begin_per_cloth[k];
				for (unsigned int i = vertex_index_begin_per_thread[thread_id]; i < vertex_index_begin_per_thread[thread_id + 1]; ++i) {
					vertex_around_vertex_for_bending_ = &mesh_struct->vertices[i].neighbor_vertex;
					p_bend = p_bending[k][i].data();
					size = vertex_around_vertex_for_bending_->size() + 1;
					for (unsigned int j = 0; j < 3; ++j) {
						q[j].resize(size);
					}
					q[0][0] = u[0].data()[i + vertex_index_start];
					q[1][0] = u[1].data()[i + vertex_index_start];
					q[2][0] = u[2].data()[i + vertex_index_start];
					for (unsigned int h = 1; h < size; h++) {
						q[0][h] = u[0].data()[vertex_around_vertex_for_bending_->data()[h - 1] + vertex_index_start];
						q[1][h] = u[1].data()[vertex_around_vertex_for_bending_->data()[h - 1] + vertex_index_start];
						q[2][h] = u[2].data()[vertex_around_vertex_for_bending_->data()[h - 1] + vertex_index_start];
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
	
	if (with_energy) {
		double delta_q[3];
		TriangleMeshStruct* mesh_struct;
		int vertex_index_start;
		for (unsigned int j = 0; j < total_cloth_num; ++j) {
			vertex_index_start = vertex_begin_per_cloth[j];
			if (!(*cloth)[j].mesh_struct.anchor_vertex.empty()) {
				mesh_struct = &(*cloth)[j].mesh_struct;
				for (unsigned int i = mesh_struct->anchor_index_begin_per_thread[thread_id]; i < mesh_struct->anchor_index_begin_per_thread[thread_id + 1]; ++i)
				{
					delta_q[0] = u[0].data()[mesh_struct->anchor_vertex[i] + vertex_index_start] - mesh_struct->anchor_position[i][0];
					delta_q[1] = u[1].data()[mesh_struct->anchor_vertex[i] + vertex_index_start] - mesh_struct->anchor_position[i][1];
					delta_q[2] = u[2].data()[mesh_struct->anchor_vertex[i] + vertex_index_start] - mesh_struct->anchor_position[i][2];
					temEnergy[thread_id] += 0.5 * (*cloth)[j].position_stiffness * dotProduct(delta_q, delta_q);
				}
			}
		}
		TetrahedronMeshStruct* mesh_struct2;
		for (unsigned int j = 0; j < total_tetrahedron_num; ++j) {
			vertex_index_start = vertex_begin_per_tetrahedron[j];
			if (!(*tetrahedron)[j].mesh_struct.anchor_vertex.empty()) {
				mesh_struct2 = &(*tetrahedron)[j].mesh_struct;
				for (unsigned int i = mesh_struct2->anchor_index_begin_per_thread[thread_id]; i < mesh_struct2->anchor_index_begin_per_thread[thread_id + 1]; ++i)
				{
					delta_q[0] = u[0].data()[mesh_struct2->anchor_vertex[i] + vertex_index_start] - mesh_struct2->anchor_position[i][0];
					delta_q[1] = u[1].data()[mesh_struct2->anchor_vertex[i] + vertex_index_start] - mesh_struct2->anchor_position[i][1];
					delta_q[2] = u[2].data()[mesh_struct2->anchor_vertex[i] + vertex_index_start] - mesh_struct2->anchor_position[i][2];
					temEnergy[thread_id] += 0.5 * (*tetrahedron)[j].position_stiffness * dotProduct(delta_q, delta_q);
				}
			}
		}
	}
}



//CONSTRUCT_B
//CONSTRUCT_B_WITHOUT_COLLISION
void ProjectDynamic::constructbPerThead(int thread_id, bool with_collision)
{
	int obj_No;
	int dimension;
	for (int i = 0; i < cloth_dimension_per_thread[thread_id].size(); ++i) {
		obj_No = cloth_dimension_per_thread[thread_id][i] / 3;
		dimension = cloth_dimension_per_thread[thread_id][i] % 3;
		constructbPerThead(b[dimension].data(), (*cloth)[obj_No].mesh_struct,
			p_edge_length[obj_No], obj_No, dimension, vertex_lbo[obj_No], p_bending[obj_No],
			thread_id, collision.obj_target_pos.b_sum[obj_No], collision.obj_target_pos.need_update[obj_No],
			with_collision, vertex_begin_per_cloth[obj_No]);

	}

	for (int i = 0; i < tetrahedron_dimension_per_thread[thread_id].size(); ++i) {
		obj_No = tetrahedron_dimension_per_thread[thread_id][i] / 3;
		dimension = tetrahedron_dimension_per_thread[thread_id][i] % 3;
		constructbTetPerThead(b[dimension].data(), (*tetrahedron)[obj_No].mesh_struct, p_ARAP_volume_preserve[obj_No], obj_No, dimension,
			collision.obj_target_pos.b_sum[obj_No+tetrahedron_begin_obj_index], collision.obj_target_pos.need_update[obj_No + tetrahedron_begin_obj_index],
			with_collision, (*tetrahedron)[obj_No].position_stiffness, vertex_begin_per_tetrahedron[obj_No],
			(*tetrahedron)[obj_No].mesh_struct.indices.data(), (*tetrahedron)[obj_No].mesh_struct.anchor_vertex,
			(*tetrahedron)[obj_No].mesh_struct.anchor_position.data());

	}
}


void ProjectDynamic::constructbTetPerThead(double* b, TetrahedronMeshStruct& mesh_struct,
	std::vector<Matrix<double,4,3>>& p_ARAP_volume_preserve, int tet_No, int dimension,
	std::vector<std::array<double, 3>>& collision_b_sum, bool* collision_b_need_update, bool with_collision,
	double position_stiffness, int vertex_index_start, std::array<int, 4>* indices, std::vector<int>& anchor_vertex,
	std::array<double, 3>* anchor_pos)
{
	unsigned int vertex_num = tetrahedron_sys_size[tet_No];
	memset(b + vertex_index_start, 0, vertex_num<<3);
	//collision
	if (with_collision) {
		for (int i = 0; i < vertex_num; ++i) {
			if (collision_b_need_update[i]) {
				b[i + vertex_index_start] += collision_b_sum[i][dimension];
			}
		}
	}

	//ARAP + volime_preserve
	int tet_size = mesh_struct.indices.size();
	for (int i = 0; i < tet_size; ++i) {
		for (int j = 0; j < 4; ++j) {
			b[indices[i][j] + vertex_index_start] += p_ARAP_volume_preserve[i].data()[4 * dimension + j];
		}		
	}
	//position
	for (int i = 0; i < anchor_vertex.size(); ++i) {
		b[anchor_vertex[i] + vertex_index_start] += position_stiffness * anchor_pos[i][dimension];
	}
}


void ProjectDynamic::constructbPerThead(double* b, TriangleMeshStruct& mesh_struct,
	std::vector<Vector3d>& p_edge_length, int cloth_No, int dimension,
	std::vector<VectorXd>& vertex_lbo, std::vector<std::vector<VectorXd>>& p_bending, int thread_id,
	std::vector<std::array<double, 3>>& collision_b_sum, bool* collision_b_need_update, bool with_collision,
	int vertex_index_start)
{	
	unsigned int vertex_num = cloth_sys_size[cloth_No];
	memset(b + vertex_index_start, 0, vertex_num<<3);
	//collision
	if (with_collision) {
		for (unsigned int i = 0; i < vertex_num; ++i) {
			if (collision_b_need_update[i]) {
				b[i + vertex_index_start] += collision_b_sum[i][dimension];
			}
		}
	}
	//edge length
	for (int i = 0; i < mesh_struct.edges.size(); ++i) {
		b[mesh_struct.edge_vertices[i << 1] + vertex_index_start] += p_edge_length[i].data()[dimension];
		b[mesh_struct.edge_vertices[(i << 1) + 1] + vertex_index_start] -= p_edge_length[i].data()[dimension];
	}
	//for (int i = 0; i < cloth_sys_size[cloth_No]; ++i) {
	//	b.data()[i] += p_edge_length[i].data()[dimension];
	//}
	//bending	
	if (!mesh_struct.faces.empty()) {
		for (unsigned int i = 0; i < cloth_sys_size[cloth_No]; ++i) {
			b[i + vertex_index_start] += p_bending[i][dimension].data()[0];
			for (unsigned int j = 1; j < vertex_lbo[i].size(); ++j) {
				b[mesh_struct.vertices[i].neighbor_vertex[j-1] + vertex_index_start] += p_bending[i][dimension].data()[j];
			}
		}
	}
	//position
	for (int i = 0; i < mesh_struct.anchor_vertex.size(); ++i) {
		b[mesh_struct.anchor_vertex[i] + vertex_index_start] += (*cloth)[cloth_No].position_stiffness * mesh_struct.anchor_position[i][dimension];
	}
	
}

//SOLVE_WITH_COLLISION
//SOLVE_WITHOUT_COLLISION
void ProjectDynamic::solveSystemPerThread(int thread_id, bool with_collision)
{
	for (unsigned int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id+1]; ++i) {
		b[i] += (1.0 / (sub_time_step * sub_time_step)) * (mass.cwiseProduct(u_prediction[i]));
		if (with_collision) {
			switch (itr_solver_method)
			{
			case DIRECT_SOLVE:
				u[i] = global_llt.solve(b[i]);
				break;
			}
		}
		else {
			switch (itr_solver_method)
			{
			case DIRECT_SOLVE:
				u[i] = collision_free_llt.solve(b[i]);
				if (i == 1) {
				//	std::cout << "=====" << std::endl;
				//	std::cout << b[i] << std::endl;
				}

				break;
			}
		}
	}
}



void ProjectDynamic::solveClothSystem2(bool compute_energy)
{
	int itr_num;
	switch (itr_solver_method)
	{
	case JACOBI: {
		//std::cout << "refer" << std::endl;
		//testJacobiMatrix();
		//testJacobi();
		//std::cout << b[0] << std::endl;
			iteration_method.solveByJacobi(u.data(), b.data(), itr_num);
		//std::cout << "Jacobi "<< itr_num << std::endl;
	}
			   break;
	case A_JACOBI: {
			iteration_method.solveByAJacobi_2(u.data(), b.data(), itr_num);
		//std::cout << "super jacobi " << itr_num << std::endl;
	}
				break;
	case CHEBYSHEV_A_JACOBI: {
			iteration_method.solveByChebyshevSemiIterativeAJacobi2(u.data(), b.data(), itr_num);
		//std::cout << "chebyshev super jacobi "<< itr_num << std::endl;
	}
							   break;
	case GAUSS_SEIDEL: {
			iteration_method.solveByGaussSeidel(u.data(), b.data(), itr_num);
		//std::cout << "gauss_seidel "<< itr_num << std::endl;
	}
					 break;
	case GAUSS_SEIDEL_CHEBYSHEV: {
			iteration_method.solveByChebyshevGaussSeidel(u.data(), b.data(), itr_num, 0.8);
		//std::cout << "gauss_seidel_chebysev "<< itr_num << std::endl;
	}
							   break;
	case CHEBYSHEV_JACOBI: {
			iteration_method.solveByChebyshevSemiIterativeJacobi(u.data(), b.data(), itr_num);
		//std::cout <<"chebyshev jacobi "<< itr_num << std::endl;
	}
						 break;
	case PCG: {
			iteration_method.solveByPCG(u.data(), b.data(), itr_num);
		//std::cout <<"PCG "<< itr_num << std::endl;
	}
			break;
	}
	if (compute_energy) {
		thread->assignTask(this, COMPUTE_ENERGY);
	}
	ave_iteration += itr_num;
}

//COMPUTE_ENERGY
void ProjectDynamic::computeEnergyPerThread(int thread_id)
{
	temEnergy[thread_id] = 0.0;
	for (unsigned int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id+1]; ++i) {
		VectorXd p0, p1;
		p0 = u[i] - u_prediction[i];
		p1 = p0.cwiseProduct(mass);
		temEnergy[thread_id] += 0.5 / (sub_time_step * sub_time_step) * p0.dot(p1);
	}
}



void ProjectDynamic::addExternalClothForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int cloth_No)
{
	if (!coe.empty()) {
		unsigned int vertex_index_start = vertex_begin_per_cloth[cloth_No];
		for (unsigned int i = 0; i < coe.size(); ++i) {
			for (unsigned int j = 0; j < 3; ++j) {
				f_ext[j].data()[neighbor_vertex[i]+ vertex_index_start] += coe[i] * neighbor_vertex_force_direction[j];
			}
		}
	}
}

void ProjectDynamic::addExternalTetForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int tet_No)
{
	if (!coe.empty()) {
		unsigned int vertex_index_start = vertex_begin_per_tetrahedron[tet_No];
		for (unsigned int i = 0; i < coe.size(); ++i) {
			for (unsigned int j = 0; j < 3; ++j) {
				f_ext[j].data()[neighbor_vertex[i] + vertex_index_start] += coe[i] * neighbor_vertex_force_direction[j];
			}
		}
	}
}

void ProjectDynamic::resetExternalForce()
{
	f_ext = total_gravity;
}

void ProjectDynamic::mainProcess()
{
	PDsetPosPredict();
	for (unsigned int i = 0; i < total_cloth_num; ++i) {
		(*cloth)[i].mesh_struct.getNormal();
	}
	for (unsigned int i = 0; i < total_tetrahedron_num; ++i) {
		(*tetrahedron)[i].mesh_struct.getNormal();
	}

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
	ave_iteration = 0;
}

void ProjectDynamic::update_ave_iteration_record(double& ave_itr)
{

	ave_iteration /=(double) (max_inner_iteration_num * (outer_iteration_num - 1));	
	ave_itr = ave_iteration;
}

Matrix4d ProjectDynamic::getARAPmatrix()
{
	Matrix4d A=Matrix4d::Ones();
	A *= -0.25;
	A += Matrix4d::Identity();
	return A;
}