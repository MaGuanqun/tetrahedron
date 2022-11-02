#include"second_order.h"
#include"FEM_relate.h"

void SecondOrderConstraint::computeEdgeLengthForce(double* vertex_0, double* vertex_1, double stiffness,
	double* potential_0, double* potential_1, double rest_length)
{
	double C = sqrt(EDGE_LENGTH(vertex_0, vertex_1));
	double grad[3];
	grad[0] = (vertex_0[0] - vertex_1[0]) / C;
	grad[1] = (vertex_0[1] - vertex_1[1]) / C;
	grad[2] = (vertex_0[2] - vertex_1[2]) / C;
	C -= rest_length;
	double coe = C * stiffness;
	MULTI(potential_0, grad, coe);
	coe *= -1.0;
	MULTI(potential_1, grad, coe);
}



void SecondOrderConstraint::computeARAPForce(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	double stiffness, Matrix<double, 3, 4>& A, double volume, Matrix<double,3,4>& force)
{
	Vector3d eigen_value;
	Matrix3d q_e;
	double determinant;
	Vector3d position;
	memcpy(q_e.data(), vertex_position_1, 24);
	memcpy(q_e.data() + 3, vertex_position_2, 24);
	memcpy(q_e.data() + 6, vertex_position_3, 24);
	//first use eigen value to store the position of first vertex
	memcpy(eigen_value.data(), vertex_position_0, 24);
	for (unsigned int i = 0; i < 3; ++i) {
		q_e.col(i) -= eigen_value;
	}
	Matrix3d deformation_gradient;
	Matrix3d P_inv;
	memcpy(P_inv.data(), A.data() + 3, 72);
	deformation_gradient = q_e * P_inv.transpose();
	JacobiSVD<Matrix3d> svd;
	svd.compute(deformation_gradient, ComputeFullU | ComputeFullV);

	eigen_value = svd.singularValues();
	determinant = eigen_value[0] * eigen_value[1] * eigen_value[2];

	Matrix3d V = svd.matrixV();
	if (determinant < 0) {
		V.col(2) *= -1.0;
	}

	P_inv = svd.matrixU() * V.transpose();

	force = (stiffness*volume)* (deformation_gradient - P_inv) * A;

}

void SecondOrderConstraint::solveSingleVertexNewton(std::array<double, 3>* vertex_position, double stiffness, double dt,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn)
{
	unsigned int tet_index;
	Matrix3d Hessian_single;
	Vector3d grad_single;
	unsigned int vertex_no;
	Matrix3d Hessian;
	Vector3d grad;
	Hessian.setZero();
	grad.setZero();
	double C;
	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		tet_index = tet_indices[i];
		vertex_no = findVertexNo(vertex_index, indices[tet_index].data(),4);
		if (getARAPGradHessianNewton(vertex_position[indices[tet_index][0]].data(), vertex_position[indices[tet_index][1]].data(),
			vertex_position[indices[tet_index][2]].data(), vertex_position[indices[tet_index][3]].data(),
			A[tet_index], Hessian_single, grad_single, C, vertex_no)) {
			Hessian += (0.5 * stiffness * volume[tet_index]) * Hessian_single;
			grad -= (0.5 * stiffness * volume[tet_index]) * grad_single;
		}
	}

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



void SecondOrderConstraint::solveNewtonCD_ARAP(std::array<double, 3>* vertex_position, double stiffness, double dt,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn)
{
	unsigned int tet_index;
	Matrix3d Hessian_single;
	Vector3d grad_single;
	unsigned int vertex_no;
	Matrix3d Hessian;
	Vector3d grad;
	Hessian.setZero();
	grad.setZero();
	double C;
	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		tet_index = tet_indices[i];
		vertex_no = findVertexNo(vertex_index, indices[tet_index].data(),4);
		if (getARAPGradHessianNewton(vertex_position[indices[tet_index][0]].data(), vertex_position[indices[tet_index][1]].data(),
			vertex_position[indices[tet_index][2]].data(), vertex_position[indices[tet_index][3]].data(),
			A[tet_index], Hessian_single, grad_single, C, vertex_no)) {
			Hessian += (0.5*stiffness*volume[tet_index])* Hessian_single;
			grad-= (0.5*stiffness * volume[tet_index]) * grad_single;
		}
	}

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

void SecondOrderConstraint::solveSingleVertexCD_ARAP(std::array<double, 3>* vertex_position, double stiffness, double dt,
	Matrix<double, 3, 4>* A, double* lambda, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn)
{
	unsigned int tet_index;
	Matrix3d Hessian_single;
	Vector3d grad_single;
	unsigned int vertex_no;
	MatrixXd Hessian;
	Hessian.resize(3 + tet_indices.size(), 3 + tet_indices.size());
	Hessian.setZero();

	VectorXd b;
	b.resize(3 + tet_indices.size());
	b.setZero();
	double C;
	double alpha;

	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		tet_index = tet_indices[i];
		vertex_no = findVertexNo(vertex_index, indices[tet_index].data(),4);
		if (getARAPGradHessian(vertex_position[indices[tet_index][0]].data(), vertex_position[indices[tet_index][1]].data(),
			vertex_position[indices[tet_index][2]].data(), vertex_position[indices[tet_index][3]].data(),
			A[tet_index], Hessian_single, grad_single, C, vertex_no)) {
			Hessian_single *= -lambda[tet_index];
			Hessian.block<3, 3>(0, 0) += Hessian_single;
			Hessian.block<3, 1>(0, 3 + i) = -grad_single;
			Hessian.block<1, 3>(3 + i, 0) = -grad_single;
			b.head(3) += lambda[tet_index] * grad_single;
		}
		alpha = 1.0 / (stiffness * volume[tet_index] * dt * dt);
		Hessian(3 + i, 3 + i) = -alpha;
		b.data()[3 + i] = C + alpha * lambda[tet_index];
	}

	Hessian.data()[0] += mass[vertex_index];
	Hessian(1, 1) += mass[vertex_index];
	Hessian(2, 2) += mass[vertex_index];
	b.data()[0] -= mass[vertex_index] * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	b.data()[1] -= mass[vertex_index] * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	b.data()[2] -= mass[vertex_index] * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(b);

	SUM_(vertex_position[vertex_index], result);
	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		lambda[tet_indices[i]] += result[3 + i];
	}
}

void SecondOrderConstraint::solveCD_ARAP_block(MatrixXd& Hessian, VectorXd& grad, std::array<double, 3>* vertex_position, double stiffness,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices, std::array<int, 4>* indices,
	double* volume, unsigned int tet_index,  unsigned int* common_vertex_in_order,
	int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num)
{

	Vector3d eigen_value;
	Matrix3d deformation_gradient;
	Matrix3d S, rotation;
	Hessian.setZero();
	grad.setZero();

	FEM::getDeformationGradient(vertex_position[tet_vertex_index[0]].data(), vertex_position[tet_vertex_index[1]].data(),
		vertex_position[tet_vertex_index[2]].data(), vertex_position[tet_vertex_index[3]].data(), A[tet_index], 
		deformation_gradient);
	FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);

	if((deformation_gradient - rotation).squaredNorm()>1e-16){
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (2.0 * stiffness * volume[tet_index]) * (deformation_gradient - rotation) * A[tet_index];
		for (int i = 0; i < unfixed_vertex_num; ++i) {
			memcpy(grad.data() + 3 * i, grad_C_transpose.data() + 3 * unfixed_tet_vertex_index[i], 24);
		}
		Matrix3d Dm = A[tet_index].block<3, 3>(0, 1).transpose();
		FEM::getHessianForSeveralVertex(Hessian, S, rotation, Dm, A[tet_index], unfixed_tet_vertex_index,
			unfixed_vertex_num);
		Hessian *= volume[tet_index] * stiffness;
	}
	for (auto i = neighbor_tet_indices.begin(); i < neighbor_tet_indices.end(); ++i) {
		solveCertainHessianForNeighborTet(vertex_position, stiffness, A[*i], common_vertex_in_order, indices[*i].data(),
			Hessian, volume[*i], grad);
	}
}



bool SecondOrderConstraint::solveCertainHessianForNeighborTet(std::array<double, 3>* vertex_position, double stiffness,
	Matrix<double, 3, 4>& A, unsigned int* &common_vertex_in_order, int* neighbor_tet_vetex_indices, 
	MatrixXd& sys_matrix, double volume, VectorXd& grad)
{// for a single neighbor tet, only compute common vertices
	if (*common_vertex_in_order == 0) {
		common_vertex_in_order++;
		return false;
	}
	int common_vertex_num = *common_vertex_in_order;
	common_vertex_in_order++;


	Matrix<double, 12, 3> Hessian_vertex;
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d S, rotation;


	FEM::getDeformationGradient(vertex_position[neighbor_tet_vetex_indices[0]].data(), vertex_position[neighbor_tet_vetex_indices[1]].data(),
		vertex_position[neighbor_tet_vetex_indices[2]].data(), vertex_position[neighbor_tet_vetex_indices[3]].data(), A,
		deformation_gradient);
	FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);

	if ((deformation_gradient - rotation).squaredNorm() < 1e-16) {
		common_vertex_in_order += (common_vertex_num * 2);
		return false;
	}
	Matrix<double, 3, 4> grad_C_transpose;
	grad_C_transpose = (2.0 * stiffness * volume) * (deformation_gradient - rotation) * A;

	for (int i = 0; i < common_vertex_num; ++i) {
		grad.segment(3 * (*(common_vertex_in_order + i + common_vertex_num)), 3)
			+= grad_C_transpose.col(*(common_vertex_in_order + i));
	}

	Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
	for (int i = 0; i < common_vertex_num; ++i) {
		FEM::getHessianForOneVertex(Hessian_vertex, S, rotation, Dm, A, *(common_vertex_in_order + i));
		Hessian_vertex *= volume * stiffness;
		for (int j = 0; j < common_vertex_num; ++j) {
			sys_matrix.block<3, 3>(3 * (*(common_vertex_in_order + j + common_vertex_num)),
				3 * (*(common_vertex_in_order + i + common_vertex_num))) +=
				Hessian_vertex.block<3, 3>(3 * (*(common_vertex_in_order + j)), 0);
		}
	}

	common_vertex_in_order += (common_vertex_num << 1);
	return true;
}

void SecondOrderConstraint::solveCD_ARAP(std::array<double, 3>* vertex_position, double stiffness, double dt,
	Matrix<double, 3, 4>* A, double* lambda, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass, 
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn)
{
	unsigned int tet_index;
	Matrix3d Hessian_single;
	Vector3d grad_single;
	unsigned int vertex_no;
	MatrixXd Hessian;
	Hessian.resize(3 + tet_indices.size(), 3 + tet_indices.size());
	Hessian.setZero();

	VectorXd b;
	b.resize(3 + tet_indices.size());
	b.setZero();
	double C;
	double alpha;

	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		tet_index = tet_indices[i];
		vertex_no = findVertexNo(vertex_index, indices[tet_index].data(),4);
		if (getARAPGradHessian(vertex_position[indices[tet_index][0]].data(), vertex_position[indices[tet_index][1]].data(),
			vertex_position[indices[tet_index][2]].data(), vertex_position[indices[tet_index][3]].data(),
			A[tet_index],  Hessian_single, grad_single, C, vertex_no)) {
			Hessian_single *= -lambda[tet_index];
			Hessian.block<3, 3>(0, 0) += Hessian_single;
			Hessian.block<3, 1>(0, 3 + i) = -grad_single;
			Hessian.block<1, 3>(3 + i, 0) = -grad_single;
			b.head(3) += lambda[tet_index] * grad_single;
		}		
		alpha = 1.0 / (stiffness * volume[tet_index] * dt * dt);
		Hessian(3 + i, 3 + i) = -alpha;
		b.data()[3 + i] = C + alpha * lambda[tet_index];
	}
	
	Hessian.data()[0] += mass[vertex_index];
	Hessian(1,1) += mass[vertex_index];
	Hessian(2,2) += mass[vertex_index];
	b.data()[0] -= mass[vertex_index] * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	b.data()[1] -= mass[vertex_index] * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	b.data()[2] -= mass[vertex_index] * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(b);

	SUM_(vertex_position[vertex_index], result);
	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		lambda[tet_indices[i]] += result[3 + i];
	}
}


bool SecondOrderConstraint::getCollisionPairHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	double ori_volume, Matrix3d& Hessian, Vector3d& grad, unsigned int vertex_no)
{
	double C;

	double volume =getTetCubeVolume(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3);
	if (abs(volume) >= ori_volume) {
		return false;
	}
	if (abs(volume) < 1e-9) {
		std::cout << "may be edge edge parallel " << std::endl;
		return false;
	}
	double e0[3], e1[3];
	switch (vertex_no)
	{
	case 0: {
		SUB(e0, vertex_position_2, vertex_position_1);
		CROSS(grad.data(), vertex_position_3,  e0);
		double temp[3];
		CROSS(temp, vertex_position_2, vertex_position_1);
		SUM_(grad, temp);
	}
		break;
	case 1:
		SUB(e0, vertex_position_2, vertex_position_0);
		SUB(e1, vertex_position_3, vertex_position_0);
		CROSS(grad.data(), e0, e1);
		break;
	case 2:
		SUB(e0, vertex_position_3, vertex_position_0);
		SUB(e1, vertex_position_1, vertex_position_0);
		CROSS(grad.data(), e0, e1);
		break;	
	case 3:
		SUB(e0, vertex_position_1, vertex_position_0);
		SUB(e1, vertex_position_2, vertex_position_0);
		CROSS(grad.data(), e0, e1);
		break;
	}
	if (volume < 0) {
		grad *= -1.0;
		volume = -volume;
		//std::cout << "collision volume negative" << std::endl;
	}
	double ln_ = log(volume / ori_volume);
	Hessian = ( - 2.0 * ln_ + (ori_volume - volume) * (ori_volume + 3.0 * volume) / (volume * volume))*grad*grad.transpose();
	grad *= (ori_volume - volume) * (2.0 * ln_ - ori_volume / volume + 1.0);
	
	return true;
}




bool SecondOrderConstraint::getARAPGradHessianNewton(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	Matrix<double, 3, 4>& A, Matrix3d& Hessian, Vector3d& grad, double& C, unsigned int vertex_no)
{
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d S, rotation;
	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);
	FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);
	C = (deformation_gradient - rotation).squaredNorm();
	if (C < 1e-16) {
		C = 0;
		return false;
	}

	grad =2.0 * (deformation_gradient - rotation) * A.col(vertex_no);
	Matrix3d Dm = A.block<3, 3>(0, 1).transpose();

	//Matrix<double, 12, 12> Hessian_global;
	//FEM::getHessian(Hessian_global, S, rotation, Dm, A);
	//Hessian = Hessian_global.block<3, 3>(3 * vertex_no, 3 * vertex_no);

	FEM::getHessianForOneVertex(Hessian, S, rotation, Dm, A, vertex_no);
	return true;
}


bool SecondOrderConstraint::getARAPGradHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	Matrix<double, 3, 4>& A, Matrix3d& Hessian, Vector3d& grad, double& C, unsigned int vertex_no)
{
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d S, rotation;
	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);
	FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);
	C = (deformation_gradient - rotation).norm();
	if (C < 1e-8) {
		C = 0;
		return false;
	}
	grad = (1.0 / C) * (deformation_gradient - rotation) * A.col(vertex_no);

	Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
	FEM::getHessianForOneVertex(Hessian, S, rotation, Dm, A,vertex_no);
	Hessian *= (0.5 / C);
	Hessian -= ((1.0 / C) * grad) * grad.transpose();
	return true;
}


void SecondOrderConstraint::solveARAPConstraint(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	double stiffness, double dt,
	Matrix<double, 3, 4>& A, double* inv_mass, double& lambda, const double damping_stiffness, double* mass,
	double volume)
{
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d S, rotation;

	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);

	FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);

	//use P_inv to record transform

	double C = (deformation_gradient - rotation).norm();

	if (C < 1e-8) {
		return;
	}

	double alpha_ = 1.0 / (stiffness * volume * dt * dt);

	Matrix<double, 3, 4> grad_C_transpose;
	grad_C_transpose = (1.0 / C) * (deformation_gradient - rotation) * A;	
	Matrix<double, 12, 1> grad;
	memcpy(grad.data(), grad_C_transpose.data(), 96);

	Matrix<double, 12, 12> Hessian;

	Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
	FEM::getHessian(Hessian, S, rotation, Dm, A);

	Hessian *= (0.5 / C);
	Hessian -= ((1.0 / C) * grad) * grad.transpose();

	for (unsigned int i = 0; i < 4; ++i) {
		if (inv_mass[i] == 0.0) {
			Hessian.block<12, 3>(0, 3 * i).setZero();
			Hessian.block<3, 12>(3 * i, 0).setZero();
			grad.segment(3 * i, 3).setZero();
		}
	}


	//Matrix<double, 12, 1> inv_mass_sys;
	//inv_mass_sys << inv_mass[0], inv_mass[0], inv_mass[0], inv_mass[1], inv_mass[1], inv_mass[1], inv_mass[2], inv_mass[2], inv_mass[2],
	//	inv_mass[3], inv_mass[3], inv_mass[3];
	////how to improve this
	//for (unsigned int j = 0; j < 12; ++j) {
	//	for (unsigned int i = 0; i < 4; ++i) {
	//		Hessian(i + i + i, j) *= inv_mass[i];
	//		Hessian(i + i + i+1, j) *= inv_mass[i];
	//		Hessian(i + i + i+2, j) *= inv_mass[i];
	//	}
	//}

	Hessian *= -lambda;

	for (unsigned int i = 0; i < 4; ++i) {
		Hessian.data()[39 * i] += mass[i];
		Hessian.data()[39 * i+13] += mass[i];
		Hessian.data()[39 * i+26] += mass[i];
	}

	//Hessian.setZero();
	//for (unsigned int i = 0; i < 144; i += 13) {
	//	Hessian.data()[i] += 1.0;
	//}

	//Matrix<double, 12, 1> M_inv_g;
	//SUB(M_inv_g.data(), vertex_position_0, sn_0);
	//M_inv_g.data()[3] = vertex_position_1[0] - sn_1[0];
	//M_inv_g.data()[4] = vertex_position_1[1] - sn_1[1];
	//M_inv_g.data()[5] = vertex_position_1[2] - sn_1[2];
	//M_inv_g.data()[6] = vertex_position_2[0] - sn_2[0];
	//M_inv_g.data()[7] = vertex_position_2[1] - sn_2[1];
	//M_inv_g.data()[8] = vertex_position_2[2] - sn_2[2];
	//M_inv_g.data()[9] = vertex_position_3[0] - sn_3[0];
	//M_inv_g.data()[10] = vertex_position_3[1] - sn_3[1];
	//M_inv_g.data()[11] = vertex_position_3[2] - sn_3[2];
	//M_inv_g -= lambda * inv_mass_sys.cwiseProduct(grad);

	ColPivHouseholderQR <Matrix<double,12,12>> linear(Hessian);

	//if (abs(Hessian.determinant()) < 1e-8) {
	//	std::cout << Hessian.determinant() << std::endl;
	//}

	double h = C + alpha_ * lambda;

	double delta_lambda = (-h) / (grad.dot(linear.solve(grad)) + alpha_);//  + grad.dot(linear.solve(M_inv_g))
	Matrix<double, 12, 1> delta_x = linear.solve(delta_lambda * grad); //-M_inv_g

	//double delta_lambda = (-h) / (grad.dot(linear.solve(inv_mass_sys.cwiseProduct(grad))) + alpha_);//  + grad.dot(linear.solve(M_inv_g))
	//Matrix<double, 12, 1> delta_x = linear.solve(delta_lambda * (inv_mass_sys.cwiseProduct(grad))); //-M_inv_g


	vertex_position_0[0] += delta_x[0];
	vertex_position_0[1] += delta_x[1];
	vertex_position_0[2] += delta_x[2];

	vertex_position_1[0] += delta_x[3];
	vertex_position_1[1] += delta_x[4];
	vertex_position_1[2] += delta_x[5];
	
	vertex_position_2[0] += delta_x[6];
	vertex_position_2[1] += delta_x[7];
	vertex_position_2[2] += delta_x[8];

	vertex_position_3[0] += delta_x[9];
	vertex_position_3[1] += delta_x[10];
	vertex_position_3[2] += delta_x[11];
	 
	lambda += delta_lambda;

}





void SecondOrderConstraint::solveEdgeLengthConstraint(Vector3f& p1, Vector3f& p2, const double d, double mass_0,
	double mass_1, Vector3f& ori_p1, Vector3f& ori_p2, bool v0_fixed, bool v1_fixed, double& lambda)
{
	float w1 =1.0/mass_0;
	float w2 = 1.0/mass_1;

	if (v0_fixed) {
		w1 = 0;
	}
	if (v1_fixed) {
		w2 = 0;
	}
	Matrix<float, 6,1> g;
	g.setZero();
	if (w1 > 0)
		g.segment(0, 3) = (p1 - ori_p1) / w1;
	if (w2 > 0)
		g.segment(3, 3) = (p2 - ori_p2) / w2;

	//std::cout <<"==="<< g.transpose() << std::endl;

	Matrix<float, 6, 1> M_1;
	M_1 << w1, w1, w1, w2, w2, w2;
	Eigen::Vector3f n = p1 - p2;
	float norm = n.norm();
	n /= norm;
	Matrix<float, 6,1> grad;
	grad << n, -n;
	Eigen::Matrix<float, 6, 6> H;
	Eigen::Matrix3f s = Eigen::Matrix3f::Identity() - n * n.transpose();
	H.topLeftCorner(3, 3) = s;
	H.topRightCorner(3, 3) = -s;
	H.bottomLeftCorner(3, 3) = -s;
	H.bottomRightCorner(3, 3) = s;
	H /= norm;

	Eigen::Matrix<float, 6, 6> M_lambda_H_1 = (Eigen::Matrix<float, 6, 6>::Identity() + lambda * M_1.asDiagonal() * H).inverse() * M_1.asDiagonal();
	Eigen::Matrix<float, 6,1> tmp = M_lambda_H_1 * grad;
	g += lambda * grad;
	// Eigen::Vector<float, 6> tmp = M_1 * grad;
	float d_lambda = (norm - d - g.dot(tmp)) / (grad.dot(tmp));
	lambda += d_lambda;
	// std::cout << lambda << std::endl;
	Eigen::Matrix<float, 6,1> d_x = d_lambda * tmp + M_lambda_H_1 * g;
	p1 -= d_x.head(3);
	p2 -= d_x.tail(3);
}

void SecondOrderConstraint::solveEdgeLengthConstraint(double* p0, double* p1, const double rest_length, double stiffness, double mass_0,
	double mass_1, double time_step, double* sn_0, double* sn_1, bool v0_fixed, bool v1_fixed, unsigned int edge_index, double& lambda, unsigned int vertex_0_index, unsigned int vertex_1_index)
{
	if (v0_fixed && v1_fixed) {
		return;
	}

	//std::cout << "true " << std::endl;
	//std::cout << edge_index << " " << mass_0 << " " << mass_1 << " " << rest_length<<" "<< sn_0[0] << " " << sn_0[1] << " " << sn_0[2] << std::endl;


	//Vector3f p_1, p_2,ori_1,ori_2;
	//for (unsigned int i = 0; i < 3; ++i) {
	//	p_1[i] = p0[i];
	//	p_2[i] = p1[i];
	//	ori_1[i] = sn_0[i];
	//	ori_2[i] = sn_1[i];
	//}

	//if (edge_index == 3) {
	//	std::cout << lambda << std::endl;

	//}

	//solveEdgeLengthConstraint(p_1, p_2, rest_length, mass_0, mass_1, ori_1, ori_2, v0_fixed, v1_fixed, lambda);

	//if (edge_index == 3) {
	//	std::cout << lambda << std::endl;
	//	std::cout << sn_0[0] << " " << sn_0[1] << " " << sn_0[2] << std::endl;
	//	std::cout << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
	//}


	//p0[0] = p_1[0];
	//p0[1] = p_1[1];
	//p0[2] = p_1[2];
	//p1[0] = p_2[0];
	//p1[1] = p_2[1];
	//p1[2] = p_2[2];
	Vector3d n;
	SUB(n.data(), p0, p1);
	double n_norm = sqrt(DOT(n, n));
	if (abs(n_norm - rest_length) < 1e-8) {
		return;
	}
	Vector3d grad_ = n / n_norm;
	double alpha =  1.0 / (time_step * time_step * stiffness);
	Matrix3d He = (1.0 / n_norm) * (Matrix3d::Identity() - grad_ * grad_.transpose());
	//Matrix<double, 6, 6>Hessian;
	//Hessian.block<3, 3>(0, 0) = He;
	//Hessian.block<3, 3>(3, 3) = He;
	//Hessian.block<3, 3>(3, 0) = -He;
	//Hessian.block<3, 3>(0, 3) = -He;
	VectorXd mass_inv(6);
	mass_inv[0] = mass_inv[1] = mass_inv[2] = 1.0 / mass_0;
	mass_inv[3] = mass_inv[4] = mass_inv[5] = 1.0 / mass_1;
	if (v0_fixed) {
		mass_inv[0] = mass_inv[1] = mass_inv[2] = 0;
	}
	if (v1_fixed) {
		mass_inv[3] = mass_inv[4] = mass_inv[5] = 0;
	}
	//Matrix<double, 6, 1> gradient;
	//gradient.segment(0, 3) = grad_;
	//gradient.segment(3, 3) = -grad_;
	/*
	Matrix<double, 6, 6>sys_matrix = Matrix<double, 6, 6>::Identity() + lambda * mass_inv.asDiagonal() * Hessian;

	//ColPivHouseholderQR <Matrix<double, 6, 6>> linear(sys_matrix);
	Matrix<double, 6, 6>  linear = sys_matrix.inverse() * mass_inv.asDiagonal();
	Matrix<double, 6, 1>g;
	g.data()[0] = mass_0 * (p0[0] - sn_0[0]);
	g.data()[1] = mass_0 * (p0[1] - sn_0[1]);
	g.data()[2] = mass_0 * (p0[2] - sn_0[2]);
	g.data()[3] = mass_1 * (p1[0] - sn_1[0]);
	g.data()[4] = mass_1 * (p1[1] - sn_1[1]);
	g.data()[5] = mass_1 * (p1[2] - sn_1[2]);
	g += lambda * gradient;

	Matrix<double, 6,1> tmp = linear * gradient;

	//double coe = linear.solve(mass_inv.asDiagonal()*gradient).dot(gradient)+alpha;
	double coe = tmp.dot(gradient)+alpha;
	double delta_lambda = (n_norm - rest_length - alpha * lambda) - (g.dot(tmp));
	delta_lambda /= coe;
	Matrix<double, 6, 1>delta_x;
	delta_x =-1.0*(delta_lambda*tmp + sys_matrix*g);
	//delta_x =-1.0*( delta_lambda *linear *  (mass_inv.asDiagonal() * gradient) + linear *  M_inv_g);
	p0[0] += delta_x[0];
	p0[1] += delta_x[1];
	p0[2] += delta_x[2];

	p1[0] += delta_x[3];
	p1[1] += delta_x[4];
	p1[2] += delta_x[5];		
	lambda += delta_lambda;
	*/
	//Matrix3d He = (((- lambda) /n_norm) * grad_) * grad_.transpose();
	//double coe = lambda / n_norm;
	//He.data()[0] += coe;
	//He.data()[4] += coe;
	//He.data()[8] += coe;
	double coe;
		if (v0_fixed) {
			/*
			Matrix<double, 3, 3> sys_matrix;
			sys_matrix = mass_1/(time_step*time_step) * Matrix3d::Identity() + stiffness * grad_*grad_.transpose() + (stiffness*(n_norm - rest_length)) * He;
			Vector3d g;
			//g.data()[0] = mass_1 * (p1[0] - sn_1[0]);
			//g.data()[1] = mass_1 * (p1[1] - sn_1[1]);
			//g.data()[2] = mass_1 * (p1[2] - sn_1[2]);
			g.setZero();
			g -= (stiffness*(n_norm - rest_length))*grad_;
			ColPivHouseholderQR <Matrix3d> linear(sys_matrix);
			Vector3d delta_x = -linear.solve(g);
			p1[0] += delta_x[0];
			p1[1] += delta_x[1];
			p1[2] += delta_x[2];
			*/
			
			Matrix<double, 3, 3> sys_matrix;
			sys_matrix = mass_1 * Matrix3d::Identity()-lambda * He;
			//sys_matrix.setZero();
			//sys_matrix -= He;
			//sys_matrix.data()[0] += mass_1;
			//sys_matrix.data()[4] += mass_1;
			//sys_matrix.data()[8] += mass_1;
			Vector3d g;
			g.setZero();
			//g.data()[0] = mass_1 * (p1[0] - sn_1[0]);
			//g.data()[1] = mass_1 * (p1[1] - sn_1[1]);
			//g.data()[2] = mass_1 * (p1[2] - sn_1[2]);
			g += lambda * grad_;
			double h = n_norm - rest_length + alpha * lambda;
			ColPivHouseholderQR <Matrix3d> linear(sys_matrix);
			//LDLT <Matrix3d> linear(sys_matrix);
			coe = (linear.solve(grad_)).dot(grad_) + alpha;
			double delta_lambda =( - h - (linear.solve(g)).dot(grad_))/coe;//
			Vector3d delta_x = linear.solve((- delta_lambda) * grad_-g);//-g
			p1[0] += delta_x[0];
			p1[1] += delta_x[1];
			p1[2] += delta_x[2];
			lambda += delta_lambda;
			
			return;
		}		
		if (v1_fixed) {
			/*
			Matrix<double, 3, 3> sys_matrix;
			sys_matrix = mass_0 / (time_step * time_step) * Matrix3d::Identity() + stiffness * grad_ * grad_.transpose() + (stiffness * (n_norm - rest_length)) * He;
			Vector3d g;
			//g.data()[0] = mass_0 * (p0[0] - sn_0[0]);
			//g.data()[1] = mass_0 * (p0[1] - sn_0[1]);
			//g.data()[2] = mass_0 * (p0[2] - sn_0[2]);
			g.setZero();
			g += (stiffness * (n_norm - rest_length)) * grad_;
			ColPivHouseholderQR <Matrix3d> linear(sys_matrix);
			Vector3d delta_x = -linear.solve(g);
			p0[0] += delta_x[0];
			p0[1] += delta_x[1];
			p0[2] += delta_x[2];
			*/
			
			Matrix<double, 3, 3> sys_matrix;
			sys_matrix = mass_0 * Matrix3d::Identity() -lambda * He;
			//sys_matrix.setZero();
			//sys_matrix -= He;
			//sys_matrix.data()[0] += mass_0;
			//sys_matrix.data()[4] += mass_0;
			//sys_matrix.data()[8] += mass_0;
			Vector3d g;
			g.setZero();
			//g.data()[0] = mass_0 * (p0[0] - sn_0[0]);
			//g.data()[1] = mass_0 * (p0[1] - sn_0[1]);
			//g.data()[2] = mass_0 * (p0[2] - sn_0[2]);
			g -= lambda * grad_;

			double h = n_norm - rest_length + alpha * lambda;
			//LDLT <Matrix3d> linear(sys_matrix);
			ColPivHouseholderQR <Matrix3d> linear(sys_matrix);
			coe = linear.solve(grad_).dot(grad_) + alpha;
			double delta_lambda = (-h + linear.solve(g).dot(grad_)) / coe;//
			Vector3d delta_x = linear.solve(delta_lambda * grad_-g);//- g
			p0[0] += delta_x[0];
			p0[1] += delta_x[1];
			p0[2] += delta_x[2];
			lambda += delta_lambda;
			
			return;
		}

		/*
		///test
		Matrix<double, 6, 1> gradient;
		gradient.segment(0, 3) = grad_;
		gradient.segment(3, 3) = -grad_;
		Matrix<double, 6, 6> sys_matrix;
		sys_matrix.setZero();
		sys_matrix.block<3, 3>(0, 0) = -He;
		sys_matrix.block<3, 3>(3, 3) = -He;
		sys_matrix.block<3, 3>(3, 0) = He;
		sys_matrix.block<3, 3>(0, 3) = He;
		sys_matrix *= (n_norm - rest_length) * stiffness;
		sys_matrix += stiffness * gradient * gradient.transpose();
		for (unsigned int i = 0; i < 18; i += 7) {
			sys_matrix.data()[i] += mass_0 / (time_step * time_step);
		}
		for (unsigned int i = 21; i < 36; i += 7) {
			sys_matrix.data()[i] += mass_1 / (time_step * time_step);
		}
		Matrix<double, 6, 1>g;
		g.setZero();
		//g.data()[0] = mass_0 / (time_step * time_step) * (p0[0] - sn_0[0]);
		//g.data()[1] = mass_0 / (time_step * time_step) * (p0[1] - sn_0[1]);
		//g.data()[2] = mass_0 / (time_step * time_step) * (p0[2] - sn_0[2]);
		//g.data()[3] = mass_1 / (time_step * time_step) * (p1[0] - sn_1[0]);
		//g.data()[4] = mass_1 / (time_step * time_step) * (p1[1] - sn_1[1]);
		//g.data()[5] = mass_1 / (time_step * time_step) * (p1[2] - sn_1[2]);
		g += (n_norm - rest_length) * stiffness * gradient;
		Matrix<double, 6, 1> result;
		ColPivHouseholderQR <Matrix<double, 6, 6>> linear_(sys_matrix);
		result = -1.0 * linear_.solve(g);
		p0[0] += result[0];
		p0[1] += result[1];
		p0[2] += result[2];
		p1[0] += result[3];
		p1[1] += result[4];
		p1[2] += result[5];

		////
		*/
		
		Matrix<double, 6, 1> gradient;
		gradient.segment(0, 3) = grad_;
		gradient.segment(3, 3) = -grad_;
		Matrix<double, 6, 6> sys_matrix;
		sys_matrix.setZero();
		sys_matrix.block<3, 3>(0, 0) =- He;
		sys_matrix.block<3, 3>(3, 3) =- He;
		sys_matrix.block<3, 3>(3, 0) = He;
		sys_matrix.block<3, 3>(0, 3) = He;
		sys_matrix *= lambda;
		for (unsigned int i = 0; i < 18; i+=7) {
			sys_matrix.data()[i] += mass_0;
		}
		for (unsigned int i = 21; i < 36; i+=7) {
			sys_matrix.data()[i] += mass_1;
		}
		//double det = sys_matrix.determinant();
		//if (det < 1e-8) {
		//	std::cout << det << std::endl;
		//}
		//ColPivHouseholderQR <Matrix<double, 6, 6>> linear(sys_matrix);
		//LDLT <Matrix<double, 6, 6>> linear(sys_matrix);
		double h = n_norm - rest_length + alpha*lambda;
		Matrix<double, 6, 1>g;
		g.setZero();
		//g.data()[0] =mass_0 *(p0[0] - sn_0[0]);
		//g.data()[1] =mass_0 *(p0[1] - sn_0[1]);
		//g.data()[2] =mass_0 *(p0[2] - sn_0[2]);
		//g.data()[3] = mass_1 * (p1[0] - sn_1[0]);
		//g.data()[4] = mass_1 * (p1[1] - sn_1[1]);
		//g.data()[5] = mass_1 * (p1[2] - sn_1[2]);
		g -= lambda*gradient;
	//	coe = linear.solve(gradient).dot(gradient) + alpha;
		double delta_lambda;
	//	delta_lambda = (-h + linear.solve(g).dot(gradient)) / coe;//+
		//if (edge_index == 0) {
		////	std::cout << g << std::endl;
		//std::cout << linear.solve(g) << std::endl;
		//}
		Matrix<double, 6, 1> delta_x; 
		// delta_x = linear.solve(delta_lambda * gradient - g);//		-g
		//if (vertex_0_index == 2451 || vertex_1_index == 2451) {
		//	std::cout << delta_lambda<<" "<<edge_index << std::endl;
		//	std::cout << delta_x.segment(0, 3).norm() << " " << delta_x.segment(3, 3).norm() << std::endl;
		//}
		Matrix<double, 7, 7> check;
		check.block<6, 6>(0, 0) = sys_matrix;
		check.block<6, 1>(0, 6) = -gradient;
		check.block<1, 6>(6, 0) = gradient.transpose();
		check.data()[7 * 7 - 1] = alpha;
		Matrix<double, 7, 1> vector;
		//vector.setZero();
		memcpy(vector.data(), g.data(), 48);
		vector.data()[6] = h;		
		Matrix<double, 7, 1> result;
		ColPivHouseholderQR <Matrix<double, 7, 7>> linear_(check);
	   result =-1.0* linear_.solve(vector);

	   //double k = check.determinant();
	  // if (k < 1e-7) {
		//   std::cout << edge_index << " " << k << std::endl;
	   //}

		//memcpy(result.data(), delta_x.data(), 48);
		//result.data()[6] = delta_lambda;
		//double error = (check * result + vector).norm();
		//if (error > 1e-10) {
		//	std::cout << "error " << error << std::endl;
		//	
			//std::cout << check << std::endl;
			//std::cout << lambda << std::endl;
			//std::cout << check.determinant() << std::endl;;
			//std::cout << g.transpose() << std::endl;;
			//std::cout << vector.transpose() << std::endl;;
			//std::cout <<  std::endl;
			//std::cout << result.transpose() << std::endl;;
		//}
		//}
		//else {
		//	//std::cout << error << std::endl;
		//}

	  


		p0[0] += result[0];
		p0[1] += result[1];
		p0[2] += result[2];
		p1[0] += result[3];
		p1[1] += result[4];
		p1[2] += result[5];
		lambda += result[6];
	
		//p0[0] += delta_x[0];
		//p0[1] += delta_x[1];
		//p0[2] += delta_x[2];
		//p1[0] += delta_x[3];
		//p1[1] += delta_x[4];
		//p1[2] += delta_x[5];		
		//lambda += delta_lambda;


 


	//Vector3d p0_current;
	//Vector3d p1_current;

	//SUM(p0_current, p0, delta_x);
	//p1_current[0]=p1[0] + delta_x[3];
	//p1_current[1]=p1[1] + delta_x[4];
	//p1_current[2]=p1[2] + delta_x[5];

	//double constraint =(p0_current - p1_current).norm() - rest_length;
	//constraint = 0.5 * constraint * constraint;
	//n_norm = sqrt(EDGE_LENGTH(p0, p1));
	//double constraint_1  = 1.0 / (2.0 * time_step * time_step) * (mass_0 * EDGE_LENGTH(p0, sn_0) + mass_1 * EDGE_LENGTH(p1, sn_1)) + 0.5 * stiffness * (n_norm - rest_length) * (n_norm - rest_length);


	//if (edge_index == 68) {
	//	//std::cout << delta_x.transpose() << std::endl;
	//	//std::cout << sn_0[0] << " " << sn_0[1] << " " << sn_0[2] << std::endl;
	//}

	//if (constraint_0 >= constraint_1) {
	//	//std::cout << "true " << std::endl;
	//}
	//else {
		//std::cout << "false "<<constraint_0<<" "<<constraint_1 << std::endl;
		//std::cout << delta_x.transpose() << std::endl;
		//std::cout << p0[0] - delta_x[0] << " " << p0[1] - delta_x[1] << " " << p0[2] - delta_x[2] << std::endl;
		//std::cout << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
		//std::cout << std::endl;
		//std::cout << p1[0] - delta_x[3] << " " << p1[1] - delta_x[4] << " " << p1[2] - delta_x[5] << std::endl;
		//std::cout << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
		//std::cout << edge_index << std::endl;
		//JacobiSVD<MatrixXd> svd;
		//svd.compute(H);
		//std::cout << H.determinant() << std::endl;
	//}
	

}



double SecondOrderConstraint::solveBendingConstraint(double* center_vertex, std::array<double, 3>* vertex_position, unsigned int* neighbor_vertex, unsigned int neighbor_vertex_size,
	double rest_Aq_norm, double lbo_weight, VectorXd& vertex_lbo)
{
	std::vector<VectorXd>q(3);
	//std::vector<VectorXd>q_initial(3);
	double aq[3];
	unsigned int size = neighbor_vertex_size + 1;
	VectorXd inv_m(size);
	for (unsigned int j = 0; j < 3; ++j) {
		q[j].resize(size);
	}
	q[0][0] = center_vertex[0];
	q[1][0] = center_vertex[1];
	q[2][0] = center_vertex[2];


	for (unsigned int h = 1; h < size; h++) {
		q[0][h] = vertex_position[*neighbor_vertex][0];
		q[1][h] = vertex_position[*neighbor_vertex][1];
		q[2][h] = vertex_position[*neighbor_vertex][2];
		neighbor_vertex++;
	}
	neighbor_vertex -= neighbor_vertex_size;

	aq[0] = vertex_lbo.dot(q[0]);
	aq[1] = vertex_lbo.dot(q[1]);
	aq[2] = vertex_lbo.dot(q[2]);

	double aq_norm = sqrt(DOT(aq, aq));
	if (aq_norm < epsilon_for_bending) {
		return 0;
	}
	//double K = aq_norm - rest_Aq_norm;
	VectorXd ATAx;
	ATAx.resize(3*size);
	for (unsigned int j = 0; j < 3; ++j) {
		ATAx.segment(j*size,size) = vertex_lbo * aq[j];
	}
	MatrixXd ATA_part = ((1.0- rest_Aq_norm/ aq_norm)*vertex_lbo) * vertex_lbo.transpose();
	MatrixXd H(3 * size, 3 * size);
	H = ((rest_Aq_norm / (aq_norm * aq_norm * aq_norm)) * ATAx) * ATAx.transpose();
	for (unsigned int i = 0; i < 3; ++i) {
		H.block(i * size, i * size, size, size) += ATA_part;
	}
	ATAx *= (aq_norm - rest_Aq_norm) / aq_norm;
	ColPivHouseholderQR <MatrixXd> linear(H);
	VectorXd delta_x = linear.solve(ATAx);

	//center_vertex[0] += delta_x[0];
	//center_vertex[1] += delta_x[size];
	//center_vertex[2] += delta_x[size+size];
	//for (unsigned int h = 1; h < size; h++) {
	//	vertex_position[*neighbor_vertex][0] += delta_x[h];
	//	vertex_position[*neighbor_vertex][1] += delta_x[size + h];
	//	vertex_position[*neighbor_vertex][2] += delta_x[size+size+h];
	//	neighbor_vertex++;
	//}


	//JacobiSVD<MatrixXd> svd;
	//svd.compute(H);
	//std::cout << svd.singularValues() << std::endl;


	q[0] -= delta_x.segment(0, size);
	q[1] -= delta_x.segment(size, size);
	q[2] -= delta_x.segment(size<<1, size);
	aq[0] = vertex_lbo.dot(q[0]);
	aq[1] = vertex_lbo.dot(q[1]);
	aq[2] = vertex_lbo.dot(q[2]);

	double constraint = sqrt(DOT(aq, aq)) - rest_Aq_norm;
	constraint = 0.5 * constraint * constraint;
	return constraint;
}





double SecondOrderConstraint::solveEdgeLengthConstraintFirstOrder(double* p0, double* p1, const double rest_length)
{
	double n[3];
	SUB(n, p0, p1);
	double n_norm = sqrt(DOT(n, n));
	DEV_(n, n_norm);

	double delta_lambda = -0.5* (n_norm - rest_length);//(1.0+gamma)*(inv_mass0 + inv_mass1) 
	MULTI_(n, delta_lambda);
	//SUM_(p0, n);
	//SUB_(p1, n);

	Vector3d p0_current;
	Vector3d p1_current;
	SUM(p0_current, p0, n);
	SUB(p1_current, p1, n);

	double constraint = (p0_current - p1_current).norm() - rest_length;
	constraint = 0.5 * constraint * constraint;

	//std::cout << p0_current[0] << " " << p0_current[1] << p0_current[2] << std::endl;
	//std::cout << p1_current[0] << " " << p1_current[1] << p1_current[2] << std::endl;


	return constraint;
}





double SecondOrderConstraint::solveBendingConstraintFirstOrder(double* center_vertex, std::array<double, 3>* vertex_position, unsigned int* neighbor_vertex, unsigned int neighbor_vertex_size,
	double rest_Aq_norm, double lbo_weight, VectorXd& vertex_lbo)
{
	//	energy = 0.0;
	std::vector<VectorXd>q(3);
	//std::vector<VectorXd>q_initial(3);
	double aq[3];
	unsigned int size = neighbor_vertex_size + 1;
	for (unsigned int j = 0; j < 3; ++j) {
		q[j].resize(size);
	}
	q[0][0] = center_vertex[0];
	q[1][0] = center_vertex[1];
	q[2][0] = center_vertex[2];


	for (unsigned int h = 1; h < size; h++) {
		q[0][h] = vertex_position[*neighbor_vertex][0];
		q[1][h] = vertex_position[*neighbor_vertex][1];
		q[2][h] = vertex_position[*neighbor_vertex][2];
		neighbor_vertex++;
	}
	neighbor_vertex -= neighbor_vertex_size;


	aq[0] = vertex_lbo.dot(q[0]);
	aq[1] = vertex_lbo.dot(q[1]);
	aq[2] = vertex_lbo.dot(q[2]);

	double aq_norm = sqrt(DOT(aq, aq));
	if (aq_norm < epsilon_for_bending) {
		return 0;
	}
	DEV_(aq, aq_norm);
	double K = aq_norm - rest_Aq_norm;
	//use grad to store delta_C

	std::vector<VectorXd>grad(3);
	//std::vector<VectorXd>q_initial(3);

	for (unsigned int j = 0; j < 3; ++j) {
		//q_initial[j] = q[j] - q_initial[j];
		grad[j] = vertex_lbo * aq[j];
	}

	
	double delta_lambda = -0.5*K
		/ (grad[0].squaredNorm() + grad[1].squaredNorm()+ grad[2].squaredNorm());

	center_vertex[0] += delta_lambda * grad[0][0];
	center_vertex[1] += delta_lambda * grad[1][0];
	center_vertex[2] += delta_lambda * grad[2][0];

	for (unsigned int h = 1; h < size; h++) {
		vertex_position[*neighbor_vertex][0] += delta_lambda * grad[0][h];
		vertex_position[*neighbor_vertex][1] += delta_lambda * grad[1][h];
		vertex_position[*neighbor_vertex][2] += delta_lambda * grad[2][h];
		neighbor_vertex++;
	}



	q[0] += delta_lambda* grad[0];
	q[1] += delta_lambda* grad[1];
	q[2] += delta_lambda* grad[2];

	aq[0] = vertex_lbo.dot(q[0]);
	aq[1] = vertex_lbo.dot(q[1]);
	aq[2] = vertex_lbo.dot(q[2]);

	double constraint = sqrt(DOT(aq, aq)) - rest_Aq_norm;
	constraint = 0.5 * constraint * constraint;

	return constraint;
}


void SecondOrderConstraint::test(MeshStruct& mesh_struct, std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo,
	std::vector<double>& rest_mean_curvature_norm)
{
	double constriant_1, constraint_2;

	//std::cout << "edge length " << std::endl;

	//for(unsigned int i = 0; i < mesh_struct.edge_length.size(); ++i) {
	//	constraint_2 = solveEdgeLengthConstraint(mesh_struct.vertex_position[mesh_struct.edge_vertices[i < 1]].data(),
	//		mesh_struct.vertex_position[mesh_struct.edge_vertices[(i < 1) + 1]].data(), mesh_struct.edge_length[i]);
	//	constriant_1 = solveEdgeLengthConstraintFirstOrder(mesh_struct.vertex_position[mesh_struct.edge_vertices[i < 1]].data(),
	//		mesh_struct.vertex_position[mesh_struct.edge_vertices[(i < 1) + 1]].data(), mesh_struct.edge_length[i]);
	//	if (constraint_2 < constriant_1) {
	//		std::cout << "true ";
	//	}
	//	else {
	//		std::cout << "false ";
	//	}
	//	std::cout << constraint_2 << " " << constriant_1 << std::endl;
	//}


	//std::cout << "bending " << std::endl;

	for (unsigned int i = 0; i < mesh_struct.vertex_position.size(); ++i) {
		constraint_2 = solveBendingConstraint(mesh_struct.vertex_position[i].data(), mesh_struct.vertex_position.data(), mesh_struct.vertices[i].neighbor_vertex.data(),
			mesh_struct.vertices[i].neighbor_vertex.size(), rest_mean_curvature_norm[i], lbo_weight[i], vertex_lbo[i]);
		constriant_1 = solveBendingConstraintFirstOrder(mesh_struct.vertex_position[i].data(), mesh_struct.vertex_position.data(), mesh_struct.vertices[i].neighbor_vertex.data(),
			mesh_struct.vertices[i].neighbor_vertex.size(), rest_mean_curvature_norm[i], lbo_weight[i], vertex_lbo[i]);
		if (constraint_2 < constriant_1) {
			std::cout << "true ";
		}
		else {
			std::cout << "false ";
		}
		std::cout << constraint_2 << " " << constriant_1 << std::endl;
	}
}