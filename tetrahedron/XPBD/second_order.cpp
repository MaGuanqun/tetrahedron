#include"second_order.h"


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


void SecondOrderConstraint::solveARAPConstraint(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	double stiffness, double dt,
	Matrix<double, 3, 4>& A, double* inv_mass, double& lambda, const double damping_stiffness, double sigma_min,
	double sigma_max, double volume)
{
	Vector3d eigen_value;
	Matrix3d q_e;
	double determinant;
	Vector3d position;
	memcpy(q_e.data(), vertex_position_1, 24);
	memcpy(q_e.data()+3, vertex_position_2, 24);
	memcpy(q_e.data()+6, vertex_position_3, 24);
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

	Matrix3d U = svd.matrixU();
	Matrix3d V = svd.matrixV();
	if (determinant < 0) {
		V.col(2) *= -1.0;
		eigen_value[2] *= -1.0;
	}

	//use P_inv to record transform
	P_inv = U *V.transpose();
	double C = (deformation_gradient - P_inv).norm();

	if (C < 1e-8) {
		return;
	}

	Matrix<double, 3, 4> grad_C_transpose;

	grad_C_transpose = (1.0 / C) * (deformation_gradient - P_inv) * A;
	double alpha_ = 1.0 / (stiffness*volume * dt * dt);

	
	Matrix<double, 12, 1> grad;
	memcpy(grad.data(), grad_C_transpose.data(), 96);

	Matrix<double, 9, 9> dPdF;

	getdPdF(U, V, eigen_value, dPdF);

	Matrix<double, 12, 12> Hessian;
	backpropagateElementHessian(Hessian, dPdF, A);
	Hessian *= (0.5 / C);
	Hessian -= ((1.0 / C) * grad) * grad.transpose();

	Matrix<double, 12, 1> inv_mass_sys;
	inv_mass_sys << inv_mass[0], inv_mass[0], inv_mass[0], inv_mass[1], inv_mass[1], inv_mass[1], inv_mass[2], inv_mass[2], inv_mass[2],
		inv_mass[3], inv_mass[3], inv_mass[3];
	//how to improve this
	for (unsigned int j = 0; j < 12; ++j) {
		for (unsigned int i = 0; i < 4; ++i) {
			Hessian(i + i + i, j) *= inv_mass[i];
			Hessian(i + i + i+1, j) *= inv_mass[i];
			Hessian(i + i + i+2, j) *= inv_mass[i];
		}
	}
	Hessian *= -lambda;
	for (unsigned int i = 0; i < 144; i += 13) {
		Hessian.data()[i] += 1.0;
	}
	ColPivHouseholderQR <Matrix<double,12,12>> linear(Hessian);
	double h = C + alpha_ * lambda;
	double delta_lambda = -h / (grad.dot(linear.solve(inv_mass_sys.cwiseProduct(grad))) + alpha_);
	Matrix<double, 12, 1> delta_x = linear.solve(delta_lambda * (inv_mass_sys.cwiseProduct(grad)));

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

//the target is to solve A^T *dPdF* A
void SecondOrderConstraint::backpropagateElementHessian(Matrix<double, 12,12>& Hessian, Matrix<double,9,9>& dPdF, Matrix<double,3,4>&A)
{
	Matrix<double, 12, 9> temp;//temp = A^T dPdF
	for (unsigned int i = 0; i < 9; ++i) {
		temp(3, i) = A(0, 1) * dPdF(0, i) + A(1, 1) * dPdF(3, i) + A(2, 1) * dPdF(6, i);
		temp(4, i) = A(0, 1) * dPdF(1, i) + A(1, 1) * dPdF(4, i) + A(2, 1) * dPdF(7, i);
		temp(5, i) = A(0, 1) * dPdF(2, i) + A(1, 1) * dPdF(5, i) + A(2, 1) * dPdF(8, i);

		temp(6, i) = A(0, 2) * dPdF(0, i) + A(1, 2) * dPdF(3, i) + A(2, 2) * dPdF(6, i);
		temp(7, i) = A(0, 2) * dPdF(1, i) + A(1, 2) * dPdF(4, i) + A(2, 2) * dPdF(7, i);
		temp(8, i) = A(0, 2) * dPdF(2, i) + A(1, 2) * dPdF(5, i) + A(2, 2) * dPdF(8, i);

		temp(9, i) = A(0, 3) * dPdF(0, i) + A(1, 3) * dPdF(3, i) + A(2, 3) * dPdF(6, i);
		temp(10, i) = A(0, 3) * dPdF(1, i) + A(1, 3) * dPdF(4, i) + A(2, 3) * dPdF(7, i);
		temp(11, i) = A(0, 3) * dPdF(2, i) + A(1, 3) * dPdF(5, i) + A(2, 3) * dPdF(8, i);

		temp(0, i) = -temp(3, i) - temp(6, i) - temp(9, i);
		temp(1, i) = -temp(4, i) - temp(7, i) - temp(10, i);
		temp(2, i) = -temp(5, i) - temp(8, i) - temp(11, i);
	}

	//temp*A
	for (int i = 0; i < 12; ++i) {
		Hessian(i, 3) = temp(i, 0) * A(0, 1) + temp(i, 3) * A(1, 1) + temp(i, 6) * A(2, 1);
		Hessian(i, 4) = temp(i, 1) * A(0, 1) + temp(i, 4) * A(1, 1) + temp(i, 7) * A(2, 1);
		Hessian(i, 5) = temp(i, 2) * A(0, 1) + temp(i, 5) * A(1, 1) + temp(i, 8) * A(2, 1);

		Hessian(i, 6) = temp(i, 0) * A(0, 2) + temp(i, 3) * A(1, 2) + temp(i, 6) * A(2, 2);
		Hessian(i, 7) = temp(i, 1) * A(0, 2) + temp(i, 4) * A(1, 2) + temp(i, 7) * A(2, 2);
		Hessian(i, 8) = temp(i, 2) * A(0, 2) + temp(i, 5) * A(1, 2) + temp(i, 8) * A(2, 2);

		Hessian(i, 9) = temp(i, 0) * A(0, 3) + temp(i, 3) * A(1, 3) + temp(i, 6) * A(2, 3);
		Hessian(i, 10) = temp(i, 1) * A(0, 3) + temp(i, 4) * A(1, 3) + temp(i, 7) * A(2, 3);
		Hessian(i, 11) = temp(i, 2) * A(0, 3) + temp(i, 5) * A(1, 3) + temp(i, 8) * A(2, 3);

		Hessian(i, 0) = -Hessian(i, 3) - Hessian(i, 6) - Hessian(i, 9);
		Hessian(i, 1) = -Hessian(i, 4) - Hessian(i, 7) - Hessian(i, 10);
		Hessian(i, 2) = -Hessian(i, 5) - Hessian(i, 8) - Hessian(i, 11);
	}

}


void SecondOrderConstraint::getdPdF(Matrix3d& U, Matrix3d& V, Vector3d& eigen_value, Matrix<double,9,9>& dPdF)
{

	
	Matrix<double, 2, 2> B0;
	double coe = 1.0;
	if (eigen_value[0] + eigen_value[1] < 1e-6) {
		coe = (eigen_value[0] + eigen_value[1]) / 1e-6;
	}
	B0 << 1.0 + coe, 1.0 - coe, 1.0 - coe, 1.0 + coe;

	Matrix<double, 2, 2> B1;
	coe = 1.0;
	if (eigen_value[1] + eigen_value[2] < 1e-6) {
		coe = (eigen_value[1] + eigen_value[2]) / 1e-6;
	}
	B1 << 1.0 + coe, 1.0 - coe, 1.0 - coe, 1.0 + coe;

	Matrix<double, 2, 2> B2;
	coe = 1.0;
	if (eigen_value[0] + eigen_value[2] < 1e-6) {
		coe = (eigen_value[0] + eigen_value[2]) / 1e-6;
	}
	B2 << 1.0 + coe, 1.0 - coe, 1.0 - coe, 1.0 + coe;


	Matrix<double, 9, 9> M;
	M.setZero();
	M(0, 0) = eigen_value[0] + eigen_value[0];
	M(4, 4) = eigen_value[1] + eigen_value[1];
	M(8, 8) = eigen_value[2] + eigen_value[2];
	M(1, 1) = B0(0, 0);
	M(1, 3) = B0(0, 1);
	M(3, 1) = B0(1, 0);
	M(3, 3) = B0(1, 1);
	M(5, 5) = B1(0, 0);
	M(5, 7) = B1(0, 1);
	M(7, 5) = B1(1, 0);
	M(7, 7) = B1(1, 1);
	M(2, 2) = B2(1, 1);
	M(2, 6) = B2(1, 0);
	M(6, 2) = B2(0, 1);
	M(6, 6) = B2(0, 0);

	for (int j = 0; j < 3; ++j)
		for (int i = 0; i < 3; ++i)
			for (int s = 0; s < 3; ++s)
				for (int r = 0; r < 3; ++r) {
					int ij = j * 3 + i;
					int rs = s * 3 + r;
					dPdF(ij, rs) = M(0, 0) * U(i, 0) * V(j, 0) * U(r, 0) * V(s, 0)
						+ M(4, 4) * U(i, 1) * V(j, 1) * U(r, 1) * V(s, 1)
						+ M(8, 8) * U(i, 2) * V(j, 2) * U(r, 2) * V(s, 2)
						+ M(1, 1) * U(i, 0) * V(j, 1) * U(r, 0) * V(s, 1)
						+ M(1, 3) * U(i, 0) * V(j, 1) * U(r, 1) * V(s, 0)
						+ M(3, 1) * U(i, 1) * V(j, 0) * U(r, 0) * V(s, 1)
						+ M(3, 3) * U(i, 1) * V(j, 0) * U(r, 1) * V(s, 0)
						+ M(5, 5) * U(i, 1) * V(j, 2) * U(r, 1) * V(s, 2)
						+ M(5, 7) * U(i, 1) * V(j, 2) * U(r, 2) * V(s, 1)
						+ M(7, 5) * U(i, 2) * V(j, 1) * U(r, 1) * V(s, 2)
						+ M(7, 7) * U(i, 2) * V(j, 1) * U(r, 2) * V(s, 1)
						+ M(2, 2) * U(i, 0) * V(j, 2) * U(r, 0) * V(s, 2)
						+ M(2, 6) * U(i, 0) * V(j, 2) * U(r, 2) * V(s, 0)
						+ M(6, 2) * U(i, 2) * V(j, 0) * U(r, 0) * V(s, 2)
						+ M(6, 6) * U(i, 2) * V(j, 0) * U(r, 2) * V(s, 0);
					//dPdF(ij, rs) = M(0, 0) * U(i, 0) * V(j, 0) * U(r, 0) * V(s, 0)
					//	+ M(0, 4) * U(i, 0) * V(j, 0) * U(r, 1) * V(s, 1)
					//	+ M(0, 8) * U(i, 0) * V(j, 0) * U(r, 2) * V(s, 2)
					//	+ M(4, 0) * U(i, 1) * V(j, 1) * U(r, 0) * V(s, 0)
					//	+ M(4, 4) * U(i, 1) * V(j, 1) * U(r, 1) * V(s, 1)
					//	+ M(4, 8) * U(i, 1) * V(j, 1) * U(r, 2) * V(s, 2)
					//	+ M(8, 0) * U(i, 2) * V(j, 2) * U(r, 0) * V(s, 0)
					//	+ M(8, 4) * U(i, 2) * V(j, 2) * U(r, 1) * V(s, 1)
					//	+ M(8, 8) * U(i, 2) * V(j, 2) * U(r, 2) * V(s, 2)
					//	+ M(1, 1) * U(i, 0) * V(j, 1) * U(r, 0) * V(s, 1)
					//	+ M(1, 3) * U(i, 0) * V(j, 1) * U(r, 1) * V(s, 0)
					//	+ M(3, 1) * U(i, 1) * V(j, 0) * U(r, 0) * V(s, 1)
					//	+ M(3, 3) * U(i, 1) * V(j, 0) * U(r, 1) * V(s, 0)
					//	+ M(5, 5) * U(i, 1) * V(j, 2) * U(r, 1) * V(s, 2)
					//	+ M(5, 7) * U(i, 1) * V(j, 2) * U(r, 2) * V(s, 1)
					//	+ M(7, 5) * U(i, 2) * V(j, 1) * U(r, 1) * V(s, 2)
					//	+ M(7, 7) * U(i, 2) * V(j, 1) * U(r, 2) * V(s, 1)
					//	+ M(2, 2) * U(i, 0) * V(j, 2) * U(r, 0) * V(s, 2)
					//	+ M(2, 6) * U(i, 0) * V(j, 2) * U(r, 2) * V(s, 0)
					//	+ M(6, 2) * U(i, 2) * V(j, 0) * U(r, 0) * V(s, 2)
					//	+ M(6, 6) * U(i, 2) * V(j, 0) * U(r, 2) * V(s, 0);
				}


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