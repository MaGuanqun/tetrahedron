#include"second_order.h"


void SecondOrderConstraint::computeForce(double* vertex_0, double* vertex_1, double stiffness, 
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

void SecondOrderConstraint::solveARAPConstraint(std::array<double, 3>* vertex_position, std::array<double, 3>* initial_vertex_position,
	double stiffness, double dt,
	Matrix<double, 3, 4>& A, int* vertex_index, double* inv_mass, double& lambda, const double damping_stiffness, double sigma_min,
	double sigma_max, double volume, double& energy, double* mass)
{
	Vector3d eigen_value;
	Matrix3d q_e;
	double determinant;
	Vector3d position;
	for (unsigned int i = 0; i < 3; ++i) {
		memcpy(q_e.data() + 3 * i, vertex_position[vertex_index[i + 1]].data(), 24);
	}
	//first use eigen value to store the position of first vertex
	memcpy(eigen_value.data(), vertex_position[vertex_index[0]].data(), 24);
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


	Matrix4d ATA_c = A.transpose() * (A*(lambda/C));

	
	Matrix<double, 12, 1> grad;
	// grad=grad_C_transpose^T
	grad.data()[0] = grad_C_transpose.data()[0]; 	grad.data()[1] = grad_C_transpose.data()[3];	grad.data()[2] = grad_C_transpose.data()[6];	grad.data()[3] = grad_C_transpose.data()[9];
	grad.data()[4] = grad_C_transpose.data()[1]; 	grad.data()[5] = grad_C_transpose.data()[4];	grad.data()[6] = grad_C_transpose.data()[7];	grad.data()[7] = grad_C_transpose.data()[10];
	grad.data()[8] = grad_C_transpose.data()[2]; 	grad.data()[9] = grad_C_transpose.data()[5];	grad.data()[10] = grad_C_transpose.data()[8];	grad.data()[11] = grad_C_transpose.data()[11];


	//	ColPivHouseholderQR <Matrix3d> linear(sys_matrix);

	

	

}

void SecondOrderConstraint::getdPdF(Matrix3d& U, Matrix3d& V, Vector3d& eigen_value, Matrix<double,9,9>& dPdF)
{

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
			Matrix<double, 3, 3> sys_matrix;
			sys_matrix = mass_1 * Matrix3d::Identity()-lambda * He;
			//sys_matrix.setZero();
			//sys_matrix -= He;
			//sys_matrix.data()[0] += mass_1;
			//sys_matrix.data()[4] += mass_1;
			//sys_matrix.data()[8] += mass_1;
			Vector3d g;
			g.data()[0] = mass_1 * (p1[0] - sn_1[0]);
			g.data()[1] = mass_1 * (p1[1] - sn_1[1]);
			g.data()[2] = mass_1 * (p1[2] - sn_1[2]);
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
			Matrix<double, 3, 3> sys_matrix;
			sys_matrix = mass_0 * Matrix3d::Identity() -lambda * He;
			//sys_matrix.setZero();
			//sys_matrix -= He;
			//sys_matrix.data()[0] += mass_0;
			//sys_matrix.data()[4] += mass_0;
			//sys_matrix.data()[8] += mass_0;
			Vector3d g;
			g.data()[0] = mass_0 * (p0[0] - sn_0[0]);
			g.data()[1] = mass_0 * (p0[1] - sn_0[1]);
			g.data()[2] = mass_0 * (p0[2] - sn_0[2]);
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
		g.data()[0] =mass_0 *(p0[0] - sn_0[0]);
		g.data()[1] =mass_0 *(p0[1] - sn_0[1]);
		g.data()[2] =mass_0 *(p0[2] - sn_0[2]);
		g.data()[3] = mass_1 * (p1[0] - sn_1[0]);
		g.data()[4] = mass_1 * (p1[1] - sn_1[1]);
		g.data()[5] = mass_1 * (p1[2] - sn_1[2]);
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



		//////check:
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