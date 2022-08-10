#include"XPBD_constraint.h"



void XPBDconstraint::PBDsolveARAPConstraint(std::array<double, 3>* vertex_position, std::array<double, 3>* initial_vertex_position,
	double stiffness, double dt,
	Matrix<double,3,4>& A, int* vertex_index, double* inv_mass, double volume, double iteration_num_inverse)
{
	Vector3d eigen_value;
	Matrix3d q_e;
	double determinant;
	Vector3d position;
	for (unsigned int i = 0; i < 3; ++i) {
		memcpy(q_e.data() + 3 * i, vertex_position[vertex_index[i+1]].data(), 24);
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

	//for (unsigned int i = 0; i < 4; ++i) {
	//	//if (vertex_index[i] == 0) {
	//		//for (unsigned int i = 0; i < 4; ++i) {
	//			std::cout << vertex_index[i]<<" "<< vertex_position[vertex_index[i]][0] << " ";
	//		//}
	//		//std::cout << q_e << std::endl;
	//		//break;
	//	//}
	//}
	//std::cout << std::endl;
	JacobiSVD<Matrix3d> svd;
	svd.compute(deformation_gradient, ComputeFullU | ComputeFullV);

	eigen_value = svd.singularValues();
	determinant = eigen_value[0] * eigen_value[1] * eigen_value[2];

	q_e = svd.matrixU();
	if (determinant < 0) {
		q_e.col(2) *= -1.0;
	}

	//use P_inv to record transform
	P_inv = q_e * svd.matrixV().transpose();


	//get delta_c
	Matrix<double, 3, 4> grad_C_transpose;


	grad_C_transpose = volume * (deformation_gradient - P_inv) * A;

	double C = 0.5 * volume * (deformation_gradient - P_inv).squaredNorm();

	if (C < 1e-20) {
		return;
	}

	Vector3d position_;

	double delta_lambda_denominator = 0.0;
	for (unsigned int k = 0; k < 4; ++k) {
		delta_lambda_denominator += inv_mass[vertex_index[k]] * grad_C_transpose.col(k).squaredNorm();
	}
	
	double s = (C / delta_lambda_denominator) * (1.0 - pow(1.0 - stiffness, iteration_num_inverse));

	double coe;
	for (unsigned int k = 0; k < 4; ++k) {
		coe = inv_mass[vertex_index[k]] * s;
		vertex_position[vertex_index[k]][0] -= coe * grad_C_transpose.data()[3 * k];
		vertex_position[vertex_index[k]][1] -= coe * grad_C_transpose.data()[3 * k + 1];
		vertex_position[vertex_index[k]][2] -= coe * grad_C_transpose.data()[3 * k + 2];
	}

	//std::cout << vertex_index[0]<<" "<< vertex_index[1]<<" "<< vertex_index[2]<<" "<< vertex_index[3]<<" "<<  coe << std::endl;
	//std::cout << deformation_gradient << std::endl;
}






void XPBDconstraint::solveARAPConstraint(std::array<double, 3>* vertex_position, std::array<double, 3>* initial_vertex_position,
	double stiffness, double dt,
	Matrix<double, 3, 4>& A, int* vertex_index, double* inv_mass, double& lambda, const double damping_stiffness, double sigma_min,
	double sigma_max, double volume, double& energy)
{
	energy = 0.0;
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

	q_e = svd.matrixU();
	if (determinant < 0) {
		q_e.col(2) *= -1.0;
	}


	//if (determinant < 0) {
	//	eigen_value[2] = -eigen_value[2];
	//}
	//for (unsigned int j = 0; j < 3; ++j) {
	//	if (eigen_value[j] > 0) {
	//		if (eigen_value[j] < sigma_min) {
	//			eigen_value[j] = sigma_min;
	//		}
	//		else if (eigen_value[j] > sigma_max) {
	//			eigen_value[j] = sigma_max;
	//		}
	//	}
	//	else {
	//		if (eigen_value[j] > -sigma_min) {
	//			eigen_value[j] = -sigma_min;
	//		}
	//		else if (eigen_value[j] < -sigma_max) {
	//			eigen_value[j] = -sigma_max;
	//		}
	//	}
	//}
	////use q_e as a temp vector
	//for (unsigned int j = 0; j < 3; ++j) {
	//	for (unsigned int k = 0; k < 3; ++k) {
	//		q_e.data()[3 * j + k] = eigen_value[j] * svd.matrixU().data()[3 * j + k];
	//	}
	//}

	//use P_inv to record transform
	P_inv = q_e * svd.matrixV().transpose();

	//if((deformation_gradient-P_inv).squaredNorm)

	//get delta_c
	Matrix<double, 3, 4> grad_C_transpose;


	grad_C_transpose = volume * (deformation_gradient - P_inv) * A;

	double C = 0.5 * volume * (deformation_gradient - P_inv).squaredNorm();






	double alpha_ = 1.0 / (stiffness * dt * dt);
	double gamma = damping_stiffness / (stiffness * dt);

	Vector3d position_;

	double delta_lambda_numerator = 0.0;
	for (unsigned int k = 0; k < 4; ++k) {
		SUB(position_, vertex_position[vertex_index[k]], initial_vertex_position[vertex_index[k]]);
		delta_lambda_numerator += grad_C_transpose.col(k).dot(position_);
	}


	double delta_lambda_denominator = 0.0;
	for (unsigned int k = 0; k < 4; ++k) {
		delta_lambda_denominator += inv_mass[vertex_index[k]] * grad_C_transpose.col(k).squaredNorm();
	}
	double delta_lambda = -(C + alpha_ * lambda + gamma * delta_lambda_numerator)
		/ ((1.0 + gamma) * delta_lambda_denominator + alpha_);
	lambda += delta_lambda;

	double coe;
	for (unsigned int k = 0; k < 4; ++k) {
		coe = inv_mass[vertex_index[k]] * delta_lambda;
		vertex_position[vertex_index[k]][0] += coe * grad_C_transpose.data()[3 * k];
		vertex_position[vertex_index[k]][1] += coe * grad_C_transpose.data()[3 * k + 1];
		vertex_position[vertex_index[k]][2] += coe * grad_C_transpose.data()[3 * k + 2];
	}


	energy = 0.5 * stiffness * C * C;
}

void XPBDconstraint::solveARAPConstraint2(std::array<double, 3>* original_vertex_pos ,std::array<double, 3>* vertex_position, std::array<double, 3>* initial_vertex_position,
	double stiffness, double dt,
	Matrix<double, 3, 4>& A, int* vertex_index, double* inv_mass, double& lambda, const double damping_stiffness, double sigma_min,
	double sigma_max, double volume)
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

	q_e = svd.matrixU();
	if (determinant < 0) {
		q_e.col(2) *= -1.0;
	}


	

	//q_e = svd.matrixU();
	//if (determinant < 0) {
	//	eigen_value[2] = -eigen_value[2];
	//}
	//for (unsigned int j = 0; j < 3; ++j) {
	//	if (eigen_value[j] > 0) {
	//		if (eigen_value[j] < sigma_min) {
	//			eigen_value[j] = sigma_min;
	//		}
	//		else if (eigen_value[j] > sigma_max) {
	//			eigen_value[j] = sigma_max;
	//		}
	//	}
	//	else {
	//		if (eigen_value[j] > -sigma_min) {
	//			eigen_value[j] = -sigma_min;
	//		}
	//		else if (eigen_value[j] < -sigma_max) {
	//			eigen_value[j] = -sigma_max;
	//		}
	//	}
	//}
	////use q_e as a temp vector
	//for (unsigned int j = 0; j < 3; ++j) {
	//	for (unsigned int k = 0; k < 3; ++k) {
	//		q_e.data()[3 * j + k] = eigen_value[j] * svd.matrixU().data()[3 * j + k];
	//	}
	//}

	//use P_inv to record transform
	P_inv = q_e * svd.matrixV().transpose();

	Matrix<double, 3, 4> p_ori;
	for (unsigned int i = 0; i < 4; ++i) {
		memcpy(p_ori.data() + 3 * i, original_vertex_pos[vertex_index[i]].data(), 24);
	}
	Vector3d center = 0.25 * (p_ori.col(0) + p_ori.col(1) + p_ori.col(2) + p_ori.col(3));
	for (unsigned int i = 0; i < 4; ++i) {
		p_ori.col(i) -= center;
	}
	Matrix<double, 3, 4> p_target;
	p_target = P_inv * p_ori;

	//use p_ori to store the current_position
	for (unsigned int i = 0; i < 4; ++i) {
		memcpy(p_ori.data() + 3 * i, vertex_position[vertex_index[i]].data(), 24);
	}
	center = 0.25 * (p_ori.col(0) + p_ori.col(1) + p_ori.col(2) + p_ori.col(3));
	for (unsigned int i = 0; i < 4; ++i) {
		p_ori.col(i) -= center;
	}

	//get delta_c
	Matrix<double, 3, 4> grad_C;
	grad_C = volume * (p_ori - p_target);

	double C = 0.5 * volume * (p_ori - p_target).squaredNorm();

	double alpha_ = 1.0 / (stiffness * dt * dt);
	double gamma = damping_stiffness / (stiffness * dt);

	Vector3d position_;

	double delta_lambda_numerator = 0.0;
	for (unsigned int k = 0; k < 4; ++k) {
		SUB(position_, vertex_position[vertex_index[k]], initial_vertex_position[vertex_index[k]]);
		delta_lambda_numerator += grad_C.col(k).dot(position_);
	}


	double delta_lambda_denominator = 0.0;
	for (unsigned int k = 0; k < 4; ++k) {
		delta_lambda_denominator += inv_mass[vertex_index[k]] * grad_C.col(k).squaredNorm();
	}
	double delta_lambda = -(C + alpha_ * lambda + gamma * delta_lambda_numerator)
		/ ((1.0 + gamma) * delta_lambda_denominator + alpha_);
	lambda += delta_lambda;

	double coe;
	for (unsigned int k = 0; k < 4; ++k) {
		coe = inv_mass[vertex_index[k]] * delta_lambda;
		vertex_position[vertex_index[k]][0] += coe * grad_C.data()[3 * k];
		vertex_position[vertex_index[k]][1] += coe * grad_C.data()[3 * k + 1];
		vertex_position[vertex_index[k]][2] += coe * grad_C.data()[3 * k + 2];
	}
}





void XPBDconstraint::solveTetStrainConstraint(std::array<double, 3>* vertex_position, std::array<double, 3>* initial_vertex_position,
	double stiffness, double dt, Matrix<double, 3, 4>& A, int* vertex_index, double* inv_mass, double& lambda, const double damping_stiffness, double volume,
	double youngs_modulus, double poisson_ratio)
{

	Matrix<double, 3, 4> grad_C;
	double C;
	double mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
	double lambda_ = youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
	computeGreenStrainAndPiolaStress(vertex_position[vertex_index[0]].data(), vertex_position[vertex_index[1]].data(),
		vertex_position[vertex_index[2]].data(), vertex_position[vertex_index[3]].data(), A, volume, mu, lambda_, grad_C, C);
	double alpha_ = 1.0 / (stiffness * dt * dt);
	double gamma = damping_stiffness / (stiffness * dt);

	double sum_normGradC = inv_mass[vertex_index[0]] * grad_C.col(0).squaredNorm() +
		inv_mass[vertex_index[1]] * grad_C.col(1).squaredNorm() +
		inv_mass[vertex_index[2]] * grad_C.col(2).squaredNorm() +
		inv_mass[vertex_index[3]] * grad_C.col(3).squaredNorm();

	Vector3d temp;
	//k = delta_C(x-x_initial)

	double k = 0;
	for (unsigned int i = 0; i < 4; ++i) {
		SUB(temp, vertex_position[vertex_index[i]], initial_vertex_position[vertex_index[i]]);
		k += grad_C.col(i).dot(temp);
	}
	double delta_lambda = -(C + alpha_ * lambda + gamma * k) / ((1.0 + gamma) * sum_normGradC + alpha_);
	lambda += delta_lambda;

	double coe;
	int vertex_index_;
	for (unsigned int i = 0; i < 4; ++i){
		vertex_index_ = vertex_index[i];
		coe = inv_mass[vertex_index_] * delta_lambda;
		ACCUMULATE_SUM_WITH_COE(vertex_position[vertex_index_], coe, (grad_C.data() + 3 * i));
	}
}


void XPBDconstraint::scondOrderStrainConstraint(std::array<double, 3>* vertex_position,
	double stiffness, double dt, Matrix<double, 3, 4>& A, int* vertex_index, double* inv_mass, double volume, double youngs_modulus, double poisson_ratio)
{
	Matrix<double, 12, 1> grad_C;
	double C;
	double mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
	Matrix<double, 12, 12>Hessian;
	double lambda_ = youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
	computeGreenStrainAndPiolaStressHessian(vertex_position[vertex_index[0]].data(), vertex_position[vertex_index[1]].data(),
		vertex_position[vertex_index[2]].data(), vertex_position[vertex_index[3]].data(), A, volume, mu, lambda_, grad_C, C, Hessian,stiffness);
	double inverse_mass;
	for (unsigned int i = 0; i < 12; i+=3) {
		inverse_mass = inv_mass[i/3];
		for (unsigned int j = 0; j < 81; j += 9) {
			Hessian.data()[j + i] *= inverse_mass;
			Hessian.data()[j + i + 1] *= inverse_mass;
			Hessian.data()[j + i + 2] *= inverse_mass;
		}
		grad_C.data()[i] *= inverse_mass;
		grad_C.data()[i+1] *= inverse_mass;
		grad_C.data()[i+2] *= inverse_mass;
	}
	double inv_t = 1.0 / (dt * dt);
	for (unsigned int j = 0; j < 81; j += 10) {
		Hessian.data()[j] += inv_t;
	}
	here need sn

}


void XPBDconstraint::test()
{
	double v0[3] = { 0.121904, -0.150605,-0.0381831 };
	double v1[3] = { 0.168903, -0.13272,-0.0500603 };
	double v2[3] = { 0.158582, -0.210298,-0.0497997 };
	double v3[3] = { 0.0990045, -0.123995,-0.0156256 };
	double k[9] = { 0.124,1.0254,0.1241,0.145,2.32,4.154,3.156,4.412,0.454 };
	Matrix3d p;
	memcpy(p.data(), k, 72);
	Matrix<double, 3, 4> A;
	for (unsigned int i = 0; i < 3; ++i) {
		A.data()[i] = -p.col(i).sum();
		A.data()[3 + i] = p.data()[3 * i];
		A.data()[6 + i] = p.data()[3 * i + 1];
		A.data()[9 + i] = p.data()[3 * i + 2];
	}
	std::cout << p << std::endl;
	std::cout << A << std::endl;
	Matrix<double, 12, 1> grad_C;
	Matrix<double, 12, 12> Hessian;;
	double C;
	computeGreenStrainAndPiolaStressHessian(v0, v1, v2, v3, A, 1.0, 3.0, 3.0, grad_C, C,  Hessian,1.0);

}

void XPBDconstraint::computeGreenStrainAndPiolaStressHessian(double* v0, double* v1, double* v2, double* v3,
	Matrix<double, 3, 4>& inv_rest_pos,
	double rest_volume,
	double mu, double lambda, Matrix<double, 12, 1>& grad_C, double& C, Matrix<double,12,12>& Hessian, double stiffness)
{
	Matrix3d F;
	Matrix3d P;
	SUB(P.data(), v1, v0);
	SUB((P.data() + 3), v2, v0);
	SUB((P.data() + 6), v3, v0);
	//first use eigen value to store the position of first vertex
	Matrix3d P_inv_transpose;
	memcpy(P_inv_transpose.data(), inv_rest_pos.data() + 3, 72);
	F = P * P_inv_transpose.transpose();

	if (F.determinant() < 0.58) {
		JacobiSVD<Matrix3d> svd;
		svd.compute(F, ComputeFullU | ComputeFullV);
		Vector3d new_singular_value = svd.singularValues();
		double min_s_value = 0.58;
		//once a critical compression
		//threshold is reached( 58 % of undeformed dimensions, when compression occurs
		//	along a single axis) the strength of the restorative force reaches a maximum.
		//Further compression will be met with decreasing resistance
		for (unsigned char j = 0; j < 3; j++) {
			if (new_singular_value[j] < min_s_value)
				new_singular_value[j] = min_s_value;
		}
		F = svd.matrixU() * new_singular_value.asDiagonal() * svd.matrixV().transpose();
	}

	Matrix3d E;

	//E= 0.5* (F^T F - I)
	E.data()[0] = 0.5 * (F.data()[0] * F.data()[0] + F.data()[1] * F.data()[1] + F.data()[2] * F.data()[2] - 1.0);		// xx
	E.data()[4] = 0.5 * (F.data()[3] * F.data()[3] + F.data()[4] * F.data()[4] + F.data()[5] * F.data()[5] - 1.0);		// yy
	E.data()[8] = 0.5 * (F.data()[6] * F.data()[6] + F.data()[7] * F.data()[7] + F.data()[8] * F.data()[8] - 1.0);		// zz
	E.data()[3] = 0.5 * (F.data()[0] * F.data()[3] + F.data()[1] * F.data()[4] + F.data()[2] * F.data()[5]);			// xy
	E.data()[6] = 0.5 * (F.data()[0] * F.data()[6] + F.data()[1] * F.data()[7] + F.data()[2] * F.data()[8]);			// xz
	E.data()[7] = 0.5 * (F.data()[3] * F.data()[6] + F.data()[4] * F.data()[7] + F.data()[5] * F.data()[8]);			// yz
	E.data()[1] = E.data()[3];
	E.data()[2] = E.data()[6];
	E.data()[5] = E.data()[7];

	//P(F)=F(2muE +lambda tr(E)I) 
	double trace = E.data()[0] + E.data()[4] + E.data()[8];
	double lambda_trace = lambda * trace;
	P = (mu + mu) * E;
	P.data()[0] += lambda_trace;
	P.data()[4] += lambda_trace;
	P.data()[8] += lambda_trace;
	P = F * P;
	//st. venant-kirchhoff model
	// psi= mu E:E + lambda/2 * tr(E)^2
	double psi = mu * E.squaredNorm() + 0.5 * lambda * trace * trace;
	C = rest_volume * psi;

	Matrix<double, 3, 4> grad = (rest_volume*stiffness) * P * inv_rest_pos;
	memcpy(grad_C.data(), grad.data(), 96);//12*8


	mu *= rest_volume* stiffness;
	lambda *= rest_volume* stiffness;
	lambda_trace *= rest_volume* stiffness;

	//first compute second derivative of mu 2FE
	Matrix<double, 9, 9>Hessian_2FE;
	Hessian_2FE.data()[0] = mu * (3.0 * F.data()[0] * F.data()[0] + F.data()[1] * F.data()[1] + F.data()[2] * F.data()[2] - 1.0 + F.data()[3] * F.data()[3] + F.data()[6] * F.data()[6]);
	std::cout << Hessian_2FE.data()[0] << std::endl;
	Hessian_2FE.data()[1] = mu * (2.0 * F.data()[0] * F.data()[1] + F.data()[3] * F.data()[4] + F.data()[6] * F.data()[7]);
	Hessian_2FE.data()[2] = mu * (2.0 * F.data()[0] * F.data()[2] + F.data()[3] * F.data()[5] + F.data()[6] * F.data()[8]);
	Hessian_2FE.data()[3] = mu * (2.0 * F.data()[0] * F.data()[3] + F.data()[1] * F.data()[4] + F.data()[2] * F.data()[5]);
	Hessian_2FE.data()[4] = mu *(F.data()[3] * F.data()[1]);
	Hessian_2FE.data()[5] = mu *(F.data()[3] * F.data()[2]);
	Hessian_2FE.data()[6] = mu * (2.0 * F.data()[6] * F.data()[0] + F.data()[1] * F.data()[7] + F.data()[2] * F.data()[8]);
	Hessian_2FE.data()[7] = mu *(F.data()[6] * F.data()[1]);
	Hessian_2FE.data()[8] = mu *(F.data()[6] * F.data()[2]);
	///
	Hessian_2FE.data()[10] = mu *(F.data()[0] * F.data()[0] + 3.0*F.data()[1] * F.data()[1] + F.data()[2] * F.data()[2] - 1.0 + F.data()[4] * F.data()[4] + F.data()[7] * F.data()[7]);
	Hessian_2FE.data()[11] = mu * (2.0 * F.data()[1] * F.data()[2] + F.data()[5] * F.data()[4] + F.data()[8] * F.data()[7]);
	Hessian_2FE.data()[12] = mu *(F.data()[4] * F.data()[0]);
	Hessian_2FE.data()[13] = mu * (2.0 * F.data()[1] * F.data()[4] + F.data()[0] * F.data()[3] + F.data()[2] * F.data()[5]);
	Hessian_2FE.data()[14] = mu *(F.data()[4] * F.data()[2]);
	Hessian_2FE.data()[15] = mu *(F.data()[7] * F.data()[0]);
	Hessian_2FE.data()[16] = mu * (2.0 * F.data()[1] * F.data()[7] + F.data()[0] * F.data()[6]  + F.data()[2] * F.data()[8]);
	Hessian_2FE.data()[17] = mu *(F.data()[7] * F.data()[2]);
	///
	Hessian_2FE.data()[20] = mu *(F.data()[0] * F.data()[0] + F.data()[1] * F.data()[1] + 3.0*F.data()[2] * F.data()[2] - 1.0 + F.data()[5] * F.data()[5] + F.data()[8] * F.data()[8]);
	Hessian_2FE.data()[21] = mu *(F.data()[5] * F.data()[0]);
	Hessian_2FE.data()[22] = mu *(F.data()[5] * F.data()[1]);
	Hessian_2FE.data()[23] = mu * (2.0 * F.data()[2] * F.data()[5] + F.data()[0] * F.data()[3] + F.data()[1] * F.data()[4]);
	Hessian_2FE.data()[24] = mu *(F.data()[8] * F.data()[0]);
	Hessian_2FE.data()[25] = mu *(F.data()[8] * F.data()[1]);
	Hessian_2FE.data()[26] = mu * (2.0 * F.data()[2] * F.data()[8] + F.data()[0] * F.data()[6] + F.data()[1] * F.data()[7]);
	///
	Hessian_2FE.data()[30] = mu *(F.data()[0] * F.data()[0] + 3.0* F.data()[3] * F.data()[3] + F.data()[4] * F.data()[4] - 1.0 + F.data()[5] * F.data()[5] + F.data()[6] * F.data()[6]);
	Hessian_2FE.data()[31] = mu * (2.0 * F.data()[3] * F.data()[4] + F.data()[0] * F.data()[1] + F.data()[6] * F.data()[7]);
	Hessian_2FE.data()[32] = mu * (2.0 * F.data()[3] * F.data()[5] + F.data()[0] * F.data()[2] + F.data()[6] * F.data()[8]);
	Hessian_2FE.data()[33] = mu * (2.0 * F.data()[3] * F.data()[6] + F.data()[4] * F.data()[7] + F.data()[5] * F.data()[8]);
	Hessian_2FE.data()[34] = mu *(F.data()[6] * F.data()[4]);
	Hessian_2FE.data()[35] = mu *(F.data()[6] * F.data()[5]);
	///
	Hessian_2FE.data()[40] = mu *(F.data()[1] * F.data()[1] + F.data()[3] * F.data()[3] + 3.0*F.data()[4] * F.data()[4] - 1.0 + F.data()[5] * F.data()[5] + F.data()[7] * F.data()[7]);
	Hessian_2FE.data()[41] = mu * (2.0 * F.data()[4] * F.data()[5] + F.data()[1] * F.data()[2] + F.data()[7] * F.data()[8]);
	Hessian_2FE.data()[42] = mu *(F.data()[7] * F.data()[3]);
	Hessian_2FE.data()[43] = mu * (2.0 * F.data()[4] * F.data()[7] + F.data()[3] * F.data()[6] + F.data()[5] * F.data()[8]);
	Hessian_2FE.data()[44] = mu *(F.data()[7] * F.data()[5]);
	///
	Hessian_2FE.data()[50] = mu *(F.data()[2] * F.data()[2] + F.data()[3] * F.data()[3] + F.data()[4] * F.data()[4] - 1.0 + 3.0*F.data()[5] * F.data()[5] + F.data()[8] * F.data()[8]);
	Hessian_2FE.data()[51] = mu *(F.data()[8] * F.data()[3]);
	Hessian_2FE.data()[52] = mu *(F.data()[8] * F.data()[4]);
	Hessian_2FE.data()[53] = mu * (2.0 * F.data()[5] * F.data()[8] + F.data()[3] * F.data()[6] + F.data()[4] * F.data()[7]);
	///
	Hessian_2FE.data()[60] = mu *(F.data()[0] * F.data()[0] + F.data()[3] * F.data()[3] + 3.0 * F.data()[6] * F.data()[6] - 1.0 + F.data()[7] * F.data()[7] + F.data()[8] * F.data()[8]);
	Hessian_2FE.data()[61] = mu * (2.0 * F.data()[6] * F.data()[7] + F.data()[0] * F.data()[1] + F.data()[4] * F.data()[3]);
	Hessian_2FE.data()[62] = mu * (2.0 * F.data()[6] * F.data()[8] + F.data()[0] * F.data()[2] + F.data()[5] * F.data()[3]);
	///
	Hessian_2FE.data()[70] = mu *(F.data()[1] * F.data()[1] + F.data()[4] * F.data()[4] + F.data()[6] * F.data()[6] - 1.0 + 3.0 * F.data()[7] * F.data()[7] + F.data()[8] * F.data()[8]);
	Hessian_2FE.data()[71] = mu * (2.0 * F.data()[8] * F.data()[7] + F.data()[2] * F.data()[1] + F.data()[4] * F.data()[5]);
	//
	Hessian_2FE.data()[80] = mu *(F.data()[2] * F.data()[2] + F.data()[5] * F.data()[5] + F.data()[6] * F.data()[6] - 1.0 + F.data()[7] * F.data()[7] + 3.0 * F.data()[8] * F.data()[8]);

	// add  first derivative of lambda*tr(E)F
	unsigned int k = 0;
	for (unsigned int i = 0; i < 9; ++i) {
		Hessian_2FE.data()[k] += lambda_trace;
		k += 10;
	}


	for (unsigned int i = 0; i < 9; ++i) {//col
		k = 9 * i;
		for (unsigned int j = i; j <9; ++j) {//row
			Hessian_2FE.data()[k+j] += lambda * F.data()[i] * F.data()[j];
		}
	}


	Hessian_2FE.data()[9] = Hessian_2FE.data()[1];
	Hessian_2FE.data()[18] = Hessian_2FE.data()[2];
	Hessian_2FE.data()[27] = Hessian_2FE.data()[3];
	Hessian_2FE.data()[36] = Hessian_2FE.data()[4];
	Hessian_2FE.data()[45] = Hessian_2FE.data()[5];
	Hessian_2FE.data()[54] = Hessian_2FE.data()[6];
	Hessian_2FE.data()[63] = Hessian_2FE.data()[7];
	Hessian_2FE.data()[72] = Hessian_2FE.data()[8];
	///
	Hessian_2FE.data()[19] = Hessian_2FE.data()[11];
	Hessian_2FE.data()[28] = Hessian_2FE.data()[12];
	Hessian_2FE.data()[37] = Hessian_2FE.data()[13];
	Hessian_2FE.data()[46] = Hessian_2FE.data()[14];
	Hessian_2FE.data()[55] = Hessian_2FE.data()[15];
	Hessian_2FE.data()[64] = Hessian_2FE.data()[16];
	Hessian_2FE.data()[73] = Hessian_2FE.data()[17];
	///
	Hessian_2FE.data()[29] = Hessian_2FE.data()[21];
	Hessian_2FE.data()[38] = Hessian_2FE.data()[22];
	Hessian_2FE.data()[47] = Hessian_2FE.data()[23];
	Hessian_2FE.data()[56] = Hessian_2FE.data()[24];
	Hessian_2FE.data()[65] = Hessian_2FE.data()[25];
	Hessian_2FE.data()[74] = Hessian_2FE.data()[26];
	///
	Hessian_2FE.data()[39] = Hessian_2FE.data()[31];
	Hessian_2FE.data()[48] = Hessian_2FE.data()[32];
	Hessian_2FE.data()[57] = Hessian_2FE.data()[33];
	Hessian_2FE.data()[66] = Hessian_2FE.data()[34];
	Hessian_2FE.data()[75] = Hessian_2FE.data()[35];
	///
	Hessian_2FE.data()[49] = Hessian_2FE.data()[41];
	Hessian_2FE.data()[58] = Hessian_2FE.data()[42];
	Hessian_2FE.data()[67] = Hessian_2FE.data()[43];
	Hessian_2FE.data()[76] = Hessian_2FE.data()[44];
	///
	Hessian_2FE.data()[59] = Hessian_2FE.data()[51];
	Hessian_2FE.data()[68] = Hessian_2FE.data()[52];
	Hessian_2FE.data()[77] = Hessian_2FE.data()[53];
	///
	Hessian_2FE.data()[69] = Hessian_2FE.data()[61];
	Hessian_2FE.data()[78] = Hessian_2FE.data()[62];
	///
	Hessian_2FE.data()[79] = Hessian_2FE.data()[71];


	Matrix<double, 9, 12> A;
	A.setZero();
	for (unsigned int i = 0; i < 12; i+=3) {//col
		for (unsigned int j = 0; j < 3; ++j) {
			A.data()[9*i +3* j] = inv_rest_pos.data()[i + j];
			A.data()[9 * i +3* j + 10] = inv_rest_pos.data()[i + j];
			A.data()[9 * i + 3*j + 20] = inv_rest_pos.data()[i + j];
		}
	}

	Hessian = A.transpose() * Hessian_2FE* A;


}



void XPBDconstraint::computeGreenStrainAndPiolaStress(
	double* v0, double* v1, double* v2, double* v3,
	Matrix<double, 3, 4>& inv_rest_pos,
	double rest_volume,
	double mu, double lambda, Matrix<double, 3, 4>& grad_C, double& C)
{
	Matrix3d F;
	Matrix3d P;
	SUB(P.data(), v1, v0);
	SUB((P.data() + 3), v2, v0);
	SUB((P.data() + 6), v3, v0);
	//first use eigen value to store the position of first vertex
	Matrix3d P_inv_transpose;
	memcpy(P_inv_transpose.data(), inv_rest_pos.data() + 3, 72);
	F = P * P_inv_transpose.transpose();

	if (F.determinant() < 0.58) {
		JacobiSVD<Matrix3d> svd;
		svd.compute(F, ComputeFullU | ComputeFullV);
		Vector3d new_singular_value = svd.singularValues();
		double min_s_value = 0.58;
		//once a critical compression
		//threshold is reached( 58 % of undeformed dimensions, when compression occurs
		//	along a single axis) the strength of the restorative force reaches a maximum.
		//Further compression will be met with decreasing resistance

		for (unsigned char j = 0; j < 3; j++) {
			if (new_singular_value[j] < min_s_value)
				new_singular_value[j] = min_s_value;
		}
		F = svd.matrixU() * new_singular_value.asDiagonal() * svd.matrixV().transpose();
	}

	Matrix3d E;

	//E= 0.5* (F^T F - I)
	E.data()[0] = 0.5 * (F.data()[0] * F.data()[0] + F.data()[1] * F.data()[1] + F.data()[2] * F.data()[2] - 1.0);		// xx
	E.data()[4] = 0.5 * (F.data()[3] * F.data()[3] + F.data()[4] * F.data()[4] + F.data()[5] * F.data()[5] - 1.0);		// yy
	E.data()[8] = 0.5 * (F.data()[6] * F.data()[6] + F.data()[7] * F.data()[7] + F.data()[8] * F.data()[8] - 1.0);		// zz
	E.data()[3] = 0.5 * (F.data()[0] * F.data()[3] + F.data()[1] * F.data()[4] + F.data()[2] * F.data()[5]);			// xy
	E.data()[6] = 0.5 * (F.data()[0] * F.data()[6] + F.data()[1] * F.data()[7] + F.data()[2] * F.data()[8]);			// xz
	E.data()[7] = 0.5 * (F.data()[3] * F.data()[6] + F.data()[4] * F.data()[7] + F.data()[5] * F.data()[8]);			// yz
	E.data()[1] = E.data()[3];
	E.data()[2] = E.data()[6];
	E.data()[5] = E.data()[7];

	//P(F)=F(2muE +lambda tr(E)I) 
	double trace = E.data()[0] + E.data()[4] + E.data()[8];
	double lambda_trace = lambda * trace;
	P = (mu + mu) * E;
	P.data()[0] += lambda_trace;
	P.data()[4] += lambda_trace;
	P.data()[8] += lambda_trace;
	P = F * P;
	//st. venant-kirchhoff model
	// psi= mu E:E + lambda/2 * tr(E)^2
	double psi = mu * E.squaredNorm() + 0.5 * lambda * trace * trace;
	C = rest_volume * psi;
	
	////Use E to temporary store P*P_inv_transpose*volume
	//E = P * P_inv_transpose * rest_volume;
	//memcpy((grad_C.data() + 3), E.data(), 72);
	//grad_C.col(0) = -grad_C.col(1) - grad_C.col(2) - grad_C.col(3);
	//Matrix<double, 3, 4> test;

	grad_C = rest_volume * P * inv_rest_pos;

}



void XPBDconstraint::solveEdgeLengthConstraint(double* p0, double* p1, const double rest_length, const double stiffness, double dt,
	double inv_mass0, double inv_mass1, double& lambda, const double damping_stiffness, double* initial_p0, double* inital_p1)
{	
	//double gamma = damping_stiffness / (stiffness*dt);
	dt = dt * dt;
	double n[3];
	SUB(n, p0, p1);
	double C = sqrt(DOT(n,n)) - rest_length;
	double n_norm = sqrt(DOT(n, n));
	DEV_(n, n_norm);

	double delta_lambda = -(C + lambda / (stiffness * dt))// + (gamma * n[0]*(p0[0] - initial_p0[0] - p1[0] + inital_p1[0] )+ gamma * n[1] * (p0[1] - initial_p0[1] - p1[1] + inital_p1[1])+ gamma * n[2] * (p0[2] - initial_p0[2] - p1[2] + inital_p1[2])))
		/ ((inv_mass0 + inv_mass1) + 1.0 / (stiffness * dt));//(1.0+gamma)*(inv_mass0 + inv_mass1) 
	lambda += delta_lambda;

	MULTI_(n, delta_lambda);
	ACCUMULATE_SUM_WITH_COE(p0, inv_mass0, n);
	inv_mass1 *= -1.0;
	ACCUMULATE_SUM_WITH_COE(p1, inv_mass1, n);

	//energy = 0.5 * stiffness * C * C;
}



void XPBDconstraint::solveBendingConstraint(double* center_vertex, double vertex_inv_mass,  std::array<double,3>* vertex_position, unsigned int* neighbor_vertex, unsigned int neighbor_vertex_size,
	double rest_Aq_norm, double lbo_weight, VectorXd& vertex_lbo,  double stiffness, double dt, double* inv_mass, double &lambda)
	//,	const double damping_stiffness, double* initial_center_vertex, std::array<double, 3>* inital_vertex_position
{
//	energy = 0.0;
	std::vector<VectorXd>q(3);
	//std::vector<VectorXd>q_initial(3);
	double aq[3];
	unsigned int size = neighbor_vertex_size + 1;
	VectorXd inv_m(size);
	for (unsigned int j = 0; j < 3; ++j) {
		q[j].resize(size);
		//q_initial[j].resize(size);
	}
	q[0][0] = center_vertex[0];
	q[1][0] = center_vertex[1];
	q[2][0] = center_vertex[2];

	//q_initial[0][0] = initial_center_vertex[0];
	//q_initial[1][0] = initial_center_vertex[1];
	//q_initial[2][0] = initial_center_vertex[2];

	inv_m[0] = vertex_inv_mass;


	for (unsigned int h = 1; h < size; h++) {
		q[0][h] = vertex_position[*neighbor_vertex][0];
		q[1][h] = vertex_position[*neighbor_vertex][1];
		q[2][h] = vertex_position[*neighbor_vertex][2];

		//q_initial[0][h] = inital_vertex_position[index][0];
		//q_initial[1][h] = inital_vertex_position[index][1];
		//q_initial[2][h] = inital_vertex_position[index][2];

		inv_m[h] = inv_mass[*neighbor_vertex];
		neighbor_vertex++;
	}
	neighbor_vertex -= neighbor_vertex_size;


	aq[0] = vertex_lbo.dot(q[0]);
	aq[1] = vertex_lbo.dot(q[1]);
	aq[2] = vertex_lbo.dot(q[2]);

	double aq_norm = sqrt(DOT(aq, aq));
	if (aq_norm < epsilon_for_bending) {
		return;
	}
	DEV_(aq, aq_norm);
	double C = aq_norm - rest_Aq_norm;
	//use q to store delta_C
	for (unsigned int j = 0; j < 3; ++j) {
		//q_initial[j] = q[j] - q_initial[j];
		q[j] = vertex_lbo * aq[j];
	}

	//aq[0] = vertex_lbo.dot(q[0])-rest_Aq[0];
	//aq[1] = vertex_lbo.dot(q[1]) - rest_Aq[1];
	//aq[2] = vertex_lbo.dot(q[2]) - rest_Aq[2];
	//double aq_norm = sqrt(DOT(aq, aq));
	//if (aq_norm < epsilon_for_bending) {
	//	return;
	//}
	//DEV_(aq, aq_norm);
	//double C = lbo_weight *aq_norm;
	//double C = 0.5 * lbo_weight * aq_norm;
	////use q to store delta_C
	//for (unsigned int j = 0; j < 3; ++j) {
	//	q_initial[j] = q[j] - q_initial[j];
	//	q[j] = vertex_lbo * (lbo_weight *aq[j]);
	//}
	//double alpha_ = 1.0 / (stiffness * dt * dt);
	//double gamma = damping_stiffness / (stiffness * dt);
	//double delta_lambda = -(C + alpha_ * lambda + gamma * (q[0].dot(q_initial[0]) + q[1].dot(q_initial[1]) + q[2].dot(q_initial[2])))
	//	/ ((1.0+gamma)*(q[0].dot(q[0].cwiseProduct(inv_m)) + q[1].dot(q[1].cwiseProduct(inv_m))
	//		+ q[2].dot(q[2].cwiseProduct(inv_m))) + alpha_);	
	//lambda += delta_lambda;
	//inv_m *= delta_lambda;

	double alpha_ = lbo_weight / (stiffness * dt * dt);
	//double gamma = lbo_weight * damping_stiffness / (stiffness * dt);
	double delta_lambda = -(C + alpha_ * lambda )//+ gamma * (q[0].dot(q_initial[0]) + q[1].dot(q_initial[1]) + q[2].dot(q_initial[2]))
		/ (((q[0].cwiseProduct(inv_m)).dot(q[0]) + (q[1].cwiseProduct(inv_m)).dot(q[1])//(1.0 + gamma) * (q[0].dot(q[0].cwiseProduct(inv_m))
			+ (q[2].cwiseProduct(inv_m)).dot(q[2])) + alpha_);
	lambda += delta_lambda;


	inv_m *= delta_lambda;
	center_vertex[0] += inv_m.data()[0] * q[0][0];
	center_vertex[1] += inv_m.data()[0] * q[1][0];
	center_vertex[2] += inv_m.data()[0] * q[2][0];

	for (unsigned int h = 1; h < size; h++) {
		vertex_position[*neighbor_vertex][0] += inv_m.data()[h] * q[0][h];
		vertex_position[*neighbor_vertex][1] += inv_m.data()[h] * q[1][h];
		vertex_position[*neighbor_vertex][2] += inv_m.data()[h] * q[2][h];
		neighbor_vertex++;
	}

	//energy = 0.5 * stiffness * C * C/ lbo_weight;
}

void XPBDconstraint::initial_LBO_EdgeCotWeight(TriangleMeshStruct& mesh_struct, std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo,
	std::vector<double>& rest_mean_curvature_norm)
{
	std::vector<double> edge_cot_weight;
	initialEdgeCotWeight(mesh_struct, edge_cot_weight);
	computeLBOWeight(lbo_weight, mesh_struct);
	computeVertexLBO(mesh_struct, vertex_lbo, edge_cot_weight);
	restBendingMeanCurvature(mesh_struct, rest_mean_curvature_norm, vertex_lbo, lbo_weight);
}

void XPBDconstraint::initialEdgeCotWeight(TriangleMeshStruct& mesh_struct, std::vector<double>& edge_cot_weight)
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
		
		edge_vertex_0 = mesh_struct.edge_vertices[i << 1];
		edge_vertex_1 = mesh_struct.edge_vertices[(i << 1) + 1];
		opposite_0 = mesh_struct.edges[i].opposite_vertex[0];
	
		cotan0 = 0;
		cotan1 = 0;
		SUB(x10, mesh_struct.vertex_position[edge_vertex_0], mesh_struct.vertex_position[opposite_0]);
		SUB(x20, mesh_struct.vertex_position[edge_vertex_1], mesh_struct.vertex_position[opposite_0]);
		len10 = sqrt(DOT(x10, x10));
		len20 = sqrt(DOT(x20, x20));
		theta0 = acos(DOT(x10, x20) / (len10 * len20));
		cotan0 = 1.0 / tan(theta0);

		if (mesh_struct.edges[i].opposite_vertex.size() > 1) {
			opposite_1 = mesh_struct.edges[i].opposite_vertex[1];

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
			edge_cot_weight[i] = -0.5 * cotan0;
		}
	}
}

void XPBDconstraint::computeLBOWeight(std::vector<double>& lbo_weight, TriangleMeshStruct& mesh_struct)
{
	double m;
	lbo_weight.resize(mesh_struct.vertex_position.size(),0.0);
	for (int i = 0; i < mesh_struct.faces.size(); ++i) {
		m = mesh_struct.faces[i].area / 3.0;
		lbo_weight[mesh_struct.triangle_indices[i][0]] += m;
		lbo_weight[mesh_struct.triangle_indices[i][1]] += m;
		lbo_weight[mesh_struct.triangle_indices[i][2]] += m;
	}

	//for (unsigned int i = 0; i < lbo_weight.size(); ++i) {
	//	lbo_weight[i] = 1.0 / sqrt(lbo_weight[i]);
	//}

}





void XPBDconstraint::computeVertexLBO(TriangleMeshStruct& mesh_struct, std::vector<VectorXd>& vertex_lbo, std::vector<double>& edge_cot_weight)
{
	double total;
	double edge_weight;
	vertex_lbo.resize(mesh_struct.vertex_position.size());
	for (int i = 0; i < mesh_struct.vertex_position.size(); ++i) {
		vertex_lbo[i].resize(mesh_struct.vertices[i].edge.size() + 1);
		vertex_lbo[i].setZero();
		total = 0.0;
		for (int j = 0; j < mesh_struct.vertices[i].edge.size(); ++j) {
			edge_weight = edge_cot_weight[mesh_struct.vertices[i].edge[j]];
			total += edge_weight;
			vertex_lbo[i].data()[j + 1] = edge_weight;
		}
		vertex_lbo[i].data()[0] = -total;
	}
}

void XPBDconstraint::restBendingMeanCurvature(TriangleMeshStruct& mesh_struct, std::vector<double>& rest_mean_curvature_norm,
//void XPBDconstraint::restBendingMeanCurvature(TriangleMeshStruct& mesh_struct, std::vector<Vector3d>& rest_Aq,
	std::vector<VectorXd>& vertex_lbo, std::vector<double>& lbo_weight)
{
	rest_mean_curvature_norm.resize(mesh_struct.vertex_position.size());
	//rest_Aq.resize(mesh_struct.vertex_position.size());
	VectorXd q;
	unsigned int size;
	unsigned int* neighbor_vertex;
	double curvature[3];
	for (unsigned int i = 0; i < rest_mean_curvature_norm.size(); ++i) {
		size = vertex_lbo[i].size();
		q.resize(size);
		neighbor_vertex = mesh_struct.vertices[i].neighbor_vertex.data();
		for (unsigned int j = 0; j < 3; ++j) {
			q[0] = mesh_struct.vertex_position[i][j];
			for (unsigned int h = 1; h < size; h++) {
				q[h] = mesh_struct.vertex_position[neighbor_vertex[h - 1]][j];
			}
			curvature[j] = dotProductX_(vertex_lbo[i], q);
		}
		rest_mean_curvature_norm[i] = sqrt(DOT(curvature, curvature));
	}

}