#pragma once
#include"XPBD/FEM_relate.h"
#include"basic/global.h"
#include"constitutiveModel.h"
#include <stdlib.h>

using namespace Eigen;

namespace NeoHookean {

	inline bool testDeterminent(double* position_0, double* position_1, double* position_2, double* position_3,
		Matrix<double, 3, 4>& A)
	{
		Matrix3d deformation_gradient;
		FEM::getDeformationGradient(position_0, position_1, position_2, position_3, A, deformation_gradient);

		double k = deformation_gradient.determinant();

		if (k <= 0) {
			std::cout << "error tet volume equals zero, error to compute ||F|| " << k << std::endl;
			k = 1e-36;
			return false;
			//system("pause");


		}
		return true;
	}

	inline double energy(double* position_0, double* position_1, double* position_2, double* position_3,
		Matrix<double, 3, 4>& A, double volume, double mu, double lambda, int tet_index, int obj_index)
	{
		Matrix3d deformation_gradient;
		FEM::getDeformationGradient(position_0, position_1, position_2, position_3, A, deformation_gradient);

		double k = deformation_gradient.determinant();

		if (k <= 0) {
			std::cout << "compute energy, error tet volume equals zero, error to compute ||F|| "<<k <<" "<< tet_index <<" "<< obj_index << std::endl;
			k = 1e-36;

			//system("pause");


		}


		double J = log(k);




		double energy = mu * 0.5 * (deformation_gradient.squaredNorm() - 3.0) - mu *J + lambda * 0.5 * J * J;
		return volume * energy;
	}


	inline void dpsi_dsigma(const Vector3d& sigma, double mu, double lambda, Vector3d& de_dsigma)
	{
		double log_sigma_prod = std::log(sigma.prod());
		double inv = 1.0 / sigma[0];
		de_dsigma[0] = mu * (sigma[0] - inv) + lambda * inv * log_sigma_prod;

		inv = 1.0 / sigma[1];
		de_dsigma[1] = mu * (sigma[1] - inv) + lambda * inv * log_sigma_prod;

		inv = 1.0 / sigma[2];
		de_dsigma[2] = mu * (sigma[2] - inv) + lambda * inv * log_sigma_prod;
	}


	inline void d2psi_dsigma2(const Vector3d& sigma, double mu, double lambda, Matrix3d& d2e_dsigma2)
	{
		double log_sigma_prod = std::log(sigma.prod());
		const double inv2_0 = 1.0 / (sigma[0] * sigma[0]);
		d2e_dsigma2(0, 0) = mu * (1.0 + inv2_0) - lambda * inv2_0 * (log_sigma_prod - 1.0);

		const double inv2_1 = 1.0 / (sigma[1] * sigma[1]);
		d2e_dsigma2(1, 1) = mu * (1.0 + inv2_1) - lambda * inv2_1 * (log_sigma_prod - 1.0);
		d2e_dsigma2(0, 1) = d2e_dsigma2(1, 0) = lambda / sigma[0] / sigma[1];

		const double inv2_2 = 1.0 / (sigma[2] * sigma[2]);
		d2e_dsigma2(2, 2) = mu * (1.0 + inv2_2) - lambda * inv2_2 * (log_sigma_prod - 1.0);
		d2e_dsigma2(1, 2) = d2e_dsigma2(2, 1) = lambda / sigma[1] / sigma[2];
		d2e_dsigma2(2, 0) = d2e_dsigma2(0, 2) = lambda / sigma[2] / sigma[0];

	}



	inline void B_left_coeff(const Vector3d& sigma, double mu, double lambda, Vector3d& left_coeff)
	{
		double sigma_prod = sigma.prod();

		if (sigma_prod <= 0) {
			std::cout << "error to compute B_left_coeff, product of eigen value equals zero " << std::endl;
			sigma_prod = 1e-36;
		}

		const double middle = mu - lambda * std::log(sigma_prod);
		left_coeff[0] = (mu + middle / (sigma[0] * sigma[1])) / 2.0;
		left_coeff[1] = (mu + middle / (sigma[1] * sigma[2])) / 2.0;
		left_coeff[2] = (mu + middle / (sigma[2] * sigma[0])) / 2.0;

	}


	//d||F|| /dF
	inline void cofactor(const Eigen::Matrix3d& A, Eigen::Matrix3d& coA)
	{
		coA(0, 0) = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
		coA(0, 1) = A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2);
		coA(0, 2) = A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0);
		coA(1, 0) = A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2);
		coA(1, 1) = A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0);
		coA(1, 2) = A(0, 1) * A(2, 0) - A(0, 0) * A(2, 1);
		coA(2, 0) = A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1);
		coA(2, 1) = A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2);
		coA(2, 2) = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
	}


	inline void firstPiola(const Matrix3d& F, double mu, double lambda, Matrix3d& P)
	{
		Matrix3d F_inv;
		cofactor(F, F_inv);
		double J = F.determinant();

		if (J <= 0) {
			std::cout << "error fist piola, tet volume equals zero, error to compute ||F|| " << std::endl;
			J = 1e-36;
		}

		F_inv /= J;
		P = mu * (F - F_inv) + lambda * std::log(J) * F_inv;
	}

	inline void firstPiolaDerivative(const Matrix3d& F, double mu, double lambda, MatrixXd& dPdF)
	{
		dPdF.setZero();
		JacobiSVD<Matrix3d> svd;
		Vector3d de_dsigma;
		Matrix3d d2e_dsigma2;

		svd.compute(F, ComputeFullU | ComputeFullV);



		dpsi_dsigma(svd.singularValues(), mu, lambda, de_dsigma);
		d2psi_dsigma2(svd.singularValues(), mu, lambda, d2e_dsigma2);

		Vector3d left_coeff;
		B_left_coeff(svd.singularValues(), mu, lambda, left_coeff);

		//std::cout << "de_dsigma " << de_dsigma.transpose() << std::endl;
		//std::cout << d2e_dsigma2 << std::endl;

		ConstitutiveModel::first_piola_derivative(svd.matrixU(), svd.singularValues(), svd.matrixV(), de_dsigma,
			left_coeff, d2e_dsigma2, dPdF);
	}

	inline void gradientMPM(double* position_0, double* position_1, double* position_2, double* position_3,
		Matrix<double, 3, 4>& A, double volume, double mu, double lambda, VectorXd& grad)
	{
		Matrix3d deformation_gradient;
		FEM::getDeformationGradient(position_0, position_1, position_2, position_3, A, deformation_gradient);

		Matrix3d P;
		firstPiola(deformation_gradient, mu, lambda, P);

		Matrix3d P_inv;
		memcpy(P_inv.data(), A.data() + 3, 72);
		P_inv.transposeInPlace();
		ConstitutiveModel::backpropagate_element_gradient(P_inv, P, grad);
		grad *= volume;
	}


	inline void gradientHessianMPM(double* position_0, double* position_1, double* position_2, double* position_3,
		Matrix<double, 3, 4>& A, double volume, double mu, double lambda, VectorXd& grad, MatrixXd& Hessian)
	{
		Matrix3d deformation_gradient;
		FEM::getDeformationGradient(position_0, position_1, position_2, position_3, A, deformation_gradient);

		Matrix3d P;
		firstPiola(deformation_gradient, mu, lambda, P);

		//std::cout << "p " << std::endl;
		//std::cout <<P << std::endl;

		Matrix3d P_inv;
		memcpy(P_inv.data(), A.data() + 3, 72);
		//std::cout << P_inv << std::endl;
		P_inv.transposeInPlace();
		ConstitutiveModel::backpropagate_element_gradient(P_inv, P, grad);
		grad *= volume;

		MatrixXd dPdF(9,9);
		firstPiolaDerivative(deformation_gradient, mu, lambda, dPdF);
		ConstitutiveModel::backpropagate_element_hessian(P_inv, dPdF, Hessian);
		Hessian *= volume;
	}


	//set cross product on a large matrix
	inline void setCP(MatrixXd& H, const double* f, bool positive, int start_row, int start_col)
	{
		if (positive) {
			H(start_row, start_col + 1) = -f[2];
			H(start_row, start_col + 2) = f[1];
			H(start_row + 1, start_col) = f[2];
			H(start_row + 1, start_col + 2) = -f[0];
			H(start_row + 2, start_col) = -f[1];
			H(start_row + 2, start_col + 1) = f[0];
		}
		else {
			H(start_row, start_col + 1) = f[2];
			H(start_row, start_col + 2) = -f[1];
			H(start_row + 1, start_col) = -f[2];
			H(start_row + 1, start_col + 2) = f[0];
			H(start_row + 2, start_col) = f[1];
			H(start_row + 2, start_col + 1) = -f[0];
		}
	}




	inline void d2JdF2(const Matrix3d& F, MatrixXd& H)
	{
		H.setZero();
		setCP(H, F.data() + 6, false, 0, 3);
		setCP(H, F.data() + 3, true, 0, 6);
		setCP(H, F.data() + 6, true, 3, 0);
		setCP(H, F.data(), false, 3, 6);
		setCP(H, F.data() + 3, false, 6, 0);
		setCP(H, F.data(), true, 6, 3);
	}


	inline void dJdF(VectorXd& dJ_dF, const Matrix3d& F)
	{
		dJ_dF.segment(0, 3) = F.col(1).cross(F.col(2));
		dJ_dF.segment(3, 3) = F.col(2).cross(F.col(0));
		dJ_dF.segment(6, 3) = F.col(0).cross(F.col(1));

	}

	inline void dFdX(MatrixXd& dF_dx, Matrix<double, 3, 4>& A)
	{

		const double m = A(0, 1); const double n = A(1, 1); const double o = A(2, 1);
		const double p = A(0, 2); const double q = A(1, 2); const double r = A(2, 2);
		const double s = A(0, 3); const double t = A(1, 3); const double u = A(2, 3);

		const double t1 = -m - p - s; const double t2 = -n - q - t; const double t3 = -o - r - u;


		dF_dx.setZero(); 
		dF_dx(0, 0) = t1;
		dF_dx(0, 3) = m;
		dF_dx(0, 6) = p;
		dF_dx(0, 9) = s;
		dF_dx(1, 1) = t1;
		dF_dx(1, 4) = m;
		dF_dx(1, 7) = p;
		dF_dx(1, 10) = s;
		dF_dx(2, 2) = t1;
		dF_dx(2, 5) = m;
		dF_dx(2, 8) = p;
		dF_dx(2, 11) = s;
		dF_dx(3, 0) = t2;
		dF_dx(3, 3) = n;
		dF_dx(3, 6) = q;
		dF_dx(3, 9) = t;
		dF_dx(4, 1) = t2;
		dF_dx(4, 4) = n;
		dF_dx(4, 7) = q;
		dF_dx(4, 10) = t;
		dF_dx(5, 2) = t2;
		dF_dx(5, 5) = n;
		dF_dx(5, 8) = q;
		dF_dx(5, 11) = t;
		dF_dx(6, 0) = t3;
		dF_dx(6, 3) = o;
		dF_dx(6, 6) = r;
		dF_dx(6, 9) = u;
		dF_dx(7, 1) = t3;
		dF_dx(7, 4) = o;
		dF_dx(7, 7) = r;
		dF_dx(7, 10) = u;
		dF_dx(8, 2) = t3;
		dF_dx(8, 5) = o;
		dF_dx(8, 8) = r;
		dF_dx(8, 11) = u;

	}

	inline void gradientHessian(double* position_0, double* position_1, double* position_2, double* position_3,
		Matrix<double, 3, 4>& A, double volume, double mu, double lambda, VectorXd& grad, MatrixXd& Hessian)
	{
		Matrix3d deformation_gradient;
		FEM::getDeformationGradient(position_0, position_1, position_2, position_3, A, deformation_gradient);

		double J =deformation_gradient.determinant();
		double log_J = log(J);

		double coe1 = (mu + lambda * (1.0 - log_J)) / (J * J);
		double coe2 = (lambda * log_J - mu) / J;

		VectorXd dJ_dF(9);
		dJdF(dJ_dF, deformation_gradient);

		MatrixXd d2J_dF2(9, 9);
		d2JdF2(deformation_gradient, d2J_dF2);

		d2J_dF2 *= coe2;
		d2J_dF2 += coe1 * dJ_dF * dJ_dF.transpose();
		for (int i = 0; i < 9; ++i) {
			d2J_dF2(i, i) += mu;
		}

		MatrixXd dF_dx(9, 12);
		dFdX(dF_dx, A);
		Hessian = dF_dx.transpose() * d2J_dF2 * dF_dx;
		

		dJ_dF *= coe2;		
		for (int i = 0; i < 9; ++i) {
			dJ_dF.data()[i] += mu * deformation_gradient.data()[i];
		}
		grad = dF_dx.transpose() * dJ_dF;

		FEM::SPDprojection(Hessian);
	}


}