#pragma once

#include"external/Eigen/Dense"
#include"basic/eigenDenseOperation.h"
#include"basic/global.h"
#include"XPBD/FEM_relate.h"

using namespace Eigen;

namespace TEST_HESSIAN {


	inline double computeC(Vector3d& v0, Vector3d& v1, Vector3d& v2, Vector3d& v3, Matrix<double, 3, 4>& A)
	{
		Matrix3d deformation_gradient;
		FEM::getDeformationGradient(v0.data(), v1.data(), v2.data(), v3.data(), A, deformation_gradient);
		JacobiSVD<Matrix3d> svd;
		svd.compute(deformation_gradient);
		double norm = (svd.singularValues()[0] - 1.0) * (svd.singularValues()[0] - 1.0) + (svd.singularValues()[1] - 1.0) * (svd.singularValues()[1] - 1.0) +
			(svd.singularValues()[2] - 1.0) * (svd.singularValues()[2] - 1.0);
		return sqrt(norm);
	}

	inline Matrix<double, 12, 1> computeGrad(Vector3d& v0, Vector3d& v1, Vector3d& v2, Vector3d& v3, Matrix<double, 3, 4>& A)
	{
		Matrix3d deformation_gradient;
		FEM::getDeformationGradient(v0.data(), v1.data(), v2.data(), v3.data(), A, deformation_gradient);
		Matrix3d U, V, rotation;
		Vector3d eigen_value;
		FEM::extractRotation(deformation_gradient, eigen_value, U, V, rotation);
		double C = (deformation_gradient - rotation).norm();
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (1.0 / C) * (deformation_gradient - rotation) * A;//
		Matrix<double, 12, 1> grad;
		memcpy(grad.data(), grad_C_transpose.data(), 96);
		return grad;
	}

	inline Matrix<double, 12, 12> computeHessianByGradNumeric(std::vector<Vector3d>& v, Matrix<double, 3, 4>& A, double step_size)
	{
		std::vector<Vector3d>temp_v = v;
		Matrix<double, 12, 1>grad,grad_forward;
		Matrix<double, 12, 12> Hessian;
		grad = computeGrad(v[0], v[1], v[2], v[3], A);
		for (unsigned int i = 0; i < 4; ++i) {
			for (unsigned int j = 0; j < 3; ++j) {
				temp_v = v;
				temp_v[i][j] += step_size;
				grad_forward= computeGrad(temp_v[0], temp_v[1], temp_v[2], temp_v[3], A);
				Hessian.col(3 * i + j) = (grad_forward - grad) / step_size;
			}
		}
		return Hessian;
	}

	inline Matrix<double, 12, 12> computeHessianByAna(std::vector<Vector3d>& v, Matrix<double, 3, 4>& A, double step_size)
	{
		Matrix3d deformation_gradient;
		FEM::getDeformationGradient(v[0].data(), v[1].data(), v[2].data(), v[3].data(), A, deformation_gradient);
		Matrix3d U, V, rotation;
		Vector3d eigen_value;
		FEM::extractRotation(deformation_gradient, eigen_value, U, V, rotation);
		double C = computeC(v[0], v[1], v[2], v[3], A);
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (1.0 / C) * (deformation_gradient - rotation) * A;//
		Matrix<double, 12, 1> grad;
		memcpy(grad.data(), grad_C_transpose.data(), 96);
		Matrix<double, 9, 9> dPdF;
		//FEM::getdPdF(U, V, eigen_value, dPdF);
		dPdF = 2.0 * Matrix<double, 9, 9>::Identity();
		Matrix<double, 12, 12> Hessian;
		FEM::backpropagateElementHessian(Hessian, dPdF, A);
		Hessian *= (0.5 / C);
		Hessian -= ((1.0 / C) * grad) * grad.transpose();
		return Hessian;
	}

	inline Matrix<double, 12, 12> computeHessianByNumeric(std::vector<Vector3d>& v, Matrix<double, 3, 4>& A, double step_size)
	{
		std::vector<Vector3d>temp_v = v;
		std::vector<Vector3d> grad_c_previous;
		double C[4], C_forward[4][3], C_previous[4][3];
		for (unsigned int i = 0; i < 4; ++i) {
			C[i] = computeC(v[0],v[1],v[2],v[3], A);
			for (unsigned int j = 0; j < 3; ++j) {
				temp_v = v;
				temp_v[i][j] += step_size;
				double C_ = computeC(temp_v[0], temp_v[1], temp_v[2], temp_v[3], A);
				C_forward[i][ j] = (C_ - C[i]) / step_size;
				temp_v = v;
				temp_v[i][j] -= step_size;
				C_ = computeC(temp_v[0], temp_v[1], temp_v[2], temp_v[3], A);
				C_previous[i][j] = (C[i]-C_) / step_size;
			}			
		}
		
		Matrix<double, 12, 12> hessian;

		//for (unsigned int i = 0; i < 4; ++i) {
		//	for (unsigned int j = 0; j < 4; ++j) {
		//		for (unsigned int m = 0; m < 3; ++m) {
		//			for (unsigned int n = 0; n < 3; ++n) {
		//				hessian(3*i+m, 3*j+n) = C_forward-C_previous
		//			}
		//		}
		//	}
		//}

	}



	inline void testARAPHessian()
	{
		Vector3d x0, x1, x2, x3;
		Vector3d ori_x0, ori_x1, ori_x2, ori_x3;

		double info_v0[3], info_v1[3], info_v2[3], info_v3[3];


		info_v0[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
		info_v0[1] = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
		info_v0[2] = 1.0;
		info_v1[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
		info_v1[1] = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
		info_v1[2] = 1.0;
		info_v2[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
		info_v2[1] = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
		info_v2[2] = 1.0;
		info_v3[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
		info_v3[1] = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
		info_v3[2] = 1.0;

		ori_x0 = Vector3d(info_v0[2] * sin(info_v0[0]), info_v0[2] * cos(info_v0[0]) * sin(info_v0[1]), info_v0[2] * cos(info_v0[0]) * cos(info_v0[1]));
		ori_x1 = Vector3d(info_v1[2] * sin(info_v1[0]), info_v1[2] * cos(info_v1[0]) * sin(info_v1[1]), info_v1[2] * cos(info_v1[0]) * cos(info_v1[1]));
		ori_x2 = Vector3d(info_v2[2] * sin(info_v2[0]), info_v2[2] * cos(info_v2[0]) * sin(info_v2[1]), info_v2[2] * cos(info_v2[0]) * cos(info_v2[1]));
		ori_x3 = Vector3d(info_v3[2] * sin(info_v3[0]), info_v3[2] * cos(info_v3[0]) * sin(info_v3[1]), info_v3[2] * cos(info_v3[0]) * cos(info_v3[1]));

		if ((ori_x3 - ori_x0).dot((ori_x1 - ori_x0).cross(ori_x2 - ori_x0)) < 0) {
			Vector3d k = ori_x2;
			ori_x2 = ori_x1;
			ori_x1 = k;
		}

		if ((ori_x3 - ori_x0).dot((ori_x1 - ori_x0).cross(ori_x2 - ori_x0)) < 1e-8) {
			return;
		}

		info_v0[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
		info_v0[1] = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
		info_v0[2] = 4.0;
		info_v1[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
		info_v1[1] = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
		info_v1[2] = 4.0;
		info_v2[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
		info_v2[1] = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
		info_v2[2] = 4.0;
		info_v3[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
		info_v3[1] = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
		info_v3[2] = 4.0;

		x0 = Vector3d(info_v0[2] * sin(info_v0[0]), info_v0[2] * cos(info_v0[0]) * sin(info_v0[1]), info_v0[2] * cos(info_v0[0]) * cos(info_v0[1]));
		x1 = Vector3d(info_v1[2] * sin(info_v1[0]), info_v1[2] * cos(info_v1[0]) * sin(info_v1[1]), info_v1[2] * cos(info_v1[0]) * cos(info_v1[1]));
		x2 = Vector3d(info_v2[2] * sin(info_v2[0]), info_v2[2] * cos(info_v2[0]) * sin(info_v2[1]), info_v2[2] * cos(info_v2[0]) * cos(info_v2[1]));
		x3 = Vector3d(info_v3[2] * sin(info_v3[0]), info_v3[2] * cos(info_v3[0]) * sin(info_v3[1]), info_v3[2] * cos(info_v3[0]) * cos(info_v3[1]));

		if ((x3 - x0).dot((x1 - x0).cross(x2 - x0)) < 0) {
			Vector3d k = x2;
			x2 = x1;
			x1 = k;
		}


		Matrix<double, 3, 3> p;
		p.col(0) = ori_x1 - ori_x0;
		p.col(1) = ori_x2 - ori_x0;
		p.col(2) = ori_x3 - ori_x0;
		Matrix3d p_ = p.inverse();
		Matrix<double, 3, 4> A;
		for (unsigned int i = 0; i < 3; ++i) {
			A.data()[i] = -p_.col(i).sum();
			A.data()[3 + i] = p_.data()[3 * i];
			A.data()[6 + i] = p_.data()[3 * i + 1];
			A.data()[9 + i] = p_.data()[3 * i + 2];
		}
		std::vector<Vector3d> v(4);
		v[0] = x0; v[1] = x1; v[2] = x2; v[3] = x3;
		double step_size = 1e-5;
		Matrix<double, 12, 12> Hessian_num = computeHessianByGradNumeric(v, A, step_size);
		Matrix<double, 12, 12> Hessian_ana = computeHessianByAna(v, A, step_size);

		std::cout << (Hessian_ana - Hessian_num).norm() << std::endl;

		std::cout << Hessian_ana << std::endl;
		std::cout << "===" << std::endl;
		std::cout << Hessian_num << std::endl;

	}

	inline void testARAPHessianMulti()
	{
		for (unsigned int i = 0; i < 1; ++i) {
			testARAPHessian();
		}
	}

}

