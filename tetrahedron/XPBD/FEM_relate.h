#pragma once
#include"../basic/global.h"
#include"../external/Eigen/Dense"
#include<array>
using namespace Eigen;

namespace FEM {
	//FEM course p31
	inline void getDeltaF(Matrix3d& delta_F, Matrix3d& Dm, unsigned int type)
	{
		delta_F.setZero();
		if (type < 3) {
			delta_F(type, 0) = -Dm.col(0).sum();
			delta_F(type, 1) = -Dm.col(1).sum();
			delta_F(type, 2) = -Dm.col(2).sum();
		}
		else if(type < 6) {
			delta_F.row(type - 3) = Dm.row(0);
		}
		else if (type < 9) {
			delta_F.row(type - 6) = Dm.row(1);
		}
		else {
			delta_F.row(type - 9) = Dm.row(2);
		}
	}
	//mpmcourse p20
	inline void getDeltaR(Matrix3d& delta_F,Matrix3d& delta_R, Matrix3d& S, Matrix3d& R)
	{
		//first compute R^T Delta R
		Matrix3d left = R.transpose() * delta_F - delta_F.transpose() * R;
		Matrix3d A;
		A << S.data()[0]+ S.data()[4], S.data()[7], -S.data()[6], 
			S.data()[7], S.data()[0]+ S.data()[8], S.data()[3],
			-S.data()[6], S.data()[3], S.data()[4]+ S.data()[8];
		Vector3d b;
		b << left.data()[3], left.data()[6], left.data()[7];
		ColPivHouseholderQR<Matrix3d> sys(A);
		Vector3d result=sys.solve(b);
		Matrix3d RT_Delta_R;
		RT_Delta_R << 0, result.data()[0], result.data()[1], 
			-result.data()[0], 0, result.data()[2],
			-result.data()[1], -result.data()[2], 0;
		delta_R = R * RT_Delta_R;
	}

	inline void getdPdF(Matrix3d& U, Matrix3d& V, Vector3d& eigen_value, Matrix<double, 9, 9>& dPdF)
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
	inline void backpropagateElementHessian(Matrix<double, 12, 12>& Hessian, Matrix<double, 9, 9>& dPdF, Matrix<double, 3, 4>& A)
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
	inline void getDeformationGradient(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
		Matrix<double, 3, 4>& A, Matrix3d& deformation_gradient)
	{
		Matrix3d q_e;
		double determinant;
		memcpy(q_e.data(), vertex_position_1, 24);
		memcpy(q_e.data() + 3, vertex_position_2, 24);
		memcpy(q_e.data() + 6, vertex_position_3, 24);
		//first use eigen value to store the position of first vertex
		for (unsigned int i = 0; i < 3; ++i) {
			q_e.data()[i] -= vertex_position_0[i];
			q_e.data()[i + 3] -= vertex_position_0[i];
			q_e.data()[i + 6] -= vertex_position_0[i];
		}

		Matrix3d P_inv;
		memcpy(P_inv.data(), A.data() + 3, 72);
		deformation_gradient = q_e * P_inv.transpose();
	}


	inline void polarDecomposition(Matrix3d& deformation_gradient, Vector3d& eigen_value, Matrix3d& S, Matrix3d& rotation)
	{
		JacobiSVD<Matrix3d> svd;
		svd.compute(deformation_gradient, ComputeFullU | ComputeFullV);
		eigen_value = svd.singularValues();
		double determinant = eigen_value[0] * eigen_value[1] * eigen_value[2];
		Matrix3d U;
		U = svd.matrixU();
		if (determinant < 0) {
			U.col(2) *= -1.0;
			eigen_value[2] *= -1.0;
		}
		rotation = U * svd.matrixV().transpose();
		S = svd.matrixV() * eigen_value.asDiagonal() * svd.matrixV().transpose();
	}

	inline void extractRotation(Matrix3d& deformation_gradient, Vector3d& eigen_value, Matrix3d& U, Matrix3d& V, Matrix3d& rotation)
	{
		JacobiSVD<Matrix3d> svd;
		svd.compute(deformation_gradient, ComputeFullU | ComputeFullV);
		eigen_value = svd.singularValues();
		double determinant = eigen_value[0] * eigen_value[1] * eigen_value[2];

		U = svd.matrixU();
		V = svd.matrixV();
		if (determinant < 0) {
			U.col(2) *= -1.0;
			eigen_value[2] *= -1.0;
		}
		rotation = U * V.transpose();
	}

	inline void getHessian(Matrix<double, 12, 12>& Hessian, Matrix3d& S, Matrix3d& R, Matrix3d& Dm, Matrix<double, 3, 4>& A)
	{
		Matrix3d delta_F, delta_R;
		Matrix<double, 3, 4> delta_H;
		for (unsigned int i = 0; i < 12; ++i) {
			getDeltaF(delta_F, Dm, i);
			getDeltaR(delta_F, delta_R, S, R);
			delta_H = 2.0 * (delta_F - delta_R) * A;
			memcpy(Hessian.data() + 12 * i, delta_H.data(), 96);
		}		
	}
}

