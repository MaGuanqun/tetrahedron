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

	inline void getDeformationGradient(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
		Matrix<double, 3, 4>& A, Matrix3d& deformation_gradient)
	{
		Matrix3d q_e;
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



	inline void extractRotation(Matrix3d& deformation_gradient, Matrix3d& rotation)
	{
		JacobiSVD<Matrix3d> svd;
		svd.compute(deformation_gradient, ComputeFullU | ComputeFullV);
		if (deformation_gradient.determinant() > 0) {
			rotation = svd.matrixU() * svd.matrixV().transpose();
		}
		else {
			Matrix3d U;
			U = svd.matrixU();
			U.col(2) *= -1.0;
			rotation = U * svd.matrixV().transpose();
		}
	}

	inline void polarDecomposition(Matrix3d& deformation_gradient, Vector3d& eigen_value, Matrix3d& S, Matrix3d& rotation)
	{
		JacobiSVD<Matrix3d> svd;
		svd.compute(deformation_gradient, ComputeFullU | ComputeFullV);
		eigen_value = svd.singularValues();
		Matrix3d U;
		U = svd.matrixU();

		//std::cout << "determinant " << determinant<<" "<< << std::endl;
		if (deformation_gradient.determinant() < 0) {
		//if ((svd.matrixU()* svd.matrixV().transpose()).determinant() < 0) {
			//std::cout << "reverse " << std::endl;
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
		double determinant = deformation_gradient.determinant();

		U = svd.matrixU();
		V = svd.matrixV();
		if (determinant < 0) {
			U.col(2) *= -1.0;
			eigen_value[2] *= -1.0;
		}
		rotation = U * V.transpose();
	}

	inline void getHessianForOneVertex(Matrix3d& Hessian, Matrix3d& S, Matrix3d& R, Matrix3d& Dm, Matrix<double, 3, 4>& A, unsigned int vertex_no)
	{
		Matrix3d delta_F, delta_R;
		Matrix<double, 3, 4> delta_H;
		for (unsigned int i = 0; i < 3; ++i) {
			getDeltaF(delta_F, Dm, 3*vertex_no+ i);
			getDeltaR(delta_F, delta_R, S, R);
			delta_H = 2.0 * (delta_F - delta_R) * A;
			memcpy(Hessian.data() + 3 * i, delta_H.data()+3*vertex_no, 24);
		}
	}


	inline void getHessianForOneVertex(Matrix<double,12,3>& Hessian, Matrix3d& S, Matrix3d& R, Matrix3d& Dm, Matrix<double, 3, 4>& A, unsigned int vertex_no)
	{
		Matrix3d delta_F, delta_R;
		Matrix<double, 3, 4> delta_H;
		for (unsigned int i = 0; i < 3; ++i) {
			getDeltaF(delta_F, Dm, 3 * vertex_no + i);
			getDeltaR(delta_F, delta_R, S, R);
			delta_H = 2.0 * (delta_F - delta_R) * A;
			memcpy(Hessian.data() + 12 * i, delta_H.data(), 96);
		}
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


	inline void getHessianForSeveralVertex(MatrixXd& Hessian, Matrix3d& S, Matrix3d& R, Matrix3d& Dm, Matrix<double, 3, 4>& A, int* tet_vertex_order,
		unsigned int unfixed_vertex_num)
	{
		Matrix3d delta_F, delta_R;
		Matrix<double, 3, 4> delta_H;

		if (unfixed_vertex_num < 4) {
			for (int j = 0; j < unfixed_vertex_num; ++j) {
				for (unsigned int i = 0; i < 3; ++i) {
					getDeltaF(delta_F, Dm, 3 * tet_vertex_order[j] + i);
					getDeltaR(delta_F, delta_R, S, R);
					delta_H = 2.0 * (delta_F - delta_R) * A;
					for (int k = 0; k < unfixed_vertex_num; ++k) {
						memcpy(Hessian.data() + 3 * (unfixed_vertex_num * (3 * j + i) + k), delta_H.data() + 3 * tet_vertex_order[k], 24);
					}
				}
			}
		}
		else {
			for (unsigned int i = 0; i < 12; ++i) {
				getDeltaF(delta_F, Dm, i);
				getDeltaR(delta_F, delta_R, S, R);
				delta_H = 2.0 * (delta_F - delta_R) * A;
				memcpy(Hessian.data() + 12 * i, delta_H.data(), 96);
			}
		}
	}


	inline void SPDprojection(Eigen::Matrix<double,12,12>& A)
	{
	

		//JacobiSVD<Matrix<double, 12, 12>>  svd;
		//svd.compute(A, ComputeFullU | ComputeFullV);


		SelfAdjointEigenSolver<Matrix<double, 12, 12>> svd;
		if (svd.eigenvalues()[0] >= 0) {
			return;
		}

		//std::cout << "not SPD matrix occured " << std::endl;

		VectorXd fixed_eigen_value = svd.eigenvalues();
		for (unsigned int i = 0; i <12; ++i) {
			if (fixed_eigen_value.data()[i] < 0) {
				fixed_eigen_value.data()[i] = 0;
			}
		}
		A = svd.eigenvectors() * fixed_eigen_value.asDiagonal() * svd.eigenvectors().transpose();

	}

	inline void mssMatrix()
	{
		Matrix<double, 12, 12> mass;
		for (unsigned int i = 0; i < 12; ++i) {

		}
	}

}

