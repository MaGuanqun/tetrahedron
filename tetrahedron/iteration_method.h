#pragma once
#include"external/Eigen/Dense"
#include"basic/eigenDenseOperation.h"
#include"external/Eigen/Sparse"
#include"thread.h"


using namespace Eigen;
using namespace denseOperation;

class IterationMethod
{
public:

	void setConvergenceRate(double conv_rate, int max_jacobi_itr_num);

	void updateConvergenceRate(double conv_rate);

	void setBasicInfo(int object_num, std::vector<int>& sys_size, Thread* thread, std::vector<int>&cloth_per_thread_begin);
	void offDiagonalSize();
	void initialJacobi();
	void updateJacobi_R(int thread_id);
	void setOffDiagonal(int obj_No, std::vector<Triplet<double>>& global_mat_nnz);
	void solveByJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num);
	void solveBySuperJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num);
	void solveByChebyshevSemiIterativeSuperJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num);

	
	void estimateJacobiEigenValue(std::vector<std::vector<VectorXd>>& u);
	void solveByGaussSeidel(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num);
	void estimateSuperJacobiEigenValue(std::vector<std::vector<VectorXd>>& u);

	void estimateGaussSeidelEigenValue(std::vector<std::vector<VectorXd>>& u, std::vector<SparseMatrix<double, RowMajor>>& system_matrix);
	void solveByChebyshevGaussSeidel(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num, double weight);

	void solveByChebyshevSemiIterativeJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num);
	void solveByPCG(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num);
	void updateGlobalDiagonalInv();
	void initialGlobalDiagonalInv(std::vector<std::vector<double*>>* cloth_global_mat_diagonal_ref_address);
	void solveByWeightedJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num, double weight);


	std::vector<SparseMatrix<double, RowMajor>> off_diagonal;

	void test();
private:
	std::vector<SparseMatrix<double, RowMajor>> R_Jacobi;
	std::vector<VectorXd> global_diagonal_inv;


	std::vector<std::vector<double*>>* global_mat_diagonal_ref_address;//the address to quickly set the diagonal of system matrix A
	int total_object_num;
	std::vector<int>sys_size;
	Thread* thread;
	std::vector<int> obj_per_thread_begin;
	void updateJacobi(int cloth_No);
	int max_itr_num;
	double convergence_rate_2;
	void superJacobiSingleIteration(VectorXd& u, VectorXd& b, int cloth_No);
	int super_jacobi_step_size = 2;
	void estimateSuperJacobiEigenValue(int cloth_No, std::vector<VectorXd>& u);
	std::vector<double> super_jacobi_spectral_radius_square;
	std::vector<double> jacobi_spectral_radius_square;
	std::vector<double> gauss_seidel_spectral_radius_square;

	void estimateGaussSeidelEigenValue(int cloth_No, std::vector<VectorXd>& u, SparseMatrix<double, RowMajor>& system_matrix);

	void estimateJacobiEigenValue(int cloth_No, std::vector<VectorXd>& u);

	template <class T>
	void superJacobiSingleIteration(Matrix<T, -1, 1>& u, Matrix<T, -1, 1>& b, Matrix<T, -1, 1>& global_diagonal_inv, SparseMatrix<T, RowMajor>& R_Jacobi)
	{
		Matrix<T, -1, 1> R_D_inv_b = b.cwiseProduct(global_diagonal_inv);
		Matrix<T, -1, 1> R_x = R_Jacobi * u;
		Matrix<T, -1, 1> sum_R_b = R_D_inv_b;
		for (int i = 1; i < super_jacobi_step_size; ++i) {
			R_D_inv_b = R_Jacobi * R_D_inv_b;
			sum_R_b += R_D_inv_b;
			R_x = R_Jacobi * R_x;
		}
		u = sum_R_b + R_x;
	}
	template <class T>
	void prepareR(SparseMatrix<T, RowMajor>& R_Jacobi, Matrix<T, -1, 1>& global_diagonal_inverse, SparseMatrix<T, RowMajor>& system_matrix)
	{
		int size = system_matrix.cols();
		R_Jacobi = system_matrix;
		global_diagonal_inverse.resize(size);

		for (int i = 0; i < size; ++i) {
			R_Jacobi.coeffRef(i, i) = 0.0;
			global_diagonal_inverse[i] = 1.0 / system_matrix.coeff(i, i);
		}
		R_Jacobi.prune((T)0.0);
		R_Jacobi *= -1.0;

		for (int j = 0; j < R_Jacobi.cols(); ++j) {
			for (int k = R_Jacobi.outerIndexPtr()[j]; k < R_Jacobi.outerIndexPtr()[j + 1]; ++k) {
				R_Jacobi.valuePtr()[k] *= global_diagonal_inverse.data()[j];
			}
		}
	}


	template <class T>
	void solveByChebyshevSemiIterativeSuperJacobi(Matrix<T, -1, 1>& u, Matrix<T, -1, 1>& b, SparseMatrix<T, RowMajor>& system_matrix, SparseMatrix<T, RowMajor>& R_Jacobi, Matrix<T, -1, 1>& global_diagonal_inv, double super_jacobi_spectral_radius_square, int& itr_num, Matrix<T, -1, 1>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		Matrix<T, -1, 1> u_last = u;
		Matrix<T, -1, 1> u_previous;
		double b_norm = b.squaredNorm();
		double omega_chebyshev = 2.0;
		double residual_norm;
		relative_error.reserve(max_itr_num);
		relative_error.push_back((u - ground_truth).squaredNorm());
		superJacobiSingleIteration(u, b, global_diagonal_inv, R_Jacobi);
		itr_num = 1;
		residual_norm = (b - system_matrix * u).squaredNorm();
		relative_error.push_back((u - ground_truth).squaredNorm());
		while (residual_norm / b_norm > conv_rate_2) {
			u_previous = u;
			omega_chebyshev = 4.0 / (4.0 - super_jacobi_spectral_radius_square * omega_chebyshev);
			superJacobiSingleIteration(u, b, global_diagonal_inv, R_Jacobi);
			u = omega_chebyshev * (u - u_last) + u_last;
			u_last = u_previous;
			itr_num++;
			residual_norm = (b - system_matrix * u).squaredNorm();
			relative_error.push_back((u - ground_truth).squaredNorm());
		}
	}

	template <class T>
	void solveBySuperJacobi(Matrix<T, -1, 1>& u, Matrix<T, -1, 1>& b, SparseMatrix<T, RowMajor>& system_matrix, SparseMatrix<T, RowMajor>& R_Jacobi, Matrix<T, -1, 1>& global_diagonal_inv, int& itr_num, Matrix<T, -1, 1>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		double residual_norm;
		double b_norm = b.squaredNorm();
		residual_norm = 2.0 * b_norm;
		itr_num = 0;
		relative_error.reserve(max_itr_num);
		relative_error.push_back((u - ground_truth).squaredNorm());
		while (residual_norm / b_norm > conv_rate_2) {
			superJacobiSingleIteration(u, b, global_diagonal_inv, R_Jacobi);
			residual_norm = (b - system_matrix * u).squaredNorm();
			itr_num++;
			relative_error.push_back((u - ground_truth).squaredNorm());
		}
	}

	template <class T>
	void solveByJacobi(Matrix<T, -1, 1>& u, Matrix<T, -1, 1>& b, SparseMatrix<T, RowMajor>& system_matrix, SparseMatrix<T, RowMajor>& R_Jacobi, Matrix<T, -1, 1>& global_diagonal_inv, int& itr_num, Matrix<T, -1, 1>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		double residual_norm;
		double b_norm = b.squaredNorm();
		residual_norm = 2.0 * b_norm;
		itr_num = 0;
		relative_error.reserve(max_itr_num);
		relative_error.push_back((u - ground_truth).squaredNorm());
		while (residual_norm / b_norm > conv_rate_2) {
			//u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + RMultiplyX(u, cloth_No, 0);
			u = b.cwiseProduct(global_diagonal_inv) + R_Jacobi * u;
			residual_norm = (b - system_matrix * u).squaredNorm();
			relative_error.push_back((u - ground_truth).squaredNorm());
			itr_num++;
		}
	}
	template <class T>
	void solveByGaussSeidel(Matrix<T, -1, 1>& u, Matrix<T, -1, 1>& b, SparseMatrix<T, RowMajor>& system_matrix, int& itr_num, Matrix<T, -1, 1>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		double residual_norm;
		double b_norm = b.squaredNorm();
		residual_norm = 2.0 * b_norm;
		itr_num = 0;
		relative_error.reserve(max_itr_num);
		relative_error.push_back((u - ground_truth).squaredNorm());
		while (residual_norm / b_norm > conv_rate_2) {
			u = b - system_matrix.triangularView<StrictlyUpper>() * u;
			u = system_matrix.triangularView<Lower>().solve(u);
			residual_norm = (b - system_matrix * u).squaredNorm();
			relative_error.push_back((u - ground_truth).squaredNorm());
			itr_num++;
		}
	}

	//weighted chebyshev gauss seidel, guarantee convergence
	template <class T>
	void solveByChebyshevGaussSeidel(Matrix<T, -1, 1>& u, Matrix<T, -1, 1>& b, SparseMatrix<T, RowMajor>& system_matrix,
		double gauss_seidel_spectral_radius_square, int& itr_num, double weight, Matrix<T, -1, 1>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		Matrix<T, -1, 1> u_last = u;
		double residual_norm;
		double omega_chebyshev = 2.0;
		double b_norm = b.squaredNorm();
		relative_error.reserve(max_itr_num);
		relative_error.push_back((u - ground_truth).squaredNorm());
		u = b - system_matrix.triangularView<StrictlyUpper>() * u;
		u = system_matrix.triangularView<Lower>().solve(u);
		itr_num = 1;
		residual_norm = (b - system_matrix * u).squaredNorm();
		relative_error.push_back((u - ground_truth).squaredNorm());
		Matrix<T, -1, 1> u_previous;
		while (residual_norm / b_norm > conv_rate_2) {
			u_previous = u;
			omega_chebyshev = 4.0 / (4.0 - gauss_seidel_spectral_radius_square * omega_chebyshev);
			u = system_matrix.triangularView<Lower>().solve(b - system_matrix.triangularView<StrictlyUpper>() * u);
			//u = omega_chebyshev * (u - u_last) + u_last;
			u = omega_chebyshev * (weight * (u - u_previous) + u_previous - u_last) + u_last;
			u_last = u_previous;
			residual_norm = (b - system_matrix * u).squaredNorm();
			itr_num++;
			relative_error.push_back((u - ground_truth).squaredNorm());
		}
	}
	template <class T>
	void solveByPCG(Matrix<T, -1, 1>& u, Matrix<T, -1, 1>& b, SparseMatrix<T, RowMajor>& system_matrix, int& itr_num,
		Matrix<T, -1, 1>& global_diagonal_inv, Matrix<T, -1, 1>& ground_truth, std::vector<double>& relative_error, double conv_rate_2)
	{
		Matrix<T, -1, 1> residual;
		double residual_norm;
		double b_norm = b.squaredNorm();
		Matrix<T, -1, 1> z, p;
		residual = b - system_matrix * u;
		z = global_diagonal_inv.cwiseProduct(residual);
		p = z;
		double alpha, beta;
		double rz_k, rz_k_1;
		rz_k = residual.dot(z);
		itr_num = 0;
		relative_error.reserve(max_itr_num);
		relative_error.push_back((u - ground_truth).squaredNorm());

		while (true)
		{
			itr_num++;
			alpha = rz_k / (system_matrix * p).dot(p);
			u += alpha * p;
			residual = b - system_matrix * u;
			relative_error.push_back((u - ground_truth).squaredNorm());
			if (residual.squaredNorm() / b_norm < conv_rate_2) {
				break;
			}
			z = global_diagonal_inv.cwiseProduct(residual);
			rz_k_1 = residual.dot(z);
			beta = rz_k_1 / rz_k;
			rz_k = rz_k_1;
			p = z + beta * p;

		}
	}
	template <class T>
	void solveByChebyshevSemiIterativeJacobi(Matrix<T, -1, 1>& u, Matrix<T, -1, 1>& b, SparseMatrix<T, RowMajor>& system_matrix,
		SparseMatrix<T, RowMajor>& R_Jacobi, int& itr_num, Matrix<T, -1, 1>& global_diagonal_inv, Matrix<T, -1, 1>& ground_truth, std::vector<double>& relative_error,
		double jacobi_spectral_radius_square, double conv_rate_2)
	{
		Matrix<T, -1, 1> u_last = u;
		Matrix<T, -1, 1> B_inv_b = global_diagonal_inv.cwiseProduct(b);
		//u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + RMultiplyX(u, cloth_No, 0);
		relative_error.reserve(max_itr_num);
		relative_error.push_back((u - ground_truth).squaredNorm());
		u = b.cwiseProduct(global_diagonal_inv) + R_Jacobi * u;
		itr_num = 1;
		double omega_chebyshev = 2.0;
		Matrix<T, -1, 1> u_previous;
		double b_norm = b.squaredNorm();
		double residual_chebyshev = 2 * b_norm;
		relative_error.push_back((u - ground_truth).squaredNorm());
		while (residual_chebyshev / b_norm > conv_rate_2)
		{
			u_previous = u;
			omega_chebyshev = 4.0 / (4.0 - jacobi_spectral_radius_square * omega_chebyshev);
			u = b.cwiseProduct(global_diagonal_inv) + R_Jacobi * u;
			u = omega_chebyshev * (u - u_last) + u_last;
			u_last = u_previous;
			itr_num++;
			residual_chebyshev = (b - system_matrix * u).squaredNorm();
			relative_error.push_back((u - ground_truth).squaredNorm());
		}
	}
	template <class T>
	double estimateGaussSeidelEigenValue(std::vector<Matrix<T, -1, 1>>& u, SparseMatrix<T, RowMajor>& system_matrix)
	{
		double vec_norm2 = 0;
		double u_norm2 = 0;
		for (int i = 0; i < 3; ++i) {
			vec_norm2 += system_matrix.triangularView<Lower>().solve(system_matrix.triangularView<StrictlyUpper>() * u[i]).squaredNorm();
			u_norm2 += u[i].squaredNorm();
		}
		return vec_norm2 / u_norm2;
	}

	//void jacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, VectorXd& ground_truth, std::vector<double>& relative_error);
	template <class T>
	void jacobi(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, RowMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		SparseMatrix<T, RowMajor> R_Jacobi;
		Matrix<T, -1, 1> global_diagonal_inverse;
		prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
		std::cout << R_Jacobi.cols() << " " << R_Jacobi.rows() << std::endl;

		int itr_num;
		std::vector<std::vector<double>> relat(3);
		for (int i = 0; i < 3; ++i) {
			solveByJacobi(u[i], b[i], system_matrix, R_Jacobi, global_diagonal_inverse, itr_num, ground_truth[i], relat[i],conv_rate_2);
		}
		computeRelativeError(relat, relative_error, ground_truth);

	}

	template <class T>
	void computeRelativeError(std::vector<std::vector<double>>& relat, std::vector<double>& relative_error, std::vector<Matrix<T, -1, 1>>& ground_truth)
	{
		double ground_truth_norm;
		//ground_truth_norm = ground_truth[0].squaredNorm() + ground_truth[1].squaredNorm() + ground_truth[2].squaredNorm();
		ground_truth_norm = ground_truth[1].squaredNorm();
		relative_error.resize(relat[1].size());
		for (int i = 0; i < relative_error.size(); ++i) {
			//relative_error[i] = sqrt((relat[0][i] + relat[1][i] + relat[2][i]) / ground_truth_norm);
			relative_error[i] = sqrt(relat[1][i] / ground_truth_norm);
		}
	}

	template <class T>
	void superJacobi(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, RowMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		SparseMatrix<T, RowMajor> R_Jacobi;
		Matrix<T, -1, 1> global_diagonal_inverse;
		std::vector<std::vector<double>> relat(3);
		computeRelativeError(relat, relative_error, ground_truth);
		prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
		int itr_num;
		for (int i = 0; i < 3; ++i) {
			solveBySuperJacobi(u[i], b[i], system_matrix, R_Jacobi, global_diagonal_inverse, itr_num, ground_truth[i], relat[i],conv_rate_2);
		}
		computeRelativeError(relat, relative_error, ground_truth);
	}
	template <class T>
	void gauss_seidel(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, RowMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		int itr_num;
		std::vector<std::vector<double>> relat(3);
		for (int i = 0; i < 3; ++i) {
			solveByGaussSeidel(u[i], b[i], system_matrix, itr_num, ground_truth[i], relat[i],conv_rate_2);
		}
		computeRelativeError(relat, relative_error, ground_truth);
	}


	template <class T>
	void chebyshevSemiIterativeJacobi(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, RowMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		SparseMatrix<T, RowMajor> R_Jacobi;
		Matrix<T, -1, 1> global_diagonal_inverse;
		std::vector<std::vector<double>> relat(3);
		prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
		int itr_num;
		double eigen_value_square = estimateSuperJacobiEigenValue(u, R_Jacobi, 1);
		std::cout << "eigen value " << eigen_value_square << std::endl;
		for (int i = 0; i < 3; ++i) {
			solveByChebyshevSemiIterativeJacobi(u[i], b[i], system_matrix, R_Jacobi, itr_num, global_diagonal_inverse, ground_truth[i], relat[i],
				eigen_value_square,conv_rate_2);
		}
		computeRelativeError(relat, relative_error, ground_truth);
	}

	template <class T>
	double estimateSuperJacobiEigenValue(std::vector<Matrix<T, -1, 1>>& u, SparseMatrix<T, RowMajor>& R_Jacobi, int super_jacobi_step_size)
	{
		Matrix<T, -1, 1> Vector_change;
		double vec_norm2 = 0;
		double u_norm2 = 0;
		for (int j = 0; j < 3; ++j) {
			Vector_change = R_Jacobi * u[j];
			for (int i = 1; i < super_jacobi_step_size; ++i) {
				Vector_change = R_Jacobi * Vector_change;
			}
			vec_norm2 += Vector_change.squaredNorm();
			u_norm2 += u[j].squaredNorm();
		}
		return vec_norm2 / u_norm2;
	}

	template <class T>
	void chebyshevSemiIterativeSuperJacobi(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, RowMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		SparseMatrix<T, RowMajor> R_Jacobi;
		Matrix<T, -1, 1> global_diagonal_inverse;
		prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
		int itr_num;
		double eigen_value_square = estimateSuperJacobiEigenValue(u, R_Jacobi, 2);
		std::vector<std::vector<double>> relat(3);
		for (int i = 0; i < 3; ++i) {
			solveByChebyshevSemiIterativeSuperJacobi(u[i], b[i], system_matrix, R_Jacobi, global_diagonal_inverse, eigen_value_square, itr_num,
				ground_truth[i], relat[i],conv_rate_2);
		}
		computeRelativeError(relat, relative_error, ground_truth);
	}



	
	template <class T>
	void chebyshev_gauss_seidel(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, RowMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		int itr_num;
		std::vector<std::vector<double>> relat(3);
		double eigen_value_square = estimateGaussSeidelEigenValue(u, system_matrix);
		for (int i = 0; i < 3; ++i) {
			solveByChebyshevGaussSeidel(u[i], b[i], system_matrix, eigen_value_square, itr_num, 0.6, ground_truth[i], relat[i],conv_rate_2);
		}
		computeRelativeError(relat, relative_error, ground_truth);
	}
	template <class T>
	void PCG(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, RowMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2)
	{
		int itr_num;
		SparseMatrix<T, RowMajor> R_Jacobi;
		Matrix<T, -1, 1> global_diagonal_inverse;
		std::vector<std::vector<double>> relat(3);
		prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
		for (int i = 0; i < 3; ++i) {
			solveByPCG(u[i], b[i], system_matrix, itr_num, global_diagonal_inverse, ground_truth[i], relat[i],conv_rate_2);
		}
		computeRelativeError(relat, relative_error, ground_truth);
	}

	//std::vector<int>dimension_per_thread_test;
};

