#pragma once
#include"external/Eigen/Dense"
#include"basic/eigenDenseOperation.h"
#include"external/Eigen/Sparse"
#include"thread.h"
#include"basic/global.h"

using namespace Eigen;
using namespace denseOperation;


class ProjectDynamic;

class IterationMethod
{
public:
	void setConvergenceRate(double conv_rate, int max_jacobi_itr_num);

	void updateConvergenceRate(double conv_rate);

	void setBasicInfo(int object_num, std::vector<int>& sys_size, Thread* thread, std::vector<int>& cloth_per_thread_begin);
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

	void RMultiXPlusb(std::vector<int>* vertex_index, std::vector<double>* coefficient, double* x, double* b, double* result,
		int vertex_index_begin, int vertex_index_end, int sys_size);

private:
	struct AJacobiOperator
	{
		std::vector<std::vector<int>> vertex_index;
		std::vector<std::vector<double>> coefficient;
	};

	//AJacobiOperator A_jacobi_operator_1;
	//AJacobiOperator A_jacobi_operator_2;
	//AJacobiOperator A_jacobi_operator_3;

	std::vector<SparseMatrix<double, RowMajor>> R_Jacobi;
	std::vector<VectorXd> global_diagonal_inv;

	

	std::vector<int>dimension_per_thread;

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


	void createSuperJacobiOperator(AJacobiOperator* A_jacobi_operator, SparseMatrix<double, ColMajor>& R_jacobi,
		AJacobiOperator* A_jacobi_basic);
	void createHighOrderSuperJacobiMethod(AJacobiOperator* A_jacobi_operator_basic, AJacobiOperator* A_jacobi_operator_need_to_multi, AJacobiOperator* A_jacobi_operator);
	
	double obtainElementValue(std::vector<int>* vertex_index_of_i, std::vector<double>* value_of_i,
		std::vector<int>* vertex_index_of_j, std::vector<double>* value_of_j, std::vector<int>& indicate_if_vertex_exists);


	void testOperator(AJacobiOperator* A_jacobi_operator, SparseMatrix<double, ColMajor>& R_jacobi);
	VectorXd RMultiX(AJacobiOperator* A_jacobi_operator, VectorXd& x);

	void RMultiXPlusb(AJacobiOperator* A_jacobi_operator, double* x, double* b, double* result);

	template <class T>
	void superJacobiSingleIteration(Matrix<T, -1, 1>& u, Matrix<T, -1, 1>& b, Matrix<T, -1, 1>& global_diagonal_inv, SparseMatrix<T, ColMajor>& R_Jacobi,
		int super_jacobi_step_size)
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
	void prepareR(SparseMatrix<T, ColMajor>& R_Jacobi, Matrix<T, -1, 1>& global_diagonal_inverse, SparseMatrix<T, ColMajor>& system_matrix)
	{
		int size = system_matrix.cols();
		R_Jacobi = system_matrix;
		global_diagonal_inverse.resize(size);

		for (int i = 0; i < size; ++i) {
			R_Jacobi.coeffRef(i, i) = 0.0;
			global_diagonal_inverse[i] = 1.0 / system_matrix.coeff(i, i);
		}
		R_Jacobi *= -1.0;

		for (int j = 0; j < R_Jacobi.cols(); ++j) {
			for (int k = R_Jacobi.outerIndexPtr()[j]; k < R_Jacobi.outerIndexPtr()[j + 1]; ++k) {
				R_Jacobi.valuePtr()[k] *= global_diagonal_inverse.data()[R_Jacobi.innerIndexPtr()[k]];
			}
		}
	}


	template <class T>
	void solveByChebyshevSemiIterativeSuperJacobi_2(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, SparseMatrix<T, ColMajor>& R_Jacobi, Matrix<T, -1, 1>& global_diagonal_inv, double super_jacobi_spectral_radius_square, int& itr_num, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, int super_jacobi_step_size, std::vector<int>& time,
		AJacobiOperator* a_jacobi_operator_1, AJacobiOperator* a_jacobi_operator_2)
	{
		std::vector<int> vertex_index_begin_thread(thread->thread_num + 1);
		arrangeIndex(thread->thread_num, u[0].size(), vertex_index_begin_thread);

		time_t t0 = clock();
		std::vector<Matrix<T, -1, 1>> u_last = u;
		std::vector<Matrix<T, -1, 1>> u_previous;
		double b_norm = b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm();
		double omega_chebyshev = 2.0;
		double residual_norm;
		relative_error.reserve(max_itr_num);
		time.reserve(max_itr_num);
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		std::vector<VectorXd> R_temp(3);
		std::vector<VectorXd> temp(3);
		int sys_size = global_diagonal_inv.size();
		for (int i = 0; i < 3; ++i)
		{
			temp[i] = b[i].cwiseProduct(global_diagonal_inv);
			R_temp[i].resize(sys_size);
			//RMultiXPlusb(a_jacobi_operator_1, temp[i].data(), temp[i].data(), R_temp[i].data());
			thread->assignTask(this, a_jacobi_operator_1->vertex_index.data(), a_jacobi_operator_1->coefficient.data(),
				temp[i].data(), temp[i].data(), R_temp[i].data(), vertex_index_begin_thread.data(), sys_size);

			//R_temp[i] += RMultiX(a_jacobi_operator_1, R_temp[i]);
		}
		time.push_back(clock() - t0);
		for (int i = 0; i < 3; ++i) {
			//u[i] = R_Jacobi * (R_Jacobi * u[i]) + R_temp[i];
			//u[i] = RMultiX(a_jacobi_operator_2, u[i]) + R_temp[i];
			//RMultiXPlusb(a_jacobi_operator_2, u[i].data(), R_temp[i].data(), temp[i].data());
			thread->assignTask(this, a_jacobi_operator_2->vertex_index.data(), a_jacobi_operator_2->coefficient.data(),
				u[i].data(), R_temp[i].data(), temp[i].data(), vertex_index_begin_thread.data(), sys_size);

			memcpy(u[i].data(), temp[i].data(), 8 * sys_size);
			//superJacobiSingleIteration(u[i], b[i], global_diagonal_inv, R_Jacobi, super_jacobi_step_size);
		}
		itr_num = 1;
		residual_norm = (b[0] - system_matrix * u[0]).squaredNorm() + (b[1] - system_matrix * u[1]).squaredNorm() + (b[2] - system_matrix * u[2]).squaredNorm();
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		time.push_back(clock() - t0);
		while (residual_norm / b_norm > conv_rate_2) {
			u_previous = u;
			omega_chebyshev = 4.0 / (4.0 - super_jacobi_spectral_radius_square * omega_chebyshev);
			for (int i = 0; i < 3; ++i) {
				//u[i] = RMultiX(a_jacobi_operator_2, u[i]) + R_temp[i];
				//RMultiXPlusb(a_jacobi_operator_2, u[i].data(), R_temp[i].data(), temp[i].data());
				thread->assignTask(this, a_jacobi_operator_2->vertex_index.data(), a_jacobi_operator_2->coefficient.data(),
					u[i].data(), R_temp[i].data(), temp[i].data(), vertex_index_begin_thread.data(), sys_size);
				memcpy(u[i].data(), temp[i].data(), 8 * sys_size);
				//u[i] = R_Jacobi * (R_Jacobi * u[i]) + R_temp[i];
				//superJacobiSingleIteration(u[i], b[i], global_diagonal_inv, R_Jacobi, super_jacobi_step_size);
				u[i] = omega_chebyshev * (u[i] - u_last[i]) + u_last[i];
			}
			u_last = u_previous;
			itr_num++;
			residual_norm = (b[0] - system_matrix * u[0]).squaredNorm() + (b[1] - system_matrix * u[1]).squaredNorm() + (b[2] - system_matrix * u[2]).squaredNorm();
			relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
			time.push_back(clock() - t0);
		}
	}

	template <class T>
	void solveByChebyshevSemiIterativeSuperJacobi_3(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, SparseMatrix<T, ColMajor>& R_Jacobi, Matrix<T, -1, 1>& global_diagonal_inv, double super_jacobi_spectral_radius_square, int& itr_num, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, int super_jacobi_step_size, std::vector<int>& time,
		AJacobiOperator* a_jacobi_operator, AJacobiOperator* a_jacobi_operator2, AJacobiOperator* a_jacobi_operator3)
	{
		std::vector<int> vertex_index_begin_thread(thread->thread_num + 1);
		arrangeIndex(thread->thread_num, u[0].size(), vertex_index_begin_thread);
		time_t t0 = clock();
		std::vector<Matrix<T, -1, 1>> u_last = u;
		std::vector<Matrix<T, -1, 1>> u_previous;
		double b_norm = b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm();
		double omega_chebyshev = 2.0;
		double residual_norm;
		relative_error.reserve(max_itr_num);
		time.reserve(max_itr_num);
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());

		std::vector<VectorXd> R_temp(3);
		std::vector<VectorXd> temp(3);
		for (int i = 0; i < 3; ++i)
		{
			R_temp[i] = b[i].cwiseProduct(global_diagonal_inv);
			R_temp[i] += RMultiX(a_jacobi_operator2, R_temp[i]) + RMultiX(a_jacobi_operator, R_temp[i]);
			temp[i].resize(global_diagonal_inv.size());
		}
		time.push_back(clock() - t0);
		//std::vector<VectorXd> temp(3);
		//temp[0] = b[0].cwiseProduct(global_diagonal_inv);
		//temp[1] = b[1].cwiseProduct(global_diagonal_inv);
		//temp[2] = b[2].cwiseProduct(global_diagonal_inv);
		//std::vector<VectorXd> R_temp(3);
		//R_temp[0] = R_Jacobi * temp[0];
		//R_temp[1] = R_Jacobi * temp[1];
		//R_temp[2] = R_Jacobi * temp[2];		
		int sys_size = global_diagonal_inv.size();
		for (int i = 0; i < 3; ++i) {
			//u[i] = RMultiX(a_jacobi_operator3, u[i]) + R_temp[i];
			//RMultiXPlusb(a_jacobi_operator3, u[i].data(), R_temp[i].data(), temp[i].data());
			thread->assignTask(this, a_jacobi_operator3->vertex_index.data(), a_jacobi_operator3->coefficient.data(),
				u[i].data(), R_temp[i].data(), temp[i].data(), vertex_index_begin_thread.data(), sys_size);
			memcpy(u[i].data(), temp[i].data(), 8 * sys_size);
		}		

		//for (int i = 0; i < 3; ++i) {
		//	u[i] = R_Jacobi * (R_Jacobi * (R_Jacobi * u[i])) + temp[i] + R_temp [i] + R_Jacobi * R_temp[i];
		//	//superJacobiSingleIteration(u[i], b[i], global_diagonal_inv, R_Jacobi, super_jacobi_step_size);
		//}
		itr_num = 1;
		residual_norm = (b[0] - system_matrix * u[0]).squaredNorm() + (b[1] - system_matrix * u[1]).squaredNorm() + (b[2] - system_matrix * u[2]).squaredNorm();
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		time.push_back(clock() - t0);
		while (residual_norm / b_norm > conv_rate_2) {
			u_previous = u;
			omega_chebyshev = 4.0 / (4.0 - super_jacobi_spectral_radius_square * omega_chebyshev);
			for (int i = 0; i < 3; ++i) {

				//u[i] = RMultiX(a_jacobi_operator3, u[i]) + R_temp[i];
				//RMultiXPlusb(a_jacobi_operator3, u[i].data(), R_temp[i].data(), temp[i].data());
				thread->assignTask(this, a_jacobi_operator3->vertex_index.data(), a_jacobi_operator3->coefficient.data(),
					u[i].data(), R_temp[i].data(), temp[i].data(), vertex_index_begin_thread.data(), sys_size);
				memcpy(u[i].data(), temp[i].data(), 8 * sys_size);
				//u[i] = R_Jacobi * (R_Jacobi * (R_Jacobi * u[i])) + R_temp[i];
				//superJacobiSingleIteration(u[i], b[i], global_diagonal_inv, R_Jacobi, super_jacobi_step_size);
				u[i] = omega_chebyshev * (u[i] - u_last[i]) + u_last[i];
			}
			u_last = u_previous;
			itr_num++;
			residual_norm = (b[0] - system_matrix * u[0]).squaredNorm() + (b[1] - system_matrix * u[1]).squaredNorm() + (b[2] - system_matrix * u[2]).squaredNorm();
			relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
			time.push_back(clock() - t0);
		}
	}

	template <class T>
	void solveBySuperJacobi_2(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, SparseMatrix<T, ColMajor>& R_Jacobi, Matrix<T, -1, 1>& global_diagonal_inv, int& itr_num, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, int super_jacobi_step_size, SparseMatrix<T, ColMajor>& R_Jacobi_2, std::vector<int>& time,
		AJacobiOperator* a_jacobi_operator_2, AJacobiOperator* a_jacobi_operator_1)
	{
		std::vector<int> vertex_index_begin_thread(thread->thread_num + 1);
		arrangeIndex(thread->thread_num, u[0].size(), vertex_index_begin_thread);
		time_t t0 = clock();
		double residual_norm;
		double b_norm = b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm();
		residual_norm = 2.0 * b_norm;
		itr_num = 0;
		relative_error.reserve(max_itr_num);
		time.reserve(max_itr_num);
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		
		std::vector<VectorXd> R_temp(3);
		std::vector<VectorXd> temp(3);

		int sys_size = global_diagonal_inv.size();

		for (int i = 0; i < 3; ++i)
		{
			temp[i] = b[i].cwiseProduct(global_diagonal_inv);
			R_temp[i].resize(sys_size);
			//RMultiXPlusb(a_jacobi_operator_1, temp[i].data(), temp[i].data(), R_temp[i].data());
			thread->assignTask(this, a_jacobi_operator_1->vertex_index.data(), a_jacobi_operator_1->coefficient.data(),
				temp[i].data(), temp[i].data(), R_temp[i].data(), vertex_index_begin_thread.data(), sys_size);
			//R_temp[i] += RMultiX(a_jacobi_operator_1, R_temp[i]);
		}

		time.push_back(clock() - t0);
		

		while (residual_norm / b_norm > conv_rate_2) {
			for (int i = 0; i < 3; ++i) {
				//u[i] = R_Jacobi_2 * u[i] + R_temp[i];
				//RMultiXPlusb(a_jacobi_operator_2, u[i].data(), R_temp[i].data(), temp[i].data());
				thread->assignTask(this, a_jacobi_operator_2->vertex_index.data(), a_jacobi_operator_2->coefficient.data(),
					u[i].data(), R_temp[i].data(), temp[i].data(), vertex_index_begin_thread.data(), sys_size);
				memcpy(u[i].data(), temp[i].data(), 8 * sys_size);
				//u[i] = RMultiX(a_jacobi_operator_2, u[i]) + R_temp[i];
				//superJacobiSingleIteration(u[i], b[i], global_diagonal_inv, R_Jacobi, super_jacobi_step_size);
				//u[i] = R_Jacobi*(R_Jacobi * u[i]) + R_temp[i];
			}
			residual_norm = (b[0] - system_matrix * u[0]).squaredNorm() + (b[1] - system_matrix * u[1]).squaredNorm() + (b[2] - system_matrix * u[2]).squaredNorm();
			itr_num++;
			relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
			time.push_back(clock() - t0);
		}
	}


	template <class T>
	void solveBySuperJacobi_3(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, SparseMatrix<T, ColMajor>& R_Jacobi, Matrix<T, -1, 1>& global_diagonal_inv, int& itr_num, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, int super_jacobi_step_size, SparseMatrix<T, ColMajor>& R_Jacobi_2, SparseMatrix<T, ColMajor>& R_Jacobi_3, std::vector<int>& time,
		AJacobiOperator* a_jacobi_operator, AJacobiOperator* a_jacobi_operator2, AJacobiOperator* a_jacobi_operator3)
	{
		std::vector<int> vertex_index_begin_thread(thread->thread_num + 1);
		arrangeIndex(thread->thread_num, u[0].size(), vertex_index_begin_thread);
		time_t t0 = clock();
		double residual_norm;
		double b_norm = b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm();
		residual_norm = 2.0 * b_norm;
		itr_num = 0;
		relative_error.reserve(max_itr_num);
		time.reserve(max_itr_num);
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		time.push_back(clock() - t0);

		std::vector<VectorXd> R_temp(3);
		std::vector<VectorXd> temp(3);
		int sys_size = global_diagonal_inv.size();
		for (int i = 0; i < 3; ++i)
		{
			temp[i].resize(sys_size);
			R_temp[i] = b[i].cwiseProduct(global_diagonal_inv);
			R_temp[i] += RMultiX(a_jacobi_operator2, R_temp[i]) + RMultiX(a_jacobi_operator, R_temp[i]);
		}

	
		while (residual_norm / b_norm > conv_rate_2) {
			for (int i = 0; i < 3; ++i) {
				//superJacobiSingleIteration(u[i], b[i], global_diagonal_inv, R_Jacobi, super_jacobi_step_size);
				//u[i] = R_Jacobi_3 * u[i] + R_temp[i];
				//u[i] = RMultiX(a_jacobi_operator3, u[i]) + R_temp[i];
				//RMultiXPlusb(a_jacobi_operator3, u[i].data(), R_temp[i].data(), temp[i].data());
				thread->assignTask(this, a_jacobi_operator3->vertex_index.data(), a_jacobi_operator3->coefficient.data(),
					u[i].data(), R_temp[i].data(), temp[i].data(), vertex_index_begin_thread.data(), sys_size);
				memcpy(u[i].data(), temp[i].data(), 8 * sys_size);
			}
			residual_norm = (b[0] - system_matrix * u[0]).squaredNorm() + (b[1] - system_matrix * u[1]).squaredNorm() + (b[2] - system_matrix * u[2]).squaredNorm();
			itr_num++;
			relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
			time.push_back(clock() - t0);
		}
	}


	template <class T>
	void solveByJacobi(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, SparseMatrix<T, ColMajor>& R_Jacobi, Matrix<T, -1, 1>& global_diagonal_inv, int& itr_num, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, std::vector<int>& time)
	{
		time_t t0 = clock();
		double residual_norm;
		double b_norm = b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm();
		residual_norm = 2.0 * b_norm;
		itr_num = 0;
		relative_error.reserve(max_itr_num);
		time.reserve(max_itr_num);
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		std::vector<VectorXd> temp(3);
		temp[0] = b[0].cwiseProduct(global_diagonal_inv);
		temp[1] = b[1].cwiseProduct(global_diagonal_inv);
		temp[2] = b[2].cwiseProduct(global_diagonal_inv);
		time.push_back(clock() - t0);

		while (residual_norm / b_norm > conv_rate_2) {
			//u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + RMultiplyX(u, cloth_No, 0);
			for (int i = 0; i < 3; ++i) {
				u[i] = temp[i] + R_Jacobi * u[i];
			}
			residual_norm = (b[0] - system_matrix * u[0]).squaredNorm() + (b[1] - system_matrix * u[1]).squaredNorm() + (b[2] - system_matrix * u[2]).squaredNorm();
			relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());		
			itr_num++;
			time.push_back(clock() - t0);
		}
	}

	//template <class T>
	//void solveByJacobiPerThread(int thread_id, std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, SparseMatrix<T, ColMajor>& R_Jacobi, std::vector<Matrix<T, -1, 1>>& global_diagonal_inv, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<std::vector<double>>& relative_error, std::vector<double>& residual_norm)
	//{
	//	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
	//		u[i] = b[i].cwiseProduct(global_diagonal_inv[i]) + R_Jacobi * u[i];
	//		residual_norm[i] = (b[i] - system_matrix * u[i]).squaredNorm();
	//		relative_error[i].push_back((u[i] - ground_truth[i]).squaredNorm());
	//	}
	//}



	template <class T>
	void solveByGaussSeidel(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, int& itr_num, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, std::vector<int>& time)
	{
		time_t t0 = clock();
		double residual_norm;
		double b_norm = b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm();
		residual_norm = 2.0 * b_norm;
		itr_num = 0;
		relative_error.reserve(max_itr_num);
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		time.push_back(clock() - t0);
		while (residual_norm / b_norm > conv_rate_2) {
			for (int i = 0; i < 3; ++i) {
				u[i] = b[i] - system_matrix.triangularView<StrictlyUpper>() * u[i];
				u[i] = system_matrix.triangularView<Lower>().solve(u[i]);
			}
			residual_norm = (b[0] - system_matrix * u[0]).squaredNorm() + (b[1] - system_matrix * u[1]).squaredNorm() + (b[2] - system_matrix * u[2]).squaredNorm();
			relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
			itr_num++;
			time.push_back(clock() - t0);
		}
	}

	//weighted chebyshev gauss seidel, guarantee convergence
	template <class T>
	void solveByChebyshevGaussSeidel(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix,
		double gauss_seidel_spectral_radius_square, int& itr_num, double weight, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, std::vector<int>& time)
	{
		time_t t0 = clock();
		std::vector<Matrix<T, -1, 1>> u_last = u;
		double residual_norm;
		double omega_chebyshev = 2.0;
		double b_norm = b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm();
		relative_error.reserve(max_itr_num);
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		time.push_back(clock() - t0);
		for (int i = 0; i < 3; ++i) {
			u[i] = system_matrix.triangularView<Lower>().solve(b[i] - system_matrix.triangularView<StrictlyUpper>() * u[i]);
		}
		itr_num = 1;
		residual_norm = (b[0] - system_matrix * u[0]).squaredNorm() + (b[1] - system_matrix * u[1]).squaredNorm() + (b[2] - system_matrix * u[2]).squaredNorm();
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		time.push_back(clock() - t0);
		std::vector<Matrix<T, -1, 1>> u_previous;
		while (residual_norm / b_norm > conv_rate_2) {
			u_previous = u;
			omega_chebyshev = 4.0 / (4.0 - gauss_seidel_spectral_radius_square * omega_chebyshev);

			for (int i = 0; i < 3; ++i) {
				u[i] = system_matrix.triangularView<Lower>().solve(b[i] - system_matrix.triangularView<StrictlyUpper>() * u[i]);
				//u = omega_chebyshev * (u - u_last) + u_last;
				u[i] = omega_chebyshev * (weight * (u[i] - u_previous[i]) + u_previous[i] - u_last[i]) + u_last[i];
			}
			u_last = u_previous;
			residual_norm = (b[0] - system_matrix * u[0]).squaredNorm() + (b[1] - system_matrix * u[1]).squaredNorm() + (b[2] - system_matrix * u[2]).squaredNorm();
			itr_num++;
			relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
			time.push_back(clock() - t0);
		}
	}
	template <class T>
	void solveByPCG(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, int& itr_num,
		Matrix<T, -1, 1>& global_diagonal_inv, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error, double conv_rate_2,
		std::vector<int>& time)
	{
		time_t t0 = clock();
		double residual_norm;
		std::vector<Matrix<T, -1, 1>> residual(3);
		double b_norm = b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm();
		std::vector<Matrix<T, -1, 1>> z(3);
		std::vector<Matrix<T, -1, 1>> p(3);
		for (int i = 0; i < 3; ++i) {
			residual[i] = b[i] - system_matrix * u[i];
			z[i] = global_diagonal_inv.cwiseProduct(residual[i]);
		}	
		p = z;
		double alpha, beta;
		double rz_k, rz_k_1;
		rz_k = residual[0].dot(z[0])+ residual[1].dot(z[1])+ residual[2].dot(z[2]);
		itr_num = 0;
		relative_error.reserve(max_itr_num);
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		time.push_back(clock() - t0);
		while (true)
		{
			itr_num++;
			alpha = rz_k / ((system_matrix * p[0]).dot(p[0]) + (system_matrix * p[1]).dot(p[1]) + (system_matrix * p[2]).dot(p[2]));
			for (int i = 0; i < 3; ++i)
			{
				u[i] += alpha * p[i];
				residual[i] = b[i] - system_matrix * u[i];
			}
			residual_norm = residual[0].squaredNorm() + residual[1].squaredNorm() + residual[2].squaredNorm();
			relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
			time.push_back(clock() - t0);
			if (residual_norm / b_norm < conv_rate_2) {
				break;
			}
			for (int i = 0; i < 3; ++i) {
				z[i] = global_diagonal_inv.cwiseProduct(residual[i]);			
			}
			rz_k_1 = residual[0].dot(z[0])+ residual[1].dot(z[1])+ residual[2].dot(z[2]);
			beta = rz_k_1 / rz_k;
			rz_k = rz_k_1;
			for (int i = 0; i < 3; ++i) {
				p[i] = z[i] + beta * p[i];
			}
		}
	}
	template <class T>
	void solveByChebyshevSemiIterativeJacobi(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix,
		SparseMatrix<T, ColMajor>& R_Jacobi, int& itr_num, Matrix<T, -1, 1>& global_diagonal_inv, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double jacobi_spectral_radius_square, double conv_rate_2, std::vector<int>& time)
	{
		time_t t0 = clock();
		std::vector<Matrix<T, -1, 1>> u_last = u;
		std::vector<Matrix<T, -1, 1>> B_inv_b(3);
		for (int i = 0; i < 3; ++i) {
			B_inv_b[i] = global_diagonal_inv.cwiseProduct(b[i]);
		}		
		//u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + RMultiplyX(u, cloth_No, 0);
		relative_error.reserve(max_itr_num);
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		time.push_back(clock() - t0);
		for (int i = 0; i < 3; ++i) {
			u[i] = b[i].cwiseProduct(global_diagonal_inv) + R_Jacobi * u[i];
		}
		itr_num = 1;
		double omega_chebyshev = 2.0;
		std::vector<Matrix<T, -1, 1>> u_previous;
		double b_norm = b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm();
		double residual_chebyshev = 2 * b_norm;
		relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
		
		
		std::vector<VectorXd> temp(3);
		temp[0] = b[0].cwiseProduct(global_diagonal_inv);
		temp[1] = b[1].cwiseProduct(global_diagonal_inv);
		temp[2] = b[2].cwiseProduct(global_diagonal_inv);
		time.push_back(clock() - t0);
		while (residual_chebyshev / b_norm > conv_rate_2)
		{
			u_previous = u;
			omega_chebyshev = 4.0 / (4.0 - jacobi_spectral_radius_square * omega_chebyshev);
			for (int i = 0; i < 3; ++i) {
				u[i] = temp[i] + R_Jacobi * u[i];
				u[i] = omega_chebyshev * (u[i] - u_last[i]) + u_last[i];
			}
			u_last = u_previous;
			itr_num++;
			residual_chebyshev = (b[0] - system_matrix * u[0]).squaredNorm() + (b[1] - system_matrix * u[1]).squaredNorm() + (b[2] - system_matrix * u[2]).squaredNorm();
			relative_error.push_back((u[0] - ground_truth[0]).squaredNorm() + (u[1] - ground_truth[1]).squaredNorm() + (u[2] - ground_truth[2]).squaredNorm());
			time.push_back(clock() - t0);
		}
	}
	template <class T>
	double estimateGaussSeidelEigenValue(std::vector<Matrix<T, -1, 1>>& u, SparseMatrix<T, ColMajor>& system_matrix)
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
	void jacobi(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, std::vector<int>& time)
	{
		SparseMatrix<T, ColMajor> R_Jacobi;
		Matrix<T, -1, 1> global_diagonal_inverse;
		prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
		std::cout << R_Jacobi.cols() << " " << R_Jacobi.rows() << std::endl;

		int itr_num;
		std::vector<double> relat;
		solveByJacobi(u, b, system_matrix, R_Jacobi, global_diagonal_inverse, itr_num, ground_truth, relat, conv_rate_2, time);
		computeRelativeError(relat, relative_error, ground_truth);

	}

	template <class T>
	void computeRelativeError(std::vector<double>& relat, std::vector<double>& relative_error, std::vector<Matrix<T, -1, 1>>& ground_truth)
	{
		double ground_truth_norm;
		ground_truth_norm = ground_truth[0].squaredNorm() + ground_truth[1].squaredNorm() + ground_truth[2].squaredNorm();
		//ground_truth_norm = ground_truth[1].squaredNorm();
		relative_error.resize(relat.size());
		for (int i = 0; i < relative_error.size(); ++i) {
			relative_error[i] = log10(sqrt(relat[i] / ground_truth_norm));
		}
	}

	template <class T>
	void superJacobi(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, double super_jacobi_step_size, std::vector<int>& time)
	{
		SparseMatrix<T, ColMajor> R_Jacobi;
		Matrix<T, -1, 1> global_diagonal_inverse;
		std::vector<double> relat;
		computeRelativeError(relat, relative_error, ground_truth);
		prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
		SparseMatrix<T, ColMajor> R_Jacobi_2 = R_Jacobi * R_Jacobi;
		SparseMatrix<T, ColMajor> R_Jacobi_3 = R_Jacobi_2* R_Jacobi;
		int itr_num;

		AJacobiOperator basic_a_jacobi;
		AJacobiOperator a_jacobi_1;
		AJacobiOperator a_jacobi_2;
		AJacobiOperator a_jacobi_3;
		createSuperJacobiOperator(&a_jacobi_1, R_Jacobi, &basic_a_jacobi);
		createHighOrderSuperJacobiMethod(&basic_a_jacobi, &a_jacobi_1, &a_jacobi_2);
		createHighOrderSuperJacobiMethod(&basic_a_jacobi, &a_jacobi_2, &a_jacobi_3);

		//for (int i = 0; i < siz; ++i)
		//{
		//	result += basic_a_jacobi;
		//}
		//std::cout << "=====" << std::endl;
		//for (int i = 0; i < basic_a_jacobi.vertex_index[0].size(); ++i)
		//{
		//	std::cout << basic_a_jacobi.vertex_index[0][i] << " ";
		//}
		//std::cout << std::endl;
		//for (int i = 0; i < basic_a_jacobi.coefficient[0].size(); ++i)
		//{
		//	std::cout << basic_a_jacobi.coefficient[0][i] << " ";
		//}
		//std::cout << std::endl;
		//for (int i = 0; i < a_jacobi_2.vertex_index[0].size(); ++i)
		//{
		//	std::cout << a_jacobi_2.vertex_index[0][i] << " ";
		//}
		//std::cout << std::endl;
		//for (int i = 0; i < a_jacobi_2.coefficient[0].size(); ++i)
		//{
		//	std::cout << a_jacobi_2.coefficient[0][i] << " ";
		//}
		//std::cout << std::endl;

		//testOperator(&a_jacobi_1, R_Jacobi);
		//testOperator(&a_jacobi_2, R_Jacobi_2);
		//testOperator(&a_jacobi_3, R_Jacobi_3);

		if (super_jacobi_step_size == 2) {
			solveBySuperJacobi_2(u, b, system_matrix, R_Jacobi, global_diagonal_inverse, itr_num, ground_truth, relat, conv_rate_2,
				super_jacobi_step_size, R_Jacobi_2, time, &a_jacobi_2, &a_jacobi_1);
		}
		else if (super_jacobi_step_size == 3) {
			solveBySuperJacobi_3(u, b, system_matrix, R_Jacobi, global_diagonal_inverse, itr_num, ground_truth, relat, conv_rate_2,
				super_jacobi_step_size, R_Jacobi_2, R_Jacobi_3, time, &a_jacobi_1, &a_jacobi_2, &a_jacobi_3);
		}
		computeRelativeError(relat, relative_error, ground_truth);
	}
	template <class T>
	void gauss_seidel(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, std::vector<int>& time)
	{
		int itr_num;
		std::vector<double> relat;
		solveByGaussSeidel(u, b, system_matrix, itr_num, ground_truth, relat, conv_rate_2,time);		
		computeRelativeError(relat, relative_error, ground_truth);
	}


	template <class T>
	void chebyshevSemiIterativeJacobi(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, std::vector<int>& time)
	{
		SparseMatrix<T, ColMajor> R_Jacobi;
		Matrix<T, -1, 1> global_diagonal_inverse;
		std::vector<double> relat;
		prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
		int itr_num;
		double eigen_value_square = estimateSuperJacobiEigenValue(u, R_Jacobi, 1);
		std::cout << "eigen value " << eigen_value_square << std::endl;
		solveByChebyshevSemiIterativeJacobi(u, b, system_matrix, R_Jacobi, itr_num, global_diagonal_inverse, ground_truth, relat,
			eigen_value_square, conv_rate_2,time);
		
		computeRelativeError(relat, relative_error, ground_truth);
	}

	template <class T>
	double estimateSuperJacobiEigenValue(std::vector<Matrix<T, -1, 1>>& u, SparseMatrix<T, ColMajor>& R_Jacobi, int super_jacobi_step_size)
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
	void chebyshevSemiIterativeSuperJacobi(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, int super_jacobi_step_size, std::vector<int>& time)
	{
		SparseMatrix<T, ColMajor> R_Jacobi;
		Matrix<T, -1, 1> global_diagonal_inverse;
		prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
		int itr_num;
		double eigen_value_square = estimateSuperJacobiEigenValue(u, R_Jacobi, 2);
		std::vector<double> relat;

		SparseMatrix<T, ColMajor> R_Jacobi_2 = R_Jacobi * R_Jacobi;
		SparseMatrix<T, ColMajor> R_Jacobi_3 = R_Jacobi_2 * R_Jacobi;

		AJacobiOperator basic_a_jacobi;
		AJacobiOperator a_jacobi_1;
		AJacobiOperator a_jacobi_2;
		AJacobiOperator a_jacobi_3;
		createSuperJacobiOperator(&a_jacobi_1, R_Jacobi, &basic_a_jacobi);
		createHighOrderSuperJacobiMethod(&basic_a_jacobi, &a_jacobi_1, &a_jacobi_2);
		createHighOrderSuperJacobiMethod(&basic_a_jacobi, &a_jacobi_2, &a_jacobi_3);


		if (super_jacobi_step_size == 2) {
			solveByChebyshevSemiIterativeSuperJacobi_2(u, b, system_matrix, R_Jacobi, global_diagonal_inverse, eigen_value_square, itr_num,
				ground_truth, relat, conv_rate_2, super_jacobi_step_size, time, &a_jacobi_1, &a_jacobi_2);
		}
		else if (super_jacobi_step_size == 3)
		{
			solveByChebyshevSemiIterativeSuperJacobi_3(u, b, system_matrix, R_Jacobi, global_diagonal_inverse, eigen_value_square, itr_num,
				ground_truth, relat, conv_rate_2, super_jacobi_step_size, time, &a_jacobi_1, &a_jacobi_2, &a_jacobi_3);
		}
		computeRelativeError(relat, relative_error, ground_truth);
	}




	template <class T>
	void chebyshev_gauss_seidel(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, std::vector<int>& time)
	{
		int itr_num;
		std::vector<double> relat;
		double eigen_value_square = estimateGaussSeidelEigenValue(u, system_matrix);	
		solveByChebyshevGaussSeidel(u, b, system_matrix, eigen_value_square, itr_num, 0.125, ground_truth, relat, conv_rate_2, time);	
		computeRelativeError(relat, relative_error, ground_truth);
	}
	template <class T>
	void PCG(std::vector<Matrix<T, -1, 1>>& u, std::vector<Matrix<T, -1, 1>>& b, SparseMatrix<T, ColMajor>& system_matrix, std::vector<Matrix<T, -1, 1>>& ground_truth, std::vector<double>& relative_error,
		double conv_rate_2, std::vector<int>& time)
	{
		int itr_num;
		SparseMatrix<T, ColMajor> R_Jacobi;
		Matrix<T, -1, 1> global_diagonal_inverse;
		std::vector<double> relat;
		prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);	
		solveByPCG(u, b, system_matrix, itr_num, global_diagonal_inverse, ground_truth, relat, conv_rate_2,time);
		computeRelativeError(relat, relative_error, ground_truth);
	}

	//std::vector<int>dimension_per_thread_test;
};

