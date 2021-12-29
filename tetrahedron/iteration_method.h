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
	void solveByChebyshevGaussSeidel(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num);

	void solveByChebyshevSemiIterativeJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num);
	void solveByPCG(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num);
	void updateGlobalDiagonalInv();
	void initialGlobalDiagonalInv(std::vector<std::vector<double*>>* cloth_global_mat_diagonal_ref_address);
private:
	std::vector<SparseMatrix<double, RowMajor>> R_Jacobi;
	std::vector<VectorXd> global_diagonal_inv;
	std::vector<SparseMatrix<double, RowMajor>> off_diagonal;

	std::vector<std::vector<double*>>* global_mat_diagonal_ref_address;//the address to quickly set the diagonal of system matrix A
	int total_object_num;
	std::vector<int>sys_size;
	Thread* thread;
	std::vector<int> obj_per_thread_begin;
	void updateJacobi(int cloth_No);
	int max_jacobi_itr_num;
	double convergence_rate_2;
	void superJacobiSingleIteration(VectorXd& u, VectorXd& b, int cloth_No);
	int super_jacobi_step_size = 2;
	void estimateSuperJacobiEigenValue(int cloth_No, std::vector<VectorXd>& u);
	std::vector<double> super_jacobi_spectral_radius_square;
	std::vector<double> jacobi_spectral_radius_square;
	std::vector<double> gauss_seidel_spectral_radius_square;

	void estimateGaussSeidelEigenValue(int cloth_No, std::vector<VectorXd>& u, SparseMatrix<double, RowMajor>& system_matrix);
	double gamma=0.6;//for Chebyshev Gauss Seidel, guarantee convergence;

	void estimateJacobiEigenValue(int cloth_No, std::vector<VectorXd>& u);
};

