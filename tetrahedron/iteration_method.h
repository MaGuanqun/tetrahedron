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

	void setBasicInfo(int object_num, std::vector<int>& sys_size, Thread* thread, std::vector<int>&cloth_per_thread_begin);
	void offDiagonalSize();
	void initialJacobi(std::vector<std::vector<double*>>* cloth_global_mat_diagonal_ref_address);
	void updateJacobi_R(int thread_id);
	void setOffDiagonal(int obj_No, std::vector<Triplet<double>>& global_mat_nnz);
	void solveByJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num);

	void initialGaussSeidel();
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
};

