#include"iteration_method.h"


void IterationMethod::setConvergenceRate(double conv_rate, int max_jacobi_itr_num)
{
	convergence_rate_2 = conv_rate* conv_rate;
	this->max_jacobi_itr_num = max_jacobi_itr_num;
}


void IterationMethod::setBasicInfo(int object_num, std::vector<int>&sys_size, Thread* thread, std::vector<int>& cloth_per_thread_begin)
{
	total_object_num= object_num;
	this->sys_size=sys_size;
	this->thread = thread;
	obj_per_thread_begin = cloth_per_thread_begin;
}


void IterationMethod::offDiagonalSize()
{
	off_diagonal.resize(total_object_num);
}


void IterationMethod::setOffDiagonal(int obj_No, std::vector<Triplet<double>>&global_mat_nnz)
{
	off_diagonal[obj_No].resize(sys_size[obj_No], sys_size[obj_No]);
	off_diagonal[obj_No].setFromTriplets(global_mat_nnz.begin(), global_mat_nnz.end());
	off_diagonal[obj_No] *= -1.0;
}

void IterationMethod::initialJacobi(std::vector<std::vector<double*>>*cloth_global_mat_diagonal_ref_address)
{
	this->global_mat_diagonal_ref_address = cloth_global_mat_diagonal_ref_address;
	R_Jacobi = off_diagonal;
	global_diagonal_inv.resize(total_object_num);
	for (int i = 0; i < total_object_num; ++i) {
		global_diagonal_inv[i].resize(sys_size[i]);
		for (int j = 0; j < sys_size[i]; ++j) {
			global_diagonal_inv[i].data()[j] = 1.0 / (*((*cloth_global_mat_diagonal_ref_address)[i][j]));
			for (int k = R_Jacobi[i].outerIndexPtr()[j]; k < R_Jacobi[i].outerIndexPtr()[j + 1]; ++k) {
				R_Jacobi[i].valuePtr()[k] *= global_diagonal_inv[i].data()[j];
			}
		}
	}
}

//UPDATE_JACOBI_R
void IterationMethod::updateJacobi_R(int thread_id)
{
	for (int i = obj_per_thread_begin[thread_id]; i < obj_per_thread_begin[thread_id + 1]; ++i) {
		updateJacobi(i);
	}
}


void IterationMethod::updateJacobi(int cloth_No)
{
	R_Jacobi[cloth_No] = off_diagonal[cloth_No];
	for (int i = 0; i < sys_size[cloth_No]; ++i) {
		global_diagonal_inv[cloth_No].data()[i] = 1.0 / (*((*global_mat_diagonal_ref_address)[cloth_No][i]));
		for (int k = R_Jacobi[cloth_No].outerIndexPtr()[i]; k < R_Jacobi[cloth_No].outerIndexPtr()[i + 1]; ++k) {
			R_Jacobi[cloth_No].valuePtr()[k] *= global_diagonal_inv[cloth_No].data()[i];
		}
	}
}


void IterationMethod::solveByJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix,  int cloth_No, int& itr_num)
{
	VectorXd residual;
	double residual_norm;
	double b_norm = b.dot(b);
	residual_norm = 2.0 * b_norm;
	itr_num = 0;
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_jacobi_itr_num) {
		//u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + RMultiplyX(u, cloth_No, 0);
		u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + R_Jacobi[cloth_No] * u;
		residual = b - system_matrix * u;
		residual_norm = residual.dot(residual);
		itr_num++;
	}
}


