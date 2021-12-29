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
	super_jacobi_spectral_radius_square.resize(object_num);
	jacobi_spectral_radius_square.resize(object_num);
	gauss_seidel_spectral_radius_square.resize(object_num);
}


void IterationMethod::updateConvergenceRate(double conv_rate)
{
	convergence_rate_2= conv_rate * conv_rate;
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


void IterationMethod::initialGlobalDiagonalInv(std::vector<std::vector<double*>>* cloth_global_mat_diagonal_ref_address)
{
	this->global_mat_diagonal_ref_address = cloth_global_mat_diagonal_ref_address;
	global_diagonal_inv.resize(total_object_num);
	for (int i = 0; i < total_object_num; ++i) {
		global_diagonal_inv[i].resize(sys_size[i]);
		for (int j = 0; j < sys_size[i]; ++j) {
			global_diagonal_inv[i].data()[j] = 1.0 / (*((*cloth_global_mat_diagonal_ref_address)[i][j]));
		}
	}
}

void IterationMethod::initialJacobi()
{	
	R_Jacobi = off_diagonal;	
	for (int i = 0; i < total_object_num; ++i) {
		for (int j = 0; j < sys_size[i]; ++j) {
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


void IterationMethod::updateGlobalDiagonalInv()
{
	for (int j = 0; j < total_object_num; ++j) {
		for (int i = 0; i < sys_size[j]; ++i) {
			global_diagonal_inv[j].data()[i] = 1.0 / (*((*global_mat_diagonal_ref_address)[j][i]));
		}
	}

}

void IterationMethod::solveByJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix,  int cloth_No, int& itr_num)
{
	double residual_norm;
	double b_norm = b.squaredNorm();
	residual_norm = 2.0 * b_norm;
	itr_num = 0;
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_jacobi_itr_num) {
		//u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + RMultiplyX(u, cloth_No, 0);
		u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + R_Jacobi[cloth_No] * u;
		residual_norm = (b - system_matrix * u).squaredNorm();
		itr_num++;
	}
}


void IterationMethod::solveBySuperJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num)
{
	double residual_norm;
	double b_norm = b.squaredNorm();
	residual_norm = 2.0 * b_norm;
	itr_num = 0;
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_jacobi_itr_num) {
		superJacobiSingleIteration(u, b, cloth_No);
		residual_norm = (b - system_matrix * u).squaredNorm();
		itr_num++;
	}

}


void IterationMethod::superJacobiSingleIteration(VectorXd& u, VectorXd& b, int cloth_No)
{
	VectorXd R_D_inv_b(u.size());
	VectorXd R_x(u.size());
	VectorXd sum_R_b(u.size());
	R_D_inv_b = b.cwiseProduct(global_diagonal_inv[cloth_No]);
	sum_R_b = R_D_inv_b;
	R_x = R_Jacobi[cloth_No] * u;
	for (int i = 1; i < super_jacobi_step_size; ++i) {
		R_D_inv_b = R_Jacobi[cloth_No] * R_D_inv_b;
		sum_R_b += R_D_inv_b;
		R_x = R_Jacobi[cloth_No] * R_x;
	}
	u = sum_R_b + R_x;
}




void IterationMethod::estimateSuperJacobiEigenValue(std::vector<std::vector<VectorXd>>& u)
{
	for (int i = 0; i < total_object_num; ++i) {
		estimateSuperJacobiEigenValue(i, u[i]);
	}
	//std::cout << super_jacobi_spectral_radius_square[0] << std::endl;
}
void IterationMethod::estimateJacobiEigenValue(std::vector<std::vector<VectorXd>>& u)
{
	for (int i = 0; i < total_object_num; ++i) {
		estimateJacobiEigenValue(i, u[i]);
	}
	//std::cout << jacobi_spectral_radius_square[0] << std::endl;
}


void IterationMethod::estimateGaussSeidelEigenValue(std::vector<std::vector<VectorXd>>& u, std::vector<SparseMatrix<double, RowMajor>>& system_matrix)
{
	for (int i = 0; i < total_object_num; ++i) {
		estimateGaussSeidelEigenValue(i, u[i], system_matrix[i]);
	}
	//std::cout << gauss_seidel_spectral_radius_square[0] << std::endl;
}

void IterationMethod::estimateGaussSeidelEigenValue(int cloth_No, std::vector<VectorXd>& u, SparseMatrix<double, RowMajor>& system_matrix)
{
	double vec_norm2 = 0;
	double u_norm2 = 0;
	for (int i = 0; i < 3; ++i) {
		vec_norm2 += system_matrix.triangularView<Lower>().solve(system_matrix.triangularView<StrictlyUpper>() * u[i]).squaredNorm();
		u_norm2 += u[i].squaredNorm();
	}
	gauss_seidel_spectral_radius_square[cloth_No] = vec_norm2 / u_norm2;
}

void IterationMethod::estimateJacobiEigenValue(int cloth_No, std::vector<VectorXd>& u)
{
	double vec_norm2 = 0;
	double u_norm2 = 0;
	for (int j = 0; j < 3; ++j) {
		vec_norm2 += (R_Jacobi[cloth_No] * u[j]).squaredNorm();
		u_norm2 += u[j].squaredNorm();
		std::cout << (R_Jacobi[cloth_No] * u[j]).squaredNorm() / u[j].squaredNorm();
	}
	jacobi_spectral_radius_square[cloth_No] = vec_norm2 / u_norm2;
}


void IterationMethod::estimateSuperJacobiEigenValue(int cloth_No, std::vector<VectorXd>& u)
{
	VectorXd Vector_change;
	double vec_norm2=0;
	double u_norm2=0;
	for (int j = 0; j < 3; ++j) {
		Vector_change = R_Jacobi[cloth_No] * u[j];
		for (int i = 1; i < super_jacobi_step_size; ++i) {
			Vector_change = R_Jacobi[cloth_No] * Vector_change;
		}
		vec_norm2 += Vector_change.squaredNorm();
		u_norm2 += u[j].squaredNorm();
	}
	super_jacobi_spectral_radius_square[cloth_No] = vec_norm2 / u_norm2;	
}



void IterationMethod::solveByChebyshevSemiIterativeSuperJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num)
{
	VectorXd u_last = u;
	VectorXd u_previous;
	double b_norm = b.squaredNorm();
	double omega_chebyshev = 2.0;
	double residual_norm;
	superJacobiSingleIteration(u, b, cloth_No);
	itr_num=1;
	residual_norm = (b - system_matrix * u).squaredNorm();
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_jacobi_itr_num) {
		u_previous = u;
		omega_chebyshev = 4.0 / (4.0 - super_jacobi_spectral_radius_square[cloth_No] * omega_chebyshev);
		superJacobiSingleIteration(u, b, cloth_No);
		u = omega_chebyshev * (u - u_last) + u_last;
		u_last = u_previous;
		itr_num++;
		residual_norm = (b - system_matrix * u).squaredNorm();
	}
}

void IterationMethod::solveByGaussSeidel(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num)
{
	double residual_norm;
	double b_norm = b.squaredNorm();
	residual_norm = 2.0 * b_norm;
	itr_num = 0;
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_jacobi_itr_num) {
		u = b - system_matrix.triangularView<StrictlyUpper>() * u;
		u = system_matrix.triangularView<Lower>().solve(u);
		residual_norm = (b - system_matrix * u).squaredNorm();
		itr_num++;
	}
}

void IterationMethod::solveByChebyshevGaussSeidel(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num)
{
	VectorXd u_last = u;
	double residual_norm;
	double omega_chebyshev = 2.0;
	double b_norm = b.squaredNorm();	
	u = b - system_matrix.triangularView<StrictlyUpper>() * u;
	u = system_matrix.triangularView<Lower>().solve(u);
	itr_num = 1;
	residual_norm = (b - system_matrix * u).squaredNorm();
	VectorXd u_previous;
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_jacobi_itr_num) {
		u_previous = u;
		omega_chebyshev = 4.0 / (4.0 - gauss_seidel_spectral_radius_square[cloth_No] * omega_chebyshev);
		u = system_matrix.triangularView<Lower>().solve(b - system_matrix.triangularView<StrictlyUpper>() * u);
		//u = omega_chebyshev * (u - u_last) + u_last;
		u = omega_chebyshev * (gamma * (u - u_previous) + u_previous - u_last) + u_last;
		u_last = u_previous;
		residual_norm = (b - system_matrix * u).squaredNorm();
		itr_num++;
	}
}



void IterationMethod::solveByChebyshevSemiIterativeJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num)
{
	VectorXd u_last = u;
	VectorXd B_inv_b = global_diagonal_inv[cloth_No].cwiseProduct(b);
	//u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + RMultiplyX(u, cloth_No, 0);
	u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + R_Jacobi[cloth_No] * u;
	itr_num = 1;
	double omega_chebyshev = 2.0;
	VectorXd u_previous;
	double b_norm = b.squaredNorm();
	double residual_chebyshev = 2 * b_norm;

	while (residual_chebyshev / b_norm > convergence_rate_2 && itr_num < max_jacobi_itr_num)
	{
		u_previous = u;
		omega_chebyshev = 4.0 / (4.0 - jacobi_spectral_radius_square[cloth_No] * omega_chebyshev);
		u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + R_Jacobi[cloth_No] * u;
		u = omega_chebyshev * (u - u_last) + u_last;
		u_last = u_previous;
		itr_num++;
		residual_chebyshev = (b - system_matrix * u).squaredNorm();
	}
}

void IterationMethod::solveByPCG(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num)
{
	VectorXd residual;
	double residual_norm;
	double b_norm = b.squaredNorm();
	VectorXd z, p;
	residual = b - system_matrix * u;
	z = global_diagonal_inv[cloth_No].cwiseProduct(residual);
	p = z;
	double alpha, beta;
	double rz_k, rz_k_1;
	rz_k = residual.dot(z);
	itr_num = 0;
	while (true)
	{
		itr_num++;
		alpha = rz_k / (system_matrix * p).dot(p);
		u += alpha * p;
		residual = b - system_matrix * u;
		if (residual.squaredNorm() / b_norm < convergence_rate_2 || itr_num >= max_jacobi_itr_num) {
			break;
		}
		z = global_diagonal_inv[cloth_No].cwiseProduct(residual);
		rz_k_1 = residual.dot(z);
		beta = rz_k_1 / rz_k;
		rz_k = rz_k_1;
		p = z + beta * p;
		std::cout << residual.squaredNorm() / b_norm << std::endl;
	}
}