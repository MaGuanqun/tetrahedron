#include"iteration_method.h"
#include"basic/EigenMatrixIO.h"

void IterationMethod::setConvergenceRate(double conv_rate, int max_itr_num)
{
	convergence_rate_2 = conv_rate* conv_rate;
	this->max_itr_num = max_itr_num;
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
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
		//u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + RMultiplyX(u, cloth_No, 0);
		u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + R_Jacobi[cloth_No] * u;
		residual_norm = (b - system_matrix * u).squaredNorm();
		itr_num++;
	}
}



//over-relaxation jacobi (weighted jacobi)
void IterationMethod::solveByWeightedJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num, double weight)
{
	double residual_norm;
	double b_norm = b.squaredNorm();
	residual_norm = 2.0 * b_norm;
	itr_num = 0;
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
		//u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + RMultiplyX(u, cloth_No, 0);
		u = weight * b.cwiseProduct(global_diagonal_inv[cloth_No]) + weight * R_Jacobi[cloth_No] * u + (1.0 - weight) * u;
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
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
		superJacobiSingleIteration(u, b, cloth_No);
		residual_norm = (b - system_matrix * u).squaredNorm();
		itr_num++;
	}

}




void IterationMethod::superJacobiSingleIteration(VectorXd& u, VectorXd& b, int cloth_No)
{
	VectorXd R_D_inv_b = b.cwiseProduct(global_diagonal_inv[cloth_No]);
	VectorXd R_x= R_Jacobi[cloth_No] * u;
	VectorXd sum_R_b = R_D_inv_b;
	for (int i = 1; i < super_jacobi_step_size; ++i) {
		R_D_inv_b = R_Jacobi[cloth_No] * R_D_inv_b;
		sum_R_b += R_D_inv_b;
		R_x = R_Jacobi[cloth_No] * R_x;
	}
	u = sum_R_b + R_x;
}


void IterationMethod::superJacobiSingleIteration(VectorXd& u, VectorXd& b, VectorXd& global_diagonal_inv, SparseMatrix<double, RowMajor>& R_Jacobi)
{
	VectorXd R_D_inv_b = b.cwiseProduct(global_diagonal_inv);
	VectorXd R_x = R_Jacobi * u;
	VectorXd sum_R_b = R_D_inv_b;
	for (int i = 1; i < super_jacobi_step_size; ++i) {
		R_D_inv_b = R_Jacobi * R_D_inv_b;
		sum_R_b += R_D_inv_b;
		R_x = R_Jacobi * R_x;
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
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
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
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
		u = b - system_matrix.triangularView<StrictlyUpper>() * u;
		u = system_matrix.triangularView<Lower>().solve(u);
		residual_norm = (b - system_matrix * u).squaredNorm();
		itr_num++;
	}
}

//weighted chebyshev gauss seidel, guarantee convergence
void IterationMethod::solveByChebyshevGaussSeidel(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num, double weight)
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
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
		u_previous = u;
		omega_chebyshev = 4.0 / (4.0 - gauss_seidel_spectral_radius_square[cloth_No] * omega_chebyshev);
		u = system_matrix.triangularView<Lower>().solve(b - system_matrix.triangularView<StrictlyUpper>() * u);
		//u = omega_chebyshev * (u - u_last) + u_last;
		u = omega_chebyshev * (weight * (u - u_previous) + u_previous - u_last) + u_last;
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

	while (residual_chebyshev / b_norm > convergence_rate_2 && itr_num < max_itr_num)
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
		if (residual.squaredNorm() / b_norm < convergence_rate_2 || itr_num >= max_itr_num) {
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







void IterationMethod::solveBySuperJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, SparseMatrix<double, RowMajor>& R_Jacobi, VectorXd& global_diagonal_inv, int& itr_num, VectorXd& ground_truth, std::vector<double>& relative_error)
{
	double residual_norm;
	double b_norm = b.squaredNorm();
	residual_norm = 2.0 * b_norm;
	itr_num = 0;
	relative_error.reserve(max_itr_num);
	relative_error.push_back((u - ground_truth).squaredNorm());
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
		superJacobiSingleIteration(u, b, global_diagonal_inv, R_Jacobi);
		residual_norm = (b - system_matrix * u).squaredNorm();
		itr_num++;
		relative_error.push_back((u - ground_truth).squaredNorm());
	}
}

void IterationMethod::solveByChebyshevSemiIterativeSuperJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, SparseMatrix<double, RowMajor>& R_Jacobi, VectorXd& global_diagonal_inv, double super_jacobi_spectral_radius_square, int& itr_num, VectorXd& ground_truth, std::vector<double>& relative_error)
{
	VectorXd u_last = u;
	VectorXd u_previous;
	double b_norm = b.squaredNorm();
	double omega_chebyshev = 2.0;
	double residual_norm;
	double ground_truth_norm = ground_truth.squaredNorm();
	relative_error.reserve(max_itr_num);
	relative_error.push_back((u - ground_truth).squaredNorm() / ground_truth_norm);
	superJacobiSingleIteration(u, b, global_diagonal_inv, R_Jacobi);
	itr_num = 1;
	residual_norm = (b - system_matrix * u).squaredNorm();
	relative_error.push_back((u - ground_truth).squaredNorm() / ground_truth_norm);
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
		u_previous = u;
		omega_chebyshev = 4.0 / (4.0 - super_jacobi_spectral_radius_square * omega_chebyshev);
		superJacobiSingleIteration(u, b, global_diagonal_inv, R_Jacobi);
		u = omega_chebyshev * (u - u_last) + u_last;
		u_last = u_previous;
		itr_num++;
		residual_norm = (b - system_matrix * u).squaredNorm();
		relative_error.push_back((u - ground_truth).squaredNorm() / ground_truth_norm);
	}
}


void IterationMethod::solveByJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, SparseMatrix<double, RowMajor>& R_Jacobi, VectorXd& global_diagonal_inv, int& itr_num, VectorXd& ground_truth, std::vector<double>& relative_error)
{
	double residual_norm;
	double b_norm = b.squaredNorm();
	residual_norm = 2.0 * b_norm;
	itr_num = 0;
	relative_error.reserve(max_itr_num);
	relative_error.push_back((u - ground_truth).squaredNorm());
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
		//u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + RMultiplyX(u, cloth_No, 0);
		u = b.cwiseProduct(global_diagonal_inv) + R_Jacobi * u;
		residual_norm = (b - system_matrix * u).squaredNorm();
		relative_error.push_back((u - ground_truth).squaredNorm());
		itr_num++;
	}
}


void IterationMethod::solveByGaussSeidel(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int& itr_num, VectorXd& ground_truth, std::vector<double>& relative_error)
{
	double residual_norm;
	double b_norm = b.squaredNorm();
	residual_norm = 2.0 * b_norm;
	itr_num = 0;
	relative_error.reserve(max_itr_num);
	relative_error.push_back((u - ground_truth).squaredNorm());
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
		u = b - system_matrix.triangularView<StrictlyUpper>() * u;
		u = system_matrix.triangularView<Lower>().solve(u);
		residual_norm = (b - system_matrix * u).squaredNorm();
		relative_error.push_back((u - ground_truth).squaredNorm());
		itr_num++;
	}
}

//weighted chebyshev gauss seidel, guarantee convergence
void IterationMethod::solveByChebyshevGaussSeidel(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, double gauss_seidel_spectral_radius_square, int& itr_num, double weight, VectorXd& ground_truth, std::vector<double>& relative_error)
{
	VectorXd u_last = u;
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
	VectorXd u_previous;
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
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


void IterationMethod::solveByPCG(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int& itr_num,
	VectorXd& global_diagonal_inv, VectorXd& ground_truth, std::vector<double>& relative_error)
{
	VectorXd residual;
	double residual_norm;
	double b_norm = b.squaredNorm();
	VectorXd z, p;
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
		if (residual.squaredNorm() / b_norm < convergence_rate_2 || itr_num >= max_itr_num) {
			break;
		}
		z = global_diagonal_inv.cwiseProduct(residual);
		rz_k_1 = residual.dot(z);
		beta = rz_k_1 / rz_k;
		rz_k = rz_k_1;
		p = z + beta * p;
		
	}
}


void IterationMethod::solveByChebyshevSemiIterativeJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix,
	SparseMatrix<double, RowMajor>& R_Jacobi, int& itr_num, VectorXd& global_diagonal_inv, VectorXd& ground_truth, std::vector<double>& relative_error,
	double jacobi_spectral_radius_square)
{
	VectorXd u_last = u;
	VectorXd B_inv_b = global_diagonal_inv.cwiseProduct(b);
	//u = b.cwiseProduct(global_diagonal_inv[cloth_No]) + RMultiplyX(u, cloth_No, 0);
	relative_error.reserve(max_itr_num);
	relative_error.push_back((u - ground_truth).squaredNorm());
	u = b.cwiseProduct(global_diagonal_inv) + R_Jacobi* u;
	itr_num = 1;
	double omega_chebyshev = 2.0;
	VectorXd u_previous;
	double b_norm = b.squaredNorm();
	double residual_chebyshev = 2 * b_norm;
	relative_error.push_back((u - ground_truth).squaredNorm());
	while (residual_chebyshev / b_norm > convergence_rate_2 && itr_num < max_itr_num)
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

void IterationMethod::test()
{
	//load matrix
	std::vector<VectorXd> u(3);
	std::vector<VectorXd> b(3);
	SparseMatrix<double, RowMajor> system_matrix;
	std::vector<VectorXd> ground_truth(3);

	//dimension_per_thread_test.resize(thread->thread_num + 1, 3);
	//for (int i = 0; i < 3; ++i) {
	//	dimension_per_thread_test[i] = i;
	//}
	std::string file_name_matrix;

	EigenMatrixIO::read_sp_binary(file_name_matrix.c_str(), system_matrix);
	std::vector<std::string> file_name_u(3);
	std::vector<std::string> file_name_b(3);
	for (int i = 0; i < 3; ++i) {
		EigenMatrixIO::read_binary(file_name_u[i].c_str(), u[i]);
		EigenMatrixIO::read_binary(file_name_b[i].c_str(), b[i]);
	}	

	//ground_truth
	SimplicialLLT<SparseMatrix<double>> collision_free_cloth_llt;
	collision_free_cloth_llt.compute(system_matrix);
	for (int i = 0; i < 3; ++i) {
		ground_truth[i] = collision_free_cloth_llt.solve(b[i]);
	}	
	//jacobi
	std::vector<double> jacobi_relative_error;
	jacobi(u, b, system_matrix, ground_truth, jacobi_relative_error);

	//super_jacobi
	std::vector<double> super_jacobi_relative_error;
	superJacobi(u, b, system_matrix, ground_truth, super_jacobi_relative_error);

	//gauss_seidel
	std::vector<double> gauss_seidel_relative_error;
	gauss_seidel(u, b, system_matrix, ground_truth, gauss_seidel_relative_error);

	//jacobi_chebyshev
	std::vector<double> jacobi_chebyshev_relative_error;
	chebyshevSemiIterativeJacobi(u, b, system_matrix, ground_truth, jacobi_chebyshev_relative_error);
	
	//super_jacobi_chebyshev
	std::vector<double> super_jacobi_chebyshev_relative_error;
	chebyshevSemiIterativeSuperJacobi(u, b, system_matrix, ground_truth, super_jacobi_chebyshev_relative_error);

	//gauss_seidel_chebyshev
	std::vector<double> gauss_seidel_chebyshev_relative_error;
	chebyshev_gauss_seidel(u, b, system_matrix, ground_truth, gauss_seidel_chebyshev_relative_error);

	//PCG
	std::vector<double> PCG_chebyshev_relative_error;
	PCG(u, b, system_matrix, ground_truth, PCG_chebyshev_relative_error);

	//
}

//void IterationMethod::jacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, VectorXd& ground_truth, std::vector<double>&relative_error)
void IterationMethod::jacobi(std::vector<VectorXd>& u, std::vector<VectorXd>& b, SparseMatrix<double, RowMajor>& system_matrix, std::vector<VectorXd>& ground_truth, std::vector<double>&relative_error)
{
	SparseMatrix<double, RowMajor> R_Jacobi;
	VectorXd global_diagonal_inverse;
	prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
	int itr_num;	
	std::vector<std::vector<double>> relat(3);
	for (int i = 0; i < 3; ++i) {
		solveByJacobi(u[i], b[i], system_matrix, R_Jacobi, global_diagonal_inverse, itr_num, ground_truth[i], relat[i]);
	}
	computeRelativeError(relat, relative_error, ground_truth);
}


void IterationMethod::prepareR(SparseMatrix<double, RowMajor>& R_Jacobi, VectorXd& global_diagonal_inverse, SparseMatrix<double, RowMajor>& system_matrix)
{
	int size = system_matrix.cols();
	R_Jacobi= system_matrix;
	global_diagonal_inverse.resize(size);

	for(int i=0;i<size;++i){
		R_Jacobi.coeffRef(i, i) = 0.0;
		global_diagonal_inverse[i] = 1.0 / system_matrix.coeff(i, i);
	}
	R_Jacobi.prune(0.0);
	R_Jacobi *= -1.0;
}

void IterationMethod::superJacobi(std::vector<VectorXd>& u, std::vector<VectorXd>& b, SparseMatrix<double, RowMajor>& system_matrix, std::vector<VectorXd>& ground_truth, std::vector<double>& relative_error)
{
	SparseMatrix<double, RowMajor> R_Jacobi;
	VectorXd global_diagonal_inverse;
	std::vector<std::vector<double>> relat(3);
	computeRelativeError(relat, relative_error, ground_truth);
	prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
	int itr_num;
	for (int i = 0; i < 3; ++i) {
		solveBySuperJacobi(u[i], b[i], system_matrix, R_Jacobi, global_diagonal_inverse, itr_num, ground_truth[i], relat[i]);
	}
	computeRelativeError(relat, relative_error, ground_truth);

}



void IterationMethod::gauss_seidel(std::vector<VectorXd>& u, std::vector<VectorXd>& b, SparseMatrix<double, RowMajor>& system_matrix, std::vector<VectorXd>& ground_truth, std::vector<double>& relative_error)
{
	int itr_num; 
	std::vector<std::vector<double>> relat(3);
	for (int i = 0; i < 3; ++i) {
		solveByGaussSeidel(u[i], b[i], system_matrix, itr_num, ground_truth[i], relat[i]);
	}
	computeRelativeError(relat, relative_error, ground_truth);
}

void IterationMethod::chebyshevSemiIterativeJacobi(std::vector<VectorXd>& u, std::vector<VectorXd>& b, SparseMatrix<double, RowMajor>& system_matrix, std::vector<VectorXd>& ground_truth, std::vector<double>& relative_error)
{
	SparseMatrix<double, RowMajor> R_Jacobi;
	VectorXd global_diagonal_inverse;
	std::vector<std::vector<double>> relat(3);
	prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
	int itr_num;
	double eigen_value_square = estimateSuperJacobiEigenValue(u, R_Jacobi, 1);
	for (int i = 0; i < 3; ++i) {
		solveByChebyshevSemiIterativeJacobi(u[i], b[i], system_matrix, R_Jacobi, itr_num, global_diagonal_inverse, ground_truth[i], relat[i],
			eigen_value_square);
	}
	computeRelativeError(relat, relative_error, ground_truth);
}

void IterationMethod::chebyshevSemiIterativeSuperJacobi(std::vector<VectorXd>& u, std::vector<VectorXd>& b, SparseMatrix<double, RowMajor>& system_matrix, std::vector<VectorXd>& ground_truth, std::vector<double>& relative_error)
{
	SparseMatrix<double, RowMajor> R_Jacobi;
	VectorXd global_diagonal_inverse;
	prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
	int itr_num;
	double eigen_value_square = estimateSuperJacobiEigenValue(u, R_Jacobi, 2);
	std::vector<std::vector<double>> relat(3);
	for (int i = 0; i < 3; ++i) {
		solveByChebyshevSemiIterativeSuperJacobi(u[i], b[i], system_matrix, R_Jacobi, global_diagonal_inverse, eigen_value_square, itr_num,
			ground_truth[i], relat[i]);
	}
	computeRelativeError(relat, relative_error, ground_truth);
}


void IterationMethod::chebyshev_gauss_seidel(std::vector<VectorXd>& u, std::vector<VectorXd>& b, SparseMatrix<double, RowMajor>& system_matrix, std::vector<VectorXd>& ground_truth, std::vector<double>& relative_error)
{
	int itr_num;
	std::vector<std::vector<double>> relat(3);
	double eigen_value_square = estimateGaussSeidelEigenValue(u, system_matrix);
	for (int i = 0; i < 3; ++i) {
		solveByChebyshevGaussSeidel(u[i], b[i], system_matrix, eigen_value_square, itr_num, 0.6, ground_truth[i], relat[i]);
	}
	computeRelativeError(relat, relative_error, ground_truth);
}

void IterationMethod::PCG(std::vector<VectorXd>& u, std::vector<VectorXd>& b, SparseMatrix<double, RowMajor>& system_matrix, std::vector<VectorXd>& ground_truth, std::vector<double>& relative_error)
{
	int itr_num;
	SparseMatrix<double, RowMajor> R_Jacobi;
	VectorXd global_diagonal_inverse;
	std::vector<std::vector<double>> relat(3);
	prepareR(R_Jacobi, global_diagonal_inverse, system_matrix);
	for (int i = 0; i < 3; ++i) {
		solveByPCG(u[i], b[i], system_matrix, itr_num, global_diagonal_inverse, ground_truth[i], relat[i]);
	}
	computeRelativeError(relat, relative_error, ground_truth);
}

double IterationMethod::estimateSuperJacobiEigenValue(std::vector<VectorXd>& u, SparseMatrix<double, RowMajor>& R_Jacobi, int super_jacobi_step_size)
{
	VectorXd Vector_change;
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

double IterationMethod::estimateGaussSeidelEigenValue(std::vector<VectorXd>& u, SparseMatrix<double, RowMajor>& system_matrix)
{
	double vec_norm2 = 0;
	double u_norm2 = 0;
	for (int i = 0; i < 3; ++i) {
		vec_norm2 += system_matrix.triangularView<Lower>().solve(system_matrix.triangularView<StrictlyUpper>() * u[i]).squaredNorm();
		u_norm2 += u[i].squaredNorm();
	}
	return vec_norm2 / u_norm2;
}

void IterationMethod::computeRelativeError(std::vector<std::vector<double>>& relat, std::vector<double>& relative_error, std::vector<VectorXd>& ground_truth)
{
	double ground_truth_norm;
	ground_truth_norm = ground_truth[0].squaredNorm() + ground_truth[1].squaredNorm() + ground_truth[2].squaredNorm();
	relative_error.resize(relat[0].size());
	for (int i = 0; i < relative_error.size(); ++i) {
		relative_error[i] = sqrt((relat[0][i] + relat[1][i] + relat[2][i]) / ground_truth_norm);
	}
}




