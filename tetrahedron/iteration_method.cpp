#include"iteration_method.h"
#include"basic/EigenMatrixIO.h"
#include<algorithm>

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
		//std::cout << residual.squaredNorm() / b_norm << std::endl;
	}
}






void IterationMethod::test()
{
	//load matrix
	std::vector<VectorXd> u_(3);
	std::vector<VectorXd> b_(3);
	SparseMatrix<double, RowMajor> system_matrix_;
	std::vector<VectorXf> ground_truth(3);

	//dimension_per_thread_test.resize(thread->thread_num + 1, 3);
	//for (int i = 0; i < 3; ++i) {
	//	dimension_per_thread_test[i] = i;
	//}
	std::string file_name_matrix="./back/global.dat";

	EigenMatrixIO::read_sp_binary(file_name_matrix.c_str(), system_matrix_);
	SparseMatrix<float, RowMajor> system_matrix;


	system_matrix = system_matrix_.cast<float>();

	std::vector<std::string> file_name_u(3);
	std::vector<std::string> file_name_b(3);

	for (int i = 0; i < 3; ++i) {
		file_name_u[i] = "./back/u" + std::to_string(i) + ".dat";
		file_name_b[i] = "./back/b" + std::to_string(i) + ".dat";
		EigenMatrixIO::read_binary(file_name_u[i].c_str(), u_[i]);
		EigenMatrixIO::read_binary(file_name_b[i].c_str(), b_[i]);
		std::cout << "i" << std::endl;
	}	


	std::vector<VectorXf> u(3);
	std::vector<VectorXf> b(3);

	for (int i = 0; i < 3; ++i) {
		u[i] = u_[i].cast<float>();
		b[i] = b_[i].cast<float>();
	}
	//ground_truth
	SimplicialLLT<SparseMatrix<float>> collision_free_cloth_llt;
	collision_free_cloth_llt.compute(system_matrix);
	for (int i = 0; i < 3; ++i) {
		ground_truth[i] = collision_free_cloth_llt.solve(b[i]);
	}	

	std::vector<VectorXf> u_use;
	//jacobi
	u_use = u;
	std::vector<double> jacobi_relative_error;
	double conv_rate_2 = 1e-7;
	conv_rate_2 *= conv_rate_2;

	jacobi(u_use, b, system_matrix, ground_truth, jacobi_relative_error, conv_rate_2);

	//super_jacobi
	u_use = u;
	std::vector<double> super_jacobi_relative_error;
	superJacobi(u_use, b, system_matrix, ground_truth, super_jacobi_relative_error, conv_rate_2);

	//gauss_seidel
	u_use = u;
	std::vector<double> gauss_seidel_relative_error;
	gauss_seidel(u_use, b, system_matrix, ground_truth, gauss_seidel_relative_error, conv_rate_2);

	//jacobi_chebyshev
	u_use = u;
	std::vector<double> jacobi_chebyshev_relative_error;
	chebyshevSemiIterativeJacobi(u_use, b, system_matrix, ground_truth, jacobi_chebyshev_relative_error, conv_rate_2);
	
	//super_jacobi_chebyshev
	u_use = u;
	std::vector<double> super_jacobi_chebyshev_relative_error;
	chebyshevSemiIterativeSuperJacobi(u_use, b, system_matrix, ground_truth, super_jacobi_chebyshev_relative_error, conv_rate_2);

	//gauss_seidel_chebyshev
	u_use = u;
	std::vector<double> gauss_seidel_chebyshev_relative_error;
	chebyshev_gauss_seidel(u_use, b, system_matrix, ground_truth, gauss_seidel_chebyshev_relative_error, conv_rate_2);

	//PCG
	u_use = u;
	std::vector<double> PCG_chebyshev_relative_error;
	PCG(u, b, system_matrix, ground_truth, PCG_chebyshev_relative_error, conv_rate_2);

	
	size_t max_itr = 0;
	max_itr = (std::max)(max_itr, jacobi_relative_error.size());
	max_itr = (std::max)(max_itr, super_jacobi_relative_error.size());
	max_itr = (std::max)(max_itr, gauss_seidel_relative_error.size());
	max_itr = (std::max)(max_itr, jacobi_chebyshev_relative_error.size());
	max_itr = (std::max)(max_itr, super_jacobi_chebyshev_relative_error.size());
	max_itr = (std::max)(max_itr, gauss_seidel_chebyshev_relative_error.size());
	max_itr = (std::max)(max_itr, PCG_chebyshev_relative_error.size());

	std::cout << "1.jacobi 2.super_jacobi 3.gauss_seidel 4.jacobi_chebyshev 5.super_jacobi_chebyshev 6.gauss_seidel_chebyshev 7.PCG" << std::endl;

	for (int i = 0; i < max_itr; ++i) {
		if (i < jacobi_relative_error.size()) {
			std::cout << log10(jacobi_relative_error[i]);
		}
		else {
			std::cout << " ++";
		}
		if (i < super_jacobi_relative_error.size()) {
			std::cout << " " << log10(super_jacobi_relative_error[i]);
		}
		else {
			std::cout << " ++";
		}
		if (i < gauss_seidel_relative_error.size()) {
			std::cout << " " << log10(gauss_seidel_relative_error[i]);
		}
		else {
			std::cout << " ++";
		}
		if (i < jacobi_chebyshev_relative_error.size()) {
			std::cout << " " << log10(jacobi_chebyshev_relative_error[i]);
		}
		else {
			std::cout << " ++";
		}
		if (i < super_jacobi_chebyshev_relative_error.size()) {
			std::cout << " " << log10(super_jacobi_chebyshev_relative_error[i]);
		}
		else {
			std::cout << " ++";
		}
		if (i < gauss_seidel_chebyshev_relative_error.size()) {
			std::cout << " " << log10(gauss_seidel_chebyshev_relative_error[i]);
		}
		else {
			std::cout << " ++";
		}
		if (i < PCG_chebyshev_relative_error.size()) {
			std::cout << " " << log10(PCG_chebyshev_relative_error[i]);
		}
		else {
			std::cout << " ++";
		}
		std::cout << std::endl;
	}
}

//void IterationMethod::jacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, VectorXd& ground_truth, std::vector<double>&relative_error)



void IterationMethod::solveByChebyshevSemiIterativeSuperJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num)
{
	VectorXd u_last = u;
	VectorXd u_previous;
	double b_norm = b.squaredNorm();
	double omega_chebyshev = 2.0;
	double residual_norm;
	superJacobiSingleIteration(u, b, cloth_No);
	itr_num = 1;
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
