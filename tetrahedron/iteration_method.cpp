#include"iteration_method.h"
#include"basic/EigenMatrixIO.h"
#include<algorithm>
#include"basic/write_txt.h"



//void IterationMethod::useOperator()
//{
//
//}


void IterationMethod::testOperator(AJacobiOperator* A_jacobi_operator, SparseMatrix<double, ColMajor>& R_jacobi)
{
	std::cout << A_jacobi_operator->coefficient[0][0] << " " << R_jacobi.valuePtr()[0] << std::endl;
	int size;
	for (int i = 0; i < A_jacobi_operator->vertex_index.size(); ++i) {
		size = R_jacobi.outerIndexPtr()[i + 1] - R_jacobi.outerIndexPtr()[i];
		if (A_jacobi_operator->vertex_index[i].size() != size) {
			std::cout << "error occur on count in " << i << " th line" << std::endl;
		}
		else {
			for (int j = 0; j < size; ++j) {
				if (A_jacobi_operator->coefficient[i][j] != R_jacobi.coeff(i, A_jacobi_operator->vertex_index[i][j])) {
					std::cout << "error occur on the element " << i << " " << j<<" "<< A_jacobi_operator->coefficient[i][j] <<" "
						<< R_jacobi.coeff(i,A_jacobi_operator->vertex_index[i][j]) << std::endl;
				}
			}
		}
	}
}

//RX+b result cannot be the same with x or b
void IterationMethod::RMultiXPlusb(AJacobiOperator* A_jacobi_operator, double* x, double* b, double* result)
{
	int* index;
	double* coeff;
	int size = A_jacobi_operator->vertex_index.size();
	memcpy(result, b, 8 * size);
	std::vector<int>* vertex_index = A_jacobi_operator->vertex_index.data();
	std::vector<double>* coefficient = A_jacobi_operator->coefficient.data();
	int col_size;
	for (int i = 0; i < size; ++i)
	{
		col_size = vertex_index[i].size();
		index = vertex_index[i].data();
		coeff = coefficient[i].data();
		for (int j = 0; j < col_size; ++j)
		{
			result[i] += coeff[j] * x[index[j]];
		}
	}
}



void IterationMethod::RMultiXPlusb(std::vector<int>* vertex_index,std::vector<double>* coefficient , double* x, double* b, double* result, 
	int vertex_index_begin, int vertex_index_end, int sys_size)
{
	int* index;
	double* coeff;
	memcpy(result + vertex_index_begin, b + vertex_index_begin, 8 * (vertex_index_end - vertex_index_begin));
	int col_size;
	for (int i = vertex_index_begin; i < vertex_index_end; ++i)
	{
		col_size = vertex_index[i].size();
		index = vertex_index[i].data();
		coeff = coefficient[i].data();
		for (int j = 0; j < col_size; ++j)
		{
			result[i] += coeff[j] * x[index[j]];
		}
	}
}




VectorXd IterationMethod::RMultiX(AJacobiOperator* A_jacobi_operator, VectorXd& x)
{
	std::vector<int>* index;
	double* coeff;
	VectorXd result(x.size());

	std::vector<int>* vertex_index = A_jacobi_operator->vertex_index.data();
	std::vector<double>* coefficient= A_jacobi_operator->coefficient.data();

	for (int i = 0; i < x.size(); ++i)
	{
		index = &vertex_index[i];
		coeff = coefficient[i].data();
		result.data()[i] = coeff[0] * x.data()[index->data()[0]];
		for (int j = 1; j < index->size(); ++j)
		{
			result.data()[i] += coeff[j] * x.data()[index->data()[j]];
		}
	}
	return result;
}


void IterationMethod::createSuperJacobiOperator(AJacobiOperator*  A_jacobi_operator, SparseMatrix<double, ColMajor>& R_jacobi,
	AJacobiOperator* A_jacobi_basic)
{
	int sys_size = R_jacobi.rows();
	A_jacobi_operator->vertex_index.resize(sys_size);
	A_jacobi_operator->coefficient.resize(sys_size);

	A_jacobi_basic->vertex_index.resize(sys_size);
	A_jacobi_basic->coefficient.resize(sys_size);

	std::vector<int>* vertex_index;
	std::vector<double>* coefficient;
	int inner_size;
	int inner_index_start;

	for (int i = 0; i < sys_size; ++i) {
		inner_size = R_jacobi.outerIndexPtr()[i + 1] - R_jacobi.outerIndexPtr()[i];
		A_jacobi_operator->vertex_index[i].reserve(inner_size);
		A_jacobi_operator->coefficient[i].reserve(inner_size);
	}

	for (int i = 0; i < sys_size; ++i) {
		vertex_index = &A_jacobi_basic->vertex_index[i];
		coefficient = &A_jacobi_basic->coefficient[i];
		inner_size = R_jacobi.outerIndexPtr()[i+1]- R_jacobi.outerIndexPtr()[i];
		vertex_index->resize(inner_size);
		coefficient->resize(inner_size);
		inner_index_start = R_jacobi.outerIndexPtr()[i];
		for (int j = 0; j < inner_size; ++j) {
			vertex_index->data()[j] = R_jacobi.innerIndexPtr()[inner_index_start + j];
			coefficient->data()[j] = R_jacobi.valuePtr()[inner_index_start + j];
			A_jacobi_operator->vertex_index[vertex_index->data()[j]].push_back(i);
			A_jacobi_operator->coefficient[vertex_index->data()[j]].push_back(R_jacobi.valuePtr()[inner_index_start + j]);
		}
	}
}


//A_jacobi_operator_need_to_multi * A_jacobi_operator_basic
void IterationMethod::createHighOrderSuperJacobiMethod(AJacobiOperator* A_jacobi_operator_basic, AJacobiOperator* A_jacobi_operator_need_to_multi, AJacobiOperator* A_jacobi_operator)
{
	int sys_size = A_jacobi_operator_basic->vertex_index.size();
	std::vector<bool> is_vertex_used(sys_size,false);
	std::vector<int> indicate_if_vertex_exists(sys_size, -1);
	A_jacobi_operator->vertex_index.resize(sys_size);
	A_jacobi_operator->coefficient.resize(sys_size);

	std::vector<int>* need_to_multi_index;
	std::vector<double>* need_to_multi_value;
	int need_to_multi_index_size;
	int basic_index_in_need_to_multi;
	
	std::vector<int>* basic_index;
	int basic_index_size;
	std::vector<double>* basic_value;

	std::vector<int>* result_index;
	std::vector<double>* result_value;

	double element_value;

	int column_index;
	for (int i = 0; i < sys_size; ++i) {
		need_to_multi_index = &A_jacobi_operator_need_to_multi->vertex_index[i];
		need_to_multi_value = &A_jacobi_operator_need_to_multi->coefficient[i];
		need_to_multi_index_size = need_to_multi_index->size();
		result_index = &A_jacobi_operator->vertex_index[i];
		result_value = &A_jacobi_operator->coefficient[i];
		result_index->reserve(3 * need_to_multi_index_size);
		result_value->reserve(3 * need_to_multi_index_size);
		
		for (int j = 0; j < need_to_multi_index_size; ++j) {
			basic_index_in_need_to_multi = need_to_multi_index->data()[j];
			basic_index = &A_jacobi_operator_basic->vertex_index[basic_index_in_need_to_multi];
			basic_value = &A_jacobi_operator_basic->coefficient[basic_index_in_need_to_multi];
			basic_index_size = basic_index->size();
			for (int k = 0; k < basic_index_size; ++k) {
				column_index = basic_index->data()[k];
				if (!is_vertex_used[column_index]) {
					is_vertex_used[column_index] = true;
					//std::cout << i << " " << column_index << std::endl;
					//if (i == 0 && column_index == 0) {
					//	std::cout << i << " " << A_jacobi_operator_basic->vertex_index[column_index].size() << " " << need_to_multi_index->size() << std::endl;
					//	for (int m = 0; m < A_jacobi_operator_basic->vertex_index[column_index].size(); ++m) {
					//		std::cout << A_jacobi_operator_basic->vertex_index[column_index][m] << " ";
					//	}
					//	std::cout << std::endl;
					//	for (int m = 0; m < A_jacobi_operator_basic->coefficient[column_index].size(); ++m) {
					//		std::cout << A_jacobi_operator_basic->coefficient[column_index][m] << " ";
					//	}
					//	std::cout << std::endl;
					//	for (int m = 0; m < need_to_multi_index->size(); ++m) {
					//		std::cout << need_to_multi_index->data()[m] << " ";
					//	}
					//	std::cout << std::endl;
					//	for (int m = 0; m < need_to_multi_value->size(); ++m) {
					//		std::cout << need_to_multi_value->data()[m] << " ";
					//	}
					//	std::cout << std::endl;
					//}
					result_index->push_back(column_index);
					element_value = obtainElementValue(need_to_multi_index, need_to_multi_value, &A_jacobi_operator_basic->vertex_index[column_index], 
						&A_jacobi_operator_basic->coefficient[column_index], indicate_if_vertex_exists);
					result_value->push_back(element_value);
						//std::cout << "value "<< element_value << std::endl;

				}
			}
		}

		for (int j = 0; j < result_index->size(); ++j) {
			is_vertex_used[result_index->data()[j]] = false;
		}
		//for (int j = 0; j < result_index->size() - 1; ++j) {
		//	if (result_index->data()[j + 1] - result_index->data()[j] < 0) {
		//		std::cout << i<<" order wrong" << std::endl;
		//	}
		//}
	}

}

//to obtain the value of the element (A_jacobi_operator->vertex_index[i],A_jacobi_operator->vertex_index[k])
double IterationMethod::obtainElementValue(std::vector<int>*vertex_index_of_i, std::vector<double>* value_of_i, 
	std::vector<int>* vertex_index_of_j, std::vector<double>* value_of_j,std::vector<int>& indicate_if_vertex_exists)
{
	double value = 0;

	for (int i = 0; i < vertex_index_of_i->size(); ++i)
	{
		indicate_if_vertex_exists[vertex_index_of_i->data()[i]] = i;
	}

	for (int j = 0; j < vertex_index_of_j->size(); ++j) 
	{
		if (indicate_if_vertex_exists[vertex_index_of_j->data()[j]] > -1)
		{
			value += value_of_j->data()[j] * value_of_i->data()[indicate_if_vertex_exists[vertex_index_of_j->data()[j]]];
		}
	}

	for (int i = 0; i < vertex_index_of_i->size(); ++i)
	{
		indicate_if_vertex_exists[vertex_index_of_i->data()[i]] = -1;
	}

	return value;
}


void IterationMethod::setConvergenceRate(double conv_rate, int max_itr_num)
{
	convergence_rate_2 = conv_rate* conv_rate;
	this->max_itr_num = max_itr_num;
}


void IterationMethod::setBasicInfo(int object_num, std::vector<int>&sys_size, Thread* thread, std::vector<int>& cloth_per_thread_begin,
	std::vector<std::vector<VectorXd>>* cloth_b, std::vector<std::vector<VectorXd>>* cloth_u, std::vector<SparseMatrix<double, RowMajor>>* cloth_global_mat)
{
	total_object_num= object_num;
	this->sys_size=sys_size;
	this->thread = thread;
	obj_per_thread_begin = cloth_per_thread_begin;
	a_jacobi_2_spectral_radius_square.resize(object_num);
	a_jacobi_3_spectral_radius_square.resize(object_num);
	jacobi_spectral_radius_square.resize(object_num);
	gauss_seidel_spectral_radius_square.resize(object_num);

	this->cloth_b = cloth_b;
	this->cloth_u = cloth_u;
	this->cloth_global_mat = cloth_global_mat;

	dimension_per_thread.resize(thread->thread_num + 1, 3);
	for (int i = 0; i < 3; ++i) {
		dimension_per_thread[i] = i;
	}
	b_global_inv.resize(3);
	R_b_global_inv.resize(3);

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


void IterationMethod::solveByJacobi(VectorXd* u, VectorXd* b, int cloth_No, int& itr_num)
{
	double residual_norm_per_thread[3];
	double b_norm_cov = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm())* convergence_rate_2;
	double residual_norm = 2.0 * b_norm_cov;
	itr_num = 0;

	for (int i = 0; i < 3; ++i) {
		b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv[cloth_No]);
	}
	while (residual_norm> b_norm_cov && itr_num < max_itr_num)
	{
		thread->assignTask(this, JACOBI_ITR, u, b, residual_norm_per_thread, cloth_No,0.0,u,u);
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;
	}

}

//JACOBI_ITR
void IterationMethod::JacobiIterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, int cloth_No)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i)
	{
		u[i] = b_global_inv[i] + R_Jacobi[cloth_No] * u[i];
		residual_norm[i] = (b[i] - cloth_global_mat->data()[cloth_No] * u[i]).squaredNorm();
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
		superJacobiSingleIteration(u, b, cloth_No,2);
		residual_norm = (b - system_matrix * u).squaredNorm();
		itr_num++;
	}

}

void IterationMethod::solveByAJacobi_2(VectorXd* u, VectorXd* b, int cloth_No, int& itr_num)
{
	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm())* convergence_rate_2;
	double residual_norm = 2.0 * b_norm_conv;
	itr_num = 0;
	std::vector<VectorXd> R_temp(3);
	std::vector<VectorXd> temp(3);
	double residual_norm_per_thread[3];

	for (int i = 0; i < 3; ++i)	{
		b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv[cloth_No]);
		R_b_global_inv[i] = R_Jacobi[cloth_No] * b_global_inv[i] + b_global_inv[i];
	}

	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		thread->assignTask(this, A_JACOBI_2_ITR, u, b, residual_norm_per_thread, cloth_No,0.0,u,u);
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;
	}

}


void IterationMethod::solveByAJacobi_3(VectorXd* u, VectorXd* b, int cloth_No, int& itr_num)
{
	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm()) * convergence_rate_2;
	double residual_norm = 2.0 * b_norm_conv;
	itr_num = 0;
	std::vector<VectorXd> R_temp(3);
	std::vector<VectorXd> temp(3);
	double residual_norm_per_thread[3];

	for (int i = 0; i < 3; ++i) {
		b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv[cloth_No]);
		R_b_global_inv[i] = R_Jacobi[cloth_No] * b_global_inv[i];
		R_b_global_inv[i] = R_Jacobi[cloth_No] * R_b_global_inv[i] + R_b_global_inv[i] + b_global_inv[i];
	}

	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		thread->assignTask(this, A_JACOBI_3_ITR, u, b, residual_norm_per_thread, cloth_No,0.0,u,u);
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;
	}

}

//SUPER_JACOBI_2_ITR
void IterationMethod::SuperJacobi2IterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, int cloth_No)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i)
	{
		u[i] = R_Jacobi[cloth_No] * (R_Jacobi[cloth_No] * u[i]) + R_b_global_inv[i];
		residual_norm[i] = (b[i] - cloth_global_mat->data()[cloth_No] * u[i]).squaredNorm();
	}
}


//SUPER_JACOBI_3_ITR
void IterationMethod::SuperJacobi3IterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, int cloth_No)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i)
	{
		u[i] = R_Jacobi[cloth_No]*(R_Jacobi[cloth_No] * (R_Jacobi[cloth_No] * u[i])) + R_b_global_inv[i];
		residual_norm[i] = (b[i] - cloth_global_mat->data()[cloth_No] * u[i]).squaredNorm();
	}
}







void IterationMethod::superJacobiSingleIteration(VectorXd& u, VectorXd& b, int cloth_No, int super_jacobi_step_size)
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



void IterationMethod::estimateSuperJacobiEigenValue(std::vector<std::vector<VectorXd>>& u, int A_jacobi_step_size)
{
	switch (A_jacobi_step_size)
	{
	case 2:
		for (int i = 0; i < total_object_num; ++i) {
			estimateAJacobi2EigenValue(i, u[i]);
		}
		break;
	case 3:
		for (int i = 0; i < total_object_num; ++i) {
			estimateAJacobi3EigenValue(i, u[i]);
		}
		break;
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
		//std::cout << (R_Jacobi[cloth_No] * u[j]).squaredNorm() / u[j].squaredNorm();
	}
	jacobi_spectral_radius_square[cloth_No] = vec_norm2 / u_norm2;
}


void IterationMethod::estimateAJacobi2EigenValue(int cloth_No, std::vector<VectorXd>& u)
{
	VectorXd Vector_change;
	double vec_norm2=0;
	double u_norm2=0;
	for (int j = 0; j < 3; ++j) {
		Vector_change = R_Jacobi[cloth_No] * (R_Jacobi[cloth_No] * u[j]);
		vec_norm2 += Vector_change.squaredNorm();
		u_norm2 += u[j].squaredNorm();
	}
	a_jacobi_2_spectral_radius_square[cloth_No] = vec_norm2 / u_norm2;
}

void IterationMethod::estimateAJacobi3EigenValue(int cloth_No, std::vector<VectorXd>& u)
{
	VectorXd Vector_change;
	double vec_norm2 = 0;
	double u_norm2 = 0;
	for (int j = 0; j < 3; ++j) {
		Vector_change = R_Jacobi[cloth_No] * u[j];
		for (int i = 1; i < 3; ++i) {
			Vector_change = R_Jacobi[cloth_No] * Vector_change;
		}
		vec_norm2 += Vector_change.squaredNorm();
		u_norm2 += u[j].squaredNorm();
	}
	a_jacobi_3_spectral_radius_square[cloth_No] = vec_norm2 / u_norm2;
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


void IterationMethod::solveByGaussSeidel(VectorXd* u, VectorXd* b, int cloth_No, int& itr_num) 
{
	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm()) * convergence_rate_2;
	double residual_norm = 2.0 * b_norm_conv;
	double residual_norm_per_thread[3];
	itr_num = 0;
	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		thread->assignTask(this, GAUSS_SEIDEL_ITR, u, b, residual_norm_per_thread, cloth_No, 0.0, u, u);
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;
	}
}


//GAUSS_SEIDEL_ITR
void IterationMethod::GaussSeidelIterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, int cloth_No, double omega_chebyshev,
	VectorXd* u_last, VectorXd* u_previous)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u[i] = cloth_global_mat->data()[cloth_No].triangularView<Lower>().solve(
			b[i] - cloth_global_mat->data()[cloth_No].triangularView<StrictlyUpper>() * u[i]);
		residual_norm[i] = (b[i] - cloth_global_mat->data()[cloth_No] * u[i]).squaredNorm();
	}
}



void IterationMethod::solveByChebyshevGaussSeidel(VectorXd* u, VectorXd* b, int cloth_No, int& itr_num, double weight)
{
	std::vector<VectorXd> u_last(3);
	for (int i = 0; i < 3; ++i) {
		u_last[i] = u[i];
	}
	double omega_chebyshev = 2.0;
	weight_for_chebyshev_gauss_seidel = weight;
	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm()) * convergence_rate_2;
	double residual_norm;
	double residual_norm_per_thread[3];
	thread->assignTask(this, GAUSS_SEIDEL_ITR, u, b, residual_norm_per_thread, cloth_No, 0.0, u, u);
	residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
	itr_num = 1;
	std::vector<VectorXd> u_previous(3);
	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		omega_chebyshev = 4.0 / (4.0 - gauss_seidel_spectral_radius_square[cloth_No] * omega_chebyshev);
		thread->assignTask(this, CHEBYSHEV_GAUSS_SEIDEL_ITR, u, b, residual_norm_per_thread, cloth_No, omega_chebyshev, u_last.data(), u_previous.data());
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;
	}

}

//CHEBYSHEV_GAUSS_SEIDEL_ITR
void IterationMethod::ChebyshevSemiIterativeGaussSeidelIterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, int cloth_No, double omega_chebyshev,
	VectorXd* u_last, VectorXd* u_previous)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u_previous[i] = u[i];
		u[i] = cloth_global_mat->data()[cloth_No].triangularView<Lower>().solve(b[i] - cloth_global_mat->data()[cloth_No].triangularView<StrictlyUpper>() * u[i]);
		u[i] = omega_chebyshev * (weight_for_chebyshev_gauss_seidel * (u[i] - u_previous[i]) + u_previous[i] - u_last[i]) + u_last[i];
		u_last[i] = u_previous[i];
		residual_norm[i] = (b[i] - cloth_global_mat->data()[cloth_No] * u[i]).squaredNorm();
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

void IterationMethod::solveByChebyshevSemiIterativeSuperJacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, int cloth_No, int& itr_num)
{
	VectorXd u_last = u;
	VectorXd u_previous;
	double b_norm = b.squaredNorm();
	double omega_chebyshev = 2.0;
	double residual_norm;
	superJacobiSingleIteration(u, b, cloth_No, 2);
	itr_num = 1;
	residual_norm = (b - system_matrix * u).squaredNorm();
	while (residual_norm / b_norm > convergence_rate_2 && itr_num < max_itr_num) {
		u_previous = u;
		omega_chebyshev = 4.0 / (4.0 - a_jacobi_2_spectral_radius_square[cloth_No] * omega_chebyshev);
		superJacobiSingleIteration(u, b, cloth_No, 2);
		u = omega_chebyshev * (u - u_last) + u_last;
		u_last = u_previous;
		itr_num++;
		residual_norm = (b - system_matrix * u).squaredNorm();
	}
}




//CHEBYSHEV_A_JACOBI_2_ITR
void IterationMethod::ChebyshevSemiIterativeAJacobi2IterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, int cloth_No, double omega_chebyshev,
	VectorXd* u_last, VectorXd* u_previous)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u_previous[i] = u[i];
		u[i] = R_b_global_inv[i] + R_Jacobi[cloth_No] * (R_Jacobi[cloth_No] * u[i]);
		u[i] = omega_chebyshev * (u[i] - u_last[i]) + u_last[i];
		u_last[i] = u_previous[i];
		residual_norm[i] = (b[i] - cloth_global_mat->data()[cloth_No] * u[i]).squaredNorm();
	}
}


//CHEBYSHEV_A_JACOBI_3_ITR
void IterationMethod::ChebyshevSemiIterativeAJacobi3IterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, int cloth_No, double omega_chebyshev,
	VectorXd* u_last, VectorXd* u_previous)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u_previous[i] = u[i];
		u[i] = R_b_global_inv[i] + R_Jacobi[cloth_No] *(R_Jacobi[cloth_No] * (R_Jacobi[cloth_No] * u[i]));
		u[i] = omega_chebyshev * (u[i] - u_last[i]) + u_last[i];
		u_last[i] = u_previous[i];
		residual_norm[i] = (b[i] - cloth_global_mat->data()[cloth_No] * u[i]).squaredNorm();
	}
}


void IterationMethod::solveByChebyshevSemiIterativeAJacobi2(VectorXd* u, VectorXd* b, int cloth_No, int& itr_num) {
	std::vector<VectorXd> u_last(3);
	std::vector<VectorXd> u_previous(3);
	for (int i = 0; i < 3; ++i) {
		u_last[i] = u[i];
		b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv[cloth_No]);
		R_b_global_inv[i] = R_Jacobi[cloth_No] * b_global_inv[i] + b_global_inv[i];
	}
	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm()) * convergence_rate_2;
	double omega_chebyshev = 2.0;
	double residual_norm_per_thread[3];
	double residual_norm;
	thread->assignTask(this, A_JACOBI_2_ITR, u, b, residual_norm_per_thread, cloth_No, 0.0, u, u);
	residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
	itr_num = 1;
	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		omega_chebyshev = 4.0 / (4.0 - a_jacobi_2_spectral_radius_square[cloth_No] * omega_chebyshev);
		thread->assignTask(this, CHEBYSHEV_A_JACOBI_2_ITR, u, b, residual_norm_per_thread, cloth_No, omega_chebyshev, u_last.data(), u_previous.data());
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;
	}
}

void IterationMethod::solveByChebyshevSemiIterativeAJacobi3(VectorXd* u, VectorXd* b, int cloth_No, int& itr_num) {
	std::vector<VectorXd> u_last(3);
	std::vector<VectorXd> u_previous(3);
	for (int i = 0; i < 3; ++i) {
		u_last[i] = u[i];
		b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv[cloth_No]);
		R_b_global_inv[i] = R_Jacobi[cloth_No] * b_global_inv[i];
		R_b_global_inv[i] = R_Jacobi[cloth_No] * b_global_inv[i] + R_b_global_inv[i] + b_global_inv[i];
	}
	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm()) * convergence_rate_2;
	double omega_chebyshev = 2.0;
	double residual_norm_per_thread[3];
	double residual_norm;
	thread->assignTask(this, A_JACOBI_3_ITR, u, b, residual_norm_per_thread, cloth_No, 0.0, u, u);
	residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
	itr_num = 1;
	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		omega_chebyshev = 4.0 / (4.0 - a_jacobi_2_spectral_radius_square[cloth_No] * omega_chebyshev);
		thread->assignTask(this, CHEBYSHEV_A_JACOBI_3_ITR, u, b, residual_norm_per_thread, cloth_No, omega_chebyshev, u_last.data(), u_previous.data());
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;
	}
}


void IterationMethod::solveByChebyshevSemiIterativeJacobi(VectorXd* u, VectorXd* b, int cloth_No, int& itr_num)
{
	std::vector<VectorXd> u_last(3);
	std::vector<VectorXd> u_previous(3);
	for (int i = 0; i < 3; ++i) {
		u_last[i] = u[i];
		b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv[cloth_No]);
	}
	double residual_norm_per_thread[3];
	double residual_norm;
	thread->assignTask(this, JACOBI_ITR, u, b, residual_norm_per_thread, cloth_No, 0.0, u, u);
	residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];

	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm()) * convergence_rate_2;
	double omega_chebyshev = 2.0;	
	itr_num = 1;
	while (residual_norm  > b_norm_conv && itr_num < max_itr_num) {
		omega_chebyshev = 4.0 / (4.0 - jacobi_spectral_radius_square[cloth_No] * omega_chebyshev);
		thread->assignTask(this, CHEBYSHEV_JACOBI_ITR, u, b, residual_norm_per_thread, cloth_No, 
			omega_chebyshev, u_last.data(), u_previous.data());
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;	
	}


}


//CHEBYSHEV_JACOBI_ITR
void IterationMethod::ChebyshevSemiIterativeJacobiIterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, int cloth_No, double omega_chebyshev,
	VectorXd* u_last, VectorXd* u_previous)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u_previous[i] = u[i];
		u[i] = b_global_inv[i] + R_Jacobi[cloth_No] * u[i];
		u[i] = omega_chebyshev * (u[i] - u_last[i]) + u_last[i];
		u_last[i] = u_previous[i];
		residual_norm[i] = (b[i] - cloth_global_mat->data()[cloth_No] * u[i]).squaredNorm();
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


void IterationMethod::solveByPCG(VectorXd* u, VectorXd* b, int cloth_No, int& itr_num)
{
	std::vector<VectorXd> residual(3);
	double residual_norm;
	double residual_norm_per_thread[3];
	double b_norm_cov=(b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm())*convergence_rate_2;
	std::vector<VectorXd> z(3);
	std::vector<VectorXd> p(3);
	for (int i = 0; i < 3; ++i) {
		residual[i] = b[i] - cloth_global_mat->data()[cloth_No] * u[i];
		z[i] = global_diagonal_inv[cloth_No].cwiseProduct(residual[i]);
	}
	p = z;
	double alpha, beta;
	double rz_k, rz_k_1;
	double rz_k_1_per_thread[3];
	rz_k = residual[0].dot(z[0]) + residual[1].dot(z[1]) + residual[2].dot(z[2]);
	itr_num = 0;
	while (true)
	{
		itr_num++;
		alpha = rz_k / ((cloth_global_mat->data()[cloth_No] * p[0]).dot(p[0]) + 
			(cloth_global_mat->data()[cloth_No] * p[1]).dot(p[1]) + (cloth_global_mat->data()[cloth_No] * p[2]).dot(p[2]));
		thread->assignTask(this, PCG_ITR1, u, b, residual_norm_per_thread, cloth_No,
			alpha, residual.data(), p.data());
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		if (residual_norm < b_norm_cov || itr_num >= max_itr_num) {
			break;
		}
		thread->assignTask(this, PCG_ITR2, z.data(), residual.data(), rz_k_1_per_thread, cloth_No,
			alpha, residual.data(), p.data());
		rz_k_1 = rz_k_1_per_thread[0] + rz_k_1_per_thread[1] + rz_k_1_per_thread[2];
		beta = rz_k_1 / rz_k;
		rz_k = rz_k_1;
		for (int i = 0; i < 3; ++i) {
			p[i] = z[i] + beta * p[i];
		}
	}
}


//PCG_ITR1
void IterationMethod::PCGIterationPerThread1(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, int cloth_No, double alpha,
	VectorXd* residual, VectorXd* p)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u[i] += alpha* p[i];
		residual[i] = b[i] - cloth_global_mat->data()[cloth_No] * u[i];
		residual_norm[i] = residual[i].squaredNorm();
	}
}

//PCG_ITR2
void IterationMethod::PCGIterationPerThread2(int thread_id, VectorXd* z, VectorXd* residual, double* rz_k_1, int cloth_No, double alpha,
	VectorXd* q, VectorXd* p)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		z[i] = global_diagonal_inv[cloth_No].cwiseProduct(residual[i]);
		rz_k_1[i] = residual[i].dot(z[i]);
	}
}


void IterationMethod::testGaussSeidel()
{
	int size = 5;
	SparseMatrix<double, ColMajor> system_matrix_(size,size);
	std::vector<Eigen::Triplet<double>> triplet;
	std::vector<std::array<double, 5>> ma(5);
	ma[0] =std::array{ 2.1012, 0.1300, -1.6081, -1.1935, 0.3851 };
	ma[1] =std::array{ 0.1300,0.6209,-0.2666,-0.3431,-0.5251 };
	ma[2] =std::array{ -1.6081,-0.2666,2.6402,1.0969,-0.5009 };
	ma[3] =std::array{ -1.1935,-0.3431,1.0969,3.3753,-0.6894 };
	ma[4] =std::array{ 0.3851,-0.5251,-0.5009,-0.6894,1.5310 };

	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			triplet.push_back(Eigen::Triplet<double>(i, j, ma[i][j]));
		}		
	}	
	system_matrix_.setFromTriplets(triplet.begin(), triplet.end());
	VectorXd f(size);
	for (int i = 0; i < size; ++i)
	{
			f[i] = 0.2 * (i+1);
	}
	VectorXd x(size);
	x.setZero();
	std::vector<VectorXd> u(3);
	std::vector<VectorXd> b(3);
	for(int i=0;i<3;++i)
	{
		u[i] = x;
		b[i] = f;
	}

	std::vector<VectorXd> ground_truth(3);
	SimplicialLLT<SparseMatrix<double>> collision_free_cloth_llt;
	collision_free_cloth_llt.compute(system_matrix_);
	for (int i = 0; i < 3; ++i) {
		ground_truth[i] = collision_free_cloth_llt.solve(b[i]);
	}
	std::vector<double> relative_error;
	double conv_rate_2 = 1e-5;
	std::vector<int> time;
	gauss_seidel(u, b, system_matrix_, ground_truth, relative_error, conv_rate_2, time);

}

void IterationMethod::testRelativeError()
{
	std::vector<VectorXd> u_(3);
	std::vector<VectorXd> b_(3);
	std::string file_name_matrix = "./save_matrix_prediction2/global_0.dat";
	SparseMatrix<double, ColMajor> system_matrix_;
	EigenMatrixIO::read_sp_binary(file_name_matrix.c_str(), system_matrix_);
	SparseMatrix<double, ColMajor> system_matrix;
	system_matrix = system_matrix_;// .cast<float>();
	std::vector<std::string> file_name_u(3);
	std::vector<std::string> file_name_b(3);
	std::vector<VectorXd> ground_truth(3);
	//ground_truth
	SimplicialLLT<SparseMatrix<double>> collision_free_cloth_llt;
	collision_free_cloth_llt.analyzePattern(system_matrix);
	collision_free_cloth_llt.factorize(system_matrix_);
	double relative_error;
	double ground_truth_norm;
	//std::string file_name_matrix1 = "./save_matrix_prediction2/global_0.dat";
	for (int k = 16; k < 17; ++k)
	{
		//std::cout << k << std::endl;
		//std::string file_name_matrix1 = "./save_matrix_prediction2/global_" + std::to_string(k) + ".dat";
		//EigenMatrixIO::read_sp_binary(file_name_matrix1.c_str(), system_matrix_);
		//collision_free_cloth_llt.factorize(system_matrix_);
		//system_matrix = system_matrix_;
		for (int i = 0; i < 3; ++i) {			
			file_name_u[i] = "./save_matrix_prediction2/u_"+ std::to_string(k) + "_" + std::to_string(i) + ".dat";
			file_name_b[i] = "./save_matrix_prediction2/b_" + std::to_string(k) + "_" + std::to_string(i) + ".dat";
			EigenMatrixIO::read_binary(file_name_u[i].c_str(), u_[i]);
			EigenMatrixIO::read_binary(file_name_b[i].c_str(), b_[i]);
			for (int i = 0; i < 3; ++i) {
				ground_truth[i] = collision_free_cloth_llt.solve(b_[i]);
			}	
		}
		//relative_error = (u_[0] - ground_truth[0]).squaredNorm() + (u_[1] - ground_truth[1]).squaredNorm() + (u_[2] - ground_truth[2]).squaredNorm();
		//ground_truth_norm = ground_truth[0].squaredNorm() + ground_truth[1].squaredNorm() + ground_truth[2].squaredNorm();
		//relative_error = sqrt(relative_error / ground_truth_norm);
		//std::cout << k << " " << relative_error << std::endl;
		std::cout << "successful read" << std::endl;
		std::vector<VectorXd> u(3);
		std::vector<VectorXd> b(3);
		for (int i = 0; i < 3; ++i) {
			u[i] = u_[i];// .cast<float>();
			b[i] = b_[i];// .cast<float>();
		}
		//ground_truth
		for (int i = 0; i < 3; ++i) {
			ground_truth[i] = collision_free_cloth_llt.solve(b[i]);
		}
		std::vector<VectorXd> u_use;

		std::string txt_file_name = "add force iteration_result" + std::to_string(k);

		//jacobi
		u_use = u;
		std::vector<double> jacobi_relative_error;
		std::vector<int> jacobi_time;
		double conv_rate_2 = 1e-7;
		conv_rate_2 *= conv_rate_2;
		jacobi(u_use, b, system_matrix, ground_truth, jacobi_relative_error, conv_rate_2, jacobi_time);
		WriteTxt::writeTxt(jacobi_relative_error, jacobi_time, 10, txt_file_name, "jacobi");
		std::cout << "jaco" << std::endl;
		////super_jacobi
		u_use = u;
		std::vector<double> super_jacobi_relative_error;
		std::vector<int> super_jacobi_time;
		superJacobi(u_use, b, system_matrix, ground_truth, super_jacobi_relative_error, conv_rate_2, 2, super_jacobi_time);
		WriteTxt::addToTxt(super_jacobi_relative_error, super_jacobi_time, 10, txt_file_name, "super jacobi");
		//WriteTxt::writeTxt(super_jacobi_relative_error, super_jacobi_time, 10, txt_file_name, "super jacobi");
		std::cout << "finished a_jacobi 2" << std::endl;
		u_use = u;
		std::vector<double> super_jacobi_relative_error_3;
		std::vector<int> super_jacobi_time_3;
		superJacobi(u_use, b, system_matrix, ground_truth, super_jacobi_relative_error_3, conv_rate_2, 3, super_jacobi_time_3);
		WriteTxt::addToTxt(super_jacobi_relative_error_3, super_jacobi_time_3, 10, txt_file_name, "super jacobi 3");

		std::cout << "finished a_jacobi 3" << std::endl;


		////gauss_seidel
		u_use = u;
		std::vector<double> gauss_seidel_relative_error;
		std::vector<int> gauss_seidel_time;
		gauss_seidel(u_use, b, system_matrix, ground_truth, gauss_seidel_relative_error, conv_rate_2, gauss_seidel_time);
		WriteTxt::addToTxt(gauss_seidel_relative_error, gauss_seidel_time, 10, txt_file_name, "gauss seidel");
		////jacobi_chebyshev
		u_use = u;
		std::vector<double> jacobi_chebyshev_relative_error;
		std::vector<int> jacobi_chebyshev_time;
		chebyshevSemiIterativeJacobi(u_use, b, system_matrix, ground_truth, jacobi_chebyshev_relative_error, conv_rate_2, jacobi_chebyshev_time);
		WriteTxt::addToTxt(jacobi_chebyshev_relative_error, jacobi_chebyshev_time, 10, txt_file_name, "cheybyshev jacobi");
		//super_jacobi_chebyshev_2
		u_use = u;
		std::vector<double> super_jacobi_chebyshev_relative_error_2;
		std::vector<int> super_jacobi_chebyshev_time_2;
		chebyshevSemiIterativeSuperJacobi(u_use, b, system_matrix, ground_truth, super_jacobi_chebyshev_relative_error_2, conv_rate_2, 2, super_jacobi_chebyshev_time_2);
		WriteTxt::addToTxt(super_jacobi_chebyshev_relative_error_2, super_jacobi_chebyshev_time_2, 10, txt_file_name, "cheybyshev super jacobi 2");
		//WriteTxt::writeTxt(super_jacobi_chebyshev_relative_error_2, super_jacobi_chebyshev_time_2, 10, txt_file_name, "cheybyshev super jacobi 2");
		//super_jacobi_chebyshev_3
		u_use = u;
		std::vector<double> super_jacobi_chebyshev_relative_error_3;
		std::vector<int> super_jacobi_chebyshev_time_3;
		chebyshevSemiIterativeSuperJacobi(u_use, b, system_matrix, ground_truth, super_jacobi_chebyshev_relative_error_3, conv_rate_2, 3, super_jacobi_chebyshev_time_3);
		WriteTxt::addToTxt(super_jacobi_chebyshev_relative_error_3, super_jacobi_chebyshev_time_3, 10, txt_file_name, "cheybyshev super jacobi 3");
		//gauss_seidel_chebyshev
		std::cout << "start" << std::endl;
		u_use = u;
		std::vector<double> gauss_seidel_chebyshev_relative_error;
		std::vector<int> gauss_seidel_chebyshev_time;
		chebyshev_gauss_seidel(u_use, b, system_matrix, ground_truth, gauss_seidel_chebyshev_relative_error, conv_rate_2, gauss_seidel_chebyshev_time);
		WriteTxt::addToTxt(gauss_seidel_chebyshev_relative_error, gauss_seidel_chebyshev_time, 10, txt_file_name, "gauss seidel chebyshev 3");
		std::cout << "end" << std::endl;
		//std::cout << gauss_seidel_chebyshev_relative_error.size() << " ";
		//std::cout << "finished GS cheby" << std::endl;
		////PCG
		u_use = u;
		std::vector<double> PCG_relative_error;
		std::vector<int> PCG_time;
		PCG(u, b, system_matrix, ground_truth, PCG_relative_error, conv_rate_2, PCG_time);
		WriteTxt::addToTxt(PCG_relative_error, PCG_time, 10, txt_file_name, "PCG");
	}
}


void IterationMethod::test()
{
	//testGaussSeidel();
	//load matrix
	std::vector<VectorXd> u_(3);
	std::vector<VectorXd> b_(3);
	SparseMatrix<double, ColMajor> system_matrix_;
	SparseMatrix<double, ColMajor> system_matrix_3;
	std::vector<VectorXd> ground_truth(3);

	//dimension_per_thread_test.resize(thread->thread_num + 1, 3);
	//for (int i = 0; i < 3; ++i) {
	//	dimension_per_thread_test[i] = i;
	//}
	std::string file_name_matrix="./back/scene2_101/global.dat";

	EigenMatrixIO::read_sp_binary(file_name_matrix.c_str(), system_matrix_);
	SparseMatrix<double, ColMajor> system_matrix;


	system_matrix = system_matrix_;// .cast<float>();

	std::vector<std::string> file_name_u(3);
	std::vector<std::string> file_name_b(3);

	for (int i = 0; i < 3; ++i) {
		file_name_u[i] = "./back/scene2_101/u" + std::to_string(i) + ".dat";
		file_name_b[i] = "./back/scene2_101/b" + std::to_string(i) + ".dat";
		EigenMatrixIO::read_binary(file_name_u[i].c_str(), u_[i]);
		EigenMatrixIO::read_binary(file_name_b[i].c_str(), b_[i]);
		std::cout << "i" << std::endl;
	}	

	std::vector<VectorXd> u(3);
	std::vector<VectorXd> b(3);
	for (int i = 0; i < 3; ++i) {
		u[i] = u_[i];// .cast<float>();
		b[i] = b_[i];// .cast<float>();
	}
	//ground_truth
	SimplicialLLT<SparseMatrix<double>> collision_free_cloth_llt;
	collision_free_cloth_llt.compute(system_matrix);
	for (int i = 0; i < 3; ++i) {
		ground_truth[i] = collision_free_cloth_llt.solve(b[i]);
	}	
	std::vector<VectorXd> u_use;
	//jacobi
	u_use = u;
	std::vector<double> jacobi_relative_error;
	std::vector<int> jacobi_time;
	double conv_rate_2 = 5e-7;
	conv_rate_2 *= conv_rate_2;
	jacobi(u_use, b, system_matrix, ground_truth, jacobi_relative_error, conv_rate_2,jacobi_time);
	WriteTxt::writeTxt(jacobi_relative_error, jacobi_time, 10, "iteration_result", "jacobi");
	////super_jacobi
	u_use = u;
	std::vector<double> super_jacobi_relative_error;
	std::vector<int> super_jacobi_time;
	superJacobi(u_use, b, system_matrix, ground_truth, super_jacobi_relative_error, conv_rate_2,2, super_jacobi_time);
	WriteTxt::addToTxt(super_jacobi_relative_error, super_jacobi_time, 10, "iteration_result", "super jacobi");
	//WriteTxt::writeTxt(super_jacobi_relative_error, super_jacobi_time, 10, "iteration_result", "super jacobi");
	//std::cout << "finished a_jacobi 2" << std::endl;
	u_use = u;
	std::vector<double> super_jacobi_relative_error_3;
	std::vector<int> super_jacobi_time_3;
	superJacobi(u_use, b, system_matrix, ground_truth, super_jacobi_relative_error_3, conv_rate_2, 3, super_jacobi_time_3);
	WriteTxt::addToTxt(super_jacobi_relative_error_3, super_jacobi_time_3, 10, "iteration_result", "super jacobi 3");

	//gauss_seidel_chebyshev
	

	////gauss_seidel
	u_use = u;
	std::vector<double> gauss_seidel_relative_error;
	std::vector<int> gauss_seidel_time;
	gauss_seidel(u_use, b, system_matrix, ground_truth, gauss_seidel_relative_error, conv_rate_2, gauss_seidel_time);
	WriteTxt::addToTxt(gauss_seidel_relative_error, gauss_seidel_time, 10, "iteration_result", "gauss seidel");
	////jacobi_chebyshev
	u_use = u;
	std::vector<double> jacobi_chebyshev_relative_error;
	std::vector<int> jacobi_chebyshev_time;
	chebyshevSemiIterativeJacobi(u_use, b, system_matrix, ground_truth, jacobi_chebyshev_relative_error, conv_rate_2, jacobi_chebyshev_time);
	WriteTxt::addToTxt(jacobi_chebyshev_relative_error, jacobi_chebyshev_time, 10, "iteration_result", "cheybyshev jacobi");
	//super_jacobi_chebyshev_2
	u_use = u;
	std::vector<double> super_jacobi_chebyshev_relative_error_2;
	std::vector<int> super_jacobi_chebyshev_time_2;
	chebyshevSemiIterativeSuperJacobi(u_use, b, system_matrix, ground_truth, super_jacobi_chebyshev_relative_error_2, conv_rate_2,2, super_jacobi_chebyshev_time_2);
	WriteTxt::addToTxt(super_jacobi_chebyshev_relative_error_2, super_jacobi_chebyshev_time_2, 10, "iteration_result", "cheybyshev super jacobi 2");
	//WriteTxt::writeTxt(super_jacobi_chebyshev_relative_error_2, super_jacobi_chebyshev_time_2, 10, "iteration_result", "cheybyshev super jacobi 2");
	//super_jacobi_chebyshev_3
	u_use = u;
	std::vector<double> super_jacobi_chebyshev_relative_error_3;
	std::vector<int> super_jacobi_chebyshev_time_3;
	chebyshevSemiIterativeSuperJacobi(u_use, b, system_matrix, ground_truth, super_jacobi_chebyshev_relative_error_3, conv_rate_2, 3, super_jacobi_chebyshev_time_3);
	WriteTxt::addToTxt(super_jacobi_chebyshev_relative_error_3, super_jacobi_chebyshev_time_3, 10, "iteration_result", "cheybyshev super jacobi 3");
	
	u_use = u;
	std::vector<double> gauss_seidel_chebyshev_relative_error;
	std::vector<int> gauss_seidel_chebyshev_time;
	chebyshev_gauss_seidel(u_use, b, system_matrix, ground_truth, gauss_seidel_chebyshev_relative_error, conv_rate_2, gauss_seidel_chebyshev_time);
	WriteTxt::addToTxt(gauss_seidel_chebyshev_relative_error, gauss_seidel_chebyshev_time, 10, "iteration_result", "gauss seidel chebyshev 3");
	std::cout << gauss_seidel_chebyshev_relative_error.size() << " ";
	std::cout << "finished GS cheby" << std::endl;
	////PCG
	u_use = u;
	std::vector<double> PCG_relative_error;
	std::vector<int> PCG_time;
	PCG(u, b, system_matrix, ground_truth, PCG_relative_error, conv_rate_2, PCG_time);
	WriteTxt::addToTxt(PCG_relative_error, PCG_time, 10, "iteration_result", "PCG");
	
	//size_t max_itr = 0;
	//max_itr = (std::max)(max_itr, jacobi_relative_error.size());
	//max_itr = (std::max)(max_itr, super_jacobi_relative_error.size());
	//max_itr = (std::max)(max_itr, super_jacobi_relative_error_3.size());
	//max_itr = (std::max)(max_itr, gauss_seidel_relative_error.size());
	//max_itr = (std::max)(max_itr, jacobi_chebyshev_relative_error.size());
	//max_itr = (std::max)(max_itr, super_jacobi_chebyshev_relative_error_2.size());
	//max_itr = (std::max)(max_itr, super_jacobi_chebyshev_relative_error_3.size());
	//max_itr = (std::max)(max_itr, gauss_seidel_chebyshev_relative_error.size());
	//max_itr = (std::max)(max_itr, PCG_relative_error.size());
	std::cout << "1.jacobi 2.super_jacobi_2 3. super_jacobi_3 4.gauss_seidel 5.jacobi_chebyshev 6.super_jacobi_chebyshev_2 7.super_jacobi_chebyshev_3 8.gauss_seidel_chebyshev 9.PCG" << std::endl;
	//for (int i = 0; i < max_itr; ++i) {
	//	if (i < jacobi_relative_error.size()) {
	//		std::cout << (jacobi_relative_error[i])<<" "<< jacobi_time[i];
	//	}
	//	else {
	//		std::cout << " ++ ++";
	//	}
	//	if (i < super_jacobi_relative_error.size()) {
	//		std::cout << " " << (super_jacobi_relative_error[i])<<" "<<super_jacobi_time[i];
	//	}
	//	else {
	//		std::cout << " ++ ++";
	//	}
	//	if (i < super_jacobi_relative_error_3.size()) {
	//		std::cout << " " << (super_jacobi_relative_error_3[i])<<" "<< super_jacobi_time_3[i];
	//	}
	//	else {
	//		std::cout << " ++ ++";
	//	}
	//	if (i < gauss_seidel_relative_error.size()) {
	//		std::cout << " " << (gauss_seidel_relative_error[i])<<" "<<gauss_seidel_time[i];
	//	}
	//	else {
	//		std::cout << " ++ ++";
	//	}
	//	if (i < jacobi_chebyshev_relative_error.size()) {
	//		std::cout << " " << (jacobi_chebyshev_relative_error[i])<<" "<<jacobi_chebyshev_time[i];
	//	}
	//	else {
	//		std::cout << " ++ ++";
	//	}
	//	if (i < super_jacobi_chebyshev_relative_error_2.size()) {
	//		std::cout << " " << (super_jacobi_chebyshev_relative_error_2[i])<<" "<<super_jacobi_chebyshev_time_2[i];
	//	}
	//	else {
	//		std::cout << " ++ ++";
	//	}
	//	if (i < super_jacobi_chebyshev_relative_error_3.size()) {
	//		std::cout << " " << (super_jacobi_chebyshev_relative_error_3[i])<<" "<<super_jacobi_chebyshev_time_3[i];
	//	}
	//	else {
	//		std::cout << " ++ ++";
	//	}
	//	if (i < gauss_seidel_chebyshev_relative_error.size()) {
	//		std::cout << " " << (gauss_seidel_chebyshev_relative_error[i])<<" "<<gauss_seidel_chebyshev_time[i];
	//	}
	//	else {
	//		std::cout << " ++ ++";
	//	}
	//	if (i < PCG_relative_error.size()) {
	//		std::cout << " " << (PCG_relative_error[i])<<" "<<PCG_time[i];
	//	}
	//	else {
	//		std::cout << " ++ ++";
	//	}
	//	std::cout << std::endl;
	//}
	std::cout << gauss_seidel_chebyshev_relative_error.size() << std::endl;
}

//void IterationMethod::jacobi(VectorXd& u, VectorXd& b, SparseMatrix<double, RowMajor>& system_matrix, VectorXd& ground_truth, std::vector<double>&relative_error)



