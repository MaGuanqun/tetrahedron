#include"iteration_method.h"
#include"basic/EigenMatrixIO.h"
#include<algorithm>
#include"basic/write_txt.h"



//void IterationMethod::useOperator()
//{
//
//}


void IterationMethod::testOperator(A_JacobiOperator* A_jacobi_operator, SparseMatrix<double, ColMajor>& R_jacobi)
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
void IterationMethod::RMultiXPlusb(A_JacobiOperator* A_jacobi_operator, double* x, double* b, double* result)
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



//COMPUTE_RESIDUAL
void IterationMethod::computeResidual(int* global_matrix_vertex_index, double* global_matrix_coefficient, int* global_matrix_vertex_index_start,
	VectorXd* x, VectorXd* b, VectorXd* b_, double* residual_norm, int vertex_index_begin, int vertex_index_end)
{
	double temp_diff;
	double temp_value;
	int vertex_end;
	(*residual_norm) = 0.0;
	double* x_dimension;
	double* b_dimension;
	for (int k = 0; k < 3; ++k) {
		x_dimension = x[k].data();
		b_dimension = b[k].data();
		for (int i = vertex_index_begin; i < vertex_index_end; ++i)	{
			vertex_end = global_matrix_vertex_index_start[i + 1];
			temp_diff = 0;
			for (int j = global_matrix_vertex_index_start[i]; j < vertex_end; ++j) {
				temp_diff += global_matrix_coefficient[j] * x_dimension[global_matrix_vertex_index[j]];
			}
			temp_diff = b_dimension[i] - temp_diff;
			(*residual_norm) += temp_diff * temp_diff;
		}
	}
	//std::cout << (*residual_norm) << std::endl;
	
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




VectorXd IterationMethod::RMultiX(A_JacobiOperator* A_jacobi_operator, VectorXd& x)
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


void IterationMethod::createSuperJacobiOperator(A_JacobiOperator* A_jacobi_operator, SparseMatrix<double, ColMajor>& R_jacobi,
	A_JacobiOperator* A_jacobi_basic)
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
		inner_size = R_jacobi.outerIndexPtr()[i + 1] - R_jacobi.outerIndexPtr()[i];
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


void IterationMethod::createSuperJacobiOperator(A_JacobiOperator*  A_jacobi_operator, SparseMatrix<double, RowMajor>& R_jacobi,
	A_JacobiOperator* A_jacobi_basic)
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


void IterationMethod::createAJacobiOperator(std::vector<std::array<int, 2>>& coeff_pos, std::vector<double>& coeff)
{
	BasicJacobiOperator A_jacobi_basic; //column major
	createAJacobiOperator(&A_jacobi_operator, coeff_pos, coeff, &A_jacobi_basic);

	setRJaocbiDiagonalInv(&A_jacobi_operator, &off_diagonal_operator, diagonal_inv, original_diagonal, original_diagonal_inv,&global_matrix_operator);
	original_initial_diagonal = original_diagonal;

	//create structure of 2-order A jacobi
	createHighOrderAJacobiMethod(&A_jacobi_basic, &A_jacobi_operator, &A_jacobi_operator_2, &A_jacobi_operator);
	arrangeIndex(thread->thread_num, A_jacobi_operator_2.vertex_index.size(), A_jacobi_2_index_per_thread);

	//create structure of 3-order A jacobi
	createHighOrderAJacobiMethod(&A_jacobi_basic, &A_jacobi_operator_2, &A_jacobi_operator_3, &A_jacobi_operator);
	arrangeIndex(thread->thread_num, A_jacobi_operator_3.vertex_index.size(), A_jacobi_3_index_per_thread);

	thread->assignTask(this, UPDATE_JACOBI_OPERATOR);
	thread->assignTask(this, UPDATE_2_A_JACOBI_ITR_MATRIX);
	thread->assignTask(this, UPDATE_3_A_JACOBI_ITR_MATRIX);
	
	ori_A_jacobi_operator_2_coefficient = A_jacobi_operator_2.coefficient;
	ori_A_jacobi_operator_3_coefficient = A_jacobi_operator_3.coefficient;
	ori_A_jacobi_operator_coefficient = A_jacobi_operator.coefficient;

	//AJacobiOperator A2_jacobi_operator;
	//BasicJacobiOperator A2_jacobi_operator_;
	////SparseMatrix<double, RowMajor> global_mat_2 = (*global_mat) * (*global_mat);
	//SparseMatrix<double, RowMajor> R_jacobi_2 = R_Jacobi * R_Jacobi * R_Jacobi;
	////std::cout << global_mat_2.nonZeros() << " " << R_jacobi_2.nonZeros() << std::endl;
	//createAJacobiOperator(&A2_jacobi_operator, R_jacobi_2, &A2_jacobi_operator_);
	//testIfOperatorIsRight(&A_jacobi_operator_3, &A2_jacobi_operator);
}


void IterationMethod::setOperatorCollisionFree()
{
	A_jacobi_operator_2.coefficient = ori_A_jacobi_operator_2_coefficient;
	A_jacobi_operator_3.coefficient = ori_A_jacobi_operator_3_coefficient;
	A_jacobi_operator.coefficient = ori_A_jacobi_operator_coefficient;
}

void IterationMethod::initialRecordDiagonal_Operator()
{
	original_diagonal = original_initial_diagonal;
	for (int i = 0; i < sys_size; ++i) {
		original_diagonal_inv[i] = 1.0 / original_diagonal[i];
	}
	diagonal_inv = original_diagonal_inv;
	thread->assignTask(this, UPDATE_JACOBI_OPERATOR);
	thread->assignTask(this, UPDATE_2_A_JACOBI_ITR_MATRIX);
	thread->assignTask(this, UPDATE_3_A_JACOBI_ITR_MATRIX);
	ori_A_jacobi_operator_2_coefficient = A_jacobi_operator_2.coefficient;
	ori_A_jacobi_operator_3_coefficient = A_jacobi_operator_3.coefficient;
	ori_A_jacobi_operator_coefficient = A_jacobi_operator.coefficient;
	for (int i = 0; i < sys_size; ++i) {
		global_matrix_operator.coefficient[global_matrix_operator.index_of_diagonal[i]] = original_diagonal[i];
	}
}


void IterationMethod::updateDiagonalWithAnchorVertices(int anchor_vertex_size, int* anchor_vertex, int system_size, int vertex_index_start, double position_stiffness)
{
	memcpy(original_diagonal.data() + vertex_index_start, original_initial_diagonal.data() + vertex_index_start, 8 * system_size);
	for (int i = 0; i < anchor_vertex_size; ++i) {
		original_diagonal[anchor_vertex[i] + vertex_index_start] += position_stiffness;
	}
	for (int i = 0; i < system_size; ++i) {
		original_diagonal_inv[i + vertex_index_start] = 1.0 / original_diagonal[i + vertex_index_start];
		global_matrix_operator.coefficient[global_matrix_operator.index_of_diagonal[i + vertex_index_start]] = original_diagonal[i + vertex_index_start];
	}
}

void IterationMethod::updateDiagonalWithAnchorVerticesTotal()
{
	diagonal_inv = original_diagonal_inv;
	thread->assignTask(this, UPDATE_JACOBI_OPERATOR);
	thread->assignTask(this, UPDATE_2_A_JACOBI_ITR_MATRIX);
	thread->assignTask(this, UPDATE_3_A_JACOBI_ITR_MATRIX);
	ori_A_jacobi_operator_2_coefficient = A_jacobi_operator_2.coefficient;
	ori_A_jacobi_operator_3_coefficient = A_jacobi_operator_3.coefficient;
	ori_A_jacobi_operator_coefficient = A_jacobi_operator.coefficient;
}

void IterationMethod::updateClothDiagonalPerThread(int index_start, int index_end, int obj_index_start, double* stiffness, bool* need_update)
{
	double* ori_diagonal = original_diagonal.data();
	double* current_diagonal_inv = diagonal_inv.data();
	double* ori_diagonal_inv = original_diagonal_inv.data();
	double* global_mat_coeff = global_matrix_operator.coefficient.data();
	int* global_mat_digonal_index = global_matrix_operator.index_of_diagonal.data();


	for (int i = index_start; i < index_end; ++i) {
		if (need_update[i]) {
			global_mat_coeff[global_mat_digonal_index[i + obj_index_start]] = ori_diagonal[i + obj_index_start] + stiffness[i];
			current_diagonal_inv[i + obj_index_start] = 1.0 / (ori_diagonal[i + obj_index_start] + stiffness[i]);

		}
		else {
			global_mat_coeff[global_mat_digonal_index[i + obj_index_start]] = ori_diagonal[i + obj_index_start];
			current_diagonal_inv[i + obj_index_start] = ori_diagonal_inv[i + obj_index_start];
		}
	}
}

void IterationMethod::updateTetDiagonalPerThread(int index_start, int index_end, int obj_index_start, double* stiffness, bool* need_update, unsigned int* vertex_index_on_surface)
{
	double* ori_diagonal = original_diagonal.data();
	double* current_diagonal_inv = diagonal_inv.data();
	double* ori_diagonal_inv = original_diagonal_inv.data();
	double* global_mat_coeff = global_matrix_operator.coefficient.data();
	int* global_mat_digonal_index = global_matrix_operator.index_of_diagonal.data();

	int vertex_index;
	for (int i = index_start; i < index_end; ++i) {
		vertex_index = vertex_index_on_surface[i] + obj_index_start;
		if (need_update[vertex_index]) {
			global_mat_coeff[global_mat_digonal_index[vertex_index]] = ori_diagonal[vertex_index] + stiffness[vertex_index];
			current_diagonal_inv[vertex_index] = 1.0 / (ori_diagonal[vertex_index] + stiffness[vertex_index]);

		}
		else {
			global_mat_coeff[global_mat_digonal_index[vertex_index]] = ori_diagonal[vertex_index];
			current_diagonal_inv[vertex_index] = ori_diagonal_inv[vertex_index];
		}
	}
}

void IterationMethod::testIfOperatorIsRight(AJacobiOperator* A_jacobi_operator, AJacobiOperator* A_jacobi_operator_)
{
	if (A_jacobi_operator->vertex_index.size() != A_jacobi_operator_->vertex_index.size()) {
		std::cout << "error size not equal" << std::endl;
		return;
	}

	for (int i = 0; i < sys_size; ++i) {
		if (A_jacobi_operator->start_index[i + 1] != A_jacobi_operator_->start_index[i + 1]) {
			std::cout << "start_index_not_equal" << std::endl;
			return;
		}
	}

	for (int i = 0; i < A_jacobi_operator->vertex_index.size(); ++i) {
		if (A_jacobi_operator->vertex_index[i] != A_jacobi_operator_->vertex_index[i]) {
			std::cout << "vertex_index_not_equal" << std::endl;
			return;
		}
	}
	for (int i = 0; i < A_jacobi_operator->vertex_index.size(); ++i) {
		if (A_jacobi_operator->coefficient[i] != A_jacobi_operator_->coefficient[i]) {
			std::cout << A_jacobi_operator->coefficient[i] << " " << A_jacobi_operator_->coefficient[i] << std::endl;
		}
	}
}


void IterationMethod::testIfOperatorIsRight(AJacobiOperator* A_jacobi_operator, BasicJacobiOperator* A_jacobi_operator_)
{
	int count = 0;
	for (int i = 0; i < A_jacobi_operator_->element.size(); ++i) {
		if (count != A_jacobi_operator->start_index[i]) {
			std::cout << "start_index_error " << std::endl;
		}
		for (int j = 0; j < A_jacobi_operator_->element[i].size(); ++j) {
			if (A_jacobi_operator_->element[i][j].vertex_index != A_jacobi_operator->vertex_index[count]) {
				std::cout << "index error" << std::endl;
			}
			if (A_jacobi_operator_->element[i][j].coeff != A_jacobi_operator->coefficient[count]) {
				std::cout << "coeff error" << std::endl;
			}
			count++;
		}
	}

}

void IterationMethod::createAJacobiOperator(AJacobiOperator* A_jacobi_operator, SparseMatrix<double, RowMajor>& R_jacobi,
	BasicJacobiOperator* A_jacobi_basic)
{
	int sys_size = R_jacobi.rows();
	A_jacobi_operator->vertex_index.reserve(R_jacobi.nonZeros());
	A_jacobi_operator->coefficient.reserve(R_jacobi.nonZeros());
	A_jacobi_operator->start_index.reserve(R_jacobi.rows()+1);
	A_jacobi_basic->element.resize(sys_size);
	int inner_size;
	int inner_index_start;

	for (int i = 0; i < sys_size; ++i) {
		inner_size = R_jacobi.outerIndexPtr()[i + 1] - R_jacobi.outerIndexPtr()[i];
		A_jacobi_basic->element[i].reserve(inner_size);
	}
	A_jacobi_operator->start_index.push_back(0);
	for (int i = 0; i < sys_size; ++i) {
		inner_size = R_jacobi.outerIndexPtr()[i + 1] - R_jacobi.outerIndexPtr()[i];
		inner_index_start = R_jacobi.outerIndexPtr()[i];
		A_jacobi_operator->start_index.push_back(A_jacobi_operator->start_index[i] + inner_size);
		for (int j = 0; j < inner_size; ++j) {
			A_jacobi_operator->vertex_index.push_back(R_jacobi.innerIndexPtr()[inner_index_start + j]);
			A_jacobi_operator->coefficient.push_back(R_jacobi.valuePtr()[inner_index_start + j]);

			A_jacobi_basic->element[R_jacobi.innerIndexPtr()[inner_index_start + j]].push_back(
				ColIndexWithCoeff(i, R_jacobi.valuePtr()[inner_index_start + j]));
		}
	}
}

//UPDATE_JACOBI_OPERATOR
void IterationMethod::updateJacobiOperator(int thread_id)
{
	int vertex_end = vertex_index_begin_thread[thread_id+1];
	double diagonal_coe;
	double* coeff = A_jacobi_operator.coefficient.data();

	int end_index;

	for (int i = vertex_index_begin_thread[thread_id]; i < vertex_end; ++i) {
		diagonal_coe = diagonal_inv.data()[i];
		end_index = A_jacobi_operator.start_index[i + 1];
		for (int j = A_jacobi_operator.start_index[i]; j < end_index; ++j) {
			A_jacobi_operator.coefficient[j] = diagonal_coe* off_diagonal_operator.coefficient[j];
		}		
	}
}


void IterationMethod::setRJaocbiDiagonalInv(AJacobiOperator* A_jacobi_operator, AJacobiOperator* off_diagonal_operator, 
	std::vector<double>& diagonal_inv, std::vector<double>& ori_diagonal, std::vector<double>& ori_diagonal_inv, AJacobiOperator* global_matrix_operator)
{
	global_matrix_operator->vertex_index = A_jacobi_operator->vertex_index;
	global_matrix_operator->start_index = A_jacobi_operator->start_index;
	global_matrix_operator->coefficient = A_jacobi_operator->coefficient;
	global_matrix_operator->index_of_diagonal.resize(sys_size);

	(*off_diagonal_operator) = (*A_jacobi_operator);
	int sys_size = off_diagonal_operator->start_index.size() - 1;
	diagonal_inv.resize(sys_size);
	ori_diagonal.resize(sys_size);
	ori_diagonal_inv.resize(sys_size);
	int col_end;
	for (int j = 0; j < sys_size; ++j) {
		col_end = off_diagonal_operator->start_index[j + 1];
		for (int i = off_diagonal_operator->start_index[j]; i < col_end; ++i) {
			if (off_diagonal_operator->vertex_index[i] == j) {
				ori_diagonal[j] = off_diagonal_operator->coefficient[i];
				diagonal_inv[j] = 1.0 / ori_diagonal[j];
				ori_diagonal_inv[j] = diagonal_inv[j];
				off_diagonal_operator->coefficient[i] = 0.0;
				global_matrix_operator->index_of_diagonal[j] = i;
			}
			else {
				off_diagonal_operator->coefficient[i] *= -1.0;
			}
		}
	}
}

void IterationMethod::createAJacobiOperator(AJacobiOperator* A_jacobi_operator, std::vector<std::array<int,2>>& coeff_pos, std::vector<double>& coeff,
	BasicJacobiOperator* A_jacobi_basic)
{
	std::map<AJacobiOperatorForConstruct, double> system;
	for (int i = 0; i < coeff_pos.size(); ++i) {
		buildMap(system, coeff_pos[i][0], coeff_pos[i][1], coeff[i]);
	}
	int j = 0;
	int index0_pre[2];
	int count = 0;
	//AJacobiOperatorForConstruct a1();
	A_jacobi_operator->vertex_index.reserve(system.size());
	A_jacobi_operator->coefficient.reserve(system.size());
	A_jacobi_operator->start_index.reserve(sys_size + 1);
	A_jacobi_operator->start_index.push_back(0);

	A_jacobi_basic->element.resize(sys_size);

	std::vector<int>* vertex_index;
	std::vector<int>* start_index;
	std::vector<double>* coefficient;


	vertex_index = &A_jacobi_operator->vertex_index;
	coefficient = &A_jacobi_operator->coefficient;
	start_index = &A_jacobi_operator->start_index;

	for (auto i = system.begin(); i != system.end(); ++i) {
		if (i->first.index[0] != j) {			
			j++;
			start_index->push_back(count);
		}
		vertex_index->push_back(i->first.index[1]);
		coefficient->push_back(i->second);
		count++;
	}
	start_index->push_back(count);

	for (int i = 0; i < sys_size; ++i) {
		A_jacobi_basic->element[j].reserve(start_index->data()[j+1] - start_index->data()[j]);
	}
	std::vector<ColIndexWithCoeff>* element;
	element = A_jacobi_basic->element.data();
	for (auto i = system.begin(); i != system.end(); ++i) {
		element[i->first.index[1]].push_back(ColIndexWithCoeff(i->first.index[0], i->second));
	}
	//for (int i = 0; i < sys_size; ++i) {
	//	std::sort(element[i].begin(), element[i].end());
	//}
}


void IterationMethod::buildMap(std::map<AJacobiOperatorForConstruct, double>& system, int v0, int v1, double coeff)
{
	AJacobiOperatorForConstruct A(v0, v1);
	auto ret = system.insert({ A,coeff });
	if (!ret.second) {
		ret.first->second+=coeff;
	}
}

//A_jacobi_operator_need_to_multi * A_jacobi_operator_basic
//A_jacobi_operator_right is the row major version of A_jacobi_operator_basic
void IterationMethod::createHighOrderAJacobiMethod(BasicJacobiOperator* A_jacobi_operator_basic, AJacobiOperator* A_jacobi_operator_need_to_multi_,
	AJacobiOperator* A_jacobi_operator_result, AJacobiOperator* A_jacobi_operator_right)
{

	BasicJacobiOperator A_jacobi_operator_need_to_multi;
	transferAJacobiOperator2BasicOperator(A_jacobi_operator_need_to_multi_, &A_jacobi_operator_need_to_multi);

	BasicJacobiOperator A_jacobi_operator;
	std::vector<bool> is_vertex_used(sys_size, false);
	std::vector<int> indicate_if_vertex_exists(sys_size, -1);
	A_jacobi_operator.element.resize(sys_size);

	std::vector<ColIndexWithCoeff>* need_to_multi_element;

	int basic_index_in_need_to_multi;

	std::vector<ColIndexWithCoeff>* basic_element;

	int basic_index_size;


	std::vector<ColIndexWithCoeff>* result_element;
	std::vector<std::vector<std::array<int, 4>>>* matrix_multiple_index;

	int end_index_of_col;

	double element_value;

	int column_index;
	int need_to_multi_index_size;

	std::vector<int> multiple_index;

	for (int i = 0; i < sys_size; ++i) {		
		need_to_multi_element = &A_jacobi_operator_need_to_multi.element[i];
		need_to_multi_index_size = need_to_multi_element->size();
		result_element = &A_jacobi_operator.element[i];
		result_element->reserve(3 * need_to_multi_index_size);
		for (int j = 0; j < need_to_multi_index_size; ++j) {
			basic_index_in_need_to_multi = need_to_multi_element->data()[j].vertex_index;
			basic_element = &A_jacobi_operator_basic->element[basic_index_in_need_to_multi];		
			basic_index_size = basic_element->size();
			for (int k = 0; k < basic_index_size; ++k) {
				column_index = basic_element->data()[k].vertex_index;
				if (!is_vertex_used[column_index]) {
					is_vertex_used[column_index] = true;
					multiple_index.clear();
					element_value = obtainElementValue(need_to_multi_element, &A_jacobi_operator_basic->element[column_index],indicate_if_vertex_exists,
						multiple_index,i, column_index);
					result_element->push_back(ColIndexWithCoeff(column_index, element_value, multiple_index));
				}
			}
		}

		for (int j = 0; j < result_element->size(); ++j) {
			is_vertex_used[result_element->data()[j].vertex_index] = false;
		}
	}

	for (int i = 0; i < sys_size; ++i)//reorder
	{
		std::sort(A_jacobi_operator.element[i].begin(), A_jacobi_operator.element[i].end());
	}

	transferBasicOperator2AJacobi(A_jacobi_operator_result, &A_jacobi_operator, A_jacobi_operator_need_to_multi_, A_jacobi_operator_right);
}


//UPDATE_2_A_JACOBI_ITR_MATRIX
//AJacobiOperator* A_jacobi_operator_result, AJacobiOperator* left_operator, AJacobiOperator* right_operator
//A_jacobi_operator_2 = A_jacobi_operator* A_jacobi_operator
void IterationMethod::update2AJaocbiIterationMatrix(int thread_id)
{
	double* coefficient = A_jacobi_operator_2.coefficient.data();
	int* left_multiplier_index = A_jacobi_operator_2.left_multiplier_index.data();
	int* right_multiplier_index = A_jacobi_operator_2.right_multiplier_index.data();
	int* multiplier_start_per_element = A_jacobi_operator_2.multiplier_start_per_element.data();


	double* left_coefficient = A_jacobi_operator.coefficient.data();
	double* right_coefficient = A_jacobi_operator.coefficient.data();

	double value;
	int multiplier_end;
	int vertex_end = A_jacobi_2_index_per_thread[thread_id + 1];
	for (int i = A_jacobi_2_index_per_thread[thread_id]; i < vertex_end; ++i) { //row index=i;
		value = 0.0;
		multiplier_end = multiplier_start_per_element[i + 1];
		for (int j = multiplier_start_per_element[i]; j < multiplier_end; ++j) {
			value += left_coefficient[left_multiplier_index[j]] * right_coefficient[right_multiplier_index[j]];
		}
		coefficient[i] = value;
	}
}


//UPDATE_3_A_JACOBI_ITR_MATRIX
//A_jacobi_operator_3 = A_jacobi_operator_2* A_jacobi_operator
void IterationMethod::update3AJaocbiIterationMatrix(int thread_id)
{
	double* coefficient = A_jacobi_operator_3.coefficient.data();
	int* left_multiplier_index = A_jacobi_operator_3.left_multiplier_index.data();
	int* right_multiplier_index = A_jacobi_operator_3.right_multiplier_index.data();
	int* multiplier_start_per_element = A_jacobi_operator_3.multiplier_start_per_element.data();


	double* left_coefficient = A_jacobi_operator_2.coefficient.data();
	double* right_coefficient = A_jacobi_operator.coefficient.data();

	double value;
	int multiplier_end;
	int vertex_end = A_jacobi_3_index_per_thread[thread_id + 1];
	for (int i = A_jacobi_3_index_per_thread[thread_id]; i < vertex_end; ++i) { //row index=i;
		value = 0.0;
		multiplier_end = multiplier_start_per_element[i + 1];
		for (int j = multiplier_start_per_element[i]; j < multiplier_end; ++j) {
			value += left_coefficient[left_multiplier_index[j]] * right_coefficient[right_multiplier_index[j]];
		}
		coefficient[i] = value;
	}
}



void IterationMethod::transferAJacobiOperator2BasicOperator(AJacobiOperator* A_jacobi_operator, BasicJacobiOperator* A_jacobi_operator_basic)
{
	int sys_size = A_jacobi_operator->start_index.size() - 1;
	A_jacobi_operator_basic->element.resize(sys_size);
	int* start_index = A_jacobi_operator->start_index.data();
	int* vertex_index = A_jacobi_operator->vertex_index.data();
	double* coefficient = A_jacobi_operator->coefficient.data();

	int end_col_index;
	for (int i = 0; i < sys_size; ++i) {
		end_col_index = start_index[i + 1];
		A_jacobi_operator_basic->element[i].reserve(end_col_index- start_index[i]);
		for (int j = start_index[i]; j < end_col_index; ++j) {
			A_jacobi_operator_basic->element[i].push_back(ColIndexWithCoeff(vertex_index[j], coefficient[j]));
		}
	}
	
}

void IterationMethod::transferBasicOperator2AJacobi(AJacobiOperator* A_jacobi_operator, BasicJacobiOperator* A_jacobi_operator_basic,
	AJacobiOperator* left_multipler_operator, AJacobiOperator* right_multipler_operator)
{
	int sys_size = A_jacobi_operator_basic->element.size();
	A_jacobi_operator->start_index.reserve(sys_size+1);
	int count=0;
	int multipler_count = 0;
	for (int i = 0; i < sys_size; ++i) {
		count += A_jacobi_operator_basic->element[i].size();
		for (int j = 0; j < A_jacobi_operator_basic->element[i].size(); ++j) {
			multipler_count += A_jacobi_operator_basic->element[i][j].index_for_matrix_multiplication.size();
		}
	}
	A_jacobi_operator->vertex_index.reserve(count);
	A_jacobi_operator->coefficient.reserve(count);
	A_jacobi_operator->multiplier_start_per_element.reserve(count + 1);
	A_jacobi_operator->left_multiplier_index.reserve(multipler_count);
	A_jacobi_operator->right_multiplier_index.reserve(multipler_count);

	std::vector<ColIndexWithCoeff>* element;

	std::vector<int>* start_index = &A_jacobi_operator->start_index;
	std::vector<int>* vertex_index = &A_jacobi_operator->vertex_index;
	std::vector<double>* coefficient = &A_jacobi_operator->coefficient;

	std::vector<int>* left_multiplier_index = &A_jacobi_operator->left_multiplier_index;
	std::vector<int>* right_multiplier_index = &A_jacobi_operator->right_multiplier_index;
	std::vector<int>* multiplier_start_per_element = &A_jacobi_operator->multiplier_start_per_element;

	std::vector<int>* index_for_matrix_multiplication;

	start_index->push_back(0);

	int right_col_index;

	multiplier_start_per_element->push_back(0);

	for (int i = 0; i < sys_size; ++i) {		
		element = &A_jacobi_operator_basic->element[i];
		start_index->push_back(start_index->data()[i] + element->size());
		for (int j = 0; j < element->size(); ++j) {
			right_col_index = element->data()[j].vertex_index;
			vertex_index->push_back(right_col_index);
			coefficient->push_back(element->data()[j].coeff);
			index_for_matrix_multiplication = &element->data()[j].index_for_matrix_multiplication;

			//[i,index_for_matrix_multiplication], [index_for_matrix_multiplication,element->vertex_index]
			for (int k = 0; k < index_for_matrix_multiplication->size(); ++k) {
				left_multiplier_index->push_back(obtainIndexInAJacobiOperator(i, index_for_matrix_multiplication->data()[k],
					left_multipler_operator));
				right_multiplier_index->push_back(obtainIndexInAJacobiOperator(index_for_matrix_multiplication->data()[k], right_col_index, 
					right_multipler_operator));
			}
			multiplier_start_per_element->push_back(left_multiplier_index->size());
		}
	}
}

int IterationMethod::obtainIndexInAJacobiOperator(int row_index, int col_index, AJacobiOperator* A_jacobi_operator)
{
	int start_index = A_jacobi_operator->start_index[row_index];
	int end_index = A_jacobi_operator->start_index[row_index + 1];
	for (int i = start_index; i < end_index; ++i) {
		if (A_jacobi_operator->vertex_index[i] == col_index) {
			return i;
		}
	}
	return -1;
}



//A_jacobi_operator_need_to_multi * A_jacobi_operator_basic
void IterationMethod::createHighOrderAJacobiMethod(A_JacobiOperator* A_jacobi_operator_basic, A_JacobiOperator* A_jacobi_operator_need_to_multi, A_JacobiOperator* A_jacobi_operator)
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
//multiple_index means the index of left right multipler to obtain an element
double IterationMethod::obtainElementValue(std::vector<ColIndexWithCoeff>* element_of_i,
	std::vector<ColIndexWithCoeff>* element_of_j, std::vector<int>& indicate_if_vertex_exists, std::vector<int>& multiple_index,
	int left_row_index,	int right_column_index)
{
	multiple_index.reserve(element_of_i->size() + element_of_j->size());
	double value = 0;

	for (int i = 0; i < element_of_i->size(); ++i)
	{
		indicate_if_vertex_exists[element_of_i->data()[i].vertex_index] = i;
	}

	for (int j = 0; j < element_of_j->size(); ++j)
	{
		if (indicate_if_vertex_exists[element_of_j->data()[j].vertex_index] > -1)
		{
			value += element_of_j->data()[j].coeff * element_of_i->data()[indicate_if_vertex_exists[element_of_j->data()[j].vertex_index]].coeff;
			multiple_index.push_back(element_of_j->data()[j].vertex_index);
		}
	}

	for (int i = 0; i < element_of_i->size(); ++i)
	{
		indicate_if_vertex_exists[element_of_i->data()[i].vertex_index] = -1;
	}

	return value;
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



void IterationMethod::setBasicInfo(int sys_size, Thread* thread, SparseMatrix<double, RowMajor>* global_mat)
{
	this->sys_size=sys_size;
	this->thread = thread;
	this->global_mat = global_mat;
	dimension_per_thread.resize(thread->thread_num + 1, 3);
	for (int i = 0; i < 3; ++i) {
		dimension_per_thread[i] = i;
	}
	b_global_inv.resize(3);
	R_b_global_inv.resize(3);

	x_temp.resize(3);
	u_last_for_chebyshev.resize(3);
	u_previous_for_chebyshev.resize(3);

	for (int i = 0; i < 3; ++i) {
		b_global_inv[i].resize(sys_size);
		R_b_global_inv[i].resize(sys_size);
		x_temp[i].resize(sys_size);
		u_previous_for_chebyshev[i].resize(sys_size);
		u_last_for_chebyshev[i].resize(sys_size);
	}

	vertex_index_begin_thread.resize(thread->thread_num + 1);
	A_jacobi_2_index_per_thread.resize(thread->thread_num + 1);
	A_jacobi_3_index_per_thread.resize(thread->thread_num + 1);
	arrangeIndex(thread->thread_num, sys_size, vertex_index_begin_thread);


}





void IterationMethod::updateConvergenceRate(double conv_rate)
{
	convergence_rate_2= conv_rate * conv_rate;
}


void IterationMethod::setOffDiagonal()
{
	off_diagonal = (*global_mat);
	for (int i = 0; i < sys_size; ++i)
	{
		off_diagonal.coeffRef(i, i) = 0.0;
	}
	off_diagonal *= -1.0;
}


void IterationMethod::initialGlobalDiagonalInv(std::vector<double*>* global_mat_diagonal_ref_address)
{
	this->global_mat_diagonal_ref_address = global_mat_diagonal_ref_address;
	global_diagonal_inv.resize(sys_size);
	for (int j = 0; j < sys_size; ++j) {
		global_diagonal_inv.data()[j] = 1.0 / (*((*global_mat_diagonal_ref_address)[j]));
	}
}

void IterationMethod::initialJacobi()
{	
	R_Jacobi = off_diagonal;	
	for (int j = 0; j < sys_size; ++j) {
		for (int k = R_Jacobi.outerIndexPtr()[j]; k < R_Jacobi.outerIndexPtr()[j + 1]; ++k) {
			R_Jacobi.valuePtr()[k] *= global_diagonal_inv.data()[j];
		}
	}
}


void IterationMethod::updateJacobi()
{
	R_Jacobi = off_diagonal;
	for (int i = 0; i < sys_size; ++i) {
		global_diagonal_inv.data()[i] = 1.0 / (*((*global_mat_diagonal_ref_address)[i]));
		for (int k = R_Jacobi.outerIndexPtr()[i]; k < R_Jacobi.outerIndexPtr()[i + 1]; ++k) {
			R_Jacobi.valuePtr()[k] *= global_diagonal_inv.data()[i];
		}
	}
}


void IterationMethod::updateGlobalDiagonalInv()
{
	
	for (int i = 0; i < sys_size; ++i) {
		global_diagonal_inv.data()[i] = 1.0 / (*((*global_mat_diagonal_ref_address)[i]));
	}
	

}

void IterationMethod::solveByJacobi(VectorXd* u, VectorXd* b, int& itr_num)
{
	double residual_norm_per_thread[3];
	double b_norm_cov = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm())* convergence_rate_2;
	double residual_norm = 2.0 * b_norm_cov;
	itr_num = 0;

	for (int i = 0; i < 3; ++i) {
		b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv);
	}
	while (residual_norm> b_norm_cov && itr_num < max_itr_num)
	{
		thread->assignTask(this, JACOBI_ITR, u, b, residual_norm_per_thread,0.0,u,u);
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;
	}

}

//JACOBI_ITR
void IterationMethod::JacobiIterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i)
	{
		u[i] = b_global_inv[i] + R_Jacobi * u[i];
		residual_norm[i] = (b[i] - (*global_mat) * u[i]).squaredNorm();
	}
}

//A_JACOBI_2_ITR
void IterationMethod::SuperJacobi2IterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i)
	{
		u[i] = R_Jacobi * (R_Jacobi * u[i]) + R_b_global_inv[i];
		residual_norm[i] = (b[i] - (*global_mat) * u[i]).squaredNorm();
	}
}

void IterationMethod::solveByAJacobi_2(VectorXd* u, VectorXd* b, int& itr_num)
{
	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm())* convergence_rate_2;
	double residual_norm = 2.0 * b_norm_conv;
	itr_num = 0;
	//std::vector<VectorXd> R_temp(3);
	//std::vector<VectorXd> temp(3);
	//double residual_norm_per_thread[3];
	//for (int i = 0; i < 3; ++i)	{
	//	b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv);
	//	R_b_global_inv[i] = R_Jacobi * b_global_inv[i] + b_global_inv[i];
	//}
	//while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
	//	thread->assignTask(this, A_JACOBI_2_ITR, u, b, residual_norm_per_thread, 0.0,u,u);
	//	residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
	//	itr_num++;
	//}

	std::vector<double> residual_norm_per_thread(thread->thread_num, 0.0);
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < sys_size; ++j) {
			b_global_inv[i][j] = b[i][j]* diagonal_inv[j];
		}	
	}
	thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator.vertex_index.data(), A_jacobi_operator.coefficient.data(),
		A_jacobi_operator.start_index.data(), b_global_inv.data(), b_global_inv.data(), R_b_global_inv.data(), residual_norm_per_thread.data(),
		vertex_index_begin_thread.data(),u,u);

	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		if (itr_num % 2 == 0) {
			thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator_2.vertex_index.data(), A_jacobi_operator_2.coefficient.data(),
				A_jacobi_operator_2.start_index.data(), u, R_b_global_inv.data(), x_temp.data(), residual_norm_per_thread.data(), vertex_index_begin_thread.data(), u, u);
			thread->assignTask(this, COMPUTE_RESIDUAL, global_matrix_operator.vertex_index.data(), global_matrix_operator.coefficient.data(),
				global_matrix_operator.start_index.data(), x_temp.data(), b, b, residual_norm_per_thread.data(), vertex_index_begin_thread.data(), u, u);

		}
		else {
			thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator_2.vertex_index.data(), A_jacobi_operator_2.coefficient.data(),
				A_jacobi_operator_2.start_index.data(), x_temp.data(), R_b_global_inv.data(), u, residual_norm_per_thread.data(), vertex_index_begin_thread.data(), u, u);
			thread->assignTask(this, COMPUTE_RESIDUAL, global_matrix_operator.vertex_index.data(), global_matrix_operator.coefficient.data(),
				global_matrix_operator.start_index.data(), u, b, b, residual_norm_per_thread.data(), vertex_index_begin_thread.data(), u, u);

		}
		residual_norm = residual_norm_per_thread[0];
		for(int i=1;i<thread->thread_num;++i){
			residual_norm += residual_norm_per_thread[i];
		}
		itr_num++;
	}
	if (itr_num % 2 == 1) {
		for (int i = 0; i < 3; ++i) {
			memcpy(u[i].data(), x_temp[i].data(), 8 * sys_size);
		}		
	}

}




//A_JACOBI_3_ITR
void IterationMethod::SuperJacobi3IterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i)
	{
		u[i] = R_Jacobi*(R_Jacobi * (R_Jacobi * u[i])) + R_b_global_inv[i];
		residual_norm[i] = (b[i] - (*global_mat) * u[i]).squaredNorm();
	}
}



void IterationMethod::estimateSuperJacobiEigenValue(VectorXd* u, int A_jacobi_step_size)
{
	switch (A_jacobi_step_size)
	{
	case 2:
		estimateAJacobi2EigenValue(u);

		break;
	case 3:
		estimateAJacobi3EigenValue(u);
		break;
	}
	
	//std::cout << super_jacobi_spectral_radius_square[0] << std::endl;
}


void IterationMethod::estimateGaussSeidelEigenValue(std::vector<VectorXd>& u, SparseMatrix<double, RowMajor>& system_matrix)
{
	double vec_norm2 = 0;
	double u_norm2 = 0;
	for (int i = 0; i < 3; ++i) {
		vec_norm2 += system_matrix.triangularView<Lower>().solve(system_matrix.triangularView<StrictlyUpper>() * u[i]).squaredNorm();
		u_norm2 += u[i].squaredNorm();
	}
	gauss_seidel_spectral_radius_square = vec_norm2 / u_norm2;
}

void IterationMethod::estimateJacobiEigenValue(std::vector<VectorXd>& u)
{
	double vec_norm2 = 0;
	double u_norm2 = 0;
	for (int j = 0; j < 3; ++j) {
		vec_norm2 += (R_Jacobi * u[j]).squaredNorm();
		u_norm2 += u[j].squaredNorm();
		//std::cout << (R_Jacobi[cloth_No] * u[j]).squaredNorm() / u[j].squaredNorm();
	}
	jacobi_spectral_radius_square = vec_norm2 / u_norm2;
}


void IterationMethod::estimateAJacobi2EigenValue(VectorXd* u)
{
	//VectorXd Vector_change;
	//double vec_norm2=0;
	//double u_norm2=0;
	//for (int j = 0; j < 3; ++j) {
	//	Vector_change = R_Jacobi * (R_Jacobi * u[j]);
	//	vec_norm2 += Vector_change.squaredNorm();
	//	u_norm2 += u[j].squaredNorm();
	//}
	//a_jacobi_2_spectral_radius_square = vec_norm2 / u_norm2;

	std::vector<double> u_norm(thread->thread_num,0.0);
	std::vector<double> Ru_norm(thread->thread_num,0.0);

	int vertex_index[1]; int vertex_index_start[1];
	VectorXd b; VectorXd result; b.resize(1); result.resize(1);
	thread->assignTask(this, ESTIMATE_A_JACOBI_2_EIGEN_VALUE, vertex_index, u_norm.data(), vertex_index_start,
		u, &b, &result, Ru_norm.data(), vertex_index_begin_thread.data(), u, u);

	double total_u_norm = u_norm[0];
	double total_Ru_norm = Ru_norm[0];
	for (int i = 1; i < thread->thread_num; ++i) {
		total_u_norm += u_norm[i];
		total_Ru_norm += Ru_norm[i];
	}
	a_jacobi_2_spectral_radius_square = total_Ru_norm / total_u_norm;

}

void IterationMethod::estimateAJacobi3EigenValue(VectorXd* u)
{
	//VectorXd Vector_change;
	//double vec_norm2 = 0;
	//double u_norm2 = 0;
	//for (int j = 0; j < 3; ++j) {
	//	Vector_change = R_Jacobi * u[j];
	//	for (int i = 1; i < 3; ++i) {
	//		Vector_change = R_Jacobi * Vector_change;
	//	}
	//	vec_norm2 += Vector_change.squaredNorm();
	//	u_norm2 += u[j].squaredNorm();
	//}
	//a_jacobi_3_spectral_radius_square = vec_norm2 / u_norm2;

	std::vector<double> u_norm(thread->thread_num, 0.0);
	std::vector<double> Ru_norm(thread->thread_num, 0.0);
	int vertex_index[1]; int vertex_index_start[1];
	VectorXd b; VectorXd result; b.resize(1); result.resize(1);
	thread->assignTask(this, ESTIMATE_A_JACOBI_3_EIGEN_VALUE, vertex_index, u_norm.data(), vertex_index_start,
		u, &b, &result, Ru_norm.data(), vertex_index_begin_thread.data(), u, u);
	double total_u_norm = u_norm[0];
	double total_Ru_norm = Ru_norm[0];
	for (int i = 1; i < thread->thread_num; ++i) {
		total_u_norm += u_norm[i];
		total_Ru_norm += Ru_norm[i];
	}
	a_jacobi_3_spectral_radius_square = total_Ru_norm / total_u_norm;
}


//ESTIMATE_A_JACOBI_2_EIGEN_VALUE
void IterationMethod::estimateAJacobi2EigenValue(VectorXd* u, double* u_norm, double* Ru_norm, int vertex_index_start, int vertex_index_end)
{
	(*u_norm) = 0;
	(*Ru_norm) = 0;
	double* coeff = A_jacobi_operator_2.coefficient.data();
	int* index = A_jacobi_operator_2.vertex_index.data();
	int end;
	int start;
	int value;
	double* u_;
	for (int j = 0; j < 3; ++j) {
		u_ = u[j].data();
		for (int i = vertex_index_start; i < vertex_index_end; ++i) {
			end = A_jacobi_operator_2.start_index[i + 1];
			start = A_jacobi_operator_2.start_index[i];
			value = 0.0;
			for (int k = start; k < end; ++k) {
				value += coeff[k] * u_[index[k]];
			}
			(*Ru_norm) += value * value;
			(*u_norm) += u_[i];
		}
	}
}

//ESTIMATE_A_JACOBI_3_EIGEN_VALUE
void IterationMethod::estimateAJacobi3EigenValue(VectorXd* u, double* u_norm, double* Ru_norm, int vertex_index_start, int vertex_index_end)
{
	(*u_norm) = 0;
	(*Ru_norm) = 0;
	double* coeff = A_jacobi_operator_3.coefficient.data();
	int* index = A_jacobi_operator_3.vertex_index.data();
	int end;
	int start;
	int value;
	double* u_;
	for (int j = 0; j < 3; ++j) {
		u_ = u[j].data();
		for (int i = vertex_index_start; i < vertex_index_end; ++i) {
			end = A_jacobi_operator_3.start_index[i + 1];
			start = A_jacobi_operator_3.start_index[i];
			value = 0.0;
			for (int k = start; k < end; ++k) {
				value += coeff[k] * u_[index[k]];
			}
			(*Ru_norm) += value * value;
			(*u_norm) += u_[i];
		}
	}
}




void IterationMethod::solveByGaussSeidel(VectorXd* u, VectorXd* b, int& itr_num) 
{
	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm()) * convergence_rate_2;
	double residual_norm = 2.0 * b_norm_conv;
	double residual_norm_per_thread[3];
	itr_num = 0;
	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		thread->assignTask(this, GAUSS_SEIDEL_ITR, u, b, residual_norm_per_thread, 0.0, u, u);
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;
	}
}


//GAUSS_SEIDEL_ITR
void IterationMethod::GaussSeidelIterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, double omega_chebyshev,
	VectorXd* u_last, VectorXd* u_previous)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u[i] = (*global_mat).triangularView<Lower>().solve(
			b[i] - (*global_mat).triangularView<StrictlyUpper>() * u[i]);
		residual_norm[i] = (b[i] - (*global_mat) * u[i]).squaredNorm();
	}
}



void IterationMethod::solveByChebyshevGaussSeidel(VectorXd* u, VectorXd* b, int& itr_num, double weight)
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
	thread->assignTask(this, GAUSS_SEIDEL_ITR, u, b, residual_norm_per_thread, 0.0, u, u);
	residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
	itr_num = 1;
	std::vector<VectorXd> u_previous(3);
	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		omega_chebyshev = 4.0 / (4.0 - gauss_seidel_spectral_radius_square * omega_chebyshev);
		thread->assignTask(this, CHEBYSHEV_GAUSS_SEIDEL_ITR, u, b, residual_norm_per_thread, omega_chebyshev, u_last.data(), u_previous.data());
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;
	}

}

//CHEBYSHEV_GAUSS_SEIDEL_ITR
void IterationMethod::ChebyshevSemiIterativeGaussSeidelIterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, double omega_chebyshev,
	VectorXd* u_last, VectorXd* u_previous)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u_previous[i] = u[i];
		u[i] = (*global_mat).triangularView<Lower>().solve(b[i] - (*global_mat).triangularView<StrictlyUpper>() * u[i]);
		u[i] = omega_chebyshev * (weight_for_chebyshev_gauss_seidel * (u[i] - u_previous[i]) + u_previous[i] - u_last[i]) + u_last[i];
		u_last[i] = u_previous[i];
		residual_norm[i] = (b[i] - (*global_mat) * u[i]).squaredNorm();
	}
}






//CHEBYSHEV_A_JACOBI_3_ITR
void IterationMethod::ChebyshevSemiIterativeAJacobi3IterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, double omega_chebyshev,
	VectorXd* u_last, VectorXd* u_previous)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u_previous[i] = u[i];
		u[i] = R_b_global_inv[i] + R_Jacobi *(R_Jacobi * (R_Jacobi * u[i]));
		u[i] = omega_chebyshev * (u[i] - u_last[i]) + u_last[i];
		u_last[i] = u_previous[i];
		residual_norm[i] = (b[i] - (*global_mat) * u[i]).squaredNorm();
	}
}

//CHEBYSHEV_A_JACOBI_2_ITR
void IterationMethod::ChebyshevSemiIterativeAJacobi2IterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, double omega_chebyshev,
	VectorXd* u_last, VectorXd* u_previous)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u_previous[i] = u[i];
		u[i] = R_b_global_inv[i] + R_Jacobi * (R_Jacobi * u[i]);
		u[i] = omega_chebyshev * (u[i] - u_last[i]) + u_last[i];
		u_last[i] = u_previous[i];
		residual_norm[i] = (b[i] - (*global_mat) * u[i]).squaredNorm();
	}
}

//CHEBYSHEV_A_JACOBI_ITERATION
void IterationMethod::ChebyshevAJacobiIterationPerThread(int* vertex_index, double* coefficient, int* vertex_index_start, VectorXd* x, VectorXd* b, VectorXd* result,
	double* omega_chebyshev, int vertex_index_begin, int vertex_index_end, VectorXd* u_last, VectorXd* u_previous)
{
	double* result_dimension;
	double* x_dimension;
	double* u_previous_;
	double* u_last_;
	for (int k = 0; k < 3; ++k) {
		result_dimension = result[k].data();
		x_dimension = x[k].data();
		memcpy(result_dimension + vertex_index_begin, b[k].data() + vertex_index_begin, 8 * (vertex_index_end - vertex_index_begin));
		int vertex_end;
		u_previous_ = u_previous[k].data();
		u_last_ = u_last[k].data();
		memcpy(u_previous_ + vertex_index_begin, x_dimension + vertex_index_begin, 8 * (vertex_index_end - vertex_index_begin));
		for (int i = vertex_index_begin; i < vertex_index_end; ++i)
		{
			vertex_end = vertex_index_start[i + 1];
			for (int j = vertex_index_start[i]; j < vertex_end; ++j) {
				result_dimension[i] += coefficient[j] * x_dimension[vertex_index[j]];
			}
			result_dimension[i] = (*omega_chebyshev) * (result_dimension[i] - u_last_[i]) + u_last_[i];
		}
		memcpy(u_last_ + vertex_index_begin, u_previous_ + vertex_index_begin, 8 * (vertex_index_end - vertex_index_begin));
	}
}


//for input variables, RX+d result cannot be the same with x or d
//to compute residual, also send in operator for global_matrix
//R_MULTIPLY_X_PLUS_B
void IterationMethod::RMultiXPlusb(int* vertex_index, double* coefficient, int* vertex_index_start, VectorXd* x, VectorXd* b, VectorXd* result,
	int vertex_index_begin, int vertex_index_end)
{
	double* result_dimension;
	double* x_dimension;
	for (int k = 0; k < 3; ++k) {
		result_dimension = result[k].data();
		x_dimension = x[k].data();
		memcpy(result_dimension + vertex_index_begin, b[k].data() + vertex_index_begin, 8 * (vertex_index_end - vertex_index_begin));
		int vertex_end;
		for (int i = vertex_index_begin; i < vertex_index_end; ++i)
		{
			vertex_end = vertex_index_start[i + 1];
			for (int j = vertex_index_start[i]; j < vertex_end; ++j) {
				result_dimension[i] += coefficient[j] * x_dimension[vertex_index[j]];
			}
		}
	}
}

void IterationMethod::solveByChebyshevSemiIterativeAJacobi2(VectorXd* u, VectorXd* b, int& itr_num) {

	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm()) * convergence_rate_2;
	double residual_norm;
	double omega_chebyshev = 2.0;
	itr_num = 1;


	//std::vector<VectorXd> u_last(3);
	//std::vector<VectorXd> u_previous(3);
	//for (int i = 0; i < 3; ++i) {
	//	u_last[i] = u[i];
	//	b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv);
	//	R_b_global_inv[i] = R_Jacobi * b_global_inv[i] + b_global_inv[i];
	//}		
	//double residual_norm_per_thread[3];	
	//thread->assignTask(this, A_JACOBI_2_ITR, u, b, residual_norm_per_thread, 0.0, u, u);
	//residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];	
	//while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
	//	omega_chebyshev = 4.0 / (4.0 - a_jacobi_2_spectral_radius_square * omega_chebyshev);
	//	thread->assignTask(this, CHEBYSHEV_A_JACOBI_2_ITR, u, b, residual_norm_per_thread,  omega_chebyshev, u_last.data(), u_previous.data());
	//	residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
	//	itr_num++;
	//}

	std::vector<double> residual_norm_per_thread(thread->thread_num, 0.0);
	for (int i = 0; i < 3; ++i) {
		u_last_for_chebyshev[i] = u[i];
		for (int j = 0; j < sys_size; ++j) {
			b_global_inv[i][j] = b[i][j] * diagonal_inv[j];
		}
	}
	thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator.vertex_index.data(), A_jacobi_operator.coefficient.data(),
		A_jacobi_operator.start_index.data(), b_global_inv.data(), b_global_inv.data(), R_b_global_inv.data(), residual_norm_per_thread.data(),
		vertex_index_begin_thread.data(), u, u);


	thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator_2.vertex_index.data(), A_jacobi_operator_2.coefficient.data(),
		A_jacobi_operator_2.start_index.data(), u, R_b_global_inv.data(), x_temp.data(), residual_norm_per_thread.data(),
		vertex_index_begin_thread.data(), u, u);
	thread->assignTask(this, COMPUTE_RESIDUAL, global_matrix_operator.vertex_index.data(), global_matrix_operator.coefficient.data(),
		global_matrix_operator.start_index.data(), x_temp.data(), b, b, residual_norm_per_thread.data(),
		vertex_index_begin_thread.data(), u, u);
	residual_norm = residual_norm_per_thread[0];
	for (int i = 1; i < thread->thread_num; ++i) {
		residual_norm += residual_norm_per_thread[i];
	}
	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		omega_chebyshev = 4.0 / (4.0 - a_jacobi_2_spectral_radius_square * omega_chebyshev);
		if (itr_num % 2 == 0) {
			thread->assignTask(this, CHEBYSHEV_A_JACOBI_ITERATION, A_jacobi_operator_2.vertex_index.data(), A_jacobi_operator_2.coefficient.data(),
				A_jacobi_operator_2.start_index.data(), u, R_b_global_inv.data(), x_temp.data(), &omega_chebyshev,
				vertex_index_begin_thread.data(), u_last_for_chebyshev.data(), u_previous_for_chebyshev.data());
			thread->assignTask(this, COMPUTE_RESIDUAL, global_matrix_operator.vertex_index.data(), global_matrix_operator.coefficient.data(),
				global_matrix_operator.start_index.data(), x_temp.data(), b, b, residual_norm_per_thread.data(),
				vertex_index_begin_thread.data(), u, u);
		}
		else {
			thread->assignTask(this, CHEBYSHEV_A_JACOBI_ITERATION, A_jacobi_operator_2.vertex_index.data(), A_jacobi_operator_2.coefficient.data(),
				A_jacobi_operator_2.start_index.data(), x_temp.data(), R_b_global_inv.data(), u, &omega_chebyshev,
				vertex_index_begin_thread.data(), u_last_for_chebyshev.data(), u_previous_for_chebyshev.data());
			thread->assignTask(this, COMPUTE_RESIDUAL, global_matrix_operator.vertex_index.data(), global_matrix_operator.coefficient.data(),
				global_matrix_operator.start_index.data(), u, b, b, residual_norm_per_thread.data(),
				vertex_index_begin_thread.data(), u, u);
		}
		residual_norm = residual_norm_per_thread[0];
		for (int i = 1; i < thread->thread_num; ++i) {
			residual_norm += residual_norm_per_thread[i];
		}
		itr_num++;
	}
	if (itr_num % 2 == 1) {
		for (int i = 0; i < 3; ++i) {
			memcpy(u[i].data(), x_temp[i].data(), 8 * sys_size);
		}
	}
}




void IterationMethod::solveByAJacobi_3(VectorXd* u, VectorXd* b, int& itr_num)
{
	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm()) * convergence_rate_2;
	double residual_norm = 2.0 * b_norm_conv;
	itr_num = 0;
	//std::vector<VectorXd> R_temp(3);
	//std::vector<VectorXd> temp(3);
	//double residual_norm_per_thread[3];
	//for (int i = 0; i < 3; ++i) {
	//	b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv);
	//	R_b_global_inv[i] = R_Jacobi * b_global_inv[i];
	//	R_b_global_inv[i] = R_Jacobi * R_b_global_inv[i] + R_b_global_inv[i] + b_global_inv[i];
	//}
	//while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
	//	thread->assignTask(this, A_JACOBI_3_ITR, u, b, residual_norm_per_thread, 0.0,u,u);
	//	residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
	//	itr_num++;
	//}

	std::vector<double> residual_norm_per_thread(thread->thread_num, 0.0);
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < sys_size; ++j) {
			b_global_inv[i][j] = b[i][j] * diagonal_inv[j];
		}
	}
	thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator.vertex_index.data(), A_jacobi_operator.coefficient.data(),
		A_jacobi_operator.start_index.data(), b_global_inv.data(), b_global_inv.data(), x_temp.data(), residual_norm_per_thread.data(),
		vertex_index_begin_thread.data(), u, u);
	thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator_2.vertex_index.data(), A_jacobi_operator_2.coefficient.data(),
		A_jacobi_operator_2.start_index.data(), b_global_inv.data(), x_temp.data(), R_b_global_inv.data(), residual_norm_per_thread.data(),
		vertex_index_begin_thread.data(), u, u);

	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		if (itr_num % 2 == 0) {
			thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator_3.vertex_index.data(), A_jacobi_operator_3.coefficient.data(),
				A_jacobi_operator_3.start_index.data(), u, R_b_global_inv.data(), x_temp.data(), residual_norm_per_thread.data(), 
				vertex_index_begin_thread.data(), u, u);
			thread->assignTask(this, COMPUTE_RESIDUAL, global_matrix_operator.vertex_index.data(), global_matrix_operator.coefficient.data(),
				global_matrix_operator.start_index.data(), x_temp.data(), b, b, residual_norm_per_thread.data(), 
				vertex_index_begin_thread.data(), u, u);
		}
		else {
			thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator_3.vertex_index.data(), A_jacobi_operator_3.coefficient.data(),
				A_jacobi_operator_3.start_index.data(), x_temp.data(), R_b_global_inv.data(), u, residual_norm_per_thread.data(), 
				vertex_index_begin_thread.data(), u, u);
			thread->assignTask(this, COMPUTE_RESIDUAL, global_matrix_operator.vertex_index.data(), global_matrix_operator.coefficient.data(),
				global_matrix_operator.start_index.data(), u, b, b, residual_norm_per_thread.data(), 
				vertex_index_begin_thread.data(), u, u);
		}
		residual_norm = residual_norm_per_thread[0];
		for (int i = 1; i < thread->thread_num; ++i) {
			residual_norm += residual_norm_per_thread[i];
		}
		itr_num++;
	}
	if (itr_num % 2 == 1) {
		for (int i = 0; i < 3; ++i) {
			memcpy(u[i].data(), x_temp[i].data(), 8 * sys_size);
		}
	}
}



void IterationMethod::solveByChebyshevSemiIterativeAJacobi3(VectorXd* u, VectorXd* b, int& itr_num) {
	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm()) * convergence_rate_2;
	double omega_chebyshev = 2.0;
	double residual_norm;
	itr_num = 1;

	//std::vector<VectorXd> u_last(3);
	//std::vector<VectorXd> u_previous(3);
	//for (int i = 0; i < 3; ++i) {
	//	u_last[i] = u[i];
	//	b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv);
	//	R_b_global_inv[i] = R_Jacobi * b_global_inv[i];
	//	R_b_global_inv[i] = R_Jacobi * b_global_inv[i] + R_b_global_inv[i] + b_global_inv[i];
	//}	
	//double residual_norm_per_thread[3];	
	//thread->assignTask(this, A_JACOBI_3_ITR, u, b, residual_norm_per_thread,  0.0, u, u);
	//residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];	
	//while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
	//	omega_chebyshev = 4.0 / (4.0 - a_jacobi_2_spectral_radius_square * omega_chebyshev);
	//	thread->assignTask(this, CHEBYSHEV_A_JACOBI_3_ITR, u, b, residual_norm_per_thread,  omega_chebyshev, u_last.data(), u_previous.data());
	//	residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
	//	itr_num++;
	//}

	std::vector<double> residual_norm_per_thread(thread->thread_num, 0.0);
	for (int i = 0; i < 3; ++i) {
		u_last_for_chebyshev[i] = u[i];
		for (int j = 0; j < sys_size; ++j) {
			b_global_inv[i][j] = b[i][j] * diagonal_inv[j];
		}
	}
	thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator.vertex_index.data(), A_jacobi_operator.coefficient.data(),
		A_jacobi_operator.start_index.data(), b_global_inv.data(), b_global_inv.data(), x_temp.data(), residual_norm_per_thread.data(),
		vertex_index_begin_thread.data(), u, u);
	thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator_2.vertex_index.data(), A_jacobi_operator_2.coefficient.data(),
		A_jacobi_operator_2.start_index.data(), b_global_inv.data(), x_temp.data(), R_b_global_inv.data(), residual_norm_per_thread.data(),
		vertex_index_begin_thread.data(), u, u);


	thread->assignTask(this, R_MULTIPLY_X_PLUS_B, A_jacobi_operator_3.vertex_index.data(), A_jacobi_operator_3.coefficient.data(),
		A_jacobi_operator_3.start_index.data(), u, R_b_global_inv.data(), x_temp.data(), residual_norm_per_thread.data(),
		vertex_index_begin_thread.data(), u, u);
	thread->assignTask(this, COMPUTE_RESIDUAL, global_matrix_operator.vertex_index.data(), global_matrix_operator.coefficient.data(),
		global_matrix_operator.start_index.data(), x_temp.data(), b, b, residual_norm_per_thread.data(),
		vertex_index_begin_thread.data(), u, u);
	residual_norm = residual_norm_per_thread[0];
	for (int i = 1; i < thread->thread_num; ++i) {
		residual_norm += residual_norm_per_thread[i];
	}
	while (residual_norm > b_norm_conv && itr_num < max_itr_num) {
		omega_chebyshev = 4.0 / (4.0 - a_jacobi_2_spectral_radius_square * omega_chebyshev);
		if (itr_num % 2 == 0) {
			thread->assignTask(this, CHEBYSHEV_A_JACOBI_ITERATION, A_jacobi_operator_3.vertex_index.data(), A_jacobi_operator_3.coefficient.data(),
				A_jacobi_operator_3.start_index.data(), u, R_b_global_inv.data(), x_temp.data(), &omega_chebyshev,
				vertex_index_begin_thread.data(), u_last_for_chebyshev.data(), u_previous_for_chebyshev.data());
			thread->assignTask(this, COMPUTE_RESIDUAL, global_matrix_operator.vertex_index.data(), global_matrix_operator.coefficient.data(),
				global_matrix_operator.start_index.data(), x_temp.data(), b, b, residual_norm_per_thread.data(),
				vertex_index_begin_thread.data(), u, u);
		}
		else {
			thread->assignTask(this, CHEBYSHEV_A_JACOBI_ITERATION, A_jacobi_operator_3.vertex_index.data(), A_jacobi_operator_3.coefficient.data(),
				A_jacobi_operator_3.start_index.data(), x_temp.data(), R_b_global_inv.data(), u, &omega_chebyshev,
				vertex_index_begin_thread.data(), u_last_for_chebyshev.data(), u_previous_for_chebyshev.data());
			thread->assignTask(this, COMPUTE_RESIDUAL, global_matrix_operator.vertex_index.data(), global_matrix_operator.coefficient.data(),
				global_matrix_operator.start_index.data(), u, b, b, residual_norm_per_thread.data(),
				vertex_index_begin_thread.data(), u, u);
		}

		residual_norm = residual_norm_per_thread[0];
		for (int i = 1; i < thread->thread_num; ++i) {
			residual_norm += residual_norm_per_thread[i];
		}
		itr_num++;
	}
	if (itr_num % 2 == 1) {
		for (int i = 0; i < 3; ++i) {
			memcpy(u[i].data(), x_temp[i].data(), 8 * sys_size);
		}
	}
}


void IterationMethod::solveByChebyshevSemiIterativeJacobi(VectorXd* u, VectorXd* b, int& itr_num)
{
	std::vector<VectorXd> u_last(3);
	std::vector<VectorXd> u_previous(3);
	for (int i = 0; i < 3; ++i) {
		u_last[i] = u[i];
		b_global_inv[i] = b[i].cwiseProduct(global_diagonal_inv);
	}
	double residual_norm_per_thread[3];
	double residual_norm;
	thread->assignTask(this, JACOBI_ITR, u, b, residual_norm_per_thread,  0.0, u, u);
	residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];

	double b_norm_conv = (b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm()) * convergence_rate_2;
	double omega_chebyshev = 2.0;	
	itr_num = 1;
	while (residual_norm  > b_norm_conv && itr_num < max_itr_num) {
		omega_chebyshev = 4.0 / (4.0 - jacobi_spectral_radius_square * omega_chebyshev);
		thread->assignTask(this, CHEBYSHEV_JACOBI_ITR, u, b, residual_norm_per_thread,  
			omega_chebyshev, u_last.data(), u_previous.data());
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		itr_num++;	
	}


}


//CHEBYSHEV_JACOBI_ITR
void IterationMethod::ChebyshevSemiIterativeJacobiIterationPerThread(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm, double omega_chebyshev,
	VectorXd* u_last, VectorXd* u_previous)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u_previous[i] = u[i];
		u[i] = b_global_inv[i] + R_Jacobi * u[i];
		u[i] = omega_chebyshev * (u[i] - u_last[i]) + u_last[i];
		u_last[i] = u_previous[i];
		residual_norm[i] = (b[i] - (*global_mat) * u[i]).squaredNorm();
	}
}


void IterationMethod::solveByPCG(VectorXd* u, VectorXd* b,  int& itr_num)
{
	std::vector<VectorXd> residual(3);
	double residual_norm;
	double residual_norm_per_thread[3];
	double b_norm_cov=(b[0].squaredNorm() + b[1].squaredNorm() + b[2].squaredNorm())*convergence_rate_2;
	std::vector<VectorXd> z(3);
	std::vector<VectorXd> p(3);
	for (int i = 0; i < 3; ++i) {
		residual[i] = b[i] - (*global_mat) * u[i];
		z[i] = global_diagonal_inv.cwiseProduct(residual[i]);
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
		alpha = rz_k / (((*global_mat) * p[0]).dot(p[0]) + 
			((*global_mat) * p[1]).dot(p[1]) + ((*global_mat) * p[2]).dot(p[2]));
		thread->assignTask(this, PCG_ITR1, u, b, residual_norm_per_thread,
			alpha, residual.data(), p.data());
		residual_norm = residual_norm_per_thread[0] + residual_norm_per_thread[1] + residual_norm_per_thread[2];
		if (residual_norm < b_norm_cov || itr_num >= max_itr_num) {
			break;
		}
		thread->assignTask(this, PCG_ITR2, z.data(), residual.data(), rz_k_1_per_thread,
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
void IterationMethod::PCGIterationPerThread1(int thread_id, VectorXd* u, VectorXd* b, double* residual_norm,  double alpha,
	VectorXd* residual, VectorXd* p)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		u[i] += alpha* p[i];
		residual[i] = b[i] - (*global_mat) * u[i];
		residual_norm[i] = residual[i].squaredNorm();
	}
}

//PCG_ITR2
void IterationMethod::PCGIterationPerThread2(int thread_id, VectorXd* z, VectorXd* residual, double* rz_k_1,  double alpha,
	VectorXd* q, VectorXd* p)
{
	for (int i = dimension_per_thread[thread_id]; i < dimension_per_thread[thread_id + 1]; ++i) {
		z[i] = global_diagonal_inv.cwiseProduct(residual[i]);
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



