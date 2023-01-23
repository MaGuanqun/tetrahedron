#include"second_order.h"
#include"FEM_relate.h"
#include"../collision/primitive_distance_gradient_hessian.h"

void SecondOrderConstraint::computeEdgeLengthForce(double* vertex_0, double* vertex_1, double stiffness,
	double* potential_0, double* potential_1, double rest_length)
{
	double C = sqrt(EDGE_LENGTH(vertex_0, vertex_1));
	double grad[3];
	grad[0] = (vertex_0[0] - vertex_1[0]) / C;
	grad[1] = (vertex_0[1] - vertex_1[1]) / C;
	grad[2] = (vertex_0[2] - vertex_1[2]) / C;
	C -= rest_length;
	double coe = C * stiffness;
	MULTI(potential_0, grad, coe);
	coe *= -1.0;
	MULTI(potential_1, grad, coe);
}



void SecondOrderConstraint::computeARAPForce(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	double stiffness, Matrix<double, 3, 4>& A, double volume, Matrix<double,3,4>& force)
{
	Vector3d eigen_value;
	Matrix3d q_e;
	double determinant;
	Vector3d position;
	memcpy(q_e.data(), vertex_position_1, 24);
	memcpy(q_e.data() + 3, vertex_position_2, 24);
	memcpy(q_e.data() + 6, vertex_position_3, 24);
	//first use eigen value to store the position of first vertex
	memcpy(eigen_value.data(), vertex_position_0, 24);
	for (unsigned int i = 0; i < 3; ++i) {
		q_e.col(i) -= eigen_value;
	}
	Matrix3d deformation_gradient;
	Matrix3d P_inv;
	memcpy(P_inv.data(), A.data() + 3, 72);
	deformation_gradient = q_e * P_inv.transpose();
	JacobiSVD<Matrix3d> svd;
	svd.compute(deformation_gradient, ComputeFullU | ComputeFullV);

	eigen_value = svd.singularValues();
	determinant = eigen_value[0] * eigen_value[1] * eigen_value[2];

	Matrix3d V = svd.matrixV();
	if (determinant < 0) {
		V.col(2) *= -1.0;
	}

	P_inv = svd.matrixU() * V.transpose();

	force = (stiffness*volume)* (deformation_gradient - P_inv) * A;

}

void SecondOrderConstraint::solveSingleVertexNewton(std::array<double, 3>* vertex_position, double stiffness, double dt,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn)
{
	unsigned int tet_index;
	Matrix3d Hessian_single;
	Vector3d grad_single;
	unsigned int vertex_no;
	Matrix3d Hessian;
	Vector3d grad;
	Hessian.setZero();
	grad.setZero();
	double C;
	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		tet_index = tet_indices[i];
		vertex_no = findVertexNo(vertex_index, indices[tet_index].data(),4);
		if (getARAPGradHessianNewton(vertex_position[indices[tet_index][0]].data(), vertex_position[indices[tet_index][1]].data(),
			vertex_position[indices[tet_index][2]].data(), vertex_position[indices[tet_index][3]].data(),
			A[tet_index], Hessian_single, grad_single, C, vertex_no)) {
			Hessian += (0.5 * stiffness * volume[tet_index]) * Hessian_single;
			grad -= (0.5 * stiffness * volume[tet_index]) * grad_single;
		}
	}

	double mass_dt_2 = mass[vertex_index] / (dt * dt);

	Hessian.data()[0] += mass_dt_2;
	Hessian(1, 1) += mass_dt_2;
	Hessian(2, 2) += mass_dt_2;
	grad.data()[0] -= mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	grad.data()[1] -= mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	grad.data()[2] -= mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <Matrix3d> linear(Hessian);
	Vector3d result = linear.solve(grad);

	SUM_(vertex_position[vertex_index], result);
}



void SecondOrderConstraint::solveNewtonCD_ARAP(std::array<double, 3>* vertex_position, double stiffness, double dt,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn)
{
	unsigned int tet_index;
	Matrix3d Hessian_single;
	Vector3d grad_single;
	unsigned int vertex_no;
	Matrix3d Hessian;
	Vector3d grad;
	Hessian.setZero();
	grad.setZero();
	double C;
	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		tet_index = tet_indices[i];
		vertex_no = findVertexNo(vertex_index, indices[tet_index].data(),4);
		if (getARAPGradHessianNewton(vertex_position[indices[tet_index][0]].data(), vertex_position[indices[tet_index][1]].data(),
			vertex_position[indices[tet_index][2]].data(), vertex_position[indices[tet_index][3]].data(),
			A[tet_index], Hessian_single, grad_single, C, vertex_no)) {
			Hessian += (0.5*stiffness*volume[tet_index])* Hessian_single;
			grad-= (0.5*stiffness * volume[tet_index]) * grad_single;
		}
	}

	double mass_dt_2 = mass[vertex_index] / (dt * dt);

	Hessian.data()[0] += mass_dt_2;
	Hessian(1, 1) += mass_dt_2;
	Hessian(2, 2) += mass_dt_2;
	grad.data()[0] -= mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	grad.data()[1] -= mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	grad.data()[2] -= mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <Matrix3d> linear(Hessian);
	Vector3d result = linear.solve(grad);

	SUM_(vertex_position[vertex_index], result);

}

void SecondOrderConstraint::solveSingleVertexCD_ARAP(std::array<double, 3>* vertex_position, double stiffness, double dt,
	Matrix<double, 3, 4>* A, double* lambda, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn)
{
	unsigned int tet_index;
	Matrix3d Hessian_single;
	Vector3d grad_single;
	unsigned int vertex_no;
	MatrixXd Hessian;
	Hessian.resize(3 + tet_indices.size(), 3 + tet_indices.size());
	Hessian.setZero();

	VectorXd b;
	b.resize(3 + tet_indices.size());
	b.setZero();
	double C;
	double alpha;

	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		tet_index = tet_indices[i];
		vertex_no = findVertexNo(vertex_index, indices[tet_index].data(),4);
		if (getARAPGradHessian(vertex_position[indices[tet_index][0]].data(), vertex_position[indices[tet_index][1]].data(),
			vertex_position[indices[tet_index][2]].data(), vertex_position[indices[tet_index][3]].data(),
			A[tet_index], Hessian_single, grad_single, C, vertex_no)) {
			Hessian_single *= -lambda[tet_index];
			Hessian.block<3, 3>(0, 0) += Hessian_single;
			Hessian.block<3, 1>(0, 3 + i) = -grad_single;
			Hessian.block<1, 3>(3 + i, 0) = -grad_single;
			b.head(3) += lambda[tet_index] * grad_single;
		}
		alpha = 1.0 / (stiffness * volume[tet_index] * dt * dt);
		Hessian(3 + i, 3 + i) = -alpha;
		b.data()[3 + i] = C + alpha * lambda[tet_index];
	}

	Hessian.data()[0] += mass[vertex_index];
	Hessian(1, 1) += mass[vertex_index];
	Hessian(2, 2) += mass[vertex_index];
	b.data()[0] -= mass[vertex_index] * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	b.data()[1] -= mass[vertex_index] * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	b.data()[2] -= mass[vertex_index] * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(b);

	SUM_(vertex_position[vertex_index], result);
	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		lambda[tet_indices[i]] += result[3 + i];
	}
}

void SecondOrderConstraint::setARAPHessian(MatrixXd& Hessian,double stiffness,
	Matrix<double, 3, 4>& A, double volume)
{
	Hessian = A.transpose() * (A * (2.0 * volume * stiffness));
}


void SecondOrderConstraint::setARAPGrad(MatrixXd& grad, std::array<double, 3>* vertex_position, double stiffness,
	Matrix<double, 3, 4>& A, double volume, int* tet_vertex_index)
{
	Matrix3d deformation_gradient;
	Matrix3d rotation;

	FEM::getDeformationGradient(vertex_position[tet_vertex_index[0]].data(), vertex_position[tet_vertex_index[1]].data(),
		vertex_position[tet_vertex_index[2]].data(), vertex_position[tet_vertex_index[3]].data(), A,
		deformation_gradient);

	FEM::extractRotation(deformation_gradient, rotation);
	grad = (2.0 * stiffness * volume) * (deformation_gradient - rotation) * A;
}

void SecondOrderConstraint::solveCD_ARAP_block_fromRecord(MatrixXd& Hessian, VectorXd& grad, int* unfixed_tet_vertex_index,
	unsigned int unfixed_vertex_num, double* hessian_record, double* grad_record)
{
	Hessian.setZero();
	grad.setZero();

	for (int i = 0; i < unfixed_vertex_num; ++i) {
		memcpy(grad.data() + 3 * i, grad_record + 3 * unfixed_tet_vertex_index[i], 24);
	}

	if (unfixed_vertex_num == 4) {
		for (unsigned int i = 0; i < 144; i += 36) {
			for (unsigned int j = 0; j < 12; j += 3) {
				Hessian.data()[i + j] += *hessian_record;
				Hessian.data()[i + j + 13] += *hessian_record;
				Hessian.data()[i + j + 26] += *hessian_record;
				hessian_record++;
			}
		}
	}
	else {
		for (int i = 0; i < unfixed_vertex_num; ++i) {
			for (int j = 0; j < unfixed_vertex_num; ++j) {
				Hessian(3 * i, 3 * j) = hessian_record[(unfixed_tet_vertex_index[j] << 2) + unfixed_tet_vertex_index[i]];
				Hessian(3 * i + 1, 3 * j + 1) = hessian_record[(unfixed_tet_vertex_index[j] << 2) + unfixed_tet_vertex_index[i]];
				Hessian(3 * i + 2, 3 * j + 2) = hessian_record[(unfixed_tet_vertex_index[j] << 2) + unfixed_tet_vertex_index[i]];
			}
		}
	}
}


void SecondOrderConstraint::solveCD_ARAP_block(MatrixXd& Hessian, VectorXd& grad, std::array<double, 3>* vertex_position, double stiffness,
	Matrix<double, 3, 4>& A,
	double volume,
	int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num, double* hessian_record)
{

	
	Matrix3d deformation_gradient;
	Matrix3d rotation;
	Hessian.setZero();
	grad.setZero();

	FEM::getDeformationGradient(vertex_position[tet_vertex_index[0]].data(), vertex_position[tet_vertex_index[1]].data(),
		vertex_position[tet_vertex_index[2]].data(), vertex_position[tet_vertex_index[3]].data(), A, 
		deformation_gradient);


	if (solve_exact_ARAP_hessian) {
		Matrix3d S; Vector3d eigen_value;
		FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (2.0 * stiffness * volume) * (deformation_gradient - rotation) * A;
		for (int i = 0; i < unfixed_vertex_num; ++i) {
			memcpy(grad.data() + 3 * i, grad_C_transpose.data() + 3 * unfixed_tet_vertex_index[i], 24);
		}
		Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
		FEM::getHessianForSeveralVertex(Hessian, S, rotation, Dm, A, unfixed_tet_vertex_index,
			unfixed_vertex_num);
		Hessian *= volume * stiffness;
	}
	else {

		FEM::extractRotation(deformation_gradient, rotation);
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (2.0 * stiffness * volume) * (deformation_gradient - rotation) * A;
		for (int i = 0; i < unfixed_vertex_num; ++i) {
			memcpy(grad.data() + 3 * i, grad_C_transpose.data() + 3 * unfixed_tet_vertex_index[i], 24);
		}

		if (unfixed_vertex_num == 4) {
			for (unsigned int i = 0; i < 144; i += 36) {
				for (unsigned int j = 0; j < 12; j += 3) {
					Hessian.data()[i + j] += *hessian_record;
					Hessian.data()[i + j + 13] += *hessian_record;
					Hessian.data()[i + j + 26] += *hessian_record;
					hessian_record++;
				}
			}
		}
		else {
			for (int i = 0; i < unfixed_vertex_num; ++i) {
				for (int j = 0; j < unfixed_vertex_num; ++j) {
					Hessian(3 * i, 3 * j) = hessian_record[(unfixed_tet_vertex_index[j] << 2) + unfixed_tet_vertex_index[i]];
					Hessian(3 * i + 1, 3 * j + 1) = hessian_record[(unfixed_tet_vertex_index[j] << 2) + unfixed_tet_vertex_index[i]];
					Hessian(3 * i + 2, 3 * j + 2) = hessian_record[(unfixed_tet_vertex_index[j] << 2) + unfixed_tet_vertex_index[i]];
				}
			}
		}
	}

}



void SecondOrderConstraint::solveCD_ARAP_blockTest(MatrixXd& Hessian, VectorXd& grad, std::array<double, 3>* vertex_position, double stiffness,
	Matrix<double, 3, 4>& A,
	double volume,
	int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num)
{


	Matrix3d deformation_gradient;
	Matrix3d rotation;
	Hessian.setZero();
	grad.setZero();

	FEM::getDeformationGradient(vertex_position[tet_vertex_index[0]].data(), vertex_position[tet_vertex_index[1]].data(),
		vertex_position[tet_vertex_index[2]].data(), vertex_position[tet_vertex_index[3]].data(), A,
		deformation_gradient);


	if (solve_exact_ARAP_hessian) {
		Matrix3d S; Vector3d eigen_value;
		FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (2.0 * stiffness * volume) * (deformation_gradient - rotation) * A;
		for (int i = 0; i < unfixed_vertex_num; ++i) {
			memcpy(grad.data() + 3 * i, grad_C_transpose.data() + 3 * unfixed_tet_vertex_index[i], 24);
		}
		Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
		FEM::getHessianForSeveralVertex(Hessian, S, rotation, Dm, A, unfixed_tet_vertex_index,
			unfixed_vertex_num);
		Hessian *= volume * stiffness;
	}
	else {

		FEM::extractRotation(deformation_gradient, rotation);
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (2.0 * stiffness * volume) * (deformation_gradient - rotation) * A;
		for (int i = 0; i < unfixed_vertex_num; ++i) {
			memcpy(grad.data() + 3 * i, grad_C_transpose.data() + 3 * unfixed_tet_vertex_index[i], 24);
		}


		Matrix4d result = A.transpose() * (A * (2.0 * volume * stiffness));

		if (unfixed_vertex_num == 4) {
			double* address = result.data();

				for (unsigned int i = 0; i < 144; i += 36) {
					for (unsigned int j = 0; j < 12; j += 3) {
						Hessian.data()[i + j] += *address;
						Hessian.data()[i + j + 13] += *address;
						Hessian.data()[i + j + 26] += *address;
						address++;
					}
				}
		}
		else {
			for (int i = 0; i < unfixed_vertex_num; ++i) {
				for (int j = 0; j < unfixed_vertex_num; ++j) {
					Hessian(3 * i, 3 * j) += result.data()[(unfixed_tet_vertex_index[j] << 2) + unfixed_tet_vertex_index[i]];
					Hessian(3 * i + 1, 3 * j + 1) += result.data()[(unfixed_tet_vertex_index[j] << 2) + unfixed_tet_vertex_index[i]];
					Hessian(3 * i + 2, 3 * j + 2) += result.data()[(unfixed_tet_vertex_index[j] << 2) + unfixed_tet_vertex_index[i]];
				}
			}
		}
	}

	//std::cout <<"this tet "<< deformation_gradient.determinant() << " " << rotation.determinant() << std::endl;

}




bool SecondOrderConstraint::solveTetCertainVertices(std::array<double, 3>* vertex_position, double stiffness,
	Matrix<double, 3, 4>& A, int* vertex_in_sys, int* tet_vetex_indices,
	MatrixXd& sys_matrix, double volume, VectorXd& grad)
{
	Matrix<double, 12, 3> Hessian_vertex;


	Matrix3d deformation_gradient;
	Matrix3d rotation;

	FEM::getDeformationGradient(vertex_position[tet_vetex_indices[0]].data(), vertex_position[tet_vetex_indices[1]].data(),
		vertex_position[tet_vetex_indices[2]].data(), vertex_position[tet_vetex_indices[3]].data(), A,
		deformation_gradient);
	stiffness *= volume;


	if (solve_exact_ARAP_hessian) {
		Matrix3d S;	Vector3d eigen_value;
		FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (2.0 * stiffness) * (deformation_gradient - rotation) * A;

		for (int i = 0; i < 4; ++i) {
			if (vertex_in_sys[i] != -1) {
				grad.segment(3 * vertex_in_sys[i], 3)
					+= grad_C_transpose.col(i);
			}
		}
		Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
		for (int i = 0; i < 4; ++i) {
			if (vertex_in_sys[i] != -1) {
				FEM::getHessianForOneVertex(Hessian_vertex, S, rotation, Dm, A, i);
				Hessian_vertex *= stiffness;
				for (int j = 0; j < 4; ++j) {
					if (vertex_in_sys[j] != -1) {
						sys_matrix.block<3, 3>(3 * vertex_in_sys[j],
							3 * vertex_in_sys[i]) +=
							Hessian_vertex.block<3, 3>(3 * j, 0);
					}
				}
			}
		}
	}
	else {
		FEM::extractRotation(deformation_gradient, rotation);
		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (2.0 * stiffness) * (deformation_gradient - rotation) * A;

		for (int i = 0; i < 4; ++i) {
			if (vertex_in_sys[i] != -1) {
				grad.segment(3 * vertex_in_sys[i], 3)
					+= grad_C_transpose.col(i);
			}
		}

		double value;
		for (int i = 0; i < 4; ++i) {
			if (vertex_in_sys[i] != -1) {
				for (int j = 0; j < 4; ++j) {
					if (vertex_in_sys[j] != -1) {
						value = 2.0 * stiffness * DOT((A.data() +3*j), (A.data() + 3*i));
						sys_matrix(3 * vertex_in_sys[j],
							3 * vertex_in_sys[i]) += value;
						sys_matrix(3 * vertex_in_sys[j]+1,
							3 * vertex_in_sys[i]+1) += value;
						sys_matrix(3 * vertex_in_sys[j]+2,
							3 * vertex_in_sys[i]+2) += value;
					}
				}
			}
		}
	}
	return true;
}

void SecondOrderConstraint::solveHessianForNeighborTet(unsigned int*& common_vertex_in_order, 
	MatrixXd& sys_matrix, VectorXd& grad, double* hessian_record, double* grad_record)
{
	if (*common_vertex_in_order == 0) {
		common_vertex_in_order++;
		return;
	}
	int common_vertex_num = *common_vertex_in_order;
	common_vertex_in_order++;

	double* grad_address;
	double* grad_record_address;


	int row, col;

	double* hessian_address;
	int move_size = sys_matrix.cols() + 1;
	double value;
	for (int i = 0; i < common_vertex_num; ++i) {
		row = 3 * (*(common_vertex_in_order + i + common_vertex_num));
		grad_address = grad.data() + row;
		grad_record_address = grad_record + 3 * (*(common_vertex_in_order + i));
		(*grad_address) += *grad_record_address;
		*(grad_address + 1) += *(grad_record_address + 1);
		*(grad_address + 2) += *(grad_record_address + 2);
		
		for (int j = i + 1; j < common_vertex_num; ++j) {
			col = 3 * (*(common_vertex_in_order + j + common_vertex_num));
			value = *(hessian_record + (((*(common_vertex_in_order + i)) << 2) + (*(common_vertex_in_order + j))));

			hessian_address = &sys_matrix(row, col);
			*hessian_address += value;
			*(hessian_address + move_size)+= value;
			*(hessian_address +(move_size<<1) )+= value;

			hessian_address = &sys_matrix(col, row);
			*hessian_address += value;
			*(hessian_address + move_size) += value;
			*(hessian_address + (move_size << 1)) += value;
		}

		hessian_address = &sys_matrix(row, row);
		value = *(hessian_record + (((*(common_vertex_in_order + i)) << 2) + (*(common_vertex_in_order + i))));
		*hessian_address += value;
		*(hessian_address + move_size) += value;
		*(hessian_address + (move_size << 1)) += value;

	}


	common_vertex_in_order += (common_vertex_num << 1);
}

void SecondOrderConstraint::solveCertainHessianForNeighborTetTest(std::array<double, 3>* vertex_position, double stiffness,
	Matrix<double, 3, 4>& A, unsigned int*& common_vertex_in_order, int* neighbor_tet_vetex_indices,
	MatrixXd& sys_matrix, double volume, VectorXd& grad, double* hessian_record, double* grad_record)
{
	if (*common_vertex_in_order == 0) {
		common_vertex_in_order++;
		return;
	}
	int common_vertex_num = *common_vertex_in_order;
	common_vertex_in_order++;


	Matrix<double, 12, 3> Hessian_vertex;

	Matrix3d deformation_gradient;
	Matrix3d rotation;


	FEM::getDeformationGradient(vertex_position[neighbor_tet_vetex_indices[0]].data(), vertex_position[neighbor_tet_vetex_indices[1]].data(),
		vertex_position[neighbor_tet_vetex_indices[2]].data(), vertex_position[neighbor_tet_vetex_indices[3]].data(), A,
		deformation_gradient);
	stiffness *= volume;
	//if ((deformation_gradient - rotation).squaredNorm() < 1e-16) {
		//common_vertex_in_order += (common_vertex_num * 2);
		//return false;
	//}

	double* grad_address;
	double* grad_record_address;


	int row, col;

	if (solve_exact_ARAP_hessian) {
		Matrix3d S;	Vector3d eigen_value;
		FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);

		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (2.0 * stiffness) * (deformation_gradient - rotation) * A;

		for (int i = 0; i < common_vertex_num; ++i) {
			grad.segment(3 * (*(common_vertex_in_order + i + common_vertex_num)), 3)
				+= grad_C_transpose.col(*(common_vertex_in_order + i));
		}

		Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
		for (int i = 0; i < common_vertex_num; ++i) {
			FEM::getHessianForOneVertex(Hessian_vertex, S, rotation, Dm, A, *(common_vertex_in_order + i));
			Hessian_vertex *= stiffness;
			for (int j = 0; j < common_vertex_num; ++j) {
				sys_matrix.block<3, 3>(3 * (*(common_vertex_in_order + j + common_vertex_num)),
					3 * (*(common_vertex_in_order + i + common_vertex_num))) +=
					Hessian_vertex.block<3, 3>(3 * (*(common_vertex_in_order + j)), 0);
			}
		}
	}
	else {
		FEM::extractRotation(deformation_gradient, rotation);

		//if (rotation.determinant() < 0.9 || rotation.determinant() > 1.1) {
		//	std::cout << "rotation error " << rotation.determinant() << std::endl;
		//}

		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (2.0 * stiffness) * (deformation_gradient - rotation) * A;

		//for (int i = 0; i < common_vertex_num; ++i) {
		//	grad.segment(3 * (*(common_vertex_in_order + i + common_vertex_num)), 3)
		//		+= grad_C_transpose.col(*(common_vertex_in_order + i));
		//}

		for (int i = 0; i < common_vertex_num; ++i) {
			row = 3 * (*(common_vertex_in_order + i + common_vertex_num));
			grad_address = grad.data() + row;
			grad_record_address = grad_record +3 * (*(common_vertex_in_order + i));
			(*grad_address) += *grad_record_address;
			*(grad_address + 1) += *(grad_record_address + 1);
			*(grad_address + 2) += *(grad_record_address + 2);
		}

		//for (unsigned int i = 0; i < 12; ++i) {
		//	std::cout << grad_record[i] << " ";
		//}
		//std::cout << std::endl;
		//std::cout << grad_C_transpose << std::endl;

		double* hessian_address;
		int move_size = sys_matrix.cols() + 1;
		double value;
		for (int i = 0; i < common_vertex_num; ++i) {
			row = 3 * (*(common_vertex_in_order + i + common_vertex_num));
			for (int j = i + 1; j < common_vertex_num; ++j) {
				col = 3 * (*(common_vertex_in_order + j + common_vertex_num));
				value = *(hessian_record + (((*(common_vertex_in_order + i)) << 2) + (*(common_vertex_in_order + j))));

				hessian_address = &sys_matrix(row, col);
				*hessian_address += value;
				*(hessian_address + move_size) += value;
				*(hessian_address + (move_size << 1)) += value;

				hessian_address = &sys_matrix(col, row);
				*hessian_address += value;
				*(hessian_address + move_size) += value;
				*(hessian_address + (move_size << 1)) += value;
			}

			hessian_address = &sys_matrix(row, row);
			value = *(hessian_record + (((*(common_vertex_in_order + i)) << 2) + (*(common_vertex_in_order + i))));
			*hessian_address += value;
			*(hessian_address + move_size) += value;
			*(hessian_address + (move_size << 1)) += value;

		}


		//double value;
		//double* A_1;
		//double* A_2;
		//for (int i = 0; i < common_vertex_num; ++i) {
		//	A_1 = A.data() + ((*(common_vertex_in_order + i)) * 3);
		//	for (int j = i + 1; j < common_vertex_num; ++j) {
		//		A_2 = A.data() + ((*(common_vertex_in_order + j)) * 3);
		//		value = 2.0 * stiffness * DOT(A_1, A_2);
		//		sys_matrix(3 * (*(common_vertex_in_order + j + common_vertex_num)),
		//			3 * (*(common_vertex_in_order + i + common_vertex_num))) += value;
		//		sys_matrix(3 * (*(common_vertex_in_order + j + common_vertex_num)) + 1,
		//			3 * (*(common_vertex_in_order + i + common_vertex_num)) + 1) += value;
		//		sys_matrix(3 * (*(common_vertex_in_order + j + common_vertex_num)) + 2,
		//			3 * (*(common_vertex_in_order + i + common_vertex_num)) + 2) += value;

		//		sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)),
		//			3 * (*(common_vertex_in_order + j + common_vertex_num))) += value;
		//		sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)) + 1,
		//			3 * (*(common_vertex_in_order + j + common_vertex_num)) + 1) += value;
		//		sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)) + 2,
		//			3 * (*(common_vertex_in_order + j + common_vertex_num)) + 2) += value;
		//	}

		//	value = 2.0 * stiffness * DOT(A_1, A_1);
		//	sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)),
		//		3 * (*(common_vertex_in_order + i + common_vertex_num))) += value;
		//	sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)) + 1,
		//		3 * (*(common_vertex_in_order + i + common_vertex_num)) + 1) += value;
		//	sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)) + 2,
		//		3 * (*(common_vertex_in_order + i + common_vertex_num)) + 2) += value;
		//}
	}

	common_vertex_in_order += (common_vertex_num << 1);
}

bool SecondOrderConstraint::solveCertainHessianForNeighborTet(std::array<double, 3>* vertex_position, double stiffness,
	Matrix<double, 3, 4>& A, unsigned int* &common_vertex_in_order, int* neighbor_tet_vetex_indices, 
	MatrixXd& sys_matrix, double volume, VectorXd& grad)
{// for a single neighbor tet, only compute common vertices
	if (*common_vertex_in_order == 0) {
		common_vertex_in_order++;
		return false;
	}
	int common_vertex_num = *common_vertex_in_order;
	common_vertex_in_order++;


	Matrix<double, 12, 3> Hessian_vertex;

	Matrix3d deformation_gradient;
	Matrix3d rotation;


	FEM::getDeformationGradient(vertex_position[neighbor_tet_vetex_indices[0]].data(), vertex_position[neighbor_tet_vetex_indices[1]].data(),
		vertex_position[neighbor_tet_vetex_indices[2]].data(), vertex_position[neighbor_tet_vetex_indices[3]].data(), A,
		deformation_gradient);


	stiffness *= volume;

	//if ((deformation_gradient - rotation).squaredNorm() < 1e-16) {
		//common_vertex_in_order += (common_vertex_num * 2);
		//return false;
	//}





	if (solve_exact_ARAP_hessian) {
		Matrix3d S;	Vector3d eigen_value;
		FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);

		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (2.0 * stiffness) * (deformation_gradient - rotation) * A;

		for (int i = 0; i < common_vertex_num; ++i) {
			grad.segment(3 * (*(common_vertex_in_order + i + common_vertex_num)), 3)
				+= grad_C_transpose.col(*(common_vertex_in_order + i));
		}

		Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
		for (int i = 0; i < common_vertex_num; ++i) {
			FEM::getHessianForOneVertex(Hessian_vertex, S, rotation, Dm, A, *(common_vertex_in_order + i));
			Hessian_vertex *= stiffness;
			for (int j = 0; j < common_vertex_num; ++j) {
				sys_matrix.block<3, 3>(3 * (*(common_vertex_in_order + j + common_vertex_num)),
					3 * (*(common_vertex_in_order + i + common_vertex_num))) +=
					Hessian_vertex.block<3, 3>(3 * (*(common_vertex_in_order + j)), 0);
			}
		}
	}
	else {
		FEM::extractRotation(deformation_gradient, rotation);

		//if (rotation.determinant() < 0.9 || rotation.determinant() > 1.1) {
		//	std::cout << "rotation error " << rotation.determinant() << std::endl;
		//}

		Matrix<double, 3, 4> grad_C_transpose;
		grad_C_transpose = (2.0 * stiffness) * (deformation_gradient - rotation) * A;

		for (int i = 0; i < common_vertex_num; ++i) {
			grad.segment(3 * (*(common_vertex_in_order + i + common_vertex_num)), 3)
				+= grad_C_transpose.col(*(common_vertex_in_order + i));
		}

		double value;
		double* A_1;
		double* A_2;
		for (int i = 0; i < common_vertex_num; ++i) {
			A_1 = A.data() + ((*(common_vertex_in_order + i)) * 3);
			for (int j = i + 1; j < common_vertex_num; ++j) {
				A_2 = A.data() + ((*(common_vertex_in_order + j)) * 3);
				value = 2.0 * stiffness * DOT(A_1, A_2);
				sys_matrix(3 * (*(common_vertex_in_order + j + common_vertex_num)),
					3 * (*(common_vertex_in_order + i + common_vertex_num))) += value;
				sys_matrix(3 * (*(common_vertex_in_order + j + common_vertex_num))+1,
					3 * (*(common_vertex_in_order + i + common_vertex_num))+1) += value;
				sys_matrix(3 * (*(common_vertex_in_order + j + common_vertex_num))+2,
					3 * (*(common_vertex_in_order + i + common_vertex_num))+2) += value;

				sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)),
					3 * (*(common_vertex_in_order + j + common_vertex_num))) += value;
				sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)) + 1,
					3 * (*(common_vertex_in_order + j + common_vertex_num)) + 1) += value;
				sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)) + 2,
					3 * (*(common_vertex_in_order + j + common_vertex_num)) + 2) += value;
			}

			value = 2.0 * stiffness * DOT(A_1, A_1);
			sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)),
				3 * (*(common_vertex_in_order + i + common_vertex_num))) += value;
			sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)) + 1,
				3 * (*(common_vertex_in_order + i + common_vertex_num)) + 1) += value;
			sys_matrix(3 * (*(common_vertex_in_order + i + common_vertex_num)) + 2,
				3 * (*(common_vertex_in_order + i + common_vertex_num)) + 2) += value;

		}
	}

	common_vertex_in_order += (common_vertex_num << 1);
	return true;
}

void SecondOrderConstraint::solveCD_ARAP(std::array<double, 3>* vertex_position, double stiffness, double dt,
	Matrix<double, 3, 4>* A, double* lambda, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass, 
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn)
{
	unsigned int tet_index;
	Matrix3d Hessian_single;
	Vector3d grad_single;
	unsigned int vertex_no;
	MatrixXd Hessian;
	Hessian.resize(3 + tet_indices.size(), 3 + tet_indices.size());
	Hessian.setZero();

	VectorXd b;
	b.resize(3 + tet_indices.size());
	b.setZero();
	double C;
	double alpha;

	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		tet_index = tet_indices[i];
		vertex_no = findVertexNo(vertex_index, indices[tet_index].data(),4);
		if (getARAPGradHessian(vertex_position[indices[tet_index][0]].data(), vertex_position[indices[tet_index][1]].data(),
			vertex_position[indices[tet_index][2]].data(), vertex_position[indices[tet_index][3]].data(),
			A[tet_index],  Hessian_single, grad_single, C, vertex_no)) {
			Hessian_single *= -lambda[tet_index];
			Hessian.block<3, 3>(0, 0) += Hessian_single;
			Hessian.block<3, 1>(0, 3 + i) = -grad_single;
			Hessian.block<1, 3>(3 + i, 0) = -grad_single;
			b.head(3) += lambda[tet_index] * grad_single;
		}		
		alpha = 1.0 / (stiffness * volume[tet_index] * dt * dt);
		Hessian(3 + i, 3 + i) = -alpha;
		b.data()[3 + i] = C + alpha * lambda[tet_index];
	}
	
	Hessian.data()[0] += mass[vertex_index];
	Hessian(1,1) += mass[vertex_index];
	Hessian(2,2) += mass[vertex_index];
	b.data()[0] -= mass[vertex_index] * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	b.data()[1] -= mass[vertex_index] * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	b.data()[2] -= mass[vertex_index] * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(b);

	SUM_(vertex_position[vertex_index], result);
	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		lambda[tet_indices[i]] += result[3 + i];
	}
}




void SecondOrderConstraint::setBarrierGHWithMollifier(double barrier_,MatrixXd& dis_h, VectorXd& dis_g,
	double* ea0, double* ea1, double* eb0, double* eb1, double eps_x, double ee_cross_norm_2,
	double mollifier, double b_grad, double b_hessian, bool compute_hessian, bool is_collider)
{
	VectorXd mollifier_grad; mollifier_grad.resize(12);	
	distance::edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1, mollifier_grad.data());
	double K = CCD::internal::edgeEdgeMollifierGradient(ee_cross_norm_2, eps_x);

	if (compute_hessian)
	{
		MatrixXd mollifier_hessian(12, 12);
		distance::edge_edge_cross_squarednorm_hessian(ea0, ea1, eb0, eb1, mollifier_hessian.data());
		mollifier_hessian *= K;
		mollifier_hessian += ((-2.0 / (eps_x * eps_x)) * mollifier_grad) * mollifier_grad.transpose();
		mollifier_grad *= K;
		dis_h = (mollifier_hessian * barrier_
			+ b_grad * (dis_g * mollifier_grad.transpose() + mollifier_grad * dis_g.transpose())
			+ mollifier * (b_hessian * dis_g * dis_g.transpose() + b_grad * dis_h)).eval();

		if (is_collider) {
			MatrixXd temp = dis_h.block<6, 6>(0, 0);
			FEM::SPDprojection(temp);
			dis_h.block<6, 6>(0, 0) = temp;
		}
		else {
			FEM::SPDprojection(dis_h);
		}
		

	}
	else {
		mollifier_grad *= K;
	}

	dis_g = mollifier_grad * barrier_ + (mollifier * b_grad) * dis_g;
}


void SecondOrderConstraint::computeEEBarrierGradientHessian(double* ea0, double* ea1, double* eb0, double* eb1, MatrixXd& Hessian_, VectorXd& grad_,
	int* vertex_order_in_system, double stiffness, double d_hat_2, double rest_length_0, double rest_length_1)
{
	double distance;
	double b_grad, b_hessian;
	MatrixXd h; VectorXd g;
	int vertex_in_pair[4];
	memset(vertex_in_pair, 0, 16);


	double eps_x = 1e-3 * rest_length_0 * rest_length_0 * rest_length_1 * rest_length_1;
	double mollifier;
	double ee_cross_norm_2;
	bool need_mollifier = CCD::internal::edgeEdgeMollifier(ea0, ea1, eb0, eb1, eps_x, mollifier, ee_cross_norm_2);

	switch (CCD::internal::edgeEdgeDistanceType(ea0, ea1, eb0, eb1)) {
	case 0:
		distance = CCD::internal::pointPointDistance(ea0, eb0);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea0, eb0, g.data());
		distance::point_point_distance_hessian(ea0, eb0, h.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 2);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian,true,false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}
		break;

	case 1:
		distance = CCD::internal::pointPointDistance(ea0, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea0, eb1, g.data());
		distance::point_point_distance_hessian(ea0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 3;


		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 2);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);
			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}
		break;
	case 2:

		distance = CCD::internal::pointEdgeDistance(ea0, eb0, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(ea0, eb0, eb1, g.data());
		distance::point_edge_distance_hessian(ea0, eb0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;
		vertex_in_pair[2] = 3;


		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 3);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}
		break;

	case 3:
		distance = CCD::internal::pointPointDistance(ea1, eb0);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea1, eb0, g.data());
		distance::point_point_distance_hessian(ea1, eb0, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 1;
		vertex_in_pair[1] = 2;

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 2);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);
			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}

		break;
	case 4:
		distance = CCD::internal::pointPointDistance(ea1, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea1, eb1, g.data());
		distance::point_point_distance_hessian(ea1, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 1;
		vertex_in_pair[1] = 3;

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 2);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);
			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}
		break;
	case 5:

		distance = CCD::internal::pointEdgeDistance(ea1, eb0, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(ea1, eb0, eb1, g.data());
		distance::point_edge_distance_hessian(ea1, eb0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 1;
		vertex_in_pair[1] = 2;
		vertex_in_pair[2] = 3;


		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 3);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}

		break;
	case 6:
		distance = CCD::internal::pointEdgeDistance(eb0, ea0, ea1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(eb0, ea0, ea1, g.data());
		distance::point_edge_distance_hessian(eb0, ea0, ea1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 3);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}
		break;
	case 7:

		distance = CCD::internal::pointEdgeDistance(eb1, ea0, ea1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(eb1, ea0, ea1, g.data());
		distance::point_edge_distance_hessian(eb1, ea0, ea1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;


		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 3);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}

		break;
	case 8:
		distance = CCD::internal::edgeEdgeDistance(ea0, ea1, eb0, eb1);

		if (distance >= d_hat_2) {
			return;
		}


		h.resize(12, 12);
		g.resize(12);

		distance::edge_edge_distance_gradient(ea0, ea1, eb0, eb1, g.data());
		distance::edge_edge_distance_hessian(ea0, ea1, eb0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;
		vertex_in_pair[3] = 3;


		if (need_mollifier) {
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, h, g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);
			setTetHessianFromHessian(Hessian_, grad_, h, g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 4);
		}
		break;

	}
}



void SecondOrderConstraint::computeEEBarrierGradientHessianTest(double* ea0, double* ea1, double* eb0, double* eb1, MatrixXd& Hessian_, VectorXd& grad_,
	int* vertex_order_in_system, double stiffness, double d_hat_2, double rest_length_0, 
	double rest_length_1, int num_record, int* vertex_in_pair_, int type,
	int obj_0, int edge_index_0, int obj_1, int edge_index_1)
{
	double distance;
	double b_grad, b_hessian;
	MatrixXd h; VectorXd g;
	int vertex_in_pair[4];
	memset(vertex_in_pair, 0, 16);


	double eps_x = 1e-3 * rest_length_0 * rest_length_0 * rest_length_1 * rest_length_1;
	double mollifier;
	double ee_cross_norm_2;
	bool need_mollifier = CCD::internal::edgeEdgeMollifier(ea0, ea1, eb0, eb1, eps_x, mollifier, ee_cross_norm_2);

	int use_pair_num = 0;

	switch (CCD::internal::edgeEdgeDistanceType(ea0, ea1, eb0, eb1)) {
	case 0:
		distance = CCD::internal::pointPointDistance(ea0, eb0);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea0, eb0, g.data());
		distance::point_point_distance_hessian(ea0, eb0, h.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 2);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}

		use_pair_num = 2;

		break;

	case 1:
		distance = CCD::internal::pointPointDistance(ea0, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea0, eb1, g.data());
		distance::point_point_distance_hessian(ea0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 3;


		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				 vertex_in_pair, 2);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);
			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}

		use_pair_num = 2;
		break;
	case 2:

		distance = CCD::internal::pointEdgeDistance(ea0, eb0, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(ea0, eb0, eb1, g.data());
		distance::point_edge_distance_hessian(ea0, eb0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;
		vertex_in_pair[2] = 3;


		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				 vertex_in_pair, 3);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}


		use_pair_num = 3;

		break;

	case 3:
		distance = CCD::internal::pointPointDistance(ea1, eb0);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea1, eb0, g.data());
		distance::point_point_distance_hessian(ea1, eb0, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 1;
		vertex_in_pair[1] = 2;

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 2);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);
			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}

		use_pair_num = 2;

		break;
	case 4:
		distance = CCD::internal::pointPointDistance(ea1, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea1, eb1, g.data());
		distance::point_point_distance_hessian(ea1, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 1;
		vertex_in_pair[1] = 3;

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 2);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);
			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}

		use_pair_num = 2;

		break;
	case 5:

		distance = CCD::internal::pointEdgeDistance(ea1, eb0, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(ea1, eb0, eb1, g.data());
		distance::point_edge_distance_hessian(ea1, eb0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 1;
		vertex_in_pair[1] = 2;
		vertex_in_pair[2] = 3;


		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 3);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}

		use_pair_num = 3;

		break;
	case 6:
		distance = CCD::internal::pointEdgeDistance(eb0, ea0, ea1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(eb0, ea0, ea1, g.data());
		distance::point_edge_distance_hessian(eb0, ea0, ea1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 3);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}

		use_pair_num = 3;

		break;
	case 7:

		distance = CCD::internal::pointEdgeDistance(eb1, ea0, ea1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(eb1, ea0, ea1, g.data());
		distance::point_edge_distance_hessian(eb1, ea0, ea1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;


		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 3);
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}

		use_pair_num = 3;

		break;
	case 8:
		distance = CCD::internal::edgeEdgeDistance(ea0, ea1, eb0, eb1);

		if (distance >= d_hat_2) {
			return;
		}


		h.resize(12, 12);
		g.resize(12);

		distance::edge_edge_distance_gradient(ea0, ea1, eb0, eb1, g.data());
		distance::edge_edge_distance_hessian(ea0, ea1, eb0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;
		vertex_in_pair[3] = 3;


		if (need_mollifier) {
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, h, g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true, false);
			setTetHessianFromHessian(Hessian_, grad_, h, g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 4);
		}

		use_pair_num = 4;

		break;

	}
	if (type ==0) {
		int num_sued=0;
		int record[4];
		for (int i = 0; i < use_pair_num; ++i) {
			if (vertex_in_pair[i] < 2) {
				record[num_sued] = vertex_in_pair[i];
				num_sued++;
			}
		}
		if (num_record != num_sued) {
			std::cout << type << " ee collider  use_pair_num wrong, right: " << num_sued << " " << num_record << std::endl;
			std::cout << obj_0 << " " << edge_index_0 << " " << obj_1 << " " << edge_index_1 << std::endl;
		}
		for (int i = 0; i < num_record; ++i) {
			if (record[i] != vertex_in_pair_[i]) {
				std::cout << "ee collider record_vertex_in_use wrong , right " << record[i] << " " << vertex_in_pair_[i] << std::endl;
			}
		}
	}
	else {
		if (num_record != use_pair_num) {
			std::cout <<type<< " ee use_pair_num wrong, right: " << use_pair_num << " " << num_record << std::endl;
			std::cout << obj_0 << " " << edge_index_0 << " " << obj_1 << " " << edge_index_1 << std::endl;
		}
		for (int i = 0; i < num_record; ++i) {
			if (vertex_in_pair[i] != vertex_in_pair_[i]) {
				std::cout << "ee record_vertex_in_use wrong , right " << vertex_in_pair[i] << " " << vertex_in_pair_[i] << std::endl;
			}
		}		
	}


}
bool SecondOrderConstraint::computeBarrierEEGradientHessian(double* ea0, double* ea1, double* eb0, double* eb1, MatrixXd& h, VectorXd& g,
	int* vertex_in_pair, double stiffness, double d_hat_2, double rest_length_0, double rest_length_1, bool compute_hessian,
	bool is_collider)
{
	double distance;
	double b_grad, b_hessian;

	double eps_x = 1e-3 * rest_length_0 * rest_length_0 * rest_length_1 * rest_length_1;
	double mollifier;
	double ee_cross_norm_2;
	bool need_mollifier = CCD::internal::edgeEdgeMollifier(ea0, ea1, eb0, eb1, eps_x, mollifier, ee_cross_norm_2);

	switch (CCD::internal::edgeEdgeDistanceType(ea0, ea1, eb0, eb1)) {
	case 0:
		distance = CCD::internal::pointPointDistance(ea0, eb0);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
		g.resize(6);
		distance::point_point_distance_gradient(ea0, eb0, g.data());
		if (compute_hessian)
		{
			h.resize(6, 6);
			distance::point_point_distance_hessian(ea0, eb0, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;
		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 2;
		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian,is_collider);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g, vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				if (is_collider) {
					MatrixXd temp = h.block<3, 3>(0, 0);
					FEM::SPDprojection(temp);
					h.block<3, 3>(0, 0) = temp;
				}
				else {
					FEM::SPDprojection(h);
				}

			}
			g *= b_grad;
		}
		break;

	case 1:
		distance = CCD::internal::pointPointDistance(ea0, eb1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}

		g.resize(6);
		distance::point_point_distance_gradient(ea0, eb1, g.data());

		if (compute_hessian) {
			h.resize(6, 6);
			distance::point_point_distance_hessian(ea0, eb1, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}

		b_grad *= stiffness;
		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 3;

		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;

		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				if (is_collider) {
					MatrixXd temp = h.block<3, 3>(0, 0);
					FEM::SPDprojection(temp);
					h.block<3, 3>(0, 0) = temp;
				}
				else {
					FEM::SPDprojection(h);
				}
			}
			g *= b_grad;
		}
		break;
	case 2:

		distance = CCD::internal::pointEdgeDistance(ea0, eb0, eb1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}

		g.resize(9);
		distance::point_edge_distance_gradient(ea0, eb0, eb1, g.data());
		if (compute_hessian)
		{
			h.resize(9, 9);
			distance::point_edge_distance_hessian(ea0, eb0, eb1, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;

		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 2;
		vertex_in_pair[3] = 3;

		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();

				if (is_collider) {
					MatrixXd temp = h.block<3, 3>(0, 0);
					FEM::SPDprojection(temp);
					h.block<3, 3>(0, 0) = temp;
				}
				else {
					FEM::SPDprojection(h);
				}

			}
			g *= b_grad;

		}
		break;

	case 3:
		distance = CCD::internal::pointPointDistance(ea1, eb0);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}

		g.resize(6);
		distance::point_point_distance_gradient(ea1, eb0, g.data());
		if (compute_hessian)
		{
			h.resize(6, 6);
			distance::point_point_distance_hessian(ea1, eb0, h.data());

			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;

		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;

		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;

		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();

				if (is_collider) {
					MatrixXd temp = h.block<3, 3>(0, 0);
					FEM::SPDprojection(temp);
					h.block<3, 3>(0, 0) = temp;
				}
				else {
					FEM::SPDprojection(h);
				}
			}
			g *= b_grad;

		}

		break;
	case 4:
		distance = CCD::internal::pointPointDistance(ea1, eb1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}

		g.resize(6);
		distance::point_point_distance_gradient(ea1, eb1, g.data());

		if (compute_hessian)
		{
			h.resize(6, 6);
			distance::point_point_distance_hessian(ea1, eb1, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;
		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 3;

		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();

				if (is_collider) {
					MatrixXd temp = h.block<3, 3>(0, 0);
					FEM::SPDprojection(temp);
					h.block<3, 3>(0, 0) = temp;
				}
				else {
					FEM::SPDprojection(h);
				}
			}
			g *= b_grad;

		}
		break;
	case 5:

		distance = CCD::internal::pointEdgeDistance(ea1, eb0, eb1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}

		g.resize(9);

		distance::point_edge_distance_gradient(ea1, eb0, eb1, g.data());

		if (compute_hessian)
		{
			h.resize(9, 9);
			distance::point_edge_distance_hessian(ea1, eb0, eb1, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;
		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;
		vertex_in_pair[3] = 3;


		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				if (is_collider) {
					MatrixXd temp = h.block<3, 3>(0, 0);
					FEM::SPDprojection(temp);
					h.block<3, 3>(0, 0) = temp;
				}
				else {
					FEM::SPDprojection(h);
				}
			}
			g *= b_grad;

		}

		break;
	case 6:
		distance = CCD::internal::pointEdgeDistance(eb0, ea0, ea1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}

		g.resize(9);
		distance::point_edge_distance_gradient(eb0, ea0, ea1, g.data());
		if (compute_hessian)
		{
			h.resize(9, 9);
			distance::point_edge_distance_hessian(eb0, ea0, ea1, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;
		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 2;
		vertex_in_pair[2] = 0;
		vertex_in_pair[3] = 1;

		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;

		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				if (is_collider) {
					MatrixXd temp = h.block<6, 6>(3, 3);
					FEM::SPDprojection(temp);
					h.block<6, 6>(3, 3)= temp;
				}
				else {
					FEM::SPDprojection(h);
				}
			}
			g *= b_grad;

		}
		break;
	case 7:

		distance = CCD::internal::pointEdgeDistance(eb1, ea0, ea1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}

		g.resize(9);

		distance::point_edge_distance_gradient(eb1, ea0, ea1, g.data());
		if (compute_hessian)
		{
			h.resize(9, 9);
			distance::point_edge_distance_hessian(eb1, ea0, ea1, h.data());

			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;

		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 3;
		vertex_in_pair[2] = 0;
		vertex_in_pair[3] = 1;


		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				if (is_collider) {
					MatrixXd temp = h.block<6, 6>(3, 3);
					FEM::SPDprojection(temp);
					h.block<6, 6>(3, 3) = temp;
				}
				else {
					FEM::SPDprojection(h);
				}
			}
			g *= b_grad;
		}

		break;
	case 8:
		distance = CCD::internal::edgeEdgeDistance(ea0, ea1, eb0, eb1);

		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
		g.resize(12);
		distance::edge_edge_distance_gradient(ea0, ea1, eb0, eb1, g.data());
		if (compute_hessian) {
			h.resize(12, 12);
			distance::edge_edge_distance_hessian(ea0, ea1, eb0, eb1, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;
		vertex_in_pair[0] = 4;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;
		vertex_in_pair[3] = 2;
		vertex_in_pair[4] = 3;

		if (need_mollifier) {
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, h, g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, is_collider);
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				if (is_collider) {
					MatrixXd temp = h.block<6, 6>(0, 0);
					FEM::SPDprojection(temp);
					h.block<6, 6>(0, 0) = temp;
				}
				else {
					FEM::SPDprojection(h);
				}
			}
			g *= b_grad;
		}
		break;
	}

	return true;
}



bool SecondOrderConstraint::computeBarrierEEGradientHessianTest(double* ea0, double* ea1, double* eb0, double* eb1, MatrixXd& h, VectorXd& g,
	int* vertex_in_pair, double stiffness, double d_hat_2, double rest_length_0, double rest_length_1, bool compute_hessian)
{
	double distance;
	double b_grad, b_hessian;

	double eps_x = 1e-3 * rest_length_0 * rest_length_0 * rest_length_1 * rest_length_1;
	double mollifier;
	double ee_cross_norm_2;
	bool need_mollifier = CCD::internal::edgeEdgeMollifier(ea0, ea1, eb0, eb1, eps_x, mollifier, ee_cross_norm_2);

	std::cout<<" type " << CCD::internal::edgeEdgeDistanceType(ea0, ea1, eb0, eb1) << std::endl;

	switch (CCD::internal::edgeEdgeDistanceType(ea0, ea1, eb0, eb1)) {
	case 0:
		distance = CCD::internal::pointPointDistance(ea0, eb0);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
		g.resize(6);
		distance::point_point_distance_gradient(ea0, eb0, g.data());
		if (compute_hessian)
		{
			h.resize(6, 6);
			distance::point_point_distance_hessian(ea0, eb0, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;
		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 2;
		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, false);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g, vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, false);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				FEM::SPDprojection(h);
			}
			g *= b_grad;
		}
		break;

	case 1:
		distance = CCD::internal::pointPointDistance(ea0, eb1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}

		g.resize(6);
		distance::point_point_distance_gradient(ea0, eb1, g.data());

		if (compute_hessian) {
			h.resize(6, 6);
			distance::point_point_distance_hessian(ea0, eb1, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}

		b_grad *= stiffness;
		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 3;

		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, false);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, false);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;	
			g = dis_g;

		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				FEM::SPDprojection(h);
			}
			g *= b_grad;
		}
		break;
	case 2:

		distance = CCD::internal::pointEdgeDistance(ea0, eb0, eb1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
		
		g.resize(9);
		distance::point_edge_distance_gradient(ea0, eb0, eb1, g.data());
		if (compute_hessian)
		{
			h.resize(9, 9);
			distance::point_edge_distance_hessian(ea0, eb0, eb1, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;

		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 2;
		vertex_in_pair[3] = 3;

		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian,compute_hessian, false);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier( dis_g, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, false);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				MatrixXd temp = h.block<3, 3>(0, 0);				
				std::cout << "ori " << std::endl;
				std::cout << h << std::endl;
				FEM::SPDprojection(temp);
				h.block<3, 3>(0, 0)=temp;
			}
			g *= b_grad;

		}
		break;

	case 3:
		distance = CCD::internal::pointPointDistance(ea1, eb0);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
		
		g.resize(6);
		distance::point_point_distance_gradient(ea1, eb0, g.data());
		if (compute_hessian)
		{
			h.resize(6, 6);
			distance::point_point_distance_hessian(ea1, eb0, h.data());

			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;

		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;

		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, false);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, false);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;

		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				FEM::SPDprojection(h);
			}
			g *= b_grad;
			
		}

		break;
	case 4:
		distance = CCD::internal::pointPointDistance(ea1, eb1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
		
		g.resize(6);
		distance::point_point_distance_gradient(ea1, eb1, g.data());

		if (compute_hessian)
		{
			h.resize(6, 6);
			distance::point_point_distance_hessian(ea1, eb1, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			 b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;
		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 3;

		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian,compute_hessian, false);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier( dis_g, g,
					vertex_in_pair + 1, 2);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, false);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;			
			g = dis_g;
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				FEM::SPDprojection(h);
			}
			g *= b_grad;

		}
		break;
	case 5:

		distance = CCD::internal::pointEdgeDistance(ea1, eb0, eb1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
		
		g.resize(9);

		distance::point_edge_distance_gradient(ea1, eb0, eb1, g.data());

		if (compute_hessian)
		{
			h.resize(9, 9);
			distance::point_edge_distance_hessian(ea1, eb0, eb1, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);			
		}
		b_grad *= stiffness;
		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;
		vertex_in_pair[3] = 3;


		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian,compute_hessian, false);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier(dis_g, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, false);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;		
			g = dis_g;
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				FEM::SPDprojection(h);
			}
			g *= b_grad;

		}

		break;
	case 6:
		distance = CCD::internal::pointEdgeDistance(eb0, ea0, ea1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
	
		g.resize(9);
		distance::point_edge_distance_gradient(eb0, ea0, ea1, g.data());
		if (compute_hessian)
		{
			h.resize(9, 9);
			distance::point_edge_distance_hessian(eb0, ea0, ea1, h.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;
		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 2;
		vertex_in_pair[2] = 0;
		vertex_in_pair[3] = 1;

		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian,compute_hessian, false);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier( dis_g, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, false);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;

		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				FEM::SPDprojection(h);
			}
			g *= b_grad;

		}
		break;
	case 7:

		distance = CCD::internal::pointEdgeDistance(eb1, ea0, ea1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
	
		g.resize(9);

		distance::point_edge_distance_gradient(eb1, ea0, ea1, g.data());
		if (compute_hessian)
		{
			h.resize(9, 9);
			distance::point_edge_distance_hessian(eb1, ea0, ea1, h.data());

			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_hessian *= stiffness;
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;

		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 3;
		vertex_in_pair[2] = 0;
		vertex_in_pair[3] = 1;


		if (need_mollifier) {
			VectorXd dis_g(12);
			dis_g.setZero();
			MatrixXd dis_h;
			if (compute_hessian) {
				dis_h.resize(12, 12);
				dis_h.setZero();
				setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian,compute_hessian, false);
				h = dis_h;
			}
			else {
				setFourVertexGradFromBarrier( dis_g, g,
					vertex_in_pair + 1, 3);
				double barrier_ = stiffness * barrier(distance, d_hat_2);
				setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
					b_grad, b_hessian, compute_hessian, false);
			}
			vertex_in_pair[0] = 4;
			vertex_in_pair[1] = 0;
			vertex_in_pair[2] = 1;
			vertex_in_pair[3] = 2;
			vertex_in_pair[4] = 3;
			g = dis_g;
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				FEM::SPDprojection(h);
			}
			g *= b_grad;
		}

		break;
	case 8:
		distance = CCD::internal::edgeEdgeDistance(ea0, ea1, eb0, eb1);

		std::cout <<"distance "<< distance <<" "<< d_hat_2 << std::endl;

		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}

		g.resize(12);
		//distance::edge_edge_distance_gradient(ea0, ea1, eb0, eb1, g.data());
		if (compute_hessian) {
			h.resize(12, 12);

			double mat[4][4][9];
				distance::edge_edge_distance_all_grad_hess(ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2], eb1[0], eb1[1], eb1[2],
					g.data(), g.data() + 3, g.data() + 6, g.data() + 9,mat[0][0], mat[0][1], mat[0][2],
					mat[0][3], mat[1][0], mat[1][1], mat[1][2], mat[1][3],
					mat[2][0], mat[2][1], mat[2][2], mat[2][3], mat[3][0], mat[3][1], mat[3][2], mat[3][3]);

				for (int k = 0; k < 4; ++k) {
					for (int l = 0; l < 4; ++l) {
						for (int p0 = 0; p0 < 3; ++p0) {
							for (int p1 = 0; p1 < 3; ++p1) {
								h.data()[12 * (3 * l + p1) + (3 * k + p0)] = mat[k][l][3 * p1 + p0];
							}
						}
					}
				}

				MatrixXd kk(12,12);
			
			distance::H_EE(ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2], eb1[0], eb1[1], eb1[2], kk.data());

			std::cout << kk - h << std::endl;

			std::cout << h << std::endl;

			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			 b_hessian *= stiffness;

			 std::cout << "b_hessian " << b_hessian<<" "<< need_mollifier << std::endl;

		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
		}
		b_grad *= stiffness;
		vertex_in_pair[0] = 4;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;
		vertex_in_pair[3] = 2;
		vertex_in_pair[4] = 3;

		if (need_mollifier) {
			double barrier_ = stiffness * barrier(distance, d_hat_2);
			setBarrierGHWithMollifier(barrier_, h, g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, compute_hessian, false);
		}
		else {
			if (compute_hessian) {
				h = (b_hessian * g * g.transpose() + b_grad * h).eval();
				FEM::SPDprojection(h);
			}
			g *= b_grad;
		}
		break;
	}

	return true;
}




void SecondOrderConstraint::computeEEBarrierGradientHessianTest(double* ea0, double* ea1, double* eb0, double* eb1, MatrixXd& Hessian_, VectorXd& grad_,
	int* vertex_order_in_system, double stiffness, double d_hat_2, double rest_length_0, double rest_length_1, double& barrier_)
{
	double distance;
	double b_grad, b_hessian;
	MatrixXd h; VectorXd g;
	int vertex_in_pair[4];
	memset(vertex_in_pair, 0, 16);


	double eps_x = 1e-3 * rest_length_0 * rest_length_0 * rest_length_1 * rest_length_1;
	double mollifier;
	double ee_cross_norm_2;
	bool need_mollifier = CCD::internal::edgeEdgeMollifier(ea0, ea1, eb0, eb1, eps_x, mollifier, ee_cross_norm_2);

	switch (CCD::internal::edgeEdgeDistanceType(ea0, ea1, eb0, eb1)) {
	case 0:
		distance = CCD::internal::pointPointDistance(ea0, eb0);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea0, eb0, g.data());
		distance::point_point_distance_hessian(ea0, eb0, h.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;

		barrier_ = stiffness * barrier(distance, d_hat_2);
		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 2);

			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian,true,false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}
		break;

	case 1:
		distance = CCD::internal::pointPointDistance(ea0, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea0, eb1, g.data());
		distance::point_point_distance_hessian(ea0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 3;

		barrier_ = stiffness * barrier(distance, d_hat_2);

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 2);

			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true,false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}
		break;
	case 2:

		distance = CCD::internal::pointEdgeDistance(ea0, eb0, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(ea0, eb0, eb1, g.data());
		distance::point_edge_distance_hessian(ea0, eb0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;
		vertex_in_pair[2] = 3;

		barrier_ = stiffness * barrier(distance, d_hat_2);

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 3);

			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true,false);
			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}
		break;

	case 3:
		distance = CCD::internal::pointPointDistance(ea1, eb0);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea1, eb0, g.data());
		distance::point_point_distance_hessian(ea1, eb0, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 1;
		vertex_in_pair[1] = 2;

		barrier_ = stiffness * barrier(distance, d_hat_2);

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 2);

			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true,false);
			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}

		break;
	case 4:
		distance = CCD::internal::pointPointDistance(ea1, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(6, 6);
		g.resize(6);

		distance::point_point_distance_gradient(ea1, eb1, g.data());
		distance::point_point_distance_hessian(ea1, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 1;
		vertex_in_pair[1] = 3;
		barrier_ = stiffness * barrier(distance, d_hat_2);
		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 2);

			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true,false);
			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 2);
		}
		break;
	case 5:

		distance = CCD::internal::pointEdgeDistance(ea1, eb0, eb1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(ea1, eb0, eb1, g.data());
		distance::point_edge_distance_hessian(ea1, eb0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 1;
		vertex_in_pair[1] = 2;
		vertex_in_pair[2] = 3;
		barrier_ = stiffness * barrier(distance, d_hat_2);

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 3);
			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true,false);
			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}

		break;
	case 6:
		distance = CCD::internal::pointEdgeDistance(eb0, ea0, ea1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(eb0, ea0, ea1, g.data());
		distance::point_edge_distance_hessian(eb0, ea0, ea1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;
		barrier_ = stiffness * barrier(distance, d_hat_2);
		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 3);

			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true,false);
			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}
		break;
	case 7:

		distance = CCD::internal::pointEdgeDistance(eb1, ea0, ea1);
		if (distance >= d_hat_2) {
			return;
		}
		h.resize(9, 9);
		g.resize(9);

		distance::point_edge_distance_gradient(eb1, ea0, ea1, g.data());
		distance::point_edge_distance_hessian(eb1, ea0, ea1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;
		barrier_ = stiffness * barrier(distance, d_hat_2);

		if (need_mollifier) {
			MatrixXd dis_h(12, 12); VectorXd dis_g(12);
			dis_h.setZero(); dis_g.setZero();
			setFourVertexHessianFromBarrierHessian(dis_h, dis_g, h, g,
				vertex_in_pair, 3);

			setBarrierGHWithMollifier(barrier_, dis_h, dis_g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true,false);

			setTetHessianFromHessian(Hessian_, grad_, dis_h, dis_g, vertex_order_in_system);
		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 3);
		}

		break;
	case 8:
		distance = CCD::internal::edgeEdgeDistance(ea0, ea1, eb0, eb1);

		if (distance >= d_hat_2) {
			return;
		}


		h.resize(12, 12);
		g.resize(12);

		distance::edge_edge_distance_gradient(ea0, ea1, eb0, eb1, g.data());
		distance::edge_edge_distance_hessian(ea0, ea1, eb0, eb1, h.data());

		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;
		vertex_in_pair[3] = 3;

		barrier_ = stiffness * barrier(distance, d_hat_2);

		if (need_mollifier) {

			setBarrierGHWithMollifier(barrier_, h, g, ea0, ea1, eb0, eb1, eps_x, ee_cross_norm_2, mollifier,
				b_grad, b_hessian, true,false);
			setTetHessianFromHessian(Hessian_, grad_, h, g, vertex_order_in_system);

		}
		else {
			h = (b_hessian * g * g.transpose() + b_grad * h).eval();
			FEM::SPDprojection(h);
			g *= b_grad;
			setTetHessianFromBarrierHessian(Hessian_, grad_, h, g,
				vertex_order_in_system, vertex_in_pair, 4);
		}
		break;

	}
}



void SecondOrderConstraint::computeVTBarrierGradientHessianTest(MatrixXd& Hessian_, VectorXd& grad_, double* p, double* t0,
	double* t1, double* t2, double d_hat_2, int* triangle_vertex_order_in_system, double stiffness, double& barrier_)
{
	double distance;
	double b_grad, b_hessian;
	MatrixXd Hessian; VectorXd grad;
	int vertex_in_pair[4];
	memset(vertex_in_pair, 0, 16);
	switch (CCD::internal::pointTriangleDistanceType(p, t0, t1, t2))
	{
	case 0:
		distance = CCD::internal::pointPointDistance(p, t0);
		if (distance >= d_hat_2) {
			return;
		}
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;

		Hessian.resize(6, 6);
		grad.resize(6);

		barrier_ = barrier(distance, d_hat_2);

		distance::point_point_distance_gradient(p, t0, grad.data());
		distance::point_point_distance_hessian(p, t0, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;

		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 2);

		break;
	case 1:
		distance = CCD::internal::pointPointDistance(p, t1);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(6, 6);
		grad.resize(6);

		barrier_ = barrier(distance, d_hat_2);

		distance::point_point_distance_gradient(p, t1, grad.data());
		distance::point_point_distance_hessian(p, t1, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 2);

		break;
	case 2:
		distance = CCD::internal::pointPointDistance(p, t2);
		if (distance >= d_hat_2) {
			return;
		}

		barrier_ = barrier(distance, d_hat_2);
		Hessian.resize(6, 6);
		grad.resize(6);
		distance::point_point_distance_gradient(p, t2, grad.data());
		distance::point_point_distance_hessian(p, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 2);

		break;

	case 3:
		distance = CCD::internal::pointEdgeDistance(p, t0, t1);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(9, 9);
		grad.resize(9);

		barrier_ = barrier(distance, d_hat_2);


		distance::point_edge_distance_gradient(p, t0, t1, grad.data());
		distance::point_edge_distance_hessian(p, t0, t1, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 3);

		break;

	case 4:
		distance = CCD::internal::pointEdgeDistance(p, t1, t2);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(9, 9);
		grad.resize(9);
		distance::point_edge_distance_gradient(p, t1, t2, grad.data());
		distance::point_edge_distance_hessian(p, t1, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		barrier_ = barrier(distance, d_hat_2);


		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);


		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;
		vertex_in_pair[2] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 3);

		break;
	case 5:
		distance = CCD::internal::pointEdgeDistance(p, t0, t2);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(9, 9);
		grad.resize(9);

		barrier_ = barrier(distance, d_hat_2);


		distance::point_edge_distance_gradient(p, t0, t2, grad.data());
		distance::point_edge_distance_hessian(p, t0, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 3);
		break;

	case 6:
		distance = CCD::internal::pointTriangleDistance(p, t0, t1, t2);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(12, 12);
		grad.resize(12);
		distance::point_triangle_distance_gradient(p, t0, t1, t2, grad.data());
		distance::point_triangle_distance_hessian(p, t0, t1, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		barrier_ = barrier(distance, d_hat_2);


		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;
		vertex_in_pair[3] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 4);

		break;
	}



}

//type: 0 vt, 1: vt_c, 2: tv_c
bool SecondOrderConstraint::computeBarrierVTGradientHessian(MatrixXd& Hessian, VectorXd& grad, double* p, double* t0,
	double* t1, double* t2, double d_hat_2, int* vertex_in_pair, double  stiffness, bool compute_hessian, int type)
{
	double distance=1.0;
	double b_grad = 0.0, b_hessian = 0.0;
	switch (CCD::internal::pointTriangleDistanceType(p, t0, t1, t2))
	{
	case 0:
		distance = CCD::internal::pointPointDistance(p, t0);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}

		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;
		grad.resize(6);
		distance::point_point_distance_gradient(p, t0, grad.data());

		if (compute_hessian) {
			Hessian.resize(6, 6);
			distance::point_point_distance_hessian(p, t0, Hessian.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_grad *= stiffness; b_hessian *= stiffness;
			Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

			switch (type)
			{
			case 0:
				FEM::SPDprojection(Hessian);
				break;
			case 1: {
				MatrixXd temp = Hessian.block<3, 3>(0, 0);
				FEM::SPDprojection(temp);
				Hessian.block<3, 3>(0, 0) = temp;
			}
				  break;
			case 2: {
				MatrixXd temp = Hessian.block<3, 3>(3, 3);
				FEM::SPDprojection(temp);
				Hessian.block<3, 3>(3, 3) = temp;
			}
				  break;
			}

			
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
			b_grad *= stiffness;
		}		
		grad *= b_grad;

		break;
	case 1:
		distance = CCD::internal::pointPointDistance(p, t1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
		grad.resize(6);
		distance::point_point_distance_gradient(p, t1, grad.data());

		if (compute_hessian) {
			Hessian.resize(6, 6);
			distance::point_point_distance_hessian(p, t1, Hessian.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_grad *= stiffness; b_hessian *= stiffness;
			Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

			switch (type)
			{
			case 0:
				FEM::SPDprojection(Hessian);
				break;
			case 1: {
				MatrixXd temp = Hessian.block<3, 3>(0, 0);
				FEM::SPDprojection(temp);
				Hessian.block<3, 3>(0, 0) = temp;
			}
				  break;
			case 2: {
				MatrixXd temp = Hessian.block<3, 3>(3, 3);
				FEM::SPDprojection(temp);
				Hessian.block<3, 3>(3, 3) = temp;
			}
				  break;
			}
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
			b_grad *= stiffness;
		}

		grad *= b_grad;
		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 2;


		break;
	case 2:
		distance = CCD::internal::pointPointDistance(p, t2);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
		grad.resize(6);
		distance::point_point_distance_gradient(p, t2, grad.data());

		if (compute_hessian) {
			Hessian.resize(6, 6);
			distance::point_point_distance_hessian(p, t2, Hessian.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_grad *= stiffness; b_hessian *= stiffness;
			Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();
			switch (type)
			{
			case 0:
				FEM::SPDprojection(Hessian);
				break;
			case 1: {
				MatrixXd temp = Hessian.block<3, 3>(0, 0);
				FEM::SPDprojection(temp);
				Hessian.block<3, 3>(0, 0) = temp;
			}
				  break;
			case 2: {
				MatrixXd temp = Hessian.block<3, 3>(3, 3);
				FEM::SPDprojection(temp);
				Hessian.block<3, 3>(3, 3) = temp;
			}
				  break;
			}
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
			b_grad *= stiffness;
		}
		grad *= b_grad;
		vertex_in_pair[0] = 2;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 3;

		break;

	case 3:
		distance = CCD::internal::pointEdgeDistance(p, t0, t1);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
		grad.resize(9);
		distance::point_edge_distance_gradient(p, t0, t1, grad.data());

		if (compute_hessian) {
			Hessian.resize(9, 9);
			distance::point_edge_distance_hessian(p, t0, t1, Hessian.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_grad *= stiffness; b_hessian *= stiffness;
			Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();
			switch (type)
			{
			case 0:
				FEM::SPDprojection(Hessian);
				break;
			case 1: {
				MatrixXd temp = Hessian.block<3, 3>(0, 0);
				FEM::SPDprojection(temp);
				Hessian.block<3, 3>(0, 0) = temp;
			}
				  break;
			case 2: {
				MatrixXd temp = Hessian.block<6, 6>(3, 3);
				FEM::SPDprojection(temp);
				Hessian.block<6, 6>(3, 3) = temp;
			}
				  break;
			}
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
			b_grad *= stiffness;
		}
		grad *= b_grad;

		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;
		vertex_in_pair[3] = 2;


		break;

	case 4:
		distance = CCD::internal::pointEdgeDistance(p, t1, t2);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
		grad.resize(9);
		distance::point_edge_distance_gradient(p, t1, t2, grad.data());

		if (compute_hessian) {
			Hessian.resize(9, 9);
			distance::point_edge_distance_hessian(p, t1, t2, Hessian.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_grad *= stiffness; b_hessian *= stiffness;
			Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();
			switch (type)
			{
			case 0:
				FEM::SPDprojection(Hessian);
				break;
			case 1: {
				MatrixXd temp = Hessian.block<3, 3>(0, 0);
				FEM::SPDprojection(temp);
				Hessian.block<3, 3>(0, 0) = temp;
			}
				  break;
			case 2: {
				MatrixXd temp = Hessian.block<6, 6>(3, 3);
				FEM::SPDprojection(temp);
				Hessian.block<6, 6>(3, 3) = temp;
			}
				  break;
			}
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
			b_grad *= stiffness;
		}
		grad *= b_grad;

		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 2;
		vertex_in_pair[3] = 3;


		break;
	case 5:
		distance = CCD::internal::pointEdgeDistance(p, t0, t2);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}

		grad.resize(9);
		distance::point_edge_distance_gradient(p, t0, t2, grad.data());

		if (compute_hessian) {
			Hessian.resize(9, 9);
			distance::point_edge_distance_hessian(p, t0, t2, Hessian.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_grad *= stiffness; b_hessian *= stiffness;
			Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();
			switch (type)
			{
			case 0:
				FEM::SPDprojection(Hessian);
				break;
			case 1: {
				MatrixXd temp = Hessian.block<3, 3>(0, 0);
				FEM::SPDprojection(temp);
				Hessian.block<3, 3>(0, 0) = temp;
			}
				  break;
			case 2: {
				MatrixXd temp = Hessian.block<6, 6>(3, 3);
				FEM::SPDprojection(temp);
				Hessian.block<6, 6>(3, 3) = temp;
			}
				  break;
			}
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
			b_grad *= stiffness;
		}
		grad *= b_grad;
		vertex_in_pair[0] = 3;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;
		vertex_in_pair[3] = 3;

		break;

	case 6:
		distance = CCD::internal::pointTriangleDistance(p, t0, t1, t2);
		if (distance >= d_hat_2) {
			vertex_in_pair[0] = 0;
			return false;
		}
	
		grad.resize(12);
		distance::point_triangle_distance_gradient(p, t0, t1, t2, grad.data());

		if (compute_hessian) {
			Hessian.resize(12, 12);
			distance::point_triangle_distance_hessian(p, t0, t1, t2, Hessian.data());
			barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
			b_grad *= stiffness; b_hessian *= stiffness;
			Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();
			switch (type)
			{
			case 0:
				FEM::SPDprojection(Hessian);
				break;
			case 1: {
				MatrixXd temp = Hessian.block<3, 3>(0, 0);
				FEM::SPDprojection(temp);
				Hessian.block<3, 3>(0, 0) = temp;
			}
				  break;
			case 2: {
				MatrixXd temp = Hessian.block<9, 9>(3, 3);
				FEM::SPDprojection(temp);
				Hessian.block<9, 9>(3, 3) = temp;
			}
				  break;
			}
		}
		else {
			barrierGrad(distance, d_hat_2, b_grad);
			b_grad *= stiffness;
		}
		grad *= b_grad;

		vertex_in_pair[0] = 4;
		vertex_in_pair[1] = 0;
		vertex_in_pair[2] = 1;
		vertex_in_pair[3] = 2;
		vertex_in_pair[4] = 3;


		break;
	}


	//if (distance == 0.0) {
	//	std::cout << "error distance is zero " << std::endl;
	//}

	return true;

}




void SecondOrderConstraint::computeVTBarrierGradientHessianTest(MatrixXd& Hessian_, VectorXd& grad_, double* p, double* t0,
	double* t1, double* t2, double d_hat_2, int* triangle_vertex_order_in_system, double stiffness, int* record_vertex_in_pair, int record_vertex_in_use)
{
	int compute_num=0;

	double distance;
	double b_grad, b_hessian;
	MatrixXd Hessian; VectorXd grad;
	int vertex_in_pair[4];
	memset(vertex_in_pair, 0, 16);
	switch (CCD::internal::pointTriangleDistanceType(p, t0, t1, t2))
	{
	case 0:
		distance = CCD::internal::pointPointDistance(p, t0);
		if (distance >= d_hat_2) {
			return;
		}
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;

		Hessian.resize(6, 6);
		grad.resize(6);


		distance::point_point_distance_gradient(p, t0, grad.data());
		distance::point_point_distance_hessian(p, t0, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;

		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 2);
		compute_num = 1;
		break;
	case 1:
		distance = CCD::internal::pointPointDistance(p, t1);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(6, 6);
		grad.resize(6);

		distance::point_point_distance_gradient(p, t1, grad.data());
		distance::point_point_distance_hessian(p, t1, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 2);
		compute_num = 1;
		break;
	case 2:
		distance = CCD::internal::pointPointDistance(p, t2);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(6, 6);
		grad.resize(6);
		distance::point_point_distance_gradient(p, t2, grad.data());
		distance::point_point_distance_hessian(p, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 2);
		compute_num = 1;
		break;

	case 3:
		distance = CCD::internal::pointEdgeDistance(p, t0, t1);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(9, 9);
		grad.resize(9);
		distance::point_edge_distance_gradient(p, t0, t1, grad.data());
		distance::point_edge_distance_hessian(p, t0, t1, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 3);
		compute_num = 2;
		break;

	case 4:
		distance = CCD::internal::pointEdgeDistance(p, t1, t2);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(9, 9);
		grad.resize(9);
		distance::point_edge_distance_gradient(p, t1, t2, grad.data());
		distance::point_edge_distance_hessian(p, t1, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);


		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;
		vertex_in_pair[2] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 3);
		compute_num = 2;
		break;
	case 5:
		distance = CCD::internal::pointEdgeDistance(p, t0, t2);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(9, 9);
		grad.resize(9);
		distance::point_edge_distance_gradient(p, t0, t2, grad.data());
		distance::point_edge_distance_hessian(p, t0, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 3);
		compute_num = 2;
		break;

	case 6:
		distance = CCD::internal::pointTriangleDistance(p, t0, t1, t2);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(12, 12);
		grad.resize(12);
		distance::point_triangle_distance_gradient(p, t0, t1, t2, grad.data());
		distance::point_triangle_distance_hessian(p, t0, t1, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;
		vertex_in_pair[3] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 4);
		compute_num = 3;
		break;
	}


	if (record_vertex_in_use != compute_num+1) {
		std::cout << "record_vertex_in_use wrong, right: "<< compute_num+1<<" "<< record_vertex_in_use << std::endl;
	}
	for (int i = 0; i < record_vertex_in_use; ++i) {
		if (vertex_in_pair[i] != record_vertex_in_pair[i]) {
			std::cout << "record_vertex_in_use wrong , right " << vertex_in_pair[i] << " " << record_vertex_in_pair[i] << std::endl;
		}
	}

}






void SecondOrderConstraint::computeVTBarrierGradientHessian(MatrixXd& Hessian_, VectorXd& grad_, double* p, double* t0, 
	double* t1, double* t2,  double d_hat_2, int* triangle_vertex_order_in_system, double stiffness)
{
	double distance;
	double b_grad, b_hessian;
	MatrixXd Hessian; VectorXd grad;
	int vertex_in_pair[4];
	memset(vertex_in_pair, 0, 16);
	switch (CCD::internal::pointTriangleDistanceType(p, t0, t1, t2))
	{
	case 0:
		distance = CCD::internal::pointPointDistance(p, t0);
		if (distance >= d_hat_2) {
			return;
		}
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;

		Hessian.resize(6, 6);
		grad.resize(6);

		
		distance::point_point_distance_gradient(p, t0, grad.data());
		distance::point_point_distance_hessian(p, t0, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);
		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;

		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 2);

		break;
	case 1:
		distance = CCD::internal::pointPointDistance(p, t1);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(6, 6);
		grad.resize(6);

		distance::point_point_distance_gradient(p, t1, grad.data());
		distance::point_point_distance_hessian(p, t1, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 2);

		break;
	case 2:
		distance = CCD::internal::pointPointDistance(p, t2);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(6, 6);
		grad.resize(6);
		distance::point_point_distance_gradient(p, t2, grad.data());
		distance::point_point_distance_hessian(p, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 2);

		break;

	case 3:
		distance = CCD::internal::pointEdgeDistance(p, t0, t1);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(9, 9);
		grad.resize(9);
		distance::point_edge_distance_gradient(p, t0, t1, grad.data());
		distance::point_edge_distance_hessian(p, t0, t1, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 3);

		break;

	case 4:
		distance = CCD::internal::pointEdgeDistance(p, t1, t2);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(9, 9);
		grad.resize(9);
		distance::point_edge_distance_gradient(p, t1, t2, grad.data());
		distance::point_edge_distance_hessian(p, t1, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);


		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 2;
		vertex_in_pair[2] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 3);

		break;
	case 5:
		distance = CCD::internal::pointEdgeDistance(p, t0, t2);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(9, 9);
		grad.resize(9);
		distance::point_edge_distance_gradient(p, t0, t2, grad.data());
		distance::point_edge_distance_hessian(p, t0, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 3);
		break;

	case 6:
		distance = CCD::internal::pointTriangleDistance(p, t0, t1, t2);
		if (distance >= d_hat_2) {
			return;
		}
		Hessian.resize(12, 12);
		grad.resize(12);
		distance::point_triangle_distance_gradient(p, t0, t1, t2, grad.data());
		distance::point_triangle_distance_hessian(p, t0, t1, t2, Hessian.data());
		barrierGradHessian(distance, d_hat_2, b_grad, b_hessian);

		b_grad *= stiffness; b_hessian *= stiffness;
		Hessian = (b_hessian * grad * grad.transpose() + b_grad * Hessian).eval();

		FEM::SPDprojection(Hessian);

		grad *= b_grad;
		vertex_in_pair[0] = 0;
		vertex_in_pair[1] = 1;
		vertex_in_pair[2] = 2;
		vertex_in_pair[3] = 3;
		setTetHessianFromBarrierHessian(Hessian_, grad_, Hessian, grad,
			triangle_vertex_order_in_system, vertex_in_pair, 4);

		break;
	} 

	//std::cout << "-=-=-==-" << std::endl;
	//std::cout << Hessian << std::endl;

}


void SecondOrderConstraint::setFourVertexGradFromBarrier(VectorXd& grad_system, VectorXd& grad_,
	int* vertex_in_pair, int vertex_in_use)
{
	for (int i = 0; i < vertex_in_use; ++i) {
		grad_system.segment(3 * vertex_in_pair[i], 3) += grad_.segment(3 * i, 3);
	}
}


void SecondOrderConstraint::setFourVertexHessianFromBarrierHessian(MatrixXd& Hessian_system, VectorXd& grad_system, MatrixXd& Hessian_, VectorXd& grad_,
	int* vertex_in_pair, int vertex_in_use)
{
	for (int i = 0; i < vertex_in_use; ++i) {
		for (int j = 0; j < vertex_in_use; ++j) {
			Hessian_system.block<3, 3>(3 * vertex_in_pair[i], 3 * vertex_in_pair[j])
				+= Hessian_.block<3, 3>(3 * i, 3 * j);
		}
		grad_system.segment(3 * vertex_in_pair[i], 3) += grad_.segment(3 * i, 3);
	}
}



void SecondOrderConstraint::setTetHessianFromHessian(MatrixXd& Hessian_system, VectorXd& grad_system, MatrixXd& Hessian_, VectorXd& grad_,
	int* triangle_vertex_order_in_system)
{
	for (int i = 0; i < 4; ++i) {
		if (triangle_vertex_order_in_system[i] != -1) {
			for (int j = 0; j < 4; ++j) {
				if (triangle_vertex_order_in_system[j] != -1) {

					Hessian_system.block<3, 3>(3 * triangle_vertex_order_in_system[i], 3 * triangle_vertex_order_in_system[j])
						+= Hessian_.block<3, 3>(3 * i, 3 * j);
				}

			}
			grad_system.segment(3 * triangle_vertex_order_in_system[i], 3) += grad_.segment(3 * i, 3);
		}
	}

}


void SecondOrderConstraint::setTetHessianFromBarrierHessian(MatrixXd& Hessian_system, VectorXd& grad_system, MatrixXd& Hessian_, VectorXd& grad_, 
	int* triangle_vertex_order_in_system, int* vertex_in_pair, int vertex_in_use)
{
	for (int i = 0; i < vertex_in_use; ++i) {
		if (triangle_vertex_order_in_system[vertex_in_pair[i]] != -1) {
			for (int j = 0; j < vertex_in_use; ++j) {
				if (triangle_vertex_order_in_system[vertex_in_pair[j]] != -1) {

					Hessian_system.block<3, 3>(3 * triangle_vertex_order_in_system[vertex_in_pair[i]], 3 * triangle_vertex_order_in_system[vertex_in_pair[j]])
						+= Hessian_.block<3, 3>(3 * i, 3 * j);
				}

			}
			grad_system.segment(3 * triangle_vertex_order_in_system[vertex_in_pair[i]], 3) += grad_.segment(3 * i, 3);
		}
	}
}



bool SecondOrderConstraint::getCollisionPairHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	double ori_volume, Matrix3d& Hessian, Vector3d& grad, unsigned int vertex_no)
{
	double C;

	double volume =getTetCubeVolume(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3);
	if (abs(volume) >= ori_volume) {
		return false;
	}
	if (abs(volume) < 1e-9) {
		std::cout << "may be edge edge parallel " << std::endl;
		return false;
	}
	double e0[3], e1[3];
	switch (vertex_no)
	{
	case 0: {
		SUB(e0, vertex_position_2, vertex_position_1);
		CROSS(grad.data(), vertex_position_3,  e0);
		double temp[3];
		CROSS(temp, vertex_position_2, vertex_position_1);
		SUM_(grad, temp);
	}
		break;
	case 1:
		SUB(e0, vertex_position_2, vertex_position_0);
		SUB(e1, vertex_position_3, vertex_position_0);
		CROSS(grad.data(), e0, e1);
		break;
	case 2:
		SUB(e0, vertex_position_3, vertex_position_0);
		SUB(e1, vertex_position_1, vertex_position_0);
		CROSS(grad.data(), e0, e1);
		break;	
	case 3:
		SUB(e0, vertex_position_1, vertex_position_0);
		SUB(e1, vertex_position_2, vertex_position_0);
		CROSS(grad.data(), e0, e1);
		break;
	}
	if (volume < 0) {
		grad *= -1.0;
		volume = -volume;
		//std::cout << "collision volume negative" << std::endl;
	}
	double ln_ = log(volume / ori_volume);
	Hessian = ( - 2.0 * ln_ + (ori_volume - volume) * (ori_volume + 3.0 * volume) / (volume * volume))*grad*grad.transpose();
	grad *= (ori_volume - volume) * (2.0 * ln_ - ori_volume / volume + 1.0);
	
	return true;
}




bool SecondOrderConstraint::getARAPGradHessianNewton(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	Matrix<double, 3, 4>& A, Matrix3d& Hessian, Vector3d& grad, double& C, unsigned int vertex_no)
{
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d S, rotation;
	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);
	FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);
	C = (deformation_gradient - rotation).squaredNorm();
	if (C < 1e-16) {
		C = 0;
		return false;
	}

	grad =2.0 * (deformation_gradient - rotation) * A.col(vertex_no);
	Matrix3d Dm = A.block<3, 3>(0, 1).transpose();

	//Matrix<double, 12, 12> Hessian_global;
	//FEM::getHessian(Hessian_global, S, rotation, Dm, A);
	//Hessian = Hessian_global.block<3, 3>(3 * vertex_no, 3 * vertex_no);

	FEM::getHessianForOneVertex(Hessian, S, rotation, Dm, A, vertex_no);
	return true;
}


bool SecondOrderConstraint::getARAPGradHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	Matrix<double, 3, 4>& A, Matrix3d& Hessian, Vector3d& grad, double& C, unsigned int vertex_no)
{
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d S, rotation;
	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);
	FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);
	C = (deformation_gradient - rotation).norm();
	if (C < 1e-8) {
		C = 0;
		return false;
	}
	grad = (1.0 / C) * (deformation_gradient - rotation) * A.col(vertex_no);

	Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
	FEM::getHessianForOneVertex(Hessian, S, rotation, Dm, A,vertex_no);
	Hessian *= (0.5 / C);
	Hessian -= ((1.0 / C) * grad) * grad.transpose();
	return true;
}


void SecondOrderConstraint::solveARAPConstraint(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	double stiffness, double dt,
	Matrix<double, 3, 4>& A, double* inv_mass, double& lambda, const double damping_stiffness, double* mass,
	double volume)
{
	Vector3d eigen_value;
	Vector3d position;
	Matrix3d deformation_gradient;
	Matrix3d S, rotation;

	FEM::getDeformationGradient(vertex_position_0, vertex_position_1, vertex_position_2, vertex_position_3, A, deformation_gradient);

	FEM::polarDecomposition(deformation_gradient, eigen_value, S, rotation);

	//use P_inv to record transform

	double C = (deformation_gradient - rotation).norm();

	if (C < 1e-8) {
		return;
	}

	double alpha_ = 1.0 / (stiffness * volume * dt * dt);

	Matrix<double, 3, 4> grad_C_transpose;
	grad_C_transpose = (1.0 / C) * (deformation_gradient - rotation) * A;	
	Matrix<double, 12, 1> grad;
	memcpy(grad.data(), grad_C_transpose.data(), 96);

	Matrix<double, 12, 12> Hessian;

	Matrix3d Dm = A.block<3, 3>(0, 1).transpose();
	FEM::getHessian(Hessian, S, rotation, Dm, A);

	Hessian *= (0.5 / C);
	Hessian -= ((1.0 / C) * grad) * grad.transpose();

	for (unsigned int i = 0; i < 4; ++i) {
		if (inv_mass[i] == 0.0) {
			Hessian.block<12, 3>(0, 3 * i).setZero();
			Hessian.block<3, 12>(3 * i, 0).setZero();
			grad.segment(3 * i, 3).setZero();
		}
	}


	//Matrix<double, 12, 1> inv_mass_sys;
	//inv_mass_sys << inv_mass[0], inv_mass[0], inv_mass[0], inv_mass[1], inv_mass[1], inv_mass[1], inv_mass[2], inv_mass[2], inv_mass[2],
	//	inv_mass[3], inv_mass[3], inv_mass[3];
	////how to improve this
	//for (unsigned int j = 0; j < 12; ++j) {
	//	for (unsigned int i = 0; i < 4; ++i) {
	//		Hessian(i + i + i, j) *= inv_mass[i];
	//		Hessian(i + i + i+1, j) *= inv_mass[i];
	//		Hessian(i + i + i+2, j) *= inv_mass[i];
	//	}
	//}

	Hessian *= -lambda;

	for (unsigned int i = 0; i < 4; ++i) {
		Hessian.data()[39 * i] += mass[i];
		Hessian.data()[39 * i+13] += mass[i];
		Hessian.data()[39 * i+26] += mass[i];
	}

	//Hessian.setZero();
	//for (unsigned int i = 0; i < 144; i += 13) {
	//	Hessian.data()[i] += 1.0;
	//}

	//Matrix<double, 12, 1> M_inv_g;
	//SUB(M_inv_g.data(), vertex_position_0, sn_0);
	//M_inv_g.data()[3] = vertex_position_1[0] - sn_1[0];
	//M_inv_g.data()[4] = vertex_position_1[1] - sn_1[1];
	//M_inv_g.data()[5] = vertex_position_1[2] - sn_1[2];
	//M_inv_g.data()[6] = vertex_position_2[0] - sn_2[0];
	//M_inv_g.data()[7] = vertex_position_2[1] - sn_2[1];
	//M_inv_g.data()[8] = vertex_position_2[2] - sn_2[2];
	//M_inv_g.data()[9] = vertex_position_3[0] - sn_3[0];
	//M_inv_g.data()[10] = vertex_position_3[1] - sn_3[1];
	//M_inv_g.data()[11] = vertex_position_3[2] - sn_3[2];
	//M_inv_g -= lambda * inv_mass_sys.cwiseProduct(grad);

	ColPivHouseholderQR <Matrix<double,12,12>> linear(Hessian);

	//if (abs(Hessian.determinant()) < 1e-8) {
	//	std::cout << Hessian.determinant() << std::endl;
	//}

	double h = C + alpha_ * lambda;

	double delta_lambda = (-h) / (grad.dot(linear.solve(grad)) + alpha_);//  + grad.dot(linear.solve(M_inv_g))
	Matrix<double, 12, 1> delta_x = linear.solve(delta_lambda * grad); //-M_inv_g

	//double delta_lambda = (-h) / (grad.dot(linear.solve(inv_mass_sys.cwiseProduct(grad))) + alpha_);//  + grad.dot(linear.solve(M_inv_g))
	//Matrix<double, 12, 1> delta_x = linear.solve(delta_lambda * (inv_mass_sys.cwiseProduct(grad))); //-M_inv_g


	vertex_position_0[0] += delta_x[0];
	vertex_position_0[1] += delta_x[1];
	vertex_position_0[2] += delta_x[2];

	vertex_position_1[0] += delta_x[3];
	vertex_position_1[1] += delta_x[4];
	vertex_position_1[2] += delta_x[5];
	
	vertex_position_2[0] += delta_x[6];
	vertex_position_2[1] += delta_x[7];
	vertex_position_2[2] += delta_x[8];

	vertex_position_3[0] += delta_x[9];
	vertex_position_3[1] += delta_x[10];
	vertex_position_3[2] += delta_x[11];
	 
	lambda += delta_lambda;

}





void SecondOrderConstraint::solveEdgeLengthConstraint(Vector3f& p1, Vector3f& p2, const double d, double mass_0,
	double mass_1, Vector3f& ori_p1, Vector3f& ori_p2, bool v0_fixed, bool v1_fixed, double& lambda)
{
	float w1 =1.0/mass_0;
	float w2 = 1.0/mass_1;

	if (v0_fixed) {
		w1 = 0;
	}
	if (v1_fixed) {
		w2 = 0;
	}
	Matrix<float, 6,1> g;
	g.setZero();
	if (w1 > 0)
		g.segment(0, 3) = (p1 - ori_p1) / w1;
	if (w2 > 0)
		g.segment(3, 3) = (p2 - ori_p2) / w2;

	//std::cout <<"==="<< g.transpose() << std::endl;

	Matrix<float, 6, 1> M_1;
	M_1 << w1, w1, w1, w2, w2, w2;
	Eigen::Vector3f n = p1 - p2;
	float norm = n.norm();
	n /= norm;
	Matrix<float, 6,1> grad;
	grad << n, -n;
	Eigen::Matrix<float, 6, 6> H;
	Eigen::Matrix3f s = Eigen::Matrix3f::Identity() - n * n.transpose();
	H.topLeftCorner(3, 3) = s;
	H.topRightCorner(3, 3) = -s;
	H.bottomLeftCorner(3, 3) = -s;
	H.bottomRightCorner(3, 3) = s;
	H /= norm;

	Eigen::Matrix<float, 6, 6> M_lambda_H_1 = (Eigen::Matrix<float, 6, 6>::Identity() + lambda * M_1.asDiagonal() * H).inverse() * M_1.asDiagonal();
	Eigen::Matrix<float, 6,1> tmp = M_lambda_H_1 * grad;
	g += lambda * grad;
	// Eigen::Vector<float, 6> tmp = M_1 * grad;
	float d_lambda = (norm - d - g.dot(tmp)) / (grad.dot(tmp));
	lambda += d_lambda;
	// std::cout << lambda << std::endl;
	Eigen::Matrix<float, 6,1> d_x = d_lambda * tmp + M_lambda_H_1 * g;
	p1 -= d_x.head(3);
	p2 -= d_x.tail(3);
}

void SecondOrderConstraint::solveEdgeLengthConstraint(double* p0, double* p1, const double rest_length, double stiffness, double mass_0,
	double mass_1, double time_step, double* sn_0, double* sn_1, bool v0_fixed, bool v1_fixed, unsigned int edge_index, double& lambda, unsigned int vertex_0_index, unsigned int vertex_1_index)
{
	if (v0_fixed && v1_fixed) {
		return;
	}

	//std::cout << "true " << std::endl;
	//std::cout << edge_index << " " << mass_0 << " " << mass_1 << " " << rest_length<<" "<< sn_0[0] << " " << sn_0[1] << " " << sn_0[2] << std::endl;


	//Vector3f p_1, p_2,ori_1,ori_2;
	//for (unsigned int i = 0; i < 3; ++i) {
	//	p_1[i] = p0[i];
	//	p_2[i] = p1[i];
	//	ori_1[i] = sn_0[i];
	//	ori_2[i] = sn_1[i];
	//}

	//if (edge_index == 3) {
	//	std::cout << lambda << std::endl;

	//}

	//solveEdgeLengthConstraint(p_1, p_2, rest_length, mass_0, mass_1, ori_1, ori_2, v0_fixed, v1_fixed, lambda);

	//if (edge_index == 3) {
	//	std::cout << lambda << std::endl;
	//	std::cout << sn_0[0] << " " << sn_0[1] << " " << sn_0[2] << std::endl;
	//	std::cout << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
	//}


	//p0[0] = p_1[0];
	//p0[1] = p_1[1];
	//p0[2] = p_1[2];
	//p1[0] = p_2[0];
	//p1[1] = p_2[1];
	//p1[2] = p_2[2];
	Vector3d n;
	SUB(n.data(), p0, p1);
	double n_norm = sqrt(DOT(n, n));
	if (abs(n_norm - rest_length) < 1e-8) {
		return;
	}
	Vector3d grad_ = n / n_norm;
	double alpha =  1.0 / (time_step * time_step * stiffness);
	Matrix3d He = (1.0 / n_norm) * (Matrix3d::Identity() - grad_ * grad_.transpose());
	//Matrix<double, 6, 6>Hessian;
	//Hessian.block<3, 3>(0, 0) = He;
	//Hessian.block<3, 3>(3, 3) = He;
	//Hessian.block<3, 3>(3, 0) = -He;
	//Hessian.block<3, 3>(0, 3) = -He;
	VectorXd mass_inv(6);
	mass_inv[0] = mass_inv[1] = mass_inv[2] = 1.0 / mass_0;
	mass_inv[3] = mass_inv[4] = mass_inv[5] = 1.0 / mass_1;
	if (v0_fixed) {
		mass_inv[0] = mass_inv[1] = mass_inv[2] = 0;
	}
	if (v1_fixed) {
		mass_inv[3] = mass_inv[4] = mass_inv[5] = 0;
	}
	//Matrix<double, 6, 1> gradient;
	//gradient.segment(0, 3) = grad_;
	//gradient.segment(3, 3) = -grad_;
	/*
	Matrix<double, 6, 6>sys_matrix = Matrix<double, 6, 6>::Identity() + lambda * mass_inv.asDiagonal() * Hessian;

	//ColPivHouseholderQR <Matrix<double, 6, 6>> linear(sys_matrix);
	Matrix<double, 6, 6>  linear = sys_matrix.inverse() * mass_inv.asDiagonal();
	Matrix<double, 6, 1>g;
	g.data()[0] = mass_0 * (p0[0] - sn_0[0]);
	g.data()[1] = mass_0 * (p0[1] - sn_0[1]);
	g.data()[2] = mass_0 * (p0[2] - sn_0[2]);
	g.data()[3] = mass_1 * (p1[0] - sn_1[0]);
	g.data()[4] = mass_1 * (p1[1] - sn_1[1]);
	g.data()[5] = mass_1 * (p1[2] - sn_1[2]);
	g += lambda * gradient;

	Matrix<double, 6,1> tmp = linear * gradient;

	//double coe = linear.solve(mass_inv.asDiagonal()*gradient).dot(gradient)+alpha;
	double coe = tmp.dot(gradient)+alpha;
	double delta_lambda = (n_norm - rest_length - alpha * lambda) - (g.dot(tmp));
	delta_lambda /= coe;
	Matrix<double, 6, 1>delta_x;
	delta_x =-1.0*(delta_lambda*tmp + sys_matrix*g);
	//delta_x =-1.0*( delta_lambda *linear *  (mass_inv.asDiagonal() * gradient) + linear *  M_inv_g);
	p0[0] += delta_x[0];
	p0[1] += delta_x[1];
	p0[2] += delta_x[2];

	p1[0] += delta_x[3];
	p1[1] += delta_x[4];
	p1[2] += delta_x[5];		
	lambda += delta_lambda;
	*/
	//Matrix3d He = (((- lambda) /n_norm) * grad_) * grad_.transpose();
	//double coe = lambda / n_norm;
	//He.data()[0] += coe;
	//He.data()[4] += coe;
	//He.data()[8] += coe;
	double coe;
		if (v0_fixed) {
			/*
			Matrix<double, 3, 3> sys_matrix;
			sys_matrix = mass_1/(time_step*time_step) * Matrix3d::Identity() + stiffness * grad_*grad_.transpose() + (stiffness*(n_norm - rest_length)) * He;
			Vector3d g;
			//g.data()[0] = mass_1 * (p1[0] - sn_1[0]);
			//g.data()[1] = mass_1 * (p1[1] - sn_1[1]);
			//g.data()[2] = mass_1 * (p1[2] - sn_1[2]);
			g.setZero();
			g -= (stiffness*(n_norm - rest_length))*grad_;
			ColPivHouseholderQR <Matrix3d> linear(sys_matrix);
			Vector3d delta_x = -linear.solve(g);
			p1[0] += delta_x[0];
			p1[1] += delta_x[1];
			p1[2] += delta_x[2];
			*/
			
			Matrix<double, 3, 3> sys_matrix;
			sys_matrix = mass_1 * Matrix3d::Identity()-lambda * He;
			//sys_matrix.setZero();
			//sys_matrix -= He;
			//sys_matrix.data()[0] += mass_1;
			//sys_matrix.data()[4] += mass_1;
			//sys_matrix.data()[8] += mass_1;
			Vector3d g;
			g.setZero();
			//g.data()[0] = mass_1 * (p1[0] - sn_1[0]);
			//g.data()[1] = mass_1 * (p1[1] - sn_1[1]);
			//g.data()[2] = mass_1 * (p1[2] - sn_1[2]);
			g += lambda * grad_;
			double h = n_norm - rest_length + alpha * lambda;
			ColPivHouseholderQR <Matrix3d> linear(sys_matrix);
			//LDLT <Matrix3d> linear(sys_matrix);
			coe = (linear.solve(grad_)).dot(grad_) + alpha;
			double delta_lambda =( - h - (linear.solve(g)).dot(grad_))/coe;//
			Vector3d delta_x = linear.solve((- delta_lambda) * grad_-g);//-g
			p1[0] += delta_x[0];
			p1[1] += delta_x[1];
			p1[2] += delta_x[2];
			lambda += delta_lambda;
			
			return;
		}		
		if (v1_fixed) {
			/*
			Matrix<double, 3, 3> sys_matrix;
			sys_matrix = mass_0 / (time_step * time_step) * Matrix3d::Identity() + stiffness * grad_ * grad_.transpose() + (stiffness * (n_norm - rest_length)) * He;
			Vector3d g;
			//g.data()[0] = mass_0 * (p0[0] - sn_0[0]);
			//g.data()[1] = mass_0 * (p0[1] - sn_0[1]);
			//g.data()[2] = mass_0 * (p0[2] - sn_0[2]);
			g.setZero();
			g += (stiffness * (n_norm - rest_length)) * grad_;
			ColPivHouseholderQR <Matrix3d> linear(sys_matrix);
			Vector3d delta_x = -linear.solve(g);
			p0[0] += delta_x[0];
			p0[1] += delta_x[1];
			p0[2] += delta_x[2];
			*/
			
			Matrix<double, 3, 3> sys_matrix;
			sys_matrix = mass_0 * Matrix3d::Identity() -lambda * He;
			//sys_matrix.setZero();
			//sys_matrix -= He;
			//sys_matrix.data()[0] += mass_0;
			//sys_matrix.data()[4] += mass_0;
			//sys_matrix.data()[8] += mass_0;
			Vector3d g;
			g.setZero();
			//g.data()[0] = mass_0 * (p0[0] - sn_0[0]);
			//g.data()[1] = mass_0 * (p0[1] - sn_0[1]);
			//g.data()[2] = mass_0 * (p0[2] - sn_0[2]);
			g -= lambda * grad_;

			double h = n_norm - rest_length + alpha * lambda;
			//LDLT <Matrix3d> linear(sys_matrix);
			ColPivHouseholderQR <Matrix3d> linear(sys_matrix);
			coe = linear.solve(grad_).dot(grad_) + alpha;
			double delta_lambda = (-h + linear.solve(g).dot(grad_)) / coe;//
			Vector3d delta_x = linear.solve(delta_lambda * grad_-g);//- g
			p0[0] += delta_x[0];
			p0[1] += delta_x[1];
			p0[2] += delta_x[2];
			lambda += delta_lambda;
			
			return;
		}

		/*
		///test
		Matrix<double, 6, 1> gradient;
		gradient.segment(0, 3) = grad_;
		gradient.segment(3, 3) = -grad_;
		Matrix<double, 6, 6> sys_matrix;
		sys_matrix.setZero();
		sys_matrix.block<3, 3>(0, 0) = -He;
		sys_matrix.block<3, 3>(3, 3) = -He;
		sys_matrix.block<3, 3>(3, 0) = He;
		sys_matrix.block<3, 3>(0, 3) = He;
		sys_matrix *= (n_norm - rest_length) * stiffness;
		sys_matrix += stiffness * gradient * gradient.transpose();
		for (unsigned int i = 0; i < 18; i += 7) {
			sys_matrix.data()[i] += mass_0 / (time_step * time_step);
		}
		for (unsigned int i = 21; i < 36; i += 7) {
			sys_matrix.data()[i] += mass_1 / (time_step * time_step);
		}
		Matrix<double, 6, 1>g;
		g.setZero();
		//g.data()[0] = mass_0 / (time_step * time_step) * (p0[0] - sn_0[0]);
		//g.data()[1] = mass_0 / (time_step * time_step) * (p0[1] - sn_0[1]);
		//g.data()[2] = mass_0 / (time_step * time_step) * (p0[2] - sn_0[2]);
		//g.data()[3] = mass_1 / (time_step * time_step) * (p1[0] - sn_1[0]);
		//g.data()[4] = mass_1 / (time_step * time_step) * (p1[1] - sn_1[1]);
		//g.data()[5] = mass_1 / (time_step * time_step) * (p1[2] - sn_1[2]);
		g += (n_norm - rest_length) * stiffness * gradient;
		Matrix<double, 6, 1> result;
		ColPivHouseholderQR <Matrix<double, 6, 6>> linear_(sys_matrix);
		result = -1.0 * linear_.solve(g);
		p0[0] += result[0];
		p0[1] += result[1];
		p0[2] += result[2];
		p1[0] += result[3];
		p1[1] += result[4];
		p1[2] += result[5];

		////
		*/
		
		Matrix<double, 6, 1> gradient;
		gradient.segment(0, 3) = grad_;
		gradient.segment(3, 3) = -grad_;
		Matrix<double, 6, 6> sys_matrix;
		sys_matrix.setZero();
		sys_matrix.block<3, 3>(0, 0) =- He;
		sys_matrix.block<3, 3>(3, 3) =- He;
		sys_matrix.block<3, 3>(3, 0) = He;
		sys_matrix.block<3, 3>(0, 3) = He;
		sys_matrix *= lambda;
		for (unsigned int i = 0; i < 18; i+=7) {
			sys_matrix.data()[i] += mass_0;
		}
		for (unsigned int i = 21; i < 36; i+=7) {
			sys_matrix.data()[i] += mass_1;
		}
		//double det = sys_matrix.determinant();
		//if (det < 1e-8) {
		//	std::cout << det << std::endl;
		//}
		//ColPivHouseholderQR <Matrix<double, 6, 6>> linear(sys_matrix);
		//LDLT <Matrix<double, 6, 6>> linear(sys_matrix);
		double h = n_norm - rest_length + alpha*lambda;
		Matrix<double, 6, 1>g;
		g.setZero();
		//g.data()[0] =mass_0 *(p0[0] - sn_0[0]);
		//g.data()[1] =mass_0 *(p0[1] - sn_0[1]);
		//g.data()[2] =mass_0 *(p0[2] - sn_0[2]);
		//g.data()[3] = mass_1 * (p1[0] - sn_1[0]);
		//g.data()[4] = mass_1 * (p1[1] - sn_1[1]);
		//g.data()[5] = mass_1 * (p1[2] - sn_1[2]);
		g -= lambda*gradient;
	//	coe = linear.solve(gradient).dot(gradient) + alpha;
		double delta_lambda;
	//	delta_lambda = (-h + linear.solve(g).dot(gradient)) / coe;//+
		//if (edge_index == 0) {
		////	std::cout << g << std::endl;
		//std::cout << linear.solve(g) << std::endl;
		//}
		Matrix<double, 6, 1> delta_x; 
		// delta_x = linear.solve(delta_lambda * gradient - g);//		-g
		//if (vertex_0_index == 2451 || vertex_1_index == 2451) {
		//	std::cout << delta_lambda<<" "<<edge_index << std::endl;
		//	std::cout << delta_x.segment(0, 3).norm() << " " << delta_x.segment(3, 3).norm() << std::endl;
		//}
		Matrix<double, 7, 7> check;
		check.block<6, 6>(0, 0) = sys_matrix;
		check.block<6, 1>(0, 6) = -gradient;
		check.block<1, 6>(6, 0) = gradient.transpose();
		check.data()[7 * 7 - 1] = alpha;
		Matrix<double, 7, 1> vector;
		//vector.setZero();
		memcpy(vector.data(), g.data(), 48);
		vector.data()[6] = h;		
		Matrix<double, 7, 1> result;
		ColPivHouseholderQR <Matrix<double, 7, 7>> linear_(check);
	   result =-1.0* linear_.solve(vector);

	   //double k = check.determinant();
	  // if (k < 1e-7) {
		//   std::cout << edge_index << " " << k << std::endl;
	   //}

		//memcpy(result.data(), delta_x.data(), 48);
		//result.data()[6] = delta_lambda;
		//double error = (check * result + vector).norm();
		//if (error > 1e-10) {
		//	std::cout << "error " << error << std::endl;
		//	
			//std::cout << check << std::endl;
			//std::cout << lambda << std::endl;
			//std::cout << check.determinant() << std::endl;;
			//std::cout << g.transpose() << std::endl;;
			//std::cout << vector.transpose() << std::endl;;
			//std::cout <<  std::endl;
			//std::cout << result.transpose() << std::endl;;
		//}
		//}
		//else {
		//	//std::cout << error << std::endl;
		//}

	  


		p0[0] += result[0];
		p0[1] += result[1];
		p0[2] += result[2];
		p1[0] += result[3];
		p1[1] += result[4];
		p1[2] += result[5];
		lambda += result[6];
	
		//p0[0] += delta_x[0];
		//p0[1] += delta_x[1];
		//p0[2] += delta_x[2];
		//p1[0] += delta_x[3];
		//p1[1] += delta_x[4];
		//p1[2] += delta_x[5];		
		//lambda += delta_lambda;


 


	//Vector3d p0_current;
	//Vector3d p1_current;

	//SUM(p0_current, p0, delta_x);
	//p1_current[0]=p1[0] + delta_x[3];
	//p1_current[1]=p1[1] + delta_x[4];
	//p1_current[2]=p1[2] + delta_x[5];

	//double constraint =(p0_current - p1_current).norm() - rest_length;
	//constraint = 0.5 * constraint * constraint;
	//n_norm = sqrt(EDGE_LENGTH(p0, p1));
	//double constraint_1  = 1.0 / (2.0 * time_step * time_step) * (mass_0 * EDGE_LENGTH(p0, sn_0) + mass_1 * EDGE_LENGTH(p1, sn_1)) + 0.5 * stiffness * (n_norm - rest_length) * (n_norm - rest_length);


	//if (edge_index == 68) {
	//	//std::cout << delta_x.transpose() << std::endl;
	//	//std::cout << sn_0[0] << " " << sn_0[1] << " " << sn_0[2] << std::endl;
	//}

	//if (constraint_0 >= constraint_1) {
	//	//std::cout << "true " << std::endl;
	//}
	//else {
		//std::cout << "false "<<constraint_0<<" "<<constraint_1 << std::endl;
		//std::cout << delta_x.transpose() << std::endl;
		//std::cout << p0[0] - delta_x[0] << " " << p0[1] - delta_x[1] << " " << p0[2] - delta_x[2] << std::endl;
		//std::cout << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
		//std::cout << std::endl;
		//std::cout << p1[0] - delta_x[3] << " " << p1[1] - delta_x[4] << " " << p1[2] - delta_x[5] << std::endl;
		//std::cout << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
		//std::cout << edge_index << std::endl;
		//JacobiSVD<MatrixXd> svd;
		//svd.compute(H);
		//std::cout << H.determinant() << std::endl;
	//}
	

}



double SecondOrderConstraint::solveBendingConstraint(double* center_vertex, std::array<double, 3>* vertex_position, unsigned int* neighbor_vertex, unsigned int neighbor_vertex_size,
	double rest_Aq_norm, double lbo_weight, VectorXd& vertex_lbo)
{
	std::vector<VectorXd>q(3);
	//std::vector<VectorXd>q_initial(3);
	double aq[3];
	unsigned int size = neighbor_vertex_size + 1;
	VectorXd inv_m(size);
	for (unsigned int j = 0; j < 3; ++j) {
		q[j].resize(size);
	}
	q[0][0] = center_vertex[0];
	q[1][0] = center_vertex[1];
	q[2][0] = center_vertex[2];


	for (unsigned int h = 1; h < size; h++) {
		q[0][h] = vertex_position[*neighbor_vertex][0];
		q[1][h] = vertex_position[*neighbor_vertex][1];
		q[2][h] = vertex_position[*neighbor_vertex][2];
		neighbor_vertex++;
	}
	neighbor_vertex -= neighbor_vertex_size;

	aq[0] = vertex_lbo.dot(q[0]);
	aq[1] = vertex_lbo.dot(q[1]);
	aq[2] = vertex_lbo.dot(q[2]);

	double aq_norm = sqrt(DOT(aq, aq));
	if (aq_norm < epsilon_for_bending) {
		return 0;
	}
	//double K = aq_norm - rest_Aq_norm;
	VectorXd ATAx;
	ATAx.resize(3*size);
	for (unsigned int j = 0; j < 3; ++j) {
		ATAx.segment(j*size,size) = vertex_lbo * aq[j];
	}
	MatrixXd ATA_part = ((1.0- rest_Aq_norm/ aq_norm)*vertex_lbo) * vertex_lbo.transpose();
	MatrixXd H(3 * size, 3 * size);
	H = ((rest_Aq_norm / (aq_norm * aq_norm * aq_norm)) * ATAx) * ATAx.transpose();
	for (unsigned int i = 0; i < 3; ++i) {
		H.block(i * size, i * size, size, size) += ATA_part;
	}
	ATAx *= (aq_norm - rest_Aq_norm) / aq_norm;
	ColPivHouseholderQR <MatrixXd> linear(H);
	VectorXd delta_x = linear.solve(ATAx);

	//center_vertex[0] += delta_x[0];
	//center_vertex[1] += delta_x[size];
	//center_vertex[2] += delta_x[size+size];
	//for (unsigned int h = 1; h < size; h++) {
	//	vertex_position[*neighbor_vertex][0] += delta_x[h];
	//	vertex_position[*neighbor_vertex][1] += delta_x[size + h];
	//	vertex_position[*neighbor_vertex][2] += delta_x[size+size+h];
	//	neighbor_vertex++;
	//}


	//JacobiSVD<MatrixXd> svd;
	//svd.compute(H);
	//std::cout << svd.singularValues() << std::endl;


	q[0] -= delta_x.segment(0, size);
	q[1] -= delta_x.segment(size, size);
	q[2] -= delta_x.segment(size<<1, size);
	aq[0] = vertex_lbo.dot(q[0]);
	aq[1] = vertex_lbo.dot(q[1]);
	aq[2] = vertex_lbo.dot(q[2]);

	double constraint = sqrt(DOT(aq, aq)) - rest_Aq_norm;
	constraint = 0.5 * constraint * constraint;
	return constraint;
}





double SecondOrderConstraint::solveEdgeLengthConstraintFirstOrder(double* p0, double* p1, const double rest_length)
{
	double n[3];
	SUB(n, p0, p1);
	double n_norm = sqrt(DOT(n, n));
	DEV_(n, n_norm);

	double delta_lambda = -0.5* (n_norm - rest_length);//(1.0+gamma)*(inv_mass0 + inv_mass1) 
	MULTI_(n, delta_lambda);
	//SUM_(p0, n);
	//SUB_(p1, n);

	Vector3d p0_current;
	Vector3d p1_current;
	SUM(p0_current, p0, n);
	SUB(p1_current, p1, n);

	double constraint = (p0_current - p1_current).norm() - rest_length;
	constraint = 0.5 * constraint * constraint;

	//std::cout << p0_current[0] << " " << p0_current[1] << p0_current[2] << std::endl;
	//std::cout << p1_current[0] << " " << p1_current[1] << p1_current[2] << std::endl;


	return constraint;
}





double SecondOrderConstraint::solveBendingConstraintFirstOrder(double* center_vertex, std::array<double, 3>* vertex_position, unsigned int* neighbor_vertex, unsigned int neighbor_vertex_size,
	double rest_Aq_norm, double lbo_weight, VectorXd& vertex_lbo)
{
	//	energy = 0.0;
	std::vector<VectorXd>q(3);
	//std::vector<VectorXd>q_initial(3);
	double aq[3];
	unsigned int size = neighbor_vertex_size + 1;
	for (unsigned int j = 0; j < 3; ++j) {
		q[j].resize(size);
	}
	q[0][0] = center_vertex[0];
	q[1][0] = center_vertex[1];
	q[2][0] = center_vertex[2];


	for (unsigned int h = 1; h < size; h++) {
		q[0][h] = vertex_position[*neighbor_vertex][0];
		q[1][h] = vertex_position[*neighbor_vertex][1];
		q[2][h] = vertex_position[*neighbor_vertex][2];
		neighbor_vertex++;
	}
	neighbor_vertex -= neighbor_vertex_size;


	aq[0] = vertex_lbo.dot(q[0]);
	aq[1] = vertex_lbo.dot(q[1]);
	aq[2] = vertex_lbo.dot(q[2]);

	double aq_norm = sqrt(DOT(aq, aq));
	if (aq_norm < epsilon_for_bending) {
		return 0;
	}
	DEV_(aq, aq_norm);
	double K = aq_norm - rest_Aq_norm;
	//use grad to store delta_C

	std::vector<VectorXd>grad(3);
	//std::vector<VectorXd>q_initial(3);

	for (unsigned int j = 0; j < 3; ++j) {
		//q_initial[j] = q[j] - q_initial[j];
		grad[j] = vertex_lbo * aq[j];
	}

	
	double delta_lambda = -0.5*K
		/ (grad[0].squaredNorm() + grad[1].squaredNorm()+ grad[2].squaredNorm());

	center_vertex[0] += delta_lambda * grad[0][0];
	center_vertex[1] += delta_lambda * grad[1][0];
	center_vertex[2] += delta_lambda * grad[2][0];

	for (unsigned int h = 1; h < size; h++) {
		vertex_position[*neighbor_vertex][0] += delta_lambda * grad[0][h];
		vertex_position[*neighbor_vertex][1] += delta_lambda * grad[1][h];
		vertex_position[*neighbor_vertex][2] += delta_lambda * grad[2][h];
		neighbor_vertex++;
	}



	q[0] += delta_lambda* grad[0];
	q[1] += delta_lambda* grad[1];
	q[2] += delta_lambda* grad[2];

	aq[0] = vertex_lbo.dot(q[0]);
	aq[1] = vertex_lbo.dot(q[1]);
	aq[2] = vertex_lbo.dot(q[2]);

	double constraint = sqrt(DOT(aq, aq)) - rest_Aq_norm;
	constraint = 0.5 * constraint * constraint;

	return constraint;
}


void SecondOrderConstraint::test(MeshStruct& mesh_struct, std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo,
	std::vector<double>& rest_mean_curvature_norm)
{
	double constriant_1, constraint_2;

	//std::cout << "edge length " << std::endl;

	//for(unsigned int i = 0; i < mesh_struct.edge_length.size(); ++i) {
	//	constraint_2 = solveEdgeLengthConstraint(mesh_struct.vertex_position[mesh_struct.edge_vertices[i < 1]].data(),
	//		mesh_struct.vertex_position[mesh_struct.edge_vertices[(i < 1) + 1]].data(), mesh_struct.edge_length[i]);
	//	constriant_1 = solveEdgeLengthConstraintFirstOrder(mesh_struct.vertex_position[mesh_struct.edge_vertices[i < 1]].data(),
	//		mesh_struct.vertex_position[mesh_struct.edge_vertices[(i < 1) + 1]].data(), mesh_struct.edge_length[i]);
	//	if (constraint_2 < constriant_1) {
	//		std::cout << "true ";
	//	}
	//	else {
	//		std::cout << "false ";
	//	}
	//	std::cout << constraint_2 << " " << constriant_1 << std::endl;
	//}


	//std::cout << "bending " << std::endl;

	for (unsigned int i = 0; i < mesh_struct.vertex_position.size(); ++i) {
		constraint_2 = solveBendingConstraint(mesh_struct.vertex_position[i].data(), mesh_struct.vertex_position.data(), mesh_struct.vertices[i].neighbor_vertex.data(),
			mesh_struct.vertices[i].neighbor_vertex.size(), rest_mean_curvature_norm[i], lbo_weight[i], vertex_lbo[i]);
		constriant_1 = solveBendingConstraintFirstOrder(mesh_struct.vertex_position[i].data(), mesh_struct.vertex_position.data(), mesh_struct.vertices[i].neighbor_vertex.data(),
			mesh_struct.vertices[i].neighbor_vertex.size(), rest_mean_curvature_norm[i], lbo_weight[i], vertex_lbo[i]);
		if (constraint_2 < constriant_1) {
			std::cout << "true ";
		}
		else {
			std::cout << "false ";
		}
		std::cout << constraint_2 << " " << constriant_1 << std::endl;
	}
}