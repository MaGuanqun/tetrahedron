#include"second_order.h"

void SecondOrderConstraint::solveEdgeLengthConstraint(double* p0, double* p1, const double rest_length, double stiffness, double mass_0,
	double mass_1, double time_step, double* sn_0, double* sn_1, bool v0_fixed, bool v1_fixed, unsigned int edge_index)
{
	Vector3d n;
	SUB(n.data(), p0, p1);
	double n_norm = sqrt(DOT(n, n));
	double coe =stiffness * (1.0 - rest_length / n_norm);
	Vector3d grad_ = coe * n;
	Matrix3d He = ((rest_length*stiffness / (n_norm * n_norm * n_norm)) * n) * n.transpose();
	He.data()[0] += coe;
	He.data()[4] += coe;
	He.data()[8] += coe;

	if (v0_fixed) {
		if (v1_fixed) {
			return;
		}
		else {
			coe = mass_1 / (time_step * time_step);;
			He.data()[0] += coe;
			He.data()[4] += coe;
			He.data()[8] += coe;
			grad_[0] -= coe * (p1[0] - sn_1[0]);
			grad_[1] -= coe * (p1[1] - sn_1[1]);
			grad_[2] -= coe * (p1[2] - sn_1[2]);
			ColPivHouseholderQR <Matrix<double, 3, 3>> linear(He);
			Matrix<double, 3, 1> delta_x = linear.solve(grad_);
			SUM_(p1, delta_x);
			return;
		}
	}
	else {
		if (v1_fixed) {
			coe = mass_0 / (time_step * time_step);
			He.data()[0] += coe;
			He.data()[4] += coe;
			He.data()[8] += coe;
			grad_[0] += coe * (p0[0] - sn_0[0]);
			grad_[1] += coe * (p0[1] - sn_0[1]);
			grad_[2] += coe * (p0[2] - sn_0[2]);
			ColPivHouseholderQR <Matrix<double, 3, 3>> linear(He);
			Matrix<double, 3, 1> delta_x = linear.solve(grad_);
			SUB_(p0, delta_x);
			return;
		}
	}

	Matrix<double, 6, 1> neg_grad;
	neg_grad.segment(0, 3) = -grad_;
	neg_grad.segment(3, 3) = grad_;

	Matrix<double,6,6> H;
	H.block<3,3>(0, 0) = He;
	H.block<3,3>(3, 3) = He;

	for (unsigned int i = 0; i < 9; ++i) {
		He.data()[i] = 0.0 - He.data()[i];
	}
	H.block<3,3>(3, 0) = He;
	H.block<3,3>(0, 3) = He;	

	coe = mass_0 / (time_step * time_step);
	H.data()[0] += coe;
	H.data()[7] += coe;
	H.data()[14] += coe;
	neg_grad.data()[0] -= coe * (p0[0] - sn_0[0]);
	neg_grad.data()[1] -= coe * (p0[1] - sn_0[1]);
	neg_grad.data()[2] -= coe * (p0[2] - sn_0[2]);
	
	coe= mass_1 / (time_step * time_step);
	H.data()[21] += coe;
	H.data()[28] += coe;
	H.data()[35] += coe;
	neg_grad.data()[3] -= coe * (p1[0] - sn_1[0]);
	neg_grad.data()[4] -= coe * (p1[1] - sn_1[1]);
	neg_grad.data()[5] -= coe * (p1[2] - sn_1[2]);
	//memcpy(H.data(), He.data(), 24);
	//memcpy(H.data() + 21, He.data(), 24);
	// 	//memcpy(H.data() + 6, He.data() + 3, 24);
	//memcpy(H.data() + 27, He.data() + 3, 24);
	//memcpy(H.data() + 12, He.data() + 6, 24);
	//memcpy(H.data() + 33, He.data() + 6, 24);

	ColPivHouseholderQR <Matrix<double, 6, 6>> linear(H);
	Matrix<double, 6, 1> delta_x = linear.solve(neg_grad);

	JacobiSVD<MatrixXd> svd;
	svd.compute(H);
	// 	std::cout << svd.singularValues() << std::endl;
	//std::cout << H.determinant() << std::endl;
	//std::cout << H << std::endl;


	double constraint_0 = 1.0 / (2.0 * time_step * time_step) * (mass_0 * EDGE_LENGTH(p0, sn_0) + mass_1 * EDGE_LENGTH(p1, sn_1)) + 0.5 *stiffness* (n_norm - rest_length) * (n_norm - rest_length);


	SUM_(p0, delta_x);
	p1[0] += delta_x[3];
	p1[1] += delta_x[4];
	p1[2] += delta_x[5];



	//Vector3d p0_current;
	//Vector3d p1_current;

	//SUM(p0_current, p0, delta_x);
	//p1_current[0]=p1[0] + delta_x[3];
	//p1_current[1]=p1[1] + delta_x[4];
	//p1_current[2]=p1[2] + delta_x[5];

	//double constraint =(p0_current - p1_current).norm() - rest_length;
	//constraint = 0.5 * constraint * constraint;
	n_norm = sqrt(EDGE_LENGTH(p0, p1));
	double constraint_1  = 1.0 / (2.0 * time_step * time_step) * (mass_0 * EDGE_LENGTH(p0, sn_0) + mass_1 * EDGE_LENGTH(p1, sn_1)) + 0.5 * stiffness * (n_norm - rest_length) * (n_norm - rest_length);


	if (edge_index == 68) {
		//std::cout << delta_x.transpose() << std::endl;
		//std::cout << sn_0[0] << " " << sn_0[1] << " " << sn_0[2] << std::endl;
	}

	if (constraint_0 >= constraint_1) {
		//std::cout << "true " << std::endl;
	}
	else {
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

	}
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