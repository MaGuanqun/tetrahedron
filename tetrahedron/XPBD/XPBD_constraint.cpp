#include"XPBD_constraint.h"

void XPBDconstraint::solveARAPConstraint(std::array<double, 3>* vertex_position, std::array<double, 3>* initial_vertex_position,
	double stiffness, double dt,
	Matrix<double,3,4>& A, int* vertex_index, double* inv_mass, double* lambda, const double damping_stiffness, double sigma_min,
	double sigma_max, double volume)
{
	Vector3d eigen_value;
	Matrix3d q_e;
	double determinant;
	Vector3d position;
	for (unsigned int i = 0; i < 3; ++i) {
		memcpy(q_e.data() + 3 * i, vertex_position[vertex_index[i+1]].data(), 24);
	}
	//first use eigen value to store the position of first vertex
	memcpy(eigen_value.data(), vertex_position[vertex_index[0]].data(), 24);
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

	if (determinant < 0) {
		eigen_value[2] = -eigen_value[2];
	}

	for (unsigned int j = 0; j < 3; ++j) {
		if (eigen_value[j] > 0) {
			if (eigen_value[j] < sigma_min) {
				eigen_value[j] = sigma_min;
			}
			else if (eigen_value[j] > sigma_max) {
				eigen_value[j] = sigma_max;
			}
		}
		else {
			if (eigen_value[j] > -sigma_min) {
				eigen_value[j] = -sigma_min;
			}
			else if (eigen_value[j] < -sigma_max) {
				eigen_value[j] = -sigma_max;
			}
		}
	}
	//std::cout << eigen_value << std::endl;
	//std::cout << "===" << std::endl;
	//use q_e as a temp vector
	for (unsigned int j = 0; j < 3; ++j) {
		for (unsigned int k = 0; k < 3; ++k) {
			q_e.data()[3 * j + k] = eigen_value[j] * svd.matrixU().data()[3 * j + k];
		}
	}
	//use P_inv to record transform
	P_inv = q_e * svd.matrixV().transpose();

	//get delta_c
	//Ax-q is deformation gradient-P_inv

	Vector4d delta_c_transpose;
	Vector4d position_;
	//deformation_gradient -= P_inv;
	double alpha_ = 1.0 / (stiffness * dt * dt);
	double gamma = 0.0;// = damping_stiffness / (stiffness * dt);
 
	// to limit the deformation on three dimensions saparately, 
	//we add connstraints on every dimension  
	//C_x = |Aq_x|-|T_x|, C_y=|Aq_y|-|T_y|, C_z=|Aq_z|-|T_z|

	double C_;
	double delta_lambda;
	for (unsigned int k = 0; k < 3; ++k) {
		C_ = deformation_gradient.row(k).norm();
		for (unsigned int i = 0; i < 4; ++i) {
			delta_c_transpose.data()[i] = (volume / C_) * A.col(i).dot(deformation_gradient.row(k));
			position_.data()[i] = vertex_position[vertex_index[i]][k]- initial_vertex_position[vertex_index[i]][k];
		}

		C_ = volume*(C_ - P_inv.row(k).norm());		
		delta_lambda = -(C_ + alpha_ * lambda[k] + gamma * delta_c_transpose.dot(position_))
			/ ((1.0 + gamma) *
				(inv_mass[vertex_index[0]] * delta_c_transpose.data()[0] * delta_c_transpose.data()[0]
					+ inv_mass[vertex_index[1]] * delta_c_transpose.data()[1]* delta_c_transpose.data()[1]
					+ inv_mass[vertex_index[2]] * delta_c_transpose.data()[2]* delta_c_transpose.data()[2]
					+ inv_mass[vertex_index[3]] * delta_c_transpose.data()[3]* delta_c_transpose.data()[3])
				+alpha_);

		//if (abs(delta_lambda) > 1e1) {
		//	std::cout << delta_lambda << " " << k << std::endl;
		//}

		lambda[k] += delta_lambda;

		for (unsigned int i = 0; i < 4; ++i) {
			delta_c_transpose.data()[i] *= (inv_mass[i] * delta_lambda);			
			vertex_position[vertex_index[i]].data()[k] += delta_c_transpose.data()[i];
		}
	}

	//double C;
	//for (unsigned int k = 0; k < 3; ++k) {
	//	C = volume * (deformation_gradient.row(k).norm() - P_inv.row(k).norm());
	//	delta_c_transpose = ((volume / C) * deformation_gradient).transpose() * A;
	//	for (unsigned int i = 0; i < 4; ++i) {
	//		position.data()[3 * i] -= initial_vertex_position[vertex_index[i]][0];
	//		position.data()[3 * i + 1] -= initial_vertex_position[vertex_index[i]][1];
	//		position.data()[3 * i + 2] -= initial_vertex_position[vertex_index[i]][2];
	//	}
	//	double delta_lambda = -(C + alpha_ * lambda + gamma * delta_c_transpose.cwiseProduct(position).sum())
	//		/ ((1.0 + gamma) *
	//			(inv_mass[vertex_index[0]] * delta_c_transpose.col(0).squaredNorm() + inv_mass[vertex_index[1]] * delta_c_transpose.col(1).squaredNorm()
	//				+ inv_mass[vertex_index[2]] * delta_c_transpose.col(2).squaredNorm() + inv_mass[vertex_index[3]] * delta_c_transpose.col(3).squaredNorm()));
	//	lambda += delta_lambda;
	//	for (unsigned int i = 0; i < 4; ++i) {
	//		delta_c_transpose.col(i) *= (inv_mass[i] * delta_lambda);
	//		SUM_(vertex_position[vertex_index[i]].data(), (delta_c_transpose.data() + 3 * i));
	//	}
	//}
}


void XPBDconstraint::solveTetStrainConstraint(std::array<double, 3>* vertex_position, std::array<double, 3>* initial_vertex_position,
	double stiffness, double dt, Matrix<double, 3, 4>& A, int* vertex_index, double* inv_mass, double& lambda, const double damping_stiffness, double volume,
	double youngs_modulus, double poisson_ratio)
{

	Matrix<double, 3, 4> grad_C;
	double C;
	double mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
	double lambda_ = youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
	computeGreenStrainAndPiolaStress(vertex_position[vertex_index[0]].data(), vertex_position[vertex_index[1]].data(),
		vertex_position[vertex_index[2]].data(), vertex_position[vertex_index[3]].data(), A, volume, mu, lambda_, grad_C, C);
	double alpha_ = 1.0 / (stiffness * dt * dt);
	double gamma = damping_stiffness / (stiffness * dt);

	double sum_normGradC = inv_mass[vertex_index[0]] * grad_C.col(0).squaredNorm() +
		inv_mass[vertex_index[1]] * grad_C.col(1).squaredNorm() +
		inv_mass[vertex_index[2]] * grad_C.col(2).squaredNorm() +
		inv_mass[vertex_index[3]] * grad_C.col(3).squaredNorm();

	Vector3d temp;
	//k = delta_C(x-x_initial)

	double k = 0;
	for (unsigned int i = 0; i < 4; ++i) {
		SUB(temp, vertex_position[vertex_index[i]], initial_vertex_position[vertex_index[i]]);
		k += grad_C.col(i).dot(temp);
	}
	double delta_lambda = -(C + alpha_ * lambda + gamma * k) / ((1.0 + gamma) * sum_normGradC + alpha_);
	lambda += delta_lambda;

	double coe;
	int vertex_index_;
	for (unsigned int i = 0; i < 4; ++i){
		vertex_index_ = vertex_index[i];
		coe = inv_mass[vertex_index_] * delta_lambda;
		ACCUMULATE_SUM_WITH_COE(vertex_position[vertex_index_], coe, (grad_C.data() + 3 * i));
	}
}






void XPBDconstraint::computeGreenStrainAndPiolaStress(
	double* v0, double* v1, double* v2, double* v3,
	Matrix<double, 3, 4>& inv_rest_pos,
	double rest_volume,
	double mu, double lambda, Matrix<double, 3, 4>& grad_C, double& C)
{
	Matrix3d F;
	Matrix3d P;
	SUB(P.data(), v1, v0);
	SUB((P.data() + 3), v2, v0);
	SUB((P.data() + 6), v3, v0);
	//first use eigen value to store the position of first vertex
	Matrix3d P_inv;
	memcpy(P_inv.data(), inv_rest_pos.data() + 3, 72);
	F = P * P_inv.transpose();

	if (F.determinant() < 0.58) {
		JacobiSVD<Matrix3d> svd;
		svd.compute(F, ComputeFullU | ComputeFullV);
		Vector3d new_singular_value = svd.singularValues();
		double min_s_value = 0.58;
		//once a critical compression
		//threshold is reached( 58 % of undeformed dimensions, when compression occurs
		//	along a single axis) the strength of the restorative force reaches a maximum.
		//Further compression will be met with decreasing resistance

		for (unsigned char j = 0; j < 3; j++) {
			if (new_singular_value[j] < min_s_value)
				new_singular_value[j] = min_s_value;
		}
		F = svd.matrixU() * new_singular_value.asDiagonal() * svd.matrixV().transpose();
	}

	Matrix3d E;
	//E= 0.5* (F^T F - I)
	E.data()[0] = 0.5 * (F.data()[0] * F.data()[0] + F.data()[1] * F.data()[1] + F.data()[2] * F.data()[2] - 1.0);		// xx
	E.data()[4] = 0.5 * (F.data()[3] * F.data()[3] + F.data()[4] * F.data()[4] + F.data()[5] * F.data()[5] - 1.0);		// yy
	E.data()[8] = 0.5 * (F.data()[6] * F.data()[6] + F.data()[7] * F.data()[7] + F.data()[8] * F.data()[8] - 1.0);		// zz
	E.data()[3] = 0.5 * (F.data()[0] * F.data()[3] + F.data()[1] * F.data()[4] + F.data()[2] * F.data()[5]);			// xy
	E.data()[6] = 0.5 * (F.data()[0] * F.data()[6] + F.data()[1] * F.data()[7] + F.data()[2] * F.data()[8]);			// xz
	E.data()[7] = 0.5 * (F.data()[3] * F.data()[6] + F.data()[4] * F.data()[7] + F.data()[5] * F.data()[8]);			// yz
	E.data()[1] = E.data()[3];
	E.data()[2] = E.data()[6];
	E.data()[5] = E.data()[7];

	//P(F)=F(2muE +lambda tr(E)I) 
	double trace = E.data()[0] + E.data()[4] + E.data()[8];
	double lambda_trace = lambda * trace;
	P = (mu + mu) * E;
	P.data()[0] += lambda_trace;
	P.data()[4] += lambda_trace;
	P.data()[8] += lambda_trace;
	P = F * P;
	//st. venant-kirchhoff model
	// psi= mu E:E + lambda/2 * tr(E)^2
	double psi = mu * E.squaredNorm() + 0.5 * lambda * trace * trace;
	C = rest_volume * psi;
	grad_C = rest_volume * P * inv_rest_pos;
}



void XPBDconstraint::solveEdgeLengthConstraint(double* p0, double* p1, const double rest_length, const double stiffness, double dt,
	double inv_mass0, double inv_mass1, double& lambda, const double damping_stiffness, double* initial_p0, double* inital_p1)
{	
	double gamma = damping_stiffness / (stiffness*dt);
	dt = dt * dt;
	double C = sqrt(EDGE_LENGTH(p0, p1)) - rest_length;
	double n[3];
	SUB(n, p0, p1);
	double n_[3];
	double n_norm = sqrt(DOT(n, n));
	double n_coeff = gamma / n_norm;
	MULTI(n_, n, n_coeff); 

	double delta_lambda = -(C + lambda / (stiffness * dt) + (n_[0]*(p0[0] - initial_p0[0] - p1[0] + inital_p1[0] )
		+ n_[1] * (p0[1] - initial_p0[1] - p1[1] + inital_p1[1])+ n_[2] * (p0[2] - initial_p0[2] - p1[2] + inital_p1[2]))) 
		/ ((1.0+gamma)*(inv_mass0 + inv_mass1) + 1.0 / (stiffness * dt));
	lambda += delta_lambda;
	
	
	n_coeff = delta_lambda / n_norm;

	MULTI_(n, n_coeff);
	ACCUMULATE_SUM_WITH_COE(p0, inv_mass0, n);
	inv_mass1 *= -1.0;
	ACCUMULATE_SUM_WITH_COE(p1, inv_mass1, n);
}



void XPBDconstraint::solveBendingConstraint(double* center_vertex, double vertex_inv_mass,  std::array<double,3>* vertex_position, std::vector<unsigned int>& neighbor_vertex,
	double rest_curvature_norm, double lbo_weight, VectorXd& vertex_lbo,  double stiffness, double dt, double* inv_mass, double &lambda,
	const double damping_stiffness, double* initial_center_vertex, std::array<double, 3>* inital_vertex_position)
{
	std::vector<VectorXd>q(3);
	std::vector<VectorXd>q_initial(3);
	double aq[3];
	unsigned int size = neighbor_vertex.size() + 1;
	VectorXd inv_m(size);
	for (unsigned int j = 0; j < 3; ++j) {
		q[j].resize(size);
		q_initial[j].resize(size);
	}
	q[0][0] = center_vertex[0];
	q[1][0] = center_vertex[1];
	q[2][0] = center_vertex[2];

	q_initial[0][0] = initial_center_vertex[0];
	q_initial[1][0] = initial_center_vertex[1];
	q_initial[2][0] = initial_center_vertex[2];

	inv_m[0] = vertex_inv_mass;

	unsigned int index;
	for (unsigned int h = 1; h < size; h++) {
		index = neighbor_vertex[h - 1];
		q[0][h] = vertex_position[index][0];
		q[1][h] = vertex_position[index][1];
		q[2][h] = vertex_position[index][2];

		q_initial[0][h] = inital_vertex_position[index][0];
		q_initial[1][h] = inital_vertex_position[index][1];
		q_initial[2][h] = inital_vertex_position[index][2];

		inv_m[h] = inv_mass[index];
	}
	aq[0] = vertex_lbo.dot(q[0]);
	aq[1] = vertex_lbo.dot(q[1]);
	aq[2] = vertex_lbo.dot(q[2]);
	double aq_norm = sqrt(DOT(aq, aq));
	DEV_(aq, aq_norm);
	double C = aq_norm - rest_curvature_norm;

	//use q to store delta_C
	for (unsigned int j = 0; j < 3; ++j) {
		q_initial[j] = q[j] - q_initial[j];
		q[j] = vertex_lbo * aq[j];
	}
	double alpha_ = lbo_weight / (stiffness * dt * dt);
	double gamma = lbo_weight * damping_stiffness / (stiffness * dt);
	double delta_lambda = -lbo_weight * (C + alpha_ * lambda + gamma * (q[0].dot(q_initial[0]) + q[1].dot(q_initial[1]) + q[2].dot(q_initial[2])))
		/ ((1.0+gamma)*(q[0].dot(q[0].cwiseProduct(inv_m)) + q[1].dot(q[1].cwiseProduct(inv_m))
			+ q[2].dot(q[2].cwiseProduct(inv_m))) + lbo_weight * alpha_);
	
	lambda += delta_lambda;

	inv_m *= (delta_lambda/lbo_weight);
	center_vertex[0] += inv_m.data()[0] * q[0][0];
	center_vertex[1] += inv_m.data()[0] * q[1][0];
	center_vertex[2] += inv_m.data()[0] * q[2][0];

	for (unsigned int h = 1; h < size; h++) {
		vertex_position[neighbor_vertex[h - 1]][0] += inv_m.data()[h] * q[0][h];
		vertex_position[neighbor_vertex[h - 1]][1] += inv_m.data()[h] * q[1][h];
		vertex_position[neighbor_vertex[h - 1]][2] += inv_m.data()[h] * q[2][h];
	}

}

void XPBDconstraint::initial_LBO_EdgeCotWeight(TriangleMeshStruct& mesh_struct, std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo,
	std::vector<double>& rest_mean_curvature_norm)
{
	std::vector<double> edge_cot_weight;
	initialEdgeCotWeight(mesh_struct, edge_cot_weight);
	computeLBOWeight(lbo_weight, mesh_struct);
	computeVertexLBO(mesh_struct, vertex_lbo, edge_cot_weight);
	restBendingMeanCurvature(mesh_struct, rest_mean_curvature_norm, vertex_lbo, lbo_weight);
}

void XPBDconstraint::initialEdgeCotWeight(TriangleMeshStruct& mesh_struct, std::vector<double>& edge_cot_weight)
{
	int edge_num;
	edge_num = mesh_struct.edges.size();
	edge_cot_weight.resize(edge_num);
	double cotan0, cotan1;
	double len10, len20, len13, len23;
	double x10[3], x20[3];
	double x13[3], x23[3];
	double theta0, theta1;
	unsigned int edge_vertex_0;
	unsigned int edge_vertex_1;
	unsigned int opposite_0;
	unsigned int opposite_1;
	for (int i = 0; i < edge_num; ++i) {
		
		edge_vertex_0 = mesh_struct.edge_vertices[i << 1];
		edge_vertex_1 = mesh_struct.edge_vertices[(i << 1) + 1];
		opposite_0 = mesh_struct.edges[i].opposite_vertex[0];
	
		cotan0 = 0;
		cotan1 = 0;
		SUB(x10, mesh_struct.vertex_position[edge_vertex_0], mesh_struct.vertex_position[opposite_0]);
		SUB(x20, mesh_struct.vertex_position[edge_vertex_1], mesh_struct.vertex_position[opposite_0]);
		len10 = sqrt(DOT(x10, x10));
		len20 = sqrt(DOT(x20, x20));
		theta0 = acos(DOT(x10, x20) / (len10 * len20));
		cotan0 = 1.0 / tan(theta0);

		if (mesh_struct.edges[i].opposite_vertex.size() > 1) {
			opposite_1 = mesh_struct.edges[i].opposite_vertex[1];

			SUB(x13, mesh_struct.vertex_position[edge_vertex_0], mesh_struct.vertex_position[opposite_1]);
			SUB(x23, mesh_struct.vertex_position[edge_vertex_1], mesh_struct.vertex_position[opposite_1]);
			len13 = sqrt(DOT(x13, x13));
			len23 = sqrt(DOT(x23, x23));
			theta1 = acos(DOT(x13, x23) / (len13 * len23));
			cotan1 = 1.0 / tan(theta1);
			edge_cot_weight[i] = -0.5 * (cotan0 + cotan1);
		}
		else
		{
			edge_cot_weight[i] = -0.5 * cotan0;
		}
	}
}

void XPBDconstraint::computeLBOWeight(std::vector<double>& lbo_weight, TriangleMeshStruct& mesh_struct)
{
	double m;
	lbo_weight.resize(mesh_struct.vertex_position.size(),0.0);
	for (int i = 0; i < mesh_struct.faces.size(); ++i) {
		m = mesh_struct.faces[i].area / 3.0;
		lbo_weight[mesh_struct.triangle_indices[i][0]] += m;
		lbo_weight[mesh_struct.triangle_indices[i][1]] += m;
		lbo_weight[mesh_struct.triangle_indices[i][2]] += m;
	}
}





void XPBDconstraint::computeVertexLBO(TriangleMeshStruct& mesh_struct, std::vector<VectorXd>& vertex_lbo, std::vector<double>& edge_cot_weight)
{
	double total;
	double edge_weight;
	vertex_lbo.resize(mesh_struct.vertex_position.size());
	for (int i = 0; i < mesh_struct.vertex_position.size(); ++i) {
		vertex_lbo[i].resize(mesh_struct.vertices[i].edge.size() + 1);
		vertex_lbo[i].setZero();
		total = 0.0;
		for (int j = 0; j < mesh_struct.vertices[i].edge.size(); ++j) {
			edge_weight = edge_cot_weight[mesh_struct.vertices[i].edge[j]];
			total += edge_weight;
			vertex_lbo[i].data()[j + 1] = edge_weight;
		}
		vertex_lbo[i].data()[0] = -total;
	}
}

void XPBDconstraint::restBendingMeanCurvature(TriangleMeshStruct& mesh_struct, std::vector<double>& rest_mean_curvature_norm,
	std::vector<VectorXd>& vertex_lbo, std::vector<double>& lbo_weight)
{
	rest_mean_curvature_norm.resize(mesh_struct.vertex_position.size());
	VectorXd q;
	unsigned int size;
	double vertex_curvature[3];
	unsigned int* neighbor_vertex;
	for (unsigned int i = 0; i < rest_mean_curvature_norm.size(); ++i) {
		size = vertex_lbo[i].size();
		q.resize(size);
		neighbor_vertex = mesh_struct.vertices[i].neighbor_vertex.data();
		for (unsigned int j = 0; j < 3; ++j) {
			q[0] = mesh_struct.vertex_position[i][j];
			for (unsigned int h = 1; h < size; h++) {
				q[h] = mesh_struct.vertex_position[neighbor_vertex[h - 1]][j];
			}
			vertex_curvature[j] = dotProductX_(vertex_lbo[i], q);
		}
		rest_mean_curvature_norm[i] = sqrt(DOT(vertex_curvature, vertex_curvature));
	}

}