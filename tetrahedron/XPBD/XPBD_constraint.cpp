#include"XPBD_constraint.h"

void XPBDconstraint::solveEdgeLengthConstraint(double* p0, double* p1, const double rest_length, const double stiffness, double dt,
	double inv_mass0, double inv_mass1, double& lambda)
{
	dt = dt *dt;
	double C = sqrt(EDGE_LENGTH(p0, p1)) - rest_length;
	//std::cout << p0[0]<<" "<<p1[0] << std::endl;
	double delta_lambda = -(C + lambda / (stiffness * dt)) / (inv_mass0 + inv_mass1 + 1.0 / (stiffness * dt));
	lambda += delta_lambda;
	double n[3];
	SUB(n, p0, p1);
	double n_coeff = delta_lambda / sqrt(DOT(n, n));
	//std::cout << "nn coef " << sqrt(DOT(n, n))<<" "<<n[0] << std::endl;

	MULTI_(n, n_coeff);
	//std::cout << "n " << n[0] << std::endl;
	ACCUMULATE_SUM_WITH_COE(p0, inv_mass0, n);
	inv_mass1 *= -1.0;
	ACCUMULATE_SUM_WITH_COE(p1, inv_mass1, n);
	//std::cout << p0[0] << " " << p1[0] << std::endl;
}

void XPBDconstraint::solveBendingConstraint(double* center_vertex, double vertex_inv_mass,  std::array<double,3>* vertex_position, std::vector<unsigned int>& neighbor_vertex,
	double rest_curvature_norm, double lbo_weight, VectorXd& vertex_lbo,  double stiffness, double dt, double* inv_mass, double &lambda)
{
	std::vector<VectorXd>q(3);
	double aq[3];
	unsigned int size = neighbor_vertex.size() + 1;
	VectorXd inv_m(size);
	for (unsigned int j = 0; j < 3; ++j) {
		q[j].resize(size);
	}
	q[0][0] = center_vertex[0];
	q[1][0] = center_vertex[1];
	q[2][0] = center_vertex[2];
	inv_m[0] = vertex_inv_mass;
	for (unsigned int h = 1; h < size; h++) {
		q[0][h] = vertex_position[neighbor_vertex[h-1]][0];
		q[1][h] = vertex_position[neighbor_vertex[h-1]][1];
		q[2][h] = vertex_position[neighbor_vertex[h-1]][2];
		inv_m[h] = inv_mass[neighbor_vertex[h - 1]];
	}
	aq[0] = vertex_lbo.dot(q[0]);
	aq[1] = vertex_lbo.dot(q[1]);
	aq[2] = vertex_lbo.dot(q[2]);
	double aq_norm = sqrt(DOT(aq, aq));
	DEV_(aq, aq_norm);
	double C = aq_norm - rest_curvature_norm;

	//use q to store delta_C
	for (unsigned int j = 0; j < 3; ++j) {
		q[j] = vertex_lbo * aq[j];
	}
	double alpha_ = lbo_weight / (stiffness * dt * dt);
	double delta_lambda = -(C + alpha_ * lambda) / (q[0].dot(q[0].cwiseProduct(inv_m)) + q[1].dot(q[1].cwiseProduct(inv_m))
		+ q[2].dot(q[2].cwiseProduct(inv_m)) + alpha_);
	
	//if (delta_lambda<5 && delta_lambda>-5) {

	//}
	//else {
	//	//std::cout << delta_lambda << " " << q[0].dot(q[0].cwiseProduct(inv_m)) + q[1].dot(q[1].cwiseProduct(inv_m))
	//	//	+ q[2].dot(q[2].cwiseProduct(inv_m)) << std::endl;
	//	std::cout << vertex_lbo << std::endl;
	//	std::cout << "=====" << std::endl;
	//}
	
	lambda += delta_lambda;

	inv_m *= delta_lambda;
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