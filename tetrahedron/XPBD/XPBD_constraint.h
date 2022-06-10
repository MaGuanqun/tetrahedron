#pragma once
#include"../basic/global.h"
#include"../mesh_struct/triangle_mesh_struct.h"
#include"../external/Eigen/Dense"
#include"../basic/eigenDenseOperation.h"
using namespace Eigen;
using namespace denseOperation;

class XPBDconstraint
{
public:
	void solveEdgeLengthConstraint(double* p0, double* p1, const double rest_length, const double stiffness, double dt,
		double inv_mass0, double inv_mass1, double& lambda, const double damping_stiffness, double* initial_p0, double* inital_p1);
	void initial_LBO_EdgeCotWeight(TriangleMeshStruct& mesh_struct, std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo,
		std::vector<double>& rest_mean_curvature_norm);
	void solveBendingConstraint(double* center_vertex, double vertex_inv_mass, std::array<double, 3>* vertex_position, std::vector<unsigned int>& neighbor_vertex,
		double rest_curvature_norm, double lbo_weight, VectorXd& vertex_lbo, double stiffness, double dt, double* inv_mass, double& lambda,
		const double damping_stiffness, double* initial_center_vertex, std::array<double, 3>* inital_vertex_position);
	void solveARAPConstraint(std::array<double, 3>* vertex_position, std::array<double, 3>* initial_vertex_position,
		double stiffness, double dt,
		Matrix<double, 3, 4>& A, int* vertex_index, double* inv_mass, double& lambda, const double damping_stiffness,
		double sigma_min, double sigma_max, double volume);
	void solveARAPConstraint2(std::array<double, 3>* original_vertex_pos, std::array<double, 3>* vertex_position, std::array<double, 3>* initial_vertex_position,
		double stiffness, double dt,
		Matrix<double, 3, 4>& A, int* vertex_index, double* inv_mass, double& lambda, const double damping_stiffness, double sigma_min,
		double sigma_max, double volume);

	void solveTetStrainConstraint(std::array<double, 3>* vertex_position, std::array<double, 3>* initial_vertex_position,
		double stiffness, double dt, Matrix<double, 3, 4>& A, int* vertex_index, double* inv_mass, double& lambda, const double damping_stiffness, double volume,
		double youngs_modulus, double poisson_ratio);

	double epsilon_for_bending;// set min value of ||Aq-Ap||. avoid  /0

private:
	void initialEdgeCotWeight(TriangleMeshStruct& mesh_struct, std::vector<double>& edge_cot_weight);
	void computeLBOWeight(std::vector<double>& lbo_weight, TriangleMeshStruct& mesh_struct);
	void computeVertexLBO(TriangleMeshStruct& mesh_struct, std::vector<VectorXd>& vertex_lbo, std::vector<double>& edge_cot_weight);
	void restBendingMeanCurvature(TriangleMeshStruct& mesh_struct, std::vector<double>& rest_mean_curvature_norm,
		std::vector<VectorXd>& vertex_lbo, std::vector<double>& lbo_weight);

	void computeGreenStrainAndPiolaStress(
		double* v0, double* v1, double* v2, double* v3,
		Matrix<double, 3, 4>& inv_rest_pos,
		double rest_volume,
		double mu, double lambda, Matrix<double, 3, 4>& grad_C, double& C);

};
