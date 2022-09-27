#ifndef SECOND_ORDER
#define SECOND_ORDER
#include"../basic/global.h"
#include"../external/Eigen/Dense"
#include<array>

#include"../mesh_struct/triangle_mesh_struct.h"

using namespace Eigen;
class SecondOrderConstraint
{
public:
	void solveEdgeLengthConstraint(double* p0, double* p1, const double rest_length, double stiffness, double mass_0, double mass_1, double time_step, double* sn_0, double* sn_1,
		bool v0_fixed, bool v1_fixed, unsigned int edge_index, double& lambda, unsigned int vertex_0_index, unsigned int vertex_1_index);
	double solveBendingConstraint(double* center_vertex, std::array<double, 3>* vertex_position, unsigned int* neighbor_vertex,
		unsigned int neighbor_vertex_size,
		double rest_Aq_norm, double lbo_weight, VectorXd& vertex_lbo);
	double solveEdgeLengthConstraintFirstOrder(double* p0, double* p1, const double rest_length);

	double solveBendingConstraintFirstOrder(double* center_vertex,std::array<double, 3>* vertex_position, unsigned int* neighbor_vertex, unsigned int neighbor_vertex_size,
		double rest_Aq_norm, double lbo_weight, VectorXd& vertex_lbo);
	void test(MeshStruct& mesh_struct, std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo,
		std::vector<double>& rest_mean_curvature_norm);

	void solveARAPConstraint(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
		double stiffness, double dt,
		Matrix<double, 3, 4>& A, double* inv_mass, double& lambda, const double damping_stiffness, double* mass,
		double volume);
		

	void solveCD_ARAP(std::array<double,3>* vertex_position, double stiffness, double dt,
		Matrix<double, 3, 4>* A, double* lambda, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, 
		double* mass, double* volume, unsigned int vertex_index, std::array<double, 3>* sn);

	unsigned int findVertexNo(int vertex_index, int* indices);
	 
	bool getARAPGradHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3, 
		Matrix<double, 3, 4>& A, Matrix3d& Hessian, Vector3d& grad, double& C, unsigned int vertex_no);

	bool getARAPGradHessianNewton(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
		Matrix<double, 3, 4>& A, Matrix3d& Hessian, Vector3d& grad, double& C, unsigned int vertex_no);

	void computeEdgeLengthForce(double* vertex_0, double* vertex_1, double stiffness,
		double* potential_0, double* potential_1, double rest_length);
	void computeARAPForce(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
		double stiffness, Matrix<double, 3, 4>& A, double volume, Matrix<double, 3, 4>& force);

	void solveNewtonCD_ARAP(std::array<double, 3>* vertex_position, double stiffness, double dt,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
		double* volume, unsigned int vertex_index, std::array<double, 3>* sn);

private:
	double epsilon_for_bending= 1e-10;

	void solveEdgeLengthConstraint(Vector3f& p1, Vector3f& p2, const double d, double mass_0,
		double mass_1, Vector3f& ori_p1, Vector3f& ori_p2, bool v0_fixed, bool v1_fixed, double& lambda);


};


#endif // !SECOND_ORDER


