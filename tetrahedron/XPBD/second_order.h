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
		bool v0_fixed, bool v1_fixed, unsigned int edge_index, double& lambda);
	double solveBendingConstraint(double* center_vertex, std::array<double, 3>* vertex_position, unsigned int* neighbor_vertex,
		unsigned int neighbor_vertex_size,
		double rest_Aq_norm, double lbo_weight, VectorXd& vertex_lbo);
	double solveEdgeLengthConstraintFirstOrder(double* p0, double* p1, const double rest_length);

	double solveBendingConstraintFirstOrder(double* center_vertex,std::array<double, 3>* vertex_position, unsigned int* neighbor_vertex, unsigned int neighbor_vertex_size,
		double rest_Aq_norm, double lbo_weight, VectorXd& vertex_lbo);
	void test(MeshStruct& mesh_struct, std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo,
		std::vector<double>& rest_mean_curvature_norm);
private:
	double epsilon_for_bending= 1e-10;
};


#endif // !SECOND_ORDER


