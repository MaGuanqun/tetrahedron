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

	bool solve_exact_ARAP_hessian=false;

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

	void solveSingleVertexNewton(std::array<double, 3>* vertex_position, double stiffness, double dt,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
		double* volume, unsigned int vertex_index, std::array<double, 3>* sn);
		
	void solveSingleVertexCD_ARAP(std::array<double, 3>* vertex_position, double stiffness, double dt,
		Matrix<double, 3, 4>* A, double* lambda, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices,
		double* mass, double* volume, unsigned int vertex_index, std::array<double, 3>* sn);

	void solveCD_ARAP(std::array<double,3>* vertex_position, double stiffness, double dt,
		Matrix<double, 3, 4>* A, double* lambda, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, 
		double* mass, double* volume, unsigned int vertex_index, std::array<double, 3>* sn);
	 
	bool getARAPGradHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3, 
		Matrix<double, 3, 4>& A, Matrix3d& Hessian, Vector3d& grad, double& C, unsigned int vertex_no);

	bool getARAPGradHessianNewton(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
		Matrix<double, 3, 4>& A, Matrix3d& Hessian, Vector3d& grad, double& C, unsigned int vertex_no);

	bool getCollisionPairHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
		double ori_volume, Matrix3d& Hessian, Vector3d& grad, unsigned int vertex_no);




	void computeEdgeLengthForce(double* vertex_0, double* vertex_1, double stiffness,
		double* potential_0, double* potential_1, double rest_length);
	void computeARAPForce(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
		double stiffness, Matrix<double, 3, 4>& A, double volume, Matrix<double, 3, 4>& force);

	void solveNewtonCD_ARAP(std::array<double, 3>* vertex_position, double stiffness, double dt,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
		double* volume, unsigned int vertex_index, std::array<double, 3>* sn);


	void solveCD_ARAP_block(MatrixXd& Hessian, VectorXd& grad, std::array<double, 3>* vertex_position, double stiffness,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices, std::array<int, 4>* indices,
		double* volume, unsigned int tet_index, unsigned int* common_vertex_in_order, int* tet_vertex_index,
		int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num);

	bool solveCertainHessianForNeighborTet(std::array<double, 3>* vertex_position, double stiffness,
		Matrix<double, 3, 4>& A, unsigned int*& common_vertex_in_order, int* neighbor_tet_vetex_indices, 
		MatrixXd& sys_matrix, double volume, VectorXd& grad);

	bool solveTetCertainVertices(std::array<double, 3>* vertex_position, double stiffness,
		Matrix<double, 3, 4>& A, int* vertex_in_sys, int* tet_vetex_indices,
		MatrixXd& sys_matrix, double volume, VectorXd& grad);

	bool computeBarrierVTGradientHessian(MatrixXd& Hessian, VectorXd& grad, double* p, double* t0,
		double* t1, double* t2, double d_hat_2, int* vertex_in_pair, double  stiffness);


	void computeVTBarrierGradientHessian(MatrixXd& Hessian_, VectorXd& grad_, double* p, double* t0, double* t1, double* t2,
		double d_hat_2, int* triangle_vertex_order_in_system, double stiffness);

	void computeEEBarrierGradientHessian(double* ea0, double* ea1, double* eb0, double* eb1, MatrixXd& Hessian_, VectorXd& grad_,
		int* vertex_order_in_system, double stiffness, double d_hat_2, double rest_length_0, double rest_length_1);

	bool computeBarrierEEGradientHessian(double* ea0, double* ea1, double* eb0, double* eb1, MatrixXd& h, VectorXd& g,
		int* vertex_in_pair, double stiffness, double d_hat_2, double rest_length_0, double rest_length_1);


	void computeVTBarrierGradientHessianTest(MatrixXd& Hessian_, VectorXd& grad_, double* p, double* t0, double* t1, double* t2,
		double d_hat_2, int* triangle_vertex_order_in_system, double stiffness, double& barrier);

	void computeEEBarrierGradientHessianTest(double* ea0, double* ea1, double* eb0, double* eb1, MatrixXd& Hessian_, VectorXd& grad_,
		int* vertex_order_in_system, double stiffness, double d_hat_2, double rest_length_0, double rest_length_1, double& barrier_);

private:
	double epsilon_for_bending= 1e-10;

	void solveEdgeLengthConstraint(Vector3f& p1, Vector3f& p2, const double d, double mass_0,
		double mass_1, Vector3f& ori_p1, Vector3f& ori_p2, bool v0_fixed, bool v1_fixed, double& lambda);

	void setTetHessianFromBarrierHessian(MatrixXd& Hessian_system, VectorXd& grad_system, MatrixXd& Hessian_, VectorXd& grad_,
		int* triangle_vertex_order_in_system, int* vertex_in_pair, int vertex_in_use);

	void 	setBarrierGHWithMollifier(double barrier_, MatrixXd& dis_h, VectorXd& dis_g,
		double* ea0, double* ea1, double* eb0, double* eb1,double eps_x, 
		double ee_cross_norm_2, double mollifier, double b_grad, double b_hessian);


	void setFourVertexHessianFromBarrierHessian(MatrixXd& Hessian_system, VectorXd& grad_system, MatrixXd& Hessian_, VectorXd& grad_,
		int* vertex_in_pair, int vertex_in_use);

	void setTetHessianFromHessian(MatrixXd& Hessian_system, VectorXd& grad_system, MatrixXd& Hessian_, VectorXd& grad_,
		int* vertex_order_in_system);


};


#endif // !SECOND_ORDER


