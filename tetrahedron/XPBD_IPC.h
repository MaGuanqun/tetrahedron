#pragma once
#include"./external/Eigen/Dense"
#include"./basic/eigenDenseOperation.h"
#include"./external/Eigen/Sparse"
#include"./external/Eigen/SparseCholesky"
//#include"basic/EigenMatrixIO.h"
#include"./thread.h"
#include"./basic/global.h"
#include"./object/cloth.h"
#include"./object/tetrahedron.h"
#include"./object/collider.h"
#include"./collision/collision.h"
#include"./XPBD/XPBD_constraint.h"
#include"./basic/move_model.h"
#include"./basic/save_scene.h"
#include"./XPBD/second_order.h"
#include"./compute_energy.h"

using namespace Eigen;
using namespace denseOperation;

class XPBD_IPC
{
public:
	XPBD_IPC();
	double time_step;
	double gravity_;
	unsigned int sub_step_num;

	void setForXPBD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
		Thread* thread, double* tolerance_ratio);
	size_t* time_stamp;
	void setPosPredict(int thread_No);
	void computeVelocity(int thread_No);
	unsigned int iteration_number;
	unsigned int inner_iteration_number;
	unsigned int outer_itr_num;
	void initial();
	void reset();
	void resetExternalForce();
	void initialDHatTolerance(double ave_edge_length);
	void updateTetrahedronAnchorVertices();
	void addExternalForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int obj_No);
	void updateItrInfo(int* iteration_num);

	void XPBD_IPCSolve();

	void XPBD_IPC_Block_Solve();

	void XPBD_IPC_Position_Solve();//solve collision as four position constraint
	Collision collision;

	unsigned int* time_indicate_for_simu;

	MoveModel* move_model;
	bool* control_parameter;

	void saveScene(double* force_direction, int obj_No, bool have_force);
	void readScene(const char* file_name, double* force_direction, int& obj_No);
	unsigned int max_iteration_number;
	double velocity_damp;
	unsigned int* sub_step_per_detection;

	bool* has_force;
	void computeCollisionFreePosition(int thread_No);


private:
	void coordinateDescent();
	Floor* floor;
	double gravity[3];

	unsigned int total_thread_num;

	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;
	double sub_time_step;
	std::vector<std::vector<std::array<double, 3>>> f_ext;
	std::vector<std::vector<std::array<double, 3>>> velocity;

	std::vector<std::vector<std::array<double, 3>>> sn;

	//std::vector<std::vector<std::array<double, 3>>> total_gravity;

	std::vector<std::array<double, 3>*> vertex_position_collider;

	std::vector<std::array<double, 3>*> vertex_position;
	std::vector<std::array<double, 3>*> initial_vertex_position;
	std::vector<MeshStruct*> mesh_struct;
	std::vector<MeshStruct*> collider_mesh_struct;
	std::vector<unsigned int*> vertex_index_begin_per_thread;
	void reorganzieDataOfObjects();
	unsigned int total_obj_num;
	void initialVariable();
	XPBDconstraint XPBD_constraint;
	std::vector<std::vector<double>> lbo_weight;
	std::vector<std::vector<VectorXd>>vertex_lbo;
	std::vector<std::vector<double>> rest_mean_curvature_norm;
	//std::vector<std::vector<Vector3d>> rest_Aq;

	std::vector<std::array<int, 3>*> triangle_indices;
	std::vector<std::array<int, 3>*> triangle_indices_collider;


	std::vector<double> lambda;

	std::vector<double> lambda_collision;
	std::vector<unsigned int>constraint_index_start;
	std::vector<std::vector<unsigned int>>collision_constraint_index_start;

	std::vector<unsigned int*> edge_vertices;
	std::vector<unsigned int*> collider_edge_vertices;

	void initialClothBending();
	//void solveBendingConstraint();
	//void solveEdgeLengthConstraint();
	//void solveConstraint(bool need_detection);
	void setConstraintIndex();

	void updatePosition();
	double damping_coe;

	void updateRenderNormal();
	void initialCollisionConstriantNum();
	bool perform_collision;

	void updateNormal();
	void updateRenderVertexNormal();

	bool convergeCondition(unsigned int iteration_num);
	bool innerConvergeCondition(unsigned int iteration_num);

	std::vector<std::array<double, 3>*>address_of_record_vertex_position;

	std::vector<std::vector<std::array<double, 3>>> record_vertex_position;
	std::vector<std::vector<std::array<double, 3>>> record_gloabl_CCD_vertex_position;
	//std::vector<std::vector<std::array<double, 3>>> record_outer_vertex_position;
	//void recordLastStepVertexPosition();

	//std::vector<std::vector<unsigned int>* >unfixed_vertex;
	double max_move_standard;//the max displacement to stop iteration
	double max_move_standard_inner_itr;


	double converge_condition_ratio;// converge_condition_ratio* edge length

	double calEdgeLength();
	void setConvergeCondition();


	unsigned int outer_max_iteration_number;

	bool use_bending_based_on_vertex = true;


	double energy;
	double previous_energy = 1e-15;
	//std::vector<double>energy_per_thread;

	double energy_converge_ratio;

	SaveScene save_scene;

	SecondOrderConstraint second_order_constraint;
	void updateSn();

	void computeCurrentEnergy();
	double computeInertialEnergy();
	double computeCurrentARAPEnergy();

	std::vector<std::vector<Vector3d>>residual;

	ComputeEnergy compute_energy;
	void newtonCDTet();
	void firstNewtonCD();
	void solveNewtonCD_tet(std::array<double, 3>* vertex_position, double stiffness, double dt,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
		double* volume, unsigned int vertex_index, std::array<double, 3>* sn, double* lambda);

	void solveNewtonCDTetWithCollision(std::array<double, 3>* vertex_position, 
		double* record_vertex_position_,
		double* last_step_vertex_position, 
		double ARAP_stiffness, double dt,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* tet_vertex_indices, double* mass,
		double* volume, unsigned int vertex_index, std::array<double, 3>* sn, double collision_stiffness, unsigned int obj_No,
		bool vertex_on_surface, unsigned int vertex_index_on_surface);
	void getARAPHessian(Matrix3d& Hessian, Vector3d& grad, std::array<double, 3>* vertex_position, double stiffness,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, 
		double* volume, unsigned int vertex_index, unsigned int obj_No);
	void getVTCollisionHessain(Matrix3d& Hessian, Vector3d& grad, double* vertex_position_, double stiffness,
		 unsigned int* VT, unsigned int num, double* ori_volume, unsigned int obj_No, unsigned int vertex_index);
	void getTVCollisionHessain(Matrix3d& Hessian, Vector3d& grad,
		double* pos_0, double* pos_1, double* pos_2,
		unsigned int vertex_no, double stiffness, unsigned int* TV, unsigned int num, double* ori_volume);

	void getEECollisionHessian(Matrix3d& Hessian, Vector3d& grad, double* pos0, double* pos1, unsigned int* EE, unsigned int num,
		double* ori_volume, double stiffness, unsigned int obj_index, unsigned int edge_index, unsigned int vertex_no);
	void getVT_ColiderCollisionHessain(Matrix3d& Hessian, Vector3d& grad, double* vertex_position_, double stiffness,
		unsigned int* VT, unsigned int num, double* ori_volume);
	void newtonCDTetWithCollision();
	void updateCollisionFreePosition();
	void getCollisionHessian(Matrix3d& Hessian, Vector3d& grad, std::array<double, 3>* vertex_position, 
		double* last_step_vertex_position,
		double collision_stiffness, unsigned int obj_No,
		unsigned int vertex_index, unsigned int vertex_index_on_surface);
	bool getFloorHessian(double& Hessian, double& grad, double* vertex_position, double floor_value,
		double* last_step_position, unsigned int dimension, double collision_stiffness, bool direction, double tolerance);
	unsigned int min_inner_iteration, min_outer_iteration;





	bool nearly_not_move;//to indicate if the move distance of current itr is far than requirement

	void firstOnlyInertialCollision();
	void solveInertialCollision(std::array<double, 3>* vertex_position,
		double* record_vertex_position_,
		double* last_step_vertex_position,
		double dt, double* mass,
		double* volume, unsigned int vertex_index, std::array<double, 3>* sn, double collision_stiffness, unsigned int obj_No,
		bool vertex_on_surface, unsigned int vertex_index_on_surface);

	void newtonCDTetBlock();

	void solveNewtonCD_tetBlock(std::array<double, 3>* vertex_position, double stiffness, double dt, double* mass,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices, std::array<int, 4>* indices,
		double* volume, unsigned int tet_index, std::array<double, 3>* sn, unsigned int* common_vertex_in_order,
		int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num,
		std::vector<unsigned int>* triangle_of_a_tet, 	std::vector<unsigned int>* edge_of_a_tet, double collision_stiffness,
		unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
		int* vertex_index_on_surface, std::array<double, 3>* record_ori_pos);

	std::vector<double>record_energy;



	void checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, int* element_indices,
		int obj_No,int* triangle_vertex_order_in_system);

	void checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, int* element_indices, 
		int* triangle_vertex_order_in_system);
	void checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, unsigned int* element_indices,
		int* triangle_vertex_order_in_system);

	void getVTCollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
		unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int obj_No,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* VT_collider, int num_collider);


	//void getVTCollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
	//	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int obj_No,
	//	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2);





	void getTVCollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, double* t0, double* t1, double* t2, double stiffness,
		unsigned int* TV, int num, unsigned int obj_No, int* triangle_indices,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, 
		unsigned int* TV_collider, int collider_num);

	bool vertexInTet(int unfixed_tet_vertex_num, int vertex_No, int* tet_unfixed_vertex_indices);

	void getEECollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_No, unsigned int* edge_vertex_index,
		unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2,
		int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, double* ea0, double* ea1, double stiffness,
		unsigned int* EE_collider, int num_collider);

	bool edgeInSameTetDuplicate(int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet,
		unsigned int compare_edge_index);


	void getCollisionHessian(MatrixXd& Hessian, VectorXd& grad,std::vector<unsigned int>* triangle_of_a_tet,
		std::vector<unsigned int>* edge_of_a_tet, 
		double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
		int unfixed_tet_vertex_num, double d_hat_2,
		int* vertex_index_on_surface, std::array<double, 3>* vertex_position);


	void getFloorHessianForTet(MatrixXd& Hessian, VectorXd& grad, double* vertex_position, double floor_value,
		unsigned int dimension, double collision_stiffness, bool direction, double d_hat, unsigned int vertex_order_in_matrix, unsigned int unfixed_vertex_num);

	double getCollisionTime(std::vector<unsigned int>* triangle_of_a_tet,
		std::vector<unsigned int>* edge_of_a_tet,
		unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
		int unfixed_tet_vertex_num,
		int* vertex_index_on_surface, std::array<double, 3>* current_vertex_position,
		std::array<double, 3>* initial_vertex_position);

	//void solveNewtonCD_collisionBlock(unsigned int vertex_obj_no, unsigned int vertex_index, 
	//	unsigned int triangle_obj_No, unsigned int triangle_index);

	//void getCollisionBlockHessian(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangle_of_first_primitve,
	//	std::vector<unsigned int>* triangle_of_second_primitve,
	//	std::vector<unsigned int>* edge_of_first_primitve,
	//	std::vector<unsigned int>* edge_of_second_primitve,
	//	double collision_stiffness, int* pair_actual_unfixed_vertex_indices,
	//	int unfixed_vertex_num, double d_hat_2, int* vertex_index_on_surface);

	bool has_collider;
};

