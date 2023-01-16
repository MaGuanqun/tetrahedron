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

#include"collision/collision_compare.h"

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

	//CollisionCompare collision_compare;

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
	void computeCollisionFreePositionFromRecord(int thread_No);



	double energy_converge_ratio;
	unsigned int min_inner_iteration, min_outer_iteration;
	double max_move_standard;//the max displacement to stop iteration
	unsigned int outer_max_iteration_number;
	double energy_converge_standard;

	//void computeTetHessianInAdvance(int thread_No, int color_No);

	void XPBD_IPC_Block_Solve_Multithread();

	void newtonCDTetBlockAGroup(int thread_No, int color);
	void newtonCDTetBlockAGroupCollision(int thread_No, int color);
	
	void warmStart();

	void tetHessian();

	void inversionTest(int thread_No);
	void tetGradForColor(int thread_No, unsigned int color_No);
	//void newtonCDTetBlockAGroupTest(int thread_No, int color);

	//void tetGradForColorCollision(unsigned int color_No);
	void tetGradForColorCollisionNeighbor(int thread_No, unsigned int color_No);
	//void tetGradForColorCollisionEECollider(int thread_No, unsigned int color_No);
	//void tetGradForColorCollisionTVCollider(int thread_No, unsigned int color_No);
	//void tetGradForColorCollisionVT(int thread_No, unsigned int color_No);
	//void tetGradForColorCollisionEE(int thread_No, unsigned int color_No);
	//void tetGradForColorCollisionTV(int thread_No, unsigned int color_No);

	std::vector<std::array<double, 3>>vertex_trace;
	void updatePositionAverage(int thread_No);
	void lastColorVertexBelongToGroup(int thread_No);
	//std::vector<std::array<double, 3>> e0_0;
	//std::vector<std::array<double, 3>> e0_1;
	//std::vector<std::array<double, 3>> e1_0;
	//std::vector<std::array<double, 3>> e1_1;


	//std::vector<std::array<double, 3>> e0_0_;
	//std::vector<std::array<double, 3>> e0_1_;
	//std::vector<std::array<double, 3>> e1_0_;
	//std::vector<std::array<double, 3>> e1_1_;
	void computePreviousColorARAPEnergy(int thread_No, unsigned int color_No);

	void inertialEnergyPerThread(int thread_No);
	void computeARAPEnergyPerThread(int thread_No);

	double min_collision_time;
	void computeBarrierEnergy(int thread_No);


	void computePreviousColorCollisionEnergy(int thread_No);

	void computePreviousColorInertialEnergy(int thread_No, unsigned int color_No);

	void maxDisplacement(int thread_No, unsigned int color_No);

	void  computeColorInertialEnergy(int thread_No);

	void computeLastColorARAPEnergy(int thread_No);

	void sumAllGrad(int thread_No);
	void sumWithInertial(int thread_No);

	void solveBlockForWarmStart(int thread_No);

	void computePreviousColorInversion(int thread_No, int color_No);
	void computeLastColorInversion(int thread_No);


	//void testMulti(int thread_No);

	int counter = 0;

	void tetGrad(int thread_No);

private:

	void computeGradient();

	std::vector<std::vector<std::vector<unsigned int>>*> vertex_index_of_a_tet_color_per_thread_start_group;

	double computeInertialEnergy();
	double computeInertialEnergyWarmStart();
	double computeCurrentARAPEnergy();
	double computeBarrierEnergy();

	double computeFloorEnergy(int type, double collision_stiffness, std::vector<unsigned int>* record_vertex_collide_with_floor, int start, int end);

	void allPairCollisionInversionTime();
	void allPairCollisionTimeWarmStart();

	void previousColorCollisionInversionTime(unsigned int color_No);

	void computeCurrentEnergy();

	double computeWarmStartEnergy();


	double computeLastColorARAPEnergy();

	void computeInversionForWarmStart();



	

	double computeFloorEnergy(int type);


	double computeEEEnergy(std::vector<std::vector<unsigned int>>* record_pair, unsigned int** edge_v_0, unsigned int** edge_v_1,
		std::array<double, 3>** e0_current_pos, std::array<double, 3>** e1_current_pos, double collision_stiffness, std::vector<std::vector<double>>* d_hat, 
		bool** belong_to_color_group, int type);


	double computeEEEnergy(std::vector<unsigned int>* record_pair, unsigned int** edge_v_0, unsigned int** edge_v_1,
		std::array<double, 3>** e0_current_pos, std::array<double, 3>** e1_current_pos, double collision_stiffness, double* d_hat,
		bool** belong_to_color_group, int type, unsigned int start, unsigned int end);


	double computeVTEnergy(std::vector<std::vector<unsigned int>>* record_vt_pair, std::array<int, 3>** triangle_indices,
		std::array<double, 3>** v_current_pos, std::array<double, 3>** t_current_pos, double collision_stiffness, std::vector<std::vector<double>>* d_hat, 
		bool** belong_to_color_group, int type);


	double computeVTEnergy(std::vector<unsigned int>* record_vt_pair, std::array<int, 3>** triangle_indices,
		std::array<double, 3>** v_current_pos, std::array<double, 3>** t_current_pos, double collision_stiffness, double* d_hat,
		bool** belong_to_color_group, int type, unsigned int start, unsigned int end);


	double computePreviousColorEnergy(unsigned int color_no);

	double computeLastColorEnergy();
	double computeColorInertialEnergy();


	double computePreviousColorInertialEnergy(unsigned int color_No);

	double computePreviousColorCollisionEnergy();


	void solveVT_BlockPerThread(std::array<double, 3>** record_vertex_position, int** record_vertex_num, unsigned int* pair, unsigned int start, unsigned int end,
		char* vt_hessian_record_index, bool only_solve_collision_pair, int thread_No);

	void solvecollider_BlockPerThread(std::array<double, 3>** record_vertex_position, int** record_vertex_num, std::vector<unsigned int>* pair,
		bool only_solve_collision_pair, int type, int start, int end, int thread_No);



	void solveEE_BlockPerThread(std::array<double, 3>** record_vertex_position, int** record_vertex_num, unsigned int* pair, unsigned int start, unsigned int end,
		char* ee_hessian_record_index_exist, bool only_solve_collision_pair, int thread_No);

	void initialHessianMap();

	void 	initialRecordHessian();

	std::unordered_map<std::array<unsigned int, 2>, StoreHessianWithOrderInConstraint, pair_hash> common_hessian;//only include collision

	//std::unordered_map<std::array<unsigned int, 2>, std::array<double, 9>, pair_hash> common_hessian_compare;//only include collision

	std::vector<std::unordered_map<std::array<int, 2>,double, pair_hash>> tet_hessian;//size is total obj num
	std::vector<double> floor_hessian;

	//std::unordered_map<unsigned int, double> compare_floor_hessian;

	std::vector<std::vector<double>>common_grad;
	//std::vector<double>common_grad_compare;


	std::vector<std::vector<std::vector<std::vector<unsigned int>>>*>tet_color_groups;//size is tet num
	//std::vector<std::vector<std::vector<char>>*>tet_color_groups_label;

	std::vector<std::vector<unsigned int>*>tet_around_tet; //size is tet num


	void setColorNum();

	void solveNewtonCD_tetBlock();

	int max_tet_color_num;//the max number of different objects tet color

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
	std::vector<std::array<double, 3>*> vertex_position_render;
	std::vector<MeshStruct*> mesh_struct;
	std::vector<MeshStruct*> collider_mesh_struct;
	std::vector<unsigned int*> vertex_index_begin_per_thread;

	std::vector<unsigned int*>tet_index_begin_per_thread;//size is total object num

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

	std::vector<std::array<int, 4>*> tet_indices; //size is total obj num

	std::vector<Matrix<double, 3, 4>*>tet_A;
	std::vector<double*>tet_volume;

	std::vector<double> lambda;

	std::vector<double> lambda_collision;
	std::vector<unsigned int>constraint_index_start;
	std::vector<std::vector<unsigned int>>collision_constraint_index_start;

	std::vector<unsigned int*> edge_vertices;
	std::vector<unsigned int*> collider_edge_vertices;

	std::vector<std::vector<bool>*> is_vertex_fixed;

	std::vector<int*>vertex_index_surface;

	std::vector<double*>mass;



	std::vector<std::vector<unsigned int>*>triangle_around_triangle;
	std::vector<std::vector<unsigned int>*>edge_around_triangle;

	std::vector<std::vector<unsigned int>*>tet_around_vertex;
	std::vector<std::vector<unsigned int>*>tet_around_triangle;


	std::vector<std::vector<unsigned int>*>triangle_around_edge;
	std::vector<std::vector<unsigned int>*>edge_around_edge;
	std::vector<std::vector<unsigned int>*>tet_around_edge;

	std::vector<MeshStruct::Vertex*> vertices;
	
	std::vector<unsigned int*>vertex_surface_to_global;


	std::vector<double*> rest_edge_length;
	std::vector<double*> rest_edge_length_collider;


	std::vector<std::vector<std::array<double, 3>*>> record_vertex_by_thread;//record_vertex_position_every_thread in thread->no, obj_no, vertex_No
	std::vector<std::vector<int*>> record_vertex_update_num_by_thread;//record_vertex_position_num_every_thread in thread->no, obj_no, vertex_No

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
	bool convCondition(unsigned int iteration_num, unsigned int min_itr, double energy, double previous_energy, unsigned int max_itr, double energy_converge_standard,
		double energy_converge_ratio);

	//std::vector<std::array<double, 3>*>address_of_record_vertex_position;

	std::vector<std::vector<std::array<double, 3>>> record_vertex_position;
	std::vector<std::vector<std::array<double, 3>>> record_collision_free_vertex_position;


	std::vector<std::vector<std::vector<std::array<double, 3>>>> record_vertex_position_every_thread; // obj_no -> thread_No ->vertex_No
	std::vector<std::vector<std::vector<int>>> record_vertex_position_num_every_thread;// obj_no -> thread_No ->vertex_No




	std::vector<std::array<double, 3>*> record_collision_free_vertex_position_address;
	std::vector<std::array<double, 3>*> record_vertex_position_address;

	//std::vector<std::vector<std::array<double, 3>>> record_outer_vertex_position;
	//void recordLastStepVertexPosition();

	//std::vector<std::vector<unsigned int>* >unfixed_vertex;

	double max_move_standard_inner_itr;



	double calEdgeLength();
	//std::vector<double> store_tet_arap_hessian; //for every 12*12, we only store 4*4 as every block is a diagonal matrix 
	//std::vector<double> store_tet_arap_grad;
	std::atomic_flag* is_tet_arap_grad_compute;

	int max_tet_size_of_a_color_group;

	std::vector<unsigned int> prefix_sum_of_every_tet_index; //size is tet_num+1

	std::vector<int> vertex_num_on_surface_prefix_sum;



	std::vector<std::array<int, 4>*> unfix_tet_index;
	std::vector<unsigned int*> unfixed_tet_vertex_num;

	bool use_bending_based_on_vertex = true;


	double energy;
	double previous_energy = 1e-15;
	//std::vector<double>energy_per_thread;



	SaveScene save_scene;

	SecondOrderConstraint second_order_constraint;
	void updateSn();

	double computeBlockCurrentEnergy(std::array<double, 3>* vertex_position, double stiffness, double dt,
		double* mass,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices,
		double* volume, unsigned int tet_index, std::array<double, 3>* sn, 
		int* tet_vertex_index, unsigned int unfixed_vertex_num,
		double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
		std::vector<unsigned int>* triangle_of_a_tet,
		std::vector<unsigned int>* edge_of_a_tet, int* vertex_index_on_surface,
		double* mass_inv_,  std::array<int, 4>*tet_vertex_indices,
		char* indicate_vertex_collide_with_floor, double* record_vertex_collide_with_floor_d_hat);

	double computeVTCollisionEnergyPerElement(double* pos0, unsigned int* triangle_index, std::array<double, 3>** pos_t, std::array<int, 3>** triangle_indices,
		unsigned int num, double* d_hat, double collision_stiffness);





	double computeEECollisionEnergyPerElement(double* pos0, double* pos1, unsigned int* edge_index,
		unsigned int num, double* d_hat, double collision_stiffness,
		int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, std::array<double, 3>** pos_e, unsigned int obj_No,
		unsigned int** edge_1_vertex);

	double computeTVColliderCollisionEnergyPerElement(double* pos0, double* pos1, double* pos2, unsigned int* vertex_index_,
		unsigned int num, double* d_hat, double collision_stiffness, std::array<double, 3>** pos_v);

	double computeEEColliderCollisionEnergyPerElement(double* pos0, double* pos1, unsigned int* edge_index,
		unsigned int num, double* d_hat, double collision_stiffness, std::array<double, 3>** pos_e, unsigned int** edge_1_vertex);

	double computeTVCollisionEnergyPerElement(double* pos0, double* pos1, double* pos2, unsigned int* vertex_index_,
		unsigned int num, double* d_hat, double collision_stiffness,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, std::array<double, 3>** pos_v, unsigned int obj_No);

	double computeVTCollisionEnergy(unsigned int** vertex_triangle_pair_by_vertex_, unsigned int** vertex_triangle_pair_num_record_,
		std::array<double, 3>** triangle_position, std::array<double, 3>** vertex_position, unsigned int close_pair_num, bool is_TV,
		std::array<int, 3>** triangle_vertex);

	double computeEECollisionEnergy(unsigned int** edge_edge_pair_by_vertex_, unsigned int** edge_edge_pair_num_record_,
		std::array<double, 3>** edge_0_position, std::array<double, 3>** edge_1_position, unsigned int close_pair_num,
		unsigned int** edge_vertex, unsigned int** edge_1_vertex, bool is_self);

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
	void recordInitialPosition();
	void getCollisionHessian(Matrix3d& Hessian, Vector3d& grad, std::array<double, 3>* vertex_position, 
		double* last_step_vertex_position,
		double collision_stiffness, unsigned int obj_No,
		unsigned int vertex_index, unsigned int vertex_index_on_surface);
	bool getFloorHessian(double& Hessian, double& grad, double* vertex_position, double floor_value,
		double* last_step_position, unsigned int dimension, double collision_stiffness, bool direction, double tolerance);






	bool nearly_not_move;//to indicate if the move distance of current itr is far than requirement

	//void firstOnlyInertialCollision();
	void solveInertialCollision(std::array<double, 3>* vertex_position,
		double* record_vertex_position_,
		double* last_step_vertex_position,
		double dt, double* mass,
		double* volume, unsigned int vertex_index, std::array<double, 3>* sn, double collision_stiffness, unsigned int obj_No,
		bool vertex_on_surface, unsigned int vertex_index_on_surface);

	void newtonCDTetBlock();
	//void newtonCDTetBlockTest(int color_No);
	void newtonVTCollisionBlock();
	void newtonEECollisionBlock();

	void newtonCDBlock();

	void newtonEEColliderCollisionBlock();
	void newtonVTColliderCollisionBlock();
	void newtonTVColliderCollisionBlock();



	void solveNewtonCD_tetBlock(std::array<double, 3>* vertex_position, double stiffness, double dt, double* mass,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices, std::array<int, 4>* indices,
		double* volume, unsigned int tet_index, std::array<double, 3>* sn, unsigned int* common_vertex_in_order,
		int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num,
		std::vector<unsigned int>* triangle_of_a_tet, 	std::vector<unsigned int>* edge_of_a_tet, double collision_stiffness,
		unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
		int* vertex_index_on_surface, std::array<double, 3>* record_ori_pos, double* hessian_record, char* indicate_collide_with_floor);

	std::vector<double>record_energy;

	void solveTetBlock(std::array<double, 3>* vertex_position, double stiffness, double dt,
		double* mass, double* mass_inv,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices,
		double* volume, unsigned int tet_index, std::array<double, 3>* sn, unsigned int* common_vertex_in_order,
		int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num,
		double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
		std::unordered_map<std::array<int, 2>, double, pair_hash>& tet_hessian,
		std::unordered_map<std::array<unsigned int, 2>, StoreHessianWithOrderInConstraint, pair_hash>& collision_hessian,
		double* common_grad, std::vector<unsigned int>* triangle_of_a_tet,
		std::vector<unsigned int>* edge_of_a_tet, int* vertex_index_on_surface, unsigned int prefix_sum_vetex_obj,
		double* floor_map, std::array<double, 3>* record_ori_pos, char* indicate_collide_with_floor, int color_No);

	void solveTetBlockCollision(std::array<double, 3>* vertex_position, double stiffness, double dt, std::array<int, 4>* indices,
		double* mass,
		Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices,
		double* volume, unsigned int tet_index, std::array<double, 3>* sn, unsigned int* common_vertex_in_order,
		int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num, std::vector<unsigned int>* triangle_of_a_tet,
		std::vector<unsigned int>* edge_of_a_tet, double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
		int* vertex_index_on_surface, std::unordered_map<std::array<int, 2>, double, pair_hash>& tet_hessian,
		std::unordered_map<std::array<unsigned int, 2>, StoreHessianWithOrderInConstraint, pair_hash>& collision_hessian,
		double* common_grad, std::array<double, 3>* record_vertex_position,
		int* record_vertex_num, unsigned int prefix_sum_vetex_obj, double* floor_map, double& max_dis,
		std::array<double, 3>* collision_free_pos);



	void checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, int* element_indices,
		int obj_No,int* triangle_vertex_order_in_system, int size_num);

	void getARAPCollisionHessianForPair(MatrixXd& Hessian, VectorXd& grad, double stiffness, int tet_obj, int tet_index, int* tet_vertex, int* tet_unfixed_vertex_indices,
		int unfixed_tet_vertex_num,  std::array<double, 3>* vertex_position, Matrix<double, 3, 4>* A, double* volume);

	void checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, int* element_indices, 
		int* triangle_vertex_order_in_system, int size_num);
	void checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, unsigned int* element_indices,
		int* triangle_vertex_order_in_system);
	void checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, unsigned int* element_indices, int obj_No,
		int* triangle_vertex_order_in_system);
	void getVTCollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
		unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int obj_No,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* VT_collider, int num_collider);
	void getVTCollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
		unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int vertex_obj_No, unsigned int tri_obj_No,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* VT_collider, int num_collider);


	void getTVCollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, double* t0, double* t1, double* t2, double stiffness,
		unsigned int* TV, int num, unsigned int tri_obj, unsigned int obj_No_0, unsigned int obj_No_1, int* triangle_indices,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* TV_collider, int collider_num);


	void getTVCollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, double* t0, double* t1, double* t2, double stiffness,
		unsigned int* TV, int num, unsigned int obj_No, int* triangle_indices,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, 
		unsigned int* TV_collider, int collider_num);

	bool vertexInPair(int unfixed_tet_vertex_num, int vertex_No, int* tet_unfixed_vertex_indices, int obj_No);
	bool vertexInTet(int unfixed_tet_vertex_num, int vertex_No, int* tet_unfixed_vertex_indices);

	void getEECollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_No, unsigned int* edge_vertex_index,
		unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2,
		int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, double* ea0, double* ea1, double stiffness,
		unsigned int* EE_collider, int num_collider, double edge_length_0);


	void getEECollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, unsigned int ee_obj_No, unsigned int* edge_vertex_index,
		unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2,
		unsigned int obj_No_0, unsigned int obj_No_1,
		int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, double* ea0, double* ea1, double stiffness,
		unsigned int* EE_collider, int num_collider, double edge_length_0);

	bool edgeInSameTetDuplicate(int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet,
		unsigned int compare_edge_index);

	bool edgeInSameTetDuplicate(int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet,
		unsigned int compare_edge_index, unsigned int compare_edge_obj);

	void getCollisionHessian(MatrixXd& Hessian, VectorXd& grad,std::vector<unsigned int>* triangle_of_a_tet,
		std::vector<unsigned int>* edge_of_a_tet, 
		double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
		int unfixed_tet_vertex_num, double d_hat_2,
		int* vertex_index_on_surface, std::array<double, 3>* vertex_position);


	void getFloorHessianForTet(MatrixXd& Hessian, VectorXd& grad, double* vertex_position, double floor_value,
		unsigned int dimension, double collision_stiffness, bool direction, double d_hat, unsigned int vertex_order_in_matrix, unsigned int unfixed_vertex_num);


	void getFloorHessian(MatrixXd& Hessian, VectorXd& grad, int* tet_actual_unfixed_vertex_indices,
		int unfixed_tet_vertex_num, double d_hat_2,
		int* vertex_index_on_surface, std::array<double, 3>* vertex_position, double collision_stiffness);


	double getInversionTime(unsigned int tet_index, std::vector<unsigned int>* neighbor_tet_index,
		std::array<int, 4>* tet_vertex, std::array<double, 3>* current_vertex_position,
		std::array<double, 3>* initial_vertex_position);


	double getCollisionTime(std::vector<unsigned int>* triangle_of_a_tet,
		std::vector<unsigned int>* edge_of_a_tet,
		unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
		int unfixed_tet_vertex_num,
		int* vertex_index_on_surface, std::array<double, 3>* current_vertex_position,
		std::array<double, 3>* initial_vertex_position, char* indicate_collide_with_floor);

	double getCollisionTime(std::vector<unsigned int>* triangle_of_a_tet,
		std::vector<unsigned int>* edge_of_a_tet,
		int* tet_actual_unfixed_vertex_indices,
		int unfixed_tet_vertex_num,
		int** vertex_index_on_surface,
		std::array<double, 3>** current_vertex_position,
		std::array<double, 3>** initial_vertex_position);



	void getCollisionPairHessian(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_0,  unsigned int obj_1, 
		double collision_stiffness, std::vector<unsigned int>* triangle_around_0, std::vector<unsigned int>* triangle_around_1,
		std::vector<unsigned int>* edge_around_0, std::vector<unsigned int>* edge_around_1, 
		std::vector<unsigned int>* tet_around_0, std::vector<unsigned int>* tet_around_1, double d_hat_2,
		int* unfixed_pair_vertex_index, int unfixed_num, std::vector<unsigned int>* around_triangle, std::vector<unsigned int>* around_edge,
		std::vector<unsigned int>*around_tet, double arap_stiffness, bool obj_0_collider, bool obj_1_collider);


	
	void solveVT_collisionBlock(unsigned int vertex_obj_no, unsigned int vertex_index, unsigned int triangle_obj_No, unsigned int triangle_index,
		double stiffness, double dt, double collision_stiffne, std::vector<unsigned int>* triangle_around_vertex, std::vector<unsigned int>* triangle_around_triangle,
		std::vector<unsigned int>* edge_around_vertex, std::vector<unsigned int>* edge_around_triangle, 
		std::vector<unsigned int>* tet_around_vertex, std::vector<unsigned int>* tet_around_triangle, double d_hat_2, bool vertex_collider, bool triangle_collider);

	void solveEE_collisionBlock(unsigned int obj_No_0, unsigned int primitive_0_index, unsigned int obj_No_1, unsigned int primitive_1_index,
		double stiffness, double dt, double collision_stiffne, std::vector<unsigned int>* triangle_around_0, std::vector<unsigned int>* triangle_around_1,
		std::vector<unsigned int>* edge_around_0, std::vector<unsigned int>* edge_around_1,
		std::vector<unsigned int>* tet_around_0, std::vector<unsigned int>* tet_around_1,
		double d_hat_2, bool edge_0_collider, bool edge_1_collider);

	void getCollisionBlockCollisionHessian(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangles,
		std::vector<unsigned int>* edges,
		double collision_stiffness, int* pair_actual_unfixed_vertex_indices,
		int unfixed_vertex_num, double d_hat_2, int** vertex_index_on_surface, unsigned int vertex_obj_No, unsigned int tri_obj_No);


	
	void getCollisionBlockTetHessian(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* tets, double stiffness, int* pair_actual_unfixed_vertex_indices, //pair_actual_unfixed_vertex_indices first obj, second primitive index
		int unfixed_vertex_num);

	bool has_collider;

	void comparePrimitiveAroundPrimitveTogether(std::vector<unsigned int>* primitive_around_1, std::vector<unsigned int>* primitive_around_2,
		unsigned int obj_1, unsigned int obj_2, std::vector<unsigned int>* primitive_together);
	bool checkMaxDisplacement();

	bool displacement_satisfied;

	void getVTCollisionHessainForPairTest(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
		unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int vertex_obj_No, unsigned int tri_obj_No,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* VT_collider, int num_collider);

	void getCollisionPairHessianTest(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_0, unsigned int obj_1,
		double collision_stiffness, std::vector<unsigned int>* triangle_around_0, std::vector<unsigned int>* triangle_around_1,
		std::vector<unsigned int>* edge_around_0, std::vector<unsigned int>* edge_around_1,
		std::vector<unsigned int>* tet_around_0, std::vector<unsigned int>* tet_around_1, double d_hat_2,
		int* unfixed_pair_vertex_index, int unfixed_num, std::vector<unsigned int>* around_triangle, std::vector<unsigned int>* around_edge,
		std::vector<unsigned int>* around_tet, double arap_stiffness, bool obj_0_collider, bool obj_1_collider);

	void getCollisionBlockCollisionHessianTest(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangles,
		std::vector<unsigned int>* edges,
		double collision_stiffness, int* pair_actual_unfixed_vertex_indices,
		int unfixed_vertex_num, double d_hat_2, int** vertex_index_on_surface, unsigned int vertex_obj_No, unsigned int tri_obj_No);
	void getTVCollisionHessainForPairTest(MatrixXd& Hessian, VectorXd& grad, double* t0, double* t1, double* t2, double stiffness,
		unsigned int* TV, int num, unsigned int tri_obj, unsigned int obj_No_0, unsigned int obj_No_1, int* triangle_indices,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* TV_collider, int collider_num);

	void getEECollisionHessainForPairTest(MatrixXd& Hessian, VectorXd& grad, unsigned int ee_obj_No, unsigned int* edge_vertex_index,
		unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2,
		unsigned int obj_No_0, unsigned int obj_No_1,
		int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, double* ea0, double* ea1, double stiffness,
		unsigned int* EE_collider, int num_collider, double edge_length_0);


	void initalARAPHessianStorages();

	void getVTCollisionHessainForTetFromRecord(MatrixXd& Hessian, VectorXd& grad,
		unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int obj_No, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num,
		int num_collider, double* vt_hessian_record, double* vt_grad_record, int* vt_hessian_record_index,
		double* vt_collider_hessian_record, double* vt_collider_grad_record, unsigned int* vt_collider);


	void setTetHessianFromBarrierHessian(MatrixXd& Hessian_system, double* grad_system, double* Hessian_, double* grad_,
		int* triangle_vertex_order_in_system, int* vertex_in_pair, int vertex_in_use);


	void getTVCollisionHessainForTetFromRecord(MatrixXd& Hessian, VectorXd& grad, 
		unsigned int* TV, int num, unsigned int obj_No, int* triangle_indices,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, unsigned int* TV_collider, int collider_num,
		int* tv_collider_hessian_record_index, double* tv_collider_hessian_record, double* tv_collider_grad_record);


	void getEECollisionHessainForTetFromRecord(MatrixXd& Hessian, VectorXd& grad, 
		unsigned int obj_No, 
		unsigned int edge_index,
		unsigned int* edge_vertex_index,
		unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, 
		int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet,
		unsigned int* EE_collider, int num_collider, 
		double* ee_hessian_record_, double* ee_grad_record_, int* ee_hessian_record_index_,
		int* ee_collider_hessian_record_index, double* ee_collider_hessian_record, double* ee_collider_grad_record);

	void getCollisionHessianFromRecord(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangle_of_a_tet,
		std::vector<unsigned int>* edge_of_a_tet,
		double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
		int unfixed_tet_vertex_num, 
		int* vertex_index_on_surface, unsigned int* vt_prefix_sum, unsigned int* ee_prefix_sum, unsigned int* tv_collider_prefix_sum, unsigned int* ee_collider_prefix_sum, std::array<double, 3>* vertex_position);

	//std::vector<std::vector<std::array<double, 3>>> temp_save_pos;
	//void tempSavePos();
	//void tempRestorePos();

	//void setCollisionPairTetGrad(int tet_obj_No, unsigned int start, unsigned int end, unsigned int* element, std::vector<unsigned int>* tet_around_an_element);

	void setCollisionPairTetNeighborGrad(std::vector<unsigned int>& tet_involved, int start, int end, double* com_grad);

	void initialRecordPositionForThread();

	void solveVT_Block(unsigned int vertex_obj_no, unsigned int vertex_index, unsigned int triangle_obj_No, unsigned int triangle_index,
		double dt, std::array<double, 3>** record_vertex_position,
		int** record_vertex_num, bool only_solve_collision_pair, int thread_No);

	void solveCollisionWithColliderBlock(double dt, unsigned int obj_no, unsigned int index, int type, std::array<double, 3>** record_vertex_position,
		int** record_vertex_num, bool only_solve_collision_pair, int thread_No);


	void getVTCollisionHessainForPairFromRecord(MatrixXd& Hessian, VectorXd& grad,
		unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int vertex_obj_No, unsigned int tri_obj_No,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num,
		double* vt_hessian_record, double* vt_grad_record, int* vt_hessian_record_index,
		double* vt_collider_hessian_record, double* vt_collider_grad_record);

	void getTVCollisionHessainForPairFromRecord(MatrixXd& Hessian, VectorXd& grad,
		unsigned int* TV, int num, unsigned int tri_obj, unsigned int obj_No_0, unsigned int obj_No_1, int* triangle_indices,
		int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, int collider_num,
		int* tv_collider_hessian_record_index, double* tv_collider_hessian_record, double* tv_collider_grad_record);


	void getCollisionBlockCollisionHessianFromRecord(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangle_of_a_tet,
		std::vector<unsigned int>* edge_of_a_tet,
		unsigned int obj_No_0, unsigned int obj_No_1, int* pair_actual_unfixed_vertex_indices,
		int unfixed_vertex_num,MatrixXd& test_m, VectorXd& grad_m);

	void getEECollisionHessainForPairFromRecord(MatrixXd& Hessian, VectorXd& grad, unsigned int ee_obj_No,
		unsigned int edge_index,
		unsigned int* edge_vertex_index,
		unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num,
		unsigned int obj_No_0, unsigned int obj_No_1,
		int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, unsigned int* EE_collider, int num_collider,
		double* ee_hessian_record_, double* ee_grad_record_, int* ee_hessian_record_index_,
		int* ee_collider_hessian_record_index, double* ee_collider_hessian_record, double* ee_collider_grad_record);


	void getHessianForCollisionBlock(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_0, unsigned int obj_1,
		std::vector<unsigned int>* triangle_around_0, std::vector<unsigned int>* triangle_around_1,
		std::vector<unsigned int>* edge_around_0, std::vector<unsigned int>* edge_around_1,
		std::vector<unsigned int>* tet_around_0, std::vector<unsigned int>* tet_around_1,
		int* unfixed_pair_vertex_index, int unfixed_num, std::vector<unsigned int>* around_triangle, std::vector<unsigned int>* around_edge,
		std::vector<unsigned int>* around_tet);

	void getARAPCollisionHessianForPairFromRecord(MatrixXd& Hessian, VectorXd& grad, int tet_obj,int* tet_vertex, int* tet_unfixed_vertex_indices,
		int unfixed_tet_vertex_num, double* Hessian_, double* grad_);

	void setHessianFromTetHessian(MatrixXd& Hessian_system, double* grad_system, double* Hessian_, double* grad_, int* vertex_in_sys);

	void solveEE_Block(unsigned int obj_No_0, unsigned int primitive_0_index, unsigned int obj_No_1, unsigned int primitive_1_index,
		double dt,
		std::array<double, 3>** record_vertex_position,
		int** record_vertex_num, bool only_solve_collision_pair, int thread_No);

	void getEEUnfixedPairIndex(int* unfixed_pair_vertex_index, int& unfixed_num, unsigned int obj_No_0,
		unsigned int primitive_0_index, unsigned int obj_No_1, unsigned int primitive_1_index);
	void getVTUnfixedPairIndex(int* unfixed_pair_vertex_index, int& unfixed_num, unsigned int vertex_obj_no,
		unsigned int vertex_index, unsigned int triangle_obj_No, unsigned int triangle_index);

	void solveBlockWithPair(
		double dt,
		std::array<double, 3>** record_vertex_position,
		int** record_vertex_num, int* unfixed_pair_vertex_index, int unfixed_num, double* common_grad,
		std::unordered_map<std::array<unsigned int, 2>, StoreHessianWithOrderInConstraint, pair_hash>& collision_hessian,
		double* floor_map, bool only_solve_collision_pair, double& max_dis);




	//void testPrintOut()
	//{
	//	for (int i = 0; i < tetrahedron->size(); ++i) {
	//		for (int j = 0; j < tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups.size(); ++j) {
	//			for (int k = 0; k < tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[j].size(); ++k) {
	//				std::cout << "print tet_in_a_group_start_per_thread_groups "<<i << std::endl;
	//				for (int m = 0; m < total_thread_num; ++m) {
	//					std::cout << tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[j][k][m] << " ";
	//				}
	//				std::cout << std::endl;
	//			}
	//		}
	//	}
	//	std::cout << "print vt_pair_compressed_record " << collision.vt_pair_compressed_record.size() << std::endl;
	//	for (int i = 0; i <= total_thread_num; ++i) {
	//		std::cout << collision.vt_per_thread_start_index[i] << " ";
	//	}
	//	std::cout << std::endl;
	//}

	std::vector<unsigned int> record_max_displace_vertex;
	//std::vector<double> record_max_displacement;

	std::vector<unsigned int>vertex_index_prefix_sum_obj;
	std::vector<unsigned int>global_vertex_index_start_per_thread;
	void checkExceedFloor();
	std::vector<double>energy_per_thread;
	std::vector< std::vector<std::vector<unsigned int>>*> vertex_index_of_a_tet_color_group;

	void lineSearchLastColor();
	void lineSearchFirstColor(unsigned int color_No, double ori_energy);

	unsigned int inner_itr_num_standard;



	//double computePreviousColorInertialEnergyTest(unsigned int color_No);
	//double computePreviousColorCollisionEnergyTest();
	//double computeColorInertialEnergyTest();
	//double computeLastColorARAPEnergyTest();
	//double computeBarrierEnergyTest();

	void test();

	void compareIfRecordHessianIsRight(int color);

	void compareVector(std::vector<unsigned int>* a, std::vector<unsigned int>* b, int type);

	int color_group_num;

	std::vector<double>grad_max_store;

	std::vector<double>max_dis_record;

	double max_displacement;

	double grad_max;

};



