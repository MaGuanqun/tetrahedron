#pragma once
#include"external/Eigen/Dense"
#include"basic/eigenDenseOperation.h"
#include"external/Eigen/Sparse"
#include"external/Eigen/SparseCholesky"
//#include"basic/EigenMatrixIO.h"
#include"thread.h"
#include"basic/global.h"
#include"object/cloth.h"
#include"object/tetrahedron.h"
#include"object/collider.h"
#include"collision/collision.h"
#include"iteration_method.h"

using namespace Eigen;
using namespace denseOperation;

class ProjectDynamic
{
public:
	ProjectDynamic();
	double time_step;
	double gravity_;
	double outer_itr_conv_rate, local_global_conv_rate;
	void setForPD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Thread* thread);
	void reset();
	void initial();
	
	void PDsolve();
	void PD_IPC_solve();
	int local_global_iteration_num, outer_iteration_num;
	void updateMatrixPerThread(int thread_No);
	void matrixDecomposition(int thread_id);
	void localProjectionPerThread(int thread_id, bool with_energy);
	
	void solveSystemPerThead(int thread_id, bool with_collision, bool compute_energy);
	void updateUVPerThread(int thread_id);
	void updateRenderPosition();
	void addExternalClothForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int cloth_No);
	void resetExternalForce();
	void mainProcess();
	void testLocalProjectionPerThread(int thread_id);
	void computeCollisionFreePosition(int thread_No);

	void initialDHatTolerance(double ave_edge_length);
	void computeDisplacement(int thread_No);

	int itr_solver_method;

	void updateIterateSolverParameter(double conv_rate);

private:
	int local_global_itr_in_single_outer;
	int total_thread_num;
	std::vector<double> temEnergy;
	int sub_step_num;
	bool use_dierct_solve_for_coarest_mesh;
	int super_jacobi_step_size;
	int max_it;
	int max_jacobi_itr_num;
	int total_cloth_num, total_tetrahedron_num, total_collider_num;
	std::vector<int>cloth_sys_size, tetrahedron_sys_size;
	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;
	double sub_time_step;
	std::vector<std::vector<double>> lbo_weight;
	std::vector<VectorXd> cloth_mass_inv;
	std::vector<VectorXd> cloth_mass;
	std::vector<std::vector<VectorXd>> cloth_f_ext;
	std::vector<std::vector<VectorXd>> cloth_gravity;
	std::vector<std::vector<VectorXd>>vertex_lbo;
	void setForClothPD(std::vector<Cloth>* cloth);
	void setForTetrahedronPD(std::vector<Tetrahedron>* tetrahedron);
	void computeEdgeCotWeight(std::vector<std::vector<double>>& edge_cot_weight);
	void computeEdgeCotWeightSingleCloth(std::vector<double>& edge_cot_weight, TriangleMeshStruct& mesh_struct);
	void computeLBOWeight(std::vector<std::vector<double>>& edge_cot_weight);
	void computeClothGravity();
	void computeVertexLBO(std::vector<std::vector<double>>& edge_cot_weight);

	void computeGlobalStepMatrix();
	void initialClothPDvariable();
	void restBendingMeanCurvature();
	std::vector<std::vector<std::vector<int>>> vertex_around_vertex_for_bending;
	void computeLBOWeightSingleCloth(std::vector<double>& edge_cot_weight, std::vector<double>& lbo_weight, TriangleMeshStruct& mesh_struct,
		VectorXd& mass_inv, double density, VectorXd& mass);
	void setAroundVertexPrimitive(TriangleMeshStruct& mesh_struct, std::vector<std::vector<int>>& edge_around_vertex_for_bending,
		std::vector<std::vector<int>>& vertex_around_vertex_for_bending);
	void computeVertexLBOSingleCloth(TriangleMeshStruct& mesh_struct, std::vector<VectorXd>& vertex_lbo, std::vector<double>& edge_cot_weight,
		std::vector<std::vector<int>>& edge_around_vertex_for_bending, std::vector<double>& lbo_weight);

	std::vector<SparseMatrix<double, RowMajor>> cloth_global_mat;
	std::vector<SparseMatrix<double, RowMajor>> initial_cloth_global_mat;
	std::vector<std::vector<double>>cloth_global_mat_diagonal_ref;
	std::vector<std::vector<double*>>cloth_global_mat_diagonal_ref_address;
	SimplicialLLT<SparseMatrix<double>>* cloth_llt;

	SimplicialLLT<SparseMatrix<double>>* collision_free_cloth_llt;

	void computeGlobalStepMatrixSingleCloth(SparseMatrix<double, RowMajor>* global_mat, std::vector<double>& global_mat_diagonal_ref,
		std::vector<double*>& global_collision_mat_diagonal_ref, SimplicialLLT<SparseMatrix<double>>* global_llt, TriangleMeshStruct& mesh_struct,
		std::vector<std::vector<int>>& vertex_around_vertex_for_bending,
		std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo, int sys_size, double bending_stiffness, std::vector<double>& length_stiffness,
		double position_stiffness, SimplicialLLT<SparseMatrix<double>>* collision_free_global_llt);
	std::vector<std::vector<double>> rest_mean_curvature_norm;
	std::vector<std::vector<Vector3d>> p_edge_length;
	std::vector<std::vector < std::vector<VectorXd>>> p_bending;
	std::vector < std::vector<VectorXd>> cloth_b;
	std::vector < std::vector<VectorXd>> cloth_u_prediction;
	std::vector < std::vector<VectorXd>> cloth_v;
	std::vector < std::vector<VectorXd>> cloth_v_;
	std::vector < std::vector<VectorXd>> cloth_u;
	std::vector < std::vector<VectorXd>> cloth_u_;
	std::vector < std::vector<VectorXd>> cloth_acceleration;

	std::vector<int> cloth_per_thread_begin; 
	std::vector<std::vector<int>> cloth_dimension_per_thread;//every sytem need to solve Ax=b, when construct b, we rearrange the xyz dimensions of every b to different thread
	std::vector<int> tetrahedron_per_thread_begin;
	std::vector<std::vector<int>> tetrahedron_dimension_per_thread;
	void restBendingMeanCurvatureSingleCloth(TriangleMeshStruct& mesh_struct, std::vector<double>& rest_mean_curvature_norm,
		std::vector<std::vector<int>>& vertex_around_vertex_for_bending, std::vector<VectorXd>& vertex_lbo);
	void PDClothPredict();
	void setClothMatrix(int thread_No);
	void localEdgeLengthProjectionPerThread(int thread_id, bool with_energy);
	void localBendingProjectionPerThread(int thread_id, bool with_energy);
	void localPositionProjectionPerThread(int thread_id, bool with_energy);
	void solveClothSystemPerThead(int thread_id, bool with_collision, bool compute_energy);
	void solveClothSystemPerThead(VectorXd& b, VectorXd& u, TriangleMeshStruct& mesh_struct, std::vector<double>& length_stiffness,
		std::vector<Vector3d>& p_edge_length, int cloth_No, int dimension, std::vector<std::vector<int>>& vertex_around_vertex_for_bending,
		std::vector<VectorXd>& vertex_lbo, std::vector<std::vector<VectorXd>>& p_bending, VectorXd& u_prediction, int thread_id,
		std::vector<std::array<double, 3>>& collision_b_sum, bool* collision_b_need_update, bool with_collision, bool compute_energy);

	void updateModelPosition();
	void setIndexPerThread();

	std::vector < std::vector<VectorXd>> tetrahedron_b;
	std::vector < std::vector<VectorXd>> tetrahedron_u_prediction;
	std::vector < std::vector<VectorXd>> tetrahedron_v;
	std::vector < std::vector<VectorXd>> tetrahedron_v_;
	std::vector < std::vector<VectorXd>> tetrahedron_u;
	std::vector < std::vector<VectorXd>> tetrahedron_u_;
	std::vector < std::vector<VectorXd>> tetrahedron_acceleration;

	std::vector<std::vector<VectorXd>> tetrahedron_f_ext;
	std::vector<std::vector<VectorXd>> tetrahedron_gravity;

	SimplicialLLT<SparseMatrix<double>>* tetrahedron_llt;
	std::vector<SparseMatrix<double>> tetrahedron_global_mat;
	std::vector<SparseMatrix<double>> initial_tetrahedron_global_mat;
	void PDTetrahedronPredict();

	double previous_itr_PD_energy, previous_itr_constraint_energy, previous_itr_collision_energy;
	double current_constraint_energy, previous_constraint_energy;
	double current_collision_energy, previous_collision_energy;
	double	current_PD_energy, previous_PD_energy;

	void PDupdateSystemMatrix();

	bool PDConvergeCondition();
	bool PDLocalGlobalConvergeCondition();
	void initialEnergy();
	void localProjection();
	void updateCollisionMatrix();
	void updateMatrix();
	void initialEnergyOuterInteration();
	void initialEnergyLocalGlobal();
	void PDsetPosPredict();
	void updateClothUV(int thread_id);
	void updateTetrahedronUV(int thread_id);

	Collision collision;
	void initialTetrahedronPDvariable();
	bool IPC_PDConvergeCondition();
	void updateRenderPositionIPC();
	void firstPDForIPC();

	std::vector<double>displacement_norm_thread;
	std::vector<std::vector<VectorXd>> cloth_u_previous_itr;

	double displacement_bound;
	void computeEnergyIPCPD();
	double PD_energy_dif;
	double collision_energy_dif;
	double displacement_ratio_dif;
	double displacement_norm;
	double previous_displacement_norm;
	bool innerIterationConvergeCondition();
	void computeInnerEnergyIPCPD();


	
	void initialJacobi();
	void setOffDiagonal(int cloth_No, std::vector<std::vector<int>>& vertex_around_vertex_for_bending,
		std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo, int sys_size, double bending_stiffness, std::vector<double>& length_stiffness,
		TriangleMeshStruct& mesh_struct);
	void computeOffDiagonal();
	IterationMethod iteration_method;
};
