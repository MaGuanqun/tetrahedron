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
//#include<iomanip>

using namespace Eigen;
using namespace denseOperation;

class ProjectDynamic
{
public:
	ProjectDynamic();
	double time_step;
	double gravity_;
	double outer_itr_conv_rate, local_global_conv_rate;
	void setForPD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor, Thread* thread,
		double* tolerance_ratio);
	void reset();
	void initial();
	
	void PDsolve();
	void PD_IPC_solve(bool& record_matrix);
	int local_global_iteration_num, outer_iteration_num;
	void updateMatrixPerThread(int thread_No);
	void updateDiagonalPerThread(int thread_No);
	void localProjectionPerThread(int thread_id, bool with_energy);
	
	void constructbPerThead(int thread_id, bool with_collision);
	void updateUVPerThread(int thread_id);
	void updateRenderPosition();
	void addExternalClothForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int cloth_No);
	void addExternalTetForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int tet_No);
	void resetExternalForce();
	void mainProcess();
	void testLocalProjectionPerThread(int thread_id);
	void computeCollisionFreePosition(int thread_No);

	void initialDHatTolerance(double ave_edge_length);
	void computeDisplacement(int thread_No);

	int itr_solver_method;

	void updateIterateSolverParameter(double conv_rate);
	void update_ave_iteration_record(double& ave_itr);
	void computeEnergyPerThread(int thread_id);
	void solveSystemPerThread(int thread_id, bool with_collision);

	void updateTetrahedronAnchorVertices();
	size_t* time_stamp;
	Collision collision;

	void updateSystemPos();
	double velocity_damp;

	std::vector<std::array<double, 3>> position_record_to_show;

private:

	unsigned int tetrahedron_begin_obj_index;
	int local_global_itr_in_single_outer;
	unsigned int total_thread_num;
	std::vector<double> temEnergy;
	int sub_step_num;
	bool use_dierct_solve_for_coarest_mesh;
	int super_jacobi_step_size;
	int max_it;
	int max_jacobi_itr_num;
	unsigned int max_single_inner_iter_num;

	unsigned int total_cloth_num, total_tetrahedron_num, total_collider_num;
	std::vector<unsigned int>cloth_sys_size, tetrahedron_sys_size;
	unsigned int sys_size;
	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;
	double sub_time_step;
	std::vector<std::vector<double>> lbo_weight;
	VectorXd mass_inv;
	VectorXd mass;

	std::vector<VectorXd> f_ext;
	std::vector<VectorXd> total_gravity;

	std::vector<std::vector<VectorXd>>vertex_lbo;
	void setForClothPD(std::vector<Cloth>* cloth);
	void setForTetrahedronPD();
	void computeEdgeCotWeight(std::vector<std::vector<double>>& edge_cot_weight);
	void computeEdgeCotWeightSingleCloth(std::vector<double>& edge_cot_weight, TriangleMeshStruct& mesh_struct);
	void computeLBOWeight(std::vector<std::vector<double>>& edge_cot_weight);
	void computeGravity();
	void computeVertexLBO(std::vector<std::vector<double>>& edge_cot_weight);

	void initialPDvariable();
	void computeGlobalStepMatrix();
	void initialClothPDvariable();
	void restBendingMeanCurvature();
	//std::vector<std::vector<std::vector<int>>> vertex_around_vertex_for_bending;
	void computeLBOWeightSingleCloth(std::vector<double>& edge_cot_weight, std::vector<double>& lbo_weight, TriangleMeshStruct& mesh_struct,
		VectorXd& mass_inv, double density, VectorXd& mass, int vertex_index_start);
	void setAroundVertexPrimitive(TriangleMeshStruct& mesh_struct, std::vector<std::vector<int>>& edge_around_vertex_for_bending,
		std::vector<std::vector<int>>& vertex_around_vertex_for_bending);
	void computeVertexLBOSingleCloth(TriangleMeshStruct& mesh_struct, std::vector<VectorXd>& vertex_lbo, std::vector<double>& edge_cot_weight,
		std::vector<double>& lbo_weight);

	SparseMatrix<double, RowMajor> global_mat;
	SparseMatrix<double, RowMajor> initial_global_mat;
	std::vector<double>global_mat_diagonal_ref;
	std::vector<double>ori_global_mat_diagonal_ref;
	std::vector<double*>global_mat_diagonal_ref_address;
	SimplicialLLT<SparseMatrix<double>> global_llt;
	SimplicialLLT<SparseMatrix<double>> collision_free_llt;

	void computeGlobalStepMatrixSingleCloth(TriangleMeshStruct& mesh_struct, std::vector<Triplet<double>>& global_mat_nnz,
		std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo, int sys_size, double bending_stiffness, double length_stiffness,
		double position_stiffness, int vertex_index_start);
	std::vector<std::vector<double>> rest_mean_curvature_norm;
	std::vector<std::vector<Vector3d>> p_edge_length;
	std::vector<std::vector < std::vector<VectorXd>>> p_bending;

	std::vector<std::vector<Vector3d>> tet_edge_length;

	std::vector<VectorXd> b;
	std::vector<VectorXd> u_prediction;
	std::vector<VectorXd> v;
	std::vector<VectorXd> v_;
	std::vector<VectorXd> u;
	std::vector<VectorXd> u_;
	std::vector<VectorXd> acceleration;

	std::vector<unsigned int>system_vertex_index_per_thread;

	std::vector<unsigned int>dimension_per_thread;

	std::vector<unsigned int> cloth_per_thread_begin; 
	std::vector<std::vector<int>> cloth_dimension_per_thread;//every sytem need to solve Ax=b, when construct b, we rearrange the xyz dimensions of every b to different thread
	std::vector<unsigned int> tetrahedron_per_thread_begin;
	std::vector<std::vector<int>> tetrahedron_dimension_per_thread;
	void restBendingMeanCurvatureSingleCloth(TriangleMeshStruct& mesh_struct, std::vector<double>& rest_mean_curvature_norm,
		std::vector<VectorXd>& vertex_lbo);
	void PDPredict();
	void localEdgeLengthProjectionPerThread(int thread_id, bool with_energy);
	void localBendingProjectionPerThread(int thread_id, bool with_energy);
	void localPositionProjectionPerThread(int thread_id, bool with_energy);
	void constructbPerThead(double* b, TriangleMeshStruct& mesh_struct,
		std::vector<Vector3d>& p_edge_length, int cloth_No, int dimension,
		std::vector<VectorXd>& vertex_lbo, std::vector<std::vector<VectorXd>>& p_bending, int thread_id,
		std::vector<std::array<double, 3>>& collision_b_sum, bool* collision_b_need_update, bool with_collision, int vertex_index_start);

	void updateModelPosition();
	void setIndexPerThread();

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


	
	void initialTetrahedronPDvariable();
	bool IPC_PDConvergeCondition();
	void updateRenderPositionIPC();
	void firstPDForIPC(bool& record_matrix);

	std::vector<double>displacement_norm_thread;
	std::vector<VectorXd> u_previous_itr;

	double displacement_bound;
	void computeEnergyIPCPD();
	double PD_energy_dif;
	double collision_energy_dif;
	double displacement_ratio_dif;
	double displacement_norm;
	double previous_displacement_norm;
	bool innerIterationConvergeCondition();
	void computeInnerEnergyIPCPD();
	void solveClothSystem2(bool compute_energy);

	IterationMethod iteration_method;
	void saveMatrix(int dimension, VectorXd& vec, std::string name);
	void saveSparseMatrix(SparseMatrix<double, RowMajor>& matrix, std::string file_name_matrix);

	double ave_iteration;//the total iteration number of iteration method in a time step
	void reset_ave_iteration_record();

	int max_inner_iteration_num;

	std::vector<unsigned int> vertex_begin_per_cloth;
	std::vector<unsigned int> vertex_begin_per_tetrahedron;

	void setSystemIndexInfo();

	void copmuteGlobalStepMatrixSingleTetrahedron(TetrahedronMeshStruct& mesh_struct, std::vector<Triplet<double>>& global_mat_nnz, int sys_size, double& ARAP_stiffness,
		double& volume_preserve_stiffness, double position_stiffness, double length_stiffness, int vertex_index_start);
	Matrix4d getARAPmatrix();
	void updateTetrahedronAnchorVertices(int tetrahedron_index, TetrahedronMeshStruct& mesh_struct, int vertex_index_start, double position_stiffness);
	void localARAPProjectionPerThread(int thread_id, bool with_energy);

	std::vector<std::vector<Matrix<double,4,3>>> p_ARAP_volume_preserve;
	void constructbTetPerThead(double* b, TetrahedronMeshStruct& mesh_struct,
		std::vector<Matrix<double,4,3>>& p_ARAP_volume_preserve, int tet_No, int dimension,
		std::vector<std::array<double, 3>>& collision_b_sum, bool* collision_b_need_update, bool with_collision,
		double position_stiffness, int vertex_index_start, std::array<int, 4>* indices, std::vector<int>& anchor_vertex,
		std::array<double, 3>* anchor_pos, Vector3d* tet_edge_length);
	Matrix4d tet_local_A;
	void computeTetMass(double* mass_, VectorXd& mass_inv, VectorXd& mass, int vertex_index_start, int vertex_num);
	bool getDigonalForVolumePreserve(Vector3d& svd_eigen, double max, double min, Vector3d& sigma);

	void computeInnerCollisionFreeEnergy();

	void testJacobiMatrix();
	void testJacobi();
	void setOverRelaxationCoefficient();

	void localTetEdgeLengthProjectionPerThread(int thread_id, bool with_energy);

	bool perform_collision;

	void upadtePositionByRender();

};
