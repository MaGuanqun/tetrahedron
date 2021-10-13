#pragma once
#include"external/Eigen/Dense"
#include"basic/eigenDenseOperation.h"
#include"external/Eigen/Sparse"
#include"external/Eigen/SparseCholesky"
//#include"basic/EigenMatrixIO.h"
#include"thread.h"
#include"basic/global.h"
#include"object/cloth.h"
#include"object/tetrohedron_object.h"

using namespace Eigen;
using namespace denseOperation;

class ProjectDynamic
{
public:
	ProjectDynamic();
	double time_step;
	double gravity_;
	double outer_itr_conv_rate, local_global_conv_rate;
	void setForPD(std::vector<Cloth>* cloth, std::vector<TetrohedronObject>* tetrohedron, Thread* thread);
	void reset();
	void initial();
	
	void PDsolve();
	int local_global_iteration_num, outer_iteration_num;
	void updateMatrixPerThread(int thread_No);
	void matrixDecomposition(int thread_id);
	void localProjectionPerThread(int thread_id);
	void solveSystemPerThead(int thread_id);
	void updateUVPerThread(int thread_id);
	void updateRenderPosition();
private:
	int total_thread_num;
	std::vector<double> temEnergy;
	int sub_step_num;
	bool use_dierct_solve_for_coarest_mesh;
	int super_jacobi_step_size;
	int max_it;
	int max_jacobi_itr_num;
	int total_cloth_num, total_tetrohedron_num;
	std::vector<int>cloth_sys_size, tetrohedron_sys_size;
	std::vector<Cloth>* cloth;
	std::vector<TetrohedronObject>* tetrohedron;
	Thread* thread;
	double sub_time_step;
	std::vector<std::vector<double>> lbo_weight;
	std::vector<VectorXd> cloth_mass_inv;
	std::vector<VectorXd> cloth_mass;
	std::vector<std::vector<VectorXd>> cloth_f_ext;
	std::vector<std::vector<VectorXd>> cloth_gravity;
	std::vector<std::vector<VectorXd>>vertex_lbo;
	void setForClothPD(std::vector<Cloth>* cloth);
	void setForTetrohedronPD(std::vector<TetrohedronObject>* tetrohedron);
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

	std::vector<SparseMatrix<double>> cloth_global_mat;
	std::vector<SparseMatrix<double>> initial_cloth_global_mat;
	std::vector<std::vector<double>>cloth_global_mat_diagonal_ref;
	std::vector<std::vector<double*>>cloth_global_mat_diagonal_ref_address;
	SimplicialLLT<SparseMatrix<double>>* cloth_llt;
	void computeGlobalStepMatrixSingleCloth(SparseMatrix<double>* global_mat, std::vector<double>& global_mat_diagonal_ref,
		std::vector<double*>& global_collision_mat_diagonal_ref, SimplicialLLT<SparseMatrix<double>>* global_llt, TriangleMeshStruct& mesh_struct,
		std::vector<std::vector<int>>& vertex_around_vertex_for_bending,
		std::vector<double>& lbo_weight, std::vector<VectorXd>& vertex_lbo, int sys_size, double bending_stiffness, std::vector<double>& length_stiffness,
		double position_stiffness);
	std::vector<std::vector<double>> rest_mean_curvature_norm;
	std::vector<std::vector<Vector3d>> p_edge_length;
	std::vector < std::vector<Vector3d>> p_bending;
	std::vector < std::vector<VectorXd>> cloth_b;
	std::vector < std::vector<VectorXd>> cloth_u_prediction;
	std::vector < std::vector<VectorXd>> cloth_v;
	std::vector < std::vector<VectorXd>> cloth_v_;
	std::vector < std::vector<VectorXd>> cloth_u;
	std::vector < std::vector<VectorXd>> cloth_u_;
	std::vector < std::vector<VectorXd>> cloth_acceleration;

	std::vector<int> cloth_per_thread_begin; 
	std::vector<std::vector<int>> cloth_dimension_per_thread;//every sytem need to solve Ax=b, when construct b, we rearrange the xyz dimensions of every b to different thread
	std::vector<int> tetrohedron_per_thread_begin;
	std::vector<std::vector<int>> tetrohedron_dimension_per_thread;
	void restBendingMeanCurvatureSingleCloth(TriangleMeshStruct& mesh_struct, std::vector<double>& rest_mean_curvature_norm,
		std::vector<std::vector<int>>& vertex_around_vertex_for_bending, std::vector<VectorXd>& vertex_lbo);
	void PDClothPredict();
	void setClothMatrix(int thread_No);
	void localEdgeLengthProjectionPerThread(int thread_id);
	void localBendingProjectionPerThread(int thread_id);
	void localPositionProjectionPerThread(int thread_id);
	void solveClothSystemPerThead(int thread_id);
	void solveClothSystemPerThead(VectorXd& b, VectorXd& u, TriangleMeshStruct& mesh_struct, std::vector<double>& length_stiffness,
		std::vector<Vector3d>& p_edge_length, int cloth_No, int dimension, std::vector<std::vector<int>>& vertex_around_vertex_for_bending,
		std::vector<VectorXd>& vertex_lbo, std::vector<Vector3d>& p_bending, std::vector<double>& lbo_weight, VectorXd& u_prediction, int thread_id);

	void updateModelPosition();
	void setIndexPerThread();

	std::vector < std::vector<VectorXd>> tetrohedron_b;
	std::vector < std::vector<VectorXd>> tetrohedron_u_prediction;
	std::vector < std::vector<VectorXd>> tetrohedron_v;
	std::vector < std::vector<VectorXd>> tetrohedron_v_;
	std::vector < std::vector<VectorXd>> tetrohedron_u;
	std::vector < std::vector<VectorXd>> tetrohedron_u_;
	std::vector < std::vector<VectorXd>> tetrohedron_acceleration;

	std::vector<std::vector<VectorXd>> tetrohedron_f_ext;
	std::vector<std::vector<VectorXd>> tetrohedron_gravity;

	SimplicialLLT<SparseMatrix<double>>* tetrohedron_llt;
	std::vector<SparseMatrix<double>> tetrohedron_global_mat;
	std::vector<SparseMatrix<double>> initial_tetrohedron_global_mat;
	void PDTetrohedronPredict();

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
	void updateTetrohedronUV(int thread_id);
};
