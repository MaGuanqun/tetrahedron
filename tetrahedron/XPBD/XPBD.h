#pragma once
#include"../external/Eigen/Dense"
#include"../basic/eigenDenseOperation.h"
#include"../external/Eigen/Sparse"
#include"../external/Eigen/SparseCholesky"
//#include"basic/EigenMatrixIO.h"
#include"../thread.h"
#include"../basic/global.h"
#include"../object/cloth.h"
#include"../object/tetrahedron.h"
#include"../object/collider.h"
#include"../collision/collision.h"
#include"XPBD_constraint.h"


using namespace Eigen;
using namespace denseOperation;

class XPBD
{
public:
	XPBD();
	double time_step;
	double gravity_;
	unsigned int sub_step_num;
	void PBDsolve();
	void setForXPBD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
		Thread* thread, double* tolerance_ratio);
	size_t* time_stamp;
	void setPosPredict(int thread_No);
	void setPosPredictSubTimeStep(int thread_No);
	void computeVelocity(int thread_No);
	unsigned int iteration_number;
	//unsigned int inner_iteration_number;
	void initial();
	void reset();
	void resetExternalForce();
	void initialDHatTolerance(double ave_edge_length);
	void updateTetrahedronAnchorVertices();
	void addExternalForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int obj_No);
	void updateItrInfo(int* iteration_num);

	void PBD_IPCSolve();
	Collision collision;

private:
	double gravity[3];

	unsigned int total_thread_num;
	
	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;
	double sub_time_step;
	std::vector<std::vector<std::array<double, 3>>> f_ext;
	std::vector<std::vector<std::array<double, 3>>> velocity;
	//std::vector<std::vector<std::array<double, 3>>> total_gravity;

	std::vector<std::array<double, 3>*> vertex_position;
	std::vector<std::array<double, 3>*> initial_vertex_position;
	std::vector<MeshStruct*> mesh_struct;
	std::vector<unsigned int*> vertex_index_begin_per_thread;
	void reorganzieDataOfObjects();
	unsigned int total_obj_num;
	void initialVariable();
	XPBDconstraint XPBD_constraint;
	std::vector<std::vector<double>> lbo_weight;
	std::vector<std::vector<VectorXd>>vertex_lbo;
	std::vector<std::vector<double>> rest_mean_curvature_norm;
	//std::vector<std::vector<Vector3d>> rest_Aq;


	std::vector<double> lambda;
	
	std::vector<double> lambda_collision;
	std::vector<unsigned int>constraint_index_start;
	std::vector<std::vector<unsigned int>>collision_constraint_index_start;

	void initialClothBending();
	void solveBendingConstraint();
	void solveEdgeLengthConstraint();
	void solveConstraint();
	void setConstraintIndex();
	void solveTetStrainConstraint();

	void updatePosition();
	double damping_coe;

	void updateRenderNormal();
	void initialCollisionConstriantNum();
	bool perform_collision;

	void updateNormal();
	void updateRenderVertexNormal();

	bool convergeCondition(unsigned int iteration_num);

	std::vector<std::vector<std::array<double, 3>>> record_vertex_position;
	//std::vector<std::vector<std::array<double, 3>>> record_outer_vertex_position;
	void recordVertexPosition();

	std::vector<std::vector<unsigned int>* >unfixed_vertex;
	double max_move_standard;//the max displacement to stop iteration
	double outer_max_move_standard;//the max displacement to stop iteration
		

	double converge_condition_ratio;// converge_condition_ratio* edge length

	double calEdgeLength();
	void setConvergeCondition();

	unsigned int max_iteration_number;
	unsigned int outer_max_iteration_number;

	bool use_bending_based_on_vertex=true;

	double velocity_damp;

	//void recordOuterVertexPosition();
	//bool outerConvergeCondition(unsigned int iteration_num);

	bool use_PBD = false;
	void solveByPBD();
	void solveByXPBD();

};

