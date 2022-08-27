#pragma once
#include"external/Eigen/Dense"
#include"basic/eigenDenseOperation.h"
#include"external/Eigen/Sparse"
#include"external/Eigen/SparseLU"
#include"thread.h"
#include"basic/global.h"
#include"object/cloth.h"
#include"object/tetrahedron.h"
#include"object/collider.h"
#include"collision/collision.h"

using namespace Eigen;
using namespace denseOperation;

class SecondOrderLargeSystem
{
public:
	SecondOrderLargeSystem();

	double time_step;
	double gravity_;

	unsigned int iteration_number;
	double conv_rate;

	void setForNewtonMethod(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
		Thread* thread, double* tolerance_ratio);
	void massSpring(int thread_No);

	void setHessianDiagonal(int thread_No);

	void hessianCoeffAddress(int thread_No);

	void updateHessianFixedStructure(int thread_No);
	void setHessianDiagonalFixedStructure(int thread_No);
	void updateInternalForce(int thread_No);

	void sumB(int thread_No);

	void setSn(int thread_No);

	void solveNewtonMethod();
	void updatePosition(int thread_No);

	void initial();
	void initialDHatTolerance(double ave_edge_length);
	void addExternalForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int cloth_No);

	void updateVelocity(int thread_No);

	void updateIndexBeginPerObj();

	size_t* time_stamp;
	void resetExternalForce();

	double* damp_stiffness;

	double* rayleigh_damp_stiffness;

	Collision collision;

	void computeEnergy(int thread_No);

	void updatePositionFromOri(int thread_No);

	void setHessianDiagonalFixedStructureInitialStiffness(int thread_No);
private:


	unsigned int max_itr_num;

	unsigned int total_thread_num;

	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;

	bool perform_collision;


	SparseMatrix<double> Hessian;

	std::vector<double*> Hessian_coeff_address;

	unsigned int sys_size;
	unsigned int total_obj_num;

	std::vector<std::array<double, 3>*> vertex_position;
	std::vector<std::array<double, 3>*> render_position;
	void reorganzieDataOfObjects();

	std::vector<Triplet<double>> hessian_nnz;
	void initialHessianNnz();
	//void recordEdgeHessian();

	std::vector<unsigned int*> edge_index_begin_per_thread_for_mass_spring; //store unfixed edges
	std::vector<unsigned int*> only_one_vertex_fixed_edge_index_begin_per_thread;

	std::vector<unsigned int>total_vertex_num;

	//std::vector<unsigned int*> vertex_index_begin_per_thread;

	std::vector<std::vector<std::array<double, 3>>*> anchor_position;

	std::vector<unsigned int*> anchor_vertex_begin_per_thread;
	std::vector<std::vector<int>*> anchor_vertex;

	std::vector<unsigned int*> unfixed_vertex_begin_per_thread;
	std::vector<std::vector<unsigned int>*> unfixed_vertex;


	VectorXd Sn; //Sn=M(x_n + delta_t * v_n) //also store f_ext
	VectorXd f_ext;

	VectorXd position_of_max_step;

	VectorXd gravity;

	VectorXd velocity;


	VectorXd delta_x; // displacement in a newton iteration step


	std::vector<std::vector<double>> hessian_coeff_diagonal;


	std::vector<unsigned int*> real_index_to_unfixed_index;

	std::vector<unsigned int> vertex_begin_per_obj;
	std::vector<unsigned int> off_diagonal_hessian_nnz_index_begin_per_thread;

	std::vector<unsigned int> unfixed_gradC_hessian_index_begin_per_thread;

	std::vector<unsigned int>fixed_gradC_hessian_index_begin_per_thread;


	//std::vector<unsigned int> edge_begin_per_obj;

	std::vector<std::vector<unsigned int>*> edge_vertices_mass_spring;
	std::vector<std::vector<unsigned int>*> only_one_vertex_fix_edge_vertices;

	std::vector<double*>edge_length_stiffness;
	std::vector<double>anchor_stiffness;

	std::vector<std::vector<double>*>unfixed_rest_length;//
	std::vector<std::vector<double>*>fixed_one_vertices_rest_length;//


	std::vector<unsigned int> unfixed_constraint_start_per_thread_in_system;
	std::vector<unsigned int> fixed_vertex_constraint_start_per_thread_in_system;

	std::vector<double*> mass;

	void computeHessian(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
		double* diagonal_coeff_1, Triplet<double>* hessian_nnz, double alpha, double rest_length, unsigned int start_index_in_system_0,
		unsigned int start_index_in_system_1, Triplet<double>* hessian_nnz_grad_C, double lambda, unsigned int constraint_start_in_sys);

	void computeHessianFixedOneVertex(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
		double alpha, double rest_length, unsigned int start_index_in_system_0,
		Triplet<double>* hessian_nnz_grad_C, double lambda, unsigned int constraint_start_in_sys);


	void setHessian();

	double time_step_square;

	void getHessianCoeffAddress(double** address, unsigned int start_index_in_system_0, unsigned int start_index_in_system_1, double** grad_C_address, unsigned int constraint_start_index_0);
	void getHessianCoeffAddressFixedOneAddress(unsigned int start_index_in_system_0, double** grad_C_address, unsigned int constraint_start_index_0);
	void computeHessianFixedStructure(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
		double* diagonal_coeff_1, double** hessian_coeff_address, double alpha, double** grad_C_address,
		double lambda);


	void updateHessianFixedStructure();

	void updateInternalForce(double* vertex_position_0, double* vertex_position_1, double* force_0,
		double* force_1, double rest_length, double* h, double lambda, double alpha);

	void computeGravity();


	//SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >global_llt;
	SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> global_llt;

	void updateRenderPosition();
	bool convergenceCondition();



	void computeHessianOnlyOneVertexFixedEdge(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
		double alpha, double** grad_C_address, double lambda);
	void updateInternalForceOnlyOneEdgeFixed(double* vertex_position_0, double* vertex_position_1, double* force_0,
		double alpha, double rest_length, double* h, double lambda);


	void computeMassSpringEnergy(int thread_No);
	void computeInertial(int thread_No);
	double computeMassSpringEnergy(double* position_0, double* position_1, double rest_length, double stiffness);
	std::vector<double> energy_per_thread;
	void computeEnergy();
	double total_energy;
	double previous_energy;
	double displacement_coe;

	bool change_direction;
	void computeResidual();

	double residual;

	std::vector<double>store_residual;
	std::vector<double>log_store_residual;

	std::vector<VectorXd> b_thread;

	double beta, gamma;

	void 	solveNewtonMethod_();

	//void setK();


	std::vector<double>previous_frame_edge_length_stiffness;
	bool edgeLengthStiffnessHasChanged();

	unsigned int sys_total_size;
	unsigned int sys_total_size_index;

	std::vector<std::vector<double>> lambda_unfixed;
	std::vector<std::vector<double>> lambda_fixed_one_vertex;

	void initialLambda();

	void resetLambda();
	void updateNormal();
	void updatePositionFromSn();
};
