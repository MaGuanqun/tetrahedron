#pragma once
#include"external/Eigen/Dense"
#include"basic/eigenDenseOperation.h"
#include"external/Eigen/Sparse"
#include"external/Eigen/SparseCholesky"
#include"thread.h"
#include"basic/global.h"
#include"object/cloth.h"
#include"object/tetrahedron.h"
#include"object/collider.h"
#include"collision/collision.h"
#include"compute_energy.h"

using namespace Eigen;
using namespace denseOperation;

class NewtonMethod
{
public:
	NewtonMethod();

	double time_step;
	double time_step_square;
	double gravity_;

	unsigned int iteration_number;
	double conv_rate;

	void setForNewtonMethod(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
		Thread* thread, double* tolerance_ratio);
	void massSpring(int thread_No);
	//void initial_hessian_coeff(int thread_No);
	void setHessianDiagonal(int thread_No);

	void hessianCoeffAddress(int thread_No);

	void updateHessianFixedStructure(int thread_No);
	void setHessianDiagonalFixedStructure(int thread_No);
	void updateInternalForce(int thread_No);

	void sumB(int thread_No);

	void setSn(int thread_No);

	void solveNewtonMethod();
	void updatePosition(int thread_No);

	void updateHessianForFixPoint(int thread_No);

	void initial();
	void initialDHatTolerance(double ave_edge_length);
	void addExternalForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int cloth_No);

	void updateVelocity(int thread_No);
	void updateVelocity2(int thread_No);
	void updateIndexBeginPerObj();

	size_t* time_stamp;
	void resetExternalForce();

	double* damp_stiffness;

	double* rayleigh_damp_stiffness;

	void updateHessianForDamp(int thread_No);
	Collision collision;

	void computeEnergy(int thread_No);

	void updatePositionFromOri(int thread_No);



	void setBn(int thread_No);
	void updateVelocityAccelerationNewMark(int thread_No);
	void setHessianDiagonalFixedStructureInitialStiffness(int thread_No);
private:


	unsigned int max_itr_num;

	unsigned int total_thread_num;

	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;

	bool perform_collision;


	VectorXd b; //solve Ax=b

	SparseMatrix<double> Hessian;
	SparseMatrix<double> ori_stiffness_matrix_record;

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

	VectorXd position_of_beginning;

	VectorXd position_of_max_step;

	VectorXd gravity;

	VectorXd velocity;
	VectorXd acceleration;

	VectorXd delta_x; // displacement in a newton iteration step


	//struct HessianCoeff
	//{
	//	unsigned int row;
	//	unsigned int col;
	//	double coeff;
	//	HessianCoeff(unsigned int row_, unsigned int col_, double coeff_)
	//	{
	//		row = row_;
	//		col = col_;
	//		coeff = coeff_;
	//	}
	//};

	//std::vector<std::vector<HessianCoeff>> hessian_coeff_off_diagonal;// only store <i,j>
	std::vector<std::vector<double>> hessian_coeff_diagonal;// only store block<i,i> coefficient, also use this variable temporally store f_int in each thread 


	std::vector<unsigned int*> real_index_to_unfixed_index;

	std::vector<unsigned int> vertex_begin_per_obj;
	std::vector<unsigned int> off_diagonal_hessian_nnz_index_begin_per_thread;
	//std::vector<unsigned int> edge_begin_per_obj;

	std::vector<std::vector<unsigned int>*> edge_vertices_mass_spring;
	std::vector<std::vector<unsigned int>*> only_one_vertex_fix_edge_vertices; //first index is unfixed (use unfixed index, need to transfer it to get the real index), second index is fixed, it is directly the real index

	std::vector<double*>edge_length_stiffness;
	std::vector<double>anchor_stiffness;

	std::vector<std::vector<double>*>unfixed_rest_length;//
	std::vector<std::vector<double>*>fixed_one_vertices_rest_length;//


	std::vector<double*> mass;

	void computeHessian(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
		double* diagonal_coeff_1, Triplet<double>* hessian_nnz, double stiffness, double rest_length, unsigned int start_index_in_system_0,
		unsigned int start_index_in_system_1, double time_step_square);

	void setHessian();



	void getHessianCoeffAddress(double** address, unsigned int start_index_in_system_0, unsigned int start_index_in_system_1);

	void computeHessianFixedStructure(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
		double* diagonal_coeff_1, double** hessian_coeff_address, double stiffness, double rest_length, double time_step_square,
		double damp_stiffness, double* damp_force_0, double* damp_force_1);


	void updateHessianFixedStructure();

	void updateInternalForce(double* vertex_position_0, double* vertex_position_1, double* force_0,
		double* force_1, double stiffness, double rest_length, double* ori_position_0, double* ori_position_1, double damp_stiffness);

	void computeGravity();


	//SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >global_llt;
	SimplicialLDLT<SparseMatrix<double>> global_llt;

	void updateRenderPosition();
	bool convergenceCondition();

	void storeInitialPosition();

	void computeHessianOnlyOneVertexFixedEdge(double* vertex_position_0, double* vertex_position_1, double* diagonal_coeff_0,
		double stiffness, double rest_length, double time_step_square, double damp_stiffness,
		double* damp_force_0);
	void updateInternalForceOnlyOneEdgeFixed(double* vertex_position_0, double* vertex_position_1, double* force_0,
		double stiffness, double rest_length, double* ori_position_0, double* ori_position_1, double damp_stiffness);

	//double damp_coe;

	//unsigned int obj_dimension = 2;

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
	bool is_newmark = false;
	void 	solveNewtonMethod_();
	void newmarkBetaMethod();
	void updateK();
	void setK();
	void setBn();

	VectorXd pos_dis;

	std::vector<double> store_ori_value;

	double previous_frame_rayleigh_damp_stiffness_beta;
	std::vector<double>previous_frame_edge_length_stiffness;
	bool edgeLengthStiffnessHasChanged();



	unsigned int tet_hessian_index[2];
	void initialCoeffDiagonal();
	void setARAP();

	std::vector<std::array<int, 4>*> tet_indices;
	std::vector<TetrahedronMeshStruct*> tet_mesh_struct;
	std::vector<double*> inv_mass;


	void computeARAPHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
		double* diagonal_coeff_0, double* diagonal_coeff_1, double* diagonal_coeff_2, double* diagonal_coeff_3, std::vector<Triplet<double>>* hessian_nnz,
		int* vertex_index, double stiffness, Matrix<double, 3, 4>& A, bool* is_unfixed);

	void getARAPCoeffAddress();

	void getARAPHessianCoeffAddress(int* start_index_in_system, std::vector<double*>* address, bool* is_unfixed);
	void updateARAPHessianFixedStructure();

	void computeARAPHessianFixedStructure(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
		double* diagonal_coeff_0, double* diagonal_coeff_1, double* diagonal_coeff_2, double* diagonal_coeff_3, double**& address,
		double stiffness, Matrix<double, 3, 4>& A, bool* is_unfixed, double* g_0, double* g_1, double* g_2, double* g_3);

	/*void testSetHessian();
	VectorXd b_for_test;
	void setARAPHessianForTest(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3, int* vertex_index,
		std::vector<Triplet<double>>* hessian_nnz, double stiffness, Matrix<double, 3, 4>& A, bool* is_unfixed, double* g_0, double* g_1, double* g_2, double* g_3);
	std::vector<Triplet<double>>test_nnz;
	SparseMatrix<double> Matrix_test;
	SparseMatrix<double> only_constraint;
	void updateTest();
	void construct_b_Test();
	*/
	void computeARAPEnergy();
	ComputeEnergy compute_energy;
};
