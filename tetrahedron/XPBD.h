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

class XPBD
{
public:
	XPBD();
	double time_step;
	double gravity_;
	void PBDsolve();
	void setForPD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider,
		Thread* thread, double* tolerance_ratio);
	Collision collision;
	size_t* time_stamp;
	void setPosPredict(int thread_No);
private:
	unsigned int sub_step_num;
	unsigned int total_thread_num;
	
	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;
	double sub_time_step;
	std::vector<std::vector<std::array<double, 3>>> f_ext;
	std::vector<std::vector<std::array<double, 3>>> velocity;
	std::vector<std::vector<std::array<double, 3>>> total_gravity;

	std::vector<std::array<double, 3>*> vertex_position;
	std::vector<MeshStruct*> mesh_struct;
	std::vector<unsigned int*> vertex_index_begin_per_thread;
	void reorganzieDataOfObjects();
	unsigned int total_obj_num;
	void initialVariable();
};

