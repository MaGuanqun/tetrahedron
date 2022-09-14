#pragma once
#include"../basic/global.h"
#include"../external/Eigen/Dense"
#include<array>

using namespace Eigen;

class ComputeEnergy
{
public:
	double computeMassSpringEnergy(double* position_0, double* position_1, double rest_length, double stiffness);
	double computeARAPEnergy(double* position_0, double* position_1, double* position_2, double* position_3, Matrix<double, 3, 4>& A, double volume, double stiffness);
	double computeInertial(double time_step, unsigned int index_start, unsigned int index_end, double* vertex_pos, double* mass_, VectorXd& Sn, unsigned int vertex_start,
		unsigned int* unfixed_index_to_normal_index);

private:

};

