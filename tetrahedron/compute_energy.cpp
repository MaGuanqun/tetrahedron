#include"compute_energy.h"
#include"XPBD/FEM_relate.h"

double ComputeEnergy::computeMassSpringEnergy(double* position_0, double* position_1, double rest_length, double stiffness)
{
	double current_length;
	current_length = sqrt((position_0[0] - position_1[0]) * (position_0[0] - position_1[0])
		+ (position_0[1] - position_1[1]) * (position_0[1] - position_1[1])
		+ (position_0[2] - position_1[2]) * (position_0[2] - position_1[2]));
	current_length -= rest_length;
	return 0.5 * stiffness * current_length * current_length;
}

double ComputeEnergy::computeMassSpringConstraint(double* position_0, double* position_1, double rest_length, double stiffness,double lambda, double time_step_square)
{
	double current_length;
	current_length = sqrt((position_0[0] - position_1[0]) * (position_0[0] - position_1[0])
		+ (position_0[1] - position_1[1]) * (position_0[1] - position_1[1])
		+ (position_0[2] - position_1[2]) * (position_0[2] - position_1[2]));
	current_length -= rest_length;
	return -current_length*lambda-0.5/(stiffness*time_step_square)*lambda*lambda;
}

double ComputeEnergy::computeARAPEnergy(double* position_0, double* position_1, double* position_2, double* position_3, Matrix<double, 3, 4>& A, double volume, double stiffness)
{
	Matrix3d deformation_gradient;
	FEM::getDeformationGradient(position_0,position_1, position_2, position_3, A, deformation_gradient);
	JacobiSVD<Matrix3d> svd;
	svd.compute(deformation_gradient);
	Vector3d eigen_value = svd.singularValues();
	if (deformation_gradient.determinant() < 0) {
		eigen_value[2] *= -1.0;
	}
	double norm = (eigen_value[0] - 1.0) * (eigen_value[0] - 1.0) + (eigen_value[1] - 1.0) * (eigen_value[1] - 1.0) +
		(eigen_value[2] - 1.0) * (eigen_value[2] - 1.0);

	return 0.5*norm * stiffness * volume;
}



double ComputeEnergy::computeARAPConstraint(double* position_0, double* position_1, double* position_2, double* position_3, Matrix<double, 3, 4>& A, double volume, double stiffness, double lambda, double time_step_square)
{
	Matrix3d deformation_gradient;
	FEM::getDeformationGradient(position_0, position_1, position_2, position_3, A, deformation_gradient);
	JacobiSVD<Matrix3d> svd;
	svd.compute(deformation_gradient);
	Vector3d eigen_value = svd.singularValues();
	if (deformation_gradient.determinant() < 0) {
		eigen_value[2] *= -1.0;
	}
	double norm = (eigen_value[0] - 1.0) * (eigen_value[0] - 1.0) + (eigen_value[1] - 1.0) * (eigen_value[1] - 1.0) +
		(eigen_value[2] - 1.0) * (eigen_value[2] - 1.0);
	return -sqrt(norm)*lambda-0.5/(stiffness*volume* time_step_square)*lambda*lambda;
}

//compute for system that has removed fixed vertices
double ComputeEnergy::computeInertial(double time_step, unsigned int index_start, unsigned int index_end, double* vertex_pos, double* mass_,VectorXd& Sn, unsigned int vertex_start,
	unsigned int* unfixed_index_to_normal_index)
{
	unsigned int j,start;
	double energy = 0.0;
	for (unsigned int l = index_start; l < index_end; ++l) {
		j = 3 * l;
		start = 3 * unfixed_index_to_normal_index[l - vertex_start];
		energy += mass_[unfixed_index_to_normal_index[l - vertex_start]] *
			((vertex_pos[start] - Sn.data()[j]) * (vertex_pos[start] - Sn.data()[j]) +
				(vertex_pos[start + 1] - Sn.data()[j + 1]) * (vertex_pos[start + 1] - Sn.data()[j + 1]) +
				(vertex_pos[start + 2] - Sn.data()[j + 2]) * (vertex_pos[start + 2] - Sn.data()[j + 2]));
	}
	return 0.5 / (time_step * time_step) * energy; // 
}
