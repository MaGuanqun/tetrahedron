#pragma once
#include"../external/Eigen/Dense"
#include"../external/Eigen/Sparse"
using namespace Eigen;

#ifndef EIGEN_DENSE_OPERATION
#define EIGEN_DENSE_OPERATION
namespace denseOperation {

	inline double dotProduct(const Vector3d& a, const Vector3d& b) {
		return a.data()[0] * b.data()[0] + a.data()[1] * b.data()[1] + a.data()[2] * b.data()[2];
	}
	inline double dotProduct(double* a, double* b) {
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}
	inline Vector3d crossProduct(const Vector3d& a, const Vector3d& b) {
		return Vector3d(a.data()[1] * b.data()[2] - a.data()[2] * b.data()[1], a.data()[2] * b.data()[0] - a.data()[0] * b.data()[2], a.data()[0] * b.data()[1] - a.data()[1] * b.data()[0]);
	}
	inline Vector3d crossProduct(double* a, double* b) {
		return Vector3d(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
	}
	inline double norm(const Vector3d& a) {
		return sqrt(a.data()[0] * a.data()[0] + a.data()[1] * a.data()[1] + a.data()[2] * a.data()[2]);
	}
	inline double norm2(const Vector3d& a) {
		return a.data()[0] * a.data()[0] + a.data()[1] * a.data()[1] + a.data()[2] * a.data()[2];
	}
	inline double norm(double* a) {
		return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	}
	inline Vector3d normalized(const Vector3d& a) {
		double norm = sqrt(a.data()[0] * a.data()[0] + a.data()[1] * a.data()[1] + a.data()[2] * a.data()[2]);
		return Vector3d(a.data()[0] / norm, a.data()[1] / norm, a.data()[2] / norm);
	}
	inline Vector3d normalized(double* a) {
		double norm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
		return Vector3d(a[0] / norm, a[1] / norm, a[2] / norm);
	}
	inline double dotProductX(const VectorXd& a, const VectorXd& b) {
		size_t size = a.size();
		double dot = 0.0;
		for (size_t i = 0; i < size; i++) {
			dot += a.data()[i] * b.data()[i];
		}
		return dot;
	}
	inline double dotProductX_(VectorXd& a, VectorXd& b) {
		size_t size = a.size();
		double dot = 0.0;
		for (size_t i = 0; i < size; i++) {
			dot += a.data()[i] * b.data()[i];
		}
		return dot;
	}
	inline VectorXd cwiseInverse(VectorXd& a) {
		size_t size = a.size();
		for (size_t it=0; it<size; it++) {
			a.data()[it]=1.0/a.data()[it];
		}
		return a;
	}
	inline VectorXd cwiseProduct(VectorXd& a, VectorXd& b) {
		size_t size = a.size();
		for (size_t it = 0; it < size; it++) {
			a.data()[it] = a.data()[it] * b.data()[it];
		}
		return a;
	}
	inline VectorXd cwiseProduct_(VectorXd& a, VectorXd& b) {
		size_t size = a.size();
		VectorXd c(size);		
		for (size_t it = 0; it < size; it++) {
			c.data()[it] = a.data()[it] * b.data()[it];
		}
		return c;
	}
	inline Matrix3d ab_transpose(Vector3d& a, Vector3d& b)//a*b^T 	
	{
		Matrix3d m;
		m.data()[0] = a.data()[0] * b.data()[0];
		m.data()[1] = a.data()[1] * b.data()[0];
		m.data()[2] = a.data()[2] * b.data()[0];
		m.data()[3] = a.data()[0] * b.data()[1];
		m.data()[4] = a.data()[1] * b.data()[1];
		m.data()[5] = a.data()[2] * b.data()[1];
		m.data()[6] = a.data()[0] * b.data()[2];
		m.data()[7] = a.data()[1] * b.data()[2];
		m.data()[8] = a.data()[2] * b.data()[2];
		return m;
	}
	inline Matrix3d ab_transpose(double* a, double* b)//a*b^T 	
	{
		Matrix3d m;
		m.data()[0] = a[0] * b[0];
		m.data()[1] = a[1] * b[0];
		m.data()[2] = a[2] * b[0];
		m.data()[3] = a[0] * b[1];
		m.data()[4] = a[1] * b[1];
		m.data()[5] = a[2] * b[1];
		m.data()[6] = a[0] * b[2];
		m.data()[7] = a[1] * b[2];
		m.data()[8] = a[2] * b[2];
		return m;
	}
	inline Matrix3d transpose(Matrix3d& a) {
		Matrix3d m;
		m.data()[0] = a.data()[0];
		m.data()[4] = a.data()[4];
		m.data()[8] = a.data()[8];
		m.data()[1] = a.data()[3];
		m.data()[2] = a.data()[6];
		m.data()[3] = a.data()[1];
		m.data()[5] = a.data()[7];
		m.data()[6] = a.data()[2];
		m.data()[7] = a.data()[5];
		return m;
	}
	inline Matrix3d transpose(double* a) {
		Matrix3d m;
		m.data()[0] = a[0];
		m.data()[4] = a[4];
		m.data()[8] = a[8];
		m.data()[1] = a[3];
		m.data()[2] = a[6];
		m.data()[3] = a[1];
		m.data()[5] = a[7];
		m.data()[6] = a[2];
		m.data()[7] = a[5];
		return m;
	}
	inline Matrix3d inversion(Matrix3d& a) {
		Matrix3d m;
		double det1 = 1.0/(a.data()[0] * (a.data()[4] * a.data()[8] - a.data()[7] * a.data()[5]) - a.data()[3] * (a.data()[1] * a.data()[8] - a.data()[7] * a.data()[2]) + a.data()[6] * (a.data()[1] * a.data()[5] - a.data()[4] * a.data()[2]));
		m.data()[0] = det1 * (a.data()[4]* a.data()[8]- a.data()[7]* a.data()[5]);
		m.data()[1] = det1 * (a.data()[7]* a.data()[2]- a.data()[1]* a.data()[8]);
		m.data()[2] = det1 * (a.data()[1]* a.data()[5]- a.data()[4]* a.data()[2]);
		m.data()[3] = det1 * (a.data()[6]* a.data()[5]- a.data()[3]* a.data()[8]);
		m.data()[4] = det1 * (a.data()[0]* a.data()[8]- a.data()[6]* a.data()[2]);
		m.data()[5] = det1 * (a.data()[3]* a.data()[2]- a.data()[0]* a.data()[5]);
		m.data()[6] = det1 * (a.data()[3]* a.data()[7]- a.data()[6]* a.data()[4]);
		m.data()[7] = det1 * (a.data()[6]* a.data()[1]- a.data()[7]* a.data()[0]);
		m.data()[8] = det1 * (a.data()[4]* a.data()[0]- a.data()[3]* a.data()[1]);
		return m;
	}
	inline Matrix3d inversion(double* a) {
		Matrix3d m;
		double det1 = 1.0 / (a[0] * (a[4] * a[8] - a[7] * a[5]) - a[3] * (a[1] * a[8] - a[7] * a[2]) + a[6] * (a[1] * a[5] - a[4] * a[2]));
		m.data()[0] = det1 * (a[4] * a[8] - a[7] * a[5]);
		m.data()[1] = det1 * (a[7] * a[2] - a[1] * a[8]);
		m.data()[2] = det1 * (a[1] * a[5] - a[4] * a[2]);
		m.data()[3] = det1 * (a[6] * a[5] - a[3] * a[8]);
		m.data()[4] = det1 * (a[0] * a[8] - a[6] * a[2]);
		m.data()[5] = det1 * (a[3] * a[2] - a[0] * a[5]);
		m.data()[6] = det1 * (a[3] * a[7] - a[6] * a[4]);
		m.data()[7] = det1 * (a[6] * a[1] - a[7] * a[0]);
		m.data()[8] = det1 * (a[4] * a[0] - a[3] * a[1]);
		return m;
	}

	inline double det(Matrix3d& a) {
		double det = a.data()[0] * (a.data()[4] * a.data()[8] - a.data()[7] * a.data()[5]) - a.data()[3] * (a.data()[1] * a.data()[8] - a.data()[7] * a.data()[2]) + a.data()[6] * (a.data()[1] * a.data()[5] - a.data()[4] * a.data()[2]);
		return det;
	}

	inline double det(double* a) {
		double det = a[0] * (a[4] * a[8] - a[7] * a[5]) - a[3] * (a[1] * a[8] - a[7] * a[2]) + a[6] * (a[1] * a[5] - a[4] * a[2]);
		return det;
	}
}


#endif // !EIGEN_DENSE_OPERATION