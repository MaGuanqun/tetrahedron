#pragma once
#include"external/Eigen/Dense"
#include"basic/eigenDenseOperation.h"
#include<iostream>
namespace test {
	void compute_f(Vector3d& x, Vector3d& x0, Vector3d& f)
	{
		double t1 = x.dot(x - x0) / x.squaredNorm();
		f = x * t1;
	}


	void compute_grad_f(Vector3d& x, Vector3d& x0, Matrix3d& grad_f)
	{
		grad_f = (1.0 - x.dot(x0) / x.squaredNorm()) * Matrix3d::Identity();
		Vector3d part = x0 / x.squaredNorm() - 2.0 * x.dot(x0) / (x.squaredNorm() * x.squaredNorm()) * x;
		grad_f -= x * part.transpose();
	}

	double TestDerivative(Vector3d& x, Vector3d& x0, Matrix3d& grad_f, Matrix3d& grac_f_approx)
	{
		double mi = 2e-8;
		Vector3d x_x = x;
		Vector3d x_y = x;
		Vector3d x_z = x;
		x_x[0] += mi;
		x_y[1] += mi;
		x_z[2] += mi;
		Vector3d f;
		compute_f(x, x0, f);
		compute_grad_f(x, x0, grad_f);

		Vector3d f1, f2, f3;
		compute_f(x_x, x0, f1);
		compute_f(x_y, x0, f2);
		compute_f(x_z, x0, f3);
		grac_f_approx.col(0) = (f1 - f) / mi;
		grac_f_approx.col(1) = (f2 - f) / mi;
		grac_f_approx.col(2) = (f3 - f) / mi;

		return (grad_f - grac_f_approx).squaredNorm();

	}

	void testDampDerivative()
	{
		int size = 1000;
		int size_half = size/2;
		int size_de = size / 3;
		Matrix3d grad_f, grac_f_approx;
		for (unsigned int i = 0; i < 100000; ++i) {
			Vector3d x = Vector3d((double)(size_half - rand() % size) / (double)size_de, (double)(size_half - rand() % size) / (double)size_de, (double)(size_half - rand() % size) / (double)size_de);
			Vector3d x0 = Vector3d((double)(size_half - rand() % size) / (double)size_de, (double)(size_half - rand() % size) / (double)size_de, (double)(size_half - rand() % size) / (double)size_de);
			if ((x - x0).squaredNorm() > 1e-16) {
				double error = TestDerivative(x, x0, grad_f, grac_f_approx);
				if (error > 1e-12) {
					std::cout << sqrt(error) << std::endl;
					//std::cout << x.transpose() << std::endl;
					//std::cout << x0.transpose() << std::endl;
					//std::cout << grad_f << std::endl;
					//std::cout << grac_f_approx << std::endl;

				}
			}
		}
	}








};

