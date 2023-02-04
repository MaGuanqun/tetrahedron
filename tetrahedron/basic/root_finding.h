#pragma once

#include"global.h"
#include"poly_root_finder.h"

namespace find_root {
	template<class T>
	inline bool quad(T a, T b, T c, T& t0, T& t1)
	{
		int sign = 1;

		if (b < 0) {
			sign = -1;
		}
		T D = b * b - 4.0 * a * c;
		if (D >= 0) {
			T q = -0.5 * (b + sign * sqrt(D)); //here we try to use sum for two positive values, avoid use subtract
			t0 = q / a;
			t1 = c / q;

			if (t0 > 1.0 || t0 < 0.0) {
				t0 = 1.0;
			}
			if (t1 > 1.0 || t1 < 0.0) {
				t1 = 1.0;
			}
			if (t0 > t1) {
				a = t1;
				t1 = t0;
				t0 = a;
			}
			if (t0 < 1.0) {
				return true;
			}
		}
		return false;
	}

	template<class T>
	inline void sortABC(T& a, T& b, T& c)
	{
		T d;
		if (a > b) { d = a; a = b; b = d; }
		if (a > c) { d = a; a = c; c = d; }
		if (b > c) { d = b; b = c; c = d; }
	}


	template <class T>
	inline bool getSmallestPositiveRealCubicRootTest(T a, T b, T c, T d, T& t, T tol)
	{
		// return negative value if no positive real root is found
		t = -1;
		std::complex<T> i(0, 1);
		std::complex<T> delta0(b * b - 3 * a * c, 0);
		std::complex<T> delta1(2 * b * b * b - 9 * a * b * c + 27 * a * a * d, 0);
		std::complex<T> C = std::pow((delta1 + std::sqrt(delta1 * delta1 - 4.0 * delta0 * delta0 * delta0)) / 2.0, 1.0 / 3.0);
		if (std::abs(C) == 0.0) {
			// a corner case listed by wikipedia found by our collaborate from another project
			C = pow((delta1 - sqrt(delta1 * delta1 - 4.0 * delta0 * delta0 * delta0)) / 2.0, 1.0 / 3.0);
		}
		std::complex<T> u2 = (-1.0 + sqrt(3.0) * i) / 2.0;
		std::complex<T> u3 = (-1.0 - sqrt(3.0) * i) / 2.0;
		std::complex<T> t1 = (b + C + delta0 / C) / (-3.0 * a);
		std::complex<T> t2 = (b + u2 * C + delta0 / (u2 * C)) / (-3.0 * a);
		std::complex<T> t3 = (b + u3 * C + delta0 / (u3 * C)) / (-3.0 * a);

		std::cout << t1 << " " << t2 << " " << t3 << std::endl;


		if ((std::abs(std::imag(t1)) < tol) && (std::real(t1) > 0))
			t = real(t1);
		if ((std::abs(imag(t2)) < tol) && (std::real(t2) > 0) && ((std::real(t2) < t) || (t < 0)))
			t = real(t2);
		if ((std::abs(imag(t3)) < tol) && (std::real(t3) > 0) && ((std::real(t3) < t) || (t < 0)))
			t = real(t3);

		if (t <= 1.0 && t > 0.0) {
			return true;
		}
		return false;
	}

	template <class T>
	inline bool getSmallestPositiveRealCubicRootWithPoly(T* op, T& t, T tol)
	{
		// return negative value if no positive real root is found
		T zeror[3] = { 2.0,2.0,2.0 }; T zeroi[3] = {0.0,0.0,0.0};
		int degree = poly_root::rpoly(op, 3, zeror, zeroi);
		t = 2.0;


		if ((std::abs(zeroi[0]) < tol) && (zeror[0] > 0))
			t = zeror[0];
		if ((std::abs(zeroi[1]) < tol) && (zeror[1] > 0) && ((zeror[1] < t) || (t < 0)))
			t = zeror[1];
		if ((std::abs(zeroi[2]) < tol) && (zeror[2] > 0) && ((zeror[2] < t) || (t < 0)))
			t = zeror[2];

		if (t <= 1.0 && t > 0.0) {
			return true;
		}
		return false;
	}

	template <class T>
	inline bool getSmallestPositiveRealCubicRoot(T a, T b, T c, T d,  T& t, T tol)
	{
		// return negative value if no positive real root is found
		t = -1;
		std::complex<T> i(0, 1);
		std::complex<T> delta0(b * b - 3 * a * c, 0);
		std::complex<T> delta1(2 * b * b * b - 9 * a * b * c + 27 * a * a * d, 0);
		std::complex<T> C = std::pow((delta1 + std::sqrt(delta1 * delta1 - 4.0 * delta0 * delta0 * delta0)) / 2.0, 1.0 / 3.0);
		if (std::abs(C) == 0.0) {
			// a corner case listed by wikipedia found by our collaborate from another project
			C = pow((delta1 - sqrt(delta1 * delta1 - 4.0 * delta0 * delta0 * delta0)) / 2.0, 1.0 / 3.0);
		}
		std::complex<T> u2 = (-1.0 + sqrt(3.0) * i) / 2.0;
		std::complex<T> u3 = (-1.0 - sqrt(3.0) * i) / 2.0;
		std::complex<T> t1 = (b + C + delta0 / C) / (-3.0 * a);
		std::complex<T> t2 = (b + u2 * C + delta0 / (u2 * C)) / (-3.0 * a);
		std::complex<T> t3 = (b + u3 * C + delta0 / (u3 * C)) / (-3.0 * a);

		if ((std::abs(std::imag(t1)) < tol) && (std::real(t1) > 0))
			t = real(t1);
		if ((std::abs(imag(t2)) < tol) && (std::real(t2) > 0) && ((std::real(t2) < t) || (t < 0)))
			t = real(t2);
		if ((std::abs(imag(t3)) < tol) && (std::real(t3) > 0) && ((std::real(t3) < t) || (t < 0)))
			t = real(t3);
		
		if (t <= 1.0 && t > 0.0) {
			return true;
		}
		return false;
	}

	//make sure the solution is right
	template <class T>
	inline void checkSolution(T a, T b, T c, T d, T& time)
	{
		while (d+time*(c+time*(b+a*time))<=0)
		{
			time *= 0.5;
		}
	}



	template<class T>
	inline bool solveCubicEquation(T a, T b, T c, T d, T& t0, T& t1, T& t2)
	{
		t0 = 1.0;
		t1 = 1.0;
		t2 = 1.0;

		b = b / a;
		c = c / a;
		d = d / a;

		T alpha = (-2.0 * b * b * b + 9.0 * b * c - 27.0 * d) / 54.0;
		T beta = (b * b - 3.0 * c) / 9.0;
		T dum1 = beta * beta * beta;
		T delta = alpha * alpha - dum1;
		b /= 3.0;

		if (delta > 0) {
			delta = sqrt(delta);
			t0 = -b + cbrt(alpha + delta) + cbrt(alpha - delta);
			if (t0 >= 0.0 && t0 <= 1.0) {
				return true;
			}
			else {
				t0 = 1.0;
			}
		}
		else if (delta == 0.0) {	//all roots are real, at least two are equal.
			alpha = cbrt(alpha);
			if ((-b + alpha + alpha) >= 0.0 && (-b + alpha + alpha) <= 1.0)
			{
				t0 = -b + alpha + alpha;
			}
			if ((alpha + b) <= 0.0 && (alpha + b) >= -1.0)
			{
				t1 = -alpha - b;
			}
			if (t1 < t0) {
				d = t1;
				t1 = t0;
				t0 = d;
			}
			if (t0 < 1.0) {
				return true;
			}
		}
		else { //all roots are real and different
			dum1 = acos(alpha / sqrt(dum1));
			alpha = 2.0 * std::sqrt(beta);
			t0 = -b + alpha * cos(dum1 / 3.0);
			t1 = -b + alpha * cos((dum1 + 2.0 * M_PI) / 3.0);
			t2 = -b + alpha * cos((dum1 + 4.0 * M_PI) / 3.0);

			if (t0 < 0.0 || t0>1.0) {
				t0 = 1.0;
			}
			if (t1 < 0.0 || t1>1.0) {
				t1 = 1.0;
			}
			if (t2 < 0.0 || t2>1.0) {
				t2 = 1.0;
			}
			sortABC(t0, t1, t2);
			if (t0 < 1.0) {
				return true;
			}
		}
		return false;
	}
}



