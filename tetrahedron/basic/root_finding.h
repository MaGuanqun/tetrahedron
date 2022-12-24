#pragma once

#include"global.h"

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



	template<class T>
	inline bool solveCubicEquation(T a, T b, T c, T d, T& t0, T& t1, T& t2)
	{
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



