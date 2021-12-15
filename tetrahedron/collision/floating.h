#pragma once
#include"../basic/global.h"
//we set this class for inside test, TightCCD

class floating {
public:
	double v;
	double sigma;
public:
	inline floating() :v(double(0)), sigma(double(0)) {}
	inline floating(double _v) : v(_v), sigma(double(0)) {}
	inline floating(double _v, double _sigma) : v(_v), sigma(_sigma) {
	}
	inline int sign() {
		if (v > sigma)
			return 1;
		else if (v < (-sigma))
			return -1;
		else
			return 0;
	}
	inline int sign() const {
		if (v > sigma)
			return 1;
		else if (v < (-sigma))
			return -1;
		else
			return 0;
	}
	inline bool operator == (const double t) {
		if (v == t)
			return true;
		return false;
	}
	inline floating& operator += (const floating& rhs) {
		v += rhs.v;
		sigma = sigma + rhs.sigma + fabs(v) * BOUND_UNIT;
		return *this;
	}
	inline floating& operator -= (const floating& rhs) {
		v -= rhs.v;
		sigma = sigma + rhs.sigma + fabs(v) * BOUND_UNIT;
		return *this;
	}
	inline floating& operator *= (const floating& rhs) {
		sigma = sigma * rhs.sigma + fabs(v) * rhs.sigma + fabs(rhs.v) * sigma;
		v *= rhs.v;
		sigma += fabs(v) * BOUND_UNIT;
		return *this;
	}
	inline floating operator + (const floating& other) const {
		double _v = v + other.v;
		double _sigma = sigma + other.sigma + fabs(_v) * BOUND_UNIT;
		return floating(_v, _sigma);
	}
	inline floating operator - (const floating& other) const {
		double _v = v - other.v;
		double _sigma = sigma + other.sigma + fabs(_v) * BOUND_UNIT;
		return floating(_v, _sigma);
	}
	inline floating operator * (const floating& other) const {
		double _v = v * other.v;
		double _sigma = sigma * other.sigma + fabs(v) * other.sigma + fabs(other.v) * sigma + fabs(_v) * BOUND_UNIT;
		return floating(_v, _sigma);
	}
	inline floating operator - () const
	{
		return floating(-v, sigma);
	}
};
