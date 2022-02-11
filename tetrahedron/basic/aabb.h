#pragma once
#include"global.h"

struct AABB
{
	double min[3];
	double max[3];

	AABB& operator =(AABB& aabb)
	{
		memcpy(min, aabb.min, 24);
		memcpy(max, aabb.max, 24);
		return *this;
	}
	bool operator==(const AABB& rhs) 
	{
		for (unsigned int i = 0; i < 3; ++i) {
			if (min[i] != rhs.min[i]) {
				return false;
			}
			if (max[i] != rhs.max[i]) {
				return false;
			}
		}
		return true;
	}

	bool AABB_intersection(AABB& a2)
	{
		if (max[0] < a2.min[0]) {
			return false;
		}
		if (min[0] > a2.max[0]) {
			return false;
		}
		if (max[1] < a2.min[1]) {
			return false;
		}
		if (min[1] > a2.max[1]) {
			return false;
		}
		if (max[2] < a2.min[2]) {
			return false;
		}
		if (min[2] > a2.max[2]) {
			return false;
		}
		return true;
	}
	void obtainAABB(double* a, double* b, double* c)
	{
		for (int i = 0; i < 3; ++i) {
			min[i] = myMin(a[i], b[i]);
			max[i] = myMax(a[i], b[i]);
			min[i] = myMin(min[i], c[i]);
			max[i] = myMin(max[i], c[i]);
		}
	}
	void obtainAABB(double* a, double* b)//
	{
		for (int i = 0; i < 3; ++i) {
			min[i] = myMin(a[i], b[i]);// -tolerance;
			max[i] = myMax(a[i], b[i]);// +tolerance;
		}
	}

	void obtainAABB(double* a, double* b, double tolerance)//
	{
		for (int i = 0; i < 3; ++i) {
			min[i] = myMin(a[i], b[i]) - tolerance;// 
			max[i] = myMax(a[i], b[i]) + tolerance;//;
		}
	}

};
