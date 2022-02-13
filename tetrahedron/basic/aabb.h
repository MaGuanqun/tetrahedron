#pragma once
#include"global.h"
namespace AABB {
	inline bool AABB_intersection(double* a, double* b)
	{
		if (a[3] < b[0]) {
			return false;
		}
		if (a[0] > b[3]) {
			return false;
		}
		if (a[4] < b[1]) {
			return false;
		}
		if (a[1] > b[4]) {
			return false;
		}
		if (a[5] < b[2]) {
			return false;
		}
		if (a[2] > b[5]) {
			return false;
		}
		return true;
	}

	//void obtainAABB(double* result, double* a, double* b, double* c)
	//{
	//	for (int i = 0; i < 3; ++i) {
	//		result[i] = myMin(a[i], b[i]);
	//		if (result[i] > c[i]) {
	//			result[i] = c[i];
	//		}
	//	}
	//	for (int i = 3; i < 6; ++i) {
	//		result[i] = myMax(a[i], b[i]);
	//		if (result[i] < c[i]) {
	//			result[i] = c[i];
	//		}
	//	}
	//}
	
	// a and b are two AABBs.
	inline void getAABB(double* result, double* a, double* b)//
	{
		for (unsigned int i = 0; i < 3; ++i) {
			if (a[i] < b[i]) {
				result[i] = a[i];
			}
			else {
				result[i] = b[i];
			}
		}
		for (unsigned int i = 3; i < 6; ++i) 
		{
			if (a[i] < b[i]) {
				result[i] = b[i];
			}
			else {
				result[i] = a[i];
			}
		}
	}

	// a and b are two positions.
	inline void obtainAABB(double* result, double* a, double* b)//
	{
		for (int i = 0; i < 3; ++i) {
			if (a[i] < b[i]) {
				result[i] = a[i];
				result[i + 3] = b[i];
			}
			else {
				result[i] = b[i];
				result[i + 3] = a[i];
			}
		}
	}

	inline void obtainAABB(double* result, double* a, double* b, double tolerance)//
	{
		for (int i = 0; i < 3; ++i) {
			if (a[i] < b[i]) {
				result[i] = a[i] - tolerance;
				result[i + 3] = b[i] + tolerance;
			}
			else {
				result[i] = b[i] - tolerance;
				result[i + 3] = a[i] + tolerance;
			}
		}
	}
}

//struct AABB
//{
//	double min[3];
//	double max[3];
//
//	AABB& operator =(AABB& aabb)
//	{
//		memcpy(min, aabb.min, 24);
//		memcpy(max, aabb.max, 24);
//		return *this;
//	}
//	bool operator==(const AABB& rhs) 
//	{
//		for (unsigned int i = 0; i < 3; ++i) {
//			if (min[i] != rhs.min[i]) {
//				return false;
//			}
//			if (max[i] != rhs.max[i]) {
//				return false;
//			}
//		}
//		return true;
//	}
//
//	bool AABB_intersection(AABB& a2)
//	{
//		if (max[0] < a2.min[0]) {
//			return false;
//		}
//		if (min[0] > a2.max[0]) {
//			return false;
//		}
//		if (max[1] < a2.min[1]) {
//			return false;
//		}
//		if (min[1] > a2.max[1]) {
//			return false;
//		}
//		if (max[2] < a2.min[2]) {
//			return false;
//		}
//		if (min[2] > a2.max[2]) {
//			return false;
//		}
//		return true;
//	}
//	void obtainAABB(double* a, double* b, double* c)
//	{
//		for (int i = 0; i < 3; ++i) {
//			min[i] = myMin(a[i], b[i]);
//			max[i] = myMax(a[i], b[i]);
//			min[i] = myMin(min[i], c[i]);
//			max[i] = myMin(max[i], c[i]);
//		}
//	}
//	void obtainAABB(double* a, double* b)//
//	{
//		for (int i = 0; i < 3; ++i) {
//			min[i] = myMin(a[i], b[i]);// -tolerance;
//			max[i] = myMax(a[i], b[i]);// +tolerance;
//		}
//	}
//
//	void obtainAABB(double* a, double* b, double tolerance)//
//	{
//		for (int i = 0; i < 3; ++i) {
//			min[i] = myMin(a[i], b[i]) - tolerance;// 
//			max[i] = myMax(a[i], b[i]) + tolerance;//;
//		}
//	}
//
//};
