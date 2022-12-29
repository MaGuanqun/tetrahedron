#pragma once
#include"basic/global.h"
#include"basic/root_finding.h"

namespace inversionTest {
	template <class T>
	inline bool notHaveRoots(T* op)
	{
		for (int i = 0; i < 3; ++i) {
			if (op[i] < 0) {
				return false;
			}
		}
		return true;
	}
	template <class T>
	inline T volume(T* a1, T* b1, T* c1, T* d1)
	{
		T bd0[3], cd0[3], BD[3], CD[3];
		T ad0[3], AD[3];
		SUB(BD, b1, a1);
		SUB(CD, c1, a1);
		SUB(AD, d1, a1);
		CROSS(ad0, BD, CD);
		return DOT(ad0, AD);
	}

	template <class T>
	inline bool TetInversionTest(T* a0, T* b0, T* c0, T* d0,
		T* a1, T* b1, T* c1, T* d1, T* time_) 
	{

		T ba0[3], ca0[3],BA[3],CA[3];
		T da0[3], DA[3];

		//SUB(BA, b1, a1);
		//SUB(CA, c1, a1);
		//SUB(DA, d1, a1);
		//CROSS(da0, BA, CA);
		//if (DOT(da0, DA) > 0) {
		//	return false;
		//}

		SUB(ba0, b0, a0);
		SUB(ca0, c0, a0);
		BA[0] = b1[0] - a1[0] - ba0[0];
		BA[1] = b1[1] - a1[1] - ba0[1];
		BA[2] = b1[2] - a1[2] - ba0[2];

		CA[0] = c1[0] - a1[0] - ca0[0];
		CA[1] = c1[1] - a1[1] - ca0[1];
		CA[2] = c1[2] - a1[2] - ca0[2];

		T m0[3], m1[3], m2[3];
		CROSS(m0, ba0, ca0);
		CROSS(m2, BA, CA);
		CROSS_SUM(m1, ba0, CA, BA, ca0);
		

		SUB(da0, d0, a0);
		DA[0] = d1[0] - a1[0] - da0[0];
		DA[1] = d1[1] - a1[1] - da0[1];
		DA[2] = d1[2] - a1[2] - da0[2];

		T op[4]; //op[0]x^3+op[1]x^2+cx+d=0
		op[0] = DOT(DA, m2);
		op[1] = DOT(da0, m2) + DOT(DA, m1);
		op[2] = DOT(da0, m1) + DOT(DA, m0);
		op[3] = DOT(da0, m0);

		int roots = 0;
		unsigned int reducedDegree = 3;
		for (unsigned int i = 0; i < 3; i++) {
			if (std::abs(op[i]) <  1e-8)
				reducedDegree--;
			else
				break;
		}
		if (reducedDegree < 3) {
			//move lower term coeff to higher term, this is for convenience
			for (int i = 0; i <= reducedDegree; i++) {
				op[i] = op[i + 3 - reducedDegree];
			}
		}

		T time[3];
		if (reducedDegree == 3) {
			if (find_root::getSmallestPositiveRealCubicRoot(op[0], op[1], op[2], op[3], time[0], 1e-8)) {
				*time_= time[0];
				return true;
			}
			return false;
		}
		else if (reducedDegree == 2) {
			if (find_root::quad(op[0], op[1], op[2], time[0], time[1])) {
				*time_=time[0];
				return true;
			}
			return false;
		}
		else if (reducedDegree == 1) {
			time[0] = -op[1] / op[0];
			if (time[0] > 0.0 && time[0]<1.0) {
				*time_ = time[0];
				return true;
			}
			return false;
		}
		return false;
	}

	template <class T>
	inline bool TetInversionTestTest(T* a0, T* b0, T* c0, T* d0,
		T* a1, T* b1, T* c1, T* d1, T* time_)
	{

		T ba0[3], ca0[3], BA[3], CA[3];
		T da0[3], DA[3];

		//SUB(BA, b1, a1);
		//SUB(CA, c1, a1);
		//SUB(DA, d1, a1);
		//CROSS(da0, BA, CA);
		//if (DOT(da0, DA) > 0) {
		//	return false;
		//}

		SUB(ba0, b0, a0);
		SUB(ca0, c0, a0);
		BA[0] = b1[0] - a1[0] - ba0[0];
		BA[1] = b1[1] - a1[1] - ba0[1];
		BA[2] = b1[2] - a1[2] - ba0[2];

		CA[0] = c1[0] - a1[0] - ca0[0];
		CA[1] = c1[1] - a1[1] - ca0[1];
		CA[2] = c1[2] - a1[2] - ca0[2];

		T m0[3], m1[3], m2[3];
		CROSS(m0, ba0, ca0);
		CROSS(m2, BA, CA);
		CROSS_SUM(m1, ba0, CA, BA, ca0);


		SUB(da0, d0, a0);
		DA[0] = d1[0] - a1[0] - da0[0];
		DA[1] = d1[1] - a1[1] - da0[1];
		DA[2] = d1[2] - a1[2] - da0[2];

		T op[4]; //op[0]x^3+op[1]x^2+cx+d=0
		op[0] = DOT(DA, m2);
		op[1] = DOT(da0, m2) + DOT(DA, m1);
		op[2] = DOT(da0, m1) + DOT(DA, m0);
		op[3] = DOT(da0, m0);

		std::cout << "coe op " << op[0] << " " << op[1] << " " << op[2] << " " << op[3] << std::endl;

		int roots = 0;
		unsigned int reducedDegree = 3;
		for (unsigned int i = 0; i < 3; i++) {
			if (op[i] < 1e-10)
				reducedDegree--;
			else
				break;
		}
		if (reducedDegree < 3) {
			//move lower term coeff to higher term, this is for convenience
			for (int i = 0; i <= reducedDegree; i++) {
				op[i] = op[i + 3 - reducedDegree];
			}
		}

		T time[3];
		if (reducedDegree == 3) {
			if (find_root::solveCubicEquation(op[0], op[1], op[2], op[3], time[0], time[1], time[2])) {
				*time_ = time[0];
				return true;
			}
			return false;
		}
		else if (reducedDegree == 2) {
			if (find_root::quad(op[0], op[1], op[2], time[0], time[1])) {
				*time_ = time[0];
				return true;
			}
			return false;
		}
		else if (reducedDegree == 1) {
			time[0] = -op[1] / op[0];
			if (time[0] > 0.0 && time[0] < 1.0) {
				*time_ = time[0];
				return true;
			}
			return false;
		}
		return false;
	}


}