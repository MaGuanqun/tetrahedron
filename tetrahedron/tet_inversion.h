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
	inline bool TetInversionTest(T* a0, T* b0, T* c0, T* d0,
		T* a1, T* b1, T* c1, T* d1, T& time_) 
	{

		T bd0[3], cd0[3],BD[3],CD[3];
		T ad0[3], AD[3];

		SUB(BD, b1, d1);
		SUB(CD, c1, d1);
		SUB(AD, a1, d1);
		CROSS(ad0, BD, CD);
		if (DOT(ad0, AD) > 0) {
			return false;
		}

		SUB(bd0, b0, d0);
		SUB(cd0, c0, d0);
		BD[0] = b1[0] - d1[0] - bd0[0];
		BD[1] = b1[1] - d1[1] - bd0[1];
		BD[2] = b1[2] - d1[2] - bd0[2];

		CD[0] = c1[0] - d1[0] - cd0[0];
		CD[1] = c1[1] - d1[1] - cd0[1];
		CD[2] = c1[2] - d1[2] - cd0[2];

		T m0[3], m1[3], m2[3];
		CROSS(m0, bd0, cd0);
		CROSS(m2, BD, CD);
		CROSS_SUM(m1, bd0, CD, BD, cd0);
		

		SUB(ad0, a0, d0);
		AD[0] = a1[0] - d1[0] - ad0[0];
		AD[1] = a1[1] - d1[1] - ad0[1];
		AD[2] = a1[2] - d1[2] - ad0[2];

		T op[4]; //op[0]x^3+op[1]x^2+cx+d=0
		op[3] = DOT(ad0, m0);
		op[2] = DOT(ad0, m1) + DOT(AD, m0);
		op[1] = DOT(ad0, m2) + DOT(AD, m1);
		op[0] = DOT(AD, m2);

		int roots = 0;
		unsigned int reducedDegree = 3;
		for (unsigned int i = 0; i < 3; i++) {
			if (op[i] == 0)
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
			if (find_root::solveCubicEquation(op[0], op[1], op[2], op[3], time[0], time[1], time[2]) {
				time_= time[0];
				return true;
			}
			return false;
		}
		else if (reducedDegree == 2) {
			if (find_root::quad(op[0], op[1], op[2], time[0], time[1])) {
				time_=time[0];
				return true;
			}
			return false;
		}
		else if (reducedDegree == 1) {
			time[0] = -op[1] / op[0];
			if (time[0] > 0.0 && time[0]<1.0) {
				time_ = time[0];
				return true;
			}
			return false;
		}
		return false;
	}



}