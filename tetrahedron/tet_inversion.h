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
		op[3] = 0.8*DOT(da0, m0);

		T ori_op[4];
		memcpy(ori_op, op, 4 * sizeof(T));

		int roots = 0;
		unsigned int reducedDegree = 3;
		for (unsigned int i = 0; i < 3; i++) {
			if (std::abs(op[i]) < 1e-20)
				reducedDegree--;
			else
				break;
		}

		if (reducedDegree < 3) {
			//move lower term coeff to higher term, this is for convenience
			for (int i = 0; i <= reducedDegree; i++) {
				op[i] = op[i +3- reducedDegree];
			}
		}

		T time[3];
		if (reducedDegree == 3) {
			if (find_root::getSmallestPositiveRealCubicRootWithPoly(op, time[0], 1e-8)) {
				*time_= time[0];
				//find_root::checkSolution(op[0], op[1], op[2], op[3], *time_);
				return true;
			}
			return false;
		}
		else if (reducedDegree == 2) {
			if (find_root::quad(op[0], op[1], op[2], time[0], time[1])) {
				*time_=time[0];
				find_root::checkSolution(ori_op[0], ori_op[1], ori_op[2], ori_op[3], *time_);
				return true;
			}
			return false;
		}
		else if (reducedDegree == 1) {
			time[0] = -op[1] / op[0];
			if (time[0] > 0.0 && time[0]<1.0) {
				*time_ = time[0];
				find_root::checkSolution(ori_op[0], ori_op[1], ori_op[2], ori_op[3], *time_);
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
		op[3] =0.8* DOT(da0, m0);


		T ori_op[4];
		memcpy(ori_op, op, 4 * sizeof(T));

		std::cout << "coe op " << op[0] << " " << op[1] << " " << op[2] << " " << op[3] << std::endl;

		T r[3]; T g[3];

		RootFinder root_finder;
		root_finder.rpoly(op, 3, r, g);
		std::cout<<"no disdegree " << r[0] << " " << r[1] << " " << r[2] << std::endl;
		std::cout << g[0] << " " << g[1] << " " << g[2] << std::endl;


		int roots = 0;
		unsigned int reducedDegree = 3;
		for (unsigned int i = 0; i < 3; i++) {
			if (std::abs(op[i]) < 1e-20)
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

		std::cout << "coe op final " << reducedDegree << " " << op[0] << " " << op[1] << " " << op[2] << " " << op[3] << std::endl;
		T time[3];
		if (reducedDegree == 3) {

			RootFinder root_finder;
			root_finder.rpoly(op, 3, r, g);
			std::cout << r[0] << " " << r[1] << " " << r[2] << std::endl;
			std::cout <<g[0] << " " << g[1] << " " << g[2] << std::endl;

			
			std::cout << time[0] << std::endl;
			if (find_root::getSmallestPositiveRealCubicRootWithPoly(op, time[0], 1e-8)) {
				*time_ = time[0];
				//find_root::checkSolution(op[0], op[1], op[2], op[3], *time_);
				return true;
			}
			return false;
		}
		else if (reducedDegree == 2) {
			if (find_root::quad(op[0], op[1], op[2], time[0], time[1])) {
				*time_ = time[0];
				find_root::checkSolution(ori_op[0], ori_op[1], ori_op[2], ori_op[3], *time_);
				return true;
			}
			return false;
		}
		else if (reducedDegree == 1) {
			time[0] = -op[1] / op[0];
			if (time[0] > 0.0 && time[0] < 1.0) {
				*time_ = time[0];
				find_root::checkSolution(ori_op[0], ori_op[1], ori_op[2], ori_op[3], *time_);
				return true;
			}
			return false;
		}
		return false;
	}

	template <class T>
	inline bool TetInversionTestTestCompare(T* a0, T* b0, T* c0, T* d0,
		T* a1, T* b1, T* c1, T* d1, T* time_)
	{
		T x1 = a0[0];
		T x2 = b0[0];
		T x3 = c0[0];
		T x4 = d0[0];

		T y1 = a0[1];
		T y2 = b0[1];
		T y3 = c0[1];
		T y4 = d0[1];

		T z1 = a0[2];
		T z2 = b0[2];
		T z3 = c0[2];
		T z4 = d0[2];

		T p1 = a1[0] - a0[0];
		T p2 = b1[0] - b0[0];
		T p3 = c1[0] - c0[0];
		T p4 = d1[0] - d0[0];


		T q1 = a1[1] - a0[1];
		T q2 = b1[1] - b0[1];
		T q3 = c1[1] - c0[1];
		T q4 = d1[1] - d0[1];

		T r1 = a1[2] - a0[2];
		T r2 = b1[2] - b0[2];
		T r3 = c1[2] - c0[2];
		T r4 = d1[2] - d0[2];



		T a = -p1 * q2 * r3 + p1 * r2 * q3 + q1 * p2 * r3 - q1 * r2 * p3 - r1 * p2 * q3 + r1 * q2 * p3 + p1 * q2 * r4 - p1 * r2 * q4 - q1 * p2 * r4 + q1 * r2 * p4 + r1 * p2 * q4 - r1 * q2 * p4 - p1 * q3 * r4 + p1 * r3 * q4 + q1 * p3 * r4 - q1 * r3 * p4 - r1 * p3 * q4 + r1 * q3 * p4 + p2 * q3 * r4 - p2 * r3 * q4 - q2 * p3 * r4 + q2 * r3 * p4 + r2 * p3 * q4 - r2 * q3 * p4;
		T b = -x1 * q2 * r3 + x1 * r2 * q3 + y1 * p2 * r3 - y1 * r2 * p3 - z1 * p2 * q3 + z1 * q2 * p3 + x2 * q1 * r3 - x2 * r1 * q3 - y2 * p1 * r3 + y2 * r1 * p3 + z2 * p1 * q3 - z2 * q1 * p3 - x3 * q1 * r2 + x3 * r1 * q2 + y3 * p1 * r2 - y3 * r1 * p2 - z3 * p1 * q2 + z3 * q1 * p2 + x1 * q2 * r4 - x1 * r2 * q4 - y1 * p2 * r4 + y1 * r2 * p4 + z1 * p2 * q4 - z1 * q2 * p4 - x2 * q1 * r4 + x2 * r1 * q4 + y2 * p1 * r4 - y2 * r1 * p4 - z2 * p1 * q4 + z2 * q1 * p4 + x4 * q1 * r2 - x4 * r1 * q2 - y4 * p1 * r2 + y4 * r1 * p2 + z4 * p1 * q2 - z4 * q1 * p2 - x1 * q3 * r4 + x1 * r3 * q4 + y1 * p3 * r4 - y1 * r3 * p4 - z1 * p3 * q4 + z1 * q3 * p4 + x3 * q1 * r4 - x3 * r1 * q4 - y3 * p1 * r4 + y3 * r1 * p4 + z3 * p1 * q4 - z3 * q1 * p4 - x4 * q1 * r3 + x4 * r1 * q3 + y4 * p1 * r3 - y4 * r1 * p3 - z4 * p1 * q3 + z4 * q1 * p3 + x2 * q3 * r4 - x2 * r3 * q4 - y2 * p3 * r4 + y2 * r3 * p4 + z2 * p3 * q4 - z2 * q3 * p4 - x3 * q2 * r4 + x3 * r2 * q4 + y3 * p2 * r4 - y3 * r2 * p4 - z3 * p2 * q4 + z3 * q2 * p4 + x4 * q2 * r3 - x4 * r2 * q3 - y4 * p2 * r3 + y4 * r2 * p3 + z4 * p2 * q3 - z4 * q2 * p3;
		T c = -x1 * y2 * r3 + x1 * z2 * q3 + x1 * y3 * r2 - x1 * z3 * q2 + y1 * x2 * r3 - y1 * z2 * p3 - y1 * x3 * r2 + y1 * z3 * p2 - z1 * x2 * q3 + z1 * y2 * p3 + z1 * x3 * q2 - z1 * y3 * p2 - x2 * y3 * r1 + x2 * z3 * q1 + y2 * x3 * r1 - y2 * z3 * p1 - z2 * x3 * q1 + z2 * y3 * p1 + x1 * y2 * r4 - x1 * z2 * q4 - x1 * y4 * r2 + x1 * z4 * q2 - y1 * x2 * r4 + y1 * z2 * p4 + y1 * x4 * r2 - y1 * z4 * p2 + z1 * x2 * q4 - z1 * y2 * p4 - z1 * x4 * q2 + z1 * y4 * p2 + x2 * y4 * r1 - x2 * z4 * q1 - y2 * x4 * r1 + y2 * z4 * p1 + z2 * x4 * q1 - z2 * y4 * p1 - x1 * y3 * r4 + x1 * z3 * q4 + x1 * y4 * r3 - x1 * z4 * q3 + y1 * x3 * r4 - y1 * z3 * p4 - y1 * x4 * r3 + y1 * z4 * p3 - z1 * x3 * q4 + z1 * y3 * p4 + z1 * x4 * q3 - z1 * y4 * p3 - x3 * y4 * r1 + x3 * z4 * q1 + y3 * x4 * r1 - y3 * z4 * p1 - z3 * x4 * q1 + z3 * y4 * p1 + x2 * y3 * r4 - x2 * z3 * q4 - x2 * y4 * r3 + x2 * z4 * q3 - y2 * x3 * r4 + y2 * z3 * p4 + y2 * x4 * r3 - y2 * z4 * p3 + z2 * x3 * q4 - z2 * y3 * p4 - z2 * x4 * q3 + z2 * y4 * p3 + x3 * y4 * r2 - x3 * z4 * q2 - y3 * x4 * r2 + y3 * z4 * p2 + z3 * x4 * q2 - z3 * y4 * p2;
		T d = (1.0 - 0.2) * (x1 * z2 * y3 - x1 * y2 * z3 + y1 * x2 * z3 - y1 * z2 * x3 - z1 * x2 * y3 + z1 * y2 * x3 + x1 * y2 * z4 - x1 * z2 * y4 - y1 * x2 * z4 + y1 * z2 * x4 + z1 * x2 * y4 - z1 * y2 * x4 - x1 * y3 * z4 + x1 * z3 * y4 + y1 * x3 * z4 - y1 * z3 * x4 - z1 * x3 * y4 + z1 * y3 * x4 + x2 * y3 * z4 - x2 * z3 * y4 - y2 * x3 * z4 + y2 * z3 * x4 + z2 * x3 * y4 - z2 * y3 * x4);


		std::cout << "compare tet inverse  coe " << a << " " << b << " " << c << " " << d << std::endl;

		return true;

	}

}