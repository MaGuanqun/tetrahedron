#pragma once
#include"../basic/global.h"
class InsideTest
{
public:
	bool insideTest(double* a0, double* b0, double* c0, double* v0,
		double* a1, double* b1, double* c1, double* v1, double* n0, double* n1, double* cross_for_CCD, bool ee_test);

private:
	struct BC {
		double k0, k1, k2, k3;
		double kk0, kk1, kk2;
		int ct;
	};
	void getBezier4(double* a0, double* b0, double* c0, double* v0,
		double* a1, double* b1, double* c1, double* v1, double& l0, double& l1, double& l2, double& l3, double& l4,
		double* n0, double* n1, double* cross_for_CCD, int which, bool ee_test);
	void getSimplifyed(double& k0, double& k1, double& k2, double& k3, double& l0, double& l1, double& l2, double& l3, double& l4,
		double& j0, double& j1, double& j2);

	void getSigns(const double& t0, const double& t1, const BC& c, double& lt0, double& lt1);

	inline void norm(double* result, double* p1, double* p2, double* p3) {
		double temp1[3], temp2[3];
		SUB(temp1, p2, p1);
		SUB(temp2, p3, p1);
		CROSS(result, temp1, temp2);
	}

	double _evaluateBezier(const double& p0, const double& p1, const double& p2, const double& p3, const double& t, const double& s);
	double _evaluateBezier2(const double& p0, const double& p1, const double& p2, const double& t, const double& s);

	bool bezierDecomposition(const double& k0, const double& k1, const double& k2, const double& k3,
		const double& j0, const double& j1, const double& j2,
		double& m0, double& m1, double& n0, double& n1);
	bool getBezier3(double* a0, double* b0, double* c0, double* v0,
		double* a1, double* b1, double* c1, double* v1, double& p0, double& p1, double& p2, double& p3, 
		double* n0, double* n1, double* cross_for_CCD);

	int bezierClassification(const double& k0, const double& k1, const double& k2, const double& k3,
		double& kk0, double& kk1, double& kk2);

};



