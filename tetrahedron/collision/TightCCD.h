#pragma once
#include"floating.h"

//#define BOUND_UNIT 0.00000000000000022204460492503136
enum ReturnValue {
	RETURN_ZERO,
	RETURN_FALSE,
	RETURN_TRUE
};


class TightCCD
{
public:
	bool insideTest(floating* a0, floating* b0, floating* c0, floating* v0,
		floating* a1, floating* b1, floating* c1, floating* v1, floating* n0, floating* n1, floating* cross_for_CCD, bool ee_test);
	bool insideTest(double* a0, double* b0, double* c0, double* v0,
		double* a1, double* b1, double* c1, double* v1, floating* n0, floating* n1, floating* cross_for_CCD, bool ee_test);

private:
	struct BC {
		floating k0, k1, k2, k3;
		floating kk0, kk1, kk2;
		int ct;
		int sign_k0, sign_k1, sign_k2, sign_k3;
		int sign_kk0, sign_kk1, sign_kk2;
		BC(const floating& _k0, const floating& _k1, const floating& _k2, const floating& _k3,
			const floating& _kk0, const floating& _kk1, const floating& _kk2, int _ct)
		{
			this->k0 = _k0;
			this->k1 = _k1;
			this->k2 = _k2;
			this->k3 = _k3;
			this->kk0 = _kk0;
			this->kk1 = _kk1;
			this->kk2 = _kk2;
			this->ct = _ct;
			sign_k0 = k0.sign();
			sign_k1 = k1.sign();
			sign_k2 = k2.sign();
			sign_k3 = k3.sign();
			sign_kk0 = kk0.sign();
			sign_kk1 = kk1.sign();
			sign_kk2 = kk2.sign();
		}
	};
	void getBezier4(floating* a0, floating* b0, floating* c0, floating* v0,
		floating* a1, floating* b1, floating* c1, floating* v1,
		floating& l0, floating& l1, floating& l2, floating& l3, floating& l4,
		floating* n0, floating* n1, floating* deltaN, floating* nX,
		int which, bool ee_test);
	void getSimplifyed(floating& k0, floating& k1, floating& k2, floating& k3, floating& l0, floating& l1, floating& l2, floating& l3, floating& l4,
		floating& j0, floating& j1, floating& j2);

	bool getSigns(const floating& t0, const floating& t1, const BC& c, floating& lt0, floating& lt1, int root_nums);

	inline void norm(floating* result, floating* p1, floating* p2, floating* p3) {
		floating temp1[3], temp2[3];
		SUB(temp1, p2, p1);
		SUB(temp2, p3, p1);
		CROSS(result, temp1, temp2);
	}

	floating _evaluateBezier(const floating& p0, const floating& p1, const floating& p2, const floating& p3, const floating& t, const floating& s);
	floating _evaluateBezier2(const floating& p0, const floating& p1, const floating& p2, const floating& t, const floating& s);

	bool bezierDecomposition(const floating& k0, const floating& k1, const floating& k2, const floating& k3,
		const floating& j0, const floating& j1, const floating& j2,
		floating& m0, floating& m1, floating& n0, floating& n1);
	bool getBezier3(floating* a0, floating* b0, floating* c0, floating* v0,
		floating* a1, floating* b1, floating* c1, floating* v1, floating& p0, floating& p1, floating& p2, floating& p3,
		floating* n0, floating* n1, floating* cross_for_CCD, floating* deltaN, floating* nX);

	int bezierClassification(const floating& k0, const floating& k1, const floating& k2, const floating& k3,
		floating& kk0, floating& kk1, floating& kk2);

	bool subdivide(floating& k0, floating& k1, floating& k2, floating& k3, floating* g);
	int sign(const double& a);
	int sameSign(const floating& a, const floating& b);
	int diffSign(const floating& a, const floating& b);
	void make_vector(const double* v, floating* out);
	floating det2x2(const floating& a, const floating& b, const floating& c, const floating& d);
	bool insideTest(BC& c, floating* g, int root_nums);
	void make_vector(double* v, floating* out);
	void make_vector(double* v, double* sigma, floating* out);

	int coplanarTest(BC& c);
	
};



