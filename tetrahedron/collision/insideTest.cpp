#include"insideTest.h"


bool InsideTest::insideTest(double* a0, double* b0, double* c0, double* v0,
	double* a1, double* b1, double* c1, double* v1, double* n0, double* n1, double* cross_for_CCD, bool ee_test)
{

	double j0, j1, j2;

	double lt0, lt1, kt0, kt1; // for signs of lt and kt
	bool bt0 = true, bt1 = true;

	double k0, k1, k2, k3;
	if (!getBezier3(a0, b0, c0, v0, a1, b1, c1, v1, k0, k1, k2, k3,n0,n1,cross_for_CCD)) {
		return false;
	}
	double kk0, kk1, kk2;
	int ct = bezierClassification(k0, k1, k2, k3, kk0, kk1, kk2);

	double kkk0 = k2 - k1 * 2.0 + k0;
	double kkk1 = k3 - k2 * 2.0 + k1;

	if (ct == 2) {
		double g[15];//coefficient for inside test
		for (int i = 0; i < 3; i++) {
			getBezier4(a0, b0, c0, v0, a1, b1, c1, v1, g[i * 5 + 0], g[i * 5 + 1], g[i * 5 + 2], g[i * 5 + 3], g[i * 5 + 4], n0, n1, cross_for_CCD, i, ee_test);
		}


	}

	for (int i = 0; i < 3; ++i) {
		getBezier4(a0, b0, c0, v0, a1, b1, c1, v1, l0, l1, l2, l3, l4, n0, n1, cross_for_CCD, i, ee_test);
		getSimplifyed(c.k0, c.k1, c.k2, c.k3, l0, l1, l2, l3, l4, j0, j1, j2);

		if (abs(j0 + j2 - j1 * 2.0) < NEAR_ZERO2) {// degenerate j0, j1, j2
												   // the first derivative of P(t) is zero, c.ct=0/1
			getSigns(j0, j2, c, lt0, lt1);
			if (c.ct == 0) {
				if (lt0 < 0) {
					return false;
				}
			}
			else {
				if (lt0 < 0) {
					bt0 = false;
				}
				if (lt1 < 0) {
					bt1 = false;
				}
				if (!bt0 && !bt1)
					return false;
			}

			continue;
		}
		double s0, s1, t0, t1;
		if (!bezierDecomposition(c.k0, c.k1, c.k2, c.k3, j0, j1, j2, s0, s1, t0, t1)) { //why?
			getSigns(j0, j2, c, lt0, lt1);
			if (c.ct == 0) {
				if (lt0 < 0) {
					return false;
				}
			}
			else {
				if (lt0 < 0) {
					bt0 = false;
				}
				if (lt1 < 0) {
					bt1 = false;
				}
				if (!bt0 && !bt1)
					return false;
			}

			continue;
		}
		getSigns(t0, t1, c, lt0, lt1);
		getSigns(s0, s1, c, kt0, kt1);

		if (c.ct == 0) {
			if (SAME_SIGN(lt0, kt0))
				return false;

			continue;
		}
		if (SAME_SIGN(lt0, kt0))
			bt0 = false;
		if (SAME_SIGN(lt1, kt1))
			bt1 = false;
		if (!bt0 && !bt1)
			return false;
	}
	return true;

}




void InsideTest::getSimplifyed(double& k0, double& k1, double& k2, double& k3, double& l0, double& l1, double& l2, double& l3, double& l4,
	double& j0, double& j1, double& j2)
{
	double kk0 = k0 * 4.0;
	double kk1 = k0 + k1 * 3.0;
	double kk2 = (k1 + k2) * 2.0;
	double kk3 = k2 * 3.0 + k3;
	double kk4 = k3 * 4.0;

	double s0= (l1 * kk0 - l0 * kk1) * 12.0;
	double s1= (l2 * kk0 - l0 * kk2) * 6.0;
	double s2 = (l3 * kk0 - l0 * kk3) * 4.0;
	double s3 = (l4 * kk0 - l0 * kk4) * 3.0;

	j0 = (s1 * k0 - s0 * k1) * 6.0;
	j1 = (s2 * k0 - s0 * k2) * 3.0;
	j2 = (s3 * k0 - s0 * k3) * 2.0;
}

void InsideTest::getBezier4(double* a0, double* b0, double* c0, double* v0,
	double* a1, double* b1, double* c1, double* v1, double& l0, double& l1, double& l2, double& l3, double& l4,
	double* n0, double* n1, double* cross_for_CCD, int which, bool ee_test)
{
	double deltaN[3];	
	double nX[3], mX[3], deltaM[3], m0[3], m1[3];
	for (int i = 0; i < 3; ++i) {
		deltaN[i] = n0[i] + n1[i] - cross_for_CCD[i];
		nX[i] = 0.5 * (n0[i] + n1[i] - deltaN[i]);
	}
	double temp0[3], temp1[3], temp2[3];
	SUB(temp0, v1, v0);
	if (which == 0) {
		norm(m0, v0, b0, c0);
		norm(m1, v1, b1, c1);	
		SUB(temp1, b1, b0);
		SUB(temp2, c1, c0);
		norm(deltaM, temp0, temp1, temp2);
	}
	else if (which == 1) {
		norm(m0, v0, c0, a0);
		norm(m1, v1, c1, a1);
		SUB(temp2, a1, a0);
		SUB(temp1, c1, c0);
		norm(deltaM, temp0, temp1, temp2);
	}
	else{
		norm(m0, v0, a0, b0);
		norm(m1, v1, a1, b1);
		SUB(temp1, a1, a0);
		SUB(temp2, b1, b0);
		norm(deltaM, temp0, temp1, temp2);
	}
	for (int i = 0; i < 3; ++i) {
		mX[i] = 0.5 * (m0[i] + m1[i] - deltaM[i]);
	}
	l0 = DOT(m0, n0) * 6.0;
	l1 = (DOT(m0, nX) + DOT(mX, n0)) * 3.0;
	l2=DOT(m0,n1)+4.0* DOT(mX, nX)+DOT(m1, n0);
	l3 = 3.0 * (DOT(mX, n1) + DOT(m1, nX));
	l4 = 6.0 * DOT(m1, n1);
	if (!ee_test && which != 0) {
		l0 = -l0, l1 = -l1, l2 = -l2, l3 = -l3, l4 = -l4;
	}

	if (which == 2 && ee_test) {
		l0 = -l0, l1 = -l1, l2 = -l2, l3 = -l3, l4 = -l4;
	}
}

void InsideTest::getSigns(const double& t0, const double& t1, const BC& c, double& lt0, double& lt1)
{
	if (SAME_SIGN(t0, t1)) {
		lt0 = t0;
		lt1 = t0;
		return;
	}
	if ((c.ct == 0) ||
		(c.ct == 1 && DIFF_SIGN(c.k0, c.k3))) {
		double ft = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, t0, -t1);
		if (t0 < 0) {
			ft = -ft;
		}
		if (SAME_SIGN(ft, c.k0)) {
			lt0 = t1;
			lt1 = t1;
		}
		else {
			lt0 = t0;
			lt1 = t0;
		}
		return;
	}
	if (c.ct == 1) {
		double ft= _evaluateBezier(c.k0, c.k1, c.k2, c.k3, t0, -t1);
		if (t0 < 0) {
			ft = -ft;
		}
		if (DIFF_SIGN(ft, c.k0)) {
			lt0 = t0;
			lt1 = t1;
			return;
		}
		double fk = _evaluateBezier2(c.kk0, c.kk1, c.kk2, t0, -t1);
		if (SAME_SIGN(fk, c.kk0))
			lt0 = lt1 = t1;
		else
			lt0 = lt1 = t0;
		return;
	}
	std::cout << "c.ct should be 0/1, should not be here" << std::endl;
}

double InsideTest::_evaluateBezier(const double& p0, const double& p1, const double& p2, const double& p3, const double& t, const double& s)
{
	double s2 = s * s;
	double s3 = s2 * s;
	double t2 = t * t;
	double t3= t2 * t;
	return p0 * s3 + p1 * 3.0 * s2 * t + p2 * 3.0 * s * t2 + p3 * t3;
}

double InsideTest::_evaluateBezier2(const double& p0, const double& p1, const double& p2, const double& t, const double& s)
{
	double s2 = s * s;
	double t2 = t * t;
	return p0 * s2 + p1 * 2.0 * s * t + p2 * t2;
}


inline bool InsideTest::bezierDecomposition(const double& k0, const double& k1, const double& k2, const double& k3,
	const double& j0, const double& j1, const double& j2,
	double& m0, double& m1, double& n0, double& n1)
{
	double A = (j1 - j2) * 2.0;
	double B = j0 - j2;
	double C = k2 * 3.0 - k3 * 2.0 - k0;
	double D = k1 * 3.0 - k0 * 2.0 - k3;
	double E = j2 - j0;
	double F = (j1 - j0) * 2.0;

	double tt = DET2X2(A, B, E, F);
	if (abs(tt) < NEAR_ZERO2) {
		return false;
	}
	m0 = DET2X2(A, B, C, D);
	m1 = DET2X2(F, E, D, C);
	n0 = k0 * tt - m0 * j0;
	n1 = k3 * tt - m1 * j2;

	return true;
}


bool InsideTest::getBezier3(double* a0, double* b0, double* c0, double* v0,
	double* a1, double* b1, double* c1, double* v1, double& p0, double& p1, double& p2, double& p3,
	double* n0, double* n1, double* cross_for_CCD)
{
	double deltaN[3];
	double nX[3], mX[3], deltaM[3], m0[3], m1[3];
	for (int i = 0; i < 3; ++i) {
		deltaN[i] = n0[i] + n1[i] - cross_for_CCD[i];
		nX[i] = 0.5 * (n0[i] + n1[i] - deltaN[i]);
	}

	double pa0[3], pa1[3];
	SUB(pa0, v0, a0);
	SUB(pa1, v1, a1);
	
	double A = DOT(n0, pa0);
	double B = DOT(n1, pa1);
	double C = DOT(nX, pa0);
	double D = DOT(nX, pa1);
	double E = DOT(n1, pa0);
	double F = DOT(n0, pa1);

	p0 = A * 3.0;
	p1 = C * 2.0 + F;
	p2 = D * 2.0 + E;
	p3 = B * 3.0;

	if (p0 > 0 && p1 > 0 && p2 > 0 && p3 > 0)
	{
		return false;
	}
	if (p0 < 0 && p1 < 0 && p2 < 0 && p3 < 0)
	{
		return false;
	}
	return  true;
}


inline int InsideTest::bezierClassification(const double& k0, const double& k1, const double& k2, const double& k3,
	double& kk0, double& kk1, double& kk2)
{
	if (k0 < 0 && k1 > 0 && k2 > 0 && k3 > 0) {
		return 0;
	}
	if (k0 > 0 && k1 < 0 && k2 < 0 && k3 < 0) {
		return 0;
	}
	if (k3 > 0 && k1 < 0 && k2 < 0 && k0 < 0) {
		return 0;
	}
	if (k3 < 0 && k1 > 0 && k2 > 0 && k0 > 0) {
		return 0;
	}

	// f'' = 6*(k2-2*k1+k0)*B^0_1 + 6*(k3-2*k2+k1)*B^1_1
	double a = k2 - k1 * 2.0 + k0;
	double b = k3 - k2 * 2.0 + k1;
	if (DIFF_SIGN(a, b)) {
		return 2;
	}
	// f' = 3*(k1-k0) B^2_0 + 3*(k2-k1)*B^2_1 + 3*(k3-k2)*B^2_2
	kk0 = k1 - k0;
	kk1 = k2 - k1;
	kk2 = k3 - k2;
	if (DIFF_SIGN(kk0, kk2))
		return 1; // no inflexion, 1 extreme
	else
		return 0;// no inflexion, no extreme
}