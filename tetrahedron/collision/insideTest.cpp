#include"insideTest.h"


bool InsideTest::insideTest(double* a0, double* b0, double* c0, double* v0,
	double* a1, double* b1, double* c1, double* v1, floating* n0, floating* n1, floating* cross_for_CCD,
	bool edge_edge)
{

	floating d_a0[3], d_b0[3], d_c0[3], d_v0[3], d_a1[3], d_b1[3], d_c1[3], d_v1[3];
	make_vector(a0, d_a0);
	make_vector(b0, d_b0);
	make_vector(c0, d_c0);
	make_vector(v0, d_v0);

	make_vector(a1, d_a1);
	make_vector(b1, d_b1);
	make_vector(c1, d_c1);
	make_vector(v1, d_v1);


	return insideRobust(d_a0, d_b0, d_c0, d_v0, d_a1, d_b1, d_c1, d_v1,
		n0, n1, cross_for_CCD, edge_edge);
}





bool InsideTest::insideRobust(floating* a0, floating* b0, floating* c0, floating* d0,
	floating* a1, floating* b1, floating* c1, floating* d1, floating* n0, floating* n1, floating* cross_for_CCD, bool ee_test)
{
	//floating j0, j1, j2;
	//floating lt0, lt1, kt0, kt1; // for signs of lt and kt
	bool bt0 = true, bt1 = true;
	floating k0, k1, k2, k3;
	floating deltaN[3], nX[3];
	if (!getBezier3(a0, b0, c0, d0, a1, b1, c1, d1, k0, k1, k2, k3, n0, n1, cross_for_CCD,deltaN,nX)) {
		return false;
	}
	floating kk0, kk1, kk2;
	int ct = bezierClassification(k0, k1, k2, k3, kk0, kk1, kk2);
	if (ct == 2) {
		floating g[15]; // for inside test
		for (int i = 0; i < 3; i++) {
			getBezier4(a0, b0, c0, d0, a1, b1, c1, d1, g[i * 5 + 0], g[i * 5 + 1], g[i * 5 + 2], g[i * 5 + 3], g[i * 5 + 4],
				n0, n1, deltaN, nX, i, ee_test);
			return subdivide(k0, k1, k2, k3, g, true, true);
		}
	}
	BC c(k0, k1, k2, k3, kk0, kk1, kk2, ct);

	//if (root_nums == 0)
	//	return false;//return -1;
	floating g[15]; // for inside test
	for (int i = 0; i < 3; i++) {
		getBezier4(a0, b0, c0, d0, a1, b1, c1, d1, g[i * 5 + 0], g[i * 5 + 1], g[i * 5 + 2], g[i * 5 + 3], g[i * 5 + 4], n0, n1, deltaN, nX, i, ee_test);
	}
	return insideTest(c, g);
}

bool InsideTest::insideTest(BC& c, floating* g)
{
	floating l0, l1, l2, l3, l4;
	floating j0, j1, j2;
	floating s0, s1; // for L(t)
	floating t0, t1; // for K(t)

	floating lt0, lt1, kt0, kt1; // for signs of lt and kt
	bool bt0[3], bt1[3];

	for (int i = 0; i < 3; i++) {
		bt0[i] = true;
		bt1[i] = true;

		l0 = g[i * 5 + 0]; l1 = g[i * 5 + 1]; l2 = g[i * 5 + 2];
		l3 = g[i * 5 + 3]; l4 = g[i * 5 + 4];

		getSimplifyed(c.k0, c.k1, c.k2, c.k3, l0, l1, l2, l3, l4, j0, j1, j2);
		if (!bezierDecomposition(c.k0, c.k1, c.k2, c.k3, j0, j1, j2, s0, s1, t0, t1)) {
			if ((j1 - j0) == 0) {
				if (j0.sign() < 0)
					return false;
				else
					continue;
			}
			getSigns(j0, j2, c, lt0, lt1);//getSigns(j0, j1*floating<T>(2.0, 0) - j0, c, lt0, lt1, root_nums);
			if (c.ct == 0 || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE)
			{
				if (lt0.sign() < 0)
					return false;
			}
			else if (c.ct == 1) { //(root_nums == 2){
				if (lt0.sign() < 0)
					bt0[i] = false;//bt0 = false;

				if (lt1.sign() < 0)
					bt1[i] = false;//bt1 = false;

				if (!bt0[i] && !bt1[i])//(!bt0 && !bt1)
					return false;
			}
			continue;
		}
		bool dg_l = getSigns(s0, s1, c, lt0, lt1);//for L(t)
		bool dg_k = getSigns(t0, t1, c, kt0, kt1);//for K(t)
		if (!dg_l || !dg_k) {
			continue;
		}
		if ((c.ct == 0 || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE)) {//(c.ct == 0 || (c.ct == 1 && diffSign(c.k0, c.k3)) == RETURN_TRUE) {//(c.ct == 0){////wzd modify 2015.1.28
			if (sameSign(lt0, kt0) == RETURN_TRUE)
				return false;
			continue;
		}
		// kill an possiblity 
		if (sameSign(lt0, kt0) == RETURN_TRUE)
			bt0[i] = false;//bt0 = false;

		// kill an possiblity 
		if (sameSign(lt1, kt1) == RETURN_TRUE)
			bt1[i] = false;//bt1 = false;

		//if no possiblity left, return false ...
		if (!bt0[i] && !bt1[i])//(!bt0 && !bt1)
			return false;
	}

	if (c.ct == 1) {
		bool bb0 = (bt0[0] && bt0[1] && bt0[2]);
		bool bb1 = (bt1[0] && bt1[1] && bt1[2]);
		if (!bb0 && !bb1)
			return false;
	}

	return true;
}

bool InsideTest::subdivide(floating& k0, floating& k1, floating& k2, floating& k3, floating* g, bool need_left, bool need_right)
{
	floating t = k0 - k1 * floating(2.0, 0) + k2;
	floating division = k0 - k1 * floating(3.0, 0) + k2 * floating(3.0, 0) - k3;

	floating l0, l1, l2, l3;
	if (division.sign() == 0)
		return true; 

	//k0, (1-t)*k0 + t*k1, (1-t)^2*k0 + 2*(1-t)*t*k1 + t^2*k2, (1-t)^3*k0 + 3*(1-t)^2*t*k1 + 3*(1-t)*t^2*k2 + t^3*k3
	if (need_left) {
		l0 = k0 * division * division * division;
		l1 = ((division - t) * k0 + t * k1) * division * division;
		l2 = ((division - t) * (division - t) * k0 + floating(2.0, 0) * t * (division - t) * k1 + t * t * k2) * division;
	}
		l3 = (division - t) * (division - t) * (division - t) * k0 + floating(3.0, 0) * t * (division - t) * (division - t) * k1
			+ floating(3.0, 0) * t * t * (division - t) * k2 + t * t * t * k3;
	
	//(1-t)^3*k0 + 3*(1-t)^2*t*k1 + 3*(1-t)*t^2*k2 + t^3*k3, (1-t)^2*k1 + 2*(1-t)*t*k2 + t^2*k3, (1-t)*k2 + t*k3, k3 
	floating r0, r1, r2, r3;
	if (need_right) {
		r0 = l3;
		r1 = ((division - t) * (division - t) * k1 + floating(2.0, 0) * t * (division - t) * k2 + t * t * k3) * division;
		r2 = ((division - t) * k2 + t * k3) * division * division;
		r3 = k3 * division * division * division;
	}
	
	if (division.sign() < 0) {
		if (need_left) {
			l0 = -l0;
			l1 = -l1;
			l2 = -l2;
			l3 = -l3;
		}
		if (need_right) {
			r0 = -r0;
			r1 = -r1;
			r2 = -r2;
			r3 = -r3;
		}
	}
	floating ll0, ll1, ll2;
	floating rr0, rr1, rr2;
	int ct_l, ct_r;

	bool left, right;
	if (need_left) {
		ct_l = bezierClassification(l0, l1, l2, l3, ll0, ll1, ll2);
		if (ct_l == 2) {
			if (diffSign(ll0, ll2)!= RETURN_FALSE)
				ct_l = 1; // no inflexion, 1 extreme
			else
				ct_l = 0; // no inflexion, no extreme
		}
		BC c_l(l0, l1, l2, l3, ll0, ll1, ll2, ct_l);
		floating gl[15];
		for (int i = 0; i < 3; i++) {
			gl[i * 5 + 0] = g[i * 5 + 0] * division * division * division * division;
			gl[i * 5 + 1] = ((division - t) * g[i * 5 + 0]
				+ t * g[i * 5 + 1]) * division * division * division;
			gl[i * 5 + 2] = ((division - t) * (division - t) * g[i * 5 + 0]
				+ floating(2.0, 0) * t * (division - t) * g[i * 5 + 1]
				+ t * t * g[i * 5 + 2]) * division * division;
			gl[i * 5 + 3] = ((division - t) * (division - t) * (division - t) * g[i * 5 + 0]
				+ floating(3.0, 0) * t * (division - t) * (division - t) * g[i * 5 + 1]
				+ floating(3.0, 0) * t * t * (division - t) * g[i * 5 + 2]
				+ t * t * t * g[i * 5 + 3]) * division;
			gl[i * 5 + 4] = ((division - t) * (division - t) * (division - t) * (division - t) * g[i * 5 + 0]
				+ floating(4.0, 0) * t * (division - t) * (division - t) * (division - t) * g[i * 5 + 1]
				+ floating(6.0, 0) * t * t * (division - t) * (division - t) * g[i * 5 + 2]
				+ floating(4.0, 0) * t * t * t * (division - t) * g[i * 5 + 3]
				+ t * t * t * t * g[i * 5 + 4]);
		}
		left = insideTest(c_l, gl);
	}
	if (need_right) {
		floating gr[15];
		ct_r = bezierClassification(r0, r1, r2, r3, rr0, rr1, rr2);
		if (ct_r == 2) {
			if (diffSign(rr0, rr2) != RETURN_FALSE)
				ct_r = 1; // no inflexion, 1 extreme
			else
				ct_r = 0; // no inflexion, no extreme
		}

		BC c_r(r0, r1, r2, r3, rr0, rr1, rr2, ct_r);
		floating gr[15];
		for (int i = 0; i < 3; i++) {
			gr[i * 5 + 4] = g[i * 5 + 4] * division * division * division * division;
			gr[i * 5 + 3] = ((division - t) * g[i * 5 + 3]
				+ t * g[i * 5 + 4]) * division * division * division;
			gr[i * 5 + 2] = ((division - t) * (division - t) * g[i * 5 + 2]
				+ floating(2.0, 0) * t * (division - t) * g[i * 5 + 3]
				+ t * t * g[i * 5 + 4]) * division * division;
			gr[i * 5 + 1] = ((division - t) * (division - t) * (division - t) * g[i * 5 + 1]
				+ floating(3.0, 0) * t * (division - t) * (division - t) * g[i * 5 + 2]
				+ floating(3.0, 0) * t * t * (division - t) * g[i * 5 + 3]
				+ t * t * t * g[i * 5 + 4]) * division;
			gr[i * 5 + 0] = ((division - t) * (division - t) * (division - t) * (division - t) * g[i * 5 + 0]
				+ floating(4.0, 0) * t * (division - t) * (division - t) * (division - t) * g[i * 5 + 1]
				+ floating(6.0, 0) * t * t * (division - t) * (division - t) * g[i * 5 + 2]
				+ floating(4.0, 0) * t * t * t * (division - t) * g[i * 5 + 3]
				+ t * t * t * t * g[i * 5 + 4]);
		}
		right = insideTest(c_r, gr);

		if (!left && !right)//
			return false;//return -1;//          
		return true;
	}
}


void InsideTest::getSimplifyed(floating& k0, floating& k1, floating& k2, floating& k3, floating& l0, floating& l1, floating& l2, 
	floating& l3, floating& l4, floating& j0, floating& j1, floating& j2)
{
	floating kk0 = k0 * floating(4.0, 0);
	floating kk1 = k0 + k1 * floating(3.0, 0);
	floating kk2 = (k1 + k2) * floating(2.0, 0);
	floating kk3 = k2 * floating(3.0, 0) + k3;
	floating kk4 = k3 * floating(4.0, 0);

	floating s0 = (l1 * kk0 - l0 * kk1) * floating(12.0, 0);
	floating s1 = (l2 * kk0 - l0 * kk2) * floating(6.0, 0);
	floating s2 = (l3 * kk0 - l0 * kk3) * floating(4.0, 0);
	floating s3 = (l4 * kk0 - l0 * kk4) * floating(3.0, 0);

	j0 = (s1 * k0 - s0 * k1) * floating(6.0, 0);
	j1 = (s2 * k0 - s0 * k2) * floating(3.0, 0);
	j2 = (s3 * k0 - s0 * k3) * floating(2.0, 0);
}

void InsideTest::getBezier4(floating* a0, floating* b0, floating* c0, floating* v0,
	floating* a1, floating* b1, floating* c1, floating* v1, 
	floating& l0, floating& l1, floating& l2, floating& l3, floating& l4,
	floating* n0, floating* n1, floating* deltaN, floating* nX,
	int which, bool ee_test)
{
	floating mX[3], deltaM[3], m0[3], m1[3];
	floating temp0[3], temp1[3], temp2[3];
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
		mX[i] = floating(0.5,0) * (m0[i] + m1[i] - deltaM[i]);
	}
	l0 = DOT(m0, n0) * floating(6.0, 0);
	l1 = (DOT(m0, nX) + DOT(mX, n0)) * floating(3.0, 0);
	l2=DOT(m0,n1)+ floating(4.0, 0) * DOT(mX, nX)+DOT(m1, n0);
	l3 = floating(3.0, 0) * (DOT(mX, n1) + DOT(m1, nX));
	l4 = floating(6.0, 0) * DOT(m1, n1);
	if (!ee_test && which != 0) {
		l0 = -l0, l1 = -l1, l2 = -l2, l3 = -l3, l4 = -l4;
	}

	if (which == 2 && ee_test) {
		l0 = -l0, l1 = -l1, l2 = -l2, l3 = -l3, l4 = -l4;
	}
}

bool InsideTest::getSigns(const floating& t0, const floating& t1, const BC& c, floating& lt0, floating& lt1)
{
	if (sameSign(t0, t1)==RETURN_TRUE) {
		lt0 = t0;
		lt1 = t0;
		return true;
	}

	if ((t0 - t1).sign() == 0) {
		return false;//here have some problem
	}

	if ((c.ct == 0) ||
		(c.ct == 1 && diffSign(c.k0, c.k3))) { //one root
		floating ft = _evaluateBezier(c.k0, c.k1, c.k2, c.k3, t0, -t1);

		if (ft.sign() == 0) {
			lt0 = floating(0, 1);//So lt0.sign()==0
			lt1 = lt0;//So lt1.sign()==0
			return true;
		}

		if ((t0 - t1).sign() < 0){//if (t0 < 0) {
			ft = -ft;
		}
		if (ft.sign()==c.sign_k0) {
			lt0 = t1;
			lt1 = t1;
		}
		else {
			lt0 = t0;
			lt1 = t0;
		}
		return true;
	}
	if (c.ct == 1) { //two roots
		floating ft= _evaluateBezier(c.k0, c.k1, c.k2, c.k3, t0, -t1);
		if (ft.sign() ==0) {
			if (c.sign_kk0==0) {
				lt0 = floating(0, 1);
				lt1 = lt0;
				return true;
			}
			floating fk = _evaluateBezier2(c.kk0, c.kk1, c.kk2, t0, -t1);
			if (fk.sign() == c.sign_kk0) {
				lt0 = floating(0, 1);
				lt1 = t1;
				return true;
			}
			if (fk.sign() == c.sign_kk2) {
				lt0 = t0;
				lt1 = floating(0, 1);
				return true;
			}
			if (fk.sign() == 0) {
				lt0 = t0;
				lt1 = t1;
				return true;
			}	
		}
		
		if ((t0 - t1).sign() < 0)//if (t0<0)//.is_certainly_negative())
			ft = -ft;

		if (ft.sign() != c.sign_k0) {//(diffSign(ft, c.k0)) {
			lt0 = t0;
			lt1 = t1;
			return true;
		}
		floating fk = _evaluateBezier2(c.kk0, c.kk1, c.kk2, t0, -t1);
		if (c.sign_kk0 == 0) {
			lt0 = floating(0, 1);//So lt0.sign()==0
			lt1 = lt0;//So lt1.sign()==0
			return true;
		}
		if (fk.sign() == c.sign_kk0) {
			lt0 = t1;
			lt1 = t1;
			return true;
		}
		if (fk.sign() == c.sign_kk2) {
			lt0 = t0;
			lt1 = t0;
			return true;
		}
		if (fk.sign() == 0) {
			lt0 = t0;
			lt1 = t1;
			//printf("Here is unreasonable!\n");
			return true;
		}
		return false;
	}
	std::cout << "c.ct should be 0/1, should not be here" << std::endl;
}

floating InsideTest::_evaluateBezier(const floating& p0, const floating& p1, const floating& p2, const floating& p3, const floating& t, const floating& s)
{
	floating s2 = s * s;
	floating s3 = s2 * s;
	floating t2 = t * t;
	floating t3 = t2 * t;
	return p0 * s3 + p1 * 3.0 * s2 * t + p2 * 3.0 * s * t2 + p3 * t3;
}

floating InsideTest::_evaluateBezier2(const floating& p0, const floating& p1, const floating& p2, const floating& t, const floating& s)
{
	floating s2 = s * s;
	floating t2 = t * t;
	return p0 * s2 + p1 * 2.0 * s * t + p2 * t2;
}


inline bool InsideTest::bezierDecomposition(const floating& k0, const floating& k1, const floating& k2, const floating& k3,
	const floating& j0, const floating& j1, const floating& j2,
	floating& m0, floating& m1, floating& n0, floating& n1)
{
	floating A = (j1 - j2) * floating(2.0, 0);
	floating B = j0 - j2;
	floating C = k2 * floating(3.0, 0) - k3 * floating(2.0, 0) - k0;
	floating D = k1 * floating(3.0, 0) - k0 * floating(2.0, 0) - k3;
	floating E = j2 - j0;
	floating F = (j1 - j0) * floating(2.0, 0);

	floating tt = det2x2(A, B, E, F);
	//Here has some problem
	if (tt==0 || (j0 + j2 - j1 * floating(2.0, 0) == 0)) {
		return false;
	}
	m0 = det2x2(A, B, C, D);
	m1 = det2x2(F, E, D, C);
	n0 = k0 * tt - m0 * j0;
	n1 = k3 * tt - m1 * j2;

	return true;
}


bool InsideTest::getBezier3(floating* a0, floating* b0, floating* c0, floating* v0,
	floating* a1, floating* b1, floating* c1, floating* v1, floating& p0, floating& p1, floating& p2, floating& p3,
	floating* n0, floating* n1, floating* cross_for_CCD,floating* deltaN, floating* nX)
{
	floating mX[3], deltaM[3], m0[3], m1[3];
	for (int i = 0; i < 3; ++i) {
		deltaN[i] = n0[i] + n1[i] - cross_for_CCD[i];
		nX[i] = floating(0.5,0) * (n0[i] + n1[i] - deltaN[i]);
	}

	floating pa0[3], pa1[3];
	SUB(pa0, v0, a0);
	SUB(pa1, v1, a1);
	
	floating A = DOT(n0, pa0);
	floating B = DOT(n1, pa1);
	floating C = DOT(nX, pa0);
	floating D = DOT(nX, pa1);
	floating E = DOT(n1, pa0);
	floating F = DOT(n0, pa1);

	p0 = A * floating(3.0,0);
	p1 = C * floating(2.0, 0) + F;
	p2 = D * floating(2.0, 0) + E;
	p3 = B * floating(3.0, 0);

	if (p0.sign() > 0 && p1.sign() > 0 && p2.sign() > 0 && p3.sign() > 0)
	{
		return false;
	}
	if (p0.sign() < 0 && p1.sign() < 0 && p2.sign() < 0 && p3.sign() < 0)
	{
		return false;
	}
	return  true;
}


inline int InsideTest::bezierClassification(const floating& k0, const floating& k1, const floating& k2, const floating& k3,
	floating& kk0, floating& kk1, floating& kk2)
{
	if (k0.sign() > 0 && k1.sign() < 0 && k2.sign() < 0 && k3.sign() < 0) {
		return 0;
	}
	if (k0.sign() < 0 && k1.sign() > 0 && k2.sign() > 0 && k3.sign() > 0) {
		return 0;
	}
	if (k3.sign() > 0 && k1.sign() < 0 && k2.sign() < 0 && k0.sign() < 0) {
		return 0;
	}
	if (k3.sign() < 0 && k1.sign() > 0 && k2.sign() > 0 && k0.sign() > 0) {
		return 0;
	}

	// f'' = 6*(k2-2*k1+k0)*B^0_1 + 6*(k3-2*k2+k1)*B^1_1
	floating a = k2 - k1 * floating(2.0, 0) + k0;
	floating b = k3 - k2 * floating(2.0, 0) + k1;

	if (diffSign(a, b) != RETURN_FALSE) {
		return 2; // 1 inflexion
	}

	// f' = 3*(k1-k0) B^2_0 + 3*(k2-k1)*B^2_1 + 3*(k3-k2)*B^2_2
	kk0 = k1 - k0;
	kk1 = k2 - k1;
	kk2 = k3 - k2;
	if (diffSign(kk0, kk2) != RETURN_FALSE)
		return 1; // no inflexion, 1 extreme
	else
		return 0;// no inflexion, no extreme
}

inline int  InsideTest::sign(const double& a)
{
	if (a > NEAR_ZERO2) {
		return 1;
	}
	else if (a < -NEAR_ZERO2) {
		return -1;
	}
	else {
		return 0;
	}
}

inline int InsideTest::sameSign(const floating& a, const floating& b)
{
	if (a.sign() == 0 || b.sign() == 0) {
		return RETURN_ZERO;
	}
	else if ((a.sign() > 0 && b.sign() > 0) || (a.sign() < 0 && b.sign() < 0)) {
		return RETURN_TRUE;
	}
	else {
		return RETURN_FALSE;
	}
}

inline int InsideTest::diffSign(const floating& a, const floating& b)
{
	if (a.sign() == 0 || b.sign() == 0) {
		return RETURN_ZERO;
	}
	else if ((a.sign() > 0 && b.sign() < 0) || (a.sign() < 0 && b.sign() > 0)) {
		return RETURN_TRUE;
	}
	else {
		return RETURN_FALSE;
	}
}

inline void InsideTest::make_vector(const double* v, floating* out)
{
	for (int i = 0; i < 3; i++) {
		out[i] = floating(v[i], 0);
	}
}

inline floating InsideTest::det2x2(const floating& a, const floating& b, const floating& c, const floating& d)
{
	return a * d - b * c;
}


inline void InsideTest::make_vector(double* v, floating* out)
{
	for (int i = 0; i < 3; i++) {
		out[i] = floating(v[i], 0);
	}
}


inline void InsideTest::make_vector(double* v, double* sigma, floating* out)
{
	for (int i = 0; i < 3; i++) {
		out[i] = floating(v[i], sigma[i]);
	}
}