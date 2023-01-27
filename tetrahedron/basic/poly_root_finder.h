#pragma once
#include <math.h>  
namespace poly_root {
	/*  Calculate the zeros of the quadratic a*z^2 + b1*z + c.
	*  The quadratic formula, modified to avoid overflow, is used
	*  to find the larger zero if the zeros are real and both
	*  are complex. The smaller real zero is found directly from
	*  the product of the zeros c/a.
	*/
	template<class T>
	inline void quad(T a, T b1, T c, T* sr, T* si,
		T* lr, T* li)
	{
		T b, d, e;

		if (a == 0.0) {         /* less than two roots */
			if (b1 != 0.0)
				*sr = -c / b1;
			else
				*sr = 0.0;
			*lr = 0.0;
			*si = 0.0;
			*li = 0.0;
			return;
		}
		if (c == 0.0) {         /* one real root, one zero root */
			*sr = 0.0;
			*lr = -b1 / a;
			*si = 0.0;
			*li = 0.0;
			return;
		}
		/* Compute discriminant avoiding overflow. */
		b = b1 / 2.0;
		if (fabs(b) < fabs(c)) {
			if (c < 0.0)
				e = -a;
			else
				e = a;
			e = b * (b / fabs(c)) - e; // if(c>0) (b^2-4ac)/c  if(c<0) (4ac-b^2)/c
			d = sqrt(fabs(e)) * sqrt(fabs(c));//sqrt(b^2-4ac)
		}
		else {//avoid overflow.
			e = 1.0 - (a / b) * (c / b);// (b^2-4ac)/c
			d = sqrt(fabs(e)) * fabs(b);//sqrt(b^2-4ac)
		}
		if (e < 0.0) {      /* complex conjugate zeros */
			*sr = -b / a;
			*lr = *sr;
			*si = fabs(d / a);
			*li = -(*si);
		}
		else {
			if (b >= 0.0)   /* real zeros. */
				d = -d;
			*lr = (-b + d) / a;//if (b>0) lr=(-b-sqrt(b^2-4ac))/a, if(b<0) lr=(-b+sqrt(b^2-4ac))/a  this is also to avoid subtract between two positive number
			*sr = 0.0;
			if (*lr != 0.0)
				*sr = (c / *lr) / a;
			*si = 0.0;
			*li = 0.0;
		}
	}

	/*  Divides p by the quadratic 1,u,v placing the quotient
	*  in q and the remainder in a,b.
	*/
	template<class T>
	inline void quadsd(int nn, T* u, T* v, T* p, T* q,
		T* a, T* b)
	{
		T c;
		int i;
		*b = p[0];
		q[0] = *b;
		*a = p[1] - (*b) * (*u);
		q[1] = *a;
		for (i = 2; i <= nn; i++) {
			c = p[i] - (*a) * (*u) - (*b) * (*v);
			q[i] = c;
			*b = *a;
			*a = c;
		}
	}


	/*  Compute new estimates of the quadratic coefficients
	*  using the scalars computed in calcsc.
	*/
	template<class T>
	inline void newest(int type, T* uu, T* vv, T* coe, T* a, T* p, T* k, int* n)
	{
		T a4, a5, b1, b2, c1, c2, c3, c4, temp;

		/* Use formulas appropriate to setting of type. */
		if (type == 3) {
			/*  If type=3 the quadratic is zeroed. */
			*uu = 0.0;
			*vv = 0.0;
			return;
		}
		if (type == 2) {
			a4 = (coe[0] + coe[6]) * coe[5] + coe[7];
			a5 = (coe[5] + a[7]) * coe[2] + a[8] * coe[3];
		}
		else {
			a4 = coe[0] + a[7] * coe[1] + coe[7] * coe[5];
			a5 = coe[2] + (a[7] + a[8] * coe[5]) * coe[3];
		}
		/*  Evaluate new quadratic coefficients. */
		b1 = -k[*n - 1] / p[*n];
		b2 = -(k[*n - 2] + b1 * p[*n - 1]) / p[*n];
		c1 = a[8] * b2 * a[0];
		c2 = b1 * a[6];
		c3 = b1 * b1 * a[2];
		c4 = c1 - c2 - c3;
		temp = a5 + b1 * a4 - c4;
		if (temp == 0.0) {
			*uu = 0.0;
			*vv = 0.0;
			return;
		}
		*uu = a[7] - (a[7] * (c3 + c2) + a[8] * (b1 * a[0] + b2 * a[6])) / temp;
		*vv = a[8] * (1.0 + c4 / temp);
		return;
	}


	/*  Computes the next k polynomials using scalars
	*  computed in calcsc.
	*/
	template<class T>
	inline void nextk(int* type, T* coe, T* a, T* machine_const, T* qp,
		T* k, T* qk, int* n)
	{
		T temp;
		int i;

		if (*type == 3) {
			/*  Use unscaled form of the recurrence if type is 3. */
			k[0] = 0.0;
			k[1] = 0.0;
			for (i = 2; i < *n; i++) {
				k[i] = qk[i - 2];
			}
			return;
		}
		temp = coe[0];
		if (*type == 1) temp = coe[1];
		if (fabs(a[0]) <= fabs(temp) * machine_const[0] * 10.0) {
			/*  If a1 is nearly zero then use a special form of the
				*  recurrence.
				*/
			k[0] = 0.0;
			k[1] = -a[6] * qp[0];
			for (i = 2; i < *n; i++) {
				k[i] = a[2] * qk[i - 2] - a[6] * qp[i - 1];
			}
			return;
		}
		/*  Use scaled form of the recurrence. */
		a[6] /= a[0];
		a[2] /= a[0];
		k[0] = qp[0];
		k[1] = qp[1] - a[6] * qp[0];
		for (i = 2; i < *n; i++) {
			k[i] = a[2] * qk[i - 2] - a[6] * qp[i - 1] + qp[i];
		}
	}
	/*  This routine calculates scalar quantities used to
	*  compute the next k polynomial and new estimates of
	*  the quadratic coefficients.
	*  type - integer variable set here indicating how the
	*  calculations are normalized to avoid overflow.
	*/
	template<class T>
	inline void calcsc(int* type, T* coe, T* a, T* machine_const, T* k, T* qk, int* n)
	{
		/*  Synthetic division of k by the quadratic 1,u,v */
		quadsd(*n - 1, a + 7, a + 8, k, qk, coe + 2, coe + 3);
		if (fabs(coe[2]) > fabs(k[*n - 1] * 100.0 * machine_const[0])) goto _10;
		if (fabs(coe[3]) > fabs(k[*n - 2] * 100.0 * machine_const[0])) goto _10;
		*type = 3;
		/*  Type=3 indicates the quadratic is almost a factor of k. */
		return;
	_10:
		if (fabs(coe[3]) < fabs(coe[2])) {
			*type = 1;
			/*  Type=1 indicates that all formulas are divided by c. */
			coe[4] = coe[0] / coe[2];
			coe[5] = coe[3] / coe[2];
			coe[6] = a[7] * coe[4];
			coe[7] = a[8] * coe[1];
			a[2] = coe[0] * coe[4] + (coe[7] / coe[2] + coe[6]) * coe[1];
			a[0] = coe[1] - coe[0] * (coe[3] / coe[2]);
			a[6] = coe[0] + coe[6] * coe[3] + coe[7] * coe[5];
			return;
		}
		*type = 2;
		/*  Type=2 indicates that all formulas are divided by d. */
		coe[4] = coe[0] / coe[3];
		coe[5] = coe[2] / coe[3];
		coe[6] = a[7] * coe[1];
		coe[7] = a[8] * coe[1];
		a[2] = (coe[0] + coe[6]) * coe[4] + coe[7] * (coe[1] / coe[3]);
		a[0] = coe[1] * coe[5] - coe[0];
		a[6] = (coe[5] + a[7]) * coe[0] + coe[7];
	}


	/*  Variable-shift H polynomial iteration for a real zero.
	*  sss - starting iterate
	*  nz  - number of zeros found
	*  iflag - flag to indicate a pair of zeros near real axis.
	*/
	template<class T>
	inline void realit(T* sss, int* nz, int* iflag, T* root_, T* machine_const,
		T* p, T* qp, T* k, T* qk, int* n)
	{
		T pv, kv, t, s;
		T ms, mp, omp, ee;
		int i, j;

		omp = 0;
		t = 0;

		*nz = 0;
		s = *sss;
		*iflag = 0;
		j = 0;
		/*  Main loop */
		while (1) {
			pv = p[0];
			/*  Evaluate p at s. */
			qp[0] = pv;
			for (i = 1; i <= *n; i++) {
				pv = pv * s + p[i];
				qp[i] = pv;
			}
			mp = fabs(pv);
			/*  Compute a rigorous bound on the error in evaluating p. */
			ms = fabs(s);
			ee = (machine_const[2] / (machine_const[1] + machine_const[2])) * fabs(qp[0]);
			for (i = 1; i <= *n; i++) {
				ee = ee * ms + fabs(qp[i]);
			}
			/*  Iteration has converged sufficiently if the polynomial
				*  value is less than 20 times this bound.
				*/
			if (mp <= 20.0 * ((machine_const[1] + machine_const[2]) * ee - machine_const[2] * mp)) {
				*nz = 1;
				root_[2] = s;
				root_[3] = 0.0;
				return;
			}
			j++;
			/*  Stop iteration after 10 steps. */
			if (j > 10) return;
			if (j < 2) goto _50;
			if (fabs(t) > 0.001 * fabs(s - t) || mp < omp) goto _50;
			/*  A cluster of zeros near the real axis has been
				*  encountered. Return with iflag set to initiate a
				*  quadratic iteration.
				*/
			*iflag = 1;
			*sss = s;
			return;
			/*  Return if the polynomial value has increased significantly. */
		_50:
			omp = mp;
			/*  Compute t, the next polynomial, and the new iterate. */
			kv = k[0];
			qk[0] = kv;
			for (i = 1; i < *n; i++) {
				kv = kv * s + k[i];
				qk[i] = kv;
			}
			if (fabs(kv) <= fabs(k[*n - 1]) * 10.0 * machine_const[0]) {
				/*  Use unscaled form. */
				k[0] = 0.0;
				for (i = 1; i < *n; i++) {
					k[i] = qk[i - 1];
				}
			}
			else {
				/*  Use the scaled form of the recurrence if the value
					*  of k at s is nonzero.
					*/
				t = -pv / kv;
				k[0] = qp[0];
				for (i = 1; i < *n; i++) {
					k[i] = t * qk[i - 1] + qp[i];
				}
			}
			kv = k[0];
			for (i = 1; i < *n; i++) {
				kv = kv * s + k[i];
			}
			t = 0.0;
			if (fabs(kv) > fabs(k[*n - 1] * 10.0 * machine_const[0])) t = -pv / kv;
			s += t;
		}
	}
	/*  Variable-shift k-polynomial iteration for a
	*  quadratic factor converges only if the zeros are
	*  equimodular or nearly so.
	*  uu, vv - coefficients of starting quadratic.
	*  nz - number of zeros found.
	*/
	template<class T>
	inline void quadit(T* uu, T* vv, int* nz, T* coe, T* root_,
		T* a, T* machine_const, T* p, T* qp, T* k, T* qk, int* n)
	{
		T ui, vi;
		T mp, omp, ee, relstp, t, zm;
		int type, i, j, tried;

		relstp = 0;
		omp = 0;

		*nz = 0;
		tried = 0;
		a[7] = *uu;
		a[8] = *vv;
		j = 0;
		/*  Main loop. */
	_10:
		quad(1.0, a[7], a[8], root_ + 2, root_ + 3, root_ + 4, root_ + 5);
		/*  Return if roots of the quadratic are real and not
			*  close to multiple or nearly equal and of opposite
			*  sign.
			*/
		if (fabs(fabs(root_[2]) - fabs(root_[4])) > 0.01 * fabs(root_[4])) return;
		/*  Evaluate polynomial by quadratic synthetic division. */
		quadsd(*n, a + 7, a + 8, p, qp, coe, coe + 1);
		mp = fabs(coe[0] - root_[2] * coe[1]) + fabs(root_[3] * coe[1]);
		/*  Compute a rigorous bound on the rounding error in
			*  evaluating p.
			*/
		zm = sqrt(fabs(a[8]));
		ee = 2.0 * fabs(qp[0]);
		t = -root_[2] * coe[1];
		for (i = 1; i < *n; i++) {
			ee = ee * zm + fabs(qp[i]);
		}
		ee = ee * zm + fabs(coe[0] + t);
		ee *= (5.0 * machine_const[2] + 4.0 * machine_const[1]);
		ee = ee - (5.0 * machine_const[2] + 2.0 * machine_const[1]) * (fabs(coe[0] + t) + fabs(coe[1]) * zm) + 2.0 * machine_const[1] * fabs(t);
		/*  Iteration has converged sufficiently if the
			*  polynomial value is less than 20 times this bound.
			*/
		if (mp <= 20.0 * ee) {
			*nz = 2;
			return;
		}
		j++;
		/*  Stop iteration after 20 steps. */
		if (j > 20) return;
		if (j < 2) goto _50;
		if (relstp > 0.01 || mp < omp || tried) goto _50;
		/*  A cluster appears to be stalling the convergence.
			*  Five fixed shift steps are taken with a u,v close
			*  to the cluster.
			*/
		if (relstp < machine_const[0]) relstp = machine_const[0];
		relstp = sqrt(relstp);
		a[7] = a[7] - a[7] * relstp;
		a[8] = a[8] + a[8] * relstp;
		quadsd(*n, a + 7, a + 8, p, qp, coe, coe + 1);
		for (i = 0; i < 5; i++) {
			calcsc(&type, coe, a, machine_const, k, qk, n);
			nextk(&type, coe, a, machine_const, qp, k, qk, n);
		}
		tried = 1;
		j = 0;
	_50:
		omp = mp;
		/*  Calculate next k polynomial and new u and v. */
		calcsc(&type, coe, a, machine_const, k, qk, n);
		nextk(&type, coe, a, machine_const, qp, k, qk, n);
		calcsc(&type, coe, a, machine_const, k, qk, n);
		newest(type, &ui, &vi, coe, a, p, k, n);
		/*  If vi is zero the iteration is not converging. */
		if (vi == 0.0) return;
		relstp = fabs((vi - a[8]) / vi);
		a[7] = ui;
		a[8] = vi;
		goto _10;
	}



	/*  Computes up to L2 fixed shift k-polynomials,
	 *  testing for convergence in the linear or quadratic
	 *  case. Initiates one of the variable shift
	 *  iterations and returns with the number of zeros
	 *  found.
	 */
	template<class T>
	inline void fxshfr(int l2, int* nz, T* coe, T* root_, T* a, T* machine_const,
		T* p, T* qp, T* k, T* qk, T* svk, int* n)
	{
		T svu, svv, ui, vi, s;
		T betas, betav, oss, ovv, ss, vv, ts, tv;
		T ots = 0, otv = 0, tvv, tss;
		int type, i, j, iflag, vpass, spass, vtry, stry;

		*nz = 0;
		betav = 0.25;
		betas = 0.25;
		oss = root_[0];
		ovv = a[8];
		/*  Evaluate polynomial by synthetic division. */
		quadsd(*n, a + 7, a + 8, p, qp, coe, coe + 1);
		calcsc(&type, coe, a, machine_const, k, qk, n);
		for (j = 0; j < l2; j++) {
			/*  Calculate next k polynomial and estimate v. */
			nextk(&type, coe, a, machine_const, qp, k, qk, n);
			calcsc(&type, coe, a, machine_const, k, qk, n);
			newest(type, &ui, &vi, coe, a, p, k, n);
			vv = vi;
			/*  Estimate s. */
			ss = 0.0;
			if (k[*n - 1] != 0.0) ss = -p[*n] / k[*n - 1];
			tv = 1.0;
			ts = 1.0;
			if (j == 0 || type == 3) goto _70;
			/*  Compute relative measures of convergence of s and v sequences. */
			if (vv != 0.0) tv = fabs((vv - ovv) / vv);
			if (ss != 0.0) ts = fabs((ss - oss) / ss);
			/*  If decreasing, multiply two most recent convergence measures. */
			tvv = 1.0;
			if (tv < otv) tvv = tv * otv;
			tss = 1.0;
			if (ts < ots) tss = ts * ots;
			/*  Compare with convergence criteria. */
			vpass = (tvv < betav);
			spass = (tss < betas);
			if (!(spass || vpass)) goto _70;
			/*  At least one sequence has passed the convergence test.
				*  Store variables before iterating.
				*/
			svu = a[7];
			svv = a[8];
			for (i = 0; i < *n; i++) {
				svk[i] = k[i];
			}
			s = ss;
			/*  Choose iteration according to the fastest converging
				*  sequence.
				*/
			vtry = 0;
			stry = 0;
			if ((spass && (!vpass)) || tss < tvv) goto _40;
		_20:
			quadit(&ui, &vi, nz, coe, root_, a, machine_const, p, qp, k, qk, n);
			if (*nz > 0) return;
			/*  Quadratic iteration has failed. Flag that it has
				*  been tried and decrease the convergence criterion.
				*/
			vtry = 1;
			betav *= 0.25;
			/*  Try linear iteration if it has not been tried and
				*  the S sequence is converging.
				*/
			if (stry || !spass) goto _50;
			for (i = 0; i < *n; i++) {
				k[i] = svk[i];
			}
		_40:
			realit(&s, nz, &iflag, root_, machine_const, p, qp, k, qk, n);
			if (*nz > 0) return;
			/*  Linear iteration has failed. Flag that it has been
				*  tried and decrease the convergence criterion.
				*/
			stry = 1;
			betas *= 0.25;
			if (iflag == 0) goto _50;
			/*  If linear iteration signals an almost T real
				*  zero attempt quadratic iteration.
				*/
			ui = -(s + s);
			vi = s * s;
			goto _20;
			/*  Restore variables. */
		_50:
			a[7] = svu;
			a[8] = svv;
			for (i = 0; i < *n; i++) {
				k[i] = svk[i];
			}
			/*  Try quadratic iteration if it has not been tried
				*  and the V sequence is convergin.
				*/
			if (vpass && !vtry) goto _20;
			/*  Recompute QP and scalar values to continue the
				*  second stage.
				*/
			quadsd(*n, a + 7, a + 8, p, qp, coe, coe + 1);
			calcsc(&type, coe, a, machine_const, k, qk, n);
		_70:
			ovv = vv;
			oss = ss;
			otv = tv;
			ots = ts;
		}
	}

	template<class T>
	inline int rpoly(T* op, int degree, T* zeror, T* zeroi)
	{
		T t, aa, bb, cc, factor, rot;
		T lo, max, min, xx, yy, cosr, sinr, xxx, x, sc, bnd;
		T xm, ff, df, dx, infin, smalno, base;
		int cnt, nz, i, j, jj, l, nm1, zerok;

		T coe[8];//a, b, c, d, e, f, g, h;	
		T root_[6];//sr, si, szr, szi, lzr, lzi;
		T a[9];// a1, a2, a3, a4, a5, a6, a7,u,v;
		T machine_const[3];//eta, are, mre;
		T p[7], qp[7], k[7], qk[7], svk[7];
		int n;
		/*  The following statements set machine constants. */
		base = 2.0;
		machine_const[0] = 2.22e-16;
		infin = 3.4e38;
		smalno = 1.2e-38;
		machine_const[1] = machine_const[0];
		machine_const[2] = machine_const[0];

		lo = smalno / machine_const[0];



		/*  Initialization of constants for shift rotation. */
		xx = sqrt(0.5);
		yy = -xx;
		rot = 94.0;
		rot *= 0.017453293;
		cosr = cos(rot);
		sinr = sin(rot);
		n = degree;
		/*  Algorithm fails of the leading coefficient is zero. */
		if (op[0] == 0.0) return -1;
		/*  Remove the zeros at the origin, if any. */
		while (op[n] == 0.0) {
			j = degree - n;
			zeror[j] = 0.0;
			zeroi[j] = 0.0;
			n--;
		}
		if (n < 1)
			return -1;

		/*
		 *  Allocate memory here
		 */
		T temp[7];
		//		temp = new T [degree+1];
		T pt[7];
		//		pt = new T [degree+1];
		//		svk = new T [degree+1];

	/*  Make a copy of the coefficients. */
		for (i = 0; i <= n; i++)
			p[i] = op[i];
		/*  Start the algorithm for one zero. */
	_40:
		if (n == 1) {
			zeror[degree - 1] = -p[1] / p[0];
			zeroi[degree - 1] = 0.0;
			n -= 1;
			goto _99;
		}
		/*  Calculate the final zero or pair of zeros. */
		if (n == 2) {
			quad(p[0], p[1], p[2], &zeror[degree - 2], &zeroi[degree - 2],
				&zeror[degree - 1], &zeroi[degree - 1]);
			n -= 2;
			goto _99;
		}
		/*  Find largest and smallest moduli of coefficients. */
		max = 0.0;
		min = infin;
		for (i = 0; i <= n; i++) {
			x = fabs(p[i]);
			if (x > max) max = x;
			if (x != 0.0 && x < min) min = x;
		}
		/*  Scale if there are large or very small coefficients.
		 *  Computes a scale factor to multiply the coefficients of the
		 *  polynomial. The scaling si done to avoid overflow and to
		 *  avoid undetected underflow interfering with the convergence
		 *  criterion. The factor is a power of the base.
		 */
		sc = lo / min;
		if (sc > 1.0 && infin / sc < max) goto _110;//when max>min*infi/low,max min differ by very many orders of magnitude
		if (sc <= 1.0) {
			if (max < 10.0) goto _110;
			if (sc == 0.0)
				sc = smalno;//????????????????
		}
		l = (int)(log(sc) / log(base) + 0.5);
		factor = pow(base * 1.0, l);//nearest 2^n
		if (factor != 1.0) {//sc>=1.0  this is when min is extremely small, smaller than lo
			for (i = 0; i <= n; i++)
				p[i] = factor * p[i];    // Scale polynomial. at least to the magnitude of lo
		}
	_110:
		/*  Compute lower bound on moduli of roots. */
		for (i = 0; i <= n; i++) {
			pt[i] = (fabs(p[i]));
		}
		pt[n] = -pt[n];
		/*  Compute upper estimate of bound. */
		x = exp((log(-pt[n]) - log(pt[0])) / (T)n);  //  (-pt[n]/pt[0])^{1/n} ????  the result of pt[n]x^n+pt[0]=0?
		//this is the upper bound of the equation pt[n]x^n+pt[n-1]x^(n-1)+...+pt[0]. 
		//if newton step at the origin is better, use it. better means smaller
		/*  If Newton step at the origin is better, use it. */
		if (pt[n - 1] != 0.0) {
			xm = -pt[n] / pt[n - 1];
			if (xm < x)  x = xm;
		}
		/*  Chop the interval (0,x) until ff <= 0 */
		while (1) {
			xm = x * 0.1;
			ff = pt[0];
			for (i = 1; i <= n; i++)
				ff = ff * xm + pt[i];
			if (ff <= 0.0) break;
			x = xm;
		}
		dx = x;
		/*  Do Newton interation until x converges to two
		 *  decimal places.
		 */
		while (fabs(dx / x) > 0.005) {
			ff = pt[0];
			df = ff;
			for (i = 1; i < n; i++) {
				ff = ff * x + pt[i];
				df = df * x + ff;//df = ((x*pt[0]+x*pt[0]+pt[1])*x+(pt[0]*x + pt[1])*x+pt[2])*x+...  this is to compute derivative
			}
			ff = ff * x + pt[n];
			dx = ff / df;
			x -= dx;
		}
		bnd = x;
		/*  Compute the derivative as the initial k polynomial
		 *  and do 5 steps with no shift.
		 */
		nm1 = n - 1;
		for (i = 1; i < n; i++)
			k[i] = (T)(n - i) * p[i] / (T)n;//here is to compute the derivative. then divide all coefficients by n
		k[0] = p[0];
		aa = p[n];
		bb = p[n - 1];
		zerok = (k[n - 1] == 0);
		for (jj = 0; jj < 5; jj++) {
			cc = k[n - 1];
			if (!zerok) {
				/*  Use a scaled form of recurrence if value of k at 0 is nonzero. */
				t = -aa / cc; //-p[n]/k[n-1] = -n*p[n]/p[n-1] ////may change
				for (i = 0; i < nm1; i++) {
					j = n - i - 1;// j=n-1;j>0;--j
					k[j] = t * k[j - 1] + p[j];
					//k[n-1]=t*k[n-2]+p[n-1] = t*(2*p[n-2]/n)+p[n-1];
					//k[n-2]=t*k[n-3]+p[n-2]= t*(3*p[n-3]/n)+p[n-2];
					//k[1]=t*k[0]+p[1]=t*p[0]+p[1]
				}
				k[0] = p[0];
				zerok = (fabs(k[n - 1]) <= fabs(bb) * machine_const[0] * 10.0);
			}
			else {
				/*  Use unscaled form of recurrence. */
				for (i = 0; i < nm1; i++) {
					j = n - i - 1;
					k[j] = k[j - 1];
				}
				k[0] = 0.0;
				zerok = (k[n - 1] == 0.0);
			}
		}
		/*  Save k for restarts with new shifts. */
		for (i = 0; i < n; i++)
			temp[i] = k[i];
		/*  Loop to select the quadratic corresponding to each new shift. */
		for (cnt = 0; cnt < 20; cnt++) {
			/*  Quadratic corresponds to a T shift to a
			 *  non-real point and its complex conjugate. The point
			 *  has modulus bnd and amplitude rotated by 94 degrees
			 *  from the previous shift.
			 */
			xxx = cosr * xx - sinr * yy;
			yy = sinr * xx + cosr * yy;
			xx = xxx;
			root_[0] = bnd * xx;
			root_[1] = bnd * yy;
			a[7] = -2.0 * root_[0];
			a[8] = bnd;
			fxshfr(20 * (cnt + 1), &nz, coe, root_, a, machine_const, p, qp, k, qk, svk, &n);
			if (nz != 0) {
				/*  The second stage jumps directly to one of the third
				 *  stage iterations and returns here if successful.
				 *  Deflate the polynomial, store the zero or zeros and
				 *  return to the main algorithm.
				 */
				j = degree - n;
				zeror[j] = root_[2];
				zeroi[j] = root_[3];
				n -= nz;
				for (i = 0; i <= n; i++)
					p[i] = qp[i];
				if (nz != 1) {
					zeror[j + 1] = root_[4];
					zeroi[j + 1] = root_[5];
				}
				goto _40;
			}
			/*  If the iteration is unsuccessful another quadratic
			 *  is chosen after restoring k.
			 */
			for (i = 0; i < n; i++) {
				k[i] = temp[i];
			}
		}
		/*  Return with failure if no convergence with 20 shifts. */
	_99:
		return degree - n;
	}


}
