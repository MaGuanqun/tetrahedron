#include"root_finder.h"

int RootFinder::root_finder(double* op, int degree, double* zeror, double* zeroi)
{
	double t, aa, bb, cc, factor, rot;
	double lo, max, min, xx, yy, cosr, sinr, xxx, x, sc, bnd;
	double xm, ff, df, dx, infin, smalno, base;
	int cnt, nz, jj, i,j, l, nm1, zerok;

	//set machine constants
	base = 2.0;
	eta = 2.22e-16;
	infin = 3.4e38;
	smalno = 1.2e-38;

	are = eta;
	mre = eta;
	lo = smalno / eta;

	// Initialization of constants for shift rotation???
	xx = sqrt(0.5);
	yy = -xx;
	rot = 94.0;
	rot *= 0.017453293;
	cosr = cos(rot);
	sinr = sin(rot);
	n = degree;

	//Algorithm fails of the leading coefficient is zero
	if (op[0] == 0.0) return -1;

	//Remove the zeros at the origin
	while (op[n] == 0.0) {
		j = degree - n;
		zeror[j] = 0.0;
		zeroi[j] = 0.0;
		n--;
	}
	if (n < 1) {
		return -1;
	}

	double temp[7];
	double pt[7];
	//make a copy of the coefficients
	for (i = 0; i <= n; ++i)
		p[i] = op[i];
_40:
	if (n == 1) {
		zeror[degree - 1] = -p[1] / p[0];
		zeroi[degree - 1] = 0.0;
		n -= 1;
		goto _99;
	}
	if (n == 2) {
		quad(p[0], p[1], p[2], &zeror[degree - 2], &zeroi[degree - 2],
			&zeror[degree - 1], &zeroi[degree - 1]);
		n -= 2;
		goto _99;
	}
	max = 0.0;
	min = infin;
	for (i = 0; i <= n; ++i) {
		x = fabs(p[i]);
		if (x > max) {
			max = x;
		}
		if (x != 0.0 && x < min) {
			min = x;
		}
	}
	/*  Scale if there are large or very small coefficients.
		 *  Computes a scale factor to multiply the coefficients of the
		 *  polynomial. The scaling si done to avoid overflow and to
		 *  avoid undetected underflow interfering with the convergence
		 *  criterion. The factor is a power of the base.
		 */
	sc = lo / min;
	if (sc > 1.0 && infin < max* sc) goto _110; //when max>min*infi/low,max min differ by very many orders of magnitude
	if (sc <= 1.0) {
		if (max < 10.0) goto _110;
		if (sc == 0.0) {
			sc = smalno; //????????????????
		}
	}
	l = (int)(log(sc) / log(base) + 0.5); 
	factor = pow(base, l);//nearest 2^n
	if (factor != 1.0) {//sc>=1.0  this is when min is extremely small, smaller than lo
		for (i = 0; i <= n; ++i)
			p[i] *= factor;     // Scale polynomial. at least to the magnitude of lo
	}
_110:
	//  Compute lower bound on moduli of roots. 
	for (i = 0; i <= n; ++i) {
		pt[i] = fabs(p[i]);
	}
	pt[n] = -pt[n];
	//  Compute upper estimate of bound. 
	x = exp((log(-pt[n]) - log(pt[0])) / (double)n);  //  (-pt[n]/pt[0])^{1/n} ????  the result of pt[n]x^n+pt[0]=0?
	//this is the upper bound of the equation pt[n]x^n+pt[n-1]x^(n-1)+...+pt[0]. 
	//if newton step at the origin is better, use it. better means smaller
	if (pt[n - 1] != 0.0) {
		xm = -pt[n] / pt[n - 1];
		if (xm < x)  x = xm;
	}
	//  Chop the interval (0,x) until ff <= 0
	while (true) {
		xm = x * 0.1;
		ff = pt[0];
		for (i = 1; i <= n; ++i) //here we compute ((pt[0]*xm + pt[1])*xm+pt[2]+.... pt[n] is negative, all other pt are positive 
			ff = ff * xm + pt[i]; //
		if (ff <= 0.0) break;
		x = xm;
	}
	dx = x; //=10*xm 
	//  Do Newton interation until x converges to two decimal places.
	while (fabs(dx / x) > 0.005) {
		ff = pt[0];
		df = ff;
		for (i = 1; i < n; ++i) { 
			ff = ff * x + pt[i];
			df = df * x + ff; //df = ((x*pt[0]+x*pt[0]+pt[1])*x+(pt[0]*x + pt[1])*x+pt[2])*x+...  this is to compute derivative
		}
		ff = ff * x + pt[n];
		dx = ff / df;
		x -= dx;
	}
	bnd = x;
	//  Compute the derivative as the initial k polynomia and do 5 steps with no shift.
	nm1 = n - 1;
	for (i = 1; i < n; ++i) {
		k[i] = (double)(n - i) * p[i] / (double)n; //here is to compute the derivative. then divide all coefficients by n
	}
	k[0] = p[0];
	aa = p[n];
	bb = p[n - 1];
	zerok = (k[n - 1] == 0);
	for (jj = 0; jj < 5; jj++) {
		cc = k[n - 1];
		if (!zerok) {
			//  Use a scaled form of recurrence if value of k at 0 is nonzero. ???????????????
			t = -aa / cc;   //-p[n]/k[n-1] = -n*p[n]/p[n-1] ////may change
			for (i = 0; i < nm1; ++i) { 
				j = n - i - 1;  // j=n-1;j>0;--j
				k[j] = t * k[j - 1] + p[j]; 
				//k[n-1]=t*k[n-2]+p[n-1] = t*(2*p[n-2]/n)+p[n-1];
				//k[n-2]=t*k[n-3]+p[n-2]= t*(3*p[n-3]/n)+p[n-2];
				//k[1]=t*k[0]+p[1]=t*p[0]+p[1]
			}
			k[0] = p[0];
			zerok = (fabs(k[n - 1]) <= fabs(bb) * eta * 10.0);
		}
		else {
			/*  Use unscaled form of recurrence. */
			for (i = 0; i < nm1; ++i) {
				j = n - i - 1;
				k[j] = k[j - 1];
			}
			k[0] = 0.0;
			zerok = (k[n - 1] == 0.0);
		}
	}
	//  Save k for restarts with new shifts. 
	for (i = 0; i < n; i++)
		temp[i] = k[i];
	for (cnt = 0; cnt < 20; cnt++) {
		/*  Quadratic corresponds to a double shift to a
		 *  non-real point and its complex conjugate. The point
		 *  has modulus bnd and amplitude rotated by 94 degrees
		 *  from the previous shift.
		 */
		xxx = cosr * xx - sinr * yy;
		yy = sinr * xx + cosr * yy;
		xx = xxx;
		sr = bnd * xx;
		si = bnd * yy;
		u = -2.0 * sr;
		v = bnd;
		fxshfr(20 * (cnt + 1), &nz);
		if (nz != 0) {
			/*  The second stage jumps directly to one of the third
			 *  stage iterations and returns here if successful.
			 *  Deflate the polynomial, store the zero or zeros and
			 *  return to the main algorithm.
			 */
			j = degree - n;
			zeror[j] = szr;
			zeroi[j] = szi;
			n -= nz;
			for (i = 0; i <= n; i++)
				p[i] = qp[i];
			if (nz != 1) {
				zeror[j + 1] = lzr;
				zeroi[j + 1] = lzi;
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
_99:
	return degree - n;

}

void RootFinder::quad(double a, double b1, double c, double* sr, double* si,
	double* lr, double* li)
{
	double b, d, e;
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
	b = b1 * 0.5;
	if (fabs(b) < fabs(c)) {
		if (c < 0.0)
			e = -a;
		else
			e = a;
		e = b * (b / fabs(c)) - e; // if(c>0) (b^2-4ac)/c  if(c<0) (4ac-b^2)/c
		d = sqrt(fabs(e)) * sqrt(fabs(c));  //sqrt(b^2-4ac)
	}
	else { //avoid overflow.
		e = 1.0 - (a / b) * (c / b); // (b^2-4ac)/c
		d = sqrt(fabs(e)) * fabs(b); //sqrt(b^2-4ac)
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
		*lr = (-b + d) / a; //if (b>0) lr=(-b-sqrt(b^2-4ac))/a, if(b<0) lr=(-b+sqrt(b^2-4ac))/a  this is also to avoid subtract between two positive number
		*sr = 0.0;
		if (*lr != 0.0)
			*sr = (c / *lr) / a;
		*si = 0.0;
		*li = 0.0;
	}
}


/*  Computes up to L2 fixed shift k-polynomials,
	 *  testing for convergence in the linear or quadratic
	 *  case. Initiates one of the variable shift
	 *  iterations and returns with the number of zeros
	 *  found.
	 */
void RootFinder::fxshfr(int l2, int* nz)
{
	double svu, svv, ui, vi, s;
	double betas, betav, oss, ovv, ss, vv, ts, tv;
	double ots = 0, otv = 0, tvv, tss;
	int type, i, j, iflag, vpass, spass, vtry, stry;

	*nz = 0;
	betav = 0.25;
	betas = 0.25;
	oss = sr;
	ovv = v;
	/*  Evaluate polynomial by synthetic division. */
	quadsd(n, &u, &v, p, qp, &a, &b);
	calcsc(&type);
	for (j = 0; j < l2; j++) {
		/*  Calculate next k polynomial and estimate v. */
		nextk(&type);
		calcsc(&type);
		newest(type, &ui, &vi);
		vv = vi;
		/*  Estimate s. */
		ss = 0.0;
		if (k[n - 1] != 0.0) ss = -p[n] / k[n - 1];
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
		svu = u;
		svv = v;
		for (i = 0; i < n; i++) {
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
		quadit(&ui, &vi, nz);
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
		for (i = 0; i < n; i++) {
			k[i] = svk[i];
		}
	_40:
		realit(&s, nz, &iflag);
		if (*nz > 0) return;
		/*  Linear iteration has failed. Flag that it has been
		 *  tried and decrease the convergence criterion.
		 */
		stry = 1;
		betas *= 0.25;
		if (iflag == 0) goto _50;
		/*  If linear iteration signals an almost double real
		 *  zero attempt quadratic iteration.
		 */
		ui = -(s + s);
		vi = s * s;
		goto _20;
		/*  Restore variables. */
	_50:
		u = svu;
		v = svv;
		for (i = 0; i < n; i++) {
			k[i] = svk[i];
		}
		/*  Try quadratic iteration if it has not been tried
		 *  and the V sequence is convergin.
		 */
		if (vpass && !vtry) goto _20;
		/*  Recompute QP and scalar values to continue the
		 *  second stage.
		 */
		quadsd(n, &u, &v, p, qp, &a, &b);
		calcsc(&type);
	_70:
		ovv = vv;
		oss = ss;
		otv = tv;
		ots = ts;
	}
}

void RootFinder::quadsd(int nn, double* u, double* v, double* p, double* q,
	double* a, double* b)
{
	double c;
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