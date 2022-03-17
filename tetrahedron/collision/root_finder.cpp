#include"root_finder.h"

int RootFinder::root_finder(double* op, int degree, double* zeror, double* zeroi)
{
	double t, aa, bb, cc, factor, rot;
	double lo, max, min, xx, yy, cosr, sinr, xxx, x, sc, bnd;
	double xm, ff, df, dx, infin, smalno, base;
	int cnt, nz, jj, l, nm1, zerok;

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
	for (unsigned int i = 0; i <= n; i++)
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
	for (unsigned int i = 0; i <= n; ++i) {
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



_110:
_99:


}