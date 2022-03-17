//real polynomial root finder.
#pragma once
#include<vector>
class RootFinder
{
public:
	int root_finder(double* op, int degree, double* zeror, double* zeroi);

private:
	double p[7], qp[7], k[7], qk[7], svk[7];
	double sr, si, u, v, a, b, c, d, a1, a2;
	double a3, a6, a7, e, f, g, h, szr, szi, lzr, lzi;
	double eta, are, mre;
	int n, nn, nmi, zerok;
};