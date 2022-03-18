//real polynomial root finder.
#pragma once
#include<vector>
class RootFinder
{
public:
	int rpoly(double* op, int degree, double* zeror, double* zeroi);

private:
	double p[7], qp[7], k[7], qk[7], svk[7];
	//double v;
	int n;
	double eta, are, mre;
	

	void quad(double a, double b1, double c, double* sr, double* si,
		double* lr, double* li);

	/*  Computes up to L2 fixed shift k-polynomials,
	 *  testing for convergence in the linear or quadratic
	 *  case. Initiates one of the variable shift
	 *  iterations and returns with the number of zeros
	 *  found.
	 */
	void fxshfr(int l2, int* nz, double* coe, double* root_, double* a);
	void quadsd(int nn, double* u, double* v, double* p, double* q,
		double* a, double* b);
	void quadit(double* uu, double* vv, int* nz, double* coe, double* root_, double* a);
	void realit(double* sss, int* nz, int* iflag, double* root_);
	void calcsc(int* type, double* coe, double* a);
	void nextk(int* type, double* coe, double* a);
	void newest(int type, double* uu, double* vv, double* coe, double* a);
};