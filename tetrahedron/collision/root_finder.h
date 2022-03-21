//real polynomial root finder.
#include <math.h>  
#pragma once
class RootFinder
{
public:
	int rpoly(double* op, int degree, double* zeror, double* zeroi);

private:
	void quad(double a, double b1, double c, double* sr, double* si,
		double* lr, double* li);

	/*  Computes up to L2 fixed shift k-polynomials,
	 *  testing for convergence in the linear or quadratic
	 *  case. Initiates one of the variable shift
	 *  iterations and returns with the number of zeros
	 *  found.
	 */
	void fxshfr(int l2, int* nz, double* coe, double* root_, double* a, double* machine_const, 
		double* p, double* qp, double* k, double* qk, double* svk, int* n);
	void quadsd(int nn, double* u, double* v, double* p, double* q,
		double* a, double* b);
	void quadit(double* uu, double* vv, int* nz, double* coe, double* root_, double* a, 
		double* machine_const, double* p, double* qp, double* k, double* qk, int* n);
	void realit(double* sss, int* nz, int* iflag, double* root_, double* machine_const, double* p, 
		double* qp, double* k, double* qk, int* n);
	void calcsc(int* type, double* coe, double* a, double* machine_const, double* k, double* qk, int* n);
	void nextk(int* type, double* coe, double* a, double* machine_const, double* qp, double* k, double* qk, int* n);
	void newest(int type, double* uu, double* vv, double* coe, double* a, double* p, double* k, int* n);
};