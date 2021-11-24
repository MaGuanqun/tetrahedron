#include"ApproxCCD.h"
#include"../basic/global.h"

bool ApproxCCD::pointTriangleCollisionTime(double& t, double* initial_position, double* current_position,
	double* initial_triangle_0, double* current_triangle_0, 
	double* initial_normal_not_normalized, double* current_normal_not_normalized, double* cross_for_CCD)
{
	double qp_0[3];
	SUB(qp_0, initial_position, initial_triangle_0);
	double d = DOT(qp_0, initial_normal_not_normalized);
	if (d<NEAR_ZERO2 && d>-NEAR_ZERO2) {
		t = 0;
		return true;
	}
	double qp_1[3];
	SUB(qp_1, current_position, current_triangle_0);
	double b3 = DOT(qp_1, current_normal_not_normalized);
	double b1 = DOT(qp_1, initial_normal_not_normalized) + DOT(qp_0, cross_for_CCD);
	double b2 = DOT(qp_1, cross_for_CCD) + DOT(qp_0, current_normal_not_normalized);
	double a3 = d + b1 - b2 + b3;
	double a2 = 3 * d - 2 * b1 + b2;
	double a1 = -3 * d + b1;
	return solveEquation(t, a3, a2, a1, d);
}

bool ApproxCCD::edgeEdgeCollisionTime(double& t, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
	double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
	double* initial_compare_edge_vertex_1)
{
	double e01[3], e01_[3], e23[3], e23_[3], e02[3], e02_[3];
	SUB(e01, initial_edge_vertex_0, initial_edge_vertex_1);
	SUB(e01_, current_edge_vertex_0, current_edge_vertex_1);
	SUB(e02, initial_edge_vertex_0, initial_compare_edge_vertex_0);
	SUB(e02_, current_edge_vertex_0, current_compare_edge_vertex_0);
	SUB(e23, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
	SUB(e23_, current_compare_edge_vertex_0, current_compare_edge_vertex_1);
	double c1[3], c2[3], c3[3], c4[3];
	CROSS(c1, e01, e23);
	CROSS(c2, e01_, e23);
	CROSS(c3, e01, e23_);
	CROSS(c4, e01_, e23_);
	SUM_(c2, c3);
	double d = DOT(c1, e02);
	double b1 = DOT(c1, e02_) + DOT(c2, e02);
	double b2 = DOT(c2, e02_) + DOT(c4, e02);
	double b3 = DOT(c4, e02_);
	double a3, a2, a1;
	a3 = -d + b1 - b2 + b3;
	a2 = 3 * d - 2 * b1 + b2;
	a1 = -3 * d + b1;
	return solveEquation(t, a3, a2, a1, d);
}


bool ApproxCCD::solveEquation(double&t, double a3, double a2, double a1, double d)
{
	if (d<NEAR_ZERO2 && d>-NEAR_ZERO2) {
		return solveQuadraticEquation(t, a3, a2, a1);
	}	
	if (a3<NEAR_ZERO2 && a3>-NEAR_ZERO2) {
		return solveQuadraticEquation(t, a2, a1, d);
	}
	double a_2 = 3 * a3;
	double a_1 = 2 * a2;
	double check_root = a_1 * a_1 - 4 * a_2 * a1;
	if (check_root > 0) {
		check_root = sqrt(check_root);
		if (d < 0) {
			a3 = -a3;
			a2 = -a2;
			a1 = -a1;
			d = -d;
			a_2 = -a_2;
			a_1 = -a_1;
		}
		double i_0, i_1;
		if (a3 > 0) {
			i_0 = (-a_1 - check_root) / (2 * a_2);
			i_1 = (-a_1 + check_root) / (2 * a_2);
			if (i_1 < 0) {
				return false;
			}
			if (i_0 > 1) {
				return false;
			}
			if (i_0 < 0) {
				if (i_1 > 1) {//interval [0, 1]
					if (d * (a3 + a2 + a1 + d) < 0) {
						double t_s = -a2 / a_2;
						t = -d / (a_2 * t_s * t_s + a_1 * t_s + a1);
						return true;
					}						
				}
				else {//interval [0, i_1]
					if (d * (a3 * i_1 * i_1 * i_1 + a2 * i_1 * i_1 + a1 * i_1 + d) < 0) {
						double t_s = -a2 / a_2;
						t = -d / (a_2 * t_s * t_s + a_1 * t_s + a1);
						return true;
					}
				}
				return false;
			}		
			double r_0 = a3 * i_0 * i_0 * i_0 + a2 * i_0 * i_0 + a1 * i_0 + d;
			if (i_1 > 1) {//interval [i_0, 1]
				if (r_0 * (a3 + a2 + a1 + d) < 0) {
					double t_s = -a2 / a_2;
					t = i_0 - r_0 / (a_2 * t_s * t_s + a_1 * t_s + a1);
					return true;
				}						
			}
			else {//interval [i_0, i_1]
				if (r_0 * (a3 * i_1 * i_1 * i_1 + a2 * i_1 * i_1 + a1 * i_1 + d) < 0) {
					double t_s = -a2 / a_2;
					t = i_0 - r_0 / (a_2 * t_s * t_s + a_1 * t_s + a1);
					return true;
				}
			}
			return false;			
		}
		else {
			i_0 = (-a_1 + check_root) / (2 * a_2);
			i_1 = (-a_1 - check_root) / (2 * a_2);
			if (i_0 > 0) {
				if (i_0 > 1) {//interval [0, 1]
					if (d * (a3 + a2 + a1 + d) < 0) {
						t = -d / a1;
						return true;
					}
				}
				else {//interval [0, i_0]
					if (d * (a3 * i_0 * i_0 * i_0 + a2 * i_0 * i_0 + a1 * i_0 + d) < 0) {
						t = -d / a1;
						return true;
					}
				}
				return false;
			}
			if (i_1 > 1) {
				return false;
			}
			if (i_1 < 0) { //interval [0,1]
				if (d * (a3 + a2 + a1 + d) < 0) {
					t = -d / (a_2 + a_1 + a1);
					return true;
				}
			}
			else { //interval [i_1,1]
				double r_i_1 = a3 * i_1 * i_1 * i_1 + a2 * i_1 * i_1 + a1 * i_1 + d;
				if (r_i_1 * (a3 + a2 + a1 + d) < 0) {
					t = i_1 - r_i_1 / (a_2 + a_1 + a1);
					return true;
				}
			}
			return false;		
		}
	}
	else {
		if (d * (a3 + a2 + a1 + d) < 0) {
			double r_1 = a_2 + a_1 + a1;
			if (abs(a1) > abs(r_1)) {
				t = -d / a1;
			}
			else {
				t = -d / r_1;
			}
			return true;
		}
		return false;
	}	
}


bool ApproxCCD::solveQuadraticEquation(double& t, double a2, double a1, double a0)
{
	if (a0<NEAR_ZERO2 && a0>-NEAR_ZERO2) {
		if (a2<NEAR_ZERO2 && a2>-NEAR_ZERO2) {
			return false;
		}
		t = -a1 / a2;
		if (t > 0.0 && t < 1.0) {
			return true;
		}
		else {
			return false;
		}
	}
	if (a2<NEAR_ZERO2 && a2>-NEAR_ZERO2) {
		if (a1<NEAR_ZERO2 && a1>-NEAR_ZERO2) {
			return false;
		}
		t = -a0 / a1;
		if (t > 0.0 && t < 1.0) {
			return true;
		}
		else {
			return false;
		}
	}
	else {
		double temp = a1 * a1 - 4 * a2 * a0;
		if (temp < 0) {
			return false;
		}
		else {
			temp = sqrt(temp);
			if (a2 < 0) {
				a2 = -a2;
				a1 = -a1;
				a0 = -a0;
			}
			t = (-a1 - temp) / (2 * a2);
			if (t > 1.0) {
				return false;
			}
			else if (t > 0 && t < 1.0) {
				return true;
			}
			else {
				t = (-a1 + temp) / (2 * a2);
				if (t > 0 && t < 1.0) {
					return true;
				}
				else {
					return false;
				}
			}
		}
	}
}

void ApproxCCD::test()
{
	double initial_pos[3] = { -1.01,0.001,1.01 };
	double current_pos[3] = { -1.01,-0.1,1.01 };
	double initial_triangle_0[3]= { 1.0,0.0,0.0 };
	double initial_triangle_1[3] = { -1.0,0.0,-1.0 };
	double initial_triangle_2[3] = { -1.0,0.0,1.0 };
	double current_triangle_0[3]= { 1.0,0.0,0.0 };
	double current_triangle_1[3]= { -1.0,0.0,-1.0 };
	double current_triangle_2[3]= { -1.0,0.0,1.0 };


	double initial_tri_normal[3]; double tri_normal[3];
	double cross_for_CCD[3]; double temp[3];
	double e1[3], e2[3];
	double e3[3], e4[3];
	SUB(e1, initial_triangle_1, initial_triangle_0);
	SUB(e2, initial_triangle_2, initial_triangle_0);
	CROSS(initial_tri_normal, e1, e2);
	SUB(e3, current_triangle_1, current_triangle_0);
	SUB(e4, current_triangle_2, current_triangle_0);
	CROSS(tri_normal, e3, e4);
	CROSS(cross_for_CCD, e1, e4);
	CROSS(temp, e3, e2);
	SUM_(cross_for_CCD, temp);
	double t;

	if (pointTriangleCollisionTime(t, initial_pos, current_pos, initial_triangle_0, current_triangle_0, initial_tri_normal, tri_normal, cross_for_CCD)) {
		std::cout << t << std::endl;
	}
	else {
		std::cout << "does not collide " << std::endl;
	}
}