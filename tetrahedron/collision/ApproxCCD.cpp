#include"ApproxCCD.h"
#include"../basic/global.h"


bool ApproxCCD::pointTriangleCollisionTime(double& t, double* initial_position, double* current_position,
	double* initial_triangle_0, double* current_triangle_0, double* initial_triangle_1, double* current_triangle_1, double* initial_triangle_2, double* current_triangle_2,
	double* initial_normal_not_normalized, double* current_normal_not_normalized, double* cross_for_CCD, double tolerance_2)
{
	double e_0[3], e_1[3], e_2[3];
	SUB(e_0, initial_position, initial_triangle_0);
	SUB(e_1, initial_position, initial_triangle_1);
	SUB(e_2, initial_position, initial_triangle_2);

	double d = DOT(e_0, initial_normal_not_normalized);

	double u[3], u_0[3], u_1[3], u_2[3];
	SUB(u, current_position, initial_position);
	SUB(u_0, current_triangle_0, initial_triangle_0);
	SUB(u_1, current_triangle_1, initial_triangle_1);
	SUB(u_2, current_triangle_2, initial_triangle_2);

	if (d<NEAR_ZERO2 && d>-NEAR_ZERO2) {		
		double e_0_1[3], e_1_2[3], e_2_0[3];
		CROSS(e_0_1, e_0, e_1);
		CROSS(e_1_2, e_1, e_2);
		CROSS(e_2_0, e_2, e_0);
		if (DOT(e_0_1, e_1_2) > 0 && DOT(e_1_2, e_2_0) > 0 && DOT(e_0_1, e_2_0) > 0) {
			t = 0;
			return true;
		}		
	}
	double qp_1[3];
	SUB(qp_1, current_position, current_triangle_0);
	double b1= DOT(qp_1, current_normal_not_normalized);
	double b2 = DOT(qp_1, initial_normal_not_normalized);
	double b3= DOT(e_0, cross_for_CCD);
	double b4 = DOT(qp_1, cross_for_CCD);
	double b5= DOT(e_0, current_normal_not_normalized);

	double a3 = -d + b1 +b2 + b3-b4-b5;
	double a2 = 3 * d + b5 - 2 * b3 + b4 - 2 * b2;
	double a1 = -3 * d + b2 + b3;
	std::vector<double>time;
	time.reserve(7);
	if (solveEquation(t, a3, a2, a1, d)) {
		time.push_back(t);
	}
	double e_[3];
	SUB(e_, initial_triangle_1, initial_triangle_0);
	if (pointEdgeCollisionTime(t, u, u_0, u_1, e_, e_0, e_1,tolerance_2)) {
		time.push_back(t);
	}
	SUB(e_, initial_triangle_2, initial_triangle_1);
	if (pointEdgeCollisionTime(t, u, u_1, u_2, e_, e_1, e_2, tolerance_2)) {
		time.push_back(t);
	}
	SUB(e_, initial_triangle_2, initial_triangle_0);
	if (pointEdgeCollisionTime(t, u, u_0, u_2, e_, e_0, e_2, tolerance_2)) {
		time.push_back(t);
	}
	if (pointPointCollisionTime(t, e_0, u_0, u, tolerance_2)) {
		time.push_back(t);
	}
	if (pointPointCollisionTime(t, e_1, u_1, u, tolerance_2)) {
		time.push_back(t);
	}
	if (pointPointCollisionTime(t, e_2, u_2, u, tolerance_2)) {
		time.push_back(t);
	}
	if (time.empty()) {
		return false;
	}
	t = time[0];
	for (int i = 1; i < time.size(); ++i) {
		if (t > time[i]) {
			t = time[i];
		}
	}

	//if (t < 1e-6) {
	//	std::cout << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
	//}


	return true;
}

bool ApproxCCD::edgeEdgeCollisionTime(double& t, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
	double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
	double* initial_compare_edge_vertex_1, double tolerance_2)
{
	double e0[3], e0_[3], e1[3], e1_[3], e2[3], e2_[3];
	SUB(e0, initial_compare_edge_vertex_0, initial_edge_vertex_0);
	SUB(e0_, current_compare_edge_vertex_0, current_edge_vertex_0);
	SUB_(e0_, e0);

	SUB(e1, initial_edge_vertex_1, initial_edge_vertex_0);
	SUB(e1_, current_edge_vertex_1, current_edge_vertex_0);
	SUB_(e1_, e1);

	SUB(e2, initial_compare_edge_vertex_1, initial_compare_edge_vertex_0);
	SUB(e2_, current_compare_edge_vertex_1, current_compare_edge_vertex_0);
	SUB_(e2_, e2);

	double a[3];
	CROSS(a, e1, e2);	
	double d = DOT(e0, a);

	if (d<NEAR_ZERO2 && d>-NEAR_ZERO2) {
		//check when t=0,the two edges are coplanar if the two edges are collide
		double e_temp[3], compare_direction0[3], compare_direction1[3];
		SUB(e_temp, initial_compare_edge_vertex_0, initial_edge_vertex_1);
		CROSS(compare_direction0, e2, e_temp);
		CROSS(compare_direction1, e2, e0);
		double d2 = DOT(compare_direction0, compare_direction1);
		SUB(e_temp, initial_compare_edge_vertex_1, initial_edge_vertex_0);
		CROSS(compare_direction0, e1, e_temp);
		CROSS(compare_direction1, e1, e0);
		double d1 = DOT(compare_direction0, compare_direction1);		
		if ( d1 < 0 && d2<0) { //this means the vertices of one edge is on different sides of the other edge
			t = 0;
			return true;
		}
		// initial_compare_edge_vertex_0  initial_edge
		if (pointEdgeIsClose(e0, e1, tolerance_2)) {
			t = 0;
			return true;
		}
		// initial_compare_edge_vertex_1  initial_edge
		if (pointEdgeIsClose(e_temp, e1, tolerance_2)) {
			t = 0;
			return true;
		}
		//initial_edge_vertex_0 initial_compare_edge
		SUB(e_temp, initial_edge_vertex_0, initial_compare_edge_vertex_0);
		if (pointEdgeIsClose(e_temp, e2, tolerance_2)) {
			t = 0;
			return true;
		}
		SUB(e_temp, initial_edge_vertex_1, initial_compare_edge_vertex_0);
		if (pointEdgeIsClose(e_temp, e2, tolerance_2)) {
			t = 0;
			return true;
		}
	}

	double b[3], c[3];
	CROSS(b, e1, e2_);
	CROSS(c, e1_, e2);
	SUM_(b, c);
	CROSS(c, e1_, e2_);
	double a3 = DOT(e0_, c);
	double a2 = DOT(e0, c) + DOT(e0_, b);
	double a1 = DOT(e0, b) + DOT(e0_, a);
	std::vector<double>time;
	time.reserve(7);
	if (solveEquation(t, a3, a2, a1, d)) {
		time.push_back(t);
	}

	double u_0_0[3], u_0_1[3], u_1_0[3], u_1_1[3];
	SUB(u_0_0, current_edge_vertex_0, initial_edge_vertex_0);
	SUB(u_0_1, current_edge_vertex_1, initial_edge_vertex_1);
	SUB(u_1_0, current_compare_edge_vertex_0, initial_compare_edge_vertex_0);
	SUB(u_1_1, current_compare_edge_vertex_1, initial_compare_edge_vertex_1);

	double v_0[3], v_1[3];
	//initial_edge_vertex_0
	SUB(v_0, initial_edge_vertex_0, initial_compare_edge_vertex_0);
	SUB(v_1, initial_edge_vertex_0, initial_compare_edge_vertex_1);
	if (pointEdgeCollisionTime(t, u_0_0, u_1_0, u_1_1, e2, v_0, v_1, tolerance_2)) {
		time.push_back(t);
	}
	if (pointPointCollisionTime(t, v_0, u_1_0, u_0_0, tolerance_2)) {
		time.push_back(t);
	}
	if (pointPointCollisionTime(t, v_1, u_1_1, u_0_0, tolerance_2)) {
		time.push_back(t);
	}

	//initial_edge_vertex_1
	SUB(v_0, initial_edge_vertex_1, initial_compare_edge_vertex_0);
	SUB(v_1, initial_edge_vertex_1, initial_compare_edge_vertex_1);
	if (pointEdgeCollisionTime(t, u_0_1, u_1_0, u_1_1, e2, v_0, v_1, tolerance_2)) {
		time.push_back(t);
	}
	if (pointPointCollisionTime(t, v_0, u_1_0, u_0_1, tolerance_2)) {
		time.push_back(t);
	}
	if (pointPointCollisionTime(t, v_1, u_1_1, u_0_1, tolerance_2)) {
		time.push_back(t);
	}

	//initial_compare_edge_vertex_0
	SUB(v_0, initial_compare_edge_vertex_0, initial_edge_vertex_0);
	SUB(v_1, initial_compare_edge_vertex_0, initial_edge_vertex_1);
	if (pointEdgeCollisionTime(t, u_1_0, u_0_0, u_0_1, e1, v_0, v_1, tolerance_2)) {
		time.push_back(t);
	}
	//initial_compare_edge_vertex_1
	SUB(v_0, initial_compare_edge_vertex_1, initial_edge_vertex_0);
	SUB(v_1, initial_compare_edge_vertex_1, initial_edge_vertex_1);
	if (pointEdgeCollisionTime(t, u_1_1, u_0_0, u_0_1, e1, v_0, v_1, tolerance_2)) {
		time.push_back(t);
	}
	if (time.empty()) {
		return false;
	}
	t = time[0];
	for (int i = 1; i < time.size(); ++i) {
		if (t > time[i]) {
			t = time[i];
		}
	}
	return true;
}


bool ApproxCCD::solveEquation(double&t, double a3, double a2, double a1, double d)
{
	//std::cout << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
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

			//std::cout << "i_0 " << i_0 << " " << i_1 << std::endl;
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
						if (t_s < 0) {
							t = -d / a1;
						}
						else {
							if (t_s < 1) {
								t = -d / (a_2 * t_s * t_s + a_1 * t_s + a1);
							}
							else {
								t = -d / (a_2 + a_1 + a1);
							}
						}						
						return true;
					}						
				}
				else {//interval [0, i_1]
					if (d * (a3 * i_1 * i_1 * i_1 + a2 * i_1 * i_1 + a1 * i_1 + d) < 0) {
						double t_s = -a2 / a_2;
						if (t_s < 0) {
							t = -d / a1;
						}
						else {
							t = -d / (a_2 * t_s * t_s + a_1 * t_s + a1);
						}						
						return true;
					}
				}
				return false;
			}		
			double r_0 = a3 * i_0 * i_0 * i_0 + a2 * i_0 * i_0 + a1 * i_0 + d;
			if (i_1 > 1) {//interval [i_0, 1]
				if (r_0 * (a3 + a2 + a1 + d) < 0) {
					double t_s = -a2 / a_2;
					if (t_s > 1) {
						t = -d / (a_2 + a_1 + a1);
					}
					else {
						t = i_0 - r_0 / (a_2 * t_s * t_s + a_1 * t_s + a1);
					}					
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
			if (fabs(a1) > fabs(r_1)) {
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


bool ApproxCCD::pointPointCollisionTime(double& t, double* e_1_0, double* e_0, double* e_1, double tolerance_2)
{
	std::vector<double> t_list;
	t_list.reserve(3);
	double t_temp;
	//x
	t_temp = e_1_0[0] / (e_0[0] - e_1[0]);
	if (t_temp > 0 && t_temp <= 1) {
		t_list.push_back(t_temp);
	}
	//y
	t_temp = e_1_0[1] / (e_0[1] - e_1[1]);
	if (t_temp > 0 && t_temp <= 1) {
		t_list.push_back(t_temp);
	}
	//z
	t_temp = e_1_0[2] / (e_0[2] - e_1[2]);
	if (t_temp > 0 && t_temp <= 1) {
		t_list.push_back(t_temp);
	}

	if (t_list.empty()) {
		return false;
	}
	if (t_list.size() == 1) {
		if (pointPointIsClose(t_list[0], e_1_0, e_0, e_1, tolerance_2)) {
			t = t_list[0];
			return true;
		}
		return false;
	}
	std::sort(t_list.begin(), t_list.end());
	for (int i = 0; i < t_list.size(); ++i) {
		if (pointPointIsClose(t_list[i], e_1_0, e_0, e_1, tolerance_2)) {
			t = t_list[i];
			return true;
		}
	}
	return false;
}


void ApproxCCD::pointEdge2D(std::vector<double>&t, double* u, double* u0, double* u1, double* e_1_0, double* e0, double* e1, 
	double* u_0, double* u_1, double* u_10)
	// u is v(t1)-v(t0), e_1_0 is e_1(t0)-e_0(t0), e0 is v(t0)-e_0(t0)	u_0=v(t1)-v(t0)-(e_0(t1)-e_0(t0))
	//u_1_0 is e_1(t1)-e_1(t0)-(e_0(t1)-e_0(t0))
{

	double a, b, c;
	double t_temp[2];
	int root_number;
	double v_0[2], v_1[2], v_0_1[2];
	double x0t_x01t;
	double x01t_x01t;
	 

	//XY
	a = u_0[0] * u_10[1] - u_0[1] * u_10[0];
	b = u_0[0] * e_1_0[1] - u_0[1] * e_1_0[0]
		+ e0[0] * u_10[1] - e0[1] * u_10[0];
	c = e0[0] * e_1_0[1] - e0[1] * e_1_0[0];



	solveQuadratic(t_temp, a, b, c, root_number);
	for (int i = 0; i < root_number; ++i) {
		v_0[0] = e0[0] + t_temp[i] * u_0[0];
		v_0[1] = e0[1] + t_temp[i] * u_0[1];
		v_0_1[0] = e_1_0[0] + t_temp[i] * u_10[0];
		v_0_1[1] = e_1_0[1] + t_temp[i] * u_10[1];
		x0t_x01t = v_0[0] * v_0_1[0] + v_0[1] * v_0_1[1];
		x01t_x01t= v_0_1[0] * v_0_1[0] + v_0_1[1] * v_0_1[1];
		if (x0t_x01t >= 0 && x0t_x01t <= x01t_x01t) {
			t.push_back(t_temp[i]);
		}
	}	

	//YZ
	a = u_0[2] * u_10[1] - u_0[1] * u_10[2];
	b = u_0[2] * e_1_0[1] - u_0[1] * e_1_0[2]
		+ e0[2] * u_10[1] - e0[1] * u_10[2];
	c = e0[2] * e_1_0[1] - e0[1] * e_1_0[2];
	solveQuadratic(t_temp, a, b, c, root_number);
	for (int i = 0; i < root_number; ++i) {
		v_0[0] = e0[1] + t_temp[i] * u_0[1];
		v_0[1] = e0[2] + t_temp[i] * u_0[2];
		v_0_1[0] = e_1_0[1] + t_temp[i] * u_10[1];
		v_0_1[1] = e_1_0[2] + t_temp[i] * u_10[2];
		x0t_x01t = v_0[0] * v_0_1[0] + v_0[1] * v_0_1[1];
		x01t_x01t = v_0_1[0] * v_0_1[0] + v_0_1[1] * v_0_1[1];
		if (x0t_x01t >= 0 && x0t_x01t <= x01t_x01t) {
			t.push_back(t_temp[i]);
		}
	}

	//XZ
	a = u_0[2] * u_10[0] - u_0[0] * u_10[2];
	b = u_0[2] * e_1_0[0] - u_0[0] * e_1_0[2]
		+ e0[2] * u_10[0] - e0[0] * u_10[2];
	c = e0[2] * e_1_0[0] - e0[0] * e_1_0[2];
	solveQuadratic(t_temp, a, b, c, root_number);
	for (int i = 0; i < root_number; ++i) {
		v_0[0] = e0[0] + t_temp[i] * u_0[0];
		v_0[1] = e0[2] + t_temp[i] * u_0[2];
		v_0_1[0] = e_1_0[0] + t_temp[i] * u_10[0];
		v_0_1[1] = e_1_0[2] + t_temp[i] * u_10[2];
		x0t_x01t = v_0[0] * v_0_1[0] + v_0[1] * v_0_1[1];
		x01t_x01t = v_0_1[0] * v_0_1[0] + v_0_1[1] * v_0_1[1];
		if (x0t_x01t >= 0 && x0t_x01t <= x01t_x01t) {
			t.push_back(t_temp[i]);
		}
	}
}

void ApproxCCD::solveQuadratic(double* t, double a, double b, double c, int& root_number)
{
	root_number = 0;
	if (a<NEAR_ZERO2 && a>-NEAR_ZERO2) {
		if (b<NEAR_ZERO2 && b>-NEAR_ZERO2) {
			return;			
		}
		double time = -c / b;
		if (0 < time && time <= 1) {
			t[root_number++] = time;
		}		
		return;
	}
	if (a < 0) { a = -a; b = -b; c = -c; }	
	double delta = b * b - 4 * a * c;
	if (delta <= 0){
		if (-b > 0 && -b < 2 * a)  t[root_number++] = -b / (2 * a);
		return;
	}
	if (b <= 0){
		double temp = -b + sqrt(delta);
		double twice_c = 2 * c;
		double twice_a = 2 * a;
		if (twice_c > 0 && twice_c < temp)	t[root_number++] = twice_c / temp;
		if (temp < twice_a)				t[root_number++] = temp / twice_a;
	}
	else{
		double temp = -b - sqrt(delta);
		double twice_c = 2 * c;
		double twice_a = 2 * a;
		if (twice_a < temp)				t[root_number++] = temp / twice_a;
		if (twice_c < 0 && temp < twice_c)	t[root_number++] = twice_c / temp;
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


//log:
//Here, we first test vertex-edge and vertex-vertex. 
//

//edge-edge
void ApproxCCD::test()
{
	////edge 1
	//double initial_pos[3] = { 1.0,1.0,1.0 };
	//double initial_triangle_0[3] = { -1.0,0.0,-1.0 };
	//
	////edge 2
	//double initial_triangle_1[3] = { -1.0,0.0,1.0 };
	//double initial_triangle_2[3] = { 1.0,0.0,-1.0 };

	//double current_pos[3] = { 1.0,0.3,1.0 };
	//double current_triangle_0[3] = { -1.0,-0.1,-1.0 };

	//double current_triangle_1[3] = { -0.5,0.5,1.0 };
	//double current_triangle_2[3] = { 0.5,0.7,-1.0 };
	//edge 1
	double initial_pos[3] = { 1.0,1.0,1.0 };
	double initial_triangle_0[3] = { -1.0,0.0,-1.0 };

	//edge 2
	double initial_triangle_1[3] = { -1,-0.5,0.3 };
	double initial_triangle_2[3] = { 0.0,-0.5,-0.8 };

	double current_pos[3] = { -0.5,0.3,1.0 };
	double current_triangle_0[3] = { -1.0,-0.1,-0.5 };

	double current_triangle_1[3] = { -0.5,0.5,0.0 };
	double current_triangle_2[3] = { 0.5,0.5,-1.0 };
	double t;
	double e_0[3], e_1[3], e_2[3];
	SUB(e_0, initial_pos, initial_triangle_0);
	SUB(e_1, initial_pos, initial_triangle_1);
	SUB(e_2, initial_pos, initial_triangle_2);
	double u[3], u_0[3], u_1[3], u_2[3];
	SUB(u, current_pos, initial_pos);
	SUB(u_0, current_triangle_0, initial_triangle_0);
	SUB(u_1, current_triangle_1, initial_triangle_1);
	SUB(u_2, current_triangle_2, initial_triangle_2);
	double e_[3];
	SUB(e_, initial_triangle_1, initial_triangle_0);
	double tolerance_2 = 1e-6;
	if (edgeEdgeCollisionTime(t, current_pos, current_triangle_0, initial_pos, initial_triangle_0,
		current_triangle_1, current_triangle_2, initial_triangle_1, initial_triangle_2, tolerance_2)) {
		std::cout << t << std::endl;
		double result_v[3];
		double result_p0[3];
		double result_p1[3];
		double result_p2[3];

		POINT_ON_EDGE(result_v, 1.0, t, initial_pos, u);
		POINT_ON_EDGE(result_p0, 1.0, t, initial_triangle_0, u_0);
		POINT_ON_EDGE(result_p1, 1.0, t, initial_triangle_1, u_1);
		POINT_ON_EDGE(result_p2, 1.0, t, initial_triangle_2, u_2);

		double e_result1[3], e_result2[3];
		SUB(e_result1, result_v, result_p0);
		SUB(e_result2, result_p1, result_p2);
		double norm_result[3];
		CROSS(norm_result, e_result1, e_result2);

		SUB(e_result1, initial_pos, initial_triangle_0);
		SUB(e_result2, initial_triangle_1, initial_triangle_2);
		double norm_initial[3];
		CROSS(norm_initial, e_result1, e_result2);

		double e_result[3];
		SUB(e_result, result_v, result_p1);
		double e_initial[3];
		SUB(e_initial, initial_pos, initial_triangle_1);
		std::cout << "side indicate " << DOT(e_initial, norm_initial) << " " << DOT(e_result, norm_result)<<" "<< DOT(norm_result, norm_result) << std::endl;
		if (DOT(e_initial, norm_initial) * DOT(e_result, norm_result) < 0) {					
			std::cout << "side has changed " << std::endl;
		}
	}
	else {
		std::cout << "does not collide " << std::endl;
	}
	Vec3d initial_pos_a = Vec3d(initial_pos[0], initial_pos[1], initial_pos_a[2]);
	Vec3d current_pos_a = Vec3d(current_pos[0], current_pos[1], current_pos[2]);
	Vec3d initial_triangle_0_a = Vec3d(initial_triangle_0[0], initial_triangle_0[1], initial_triangle_0[2]);
	Vec3d initial_triangle_1_a = Vec3d(initial_triangle_1[0], initial_triangle_1[1], initial_triangle_1[2]);
	Vec3d initial_triangle_2_a = Vec3d(initial_triangle_2[0], initial_triangle_2[1], initial_triangle_2[2]);
	Vec3d current_triangle_0_a = Vec3d(current_triangle_0[0], current_triangle_0[1], current_triangle_0[2]);
	Vec3d current_triangle_1_a = Vec3d(current_triangle_1[0], current_triangle_1[1], current_triangle_1[2]);
	Vec3d current_triangle_2_a = Vec3d(current_triangle_2[0], current_triangle_2[1], current_triangle_2[2]);

	rootparity::RootParityCollisionTest root_parity(initial_pos_a, initial_triangle_0_a, initial_triangle_1_a, initial_triangle_2_a,
		current_pos_a, current_triangle_0_a, current_triangle_1_a, current_triangle_2_a, true);
	if (root_parity.edge_edge_collision())
	{
		std::cout << "actually collide " << std::endl;
	}
}
//vertex-triangle
//void ApproxCCD::test()
//{
//	double initial_pos[3] = { -0.5,0.1,-0.3 };
//	double current_pos[3] = { 0.5,-0.1,0.3 };
//	double initial_triangle_0[3] = { 1.0,0.0,0.0 };
//	double initial_triangle_1[3] = { -1.0,-0.1,-1.0 };
//	double initial_triangle_2[3] = { -1.0,0.1,1.0 };
//	double current_triangle_0[3] = { 1.0,0.0,0.0 };
//	double current_triangle_1[3] = { -1.0,0.1,-0.8 };
//	double current_triangle_2[3] = { -1.0,-0.1,0.8 };
//	double initial_tri_normal[3]; double tri_normal[3];
//	double cross_for_CCD[3]; double temp[3];
//	double e1[3], e2[3];
//	double e3[3], e4[3];
//	SUB(e1, initial_triangle_1, initial_triangle_0);
//	SUB(e2, initial_triangle_2, initial_triangle_0);
//	CROSS(initial_tri_normal, e1, e2);
//	SUB(e3, current_triangle_1, current_triangle_0);
//	SUB(e4, current_triangle_2, current_triangle_0);
//	CROSS(tri_normal, e3, e4);
//	CROSS(cross_for_CCD, e1, e4);
//	CROSS(temp, e3, e2);
//	SUM_(cross_for_CCD, temp);
//	double t;
//	double e_0[3], e_1[3], e_2[3];
//	SUB(e_0, initial_pos, initial_triangle_0);
//	SUB(e_1, initial_pos, initial_triangle_1);
//	SUB(e_2, initial_pos, initial_triangle_2);
//	double u[3], u_0[3], u_1[3], u_2[3];
//	SUB(u, current_pos, initial_pos);
//	SUB(u_0, current_triangle_0, initial_triangle_0);
//	SUB(u_1, current_triangle_1, initial_triangle_1);
//	SUB(u_2, current_triangle_2, initial_triangle_2);
//	double e_[3];
//	SUB(e_, initial_triangle_1, initial_triangle_0);
//	double tolerance_2 = 1e-6;
//	if (pointTriangleCollisionTime(t, initial_pos, current_pos, initial_triangle_0, current_triangle_0, 
//		initial_triangle_1, current_triangle_1, initial_triangle_2, current_triangle_2,
//		initial_tri_normal, tri_normal, cross_for_CCD,tolerance_2)) {
//		std::cout << t << std::endl;
//
//		double result_v[3];
//		double result_p0[3];
//		double result_p1[3];
//		double result_p2[3];
//
//		POINT_ON_EDGE(result_v, 1.0, t, initial_pos, u);
//		POINT_ON_EDGE(result_p0, 1.0, t, initial_triangle_0, u_0);
//		POINT_ON_EDGE(result_p1, 1.0, t, initial_triangle_1, u_1);
//		POINT_ON_EDGE(result_p2, 1.0, t, initial_triangle_2, u_2);
//
//		double e_result1[3], e_result2[3];
//		SUB(e_result1, result_p1, result_p0);
//		SUB(e_result2, result_p2, result_p0);
//		double norm_result[3];
//		CROSS(norm_result, e_result1, e_result2);
//
//		double e_result[3];
//		SUB(e_result, result_v, result_p0);
//		std::cout <<"side indicate "<< DOT(e_0, initial_tri_normal) << " " << DOT(e_result, norm_result) << std::endl;
//		if (DOT(e_0, initial_tri_normal) * DOT(e_result, norm_result) < 0) {
//			
//			std::cout << "side has changed " << std::endl;
//		}
//	}
//	else {
//		std::cout << "does not collide " << std::endl;
//	}
//	Vec3d initial_pos_a = Vec3d(initial_pos[0], initial_pos[1], initial_pos_a[2]);
//	Vec3d current_pos_a = Vec3d(current_pos[0], current_pos[1], current_pos[2]);
//	Vec3d initial_triangle_0_a = Vec3d(initial_triangle_0[0], initial_triangle_0[1], initial_triangle_0[2]);
//	Vec3d initial_triangle_1_a = Vec3d(initial_triangle_1[0], initial_triangle_1[1], initial_triangle_1[2]);
//	Vec3d initial_triangle_2_a = Vec3d(initial_triangle_2[0], initial_triangle_2[1], initial_triangle_2[2]);
//	Vec3d current_triangle_0_a = Vec3d(current_triangle_0[0], current_triangle_0[1], current_triangle_0[2]);
//	Vec3d current_triangle_1_a = Vec3d(current_triangle_1[0], current_triangle_1[1], current_triangle_1[2]);
//	Vec3d current_triangle_2_a = Vec3d(current_triangle_2[0], current_triangle_2[1], current_triangle_2[2]);
//
//	rootparity::RootParityCollisionTest root_parity(initial_pos_a, initial_triangle_0_a, initial_triangle_1_a, initial_triangle_2_a,
//		current_pos_a, current_triangle_0_a, current_triangle_1_a, current_triangle_2_a,false);
//	if (root_parity.point_triangle_collision())
//	{
//		std::cout << "actually collide " << std::endl;
//	}
//}

//void ApproxCCD::test()
//{
//	double initial_pos[3] = { -1.0001,0.1,0.0 };
//	double current_pos[3] = { -1.0001,-0.1,0.0 };
//	double initial_triangle_0[3]= { 1.0,0.0,0.0 };
//	double initial_triangle_1[3] = { -1.0,-0.1,-1.0 };
//	double initial_triangle_2[3] = { -1.0,0.1,1.0 };
//	double current_triangle_0[3]= { 1.0,0.0,0.0 };
//	double current_triangle_1[3]= { -1.0,0.1,-0.8 };
//	double current_triangle_2[3]= { -1.0,-0.1,0.8 };
//	double initial_tri_normal[3]; double tri_normal[3];
//	double cross_for_CCD[3]; double temp[3];
//	double e1[3], e2[3];
//	double e3[3], e4[3];
//	SUB(e1, initial_triangle_1, initial_triangle_0);
//	SUB(e2, initial_triangle_2, initial_triangle_0);
//	CROSS(initial_tri_normal, e1, e2);
//	SUB(e3, current_triangle_1, current_triangle_0);
//	SUB(e4, current_triangle_2, current_triangle_0);
//	CROSS(tri_normal, e3, e4);
//	CROSS(cross_for_CCD, e1, e4);
//	CROSS(temp, e3, e2);
//	SUM_(cross_for_CCD, temp);
//	double t;
//	double e_0[3], e_1[3], e_2[3];
//	SUB(e_0, initial_pos, initial_triangle_0);
//	SUB(e_1, initial_pos, initial_triangle_1);
//	SUB(e_2, initial_pos, initial_triangle_2);
//	double u[3], u_0[3], u_1[3], u_2[3];
//	SUB(u, current_pos, initial_pos);
//	SUB(u_0, current_triangle_0, initial_triangle_0);
//	SUB(u_1, current_triangle_1, initial_triangle_1);
//	SUB(u_2, current_triangle_2, initial_triangle_2);
//	double e_[3];
//	SUB(e_, initial_triangle_1, initial_triangle_0);
//	tolerance_2 = 1e-6;
//	if (pointEdgeCollisionTime(t, u, u_0, u_1, e_, e_0, e_1)) {
//		std::cout << "v e01 " << t<<std::endl;
//	}
//	SUB(e_, initial_triangle_2, initial_triangle_1);
//	if (pointEdgeCollisionTime(t, u, u_1, u_2, e_, e_1, e_2)) {
//		std::cout << "v e12 " << t << std::endl;
//	}
//	SUB(e_, initial_triangle_2, initial_triangle_0);
//	if (pointEdgeCollisionTime(t, u, u_0, u_2, e_, e_0, e_2)) {
//		std::cout << "v e02 " << t << std::endl;
//	}
//	if (pointPointCollisionTime(t, e_0, u_0, u)) {
//		std::cout << "v v0 " << t << std::endl;
//	}
//	if (pointPointCollisionTime(t, e_1, u_1, u)) {
//		std::cout << "v v1 " << t << std::endl;
//	}
//	if (pointPointCollisionTime(t, e_2, u_2, u)) {
//		std::cout << "v v2 " << t << std::endl;
//	}
//}

bool ApproxCCD::pointEdgeCollisionTime(double& t, double* u, double* u0, double* u1, double* e_1_0,
	double* e0, double* e1, double tolerance_2)
	// u is v(t1)-v(t0), e_1_0 is e_1(t0)-e_0(t0), e0 is v(t0)-e_0(t0)	//u_0=v(t1)-v(t0)-(e_0(t1)-e_0(t0))

{
	double u_0[3], u_1[3], u_10[3];
	SUB(u_0, u, u0);
	SUB(u_1, u, u1);
	SUB(u_10, u1, u0);


	std::vector<double>time;
	pointEdge2D(time, u, u0, u1, e_1_0, e0, e1, u_0, u_1, u_10);
	if (time.empty()) {
		return false;
	}
	std::sort(time.begin(), time.end());
	double v_0[3], v_1_0[3];
	for (int i = 0; i < time.size(); ++i) {
		t = time[i];
		v_0[0] = e0[0] + t * u_0[0];
		v_0[1] = e0[1] + t * u_0[1];
		v_0[2] = e0[2] + t * u_0[2];	
		v_1_0[0] = e_1_0[0] + t * u_10[0];
		v_1_0[1] = e_1_0[1] + t * u_10[1];
		v_1_0[2] = e_1_0[2] + t * u_10[2];		
		if (pointEdgeIsClose(v_0,v_1_0,tolerance_2)) {
			return true;
		}
	}
	return false;
}


bool ApproxCCD::pointEdgeIsClose(double* v_0, double* v_1_0, double tolerance_2)
//v_0 is v-e_0, v_1_0 is e_1-e_0
{
	double l_v_0 = DOT(v_0, v_0);
	double l_v_0_ = DOT(v_0, v_1_0);
	double l_v_1_0 = DOT(v_1_0, v_1_0);
	if (l_v_0 - l_v_0_ * l_v_0_ / l_v_1_0 > tolerance_2) {//use Pythagorean Theorem d^2=v_0^2-(v_0*v_1_0/||v_1_0||)^2
		return false;
	}	
	if (l_v_0_ > l_v_1_0 || l_v_0_ < 0)
		return false;
	return true;
}

bool ApproxCCD::pointPointIsClose(double t, double* e10_0, double* e_0, double* e_1, double tolerance_2)
{
	double distance[3];
	for (int i = 0; i < 3; ++i) {
		distance[i] = e10_0[i] + t * (e_1[i] - e_0[i]);
	}
	if (DOT(distance, distance) > tolerance_2) {
		return false;
	}
	return true;
}