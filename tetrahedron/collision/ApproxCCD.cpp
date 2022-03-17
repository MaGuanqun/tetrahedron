#include"ApproxCCD.h"
#include"../basic/global.h"
#include <iostream>
#include <iomanip>

#undef OBTAIN_CURRENT_POS
#define OBTAIN_CURRENT_POS(dest, v0,v1,e0,e1,t)\
dest[0] = v1[0] - v0[0] + t * (e1[0] - e0[0]);\
dest[1] = v1[1] - v0[1] + t * (e1[1] - e0[1]);\
dest[2] = v1[2] - v0[2] + t * (e1[2] - e0[2]);


bool ApproxCCD::VertexTriangleCollisionTime(double& t, double* initial_position, double* current_position,
	double* initial_triangle_1, double* current_triangle_1, double* initial_triangle_2, double* current_triangle_2, double* initial_triangle_3, double* current_triangle_3,
	double eta)
{
	eta = eta * eta;
	double v0[3], v1[3], v2[3], v3[3];
	SUB(v0, current_position, initial_position);
	SUB(v1, current_triangle_1, initial_triangle_1);
	SUB(v2, current_triangle_2, initial_triangle_2);
	SUB(v3, current_triangle_3, initial_triangle_3);

	double x10[3], x20[3], x30[3];
	double v10[3], v20[3], v30[3];

	std::vector<TimeInterval> coplane, e1, e2, e3;

	double temp[3];
	//here, the plane is constructed by v1v3 and normal,
	//if there is a collision, the collision point must on the positve (front) side of the plane,
	SUB(x10, initial_position, initial_triangle_1);
	SUB(v10, v0, v1);
	SUB(x30, initial_triangle_3, initial_triangle_1);
	SUB(v30, v3, v1);
	SUB(temp, initial_triangle_2, initial_triangle_1);
	CROSS(x20, x30, temp);
	SUB(temp, v2, v1);
	CROSS(v20, v30, temp);
	
	planePoly3D(x10, x20, x30, v10, v20, v30, e1);

	if (e1.empty()) {
		return false;
	}

	SUB(x10, initial_position, initial_triangle_2);
	SUB(v10, v0, v2);
	SUB(x30, initial_triangle_1, initial_triangle_2);
	SUB(v30, v1, v2);
	SUB(temp, initial_triangle_3, initial_triangle_2);
	CROSS(x20, x30, temp);
	SUB(temp, v3, v2);
	CROSS(v20, v30, temp);
	planePoly3D(x10, x20, x30, v10, v20, v30, e2);

	if (e2.empty())
		return false;

	SUB(x10, initial_position, initial_triangle_3);
	SUB(v10, v0, v3);
	SUB(x30, initial_triangle_2, initial_triangle_3);
	SUB(v30, v2, v3);
	SUB(temp, initial_triangle_1, initial_triangle_3);
	CROSS(x20, x30, temp);
	SUB(temp, v1, v3);
	CROSS(v20, v30, temp);

	planePoly3D(x10, x20, x30, v10, v20, v30, e3);

	if (e3.empty())
		return false;


}


bool ApproxCCD::TimeInterval::overlap(const TimeInterval& t1, const TimeInterval& t2)
{
	return !(t1.l > t2.u || t2.l > t1.u);
}


bool ApproxCCD::TimeInterval::overlap(const std::vector<TimeInterval>& intervals)
{
	for (std::vector<TimeInterval>::const_iterator it1 = intervals.begin(); it1 != intervals.end(); ++it1)
	{
		std::vector<TimeInterval>::const_iterator it2 = it1;
		for (++it2; it2 != intervals.end(); ++it2)
			if (!overlap(*it1, *it2))
				return false;
	}
	return true;
}

ApproxCCD::TimeInterval ApproxCCD::TimeInterval::intersect(const std::vector<TimeInterval>& intervals)
{
	TimeInterval isect(0.0, 1.0);
	for (std::vector<TimeInterval>::const_iterator it = intervals.begin(); it != intervals.end(); ++it)
	{
		isect.l = myMax(it->l, isect.l);
		isect.u = myMin(it->u, isect.u);
	}
	return isect;
}


void ApproxCCD::planePoly3D(double* x10, double* x20, double* x30, double* v10, double* v20, double* v30,
	std::vector<TimeInterval>& result)
{
	double op[4];
	double temp_2030[3], temp_v20_x30[3], temp_x20_v30[3];
	CROSS(temp_2030, v20, v30);
	CROSS(temp_v20_x30, v20, x30);
	CROSS(temp_x20_v30, x20, v30);
	op[0] = DOT(v10, temp_2030);
	op[1] = DOT(x10, temp_2030) + DOT(v10, temp_x20_v30) + DOT(v10,temp_v20_x30);
	CROSS(temp_2030, x20, x30);
	op[2] = DOT(x10, temp_x20_v30) + DOT(x10, temp_v20_x30) + DOT(v10, temp_2030);
	op[3] = DOT(x10, temp_2030);

	//std::cout << op[0] << " " << op[1] << " " << op[2] << " " << op[3] << std::endl;
	findIntervals(op, 3, result, true);
}


bool ApproxCCD::pointTriangleCollisionTime(double& t, double* initial_position, double* current_position,
	double* initial_triangle_0, double* current_triangle_0, double* initial_triangle_1, double* current_triangle_1, double* initial_triangle_2, double* current_triangle_2,
	double* initial_normal_not_normalized, double* current_normal_not_normalized, double* cross_for_CCD, double tolerance_2)//floating* f_initial_normal, floating* f_current_normal, floating* f_cross_for_CCD, 
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

	double qp_1[3];
	SUB(qp_1, current_position, current_triangle_0);
	double b1 = DOT(qp_1, current_normal_not_normalized);
	double b2 = DOT(qp_1, initial_normal_not_normalized);
	double b3 = DOT(e_0, cross_for_CCD);
	double b4 = DOT(qp_1, cross_for_CCD);
	double b5 = DOT(e_0, current_normal_not_normalized);

	double a3 = -d + b1 + b2 + b3 - b4 - b5;
	double a2 = 3 * d + b5 - 2 * b3 + b4 - 2 * b2;
	double a1 = -3 * d + b2 + b3;

	std::vector<double>time;
	time.reserve(7);
	double t0 = 2.0; double t1 = 2.0; double t2 = 2.0;
	bool is_collide = false;

	if(a3 > -NEAR_ZERO && a3<NEAR_ZERO
		&& a2 > -NEAR_ZERO && a2<NEAR_ZERO
		&& a1 > -NEAR_ZERO && a1 < NEAR_ZERO){
		if (d<NEAR_ZERO && d>-NEAR_ZERO) {
			double e_[3];
			SUB(e_, initial_triangle_1, initial_triangle_0);
			if (pointEdgeCollisionTime(t, u, u_0, u_1, e_, e_0, e_1, tolerance_2)) {
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
			//double e_0_1[3], e_1_2[3], e_2_0[3];
			//CROSS(e_0_1, e_0, e_1);
			//CROSS(e_1_2, e_1, e_2);
			//CROSS(e_2_0, e_2, e_0);
			//if (DOT(e_0_1, e_1_2) > 0 && DOT(e_1_2, e_2_0) > 0 && DOT(e_0_1, e_2_0) > 0) {
			//	t = 0;
			//	return true;
			//}		
		}
		else {
			is_collide = false;
		}
	}
	else {
		if (solveEquation(t0, t1, t2, a3, a2, a1, d)) {
			if (checkInside(t0, initial_position, initial_triangle_0, initial_triangle_1, initial_triangle_2,
				u, u_0, u_1, u_2)) {
				time.push_back(t0);
				/*if (t < 1e-11) {
					std::cout << "point-triangle-inside "<<t << std::endl;
				}*/
			}
			else {
				if (t1 < 1.5) {
					if (checkInside(t1, initial_position, initial_triangle_0, initial_triangle_1, initial_triangle_2,
						u, u_0, u_1, u_2)) {
						time.push_back(t1);
					}
					else {
						if (t2 < 1.5) {
							if (checkInside(t2, initial_position, initial_triangle_0, initial_triangle_1, initial_triangle_2,
								u, u_0, u_1, u_2)) {
								time.push_back(t2);
							}
						}
					}
				}
			}
		}
	}



	if (!time.empty()) {
		t = time[0];
		for (int i = 1; i < time.size(); ++i) {
			if (t > time[i]) {
				t = time[i];
			}
		}
		t *= 0.8;

		is_collide = true;
	}
	

	//We add an extral collision test for CCD 
	Vec3d initial_pos_a = Vec3d(initial_position[0], initial_position[1], initial_position[2]);
	Vec3d current_pos_a = Vec3d(current_position[0], current_position[1], current_position[2]);
	Vec3d initial_triangle_0_a = Vec3d(initial_triangle_0[0], initial_triangle_0[1], initial_triangle_0[2]);
	Vec3d initial_triangle_1_a = Vec3d(initial_triangle_1[0], initial_triangle_1[1], initial_triangle_1[2]);
	Vec3d initial_triangle_2_a = Vec3d(initial_triangle_2[0], initial_triangle_2[1], initial_triangle_2[2]);
	Vec3d current_triangle_0_a = Vec3d(current_triangle_0[0], current_triangle_0[1], current_triangle_0[2]);
	Vec3d current_triangle_1_a = Vec3d(current_triangle_1[0], current_triangle_1[1], current_triangle_1[2]);
	Vec3d current_triangle_2_a = Vec3d(current_triangle_2[0], current_triangle_2[1], current_triangle_2[2]);
	rootparity::RootParityCollisionTest root_parity(initial_pos_a, initial_triangle_0_a, initial_triangle_1_a, initial_triangle_2_a,
				current_pos_a, current_triangle_0_a, current_triangle_1_a, current_triangle_2_a,false);
	if (root_parity.point_triangle_collision())
	{
		if (!is_collide) {
			std::cout << "PT actually collide but does not test out " << std::endl;
			std::cout << "vt not equation coeff " << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
		}
	}	
	testInsideOutside(t, initial_position, u, initial_triangle_0, u_0, initial_triangle_1, u_1, initial_triangle_2, u_2);


	if (is_collide) {
		std::cout << "vt equation coeff " << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
	}

	return is_collide;
}




bool ApproxCCD::edgeEdgeCollisionTime(double& t, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
	double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
	double* initial_compare_edge_vertex_1, double tolerance_2)
{
	double e0[3], e0_[3], e1[3], e1_[3], e2[3], e2_[3];
	//double n0[3];
	//floating f_e0[3], f_e0_[3], f_e1[3], f_e1_[3], f_e2[3], f_e2_[3];
	//floating f_current_edge_vertex_0[3], f_current_edge_vertex_1[3], f_initial_edge_vertex_0[3], f_initial_edge_vertex_1[3],
	//	f_current_compare_edge_vertex_0[3], f_current_compare_edge_vertex_1[3], f_initial_compare_edge_vertex_0[3],
	//	f_initial_compare_edge_vertex_1[3];
	//make_vector(current_edge_vertex_0, f_current_edge_vertex_0);
	//make_vector(current_edge_vertex_1, f_current_edge_vertex_1);
	//make_vector(initial_edge_vertex_0, f_initial_edge_vertex_0);
	//make_vector(initial_edge_vertex_1, f_initial_edge_vertex_1);
	//make_vector(initial_compare_edge_vertex_0, f_initial_compare_edge_vertex_0);
	//make_vector(initial_compare_edge_vertex_1, f_initial_compare_edge_vertex_1);
	//make_vector(current_compare_edge_vertex_0, f_current_compare_edge_vertex_0);
	//make_vector(current_compare_edge_vertex_1, f_current_compare_edge_vertex_1);


	SUB(e0, initial_compare_edge_vertex_1, initial_edge_vertex_0);//d-a
	SUB(e0_, current_compare_edge_vertex_1, current_edge_vertex_0);
	SUB(e1, initial_edge_vertex_1, initial_edge_vertex_0);//b-a
	SUB(e1_, current_edge_vertex_1, current_edge_vertex_0);
	SUB(e2, initial_compare_edge_vertex_0, initial_edge_vertex_0);//c-a
	SUB(e2_, current_compare_edge_vertex_0, current_edge_vertex_0);

	double a[3];
	CROSS(a, e1, e2);
	double d = DOT(e0, a);

	double b[3], c[3];
	CROSS(c, e1, e2_);
	CROSS(b, e1_, e2);
	SUM_(b, c);
	CROSS(c, e1_, e2_);
	double b1 = DOT(e0_, c);
	double b2 = DOT(e0_, a);
	double b3 = DOT(e0, b);
	double b4 = DOT(e0_, b);
	double b5 = DOT(e0, c);


	double a3 = -d + b1 + b2 + b3 - b4 - b5;
	double a2 = 3 * d + b5 - 2 * b3 + b4 - 2 * b2;
	double a1 = -3 * d + b2 + b3;

	std::vector<double>time;
	time.reserve(7);

	double d_c[3];
	SUB(d_c, initial_compare_edge_vertex_1, initial_compare_edge_vertex_0);//d-c

	double u_0_0[3], u_0_1[3], u_1_0[3], u_1_1[3];
	SUB(u_0_0, current_edge_vertex_0, initial_edge_vertex_0);
	SUB(u_0_1, current_edge_vertex_1, initial_edge_vertex_1);
	SUB(u_1_0, current_compare_edge_vertex_0, initial_compare_edge_vertex_0);
	SUB(u_1_1, current_compare_edge_vertex_1, initial_compare_edge_vertex_1);

	bool is_collide = false;


	if (a3 > -NEAR_ZERO && a3<NEAR_ZERO
		&& a2 > -NEAR_ZERO && a2<NEAR_ZERO
		&& a1 > -NEAR_ZERO && a1 < NEAR_ZERO) {
		if (d<NEAR_ZERO && d>-NEAR_ZERO) { //the twos edge always coplanar
			double v_0[3], v_1[3];
			//initial_edge_vertex_0
			SUB(v_0, initial_edge_vertex_0, initial_compare_edge_vertex_0);
			SUB(v_1, initial_edge_vertex_0, initial_compare_edge_vertex_1);
			if (pointEdgeCollisionTime(t, u_0_0, u_1_0, u_1_1, d_c, v_0, v_1, tolerance_2)) {
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
			if (pointEdgeCollisionTime(t, u_0_1, u_1_0, u_1_1, d_c, v_0, v_1, tolerance_2)) {
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
			////check when t=0,the two edges are coplanar if the two edges are collide
			//double e_temp0[3], compare_direction0[3], compare_direction1[3];
			//SUB(e_temp0, initial_compare_edge_vertex_0, initial_edge_vertex_1);//c-b		
			//CROSS(compare_direction0, d_c, e_temp0);
			//CROSS(compare_direction1, d_c, e2);
			//double d2 = DOT(compare_direction0, compare_direction1);
			////SUB(e_temp0, initial_compare_edge_vertex_1, initial_edge_vertex_0);//d-a
			//CROSS(compare_direction0, e1, e0);
			//CROSS(compare_direction1, e1, e2);
			//double d1 = DOT(compare_direction0, compare_direction1);		
			//if ( d1 < 0 && d2<0) { //this means the vertices of one edge is on different sides of the other edge
			//	t = 0;
			//	return true;
			//}
			//// initial_compare_edge_vertex_0  initial_edge
			//if (pointEdgeIsClose(e2, e1, tolerance_2)) {
			//	t = 0;
			//	return true;
			//}
			//// initial_compare_edge_vertex_1  initial_edge
			//if (pointEdgeIsClose(e0, e1, tolerance_2)) {
			//	t = 0;
			//	return true;
			//}
			////initial_edge_vertex_0 initial_compare_edge
			//SUB(e_temp0, initial_edge_vertex_0, initial_compare_edge_vertex_0);//a-c
			//if (pointEdgeIsClose(e_temp0, d_c, tolerance_2)) {
			//	t = 0;
			//	return true;
			//}
			//SUB(e_temp0, initial_edge_vertex_1, initial_compare_edge_vertex_0);//b-c
			//if (pointEdgeIsClose(e_temp0, d_c, tolerance_2)) {
			//	t = 0;
			//	return true;
			//}
		}
		else {
			is_collide = false;
		}

	}
	else {
		double t0 = 2.0; double t1 = 2.0; double t2 = 2.0;
		if (solveEquation(t0,t1,t2, a3, a2, a1, d)) {
			//if(estimate){
			//	if (tight_CCD.insideTest(f_initial_edge_vertex_0, f_initial_edge_vertex_1, f_initial_compare_edge_vertex_0,
			//		f_initial_compare_edge_vertex_1, f_current_edge_vertex_0, f_current_edge_vertex_1, f_current_compare_edge_vertex_0,
			//		f_current_compare_edge_vertex_1, f_n0, f_n1, f_cross_for_CCD, true)) {
			//		time.push_back(t);
			//		/*if (t < 1e-11) {
			//			std::cout << "edge-edge approx " << t << std::endl;
			//			std::cout << std::setprecision(20) << initial_edge_vertex_0[0] << ", " << initial_edge_vertex_0[1] << ", " << initial_edge_vertex_0[2] << std::endl;
			//			std::cout << std::setprecision(20) << initial_edge_vertex_1[0] << ", " << initial_edge_vertex_1[1] << ", " << initial_edge_vertex_1[2] << std::endl;
			//			std::cout << std::setprecision(20) << initial_compare_edge_vertex_0[0] << ", " << initial_compare_edge_vertex_0[1] << ", " << initial_compare_edge_vertex_0[2] << std::endl;
			//			std::cout << std::setprecision(20) << initial_compare_edge_vertex_1[0] << ", " << initial_compare_edge_vertex_1[1] << ", " << initial_compare_edge_vertex_1[2] << std::endl;
			//			std::cout << std::setprecision(20) << current_edge_vertex_0[0] << ", " << current_edge_vertex_0[1] << ", " << current_edge_vertex_0[2] << std::endl;
			//			std::cout << std::setprecision(20) << current_edge_vertex_1[0] << ", " << current_edge_vertex_1[1] << ", " << current_edge_vertex_1[2] << std::endl;
			//			std::cout << std::setprecision(20) << current_compare_edge_vertex_0[0] << ", " << current_compare_edge_vertex_0[1] << ", " << current_compare_edge_vertex_0[2] << std::endl;
			//			std::cout << std::setprecision(20) << current_compare_edge_vertex_1[0] << ", " << current_compare_edge_vertex_1[1] << ", " << current_compare_edge_vertex_1[2] << std::endl;
			//		}*/
			//	}
			//}
			//else {
			if (edgeEdgeCheckInside(t0, initial_edge_vertex_0, initial_edge_vertex_1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1,
				u_0_0, u_0_1, u_1_0, u_1_1)) {
				time.push_back(t0);
				/*if (t < 1e-11) {
					std::cout << "edge-edge inside " << t << std::endl;
				}*/
			}
			else {
				if (t1 < 1.5) {
					if (edgeEdgeCheckInside(t1, initial_edge_vertex_0, initial_edge_vertex_1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1,
						u_0_0, u_0_1, u_1_0, u_1_1)) {
						time.push_back(t1);
					}
					else {
						if (t2 < 1.5) {
							if (edgeEdgeCheckInside(t2, initial_edge_vertex_0, initial_edge_vertex_1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1,
								u_0_0, u_0_1, u_1_0, u_1_1)) {
								time.push_back(t2);
							}
						}
					}
				}
			}
			//}
		}
	}

	

	if (!time.empty()) {
		t = time[0];
		for (int i = 1; i < time.size(); ++i) {
			if (t > time[i]) {
				t = time[i];
			}
		}
		t *= 0.8;
		is_collide = true;
	}

	Vec3d initial_edge_0 = Vec3d(initial_edge_vertex_0[0], initial_edge_vertex_0[1], initial_edge_vertex_0[2]);
	Vec3d current_edge_0 = Vec3d(current_edge_vertex_0[0], current_edge_vertex_0[1], current_edge_vertex_0[2]);
	Vec3d initial_edge_1 = Vec3d(initial_edge_vertex_1[0], initial_edge_vertex_1[1], initial_edge_vertex_1[2]);
	Vec3d current_edge_1 = Vec3d(current_edge_vertex_1[0], current_edge_vertex_1[1], current_edge_vertex_1[2]);
	Vec3d initial_comapre_edge_0 = Vec3d(initial_compare_edge_vertex_0[0], initial_compare_edge_vertex_0[1], initial_compare_edge_vertex_0[2]);
	Vec3d current_comapre_edge_0 = Vec3d(current_compare_edge_vertex_0[0], current_compare_edge_vertex_0[1], current_compare_edge_vertex_0[2]);
	Vec3d initial_comapre_edge_1 = Vec3d(initial_compare_edge_vertex_1[0], initial_compare_edge_vertex_1[1], initial_compare_edge_vertex_1[2]);
	Vec3d current_comapre_edge_1 = Vec3d(current_compare_edge_vertex_1[0], current_compare_edge_vertex_1[1], current_compare_edge_vertex_1[2]);
	rootparity::RootParityCollisionTest root_parity(initial_edge_0, initial_edge_1, initial_comapre_edge_0, initial_comapre_edge_1,
		current_edge_0, current_edge_1, current_comapre_edge_0, current_comapre_edge_1, true);
	if (root_parity.edge_edge_collision())
	{
		if (!is_collide) {
			std::cout << "EE actually collide but does not test out " << std::endl;
			std::cout << "ee not equation coeff " << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
		}
	}

	testEdgeInsideOutside(t, initial_edge_vertex_0, u_0_0, initial_edge_vertex_1, u_0_1, initial_compare_edge_vertex_0, u_1_0, initial_compare_edge_vertex_1, u_1_1);

	if (is_collide) {
		std::cout << "ee equation coeff " << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
	}

	return is_collide;
}


bool ApproxCCD::checkInside(double t, double* v0, double* v1, double* v2, double* v3,
	double* e0, double* e1, double* e2, double* e3)
{
	double current_v0[3], current_v1[3], current_v2[3];

	double e0_1[3], e0_2[3], e0_3[3];
	//for (unsigned int i = 0; i < 3; ++i) {
	//	e0_1[i] = v1[i] - v0[i] + t * (e1[i] - e0[i]);
	//	e0_2[i] = v2[i] - v0[i] + t * (e2[i] - e0[i]);
	//	e0_3[i] = v3[i] - v0[i] + t * (e3[i] - e0[i]);
	//}

	OBTAIN_CURRENT_POS(e0_1, v0, v1, e0, e1, t);
	OBTAIN_CURRENT_POS(e0_2, v0, v2, e0, e2, t);
	OBTAIN_CURRENT_POS(e0_3, v0, v3, e0, e3, t);

	CROSS(current_v0, e0_1, e0_2);//here we use current_vi to store the norm to save storage
	CROSS(current_v1, e0_2, e0_3);

	if (DOT(current_v0, current_v1) < 0) {
		return false;
	}

	CROSS(current_v2, e0_3, e0_1);
	
	if (DOT(current_v1, current_v2) < 0) {
		return false;
	}
	if (DOT(current_v2, current_v0) < 0) {
		return false;
	}
	return true;
}



bool ApproxCCD::edgeEdgeCheckInside(double t, double* v0, double* v1, double* v2, double* v3,
	double* e0, double* e1, double* e2, double* e3)
{
	//double current_v0[3], current_v1[3], current_v2[3], current_v3[3];
	//for (int i = 0; i < 3; ++i) {
	//	current_v0[i] = v0[i] + t * e0[i];
	//	current_v1[i] = v1[i] + t * e1[i];
	//	current_v2[i] = v2[i] + t * e2[i];
	//	current_v3[i] = v3[i] + t * e3[i];
	//}
	double e_0_2[3], e_1_2[3], e_0_3[3], e_1_3[3], direct0[3], direct1[3];

	OBTAIN_CURRENT_POS(e_0_2, v2, v0, e2, e0, t);
	OBTAIN_CURRENT_POS(e_1_2, v2, v1, e2, e1, t);
	OBTAIN_CURRENT_POS(e_0_3, v3, v0, e3, e0, t);
	OBTAIN_CURRENT_POS(e_1_3, v3, v1, e3, e1, t);
	//SUB(e_0, current_v0, current_v2);
	//SUB(e_1, current_v1, current_v2);
	CROSS(direct0, e_0_2, e_1_2);
	//SUB(e_0, current_v0, current_v3);
	//SUB(e_1, current_v1, current_v3);
	CROSS(direct1, e_0_3, e_1_3);

	if (DOT(direct0, direct1) > 0) {
		return false;
	}
	//SUB(e_0, current_v2, current_v0);
	//SUB(e_1, current_v3, current_v0);
	CROSS(direct0, e_0_2, e_0_3);
	//SUB(e_0, current_v2, current_v1);
	//SUB(e_1, current_v3, current_v1);
	CROSS(direct1, e_1_2, e_1_3);
	if (DOT(direct0, direct1) > 0) {
		return false;
	}
	return true;
}


bool ApproxCCD::solveCubicEquation(double a, double b, double c, double d, double& t0, double& t1, double& t2)
{
	b = b / a;
	c = c / a;
	d = d / a;

	double alpha = (-2.0 * b * b * b + 9.0 * b * c - 27.0 * d) / 54.0;
	double beta = (b * b - 3.0 * c) / 9.0;
	double dum1 = beta * beta * beta;
	double delta = alpha * alpha - dum1;
	b /= 3.0;

	std::cout << "delta " << delta<< std::endl;

	if (delta > 0) {
		delta = sqrt(delta);
		//double R1 = ;
		//double R2 = ;
		t0 = -b + cbrt(alpha + delta) + cbrt(alpha - delta);
		std::cout << "t " << t0 << std::endl;
		if (t0 >= 0.0 && t0 <= 1.0) {
			return true;
		}
		else {
			t0 = 2.0;
		}
	}
	else if (delta == 0.0) {	//all roots are real, at least two are equal.
		alpha = cbrt(alpha);
		if ((-b + alpha + alpha) >= 0.0 && (-b + alpha + alpha) <= 1.0)
		{
			t0 = -b + alpha + alpha;
		}
		if ((alpha + b) <= 0.0 && (alpha + b) >= -1.0)
		{
			t1 = -alpha - b;
		}
		if (t1 < t0) {
			d = t1;
			t1 = t0;
			t0 = d;
		}
		if (t0 < 1.5) {
			return true;
		}
	}
	else { //all roots are real and different
		dum1 = acos(alpha / sqrt(dum1));
		alpha = 2.0 * std::sqrt(beta);
		t0 = -b + alpha * cos(dum1 / 3.0);
		t1 = -b + alpha * cos((dum1 + 2.0 * M_PI) / 3.0);
		t2 = -b + alpha * cos((dum1 + 4.0 * M_PI) / 3.0);

		if (t0 < 0.0 || t0>1.0) {
			t0 = 2.0;
		}
		if (t1 < 0.0 || t1>1.0) {
			t1 = 2.0;
		}
		if (t2 < 0.0 || t2>1.0) {
			t2 = 2.0;
		}
		sortABC(t0, t1, t2);
		if (t0 < 1.5) {
			return true;
		}
	}
	return false;

}


inline void ApproxCCD::sortABC(double& a, double& b, double& c)
{
	double d;
	if (a > b) { d = a; a = b; b = d; }
	if (a > c) { d = a; a = c; c = d; }
	if (b > c) { d = b; b = c; c = d; }
}

//t
bool ApproxCCD::solveEquation(double& t0, double& t1, double& t2, double a3, double a2, double a1, double d)
{
	//std::cout << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
	if (d<NEAR_ZERO && d>-NEAR_ZERO) {
		return solveQuadraticEquation(t0,t1, a3, a2, a1);
	}
	if (a3<NEAR_ZERO && a3>-NEAR_ZERO) {
		return solveQuadraticEquation(t0,t1, a2, a1, d);
	}
	return solveCubicEquation(a3, a2, a1, d, t0, t1, t2);
}

//bool ApproxCCD::solveEquation(double&t, double a3, double a2, double a1, double d)
//{
//	//std::cout << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
//	if (d<NEAR_ZERO2 && d>-NEAR_ZERO2) {
//		return solveQuadraticEquation(t, a3, a2, a1);
//	}	
//	if (a3<NEAR_ZERO2 && a3>-NEAR_ZERO2) {
//		return solveQuadraticEquation(t, a2, a1, d);
//	}
//	double a_2 = 3 * a3;
//	double a_1 = 2 * a2;
//	double check_root = a_1 * a_1 - 4 * a_2 * a1;
//	if (check_root > 0) {
//		check_root = sqrt(check_root);
//		if (d < 0) {
//			a3 = -a3;
//			a2 = -a2;
//			a1 = -a1;
//			d = -d;
//			a_2 = -a_2;
//			a_1 = -a_1;
//		}
//		double i_0, i_1;
//		if (a3 > 0) {
//			i_0 = (-a_1 - check_root) / (2 * a_2);
//			i_1 = (-a_1 + check_root) / (2 * a_2);
//			//std::cout << "i_0 " << i_0 << " " << i_1 << std::endl;
//			if (i_1 < 0) {
//				return false;
//			}
//			if (i_0 > 1) {
//				return false;
//			}
//			if (i_0 < 0) {
//				if (i_1 > 1) {//interval [0, 1]
//					if (d * (a3 + a2 + a1 + d) < 0) {
//						double t_s = -a2 / a_2;
//						if (t_s < 0) {
//							t = -d / a1;
//						}
//						else {
//							if (t_s < 1) {
//								t = -d / (a_2 * t_s * t_s + a_1 * t_s + a1);
//							}
//							else {
//								t = -d / (a_2 + a_1 + a1);
//							}
//						}						
//						return true;
//					}						
//				}
//				else {//interval [0, i_1]
//					if (d * (a3 * i_1 * i_1 * i_1 + a2 * i_1 * i_1 + a1 * i_1 + d) < 0) {
//						double t_s = -a2 / a_2;
//						if (t_s < 0) {
//							t = -d / a1;
//						}
//						else {
//							t = -d / (a_2 * t_s * t_s + a_1 * t_s + a1);
//						}						
//						return true;
//					}
//				}
//				return false;
//			}		
//			double r_0 = a3 * i_0 * i_0 * i_0 + a2 * i_0 * i_0 + a1 * i_0 + d;
//			if (i_1 > 1) {//interval [i_0, 1]
//				if (r_0 * (a3 + a2 + a1 + d) < 0) {
//					double t_s = -a2 / a_2;
//					if (t_s > 1) {
//						t = -d / (a_2 + a_1 + a1);
//					}
//					else {
//						t = i_0 - r_0 / (a_2 * t_s * t_s + a_1 * t_s + a1);
//					}					
//					return true;
//				}						
//			}
//			else {//interval [i_0, i_1]
//				if (r_0 * (a3 * i_1 * i_1 * i_1 + a2 * i_1 * i_1 + a1 * i_1 + d) < 0) {
//					double t_s = -a2 / a_2;
//					t = i_0 - r_0 / (a_2 * t_s * t_s + a_1 * t_s + a1);
//					return true;
//				}
//			}
//			return false;			
//		}
//		else {
//			i_0 = (-a_1 + check_root) / (2 * a_2);
//			i_1 = (-a_1 - check_root) / (2 * a_2);
//			if (i_0 > 0) {
//				if (i_0 > 1) {//interval [0, 1]
//					if (d * (a3 + a2 + a1 + d) < 0) {
//						t = -d / a1;
//						return true;
//					}
//				}
//				else {//interval [0, i_0]
//					if (d * (a3 * i_0 * i_0 * i_0 + a2 * i_0 * i_0 + a1 * i_0 + d) < 0) {
//						t = -d / a1;
//						return true;
//					}
//				}
//				return false;
//			}
//			if (i_1 > 1) {
//				return false;
//			}
//			if (i_1 < 0) { //interval [0,1]
//				if (d * (a3 + a2 + a1 + d) < 0) {
//					t = -d / (a_2 + a_1 + a1);
//					return true;
//				}
//			}
//			else { //interval [i_1,1]
//				double r_i_1 = a3 * i_1 * i_1 * i_1 + a2 * i_1 * i_1 + a1 * i_1 + d;
//				if (r_i_1 * (a3 + a2 + a1 + d) < 0) {
//					t = i_1 - r_i_1 / (a_2 + a_1 + a1);
//					return true;
//				}
//			}
//			return false;		
//		}
//	}
//	else {
//		if (d * (a3 + a2 + a1 + d) < 0) {
//			double r_1 = a_2 + a_1 + a1;
//			if (fabs(a1) > fabs(r_1)) {
//				t = -d / a1;
//			}
//			else {
//				t = -d / r_1;
//			}
//			return true;
//		}
//		return false;
//	}	
//}


bool ApproxCCD::pointPointCollisionTime(double& t, double* e_1_0, double* e_0, double* e_1, double tolerance_2)
{
	std::vector<double> t_list;
	t_list.reserve(3);
	double t_temp;
	//x
	t_temp = e_1_0[0] / (e_0[0] - e_1[0]);
	if (t_temp > 0 && t_temp <= 1.0) {
		t_list.emplace_back(t_temp);
	}
	//y
	t_temp = e_1_0[1] / (e_0[1] - e_1[1]);
	if (t_temp > 0 && t_temp <= 1.0) {
		t_list.emplace_back(t_temp);
	}
	//z
	t_temp = e_1_0[2] / (e_0[2] - e_1[2]);
	if (t_temp > 0 && t_temp <= 1.0) {
		t_list.emplace_back(t_temp);
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


int ApproxCCD::getQuadRoots(double a, double b, double c, double& t0, double& t1)
{
	int roots = 0;
	int sign = 1;

	if (b < 0) {
		sign = -1;
	}
	double D= b * b - 4.0 * a * c;
	if (D >= 0) {
		roots = 2;
		double q = -0.5 * (b + sign * sqrt(D)); //here we try to use sum for two positive values, avoid use subtract
		t0 = q / a;
		t1 = c / q;
		if (t0 > t1) {
			std::swap(t0, t1);
		}
	}
	return roots;
}

bool ApproxCCD::solveQuadraticEquation(double& t0, double& t1, double a2, double a1, double a0)
{
	if (a0<NEAR_ZERO2 && a0>-NEAR_ZERO2) {
		if (a2<NEAR_ZERO2 && a2>-NEAR_ZERO2) {
			return false;
			//std::cout << "all zero "<<a2<<" "<<a1<<" "<<a0 << std::endl;
		}
		t0 = -a1 / a2;
		if (t0 > 0.0 && t0 <= 1.0) {
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
		t0 = -a0 / a1;
		if (t0 > 0.0 && t0 <= 1.0) {
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
			t0 = (-a1 - temp) / (a2 + a2);
			if (t0 > 1.0) {
				return false;
			}
			else if (t0 < 0.0) {
				t0 = 2.0;
			}
			t1 = (-a1 + temp) / (a2 + a2);
			if (t1 < 0.0) {
				return false;
			}
			if (t1 > 1.0) {
				t1 = 2.0;
			}
			if (t1 < t0) {
				temp = t1;
				t1 = t0;
				t0 = temp;
			}
			if (t0 <= 1.0 && t0 >= 0.0) {
				return true;
			}
			return false;
		}
	}
}


//log:
//Here, we first test vertex-edge and vertex-vertex. 
//



void ApproxCCD::testInsideOutside(double t, double* initial_pos, double* u, double* initial_triangle_0,
	double* u_0, double* initial_triangle_1, double* u_1, double* initial_triangle_2, double* u_2)
{
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
	if (DOT(e_initial, norm_initial) * DOT(e_result, norm_result) < 0) {
		std::cout << "vt side indicate " << DOT(e_initial, norm_initial) << " " << DOT(e_result, norm_result) << " " << DOT(norm_result, norm_result) << std::endl;
		std::cout << "vt side has changed " << std::endl;
	}
}

void ApproxCCD::testEdgeInsideOutside(double t, double* initial_pos_0_0, double* u_0_0, double* initial_pos_0_1,
	double* u_0_1, double* initial_triangle_1_0, double* u_1_0, double* initial_triangle_1_1, double* u_1_1)
{
	double result_v[3];
	double result_p0[3];
	double result_p1[3];
	double result_p2[3];
	POINT_ON_EDGE(result_v, 1.0, t, initial_pos_0_0, u_0_0);
	POINT_ON_EDGE(result_p0, 1.0, t, initial_pos_0_1, u_0_1);
	POINT_ON_EDGE(result_p1, 1.0, t, initial_triangle_1_0, u_1_0);
	POINT_ON_EDGE(result_p2, 1.0, t, initial_triangle_1_1, u_1_1);
	double e_result1[3], e_result2[3];
	SUB(e_result1, result_v, result_p0);
	SUB(e_result2, result_p1, result_p2);
	double norm_result[3];
	CROSS(norm_result, e_result1, e_result2);
	SUB(e_result1, initial_pos_0_0, initial_pos_0_1);
	SUB(e_result2, initial_triangle_1_0, initial_triangle_1_1);
	double norm_initial[3];
	CROSS(norm_initial, e_result1, e_result2);
	double e_result[3];
	SUB(e_result, result_v, result_p1);
	double e_initial[3];
	SUB(e_initial, initial_pos_0_0, initial_triangle_1_0);	
	if (DOT(e_initial, norm_initial) * DOT(e_result, norm_result) < 0) {
		std::cout << "ee side indicate " << DOT(e_initial, norm_initial) << " " << DOT(e_result, norm_result) << " " << DOT(norm_result, norm_result) << std::endl;
		std::cout << "ee side has changed " << std::endl;
	}
}


//edge-edge
//void ApproxCCD::test()
//{
//	//double a = (double)rand() / RAND_MAX * 10.0;
//	//double b = (double)rand() / RAND_MAX;
//	//double c = (double)rand() / RAND_MAX;
//	//double d = (double)rand() / RAND_MAX * 0.02;
//	////double a = 10.0;
//	////double b = 35.0;
//	////double c = 25.0;
//	////double d = 4.0;
//	//double x1, x2, x3;
//	//double x11, x12, x13;
//	//for (unsigned int i = 0; i < 10; ++i) {
//	//	a = (double)rand() / RAND_MAX * 10.0;
//	//	b = (double)rand() / RAND_MAX;
//	//	c = (double)rand() / RAND_MAX;
//	//	d = (double)rand() / RAND_MAX * 0.02;
//	//	std::cout << "a: " << a << ", b: " << b << ", c: " << c << ", d: " << d << std::endl;
//	//	solveCubicEquation(a, b, c, d, x1, x2, x3);
//	//	std::cout << "x1: " << x1 << ", x2: " << x2 << ", x3: " << x3 << std::endl;
//	//	cubicsolve(a, b, c, d, x11, x12, x13);
//	//	std::cout << "x11: " << x11 << ", x12: " << x12 << ", x13: " << x13 << std::endl;
//	//	std::cout << "===" << std::endl;
//	//}
//	//double xx1 = a * x1 * x1 * x1 + b * x1 * x1 + c * x1 + d;
//	//double xx2 = a * x2 * x2 * x2 + b * x2 * x2 + c * x2 + d;
//	//double xx3 = a * x3 * x3 * x3 + b * x3 * x3 + c * x3 + d;
//	//std::cout << "error1: " << std::abs(xx1) << ", error2: " << std::abs(xx2) << ", error3: " << std::abs(xx3) << std::endl;
//	//tight_CCD.testInside();
//	////edge 1
//	//double initial_pos[3] = { 1.0,1.0,1.0 };
//	//double initial_triangle_0[3] = { -1.0,0.0,-1.0 };
//	//
//	////edge 2
//	//double initial_triangle_1[3] = { -1.0,0.0,1.0 };
//	//double initial_triangle_2[3] = { 1.0,0.0,-1.0 };
//	//double current_pos[3] = { 1.0,0.3,1.0 };
//	//double current_triangle_0[3] = { -1.0,-0.1,-1.0 };
//	//double current_triangle_1[3] = { -0.5,0.5,1.0 };
//	//double current_triangle_2[3] = { 0.5,0.7,-1.0 };
//	//edge 1
//	double initial_pos[3] = { 1.0,1.0,1.0 };
//	double initial_triangle_0[3] = { -1.0,0.0,-1.0 };
//	//edge 2
//	double initial_triangle_1[3] = { -1,-0.5,0.3 };
//	double initial_triangle_2[3] = { 0.0,-0.5,-0.8 };
//	double current_pos[3] = { -0.5,0.3,1.0 };
//	double current_triangle_0[3] = { -1.0,-0.1,-0.5 };
//	double current_triangle_1[3] = { -0.5,0.5,0.0 };
//	double current_triangle_2[3] = { 0.5,0.5,-1.0 };
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
//	if (edgeEdgeCollisionTime(t, current_pos, current_triangle_0, initial_pos, initial_triangle_0,
//		current_triangle_1, current_triangle_2, initial_triangle_1, initial_triangle_2, tolerance_2)) {
//		std::cout << t << std::endl;
//		t *= 0.9;
//		double result_v[3];
//		double result_p0[3];
//		double result_p1[3];
//		double result_p2[3];
//		POINT_ON_EDGE(result_v, 1.0, t, initial_pos, u);
//		POINT_ON_EDGE(result_p0, 1.0, t, initial_triangle_0, u_0);
//		POINT_ON_EDGE(result_p1, 1.0, t, initial_triangle_1, u_1);
//		POINT_ON_EDGE(result_p2, 1.0, t, initial_triangle_2, u_2);
//		double e_result1[3], e_result2[3];
//		SUB(e_result1, result_v, result_p0);
//		SUB(e_result2, result_p1, result_p2);
//		double norm_result[3];
//		CROSS(norm_result, e_result1, e_result2);
//		SUB(e_result1, initial_pos, initial_triangle_0);
//		SUB(e_result2, initial_triangle_1, initial_triangle_2);
//		double norm_initial[3];
//		CROSS(norm_initial, e_result1, e_result2);
//		double e_result[3];
//		SUB(e_result, result_v, result_p1);
//		double e_initial[3];
//		SUB(e_initial, initial_pos, initial_triangle_1);
//		std::cout << "side indicate " << DOT(e_initial, norm_initial) << " " << DOT(e_result, norm_result)<<" "<< DOT(norm_result, norm_result) << std::endl;
//		if (DOT(e_initial, norm_initial) * DOT(e_result, norm_result) < 0) {					
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
//	rootparity::RootParityCollisionTest root_parity(initial_pos_a, initial_triangle_0_a, initial_triangle_1_a, initial_triangle_2_a,
//		current_pos_a, current_triangle_0_a, current_triangle_1_a, current_triangle_2_a, true);
//	if (root_parity.edge_edge_collision())
//	{
//		std::cout << "actually collide " << std::endl;
//	}
//}

//vertex-triangle
void ApproxCCD::test()
{
	double initial_pos[3] = { 0.5,0.0,0.3 };
	double current_pos[3] = { 0.5,0.0,-0.3 };
	double initial_triangle_0[3] = { 1.0,0.0,0.0 };
	double initial_triangle_1[3] = { 0.0,1.0,0.0 };
	double initial_triangle_2[3] = { 0.0,-1.0,0.0 };
	double current_triangle_0[3] = { 1.0,0.0,0.0 };
	double current_triangle_1[3] = { 0.0,1.0,0.0 };
	double current_triangle_2[3] = { 0.0,-1.0,0.0 };
	//double current_triangle_1[3] = { -1.0,0.1,-0.8 };
	//double current_triangle_2[3] = { -1.0,-0.1,0.8 };
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

	VertexTriangleCollisionTime(t, initial_pos, current_pos, initial_triangle_0, current_triangle_0,
		initial_triangle_1, current_triangle_1, initial_triangle_2, current_triangle_2, tolerance_2);

	//if (pointTriangleCollisionTime(t, initial_pos, current_pos, initial_triangle_0, current_triangle_0, 
	//	initial_triangle_1, current_triangle_1, initial_triangle_2, current_triangle_2,
	//	initial_tri_normal, tri_normal, cross_for_CCD,tolerance_2)) {
	//	std::cout << t << std::endl;
	//	t *= 0.9;
	//	double result_v[3];
	//	double result_p0[3];
	//	double result_p1[3];
	//	double result_p2[3];

	//	POINT_ON_EDGE(result_v, 1.0, t, initial_pos, u);
	//	POINT_ON_EDGE(result_p0, 1.0, t, initial_triangle_0, u_0);
	//	POINT_ON_EDGE(result_p1, 1.0, t, initial_triangle_1, u_1);
	//	POINT_ON_EDGE(result_p2, 1.0, t, initial_triangle_2, u_2);

	//	double e_result1[3], e_result2[3];
	//	SUB(e_result1, result_p1, result_p0);
	//	SUB(e_result2, result_p2, result_p0);
	//	double norm_result[3];
	//	CROSS(norm_result, e_result1, e_result2);

	//	double e_result[3];
	//	SUB(e_result, result_v, result_p0);
	//	std::cout <<"side indicate "<< DOT(e_0, initial_tri_normal) << " " << DOT(e_result, norm_result) << std::endl;
	//	if (DOT(e_0, initial_tri_normal) * DOT(e_result, norm_result) < 0) {			
	//		std::cout << "side has changed " << std::endl;
	//	}
	//}
	//else {
	//	std::cout << "does not collide " << std::endl;
	//}
	//Vec3d initial_pos_a = Vec3d(initial_pos[0], initial_pos[1], initial_pos_a[2]);
	//Vec3d current_pos_a = Vec3d(current_pos[0], current_pos[1], current_pos[2]);
	//Vec3d initial_triangle_0_a = Vec3d(initial_triangle_0[0], initial_triangle_0[1], initial_triangle_0[2]);
	//Vec3d initial_triangle_1_a = Vec3d(initial_triangle_1[0], initial_triangle_1[1], initial_triangle_1[2]);
	//Vec3d initial_triangle_2_a = Vec3d(initial_triangle_2[0], initial_triangle_2[1], initial_triangle_2[2]);
	//Vec3d current_triangle_0_a = Vec3d(current_triangle_0[0], current_triangle_0[1], current_triangle_0[2]);
	//Vec3d current_triangle_1_a = Vec3d(current_triangle_1[0], current_triangle_1[1], current_triangle_1[2]);
	//Vec3d current_triangle_2_a = Vec3d(current_triangle_2[0], current_triangle_2[1], current_triangle_2[2]);

	//rootparity::RootParityCollisionTest root_parity(initial_pos_a, initial_triangle_0_a, initial_triangle_1_a, initial_triangle_2_a,
	//	current_pos_a, current_triangle_0_a, current_triangle_1_a, current_triangle_2_a,false);
	//if (root_parity.point_triangle_collision())
	//{
	//	std::cout << "actually collide " << std::endl;
	//}
}

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

inline void ApproxCCD::make_vector(double* v, floating* out)
{
	for (int i = 0; i < 3; i++) {
		out[i] = floating(v[i], 0);
	}
}

void ApproxCCD::cubicsolve(const double& a, const double& b, const double& c, const double& d, double& x1, double& x2, double& x3)
{
	if (a < 0.000001) { // end if a == 0 
		assert("The coefficient of the cube of x is 0. Please use the utility for a SECOND degree quadratic. No further action taken.");
		return;
	}
	if (d < 0.000001) { // End if d == 0, d is definitely a root of the function
		assert("One root is 0. Now divide through by x and use the utility for a SECOND degree quadratic to solve the resulting equation for the other two roots. No further action taken.");
		return;
	}
	double bb = b / a;
	double cc = c / a;
	double dd = d / a;
	double disc, q, r, dum1, s, t, term1, r13;
	q = (3.0 * cc - (bb * bb)) / 9.0;
	r = -(27.0 * dd) + bb * (9.0 * cc - 2.0 * (bb * bb));
	r /= 54.0;
	std::cout << "q: " << q << ", r: " << r << std::endl;
	disc = q * q * q + r * r;
	x1 = 0.0;			// The first root is always real.
	term1 = (bb / 3.0);
	if (disc > 0)       // one root real, two are complex
	{
		std::cout << "one root real, two are complex" << std::endl;
		s = r + std::sqrt(disc);
		s = ((s < 0) ? -std::pow(-s, (1.0 / 3.0)) : std::pow(s, (1.0 / 3.0)));
		t = r - std::sqrt(disc);
		t = ((t < 0) ? -std::pow(-t, (1.0 / 3.0)) : std::pow(t, (1.0 / 3.0)));
		x1 = -term1 + s + t;
		term1 += (s + t) / 2.0;
		x3 = x2 = -term1;
		term1 = std::sqrt(3.0) * (-t + s) / 2;
		x2 = term1;
		x3 = -term1;
		return;
	}
	// End if (disc > 0)
	// The remaining options are all real
	x3 = x2 = 0.0;
	if (disc == 0) // All roots real, at least two are equal.
	{
		std::cout << "All roots real, at least two are equal" << std::endl;
		r13 = ((r < 0) ? -std::pow(-r, (1.0 / 3.0)) : std::pow(r, (1.0 / 3.0)));
		x1 = -term1 + 2.0 * r13;
		x3 = x2 = -(r13 + term1);
		return;
	} // End if (disc == 0)
	// Only option left is that all roots are real and unequal (to get here, q < 0)
	q = -q;
	dum1 = q * q * q;
	dum1 = std::acos(r / std::sqrt(dum1));
	r13 = 2.0 * std::sqrt(q);
	x1 = -term1 + r13 * std::cos(dum1 / 3.0);
	x2 = -term1 + r13 * std::cos((dum1 + 2.0 * M_PI) / 3.0);
	x3 = -term1 + r13 * std::cos((dum1 + 4.0 * M_PI) / 3.0);
	std::cout << "All things are ok" << std::endl;
	return;
}


void ApproxCCD::checkInterval(double t1, double t2, double* op, int degree, std::vector<TimeInterval>& intervals, bool pos)
{
	// clamp values
	if (t1 < 0.0) {
		t1 = 0.0;
	}
	if (t2 < 0.0) {
		t2 = 0.0;
	}
	if (t1 > 1.0) {
		t1 = 1.0;
	}
	if (t2 > 1.0) {
		t2 = 1.0;
	}
	//this is to compute ((ax+b)x+c)x+....
	double tmid = (t2 + t1) / 2;
	double f = op[0];
	for (int i = 1; i <= degree; i++)
	{
		f *= tmid;
		f += op[i];
	}
	//for pos =true, that is to check the intersection with plane constructed by edge and triangle normal. 
	//the sign of f will be consistent in [t1,t2], thus we use midpoint.
	//if the vertex actually intersect with triangle, the sign must be positive.
	//here is to exclude when the sign is negative,
	if (pos && f >= 0)
		intervals.push_back(TimeInterval(t1, t2));
	else if (!pos && f <= 0)
		intervals.push_back(TimeInterval(t1, t2));
}

void ApproxCCD::findIntervals(double* op, unsigned int n, std::vector<TimeInterval>& intervals, bool pos)
{
	unsigned int roots = 0;
	unsigned int reducedDegree = n;
	double time[6];
	std::uninitialized_fill(time, time + 6, 2.0);
	//normalized coeffs
	double maxval = 0;
	for (int i = 0; i <= n; i++) {
		if (maxval < fabs(op[i])) {
			maxval = fabs(op[i]);
		}
	}
	if (maxval != 0) {
		for (int i = 0; i <= n; i++) {
			op[i] /= maxval;
		}			
	}
	for (int i = 0; i < n; i++){
		if (op[i] == 0)
			reducedDegree--;
		else
			break;
	}	
	if (reducedDegree < n) {
		//move lower term coeff to higher term, this is for convenience
		for (int i = 0; i <= reducedDegree; i++) {
			op[i] = op[i + n - reducedDegree];
		}			
	}
	if (reducedDegree > 2) {
		xxxxxxxx;
	}
	else if (reducedDegree == 2) {
		roots = getQuadRoots(op[0], op[1], op[2], time[0], time[1]);
	}
	else if (reducedDegree == 1) {
		time[0] = -op[1] / op[0];
		roots = 1;
	}
	else { //reduced degree =0, for example: the vertex moving paralle to the plane
		// if  pos=true, the only coeff, coeff[0] must be positive (on the front side of the plane)
		if (!pos && op[0] <= 0 || (pos && op[0] >= 0)) {
			intervals.emplace_back(TimeInterval(0, 1.0));
		}
		return;
	}
	// need to check intervals
	//because time is the root, when pos =true, the vertex position in a time range won't change sideness of the plane
	if (roots > 0) {
		std::sort(time, time + roots);
		if (time[0] >= 0) {
			checkInterval(0, time[0], op, reducedDegree, intervals, pos);
		}
		for (unsigned int i = 0; i < roots - 1; ++i) {
			if (!((time[i] < 0 && time[i + 1] < 0) || (time[i] > 1.0 && time[i + 1] > 1.0))) {
				checkInterval(time[i], time[i + 1], op, reducedDegree, intervals, pos);
			}
		}
		if (time[roots - 1] <= 1.0)
			checkInterval(time[roots - 1], 1.0, op, reducedDegree, intervals, pos);
	}
	else {
		checkInterval(0.0, 1.0, op, reducedDegree, intervals, pos);
	}
}

