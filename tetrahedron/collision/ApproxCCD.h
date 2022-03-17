#pragma once
#define NOMINMAX
#include"root_finder.h"
#include"rootparitycollisiontest.h"
#include"TightCCD.h"

class ApproxCCD
{
public:	
	bool pointTriangleCollisionTime(double& t, double* initial_position, double* current_position,
		double* initial_triangle_0, double* current_triangle_0, double* initial_triangle_1, double* current_triangle_1, double* initial_triangle_2, double* current_triangle_2,
		double* initial_normal_not_normalized, double* current_normal_not_normalized, double* cross_for_CCD, double tolerance_2);//floating* f_initial_normal, floating* f_current_normal, floating* f_cross_for_CCD, 
	bool edgeEdgeCollisionTime(double& t, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
		double* initial_compare_edge_vertex_1, double tolerance_2);
	bool VertexTriangleCollisionTime(double& t, double* initial_position, double* current_position,
		double* initial_triangle_1, double* current_triangle_1, double* initial_triangle_2, double* current_triangle_2, double* initial_triangle_3, double* current_triangle_3,
		double eta);

	void test();
private:

	struct TimeInterval
	{
		double l, u;
		TimeInterval(double tl, double tu) : l(tl), u(tu)
		{
			if (l > u) std::swap(l, u);
			l = myMax(1, 0.0);
			u = myMin(u, 1.0);
		}
		TimeInterval() : l(0), u(0) {}

		// Returns whether or not the intersection of the intervals is nonempty
		static bool overlap(const TimeInterval& t1, const TimeInterval& t2);
		static bool overlap(const std::vector<TimeInterval>& intervals);

		// Returns the intersection of the intervals **asuming the intersection is nonempty**
		static TimeInterval intersect(const std::vector<TimeInterval>& intervals);
	};


	bool solveEquation(double& t0, double& t1, double& t2, double a3, double a2, double a1, double d);
	bool solveQuadraticEquation(double& t0, double& t1, double a2, double a1, double a0);
	bool pointEdgeCollisionTime(double& t, double* u, double* u0, double* u1, double* e_1_0,
		double* e0, double* e1, double tolerance_2);
	void pointEdge2D(std::vector<double>& t, double* u, double* u0, double* u1, double* e_1_0, double* e0, double* e1,
		double* u_0, double* u_1, double* u_10);
	void solveQuadratic(double* t, double a, double b, double c, int& root_number);
	bool pointEdgeIsClose(double* v_0, double* v_1_0, double tolerance_2);
	bool pointPointCollisionTime(double& t, double* e_1_0, double* e_0, double* e_1, double tolerance_2);
	bool pointPointIsClose(double t, double* e10_0, double* e_0, double* e_1, double tolerance_2);
	void make_vector(double* v, floating* out);
	TightCCD tight_CCD;
	bool checkInside(double t, double* v0, double* v1, double* v2, double* v3,
		double* e0, double* e1, double* e2, double* e3);
	bool edgeEdgeCheckInside(double t, double* v0, double* v1, double* v2, double* v3,
		double* e0, double* e1, double* e2, double* e3);
	bool solveCubicEquation(double a, double b, double c, double d, double& t0, double& t1, double& t2);
	void sortABC(double& a, double& b, double& c);
	void cubicsolve(const double& a, const double& b, const double& c, const double& d, double& x1, double& x2, double& x3);

	void testInsideOutside(double t, double* initial_pos, double* u, double* initial_triangle_0,
		double* u_0, double* initial_triangle_1, double* u_1, double* initial_triangle_2, double* u_2);
	void testEdgeInsideOutside(double t, double* initial_pos_0_0, double* u_0_0, double* initial_pos_0_1,
		double* u_0_1, double* initial_triangle_1_0, double* u_1_0, double* initial_triangle_1_1, double* u_1_1);

	void planePoly3D(double* x10, double* x20, double* x30, double* v10, double* v20, double* v30,
		std::vector<TimeInterval>& result);
	void checkInterval(double t1, double t2, double* op, int degree, std::vector<TimeInterval>& intervals, bool pos);
	void findIntervals(double* op, unsigned int n, std::vector<TimeInterval>& intervals, bool pos);
	int getQuadRoots(double a, double b, double c, double& t0, double& t1);
};


