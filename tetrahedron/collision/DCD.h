#pragma once

#include"../basic/global.h"
#include"primitive_distance.h"
#include"../external/Eigen/Dense"
#include"../basic/eigenDenseOperation.h"
#include "../external/Eigen/Eigenvalues"

using namespace Eigen;
using namespace denseOperation;

#undef MAGNITUDE_CROSS
#define MAGNITUDE_CROSS(dest,v1,v2, t) \
dest[0] =t*(v1[1] * v2[2] - v1[2] * v2[1]); \
dest[1] =t*(v1[2] * v2[0] - v1[0] * v2[2]); \
dest[2] =t*(v1[0] * v2[1] - v1[1] * v2[0]);

#undef MULTI_SUM
#define MULTI_SUM(dest, t,v0)\
dest[0]=t*dest[0]+v0[0];\
dest[1]=t*dest[1]+v0[1];\
dest[2]=t*dest[2]+v0[2];\

#undef MULTI_SUM2
#define MULTI_SUM2(dest,t,u0,v0)\
dest[0]=t*u0[0]+v0[0];\
dest[1]=t*u0[1]+v0[1];\
dest[2]=t*u0[2]+v0[2];\

class DCD
{
public:
	bool pointSelfTriangle(double* initial_position, double* current_position,
		double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
		double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* initial_triangle_normal, double* current_triangle_normal, double* vertex_target_pos,
		double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
		double tolerance, double mass_point, double mass_t0, double mass_t1, double mass_t2,
		double triangle_normal_magnitude_reciprocal);

	bool pointColliderTriangle(double* initial_position, double* current_position,
		double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
		double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* initial_triangle_normal, double* current_triangle_normal,
		double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
		double tolerance, double mass_t0, double mass_t1, double mass_t2,
		double triangle_normal_magnitude_reciprocal);


	bool pointTriangleCollider(double* initial_position, double* current_position,
		double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
		double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* initial_triangle_normal, double* current_triangle_normal, double* vertex_target_pos,
		double tolerance);


	bool edgeEdgeCollider(double* edge_target_pos_0, double* edge_target_pos_1,
		double* current_edge_vertex_0, double* current_edge_vertex_1,
		double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
		double* initial_compare_edge_vertex_1,
		double tolerance, double mass_e_0_0, double mass_e_0_1);

	bool edgeEdge(double* edge_target_pos_0, double* edge_target_pos_1,
		double* compare_target_pos_0, double* compare_target_pos_1, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
		double* initial_compare_edge_vertex_1, double tolerance, double mass_e_0_0, double mass_e_0_1, double mass_e_1_0, double mass_e_1_1);


	bool checkPointTriangle(double* initial_position, double* current_position,
		double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
		double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* initial_triangle_normal, double tolerance);

	bool checkEdgeEdge(double* current_edge_vertex_0, double* current_edge_vertex_1,
		double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
		double* initial_compare_edge_vertex_1, double tolerance);

	bool XPBDpointSelfTriangle(double* initial_position, double* current_position,
		double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
		double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* initial_triangle_normal, double* current_triangle_normal,
		double tolerance, double mass_inv_point, double mass_inv_t0, double mass_inv_t1, double mass_inv_t2,
		double triangle_normal_magnitude_reciprocal, double* delta_vertex, double* delta_triangle0, double* delta_triangle1, 
		double* delta_triangle2);
	bool XPBDedgeEdge(double* current_edge_vertex_0, double* current_edge_vertex_1,
		double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
		double* initial_compare_edge_vertex_1, double tolerance, double mas_inv_e_0_0, double mass_inv_e_0_1, double mass_inv_e_1_0, double mass_inv_e_1_1,
		double* delta_edge_0, double* delta_edge_1, double* delta_compare_edge_0, double* delta_compare_edge_1);

	bool XPBDpointTriangleCollider(double* initial_position, double* current_position,
		double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
		double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* initial_triangle_normal, double* current_triangle_normal, double mass_inv_v,
		double tolerance, double* delta_x);

	
	void XPBDFloor(double* initial_position, double* current_position, unsigned int dimension, bool normal_direction, double mass_inv_v,
		double tolerance, double& lambda, double stiffness, double damping_stiffness, double dt, double floor_value, double& energy);


	bool accuratePointSelfTriangle(double* initial_position, double* current_position,
		double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
		double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* initial_triangle_normal, double* current_triangle_normal, double* vertex_target_pos,
		double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
		double tolerance, double mass_point, double mass_t0, double mass_t1, double mass_t2,
		double current_triangle_area);

	bool PDFloor(double* target_position, double* current_position, unsigned int dimension, bool normal_direction,
		double tolerance, double floor_value);

	void test();

private:
	bool pointProjectOnTriangle(
		const double* p,
		const double* t0,
		const double* t1,
		const double* t2,
		const double* triangle_normal, double* barycentric);

	bool checkEdgeEdgeCollision(double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0, double* initial_compare_edge_vertex_1,
		double* alpha, double* norm, double& distance2, double tolerance, bool collider);
	void calDistanceEdgeEdge(double* edge_target_pos_0, double* edge_target_pos_1,
		double* compare_target_pos_0, double* compare_target_pos_1,
		double* norm, double distance, double* alpha, double* current_edge_vertex_0, double* current_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1,
		double mass_e_0_0, double mass_e_0_1, double mass_e_1_0, double mass_e_1_1);
	void calDistanceEdgeEdgeCollider(double* edge_target_pos_0, double* edge_target_pos_1,
		double* norm, double distance, double* alpha, double* current_edge_vertex_0, double* current_edge_vertex_1,
		double mass_e_0_0, double mass_e_0_1);
	bool checkIfCollidePointTriangle(double* initial_point_position, double* current_point_position,
		double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
		double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* initial_triangle_normal, double* current_triangle_normal, double* barycentric,
		double tolerance, double& triangle_side2, bool& should_be_front);
	bool checkIfCollidePointTriangle(double* initial_point_position, double* current_point_position,
		double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
		double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* initial_triangle_normal, double* barycentric, double tolerance);
	void calDistancePointTriangle(double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
		double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* current_triangle_normal, double constraint, double tolerance, bool is_front, double triangle_normal_magnitude_reciprocal,
		double mass_point, double mass_t0, double mass_t1, double mass_t2);
	bool checkIfCollideEdgeEdge(double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0, double* initial_compare_edge_vertex_1,
		double* alpha, double tolerance);
	bool XPBDcalDistancePointTriangle(
		double* initial_position,
		double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
		double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* current_triangle_normal, double constraint, double tolerance, bool is_front, double triangle_normal_magnitude_reciprocal,
		double mass_inv_point, double mass_inv_t0, double mass_inv_t1, double mass_inv_t2, 
		double* delta_vertex, double* delta_triangle0, double* delta_triangle1, double* delta_triangle2);
	void XPBDcalDistanceEdgeEdge(double* norm, double distance, double* alpha, double* current_edge_vertex_0, double* current_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1,
		double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* initial_compare_edge_vertex_0, double* initial_compare_edge_vertex_1,
		double mass_inv_e_0_0, double mass_inv_e_0_1, double mass_inv_e_1_0, double mass_inv_e_1_1,
		double* delta_edge_0, double* delta_edge_1, double* delta_compare_edge_0, double* delta_compare_edge_1);
	void XPBDcalDistancePointTriangleCollider(double* initial_position,
		double* current_position, double mass_inv_vertex,
		double* current_triangle_normal, double constraint, double tolerance,
		double* delta_x);
	bool checkIfCollidePointTriangleCollider(double* initial_point_position, double* current_point_position,
		double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
		double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* initial_triangle_normal, double* barycentric,
		double tolerance);

	void calAccurateDistancePointTriangle(double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
		double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double tolerance, double mass_point, double mass_t0, double mass_t1, double mass_t2);

	double distanceFromTriangle(double* vertex_target_pos,
		double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2);
	void test(double* info_p, double* info_v0, double* info_v1, double* info_v2);

	void iterationToGetResult(double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
		double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* triangle_normal_ori, double constraint, double tolerance, bool is_front, double triangle_normal_magnitude_reciprocal,
		double mass_point, double mass_t0, double mass_t1, double mass_t2);
	void calPositionMove(double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
		double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* current_triangle_normal, double constraint, double tolerance, bool is_front, double triangle_normal_magnitude_reciprocal,
		double mass_point, double mass_t0, double mass_t1, double mass_t2, double& lambda);
	void computeNormalAndConstraint(double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
		double* normal, double& area, double& constraint, bool is_front, double tolerance);

	bool convergeCondition(int& itr_num,double lambda, double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
		double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double mass_point, double mass_t0, double mass_t1, double mass_t2, double constraint, double is_front,
		double& L_);
	void decideTolerance(double& tolerance, double ratio_tolerance, double ratio_max_move, double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2);
	void calDistancePointColliderTriangle(double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
		double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
		double* current_triangle_normal, double constraint, double tolerance, bool is_front, double triangle_normal_magnitude_reciprocal,
		double mass_t0, double mass_t1, double mass_t2);

};



