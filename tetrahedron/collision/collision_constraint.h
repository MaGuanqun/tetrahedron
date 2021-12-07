#pragma once

#include"../basic/global.h"
#include<array>

class CollisionConstraint
{
public:
	bool pointSelfTriangle(double* initial_position, double* current_position,
		std::vector<double*>& initial_triangle_position, std::vector<double*>& current_triangle_position,
		double* initial_triangle_normal, double* vertex_target_pos,
		std::array<double, 3>* triangle_target_pos, double d_hat_2, double& stiffness,
		double vertex_mass, double* triangle_mass);
	bool pointColliderTriangle(double* initial_position, double* current_position,
		std::vector<double*>& current_triangle_position, double* triangle_normal, double* vertex_target_pos,
		double d_hat_2, double& stiffness);
	bool edgeEdgeCollision(std::vector<std::array<double, 3>>& target_pos,
		std::vector<std::array<double, 3>>& compare_target_pos, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
		double* initial_compare_edge_vertex_1, double* mass, double d_hat_2, double& stiffness);

	void testPT();
	void testEE();

private:
	struct NearestPointInfo
	{
		double distance;
		bool is_edge;//the closest point on triangle is 0 on triangle 1 edge 2 point
		int index[2];
		NearestPointInfo() {
			distance = DBL_MAX;
		}
	};
	bool vertexTriangleDistance(double* vertex, double& d_2, std::vector<double*>& triangle_position, double* triangle_normal, double* barycentric);
	void calNearestPoint(double* barycentric, double* v0, double* v1, double* v2, double* center);
	void vertexLineSegmentDistance(int edge_vertex_0, int edge_vertex_1, double* vertex, std::vector<double*>& triangle_position,
		NearestPointInfo& nearest_point_info, double* barycentric);
	void calNearestPoint(NearestPointInfo& nearest_point_info, double* barycentric, std::vector<double*>& triangle_position, double* center);
	void pointOnEdge(double* center, double coe0, double coe1, double* vertex0, double* vertex1);
	double barrier(double d_2_minus_d_hat_2_over_d_hat2, double d_2_div_d_hat_2);
	void moveDistance(double vertex_mass, double* triangle_mass, double* vertex_target_pos, std::array<double, 3>* triangle_target_pos,
		double d_move, double* barycentric, double* direction, double* initial_pos, std::vector<double*>& initial_triangle_position);
	void setBarycentric(NearestPointInfo& nearest_point_info, double* barycentric, double* edge_barycentric);
	bool getClosestPoint(double* alpha, double* p1, double* p2, double* p3, double* p4, double d_hat_2, double* c1, double* c2, double& distance);
	bool getClosestDistanceCornerCornerCase(double* alpha, double* p1, double* p2, double* p3, double* p4, double d_hat_2, double& distance);
	bool getClosestPointBetweenPointSegement(double* alpha, double* p0, double* e0, double* e1, double& distance);
	void moveDistance(double* mass, std::vector<std::array<double, 3>>& target_pos,
		std::vector<std::array<double, 3>>& compare_target_pos, double d_move, double* barycentric, double* direction, double* initial_edge_vertex_0,
		double* initial_edge_vertex_1, double* initial_compare_edge_vertex_0, double* initial_compare_edge_vertex_1);
};


