#pragma once
#include"../basic/global.h"
#include<array>

class PredictiveContact
{
public:

	bool pointTriangleCollision(double* initial_position, double* current_position, std::vector<double*>& initial_triangle_position, std::vector<double*>& current_triangle_position,
		double* initial_triangle_normal, double* current_triangle_normal,
		double* vertex_target_pos, std::vector<std::array<double, 3>>& triangle_target_pos, double radius,
		double normal_magnitude_reciprocal, double vertex_mass, double* triangle_mass);
	bool pointBodyTriangleCollision(double* initial_position, double* current_position, std::vector<double*>& initial_triangle_position, std::vector<double*>& current_triangle_position,
		double* initial_triangle_normal, double* current_triangle_normal,
		double* vertex_target_pos, double radius);
	bool edgeEdgeCollision(std::vector<std::array<double, 3>>& target_pos,
		std::vector<std::array<double, 3>>& compare_target_pos, double radius, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
		double* initial_compare_edge_vertex_1, double* mass);


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
	bool vertexTriangleDistance(double* vertex, std::vector<double*>& triangle_position, double* triangle_normal, double* barycentric);
	void calNearestPoint(double* barycentric, double* v0, double* v1, double* v2, double* center);
	void vertexLineSegmentDistance(int edge_vertex_0, int edge_vertex_1, double* vertex, std::vector<double*>& triangle_position, 
		NearestPointInfo& nearest_point_info, double* barycentric);
	void calNearestPoint(NearestPointInfo& nearest_point_info, double* barycentric, std::vector<double*>& triangle_position, double* center);
	void pointOnEdge(double* center, double coe0, double coe1, double* vertex0, double* vertex1);
	bool checkIfCollidePointTriangle(double* initial_point_position, double* current_point_position, double* initial_nearest_point, double* current_nearest_point, double* initial_triangle_normal,
		double* initial_triangle_v0, double radius, bool& vertex_on_front);
	void obtainPointTriangleTargetPosition(double* point_position, std::vector<double*>& triangle_position,
		double* vertex_target_pos, std::vector<std::array<double, 3>>& triangle_target_pos, double radius, double* triangle_normal, bool vertex_on_front,
	double normal_magnitude_reciprocal, double vertex_mass, double* triangle_mass);
	bool checkIfCollidePointBodyTriangle(double* initial_point_position, double* current_point_position,
		double* initial_nearest_point, double* current_nearest_point, double* initial_triangle_normal,
		double* current_triangle_normal, double& current_vertex_triangle_dis, double radius);
	void obtainPointBodyTriangleTargetPosition(double* point_position,double* vertex_target_pos, double radius,
		double* current_triangle_normal, double* initial_triangle_normal, double current_vertex_triangle_dis);
	bool getClosestPoint(double* alpha, double* p1, double* p2, double* p3, double* p4);
	void getClosestDistanceCornerCornerCase(double* alpha, double* p1, double* p2, double* p3, double* p4);
	bool getClosestPointBetweenPointSegement(double* alpha, double* p0, double* e0, double* e1, double& distance);
	bool checkEdgeEdgeCollision(double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0, double* initial_compare_edge_vertex_1,
		double* alpha, double* norm, double radius);
	void obtaindgeEdgeTargetPosition(std::vector<std::array<double, 3>>& target_pos, std::vector<std::array<double, 3>>& compare_target_pos,
		double norm[3], double radius, double* alpha, double* current_edge_vertex_0, double* current_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* mass);
};


