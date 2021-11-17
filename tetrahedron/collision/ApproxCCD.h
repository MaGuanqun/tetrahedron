#pragma once
#include<vector>

class ApproxCCD
{
public:
	bool pointTriangleCollisionTime(double& alpha, double* initial_position, double* current_position,
		double* initial_triangle_0, double* current_triangle_0, double* initial_triangle_1, double* current_triangle_1,
		double* initial_triangle_2, double* current_triangle_2,
		double* initial_normal_not_normalized, double* current_normal_not_normalized, double* cross_for_CCD);
	bool edgeEdgeCollisionTime(double& alpha, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
		double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
		double* initial_compare_edge_vertex_1);
private:

};


