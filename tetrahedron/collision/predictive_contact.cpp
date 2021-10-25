#include"predictive_contact.h"

bool PredictiveContact::pointTriangleCollision(double* initial_position, double* current_position,
	std::vector<double*>& initial_triangle_position, std::vector<double*>& current_triangle_position, 
	double* initial_triangle_normal, double* current_triangle_normal,
	double* vertex_target_pos, std::vector<std::array<double, 3>>& triangle_target_pos, double radius,
	double normal_magnitude_reciprocal, double vertex_mass, double* triangle_mass)
{
	double nearest_point_barycentric[3];
	double initial_nearest_point[3], current_nearest_point[3];
	bool vertex_on_front;
	if (!vertexTriangleDistance(initial_position, initial_triangle_position, initial_triangle_normal, nearest_point_barycentric)) {
		NearestPointInfo nearest_point_info;
		vertexLineSegmentDistance(0, 1, initial_position, initial_triangle_position, nearest_point_info, nearest_point_barycentric);
		vertexLineSegmentDistance(1, 2, initial_position, initial_triangle_position, nearest_point_info, nearest_point_barycentric);
		vertexLineSegmentDistance(2, 0, initial_position, initial_triangle_position, nearest_point_info, nearest_point_barycentric);
		calNearestPoint(nearest_point_info, nearest_point_barycentric, initial_triangle_position, initial_nearest_point);
		calNearestPoint(nearest_point_info, nearest_point_barycentric, current_triangle_position, current_nearest_point);
	}
	else {
		calNearestPoint(nearest_point_barycentric, initial_triangle_position[0], initial_triangle_position[1], initial_triangle_position[2], initial_nearest_point);
		calNearestPoint(nearest_point_barycentric, current_triangle_position[0], current_triangle_position[1], current_triangle_position[2], current_nearest_point);
	}
	if (checkIfCollidePointTriangle(initial_position, current_position, initial_nearest_point, current_nearest_point,
		initial_triangle_normal, initial_triangle_position[0], radius, vertex_on_front))
	{
		obtainPointTriangleTargetPosition(current_position, current_triangle_position, vertex_target_pos, triangle_target_pos, radius, current_triangle_normal,
			vertex_on_front, normal_magnitude_reciprocal, vertex_mass, triangle_mass);
		return true;
	}
	return false;
}

bool PredictiveContact::pointBodyTriangleCollision(double* initial_position, double* current_position, std::vector<double*>& initial_triangle_position, std::vector<double*>& current_triangle_position,
	double* initial_triangle_normal, double* current_triangle_normal,
	double* vertex_target_pos, double radius)
{
	double nearest_point_barycentric[3];
	double initial_nearest_point[3], current_nearest_point[3];	
	if (!vertexTriangleDistance(initial_position, initial_triangle_position, initial_triangle_normal, nearest_point_barycentric)) {
		NearestPointInfo nearest_point_info;
		vertexLineSegmentDistance(0, 1, initial_position, initial_triangle_position, nearest_point_info, nearest_point_barycentric);
		vertexLineSegmentDistance(1, 2, initial_position, initial_triangle_position, nearest_point_info, nearest_point_barycentric);
		vertexLineSegmentDistance(2, 0, initial_position, initial_triangle_position, nearest_point_info, nearest_point_barycentric);
		calNearestPoint(nearest_point_info, nearest_point_barycentric, initial_triangle_position, initial_nearest_point);
		calNearestPoint(nearest_point_info, nearest_point_barycentric, current_triangle_position, current_nearest_point);
	}
	else {
		calNearestPoint(nearest_point_barycentric, initial_triangle_position[0], initial_triangle_position[1], initial_triangle_position[2], initial_nearest_point);
		calNearestPoint(nearest_point_barycentric, current_triangle_position[0], current_triangle_position[1], current_triangle_position[2], current_nearest_point);
	}
	double current_vertex_triangle_dis;
	if (checkIfCollidePointBodyTriangle(initial_position, current_position, initial_nearest_point, current_nearest_point,
		initial_triangle_normal,current_triangle_normal, current_vertex_triangle_dis,  radius)) {
		obtainPointBodyTriangleTargetPosition(current_position, vertex_target_pos, radius, current_triangle_normal, initial_triangle_normal, current_vertex_triangle_dis);
		return true;
	}
	return false;
}

bool PredictiveContact::edgeEdgeCollision(std::vector<std::array<double, 3>>& target_pos,
	std::vector<std::array<double, 3>>& compare_target_pos, double radius, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
	double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
	double* initial_compare_edge_vertex_1, double* mass)
{
	double alpha[4];
	double norm[3];
	if (!getClosestPoint(alpha, initial_edge_vertex_0, initial_edge_vertex_1,
		initial_compare_edge_vertex_0, initial_compare_edge_vertex_1)) {
		return false;
	}
	if (checkEdgeEdgeCollision(current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0, initial_edge_vertex_1, current_compare_edge_vertex_0, current_compare_edge_vertex_1,
		initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, alpha, norm, radius)) {
		obtaindgeEdgeTargetPosition(target_pos, compare_target_pos, norm, radius, alpha, current_edge_vertex_0, current_edge_vertex_1, current_compare_edge_vertex_0,
			current_compare_edge_vertex_1, mass);
		return true;
	}
	return false;
}

bool PredictiveContact::checkEdgeEdgeCollision(double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
	double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0, double* initial_compare_edge_vertex_1,
	double* alpha, double* norm, double radius)
{
	double e0[3];
	double e1[3];

	POINT_ON_EDGE(e0, alpha[0], alpha[1], initial_edge_vertex_0, initial_edge_vertex_1);
	POINT_ON_EDGE(e1, alpha[2], alpha[3], initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
	SUB(norm, e1, e0);
	double distance = sqrt(DOT(norm, norm));
	if (distance > NORM_NEAR_ZERO) {
		DEV(norm,norm,distance);
	}
	else {
		double dis2;
		double edge0[3], edge1[3];
		SUB(edge0, initial_edge_vertex_0, initial_edge_vertex_1);
		SUB(edge1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
		CROSS(norm, edge0, edge1);
		dis2 = sqrt(DOT(norm, norm));
		if (dis2 > NORM_NEAR_ZERO) {
			//std::cout << "One edge is getting too close to the other edge" << std::endl;
			DEV(norm, norm, dis2);
			SUB(edge0, initial_compare_edge_vertex_1, initial_edge_vertex_0);
			if (DOT(norm, edge0) < 0) {
				MULTI(norm, norm, -1.0);
			}
		}
		else {
			//std::cout << "One edge is getting too close to the other edge and they are parallel" << std::endl;
			return false;
		}
	}	
	double displacement0[3], displacement1[3];
	SUB(displacement0, current_edge_vertex_0, initial_edge_vertex_0);
	SUB(displacement1, current_edge_vertex_1, initial_edge_vertex_1);
	POINT_ON_EDGE(e0, alpha[0], alpha[1], displacement0, displacement1);
	SUB(displacement0, current_compare_edge_vertex_0, initial_compare_edge_vertex_0);
	SUB(displacement1, current_compare_edge_vertex_1, initial_compare_edge_vertex_1);
	POINT_ON_EDGE(e1, alpha[2], alpha[3], displacement0, displacement1);
	double dp_dc_project = DOT(norm, e0) - DOT(norm, e1);
	if (dp_dc_project > distance - radius) {
		return true;
	}
	return false;

}

bool PredictiveContact::vertexTriangleDistance(double* vertex, std::vector<double*>& triangle_position, double* triangle_normal, double* barycentric)
{
	double S[3];
	SUB(S, vertex, triangle_position[0]);
	double E1[3], E2[3], S1[3], S2[3];
	//double tnear;
	SUB(E1, triangle_position[1], triangle_position[0]);
	SUB(E2, triangle_position[2], triangle_position[0]);
	CROSS(S1, triangle_normal, E2);
	CROSS(S2, S, E1);
	double temp = 1.0 / DOT(S1, E1);
	//tnear = temp * DOT(S2, E2);
	barycentric[1] = temp * DOT(S1, S);
	barycentric[2] = temp * DOT(S2, triangle_normal);
	barycentric[0] = 1.0 - barycentric[1] - barycentric[2];
	if (barycentric[0] > EPSILON && barycentric[1] > EPSILON && barycentric[2] > EPSILON) {
		return true;
	}
	return false;
}

void PredictiveContact::calNearestPoint(double* barycentric, double* v0, double*v1, double*v2, double* center)
{
	BARYCENTRIC(center,barycentric, v0,v1,v2);
}

void PredictiveContact::vertexLineSegmentDistance(int edge_vertex_0, int edge_vertex_1, double* vertex, std::vector<double*>& triangle_position, 
	NearestPointInfo& nearest_point_info, double* barycentric)
{

	double AP[3], AB[3];
	double r, dot_, ab_2, ap, ab;

	SUB(AP, vertex, triangle_position[edge_vertex_0]);
	SUB(AB, triangle_position[edge_vertex_1],
		triangle_position[edge_vertex_0]);
	dot_ = DOT(AP, AB);
	ab_2 = DOT(AB, AB);
	r = dot_ / ab_2;
	if (r < 0) {
		ap = sqrt(DOT(AP, AP));
		if (ap < nearest_point_info.distance) {
			nearest_point_info.is_edge = false;
			nearest_point_info.index[0] = edge_vertex_0;
			nearest_point_info.distance = ap;
		}
	}
	else if (r > 1) {
		double BP[3], bp;
		SUB(BP, vertex, triangle_position[edge_vertex_1]);
		bp = sqrt(DOT(BP, BP));
		if (bp < nearest_point_info.distance) {
			nearest_point_info.is_edge = false;
			nearest_point_info.index[0] = edge_vertex_1;
			nearest_point_info.distance = bp;
		}
	}
	else {
		double ap_2 = DOT(AP, AP);
		double dis = sqrt(ap_2 - ab_2 * r * r);
		if (dis < nearest_point_info.distance) {
			nearest_point_info.is_edge = true;
			nearest_point_info.index[0] = edge_vertex_0;
			nearest_point_info.index[1] = edge_vertex_1;
			barycentric[0] = 1.0 - r;
			barycentric[1] = r;
			nearest_point_info.distance = dis;
		}
	}
}

void PredictiveContact::calNearestPoint(NearestPointInfo& nearest_point_info, double* barycentric, std::vector<double*>& triangle_position, double* center)
{
	if (nearest_point_info.is_edge) {
		pointOnEdge(center, barycentric[0], barycentric[1], triangle_position[nearest_point_info.index[0]],
			triangle_position[nearest_point_info.index[1]]);
	}
	else {
		memcpy(center, triangle_position[nearest_point_info.index[0]], 24);		
	}
}

void PredictiveContact::pointOnEdge(double* center, double coe0, double coe1, double* vertex0, double* vertex1)
{
	center[0]= coe0 * vertex0[0] + coe1 * vertex1[0];
	center[1]= coe0 * vertex0[1] + coe1 * vertex1[1];
	center[2]= coe0 * vertex0[2] + coe1 * vertex1[2];
}

bool PredictiveContact::checkIfCollidePointTriangle(double* initial_point_position, double* current_point_position, 
	double* initial_nearest_point, double* current_nearest_point, double* initial_triangle_normal,
	double* initial_triangle_v0, double radius, bool& vertex_on_front) {
	double p_c[3];
	double distance_p_c;

	SUB(p_c, initial_point_position, initial_triangle_v0);
	if (DOT(p_c, initial_triangle_normal) > 0) {
		vertex_on_front = true;
	}
	else {
		vertex_on_front = false;
	}
	SUB(p_c, initial_point_position, initial_nearest_point);
	distance_p_c = sqrt(DOT(p_c, p_c));
	
	if (distance_p_c > NORM_NEAR_ZERO) {
		DEV(p_c, p_c, distance_p_c);
	}
	else {
		//std::cout << "the distance between point and its closest point on the triangle is too close: "<<distance_p_c<< std::endl;
		if (vertex_on_front) {
			memcpy(p_c, initial_triangle_normal, 24);
		}
		else {
			MULTI(p_c, initial_triangle_normal, -1.0);
		}
	}

	double displacement_c[3], displacement_p[3];
	SUB(displacement_c, current_nearest_point, initial_nearest_point);
	SUB(displacement_p, current_point_position, initial_point_position);

	double dp_dc_project = DOT(displacement_c, p_c) -DOT(displacement_p, p_c);

	if (dp_dc_project > distance_p_c - radius) {
		return true;
	}
	return false;
}
bool PredictiveContact::checkIfCollidePointBodyTriangle(double* initial_point_position, double* current_point_position,
	double* initial_nearest_point, double* current_nearest_point, double* initial_triangle_normal, double* current_triangle_normal, double& current_vertex_triangle_dis,double radius)
{
	double p_c[3];

	SUB(p_c, current_point_position, current_nearest_point);
	current_vertex_triangle_dis = DOT(p_c, current_triangle_normal);
	if (current_vertex_triangle_dis > radius) {
		return false;
	}

	double distance_p_c;

	SUB(p_c, initial_point_position, initial_nearest_point);
	distance_p_c = sqrt(DOT(p_c, p_c));
	if (distance_p_c > NORM_NEAR_ZERO) {
		DEV(p_c, p_c, distance_p_c);
	}
	else {
		memcpy(p_c, initial_triangle_normal, 24);
	}
	double displacement_c[3], displacement_p[3];
	SUB(displacement_c, current_nearest_point, initial_nearest_point);
	SUB(displacement_p, current_point_position, initial_point_position);
	double dp_dc_project = DOT(displacement_c, p_c) - DOT(displacement_p, p_c);

	if (dp_dc_project > distance_p_c - radius) {
		return true;
	}
	return false;
}

void PredictiveContact::obtainPointTriangleTargetPosition(double* point_position, std::vector<double*>& triangle_position,
	double* vertex_target_pos, std::vector<std::array<double, 3>>& triangle_target_pos, double radius, double* triangle_normal, bool vertex_on_front,
	double normal_magnitude_reciprocal, double vertex_mass, double* triangle_mass)
{
	triangle_target_pos.resize(3);
	double temp_vec[3];
	SUB(temp_vec, point_position, triangle_position[0]);
	double constraint = DOT(temp_vec, triangle_normal);
	double in_triangle[3], scale_norm[3];
	SUB(in_triangle, point_position, triangle_position[0]);
	MULTI(scale_norm, triangle_normal, constraint);
	SUB(in_triangle, in_triangle, scale_norm);
	if (vertex_on_front) {
		constraint -= radius;
	}
	else {
		constraint += radius;
	}
	memcpy(vertex_target_pos, triangle_normal, 24);
	SUB(triangle_target_pos[0], triangle_position[1], triangle_position[2]);
	SUB(triangle_target_pos[1], triangle_position[2], triangle_position[0]);
	SUB(triangle_target_pos[2], triangle_position[0], triangle_position[1]);
	double temp_vec[3];
	for (int i = 0; i < 3; ++i) {
		CROSS(temp_vec, triangle_target_pos[i], in_triangle);
		MULTI(triangle_target_pos[i], temp_vec, normal_magnitude_reciprocal);
	}
	SUB(triangle_target_pos[0], triangle_target_pos[0], vertex_target_pos);
	double s = DOT(vertex_target_pos, vertex_target_pos) / vertex_mass;
	for (int i = 0; i < 3; ++i) {
		s += DOT(triangle_target_pos[i], triangle_target_pos[i]) / triangle_mass[i];
	}
	s = -constraint / s;
	double temp_value;
	temp_value = s / vertex_mass;
	MULTI(vertex_target_pos, vertex_target_pos, temp_value);
	SUM(vertex_target_pos, vertex_target_pos, point_position);
	for (int i = 0; i < 3; ++i) {
		temp_value = s / triangle_mass[i];
		MULTI(triangle_target_pos[i], triangle_target_pos[i], temp_value);
		SUM(triangle_target_pos[i], triangle_target_pos[i], triangle_position[i]);
	}
}

void PredictiveContact::obtainPointBodyTriangleTargetPosition(double* point_position,double* vertex_target_pos, 
	double radius, double* current_triangle_normal, double* initial_triangle_normal, double current_vertex_triangle_dis)
{
	current_vertex_triangle_dis -= radius;
	double angle= DOT(initial_triangle_normal, current_triangle_normal);
	if (abs(angle) > NORM_NEAR_ZERO) {
		current_vertex_triangle_dis /= angle;
		double displacement[3];
		MULTI(displacement, initial_triangle_normal, current_vertex_triangle_dis);
		SUB(vertex_target_pos, vertex_target_pos, displacement);
	}
	else {
		double displacement[3];
		MULTI(displacement, initial_triangle_normal, current_vertex_triangle_dis);
		SUB(vertex_target_pos, vertex_target_pos, displacement);
	}
}

bool PredictiveContact::getClosestPoint(double* alpha, double* p1, double* p2, double* p3, double* p4)
{
	double s1[3];
	SUB(s1, p2, p1);
	double s2[3];
	SUB(s2, p4, p3);
	double s1s2 = DOT(s1, s2);
	double s1_2 = DOT(s1, s1);
	double s2_2 = DOT(s2, s2);
	double p3p1[3];
	SUB(p3p1, p1, p3);
	double sub_ = s1_2 * s2_2 - s1s2 * s1s2;
	if (sub_ < NEAR_ZERO2) {   //NEAR_ZERO2 is too small
		return false;
	}
	else {
		double a0 = DOT(p3p1, s1);
		double a1 = DOT(p3p1, s2);

		alpha[1] = (s1s2 * a1 - s2_2 * a0) / sub_;
		alpha[3] = -(s1s2 * a0 - s1_2 * a1) / sub_;
		if (alpha[1] < 0) {
			getClosestDistanceCornerCornerCase(alpha, p1, p2, p3, p4);
			return true;
		}
		else if (alpha[1] > 1.0) {
			getClosestDistanceCornerCornerCase(alpha, p1, p2, p3, p4);
			return true;
		}
		if (alpha[3] < 0) {
			getClosestDistanceCornerCornerCase(alpha, p1, p2, p3, p4);
			return true;
		}
		else if (alpha[3] > 1.0) {
			getClosestDistanceCornerCornerCase(alpha, p1, p2, p3, p4);
			return true;
		}
		alpha[0] = 1.0 - alpha[1];
		alpha[2] = 1.0 - alpha[3];
	}
	return true;
}

void PredictiveContact::getClosestDistanceCornerCornerCase(double* alpha, double* p1, double* p2, double* p3, double* p4)
{
	double alpha_edge[2];
	double distance = DBL_MAX;
	alpha[0] = 1.0;
	alpha[1] = 0.0;
	getClosestPointBetweenPointSegement(alpha_edge, p1, p3, p4, distance);
	alpha[2] = alpha_edge[0];
	alpha[3] = alpha_edge[1];

	if (getClosestPointBetweenPointSegement(alpha_edge, p2, p3, p4, distance)) {
		alpha[0] = 0.0;
		alpha[1] = 1.0;
		alpha[2] = alpha_edge[0];
		alpha[3] = alpha_edge[1];
	}
	if (getClosestPointBetweenPointSegement(alpha_edge, p3, p1, p2, distance)) {
		alpha[0] = alpha_edge[0];
		alpha[1] = alpha_edge[1];
		alpha[2] = 1.0;
		alpha[3] = 0.0;
	}
	if (getClosestPointBetweenPointSegement(alpha_edge, p4, p1, p2, distance)) {
		alpha[0] = alpha_edge[0];
		alpha[1] = alpha_edge[1];
		alpha[2] = 0.0;
		alpha[3] = 1.0;
	}
}

bool PredictiveContact::getClosestPointBetweenPointSegement(double* alpha, double* p0, double* e0, double* e1, double& distance)
{
	double AP[3];
	double AB[3];
	SUB(AP, p0, e0);
	SUB(AB, e1, e0);
	double dot_ = DOT(AP, AB);
	double ab_2 = DOT(AB, AB);
	double r = dot_ / ab_2;
	if (r < 0) {
		double ap = sqrt(DOT(AP, AP));
		if (ap < distance) {
			alpha[0] = 1.0;
			alpha[1] = 0.0;
			distance = ap;
			return true;
		}
	}
	else if (r > 1) {
		double BP[3], bp;
		SUB(BP, p0, e1);
		bp = sqrt(DOT(BP, BP));
		if (bp < distance) {
			alpha[0] = 0.0;
			alpha[1] = 1.0;
			distance = bp;
			return true;
		}
	}
	else {
		double ap_2 = DOT(AP, AP);
		double dis = sqrt(ap_2 - ab_2 * r * r);
		if (dis < distance) {
			alpha[0] = 1.0 - r;
			alpha[1] = r;
			distance = dis;
			return true;
		}
	}
	return false;
}


void PredictiveContact::obtaindgeEdgeTargetPosition(std::vector<std::array<double, 3>>& target_pos, std::vector<std::array<double, 3>>& compare_target_pos,
	double norm[3], double radius, double* alpha, double* current_edge_vertex_0, double* current_edge_vertex_1,
	double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* mass)
{
	target_pos.resize(2);
	compare_target_pos.resize(2);

	double distance;
	double e0[3], e1[3];
	POINT_ON_EDGE(e0, alpha[0], alpha[1], current_edge_vertex_0, current_edge_vertex_1);
	POINT_ON_EDGE(e1, alpha[2], alpha[3], current_compare_edge_vertex_0, current_compare_edge_vertex_1);
	SUB(e0, e1, e0);
	distance = DOT(norm, e0)-radius;

	double coe = 0.0;
	for (int i = 0; i < 4; ++i) {
		coe += alpha[i] * alpha[i] / mass[i];
	}
	double temp;
	distance /= coe;
	temp = distance * alpha[0] / mass[0];
	MULTI(target_pos[0], norm, temp);
	SUM(target_pos[0], current_edge_vertex_0, target_pos[0]);

	temp = distance * alpha[1] / mass[1];
	MULTI(target_pos[1], norm, temp);
	SUM(target_pos[1], current_edge_vertex_1, target_pos[1]);

	temp = -distance * alpha[2] / mass[2];
	MULTI(compare_target_pos[0], norm, temp);
	SUM(compare_target_pos[0], current_compare_edge_vertex_0, compare_target_pos[0]);

	temp = -distance * alpha[3] / mass[3];
	MULTI(compare_target_pos[1], norm, temp);
	SUM(compare_target_pos[1], current_compare_edge_vertex_1, compare_target_pos[1]);
}