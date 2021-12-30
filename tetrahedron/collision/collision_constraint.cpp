#include"collision_constraint.h"

bool CollisionConstraint::pointSelfTriangle(double* initial_position, double* current_position,
	std::vector<double*>& initial_triangle_position, std::vector<double*>& current_triangle_position,
	double* initial_triangle_normal, double* vertex_target_pos,
	std::array<double, 3>* triangle_target_pos, double d_hat, double& stiffness,
	double vertex_mass, double* triangle_mass)
{
	double d_2;
	double nearest_point_barycentric[3];
	double initial_nearest_point[3], current_nearest_point[3];
	double d_hat_2 = d_hat * d_hat;
	if (!vertexTriangleDistance(initial_position, d_2, initial_triangle_position, initial_triangle_normal, nearest_point_barycentric)) {
		NearestPointInfo nearest_point_info;
		double barycentric_edge[2];
		vertexLineSegmentDistance(0, 1, initial_position, initial_triangle_position, nearest_point_info, barycentric_edge);
		vertexLineSegmentDistance(1, 2, initial_position, initial_triangle_position, nearest_point_info, barycentric_edge);
		vertexLineSegmentDistance(2, 0, initial_position, initial_triangle_position, nearest_point_info, barycentric_edge);
		if (nearest_point_info.distance >= d_hat_2) {
			return false;
		}
		d_2 = nearest_point_info.distance;
		calNearestPoint(nearest_point_info, barycentric_edge, initial_triangle_position, initial_nearest_point);
		calNearestPoint(nearest_point_info, barycentric_edge, current_triangle_position, current_nearest_point);
		setBarycentric(nearest_point_info, nearest_point_barycentric, barycentric_edge);

	}
	else {
		if (d_2 >= d_hat_2) {
			return false;
		}
		calNearestPoint(nearest_point_barycentric, initial_triangle_position[0], initial_triangle_position[1], initial_triangle_position[2], initial_nearest_point);
		calNearestPoint(nearest_point_barycentric, current_triangle_position[0], current_triangle_position[1], current_triangle_position[2], current_nearest_point);
	}


	stiffness *= barrier((d_hat_2 - d_2)/d_hat_2, d_2 / d_hat_2);
	double d[3];
	for (int i = 0; i < 3; ++i) {
		d[i] = initial_position[i] - current_position[i] - (initial_nearest_point[i] - current_nearest_point[i]);
	}
	double d_move = sqrt(DOT(d, d));
	double d_min= d_hat - sqrt(d_2);
	if (d_move < d_min) {
		d_move = d_min;
	}
	//else {
	//	std::cout << d_move<<" "<< d_min << std::endl;
	//}
	//decide the direction
	double direction[3];
	SUB(direction, initial_nearest_point, initial_position);
	double t1 = sqrt(DOT(direction, direction));
	if (t1 > NEAR_ZERO) {
		DEV_(direction, t1);
	}
	else {
		if (DOT(initial_triangle_normal, direction) < 0) {
			MULTI(direction, initial_triangle_normal, -1.0);
		}
		else {
			memcpy(direction, initial_triangle_normal, 24);
		}
	}	
	//std::cout << d_move<<" "<<stiffness << std::endl;
	moveDistance(vertex_mass, triangle_mass, vertex_target_pos, triangle_target_pos,d_move, nearest_point_barycentric, direction,initial_position, initial_triangle_position);
	//std::cout << initial_triangle_position[0][0] << " " << initial_triangle_position[0][1] << " " << initial_triangle_position[0][2] << std::endl;
	//std::cout << triangle_target_pos[0][0] << " " << triangle_target_pos[0][1] << " " << triangle_target_pos[0][2] << std::endl;
	//std::cout << current_triangle_position[0][0] << " " << current_triangle_position[0][1] << " " << current_triangle_position[0][2] << std::endl;
	//std::cout << initial_triangle_position[1][0] << " " << initial_triangle_position[1][1] << " " << initial_triangle_position[1][2] << std::endl;
	//std::cout << triangle_target_pos[1][0] << " " << triangle_target_pos[1][1] << " " << triangle_target_pos[1][2] << std::endl;
	//std::cout << current_triangle_position[1][0] << " " << current_triangle_position[1][1] << " " << current_triangle_position[1][2] << std::endl;
	//std::cout << initial_triangle_position[2][0] << " " << initial_triangle_position[2][1] << " " << initial_triangle_position[2][2] << std::endl;
	//std::cout << triangle_target_pos[2][0] << " " << triangle_target_pos[2][1] << " " << triangle_target_pos[2][2] << std::endl;
	//std::cout << current_triangle_position[2][0] << " " << current_triangle_position[2][1] << " " << current_triangle_position[2][2] << std::endl;
	//std::cout << "===" << std::endl;


	return true;
}

bool CollisionConstraint::pointColliderTriangle(double* initial_position, double* current_position,
	std::vector<double*>& current_triangle_position, double* triangle_normal, double* vertex_target_pos,
	double d_hat, double& stiffness)
{
	double d_2;
	double nearest_point_barycentric[3];
	double current_nearest_point[3];
	double d_hat_2 = d_hat * d_hat;
	if (!vertexTriangleDistance(initial_position, d_2, current_triangle_position, triangle_normal, nearest_point_barycentric)) {
		NearestPointInfo nearest_point_info;
		double barycentric_edge[2];
		vertexLineSegmentDistance(0, 1, initial_position, current_triangle_position, nearest_point_info, barycentric_edge);
		vertexLineSegmentDistance(1, 2, initial_position, current_triangle_position, nearest_point_info, barycentric_edge);
		vertexLineSegmentDistance(2, 0, initial_position, current_triangle_position, nearest_point_info, barycentric_edge);		

		//if (initial_position[1] < 0.769645) {
		//	//std::cout << "distance when <floor "<< initial_position[0]<<" "<< initial_position[1]<<" "<< initial_position[2]<<" "
		//		<< " " << current_triangle_position[0][0]<<" " << current_triangle_position[0][1] << " " << current_triangle_position[0][2] << " "
		//		<< current_triangle_position[1][0] << " " << current_triangle_position[1][1] << " " << current_triangle_position[1][2] << " "
		//		<< current_triangle_position[2][0] << " " << current_triangle_position[2][1] << " " << current_triangle_position[2][2] << " "
		//		<< " " << nearest_point_info.distance << " " << d_hat_2 << std::endl;
		//}
		if (nearest_point_info.distance >= d_hat_2) {
			return false;
		}
		d_2 = nearest_point_info.distance;
		calNearestPoint(nearest_point_info, barycentric_edge, current_triangle_position, current_nearest_point);
	}
	else {
		//if (initial_position[1] < 0.769645) {
		//	//std::cout << "distance when <floor" << initial_position[1] << " " << d_2 << " " << d_hat_2 << std::endl;
		//}
		if (d_2 >= d_hat_2) {
			return false;
		}
		calNearestPoint(nearest_point_barycentric, current_triangle_position[0], current_triangle_position[1], current_triangle_position[2], current_nearest_point);
	}

	stiffness *= barrier((d_hat_2- d_2)/d_hat_2, d_2 / d_hat_2);
	//std::cout << stiffness<<" "<< barrier(d_2 - d_hat_2, d_2 / d_hat_2) << std::endl;
	

	double d_move;// = sqrt(DOT(d, d));
	d_move = d_hat - sqrt(d_2);
	////////////////////////////////////////////
	//double d[3];
	//for (int i = 0; i < 3; ++i) {
	//	d[i] = initial_position[i] - current_position[i];
	//}
	//d_move = sqrt(DOT(d, d));
	//double d_min = d_hat - sqrt(d_2);
	//if (d_move < d_min) {
	//	d_move = d_min;
	//}
	///////////////////////////////////////////////
	//decide the direction
	double direction[3];
	SUB(direction, initial_position, current_nearest_point);
	double t1 = sqrt(DOT(direction, direction));
	if (t1 > NEAR_ZERO) {
		DEV_(direction, t1);
	}
	else {
		if (DOT(triangle_normal, direction) > 0) {
			memcpy(direction, triangle_normal, 24);
		}
		else {
			MULTI(direction, triangle_normal, -1.0);
		}
	}
	MULTI(vertex_target_pos, direction, d_move);
	SUM_(vertex_target_pos, initial_position);

	//if (vertex_target_pos[1] - current_triangle_position[0][1] < sqrt(d_hat_2)) {
	//	//std::cout << vertex_target_pos[1] - current_triangle_position[0][1] << std::endl;
	//}
	////std::cout << "build constraint " << std::endl;

	return true;
}


bool CollisionConstraint::edgeEdgeCollision(std::vector<std::array<double, 3>>& target_pos,
	std::vector<std::array<double, 3>>& compare_target_pos, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
	double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
	double* initial_compare_edge_vertex_1, double* mass, double d_hat, double& stiffness)
{
	double d_hat_2 = d_hat * d_hat;
	double alpha[4]; double distance;
	double initial_close_point_1[3]; double initial_close_point_2[3];
	if (!getClosestPoint(alpha, initial_edge_vertex_0, initial_edge_vertex_1,
		initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, d_hat_2, initial_close_point_1, initial_close_point_2, distance)) {
		return false;
	}
	stiffness *= barrier((d_hat_2- distance)/d_hat_2, distance / d_hat_2);
	double d[3];
	double current_close_point_1[3]; double current_close_point_2[3];
	POINT_ON_EDGE(current_close_point_1, alpha[0], alpha[1], current_edge_vertex_0, current_edge_vertex_1);
	POINT_ON_EDGE(current_close_point_2, alpha[2], alpha[3], current_compare_edge_vertex_0, current_compare_edge_vertex_1);

	////std::cout << "initial close point " << initial_close_point_1[0] << " " << initial_close_point_1[1] << " " << initial_close_point_1[2] << std::endl;
	////std::cout << "initial close point " << initial_close_point_2[0] << " " << initial_close_point_2[1] << " " << initial_close_point_2[2] << std::endl;
	////std::cout << "current close point " << current_close_point_1[0] << " " << current_close_point_1[1] << " " << current_close_point_1[2] << std::endl;
	////std::cout << "current close point " << current_close_point_2[0] << " " << current_close_point_2[1] << " " << current_close_point_2[2] << std::endl;

	for (int i = 0; i < 3; ++i) {
		d[i] = initial_close_point_2[i] - current_close_point_2[i] - (initial_close_point_1[i] - current_close_point_1[i]);
	}

	double d_move = sqrt(DOT(d, d));
	double d_min = d_hat + d_hat - sqrt(distance);
	if (d_move < d_min) {
		d_move = d_min;
	}
	//if (d_move < d_hat) {
	//	d_move = d_hat;
	//}

	double direction[3];
	SUB(direction, initial_close_point_2, initial_close_point_1);
	distance = sqrt(distance);
	if (distance > NEAR_ZERO) {
		DEV_(direction, distance);
	}
	else {
		double e1[3]; double e2[3];
		SUB(e1, initial_edge_vertex_1, initial_edge_vertex_0);
		SUB(e2, initial_compare_edge_vertex_1, initial_compare_edge_vertex_0);
		CROSS(direction, e1, e2);
		normalize(direction);
		SUB(e1, current_edge_vertex_0, current_compare_edge_vertex_0);
		if (DOT(e1, direction) > 0) {
			MULTI_(direction, -1.0);
		}
	}
	moveDistance(mass, target_pos, compare_target_pos, d_move, alpha, direction, initial_edge_vertex_0, initial_edge_vertex_1,
		initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
	return true;
}


void CollisionConstraint::testPT()
{
	double initial_pos[3] = { -1.01,0.001,1.01 };
	double current_pos[3] = { -1.01,-0.1,1.01 };
	std::vector<std::array<double, 3>> initial_triangle_pos(3);
	std::vector<std::array<double, 3>> current_triangle_pos(3);
	//initial_triangle_pos[0] = { 1.0,0.0,0.0 };
	//initial_triangle_pos[1] = { -1.0,0.0,-1.0 };
	//initial_triangle_pos[2] = { -1.0,0.0,1.0 };
	current_triangle_pos[0] = { 1.0,0.0,0.0 };
	current_triangle_pos[1] = { -1.0,0.0,-1.0 };
	current_triangle_pos[2] = { -1.0,0.0,1.0 };
	std::vector<double*> initial_tri(3);
	std::vector<double*> current_tri(3);
	for (int i = 0; i < 3; ++i) {
		initial_tri[i] = initial_triangle_pos[i].data();
		current_tri[i] = current_triangle_pos[i].data();
	}
	double initial_tri_normal[3]; double tri_normal[3];
	double e1[3], e2[3];
	SUB(e1, initial_triangle_pos[2], initial_triangle_pos[1]);
	SUB(e2, initial_triangle_pos[0], initial_triangle_pos[1]);
	CROSS(initial_tri_normal, e1, e2);
	normalize(initial_tri_normal);
	SUB(e1, current_triangle_pos[2], current_triangle_pos[1]);
	SUB(e2, current_triangle_pos[0], current_triangle_pos[1]);
	CROSS(tri_normal, e1, e2);
	normalize(tri_normal);

	double triangle_mass[3] = {1.0,1.0,1.0};
	double vertex_mass = 1.0;
	double d_hat_2 = 0.03 * 0.03;
	std::array<double, 3> vertex_target_pos;
	std::vector<std::array<double, 3>>triangle_target_pos(3);
	double stiffness;
	//if (pointSelfTriangle(initial_pos, current_pos, initial_tri, current_tri, initial_tri_normal, vertex_target_pos.data(),
	//	triangle_target_pos.data(), d_hat_2, stiffness, vertex_mass, triangle_mass)) {
	if (pointColliderTriangle(initial_pos, current_pos, current_tri, tri_normal, vertex_target_pos.data(),
		d_hat_2, stiffness)) {
		//std::cout <<"vertex "<< vertex_target_pos[0] << " " << vertex_target_pos[1] << " " << vertex_target_pos[2] << std::endl;
		////std::cout << "triangle vertex" << std::endl;
		//for (int i = 0; i < 3; ++i) {
		//	//std::cout << triangle_target_pos[i][0] << " " << triangle_target_pos[i][1] << " " << triangle_target_pos[i][2] << std::endl;
		//}
		//std::cout << stiffness << std::endl;
	}
	else {
		//std::cout << "does not collide " << std::endl;
	}
}

void CollisionConstraint::testEE()
{
	double initial_edge_0[3] = { 0.0,0.001,1.0 };
	double current_edge_0[3] = { 0.0,0.0,1.0 };
	double initial_edge_1[3] = { 0.0,0.001,0.0 };
	double current_edge_1[3] = { 0.0,0.0,0.0 };

	double initial_compare_edge_0[3] = { -0.01,0.0,-0.01 };
	double current_compare_edge_0[3] = { -0.01,0.3,-0.01 };
	double initial_compare_edge_1[3] = { -1.01,0.0,-0.01 };
	double current_compare_edge_1[3] = { -1.01,0.3,-0.01 };

	double mass[4] = { 0.5,1.5,0.5,1.5 };
	double d_hat_2 = 0.03 * 0.03;
	std::vector<std::array<double, 3>> target_pos(2);
	std::vector<std::array<double, 3>> compare_pos(2);
	double stiffness;
	if (edgeEdgeCollision(target_pos,compare_pos,current_edge_0, current_edge_1,initial_edge_0,initial_edge_1,
		current_compare_edge_0,current_compare_edge_1,initial_compare_edge_0,initial_compare_edge_1,mass,
		d_hat_2,stiffness)) {
		//std::cout << "edge 0 " << target_pos[0][0] << " " << target_pos[0][1] << " " << target_pos[0][2] << std::endl;
		//std::cout << "edge 1 " << target_pos[1][0] << " " << target_pos[1][1] << " " << target_pos[1][2] << std::endl;
		//std::cout << "compare 0 " << compare_pos[0][0] << " " << compare_pos[0][1] << " " << compare_pos[0][2] << std::endl;
		//std::cout << "compare 1 " << compare_pos[1][0] << " " << compare_pos[1][1] << " " << compare_pos[1][2] << std::endl;
		//std::cout << stiffness << std::endl;
		SUB_(target_pos[0], target_pos[1]);
		//std::cout << "dis " << DOT(target_pos[0], target_pos[0]) << std::endl;
	}
	else {
		//std::cout << "does not collide " << std::endl;
	}
}

bool CollisionConstraint::getClosestPoint(double* alpha, double* p1, double* p2, double* p3, double* p4, 
	double d_hat_2, double* c1, double* c2, double& distance)
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
	if (sub_ < NEAR_ZERO2) {   //parallel NEAR_ZERO2 is too small		
		if (getClosestDistanceCornerCornerCase(alpha, p1, p2, p3, p4, d_hat_2, distance)) {
			POINT_ON_EDGE(c1, alpha[0], alpha[1], p1, p2);
			POINT_ON_EDGE(c2, alpha[2], alpha[3], p3, p4);
			return true;
		}
	}
	else {
		double a0 = DOT(p3p1, s1);
		double a1 = DOT(p3p1, s2);
		alpha[1] = (s1s2 * a1 - s2_2 * a0) / sub_;
		alpha[3] = -(s1s2 * a0 - s1_2 * a1) / sub_;
		if (alpha[1] < 0 || alpha[1] > 1.0 || alpha[3] < 0 || alpha[3] > 1.0) {
			if (getClosestDistanceCornerCornerCase(alpha, p1, p2, p3, p4, d_hat_2, distance)) {
				POINT_ON_EDGE(c1, alpha[0], alpha[1], p1, p2);
				POINT_ON_EDGE(c2, alpha[2], alpha[3], p3, p4);
				return true;
			}
			return false;
		}
		alpha[0] = 1.0 - alpha[1];
		alpha[2] = 1.0 - alpha[3];
		POINT_ON_EDGE(c1, alpha[0], alpha[1], p1, p2);
		POINT_ON_EDGE(c2, alpha[2], alpha[3], p3, p4);
		SUB(s1,c1, c2);
		distance = DOT(s1, s1);
		if (distance<d_hat_2) {
			return true;
		}		
	}
	return false;
}


bool CollisionConstraint::getClosestDistanceCornerCornerCase(double* alpha, double* p1, double* p2, double* p3, double* p4, double d_hat_2, double& distance)
{
	double alpha_edge[2];
	distance = DBL_MAX;
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
	if (distance >= d_hat_2) {
		return false;
	}
	return true;
}

bool CollisionConstraint::getClosestPointBetweenPointSegement(double* alpha, double* p0, double* e0, double* e1, double& distance)
{
	double AP[3];
	double AB[3];
	SUB(AP, p0, e0);
	SUB(AB, e1, e0);
	double dot_ = DOT(AP, AB);
	double ab_2 = DOT(AB, AB);
	double r = dot_ / ab_2;
	if (r < 0) {
		double ap = DOT(AP, AP);
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
		bp = DOT(BP, BP);
		if (bp < distance) {
			alpha[0] = 0.0;
			alpha[1] = 1.0;
			distance = bp;
			return true;
		}
	}
	else {
		double ap_2 = DOT(AP, AP);
		double dis = ap_2 - ab_2 * r * r;
		if (dis < distance) {
			alpha[0] = 1.0 - r;
			alpha[1] = r;
			distance = dis;
			return true;
		}
	}
	return false;
}


void CollisionConstraint::moveDistance(double* mass, std::vector<std::array<double, 3>>& target_pos,
	std::vector<std::array<double, 3>>& compare_target_pos,	double d_move, double* barycentric, double* direction, double* initial_edge_vertex_0, 
	double* initial_edge_vertex_1, double* initial_compare_edge_vertex_0,double* initial_compare_edge_vertex_1)
{

	double mass_total = mass[0] + mass[1] + mass[2]+ mass[3];
	double d_close[2];
	d_close[1] = -mass[0]*(mass[2]+mass[3])* d_move / (mass_total*(barycentric[0]*mass[1]+barycentric[1]*mass[0]));
	d_close[0] = mass[1] * d_close[1] / mass[0];
	MULTI(target_pos[0], direction, d_close[0]);
	SUM_(target_pos[0], initial_edge_vertex_0);
	MULTI(target_pos[1], direction, d_close[1]);
	SUM_(target_pos[1], initial_edge_vertex_1);

	////std::cout << "moving distance " << d_close[0] << " " << d_close[1]<<" ";
	d_close[1] = mass[2] * (mass[0] + mass[1]) * d_move / (mass_total * (barycentric[2] * mass[3] + barycentric[3] * mass[2]));
	d_close[0] = mass[3] * d_close[1] / mass[2];
	MULTI(compare_target_pos[0], direction, d_close[0]);
	SUM_(compare_target_pos[0], initial_compare_edge_vertex_0);
	MULTI(compare_target_pos[1], direction, d_close[1]);
	SUM_(compare_target_pos[1], initial_compare_edge_vertex_1);
	////std::cout << d_close[0] << " " << d_close[1] << std::endl;

}


void CollisionConstraint::moveDistance(double vertex_mass, double* triangle_mass,double* vertex_target_pos, std::array<double, 3>* triangle_target_pos, 
	double d_move, double* barycentric, double* direction, double* initial_pos, std::vector<double*>& initial_triangle_position)
{
	double vertex_d;
	double triangle_d[3];
	double mass_triangle = triangle_mass[0] + triangle_mass[1] + triangle_mass[2];
	double d_close = d_move / (vertex_mass + mass_triangle);
	vertex_d = -mass_triangle * d_close;
	d_close *= vertex_mass;
	//std::cout << triangle_mass[0] << " " << triangle_mass[1] << " " << triangle_mass[2] << " " << barycentric[0]<<" "<< barycentric[1]<<" "<< barycentric[2] << std::endl;

	triangle_d[1] = triangle_mass[0] * triangle_mass[2] * d_close / (barycentric[0] * triangle_mass[1] * triangle_mass[2] + barycentric[1] * triangle_mass[0] * triangle_mass[2]
		+ barycentric[2] * triangle_mass[0] * triangle_mass[1]);
	triangle_d[0] = triangle_mass[1] * triangle_d[1] / triangle_mass[0];
	triangle_d[2] = triangle_mass[1] * triangle_d[1] / triangle_mass[2];
	MULTI(vertex_target_pos, direction, vertex_d);
	SUM_(vertex_target_pos, initial_pos);
	
	for (int i = 0; i < 3; ++i) {
		MULTI(triangle_target_pos[i], direction, triangle_d[i]);
		SUM_(triangle_target_pos[i], initial_triangle_position[i]);
		////std::cout << "direction " << triangle_d[i] << " " << triangle_d[i] << " " << triangle_d[i] << std::endl;
		////std::cout << "direction " << triangle_d[i][0] << " " << triangle_d[i][1] << " " << triangle_d[i][2] << std::endl;
	}
	
}


bool CollisionConstraint::vertexTriangleDistance(double* vertex, double& d_2, std::vector<double*>& triangle_position, double* triangle_normal, double* barycentric)
{
	double S[3];
	SUB(S, vertex, triangle_position[0]);
	double E1[3], E2[3], S1[3], S2[3];
	SUB(E1, triangle_position[1], triangle_position[0]);
	SUB(E2, triangle_position[2], triangle_position[0]);
	CROSS(S1, triangle_normal, E2);
	CROSS(S2, S, E1);
	double temp = 1.0 / DOT(S1, E1);
	d_2 = temp * DOT(S2, E2);
	d_2 *= d_2;
	barycentric[1] = temp * DOT(S1, S);
	barycentric[2] = temp * DOT(S2, triangle_normal);
	barycentric[0] = 1.0 - barycentric[1] - barycentric[2];
	//if (vertex[1] < 0.769645) {
	//	//std::cout << "barycentric "<< barycentric[0]<<" "<< barycentric[1]<<" "<< barycentric[2] << std::endl;
	//}

	if (barycentric[0] > EPSILON && barycentric[1] > EPSILON && barycentric[2] > EPSILON) {
		return true;
	}
	return false;
}

void CollisionConstraint::calNearestPoint(double* barycentric, double* v0, double* v1, double* v2, double* center)
{
	BARYCENTRIC(center, barycentric, v0, v1, v2);
}

void CollisionConstraint::vertexLineSegmentDistance(int edge_vertex_0, int edge_vertex_1, double* vertex, std::vector<double*>& triangle_position,
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
		ap = DOT(AP, AP);
		if (ap < nearest_point_info.distance) {
			nearest_point_info.is_edge = false;
			nearest_point_info.index[0] = edge_vertex_0;
			nearest_point_info.distance = ap;
		}
	}
	else if (r > 1) {
		double BP[3], bp;
		SUB(BP, vertex, triangle_position[edge_vertex_1]);
		bp = DOT(BP, BP);
		if (bp < nearest_point_info.distance) {
			nearest_point_info.is_edge = false;
			nearest_point_info.index[0] = edge_vertex_1;
			nearest_point_info.distance = bp;
		}
	}
	else {
		double ap_2 = DOT(AP, AP);
		double dis = ap_2 - ab_2 * r * r;
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


void CollisionConstraint::calNearestPoint(NearestPointInfo& nearest_point_info, double* barycentric, std::vector<double*>& triangle_position, double* center)
{
	if (nearest_point_info.is_edge) {
		pointOnEdge(center, barycentric[0], barycentric[1], triangle_position[nearest_point_info.index[0]],
			triangle_position[nearest_point_info.index[1]]);
	}
	else {
		memcpy(center, triangle_position[nearest_point_info.index[0]], 24);
	}
}


void CollisionConstraint::setBarycentric(NearestPointInfo& nearest_point_info, double* barycentric, double* edge_barycentric)
{
	memset(barycentric, 0, 24);
	if (nearest_point_info.is_edge) {
		barycentric[nearest_point_info.index[0]] = edge_barycentric[0];
		barycentric[nearest_point_info.index[1]] = edge_barycentric[1];
	}
	else {
		barycentric[nearest_point_info.index[0]] = 1;
	}
}


void CollisionConstraint::pointOnEdge(double* center, double coe0, double coe1, double* vertex0, double* vertex1)
{
	POINT_ON_EDGE(center, coe0, coe1, vertex0, vertex1);
}

inline double CollisionConstraint::barrier(double d_2_minus_d_hat_2_over_d_hat2, double d_2_div_d_hat_2) 
{
	return -d_2_minus_d_hat_2_over_d_hat2 * d_2_minus_d_hat_2_over_d_hat2 * log(d_2_div_d_hat_2);
}