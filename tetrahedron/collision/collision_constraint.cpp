#include"collision_constraint.h"
#undef PROJECT_VELOCITY_COE_ON_NORMAL
#define PROJECT_VELOCITY_COE_ON_NORMAL(dest, p0, p1, normal)\
dest=(p0[0]-p1[0])*normal[0]+(p0[1]-p1[1])*normal[1]+(p0[2]-p1[2])*normal[2];

#undef SUBTRACT_A_WITH_COE_B
#define SUBTRACT_A_WITH_COE_B(dest, A, B, coe)\
dest[0]=A[0]-coe*B[0];\
dest[1]=A[1]-coe*B[1];\
dest[2]=A[2]-coe*B[2];



bool CollisionConstraint::pointTriangleResponse(double* initial_position, double* current_position,
	double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
	double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
	double* previous_free_vertex_position, double* previous_free_triangle_position_0,
	double* previous_free_triangle_position_1, double* previous_free_triangle_position_2,
	double* initial_triangle_normal, double* vertex_target_pos,
	double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
	double d_hat, double& stiffness, double epsilon, double mass_point, double mass_t0, double mass_t1, double mass_t2,unsigned int* obj_index,
	double collision_time)
{
	double d_hat_2 = d_hat * d_hat;
	double barycentric[3];
	double d_2 = CCD::internal::pointTriangleNearestDistance(initial_position, initial_triangle_position_0, initial_triangle_position_1,
		initial_triangle_position_2, initial_triangle_normal, barycentric,d_hat_2);

	if (d_2 == 0.0) {
		d_2 += 1e-6;
		//std::cout <<  "yy "<< obj_index[0]<<" "<<obj_index[1]<<" "<< obj_index[2]<<" "<< obj_index[3] << std::endl;
		//std::cout << initial_position[0] << " " << initial_position[1] << " " << initial_position[2] << std::endl;
		//std::cout << initial_triangle_position_0[0] << " " << initial_triangle_position_0[1] << " " << initial_triangle_position_0[2] << std::endl;
		//std::cout << initial_triangle_position_1[0] << " " << initial_triangle_position_1[1] << " " << initial_triangle_position_1[2] << std::endl;
		//std::cout << initial_triangle_position_2[0] << " " << initial_triangle_position_2[1] << " " << initial_triangle_position_2[2] << std::endl;
	}


	if (d_2 >= d_hat_2) {
		return false;
	}

	std::cout << "collision time " << collision_time << std::endl;

	////std::cout << d_2<<" "<< barycentric[0]<<" "<< barycentric[1]<<" "<< barycentric[2] << std::endl;
	stiffness *= barrier((d_hat_2 - d_2) / d_hat_2, d_2 / d_hat_2);


	double collision_vertex[3];
	double collision_tri_0[3]; double collision_tri_1[3]; double collision_tri_2[3];

	COLLISION_POS(collision_vertex, collision_time, previous_free_vertex_position, current_position);
	COLLISION_POS(collision_tri_0, collision_time, previous_free_triangle_position_0, current_triangle_position_0);
	COLLISION_POS(collision_tri_1, collision_time, previous_free_triangle_position_1, current_triangle_position_1);
	COLLISION_POS(collision_tri_2, collision_time, previous_free_triangle_position_2, current_triangle_position_2);

	double normal[3];
	getTriangleNormal(collision_tri_0, collision_tri_1, collision_tri_2, normal);
	CCD::internal::pointTriangleNearestPoint(collision_vertex, collision_tri_0, collision_tri_1, collision_tri_2, normal, barycentric);

	double initial_nearest_point[3];
	double current_nearest_point[3];
	BARYCENTRIC(initial_nearest_point, barycentric, collision_tri_0, collision_tri_1, collision_tri_2);
	BARYCENTRIC(current_nearest_point, barycentric, current_triangle_position_0, current_triangle_position_1, current_triangle_position_2);

	////std::cout << "barycentric " << barycentric[0] << " " << barycentric[1] << " " << barycentric[2] << std::endl;
	////std::cout << "move p " << initial_nearest_point[0] << " " << initial_nearest_point[1]<<" "<< initial_nearest_point[2] << std::endl;


	double relative_velocity[3];
	relative_velocity[0] = current_position[0] - collision_vertex[0] - current_nearest_point[0] + initial_nearest_point[0];
	relative_velocity[1] = current_position[1] - collision_vertex[1] - current_nearest_point[1] + initial_nearest_point[1];
	relative_velocity[2] = current_position[2] - collision_vertex[2] - current_nearest_point[2] + initial_nearest_point[2];
	
	double total_distance;
	epsilon += 1.0;

	total_distance = epsilon * DOT(relative_velocity, normal);
	////std::cout << "epsilon " << epsilon << " " << total_distance << std::endl;

	double sub[3];
	SUB(sub, initial_position, initial_triangle_position_0);

	double sideness = DOT(initial_triangle_normal, sub);


	SUB(sub, collision_vertex, initial_nearest_point);
	double collision_dis = DOT(sub, normal);

	if (sideness == 0.0) {
		double direct_velocity = DOT(initial_triangle_normal, relative_velocity);
		if (direct_velocity > 0) {
			total_distance = abs(total_distance);
			if (total_distance < d_hat) {
				total_distance = d_hat;
			}
		}
		else {
			total_distance = -abs(total_distance);
			if (total_distance > -d_hat) {
				total_distance = -d_hat;
			}
		}
	}
	else if(sideness < 0) {
		if (collision_dis < -d_hat) {
			memcpy(vertex_target_pos, collision_vertex, 24);
			memcpy(triangle_target_pos_0, collision_tri_0, 24);
			memcpy(triangle_target_pos_1, collision_tri_1, 24);
			memcpy(triangle_target_pos_2, collision_tri_2, 24);
			return true;
		}
		if (total_distance  - collision_dis < d_hat) {
			total_distance = d_hat + collision_dis;
		}
	}
	else {
		if (collision_dis > d_hat) {
			memcpy(vertex_target_pos, collision_vertex, 24);
			memcpy(triangle_target_pos_0, collision_tri_0, 24);
			memcpy(triangle_target_pos_1, collision_tri_1, 24);
			memcpy(triangle_target_pos_2, collision_tri_2, 24);
			return true;
		}
		if (collision_dis - total_distance < d_hat) {
			total_distance = -d_hat + collision_dis;
		}
	}



	double mass_nearest_point = mass_t0 + mass_t1 + mass_t2;
	double move_d_point = mass_nearest_point * total_distance / (mass_point + mass_nearest_point);
	double move_d_t1 = - mass_t0 * mass_t2 * mass_point * total_distance /
		((mass_point + mass_nearest_point)*(barycentric[0]*mass_t1*mass_t2 + barycentric[1]*mass_t0*mass_t2
			+barycentric[2]*mass_t0*mass_t1));

	

	double move_d_t0 = mass_t1 / mass_t0 * move_d_t1;
	double move_d_t2= mass_t1 / mass_t2 * move_d_t1;

	SUBTRACT_A_WITH_COE_B(vertex_target_pos, collision_vertex, normal, move_d_point);
	SUBTRACT_A_WITH_COE_B(triangle_target_pos_0, collision_tri_0, normal, move_d_t0);
	SUBTRACT_A_WITH_COE_B(triangle_target_pos_1, collision_tri_1, normal, move_d_t1);
	SUBTRACT_A_WITH_COE_B(triangle_target_pos_2, collision_tri_2, normal, move_d_t2);


	double test[3];
	BARYCENTRIC(test, barycentric, triangle_target_pos_0, triangle_target_pos_1, triangle_target_pos_2);
	SUB_(test, vertex_target_pos);
	std::cout << "dis " << sqrt(DOT(test, test)) << " " << d_hat<<" "<<stiffness << std::endl;

	
	//double e1[3], e2[3];
	//SUB(e1, triangle_target_pos_2, triangle_target_pos_1);
	//SUB(e2, triangle_target_pos_0, triangle_target_pos_1);
	//double norm[3];
	//CROSS(norm, e1, e2);
	//normalize(norm);
	//double d_3 = CCD::internal::pointTriangleNearestDistance(vertex_target_pos, triangle_target_pos_0, triangle_target_pos_1,
	//	triangle_target_pos_2, norm, barycentric, d_hat_2);

	//std::cout << "d compare " << d_2 << " " << d_3<<" "<<d_hat_2 << std::endl;

	//PROJECT_VELOCITY_COE_ON_NORMAL(coe, current_position, initial_position, initial_triangle_normal);
	

	////std::cout <<"coe "<< coe << std::endl;
	////std::cout <<"coe "<< initial_triangle_normal[0]<<" "<< initial_triangle_normal[1]<<" "<<initial_triangle_normal[2] << std::endl;

	//SUBTRACT_A_WITH_COE_B(vertex_target_pos, initial_position, initial_triangle_normal, coe);	
	//PROJECT_VELOCITY_COE_ON_NORMAL(coe, current_triangle_position_0, initial_triangle_position_0, initial_triangle_normal);
	//coe *= epsilon;
	//SUBTRACT_A_WITH_COE_B(triangle_target_pos_0, initial_triangle_position_0, initial_triangle_normal, coe);
	//PROJECT_VELOCITY_COE_ON_NORMAL(coe, current_triangle_position_1, initial_triangle_position_1, initial_triangle_normal);
	//coe *= epsilon;
	//SUBTRACT_A_WITH_COE_B(triangle_target_pos_1, initial_triangle_position_1, initial_triangle_normal, coe);
	//PROJECT_VELOCITY_COE_ON_NORMAL(coe, current_triangle_position_2, initial_triangle_position_2, initial_triangle_normal);
	//coe *= epsilon;
	//SUBTRACT_A_WITH_COE_B(triangle_target_pos_2, initial_triangle_position_2, initial_triangle_normal, coe);

	return true;
}


bool CollisionConstraint::floorResponse(double* target_position, double* current_position, double* initial_position,
	unsigned int dimension, bool normal_direction,
	double floor_value, double d_hat, double& stiffness, double epsilon, double* previous_free_pos)
{
	double d_2;
	d_2 = (initial_position[dimension] - floor_value) * (initial_position[dimension] - floor_value);
	if (normal_direction) {
		if (initial_position[dimension] > floor_value + d_hat) {
			return false;
		}
	}
	else {
		if (initial_position[dimension] < floor_value - d_hat) {
			return false;
		}
	}

	double d_hat_2 = d_hat * d_hat;



	stiffness *= barrier((d_hat_2 - d_2) / d_hat_2, d_2 / d_hat_2);

	epsilon += 1.0;


	double collision_pos[3];



	double collision_time;
	if (abs(current_position[dimension] - initial_position[dimension]) < 1e-10) {
		collision_time = 1.0;
	}
	else {
		collision_time = (floor_value - initial_position[dimension]) / (current_position[dimension] - initial_position[dimension]);
	}	

	if (collision_time > 1.0) {
		collision_time = 1.0;
	}

	COLLISION_POS(collision_pos, collision_time, previous_free_pos, current_position);

	memcpy(target_position, collision_pos, 24);

	if (normal_direction) {
		if (floor_value-current_position[dimension] < d_hat) {
			target_position[dimension] = floor_value + d_hat;
		}
		else if (current_position[dimension] > d_hat + floor_value) {
			memcpy(target_position, current_position, 24);
		}
		else {
			target_position[dimension] += floor_value - current_position[dimension];
		}
	}
	else {
		if (current_position[dimension]- floor_value < d_hat) {
			target_position[dimension] = floor_value - d_hat;
		}
		else  if (current_position[dimension] < d_hat - floor_value) {
			memcpy(target_position, current_position, 24);
		}
		else {
			target_position[dimension] -= current_position[dimension]- floor_value;
		}
	}
	return true;
}




bool CollisionConstraint::pointTriangleColliderResponse(double* initial_position, double* current_position,
	double* previous_free_vertex_position,
	double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
	double* initial_triangle_normal, double* vertex_target_pos,
	double d_hat, double& stiffness, double epsilon, unsigned int vertex_index, double collision_time)
{
	double d_hat_2 = d_hat * d_hat;
	double barycentric[3];

	double d_2 = CCD::internal::pointTriangleNearestDistance(initial_position, initial_triangle_position_0, initial_triangle_position_1,
		initial_triangle_position_2, initial_triangle_normal, barycentric,d_hat_2);
	if (d_2 >= d_hat_2) {
		return false;
	}
	////std::cout << d_2 << std::endl;

	stiffness *= barrier((d_hat_2 - d_2) / d_hat_2, d_2 / d_hat_2);

	epsilon += 1.0;

	double collision_vertex[3];

	COLLISION_POS(collision_vertex, collision_time, previous_free_vertex_position, current_position);


	double coe;
	PROJECT_VELOCITY_COE_ON_NORMAL(coe, current_position, collision_vertex, initial_triangle_normal);
	coe *= epsilon;

	double sub[3];
	SUB(sub, initial_position, initial_triangle_position_0);
	double relative_velocity[3];
	SUB(relative_velocity, current_position, collision_vertex);

	double sideness = DOT(initial_triangle_normal, sub);

	SUB(sub, collision_vertex, initial_triangle_position_0);
	double collision_dis = DOT(initial_triangle_normal, sub);
	//double indicate_move=0.0;


	if (sideness == 0.0) {
		double direct_velocity = DOT(initial_triangle_normal, relative_velocity);
		if (direct_velocity > 0) {
			coe = abs(coe);
			if (coe < d_hat) {
				coe = d_hat;
			}
		}
		else {
			coe = -abs(coe);
			if (coe > -d_hat) {
				coe = -d_hat;
			}
		}
	}
	else if (sideness < 0) {
		if (collision_dis < -d_hat) {
			memcpy(vertex_target_pos, collision_vertex, 24);
			return true;
		}
		
		if (coe - collision_dis < d_hat) {
			coe = d_hat + collision_dis;
		}
	}
	else {
		if (collision_dis > d_hat) {
			memcpy(vertex_target_pos, collision_vertex, 24);
			return true;
		}
		if (collision_dis -coe < d_hat) {
			coe = -d_hat + collision_dis;
		}
	}

	//if (vertex_index == 350) {
	//	std::cout << coe << " " << indicate_move <<" "<< epsilon << std::endl;
	//}
	//if (coe == 0.0) {
	//	coe = sqrt(d_hat_2);
	//}
	//double direct_velocity = DOT(initial_triangle_normal, relative_velocity);
	//if (sideness * direct_velocity > 0) {
	//	coe *= -1.0;
	//}
	////std::cout <<"coe "<< coe << std::endl;
	////std::cout <<"coe "<< initial_triangle_normal[0]<<" "<< initial_triangle_normal[1]<<" "<<initial_triangle_normal[2] << std::endl;

	SUBTRACT_A_WITH_COE_B(vertex_target_pos, collision_vertex, initial_triangle_normal, coe);
	return true;
}


bool CollisionConstraint::pointColliderTriangleResponse(double* initial_position,
	double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
	double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
	double* initial_triangle_normal,
	double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
	double d_hat, double& stiffness, double epsilon, double mass_t0, double mass_t1, double mass_t2)
{
	double d_hat_2 = d_hat * d_hat;
	double barycentric[3];
	double d_2 = CCD::internal::pointTriangleNearestDistance(initial_position, initial_triangle_position_0, initial_triangle_position_1,
		initial_triangle_position_2, initial_triangle_normal, barycentric, d_hat_2);
	if (d_2 >= d_hat_2) {
		return false;
	}

	////std::cout << d_2 << std::endl;

	stiffness *= barrier((d_hat_2 - d_2) / d_hat_2, d_2 / d_hat_2);

	double initial_nearest_point[3];
	double current_nearest_point[3];
	BARYCENTRIC(initial_nearest_point, barycentric, initial_triangle_position_0, initial_triangle_position_1, initial_triangle_position_2);
	BARYCENTRIC(current_nearest_point, barycentric, current_triangle_position_0, current_triangle_position_1, current_triangle_position_2);

	double relative_velocity[3];
	SUB(relative_velocity, current_nearest_point, initial_nearest_point);

	double total_distance;
	epsilon += 1.0;
	total_distance = epsilon * DOT(relative_velocity, initial_triangle_normal);

	double sub[3];
	SUB(sub, initial_position, initial_triangle_position_0);
	double sideness = DOT(initial_triangle_normal, sub);
	//double direct_velocity = DOT(initial_triangle_normal, relative_velocity);

	if (sideness == 0.0) {
		double direct_velocity = DOT(initial_triangle_normal, relative_velocity);
		if (direct_velocity > 0) {
			total_distance = abs(total_distance);
			if (total_distance < d_hat) {
				total_distance = d_hat;
			}
		}
		else {
			total_distance = -abs(total_distance);
			if (total_distance > -d_hat) {
				total_distance = -d_hat;
			}
		}
	}
	else if (sideness < 0) {
		total_distance = abs(total_distance);
		if (total_distance < d_hat + sideness) {
			total_distance = d_hat + sideness;
		}

	}
	else {
		total_distance = -abs(total_distance);
		if (total_distance > -d_hat + sideness) {
			total_distance = -d_hat + sideness;
		}
	}

	////std::cout << total_distance << std::endl;

	//if (sideness * direct_velocity > 0) {
	//	total_distance *= -1.0;
	//}
	//if (total_distance == 0.0) {
	//	total_distance = sqrt(d_hat_2);
	//}

	double move_d_t1 = mass_t0 * mass_t2 * total_distance /
		(barycentric[0] * mass_t1 * mass_t2 + barycentric[1] * mass_t0 * mass_t2
			+ barycentric[2] * mass_t0 * mass_t1);
	double move_d_t0 = mass_t1 / mass_t0 * move_d_t1;
	double move_d_t2 = mass_t1 / mass_t2 * move_d_t1;

	SUBTRACT_A_WITH_COE_B(triangle_target_pos_0, initial_triangle_position_0, initial_triangle_normal, move_d_t0);

	SUBTRACT_A_WITH_COE_B(triangle_target_pos_1, initial_triangle_position_1, initial_triangle_normal, move_d_t1);

	SUBTRACT_A_WITH_COE_B(triangle_target_pos_2, initial_triangle_position_2, initial_triangle_normal, move_d_t2);

	return true;
}


bool CollisionConstraint::edgeEdgeColliderResponse(double* edge_target_pos_0, double* edge_target_pos_1,
	double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
	double* initial_compare_edge_vertex_0, double* initial_compare_edge_vertex_1,
	double d_hat, double& stiffness, double epsilon, double mass_e_0_0, double mass_e_0_1)
{
	double d_hat_2 = d_hat * d_hat;
	double bary_centric[4];
	double d_2 = CCD::internal::edgeEdgeNearestPoint(initial_edge_vertex_0, initial_edge_vertex_1, initial_compare_edge_vertex_0,
		initial_compare_edge_vertex_1, bary_centric);
	if (d_2 >= d_hat_2) {
		return false;
	}
	////std::cout << d_2 << std::endl;
	stiffness *= barrier((d_hat_2 - d_2) / d_hat_2, d_2 / d_hat_2);
	double coe;

	double relative_velocity[3];
	for (unsigned int i = 0; i < 3; ++i) {
		relative_velocity[i] = bary_centric[0] * (current_edge_vertex_0[i] - initial_edge_vertex_0[i])
			+ bary_centric[1] * (current_edge_vertex_1[i] - initial_edge_vertex_1[i]);
	}

	double collision_normal[3];
	CROSS_FOUR_POINTS(collision_normal, initial_edge_vertex_0, initial_edge_vertex_1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
	if (DOT(collision_normal, collision_normal) < 1e-30) {
		double nearest_vec[3];
		for (unsigned int i = 0; i < 3; ++i) {
			nearest_vec[i] = bary_centric[0] * initial_edge_vertex_0[i] + bary_centric[1] * initial_edge_vertex_1[i]
				- bary_centric[2] * initial_compare_edge_vertex_0[i] - bary_centric[3] * initial_compare_edge_vertex_1[i];
		}
		if (DOT(nearest_vec, nearest_vec) < 1e-30) {
			memcpy(collision_normal, relative_velocity, 24);
		}
		else {
			memcpy(collision_normal, nearest_vec, 24);
		}
	}
	normalize(collision_normal);
	epsilon += 1.0;

	double total_distance = epsilon * DOT(relative_velocity, collision_normal);


	double sub[3];
	SUB(sub, initial_edge_vertex_0, initial_compare_edge_vertex_0);
	double sideness = DOT(collision_normal, sub);
	//double direct_velocity = DOT(collision_normal, relative_velocity);

	if (sideness == 0.0) {
		double direct_velocity = DOT(collision_normal, relative_velocity);
		if (direct_velocity > 0) {
			total_distance = abs(total_distance);
			if (total_distance < d_hat) {
				total_distance = d_hat;
			}
		}
		else {
			total_distance = -abs(total_distance);
			if (total_distance > -d_hat) {
				total_distance = -d_hat;
			}
		}
	}
	else if (sideness < 0) {
		total_distance = abs(total_distance);
		if (total_distance < d_hat + sideness) {
			total_distance = d_hat + sideness;
		}
	}
	else {
		total_distance = -abs(total_distance);
		if (total_distance > -d_hat + sideness) {
			total_distance = -d_hat + sideness;
		}
	}


	//if (sideness * direct_velocity > 0) {
	//	total_distance *= -1.0;
	//}



	double d_0_1 = (mass_e_0_0 / (mass_e_0_1 * bary_centric[0] + mass_e_0_0 * bary_centric[1])) * total_distance;
	double d_0_0 = (mass_e_0_1 / mass_e_0_0) * d_0_1;
	SUBTRACT_A_WITH_COE_B(edge_target_pos_0, initial_edge_vertex_0, collision_normal, d_0_0);
	SUBTRACT_A_WITH_COE_B(edge_target_pos_1, initial_edge_vertex_1, collision_normal, d_0_1);

	return true;
}

bool CollisionConstraint::edgeEdgeResponse(double* edge_target_pos_0, double* edge_target_pos_1,
	double* compare_target_pos_0, double* compare_target_pos_1,
	double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
	double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
	double* initial_compare_edge_vertex_1, 
	double* previous_free_edge_v0, double* previous_free_edge_v1,
	double* previous_free_compare_edge_v0, double* previous_free_compare_edge_v1,
	double d_hat, double& stiffness, double epsilon,
	double mass_e_0_0, double mass_e_0_1,double mass_e_1_0,double mass_e_1_1, double collision_time)
{
	//std::cout << collision_time << std::endl;

	double d_hat_2 = d_hat * d_hat;

	double bary_centric[4];
	double d_2 = CCD::internal::edgeEdgeNearestPoint(initial_edge_vertex_0, initial_edge_vertex_1, initial_compare_edge_vertex_0,
		initial_compare_edge_vertex_1, bary_centric);
	if (d_2 >= d_hat_2) {
		return false;
	}
	////std::cout << d_2 << std::endl;
	stiffness *= barrier((d_hat_2 - d_2) / d_hat_2, d_2 / d_hat_2);

	double collision_normal[3];
	epsilon += 1.0;

	double collision_edge_v0[3];	double collision_edge_v1[3];	double collision_edge_compare_v0[3]; double collision_edge_compare_v1[3];

	COLLISION_POS(collision_edge_v0, collision_time, previous_free_edge_v0, current_edge_vertex_0);
	COLLISION_POS(collision_edge_v1, collision_time, previous_free_edge_v1, current_edge_vertex_1);
	COLLISION_POS(collision_edge_compare_v0, collision_time, previous_free_compare_edge_v0, current_compare_edge_vertex_0);
	COLLISION_POS(collision_edge_compare_v1, collision_time, previous_free_compare_edge_v1, current_compare_edge_vertex_1);

	CCD::internal::edgeEdgeDistanceType(previous_free_edge_v0, previous_free_edge_v1, previous_free_compare_edge_v0, previous_free_compare_edge_v1, bary_centric);


	double relative_velocity[3];
	for (unsigned int i = 0; i < 3; ++i) {
		relative_velocity[i] = bary_centric[0] * (current_edge_vertex_0[i] - collision_edge_v0[i])
			+ bary_centric[1] * (current_edge_vertex_1[i] - collision_edge_v1[i])
			- bary_centric[2] * (current_compare_edge_vertex_0[i] - collision_edge_compare_v0[i])
			- bary_centric[3] * (current_compare_edge_vertex_1[i] - collision_edge_compare_v1[i]);
	}

	EDGE_OF_TWO_EDGE(collision_normal, bary_centric, previous_free_edge_v0, previous_free_edge_v1, previous_free_compare_edge_v0, previous_free_compare_edge_v1);

	//if (DOT(collision_normal, collision_normal) < 1e-30) {
	//	double nearest_vec[3];
	//	for (unsigned int i = 0; i < 3; ++i) {
	//		nearest_vec[i] = bary_centric[0] * collision_edge_v0[i] + bary_centric[1] * collision_edge_v1[i]
	//			- bary_centric[2] * collision_edge_compare_v0[i] - bary_centric[3] * collision_edge_compare_v1[i];
	//	}
	//	if (DOT(nearest_vec, nearest_vec) < 1e-30) {
	//		memcpy(collision_normal, relative_velocity, 24);
	//	}
	//	else {
	//		memcpy(collision_normal, nearest_vec, 24);
	//	}
	//}
	normalize(collision_normal);


	double total_distance = epsilon * DOT(relative_velocity, collision_normal);

	double current_edge[3];
	EDGE_OF_TWO_EDGE(current_edge, bary_centric, current_edge_vertex_0, current_edge_vertex_1, current_compare_edge_vertex_0, current_compare_edge_vertex_1);
	double current_distance = DOT(current_edge, collision_normal);
	if (current_distance > d_hat) {
		memcpy(edge_target_pos_0, current_edge_vertex_0, 24);
		memcpy(edge_target_pos_1, current_edge_vertex_1, 24);
		memcpy(compare_target_pos_0, current_compare_edge_vertex_0, 24);
		memcpy(compare_target_pos_1, current_compare_edge_vertex_1, 24);
		return true;
	}
	if (total_distance > -d_hat) {
		total_distance = -d_hat;
	}

	double move_d_0_1 = mass_e_0_0 * (mass_e_1_0 + mass_e_1_1) * total_distance
		/ ((mass_e_0_0 + mass_e_0_1 + mass_e_1_0 + mass_e_1_1) * (bary_centric[0] * mass_e_0_1 + bary_centric[1] * mass_e_0_0));
	double move_d_0_0 = (mass_e_0_1 / mass_e_0_0) * move_d_0_1;
	double move_d_1_1 = -mass_e_1_0 * (mass_e_0_0 + mass_e_0_1) * total_distance
		/ ((mass_e_0_0 + mass_e_0_1 + mass_e_1_0 + mass_e_1_1) * (bary_centric[2] * mass_e_1_1 + bary_centric[3] * mass_e_1_0));
	double move_d_1_0 = (mass_e_1_1 / mass_e_1_0) * move_d_1_1;
	////std::cout <<"coe "<< coe << std::endl;
	////std::cout <<"coe "<< collision_normal[0]<<" "<< collision_normal[1]<<" "<<collision_normal[2] << std::endl;
	SUBTRACT_A_WITH_COE_B(edge_target_pos_0, initial_edge_vertex_0, collision_normal, move_d_0_0);
	SUBTRACT_A_WITH_COE_B(edge_target_pos_1, initial_edge_vertex_1, collision_normal, move_d_0_1);
	SUBTRACT_A_WITH_COE_B(compare_target_pos_0, initial_compare_edge_vertex_0, collision_normal, move_d_1_0);
	SUBTRACT_A_WITH_COE_B(compare_target_pos_1, initial_compare_edge_vertex_1, collision_normal, move_d_1_1);
	return true;
}




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
	double d_move =4.0 * sqrt(DOT(d, d));
	double d_min= d_hat - sqrt(d_2);
	if (d_move < d_min) {
		d_move =4.0 * d_min;
	}
	//else {
	//	//std::cout << d_move<<" "<< d_min << std::endl;
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
	////std::cout << d_move<<" "<<stiffness << std::endl;
	moveDistance(vertex_mass, triangle_mass, vertex_target_pos, triangle_target_pos,d_move, nearest_point_barycentric, direction,initial_position, initial_triangle_position);
	////std::cout << initial_triangle_position[0][0] << " " << initial_triangle_position[0][1] << " " << initial_triangle_position[0][2] << std::endl;
	////std::cout << triangle_target_pos[0][0] << " " << triangle_target_pos[0][1] << " " << triangle_target_pos[0][2] << std::endl;
	////std::cout << current_triangle_position[0][0] << " " << current_triangle_position[0][1] << " " << current_triangle_position[0][2] << std::endl;
	////std::cout << initial_triangle_position[1][0] << " " << initial_triangle_position[1][1] << " " << initial_triangle_position[1][2] << std::endl;
	////std::cout << triangle_target_pos[1][0] << " " << triangle_target_pos[1][1] << " " << triangle_target_pos[1][2] << std::endl;
	////std::cout << current_triangle_position[1][0] << " " << current_triangle_position[1][1] << " " << current_triangle_position[1][2] << std::endl;
	////std::cout << initial_triangle_position[2][0] << " " << initial_triangle_position[2][1] << " " << initial_triangle_position[2][2] << std::endl;
	////std::cout << triangle_target_pos[2][0] << " " << triangle_target_pos[2][1] << " " << triangle_target_pos[2][2] << std::endl;
	////std::cout << current_triangle_position[2][0] << " " << current_triangle_position[2][1] << " " << current_triangle_position[2][2] << std::endl;
	////std::cout << "===" << std::endl;


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
		//	////std::cout << "distance when <floor "<< initial_position[0]<<" "<< initial_position[1]<<" "<< initial_position[2]<<" "
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
		//	////std::cout << "distance when <floor" << initial_position[1] << " " << d_2 << " " << d_hat_2 << std::endl;
		//}
		if (d_2 >= d_hat_2) {
			return false;
		}
		calNearestPoint(nearest_point_barycentric, current_triangle_position[0], current_triangle_position[1], current_triangle_position[2], current_nearest_point);
	}

	stiffness *= barrier((d_hat_2- d_2)/d_hat_2, d_2 / d_hat_2);
	////std::cout << stiffness<<" "<< barrier(d_2 - d_hat_2, d_2 / d_hat_2) << std::endl;
	

	double d_move;// = sqrt(DOT(d, d));
	//d_move = d_hat - sqrt(d_2);
	////////////////////////////////////////////
	double d[3];
	for (int i = 0; i < 3; ++i) {
		d[i] = initial_position[i] - current_position[i];
	}
	d_move = sqrt(DOT(d, d));
	double d_min = d_hat - sqrt(d_2);
	if (d_move < d_min) {
		d_move = d_min;
	}
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
	//	////std::cout << vertex_target_pos[1] - current_triangle_position[0][1] << std::endl;
	//}
	//////std::cout << "build constraint " << std::endl;

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

	//////std::cout << "initial close point " << initial_close_point_1[0] << " " << initial_close_point_1[1] << " " << initial_close_point_1[2] << std::endl;
	//////std::cout << "initial close point " << initial_close_point_2[0] << " " << initial_close_point_2[1] << " " << initial_close_point_2[2] << std::endl;
	//////std::cout << "current close point " << current_close_point_1[0] << " " << current_close_point_1[1] << " " << current_close_point_1[2] << std::endl;
	//////std::cout << "current close point " << current_close_point_2[0] << " " << current_close_point_2[1] << " " << current_close_point_2[2] << std::endl;

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



void CollisionConstraint::testBarycentric()
{
	std::vector<std::array<double, 3>> initial_pos(5);// = ;
	std::vector<std::array<double, 3>> current_pos(5);// = { 0.01,-1.0,1.01 };

	initial_pos[0] = { 1.0,0.0,-1.0 };
	current_pos[0] = { 1.0,0.0,1.0 };

	initial_pos[1] = {1.0,0.01,0.0 };
	current_pos[1] = { -2.0,0.01,0.0 };

	initial_pos[2] = { 1.0,0.01,-2.0 };
	current_pos[2] = { -2.0,0.01,-2.0 };

	initial_pos[3] = { 1.0,0.01,2.0 };
	current_pos[3] = { -2.0,0.01, 2.0 };

	initial_pos[4] = { -2.0,0.0,-1.0 };
	current_pos[4] = { -3.0,0.0,1.0 };


	std::vector<std::array<double, 3>> initial_triangle_pos(3);
	std::vector<std::array<double, 3>> current_triangle_pos(3);
	initial_triangle_pos[0] = { 1.0,0.0,0.0 };
	initial_triangle_pos[1] = { -1.0,0.0,-1.0 };
	initial_triangle_pos[2] = { -1.0,0.0,1.0 };

	current_triangle_pos[0] = { 1.0,0.0,0.0 };
	current_triangle_pos[1] = { -1.0,0.0,-1.0 };
	current_triangle_pos[2] = { -1.0,0.0,1.0 };

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
	double barycentric[4];
	double d_hat_2 = 0.03 * 0.03;

	for (unsigned int i = 0; i < 5; ++i) {
		double d_2 = CCD::internal::edgeEdgeNearestPoint(initial_pos[i].data(), current_pos[i].data(),
			initial_triangle_pos[1].data(), initial_triangle_pos[2].data(), barycentric);
		//std::cout << barycentric[0] << " " << barycentric[1] << " " << barycentric[2]<<" "<<barycentric[3] << std::endl;
	}

}

void CollisionConstraint::testPT()
{
	//testBarycentric();

	double initial_pos[3] = { 0.0,0.01,0.0 };
	double current_pos[3] = { 0.0,0.01,0.1 };
	std::vector<std::array<double, 3>> initial_triangle_pos(3);
	std::vector<std::array<double, 3>> current_triangle_pos(3);
	initial_triangle_pos[0] = { 1.0,0.0,0.0 };
	initial_triangle_pos[1] = { -1.0,0.0,-1.0 };
	initial_triangle_pos[2] = { -1.0,0.0,1.0 };
	current_triangle_pos[0] = { 1.0,0.0,0.0 };
	current_triangle_pos[1] = { -1.0,0.0,-1.0 };
	current_triangle_pos[2] = { -1.0,0.0,1.0 };
	//std::vector<double*> initial_tri(3);
	//std::vector<double*> current_tri(3);
	//for (int i = 0; i < 3; ++i) {
	//	initial_tri[i] = initial_triangle_pos[i].data();
	//	current_tri[i] = current_triangle_pos[i].data();
	//}
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
	double d_hat = 0.03;
	std::array<double, 3> vertex_target_pos;
	std::vector<std::array<double, 3>>triangle_target_pos(3);
	double stiffness;
	double epsilon = 0.0;

	double mass_0_0 = 1.0;
	double mass_0_1 = 1.0;
	double mass_1_0 = 2.0;
	double mass_1_1 = 1.0;
	//if (pointSelfTriangle(initial_pos, current_pos, initial_tri, current_tri, initial_tri_normal, vertex_target_pos.data(),
	//	triangle_target_pos.data(), d_hat_2, stiffness, vertex_mass, triangle_mass)) {
	/*if (pointTriangleResponse(initial_pos, current_pos, initial_triangle_pos[0].data(), initial_triangle_pos[1].data(),
		initial_triangle_pos[2].data(), current_triangle_pos[0].data(), current_triangle_pos[1].data(), current_triangle_pos[2].data(),
		initial_tri_normal, vertex_target_pos.data(), triangle_target_pos[0].data(), triangle_target_pos[1].data(), triangle_target_pos[2].data(),
		d_hat, stiffness, epsilon, mass_0_0, mass_0_1, mass_1_0, mass_1_1)) {*/
	//if (pointTriangleColliderResponse(initial_pos, current_pos, initial_triangle_pos[0].data(), initial_triangle_pos[1].data(),
	//	initial_triangle_pos[2].data(),
	//	initial_tri_normal, vertex_target_pos.data(),
	//	d_hat, stiffness, epsilon)) {
	if (pointColliderTriangleResponse(initial_pos, initial_triangle_pos[0].data(), initial_triangle_pos[1].data(),
		initial_triangle_pos[2].data(), current_triangle_pos[0].data(), current_triangle_pos[1].data(), current_triangle_pos[2].data(),
		initial_tri_normal, triangle_target_pos[0].data(), triangle_target_pos[1].data(), triangle_target_pos[2].data(),
		d_hat, stiffness, epsilon, mass_0_1, mass_1_0, mass_1_1)) {
		////std::cout <<"vertex "<< vertex_target_pos[0] << " " << vertex_target_pos[1] << " " << vertex_target_pos[2] << std::endl;
		//////std::cout << "triangle vertex" << std::endl;
		for (int i = 0; i < 3; ++i) {
			//std::cout << triangle_target_pos[i][0] << " " << triangle_target_pos[i][1] << " " << triangle_target_pos[i][2] << std::endl;
		}
		////std::cout << stiffness << std::endl;
	}
	else {
		////std::cout << "does not collide " << std::endl;
	}
}

void CollisionConstraint::testEE()
{
	double initial_edge_0[3] = { 0.0,0.001,1.0 };
	double current_edge_0[3] = { 0.1,0.001,1.0 };
	double initial_edge_1[3] = { 0.0,0.001,0.0 };
	double current_edge_1[3] = { 0.1,0.001,0.0 };

	double initial_compare_edge_0[3] = { 0.0,0.0,0.01 };
	double current_compare_edge_0[3] = { 0.0,0.0,0.01 };
	double initial_compare_edge_1[3] = { 1.0,0.0,0.01 };
	double current_compare_edge_1[3] = { 1.0,0.0,0.01 };

	double mass[4] = { 0.5,1.5,0.5,1.5 };
	double d_hat = 0.03;
	std::vector<std::array<double, 3>> target_pos(2);
	std::vector<std::array<double, 3>> compare_pos(2);
	double stiffness;
	double epsilon = 0.0;
	double mass_0_0 = 1.0;
	double mass_0_1 = 2.0;
	double mass_1_0 = 1.0;
	double mass_1_1 = 1.0;
	//if (edgeEdgeResponse(target_pos[0].data(), target_pos[1].data(),compare_pos[0].data(), compare_pos[1].data(),current_edge_0, current_edge_1,initial_edge_0,initial_edge_1,
	//	current_compare_edge_0,current_compare_edge_1,initial_compare_edge_0,initial_compare_edge_1, d_hat,stiffness, epsilon,
	//	mass_0_0, mass_0_1, mass_1_0, mass_1_1)) {
	////if (edgeEdgeColliderResponse(target_pos[0].data(), target_pos[1].data(), current_edge_0, current_edge_1, initial_edge_0, initial_edge_1,
	////	initial_compare_edge_0, initial_compare_edge_1, d_hat, stiffness, epsilon,
	////	mass_0_0, mass_0_1)) {
	//	//std::cout << "edge 0 " << target_pos[0][0] << " " << target_pos[0][1] << " " << target_pos[0][2] << std::endl;
	//	//std::cout << "edge 1 " << target_pos[1][0] << " " << target_pos[1][1] << " " << target_pos[1][2] << std::endl;
	//	//std::cout << "compare 0 " << compare_pos[0][0] << " " << compare_pos[0][1] << " " << compare_pos[0][2] << std::endl;
	//	//std::cout << "compare 1 " << compare_pos[1][0] << " " << compare_pos[1][1] << " " << compare_pos[1][2] << std::endl;
	//	//////std::cout << stiffness << std::endl;
	//	SUB_(target_pos[0], target_pos[1]);
	//	//std::cout << "dis " << DOT(target_pos[0], target_pos[0]) << std::endl;
	//}
	//else {
	//	////std::cout << "does not collide " << std::endl;
	//}
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

	//////std::cout << "moving distance " << d_close[0] << " " << d_close[1]<<" ";
	d_close[1] = mass[2] * (mass[0] + mass[1]) * d_move / (mass_total * (barycentric[2] * mass[3] + barycentric[3] * mass[2]));
	d_close[0] = mass[3] * d_close[1] / mass[2];
	MULTI(compare_target_pos[0], direction, d_close[0]);
	SUM_(compare_target_pos[0], initial_compare_edge_vertex_0);
	MULTI(compare_target_pos[1], direction, d_close[1]);
	SUM_(compare_target_pos[1], initial_compare_edge_vertex_1);
	//////std::cout << d_close[0] << " " << d_close[1] << std::endl;

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
	////std::cout << triangle_mass[0] << " " << triangle_mass[1] << " " << triangle_mass[2] << " " << barycentric[0]<<" "<< barycentric[1]<<" "<< barycentric[2] << std::endl;

	triangle_d[1] = triangle_mass[0] * triangle_mass[2] * d_close / (barycentric[0] * triangle_mass[1] * triangle_mass[2] + barycentric[1] * triangle_mass[0] * triangle_mass[2]
		+ barycentric[2] * triangle_mass[0] * triangle_mass[1]);
	triangle_d[0] = triangle_mass[1] * triangle_d[1] / triangle_mass[0];
	triangle_d[2] = triangle_mass[1] * triangle_d[1] / triangle_mass[2];
	MULTI(vertex_target_pos, direction, vertex_d);
	SUM_(vertex_target_pos, initial_pos);
	
	for (int i = 0; i < 3; ++i) {
		MULTI(triangle_target_pos[i], direction, triangle_d[i]);
		SUM_(triangle_target_pos[i], initial_triangle_position[i]);
		//////std::cout << "direction " << triangle_d[i] << " " << triangle_d[i] << " " << triangle_d[i] << std::endl;
		//////std::cout << "direction " << triangle_d[i][0] << " " << triangle_d[i][1] << " " << triangle_d[i][2] << std::endl;
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
	//	////std::cout << "barycentric "<< barycentric[0]<<" "<< barycentric[1]<<" "<< barycentric[2] << std::endl;
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