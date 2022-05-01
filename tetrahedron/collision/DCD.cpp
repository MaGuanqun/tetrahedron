#include"DCD.h"


bool DCD::checkPointTriangle(double* initial_position, double* current_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double tolerance)
{
    double barycentric[3];
    CCD::internal::pointTriangleNearestPoint(initial_position, initial_triangle_position_0, initial_triangle_position_1,
        initial_triangle_position_2, initial_triangle_normal, barycentric);
    if (checkIfCollidePointTriangle(initial_position, current_position,
        initial_triangle_position_0, initial_triangle_position_1, initial_triangle_position_2,
        current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
        initial_triangle_normal, barycentric,tolerance)) {
        return true;
    }
    return false;
}


void DCD::XPBDpointSelfTriangle(double* initial_position, double* current_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double* current_triangle_normal, 
    double tolerance, double mass_inv_point, double mass_inv_t0, double mass_inv_t1, double mass_inv_t2,
    double current_triangle_area, double& lambda, double stiffness, double damping_stiffness, double dt)
{
    double barycentric[3];
    double current_side;
    bool should_be_front;

    CCD::internal::pointTriangleNearestPoint(initial_position, initial_triangle_position_0, initial_triangle_position_1,
        initial_triangle_position_2, initial_triangle_normal, barycentric);
    if (checkIfCollidePointTriangle(initial_position, current_position,
        initial_triangle_position_0, initial_triangle_position_1, initial_triangle_position_2,
        current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
        initial_triangle_normal, current_triangle_normal, barycentric,
        tolerance, current_side, should_be_front)) {
        XPBDcalDistancePointTriangle(initial_position, initial_triangle_position_0, initial_triangle_position_1,
            initial_triangle_position_2,
            current_position, current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
            current_triangle_normal, current_side, tolerance, should_be_front, current_triangle_area,
            mass_inv_point, mass_inv_t0, mass_inv_t1, mass_inv_t2,lambda, stiffness, damping_stiffness,dt);
    }
}






bool DCD::pointSelfTriangle(double* initial_position, double* current_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double* current_triangle_normal, double* vertex_target_pos,
    double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
    double tolerance, double mass_point, double mass_t0, double mass_t1, double mass_t2,
    double current_triangle_area)
{
    double barycentric[3];
	double c_p[3];
    double current_side;
    bool should_be_front;
	//SUB(c_p, current_position, current_triangle_position_0);
	//current_side = DOT(c_p, current_triangle_normal);
	//SUB(c_p, initial_position, initial_triangle_position_0);
	//double initial_side = DOT(c_p, initial_triangle_normal);
	//if (initial_side * current_side > 0) {//abs(current_side) > tolerance && 
	//	return false;
	//}
	//else {
 //       if (!pointProjectOnTriangle(initial_position, initial_triangle_position_0, initial_triangle_position_1,
 //           initial_triangle_position_2, initial_triangle_normal, barycentric)) {
 //           return false;
 //       }
	//}
 //   if (initial_side > 0) {
 //       should_be_front = true;
 //   }
 //   else {
 //       should_be_front = false;
 //   }
 //   calDistancePointTriangle(vertex_target_pos, triangle_target_pos_0, triangle_target_pos_1, triangle_target_pos_2,
 //       current_position, current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
 //       current_triangle_normal, current_side, tolerance, should_be_front, current_triangle_area,
 //       mass_point, mass_t0, mass_t1, mass_t2);
 //   return true;


    CCD::internal::pointTriangleNearestPoint(initial_position, initial_triangle_position_0, initial_triangle_position_1,
        initial_triangle_position_2, initial_triangle_normal, barycentric);
    if (checkIfCollidePointTriangle(initial_position, current_position,
        initial_triangle_position_0, initial_triangle_position_1, initial_triangle_position_2,
        current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
        initial_triangle_normal, current_triangle_normal, barycentric,
        tolerance, current_side, should_be_front)){
        calDistancePointTriangle(vertex_target_pos, triangle_target_pos_0, triangle_target_pos_1, triangle_target_pos_2,
            current_position, current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
            current_triangle_normal, current_side, tolerance, should_be_front, current_triangle_area,
            mass_point, mass_t0, mass_t1, mass_t2);
        return true;
    }
    return false;
}

void DCD::XPBDcalDistancePointTriangle(
    double* initial_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* current_triangle_normal, double constraint, double tolerance, bool is_front, double current_triangle_area,
    double mass_inv_point, double mass_inv_t0, double mass_inv_t1, double mass_inv_t2, double& lambda, double stiffness, double damping_stiffness, double dt)
{
    double grad_c_vertex[3];
    double grad_c_vertex_0[3];
    double grad_c_vertex_1[3];
    double grad_c_vertex_2[3];
    //C, grac_C
    double in_triangle[3], scale_norm[3];
    SUB(in_triangle, current_position, current_triangle_position_0);
    MULTI(scale_norm, current_triangle_normal, constraint);
    SUB_(in_triangle, scale_norm);
    if (is_front) {
        if (constraint > tolerance) {
            return;
        }
        constraint -= tolerance;
    }
    else {
        if (constraint > -tolerance) {
            return;
        }
        constraint += tolerance;
    }
    memcpy(grad_c_vertex, current_triangle_normal, 24);
    current_triangle_area = 1.0 / current_triangle_area;

    double temp_vec[3];
    SUB(temp_vec, current_triangle_position_1, current_triangle_position_2);
    MAGNITUDE_CROSS(grad_c_vertex_0, temp_vec, in_triangle, current_triangle_area);
    SUB_(grad_c_vertex_0, grad_c_vertex);

    SUB(temp_vec, current_triangle_position_2, current_triangle_position_0);
    MAGNITUDE_CROSS(grad_c_vertex_1, temp_vec, in_triangle, current_triangle_area);

    SUB(temp_vec, current_triangle_position_0, current_triangle_position_1);
    MAGNITUDE_CROSS(grad_c_vertex_2, temp_vec, in_triangle, current_triangle_area);

    

    // lambda
    double alpha_ = 1.0 / (stiffness * dt * dt);
    double gamma = damping_stiffness / (stiffness * dt);

    //use in_triangle[3] to store x_i-x_ori
    double lambda_numerator = 0.0;
    SUB(in_triangle, current_position, initial_position);
    lambda_numerator += DOT(grad_c_vertex, in_triangle);
    SUB(in_triangle, current_triangle_position_0, initial_triangle_position_0);
    lambda_numerator += DOT(grad_c_vertex_0, in_triangle);
    SUB(in_triangle, current_triangle_position_1, initial_triangle_position_1);
    lambda_numerator += DOT(grad_c_vertex_1, in_triangle);
    SUB(in_triangle, current_triangle_position_2, initial_triangle_position_2);
    lambda_numerator += DOT(grad_c_vertex_2, in_triangle);
    
    double delta_lambda = -(constraint + alpha_ * lambda + gamma * lambda_numerator)
        / ((1 + gamma) * (DOT(grad_c_vertex, grad_c_vertex) * mass_inv_point + DOT(grad_c_vertex_0, grad_c_vertex_0) * mass_inv_t0
            + DOT(grad_c_vertex_1, grad_c_vertex_1) * mass_inv_t1 + DOT(grad_c_vertex_2, grad_c_vertex_2) * mass_inv_t2)
            + alpha_);
    lambda += delta_lambda;

    double coe = delta_lambda * mass_inv_point;
    ACCUMULATE_SUM_WITH_COE(current_position, coe, grad_c_vertex);
    coe = delta_lambda * mass_inv_t0;
    ACCUMULATE_SUM_WITH_COE(current_triangle_position_0, coe, grad_c_vertex_0);
    coe = delta_lambda * mass_inv_t1;
    ACCUMULATE_SUM_WITH_COE(current_triangle_position_1, coe, grad_c_vertex_1);
    coe = delta_lambda * mass_inv_t2;
    ACCUMULATE_SUM_WITH_COE(current_triangle_position_2, coe, grad_c_vertex_2);


}

void DCD::calDistancePointTriangle(double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
    double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* current_triangle_normal, double constraint, double tolerance, bool is_front, double current_triangle_area,
    double mass_point, double mass_t0, double mass_t1, double mass_t2)
{
    double in_triangle[3], scale_norm[3];
    SUB(in_triangle, current_position, current_triangle_position_0);
    MULTI(scale_norm, current_triangle_normal, constraint);
    SUB_(in_triangle, scale_norm);
    if (is_front) {
        constraint -= tolerance;
    }
    else {
        constraint += tolerance;
    }
    memcpy(vertex_target_pos, current_triangle_normal, 24);
    current_triangle_area = 1.0 / current_triangle_area;

    double temp_vec[3];
    SUB(temp_vec, current_triangle_position_1, current_triangle_position_2);
    MAGNITUDE_CROSS(triangle_target_pos_0, temp_vec, in_triangle, current_triangle_area);

    SUB(temp_vec, current_triangle_position_2, current_triangle_position_0);
    MAGNITUDE_CROSS(triangle_target_pos_1, temp_vec, in_triangle, current_triangle_area);

    SUB(temp_vec, current_triangle_position_0, current_triangle_position_1);
    MAGNITUDE_CROSS(triangle_target_pos_2, temp_vec, in_triangle, current_triangle_area);

    SUB_(triangle_target_pos_0, vertex_target_pos);

    double s = -constraint / (DOT(vertex_target_pos, vertex_target_pos) / mass_point + DOT(triangle_target_pos_0, triangle_target_pos_0) / mass_t0
        + DOT(triangle_target_pos_1, triangle_target_pos_1) / mass_t1 + DOT(triangle_target_pos_2, triangle_target_pos_2) / mass_t2);

    double tem_value;

    tem_value = s / mass_point;
    MULTI_SUM(vertex_target_pos, tem_value, current_position);
    tem_value = s / mass_t0;
    MULTI_SUM(triangle_target_pos_0, tem_value, current_triangle_position_0);
    tem_value = s / mass_t1;
    MULTI_SUM(triangle_target_pos_1, tem_value, current_triangle_position_1);
    tem_value = s / mass_t2;
    MULTI_SUM(triangle_target_pos_2, tem_value, current_triangle_position_2);
}



void DCD::XPBDedgeEdge(double* current_edge_vertex_0, double* current_edge_vertex_1,
    double* initial_edge_vertex_0, double* initial_edge_vertex_1,
    double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
    double* initial_compare_edge_vertex_1, double tolerance, double mass_inv_e_0_0, double mass_inv_e_0_1, double mass_inv_e_1_0, double mass_inv_e_1_1,
    double& lambda, double stiffness, double damping_stiffness, double dt)
{
    double barycentric[4];
    if (CCD::internal::edgeEdgeDistanceType(initial_edge_vertex_0, initial_edge_vertex_1,
        initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, barycentric) != 8) {
        return;
    }
    double norm[3];
    double distance2;
    if (checkEdgeEdgeCollision(current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0, initial_edge_vertex_1, current_compare_edge_vertex_0, current_compare_edge_vertex_1,
        initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, barycentric, norm, distance2, tolerance)) {
        if (distance2 > 0) {
            return;
        }
        XPBDcalDistanceEdgeEdge(
            norm, distance2, barycentric, current_edge_vertex_0, current_edge_vertex_1, current_compare_edge_vertex_0, current_compare_edge_vertex_1,
            initial_edge_vertex_0, initial_edge_vertex_1,
            initial_compare_edge_vertex_0, initial_compare_edge_vertex_1,
            mass_inv_e_0_0, mass_inv_e_0_1, mass_inv_e_1_0, mass_inv_e_1_1,
            lambda, stiffness, damping_stiffness, dt);
    }
}



bool DCD::edgeEdge(double* edge_target_pos_0, double* edge_target_pos_1,
    double* compare_target_pos_0, double* compare_target_pos_1, double* current_edge_vertex_0, double* current_edge_vertex_1, 
    double* initial_edge_vertex_0, double* initial_edge_vertex_1,
    double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
    double* initial_compare_edge_vertex_1, double tolerance, double mass_e_0_0, double mass_e_0_1, double mass_e_1_0, double mass_e_1_1)
{
    double barycentric[4];
    if (CCD::internal::edgeEdgeDistanceType(initial_edge_vertex_0, initial_edge_vertex_1,
        initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, barycentric) != 8) {
        return false;
    }
    double norm[3];
    double distance2;
    if (checkEdgeEdgeCollision(current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0, initial_edge_vertex_1, current_compare_edge_vertex_0, current_compare_edge_vertex_1,
        initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, barycentric, norm, distance2, tolerance)) {
        calDistanceEdgeEdge(edge_target_pos_0, edge_target_pos_1, compare_target_pos_0, compare_target_pos_1,
            norm, distance2, barycentric, current_edge_vertex_0, current_edge_vertex_1, current_compare_edge_vertex_0,
            current_compare_edge_vertex_1, mass_e_0_0, mass_e_0_1, mass_e_1_0, mass_e_1_1);
        return true;
    }
    return false;

}



bool DCD::checkEdgeEdge(double* current_edge_vertex_0, double* current_edge_vertex_1,
    double* initial_edge_vertex_0, double* initial_edge_vertex_1,
    double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
    double* initial_compare_edge_vertex_1, double tolerance)
{
    double barycentric[4];
    if (CCD::internal::edgeEdgeDistanceType(initial_edge_vertex_0, initial_edge_vertex_1,
        initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, barycentric) != 8) {
        return false;
    }
    if (checkIfCollideEdgeEdge(current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0, initial_edge_vertex_1, current_compare_edge_vertex_0, current_compare_edge_vertex_1,
        initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, barycentric, tolerance)) {
        return true;
    }
    return false;

}


void DCD::XPBDcalDistanceEdgeEdge(double* norm, double distance, double* alpha, double* current_edge_vertex_0, double* current_edge_vertex_1,
    double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1,
    double* initial_edge_vertex_0, double* initial_edge_vertex_1,
    double* initial_compare_edge_vertex_0, double* initial_compare_edge_vertex_1,
    double mass_inv_e_0_0, double mass_inv_e_0_1, double mass_inv_e_1_0, double mass_inv_e_1_1,
    double& lambda, double stiffness, double damping_stiffness, double dt)
{
    double coe = alpha[0] * alpha[0] * mass_inv_e_0_0 + alpha[1] * alpha[1] * mass_inv_e_0_1 + alpha[2] * alpha[2] * mass_inv_e_1_0 
        + alpha[3] * alpha[3] * mass_inv_e_1_1;
    
    //distance is C
    //grad_c
    //grad_c_p0 = delta_p_0_coe * norm and so on
    double delta_p_0_coe = alpha[0] * distance * mass_inv_e_0_0 / coe;
    double delta_p_1_coe = alpha[1] * distance * mass_inv_e_0_1 / coe;
    double delta_compare_p_0_coe = alpha[2] * distance * mass_inv_e_1_0 / coe;
    double delta_compare_p_1_coe = alpha[3] * distance * mass_inv_e_1_1 / coe;

    double e[3];
    //use e[3] to store x_i-x_ori
    double lambda_numerator = 0.0;
    SUB(e, current_edge_vertex_0, initial_edge_vertex_0);
    lambda_numerator += delta_p_0_coe * DOT(e, norm);
    SUB(e, current_edge_vertex_1, initial_edge_vertex_1);
    lambda_numerator += delta_p_1_coe * DOT(e, norm);
    SUB(e, current_compare_edge_vertex_0, initial_compare_edge_vertex_0);
    lambda_numerator += delta_compare_p_0_coe * DOT(e, norm);
    SUB(e, current_compare_edge_vertex_1, initial_compare_edge_vertex_1);
    lambda_numerator += delta_compare_p_1_coe * DOT(e, norm);


    // lambda
    double alpha_ = 1.0 / (stiffness * dt * dt);
    double gamma = damping_stiffness / (stiffness * dt);
    double delta_lambda = -(distance + alpha_ * lambda + gamma * lambda_numerator) / ((1 + gamma) * (mass_inv_e_0_0 * delta_p_0_coe * delta_p_0_coe
        + mass_inv_e_0_1 * delta_p_1_coe * delta_p_1_coe + mass_inv_e_1_0 * delta_compare_p_0_coe * delta_compare_p_0_coe +
        mass_inv_e_1_1 * delta_compare_p_1_coe * delta_compare_p_1_coe) + alpha_);
    lambda += delta_lambda;

    coe = mass_inv_e_0_0 * delta_lambda * delta_p_0_coe;
    ACCUMULATE_SUM_WITH_COE(current_edge_vertex_0, coe, norm);
    coe = mass_inv_e_0_1 * delta_lambda * delta_p_1_coe;
    ACCUMULATE_SUM_WITH_COE(current_edge_vertex_1, coe, norm);
    coe = mass_inv_e_1_0 * delta_lambda * delta_compare_p_0_coe;
    ACCUMULATE_SUM_WITH_COE(current_compare_edge_vertex_0, coe, norm);
    coe = mass_inv_e_1_1 * delta_lambda * delta_compare_p_1_coe;
    ACCUMULATE_SUM_WITH_COE(current_compare_edge_vertex_1, coe, norm);


}

void DCD::calDistanceEdgeEdge(double* edge_target_pos_0, double* edge_target_pos_1,
    double* compare_target_pos_0, double* compare_target_pos_1,
    double* norm, double distance, double* alpha, double* current_edge_vertex_0, double* current_edge_vertex_1,
    double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, 
    double mass_e_0_0, double mass_e_0_1, double mass_e_1_0, double mass_e_1_1)
{  
    double coe = alpha[0] * alpha[0] / mass_e_0_0 + alpha[1] * alpha[1] / mass_e_0_1 + alpha[2] * alpha[2] / mass_e_1_0+ alpha[3] * alpha[3] / mass_e_1_1;
    double temp;
    distance /= coe;

    temp = distance * alpha[0] / mass_e_0_0;
    MULTI_SUM2(edge_target_pos_0, temp, norm, current_edge_vertex_0);

    temp = distance * alpha[1] / mass_e_0_1;
    MULTI_SUM2(edge_target_pos_1, temp, norm, current_edge_vertex_1);

    temp = -distance * alpha[2] / mass_e_1_0;
    MULTI_SUM2(compare_target_pos_0, temp, norm, current_compare_edge_vertex_0);

    temp = -distance * alpha[3] / mass_e_1_1;
    MULTI_SUM2(compare_target_pos_1, temp, norm, current_compare_edge_vertex_1);
}


bool DCD::checkIfCollidePointTriangle(double* initial_point_position, double* current_point_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double* barycentric,double tolerance)
{
    bool should_be_front;
    double c_p[3];
    SUB(c_p, initial_point_position, initial_triangle_position_0);
    if (DOT(c_p, initial_triangle_normal) > 0) {
        should_be_front = true;
    }
    else {
        should_be_front = false;
    }
    double tc[3];
    double tc_new[3];
    double p_c[3];
    BARYCENTRIC(tc_new, barycentric, current_triangle_position_0, current_triangle_position_1, current_triangle_position_2);
    BARYCENTRIC(tc, barycentric, initial_triangle_position_0, initial_triangle_position_1, initial_triangle_position_2);
    SUB(p_c, initial_point_position, tc);
    double distance_p_c = sqrt(DOT(p_c, p_c));
    if (distance_p_c > NORM_NEAR_ZERO) {
        normalize(p_c);
    }
    else {
        //std::cout << "the distance between point and its closest point on the triangle is too close: "<<distance_p_c<< std::endl;
        if (should_be_front) {
            memcpy(p_c, initial_triangle_normal, 24);
        }
        else {
            MULTI(p_c, initial_triangle_normal, -1.0);
        }
    }
    double displacement_c_p[3];
    displacement_c_p[0] = tc_new[0] - tc[0] - current_point_position[0] + initial_point_position[0];
    displacement_c_p[1] = tc_new[1] - tc[1] - current_point_position[1] + initial_point_position[1];
    displacement_c_p[2] = tc_new[2] - tc[2] - current_point_position[2] + initial_point_position[2];
    double distance;
    distance = distance_p_c - tolerance;
    double dp_dc_project = DOT(displacement_c_p, p_c);

    if (dp_dc_project > distance) {
        return true;
    }
    return false;
}


bool DCD::checkIfCollidePointTriangleCollider(double* initial_point_position, double* current_point_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal,  double* barycentric,
    double tolerance)
{
    double tc[3];
    double tc_new[3];
    double p_c[3];
    BARYCENTRIC(tc_new, barycentric, current_triangle_position_0, current_triangle_position_1, current_triangle_position_2);
    BARYCENTRIC(tc, barycentric, initial_triangle_position_0, initial_triangle_position_1, initial_triangle_position_2);
    SUB(p_c, initial_point_position, tc);
    double distance_p_c = sqrt(DOT(p_c, p_c));
    if (distance_p_c > NORM_NEAR_ZERO) {
        normalize(p_c);
    }
    else {
        //std::cout << "the distance between point and its closest point on the triangle is too close: "<<distance_p_c<< std::endl;
        memcpy(p_c, initial_triangle_normal, 24);
    }

    double displacement_c_p[3];
    displacement_c_p[0] = tc_new[0] - tc[0] - current_point_position[0] + initial_point_position[0];
    displacement_c_p[1] = tc_new[1] - tc[1] - current_point_position[1] + initial_point_position[1];
    displacement_c_p[2] = tc_new[2] - tc[2] - current_point_position[2] + initial_point_position[2];
    double distance;
    distance = distance_p_c - tolerance;
    double dp_dc_project = DOT(displacement_c_p, p_c);

    if (dp_dc_project > distance) {
        return true;
    }
    return false;
}


bool DCD::checkIfCollidePointTriangle(double* initial_point_position, double* current_point_position, 
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double* current_triangle_normal, double* barycentric,
    double tolerance, double& triangle_side2, bool& should_be_front) 
{

    double c_p[3];
    SUB(c_p, initial_point_position, initial_triangle_position_0);
    if (DOT(c_p, initial_triangle_normal)>0) {
        should_be_front = true;
    }
    else {
        should_be_front = false;
    }
    SUB(c_p, current_point_position, current_triangle_position_0);
    triangle_side2 = DOT(c_p, current_triangle_normal);

    double tc[3];
    double tc_new[3];

    double p_c[3];

    BARYCENTRIC(tc_new, barycentric, current_triangle_position_0,current_triangle_position_1, current_triangle_position_2);
    BARYCENTRIC(tc, barycentric, initial_triangle_position_0, initial_triangle_position_1, initial_triangle_position_2);
    SUB(p_c, initial_point_position, tc);
    double distance_p_c = sqrt(DOT(p_c, p_c));
    if (distance_p_c > NORM_NEAR_ZERO) {
        normalize(p_c);
    }
    else {
        //std::cout << "the distance between point and its closest point on the triangle is too close: "<<distance_p_c<< std::endl;
        if (should_be_front) {
            memcpy(p_c, initial_triangle_normal, 24);
        }
        else {
            MULTI(p_c, initial_triangle_normal, -1.0);
        }
    }
  
    double displacement_c_p[3];
    displacement_c_p[0] = tc_new[0] - tc[0] - current_point_position[0] + initial_point_position[0];
    displacement_c_p[1] = tc_new[1] - tc[1] - current_point_position[1] + initial_point_position[1];
    displacement_c_p[2] = tc_new[2] - tc[2] - current_point_position[2] + initial_point_position[2];
    double distance;
    distance = distance_p_c - tolerance;
    double dp_dc_project = DOT(displacement_c_p, p_c);

    if (dp_dc_project > distance) {
        //std::cout << "===========" << std::endl;
        //std::cout << barycentric[0] << " " << barycentric[1] << " " << barycentric[2]<<" "<<DOT(initial_triangle_normal, initial_triangle_normal)<<" "<< triangle_side2 << std::endl;
        //std::cout << initial_triangle_normal[0] << " " << initial_triangle_normal[1] << " " << initial_triangle_normal[2] << std::endl;
        //SUB(c_p, initial_triangle_position_1, initial_triangle_position_0);
        //std::cout << DOT(c_p, initial_triangle_normal) << std::endl;
        //std::cout << distance_p_c << " " << distance << std::endl;
        //std::cout << "==" << std::endl;
        //std::cout << initial_point_position[0] << " " << initial_point_position[1] << " "
        //    << initial_point_position[2] << std::endl;
        //std::cout << current_point_position[0] << " " << current_point_position[1] << " "
        //    << current_point_position[2] << std::endl;
        //std::cout << initial_triangle_position_0[0] << " " << initial_triangle_position_0[1] << " "
        //    << initial_triangle_position_0[2] << std::endl;
        //std::cout << current_triangle_position_0[0] << " " << current_triangle_position_0[1] << " "
        //    << current_triangle_position_0[2] << std::endl;
        //std::cout << initial_triangle_position_1[0] << " " << initial_triangle_position_1[1] << " "
        //    << initial_triangle_position_1[2] << std::endl;
        //std::cout << current_triangle_position_1[0] << " " << current_triangle_position_1[1] << " "
        //    << current_triangle_position_1[2] << std::endl;
        //std::cout << initial_triangle_position_2[0] << " " << initial_triangle_position_2[1] << " "
        //    << initial_triangle_position_2[2] << std::endl;
        //std::cout << current_triangle_position_2[0] << " " << current_triangle_position_2[1] << " "
        //    << current_triangle_position_2[2] << std::endl;
        return true;
    }
    return false;
}



bool DCD::checkIfCollideEdgeEdge(double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
    double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0, double* initial_compare_edge_vertex_1,
    double* alpha, double tolerance)
{
    double norm[3];
    double e0[3];
    double e1[3];
    POINT_ON_EDGE(e0, alpha[0], alpha[1], initial_edge_vertex_0, initial_edge_vertex_1);
    POINT_ON_EDGE(e1, alpha[2], alpha[3], initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
    SUB(norm, e1, e0);
    double distance = sqrt(DOT(norm, norm));

    if (distance > NORM_NEAR_ZERO) {
        DEV_(norm,distance);
    }
    else {
        double edge0[3], edge1[3];
        SUB(edge0, initial_edge_vertex_0, initial_edge_vertex_1);
        SUB(edge1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
        CROSS(norm, edge0, edge1);
        if (DOT(norm, norm) > NEAR_ZERO2) {
            //std::cout << "One edge is getting too close to the other edge" << std::endl;
            normalize(norm);
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
    SUB_(e0, e1);
    double dp_dc_project = DOT(norm, e0);
    if (dp_dc_project > distance - tolerance) {
        return true;
    }
    else {
        return false;
    }
}



bool DCD::checkEdgeEdgeCollision(double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
    double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0, double* initial_compare_edge_vertex_1,
    double* alpha, double* norm, double& distance2, double tolerance)
{
    double e0[3];
    double e1[3];
    POINT_ON_EDGE(e0, alpha[0], alpha[1], initial_edge_vertex_0, initial_edge_vertex_1);
    POINT_ON_EDGE(e1, alpha[2], alpha[3], initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
    SUB(norm, e1, e0);
    double distance = sqrt(DOT(norm, norm));

    if (distance > NORM_NEAR_ZERO) {
        normalize(norm);
    }
    else {
        double edge0[3], edge1[3];
        SUB(edge0, initial_edge_vertex_0, initial_edge_vertex_1);
        SUB(edge1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
        CROSS(norm, edge0, edge1);
        if (DOT(norm, norm) > NEAR_ZERO2) {
            //std::cout << "One edge is getting too close to the other edge" << std::endl;
            normalize(norm);
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

    POINT_ON_EDGE(e0, alpha[0], alpha[1], current_edge_vertex_0, current_edge_vertex_1);
    POINT_ON_EDGE(e1, alpha[2], alpha[3], current_compare_edge_vertex_0, current_compare_edge_vertex_1);
    SUB(e0, e1, e0);
    distance2 = DOT(norm, e0)-tolerance;

    double displacement0[3], displacement1[3];
    SUB(displacement0, current_edge_vertex_0, initial_edge_vertex_0);
    SUB(displacement1, current_edge_vertex_1, initial_edge_vertex_1);
    POINT_ON_EDGE(e0, alpha[0], alpha[1], displacement0, displacement1);
    SUB(displacement0, current_compare_edge_vertex_0, initial_compare_edge_vertex_0);
    SUB(displacement1, current_compare_edge_vertex_1, initial_compare_edge_vertex_1);
    POINT_ON_EDGE(e1, alpha[2], alpha[3], displacement0, displacement1);
    SUB_(e0, e1);
    double dp_dc_project = DOT(norm, e0);
    if (dp_dc_project > distance - tolerance) {
        return true;
    }
    else {
        return false;
    }
}


void DCD::XPBDpointTriangleCollider(double* initial_position, double* current_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double* current_triangle_normal, double mass_inv_v,
    double tolerance, double& lambda, double stiffness, double damping_stiffness, double dt)
{
    double barycentric[3];
    double current_side;

    double c_p[3];
    SUB(c_p, current_position, current_triangle_position_0);
    current_side = DOT(c_p, current_triangle_normal);
    if (current_side > tolerance) {
        return;
    }
    CCD::internal::pointTriangleNearestPoint(initial_position, initial_triangle_position_0, initial_triangle_position_1,
        initial_triangle_position_2, initial_triangle_normal, barycentric);
    if (checkIfCollidePointTriangleCollider(initial_position, current_position,
        initial_triangle_position_0, initial_triangle_position_1, initial_triangle_position_2,
        current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
        initial_triangle_normal, barycentric, tolerance)) {
        XPBDcalDistancePointTriangleCollider(initial_position,
            current_position, mass_inv_v,
            current_triangle_normal, current_side, tolerance, lambda,stiffness,damping_stiffness,dt);
    }
}


bool DCD::pointTriangleCollider(double* initial_position, double* current_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double* current_triangle_normal, double* vertex_target_pos,
    double tolerance)
{
    double barycentric[2];
    double c_p[3];
    SUB(c_p, initial_position, initial_triangle_position_0);
    double initial_side = DOT(c_p, initial_triangle_normal);
    SUB(c_p, current_position, current_triangle_position_0);
    double current_side = DOT(c_p, current_triangle_normal);
  
    if (current_side > tolerance) {
    //if (abs(current_side) > tolerance && initial_side * current_side > 0) {
        return false;
    }
    else {
        if (!pointProjectOnTriangle(initial_position, initial_triangle_position_0, initial_triangle_position_1,
            initial_triangle_position_2, initial_triangle_normal, barycentric)) {
            return false;
        }
    }


    //if (initial_side * current_side < 0) {
        //if (initial_side > 0) {
        //    current_side += tolerance;
        //}
        //else {
            current_side -= tolerance;
        //}      
    //}
    //else {
    //    if (current_side > 0) {
    //        current_side -= tolerance;
    //    }
    //    else {
    //        current_side += tolerance;
    //    }
    //}  
    vertex_target_pos[0] = current_position[0] - current_side * current_triangle_normal[0];
    vertex_target_pos[1] = current_position[1] - current_side * current_triangle_normal[1];
    vertex_target_pos[2] = current_position[2] - current_side * current_triangle_normal[2];
    return true;
}


void DCD::XPBDcalDistancePointTriangleCollider(double* initial_position,
    double* current_position,double mass_inv_vertex,
    double* current_triangle_normal, double constraint, double tolerance,
    double& lambda, double stiffness, double damping_stiffness, double dt)
{
    constraint -= tolerance;
    // lambda
    double alpha_ = 1.0 / (stiffness * dt * dt);
    double gamma = damping_stiffness / (stiffness * dt);

    double e[3];
    SUB(e, current_position, initial_position);
    double delta_lambda = -(constraint + alpha_ * lambda + gamma * DOT(current_triangle_normal, e)) /
        ((1 + gamma) * mass_inv_vertex + alpha_);
    lambda += delta_lambda;
    double coe = mass_inv_vertex * delta_lambda;
    ACCUMULATE_SUM_WITH_COE(current_position, coe, current_triangle_normal);
}

bool DCD::pointProjectOnTriangle(
    const double* p,
    const double* t0,
    const double* t1,
    const double* t2,
    const double* triangle_normal, double* barycentric)
{
    /*double barycentric[2];*/
    double S[3];
    SUB(S, p, t0);
    double E1[3], E2[3], S1[3], S2[3];
    SUB(E1, t1, t0);
    SUB(E2, t2, t0);
    CROSS(S1, triangle_normal, E2);
    CROSS(S2, S, E1);
    double temp = 1.0 / DOT(S1, E1);
    barycentric[0] = temp * DOT(S1, S);
    barycentric[1] = temp * DOT(S2, triangle_normal);

    if (barycentric[0] >= 0.0 && barycentric[1] >= 0.0 && 1.0 - barycentric[0] - barycentric[1] >= 0.0) {
        return true;
    }
    else {
        return false;
    }
  
}