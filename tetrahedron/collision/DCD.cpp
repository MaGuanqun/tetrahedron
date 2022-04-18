#include"DCD.h"




bool DCD::pointSelfTriangle(double* initial_position, double* current_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double* current_triangle_normal, double* vertex_target_pos,
    double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
    double tolerance, double mass_point, double mass_t0, double mass_t1, double mass_t2,
    double current_triangle_area)
{
	double c_p[3];
	SUB(c_p, current_position, current_triangle_position_0);
	double current_side = DOT(c_p, current_triangle_normal);
	SUB(c_p, initial_position, initial_triangle_position_0);
	double initial_side = DOT(c_p, initial_triangle_normal);

	if (initial_side * current_side > 0) {//abs(current_side) > tolerance && 
		return false;
	}
	else {
        if (!pointProjectOnTriangle(initial_position, initial_triangle_position_0, initial_triangle_position_1,
            initial_triangle_position_2, initial_triangle_normal)) {
            return false;
        }
	}

    //compute moving distance
    
    double in_triangle[3], scale_norm[3];
    SUB(in_triangle, current_position, current_triangle_position_0);
    MULTI(scale_norm, current_triangle_normal, current_side);
    SUB_(in_triangle, scale_norm);
    if (initial_side > 0) {
        current_side -= tolerance;
    }
    else {
        current_side += tolerance;
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

    double s =-current_side /( DOT(vertex_target_pos, vertex_target_pos) / mass_point + DOT(triangle_target_pos_0, triangle_target_pos_0) / mass_t0
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

    return true;


}

bool DCD::pointTriangleCollider(double* initial_position, double* current_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double* current_triangle_normal, double* vertex_target_pos,
    double tolerance)
{
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
            initial_triangle_position_2, initial_triangle_normal)) {
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


bool DCD::pointProjectOnTriangle(
    const double* p,
    const double* t0,
    const double* t1,
    const double* t2,
    const double* triangle_normal)
{
    double barycentric[2];
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

    if (barycentric[0] > 0.0 && barycentric[1] > 0.0 && 1.0 - barycentric[0] - barycentric[1] > 0.0) {
        return true;
    }
    else {
        return false;
    }
  
}