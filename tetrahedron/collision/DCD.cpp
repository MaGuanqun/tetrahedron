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
    double triangle_normal_magnitude_reciprocal, double& lambda, double stiffness, double damping_stiffness, double dt,
    double& energy)
{
    double barycentric[3];
    double current_side;
    bool should_be_front;
    energy = 0.0;

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
            current_triangle_normal, current_side, tolerance, should_be_front, triangle_normal_magnitude_reciprocal,
            mass_inv_point, mass_inv_t0, mass_inv_t1, mass_inv_t2,lambda, stiffness, damping_stiffness,dt, energy);//
    }
}


void DCD::test()
{
    double info_p[3], info_v0[3], info_v1[3], info_v2[3];
    for (unsigned int k = 0; k < 100; ++k) {
        info_p[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
        info_p[1]= ((double)rand() / (double)RAND_MAX)*2.0 * M_PI;
        info_p[2] = 1.0;

        info_v0[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
        info_v0[1] = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
        info_v0[2] = 1.0;

        info_v1[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
        info_v1[1] = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
        info_v1[2] = 1.0;

        info_v2[0] = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI;
        info_v2[1] = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
        info_v2[2] = 1.0;

        test(info_p, info_v0, info_v1, info_v2);

	}
}




//vertex-triangle
void DCD::test(double* info_p, double* info_v0, double* info_v1, double* info_v2) //[0] angle0, [1] angle1, [2]radius
{



	double initial_pos[3] = { info_p[2]*sin(info_p[0]),info_p[2] * cos(info_p[0]) * sin(info_p[1]),info_p[2] * cos(info_p[0])* cos(info_p[1]) };
    double current_pos[3];
    memcpy(current_pos, initial_pos, 24);
	double initial_triangle_0[3] = { info_v0[2] * sin(info_v0[0]),info_v0[2] * cos(info_v0[0]) * sin(info_v0[1]),info_v0[2] * cos(info_v0[0]) * cos(info_v0[1]) };
	double initial_triangle_1[3] = { info_v1[2] * sin(info_v1[0]),info_v1[2] * cos(info_v1[0]) * sin(info_v1[1]),info_v1[2] * cos(info_v1[0]) * cos(info_v1[1]) };
	double initial_triangle_2[3] = { info_v2[2] * sin(info_v2[0]),info_v2[2] * cos(info_v2[0]) * sin(info_v2[1]),info_v2[2] * cos(info_v2[0]) * cos(info_v2[1]) };
	double current_triangle_0[3];
	double current_triangle_1[3];
	double current_triangle_2[3];
    memcpy(current_triangle_0, initial_triangle_0, 24);
    memcpy(current_triangle_1, initial_triangle_1, 24);
    memcpy(current_triangle_2, initial_triangle_2, 24);
	//double current_triangle_1[3] = { -1.0,0.1,-0.8 };
	//double current_triangle_2[3] = { -1.0,-0.1,0.8 };
	double initial_tri_normal[3]; double tri_normal[3];
    double temp[3];
	double e1[3], e2[3];
	double e3[3], e4[3];
	SUB(e1, initial_triangle_1, initial_triangle_0);
	SUB(e2, initial_triangle_2, initial_triangle_0);
	CROSS(initial_tri_normal, e1, e2);
    normalize(initial_tri_normal);
	SUB(e3, current_triangle_1, current_triangle_0);
	SUB(e4, current_triangle_2, current_triangle_0);
	CROSS(tri_normal, e3, e4);
    double current_area = 0.5 * sqrt(DOT(tri_normal, tri_normal));
    normalize(tri_normal);
	double t;
	double e_0[3], e_1[3], e_2[3];
	SUB(e_0, initial_pos, initial_triangle_0);
	SUB(e_1, initial_pos, initial_triangle_1);
	SUB(e_2, initial_pos, initial_triangle_2);
	double u[3], u_0[3], u_1[3], u_2[3];
	SUB(u, current_pos, initial_pos);
	SUB(u_0, current_triangle_0, initial_triangle_0);
	SUB(u_1, current_triangle_1, initial_triangle_1);
	SUB(u_2, current_triangle_2, initial_triangle_2);
	double e_[3];
	SUB(e_, initial_triangle_1, initial_triangle_0);
	double tolerance_2 = 0.9;
    double massv=1.0;
    double mass0=2.0;
    double mass1=1.5;
    double mass2=1.0;
    double target_pos_v[3]; double target_pos_t0[3]; double target_pos_t1[3]; double target_pos_t2[3];
    double distance_current = distanceFromTriangle(current_pos, current_triangle_0, current_triangle_1, current_triangle_2);

    //if (distance_current > 1e-15) {
    //    calAccurateDistancePointTriangle(target_pos_v, target_pos_t0, target_pos_t1, target_pos_t2,
    //        current_pos, current_triangle_0, current_triangle_1, current_triangle_2,
    //        tolerance_2, massv, mass0, mass1, mass2);
    //    double distance_result = distanceFromTriangle(target_pos_v, target_pos_t0, target_pos_t1, target_pos_t2);
    //    std::cout << "distance " << distance_current << " " << distance_result << std::endl;
    //}

    bool is_front;
    if (abs(distance_current) > 1e-15) {
        if (distance_current > 0) {
            is_front = false;
        }
        else {
            is_front = true;
        }
        double tolerance;

        double triangle_normal[3];
        double triangle_normal_magnitude_reciprocal, constraint;
        decideTolerance(tolerance, 1e-1, 1e-2, current_pos, current_triangle_0, current_triangle_1, current_triangle_2);
        computeNormalAndConstraint(current_pos, current_triangle_0, current_triangle_1, current_triangle_2,
        triangle_normal, triangle_normal_magnitude_reciprocal, constraint,is_front,tolerance);
        iterationToGetResult(target_pos_v, target_pos_t0, target_pos_t1, target_pos_t2, current_pos, current_triangle_0, current_triangle_1, current_triangle_2,
            triangle_normal, constraint, tolerance, is_front, triangle_normal_magnitude_reciprocal, massv, mass0, mass1, mass2);
    }
}



double DCD::distanceFromTriangle(double* vertex_target_pos,
    double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2)
{
    double normal[3];
    double e0[3], e1[3], e2[3];
    SUB(e0, vertex_target_pos, triangle_target_pos_0);
    SUB(e1, triangle_target_pos_1, triangle_target_pos_0);
    SUB(e2, triangle_target_pos_2, triangle_target_pos_0);
    CROSS(normal, e1, e2);
    if (sqrt(DOT(normal, normal)) < 1e-15) {
        return 0.0;
    }
    normalized(normal);
    return DOT(e0, normal);
}


bool DCD::accuratePointSelfTriangle(double* initial_position, double* current_position,
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

    CCD::internal::pointTriangleNearestPoint(initial_position, initial_triangle_position_0, initial_triangle_position_1,
        initial_triangle_position_2, initial_triangle_normal, barycentric);
    if (checkIfCollidePointTriangle(initial_position, current_position,
        initial_triangle_position_0, initial_triangle_position_1, initial_triangle_position_2,
        current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
        initial_triangle_normal, current_triangle_normal, barycentric,
        tolerance, current_side, should_be_front)) {
        calAccurateDistancePointTriangle(vertex_target_pos, triangle_target_pos_0, triangle_target_pos_1, triangle_target_pos_2,
            current_position, current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
           tolerance, mass_point, mass_t0, mass_t1, mass_t2);
        return true;
    }
    return false;

}




bool DCD::pointSelfTriangle(double* initial_position, double* current_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double* current_triangle_normal, double* vertex_target_pos,
    double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
    double tolerance, double mass_point, double mass_t0, double mass_t1, double mass_t2,
    double triangle_normal_magnitude_reciprocal)
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
        //iterationToGetResult(vertex_target_pos, triangle_target_pos_0, triangle_target_pos_1, triangle_target_pos_2,
        //    current_position, current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
        //    current_triangle_normal, current_side, tolerance, should_be_front, triangle_normal_magnitude_reciprocal,
        //    mass_point, mass_t0, mass_t1, mass_t2);
        calDistancePointTriangle(vertex_target_pos, triangle_target_pos_0, triangle_target_pos_1, triangle_target_pos_2,
            current_position, current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
            current_triangle_normal, current_side, tolerance, should_be_front, triangle_normal_magnitude_reciprocal,
            mass_point, mass_t0, mass_t1, mass_t2);
        return true;
    }
    return false;
}

void DCD::XPBDcalDistancePointTriangle(
    double* initial_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* current_triangle_normal, double constraint, double tolerance, bool is_front, double triangle_normal_magnitude_reciprocal,
    double mass_inv_point, double mass_inv_t0, double mass_inv_t1, double mass_inv_t2, double& lambda, double stiffness, double damping_stiffness, double dt,
    double& energy)
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
        if (constraint < -tolerance) {
            return;
        }
        constraint += tolerance;
    }

    energy = 0.5 * stiffness * constraint * constraint;

    memcpy(grad_c_vertex, current_triangle_normal, 24);

    double temp_vec[3];
    SUB(temp_vec, current_triangle_position_1, current_triangle_position_2);
    MAGNITUDE_CROSS(grad_c_vertex_0, temp_vec, in_triangle, triangle_normal_magnitude_reciprocal);
    SUB_(grad_c_vertex_0, grad_c_vertex);

    SUB(temp_vec, current_triangle_position_2, current_triangle_position_0);
    MAGNITUDE_CROSS(grad_c_vertex_1, temp_vec, in_triangle, triangle_normal_magnitude_reciprocal);

    SUB(temp_vec, current_triangle_position_0, current_triangle_position_1);
    MAGNITUDE_CROSS(grad_c_vertex_2, temp_vec, in_triangle, triangle_normal_magnitude_reciprocal);

    

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

    //double delta_lambda = -constraint / (DOT(grad_c_vertex, grad_c_vertex) * mass_inv_point + DOT(grad_c_vertex_0, grad_c_vertex_0) * mass_inv_t0
    //    + DOT(grad_c_vertex_1, grad_c_vertex_1) * mass_inv_t1 + DOT(grad_c_vertex_2, grad_c_vertex_2) * mass_inv_t2);

    double coe = delta_lambda * mass_inv_point;
    ACCUMULATE_SUM_WITH_COE(current_position, coe, grad_c_vertex);
    coe = delta_lambda * mass_inv_t0;
    ACCUMULATE_SUM_WITH_COE(current_triangle_position_0, coe, grad_c_vertex_0);
    coe = delta_lambda * mass_inv_t1;
    ACCUMULATE_SUM_WITH_COE(current_triangle_position_1, coe, grad_c_vertex_1);
    coe = delta_lambda * mass_inv_t2;
    ACCUMULATE_SUM_WITH_COE(current_triangle_position_2, coe, grad_c_vertex_2);


}





void DCD::calAccurateDistancePointTriangle(double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
    double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double tolerance, double mass_point, double mass_t0, double mass_t1, double mass_t2)
{
    double center[3];
    double mass_total = mass_point + mass_t0 + mass_t1 + mass_t2;
    center[0] = (mass_point * current_position[0] + mass_t0 * current_triangle_position_0[0] + mass_t1 * current_triangle_position_1[0] + 
        mass_t2 * current_triangle_position_2[0])/mass_total;
    center[1] = (mass_point * current_position[1] + mass_t0 * current_triangle_position_0[1] + mass_t1 * current_triangle_position_1[1] +
        mass_t2 * current_triangle_position_2[1]) / mass_total;
    center[2] = (mass_point * current_position[2] + mass_t0 * current_triangle_position_0[2] + mass_t1 * current_triangle_position_1[2] +
        mass_t2 * current_triangle_position_2[2]) / mass_total;
    
    Matrix<double,3,4> p;
    for (unsigned int i = 0; i < 3; ++i) {
        p.data()[i] = current_position[i] - center[i];
        p.data()[3+i] = current_triangle_position_0[i] - center[i];
        p.data()[6+i] = current_triangle_position_1[i] - center[i];
        p.data()[9+i] = current_triangle_position_2[i] - center[i];
    }
    Matrix3d P_;
    Vector4d mass = Vector4d(mass_point, mass_t0, mass_t1, mass_t2);
    P_ = p * mass.asDiagonal() * p.transpose();
    SelfAdjointEigenSolver<Matrix3d> es(P_);
    Vector3d eigen_value = es.eigenvectors().col(0);
    
    Vector3d delta;
    delta = eigen_value * (- dotProduct(eigen_value, p.col(0)));
    SUM(vertex_target_pos, current_position, delta);
    delta = eigen_value * (-dotProduct(eigen_value, p.col(1)));
    SUM(triangle_target_pos_0, current_triangle_position_0, delta);
    delta = eigen_value * (-dotProduct(eigen_value, p.col(2)));
    SUM(triangle_target_pos_1, current_triangle_position_1, delta);
    delta = eigen_value *(-dotProduct(eigen_value, p.col(3)));
    SUM(triangle_target_pos_2, current_triangle_position_2, delta);   
}



void DCD::iterationToGetResult(double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
    double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double*  triangle_normal_ori, double constraint, double tolerance, bool is_front, double triangle_normal_magnitude_reciprocal,
    double mass_point, double mass_t0, double mass_t1, double mass_t2)
{
    double delta_p[3];
    double delta_v0[3]; double delta_v1[3]; double delta_v2[3];
    double triangle_normal[3];
    memcpy(triangle_normal, triangle_normal_ori, 24);
    memcpy(vertex_target_pos, current_position, 24);
    memcpy(triangle_target_pos_0, current_triangle_position_0, 24);
    memcpy(triangle_target_pos_1, current_triangle_position_1, 24);
    memcpy(triangle_target_pos_2, current_triangle_position_2, 24);

    //decideTolerance(tolerance, 1e-1, 1e-2, vertex_target_pos, triangle_target_pos_0, triangle_target_pos_1, triangle_target_pos_2);
    //computeNormalAndConstraint(vertex_target_pos, triangle_target_pos_0, triangle_target_pos_1, triangle_target_pos_2,
    //    triangle_normal, triangle_normal_magnitude_reciprocal, constraint,is_front,tolerance);
    double lambda = 0.0;
    double L = 0;
    int itr_num = 0;

    while (!convergeCondition(itr_num, lambda, vertex_target_pos, triangle_target_pos_0, triangle_target_pos_1, triangle_target_pos_2,
        current_position, current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
        mass_point, mass_t0, mass_t1, mass_t2,constraint,is_front, L))
    {
        calPositionMove(delta_p, delta_v0, delta_v1, delta_v2,
            vertex_target_pos, triangle_target_pos_0, triangle_target_pos_1, triangle_target_pos_2,
            triangle_normal, constraint, tolerance, is_front, triangle_normal_magnitude_reciprocal,
            mass_point, mass_t0, mass_t1, mass_t2,lambda);
        SUM_(vertex_target_pos, delta_p);
        SUM_(triangle_target_pos_0, delta_v0);
        SUM_(triangle_target_pos_1, delta_v1);
        SUM_(triangle_target_pos_2, delta_v2);

        computeNormalAndConstraint(vertex_target_pos, triangle_target_pos_0, triangle_target_pos_1, triangle_target_pos_2,
            triangle_normal, triangle_normal_magnitude_reciprocal, constraint, is_front, tolerance);
    }

}


bool DCD::convergeCondition(int& itr_num, double lambda, double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
    double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double mass_point, double mass_t0, double mass_t1, double mass_t2, double constraint, double is_front,
    double& L_)
{
    double delta_v[3];
    double delta_t0[3]; double delta_t1[3]; double delta_t2[3];

    SUB(delta_v, vertex_target_pos, current_position);
    SUB(delta_t0, triangle_target_pos_0, current_triangle_position_0);
    SUB(delta_t1, triangle_target_pos_1, current_triangle_position_1);
    SUB(delta_t2, triangle_target_pos_2, current_triangle_position_2);

    double L = 0.5 * (mass_point * DOT(delta_v, delta_v) + mass_t0 * DOT(delta_t0, delta_t0) + mass_t1 * DOT(delta_t1, delta_t1)
        + mass_t2 * DOT(delta_t2, delta_t2));

    if (is_front) {
        L += lambda * constraint;
    }
    else {
        L -= lambda * constraint;
    }
    

    double max = abs(delta_v[0]);
    for (unsigned int i = 0; i < 3; ++i) {
        if (max > abs(delta_v[i])) {
            max = abs(delta_v[i]);
        }
        if (max > abs(delta_t0[i])) {
            max = abs(delta_t0[i]);
        }
        if (max > abs(delta_t1[i])) {
            max = abs(delta_t1[i]);
        }
        if (max > abs(delta_t2[i])) {
            max = abs(delta_t2[i]);
        }
    }

    //if (is_front) {
    //    std::cout << "itr num " << itr_num << " " << is_front << " " << L << " " << " " << constraint<< std::endl;
    //}
    //else {
    //    std::cout << "itr num " << itr_num << " " << is_front << " " << L << " " << " " << constraint << std::endl;
    //}
    itr_num++;

    if (itr_num >= 10) {
        L_ = L;
        //std::cout << itr_num<<" " << std::endl;
        return true;
    }

    //if (max > max_move  && itr_num >1) {
    if (abs(L-L_)/L_ < 1e-3  && itr_num >1) {
        L_ = L;
        return true;
    }
    L_ = L;
    return false;

}

void DCD::decideTolerance(double& tolerance, double ratio_tolerance, double ratio_max_move, double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2)
{
    double e0[3], e1[3], e2[3];
    SUB(e0, triangle_target_pos_1, triangle_target_pos_2);
    SUB(e1, triangle_target_pos_1, triangle_target_pos_0);
    SUB(e2, triangle_target_pos_2, triangle_target_pos_0);
    tolerance = ratio_tolerance * (sqrt(DOT(e0, e0)) + sqrt(DOT(e1, e1)) + sqrt(DOT(e2, e2))) / 3.0;
    //max_move = ratio_max_move * (sqrt(DOT(e0, e0)) + sqrt(DOT(e1, e1)) + sqrt(DOT(e2, e2))) / 3.0;

}


void DCD::computeNormalAndConstraint(double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
    double* normal, double& triangle_normal_magnitude_reciprocal, double& constraint, bool is_front, double tolerance)
{
    double e0[3], e1[3], e2[3];
    SUB(e0, vertex_target_pos, triangle_target_pos_0);
    SUB(e1, triangle_target_pos_1, triangle_target_pos_0);
    SUB(e2, triangle_target_pos_2, triangle_target_pos_0);
    CROSS(normal, e1, e2);
    triangle_normal_magnitude_reciprocal = 1.0 / sqrt(DOT(normal, normal));
    MULTI_(normal, triangle_normal_magnitude_reciprocal);
    constraint = DOT(e0, normal);
    if (is_front) {
        constraint -= tolerance;
    }
    else {
        constraint += tolerance;
    }
}


void DCD::calPositionMove(double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
    double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* current_triangle_normal, double constraint, double tolerance, bool is_front, double triangle_normal_magnitude_reciprocal,
    double mass_point, double mass_t0, double mass_t1, double mass_t2, double& lambda)
{
    double in_triangle[3], scale_norm[3];
    SUB(in_triangle, current_position, current_triangle_position_0);
    MULTI(scale_norm, current_triangle_normal, constraint);
    SUB_(in_triangle, scale_norm);
    memcpy(vertex_target_pos, current_triangle_normal, 24);

    double temp_vec[3];
    SUB(temp_vec, current_triangle_position_1, current_triangle_position_2);
    MAGNITUDE_CROSS(triangle_target_pos_0, temp_vec, in_triangle, triangle_normal_magnitude_reciprocal);

    SUB(temp_vec, current_triangle_position_2, current_triangle_position_0);
    MAGNITUDE_CROSS(triangle_target_pos_1, temp_vec, in_triangle, triangle_normal_magnitude_reciprocal);

    SUB(temp_vec, current_triangle_position_0, current_triangle_position_1);
    MAGNITUDE_CROSS(triangle_target_pos_2, temp_vec, in_triangle, triangle_normal_magnitude_reciprocal);

    SUB_(triangle_target_pos_0, vertex_target_pos);

    double s = -constraint / (DOT(vertex_target_pos, vertex_target_pos) / mass_point + DOT(triangle_target_pos_0, triangle_target_pos_0) / mass_t0
        + DOT(triangle_target_pos_1, triangle_target_pos_1) / mass_t1 + DOT(triangle_target_pos_2, triangle_target_pos_2) / mass_t2);
    lambda = -s;

    double tem_value;

    tem_value = s / mass_point;
    MULTI_(vertex_target_pos, tem_value);
    tem_value = s / mass_t0;
    MULTI_(triangle_target_pos_0, tem_value);
    tem_value = s / mass_t1;
    MULTI_(triangle_target_pos_1, tem_value);
    tem_value = s / mass_t2;
    MULTI_(triangle_target_pos_2, tem_value);
}



void DCD::calDistancePointTriangle(double* vertex_target_pos, double* triangle_target_pos_0, double* triangle_target_pos_1, double* triangle_target_pos_2,
    double* current_position, double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* current_triangle_normal, double constraint, double tolerance, bool is_front, double triangle_normal_magnitude_reciprocal,
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

    double temp_vec[3];
    SUB(temp_vec, current_triangle_position_1, current_triangle_position_2);
    MAGNITUDE_CROSS(triangle_target_pos_0, temp_vec, in_triangle, triangle_normal_magnitude_reciprocal);

    SUB(temp_vec, current_triangle_position_2, current_triangle_position_0);
    MAGNITUDE_CROSS(triangle_target_pos_1, temp_vec, in_triangle, triangle_normal_magnitude_reciprocal);

    SUB(temp_vec, current_triangle_position_0, current_triangle_position_1);
    MAGNITUDE_CROSS(triangle_target_pos_2, temp_vec, in_triangle, triangle_normal_magnitude_reciprocal);

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
    double& lambda, double stiffness, double damping_stiffness, double dt, double& energy)
{
    energy = 0.0;
    double barycentric[4];
    CCD::internal::edgeEdgeDistanceType(initial_edge_vertex_0, initial_edge_vertex_1,
        initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, barycentric);
    //if (CCD::internal::edgeEdgeDistanceType(initial_edge_vertex_0, initial_edge_vertex_1,
    //    initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, barycentric) != 8) {
    //    return;
    //}
    double norm[3];
    double distance2;
    if (checkEdgeEdgeCollision(current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0, initial_edge_vertex_1, current_compare_edge_vertex_0, current_compare_edge_vertex_1,
        initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, barycentric, norm, distance2, tolerance)) {
        if (distance2 > 0) {
            return;
        }

        energy = 0.5 * stiffness * distance2 * distance2;

        XPBDcalDistanceEdgeEdge(
            norm, distance2, barycentric, current_edge_vertex_0, current_edge_vertex_1, current_compare_edge_vertex_0, current_compare_edge_vertex_1,
            initial_edge_vertex_0, initial_edge_vertex_1,
            initial_compare_edge_vertex_0, initial_compare_edge_vertex_1,
            //1.0,1.0,1.0,1.0,
            mass_inv_e_0_0, mass_inv_e_0_1, mass_inv_e_1_0, mass_inv_e_1_1,
            lambda, stiffness, damping_stiffness, dt);
    }
}
bool DCD::edgeEdgeCollider(double* edge_target_pos_0, double* edge_target_pos_1,
    double* current_edge_vertex_0, double* current_edge_vertex_1,
    double* initial_edge_vertex_0, double* initial_edge_vertex_1,
    double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
    double* initial_compare_edge_vertex_1,
    double tolerance, double mass_e_0_0, double mass_e_0_1)
{
    double barycentric[4];
    if (CCD::internal::edgeEdgeDistanceType(initial_edge_vertex_0, initial_edge_vertex_1,
        initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, barycentric) != 8) {
        // return false;
    }

    double norm[3];
    double distance2;
    if (checkEdgeEdgeCollision(current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0, initial_edge_vertex_1, current_compare_edge_vertex_0, current_compare_edge_vertex_1,
        initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, barycentric, norm, distance2, tolerance)) {
        calDistanceEdgeEdgeCollider(edge_target_pos_0, edge_target_pos_1,
            norm, distance2, barycentric, current_edge_vertex_0, current_edge_vertex_1,
            mass_e_0_0, mass_e_0_1);
        return true;
    }
    return false;
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
       // return false;
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
       // return false;
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
    //grad_c_p0 = (+-) alpha * norm and so on

    double e[3];
    //use e[3] to store x_i-x_ori
    double lambda_numerator = 0.0;
    SUB(e, current_edge_vertex_0, initial_edge_vertex_0);
    lambda_numerator -= alpha[0] * DOT(e, norm);
    SUB(e, current_edge_vertex_1, initial_edge_vertex_1);
    lambda_numerator -= alpha[1] * DOT(e, norm);
    SUB(e, current_compare_edge_vertex_0, initial_compare_edge_vertex_0);
    lambda_numerator += alpha[2] * DOT(e, norm);
    SUB(e, current_compare_edge_vertex_1, initial_compare_edge_vertex_1);
    lambda_numerator += alpha[3] * DOT(e, norm);


    // lambda
    double alpha_ = 1.0 / (stiffness * dt * dt);
    double gamma = damping_stiffness / (stiffness * dt);
    double delta_lambda = -(distance + alpha_ * lambda + gamma * lambda_numerator) / ((1 + gamma) * (mass_inv_e_0_0 * alpha[0] * alpha[0]
        + mass_inv_e_0_1 * alpha[1] * alpha[1] + mass_inv_e_1_0 * alpha[2] * alpha[2] +
        mass_inv_e_1_1 * alpha[3] * alpha[3]) + alpha_);
    lambda += delta_lambda;

    //double delta_lambda= -distance/(mass_inv_e_0_0 * alpha[0] * alpha[0]
    //    + mass_inv_e_0_1 * alpha[1] * alpha[1] + mass_inv_e_1_0 * alpha[2] * alpha[2] +
    //    mass_inv_e_1_1 * alpha[3] * alpha[3]);

    coe = -mass_inv_e_0_0 * delta_lambda * alpha[0];
    ACCUMULATE_SUM_WITH_COE(current_edge_vertex_0, coe, norm);
    coe = -mass_inv_e_0_1 * delta_lambda * alpha[1];
    ACCUMULATE_SUM_WITH_COE(current_edge_vertex_1, coe, norm);
    coe = mass_inv_e_1_0 * delta_lambda * alpha[2];
    ACCUMULATE_SUM_WITH_COE(current_compare_edge_vertex_0, coe, norm);
    coe = mass_inv_e_1_1 * delta_lambda * alpha[3];
    ACCUMULATE_SUM_WITH_COE(current_compare_edge_vertex_1, coe, norm);


}




void DCD::calDistanceEdgeEdgeCollider(double* edge_target_pos_0, double* edge_target_pos_1,
    double* norm, double distance, double* alpha, double* current_edge_vertex_0, double* current_edge_vertex_1,
    double mass_e_0_0, double mass_e_0_1)
{
    double coe = alpha[0] * alpha[0] / mass_e_0_0 + alpha[1] * alpha[1] / mass_e_0_1;
    double temp;
    distance /= coe;

    temp = distance * alpha[0] / mass_e_0_0;
    MULTI_SUM2(edge_target_pos_0, temp, norm, current_edge_vertex_0);

    temp = distance * alpha[1] / mass_e_0_1;
    MULTI_SUM2(edge_target_pos_1, temp, norm, current_edge_vertex_1);

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



void DCD::XPBDFloor(double* initial_position, double* current_position,unsigned int dimension, bool normal_direction, double mass_inv_v,
    double tolerance, double& lambda, double stiffness, double damping_stiffness, double dt, double floor_value, double& energy)
{
    energy = 0.0;
    double  constraint;
    if (normal_direction) {
        if (current_position[dimension] >= floor_value + tolerance) {
            return;
        }
        else {
            constraint = current_position[dimension] - floor_value - tolerance;
        }
    }
    else {
        if (current_position[dimension] <= floor_value - tolerance) {
            return;
        }
        else {
            constraint = floor_value - tolerance - current_position[dimension];
        }
    }
    // lambda
    double alpha_ = 1.0 / (stiffness * dt * dt);
    double gamma = damping_stiffness / (stiffness * dt);
    double coe_for_direction;

    if (normal_direction) {
        coe_for_direction = 1.0;
    }
    else {
        coe_for_direction = -1.0;
    }

    double e[3];
    SUB(e, current_position, initial_position);
    double delta_lambda = -(constraint + alpha_ * lambda + gamma * e[dimension]* coe_for_direction) /
        ((1 + gamma) * mass_inv_v + alpha_);
    lambda += delta_lambda;
    double coe = mass_inv_v * delta_lambda;
    current_position[dimension] += coe * coe_for_direction;
    energy = 0.5 * stiffness * constraint * constraint;
}


void DCD::XPBDpointTriangleCollider(double* initial_position, double* current_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double* current_triangle_normal, double mass_inv_v,
    double tolerance, double& lambda, double stiffness, double damping_stiffness, double dt, double& energy)
{
    energy = 0.0;
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
            current_triangle_normal, current_side, tolerance, lambda,stiffness,damping_stiffness,dt, energy);
    }
}


bool DCD::pointTriangleCollider(double* initial_position, double* current_position,
    double* initial_triangle_position_0, double* initial_triangle_position_1, double* initial_triangle_position_2,
    double* current_triangle_position_0, double* current_triangle_position_1, double* current_triangle_position_2,
    double* initial_triangle_normal, double* current_triangle_normal, double* vertex_target_pos,
    double tolerance)
{
    double barycentric[3];
    double c_p[3];

    double current_side;
    bool should_be_front;

    SUB(c_p, current_position, current_triangle_position_0);
    current_side = DOT(c_p, current_triangle_normal);  
    if (current_side > tolerance) {
        return false;
    }
    else {
        //if (!pointProjectOnTriangle(initial_position, initial_triangle_position_0, initial_triangle_position_1,
        //    initial_triangle_position_2, initial_triangle_normal, barycentric)) {
        //    return false;
        //}
        CCD::internal::pointTriangleNearestPoint(initial_position, initial_triangle_position_0, initial_triangle_position_1,
            initial_triangle_position_2, initial_triangle_normal, barycentric);

        if (checkIfCollidePointTriangleCollider(initial_position, current_position,
            initial_triangle_position_0, initial_triangle_position_1, initial_triangle_position_2,
            current_triangle_position_0, current_triangle_position_1, current_triangle_position_2,
            initial_triangle_normal, barycentric, tolerance)) {

            current_side -= tolerance;
            vertex_target_pos[0] = current_position[0] - current_side * current_triangle_normal[0];
            vertex_target_pos[1] = current_position[1] - current_side * current_triangle_normal[1];
            vertex_target_pos[2] = current_position[2] - current_side * current_triangle_normal[2];

            //double in_triangle[3], scale_norm[3];
            //SUB(in_triangle, current_position, current_triangle_position_0);
            //MULTI(scale_norm, current_triangle_normal, current_side);
            //SUB_(in_triangle, scale_norm);
            //if (should_be_front) {
            //    current_side -= tolerance;
            //}
            //else {
            //    current_side += tolerance;
            //}
            //vertex_target_pos[0] = current_position[0] - current_side * current_triangle_normal[0];
            //vertex_target_pos[1] = current_position[1] - current_side * current_triangle_normal[1];
            //vertex_target_pos[2] = current_position[2] - current_side * current_triangle_normal[2];


            return true;
        }
    }
 
    return false;
}


bool DCD::PDFloor(double* target_position, double* current_position, unsigned int dimension, bool normal_direction,
    double tolerance, double floor_value)
{
    if (normal_direction) {
        if (current_position[dimension] > floor_value + tolerance) {
            return false;
        }
        else {
            memcpy(target_position, current_position, 24);
            target_position[dimension] = floor_value + tolerance;
            return true;
        }
    }
    else {
        if (current_position[dimension] < floor_value - tolerance) {
            return false;
        }
        else {
            memcpy(target_position, current_position, 24);
            target_position[dimension] = floor_value - tolerance;
            return true;
        }
    }
}





void DCD::XPBDcalDistancePointTriangleCollider(double* initial_position,
    double* current_position,double mass_inv_vertex,
    double* current_triangle_normal, double constraint, double tolerance,
    double& lambda, double stiffness, double damping_stiffness, double dt, double& energy)
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
    energy = 0.5 * stiffness * constraint * constraint;
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



