#pragma once
#include"basic/global.h"
#include"../external/Eigen/Dense"
#include"../basic/eigenDenseOperation.h"
#include<array>
#include <iostream>


void d2V_spring_particle_particle_dq2(Eigen::MatrixXd& H, double* q0, double* q1, double l0, double stiffness) {

    Vector3d Ax;
    SUB(Ax, q0, q1);
    VectorXd ATAx(6);
    ATAx.block(0, 0, 3, 1) = Ax;
    ATAx.block(3, 0, 3, 1) = -Ax;

    double coe1 = stiffness * l0 / (Ax.norm() * Ax.squaredNorm());
    H = coe1 * ATAx * ATAx.transpose();
    double coe0 = stiffness - stiffness * l0 / Ax.norm();
    H.block(0, 0, 3, 3) += coe0 * Matrix3d::Identity();
    H.block(0, 3, 3, 3) -= coe0 * Matrix3d::Identity();
    H.block(3, 0, 3, 3) -= coe0 * Matrix3d::Identity();
    H.block(3, 3, 3, 3) += coe0 * Matrix3d::Identity();


    // Test for different output from equation
    // Both works
    //Eigen::Matrix66d test;
   /* H(0, 0) = -stiffness * (-pow(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2], 3.0 / 2.0) + \
        l0 * (q0[1] * q0[1]) + l0 * (q0[2] * q0[2]) + l0 * (q1[1] * q1[1]) + l0 * (q1[2] * q1[2]) - l0 * q0[1] * q1[1] * 2.0 - l0 * q0[2] * q1[2] * 2.0) * \
        1.0 / pow(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2], 3.0 / 2.0);
    H(0, 1) = stiffness * l0 * (q0[0] - q1[0]) * (q0[1] - q1[1]) * 1.0 / pow(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2], 3.0 / 2.0);
    H(0, 2) = stiffness * l0 * (q0[0] - q1[0]) * (q0[2] - q1[2]) * 1.0 / pow(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2], 3.0 / 2.0);
    H(0,3) = stiffness*(-pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0)+l0*(q0[1]*q0[1])+l0*(q0[2]*q0[2])+l0*(q1[1]*q1[1])+l0*(q1[2]*q1[2])-l0*q0[1]*q1[1]*2.0-l0*q0[2]*q1[2]*2.0)*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(0,4) = -stiffness*l0*(q0[0]-q1[0])*(q0[1]-q1[1])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(0,5) = -stiffness*l0*(q0[0]-q1[0])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(1,0) = stiffness*l0*(q0[0]-q1[0])*(q0[1]-q1[1])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(1,1) = -stiffness*(-pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0)+l0*(q0[0]*q0[0])+l0*(q0[2]*q0[2])+l0*(q1[0]*q1[0])+l0*(q1[2]*q1[2])-l0*q0[0]*q1[0]*2.0-l0*q0[2]*q1[2]*2.0)*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(1,2) = stiffness*l0*(q0[1]-q1[1])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(1,3) = -stiffness*l0*(q0[0]-q1[0])*(q0[1]-q1[1])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(1,4) = stiffness*(-pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0)+l0*(q0[0]*q0[0])+l0*(q0[2]*q0[2])+l0*(q1[0]*q1[0])+l0*(q1[2]*q1[2])-l0*q0[0]*q1[0]*2.0-l0*q0[2]*q1[2]*2.0)*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(1,5) = -stiffness*l0*(q0[1]-q1[1])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(2,0) = stiffness*l0*(q0[0]-q1[0])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(2,1) = stiffness*l0*(q0[1]-q1[1])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(2,2) = -stiffness*(-pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0)+l0*(q0[0]*q0[0])+l0*(q0[1]*q0[1])+l0*(q1[0]*q1[0])+l0*(q1[1]*q1[1])-l0*q0[0]*q1[0]*2.0-l0*q0[1]*q1[1]*2.0)*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(2,3) = -stiffness*l0*(q0[0]-q1[0])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(2,4) = -stiffness*l0*(q0[1]-q1[1])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(2,5) = stiffness*(-pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0)+l0*(q0[0]*q0[0])+l0*(q0[1]*q0[1])+l0*(q1[0]*q1[0])+l0*(q1[1]*q1[1])-l0*q0[0]*q1[0]*2.0-l0*q0[1]*q1[1]*2.0)*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(3,0) = stiffness*(-pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0)+l0*(q0[1]*q0[1])+l0*(q0[2]*q0[2])+l0*(q1[1]*q1[1])+l0*(q1[2]*q1[2])-l0*q0[1]*q1[1]*2.0-l0*q0[2]*q1[2]*2.0)*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(3,1) = -stiffness*l0*(q0[0]-q1[0])*(q0[1]-q1[1])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(3,2) = -stiffness*l0*(q0[0]-q1[0])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(3,3) = -stiffness*(-pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0)+l0*(q0[1]*q0[1])+l0*(q0[2]*q0[2])+l0*(q1[1]*q1[1])+l0*(q1[2]*q1[2])-l0*q0[1]*q1[1]*2.0-l0*q0[2]*q1[2]*2.0)*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(3,4) = stiffness*l0*(q0[0]-q1[0])*(q0[1]-q1[1])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(3,5) = stiffness*l0*(q0[0]-q1[0])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(4,0) = -stiffness*l0*(q0[0]-q1[0])*(q0[1]-q1[1])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(4,1) = stiffness*(-pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0)+l0*(q0[0]*q0[0])+l0*(q0[2]*q0[2])+l0*(q1[0]*q1[0])+l0*(q1[2]*q1[2])-l0*q0[0]*q1[0]*2.0-l0*q0[2]*q1[2]*2.0)*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(4,2) = -stiffness*l0*(q0[1]-q1[1])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(4,3) = stiffness*l0*(q0[0]-q1[0])*(q0[1]-q1[1])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(4,4) = -stiffness*(-pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0)+l0*(q0[0]*q0[0])+l0*(q0[2]*q0[2])+l0*(q1[0]*q1[0])+l0*(q1[2]*q1[2])-l0*q0[0]*q1[0]*2.0-l0*q0[2]*q1[2]*2.0)*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(4,5) = stiffness*l0*(q0[1]-q1[1])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(5,0) = -stiffness*l0*(q0[0]-q1[0])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(5,1) = -stiffness*l0*(q0[1]-q1[1])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(5,2) = stiffness*(-pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0)+l0*(q0[0]*q0[0])+l0*(q0[1]*q0[1])+l0*(q1[0]*q1[0])+l0*(q1[1]*q1[1])-l0*q0[0]*q1[0]*2.0-l0*q0[1]*q1[1]*2.0)*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(5,3) = stiffness*l0*(q0[0]-q1[0])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(5,4) = stiffness*l0*(q0[1]-q1[1])*(q0[2]-q1[2])*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);
    H(5,5) = -stiffness*(-pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0)+l0*(q0[0]*q0[0])+l0*(q0[1]*q0[1])+l0*(q1[0]*q1[0])+l0*(q1[1]*q1[1])-l0*q0[0]*q1[0]*2.0-l0*q0[1]*q1[1]*2.0)*1.0/pow(q0[0]*q1[0]*-2.0-q0[1]*q1[1]*2.0-q0[2]*q1[2]*2.0+q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q1[0]*q1[0]+q1[1]*q1[1]+q1[2]*q1[2],3.0/2.0);*/

    /*   H(0, 0) = -stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) *
           (l0 * pow((q0[1] - q1[1]), 2.0) + l0 * pow((q0[2] - q1[2]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));
       H(0, 1) = stiffness * l0 *(q0[0] - q1[0]) *(q0[1] - q1[1]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(0, 2) = stiffness * l0 *(q0[0] - q1[0]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(0, 3) = stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) * (l0 * pow((q0[1] - q1[1]), 2.0) + l0 * pow((q0[2] - q1[2]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));
       H(0, 4) = -stiffness * l0 *(q0[0] - q1[0]) *(q0[1] - q1[1]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(0, 5) = -stiffness * l0 *(q0[0] - q1[0]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(1, 0) = stiffness * l0 *(q0[0] - q1[0]) *(q0[1] - q1[1]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(1, 1) = -stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) * (l0 * pow((q0[0] - q1[0]), 2.0) + l0 * pow((q0[2] - q1[2]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));
       H(1, 2) = stiffness * l0 *(q0[1] - q1[1]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(1, 3) = -stiffness * l0 *(q0[0] - q1[0]) *(q0[1] - q1[1]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(1, 4) = stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) * (l0 * pow((q0[0] - q1[0]), 2.0) + l0 * pow((q0[2] - q1[2]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));
       H(1, 5) = -stiffness * l0 *(q0[1] - q1[1]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(2, 0) = stiffness * l0 *(q0[0] - q1[0]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(2, 1) = stiffness * l0 *(q0[1] - q1[1]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(2, 2) = -stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) * (l0 * pow((q0[0] - q1[0]), 2.0) + l0 * pow((q0[1] - q1[1]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));
       H(2, 3) = -stiffness * l0 *(q0[0] - q1[0]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(2, 4) = -stiffness * l0 *(q0[1] - q1[1]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(2, 5) = stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) * (l0 * pow((q0[0] - q1[0]), 2.0) + l0 * pow((q0[1] - q1[1]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));
       H(3, 0) = stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) * (l0 * pow((q0[1] - q1[1]), 2.0) + l0 * pow((q0[2] - q1[2]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));
       H(3, 1) = -stiffness * l0 *(q0[0] - q1[0]) *(q0[1] - q1[1]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(3, 2) = -stiffness * l0 *(q0[0] - q1[0]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(3, 3) = -stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) * (l0 * pow((q0[1] - q1[1]), 2.0) + l0 * pow((q0[2] - q1[2]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));
       H(3, 4) = stiffness * l0 *(q0[0] - q1[0]) *(q0[1] - q1[1]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(3, 5) = stiffness * l0 *(q0[0] - q1[0]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(4, 0) = -stiffness * l0 *(q0[0] - q1[0]) *(q0[1] - q1[1]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(4, 1) = stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) * (l0 * pow((q0[0] - q1[0]), 2.0) + l0 * pow((q0[2] - q1[2]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));
       H(4, 2) = -stiffness * l0 *(q0[1] - q1[1]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(4, 3) = stiffness * l0 *(q0[0] - q1[0]) *(q0[1] - q1[1]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(4, 4) = -stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) * (l0 * pow((q0[0] - q1[0]), 2.0) + l0 * pow((q0[2] - q1[2]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));
       H(4, 5) = stiffness * l0 *(q0[1] - q1[1]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(5, 0) = -stiffness * l0 *(q0[0] - q1[0]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(5, 1) = -stiffness * l0 *(q0[1] - q1[1]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(5, 2) = stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) * (l0 * pow((q0[0] - q1[0]), 2.0) + l0 * pow((q0[1] - q1[1]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));
       H(5, 3) = stiffness * l0 *(q0[0] - q1[0]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(5, 4) = stiffness * l0 *(q0[1] - q1[1]) *(q0[2] - q1[2]) *   1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0);
       H(5, 5) = -stiffness * 1.0 / pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0) * (l0 * pow((q0[0] - q1[0]), 2.0) + l0 * pow((q0[1] - q1[1]), 2.0) - pow(pow((q0[0] - q1[0]), 2.0) + pow((q0[1] - q1[1]), 2.0) + pow((q0[2] - q1[2]), 2.0), 3.0 / 2.0));*/



       //std::cout << "test: \n" << test << std::endl;
       //std::cout << "H: \n" << H << std::endl;

}



void setForce(VectorXd& f, double* q0, double* q1, double l0, double stiffness)
{
    //Vector3d Ax;
    //SUB(Ax, q0, q1);
    //VectorXd ATAx(6);
    //ATAx.block(0, 0, 3, 1) = Ax;
    //ATAx.block(3, 0, 3, 1) = -Ax;
    //double coe = stiffness * (l0 - Ax.norm()) / Ax.norm();
    //f = ATAx * coe;

    f[0] = stiffness * (l0 - sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2])) * (q0[0] - q1[0]) * 1.0 / sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2]);
    f[1] = stiffness * (l0 - sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2])) * (q0[1] - q1[1]) * 1.0 / sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2]);
    f[2] = stiffness * (l0 - sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2])) * (q0[2] - q1[2]) * 1.0 / sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2]);
    f[3] = -(stiffness * (q0[0] * 2.0 - q1[0] * 2.0) * (l0 - sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2])) * 1.0 / sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2])) / 2.0;
    f[4] = -(stiffness * (q0[1] * 2.0 - q1[1] * 2.0) * (l0 - sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2])) * 1.0 / sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2])) / 2.0;
    f[5] = -(stiffness * (q0[2] * 2.0 - q1[2] * 2.0) * (l0 - sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2])) * 1.0 / sqrt(q0[0] * q1[0] * -2.0 - q0[1] * q1[1] * 2.0 - q0[2] * q1[2] * 2.0 + q0[0] * q0[0] + q0[1] * q0[1] + q0[2] * q0[2] + q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2])) / 2.0;
}

void assemble_force(VectorXd& force, double* q0, double* q1, double l0, double stiffness,
    unsigned int vertex_index_0, unsigned int vertex_index_1, bool all_unfixed_vertex)
{
    VectorXd f(6);
    setForce(f, q0, q1, l0, stiffness);
    if (all_unfixed_vertex) {
        force.block(3 * vertex_index_0, 0, 3, 1) += f.block(0, 0, 3, 1);
        force.block(3 * vertex_index_1, 0, 3, 1) += f.block(3, 0, 3, 1);
    }
    else {
        force.block(3 * vertex_index_0, 0, 3, 1) += f.block(0, 0, 3, 1);
    }
}

void assemble_stiffness(std::vector<Triplet<double>>& tripleList, double* q0, double* q1, double l0, double stiffness,
    unsigned int vertex_index_0, unsigned int vertex_index_1, bool all_unfixed_vertex)
{
    MatrixXd H_i;
    H_i.resize(6, 6);
    d2V_spring_particle_particle_dq2(H_i, q0, q1, l0, stiffness);

    // std::cout << H_i << std::endl;

    Eigen::Matrix3d H_iAA = H_i.block(0, 0, 3, 3);
    if (all_unfixed_vertex) {
        Eigen::Matrix3d H_iAB = H_i.block(0, 3, 3, 3);
        Eigen::Matrix3d H_iABT = H_i.block(3, 0, 3, 3);
        Eigen::Matrix3d H_iBB = H_i.block(3, 3, 3, 3);
        for (int ii = 0; ii < 3; ii++)
        {
            for (int jj = 0; jj < 3; jj++)
            {
                tripleList.push_back(Triplet<double>(3 * vertex_index_0 + ii, 3 * vertex_index_0 + jj, H_iAA.coeff(ii, jj)));
                tripleList.push_back(Triplet<double>(3 * vertex_index_1 + ii, 3 * vertex_index_1 + jj, H_iBB.coeff(ii, jj)));
                tripleList.push_back(Triplet<double>(3 * vertex_index_0 + ii, 3 * vertex_index_1 + jj, H_iAB.coeff(ii, jj)));
                tripleList.push_back(Triplet<double>(3 * vertex_index_1 + ii, 3 * vertex_index_0 + jj, H_iABT.coeff(ii, jj)));
            }
        }
    }
    else {
        for (int ii = 0; ii < 3; ii++)
        {
            for (int jj = 0; jj < 3; jj++)
            {
                tripleList.push_back(Triplet<double>(3 * vertex_index_0 + ii, 3 * vertex_index_0 + jj, H_iAA.coeff(ii, jj)));
            }
        }
    }

}


void setForce(VectorXd& force, std::array<double, 3>* vertex_position, unsigned int* edge_vertices_mass_spring,
    unsigned int* only_one_vertex_fix_edge_vertices, double stiffness, double* unfixed_rest_length,
    double* fixed_one_vertices_rest_length, unsigned int* unfixed_vertex,
    unsigned int unfixed_edge_num, unsigned int only_one_vertex_fixed_edge_num, unsigned int vertex_start)
{
    unfixed_edge_num <<= 1;
    for (unsigned int i = 0; i < unfixed_edge_num; i += 2) {
        assemble_force(force, vertex_position[unfixed_vertex[edge_vertices_mass_spring[i]]].data(), vertex_position[unfixed_vertex[edge_vertices_mass_spring[i + 1]]].data(),
            unfixed_rest_length[i >> 1], stiffness, vertex_start + edge_vertices_mass_spring[i], vertex_start + edge_vertices_mass_spring[i + 1], true);
    }
    // std::cout << "end" << std::endl;
    for (unsigned int i = 0; i < only_one_vertex_fixed_edge_num; i += 2) {
        assemble_force(force, vertex_position[unfixed_vertex[only_one_vertex_fix_edge_vertices[i]]].data(), vertex_position[only_one_vertex_fix_edge_vertices[i + 1]].data(),
            fixed_one_vertices_rest_length[i >> 1], stiffness, vertex_start + only_one_vertex_fix_edge_vertices[i], vertex_start + only_one_vertex_fix_edge_vertices[i], false);
    }
}


void setMatrix(std::vector<Triplet<double>>& tripleList, std::array<double, 3>* vertex_position, unsigned int* edge_vertices_mass_spring,
    unsigned int* only_one_vertex_fix_edge_vertices, double stiffness, double* unfixed_rest_length,
    double* fixed_one_vertices_rest_length, unsigned int* unfixed_vertex,
    unsigned int unfixed_edge_num, unsigned int only_one_vertex_fixed_edge_num, unsigned int vertex_start)
{
    unfixed_edge_num <<= 1;
    for (unsigned int i = 0; i < unfixed_edge_num; i += 2) {
        assemble_stiffness(tripleList, vertex_position[unfixed_vertex[edge_vertices_mass_spring[i]]].data(), vertex_position[unfixed_vertex[edge_vertices_mass_spring[i + 1]]].data(),
            unfixed_rest_length[i >> 1], stiffness, vertex_start + edge_vertices_mass_spring[i], vertex_start + edge_vertices_mass_spring[i + 1], true);
    }
    // std::cout << "end" << std::endl;
    for (unsigned int i = 0; i < only_one_vertex_fixed_edge_num; i += 2) {
        assemble_stiffness(tripleList, vertex_position[unfixed_vertex[only_one_vertex_fix_edge_vertices[i]]].data(), vertex_position[only_one_vertex_fix_edge_vertices[i + 1]].data(),
            fixed_one_vertices_rest_length[i >> 1], stiffness, vertex_start + only_one_vertex_fix_edge_vertices[i], vertex_start + only_one_vertex_fix_edge_vertices[i], false);
    }
    //std::cout << "end_2" << std::endl;
}




void setForce(VectorXd& Sn, VectorXd& force, std::vector<std::array<double, 3>*>& vertex_position, std::vector<std::vector<unsigned int>*>& edge_vertices_mass_spring,
    std::vector<std::vector<unsigned int>*>& only_one_vertex_fix_edge_vertices, double stiffness, std::vector<std::vector<double>*>& unfixed_rest_length,
    std::vector<std::vector<double>*>& fixed_one_vertices_rest_length, std::vector<std::vector<unsigned int>*>& unfixed_vertex,
    std::vector<double*>& mass, double time_step)
{
    unsigned int total_num = 0;
    for (unsigned int i = 0; i < unfixed_vertex.size(); ++i) {
        total_num += unfixed_vertex[i]->size();
    }
    force.resize(3 * total_num);
    force.setZero();
    unsigned int vertex_start = 0;
    for (unsigned int i = 0; i < vertex_position.size(); ++i) {
        setForce(force, vertex_position[i], edge_vertices_mass_spring[i]->data(), only_one_vertex_fix_edge_vertices[i]->data(),
            stiffness, unfixed_rest_length[i]->data(), fixed_one_vertices_rest_length[i]->data(), unfixed_vertex[i]->data(),
            unfixed_rest_length[i]->size(), only_one_vertex_fix_edge_vertices[i]->size(), vertex_start);
        vertex_start += unfixed_vertex[i]->size();
    }

    unsigned int* unfixed_vertex_obj;
    vertex_start = 0;
    time_step *= time_step;
    std::array<double, 3>* vertex_pos;
    for (unsigned int i = 0; i < mass.size(); ++i) {
        unfixed_vertex_obj = unfixed_vertex[i]->data();
        vertex_pos = vertex_position[i];
        for (unsigned int j = 0; j < unfixed_vertex[i]->size(); ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                force[3 * (vertex_start + j) + k] = time_step * force[3 * (vertex_start + j) + k] + mass[i][unfixed_vertex_obj[j]] * (Sn[3 * (vertex_start + j) + k] -
                    vertex_pos[unfixed_vertex_obj[j]][k]);
            }
        }
    }

}

void setMatrix(Eigen::SparseMatrix<double>& K, std::vector<std::array<double, 3>*>& vertex_position, std::vector<std::vector<unsigned int>*>& edge_vertices_mass_spring,
    std::vector<std::vector<unsigned int>*>& only_one_vertex_fix_edge_vertices, double stiffness, std::vector<std::vector<double>*>& unfixed_rest_length,
    std::vector<std::vector<double>*>& fixed_one_vertices_rest_length, std::vector<std::vector<unsigned int>*>& unfixed_vertex,
    std::vector<double*>& mass, double time_step)
{
    std::vector<Triplet<double>> tripleList;
    tripleList.reserve(36 * unfixed_rest_length[0]->size());
    unsigned int vertex_start = 0;
    for (unsigned int i = 0; i < vertex_position.size(); ++i) {
        setMatrix(tripleList, vertex_position[i], edge_vertices_mass_spring[i]->data(), only_one_vertex_fix_edge_vertices[i]->data(),
            stiffness, unfixed_rest_length[i]->data(), fixed_one_vertices_rest_length[i]->data(), unfixed_vertex[i]->data(),
            unfixed_rest_length[i]->size(), only_one_vertex_fix_edge_vertices[i]->size(), vertex_start);
        vertex_start += unfixed_vertex[i]->size();
    }
    K.resize(3 * vertex_start, 3 * vertex_start);
    K.setFromTriplets(tripleList.begin(), tripleList.end());
    K *= (time_step * time_step);
    vertex_start = 0;
    unsigned int* unfixed_vertex_obj;
    for (unsigned int i = 0; i < mass.size(); ++i) {
        unfixed_vertex_obj = unfixed_vertex[i]->data();
        for (unsigned int j = 0; j < unfixed_vertex[i]->size(); ++j) {
            K.coeffRef(3 * (vertex_start + j), 3 * (vertex_start + j)) += mass[i][unfixed_vertex_obj[j]];
            K.coeffRef(3 * (vertex_start + j) + 1, 3 * (vertex_start + j) + 1) += mass[i][unfixed_vertex_obj[j]];
            K.coeffRef(3 * (vertex_start + j) + 2, 3 * (vertex_start + j) + 2) += mass[i][unfixed_vertex_obj[j]];
        }
        vertex_start += unfixed_vertex[i]->size();
    }

}


void compareTwoForce(VectorXd& Sn, VectorXd& ori, std::vector<std::array<double, 3>*>& vertex_position, std::vector<std::vector<unsigned int>*>& edge_vertices_mass_spring,
    std::vector<std::vector<unsigned int>*>& only_one_vertex_fix_edge_vertices, double stiffness, std::vector<std::vector<double>*>& unfixed_rest_length,
    std::vector<std::vector<double>*>& fixed_one_vertices_rest_length, std::vector<std::vector<unsigned int>*>& unfixed_vertex,
    std::vector<double*>& mass, double time_step)
{
    VectorXd force;
    setForce(Sn, force, vertex_position, edge_vertices_mass_spring, only_one_vertex_fix_edge_vertices, stiffness, unfixed_rest_length, fixed_one_vertices_rest_length,
        unfixed_vertex, mass, time_step);
    std::cout << "force " << (force - ori).norm() << std::endl;
    if ((force - ori).norm() > 1e-7) {
        std::cout << "force error " << std::endl;
    }

    //std::cout << force << std::endl;
    //std::cout << "==== " << std::endl;
    //std::cout << ori << std::endl;
}


void compareTwoMatrix(Eigen::SparseMatrix<double>& ori, std::vector<std::array<double, 3>*>& vertex_position, std::vector<std::vector<unsigned int>*>& edge_vertices_mass_spring,
    std::vector<std::vector<unsigned int>*>& only_one_vertex_fix_edge_vertices, double stiffness, std::vector<std::vector<double>*>& unfixed_rest_length,
    std::vector<std::vector<double>*>& fixed_one_vertices_rest_length, std::vector<std::vector<unsigned int>*>& unfixed_vertex,
    std::vector<double*>& mass, double time_step)
{
    Eigen::SparseMatrix<double> K;
    setMatrix(K, vertex_position, edge_vertices_mass_spring, only_one_vertex_fix_edge_vertices, stiffness, unfixed_rest_length, fixed_one_vertices_rest_length,
        unfixed_vertex, mass, time_step);
    std::cout << "matrix " << (K - ori).norm() << std::endl;
    if ((K - ori).norm() > 1e-7) {
        std::cout << "matrix error " << std::endl;
    }
    //std::cout << K << std::endl;
    //std::cout << "==== "  << std::endl;
    //std::cout << ori << std::endl;
}


