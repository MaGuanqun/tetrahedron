#pragma once

#include"XPBD/FEM_relate.h"

using namespace Eigen;

namespace ConstitutiveModel {

    inline void first_piola_derivative(const Matrix3d& U, const Vector3d& sigma, const Matrix3d& V, const Vector3d& de_dsigma,
        const Vector3d& B_left_coeff, const Matrix3d& _d2e_dsigma2, MatrixXd& dPdF)
    {

        Matrix3d d2e_dsigma2 = _d2e_dsigma2;
        FEM::SPDprojection(d2e_dsigma2);

        double left_coef = B_left_coeff[0];
        double right_coef = de_dsigma[0] + de_dsigma[1];
        double sum_sigma = (std::max)(sigma[0] + sigma[1], 0.000001);
        right_coef /= 2.0 * sum_sigma;
        MatrixXd B0(2, 2);
        B0 << left_coef + right_coef, left_coef - right_coef, left_coef - right_coef, left_coef + right_coef;

        FEM::SPDprojection(B0);

        left_coef = B_left_coeff[1];
        right_coef = de_dsigma[1] + de_dsigma[2];
        sum_sigma = (std::max)(sigma[1] + sigma[2], 0.000001);
        right_coef /= 2.0 * sum_sigma;
        MatrixXd B1(2, 2);
        B1 << left_coef + right_coef, left_coef - right_coef, left_coef - right_coef, left_coef + right_coef;
        FEM::SPDprojection(B1);

        left_coef = B_left_coeff[2];
        right_coef = de_dsigma[2] + de_dsigma[0];
        sum_sigma = (std::max)(sigma[2] + sigma[0], 0.000001);
        right_coef /= 2.0 * sum_sigma;
        MatrixXd B2(2, 2);
        B2 << left_coef + right_coef, left_coef - right_coef, left_coef - right_coef, left_coef + right_coef;
        FEM::SPDprojection(B2);

        Matrix<double, 9, 9> M;
        M.setZero();
        M(0, 0) = d2e_dsigma2(0, 0);
        M(0, 4) = d2e_dsigma2(0, 1);
        M(0, 8) = d2e_dsigma2(0, 2);
        M(4, 0) = d2e_dsigma2(1, 0);
        M(4, 4) = d2e_dsigma2(1, 1);
        M(4, 8) = d2e_dsigma2(1, 2);
        M(8, 0) = d2e_dsigma2(2, 0);
        M(8, 4) = d2e_dsigma2(2, 1);
        M(8, 8) = d2e_dsigma2(2, 2);
        M(1, 1) = B0(0, 0);
        M(1, 3) = B0(0, 1);
        M(3, 1) = B0(1, 0);
        M(3, 3) = B0(1, 1);
        M(5, 5) = B1(0, 0);
        M(5, 7) = B1(0, 1);
        M(7, 5) = B1(1, 0);
        M(7, 7) = B1(1, 1);
        M(2, 2) = B2(1, 1);
        M(2, 6) = B2(1, 0);
        M(6, 2) = B2(0, 1);
        M(6, 6) = B2(0, 0);

        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                for (int s = 0; s < 3; ++s)
                    for (int r = 0; r < 3; ++r) {
                        int ij = j * 3 + i;
                        int rs = s * 3 + r;
                        dPdF(ij, rs) = M(0, 0) * U(i, 0) * V(j, 0) * U(r, 0) * V(s, 0)
                            + M(0, 4) * U(i, 0) * V(j, 0) * U(r, 1) * V(s, 1)
                            + M(0, 8) * U(i, 0) * V(j, 0) * U(r, 2) * V(s, 2)
                            + M(4, 0) * U(i, 1) * V(j, 1) * U(r, 0) * V(s, 0)
                            + M(4, 4) * U(i, 1) * V(j, 1) * U(r, 1) * V(s, 1)
                            + M(4, 8) * U(i, 1) * V(j, 1) * U(r, 2) * V(s, 2)
                            + M(8, 0) * U(i, 2) * V(j, 2) * U(r, 0) * V(s, 0)
                            + M(8, 4) * U(i, 2) * V(j, 2) * U(r, 1) * V(s, 1)
                            + M(8, 8) * U(i, 2) * V(j, 2) * U(r, 2) * V(s, 2)
                            + M(1, 1) * U(i, 0) * V(j, 1) * U(r, 0) * V(s, 1)
                            + M(1, 3) * U(i, 0) * V(j, 1) * U(r, 1) * V(s, 0)
                            + M(3, 1) * U(i, 1) * V(j, 0) * U(r, 0) * V(s, 1)
                            + M(3, 3) * U(i, 1) * V(j, 0) * U(r, 1) * V(s, 0)
                            + M(5, 5) * U(i, 1) * V(j, 2) * U(r, 1) * V(s, 2)
                            + M(5, 7) * U(i, 1) * V(j, 2) * U(r, 2) * V(s, 1)
                            + M(7, 5) * U(i, 2) * V(j, 1) * U(r, 1) * V(s, 2)
                            + M(7, 7) * U(i, 2) * V(j, 1) * U(r, 2) * V(s, 1)
                            + M(2, 2) * U(i, 0) * V(j, 2) * U(r, 0) * V(s, 2)
                            + M(2, 6) * U(i, 0) * V(j, 2) * U(r, 2) * V(s, 0)
                            + M(6, 2) * U(i, 2) * V(j, 0) * U(r, 0) * V(s, 2)
                            + M(6, 6) * U(i, 2) * V(j, 0) * U(r, 2) * V(s, 0);
                    }
    }


    inline void backpropagate_element_gradient(const Matrix3d& IB, const Matrix3d& de_dF, VectorXd& de_dX)
    {
        de_dX.setZero();

        double R10 = IB(0, 0) * de_dF(0, 0) + IB(0, 1) * de_dF(0, 1) + IB(0, 2) * de_dF(0, 2);
        double R11 = IB(0, 0) * de_dF(1, 0) + IB(0, 1) * de_dF(1, 1) + IB(0, 2) * de_dF(1, 2);
        double R12 = IB(0, 0) * de_dF(2, 0) + IB(0, 1) * de_dF(2, 1) + IB(0, 2) * de_dF(2, 2);
        double R20 = IB(1, 0) * de_dF(0, 0) + IB(1, 1) * de_dF(0, 1) + IB(1, 2) * de_dF(0, 2);
        double R21 = IB(1, 0) * de_dF(1, 0) + IB(1, 1) * de_dF(1, 1) + IB(1, 2) * de_dF(1, 2);
        double R22 = IB(1, 0) * de_dF(2, 0) + IB(1, 1) * de_dF(2, 1) + IB(1, 2) * de_dF(2, 2);
        double R30 = IB(2, 0) * de_dF(0, 0) + IB(2, 1) * de_dF(0, 1) + IB(2, 2) * de_dF(0, 2);
        double R31 = IB(2, 0) * de_dF(1, 0) + IB(2, 1) * de_dF(1, 1) + IB(2, 2) * de_dF(1, 2);
        double R32 = IB(2, 0) * de_dF(2, 0) + IB(2, 1) * de_dF(2, 1) + IB(2, 2) * de_dF(2, 2);
        de_dX(1 * 3 + 0) = R10;
        de_dX(1 * 3 + 1) = R11;
        de_dX(1 * 3 + 2) = R12;
        de_dX(2 * 3 + 0) = R20;
        de_dX(2 * 3 + 1) = R21;
        de_dX(2 * 3 + 2) = R22;
        de_dX(3 * 3 + 0) = R30;
        de_dX(3 * 3 + 1) = R31;
        de_dX(3 * 3 + 2) = R32;
        de_dX(0 * 3 + 0) = -R10 - R20 - R30;
        de_dX(0 * 3 + 1) = -R11 - R21 - R31;
        de_dX(0 * 3 + 2) = -R12 - R22 - R32;
    }


    inline void backpropagate_element_hessian(const Matrix3d& IB, const Eigen::MatrixXd& d2e_dF2, MatrixXd& d2e_dX2)
    {
        d2e_dX2.setZero();
        Matrix<double, 12, 9>intermediate;
        intermediate.setZero();      
        for (int colI = 0; colI < 9; ++colI) {
            intermediate(3, colI) = IB(0, 0) * d2e_dF2(0, colI) + IB(0, 1) * d2e_dF2(3, colI) + IB(0, 2) * d2e_dF2(6, colI);
            intermediate(4, colI) = IB(0, 0) * d2e_dF2(1, colI) + IB(0, 1) * d2e_dF2(4, colI) + IB(0, 2) * d2e_dF2(7, colI);
            intermediate(5, colI) = IB(0, 0) * d2e_dF2(2, colI) + IB(0, 1) * d2e_dF2(5, colI) + IB(0, 2) * d2e_dF2(8, colI);
            intermediate(6, colI) = IB(1, 0) * d2e_dF2(0, colI) + IB(1, 1) * d2e_dF2(3, colI) + IB(1, 2) * d2e_dF2(6, colI);
            intermediate(7, colI) = IB(1, 0) * d2e_dF2(1, colI) + IB(1, 1) * d2e_dF2(4, colI) + IB(1, 2) * d2e_dF2(7, colI);
            intermediate(8, colI) = IB(1, 0) * d2e_dF2(2, colI) + IB(1, 1) * d2e_dF2(5, colI) + IB(1, 2) * d2e_dF2(8, colI);
            intermediate(9, colI) = IB(2, 0) * d2e_dF2(0, colI) + IB(2, 1) * d2e_dF2(3, colI) + IB(2, 2) * d2e_dF2(6, colI);
            intermediate(10, colI) = IB(2, 0) * d2e_dF2(1, colI) + IB(2, 1) * d2e_dF2(4, colI) + IB(2, 2) * d2e_dF2(7, colI);
            intermediate(11, colI) = IB(2, 0) * d2e_dF2(2, colI) + IB(2, 1) * d2e_dF2(5, colI) + IB(2, 2) * d2e_dF2(8, colI);
            intermediate(0, colI) = -intermediate(3, colI) - intermediate(6, colI) - intermediate(9, colI);
            intermediate(1, colI) = -intermediate(4, colI) - intermediate(7, colI) - intermediate(10, colI);
            intermediate(2, colI) = -intermediate(5, colI) - intermediate(8, colI) - intermediate(11, colI);
        }

            for (int rowI = 0; rowI < 12; ++rowI) {
                double _000 = IB(0, 0) * intermediate(rowI, 0);
                double _013 = IB(0, 1) * intermediate(rowI, 3);
                double _026 = IB(0, 2) * intermediate(rowI, 6);
                double _001 = IB(0, 0) * intermediate(rowI, 1);
                double _014 = IB(0, 1) * intermediate(rowI, 4);
                double _027 = IB(0, 2) * intermediate(rowI, 7);
                double _002 = IB(0, 0) * intermediate(rowI, 2);
                double _015 = IB(0, 1) * intermediate(rowI, 5);
                double _028 = IB(0, 2) * intermediate(rowI, 8);
                double _100 = IB(1, 0) * intermediate(rowI, 0);
                double _113 = IB(1, 1) * intermediate(rowI, 3);
                double _126 = IB(1, 2) * intermediate(rowI, 6);
                double _101 = IB(1, 0) * intermediate(rowI, 1);
                double _114 = IB(1, 1) * intermediate(rowI, 4);
                double _127 = IB(1, 2) * intermediate(rowI, 7);
                double _102 = IB(1, 0) * intermediate(rowI, 2);
                double _115 = IB(1, 1) * intermediate(rowI, 5);
                double _128 = IB(1, 2) * intermediate(rowI, 8);
                double _200 = IB(2, 0) * intermediate(rowI, 0);
                double _213 = IB(2, 1) * intermediate(rowI, 3);
                double _226 = IB(2, 2) * intermediate(rowI, 6);
                double _201 = IB(2, 0) * intermediate(rowI, 1);
                double _214 = IB(2, 1) * intermediate(rowI, 4);
                double _227 = IB(2, 2) * intermediate(rowI, 7);
                double _202 = IB(2, 0) * intermediate(rowI, 2);
                double _215 = IB(2, 1) * intermediate(rowI, 5);
                double _228 = IB(2, 2) * intermediate(rowI, 8);
                d2e_dX2(rowI, 3) = _000 + _013 + _026; // 1
                d2e_dX2(rowI, 4) = _001 + _014 + _027; // 2
                d2e_dX2(rowI, 5) = _002 + _015 + _028; // 3
                d2e_dX2(rowI, 6) = _100 + _113 + _126; // 1
                d2e_dX2(rowI, 7) = _101 + _114 + _127; // 2
                d2e_dX2(rowI, 8) = _102 + _115 + _128; // 3
                d2e_dX2(rowI, 9) = _200 + _213 + _226; // 1
                d2e_dX2(rowI, 10) = _201 + _214 + _227; // 2
                d2e_dX2(rowI, 11) = _202 + _215 + _228; // 3
                d2e_dX2(rowI, 0) = -_200 - _213 - _226 - _100 - _113 - _126 - _000 - _013 - _026;
                d2e_dX2(rowI, 1) = -_001 - _014 - _027 - _101 - _114 - _127 - _201 - _214 - _227;
                d2e_dX2(rowI, 2) = -_002 - _015 - _028 - _102 - _115 - _128 - _202 - _215 - _228;
            }
        
    }

}