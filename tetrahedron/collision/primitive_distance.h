#pragma once
#include"../basic/global.h"
namespace CCD {
    namespace internal {
        //For a matrix M(2X3)=(basic0,basis1)^T, here to solve MM^Tx=Mvec    
        template <class T>
        void solveSyStem(const T* basis0, const T* basis1, const T* vec, T* result)
        {
            T MMT[3];
            MMT[0] = DOT(basis0, basis0);//MMT(0,0)
            MMT[1] = DOT(basis1, basis1);//MMT(1,1)
            MMT[2] = DOT(basis0, basis1);//MMT(0,1), MMT(1,0)

            T mVec[2];
            mVec[0] = DOT(basis0, vec);
            mVec[1] = DOT(basis1, vec);

            T dev = MMT[0] * MMT[1] - MMT[2] * MMT[2];
            result[0] = (MMT[1] * mVec[0] - MMT[2] * mVec[1]) / dev;
            result[1] = (MMT[0] * mVec[1] - MMT[2] * mVec[0]) / dev;
        }



        template <class T>
        int pointTriangleDistanceType(
            const T* p,
            const T* t0,
            const T* t1,
            const T* t2)
        {
            T basis0[3], basis1[3];
            SUB(basis0, t1, t0);
            SUB(basis1, t2, t0);

            T n[3];
            CROSS(n, basis0, basis1);
            CROSS(basis1, basis0, n);
            T pt[3];
            SUB(pt, p, t0);
            T para0[2];
            solveSyStem(basis0, basis1, pt, para0);
            if (para0[0] > 0.0 && para0[0] < 1.0 && para0[1] >= 0.0) {
                return 3; // PE t0t1
            }
            else {
                SUB(basis0, t2, t1);
                CROSS(basis1, basis0, n);
                T para1[2];
                SUB(pt, p, t1);
                solveSyStem(basis0, basis1, pt, para1);

                if (para1[0] > 0.0 && para1[0] < 1.0 && para1[1] >= 0.0) {
                    return 4; // PE t1t2
                }
                else {
                    SUB(basis0, t0, t2);
                    CROSS(basis1, basis0, n);
                    T para2[2];
                    SUB(pt, p, t2);
                    solveSyStem(basis0, basis1, pt, para2);
                    if (para2[0] > 0.0 && para2[0] < 1.0 && para2[1] >= 0.0) {
                        return 5; // PE t2t0
                    }
                    else {
                        if (para0[0] <= 0.0 && para2[0] >= 1.0) {
                            return 0; // PP t0
                        }
                        else if (para1[0] <= 0.0 && para0[0] >= 1.0) {
                            return 1; // PP t1
                        }
                        else if (para2[0] <= 0.0 && para1[0] >= 1.0) {
                            return 2; // PP t2
                        }
                        else {
                            return 6; // PT
                        }
                    }
                }
            }
        }

        template <class T>
        int edgeEdgeDistanceType(
            const T* ea0,
            const T* ea1,
            const T* eb0,
            const T* eb1)
        {
            T u[3], v[3], w[3];
            SUB(u, ea1, ea0);
            SUB(v, eb1, eb0);
            SUB(w, ea0, eb0);
            T a = DOT(u, u); // always >= 0
            T b = DOT(u, v);
            T c = DOT(v, v); // always >= 0
            T d = DOT(u, w);
            T e = DOT(v, w);
            T D = a * c - b * b; // always >= 0
            T tD = D; // tc = tN / tD, default tD = D >= 0
            T sN, tN;

            int defaultCase = 8;
            // compute the line parameters of the two closest points
            sN = (b * e - c * d);
            if (sN <= 0.0) { // sc < 0 => the s=0 edge is visible
                tN = e;
                tD = c;
                defaultCase = 2;
            }
            else if (sN >= D) { // sc > 1  => the s=1 edge is visible
                tN = e + b;
                tD = c;
                defaultCase = 5;
            }
            else {
                tN = (a * e - b * d);
                T n[3];
                CROSS(n, u, v);
                if (tN > 0.0 && tN < tD && (DOT(n, w) == 0.0 || DOT(n, n) < NEAR_ZERO2 * a * c)) {
                    // if (tN > 0.0 && tN < tD && (u.cross(v).dot(w) == 0.0 || u.cross(v).squaredNorm() == 0.0)) {
                    // std::cout << u.cross(v).squaredNorm() / (a * c) << ": " << sN << " " << D << ", " << tN << " " << tD << std::endl;
                    // avoid coplanar or nearly parallel EE
                    if (sN < 0.5 * D) {
                        tN = e;
                        tD = c;
                        defaultCase = 2;
                    }
                    else {
                        tN = e + b;
                        tD = c;
                        defaultCase = 5;
                    }
                }
                // else defaultCase stays as 8
            }

            if (tN <= 0.0) { // tc < 0 => the t=0 edge is visible
    // recompute sc for this edge
                if (-d <= 0.0) {
                    return 0;
                }
                else if (-d >= a) {
                    return 3;
                }
                else {
                    return 6;
                }
            }
            else if (tN >= tD) {
                if ((-d + b) <= 0.0) {
                    return 1;
                }
                else if ((-d + b) >= a) {
                    return 4;
                }
                else {
                    return 7;
                }
            }
            return defaultCase;
        }

        template <class T>
        int edgeEdgeDistanceType(
            const T* ea0,
            const T* ea1,
            const T* eb0,
            const T* eb1,
            T* barycentric)
        {
            T u[3], v[3], w[3];
            SUB(u, ea1, ea0);
            SUB(v, eb1, eb0);
            SUB(w, ea0, eb0);
            T a = DOT(u, u); // always >= 0
            T b = DOT(u, v);
            T c = DOT(v, v); // always >= 0
            T d = DOT(u, w);
            T e = DOT(v, w);
            T D = a * c - b * b; // always >= 0
            T tD = D; // tc = tN / tD, default tD = D >= 0
            T sN, tN;

            int defaultCase = 8;
            // compute the line parameters of the two closest points
            sN = (b * e - c * d);
            if (sN <= 0.0) { // sc < 0 => the s=0 edge is visible
                tN = e;
                tD = c;
                barycentric[0] = 1.0;
                barycentric[1] = 0.0;
                barycentric[2] = 1.0-tN/tD;
                barycentric[3] = tN / tD;
                defaultCase = 2;
            }
            else if (sN >= D) { // sc > 1  => the s=1 edge is visible
                tN = e + b;
                tD = c;

                barycentric[0] = 0.0;
                barycentric[1] = 1.0;
                barycentric[2] = 1.0 - tN / tD;
                barycentric[3] = tN / tD;
                defaultCase = 5;
            }
            else {
                tN = (a * e - b * d);
                T n[3];
                CROSS(n, u, v);
                if (tN > 0.0 && tN < tD && (DOT(n, w) == 0.0 || DOT(n, n) < NEAR_ZERO2 * a * c)) {
                    // if (tN > 0.0 && tN < tD && (u.cross(v).dot(w) == 0.0 || u.cross(v).squaredNorm() == 0.0)) {
                    // std::cout << u.cross(v).squaredNorm() / (a * c) << ": " << sN << " " << D << ", " << tN << " " << tD << std::endl;
                    // avoid coplanar or nearly parallel EE
                    if (sN < 0.5* D) {
                        tN = e;
                        tD = c;
                        defaultCase = 2;
                        barycentric[0] = 1.0;
                        barycentric[1] = 0.0;
                        barycentric[2] = 1.0 - tN / tD;
                        barycentric[3] = tN / tD;
                    }
                    else {
                        tN = e + b;
                        tD = c;
                        defaultCase = 5;
                        barycentric[0] = 0.0;
                        barycentric[1] = 1.0;
                        barycentric[2] = 1.0 - tN / tD;
                        barycentric[3] = tN / tD;

                    }
                }
                // else defaultCase stays as 8
                else {
                    barycentric[0] = 1.0 - sN / D;
                    barycentric[1] = sN / D;
                    barycentric[2] = 1.0 - tN / D;
                    barycentric[3] = tN / D;
                }

            }

            if (tN <= 0.0) { // tc < 0 => the t=0 edge is visible
    // recompute sc for this edge
                if (-d <= 0.0) {
                    barycentric[0] = 1.0;
                    barycentric[1] = 0.0;
                    barycentric[2] = 1.0;
                    barycentric[3] = 0.0;
                    return 0;
                }
                else if (-d >= a) {
                    barycentric[0] = 0.0;
                    barycentric[1] = 1.0;
                    barycentric[2] = 1.0;
                    barycentric[3] = 0.0;
                    return 3;
                }
                else {
                    barycentric[0] = 1.0+d/a;
                    barycentric[1] =-d/a;
                    barycentric[2] = 1.0;
                    barycentric[3] = 0.0;
                    return 6;
                }
            }
            else if (tN >= tD) {
                if ((-d + b) <= 0.0) {
                    barycentric[0] = 1.0;
                    barycentric[1] = 0.0;
                    barycentric[2] = 0.0;
                    barycentric[3] = 1.0;
                    return 1;
                }
                else if ((-d + b) >= a) {
                    barycentric[0] = 0.0;
                    barycentric[1] = 1.0;
                    barycentric[2] = 0.0;
                    barycentric[3] = 1.0;
                    return 4;
                }
                else {
                    barycentric[0] = 1.0-(b-d)/a;
                    barycentric[1] = (b - d) / a;
                    barycentric[2] = 0.0;
                    barycentric[3] = 1.0;
                    return 7;
                }
            }
            return defaultCase;
        }

        template <class T> //squared distance
        inline T pointPointDistance(const T* a, const T* b)
        {
            T vec[3];
            SUB(vec, a, b);
            return DOT(vec, vec);
        }

        template <class T> //squared distance
        inline T pointEdgeDistance(const T* p, const T* e0, const T* e1)
        {
            T temp0[3], temp1[3], temp2[3];
            SUB(temp0, e0, p);
            SUB(temp1, e1, p);
            CROSS(temp2, temp0, temp1);
            SUB(temp0, e1, e0);
            return DOT(temp2, temp2) / DOT(temp0, temp0);
        }

        template <class T> //squared distance
        inline T pointEdgeDistance(const T* p, const T* e0, const T* e1, T* barycentric)
        {
            T temp0[3], temp1[3], temp2[3];
            SUB(temp0, e0, p);
            SUB(temp1, e1, p);
            CROSS(temp2, temp0, temp1);
            SUB(temp0, e1, e0);
            return DOT(temp2, temp2) / DOT(temp0, temp0);
        }

        template <class T>
        int pointEdgeDistanceType(const T* p,
            const T* e0,
            const T* e1)
        {
            T e[3], temp[3];
            SUB(e, e1, e0);
            SUB(temp, p, e0);
            T ratio = DOT(e, temp) / DOT(e, e);
            if (ratio < 0) {
                return 0; // PP (p-e0)
            }
            else if (ratio > 1) {
                return 1; // PP (p-e1)
            }
            else {
                return 2; // PE
            }
        }
        template <class T> //squared distance
        inline T pointTriangleDistance(const T* p, const T* t0, const T* t1, const T* t2)
        {
            T temp0[3], temp1[3], temp2[3];
            SUB(temp0, t1, t0);
            SUB(temp1, t2, t0);
            CROSS(temp2, temp0, temp1);
            SUB(temp0, p, t0);
            T aTb = DOT(temp0, temp2);
            return aTb * aTb / DOT(temp2, temp2);
        }
        template <class T>
        inline T edgeEdgeDistance(const T* ea0, const T* ea1, const T* eb0, const T* eb1)
        {
            T temp0[3], temp1[3], temp2[3];
            SUB(temp0, ea1, ea0);
            SUB(temp1, eb1, eb0);
            CROSS(temp2, temp0, temp1);
            SUB(temp0, eb0, ea0);
            T aTb = DOT(temp0, temp2);
            return aTb * aTb / DOT(temp2, temp2);
        }

        template <class T>
        inline T edgeEdgeDistance(const T* ea0, const T* ea1, const T* eb0, const T* eb1, const T* barycentric)
        {
            T temp0[3];
            temp0[0] = barycentric[0] * ea0[0] + barycentric[1] * ea1[0] + barycentric[2] * eb0[0] + barycentric[3] * eb1[0];
            temp0[1] = barycentric[0] * ea0[1] + barycentric[1] * ea1[1] + barycentric[2] * eb0[1] + barycentric[3] * eb1[1];
            temp0[2] = barycentric[0] * ea0[2] + barycentric[1] * ea1[2] + barycentric[2] * eb0[2] + barycentric[3] * eb1[2];
            return DOT(temp0, temp0);
        }

        template <class T>
        T pointEdgeDistanceUnclassified(
            const T* p,
            const T* e0,
            const T* e1)
        {
            T e[3], temp[3];
            SUB(e, e1, e0);
            SUB(temp, p, e0);
            T ratio = DOT(e, temp);
            if (ratio <= 0.0) {
                return pointPointDistance(p, e0); // PP (p-e0)
            }
            else
            {
                T c2 = DOT(e, e);
                if (ratio >= c2) {
                    return pointPointDistance(p, e1); // PP (p-e1)
                }
                else {
                    ratio /= c2;
                    e[0] = e0[0] + ratio * e[0];
                    e[1] = e0[1] + ratio * e[1];
                    e[2] = e0[2] + ratio * e[2];
                    return pointPointDistance(p, e); // PE
                }
            }
        }
        template <class T>
        T pointEdgeDistanceUnclassified(
            const T* p,
            const T* e0,
            const T* e1,
            T* barycentric)
        {
            T e[3], temp[3];
            SUB(e, e1, e0);
            SUB(temp, p, e0);
            T ratio = DOT(e, temp);
            if (ratio <= 0.0) {
                barycentric[0] = 1.0;
                barycentric[1] = 0.0;
                return pointPointDistance(p, e0); // PP (p-e0)
            }
            else
            {
                T c2 = DOT(e, e);
                if (ratio >= c2) {
                    barycentric[0] = 0.0;
                    barycentric[1] = 1.0;
                    return pointPointDistance(p, e1); // PP (p-e1)
                }
                else {
                    ratio /= c2;
                    e[0] = e0[0] + ratio * e[0];
                    e[1] = e0[1] + ratio * e[1];
                    e[2] = e0[2] + ratio * e[2];
                    barycentric[0] = 1.0 - ratio;
                    barycentric[1] = ratio;
                    return pointPointDistance(p, e); // PE
                }
            }
        }

        template <class T>
        T pointTriangleDistanceUnclassified(
            const T* p,
            const T* t0,
            const T* t1,
            const T* t2) //return squared distance
        {
            switch (pointTriangleDistanceType(p, t0, t1, t2)) {
            case 0:
                return pointPointDistance(p, t0);
            case 1:
                return pointPointDistance(p, t1);
            case 2:
                return pointPointDistance(p, t2);
            case 3:
                return pointEdgeDistance(p, t0, t1);
            case 4:
                return pointEdgeDistance(p, t1, t2);
            case 5:
                return pointEdgeDistance(p, t2, t0);
            case 6:
                return pointTriangleDistance(p, t0, t1, t2);
            default:
                return (std::numeric_limits<T>::max)();
            }
        }


        template <class T>
        T edgeEdgeDistanceUnclassified(
            const T* ea0,
            const T* ea1,
            const T* eb0,
            const T* eb1)
        {
            switch (edgeEdgeDistanceType(ea0, ea1, eb0, eb1)) {
            case 0:
                return pointPointDistance(ea0, eb0);
            case 1:
                return pointPointDistance(ea0, eb1);
            case 2:
                return pointEdgeDistance(ea0, eb0, eb1);
            case 3:
                return pointPointDistance(ea1, eb0);
            case 4:
                return pointPointDistance(ea1, eb1);
            case 5:
                return pointEdgeDistance(ea1, eb0, eb1);
            case 6:
                return pointEdgeDistance(eb0, ea0, ea1);
            case 7:
                return pointEdgeDistance(eb1, ea0, ea1);
            case 8:
                return edgeEdgeDistance(ea0, ea1, eb0, eb1);
            default:
                return (std::numeric_limits<T>::max)();
            }
        }

        template <class T>
        T edgeEdgeNearestPoint(
            const T* ea0,
            const T* ea1,
            const T* eb0,
            const T* eb1,
            T* barycentrc)
        {
            switch (edgeEdgeDistanceType(ea0, ea1, eb0, eb1,barycentrc)) {
            case 0:
                return pointPointDistance(ea0, eb0);
            case 1:
                return pointPointDistance(ea0, eb1);
            case 2:
                return pointEdgeDistance(ea0, eb0, eb1);
            case 3:
                return pointPointDistance(ea1, eb0);
            case 4:
                return pointPointDistance(ea1, eb1);
            case 5:
                return pointEdgeDistance(ea1, eb0, eb1);
            case 6:
                return pointEdgeDistance(eb0, ea0, ea1);
            case 7:
                return pointEdgeDistance(eb1, ea0, ea1);
            case 8:
                return edgeEdgeDistance(ea0, ea1, eb0, eb1,barycentrc);
            default:
                return (std::numeric_limits<T>::max)();
            }
        }

        template <class T>
        T pointTriangleNearestDistance(
            const T* p,
            const T* t0,
            const T* t1,
            const T* t2)
        {
            T d_2;
            T barycentric[3];
            T S[3];
            SUB(S, p, t0);
            T E1[3], E2[3], S1[3], S2[3];
            SUB(E1, t1, t0);
            SUB(E2, t2, t0);
            T triangle_normal[3];
            CROSS(triangle_normal, E1, E2);
            normalize(triangle_normal);

            CROSS(S1, triangle_normal, E2);
            CROSS(S2, S, E1);
            T temp = 1.0 / DOT(S1, E1);
            d_2 = temp * DOT(S2, E2);
            d_2 *= d_2;
            barycentric[1] = temp * DOT(S1, S);
            barycentric[2] = temp * DOT(S2, triangle_normal);
            barycentric[0] = 1.0 - barycentric[1] - barycentric[2];

            if (barycentric[0] >= 0.0 && barycentric[1] >= 0.0 && barycentric[2] >= 0.0) {
                return d_2;
            }
            else {
                d_2 = pointEdgeDistanceUnclassified(p, t0, t1);
                T d_;
                d_ = pointEdgeDistanceUnclassified(p, t0, t2);
                if (d_ < d_2) {
                    d_2 = d_;
                }
                d_ = pointEdgeDistanceUnclassified(p, t1, t2);
                if (d_ < d_2) {
                    d_2 = d_;
                }
            }
            return d_2;
        }


        template <class T>
        T pointTriangleNearestDistance(
            const T* p,
            const T* t0,
            const T* t1,
            const T* t2,
            const T* triangle_normal,
            T* barycentric, T d_hat_2)
        {
            T d_2;
           // T barycentric[3];
            T S[3];
            SUB(S, p, t0);
            T E1[3], E2[3], S1[3], S2[3];
            SUB(E1, t1, t0);
            SUB(E2, t2, t0);
            CROSS(S1, triangle_normal, E2);
            CROSS(S2, S, E1);
            T temp = 1.0 / DOT(S1, E1);
            d_2 = temp * DOT(S2, E2);
            d_2 *= d_2;

            if (d_2 >= d_hat_2) {
                return d_2;
            }

            barycentric[1] = temp * DOT(S1, S);
            barycentric[2] = temp * DOT(S2, triangle_normal);
            barycentric[0] = 1.0 - barycentric[1] - barycentric[2];

            if (barycentric[0] > 0.0 && barycentric[1] > 0.0 && barycentric[2] > 0.0) {
                return d_2;
            }
            else {
                T bary_centric2[2];

                d_2 = pointEdgeDistanceUnclassified(p, t0, t1, bary_centric2);
                barycentric[0] = bary_centric2[0];
                barycentric[1] = bary_centric2[1];
                barycentric[2] = 0.0;
                T d_;
                d_ = pointEdgeDistanceUnclassified(p, t0, t2, bary_centric2);
                if (d_ < d_2) {
                    d_2 = d_;
                    barycentric[0] = bary_centric2[0];
                    barycentric[2] = bary_centric2[1];
                    barycentric[1] = 0.0;
                }
                d_ = pointEdgeDistanceUnclassified(p, t1, t2, bary_centric2);
                if (d_ < d_2) {
                    d_2 = d_;
                    barycentric[1] = bary_centric2[0];
                    barycentric[2] = bary_centric2[1];
                    barycentric[0] = 0.0;
                }
            }
            return d_2;
        }
    }
}
