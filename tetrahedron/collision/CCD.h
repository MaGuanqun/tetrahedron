#pragma once
#include <algorithm>
#include"../basic/global.h"
#include"primitive_distance.h"
#include<vector>
namespace CCD {

    template <class T>
    T pointTriangleCcdTest(const T* _p,
        const T* _t0,
        const T* _t1,
        const T* _t2,
        const T* p_new,
        const T* t0_new,
        const T* t1_new,
        const T* t2_new,
        T eta, T thickness)
    {
        T p[3], t0[3], t1[3], t2[3], dp[3], dt0[3], dt1[3], dt2[3];
        int size = 3 * sizeof(T);
        memcpy(p, _p, size);
        memcpy(t0, _t0, size);
        memcpy(t1, _t1, size);
        memcpy(t2, _t2, size);
        SUB(dp, p_new, p);
        SUB(dt0, t0_new, t0);
        SUB(dt1, t1_new, t1);
        SUB(dt2, t2_new, t2);
        T move[3];
        for (int i = 0; i < 3; ++i) {
            move[i] = 0.25 * (dp[i] + dt0[i] + dt1[i] + dt2[i]);
        }
        SUB_(dt0, move);
        SUB_(dt1, move);
        SUB_(dt2, move);
        SUB_(dp, move);
        T disp_mag2_vec[3] = { DOT(dt0,dt0),DOT(dt1,dt1), DOT(dt2,dt2) };
        T max_disp_mag = sqrt(DOT(dp, dp)) + std::sqrt((std::max)((std::max)(disp_mag2_vec[0], disp_mag2_vec[1]), disp_mag2_vec[2]));
        if (max_disp_mag == 0)
            return 1.0;




        T dist2_cur = internal::pointTriangleDistanceUnclassified(p, t0, t1, t2);
        T dist_cur = std::sqrt(dist2_cur);


    


        T gap = eta * (dist2_cur - thickness * thickness) / (dist_cur + thickness);
        T toc = 0.0;
        int itr = 0;

        T ori_dist2_cur = dist2_cur;

        double dis_recrod;




        while (true) {
            T toc_lower_bound;
            if (itr == 0) {
                toc_lower_bound = (1.0 - eta) * (dist2_cur - thickness * thickness) / ((dist_cur + thickness) * max_disp_mag);
            }
            std::cout << toc_lower_bound << " " << (1.0 - eta) << " " << (dist2_cur - thickness * thickness) / ((dist_cur + thickness) * max_disp_mag) << " " 
                << dist2_cur << " " << dist_cur << " " << max_disp_mag << std::endl;

            T p6[3], t60[3], t61[3], t62[3];
            int size = 3 * sizeof(T);
            memcpy(p6, p, size);
            memcpy(t60, t0, size);
            memcpy(t61, t1, size);
            memcpy(t62, t2, size);

            ACCUMULATE_SUM_WITH_COE(p, toc_lower_bound, dp);
            ACCUMULATE_SUM_WITH_COE(t0, toc_lower_bound, dt0);
            ACCUMULATE_SUM_WITH_COE(t1, toc_lower_bound, dt1);
            ACCUMULATE_SUM_WITH_COE(t2, toc_lower_bound, dt2);

            std::cout << "dis " << dp[1] << " " << toc_lower_bound << " " << p[1] - t0[1] << std::endl;

            std::cout << -0.000735 * 1.01959e-13 << " " << -0.000735 * 1.01959e-13 / 0.9 << std::endl;
            std::cout << p[1] -0.000735 * 1.01959e-13 -t0[1] << " " <<p[1] - 0.000735 * 1.01959e-13 / 0.9 -t0[1] << std::endl;

            dis_recrod = dist2_cur;

            dist2_cur = internal::pointTriangleDistanceUnclassified(p, t0, t1, t2);
            dist_cur = std::sqrt(dist2_cur);

            T temp0[3], temp1[3], temp2[3];
            SUB(temp0, t1, t0);
            SUB(temp1, t2, t0);
            CROSS(temp2, temp0, temp1);
            SUB(temp0, p, t0);
            T aTb = DOT(temp0, temp2);
            std::cout <<"temp "<< temp0[0] << " " << temp0[1] << " " << temp0[2] << " " << temp2[0] << " " << temp2[1] << " " << temp2[2] << std::endl;
            std::cout <<DOT(temp0, temp0)<<" "<<aTb << " " << aTb * aTb << " " << DOT(temp2, temp2) << std::endl;

            std::cout << "toc_lower_bound " << max_disp_mag << " " << toc_lower_bound << " " << dist2_cur << std::endl;
            std::cout << "type2 " << internal::pointTriangleDistanceType(p, t0, t1, t2) << " " << dist_cur << " " << dist2_cur << " " << max_disp_mag << " " << toc_lower_bound << std::endl;

            if (dist2_cur == 0.0) {

                std::cout << internal::pointTriangleDistanceUnclassified(p, t0, t1, t2) << " " << internal::pointTriangleDistanceType(p, t0, t1, t2) << std::endl;
                std::cout << "point triangle distance " << internal::pointTriangleDistance(p, t0, t1, t2) << std::endl;;
                std::cout << "PT CCD distance equals zero, warning " << itr << std::endl;
                //std::cout<<
                //std::cout << toc_lower_bound << std::endl;

                break;
            }

            if ((toc && ((dist2_cur - thickness * thickness) / (dist_cur + thickness) < gap))) {// || (
                // std::cout << toc<<" "<< toc_lower_bound<<" "<< dist_cur - thickness<<" "<<gap << std::endl;
                break;
            }
            toc += toc_lower_bound;

            if (toc < 0.0) {
                std::cout << "warning: the vt distance is less than the thickness " << dis_recrod << " " << dist2_cur - thickness * thickness << " " << thickness << std::endl;
                return 0.0;
            }

            if (toc > 1.0) {
                return 1.0;
            }
            itr++;
        }
        return toc;
    }



    template <class T>
    T pointTriangleCcd(
        const T* _p,
        const T* _t0,
        const T* _t1,
        const T* _t2,
        const T* p_new,
        const T* t0_new,
        const T* t1_new,
        const T* t2_new,
        T eta, T thickness)
    {
        T p[3], t0[3], t1[3], t2[3], dp[3], dt0[3], dt1[3], dt2[3];
        int size = 3 * sizeof(T);
        memcpy(p, _p, size);
        memcpy(t0, _t0, size);
        memcpy(t1, _t1, size);
        memcpy(t2, _t2, size);
        SUB(dp, p_new, p);
        SUB(dt0, t0_new, t0);
        SUB(dt1, t1_new, t1);
        SUB(dt2, t2_new, t2);
        T move[3];
        for (int i = 0; i < 3; ++i) {
            move[i] = 0.25 * (dp[i] + dt0[i] + dt1[i] + dt2[i]);
        }

        //double k = sqrt(DOT(dp, dp));
        //std::cout << "real div " <<k << std::endl;

        SUB_(dt0, move);
        SUB_(dt1, move);
        SUB_(dt2, move);
        SUB_(dp, move);
        T disp_mag2_vec[3] = { DOT(dt0,dt0),DOT(dt1,dt1), DOT(dt2,dt2) };
        T max_disp_mag = sqrt(DOT(dp, dp)) + std::sqrt((std::max)((std::max)(disp_mag2_vec[0], disp_mag2_vec[1]), disp_mag2_vec[2]));
        if (max_disp_mag == 0)
            return 1.0;
        T dist2_cur = internal::pointTriangleDistanceUnclassified(p, t0, t1, t2);
        T dist_cur = std::sqrt(dist2_cur);
        T gap = eta * (dist2_cur - thickness * thickness) / (dist_cur + thickness);

        T toc = 0.0;
        int itr = 0;

        T ori_dist2_cur = dist2_cur;

        double dis_recrod;
        //std::cout << "=== "<<k-max_disp_mag << std::endl;

        while (true) {
            T toc_lower_bound = (1.0 - eta) * (dist2_cur - thickness * thickness) / ((dist_cur + thickness) * max_disp_mag);
            ACCUMULATE_SUM_WITH_COE(p, toc_lower_bound, dp);
            ACCUMULATE_SUM_WITH_COE(t0, toc_lower_bound, dt0);
            ACCUMULATE_SUM_WITH_COE(t1, toc_lower_bound, dt1);
            ACCUMULATE_SUM_WITH_COE(t2, toc_lower_bound, dt2);

            dis_recrod = dist2_cur;

            dist2_cur = internal::pointTriangleDistanceUnclassified(p, t0, t1, t2);
            dist_cur = std::sqrt(dist2_cur);


            //std::cout << "toc_lower_bound " <<max_disp_mag<<" "<<  toc_lower_bound << " " << dist2_cur << std::endl;
            //std::cout << _p[1] - t0_new[1] << " " << p_new[1] << " " << p[1] - t0[1] << std::endl;

            if (dist2_cur == 0.0) {
               
                //std::cout << internal::pointTriangleDistanceUnclassified(p, t0, t1, t2)<<" "<< internal::pointTriangleDistanceType(p, t0, t1, t2) << std::endl;
                std::cout << "point triangle distance " <<std::endl;// internal::pointTriangleDistance(p, t0, t1, t2) << 
                std::cout << "PT CCD distance equals zero, warning " << itr << std::endl;
                //pointTriangleCcdTest(_p, _t0, _t1, _t2, p_new, t0_new, t1_new, t2_new, eta, thickness);
                //std::cout<<
                //std::cout << toc_lower_bound << std::endl;

                break;
            }

            if ((toc && ((dist2_cur - thickness * thickness) / (dist_cur + thickness) < gap)) ) {// || itr > 500
                //if (itr > 500) {
                    //std::cout <<"itr "<< itr << std::endl;
                //}
                // std::cout << toc<<" "<< toc_lower_bound<<" "<< dist_cur - thickness<<" "<<gap << std::endl;
                break;
            }
            toc += toc_lower_bound;

            if (toc < 0.0) {
                std::cout << "warning: the vt distance is less than the thickness " << dis_recrod << " " << dist2_cur - thickness * thickness << " " << thickness << std::endl;
                return 0.0;
            }

            if (toc > 1.0) {
                if (itr > 500) {
                    std::cout << "itr " << itr<<" "<< max_disp_mag<<" "<< gap << std::endl;
                }
                return 1.0;
            }
            itr++;
        }
        return toc;
    }

   

    //template <class T>
    //T pointTriangleCcdPrint(
    //    const T* _p,
    //    const T* _t0,
    //    const T* _t1,
    //    const T* _t2,
    //    const T* p_new,
    //    const T* t0_new,
    //    const T* t1_new,
    //    const T* t2_new,
    //    T eta, T thickness)
    //{
    //    T p[3], t0[3], t1[3], t2[3], dp[3], dt0[3], dt1[3], dt2[3];
    //    int size = 3 * sizeof(T);
    //    memcpy(p, _p, size);
    //    memcpy(t0, _t0, size);
    //    memcpy(t1, _t1, size);
    //    memcpy(t2, _t2, size);
    //    SUB(dp, p_new, p);
    //    SUB(dt0, t0_new, t0);
    //    SUB(dt1, t1_new, t1);
    //    SUB(dt2, t2_new, t2);
    //    T move[3];
    //    for (int i = 0; i < 3; ++i) {
    //        move[i] = 0.25 * (dp[i] + dt0[i] + dt1[i] + dt2[i]);
    //    }
    //    SUB_(dt0, move);
    //    SUB_(dt1, move);
    //    SUB_(dt2, move);
    //    SUB_(dp, move);
    //    T disp_mag2_vec[3] = { DOT(dt0,dt0),DOT(dt1,dt1), DOT(dt2,dt2) };
    //    T max_disp_mag = sqrt(DOT(dp, dp)) + std::sqrt((std::max)((std::max)(disp_mag2_vec[0], disp_mag2_vec[1]), disp_mag2_vec[2]));
    //    if (max_disp_mag == 0)
    //        return 1.0;
    //    T dist2_cur = internal::pointTriangleDistanceUnclassified(p, t0, t1, t2);

    //    //std::cout << internal::pointTriangleDistanceType(p, t0, t1, t2) << std::endl;

    //    T dist_cur = std::sqrt(dist2_cur);
    //    //T gap = eta * (dist_cur - thickness);
    //    T gap = eta * (dist2_cur - thickness * thickness) / (dist_cur + thickness);

    //    T toc = 0.0;
    //    int itr = 0;

    //    T ori_dist2_cur = dist2_cur;

    //    double dis_recrod;

    //    //std::cout << "max_disp_mag " << max_disp_mag << std::endl;
    //    //std::cout << "gap " << gap<<" "<< dist2_cur<<" "<<(p[1]-t0[1])*(p[1]-t0[1]) << std::endl;

    //    while (true) {
    //        // T toc_lower_bound = (1 - eta) * (dist_cur - thickness) / max_disp_mag;
    //        T toc_lower_bound = (1.0 - eta) * (dist_cur - thickness) / max_disp_mag;

    //     
    //        ACCUMULATE_SUM_WITH_COE(p, toc_lower_bound, dp);
    //        ACCUMULATE_SUM_WITH_COE(t0, toc_lower_bound, dt0);
    //        ACCUMULATE_SUM_WITH_COE(t1, toc_lower_bound, dt1);
    //        ACCUMULATE_SUM_WITH_COE(t2, toc_lower_bound, dt2);

    //        dis_recrod = dist2_cur;

    //        dist2_cur = internal::pointTriangleDistanceUnclassified(p, t0, t1, t2);
    //        dist_cur = std::sqrt(dist2_cur);
    //        std::cout << internal::pointTriangleDistanceType(p, t0, t1, t2) << std::endl;
    //        std::cout << toc<<" "<< toc_lower_bound << " " << dist_cur << " " << dist_cur - thickness - gap <<" "<<p[1]-t0[1]
    //            << " " << thickness << std::endl;


    //        if ((toc>0.0 && (dist_cur -thickness< gap)) || itr > 300) {// || 
    //            // std::cout << toc<<" "<< toc_lower_bound<<" "<< dist_cur - thickness<<" "<<gap << std::endl;

    //            break;
    //        }
    //        toc += toc_lower_bound;

    //        if (toc < 0.0) {
    //            std::cout << "warning: the vt distance is less than the thickness " << dis_recrod << " " << dist2_cur - thickness * thickness << " " << thickness << std::endl;
    //            return 0.0;
    //        }

    //        if (toc > 1.0) {
    //            return 1.0;
    //        }
    //        itr++;
    //    }
    //    return toc;
    //}


    template <class T>
    T edgeEdgeCcd(
        const T* _ea0, const T* _ea1, const T* _eb0, const T* _eb1,
        const T* ea0_new, const T* ea1_new, const T* eb0_new, const T* eb1_new,
        T eta, T thickness)
    {
        T ea0[3], ea1[3], eb0[3], eb1[3], dea0[3], dea1[3], deb0[3], deb1[3];
        unsigned int size = 3 * sizeof(T);
        memcpy(ea0, _ea0, size);
        memcpy(ea1, _ea1, size);
        memcpy(eb0, _eb0, size);
        memcpy(eb1, _eb1, size);
        SUB(dea0, ea0_new, ea0);
        SUB(dea1, ea1_new, ea1);
        SUB(deb0, eb0_new, eb0);
        SUB(deb1, eb1_new, eb1);
        T move[3];
        for (int i = 0; i < 3; ++i) {
            move[i] = 0.25 * (dea0[i] + dea1[i] + deb0[i] + deb1[i]);
        }
        SUB_(dea0, move);
        SUB_(dea1, move);
        SUB_(deb0, move);
        SUB_(deb1, move);
        T max_disp_mag = std::sqrt((std::max)(DOT(dea0, dea0), DOT(dea1, dea1))) + std::sqrt((std::max)(DOT(deb0, deb0), DOT(deb1, deb1)));
        if (max_disp_mag == 0)
            return 1.0;
        T dist2_cur = internal::edgeEdgeDistanceUnclassified(ea0, ea1, eb0, eb1);
        T dFunc = dist2_cur - thickness * thickness;
        // since we ensured other place that all dist smaller than dHat are positive,
       // this must be some far away nearly parallel edges
        std::vector<T> dists(4);
        if (dFunc <= 0) {
            dists[0] = SQUARED_LENGTH(ea0, eb0);
            dists[1] = SQUARED_LENGTH(ea0, eb1);
            dists[2] = SQUARED_LENGTH(ea1, eb0);
            dists[3] = SQUARED_LENGTH(ea1, eb1);
            dist2_cur = *std::min_element(dists.begin(), dists.end());
            dFunc = dist2_cur - thickness * thickness;
        }
        T dist_cur = std::sqrt(dist2_cur);
        T gap = eta * dFunc / (dist_cur + thickness);
        T toc = 0.0;
        int itr = 0;
        while (true)
        {
            T toc_lower_bound = (1.0 - eta) * dFunc / ((dist_cur + thickness) * max_disp_mag); 
            ACCUMULATE_SUM_WITH_COE(ea0, toc_lower_bound, dea0);
            ACCUMULATE_SUM_WITH_COE(ea1, toc_lower_bound, dea1);
            ACCUMULATE_SUM_WITH_COE(eb0, toc_lower_bound, deb0);
            ACCUMULATE_SUM_WITH_COE(eb1, toc_lower_bound, deb1);

            dist2_cur = internal::edgeEdgeDistanceUnclassified(ea0, ea1, eb0, eb1);

            if (dist2_cur == 0.0) {
                std::cout << "EE CCD distance equals zero, warning " << std::endl;
                break;
            }
            dFunc = dist2_cur - thickness * thickness;
            if (dFunc <= 0) {
                // since we ensured other place that all dist smaller than dHat are positive,
                // this must be some far away nearly parallel edges
                dists[0] = SQUARED_LENGTH(ea0, eb0);
                dists[1] = SQUARED_LENGTH(ea0, eb1);
                dists[2] = SQUARED_LENGTH(ea1, eb0);
                dists[3] = SQUARED_LENGTH(ea1, eb1);
                dist2_cur = *std::min_element(dists.begin(), dists.end());
                dFunc = dist2_cur - thickness * thickness;
            }
            dist_cur = std::sqrt(dist2_cur);

            if ((toc && (dFunc / (dist_cur + thickness) < gap)) ) {//(|| itr > 500
                break;
            }
            toc += toc_lower_bound;
            if (toc > 1.0)
                return 1.0;
            if (toc < 0.0) {
                std::cout << "warning: the ee distance is less than the thickness " << std::endl;
                return 0.0;
            }
            itr++;
        }
        return toc;
    }
    inline void test()
    {
        //double initial_pos[3] = { 0.657142857142857, 0.7, 0.785714285714286 };
    //double current_pos[3] = { 0.657142857142857, 0.699755000142151, 0.785714285714286 };
    //double initial_triangle_0[3] = { 0.657142857142857, 0.633333333333333, 0.776190476190476 };
    //double initial_triangle_1[3] = { 0.692857142857143, 0.666666666666667, 0.780952380952381 };
    //double initial_triangle_2[3] = { 0.692857142857143, 0.633333333333333, 0.776190476190476 };
    //double current_triangle_0[3] = { 0.657142857142857, 0.633088333368321, 0.776190476190476 };
    //double current_triangle_1[3] = { 0.692857142857143, 0.666421666742275, 0.780952380952381 };
    //double current_triangle_2[3] = { 0.692857142857143, 0.633088333369544, 0.776190476190476 };
/*    double initial_pos[3] = { 1.09,0.00001,0.0 };
    double current_pos[3] = { -1.09,0.00001,0.0 };
    double initial_triangle_0[3] = { 1.0,0.0,0.0 };
    double initial_triangle_1[3] = { -1.0,0.0,-1.0 };
    double initial_triangle_2[3] = { -1.0,0.0,1.0 };
    double current_triangle_0[3] = { 1.0,0.0,0.0 };
    double current_triangle_1[3] = { -1.0,0.0,-1.0 };
    double current_triangle_2[3] = { -1.0,0.0,1.0 };*/

    ////edge 1
    //    double initial_pos[3] = { 0.0,1.0,1.0 };
    //    double initial_triangle_0[3] = { 0.0,1.0,3.0 };
    //    double current_pos[3] = { 0.0,0.0,2.0 };
    //    double current_triangle_0[3] = { 0.0,0.0,4.0 };

    //    //edge 2
    //    double initial_triangle_1[3] = { 0.0,0.0,1.0 };
    //    double initial_triangle_2[3] = { 0.0,0.0,3.0 };

    //    double current_triangle_1[3] = { 0.0,0.0,1.0 };
    //    double current_triangle_2[3] = { 0.0,0.0,3.0 };


        //double* initial_pos_a; double* initial_triangle_0_a; double* initial_triangle_1_a; double* initial_triangle_2_a;
        //double* current_pos_a; double* current_triangle_0_a; double* current_triangle_1_a; double* current_triangle_2_a;
        //initial_pos_a = initial_pos;
        //initial_triangle_0_a = initial_triangle_0;
        //initial_triangle_1_a = initial_triangle_1;
        //initial_triangle_2_a = initial_triangle_2;
        //current_pos_a = current_pos;
        //current_triangle_0_a = current_triangle_0;
        //current_triangle_1_a = current_triangle_1;
        //current_triangle_2_a = current_triangle_2;
        //double edge[3];
        //SUB(edge, initial_triangle_0, initial_triangle_1);
        //double tolerance = sqrt(DOT(edge, edge)) * 1e-4;
        //double time = edgeEdgeCcd(initial_pos_a, initial_triangle_0_a, initial_triangle_1_a, initial_triangle_2_a,
        //    current_pos_a, current_triangle_0_a, current_triangle_1_a, current_triangle_2_a, 0.1, tolerance);
        //std::cout << "ccd " << time << std::endl;


        double initial_pos[3] = { 0.6, - 0.530341, - 0.5 };
        double initial_triangle_0[3] = { 0.1, - 0.530341, 0.2 };

        double current_pos[3] = { -0.600375, - 1.9947, - 0.198086 };
        double current_triangle_0[3] = { 0.401501, - 0.519231, - 0.207114 };

        double bary[4];
        double d_2 =internal::edgeEdgeNearestPoint(initial_pos, initial_triangle_0, current_pos,
            current_triangle_0, bary);
        std::cout << internal::edgeEdgeDistanceType(initial_pos, initial_triangle_0, current_pos,
            current_triangle_0) << std::endl;

        double point_0[3];
        double point_1[3];
        POINT_ON_EDGE(point_0, bary[0], bary[1], initial_pos, initial_triangle_0);
        POINT_ON_EDGE(point_1, bary[2], bary[3], current_pos, current_triangle_0);

        std::cout << point_0[0] << " " << point_0[1] << " " << point_0[2] << std::endl;
        std::cout << point_1[0] << " " << point_1[1] << " " << point_1[2] << std::endl;

        std::cout << internal::edgeEdgeDistance(initial_pos, initial_triangle_0, current_pos, current_triangle_0, bary) << std::endl;;

        std::cout << "initial d " << d_2<<" "<<bary[0]<<" "<<bary[1]<<" "<<bary[2]<<" "<<bary[3] << std::endl;

    }
} //


