#pragma once
#include"primitive_distance.h"
#include"primitive_distance_gradient_hessian.h"

namespace distance_grad_hessian {
	template <class T>
	int getVTDistanceGradHessian(T* p, T* t0, T* t1, T* t2, T* g, T* h)
	{
		switch (CCD::internal::pointTriangleDistanceType(p,t0,t1,t2))
		{
        case 0:
      
            return 0;
        case 1:
         
            return 1;
        case 2:
            distance::point_point_distance_gradient(p, t2, g);
            distance::point_point_distance_hessian(p, t2, h);
            return 2;
        case 3:
            distance::point_edge_distance_gradient(p, t0, t1, g);
            distance::point_edge_distance_hessian(p, t0, t1, h);
            return 3;
        case 4:
            distance::point_edge_distance_gradient(p, t1, t2, g);
            distance::point_edge_distance_hessian(p, t1, t2, h);
            return 4;
        case 5:
            distance::point_edge_distance_gradient(p, t0, t2, g);
            distance::point_edge_distance_hessian(p, t0, t2, h);
            return 5;
        case 6:
            distance::point_triangle_distance_gradient(p, t0, t1, t2, g);
            distance::point_triangle_distance_hessian(p, t0, t1, t2, h);
            return 6;
		} 
	}

    template <class T>
    void getEEDistanceGradHessian(T* ea0, T* ea1, T* eb0, T* eb1, T* g, T* h)
    {
        switch (CCD::internal::edgeEdgeDistanceType(ea0, ea1, eb0, eb1)) {
        case 0:
            distance::point_point_distance_gradient(ea0, eb0, g);
            distance::point_point_distance_hessian(ea0, eb0, h);
            break;
        case 1:
            distance::point_point_distance_gradient(ea0, eb1, g);
            distance::point_point_distance_hessian(ea0, eb1, h);
            break;
        case 2:
            distance::point_edge_distance_gradient(ea0, eb0, eb1, g);
            distance::point_edge_distance_hessian(ea0, eb0, eb1, h);
            break;
        case 3:
            distance::point_point_distance_gradient(ea1, eb0, g);
            distance::point_point_distance_hessian(ea1, eb0, h);
            break;
        case 4:
            distance::point_point_distance_gradient(ea1, eb1, g);
            distance::point_point_distance_hessian(ea1, eb1, h);
            break;
        case 5:
            distance::point_edge_distance_gradient(ea1, eb0, eb1, g);
            distance::point_edge_distance_hessian(ea1, eb0, eb1, h);
            break;
        case 6:
            distance::point_edge_distance_gradient(eb0, ea0, ea1, g);
            distance::point_edge_distance_hessian(eb0, ea0, ea1, h);
            break;
        case 7:
            distance::point_edge_distance_gradient(eb1, ea0, ea1, g);
            distance::point_edge_distance_hessian(eb1, ea0, ea1, h);
            break;
        case 8:
            distance::edge_edge_distance_gradient(ea0, ea1, eb0, eb1,g);
            distance::edge_edge_distance_hessian(ea0, ea1, eb0, eb1,h);
            break;

        }
    }

}
