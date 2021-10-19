#pragma once
#include <iostream>
#include <cmath>
#include <random>

#undef BUFSIZ
#define BUFSIZ 256

#undef EPSILON
#define EPSILON -DBL_MIN

#undef NEAR_ZERO
#define NEAR_ZERO 1.0e-10

#undef NORM_NEAR_ZERO
#define NORM_NEAR_ZERO 1.0e-6

#undef NEAR_ZERO2
#define NEAR_ZERO2 1.0e-20

static const int HeightOffset = 0;
static const int WidthOffset = 0;

const int SCR_WIDTH = 1920;
const int SCR_HEIGHT = 1080;

#undef M_PI
#define M_PI 3.141592653589793

const double sphere_radius = 0.25;
#undef M_PI
#define M_PI 3.141592653589793

#undef myMax
#define myMax(a,b) (a>b?a:b)

#undef myMin
#define myMin(a,b) (a<b?a:b)

#undef SUB
#define SUB(dest,v1,v2) \
dest[0] = v1[0] - v2[0];\
dest[1] = v1[1] - v2[1];\
dest[2] = v1[2] - v2[2];


#undef MID
#define MID(dest,v1,v2)\
dest[0] =  0.5*(v1[0] + v2[0]);\
dest[1] =  0.5*(v1[1] + v2[1]);\
dest[2] =  0.5*(v1[2] + v2[2]);

#undef SUM
#define SUM(dest,v1,v2)\
dest[0] = v1[0] + v2[0];\
dest[1] = v1[1] + v2[1];\
dest[2] = v1[2] + v2[2];

#undef SUM_2
#define SUM_2(dest,v1,v2)\
dest[0] = v1[0] + v2[0];\
dest[1] = v1[1] + v2[1];

#undef SUB_2
#define SUB_2(dest,v1,v2) \
dest[0] = v1[0] - v2[0];\
dest[1] = v1[1] - v2[1];

#undef DEV
#define DEV(dest,v1,value)\
dest[0]=v1[0]/value;\
dest[1]=v1[1]/value;\
dest[2]=v1[2]/value;

#undef DOT
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#undef CROSS
#define CROSS(dest,v1,v2) \
dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#undef SIGN
#define SIGN(a) (a) > 1e-7 ? 1 : ((a) < -1e-7 ? -1 : 0);

#undef MULTI
#define MULTI(dest,v1,value)\
dest[0]=value*v1[0];\
dest[1]=value*v1[1];\
dest[2]=value*v1[2];

#undef BARYCENTRIC
#define BARYCENTRIC(dest,barycoe,v0,v1,v2)\
dest[0]=barycoe[0]*v0[0]+barycoe[1]*v1[0]+barycoe[2]*v2[0];\
dest[1]=barycoe[0]*v0[1]+barycoe[1]*v1[1]+barycoe[2]*v2[1];\
dest[2]=barycoe[0]*v0[2]+barycoe[1]*v1[2]+barycoe[2]*v2[2];

#undef POINT_ON_EDGE
#define POINT_ON_EDGE(dest,coe0,coe1,v0,v1)\
dest[0]=coe0*v0[0]+coe1*v1[0];\
dest[1]=coe0*v0[1]+coe1*v1[1];\
dest[2]=coe0*v0[2]+coe1*v1[2];

#undef POINT_ON_EDGE_MINUS
#define POINT_ON_EDGE_MINUS(dest,coe0,coe1,v0,v1)\
dest[0]=coe0*v0[0]-coe1*v1[0];\
dest[1]=coe0*v0[1]-coe1*v1[1];\
dest[2]=coe0*v0[2]-coe1*v1[2];

#undef TWO_POINTS_CENTER
#define TWO_POINTS_CENTER(dest,v0,v1)\
dest[0]=0.5*(v0[0]+v1[0]);\
dest[1]=0.5*(v0[1]+v1[1]);\
dest[2]=0.5*(v0[2]+v1[2]);



inline double gaussian(double x, double sigma) {
    return std::exp(-x * x / (2.0 * sigma * sigma));

}

inline void normalize(double x[3]) {
    double norm = sqrt(DOT(x, x));
    if (norm > 1e-10) {
        x[0] /= norm;
        x[1] /= norm;
        x[2] /= norm;
    }
}

#undef DEG_RADIANS
#define DEG_RADIANS(a) (a/180.0*M_PI)


