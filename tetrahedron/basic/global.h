#ifndef GLOBAL_H
#define GLOBAL_H


#include <iostream>
#include <cmath>
#include <random>

#undef BUFSIZ
#define BUFSIZ 256

#undef BOUND_UNIT
#define BOUND_UNIT 0.00000000000000022204460492503136

#undef EPSILON
#define EPSILON -DBL_MIN

#undef NEAR_ZERO
#define NEAR_ZERO 1.0e-10

#undef NORM_NEAR_ZERO
#define NORM_NEAR_ZERO 1.0e-7

#undef NEAR_ZERO2
#define NEAR_ZERO2 1.0e-20

static const int HeightOffset = 0;
static const int WidthOffset = 0;

const int SCR_WIDTH = 1920;
const int SCR_HEIGHT = 1080;

#undef M_PI
#define M_PI 3.14159265358979323846

const double sphere_radius = 0.25;

#undef myMax
#define myMax(a,b) (a>b?a:b)

#undef myMin
#define myMin(a,b) (a<b?a:b)

#undef SUB
#define SUB(dest,v1,v2) \
dest[0] = v1[0] - v2[0];\
dest[1] = v1[1] - v2[1];\
dest[2] = v1[2] - v2[2];



#undef SUB_
#define SUB_(dest,v1) \
dest[0] -= v1[0];\
dest[1] -= v1[1];\
dest[2] -= v1[2];

#undef SUB_2_
#define SUB_2_(dest,v1) \
dest[0] -= v1[0];\
dest[1] -= v1[1];

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


#undef SUM_
#define SUM_(dest,v1)\
dest[0] += v1[0];\
dest[1] += v1[1];\
dest[2] += v1[2];

#undef SUM_2_
#define SUM_2_(dest,v1)\
dest[0] += v1[0];\
dest[1] += v1[1];

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

#undef DEV_
#define DEV_(dest,value)\
dest[0]/=value;\
dest[1]/=value;\
dest[2]/=value;

#undef DOT
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#undef DOT_2
#define DOT_2(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1])

#undef EDGE_LENGTH
#define EDGE_LENGTH(v1,v2) ((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2]))



#undef CROSS
#define CROSS(dest,v1,v2) \
dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
dest[2] = v1[0] * v2[1] - v1[1] * v2[0];


#undef CROSS_FOUR_POINTS
#define CROSS_FOUR_POINTS(dest, p00,p01,p10,p11) \
dest[0] = (p01[1]-p00[1])*(p11[2]-p10[2])- (p01[2]-p00[2])*(p11[1]-p10[1]);\
dest[1] = (p01[2]-p00[2])*(p11[0]-p10[0])- (p01[0]-p00[0])*(p11[2]-p10[2]);\
dest[2] = (p01[0]-p00[0])*(p11[1]-p10[1])- (p01[1]-p00[1])*(p11[0]-p10[0]);




#undef SIGN
#define SIGN(a) (a) > 1e-7 ? 1 : ((a) < -1e-7 ? -1 : 0);

#undef MULTI
#define MULTI(dest,v1,value)\
dest[0]=value*v1[0];\
dest[1]=value*v1[1];\
dest[2]=value*v1[2];


#undef MULTI_
#define MULTI_(dest,value)\
dest[0]*=value;\
dest[1]*=value;\
dest[2]*=value;

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

#undef AABB_CENTER
#define AABB_CENTER(dest,v0)\
dest[0]=0.5*(v0[0]+v0[3]);\
dest[1]=0.5*(v0[1]+v0[4]);\
dest[2]=0.5*(v0[2]+v0[5]);

#undef DIFF_SIGN
#define DIFF_SIGN(v0,v1) ((v0>0&&v1<0)||(v0<0||v1>0))

#undef SAME_SIGN
#define SAME_SIGN(v0,v1) ((v0>0&&v1>0)||(v0<0||v1<0))

#undef ACCUMULATE_SUM_WITH_COE
#define ACCUMULATE_SUM_WITH_COE(dest,coe,v0)\
dest[0] += coe* v0[0];\
dest[1] += coe* v0[1];\
dest[2] += coe* v0[2];


#undef ACCUMULATE_SUM_WITH_COE_DEV
#define ACCUMULATE_SUM_WITH_COE_DEV(dest,coe,v0)\
dest[0] += v0[0]/coe;\
dest[1] += v0[1]/coe;\
dest[2] += v0[2]/coe;


#undef SQUARED_LENGTH
#define SQUARED_LENGTH(v0,v1) ((v0[0]-v1[0])*(v0[0]-v1[0])+ (v0[1]-v1[1])*(v0[1]-v1[1]) + (v0[2]-v1[2])*(v0[2]-v1[2]))


inline double gaussian(double x, double sigma) {
	return std::exp(-x * x / (2.0 * sigma * sigma));

}

inline void normalize(double* x) {
	double norm = sqrt(DOT(x, x));
	if (norm > 1e-15) {
		x[0] /= norm;
		x[1] /= norm;
		x[2] /= norm;
	}
}


#undef DEG_RADIANS
#define DEG_RADIANS(a) (a/180.0*M_PI)


#undef DET2X2
#define DET2X2(a,b,c,d) (a*d-b*c)

inline void countInEveryThread(unsigned int total_thread_num, unsigned int total_num, unsigned int* begin)
{
	unsigned int interval1 = total_num / total_thread_num;
	unsigned int resi1 = total_num % total_thread_num;

	for (unsigned int i = 0; i < resi1; ++i) {
		begin[i] = interval1 + 1;
	}
	for (unsigned int i = resi1; i < total_thread_num; ++i) {
		begin[i] = interval1;
	}
}

inline void arrangeIndex(unsigned int total_thread_num, unsigned int total_num, unsigned int* begin) {
	unsigned int interval1 = total_num / total_thread_num;
	unsigned int resi1 = total_num % total_thread_num;
	if (resi1 == 0) {
		for (unsigned int i = 0; i < total_thread_num; ++i) {
			begin[i] = interval1 * i;
		}
	}
	else {
		for (unsigned int i = 0; i < total_thread_num; ++i) {
			if (i < resi1) {
				begin[i] = (interval1 + 1) * i;
			}
			else {
				begin[i] = (interval1 + 1) * resi1 + interval1 * (i - resi1);
			}
		}
	}
	begin[total_thread_num] = total_num;
}


inline void arrangeIndex(int total_thread_num, int total_num, std::vector<int>& begin) {
	int interval1 = total_num / total_thread_num;
	int resi1 = total_num % total_thread_num;
	if (resi1 == 0) {
		for (int i = 0; i < total_thread_num; ++i) {
			begin[i] = interval1 * i;
		}
	}
	else {
		for (int i = 0; i < total_thread_num; ++i) {
			if (i < resi1) {
				begin[i] = (interval1 + 1) * i;
			}
			else {
				begin[i] = (interval1 + 1) * resi1 + interval1 * (i - resi1);
			}
		}
	}
	begin[total_thread_num] = total_num;
}

//rotation matrix store in column major
inline void rotateAroundVector(double* rotation_matrix, double* rotate_axis, double angle)
{
	double cos_ = cos(angle);
	double sin_ = sin(angle);
	rotation_matrix[0] = rotate_axis[0] * rotate_axis[0] * (1.0 - cos_) + cos_;
	rotation_matrix[1] = rotate_axis[0] * rotate_axis[1] * (1.0 - cos_) + rotate_axis[2] * sin_;
	rotation_matrix[2] = rotate_axis[0] * rotate_axis[2] * (1.0 - cos_) - rotate_axis[1] * sin_;
	rotation_matrix[3] = rotate_axis[0] * rotate_axis[1] * (1.0 - cos_) - rotate_axis[2] * sin_;
	rotation_matrix[4] = rotate_axis[1] * rotate_axis[1] * (1.0 - cos_) + cos_;
	rotation_matrix[5] = rotate_axis[2] * rotate_axis[1] * (1.0 - cos_) + rotate_axis[0] * sin_;
	rotation_matrix[6] = rotate_axis[2] * rotate_axis[0] * (1.0 - cos_) + rotate_axis[1] * sin_;
	rotation_matrix[7] = rotate_axis[2] * rotate_axis[1] * (1.0 - cos_) - rotate_axis[0] * sin_;
	rotation_matrix[8] = rotate_axis[2] * rotate_axis[2] * (1.0 - cos_) + cos_;
	normalize(rotation_matrix);
	normalize(rotation_matrix + 3);
	normalize(rotation_matrix + 6);
}

inline void inverse3X3(double* matrix, double* inverse)
{
	double det = matrix[0] * (matrix[8] * matrix[4] - matrix[7] * matrix[5])
		- matrix[1] * (matrix[8] * matrix[3] - matrix[5] * matrix[6])
		+ matrix[2] * (matrix[7] * matrix[3] - matrix[4] * matrix[6]);

	//if (det < 1e-15) {
	//	std::cout << det << std::endl;
	//}

	//std::cout << "my determi " << det << std::endl;
	inverse[0] = (matrix[8] * matrix[4] - matrix[7] * matrix[5]) / det;
	inverse[1] = (matrix[2] * matrix[7] - matrix[8] * matrix[1]) / det;
	inverse[2] = (matrix[5] * matrix[1] - matrix[2] * matrix[4]) / det;
	inverse[3] = (matrix[5] * matrix[6] - matrix[8] * matrix[3]) / det;
	inverse[4] = (matrix[8] * matrix[0] - matrix[2] * matrix[6]) / det;
	inverse[5] = (matrix[2] * matrix[3] - matrix[5] * matrix[0]) / det;
	inverse[6] = (matrix[7] * matrix[3] - matrix[4] * matrix[6]) / det;
	inverse[7] = (matrix[1] * matrix[6] - matrix[7] * matrix[0]) / det;
	inverse[8] = (matrix[4] * matrix[0] - matrix[1] * matrix[3]) / det;
}

inline bool edgeEdgeconnected(unsigned int* edge1, unsigned int* edge2)
{
	if (edge1[0] == edge2[0])
		return true;
	if (edge1[1] == edge2[1])
		return true;
	if (edge1[0] == edge2[1])
		return true;
	if (edge1[1] == edge2[0])
		return true;
	return false;
}

inline bool vertexInTriangle(int* face_index, int vertex_index)
{
	if (face_index[0] == vertex_index)
		return true;
	if (face_index[1] == vertex_index)
		return true;
	if (face_index[2] == vertex_index)
		return true;
	return false;
}

#endif // !GLOBAL_H