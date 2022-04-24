#pragma once
#include"mesh_struct.h"
#include"../external/Eigen/Dense"
#include<map>
using namespace Eigen;

class TetrahedronMeshStruct :public MeshStruct
{
public:
	void findSurface();
	void setThreadIndex(int total_thread_num_);
	void getRenderNormal();
	//void getRenderFaceNormalPerThread(int thread_id);
	std::vector<std::array<int, 4>> indices;
	std::vector<double>volume;
	std::vector<int>tet_index_of_surface_face;

	std::vector<unsigned int> tetrahedron_index_begin_per_thread;
	void setVolume(int thread_No);
	double setMass(double density);
	std::vector<bool>vertex_on_surface;
	std::vector<unsigned int>vertex_index_on_sureface; //size is the surface vertex size, surface index -> vertex index
	std::vector<int> vertex_surface_index;//size is the global vertex size, verted index -> surface index
	std::vector<unsigned int> vertex_index_on_surface_begin_per_thread;

	std::vector<Matrix3d> P_inv;
	std::vector<Matrix<double,4,3>> AT;

	void prepareForDeformationGradient();
	//std::vector<Matrix<double, 3, 3>>PT;//to record the restshape (PP^T)^-1 for deformation gradient
	//std::vector<Matrix<double, 4, 3>>PT_position;//to record the restshape (PP^T)^-1 for deformation gradien
	//std::vector<Matrix<double, 4, 3>>PT_PPT_inv;//to record the restshape (PP^T)^-1 for deformation gradien
	//std::vector<double>PPT_determinant;//to record the restshape (PP^T)^-1 for deformation gradient
	//std::vector<Matrix<double, 4, 3>> PT;//to record the restshape P^T for deformation gradient

	//void getFaceNormalPerThread(int thread_id);
	void getNormal();
	void recordTetIndexForSurfaceIndex();

	void getVertexNormalPerThread(int thread_id);
	void getRenderVertexNormalPerThread(int thread_id);

	//std::vector<std::array<int, 2>> edge_vertex_index_on_surface;


	//void setVertexIndexOnSurfaceEdgeTriangle();

private:
	struct TetrahedronFace {
		std::array<int, 3> index;
		int sorted_index[3];
		TetrahedronFace(int v0, int v1, int v2) {
			index.data()[0] = v0;
			index.data()[1] = v1;
			index.data()[2] = v2;
			memcpy(sorted_index, index.data(), 12);
			std::sort(sorted_index, sorted_index + 3);
		}
		bool operator<(const TetrahedronFace& t1) const
		{
			if (sorted_index[0] < t1.sorted_index[0])
				return true;
			else if (sorted_index[0] == t1.sorted_index[0]) {
				if (sorted_index[1] < t1.sorted_index[1])
					return true;
				else if (sorted_index[1] == t1.sorted_index[1]) {
					if (sorted_index[2] < t1.sorted_index[2])
						return true;
				}
			}
			return false;
		}

	};
	void buildMap(std::map<TetrahedronFace, int>& face_in_tet, int v0, int v1, int v2, std::vector<int>& face_tet_index, int tet_index);
	double getTetrahedronVolume(double* v1, double* v2, double* v3, double* v4);
	Matrix<double, 3, 3> constructMatrixP(int tetra_index);
	Matrix<double, 3, 4> constructMatrixP_pos(int tetra_index);
	Matrix<double, 3, 4> constructDeformGradientA(Matrix3d& p);
};

