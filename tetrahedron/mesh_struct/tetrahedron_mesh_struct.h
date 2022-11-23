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

	std::vector<std::vector<unsigned int>>triangle_index_of_a_tet;//record all surface triangle indices of a tet
	std::vector<std::vector<unsigned int>>edge_index_of_a_tet;//record all surface edge indices of a tet

	//std::vector<std::vector<unsigned int>> edge_index_in_a_tet;//record all edge whose vertices all contained in a tet

	std::vector<std::array<int, 4>> unfixied_indices; //put the unfixed index at first in order between [0,3], the index in that tet
	std::vector<std::array<int, 4>> unfixied_actual_indices; //unfixed index in a tet 



	std::vector<std::vector<unsigned int>> triangle_index_of_a_tet_color;// record surface triangle index which contains by a tet color group
	std::vector<std::vector<unsigned int>> edge_index_of_a_tet_color;// record edge index which contains by a tet color group
	std::vector<std::vector<unsigned int>> surface_vertex_index_of_a_tet_color;// record vertex index which contains by a tet color group

	std::vector<std::vector<unsigned int>> triangle_index_of_a_tet_color_per_thread_start;
	std::vector<std::vector<unsigned int>> edge_index_of_a_tet_color_per_thread_start;
	std::vector<std::vector<unsigned int>> surface_vertex_index_of_a_tet_color_per_thread_start;


	std::vector<std::vector<unsigned int>>tet_around_tet_color_group;//excluede tet in the group

	std::vector<unsigned int> tetrahedron_index_begin_per_thread;
	void setVolume(int thread_No);
	double setMass(double density);
	std::vector<bool>vertex_on_surface;
	std::vector<unsigned int>vertex_index_on_sureface; //size is the surface vertex size, surface index -> vertex index

	std::vector<unsigned int> vertex_index_on_surface_begin_per_thread;


	void recordTetIndexForVertex();
	//std::vector<Matrix3d> P_inv;
	std::vector<Matrix<double,3,4>> A;

	void prepareForDeformationGradient();
	//std::vector<Matrix<double, 3, 3>>PT;//to record the restshape (PP^T)^-1 for deformation gradient
	//std::vector<Matrix<double, 4, 3>>PT_position;//to record the restshape (PP^T)^-1 for deformation gradien
	//std::vector<Matrix<double, 4, 3>>PT_PPT_inv;//to record the restshape (PP^T)^-1 for deformation gradien
	//std::vector<double>PPT_determinant;//to record the restshape (PP^T)^-1 for deformation gradient
	//std::vector<Matrix<double, 4, 3>> PT;//to record the restshape P^T for deformation gradient

	//void getFaceNormalPerThread(int thread_id);
	void getNormal();
	void recordTetIndexForSurfaceIndex();
	void updateTetNeighborInfo();
	void recordTetIndexForTet();

	void getVertexNormalPerThread(int thread_id);
	void getRenderVertexNormalPerThread(int thread_id);
	void getVertexNormalFromRenderPerThread(int thread_id);

	//std::vector<std::array<int, 2>> edge_vertex_index_on_surface;

	std::vector<unsigned int> tet_edge_vertices;

	std::vector<double>tet_edge_rest_length;

	void initialUnfixedIndex();


	void updateUnfixedPointData();

	//void setVertexIndexOnSurfaceEdgeTriangle();
	void setTetEdges();
	void updateTetNeighborTetVertexIndex(int thread_id);

	void recordPrimitiveIndexOfATet();

	void sortTriangleAroundElement();
	void sortTetAroundVertexEdge(int thread_id);
	void setTetAroundFace(int thread_id);

	void testTetAroundFaceEdge();

	void obtainVETofColors();

	void setTetAroundTetColor(int thread_No);

	std::vector<unsigned int> tet_color_index_start_per_thread;
	std::vector<std::vector<unsigned int>> tet_around_tet_color_group_start_per_thread; //with tet_around_tet_color_group
	

	std::vector<std::vector<unsigned int>>tet_in_a_group_start_per_thread;

	void setTetColorStartPerThread();
private:


	void obtainVETofAColor(int color);
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

	Matrix<double, 3, 3> constructMatrixP(int tetra_index);
	Matrix<double, 3, 4> constructMatrixP_pos(int tetra_index);
	Matrix<double, 3, 4> constructDeformGradientA(Matrix3d& p);


	void addTetEdges(unsigned int p0, unsigned int p1, std::vector<std::vector<unsigned int>>& edge_vertex);

	void updateTetNeighborTetVertexIndex();
	void findCommonVertexInOrder(int* tet_0_index, int* tet_1_index, std::vector<unsigned int>* vertex_in_order, double* inv_mass);
	unsigned int unfixedVertexIndexInATet(int* tet, int index, double* inv_mass);
	void updateUnfixedTetVertexIndexInfo();
	void recordTriangleIndexOfATet();
	void recordEdgeIndexOfATet();

	//void recordEdgeIndexInATet();
	bool checkEdgeInATet(unsigned int* edge_vertices, int* tet_indices);

	void setTetAroundTetColor(std::vector<bool>& is_used, std::vector<unsigned int>* unconnected_tet_index,
		std::vector<unsigned int>* tet_around_tet_color_group);

};

