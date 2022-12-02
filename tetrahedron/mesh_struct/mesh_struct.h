#pragma once
#include<vector>
#include"../thread.h"
#include<array>
#include"../basic/global.h"
#include"../global_struct.h"
#include"../collision/floating.h"
class MeshStruct
{
public:
	struct Face {
		//int vertex[3];
		//std::vector<int> edge;
		double area;
		std::vector<int>vertex_around;
		double initial_area;
	};
	struct Vertex {
		std::vector<unsigned int>face;
		std::vector<unsigned int>edge;
		std::vector<unsigned int>neighbor_vertex;
		std::vector<unsigned int>around_face;
		std::vector<unsigned int>around_vertex;
		bool on_border = false;
		std::vector<int>tetrahedron; //only for surface vertex

	};
	struct Edge {
		//int vertex[2];
		std::vector<int> face;
		std::vector<int> opposite_vertex;
		bool is_border = false;
		//bool isSame(int v0, int v1)
		//{
		//	if ((vertex[0] == v0) && (vertex[1] == v1))	{
		//		return true;
		//	}
		//	if ((vertex[0] == v1) && (vertex[1] == v0)) {
		//		return true;
		//	}
		//	return false;
		//}
	};

	std::vector<Vertex> vertices;
	std::vector<Face> faces;
	std::vector<Edge> edges;


	std::vector<std::vector<unsigned int>>face_around_face;//include itself
	std::vector<std::vector<unsigned int>>edge_around_face;
	std::vector<std::vector<unsigned int>>tet_around_face;

	std::vector<std::vector<unsigned int>>face_around_edge;
	std::vector<std::vector<unsigned int>>edge_around_edge;//include itself
	std::vector<std::vector<unsigned int>>tet_around_edge;


	std::vector<unsigned int>face_edges;// edge indices of every triangle
	std::vector<double>edge_length;// edge length of every triangle
	std::vector<unsigned int>edge_vertices;//vertex indices of every edge

	std::vector<std::vector<unsigned int>>unconnected_vertex_index; //for bending
	std::vector<std::vector<unsigned int>>unconnected_edge_index; //for edge_length
	//std::vector<std::vector<unsigned int>>unconnected_tet_index; 

	std::vector<std::vector<std::vector<unsigned int>>>tet_color_group; //for tet (ARAP)
	std::vector<std::vector<std::vector<char>>> tet_in_collision;//record if the tet is involved in collision
	std::vector<int>tet_order_in_color_group;// for every tet, we store group_0_color, index_in_this_color, group_1_color, index_in_this_color, ...


	std::vector<std::array<double, 3>> vertex_position;
	std::vector<std::array<int, 3>> triangle_indices;//if for tetrahedron, store the surface triangle	
	std::vector<std::array<double, 3>> vertex_for_render;
	std::vector<std::array<double, 3>> vertex_normal_for_render;
	std::vector<std::array<double, 3>> face_normal_for_render;
	std::vector<std::array<double, 3>> face_normal;
	std::vector<std::array<double, 3>> vertex_normal;
	std::vector<int>anchor_vertex;
	std::vector<std::array<double, 3>>anchor_position;

	std::vector<std::array<double, 3>> ori_face_normal_for_render;//not normalized
	std::vector<std::array<double, 3>> ori_face_normal;//not normalized
	std::vector<std::array<double, 3>> cross_for_approx_CCD;// (p1(t1)-p0(t1))?p2(t0)-p0(t0))+(p1(t0)-p0(t0))x(p2(t1)-p0(t1)) t0 is the start of the time step, t1 is the current

	std::vector<unsigned int> unfixed_point_index;

	std::vector<unsigned int>real_index_to_unfixed_index;

	std::vector<unsigned int> unfixed_edge_vertex_index; // the index is for the unfixed index, need to transfer it to get the real index
	std::vector<unsigned int> only_one_vertex_fix_edge; //first index is unfixed (use unfixed index, need to transfer it to get the real index), second index is fixed, it is directly the real index

	
	Thread* thread;
	std::vector<double>mass;
	std::vector<double>mass_inv;
	//std::vector<double> initial_mass_inv;
	std::vector<double>triangle_normal_magnitude_reciprocal;

	std::vector<unsigned int> face_index_begin_per_thread;
	std::vector<unsigned int> anchor_index_begin_per_thread;

	std::vector<unsigned int> edge_index_begin_per_thread;

	std::vector<unsigned int> tet_edge_index_begin_per_thread;

	std::vector<unsigned int> only_one_vertex_fixed_edge_index_begin_per_thread;
	std::vector<unsigned int> unfixed_vertex_index_begin_per_thread;
	std::vector<unsigned int> unfixed_edge_index_begin_per_thread;

	std::vector<double> unfixed_rest_edge_length;
	std::vector<double> fixed_one_vertex_rest_edge_length;


	std::vector<std::array<floating, 3>> f_face_normal_for_render;
	std::vector<std::array<floating, 3>> f_face_normal;
	std::vector<std::array<floating, 3>> f_cross_for_approx_CCD;
	void initialNormalSize();
	void setAnchorPosition();
	void setVertex();
	void setFace();
	void setEdge();
	void setEdgeForSpring();
	void addArounVertex();

	std::vector<unsigned int> vertex_index_begin_per_thread;

	void getRenderFaceNormalPerThread(int thread_id);
	void getFaceNormalPerThread(int thread_id);


	void getFaceNormal();
	void getRenderFaceNormal();

	std::vector<std::array<int, 3>> surface_triangle_index_in_order;//this is for representative triangle

	void updateAnchorPerThread(int total_thread_num);
	void resetMassInv();

	std::vector<std::vector<unsigned int>> vertex_tet_index;//tet indices that contain the vertex

	std::vector<std::vector<unsigned int>> tet_tet_index;// tet indices that around a tet, does not include itself

	std::vector<unsigned int>tet_unfixed_vertex_num;
	//store common vertex order of tet's neighbor tets. For every tet, 
	//store  2,2,1,1,2,1,0,1 means for the first neighbor tet, there are 2 vertices in common, 
	// 2,1 means the two common vertices in the neighbor tet's vertex index's order is 2,1
	// 1,2 means the two common vertices in this tet's vertex index's order is 1,2 (Here exclude fixed vertex)
	//1 means for the second neighbor tet, there is 1 vertex in common
	// 0 means the two common vertices in the neighbor tet's vertex index's order is 0
	// 1 means the two common vertices in this tet's vertex index's order is 1 (Here exclude fixed vertex)
	std::vector<std::vector<unsigned int>> tet_neighbor_tet_vertex_order; 

	std::vector<bool>is_vertex_fixed;


	void sortTriangleAroundVertexEdge(int thread_id);
	void setFaceEdgeAroundFace(int thread_id);

	std::vector<int> vertex_surface_index;//size is the global vertex size, verted index -> surface index


	void testIfRIght(std::vector<std::vector<unsigned int>>& k);

	void testFaceEdgeAroundFaceEdge();


protected:
	int type;
	bool isEdgeExist(unsigned int v0, unsigned int v1, unsigned int& edge_index);
	void addEdge(int v0, int v1, int face, int opposite_vertex, int edge_index_indicator);
	void addNeighborVertex();

	bool isCommonUsed(unsigned int tet_num, std::vector<unsigned int>* tet_index);

private:

};

