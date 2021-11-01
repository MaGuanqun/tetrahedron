#pragma once
#include<vector>
#include"../thread.h"
#include<array>
#include"../basic/global.h"
#include"../global_struct.h"
class MeshStruct
{
public:
	struct Face {
		int vertex[3];
		std::vector<int> edge;
		double area;
		std::vector<int>vertex_around;
		double initial_area;
	};
	struct Vertex {
		std::vector<int>face;
		std::vector<int>edge;
		std::vector<int>neighbor_vertex;
		std::vector<int>around_face;
		bool on_border = false;
	
	};
	struct Edge {
		int vertex[2];
		std::vector<int> face;
		std::vector<int> opposite_vertex;
		double length;
		bool is_border = false;
		bool isSame(int v0, int v1)
		{
			bool b0 = (vertex[0] == v0) && (vertex[1] == v1);
			bool b1 = (vertex[0] == v1) && (vertex[1] == v0);
			return b0 || b1;
		}
	};

	std::vector<std::array<double, 3>> vertex_position;
	std::vector<std::array<int,3>> triangle_indices;//if for tetrahedron, store the surface triangle	
	std::vector<std::array<double, 3>> vertex_for_render;
	std::vector<std::array<double, 3>> vertex_norm_for_render;
	std::vector<std::array<double, 3>> face_norm_for_render;
	std::vector<std::array<double, 3>> face_norm;
	std::vector<int>anchor_vertex;
	std::vector<std::array<double, 3>>anchor_position;

	Thread* thread;
	std::vector<double>mass;
	std::vector<double>triangle_normal_magnitude_reciprocal;

	std::vector<int> face_index_begin_per_thread;
	std::vector<int> anchor_index_begin_per_thread;
	std::vector<int> vertex_index_begin_per_thread;

	void initialNormalSize();


	void setAnchorPosition();
protected:
	int type;
private:

};

