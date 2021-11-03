#pragma once
#include"mesh_struct.h"

class TriangleMeshStruct:public MeshStruct
{
public:
	TriangleMeshStruct();

	std::vector<Vertex> vertices;
	std::vector<Face> faces;
	std::vector<Edge> edges;

	std::vector<std::array<double, 3>> vertex_normal;

	std::vector<int> edge_index_begin_per_thread;

	void setVertex();
	void setFace();
	void setEdge();
	void getRenderVertexNormalPerThread(int thread_id);
	void setThreadIndex(int total_thread_num_);
	void getNormal();
	void getRenderNormal();
	void getRenderFaceNormalPerThread(int thread_id);
	void getFaceNormalPerThread(int thread_id);
	void getVertexNormalPerThread(int thread_id);
	void initialInfo();

	//void getVertexNormal();
private:
	bool isEdgeExist(int v0, int v1, int& edge_index);
	void addEdge(int v0, int v1, int face, int opposite_vertex);

};
