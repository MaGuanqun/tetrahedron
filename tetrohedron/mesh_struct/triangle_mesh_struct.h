#pragma once
#include"mesh_struct.h"

class TriangleMeshStruct:public MeshStruct
{
public:
	TriangleMeshStruct();

	std::vector<Vertex> vertices;
	std::vector<Face> faces;
	std::vector<Edge> edges;

	std::vector<std::array<double, 3>> face_normal;//if for tetrohedron, store the surface triangle
	std::vector<std::array<double, 3>> vertex_normal;//if for tetrohedron, store the surface triangle

	std::vector<int> vertex_index_begin_per_thread;
	std::vector<int> edge_index_begin_per_thread;



	void setVertex() override;
	void setFace();
	void setEdge();
	void getRenderVertexNormalPerThread(int thread_id) override;
	void setThreadIndex(int total_thread_num_) override;
	void getNormal();

	void getFaceNormalPerThread(int thread_id) override;
	void getVertexNormalPerThread(int thread_id) override;
	void initialInfo();
private:
	bool isEdgeExist(int v0, int v1, int& edge_index);
	void addEdge(int v0, int v1, int face, int opposite_vertex);

};
