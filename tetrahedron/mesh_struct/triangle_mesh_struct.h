#pragma once
#include"mesh_struct.h"


class TriangleMeshStruct:public MeshStruct
{
public:
	TriangleMeshStruct();
	
	void setThreadIndex(int total_thread_num_);
	void getNormal();
	void getRenderNormal();
	
	
	void initialInfo();
	
	std::vector<int> vertex_index_begin_per_thread;

	void getVertexNormalPerThread(int thread_id);
	void getRenderVertexNormalPerThread(int thread_id);
	//void getVertexNormal();
private:


};
