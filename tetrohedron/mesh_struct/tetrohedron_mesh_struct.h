#pragma once
#include"mesh_struct.h"
class TetrohedronMeshStruct:public MeshStruct
{
public:

	void getRenderVertexNormalPerThread(int thread_id);
	void setThreadIndex(int total_thread_num_) override {};

	void getFaceNormalPerThread(int thread_id);
	void getVertexNormalPerThread(int thread_id);

	void setVertex() override {};
private:

};

