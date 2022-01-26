#pragma once
#include"mesh_struct.h"


class TriangleMeshStruct:public MeshStruct
{
public:
	TriangleMeshStruct();
	
	void setThreadIndex(int total_thread_num_);
	void getNormal();
	void getRenderNormal();
	void getRenderFaceNormalPerThread(int thread_id);
	void getFaceNormalPerThread(int thread_id);
	
	void initialInfo();
	
	//void getVertexNormal();
private:


};
