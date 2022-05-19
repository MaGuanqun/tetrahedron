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
	
	
	void getVertexNormalPerThread(int thread_id);
	void getRenderVertexNormalPerThread(int thread_id);
	//void getVertexNormal();
	double setMass(double density);

	void updateUnfixedPointData();
	void initialUnfixedIndex();

private:


};
