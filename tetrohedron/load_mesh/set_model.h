#pragma once

#include"readObj.h"
#include"../thread.h"
#include"../not used/createMesh.h"

class SetModel
{
public:
	struct RegularizationInfo
	{
		double body_center[3];
		double max_dis_from_center;
		double move_info[3];
		double scaler;
	};

	AABB aabb;
	OriMesh ori_mesh;
	RegularizationInfo regularization_info;
	void load_getAABB(std::string& path);
	void regularization(RegularizationInfo& regularization_info);

private:

	void setBackMaterial(OriMesh& ori_mesh);
	void moveBodyCapsule();
	void getAABB();
};

