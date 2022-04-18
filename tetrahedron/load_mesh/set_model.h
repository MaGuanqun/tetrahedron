#pragma once

#include"../basic/aabb.h"
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

	std::array<double,6> aabb;
	OriMesh ori_mesh;
	RegularizationInfo regularization_info;
	void load_getAABB(std::string& path, int& index, int obj_index);
	void regularization(RegularizationInfo& regularization_info);

private:

	void setBackMaterial(OriMesh& ori_mesh);
	void moveBodyCapsule(OriMesh& ori_mesh);
	void getAABB();
	void setTetFrontMaterial(OriMesh& ori_mesh, int& index);
	void splitPath(std::string& path, std::string& name);
};

