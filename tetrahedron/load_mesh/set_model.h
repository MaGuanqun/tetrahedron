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
	void load(std::string& path, int index, bool collider);
	void regularization(int obj_index,
		std::array<double, 3>& translation_info, double& resize_scaler, std::array<double, 3>& rotate_angle);
	void getAABB();
private:

	void setBackMaterial(OriMesh& ori_mesh);
	void moveBodyCapsule(OriMesh& ori_mesh, unsigned int obj_No, bool collider);

	void setTetFrontMaterial(OriMesh& ori_mesh, int index);
	void splitPath(std::string& path, std::string& name);
};

