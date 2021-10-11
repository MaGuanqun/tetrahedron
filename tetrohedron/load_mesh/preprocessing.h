#pragma once
#include"./set_model.h"

class Preprocessing {
private:
	struct SceneInfo
	{
		double max_dis_from_center;
		double camera_center[3];
	};

	std::vector<SetModel>collider_model;
	std::vector<SetModel>simulation_model;
	void getFinalMesh();
	SetModel::RegularizationInfo regularization_info;
	void getRegularizationInfo();
	void getPresetRegularizationInfo();

public:
	SceneInfo scene_info;
	std::vector<OriMesh> ori_simulation_mesh;
	std::vector<OriMesh> ori_collider_mesh;
	void load_all_model(std::vector<std::string>& body_path, std::vector<std::string>& cloth_path);
};