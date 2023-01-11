#pragma once
#include"./set_model.h"

class Preprocessing {
private:
	struct SceneInfo
	{
		double max_dis_from_center;
		double camera_center[3];
	};

	void getFinalMesh();
	SetModel::RegularizationInfo regularization_info;
	void getRegularizationInfo();
	void getPresetRegularizationInfo();

public:
	SceneInfo scene_info;

	std::vector<SetModel>collider_model;
	std::vector<SetModel>simulation_model;

	//std::vector<OriMesh> ori_simulation_mesh;
	//std::vector<OriMesh> ori_collider_mesh;
	void load_all_model(std::vector<std::string>& body_path, std::vector<std::string>& object_path,
		std::vector<std::array<double, 3>>& translation_info, std::vector<double>& resize_scaler, std::vector<std::array<double, 3>>& rotate_angle,
		std::vector<std::array<double, 3>>& collider_translation_info, std::vector<double>& collider_resize_scaler, std::vector<std::array<double, 3>>& collider_rotate_angle);
};