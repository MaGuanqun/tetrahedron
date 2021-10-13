#include"preprocessing.h"

void Preprocessing::load_all_model(std::vector<std::string>& body_path, std::vector<std::string>& cloth_path)
{
	simulation_model.resize(cloth_path.size());
	for (int i = 0; i < cloth_path.size(); ++i) {
		simulation_model[i].load_getAABB(cloth_path[i]);
	}
	if (!body_path.empty()) {
		collider_model.resize(body_path.size());
		for (int i = 0; i < body_path.size(); ++i) {
			collider_model[i].load_getAABB(body_path[i]);
		}
	}
	//getRegularizationInfo();
	getPresetRegularizationInfo();
	for (int i = 0; i < cloth_path.size(); ++i) {
		simulation_model[i].regularization(regularization_info);
	}
	if (!body_path.empty()) {
		for (int i = 0; i < body_path.size(); ++i) {
			collider_model[i].regularization(regularization_info);
		}
	}
	scene_info.max_dis_from_center = regularization_info.max_dis_from_center;
	memcpy(scene_info.camera_center, regularization_info.move_info, 24);
	getFinalMesh();
}

void Preprocessing::getFinalMesh()
{
	if (!collider_model.empty()) {
		ori_collider_mesh.resize(collider_model.size());
		for (int i = 0; i < ori_collider_mesh.size(); ++i) {
			ori_collider_mesh[i] = collider_model[i].ori_mesh;
		}
	}
	ori_simulation_mesh.resize(simulation_model.size());
	for (int i = 0; i < simulation_model.size(); ++i) {
		ori_simulation_mesh[i] = simulation_model[i].ori_mesh;
	}
	//write_obj.write(ori_cloth_mesh);
}

void Preprocessing::getRegularizationInfo()
{
	double max_pos[3];
	double min_pos[3];
	memcpy(max_pos, simulation_model[0].aabb.max, 24);
	memcpy(min_pos, simulation_model[0].aabb.min, 24);
	for (int i = 1; i < simulation_model.size(); ++i) {
		for (int j = 0; j < 3; ++j) {
			max_pos[j] = myMax(max_pos[j], simulation_model[i].aabb.max[j]);
			min_pos[j] = myMin(min_pos[j], simulation_model[i].aabb.min[j]);
		}
	}
	//std::cout << max_pos[0] << " " << max_pos[1] << std::endl;
	if (!collider_model.empty()) {
		for (int i = 0; i < collider_model.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				max_pos[j] = myMax(max_pos[j], collider_model[i].aabb.max[j]);
				min_pos[j] = myMin(min_pos[j], collider_model[i].aabb.min[j]);
			}
		}
		//std::cout << max_pos[0] << " " << max_pos[1] << std::endl;
	}

	for (int i = 0; i < 3; ++i) {
		regularization_info.body_center[i] = 0.5 * (max_pos[i] + min_pos[i]);
	}
	double max_cen = myMax(myMax((max_pos[0] - min_pos[0]), (max_pos[1] - min_pos[1])), (max_pos[2] - min_pos[2]));
	regularization_info.scaler = 1.0 / max_cen;
	regularization_info.move_info[0] = 0.8;
	regularization_info.move_info[1] = 0.8;
	regularization_info.move_info[2] = 0.8;
	regularization_info.max_dis_from_center = 1.2 * sqrt(3.0);
}

void Preprocessing::getPresetRegularizationInfo()
{
	//for (int i = 0; i < 3; ++i) {
	regularization_info.body_center[0] = 0.0;
	regularization_info.body_center[1] = 0.0;
	regularization_info.body_center[2] = 0.0;
	//}
	regularization_info.scaler = 1.0;
	regularization_info.move_info[0] = 1.5;
	regularization_info.move_info[1] = 1.5;
	regularization_info.move_info[2] = 1.5;
	regularization_info.max_dis_from_center = 3.0 * sqrt(3.0);
}