#include"preprocessing.h"

void Preprocessing::load_all_model(std::vector<std::string>& body_path, std::vector<std::string>& object_path)
{
	int tet_index = 0;
	simulation_model.resize(object_path.size());
	for (int i = 0; i < object_path.size(); ++i) {
		simulation_model[i].load_getAABB(object_path[i], tet_index,i,false);
	}
	if (!body_path.empty()) {
		collider_model.resize(body_path.size());
		for (int i = 0; i < body_path.size(); ++i) {
			collider_model[i].load_getAABB(body_path[i], tet_index,i,true);
		}
	}
	getRegularizationInfo();
	//getPresetRegularizationInfo();
	for (int i = 0; i < object_path.size(); ++i) {
		simulation_model[i].regularization(regularization_info);
	}
	//if (!body_path.empty()) {
	//	for (int i = 0; i < body_path.size(); ++i) {
	//		collider_model[i].regularization(regularization_info);
	//	}
	//}
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
	double aabb_pos[6];
	memcpy(aabb_pos, simulation_model[0].aabb.data(), 48);
	for (int i = 1; i < simulation_model.size(); ++i) {
		for (int j = 0; j < 3; ++j) {
			aabb_pos[j] = myMin(aabb_pos[j], simulation_model[i].aabb[j]);
		}
		for (int j = 3; j < 6; ++j) {
			aabb_pos[j] = myMax(aabb_pos[j], simulation_model[i].aabb[j]);
		}
	}
	if (!collider_model.empty()) {
		for (int i = 0; i < collider_model.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				aabb_pos[j] = myMin(aabb_pos[j], collider_model[i].aabb[j]);
			}
			for (int j = 3; j < 6; ++j) {
				aabb_pos[j] = myMax(aabb_pos[j], collider_model[i].aabb[j]);
			}
		}
	}

	for (int i = 0; i < 3; ++i) {
		regularization_info.body_center[i] = 0.5 * (aabb_pos[i] + aabb_pos[i+3]);
	}
	double max_cen = myMax(myMax((aabb_pos[3] - aabb_pos[0]), (aabb_pos[4] - aabb_pos[1])), (aabb_pos[5] - aabb_pos[2]));
	regularization_info.scaler = 1.0;// 1.0 / max_cen;//;
	regularization_info.move_info[0] = 0.0;
	regularization_info.move_info[1] = 0.0;
	regularization_info.move_info[2] = 0.0;
	regularization_info.max_dis_from_center = 1.2 * sqrt(3.0* max_cen);

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