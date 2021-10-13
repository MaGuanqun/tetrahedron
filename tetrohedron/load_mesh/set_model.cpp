#include"set_model.h"
#include"../basic/enum_setting.h"

void SetModel::load_getAABB(std::string& path)
{
	std::string extension_name;
	extension_name = path.substr(path.length() - 3, path.length());
	if (extension_name == "obj") {
		//ReadObj read_obj;
		//read_obj.load(path.c_str(), ori_mesh);
		ori_mesh.type == TRIANGLE;
	CreateMesh create_mesh(ori_mesh);
	//if (path == "./model/Avatar_low.obj") {
	//	create_mesh.setSphere(ori_mesh);
	//	//create_mesh.setFloor(ori_mesh);
	//}
	//else if (path == "./model/floor.obj") {
	//	//create_mesh.setCapsule(ori_mesh); 
	//	create_mesh.setFloor(ori_mesh);
	//}
	//else if (path == "./model/Tshirts_simpleModel_high.obj") {
		create_mesh.setMaterial1(ori_mesh);
	//}
	//else if (path == "./model/Outer_simpleModel_high.obj") {
	//	create_mesh.setMaterial2(ori_mesh);
	//}
	//else if (path == "./model/Pants_simpleModel_high.obj") {
	//	create_mesh.setMaterial3(ori_mesh);
	//}
	}
	else {
		ori_mesh.type == TETROHEDRON;
	}
	
	getAABB();
	setBackMaterial(ori_mesh);
	if (path == "./model/floor.obj") {
		//moveBodyCapsule();
	}
	if (path == "./model/Avatar_low.obj") {
		//moveBodyCapsule();
		//moveSphere();
		//moveSphere();
	}
}

void SetModel::setBackMaterial(OriMesh& ori_mesh){
	ori_mesh.back_material.Ka[0] = ori_mesh.front_material.Ka[1] * 0.98;
	ori_mesh.back_material.Ka[1] = ori_mesh.front_material.Ka[2] * 0.98;
	ori_mesh.back_material.Ka[2] = ori_mesh.front_material.Ka[0] * 0.98;

	ori_mesh.back_material.Kd[0] = ori_mesh.front_material.Kd[1] * 0.98;
	ori_mesh.back_material.Kd[1] = ori_mesh.front_material.Kd[2] * 0.98;
	ori_mesh.back_material.Kd[2] = ori_mesh.front_material.Kd[0] * 0.98;

	ori_mesh.back_material.Ks[0] = ori_mesh.front_material.Ks[1] * 0.98;
	ori_mesh.back_material.Ks[1] = ori_mesh.front_material.Ks[2] * 0.98;
	ori_mesh.back_material.Ks[2] = ori_mesh.front_material.Ks[0] * 0.98;

	ori_mesh.back_material.illum = ori_mesh.front_material.illum;
	memcpy(ori_mesh.back_material.Tf, ori_mesh.front_material.Tf, 12);
	ori_mesh.back_material.Ni = ori_mesh.front_material.Ni;
	ori_mesh.back_material.Ns = ori_mesh.front_material.Ns;
}

void SetModel::getAABB()
{
	memcpy(aabb.max, &ori_mesh.vertices[0], 24);
	memcpy(aabb.min, &ori_mesh.vertices[0], 24);
	for (int i = 1; i < ori_mesh.vertices.size(); ++i) {
		for (int j = 0; j < 3; ++j) {
			aabb.max[j] = myMax(aabb.max[j], ori_mesh.vertices[i][j]);
			aabb.min[j] = myMin(aabb.min[j], ori_mesh.vertices[i][j]);
		}
	}


}

void SetModel::regularization(RegularizationInfo& regularization_info)
{
	for (int i = 0; i < ori_mesh.vertices.size(); ++i) {
		SUB(ori_mesh.vertices[i], ori_mesh.vertices[i], regularization_info.body_center);
	}
	for (int i = 0; i < ori_mesh.vertices.size(); ++i) {
		MULTI(ori_mesh.vertices[i], ori_mesh.vertices[i], regularization_info.scaler);
	}
	for (int i = 0; i < ori_mesh.vertices.size(); ++i) {
		SUM(ori_mesh.vertices[i], ori_mesh.vertices[i], regularization_info.move_info);
	}
}