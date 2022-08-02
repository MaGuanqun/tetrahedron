#include"set_model.h"
#include"../basic/enum_setting.h"
#include"readObj.h"
#include"read_ele.h"

void SetModel::load_getAABB(std::string& path, int& index, int obj_index)
{
	std::string extension_name;
	extension_name = path.substr(path.length() - 3, 3);
	if (extension_name == "obj") {
	

		ori_mesh.type = TRIANGLE;
		CreateMesh create_mesh(ori_mesh);
		std::string name;
		splitPath(path, name);
		std::cout << name << std::endl;
		//if (name == "Avatar_low.obj") {
		//	//create_mesh.setSphere(ori_mesh);
			//create_mesh.setFloor(ori_mesh);
		//	create_mesh.setMaterial1(ori_mesh);
		//}
		if (name == "floor.obj") {
			create_mesh.setFloor(ori_mesh);
		}
		else if (name == "test_cloth.obj") {
			create_mesh.setMaterial1(ori_mesh);
		}
		else if (name == "capsule.obj") {
			create_mesh.setCapsule(ori_mesh);
		}
		else if (name == "sphere.obj") {
			create_mesh.setSphere(ori_mesh);
		}
		else if (name == "skirt_outer.obj") {
			create_mesh.setCylinder(ori_mesh, 0.18, 0.36, 20, 7, 1);
		}
		else if (name == "skirt_inner.obj") {
			create_mesh.setMaterial2(ori_mesh);
		}
		else {
			ReadObj read_obj;
			read_obj.load(path.c_str(), ori_mesh);
		}
	//else if (path == "./model/Outer_simpleModel_high.obj") {
	//	create_mesh.setMaterial2(ori_mesh);
	//}
	//else if (path == "./model/Pants_simpleModel_high.obj") {
	//	create_mesh.setMaterial3(ori_mesh);
	//}
	//		moveBodyCapsule(ori_mesh);
	}
	else {
		ori_mesh.type = TETRAHEDRON;
		ReadEle read_ele;
		read_ele.load(path.c_str(), ori_mesh);
		setTetFrontMaterial(ori_mesh, index);

		//if (obj_index == 1) {
		//	moveBodyCapsule(ori_mesh);
		//}

	}
	
	getAABB();
	setBackMaterial(ori_mesh);
	if (path == "./model/floor.obj") {
		//moveBodyCapsule();
	}
	if (path == "./model/Avatar_low.obj") {
		//
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
	memcpy(aabb.data(), &ori_mesh.vertices[0], 24);
	memcpy(aabb.data()+3, &ori_mesh.vertices[0], 24);
	for (int i = 1; i < ori_mesh.vertices.size(); ++i) {
		for (int j = 0; j < 3; ++j) {
			if (aabb[j] > ori_mesh.vertices[i][j]) {
				aabb[j] = ori_mesh.vertices[i][j];
			}
			if (aabb[j + 3] < ori_mesh.vertices[i][j]) {
				aabb[j + 3] = ori_mesh.vertices[i][j];
			}
		}
	}

	//std::cout << aabb[3]-aabb[0] << " " << aabb[4] - aabb[1] << " " << aabb[5] - aabb[2] << std::endl;
	//std::cout << 0.5*(aabb[3]+aabb[0]) << " " << 0.5*(aabb[4] + aabb[1]) << " " << 0.5*(aabb[5] + aabb[2]) << std::endl;
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
	//std::cout <<"scaler "<< regularization_info.scaler << std::endl;
}

void SetModel::setTetFrontMaterial(OriMesh& ori_mesh, int& index)
{
	std::array<float, 3> Kd, Ka, Ks, Tf;
	int illum;
	float Ni, Ns;
	Tf = { 1.00, 1.00, 1.00 };
	Ni = 1.00;
	Ns = 28.76;
	illum = 4;
	switch (index)
	{
	case 0:
		Kd = { 0.235, 0.5, 0.605 };		
		Ka = { 0.1175, 0.25, 0.3025 };
		Ks = { 0.094, 0.2, 0.242 };
		break;
	case 1:
		Kd = { 0.3, 0.1, 0.6 };
		Ka = { 0.15, 0.05, 0.3 };
		Ks = { 0.12, 0.04, 0.24 };
		break;
	default:
		Kd = { 0.5, 0.5, 0.5 };
		Ka = { 0.02, 0.02, 0.02 };
		Ks = { 0.24, 0.24, 0.24 };
	}
	memcpy(ori_mesh.front_material.Ka, Ka.data(), 12);
	memcpy(ori_mesh.front_material.Kd, Kd.data(), 12);
	memcpy(ori_mesh.front_material.Ks, Ks.data(), 12);
	memcpy(ori_mesh.front_material.Tf, Tf.data(), 12);
	ori_mesh.front_material.Ns = Ns;
	ori_mesh.front_material.Ni = Ni;
	ori_mesh.front_material.illum = illum;
	index++;
}


void SetModel::splitPath(std::string& path, std::string& name)
{
	int index = path.find_last_of("\\");
	name = path.substr(index+1, path.length()-index);
}

void SetModel::moveBodyCapsule(OriMesh& ori_mesh)
{
	//band capsule
	double move[3] = { 1.5,-2.8,-0.3 };//
	//double move[3] = { 0.0,-0.9,-0.3 };//this is for two prisms
	//double move[3] = { -60, -130,-30 };//this is for two dragons
	// move capsule
	//double move[3] = { 0.0,-0.3,-0.35 };
	//sphere
	//double move[3] = { 0.0,-0.15,0.0 };
	//skirt
	//double move[3] = { 0.0,-0.3,0.0 };
	//capsule collide two clothes
	//rotate(0.5*M_PI, mesh_struct.mesh.vertices);
	//double move[3] = { 0.15,0.0,0.0 };
	for (int i = 0; i < ori_mesh.vertices.size(); ++i) {
		SUM(ori_mesh.vertices[i], ori_mesh.vertices[i], move);
	}



}