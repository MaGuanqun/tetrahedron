#include"scene.h"


Scene::Scene()
{
	light.ambient = glm::vec3(1.0, 1.0, 1.0);
	light.diffuse = glm::vec3(0.8, 0.8, 0.8);
	light.specular = glm::vec3(0.95, 0.95, 0.95);

	
	time_step = 1.0 / 50.0;
	//project_dynamic.time_step = time_step;
	//spatial_hashing.time_step = time_step;
	//max_force_magnitude = 12.0;

	last_output_obj_stamp = -1;
	//project_dynamic.setSceneAddress(this);
	time_stamp = 0;
}


void Scene::loadMesh(std::vector<std::string>& collider_path, std::vector<std::string>& object_path)
{
	Preprocessing preprocessing;
	preprocessing.load_all_model(collider_path, object_path);
	initialSceneSetting(preprocessing);
	collider_num = collider_path.size();
	collider.resize(collider_num);
	for (int i = 0; i < collider_num; ++i) {
		collider[i].loadMesh(preprocessing.ori_collider_mesh[i], &thread);
	}
	for (int i = 0; i < object_path.size(); ++i) {
		if (preprocessing.ori_simulation_mesh[i].type == TRIANGLE) {
			cloth_index_in_object.push_back(i);
		}
		else {
			tetrohedron_index_in_object.push_back(i);
		}
	}
	cloth_num = cloth_index_in_object.size();
	tetrohedron_num = tetrohedron_index_in_object.size();
	cloth.resize(cloth_num);
	tetrohedron.resize(tetrohedron_num);
	double cloth_density=25.0;
	for (int i = 0; i < cloth_num; ++i) {
		cloth[i].loadMesh(preprocessing.ori_simulation_mesh[cloth_index_in_object[i]], cloth_density, &thread);
	}
	//for (int i = 0; i < tetrohedron_num; ++i) {
	//	tetrohedron[i].loadMesh(preprocessing.ori_simulation_mesh[cloth_index[i]], 25, &thread);
	//}
	setWireframwColor();
	std::vector<SingleClothInfo> single_cloth_info;
	std::array<double, 4>collision_stiffness_per = { 5e3,5e3,5e3,5e3 };// stiffness of collision constraint //=0 body point triangle, =1 point-triangle =2 edge-edge =3 point-point,
	std::vector<std::array<double, 4>>collision_stiffness(cloth_num, collision_stiffness_per);
	for (int i = 0; i < cloth_num; ++i) {
		single_cloth_info.push_back(SingleClothInfo(cloth_density, 5e3, 1e5, 1e-4, collision_stiffness[i].data(), 0.5, 0.4, collision_stiffness_per[1]));
	}
	for (int i = 0; i < cloth_num; ++i) {
		cloth[i].recordInitialMesh(single_cloth_info[i]);
	}

}

void Scene::setWireframwColor()
{
	double color3[4][3] = { {1.0,1.0,0.0},{0.0,0.0,1.0},{0.0,0.3,0.0},{0.0, 1.0,1.0} };
	double** color0;
	color0 = new double* [collider_num];
	for (int i = 0; i < collider_num; ++i)
	{
		color0[i] = new double[3];
		memcpy(color0[i], color3[0], 24);
	}
	double** color1;
	color1 = new double* [cloth_num+ tetrohedron_num];
	for (int i = 0; i < cloth_num + tetrohedron_num; ++i)
	{
		color1[i] = new double[3];
		memcpy(color1[i], color3[(i + 1) % 4], 24);
	}
	for (int i = 0; i < cloth_num; ++i) {
		cloth[i].setWireframwColor(color1[i]);//
	}
	for (int i = 0; i < collider_num; ++i) {
		collider[i].setWireframwColor(color0[i]);
	}
	for (int i = 0; i < tetrohedron_num; ++i) {
		tetrohedron[i + cloth_num].setWireframwColor(color1[i + cloth_num]);
	}

}

void Scene::initialSceneSetting(Preprocessing& preprocessing)
{
	shadow.camera_from_origin = 0.8 * preprocessing.scene_info.max_dis_from_center;//0.8
	shadow.far_plane = 2.0 * preprocessing.scene_info.max_dis_from_center;
	memcpy(camera_center, preprocessing.scene_info.camera_center, 24);
	shadow.setBasic();

}


void Scene::getClothInfo(std::vector<std::array<int, 3>>& mesh_info, std::vector<double>& mass, std::vector<std::array<double, 3>>& mesh_stiffness, double* simulation_parameter, std::vector<std::array<double, 4>>& collision_stiffness)
{
	mesh_info.resize(2);
	mass.resize(cloth_num);
	mesh_info.resize(cloth_num);
	for (int i = 0; i < cloth_num; ++i) {
		mesh_info[i][0] = cloth[i].mesh_struct.vertices.size();
		mesh_info[i][1] = cloth[i].mesh_struct.edges.size();
		mesh_info[i][2] = cloth[i].mesh_struct.faces.size();
		mass[i] = cloth[i].mass;
	}
	mesh_stiffness.resize(cloth_num);
	for (int i = 0; i < cloth_num; ++i) {
		mesh_stiffness[i][0] = cloth[i].single_cloth_info_ref.length_stiffness;
		mesh_stiffness[i][1] = cloth[i].single_cloth_info_ref.bending_stiffness;
		mesh_stiffness[i][2] = cloth[i].single_cloth_info_ref.position_stiffness;
	}

	simulation_parameter[0] = time_step;
	//simulation_parameter[1] = project_dynamic.gravity_;
	collision_stiffness.resize(cloth_num);
	for (int i = 0; i < cloth_num; ++i) {
		memcpy(collision_stiffness[i].data(), cloth[i].single_cloth_info_ref.collision_stiffness, 32);
	}

}

void Scene::drawScene(Camera* camera, std::vector<std::vector<bool>>& wireframe, std::vector<std::vector<bool>>& hide, bool start_save_obj)
{
	shadow.drawShadow(camera, hide,cloth,collider,tetrohedron,cloth_index_in_object,tetrohedron_index_in_object);
	glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_CUBE_MAP, shadow.depth_map);
	glEnable(GL_CULL_FACE);
	for (int j = 0; j < cloth_num; ++j) {
		if (!hide[1][cloth_index_in_object[j]]) {
			cloth[j].setSceneShader(light, camera, shadow.far_plane);
			cloth[j].draw(camera);
		}
	}
	for (int j = 0; j < tetrohedron_num; ++j) {
		if (!hide[1][tetrohedron_index_in_object[j]]) {
			//tetrohedron[j].setSceneShader(light, camera, shadow.far_plane);
			//tetrohedron[j].draw(camera);
		}
	}
	glDisable(GL_CULL_FACE);


	for (int j = 0; j < collider.size(); ++j) {
		if (!hide[0][j]) {
			collider[j].setSceneShader(light, camera, shadow.far_plane);
			collider[j].draw(camera);
		}
	}

	for (int i = 0; i < collider.size(); ++i) {
		if (wireframe[0][i]) {
			collider[i].drawWireframe(camera);
		}
	}
	for (int j = 0; j < cloth_num; ++j) {
		if (wireframe[1][cloth_index_in_object[j]]) {
			cloth[j].drawWireframe(camera);
		}
	}
	for (int j = 0; j < tetrohedron_num; ++j) {
		if (wireframe[1][tetrohedron_index_in_object[j]]) {
			tetrohedron[j].drawWireframe(camera);
		}
	}

	//drawTriangle1.draw(camera, glm::vec3(0.0, 1.0, 0.0));
	if (intersection.happened) {
		cursor.draw(camera);
	}

	if (start_save_obj) {
		saveObj();
	}

}

void Scene::saveObj()
{
	if (time_stamp != last_output_obj_stamp) {
		save_obj.write(cloth, 10, time_stamp);
		last_output_obj_stamp = time_stamp;
	}
}

