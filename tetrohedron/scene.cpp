#include"scene.h"
#include<bitset>

Scene::Scene()
{
	light.ambient = glm::vec3(1.0, 1.0, 1.0);
	light.diffuse = glm::vec3(0.8, 0.8, 0.8);
	light.specular = glm::vec3(0.95, 0.95, 0.95);

	
	time_step = 1.0 / 50.0;
	project_dynamic.time_step = time_step;

	max_force_magnitude = 12.0;

	last_output_obj_stamp = -1;
	time_stamp = 0;

}


void Scene::obtainConvergenceRate(double* convergence_rate)
{
	convergence_rate[LOCAL_GLOBAL] = project_dynamic.local_global_conv_rate;
	convergence_rate[OUTER] = project_dynamic.outer_itr_conv_rate;
}

void Scene::updateConvRate(double* convergence_rate)
{
	project_dynamic.outer_itr_conv_rate = convergence_rate[OUTER];
	project_dynamic.local_global_conv_rate = convergence_rate[LOCAL_GLOBAL];
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
	project_dynamic.setForPD(&cloth, &tetrohedron, &thread);
	collision.initial(&cloth, &collider, &tetrohedron, &thread);
	setAveEdgeLength();
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

void Scene::initialCloth()
{
	project_dynamic.initial();
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].initial();
	}
	intersection.initialIntersection();
	//spatial_hashing.body_time_stamp++;
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].reset();
	}
	time_stamp = 0;
}

void Scene::resetCloth()
{
	project_dynamic.reset();
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].reset();
	}
	intersection.initialIntersection();
	//spatial_hashing.time_stamp++;
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].reset();
	}
	time_stamp = 0;
}

void Scene::initialIntersection()
{
	intersection.initialIntersection();
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].initialMouseChosenVertex();
	}
}

void Scene::updateCloth(Camera* camera, double* cursor_screen, bool* control_parameter, float force_coe)
{


	if (intersection.happened && !control_parameter[START_TEST]) {
	
		setCursorForce(camera,cursor_screen,force_coe);
	}

	if (control_parameter[START_SIMULATION] || control_parameter[ONE_FRAME]) {
		project_dynamic.PDsolve();
		time_stamp++;
	}

	updateBuffer();
}

void Scene::updateBuffer()
{
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].setBuffer();
	}
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].setBuffer();
	}
	for (int i = 0; i < tetrohedron.size(); ++i) {
		tetrohedron[i].setBuffer();
	}
}


void Scene::setAveEdgeLength()
{
	double edge_length_temp = 0.0;
	int edge_size=0;
	for (int i = 0; i < cloth.size(); ++i) {
		for (int j = 0; j < cloth[i].mesh_struct.edges.size(); ++j) {
			edge_length_temp += cloth[i].mesh_struct.edges[j].length;
		}
		edge_size += cloth[i].mesh_struct.edges.size();
	}

	ave_edge_length = edge_length_temp / (double)edge_size;
}

void Scene::setTolerance(double* tolerance_ratio)
{
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].setTolerance(tolerance_ratio, ave_edge_length);
	}
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].setTolerance(tolerance_ratio, ave_edge_length);
	}
}


void Scene::testBVH()
{
/*	collision.buildBVH();
	std::vector<std::vector<int>>cloth_neighbor_index;
	std::vector<std::vector<int>>collider_neighbor_index;
	cloth_neighbor_index.resize(cloth.size());
	collision.searchTriangle(cloth[0].triangle_AABB[60], 60, &cloth_neighbor_index, &collider_neighbor_index);
	for (int i = 0; i < cloth_neighbor_index[0].size(); ++i) {
		std::cout << cloth_neighbor_index[0][i] << std::endl;
	}*/	


	std::bitset<64>t1(0x1f00000000ffff);
	std::cout << t1 << std::endl;
	std::bitset<64>t2(0x1f0000ff0000ff);
	std::cout << t2 << std::endl;
	std::bitset<64>t3(0x100f00f00f00f00f);
	std::cout << t3 << std::endl;
	std::bitset<64>t4(0x10c30c30c30c30c3);
	std::cout << t4 << std::endl;
	std::bitset<64>t5(0x1249249249249249);
	std::cout << t5 << std::endl;
}

void Scene::obtainCursorPosition(double* pos, Camera* camera, std::vector<std::vector<bool>>& hide)
{
	int mouse_pos[2];
	mouse_pos[0] = pos[0];
	mouse_pos[1] = pos[1];
	int chosen_index[2];
	pick_triangle.pickTriangle(&cloth, &collider, camera, hide, chosen_index, mouse_pos);
	intersection.initialIntersection();
	if (chosen_index[0] > -1) {
		double cursor_pos[3];
		intersection.setIntersection(chosen_index);
		getCursorPos(cursor_pos, cloth[intersection.cloth_No].mesh_struct.vertex_position,
			cloth[intersection.cloth_No].mesh_struct.faces[intersection.face_index].vertex);
		cloth[chosen_index[1]].findAllNeighborVertex(chosen_index[0], cursor_pos, ave_edge_length);
	}
}

void Scene::getCursorPos(double* cursor_pos, std::vector<std::array<double,3>>& vertex, int* vertex_index)
{
	for (int i = 0; i < 3; ++i) {
		cursor_pos[i] = (vertex[vertex_index[0]][i] + vertex[vertex_index[1]][i] + vertex[vertex_index[2]][i]) / 3.0;
	}
}

void Scene::setCursorForce(Camera* camera, double* cursor_screen, float force_coe)
{
	double cursor_pos[3];
	double force_direction[3];
	getCursorPos(cursor_pos, cloth[intersection.cloth_No].mesh_struct.vertex_position,
		cloth[intersection.cloth_No].mesh_struct.faces[intersection.face_index].vertex);
	cursor.translate(cursor_pos);
	cursorMovement(camera, cursor_screen, force_direction, force_coe, cursor_pos);
	//project_dynamic.add
}

void Scene::cursorMovement(Camera* camera, double* cursor_screen, double* force_direction, float force_coe, double* object_position)
{
	double cursor_pos_in_space[3];
	camera->getCursorPosInSpace(cursor_pos_in_space, cursor_screen, object_position);
	SUB(force_direction, cursor_pos_in_space, object_position);
	MULTI(force_direction, force_direction, force_coe);
	double force_magnitude = sqrt(DOT(force_direction, force_direction));
	if (force_magnitude > max_force_magnitude) {
		force_magnitude = max_force_magnitude / force_magnitude;
		MULTI(force_direction, force_direction, force_magnitude);
	}
	//std::cout << "force mag " << force_direction[0] << " " << force_direction[1] << " " << force_direction[2] << std::endl;
}