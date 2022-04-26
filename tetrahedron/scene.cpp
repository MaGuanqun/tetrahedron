#include"scene.h"
#include<bitset>

Scene::Scene()
{
	light.ambient = glm::vec3(1.0, 1.0, 1.0);
	light.diffuse = glm::vec3(0.8, 0.8, 0.8);
	light.specular = glm::vec3(0.95, 0.95, 0.95);


	time_step = 1.0 / 100.0;
	project_dynamic.time_step = time_step;

	max_force_magnitude = 2.0;

	last_output_obj_stamp = -1;
	time_stamp = 0;
	genShader();
	project_dynamic.collision.time_stamp = &time_stamp;
	project_dynamic.time_stamp = &time_stamp;

	//test_array_size = 10000000;
	//iteration_num_for_test_array = 1000;
	//test_pair = new unsigned int* [thread.thread_num];
	//test_pair_single = new unsigned int* [thread.thread_num];
	//for (unsigned int i = 0; i < thread.thread_num; ++i)
	//{
	//	test_pair[i] = new unsigned int[test_array_size];
	//	test_pair_single[i] = new unsigned int[test_array_size];
	//}
	//voidForWritetToArraySingle(thread.thread_num);

	//thread.assignTask(this, TEST_ARRAY);
	//time_t t1 = clock();
	//time_t t = clock();
	//for (unsigned int i = 0; i < 2; ++i) {
	//testForWritetToArraySingle(thread.thread_num);
	//thread.assignTask(this, TEST_ARRAY);


	//t = clock();
	//for (unsigned int i = 0; i < 10; ++i) {
	//	testForWritetToArraySingle(thread.thread_num);
	//}
	//t1 = clock();
	//std::cout<< "single thread " << (t1 - t)/1 << std::endl;
	//t = clock();
	//for (unsigned int i = 0; i < 10; ++i) {
	//	thread.assignTask(this, TEST_ARRAY);
	//}
	//t1 = clock();
	//std::cout << "multi thread " << (t1 - t) / 1 << std::endl;
	//compareArray();

	//std::cout << cbrt(-27.0) << " " << cbrt(27.0) << std::endl;
}


void Scene::compareArray()
{
	for (unsigned int i = 0; i < thread.thread_num; ++i) {
		for (unsigned int j = 0; j < test_array_size; ++j) {
			if (test_pair[i][j] != test_pair_single[i][j]) {
				std::cout << "not equal " << test_pair[i][j] << " " << test_pair_single[i][j] << std::endl;
			}
		}
	}
}


void Scene::obtainConvergenceInfo(double* convergence_rate, int* iteration_num)
{
	convergence_rate[LOCAL_GLOBAL] = project_dynamic.local_global_conv_rate;
	convergence_rate[OUTER] = project_dynamic.outer_itr_conv_rate;
	iteration_num[LOCAL_GLOBAL] = project_dynamic.local_global_iteration_num;
	iteration_num[OUTER] = project_dynamic.outer_iteration_num;
}

void Scene::updateConvRate(double* convergence_rate)
{
	project_dynamic.outer_itr_conv_rate = convergence_rate[OUTER];
	project_dynamic.local_global_conv_rate = convergence_rate[LOCAL_GLOBAL];
}

void Scene::loadMesh(std::vector<std::string>& collider_path, std::vector<std::string>& object_path, double* tolerance_ratio)
{
	Preprocessing preprocessing;
	preprocessing.load_all_model(collider_path, object_path);

	initialSceneSetting(preprocessing);
	collider_num = collider_path.size();
	collider.resize(collider_num);
	for (int i = 0; i < collider_num; ++i) {
		collider[i].loadMesh(preprocessing.ori_collider_mesh[i], &thread);
	}
	std::vector<int>cloth_index_in_object, tetrahedron_index_in_object;
	for (int i = 0; i < object_path.size(); ++i) {
		if (preprocessing.ori_simulation_mesh[i].type == TRIANGLE) {
			cloth_index_in_object.push_back(i);
		}
		else {
			tetrahedron_index_in_object.push_back(i);
		}
	}

	total_obj_num = collider_path.size() + object_path.size();
	obj_num_except_collider = object_path.size();

	cloth_num = cloth_index_in_object.size();
	tetrahedron_num = tetrahedron_index_in_object.size();
	cloth.resize(cloth_num);
	tetrahedron.resize(tetrahedron_num);
	setTolerance(tolerance_ratio);
	double cloth_density = 15.0;
	double tetrahedron_density = 0.1;
	for (int i = 0; i < cloth_num; ++i) {
		cloth[i].loadMesh(preprocessing.ori_simulation_mesh[cloth_index_in_object[i]], cloth_density, &thread);
	}
	for (int i = 0; i < tetrahedron_num; ++i) {
		tetrahedron[i].loadMesh(preprocessing.ori_simulation_mesh[tetrahedron_index_in_object[i]], tetrahedron_density, &thread);
	}
	setWireframwColor();
	std::vector<SingleClothInfo> single_cloth_info;
	std::array<double, 4>collision_stiffness_per = { 2e5,2e5,2e5, 2e5 };// stiffness of collision constraint //=0 body point triangle, =1 point-triangle =2 edge-edge =3 point-point,
	std::vector<std::array<double, 4>>collision_stiffness(cloth_num, collision_stiffness_per);
	for (int i = 0; i < cloth_num; ++i) {
		single_cloth_info.push_back(SingleClothInfo(cloth_density, 1e3, 1e6, 3e-4, collision_stiffness[i].data(), 0.5, 0.4, collision_stiffness_per[1]));
	}
	for (int i = 0; i < cloth_num; ++i) {
		cloth[i].recordInitialMesh(single_cloth_info[i]);
	}

	std::array<double, 4>tetrahedron_collision_stiffness_per = {7,8, 8,2e1 };
	double sigma_limit[2] = { 0.99,1.01 };
	SingleTetrahedronInfo single_tetrahedron_info(tetrahedron_density, 2e3, 1e2, 0.0, tetrahedron_collision_stiffness_per.data(), sigma_limit);
	for (int i = 0; i < tetrahedron_num; ++i) {
		tetrahedron[i].recordInitialMesh(single_tetrahedron_info);
	}


	project_dynamic.setForPD(&cloth, &tetrahedron, &collider, &thread, tolerance_ratio);
	setAveEdgeLength();
	cursor.createVertices(0.03, camera_center);
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
	color1 = new double* [cloth_num + tetrahedron_num];
	for (int i = 0; i < cloth_num + tetrahedron_num; ++i)
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
	for (int i = 0; i < tetrahedron_num; ++i) {
		tetrahedron[i + cloth_num].setWireframwColor(color1[i + cloth_num]);
	}

}

void Scene::initialSceneSetting(Preprocessing& preprocessing)
{
	shadow.camera_from_origin = 0.8 * preprocessing.scene_info.max_dis_from_center;//0.8
	shadow.far_plane = 10.0 * preprocessing.scene_info.max_dis_from_center;
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
	simulation_parameter[1] = project_dynamic.gravity_;
	collision_stiffness.resize(cloth_num);
	for (int i = 0; i < cloth_num; ++i) {
		memcpy(collision_stiffness[i].data(), cloth[i].single_cloth_info_ref.collision_stiffness, 32);
	}
}


void Scene::updateIterateSolverParameter(double rate, int itr_solver_method)
{
	project_dynamic.updateIterateSolverParameter(rate);
	project_dynamic.itr_solver_method = itr_solver_method;
}

void Scene::getTetrahedronInfo(std::vector<std::array<int, 3>>& mesh_info, std::vector<double>& mass, std::vector<std::array<double, 2>>& mesh_stiffness, double* simulation_parameter, std::vector<std::array<double, 4>>& collision_stiffness)
{
	mesh_info.resize(2);
	mass.resize(tetrahedron_num);
	mesh_info.resize(tetrahedron_num);
	for (int i = 0; i < tetrahedron_num; ++i) {
		mesh_info[i][0] = tetrahedron[i].mesh_struct.vertex_position.size();
		mesh_info[i][1] = tetrahedron[i].mesh_struct.indices.size();
		mesh_info[i][2] = tetrahedron[i].mesh_struct.triangle_indices.size();
		mass[i] = tetrahedron[i].mass;
	}
	mesh_stiffness.resize(tetrahedron_num);
	for (int i = 0; i < tetrahedron_num; ++i) {
		mesh_stiffness[i][0] = tetrahedron[i].single_tetrahedron_info_ref.ARAP_stiffness;
		mesh_stiffness[i][1] = tetrahedron[i].single_tetrahedron_info_ref.position_stiffness;
	}
	simulation_parameter[0] = time_step;
	simulation_parameter[1] = project_dynamic.gravity_;
	collision_stiffness.resize(tetrahedron_num);
	for (int i = 0; i < tetrahedron_num; ++i) {
		memcpy(collision_stiffness[i].data(), tetrahedron[i].single_tetrahedron_info_ref.collision_stiffness, 32);
	}
}

void Scene::drawScene(Camera* camera, std::vector<std::vector<bool>>& wireframe, std::vector<std::vector<bool>>& hide,
	bool* control_parameter)
{
	shadow.drawShadow(camera, hide, cloth, collider, tetrahedron);
	glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_CUBE_MAP, shadow.depth_map);
	glEnable(GL_CULL_FACE);
	for (int j = 0; j < cloth_num; ++j) {
		if (!hide[CLOTH_][j]) {
			cloth[j].setSceneShader(light, camera, shadow.far_plane, object_shader_front);
			cloth[j].draw(camera, object_shader_front);
		}
	}
	glDisable(GL_CULL_FACE);
	for (int j = 0; j < tetrahedron_num; ++j) {
		if (!hide[TETRAHEDRON_][j]) {
			tetrahedron[j].setSceneShader(light, camera, shadow.far_plane, object_shader_front);
			tetrahedron[j].draw(camera, object_shader_front);
		}
	}
	for (int j = 0; j < collider.size(); ++j) {
		if (!hide[COLLIDER_][j]) {
			collider[j].setSceneShader(light, camera, shadow.far_plane, object_shader_front);
			collider[j].draw(camera, object_shader_front);
		}
	}

	for (int i = 0; i < collider.size(); ++i) {
		if (wireframe[COLLIDER_][i]) {
			collider[i].drawWireframe(camera, wireframe_shader);
		}
	}
	for (int j = 0; j < cloth_num; ++j) {
		if (wireframe[CLOTH_][j]) {
			cloth[j].drawWireframe(camera, wireframe_shader);
		}
	}
	for (int j = 0; j < tetrahedron_num; ++j) {
		if (wireframe[TETRAHEDRON_][j]) {
			tetrahedron[j].drawWireframe(camera, wireframe_shader);
		}
	}
	//project_dynamic.collision.draw_culling.drawTetTriangle(camera, wireframe_shader, light, shadow.far_plane);
	//project_dynamic.collision.draw_culling.draw(camera, shadow.far_plane);


	if (control_parameter[MOVE_OBJ]) {
		if (intersection.happened) {
			object_chosen_indicator.draw(wireframe_shader, camera);
		}

	}
	//project_dynamic.collision.draw_culling.drawCell(camera, wireframe_shader);

	//if (!draw_mesh_patch) {
	//	//drawTriangle1.draw(camera, glm::vec3(0.0, 1.0, 0.0));
		if (intersection.happened && (!control_parameter[MOVE_OBJ])) {
			cursor.draw(camera);
		}
	//}
	//else {
	//	//project_dynamic.collision.drawMeshPatch(camera);
	//}

	if (control_parameter[SAVE_OBJ]) {
		saveObj();
	}
}

void Scene::saveObj()
{
	if (time_stamp != last_output_obj_stamp) {
		for (int i = 0; i < tetrahedron.size(); ++i) {
			save_obj.write(tetrahedron[i].mesh_struct.vertex_for_render, tetrahedron[i].mesh_struct.triangle_indices, 10, time_stamp, i);
		}
		last_output_obj_stamp = time_stamp;
	}
}

void Scene::initial()
{
	project_dynamic.initial();
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].initial();
	}
	for (int i = 0; i < tetrahedron.size(); ++i) {
		tetrahedron[i].initial();
	}
	intersection.initialIntersection();
	//spatial_hashing.body_time_stamp++;
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].reset();
	}
	time_stamp = 0;
}

void Scene::reset()
{
	project_dynamic.reset();
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].reset();
	}
	for (int i = 0; i < tetrahedron.size(); ++i) {
		tetrahedron[i].reset();
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

void Scene::resetIntersectionState()
{
	intersection.happened = false;
	intersection.happened_include_collider = false;
}


void Scene::updateCloth(Camera* camera, double* cursor_screen, bool* control_parameter, float force_coe, bool& record_matrix,
	double& ave_iteration, bool mouse_is_pressed_previous_current_frame)
{
	if (control_parameter[MOVE_OBJ_SCRIPT]) {
		project_dynamic.collision.draw_culling.moveScript(2);
		control_parameter[MOVE_OBJ_SCRIPT] = false;
	}

	if (control_parameter[MOVE_OBJ]) {
		if (intersection.happened_include_collider) {
			setObjMoveInfo(camera, cursor_screen);
			project_dynamic.collision.getAABB();
			project_dynamic.updateSystemPos();
		}
		setChosenIndicator();
	}

	if (control_parameter[START_SIMULATION] || control_parameter[ONE_FRAME]) {
		project_dynamic.resetExternalForce();
		//std::cout << intersection.happened << " " << control_parameter[START_TEST] << std::endl;
		if (intersection.happened && !control_parameter[START_TEST]) {
			//std::cout << cursor_screen[0] << " " << cursor_screen[1] << " " << force_coe << std::end
			setCursorForce(camera, cursor_screen, force_coe);
		}
		project_dynamic.PDsolve();
		//project_dynamic.PD_IPC_solve(record_matrix);
		

		//project_dynamic.update_ave_iteration_record(ave_iteration);
		time_stamp++;
	}

	updateBuffer();
}


void Scene::setChosenIndicator()
{
	if (intersection.happened) {
		//std::cout << cursor_screen[0] << " " << cursor_screen[1] << " " << force_coe << std::end
		if (intersection.obj_No < cloth.size()) {
			object_chosen_indicator.updatePosition(cloth[intersection.obj_No].obj_aabb);
		}
		else if (intersection.obj_No < obj_num_except_collider) {
			object_chosen_indicator.updatePosition(tetrahedron[intersection.obj_No - cloth.size()].obj_aabb);
		}
		else {
			object_chosen_indicator.updatePosition(collider[intersection.obj_No - obj_num_except_collider].obj_aabb);
		}
	}
}

void Scene::updateBuffer()
{
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].setBuffer();
	}
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].setBuffer();
	}
	for (int i = 0; i < tetrahedron.size(); ++i) {
		tetrahedron[i].setBuffer();
	}
}


void Scene::setAveEdgeLength()
{
	double edge_length_temp = 0.0;
	int edge_size = 0;
	for (int i = 0; i < cloth.size(); ++i) {
		for (int j = 0; j < cloth[i].mesh_struct.edges.size(); ++j) {
			edge_length_temp += cloth[i].mesh_struct.edge_length[j];
		}
		edge_size += cloth[i].mesh_struct.edges.size();
	}

	for (int i = 0; i < tetrahedron.size(); ++i) {
		for (int j = 0; j < tetrahedron[i].mesh_struct.edges.size(); ++j) {
			edge_length_temp += tetrahedron[i].mesh_struct.edge_length[j];
		}
		edge_size += tetrahedron[i].mesh_struct.edges.size();
	}

	ave_edge_length = edge_length_temp / (double)edge_size;
	project_dynamic.initialDHatTolerance(ave_edge_length);
}

void Scene::setTolerance(double* tolerance_ratio)
{
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].setTolerance(tolerance_ratio, ave_edge_length);
	}
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].setTolerance(tolerance_ratio, ave_edge_length);
	}
	for (int i = 0; i < tetrahedron.size(); ++i) {
		tetrahedron[i].setTolerance(tolerance_ratio, ave_edge_length);
	}

}


void Scene::selectAnchor(bool* control_parameter, bool* select_anchor, double* screen_pos, bool press_state, bool pre_press_state, Camera* camera,
	std::vector<bool>& hide)
{
	if (select_anchor[0]) {
		control_parameter[START_SIMULATION] = false;
		if (press_state) {
			set_tetrahedron_anchor.setCorner(screen_pos, pre_press_state, tetrahedron, camera, hide);
		}
	}
	if (select_anchor[1]) {
		for (int i = 0; i < tetrahedron.size(); ++i) {
			tetrahedron[i].mesh_struct.setAnchorPosition();

		}
		select_anchor[1] = false;
		project_dynamic.updateTetrahedronAnchorVertices();
	}



}

void Scene::drawSelectRange(bool* select_anchor, bool press_state, bool pre_press_state)
{
	if (select_anchor[0]) {
		if (press_state && pre_press_state) {
			set_tetrahedron_anchor.draw();
		}
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

void Scene::obtainCursorIntersection(double* pos, Camera* camera, std::vector<std::vector<bool>>& hide)
{
	int mouse_pos[2];
	mouse_pos[0] = pos[0];
	mouse_pos[1] = SCR_HEIGHT - pos[1];
	int chosen_index[2];
	pick_triangle.pickTriangle(&cloth, &collider, &tetrahedron, camera, hide, chosen_index, mouse_pos);
	intersection.initialIntersection();
	//std::cout << chosen_index[0] << std::endl;
	if (chosen_index[0] > -1) {
		//double cursor_pos[3];
		intersection.setIntersection(chosen_index);
		if (chosen_index[1] < obj_num_except_collider) {
			intersection.happened = true;
		}
	}
}

void Scene::getCursorPos(double* cursor_pos, std::vector<std::array<double, 3>>& vertex, int* vertex_index)
{
	for (int i = 0; i < 3; ++i) {
		cursor_pos[i] = (vertex[vertex_index[0]][i] + vertex[vertex_index[1]][i] + vertex[vertex_index[2]][i]) / 3.0;
	}
}




void Scene::setObjMoveInfo(Camera* camera, double* cursor_screen)
{
	double cursor_pos[3];
	double force_direction[3];
	if (intersection.obj_No < cloth.size()) {
		getCursorPos(cursor_pos, cloth[intersection.obj_No].mesh_struct.vertex_for_render,
			cloth[intersection.obj_No].mesh_struct.triangle_indices[intersection.face_index].data());

	}
	else if (intersection.obj_No < obj_num_except_collider) {
		getCursorPos(cursor_pos, tetrahedron[intersection.obj_No - cloth.size()].mesh_struct.vertex_for_render,
			tetrahedron[intersection.obj_No - cloth.size()].mesh_struct.triangle_indices[intersection.face_index].data());
	}
	else {
		getCursorPos(cursor_pos, collider[intersection.obj_No - obj_num_except_collider].mesh_struct.vertex_for_render,
			collider[intersection.obj_No - obj_num_except_collider].mesh_struct.triangle_indices[intersection.face_index].data());
	}
	double cursor_pos_in_space[3];
	camera->getCursorPosInSpace(cursor_pos_in_space, cursor_screen, cursor_pos);
	double displacement[3];
	SUB(displacement, cursor_pos_in_space, cursor_pos);
	project_dynamic.collision.draw_culling.move(intersection.obj_No, displacement);

	/*std::cout << "original after move " << std::endl;
	std::cout << tetrahedron[1].mesh_struct.vertex_for_render[0][0] << " " << tetrahedron[1].mesh_struct.vertex_for_render[0][1] << " " << tetrahedron[1].mesh_struct.vertex_for_render[0][2] << std::endl;
	std::cout << tetrahedron[1].mesh_struct.vertex_position[0][0] << " " << tetrahedron[1].mesh_struct.vertex_position[0][1] << " " << tetrahedron[1].mesh_struct.vertex_position[0][2] << std::endl;
	*/

	for (unsigned int j = 0; j < cloth.size(); ++j) {
		cloth[j].ori_vertices= cloth[j].mesh_struct.vertex_position;
	}
	for (unsigned int j = 0; j < tetrahedron.size(); ++j) {
		tetrahedron[j].ori_vertices= tetrahedron[j].mesh_struct.vertex_position;
	}
	for (unsigned int j = 0; j < collider.size(); ++j) {
		collider[j].ori_vertices = collider[j].mesh_struct.vertex_position;
	}
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].initial();
	}
	for (int i = 0; i < tetrahedron.size(); ++i) {
		tetrahedron[i].initial();
	}
	//intersection.initialIntersection();
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].reset();
	}
}

void Scene::setCursorForce(Camera* camera, double* cursor_screen, float force_coe)
{
	double cursor_pos[3];
	double force_direction[3];

	if (intersection.obj_No < cloth.size()) {
		getCursorPos(cursor_pos, cloth[intersection.obj_No].mesh_struct.vertex_for_render,
			cloth[intersection.obj_No].mesh_struct.triangle_indices[intersection.face_index].data());
		cloth[intersection.obj_No].findAllNeighborVertex(intersection.face_index, cursor_pos, ave_edge_length);
	}
	else{
		getCursorPos(cursor_pos, tetrahedron[intersection.obj_No - cloth.size()].mesh_struct.vertex_for_render,
			tetrahedron[intersection.obj_No - cloth.size()].mesh_struct.triangle_indices[intersection.face_index].data());
		tetrahedron[intersection.obj_No-cloth.size()].findAllNeighborVertex(intersection.face_index, cursor_pos, ave_edge_length);
	}
	double cursor_pos_in_space[3];
	cursorMovement(camera, cursor_screen, force_direction, force_coe, cursor_pos, cursor_pos_in_space);
	cursor.translate(cursor_pos, cursor_pos_in_space);
	if (intersection.obj_No < cloth.size()) {
		project_dynamic.addExternalClothForce(force_direction, cloth[intersection.obj_No].coe_neighbor_vertex_force, cloth[intersection.obj_No].neighbor_vertex, intersection.obj_No);
	}
	else {
		//std::cout << force_direction[0] << " " << force_direction[1] << " " << force_direction[2] << std::endl;
		project_dynamic.addExternalTetForce(force_direction, tetrahedron[intersection.obj_No].coe_neighbor_vertex_force, tetrahedron[intersection.obj_No].neighbor_vertex, intersection.obj_No);
	}

}



void Scene::cursorMovement(Camera* camera, double* cursor_screen, double* force_direction, float force_coe, double* object_position, double* cursor_pos_in_space)
{
	camera->getCursorPosInSpace(cursor_pos_in_space, cursor_screen, object_position);
	SUB(force_direction, cursor_pos_in_space, object_position);
	MULTI(force_direction, force_direction, force_coe);
	double force_magnitude = sqrt(DOT(force_direction, force_direction));
	if (force_magnitude > max_force_magnitude) {
		force_magnitude = max_force_magnitude / force_magnitude;
		MULTI(force_direction, force_direction, force_magnitude);
	}
	//std::cout << cursor_pos_in_space[0] << " " << cursor_pos_in_space[1] << " " << cursor_pos_in_space[2] << std::endl;
	//std::cout << object_position[0] << " " << object_position[1] << " " << object_position[2] << std::endl;

	//std::cout << "force mag " << force_direction[0] << " " << force_direction[1] << " " << force_direction[2] << std::endl;
}

void Scene::genShader()
{
	object_shader_front = new Shader("./shader/object_triangle.vs", "./shader/object_triangle.fs");
	*object_shader_front = Shader("./shader/object_triangle.vs", "./shader/object_triangle.fs");
	wireframe_shader = new Shader("./shader/wireframe.vs", "./shader/wireframe.fs");
	*wireframe_shader = Shader("./shader/wireframe.vs", "./shader/wireframe.fs");
}

////TEST_ARRAY
//void Scene::testForWritetToArray(int thread_No)
//{
//	unsigned int* test_pair_ = test_pair[thread_No];
//	for (unsigned int j = 0; j < 100000; ++j) {
//		for (unsigned int i = 0; i < 10000; ++i) {
//			test_pair_[i] = i + j;
//		}
//	}
//}

//TEST_ARRAY
void Scene::testForWritetToArray(int thread_No)
{
	unsigned int array_size = test_array_size / 2;
	unsigned int* test_pair_ = test_pair[thread_No];
	unsigned int test = 0;
	unsigned int kk = 0;
	for (unsigned int i = 0; i < array_size; i++) {
		//kk[0] = i;
		//kk[1] = i + 1;
		//kk[2] = i + 2;
		//kk[3] = i + 3;
		//test_pair_[i] = i;
		//_mm_stream_ps((float*)&test_pair_[i], _mm_loadu_ps((float*)kk));
		kk = i * 2;
		_mm_stream_si32((int*)test_pair_++, *(int*)&kk);

		//test += test_pair_[i];
		//test += i;
	}
	//thread_test[10000 * thread_No] = test;
}

//for (unsigned int m = 0; m < 100; ++m) {
//	test_pair_[i] -= 2;
//	test_pair_[i] *= 2;
//}
//std::cout << thread_test[10000 * thread_No] << std::endl;


void Scene::testForWritetToArraySingle(int total_thread_num)
{
	unsigned int array_size = test_array_size / 2;
	unsigned int* test_pair_;
	unsigned int test = 0;
	unsigned int kk = 0;
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		test_pair_ = test_pair_single[i];
		for (unsigned int j = 0; j < array_size; j++) {
			kk = j * 2;
			//kk[0] = j;
			//kk[1] = j+1;
			//kk[2] = j+2;
			//kk[3] = j+3;
			//_mm_stream_ps((float*)&test_pair_[j], _mm_loadu_ps((float*)kk));
			_mm_stream_si32((int*)&test_pair_[j], *(int*)&kk);

		}
		//thread_test[10000 * i] = test;
	}

	//std::cout << test_pair_single[0][4] << " " << test_pair_single[0][5] << " " << test_pair_single[0][6] << " " << test_pair_single[0][7] << std::endl;
}


//for (unsigned int m = 0; m < 100; ++m) {
//	test_pair_[j] -= 2;
//	test_pair_[j] *= 2;
//}

//thread_test[10000 * k] = test;
//std::cout << thread_test[10000 * k] << std::endl;

void Scene::voidForWritetToArraySingle(int total_thread_num)
{
	unsigned int* test_pair_;
	unsigned int* test_pair_2_;
	for (unsigned int k = 0; k < total_thread_num; ++k) {
		test_pair_ = test_pair_single[k];
		test_pair_2_ = test_pair[k];
		for (unsigned int i = 0; i < test_array_size; ++i) {
			test_pair_[i] = 0;
			test_pair_2_[i] = 0;
		}
	}


}