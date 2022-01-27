#include"scene.h"
#include<bitset>

Scene::Scene()
{
	light.ambient = glm::vec3(1.0, 1.0, 1.0);
	light.diffuse = glm::vec3(0.8, 0.8, 0.8);
	light.specular = glm::vec3(0.95, 0.95, 0.95);

	
	time_step = 1.0 / 50.0;
	project_dynamic.time_step = time_step;

	max_force_magnitude = 2.0;

	last_output_obj_stamp = -1;
	time_stamp = 0;
	genShader();	
	project_dynamic.collision.time_stamp = &time_stamp;
	project_dynamic.time_stamp = &time_stamp;
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
	std::vector<int>cloth_index_in_object, tetrahedron_index_in_object;
	for (int i = 0; i < object_path.size(); ++i) {
		if (preprocessing.ori_simulation_mesh[i].type == TRIANGLE) {
			cloth_index_in_object.push_back(i);
		}
		else {
			tetrahedron_index_in_object.push_back(i);
		}
	}
	cloth_num = cloth_index_in_object.size();
	tetrahedron_num = tetrahedron_index_in_object.size();
	cloth.resize(cloth_num);
	tetrahedron.resize(tetrahedron_num);
	double cloth_density=15.0;
	double tetrahedron_density = 10.0;
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
	
	std::array<double, 4>tetrahedron_collision_stiffness_per = { 2e5,2e5,2e5, 2e5 };
	double sigma_limit[2] = { 0.95,1.05 };
	SingleTetrahedronInfo single_tetrahedron_info(tetrahedron_density,1e5,1e3,1e3, tetrahedron_collision_stiffness_per.data(), sigma_limit);
	for (int i = 0; i < tetrahedron_num; ++i) {
		tetrahedron[i].recordInitialMesh(single_tetrahedron_info);
	}


	project_dynamic.setForPD(&cloth, &tetrahedron,&collider, &thread);
	setAveEdgeLength();
	cursor.createVertices(4.0 * ave_edge_length, camera_center);
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
	color1 = new double* [cloth_num+ tetrahedron_num];
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

void Scene::drawScene(Camera* camera, std::vector<std::vector<bool>>& wireframe, std::vector<std::vector<bool>>& hide, bool start_save_obj)
{
	shadow.drawShadow(camera, hide,cloth,collider,tetrahedron);
	glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_CUBE_MAP, shadow.depth_map);
	glEnable(GL_CULL_FACE);
	for (int j = 0; j < cloth_num; ++j) {
		if (!hide[CLOTH_][j]) {
			cloth[j].setSceneShader(light, camera, shadow.far_plane, object_shader_front);
			cloth[j].draw(camera,object_shader_front);
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
			collider[i].drawWireframe(camera,wireframe_shader);
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
		for (int i = 0; i < cloth.size(); ++i) {
			save_obj.write(cloth[i].mesh_struct.vertex_for_render, cloth[i].mesh_struct.triangle_indices, 10, time_stamp,i);
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

void Scene::updateCloth(Camera* camera, double* cursor_screen, bool* control_parameter, float force_coe, bool& record_matrix,
	double& ave_iteration)
{
	if (control_parameter[START_SIMULATION] || control_parameter[ONE_FRAME]) {
		project_dynamic.resetExternalForce();
		//std::cout << intersection.happened << " " << control_parameter[START_TEST] << std::endl;
		if (intersection.happened && !control_parameter[START_TEST]) {
			//std::cout << cursor_screen[0] << " " << cursor_screen[1] << " " << force_coe << std::endl;
			setCursorForce(camera, cursor_screen, force_coe);
		}
		project_dynamic.PDsolve();
		//project_dynamic.PD_IPC_solve(record_matrix);
		project_dynamic.update_ave_iteration_record(ave_iteration);
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
	for (int i = 0; i < tetrahedron.size(); ++i) {
		tetrahedron[i].setBuffer();
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
	
}


void Scene::selectAnchor(bool* control_parameter, bool* select_anchor, double* screen_pos, bool press_state, bool pre_press_state, Camera* camera, 
	std::vector<bool>& hide)
{	
	if (select_anchor[0]) {
		control_parameter[START_SIMULATION] = false;
		if (press_state) {
			set_tetrahedron_anchor.setCorner(screen_pos, pre_press_state,tetrahedron,camera,hide);
		}
	}
	if (select_anchor[1]) {
		for (int i = 0; i < tetrahedron.size(); ++i) {
			tetrahedron[i].mesh_struct.setAnchorPosition();

		}		
		select_anchor[1] = false;
		project_dynamic.updateTetrohedronAnchorVertices();
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
	mouse_pos[1] =SCR_HEIGHT-pos[1];
	int chosen_index[2];
	bool is_cloth;
	pick_triangle.pickTriangle(&cloth, &collider, &tetrahedron, camera, hide, chosen_index, is_cloth, mouse_pos);
	intersection.initialIntersection();
	//std::cout << chosen_index[0] << std::endl;
	if (chosen_index[0] > -1) {
		double cursor_pos[3];
		intersection.setIntersection(chosen_index, is_cloth);
		if (is_cloth) {
			getCursorPos(cursor_pos, cloth[intersection.obj_No].mesh_struct.vertex_for_render,
				cloth[intersection.obj_No].mesh_struct.triangle_indices[intersection.face_index].data());
			cloth[chosen_index[1]].findAllNeighborVertex(chosen_index[0], cursor_pos, ave_edge_length);
		}
		else {
			getCursorPos(cursor_pos, tetrahedron[intersection.obj_No].mesh_struct.vertex_for_render,
				tetrahedron[intersection.obj_No].mesh_struct.triangle_indices[intersection.face_index].data());
		}
	
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

	if (intersection.is_cloth) {
		getCursorPos(cursor_pos, cloth[intersection.obj_No].mesh_struct.vertex_for_render,
			cloth[intersection.obj_No].mesh_struct.triangle_indices[intersection.face_index].data());
	}
	else {
		getCursorPos(cursor_pos, tetrahedron[intersection.obj_No].mesh_struct.vertex_for_render,
			tetrahedron[intersection.obj_No].mesh_struct.triangle_indices[intersection.face_index].data());
	}
	
	cursor.translate(cursor_pos);
	
	cursorMovement(camera, cursor_screen, force_direction, force_coe, cursor_pos);
	if (intersection.is_cloth) {
		project_dynamic.addExternalClothForce(force_direction, cloth[intersection.obj_No].coe_neighbor_vertex_force, cloth[intersection.obj_No].neighbor_vertex, intersection.obj_No);
	}
	else {
		project_dynamic.addExternalTetForce(force_direction, tetrahedron[intersection.obj_No].coe_neighbor_vertex_force, tetrahedron[intersection.obj_No].neighbor_vertex, intersection.obj_No);
	}
	
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