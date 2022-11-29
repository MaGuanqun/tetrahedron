#include"scene.h"
#include"basic/save_simulate_parameter.h"
#include<bitset>
#include<algorithm>

Scene::Scene()
{
	time_record_interval = 10;

	light.ambient = glm::vec3(1.0, 1.0, 1.0);
	light.diffuse = glm::vec3(0.9, 0.9, 0.9);
	light.specular = glm::vec3(0.85, 0.85, 0.85);

	time_step = 1.0 / 100.0;


	max_force_magnitude = 200.0;

	last_output_obj_stamp = -1;
	time_stamp = 0;
	time_indicate_for_simu = 0;

	xpbd.time_indicate_for_simu = &time_indicate_for_simu;
	newton_method.time_indicate_for_simu = &time_indicate_for_simu;
	xpbd.move_model = &move_model;
	xpbd_ipc.time_indicate_for_simu = &time_indicate_for_simu;
	xpbd_ipc.move_model = &move_model;
	
	genShader();


	project_dynamic.time_step = time_step;
	project_dynamic.collision.time_stamp = &time_stamp;
	project_dynamic.time_stamp = &time_stamp;
	xpbd.time_step = time_step;
	xpbd.time_stamp= &time_stamp;
	xpbd_ipc.time_step = time_step;
	xpbd_ipc.time_stamp = &time_stamp;

	newton_method.time_step = time_step;
	newton_method.time_step_square = time_step * time_step;
	newton_method.time_stamp = &time_stamp;

	second_order_xpbd_large.time_step = time_step;
	second_order_xpbd_large.time_stamp = &time_stamp;

	select_dimension_index = 4;
	intersect_when_rotation = false;
	start_rotation = false;
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
	xpbd.has_force = &read_force;
	xpbd_ipc.has_force = &read_force;
	newton_method.has_force = &read_force;
}


void Scene::reorganizeData()
{
	anchor_vertex.resize(cloth.size() + tetrahedron.size());
	for (unsigned int i = 0; i < cloth.size(); ++i) {
		anchor_vertex[i] = &cloth[i].mesh_struct.anchor_vertex;
	}
	for (unsigned int i = 0; i < tetrahedron.size(); ++i) {
		anchor_vertex[i] = &tetrahedron[i].mesh_struct.anchor_vertex;
	}

	mesh_struct.resize(cloth.size() + tetrahedron.size());
	for (unsigned int i = 0; i < cloth.size(); ++i) {
		mesh_struct[i] = &cloth[i].mesh_struct;
	}
	for (unsigned int i = 0; i < tetrahedron.size(); ++i) {
		mesh_struct[i] = &tetrahedron[i].mesh_struct;
	}
}


void Scene::saveParameter(std::vector<std::string>& path, std::vector<std::string>& collider_path, std::vector<std::array<double, 6>>& cloth_stiffness, std::vector<std::array<double, 6>>& tet_stiffness,
	std::vector<std::array<double, 8>>& cloth_collision_stiffness, std::vector<std::array<double, 8>>& tet_collision_stiffness, double* tolerance_ratio, double* friction_coe)
{
	double velocity_damp=1.0;
	unsigned int sub_step_per_detection = 1;
	double d_hat = 0.0;

	double local_conv_rate = 1e-2;


	unsigned int min_inner_itr = 1;
	unsigned int min_outer_itr = 1;
	double displacement_standard = 1e-4;

	double distance_record_as_collision = 1e-6;


	unsigned int max_outer_iteration_num=10;
	unsigned int max_inner_iteration_num = 10;

	double energy_converge_standard = 1e-7;

	switch (use_method)
	{
	case PD_:
		velocity_damp = project_dynamic.velocity_damp;
		local_conv_rate = project_dynamic.local_global_conv_rate;
		break;
	case XPBD_:
		velocity_damp = xpbd.velocity_damp;
		sub_step_per_detection = *xpbd.sub_step_per_detection;
		max_outer_iteration_num = xpbd.max_iteration_number;
		break;
	case NEWTON_:
		velocity_damp = 1.0;
		break;
	case XPBD_SECOND_ORDER_LARGE_:
		velocity_damp = second_order_xpbd_large.velocity_damp;
		break;
	case XPBD_IPC_:
		velocity_damp = xpbd_ipc.velocity_damp;
		sub_step_per_detection = *xpbd_ipc.sub_step_per_detection;
		d_hat = xpbd_ipc.collision.d_hat;
		local_conv_rate = xpbd_ipc.energy_converge_ratio;
		min_inner_itr = xpbd_ipc.min_inner_iteration;
		min_outer_itr = xpbd_ipc.min_outer_iteration;
		displacement_standard = xpbd_ipc.max_move_standard;
		distance_record_as_collision = xpbd_ipc.collision.tolerance;
		max_outer_iteration_num = xpbd_ipc.outer_max_iteration_number;
		max_inner_iteration_num = xpbd_ipc.max_iteration_number;
		energy_converge_standard = xpbd_ipc.energy_converge_standard;

		break;
	}

	SaveParameter::writeParameter(path, collider_path, cloth_stiffness, tet_stiffness, cloth_collision_stiffness, tet_collision_stiffness, use_method,	
		anchor_vertex, time_step,project_dynamic.outer_itr_conv_rate, local_conv_rate,
		xpbd.sub_step_num, max_outer_iteration_num, cloth_density,tetrahedron_density,velocity_damp,
		friction_coe, sub_step_per_detection,floor.exist,floor.dimension,floor.normal_direction,floor.value, d_hat,min_inner_itr,min_outer_itr,
		displacement_standard, distance_record_as_collision, camera->position, camera->up,camera->center, max_inner_iteration_num, energy_converge_standard);



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

void Scene::setFloorInfo(bool exist, bool show, bool normal_direction, unsigned int dimension, double value, bool& eidit,bool* control_parameter)
{
	if (eidit) {
		control_parameter[START_SIMULATION] = false;
		eidit = false;
	}
	floor.show = show;
	floor.exist = exist;
	floor.setFloor(dimension, value, normal_direction);
}

void Scene::updateItrInfo(int* iteration_num)
{
	//if (!only_test_collision) {
		switch (use_method)
		{
		case XPBD_:
			xpbd.updateItrInfo(iteration_num);
			break;
		case XPBD_IPC_:
			xpbd_ipc.updateItrInfo(iteration_num);
			break;
		}
	//}
}

void Scene::obtainConvergenceInfo(double* convergence_rate, int* iteration_num)
{
	//if (!only_test_collision) {
		switch (use_method)
		{
		case PD_:
			convergence_rate[LOCAL_GLOBAL] = project_dynamic.local_global_conv_rate;
			convergence_rate[OUTER] = project_dynamic.outer_itr_conv_rate;
			iteration_num[LOCAL_GLOBAL] = project_dynamic.local_global_iteration_num;
			iteration_num[OUTER] = project_dynamic.outer_iteration_num;
			break;
		case XPBD_:
			iteration_num[LOCAL_GLOBAL] = xpbd.iteration_number;
			iteration_num[OUTER] = xpbd.sub_step_num;
			break;
		case NEWTON_:
			iteration_num[LOCAL_GLOBAL] = newton_method.iteration_number;
			iteration_num[OUTER] = newton_method.iteration_number;
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			iteration_num[LOCAL_GLOBAL] = second_order_xpbd_large.iteration_number;
			iteration_num[OUTER] = second_order_xpbd_large.iteration_number;
			break;
		case XPBD_IPC_:
			iteration_num[LOCAL_GLOBAL] = xpbd_ipc.iteration_number;
			iteration_num[OUTER] = xpbd_ipc.outer_itr_num;
			break;
		}
	//}
}

void Scene::updateConvRate(double* convergence_rate)
{
	// if (!only_test_collision) {
		switch (use_method)
		{
		case PD_:
			project_dynamic.outer_itr_conv_rate = convergence_rate[OUTER];
			project_dynamic.local_global_conv_rate = convergence_rate[LOCAL_GLOBAL];
			break;
		case NEWTON_:
			newton_method.conv_rate = time_step * convergence_rate[OUTER];
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			second_order_xpbd_large.conv_rate = time_step * convergence_rate[OUTER];
			break;
		}
	// }
}

bool Scene::loadMesh(std::string& scene_path, std::vector<std::string>& collider_path, std::vector<std::string>& object_path, double* tolerance_ratio, bool* control_parameter,
	double* initial_stiffness, double* friction_coe, unsigned int* sub_step_per_detection, bool* floor_indicator, unsigned int& floor_dimension, double& floor_value)
{

	use_method =100;
	this->control_parameter = control_parameter;
	xpbd.control_parameter = control_parameter;
	xpbd_ipc.control_parameter = control_parameter;
	if (control_parameter[USE_XPBD]) {
		use_method = XPBD_;
	}
	else if (control_parameter[USE_PD_]) {
		use_method = PD_;
	}
	else if(control_parameter[USE_NEWTON_]) {
		use_method = NEWTON_;
	}
	else if (control_parameter[USE_XPBD_LARGE]) {
		use_method = XPBD_SECOND_ORDER_LARGE_;
	}
	else if (control_parameter[USE_XPBD_IPC]) {
		use_method = XPBD_IPC_;
	}
	only_test_collision = control_parameter[ONLY_COLLISION_TEST];
	cloth_density = 1;
	tetrahedron_density = 1;


	std::vector<std::vector<double>>obj_stiffness,collide_stiffness;
	std::vector<std::vector<int>>anchor_vertex;
	bool load_by_scene_file=false;
	if (!scene_path.empty()) {
		load_by_scene_file = true;
	}

	double d_hat=5e-3;

	double velocity_damp = 1.0;
	double local_global_conv_rate = 1e-2;
	unsigned int min_inner_itr = 1; 
	unsigned int min_outer_itr = 1 ;
	double displacement_standard = 1e-4;
	double distance_record_as_collision = 1e-6;
	unsigned int max_outer_itr_num=10;

	unsigned int max_inner_itr_num = 10;
	double energy_converge_standard = 1e-7;

	if (load_by_scene_file) {
		SaveParameter::readFile(scene_path, object_path, collider_path, obj_stiffness, collide_stiffness, anchor_vertex, time_step, use_method, xpbd.sub_step_num, max_outer_itr_num,
			local_global_conv_rate, project_dynamic.outer_itr_conv_rate,cloth_density,tetrahedron_density, velocity_damp,friction_coe,
			*sub_step_per_detection, floor_indicator[0], floor_dimension, floor_indicator[2], floor_value,d_hat, min_inner_itr, min_outer_itr,
			displacement_standard, distance_record_as_collision, camera->ori_position,camera->ori_up,camera->ori_center, max_inner_itr_num,
			energy_converge_standard	);

		camera->position = camera->ori_position;
		camera->up = camera->ori_up;
		camera->center = camera->ori_center;

		control_parameter[USE_XPBD] = false;
		control_parameter[USE_PD_] = false;
		control_parameter[USE_NEWTON_] = false;
		control_parameter[USE_XPBD_LARGE] = false;
		control_parameter[USE_XPBD_IPC] = false;
		switch (use_method)
		{
		case XPBD_:
			control_parameter[USE_XPBD] = true;
			xpbd.max_iteration_number = max_outer_itr_num;
			break;
		case PD_:
			control_parameter[USE_PD_] = true;
			project_dynamic.local_global_conv_rate = local_global_conv_rate;
			break;
		case NEWTON_:
			control_parameter[USE_NEWTON_] = true;
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			control_parameter[USE_XPBD_LARGE] = true;
			break;
		case XPBD_IPC_:
			control_parameter[USE_XPBD_IPC] = true;
			xpbd_ipc.collision.d_hat = d_hat;
			xpbd_ipc.energy_converge_ratio = local_global_conv_rate;
			xpbd_ipc.min_inner_iteration = min_inner_itr;
			xpbd_ipc.min_outer_iteration = min_outer_itr;
			xpbd_ipc.max_move_standard = displacement_standard;
			xpbd_ipc.collision.tolerance = distance_record_as_collision;
			xpbd_ipc.outer_max_iteration_number = max_outer_itr_num;
			xpbd_ipc.max_iteration_number = max_inner_itr_num;
			xpbd_ipc.energy_converge_standard = energy_converge_standard;

			std::cout << "d_hat " << d_hat << std::endl;
			std::cout << "energy_converge_ratio " << local_global_conv_rate << std::endl;
			std::cout << "min_inner_iteration " << min_inner_itr << std::endl;
			std::cout << "min_outer_itr " << min_outer_itr << std::endl;
			std::cout << "max_move_standard " << displacement_standard << std::endl;
			std::cout << "distance_record_as_collision " << distance_record_as_collision << std::endl;
		
			std::cout << "camera position " << camera->position.x << " " << camera->position.y << " " << camera->position.z << std::endl;
			std::cout << "camera up " << camera->up.x << " " << camera->up.y << " " << camera->up.z << std::endl;
			std::cout << "camera center " << camera->center.x << " " << camera->center.y << " " << camera->center.z << std::endl;

			std::cout << "outer_max_iteration_number " << xpbd_ipc.outer_max_iteration_number << std::endl;
			std::cout << "max_inner_itr_num " << xpbd_ipc.max_iteration_number << std::endl;
			std::cout << "energy_converge_standard  " << xpbd_ipc.energy_converge_standard << std::endl;

			break;
		}


		if (floor_indicator[0]) {
			floor_indicator[1] = true;
		}
	}


	if (use_method== XPBD_) {
		xpbd.collision.friction_coe = friction_coe;
		xpbd.sub_step_per_detection = sub_step_per_detection;
	}
	else if (use_method == PD_) {
		project_dynamic.collision.friction_coe = friction_coe;
	}
	else if (use_method == XPBD_IPC_) {
		xpbd_ipc.collision.friction_coe = friction_coe;
		xpbd_ipc.sub_step_per_detection = sub_step_per_detection;
	}

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

	for (int i = 0; i < cloth_num; ++i) {
		cloth[i].loadMesh(preprocessing.ori_simulation_mesh[cloth_index_in_object[i]], cloth_density, &thread);
		cloth[i].mesh_struct.initialUnfixedIndex();
	}
	for (int i = 0; i < tetrahedron_num; ++i) {
		tetrahedron[i].loadMesh(preprocessing.ori_simulation_mesh[tetrahedron_index_in_object[i]], tetrahedron_density, &thread);
		tetrahedron[i].mesh_struct.initialUnfixedIndex();
		
	}
	setWireframwColor();


	double cloth_position_stiffness = 1e4;
	double tet_position_stiffness = 1e4;

	if (load_by_scene_file) {
		switch (use_method)
		{
		case XPBD_:
			xpbd.velocity_damp = velocity_damp;
			break;
		case PD_:
			project_dynamic.velocity_damp = velocity_damp;
			break;
		case XPBD_IPC_:
			xpbd_ipc.velocity_damp = velocity_damp;
			break;
		}
		if (!cloth.empty()) {
			initial_stiffness[LENGTH] = obj_stiffness[0][0];
			initial_stiffness[BENDING]= obj_stiffness[0][1];
			cloth_position_stiffness = obj_stiffness[0][2];
			initial_stiffness[DAMP_LENGTH] = obj_stiffness[0][3];
			initial_stiffness[DAMP_BENDING] = obj_stiffness[0][4];
			memcpy(initial_stiffness, collide_stiffness[0].data(), 32);
			memcpy(initial_stiffness+ DAMP_BODY_POINT_TRIANGLE, collide_stiffness[0].data()+4, 32);			
			for (unsigned int i = 0; i < cloth.size(); ++i) {
				if (!anchor_vertex[i].empty()) {
					cloth[i].mesh_struct.anchor_vertex = anchor_vertex[i];
					cloth[i].mesh_struct.setAnchorPosition();
					cloth[i].mesh_struct.resetMassInv();
				}
				else {
					cloth[i].mesh_struct.anchor_vertex.clear();
					cloth[i].mesh_struct.anchor_position.clear();
					cloth[i].mesh_struct.resetMassInv();
				}
			}
		}
		if (!tetrahedron.empty()) {
			initial_stiffness[ARAP]= obj_stiffness[cloth.size()][0];
			tet_position_stiffness = obj_stiffness[cloth.size()][1];
			initial_stiffness[TET_EDGE_LENGTH] = obj_stiffness[cloth.size()][2];
			initial_stiffness[DAMP_ARAP] = obj_stiffness[cloth.size()][3];
			memcpy(initial_stiffness, collide_stiffness[cloth.size()].data(), 32);
			memcpy(initial_stiffness + DAMP_BODY_POINT_TRIANGLE, collide_stiffness[0].data() + 4, 32);
			for (unsigned int i = 0; i < tetrahedron.size(); ++i) {
				if (!anchor_vertex[i+cloth.size()].empty()) {
					tetrahedron[i].mesh_struct.anchor_vertex = anchor_vertex[i + cloth.size()];
					tetrahedron[i].mesh_struct.setAnchorPosition();
					tetrahedron[i].mesh_struct.resetMassInv();
				}
				else {
					tetrahedron[i].mesh_struct.anchor_vertex.clear();
					tetrahedron[i].mesh_struct.anchor_position.clear();
					tetrahedron[i].mesh_struct.resetMassInv();
				}
			}
		}
	}


	std::vector<SingleClothInfo> single_cloth_info;
	//std::array<double, 4>collision_stiffness_per = { 2e5,2e5,2e5, 2e5 };// stiffness of collision constraint //=0 body point triangle, =1 point-triangle =2 edge-edge =3 point-point,
	//std::vector<std::array<double, 4>>collision_stiffness(cloth_num, collision_stiffness_per);
	for (int i = 0; i < cloth_num; ++i) {
		single_cloth_info.push_back(SingleClothInfo(cloth_density, initial_stiffness[LENGTH], cloth_position_stiffness, initial_stiffness[BENDING], initial_stiffness, 0.5, 0.4, initial_stiffness[LENGTH],
			initial_stiffness[DAMP_LENGTH], initial_stiffness[DAMP_BENDING]));
	}
	for (int i = 0; i < cloth_num; ++i) {
		cloth[i].recordInitialMesh(single_cloth_info[i]);
	}
	//std::array<double, 4>tetrahedron_collision_stiffness_per = {1e1,1e1, 1e1,1e1 };
	double sigma_limit[2] = { 0.99,1.01 };
	SingleTetrahedronInfo single_tetrahedron_info(tetrahedron_density, 5e4, initial_stiffness[ARAP], 0.0, initial_stiffness, sigma_limit,
		5e4,0.45, initial_stiffness[TET_EDGE_LENGTH],initial_stiffness[DAMP_ARAP], 0.0);
	for (int i = 0; i < tetrahedron_num; ++i) {
		tetrahedron[i].recordInitialMesh(single_tetrahedron_info);
	}

	if (use_method == XPBD_ || use_method == XPBD_IPC_) {
		setGroup();
	}


	move_object.initial(&cloth, &collider, &tetrahedron, &thread);

	switch (use_method)
	{
	case PD_:
		project_dynamic.setForPD(&cloth, &tetrahedron, &collider, &floor, &thread, tolerance_ratio);
		break;
	case XPBD_:
		xpbd.setForXPBD(&cloth, &tetrahedron, &collider, &floor, &thread, tolerance_ratio);
		break;
	case NEWTON_:
		newton_method.setForNewtonMethod(&cloth, &tetrahedron, &collider, &floor, &thread, tolerance_ratio);
		break;
	case XPBD_SECOND_ORDER_LARGE_:
		second_order_xpbd_large.setForNewtonMethod(&cloth, &tetrahedron, &collider, &floor, &thread, tolerance_ratio);
		break;
	case XPBD_IPC_:
		xpbd_ipc.setForXPBD(&cloth, &tetrahedron, &collider, &floor, &thread, tolerance_ratio);
		break;
	}

	if (load_by_scene_file) {
		for (unsigned int i = 0; i < tetrahedron.size(); ++i) {
			if (!anchor_vertex[i + cloth.size()].empty()) {
				tetrahedron[i].mesh_struct.updateAnchorPerThread(thread.thread_num);
				tetrahedron[i].mesh_struct.updateUnfixedPointData();
				if (use_method == NEWTON_) {
					newton_method.updateIndexBeginPerObj();
				}
			}
			else {
				tetrahedron[i].mesh_struct.updateAnchorPerThread(thread.thread_num);
				tetrahedron[i].mesh_struct.updateUnfixedPointData();
				if (use_method == NEWTON_) {
					newton_method.updateIndexBeginPerObj();
				}
			}
		}
	}


	if (control_parameter[ONLY_COLLISION_TEST]) {
		test_draw_collision.initial(&cloth, &collider, &tetrahedron, &thread, &floor, tolerance_ratio, &project_dynamic.collision,
			&xpbd.collision, &newton_method.collision,&xpbd_ipc.collision, use_method);
	}

	setAveEdgeLength();
	cursor.createVertices(0.03, camera_center);

	reorganizeData();


	if (load_by_scene_file) {
		updateAnchorTet();
	}

	//unsigned int  k = cloth[0].mesh_struct.vertices[49 * 50].edge[0];
	//unsigned int v;
	//if (cloth[0].mesh_struct.edge_vertices[k << 1] == 49 * 50) {
	//	v = cloth[0].mesh_struct.edge_vertices[(k << 1) + 1];
	//}
	//else {
	//	v = cloth'[0].mesh_struct.edge_vertices[(k << 1)];
	//}
	// 

	//reflectModel();

	return load_by_scene_file;
}

void Scene::reflectModel()
{
	//for (unsigned int i = 0; i < tetrahedron[0].mesh_struct.vertex_position.size(); ++i) {
	//	tetrahedron[0].mesh_struct.vertex_position[i][1] = -0.32 - tetrahedron[0].mesh_struct.vertex_position[i][1];
	//	tetrahedron[0].mesh_struct.vertex_for_render[i][1] = -0.32 - tetrahedron[0].mesh_struct.vertex_for_render[i][1];
	//}


	std::uniform_real_distribution<double> distribution(-2, 2);
	std::default_random_engine engine;
	std::vector<double>x(tetrahedron[0].mesh_struct.vertex_position.size() * 3);
	auto gen = [&distribution, &engine]() {
		return distribution(engine);
	};
	std::generate(std::begin(x), std::end(x),gen);
	memcpy(tetrahedron[0].mesh_struct.vertex_position[0].data(), x.data(),8 * x.size());
	memcpy(tetrahedron[0].mesh_struct.vertex_for_render[0].data(), x.data(),8 * x.size());


}

void Scene::checkVolume()
{
	for (int i = 0; i < tetrahedron[0].mesh_struct.indices.size(); ++i) {
		double volume = getTetrahedronVolume(tetrahedron[0].mesh_struct.vertex_position[tetrahedron[0].mesh_struct.indices[i][0]].data(), tetrahedron[0].mesh_struct.vertex_position[tetrahedron[0].mesh_struct.indices[i][1]].data(), tetrahedron[0].mesh_struct.vertex_position[tetrahedron[0].mesh_struct.indices[i][2]].data(),
			tetrahedron[0].mesh_struct.vertex_position[tetrahedron[0].mesh_struct.indices[i][3]].data());
		if (volume < 0) {
			std::cout << "error volume negative " << std::endl;
			std::cout << volume << std::endl;
		}
	}
}

void Scene::setGroup()
{
	//time_t t = clock();
	//for (unsigned int k = 0; k < 10; ++k) {
		for (int i = 0; i < cloth_num; ++i) {
			for (unsigned int j = 0; j < cloth[i].mesh_struct.unconnected_edge_index.size(); j++) {
				cloth[i].mesh_struct.unconnected_edge_index[j].clear();
			}
			for (unsigned int j = 0; j < cloth[i].mesh_struct.unconnected_vertex_index.size(); j++) {
				cloth[i].mesh_struct.unconnected_vertex_index.clear();
			}
			if (cloth[i].mesh_struct.edge_length.size() > 3) {
				graph_color.graphColorEdgeLength(cloth[i].mesh_struct);
				if (cloth[i].mesh_struct.vertex_position.size() > 2) {
					graph_color.graphColorBending(cloth[i].mesh_struct);
				}
			}
			else {
				cloth[i].mesh_struct.unconnected_edge_index.resize(cloth[i].mesh_struct.edge_length.size());
				for (unsigned int j = 0; j < cloth[i].mesh_struct.unconnected_edge_index.size(); ++j) {
					cloth[i].mesh_struct.unconnected_edge_index[j].emplace_back(j);
				}
			}			
			//graph_color.testBend(cloth[i].mesh_struct.unconnected_vertex_index, cloth[i].mesh_struct);
			//graph_color.testEdge(cloth[i].mesh_struct.unconnected_edge_index, cloth[i].mesh_struct,
			//	cloth[i].mesh_struct.edge_vertices);
		}

		for (int i = 0; i < tetrahedron.size(); ++i) {
			if (tetrahedron[i].mesh_struct.indices.size() > 1) {
				graph_color.graphColor(tetrahedron[i].mesh_struct.tet_tet_index, tetrahedron[i].mesh_struct.unconnected_tet_index);
				tetrahedron[i].mesh_struct.obtainVETofColors();
				tetrahedron[i].mesh_struct.setTetColorStartPerThread();
			}
			else {
				tetrahedron[i].mesh_struct.unconnected_tet_index.resize(1);
				tetrahedron[i].mesh_struct.unconnected_tet_index[0].emplace_back(0);
				tetrahedron[i].mesh_struct.obtainVETofColors();
				tetrahedron[i].mesh_struct.setTetColorStartPerThread();
			}
		}

}

void Scene::setWireframwColor()
{
	double color3[4][3] = { {1.0,0.0,0.0},{0.0,0.0,1.0},{0.0,0.3,0.0},{0.0, 1.0,1.0} };
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

	for (int i = 0; i < cloth_num + tetrahedron_num; ++i)
	{
		delete[] color1[i];
	}
}

void Scene::initialSceneSetting(Preprocessing& preprocessing)
{
	shadow.camera_from_origin = 0.8 * preprocessing.scene_info.max_dis_from_center;//0.8
	shadow.far_plane = 10.0 * preprocessing.scene_info.max_dis_from_center;
	memcpy(camera_center, preprocessing.scene_info.camera_center, 24);
	shadow.setBasic();

}


void Scene::getClothInfo(std::vector<std::array<int, 3>>& mesh_info, std::vector<double>& mass, std::vector<std::array<double, 6>>& mesh_stiffness, double* simulation_parameter, std::vector<std::array<double, 8>>& collision_stiffness)
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
		mesh_stiffness[i][0] = cloth[i].single_cloth_info_ref.length_stiffness[0];
		mesh_stiffness[i][1] = cloth[i].single_cloth_info_ref.bending_stiffness[0];
		mesh_stiffness[i][2] = cloth[i].single_cloth_info_ref.position_stiffness;
		mesh_stiffness[i][3] = cloth[i].single_cloth_info_ref.length_stiffness[1];
		mesh_stiffness[i][4] = cloth[i].single_cloth_info_ref.bending_stiffness[1];
	}

	simulation_parameter[0] = time_step;
	switch (use_method)
	{
	case PD_:
		simulation_parameter[1] = project_dynamic.gravity_;
		break;
	case XPBD_:
		simulation_parameter[1] = xpbd.gravity_;
		break;
	case NEWTON_:
		simulation_parameter[1] = newton_method.gravity_;
		break;
	case XPBD_SECOND_ORDER_LARGE_:
		simulation_parameter[1] = second_order_xpbd_large.gravity_;
		break;
	case XPBD_IPC_:
		simulation_parameter[1] = xpbd_ipc.gravity_;
		break;
	}
	collision_stiffness.resize(cloth_num);
	for (int i = 0; i < cloth_num; ++i) {
		memcpy(collision_stiffness[i].data(), cloth[i].single_cloth_info_ref.collision_stiffness, 32);
	}
}


void Scene::updateIterateSolverParameter(double rate, int itr_solver_method)
{
	if (use_method==PD_) {
		project_dynamic.updateIterateSolverParameter(rate);
		project_dynamic.itr_solver_method = itr_solver_method;
	}
}

void Scene::getTetrahedronInfo(std::vector<std::array<int, 3>>& mesh_info, std::vector<double>& mass, std::vector<std::array<double, 6>>& mesh_stiffness, double* simulation_parameter, std::vector<std::array<double, 8>>& collision_stiffness)
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
		mesh_stiffness[i][0] = tetrahedron[i].single_tetrahedron_info_ref.ARAP_stiffness[0];
		mesh_stiffness[i][1] = tetrahedron[i].single_tetrahedron_info_ref.position_stiffness;
		mesh_stiffness[i][2] = tetrahedron[i].single_tetrahedron_info_ref.edge_length_stiffness;
		mesh_stiffness[i][3] = tetrahedron[i].single_tetrahedron_info_ref.ARAP_stiffness[1];
	}
	simulation_parameter[0] = time_step;
	switch (use_method)
	{
	case PD_:
		simulation_parameter[1] = project_dynamic.gravity_;
		break;
	case XPBD_:
		simulation_parameter[1] = xpbd.gravity_;
		break;
	case NEWTON_:
		simulation_parameter[1] = newton_method.gravity_;
		break;
	case XPBD_SECOND_ORDER_LARGE_:
		simulation_parameter[1] = second_order_xpbd_large.gravity_;
		break;
	case XPBD_IPC_:
		simulation_parameter[1] = xpbd_ipc.gravity_;
		break;
	}
	collision_stiffness.resize(tetrahedron_num);
	for (int i = 0; i < tetrahedron_num; ++i) {
		memcpy(collision_stiffness[i].data(), tetrahedron[i].single_tetrahedron_info_ref.collision_stiffness, 32);
	}
}

void Scene::drawScene(Camera* camera, std::vector<std::vector<bool>>& show_element,
	bool* control_parameter)
{
	shadow.drawShadow(camera, show_element, cloth, collider, tetrahedron);
	glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_CUBE_MAP, shadow.depth_map);
	glEnable(GL_CULL_FACE);

	if (control_parameter[SHARP_EDGE_SHADING]) {
		object_shader_front = object_shader_front_sharp_edge;
	}
	else {
		object_shader_front = object_shader_front_soft_edge;
	}

	for (int j = 0; j < cloth_num; ++j) {
		if (!show_element[CLOTH_][j]) {
			if (!cloth[j].mesh_struct.triangle_indices.empty()) {
				cloth[j].setSceneShader(light, camera, shadow.far_plane, object_shader_front);
				cloth[j].draw(camera, object_shader_front);
			}
			else {
				cloth[j].drawWireframe(camera, wireframe_shader);
			}
		}
	}
	glDisable(GL_CULL_FACE);
	for (int j = 0; j < tetrahedron_num; ++j) {
		if (!show_element[TETRAHEDRON_][j]) {
			tetrahedron[j].setSceneShader(light, camera, shadow.far_plane, object_shader_front);
			tetrahedron[j].draw(camera, object_shader_front);
		}
	}
	for (int j = 0; j < collider.size(); ++j) {
		if (!show_element[COLLIDER_][j]) {
			collider[j].setSceneShader(light, camera, shadow.far_plane, object_shader_front);
			collider[j].draw(camera, object_shader_front);
		}
	}
	for (int i = 0; i < collider.size(); ++i) {
		if (show_element[3 + COLLIDER_][i]) {
			collider[i].drawWireframe(camera, wireframe_shader);
		}
	}
	for (int j = 0; j < cloth_num; ++j) {
		if (show_element[3 + CLOTH_][j]) {
			cloth[j].drawWireframe(camera, wireframe_shader);
		}
	}
	for (int j = 0; j < tetrahedron_num; ++j) {
		if (show_element[3 + TETRAHEDRON_][j]) {
			tetrahedron[j].drawWireframe(camera, wireframe_shader);
		}
	}


	//project_dynamic.collision.draw_culling.drawTetTriangle(camera, wireframe_shader, light, shadow.far_plane);
	//project_dynamic.collision.draw_culling.draw(camera, shadow.far_plane);

	if (floor.show) {
		floor.draw(camera, object_shader_texture, &shadow, light, shadow.far_plane);
	}
	if (control_parameter[MOVE_OBJ] || control_parameter[ONLY_MOVE_CURRENT_POSITION]) {
		if (intersection.happened) {
			object_chosen_indicator.draw(wireframe_shader, camera, select_dimension_index,3.0);
		}
	}
	else {
		if (intersect_when_rotation) {
			object_chosen_indicator.draw(wireframe_shader, camera, select_dimension_index,3.0);
		}
	}
	//draw_culling.drawCell(camera, wireframe_shader);

	if (control_parameter[ONLY_COLLISION_TEST]) {
		test_draw_collision.drawCollision(control_parameter[DRAW_VT], light, camera, object_shader_front, show_element, &shadow, wireframe_shader,
			!control_parameter[DRAW_SPATIAL_HASHING], control_parameter[DRAW_ALL_PAIRS_IN_A_CELL]);
	}
	if (intersection.happened && (!(control_parameter[MOVE_OBJ]|| control_parameter[ONLY_MOVE_CURRENT_POSITION]
		|| control_parameter[ROTATION] || control_parameter[ONLY_ROTATE_CURRENT]))) {
		cursor.draw(camera);
	}

	if (control_parameter[ONLY_COLLISION_TEST]) {
		if (control_parameter[DRAW_SPATIAL_HASHING]) {
			test_draw_collision.draw_spatial_hashing.drawCell(camera, wireframe_shader);
			test_draw_collision.draw_spatial_hashing.drawCellSelect(camera, wireframe_shader);
			test_draw_collision.draw_spatial_hashing.drawCellSelectOne(camera, wireframe_shader);
		}
	}


	if (control_parameter[SAVE_OBJ]) {
		saveObj();
	}

	if (control_parameter[SAVE_SIMULATION_DATA]) {

	}

	//std::vector<std::array<double, 3>> pos;
	//pos.push_back(tetrahedron[0].mesh_struct.vertex_position[tetrahedron[0].mesh_struct.edge_vertices[6302 * 2]]);
	//pos.push_back(tetrahedron[0].mesh_struct.vertex_position[tetrahedron[0].mesh_struct.edge_vertices[6302 * 2+1]]);

	//pos.push_back(tetrahedron[0].mesh_struct.vertex_position[tetrahedron[0].mesh_struct.edge_vertices[3533 * 2]]);
	//pos.push_back(tetrahedron[0].mesh_struct.vertex_position[tetrahedron[0].mesh_struct.edge_vertices[3533 * 2 + 1]]);



	//draw_vertex.setVertex(pos,//collision.draw_target_position,
	//	0.005);
	////draw_vertex.setVertex(tetrahedron[0].mesh_struct.vertex_position[project_dynamic.collision.chosen_show_vertex].data(),
	////	0.005);
	//draw_vertex.draw(camera, glm::vec3(1.0, 0.0, 0.0));

	//draw_vertex.setVertex(project_dynamic.collision.draw_target_position, //
	//	0.005);
	//draw_vertex.draw(camera, glm::vec3(0.0, 1.0, 0.0));
	//draw_triangle.drawTriangle(camera, object_shader_front, collider[0].mesh_struct.vertex_for_render,
	//	collider[0].mesh_struct.triangle_indices, collider[0].mesh_struct.face_normal_for_render,
	//	project_dynamic.collision.test_triangle_index, glm::vec3(0.0, 1.0, 0.0));



	//draw_edge_.drawEdge(camera, wireframe_shader, cloth[0].mesh_struct.vertex_position, cloth[0].mesh_struct.unconnected_edge_index[0], glm::vec3(1.0, 0.0, 0.0), cloth[0].mesh_struct.edge_vertices);
}

void Scene::saveScene()
{
	switch (use_method)
	{
	case PD_:
		break;
	case XPBD_:
		xpbd.saveScene(force_direction, intersection.obj_No,intersection.happened);
		break;
	case NEWTON_:
		newton_method.saveScene(force_direction, intersection.obj_No, intersection.happened);
		break;
	case XPBD_IPC_:
		xpbd_ipc.saveScene(force_direction, intersection.obj_No, intersection.happened);
		break;
	}
}


void Scene::readScene(std::string& path)
{
	switch (use_method)
	{
	case PD_:
		break;
	case XPBD_:
		xpbd.readScene(path.c_str(),force_direction,intersection.obj_No);
		break;
	case NEWTON_:
		newton_method.readScene(path.c_str(), force_direction, intersection.obj_No);
		break;
	case XPBD_IPC_:
		xpbd_ipc.readScene(path.c_str(), force_direction, intersection.obj_No);
		break;
	}
}


void Scene::saveObj()
{
	if (time_stamp != last_output_obj_stamp) {
		for (unsigned int i = 0; i < cloth.size(); ++i) {
			save_obj.write(cloth[i].mesh_struct.vertex_position, cloth[i].mesh_struct.triangle_indices, 10, time_stamp, i,"object",cloth[i].material.front_material);
		}
		for (unsigned int i = 0; i < tetrahedron.size(); ++i) {
			save_obj.write(tetrahedron[i].mesh_struct.vertex_position, tetrahedron[i].mesh_struct.triangle_indices, 10, time_stamp, i+cloth.size(), "object",tetrahedron[i].material);
		}
		for (unsigned int i = 0; i < collider.size(); ++i) {
			save_obj.write(collider[i].mesh_struct.vertex_position, collider[i].mesh_struct.triangle_indices, 10, time_stamp, i, "collider",collider[i].material.front_material);
		}
		last_output_obj_stamp = time_stamp;
	}
}

void Scene::initial()
{
	//if (!only_test_collision) {
		switch (use_method)
		{
		case PD_:
			project_dynamic.initial();
			break;
		case XPBD_:
			xpbd.initial();
			break;
		case NEWTON_:
			newton_method.initial();
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			second_order_xpbd_large.initial();
			break;
		case XPBD_IPC_:
			xpbd_ipc.initial();
			break;
		}
	//}
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
	time_indicate_for_simu = 0;
}

void Scene::reset()
{

	//if (!only_test_collision) {
		switch (use_method)
		{
		case PD_:
			project_dynamic.reset();
			break;
		case XPBD_:
			xpbd.reset();
			break;
		case NEWTON_:
			newton_method.initial();
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			second_order_xpbd_large.initial();
			break;
		case XPBD_IPC_:
			xpbd_ipc.reset();
			break;
		}
	//}
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].reset();
	}



	for (int i = 0; i < tetrahedron.size(); ++i) {
		tetrahedron[i].reset(use_method);
	}
	intersection.initialIntersection();
	//spatial_hashing.time_stamp++;
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].reset();
	}
	time_stamp = 0;
	time_indicate_for_simu = 0;
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


void Scene::updateObjSimulation(Camera* camera, double* cursor_screen, bool* control_parameter, float force_coe, bool& record_matrix,
	double& ave_iteration)
{

	if (control_parameter[START_SIMULATION] || control_parameter[ONE_FRAME]) {
		switch (use_method)
		{
		case PD_:
			project_dynamic.resetExternalForce();
			break;
		case XPBD_:
			xpbd.resetExternalForce();
			break;
		case NEWTON_:
			newton_method.resetExternalForce();
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			second_order_xpbd_large.resetExternalForce();
			break;
		case XPBD_IPC_:
			xpbd_ipc.resetExternalForce();
			break;
		}
	}

		if (read_force) {
			addExternalForce();
			read_force = false;
		}
		if (intersection.happened && !control_parameter[START_TEST]) {
			setCursorForce(camera, cursor_screen, force_coe);
		}

	if (control_parameter[START_SIMULATION] || control_parameter[ONE_FRAME]) {
		if (control_parameter[START_TEST]) {
			time_indicate_for_simu++;
			//if (!collider.empty()) {
			//	move_model.moveSphere(time_indicate_for_simu, collider[0].mesh_struct.vertex_for_render, collider[0].mesh_struct.vertex_position, 1.0);
			//}
			//move_model.moveSkirt(time_indicate_for_simu, mesh_struct, use_method == PD_, 1.0);

			//rorate band capsule
			if (!collider.empty()) {
				move_model.sceneRotateCapsule(time_indicate_for_simu, collider[0].mesh_struct.vertex_for_render, collider[0].mesh_struct.vertex_position, &cloth[0].mesh_struct, use_method == PD_,1.0);
			}
		}
		auto t0 = std::chrono::system_clock::now();

		

		switch (use_method)
		{
		case PD_:
			project_dynamic.PDsolve();
			//project_dynamic.PD_IPC_solve(record_matrix);
			break;
		case XPBD_:
			xpbd.PBDsolve();
			break;
		case NEWTON_:
			newton_method.solveNewtonMethod();
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			second_order_xpbd_large.solveNewtonMethod();
			break;
		case XPBD_IPC_:
			//xpbd_ipc.XPBD_IPCSolve();
			//xpbd_ipc.XPBD_IPC_Block_Solve();
			xpbd_ipc.XPBD_IPC_Block_Solve_Multithread();
			break;
		}

		if (control_parameter[START_TEST]) {
			move_model.updateColliderPosition(collider);
		}
		auto t1 = std::chrono::system_clock::now();
		//if (time_stamp % time_record_interval == 0) {
		//	time_accumulation = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0);
		//}
		//else if (time_stamp % time_record_interval < (time_record_interval - 1)) {
		//	time_accumulation += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);;
		//}
		//else {
		//	time_accumulation += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
		//	*time_per_frame =double(time_accumulation.count())* std::chrono::microseconds::period::num/ std::chrono::microseconds::period::den*1000.0 / (double)time_record_interval;
		//}

		*time_per_frame = double(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den * 1000.0;

		std::cout << "time " << time_stamp << " " << *time_per_frame << std::endl;

		//project_dynamic.PD_IPC_solve(record_matrix);
		//project_dynamic.update_ave_iteration_record(ave_iteration);
		time_stamp++;

		if (control_parameter[ONLY_COLLISION_TEST]) {
			test_draw_collision.draw_collision.setElementIndices();
		}
		//checkVolume();

	}
}

void Scene::setMove(Camera* camera, Input* input, bool* control_parameter)
{
	if (!(control_parameter[ONLY_ROTATE_CURRENT] || control_parameter[ROTATION])) {
		intersect_when_rotation = false;
	}

	if (control_parameter[MOVE_OBJ]) {
		if (intersection.happened_include_collider) {
			moveObj(camera, input->mouse.screen_pos, false);
			switch (use_method)
			{
			case PD_:
				project_dynamic.updateSystemPos();
				break;
			}
		}
		setChosenIndicator();
	}
	if (control_parameter[ROTATION]) {
		if (intersect_when_rotation) {

			std::cout << input->mouse.angle[0] << " " << input->mouse.angle[1] << std::endl;

			rotate(camera, input->mouse.angle, false);
			switch (use_method)
			{
			case PD_:
				project_dynamic.updateSystemPos();
				break;
			}
		}
		setChosenIndicator();
	}
	if (control_parameter[ONLY_MOVE_CURRENT_POSITION]) {
		if (intersection.happened_include_collider) {
			moveObj(camera, input->mouse.screen_pos, true);
			switch (use_method)
			{
			case PD_:
				project_dynamic.updateSystemPos();
				break;
			}
		}
		setChosenIndicator();
	}
	if (control_parameter[ONLY_ROTATE_CURRENT]) {
		if (intersect_when_rotation) {
			rotate(camera, input->mouse.angle, true);
			switch (use_method)
			{
			case PD_:
				project_dynamic.updateSystemPos();
				break;
			}
		}
		setChosenIndicator();
	}
}


void Scene::drawSpatialHashing(Camera* camera, Input* input, int& select_hash_cell)
{
	if (control_parameter[SPATIAL_HASHING_UPDATE]) {
		getAABB();
		test_draw_collision.setForOriSpatialHashing();

		control_parameter[SPATIAL_HASHING_UPDATE] = false;
	}

	if (control_parameter[DRAW_SPATIAL_HASHING]) {
		if (input->mouse.leftButtonIsPressed() && !control_parameter[SEARCH_LEFT_SH_CELL]
			&& !control_parameter[SEARCH_RIGHT_SH_CELL] && !input->keyboard.keyIsPressed(GLFW_KEY_LEFT_CONTROL)
			&& !input->keyboard.keyIsPressed(GLFW_KEY_RIGHT_CONTROL)) {
			select_hash_cell = 0;
			test_draw_collision.obtianSpatialHashingCell(camera, input->mouse.screen_pos);
			test_draw_collision.setForSelectCell();
			test_draw_collision.obtainElementsInOneCell(select_hash_cell);
			test_draw_collision.setForSelectCellOne(select_hash_cell);
		}
		if (control_parameter[SEARCH_LEFT_SH_CELL]) {
			select_hash_cell--;
			test_draw_collision.obtainElementsInOneCell(select_hash_cell);
			test_draw_collision.setForSelectCellOne(select_hash_cell);
			control_parameter[SEARCH_LEFT_SH_CELL] = false;
		}
		if (control_parameter[SEARCH_RIGHT_SH_CELL]) {
			select_hash_cell++;
			test_draw_collision.obtainElementsInOneCell(select_hash_cell);
			test_draw_collision.setForSelectCellOne(select_hash_cell);
			control_parameter[SEARCH_RIGHT_SH_CELL] = false;
		}

	}
}


void Scene::updateSceneCollisionTest(bool* control_parameter)
{
	if (control_parameter[ONE_FRAME]) {
		test_draw_collision.setCollisionData();
	}

}


void Scene::updateCloth(Camera* camera, Input* input, bool* control_parameter, float force_coe, bool& record_matrix,
	double& ave_iteration, int& select_hash_cell_index)
{
	setMove(camera, input, control_parameter);
	if (use_method < 100) {
		updateObjSimulation(camera, input->mouse.screen_pos, control_parameter, force_coe, record_matrix, ave_iteration);

	}
	else {
		if (control_parameter[ONLY_COLLISION_TEST]) {
			updateSceneCollisionTest(control_parameter);
		}
	}
	drawSpatialHashing(camera, input, select_hash_cell_index);
	updateBufferOriPos();

	updateBuffer();
}


void Scene::setChosenIndicator()
{
	if (intersect_when_rotation){
		if (input->mouse.leftButtonWasPressedPreviousAndThisFrame()) {


			if (move_object.select_object_index < cloth.size()) {
				object_chosen_indicator.updatePosition(cloth[move_object.select_object_index].center, cloth[move_object.select_object_index].move_circle_radius, cloth[move_object.select_object_index].rotation_matrix);
				cloth[move_object.select_object_index].mesh_struct.setAnchorPosition();
			}
			else if (move_object.select_object_index < obj_num_except_collider) {
				object_chosen_indicator.updatePosition(tetrahedron[move_object.select_object_index - cloth.size()].center, tetrahedron[move_object.select_object_index - cloth.size()].move_circle_radius, tetrahedron[move_object.select_object_index - cloth.size()].rotation_matrix);
				tetrahedron[move_object.select_object_index - cloth.size()].mesh_struct.setAnchorPosition();
			}
			else {
				object_chosen_indicator.updatePosition(collider[move_object.select_object_index - obj_num_except_collider].center, collider[move_object.select_object_index - obj_num_except_collider].move_circle_radius, collider[move_object.select_object_index - obj_num_except_collider].rotation_matrix);
			}
		}
	}
	else if (intersection.happened) {
		if (intersection.obj_No < cloth.size()) {
			object_chosen_indicator.updatePosition(cloth[intersection.obj_No].center, cloth[intersection.obj_No].move_circle_radius, cloth[intersection.obj_No].rotation_matrix);
			cloth[move_object.select_object_index].mesh_struct.setAnchorPosition();
		}
		else if (intersection.obj_No < obj_num_except_collider) {
			object_chosen_indicator.updatePosition(tetrahedron[intersection.obj_No - cloth.size()].center, tetrahedron[intersection.obj_No - cloth.size()].move_circle_radius, tetrahedron[intersection.obj_No - cloth.size()].rotation_matrix);
			tetrahedron[move_object.select_object_index - cloth.size()].mesh_struct.setAnchorPosition();
		}
		else {
			object_chosen_indicator.updatePosition(collider[intersection.obj_No - obj_num_except_collider].center, collider[intersection.obj_No - obj_num_except_collider].move_circle_radius, collider[intersection.obj_No - obj_num_except_collider].rotation_matrix);
		}
	}
}

void Scene::updateBufferOriPos()
{
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].setBufferOriPos();
	}
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].setBufferOriPos();
	}
	for (int i = 0; i < tetrahedron.size(); ++i) {
		tetrahedron[i].setBufferOriPos();
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
	//if (!only_test_collision) {
		switch (use_method)
		{
		case PD_:
			project_dynamic.initialDHatTolerance(ave_edge_length);
			break;
		case XPBD_:
			xpbd.initialDHatTolerance(ave_edge_length);
			break;
		case NEWTON_:
			newton_method.initialDHatTolerance(ave_edge_length);
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			second_order_xpbd_large.initialDHatTolerance(ave_edge_length);
			break;
		case XPBD_IPC_:
			xpbd_ipc.initialDHatTolerance(ave_edge_length);
			break;
		}
	//}
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

void Scene::updateAnchorTet()
{
	for (int i = 0; i < tetrahedron.size(); ++i) {
		tetrahedron[i].mesh_struct.setAnchorPosition();
		tetrahedron[i].mesh_struct.updateAnchorPerThread(thread.thread_num);
		tetrahedron[i].mesh_struct.updateUnfixedPointData();
		tetrahedron[i].mesh_struct.resetMassInv();
		tetrahedron[i].mesh_struct.updateTetNeighborInfo();
	}
	if (use_method == XPBD_SECOND_ORDER_LARGE_) {
		second_order_xpbd_large.updateIndexBeginPerObj();
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
		select_anchor[1] = false;
		switch (use_method)
		{
		case PD_:
			for (int i = 0; i < tetrahedron.size(); ++i) {
				tetrahedron[i].mesh_struct.setAnchorPosition();
				tetrahedron[i].mesh_struct.updateAnchorPerThread(thread.thread_num);
				//tetrahedron[i].mesh_struct.updateUnfixedPointData();
			}
			project_dynamic.updateTetrahedronAnchorVertices();
			break;
		case NEWTON_:
			for (int i = 0; i < tetrahedron.size(); ++i) {
				tetrahedron[i].mesh_struct.setAnchorPosition();
				tetrahedron[i].mesh_struct.updateAnchorPerThread(thread.thread_num);
				tetrahedron[i].mesh_struct.updateUnfixedPointData();
				tetrahedron[i].mesh_struct.resetMassInv();
				newton_method.updateIndexBeginPerObj();
			}
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			updateAnchorTet();
			break;
		case XPBD_:
			xpbd.updateTetrahedronAnchorVertices();
			for (unsigned int i = 0; i < tetrahedron.size(); ++i) {
				tetrahedron[i].mesh_struct.updateTetNeighborInfo();
			}
			
			break;
		case XPBD_IPC_:
			xpbd_ipc.updateTetrahedronAnchorVertices();
			for (unsigned int i = 0; i < tetrahedron.size(); ++i) {
				tetrahedron[i].mesh_struct.updateTetNeighborInfo();
			}
			break;
		}
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


void Scene::pickAxes(double* pos, Camera* camera)
{
	if (intersect_when_rotation) {
		int mouse_pos[2];
		mouse_pos[0] = pos[0];
		mouse_pos[1] = SCR_HEIGHT - pos[1];
		if (object_chosen_indicator.pickAxes(wireframe_shader, camera, select_dimension_index, mouse_pos)) {
			std::cout <<"pick axe "<<  select_dimension_index << std::endl;
		}
	}
}


void Scene::obtainCursorIntersection(double* pos, Camera* camera, std::vector<std::vector<bool>>& hide)
{
	int mouse_pos[2];
	mouse_pos[0] = pos[0];
	mouse_pos[1] = SCR_HEIGHT - pos[1];
	int chosen_index[2];
	pick_triangle.pickTriangle(&cloth, &collider, &tetrahedron, camera, hide, chosen_index, mouse_pos);
	intersection.initialIntersection();
	if (chosen_index[0] > -1) {
		//double cursor_pos[3];
		intersection.setIntersection(chosen_index);
		if (chosen_index[1] < obj_num_except_collider) {
			intersection.happened = true;
			intersection.first_intersection_frame = true;
			if (control_parameter[ROTATION] || control_parameter[ONLY_ROTATE_CURRENT]) {
				intersect_when_rotation = true;
				move_object.select_object_index = intersection.obj_No;
			}
		}
	}
}

void Scene::getCursorPos(double* cursor_pos, std::vector<std::array<double, 3>>& vertex, int* vertex_index)
{
	for (int i = 0; i < 3; ++i) {
		cursor_pos[i] = (vertex[vertex_index[0]][i] + vertex[vertex_index[1]][i] + vertex[vertex_index[2]][i]) / 3.0;
	}
}


void Scene::rotate(Camera* camera, float* angle, bool only_move_vertex_pos)
{
	if (input->mouse.leftButtonWasPressedPreviousAndThisFrame()) {
		float angle_move;
		if (abs(angle[0]) > abs(angle[1])) {
			angle_move = -angle[0];
			//if (sameDirection(camera, move_object.select_object_index)) {				
			//}
			//else {
			//	angle_move = angle[0];
			//}
		}
		else {
			angle_move = -angle[1];
			//if (sameDirection(camera, move_object.select_object_index)) {				
			//}
			//else {
			//	angle_move = angle[1];
			//}
		}
		move_object.rotation(angle_move, select_dimension_index, only_move_vertex_pos);
	}
}

bool Scene::sameDirection(Camera* camera, unsigned int obj_index)
{
	double camera_y[3] = { camera->up.x, camera->up.y, camera->up.z };
	double rotation_axis[3];
	if (obj_index < cloth.size()) {
		rotation_axis[0] = cloth[obj_index].rotation_matrix[0];
		rotation_axis[1] = cloth[obj_index].rotation_matrix[3];
		rotation_axis[2] = cloth[obj_index].rotation_matrix[6];
	}
	else if (obj_index < obj_num_except_collider) {
		rotation_axis[0] = tetrahedron[obj_index - cloth.size()].rotation_matrix[0];
		rotation_axis[1] = tetrahedron[obj_index - cloth.size()].rotation_matrix[3];
		rotation_axis[2] = tetrahedron[obj_index - cloth.size()].rotation_matrix[6];
	}
	else {
		rotation_axis[0] = collider[obj_index - obj_num_except_collider].rotation_matrix[0];
		rotation_axis[1] = collider[obj_index - obj_num_except_collider].rotation_matrix[3];
		rotation_axis[2] = collider[obj_index - obj_num_except_collider].rotation_matrix[6];
	}
	if (DOT(camera_y, rotation_axis) < 0) {
		return false;
	}
	return true;
}





void Scene::moveObj(Camera* camera, double* cursor_screen, bool only_move_vertex_pos)
{
	double cursor_pos[3];	
	if (intersection.obj_No < cloth.size()) {
		getCursorPos(cursor_pos, cloth[intersection.obj_No].mesh_struct.vertex_position,
			cloth[intersection.obj_No].mesh_struct.triangle_indices[intersection.face_index].data());

	}
	else if (intersection.obj_No < obj_num_except_collider) {
		getCursorPos(cursor_pos, tetrahedron[intersection.obj_No - cloth.size()].mesh_struct.vertex_position,
			tetrahedron[intersection.obj_No - cloth.size()].mesh_struct.triangle_indices[intersection.face_index].data());
	}
	else {
		getCursorPos(cursor_pos, collider[intersection.obj_No - obj_num_except_collider].mesh_struct.vertex_position,
			collider[intersection.obj_No - obj_num_except_collider].mesh_struct.triangle_indices[intersection.face_index].data());
	}

	double cursor_pos_in_space[3];
	camera->getCursorPosInSpace(cursor_pos_in_space, cursor_screen, cursor_pos);
	double displacement[3];
	SUB(displacement, cursor_pos_in_space, cursor_pos);
	move_object.move(intersection.obj_No, displacement, only_move_vertex_pos);
}



void Scene::setObjMoveInfo(Camera* camera, double* cursor_screen)
{
	moveObj(camera, cursor_screen, false);
	for (unsigned int j = 0; j < cloth.size(); ++j) {
		cloth[j].ori_vertices = cloth[j].mesh_struct.vertex_position;
	}
	for (unsigned int j = 0; j < tetrahedron.size(); ++j) {
		tetrahedron[j].ori_vertices = tetrahedron[j].mesh_struct.vertex_position;
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
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].reset();
	}
}

void Scene::setCursorForce(Camera* camera, double* cursor_screen, float force_coe)
{
	double cursor_pos[3];

	if (intersection.obj_No < cloth.size()) {
		getCursorPos(cursor_pos, cloth[intersection.obj_No].mesh_struct.vertex_position,
			cloth[intersection.obj_No].mesh_struct.triangle_indices[intersection.face_index].data());
		if (intersection.first_intersection_frame) {
			cloth[intersection.obj_No].findAllNeighborVertex(intersection.face_index, cursor_pos, ave_edge_length);
			intersection.first_intersection_frame = false;
		}
	}
	else {
		getCursorPos(cursor_pos, tetrahedron[intersection.obj_No - cloth.size()].mesh_struct.vertex_position,
			tetrahedron[intersection.obj_No - cloth.size()].mesh_struct.triangle_indices[intersection.face_index].data());

		int* indices = tetrahedron[intersection.obj_No - cloth.size()].mesh_struct.triangle_indices[intersection.face_index].data();
		std::cout << "chosen vertex " << indices[0] << " " << indices[1] << " " << indices[2]<<" "<< intersection.face_index <<" "<< intersection.obj_No << std::endl;

		if (intersection.first_intersection_frame) {
			tetrahedron[intersection.obj_No - cloth.size()].findAllNeighborVertex(intersection.face_index, cursor_pos, ave_edge_length);
			intersection.first_intersection_frame = false;
		}
	}


	double cursor_pos_in_space[3];
	cursorMovement(camera, cursor_screen, force_direction, force_coe, cursor_pos, cursor_pos_in_space);
	cursor.translate(cursor_pos, cursor_pos_in_space);

	addExternalForce();
}




void  Scene::addExternalForce()
{
	if (intersection.obj_No < cloth.size()) {
		switch (use_method)
		{
		case PD_:
			project_dynamic.addExternalClothForce(force_direction, cloth[intersection.obj_No].coe_neighbor_vertex_force, cloth[intersection.obj_No].neighbor_vertex, intersection.obj_No);
			break;
		case XPBD_:
			xpbd.addExternalForce(force_direction, cloth[intersection.obj_No].coe_neighbor_vertex_force, cloth[intersection.obj_No].neighbor_vertex, intersection.obj_No);
			break;
		case NEWTON_:
			newton_method.addExternalForce(force_direction, cloth[intersection.obj_No].coe_neighbor_vertex_force, cloth[intersection.obj_No].neighbor_vertex, intersection.obj_No);
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			second_order_xpbd_large.addExternalForce(force_direction, cloth[intersection.obj_No].coe_neighbor_vertex_force, cloth[intersection.obj_No].neighbor_vertex, intersection.obj_No);
			break;
		case XPBD_IPC_:
			xpbd_ipc.addExternalForce(force_direction, cloth[intersection.obj_No].coe_neighbor_vertex_force, cloth[intersection.obj_No].neighbor_vertex, intersection.obj_No);
			break;
		}
	}
	else {
		switch (use_method)
		{
		case PD_:
			project_dynamic.addExternalTetForce(force_direction, tetrahedron[intersection.obj_No - cloth.size()].coe_neighbor_vertex_force, tetrahedron[intersection.obj_No - cloth.size()].neighbor_vertex, intersection.obj_No - cloth.size());
			break;
		case XPBD_:
			xpbd.addExternalForce(force_direction, tetrahedron[intersection.obj_No - cloth.size()].coe_neighbor_vertex_force, tetrahedron[intersection.obj_No - cloth.size()].neighbor_vertex, intersection.obj_No);
			break;
		case NEWTON_:
			newton_method.addExternalForce(force_direction, tetrahedron[intersection.obj_No - cloth.size()].coe_neighbor_vertex_force, tetrahedron[intersection.obj_No - cloth.size()].neighbor_vertex, intersection.obj_No);
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			second_order_xpbd_large.addExternalForce(force_direction, tetrahedron[intersection.obj_No - cloth.size()].coe_neighbor_vertex_force, tetrahedron[intersection.obj_No - cloth.size()].neighbor_vertex, intersection.obj_No);
			break;
		case XPBD_IPC_:
			xpbd_ipc.addExternalForce(force_direction, tetrahedron[intersection.obj_No - cloth.size()].coe_neighbor_vertex_force, tetrahedron[intersection.obj_No - cloth.size()].neighbor_vertex, intersection.obj_No);
			break;
		}
	}
}


void Scene::cursorMovement(Camera* camera, double* cursor_screen, double* force_direction, float force_coe, double* object_position, double* cursor_pos_in_space)
{
	camera->getCursorPosInSpace(cursor_pos_in_space, cursor_screen, object_position);
	SUB(force_direction, cursor_pos_in_space, object_position);
	force_coe *= 100.0;
	MULTI(force_direction, force_direction, force_coe);
	double force_magnitude =sqrt(DOT(force_direction, force_direction));
	if (force_magnitude > max_force_magnitude) {
		force_magnitude = max_force_magnitude / force_magnitude;
		MULTI(force_direction, force_direction, force_magnitude);
	}
}

void Scene::genShader()
{
	object_shader_front_soft_edge = new Shader("./shader/object_triangle.vs", "./shader/object_triangle.fs");
	*object_shader_front_soft_edge = Shader("./shader/object_triangle.vs", "./shader/object_triangle.fs");

	object_shader_front_sharp_edge = new Shader("./shader/object_tetrahedron.vs", "./shader/object_tetrahedron.fs", "./shader/object_tetrahedron.gs");
	*object_shader_front_sharp_edge = Shader("./shader/object_tetrahedron.vs", "./shader/object_tetrahedron.fs", "./shader/object_tetrahedron.gs");


	object_shader_texture = new Shader("./shader/object_triangle_texture.vs", "./shader/object_triangle_texture.fs");
	*object_shader_texture = Shader("./shader/object_triangle_texture.vs", "./shader/object_triangle_texture.fs");

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
}



void Scene::updateStiffness(UpdateObjStiffness& update_obj_stiffness, std::vector<std::array<double, 6>>& cloth_stiffness, std::vector<std::array<double, 6>>& tet_stiffness,
	std::vector<std::array<double, 8>>& cloth_collision_stiffness,
	std::vector<std::array<double, 8>>& tet_collision_stiffness)
{
	if (update_obj_stiffness.update_length) {
		for (unsigned int i = 0; i < cloth_num; ++i) {
			cloth[i].length_stiffness = update_obj_stiffness.length_stiffness[0];
			cloth[i].damp_length_stiffness = update_obj_stiffness.length_stiffness[1];
			//std::fill(cloth[i].length_stiffness.begin(), cloth[i].length_stiffness.end(), update_obj_stiffness.length_stiffness);
			cloth_stiffness[i][0] = update_obj_stiffness.length_stiffness[0];
			cloth_stiffness[i][3] = update_obj_stiffness.length_stiffness[1];
		}
		update_obj_stiffness.update_length = false;

	}
	if (update_obj_stiffness.update_bend) {
		for (unsigned int i = 0; i < cloth_num; ++i) {
			cloth[i].bend_stiffness = update_obj_stiffness.bend_stiffness[0];
			cloth[i].damp_bend_stiffness = update_obj_stiffness.bend_stiffness[1];
			cloth_stiffness[i][1] = update_obj_stiffness.bend_stiffness[0];
			cloth_stiffness[i][4] = update_obj_stiffness.bend_stiffness[1];
		}
		update_obj_stiffness.update_bend = false;
	}
	if (update_obj_stiffness.update_ARAP) {
		for (unsigned int i = 0; i < tetrahedron.size(); ++i) {
			tetrahedron[i].ARAP_stiffness = update_obj_stiffness.ARAP_stiffness[0];
			tetrahedron[i].damp_ARAP_stiffness = update_obj_stiffness.ARAP_stiffness[1];
			tet_stiffness[i][0]= update_obj_stiffness.ARAP_stiffness[0];
			tet_stiffness[i][3]= update_obj_stiffness.ARAP_stiffness[1];
		}
		update_obj_stiffness.update_ARAP = false;
	}
	if (update_obj_stiffness.update_tet_edge_length) {
		for (unsigned int i = 0; i < tetrahedron.size(); ++i) {
			tetrahedron[i].edge_length_stiffness = update_obj_stiffness.tet_edge_length_stiffness;
			tet_stiffness[i][2] = update_obj_stiffness.tet_edge_length_stiffness;
		}
	}

	for (unsigned int j = 0; j < 4; ++j) {
		if (update_obj_stiffness.update_collision[j]) {
			for (unsigned int i = 0; i < tetrahedron.size(); ++i) {
				tetrahedron[i].collision_stiffness[j] = update_obj_stiffness.collision_stiffness[j];
				tet_collision_stiffness[i][j] = update_obj_stiffness.collision_stiffness[j];
			}
			for (unsigned int i = 0; i < cloth.size(); ++i) {
				cloth[i].collision_stiffness[j] = update_obj_stiffness.collision_stiffness[j];
				cloth_collision_stiffness[i][j] = update_obj_stiffness.collision_stiffness[j];
			}
			update_obj_stiffness.update_collision[j] = false;
		}
	}
}


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



void Scene::getAABB()
{
	for (int i = 0; i < cloth.size(); ++i) {
		cloth[i].obtainAABB(true);
	}
	for (int i = 0; i < collider.size(); ++i) {
		collider[i].obtainAABB(true);
	}
	for (int i = 0; i < tetrahedron.size(); ++i) {
		tetrahedron[i].obtainAABB(true);
	}
	//thread->assignTask(&mesh_patch, PATCH_AABB);
}

void Scene::setDampStiffness(double* damp_stiffness, double* rayleigh_damp_stiffness)
{
	newton_method.damp_stiffness = damp_stiffness;
	newton_method.rayleigh_damp_stiffness = rayleigh_damp_stiffness;

	second_order_xpbd_large.damp_stiffness = damp_stiffness;
	second_order_xpbd_large.rayleigh_damp_stiffness = rayleigh_damp_stiffness;

}


void Scene::readScene()
{

}

//void Scene::getCurrentAABB()
//{
//	for (int i = 0; i < cloth.size(); ++i) {
//		cloth[i].obtainCurrentAABB();
//	}
//	for (int i = 0; i < collider.size(); ++i) {
//		collider[i].obtainCurrentAABB();
//	}
//	for (int i = 0; i < tetrahedron.size(); ++i) {
//		tetrahedron[i].obtainCurrentAABB();
//	}
//
//}