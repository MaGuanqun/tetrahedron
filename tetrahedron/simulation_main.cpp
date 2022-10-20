#define _CRT_SECURE_NO_DEPRECATE
#include "simulation_main.h"
#include "./basic/Coordinate.h"
#include"basic/enum_setting.h"
#include"render/imgui_windows.h"
#include"scene.h"
#include"render/saveImage.h"
#include"test_simu.h"

void simu_main(GLFWwindow* window, Input* input) {
	//test::testDampDerivative();
	// Set up GUI
	BasicImGui basic_imgui;
	basic_imgui.initialImGui(window);
	ImGuiIO& io = ImGui::GetIO();
	(void)io;
	ImFontConfig config;
	io.Fonts->AddFontDefault(&config)->Scale = 1.2;

	glm::vec3 cameraPos = glm::vec3(0.0, 0.0, 1.0);
	double camera_from_origin = 1.0;
	Camera camera(cameraPos, normalize(glm::vec3(0.0f, 1.0f, 0.0f)));
	float zoom_value = 1.0;
	CoordinateSystem coordinateSystem;
	bool control_parameter[32];
	memset(control_parameter, 0, 31);
	control_parameter[ONLY_COLLISION_TEST] =false;
	control_parameter[USE_XPBD] = false;
	control_parameter[USE_PD_] = true;
	control_parameter[USE_NEWTON_] = false;
	control_parameter[USE_XPBD_LARGE] = false;
	control_parameter[USE_XPBD_IPC] = false;
	control_parameter[DRAW_VT] = true;


	ImGuiWindows imgui_windows;
	float force_coe = 0.1;

	std::vector<std::vector<bool>> show_element(15); //0~2 hide, 3~5 wireframe, 6~8 show_colision_element, 9~11 show ori position 12~14 show ori position wireframe
	//std::vector<std::vector<bool>> wireframe(3);
	//std::vector<std::vector<bool>> hide(3);
	//std::vector<std::vector<bool>> show_collision_element(3);
	std::vector<std::string> collider_path;
	std::vector<std::string> object_path;
	std::string scene_path;
	bool already_load_model = false;
	Scene scene;
	scene.control_parameter = control_parameter;
	std::vector<std::array<int, 3>>cloth_info;//vertices,edges
	std::vector<std::array<int, 3>>tetrahedron_info;//vertices,edges
	std::vector<double>cloth_mass;//vertices,edges
	std::vector<double>tetrahedron_mass;//vertices,edges
	double time = 1.0;//fps, mass
	bool reset_camera = false;
	std::vector<std::array<double, 6>> cloth_stiffness;//stretch, bending, position, ARAP
	std::vector<std::array<double, 8>> cloth_collision_stiffness;//
	double simulation_parameter[2] = { 0.0,0.0 };//timestep, gravity
	SaveImage save_image;
	int iteration_number[2] = { 0,0 };
	int set_iteration_num[2] = { 100,1 }; //0:itr in a substep, 1:number of substep
	double convergence_rate[2] = { 0.1,0.1 };
	bool edit_PD_conv_rate = false;
	time_t start_time;

	std::vector<std::array<double, 6>> tetrahedron_stiffness;//ARAP, position, edge_length
	std::vector<std::array<double, 8>> tetrahedron_collision_stiffness;
	bool set_stiffness[13];
	memset(set_stiffness, 0, 13);
	std::vector<double> temp_stiffness(18);
	double temp_data[18] = {1e3,2e3,2e3,2e1,1e1,3e-5,2e2,1.0,0.0,
		0.0,0.0,
	//1e-3, 2e-3,2e-3,2e-3, 1e-3,1e-9,1e-2};
	0.0, 0.0,0.0,0.0, 0.0,0.0,0.0 };
	memcpy(temp_stiffness.data(), temp_data, 18 * 8);
	//memset(temp_stiffness, 0, 64);
	UpdateObjStiffness update_obj_stiffness;
	double tolerance_ratio[7] = { 5e-2,5e-2,5e-2,5e-2, 5e-2, 5e-2,5e-2 };

	double damp_stiffness = temp_stiffness[DAMP_STIFFNESS];
	double rayleigh_damp_stiffness[2] = { temp_stiffness[RAYLEIGH_DAMP_STIFFNESS_ALPHA], temp_stiffness[RAYLEIGH_DAMP_STIFFNESS_BETA] };
	scene.setDampStiffness(&damp_stiffness, rayleigh_damp_stiffness);



	bool set_anchor[2] = { false,false };

	const char* itr_solver_items_[] = { "Direct", "Jacobi","Chebyshev jacobi","Super Jacobi","Chebyshev super Jacobi", "Gauss Seidel", "PCG", "Chebyshev Gauss Seidel" };
	char* itr_solver_item = (char*)"Direct";
	char* itr_solver_items[IM_ARRAYSIZE(itr_solver_items_)];
	for (int i = 0; i < IM_ARRAYSIZE(itr_solver_items_); ++i) {
		itr_solver_items[i] = _strdup(itr_solver_items_[i]);
	}
	int use_itr_solver_method = DIRECT_SOLVE;

	double iteration_solver_conve_rate = 1e-7;
	bool record_matrix;

	double iteration_solver_iteration_num;


	unsigned int floor_dimension = 1;
	bool floor_control[4] = { false,false,true,false };
	double floor_value = -0.55;

	scene.input = input;

	int cell_index_chose = 0;

	double time_per_frame = 0;
	scene.time_per_frame = &time_per_frame;

	int time_record_inverval = 10;

	//glfwSwapInterval(0);

	std::string load_scene_path;


	double t1, t2;


	double friction_coe[3] = { 0.9,0.8,0.9 };//self, collider, floor
	unsigned int sub_step_per_detection = 1;

	while (!glfwWindowShouldClose(window))
	{

		//if (scene.time_stamp%time_record_inverval == 0) 
		//{
		//	t1 = clock();
		//}

		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
		basic_imgui.imguiNewFrame();

		input->guiCaptureMouse = io.WantCaptureMouse;
		input->guiCaptureKeyboard = io.WantCaptureKeyboard;

		imgui_windows.controlWindow(control_parameter, &force_coe);

		if (couldMoveCamera(input, control_parameter)) {
			if (input->mouse.scroll_callback) {
				if (!scene.intersection.happened && !set_anchor[0]) {
					if (zoom_value > 0.025) {
						zoom_value -= 0.025 * input->mouse.scroll;
						zoom_value = myMax(zoom_value, 0.026);
						zoom_value = myMin(zoom_value, 30.0);
					}
					camera.zoomInOut((float)zoom_value * camera_from_origin);
				}
			}
		}
		if (input->mouse.mouse_callback) {
			if (!scene.intersection.happened && !set_anchor[0]) {
				if (couldMoveCamera(input, control_parameter)) {
					if (input->mouse.leftButtonIsPressed()) {
						if (!input->mouse.rightButtonIsPressed()) {
							camera.rotation(input->mouse.angle[0], input->mouse.angle[1], 1);
						}
						if (input->mouse.rightButtonIsPressed()) {
							camera.move(-input->mouse.move_direction[0], -input->mouse.move_direction[1]);
						}
					}
				}
			}
			if (control_parameter[INITIAL_CAMERA]) {
				camera.resetCam();
				zoom_value = 1.0;
				control_parameter[INITIAL_CAMERA] = false;
			}
		}

		if (control_parameter[RESET_SIMULATION]) {
			for (int i = 3; i < 15; ++i) {
				std::fill(show_element[i].begin(), show_element[i].end(), false);
			}
			control_parameter[START_SIMULATION] = false;
			scene.reset();
			control_parameter[RESET_SIMULATION] = false;
			control_parameter[START_TEST] = false;
		}
		if (control_parameter[INITIAL_SIMULATION]) {
			for (int i = 3; i < 15; ++i) {
				std::fill(show_element[i].begin(), show_element[i].end(), false);
			}
			control_parameter[START_SIMULATION] = false;
			scene.initial();
			control_parameter[INITIAL_SIMULATION] = false;
			control_parameter[START_TEST] = false;
		}
		if (scene.intersection.happened && control_parameter[STOP_AFTER_RELEASE_MOUSE] && input->mouse.leftButtonWasReleasedThisFrame()) {
			control_parameter[START_SIMULATION] = false;
		}
		if (already_load_model) {
			scene.updateIterateSolverParameter(iteration_solver_conve_rate, use_itr_solver_method);


			if (input->mouse.leftButtonWasReleasedThisFrame() && !control_parameter[START_TEST]) {// 

				scene.resetIntersectionState();
			}
			if (input->mouse.leftButtonWasPressedThisFrame() && !input->mouse.rightButtonIsPressed()
				&& !control_parameter[START_TEST]) {
				scene.obtainCursorIntersection(input->mouse.screen_pos, &camera, show_element);
				scene.pickAxes(input->mouse.screen_pos, &camera);
			}

		}
		start_time = clock();
		if (control_parameter[USE_PD_] || control_parameter[USE_XPBD] || control_parameter[USE_NEWTON_]|| control_parameter[USE_XPBD_LARGE] || control_parameter[USE_XPBD_IPC]) {
			imgui_windows.operationWindow(cloth_stiffness, tetrahedron_stiffness, simulation_parameter, cloth_collision_stiffness, tetrahedron_collision_stiffness, set_stiffness, temp_stiffness.data(),
				update_obj_stiffness, set_anchor, !scene.tetrahedron.empty(), &damp_stiffness, rayleigh_damp_stiffness);
			if (update_obj_stiffness.update_length) {
				control_parameter[RESET_SIMULATION] = true;
			}
		}
		if (!already_load_model) {
			if (imgui_windows.loadModel(scene_path, collider_path, object_path)) {
				already_load_model = true;
				scene.loadMesh(scene_path, collider_path, object_path, tolerance_ratio, control_parameter,temp_stiffness.data(),friction_coe,&sub_step_per_detection,
					floor_control,floor_dimension,floor_value);
				//glm::vec3 camera_pos = glm::vec3(0.6 * scene.shadow.camera_from_origin + scene.camera_center[0], scene.camera_center[1], -0.8 * scene.shadow.camera_from_origin + scene.camera_center[2]);
				glm::vec3 camera_pos = glm::vec3(scene.camera_center[0], scene.camera_center[1], -scene.shadow.camera_from_origin + scene.camera_center[2]);
				camera.updateCamera(camera_pos, glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(scene.camera_center[0], scene.camera_center[1], scene.camera_center[2]));
				scene.getClothInfo(cloth_info, cloth_mass, cloth_stiffness, simulation_parameter, cloth_collision_stiffness);
				scene.getTetrahedronInfo(tetrahedron_info, tetrahedron_mass, tetrahedron_stiffness, simulation_parameter, tetrahedron_collision_stiffness);
				camera_from_origin = scene.shadow.camera_from_origin;
				setHideWireframe(show_element, scene.collider.size(), scene.cloth.size(), scene.tetrahedron.size());

				//scene.testBVH();
			}
		}
		else {
			//std::cout << camera.position.x<<" "<< camera.position.y<<" "<< camera.position.z << std::endl;
			if (input->keyboard.keyWasPressedThisFrame(GLFW_KEY_W))
			{
				record_matrix = true;
				control_parameter[SAVE_OBJ] = true;
			}

			scene.setFloorInfo(floor_control[0], floor_control[1], floor_control[2], floor_dimension, floor_value, floor_control[3], control_parameter);

			//if (scene.time_stamp == 20 || scene.time_stamp == 50 || scene.time_stamp == 100) {
			//	record_matrix = true;
			//	control_parameter[SAVE_OBJ] = true;
			//}
			//if (scene.time_stamp == 100) {
			//	record_matrix = false;
			//	control_parameter[SAVE_OBJ] = false;
			//}
			scene.updateStiffness(update_obj_stiffness, cloth_stiffness, tetrahedron_stiffness, cloth_collision_stiffness, tetrahedron_collision_stiffness);
			scene.updateItrInfo(set_iteration_num);
			scene.setTolerance(tolerance_ratio);
			scene.updateCloth(&camera, input, control_parameter, force_coe, record_matrix,
				iteration_solver_iteration_num, cell_index_chose);
			start_time = clock();
			scene.drawScene(&camera, show_element, control_parameter);
			scene.selectAnchor(control_parameter, set_anchor, input->mouse.screen_pos, input->mouse.left_press, input->mouse.prev_left_press, &camera, show_element[TETRAHEDRON_]);
			scene.obtainConvergenceInfo(convergence_rate, iteration_number);


			if (imgui_windows.loadScene(load_scene_path)) {
				scene.readScene(load_scene_path);
			}

			if (control_parameter[SAVE_SIMULATION_DATA]) {
				scene.saveScene();
			}

			if (control_parameter[SAVE_SCENE_DATA]) {
				scene.saveParameter(object_path, collider_path,cloth_stiffness,tetrahedron_stiffness,cloth_collision_stiffness,tetrahedron_collision_stiffness,
					tolerance_ratio,friction_coe);
				control_parameter[SAVE_SCENE_DATA] = false;
			}


		}

		if (control_parameter[ONE_FRAME]) {
			control_parameter[ONE_FRAME] = false;
		}

		imgui_windows.floorInfo(floor_control[0], floor_control[1], floor_control[2], floor_dimension, &floor_value, floor_control[3]);
	

		imgui_windows.visualizationControlPanel(control_parameter[INITIAL_CAMERA], show_element, control_parameter[ONLY_COLLISION_TEST], control_parameter);

		start_time = clock();
		coordinateSystem.draw(&camera, cameraPos);



		scene.drawSelectRange(set_anchor, input->mouse.left_press, input->mouse.prev_left_press);
		time = (double)(clock() - start_time);
		imgui_windows.infoWindow(cloth_info, cloth_mass, tetrahedron_info, tetrahedron_mass, time_per_frame, iteration_number, set_iteration_num, convergence_rate, scene.time_stamp, edit_PD_conv_rate, control_parameter[START_SIMULATION]);
		if (!control_parameter[ONLY_COLLISION_TEST]) {
			if (control_parameter[USE_PD_]) {
				imgui_windows.iterationSolverInfoWindow(iteration_solver_iteration_num, use_itr_solver_method, itr_solver_items, itr_solver_item,
					IM_ARRAYSIZE(itr_solver_items), &iteration_solver_conve_rate, record_matrix);
			}
		}

		if (control_parameter[ONLY_COLLISION_TEST]) {
			switchObjMoveMode(input, control_parameter);
		}
		if (input->keyboard.keyIsPressed(GLFW_KEY_S)) {
			control_parameter[ONE_FRAME] = true;
		}



		basic_imgui.imguiEndFrame();

		glfwSwapBuffers(window);
		input->beginFrame();
		glfwPollEvents();

		if (control_parameter[OUTPUT_IMAGE]) {
			save_image.outputImage(scene.time_stamp);
		}


		//if (scene.time_stamp > 0) {
		//	if ((scene.time_stamp-1) % time_record_inverval == (time_record_inverval - 1)) {
		//		t2 = (double)(clock() - t1) / (double)time_record_inverval;
		//		std::cout << "interval " << t2 << std::endl;
		//	}
		//}

	}
	basic_imgui.imguiShutdown();
}


void setHideWireframe(std::vector<std::vector<bool>>& show_element, int collider_num, int cloth_num, int tetrahedron_num)
{
	if (collider_num > 0) {
		for (unsigned int i = 0; i < 5; ++i) {
			show_element[COLLIDER_ + 3 * i].resize(collider_num, false);
		}
	}
	if (tetrahedron_num > 0) {
		for (unsigned int i = 0; i < 5; ++i) {
			show_element[TETRAHEDRON_ + 3 * i].resize(tetrahedron_num, false);
		}
	}
	if (cloth_num > 0) {
		for (unsigned int i = 0; i < 5; ++i) {
			show_element[CLOTH_ + 3 * i].resize(cloth_num, false);
		}
	}

}

bool couldMoveCamera(Input* input, bool* control_parameter)
{
	if (!(control_parameter[ROTATION] || control_parameter[ONLY_ROTATE_CURRENT])) {
		return true;
	}
	else {
		if (input->keyboard.keyIsPressed(GLFW_KEY_M)) {
			return true;
		}
	}
	return false;
}

void switchObjMoveMode(Input* input, bool* control_parameter)
{
	if (input->keyboard.keyIsPressed(GLFW_KEY_R) && !input->keyboard.keyIsPressed(GLFW_KEY_LEFT_CONTROL) && !input->keyboard.keyIsPressed(GLFW_KEY_RIGHT_CONTROL)) {
		control_parameter[ROTATION] = true;
		control_parameter[ONLY_MOVE_CURRENT_POSITION] = false;
		control_parameter[MOVE_OBJ] = false;
		control_parameter[ONLY_ROTATE_CURRENT] = false;
	}
	if ((input->keyboard.keyIsPressed(GLFW_KEY_LEFT_CONTROL) || input->keyboard.keyIsPressed(GLFW_KEY_RIGHT_CONTROL)) &&
		input->keyboard.keyIsPressed(GLFW_KEY_R)) {
		control_parameter[ONLY_ROTATE_CURRENT] = true;
		control_parameter[ONLY_MOVE_CURRENT_POSITION] = false;
		control_parameter[MOVE_OBJ] = false;
		control_parameter[ROTATION] = false;
	}
	if (input->keyboard.keyIsPressed(GLFW_KEY_T) && !input->keyboard.keyIsPressed(GLFW_KEY_LEFT_CONTROL) && !input->keyboard.keyIsPressed(GLFW_KEY_RIGHT_CONTROL)) {
		control_parameter[MOVE_OBJ] = true;
		control_parameter[ONLY_MOVE_CURRENT_POSITION] = false;
		control_parameter[ROTATION] = false;
		control_parameter[ONLY_ROTATE_CURRENT] = false;
	}
	if ((input->keyboard.keyIsPressed(GLFW_KEY_LEFT_CONTROL) || input->keyboard.keyIsPressed(GLFW_KEY_RIGHT_CONTROL)) &&
		input->keyboard.keyIsPressed(GLFW_KEY_T)) {
		control_parameter[ONLY_MOVE_CURRENT_POSITION] = true;
		control_parameter[ONLY_ROTATE_CURRENT] = false;
		control_parameter[MOVE_OBJ] = false;
		control_parameter[ROTATION] = false;
	}
}