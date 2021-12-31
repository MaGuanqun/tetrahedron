#define _CRT_SECURE_NO_DEPRECATE
#include "simulation_main.h"
#include "./basic/Coordinate.h"
#include"basic/enum_setting.h"
#include"render/imgui_windows.h"
#include"scene.h"
#include"render/saveImage.h"


void simu_main(GLFWwindow* window, Input* input) {
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
	bool control_parameter[11];
	memset(control_parameter, 0, 11);
	ImGuiWindows imgui_windows;
	float force_coe=5.0;
	std::vector<std::vector<bool>> wireframe(3);
	std::vector<std::vector<bool>> hide(3);
	std::vector<std::string> collider_path;
	std::vector<std::string> object_path;
	bool already_load_model = false;
	Scene scene;
	scene.control_parameter = control_parameter;
	std::vector<std::array<int, 3>>cloth_info;//vertices,edges
	std::vector<std::array<int, 3>>tetrahedron_info;//vertices,edges
	std::vector<double>cloth_mass;//vertices,edges
	std::vector<double>tetrahedron_mass;//vertices,edges
	double time = 1.0;//fps, mass
	bool reset_camera = false;
	std::vector<std::array<double, 3>> cloth_stiffness;//stretch, bending, position
	std::vector<std::array<double, 4>> cloth_collision_stiffness;//stretch, bending, position, collision, fricition
	double simulation_parameter[2] = {0.0,0.0};//timestep, gravity
	SaveImage save_image;
	int iteration_number[2] = { 0,0 };
	double convergence_rate[2] = { 0.1,0.1 };
	bool edit_PD_conv_rate = false;
	time_t start_time;

	std::vector<std::array<double, 2>> tetrahedron_stiffness;//ARAP, position
	std::vector<std::array<double, 4>> tetrahedron_collision_stiffness;//stretch, bending, position, collision, fricition
	bool set_stiffness[10];
	memset(set_stiffness, 0, 10);
	double temp_stiffness[6] = { 0.0,0.0,0.0,0.0,0.0,0.0 };
	UpdateClothStiffness update_cloth_stiffness;
	double tolerance_ratio[4] = {1e-1,1e-1,1e-1,1e-1 };

	bool set_anchor[2] = { false,false };

	const char* itr_solver_items_[] = { "Direct", "Jacobi","Chebyshev jacobi","Super Jacobi","Chebyshev super Jacobi", "Gauss Seidel", "PCG", "Chebyshev Gauss Seidel", "Weighted Jacobi"};
	char* itr_solver_item = (char*)"Direct";
	char* itr_solver_items[IM_ARRAYSIZE(itr_solver_items_)];
	for (int i = 0; i < IM_ARRAYSIZE(itr_solver_items_); ++i) {
		itr_solver_items[i] = _strdup(itr_solver_items_[i]);
	}
	int use_itr_solver_method = DIRECT_SOLVE;

	double iteration_solver_conve_rate = 1e-7;
	bool record_matrix;

	std::vector<std::vector<double>> iteration_solver_iteration_num;

	while (!glfwWindowShouldClose(window))
	{
		start_time = clock();
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
		basic_imgui.imguiNewFrame();
	
		input->guiCaptureMouse = io.WantCaptureMouse;
		input->guiCaptureKeyboard = io.WantCaptureKeyboard;
		
		if (input->mouse.scroll_callback) {
			if (!scene.intersection.happened && !set_anchor[0]) {
				if (zoom_value > 0.025) {
					zoom_value -= 0.025 * input->mouse.scroll;
					zoom_value = myMax(zoom_value, 0.026);
					zoom_value = myMin(zoom_value, 1.3);
				}
				camera.zoomInOut((float)zoom_value * camera_from_origin);
			}
		}

		if (input->mouse.mouse_callback) {
			if (!scene.intersection.happened && !set_anchor[0]) {
				if (input->mouse.leftButtonIsPressed()) {
					if (!input->mouse.rightButtonIsPressed()) {
						camera.rotation(input->mouse.angle[0], input->mouse.angle[1], 1);
					}
					if (input->mouse.rightButtonIsPressed()) {
						camera.move(-input->mouse.move_direction[0], -input->mouse.move_direction[1]);
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
			for (int i = 0; i < wireframe.size(); ++i) {
				std::fill(wireframe[i].begin(), wireframe[i].end(), false);
			}
			control_parameter[START_SIMULATION] = false;
			scene.resetCloth();
			control_parameter[RESET_SIMULATION] = false;
		}
		if (control_parameter[INITIAL_SIMULATION]) {
			for (int i = 0; i < wireframe.size(); ++i) {
				std::fill(wireframe[i].begin(), wireframe[i].end(), false);
			}
			control_parameter[START_SIMULATION] = false;
			scene.initialCloth();
			control_parameter[INITIAL_SIMULATION] = false;
		}
		if (scene.intersection.happened && control_parameter[STOP_AFTER_RELEASE_MOUSE] && input->mouse.leftButtonWasReleasedThisFrame()) {
			control_parameter[START_SIMULATION] = false;
		}
		if (already_load_model) {
			scene.updateIterateSolverParameter(iteration_solver_conve_rate, use_itr_solver_method);


			if (input->mouse.leftButtonWasReleasedThisFrame() && !control_parameter[START_TEST]) {// 

				scene.initialIntersection();
			}
			if (input->mouse.leftButtonWasPressedThisFrame() && !input->mouse.rightButtonIsPressed()
				&& !control_parameter[START_TEST]) {
				scene.obtainCursorIntersection(input->mouse.screen_pos, &camera, hide);
			}
		}

		imgui_windows.operationWindow(cloth_stiffness, simulation_parameter, cloth_collision_stiffness, set_stiffness, temp_stiffness, update_cloth_stiffness, set_anchor);
		if (!already_load_model) {
			if (imgui_windows.loadModel(collider_path, object_path)) {
				already_load_model = true;
				scene.loadMesh(collider_path, object_path);
				glm::vec3 camera_pos = glm::vec3(0.6 * scene.shadow.camera_from_origin + scene.camera_center[0], scene.camera_center[1], -0.8 * scene.shadow.camera_from_origin + scene.camera_center[2]);
				camera.updateCamera(camera_pos, glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(scene.camera_center[0], scene.camera_center[1], scene.camera_center[2]));
				scene.getClothInfo(cloth_info, cloth_mass, cloth_stiffness, simulation_parameter, cloth_collision_stiffness);
				scene.getTetrahedronInfo(tetrahedron_info, tetrahedron_mass, tetrahedron_stiffness, simulation_parameter, tetrahedron_collision_stiffness);
				camera_from_origin = scene.shadow.camera_from_origin;				
				setHideWireframe(hide, wireframe, scene.collider.size(), scene.cloth.size(), scene.tetrahedron.size());
				scene.setTolerance(tolerance_ratio);
				
				//scene.testBVH();
			}
		}
		else {		
			scene.setTolerance(tolerance_ratio);			
			scene.updateCloth(&camera, input->mouse.screen_pos, control_parameter, force_coe,record_matrix, iteration_solver_iteration_num);
			scene.drawScene(&camera, wireframe, hide, control_parameter[SAVE_OBJ]);
			scene.selectAnchor(control_parameter, set_anchor, input->mouse.screen_pos, input->mouse.left_press, input->mouse.prev_left_press, &camera, hide[TETROHEDRON_]);
			scene.obtainConvergenceInfo(convergence_rate, iteration_number);

		}


		if (control_parameter[ONE_FRAME]) {
			control_parameter[ONE_FRAME] = false;
		}

		imgui_windows.controlWindow(control_parameter, &force_coe);
		imgui_windows.visualizationControlPanel(control_parameter[INITIAL_CAMERA], wireframe, hide);

		coordinateSystem.draw(&camera, cameraPos);
		scene.drawSelectRange(set_anchor, input->mouse.left_press, input->mouse.prev_left_press);
		time = (double)(clock() - start_time);
		imgui_windows.infoWindow(cloth_info, cloth_mass, tetrahedron_info, tetrahedron_mass, time, iteration_number, convergence_rate, scene.time_stamp, edit_PD_conv_rate, control_parameter[START_SIMULATION]);
		imgui_windows.iterationSolverInfoWindow(iteration_solver_iteration_num, use_itr_solver_method, itr_solver_items, itr_solver_item,
			IM_ARRAYSIZE(itr_solver_items), &iteration_solver_conve_rate, record_matrix);

		basic_imgui.imguiEndFrame();
		glfwSwapBuffers(window);
		input->beginFrame();
		glfwPollEvents();

		if (control_parameter[OUTPUT_IMAGE]) {
			save_image.outputImage(scene.time_stamp);
		}
	}
	basic_imgui.imguiShutdown();
}


void setHideWireframe(std::vector<std::vector<bool>>& hide, std::vector<std::vector<bool>>& wireframe, int collider_num, int cloth_num, int tetrahedron_num)
{
	if (collider_num > 0) {
		hide[COLLIDER_].resize(collider_num, false);
		wireframe[COLLIDER_].resize(collider_num, false);
	}
	if (tetrahedron_num > 0) {
		hide[TETROHEDRON_].resize(tetrahedron_num, false);
		wireframe[TETROHEDRON_].resize(tetrahedron_num, false);
	}
	if (cloth_num > 0) {
		hide[CLOTH_].resize(cloth_num, false);
		wireframe[CLOTH_].resize(cloth_num, false);
	}
	
}