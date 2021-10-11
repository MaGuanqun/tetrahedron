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
	float force_coe=1.0;
	std::vector<std::vector<bool>> wireframe(2);
	std::vector<std::vector<bool>> hide(2);
	std::vector<std::string> collider_path;
	std::vector<std::string> object_path;
	bool already_load_model = false;
	Scene scene;
	scene.control_parameter = control_parameter;
	std::vector<std::array<int, 3>>cloth_info;//vertices,edges
	std::vector<std::array<int, 2>>tetrohedron_info;//vertices,edges
	std::vector<double>cloth_mass;//vertices,edges
	std::vector<double>tetrohedron_mass;//vertices,edges
	double time = 1.0;//fps, mass
	bool reset_camera = false;
	std::vector<std::array<double, 3>> cloth_stiffness;//stretch, bending, position, collision, fricition
	std::vector<std::array<double, 4>> collision_stiffness;//stretch, bending, position, collision, fricition
	double simulation_parameter[2] = {0.0,0.0};//timestep, gravity
	SaveImage save_image;
	int iteration_number[2] = { 0,0 };
	double convergence_rate[2] = { 0.1,0.1 };
	bool edit_PD_conv_rate = false;
	time_t start_time;
	
	bool set_stiffness[10];
	memset(set_stiffness, 0, 10);
	double temp_stiffness[6] = { 0.0,0.0,0.0,0.0,0.0,0.0 };
	UpdateStiffness update_stiffness;
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
			if (zoom_value > 0.025) {
				zoom_value -= 0.025 * input->mouse.scroll;
				zoom_value = myMax(zoom_value, 0.026);
				zoom_value = myMin(zoom_value, 1.3);
			}
			camera.zoomInOut((float)zoom_value * camera_from_origin);
		}

		if (input->mouse.mouse_callback) {
			if (input->mouse.leftButtonIsPressed()) {
				if (!input->mouse.rightButtonIsPressed()) {
					camera.rotation(input->mouse.angle[0], input->mouse.angle[1], 1);
				}
				if (input->mouse.rightButtonIsPressed()) {
					camera.move(-input->mouse.move_direction[0], -input->mouse.move_direction[1]);
				}
			}

			if (control_parameter[INITIAL_CAMERA]) {
				camera.resetCam();
				zoom_value = 1.0;
				control_parameter[INITIAL_CAMERA] = false;
			}
		}
		imgui_windows.operationWindow(cloth_stiffness, simulation_parameter, collision_stiffness, set_stiffness, temp_stiffness, update_stiffness);
		if (!already_load_model) {
			if (imgui_windows.loadModel(collider_path, object_path)) {
				already_load_model = true;
				scene.loadMesh(collider_path, object_path);
				glm::vec3 camera_pos = glm::vec3(0.6 * scene.shadow.camera_from_origin + scene.camera_center[0], scene.camera_center[1], -0.8 * scene.shadow.camera_from_origin + scene.camera_center[2]);
				camera.updateCamera(camera_pos, glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(scene.camera_center[0], scene.camera_center[1], scene.camera_center[2]));
				scene.getClothInfo(cloth_info, cloth_mass, cloth_stiffness, simulation_parameter, collision_stiffness);
				camera_from_origin = scene.shadow.camera_from_origin;				
				setHideWireframe(hide, wireframe, scene.collider.size(), scene.cloth.size() + scene.tetrohedron.size());
			}
		}
		else {
			scene.drawScene(&camera, wireframe, hide, control_parameter[SAVE_OBJ]);
		}
		imgui_windows.controlWindow(control_parameter, &force_coe);
		imgui_windows.visualizationControlPanel(control_parameter[INITIAL_CAMERA], wireframe, hide);

		coordinateSystem.draw(&camera, cameraPos);
		time = (double)(clock() - start_time);
		imgui_windows.infoWindow(cloth_info, cloth_mass, tetrohedron_info, tetrohedron_mass, time, iteration_number, convergence_rate, scene.time_stamp, edit_PD_conv_rate, control_parameter[START_SIMULATION]);


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


void setHideWireframe(std::vector<std::vector<bool>>& hide, std::vector<std::vector<bool>>& wireframe, int collider_num, int object_num)
{
	if (collider_num > 0) {
		hide[COLLIDER].resize(collider_num, false);
		wireframe[COLLIDER].resize(collider_num, false);
	}
	hide[OBJECT].resize(object_num, false);
	wireframe[OBJECT].resize(object_num, false);
}