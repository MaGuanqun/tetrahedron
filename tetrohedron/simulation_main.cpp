#include "simulation_main.h"
#include "./basic/Coordinate.h"
#include"basic/enum_setting.h"
#include"render/imgui_windows.h"

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
	bool control_paramenter[11];
	memset(control_paramenter, 0, 11);
	ImGuiWindows imgui_windows;
	float force_coe=1.0;
	std::vector<std::vector<bool>> wireframe(2);
	std::vector<std::vector<bool>> hide(2);
	std::vector<std::string> collider_path;
	std::vector<std::string> tetrohedron_path;
	bool already_load_model = false;
	while (!glfwWindowShouldClose(window))
	{
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

			if (control_paramenter[INITIAL_CAMERA]) {
				camera.resetCam();
				zoom_value = 1.0;
				control_paramenter[INITIAL_CAMERA] = false;
			}
		}
		if (!already_load_model) {
			if (imgui_windows.loadModel(collider_path, tetrohedron_path)) {
				already_load_model = true;

			}
		}
		imgui_windows.controlWindow(control_paramenter, &force_coe);
		imgui_windows.visualizationControlPanel(control_paramenter[INITIAL_CAMERA], wireframe, hide);

		coordinateSystem.draw(&camera, cameraPos);

		basic_imgui.imguiEndFrame();
		glfwSwapBuffers(window);
		input->beginFrame();
		glfwPollEvents();
	}
	basic_imgui.imguiShutdown();
}