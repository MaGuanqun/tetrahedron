#include "simulation_main.h"
#include "./basic/Coordinate.h"
void simu_main(GLFWwindow* window, Input* input) {
	// Set up GUI
	BasicImGui basic_imgui;
	basic_imgui.initialImGui(window);
	ImGuiIO& io = ImGui::GetIO();
	(void)io;
	ImFontConfig config;
	io.Fonts->AddFontDefault(&config)->Scale = 1.2;

	glm::vec3 cameraPos = glm::vec3(0.0, 0.0, 1.0);
	Camera camera(cameraPos, normalize(glm::vec3(0.0f, 1.0f, 0.0f)));
	float zoom_value = 1.0;
	CoordinateSystem coordinateSystem;

	
	while (!glfwWindowShouldClose(window))
	{
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
		basic_imgui.imguiNewFrame();
	
		input->guiCaptureMouse = io.WantCaptureMouse;
		input->guiCaptureKeyboard = io.WantCaptureKeyboard;
		
		coordinateSystem.draw(&camera, cameraPos);
		basic_imgui.imguiEndFrame();
		glfwSwapBuffers(window);
		input->beginFrame();
		glfwPollEvents();
	}
	basic_imgui.imguiShutdown();
}