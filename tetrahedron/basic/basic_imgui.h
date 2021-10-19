#pragma once
#include"../external/imgui/imgui.h"
#include"../external/imgui/imgui_impl_glfw.h"
#include"../external/imgui/imgui_impl_opengl3.h"

class BasicImGui
{
public:
	void initialImGui(GLFWwindow* window);
	void imguiNewFrame();
	void imguiEndFrame();
	void imguiShutdown();
private:

};
