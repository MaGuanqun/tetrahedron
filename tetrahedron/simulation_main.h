#pragma once
#include "./external/glad.h"
#include "./external/glfw3.h"
#include "./basic/opengl_input.h"
#include"basic/basic_imgui.h"

void simu_main(GLFWwindow* window, Input* input);

void setHideWireframe(std::vector<std::vector<bool>>& hide, std::vector<std::vector<bool>>& wireframe, std::vector<std::vector<bool>>& show_collision_element,
	int collider_num,  int cloth_num, int tetrahedron_num);