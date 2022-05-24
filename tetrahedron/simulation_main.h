#pragma once
#include "./external/glad.h"
#include "./external/glfw3.h"
#include "./basic/opengl_input.h"
#include"basic/basic_imgui.h"


void switchObjMoveMode(Input* input, bool* control_parameter);
bool couldMoveCamera(Input* input, bool* control_parameter);
void setHideWireframe(std::vector<std::vector<bool>>& show_element,
	int collider_num, int cloth_num, int tetrahedron_num);

void simu_main(GLFWwindow* window, Input* input);



