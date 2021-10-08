#pragma once
#include<vector>
#include"../external/imgui/imfilebrowser.h"

class ImGuiWindows
{
public:
	void controlWindow(bool* control_parameter, float* force_coe);
	void visualizationControlPanel(bool& reset_camera, std::vector<std::vector<bool>>& wireframe, std::vector<std::vector<bool>>& hide);
	bool loadModel(std::vector<std::string>& collider_path, std::vector<std::string>& tetrohedron_path);

private:
	bool open_load_model;
	void helpMarker(const char* desc);
	ImGui::FileBrowser m_file_dialog_info;
	bool load_collider_is_open = false;
};

