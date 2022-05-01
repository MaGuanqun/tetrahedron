#pragma once
#include<vector>
#include"../external/imgui/imfilebrowser.h"
#include"../global_struct.h"

class ImGuiWindows
{
public:
	void controlWindow(bool* control_parameter, float* force_coe);
	void visualizationControlPanel(bool& reset_camera, std::vector<std::vector<bool>>& wireframe, std::vector<std::vector<bool>>& hide);
	bool loadModel(std::vector<std::string>& collider_path, std::vector<std::string>& object_path);
	void infoWindow(std::vector<std::array<int, 3>>& cloth_info, std::vector<double>& cloth_mass,
		std::vector<std::array<int, 3>>& tetrahedron_info, std::vector<double>& tetrahedron_mass,
		double time, int* iteration_num, int* set_itr_num, double* convergence_rate, int time_stamp, bool& start_edit, bool& start_simulation);
	void operationWindow(std::vector<std::array<double, 3>>& cloth_stiffness, std::vector<std::array<double, 2>>& tet_stiffness,
		double* simulation_parameters,
		std::vector<std::array<double, 4>>& cloth_collision_stiffness,
		std::vector<std::array<double, 4>>& tet_collision_stiffness,
		bool* set_stiffness, double* temp_stiffness, UpdateObjStiffness& update_obj_stiffness, bool* set_anchor_point, bool tetrahedron_exist);
	void iterationSolverInfoWindow(double& solver_iteration_num, int& itr_sovler_method,
		char** itr_solver_items, char*& itr_solver_item, int item_num, double* conv_rate, bool& record_matrix);

	bool open_load_model;
	void helpMarker(const char* desc);
	ImGui::FileBrowser m_file_dialog_info;
	bool load_collider_is_open = false;
};

