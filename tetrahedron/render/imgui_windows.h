#pragma once
#include<vector>
#include"../external/imgui/imfilebrowser.h"
#include"../global_struct.h"

class ImGuiWindows
{
public:
	void controlWindow(bool* control_parameter, float* force_coe);
	void visualizationControlPanel(bool& reset_camera, std::vector<std::vector<bool>>& show_element,
		bool only_collision_test, bool* control_parameter);
	bool loadModel(std::string& scene_path, std::vector<std::string>& collider_path, std::vector<std::string>& object_path);
	void infoWindow(std::vector<std::array<int, 3>>& cloth_info, std::vector<double>& cloth_mass,
		std::vector<std::array<int, 3>>& tetrahedron_info, std::vector<double>& tetrahedron_mass,
		double time, int* iteration_num, int* set_itr_num, double* convergence_rate, int time_stamp, bool& start_edit, bool& start_simulation);
	void operationWindow(std::vector<std::array<double, 6>>& cloth_stiffness, std::vector<std::array<double, 8>>& tet_stiffness,
		double* simulation_parameters,
		std::vector<std::array<double, 8>>& cloth_collision_stiffness,
		std::vector<std::array<double, 8>>& tet_collision_stiffness,
		bool* set_stiffness, double* temp_stiffness, UpdateObjStiffness& update_obj_stiffness, bool* set_anchor_point, bool tetrahedron_exist,
		double* damp_stiffness, double* rayleigh_damp_stiffness);
	void iterationSolverInfoWindow(double& solver_iteration_num, int& itr_sovler_method,
		char** itr_solver_items, char*& itr_solver_item, int item_num, double* conv_rate, bool& record_matrix);

	bool open_load_model;

	void helpMarker(const char* desc);
	ImGui::FileBrowser m_file_dialog_info;
	int load_collider_is_open = -2; //-1:scene, 0:collider, 1:object
	bool load_simulation_is_open = false;
	//bool load_scene_is_open = false;
	
	void floorInfo(bool& exist, bool& show, bool& normal_direction, unsigned int& dimension, double* value, bool& eidit);

	bool loadScene(std::string& path);


};

