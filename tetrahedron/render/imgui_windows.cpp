#include"../external/imgui/imgui.h"
#include"../external/imgui/imgui_impl_glfw.h"
#include"../external/imgui/imgui_impl_opengl3.h"
#include"imgui_windows.h"
#include"../basic/global.h"
#include"../basic/enum_setting.h"

void ImGuiWindows::controlWindow(bool* control_parameter, float* force_coe)
{
	ImGui::SetNextWindowPos(ImVec2(SCR_WIDTH - 270, 70));
	ImGui::SetNextWindowSize(ImVec2(270, 850));
	ImGui::Begin("Control Panel");
	ImGui::SetNextItemOpen(true);
	if (ImGui::TreeNode("State Control")) {
		if (control_parameter[USE_XPBD] || control_parameter[USE_PD_] || control_parameter[USE_NEWTON_]) {
			if (ImGui::Button("1 Frame", ImVec2(160, 25)))
			{
				control_parameter[ONE_FRAME] = true;
				control_parameter[MOVE_OBJ] = false;
				control_parameter[ONLY_MOVE_CURRENT_POSITION] = false;
				control_parameter[ROTATION] = false;
				control_parameter[ONLY_ROTATE_CURRENT] = false;
			}
			if (control_parameter[START_SIMULATION]) {
				if (ImGui::Button("Pause Simulation", ImVec2(160, 25)))
				{
					control_parameter[START_SIMULATION] = false;
				}
			}
			else {
				if (ImGui::Button("Continue Simulation", ImVec2(160, 25)))
				{
					control_parameter[START_SIMULATION] = true;
					control_parameter[MOVE_OBJ] = false;
					control_parameter[ONLY_MOVE_CURRENT_POSITION] = false;
					control_parameter[ROTATION] = false;
					control_parameter[ONLY_ROTATE_CURRENT] = false;
				}
			}
			if (ImGui::Button("Reset Simulation", ImVec2(160, 25))) {
				control_parameter[RESET_SIMULATION] = true;
			}
			ImGui::SameLine();
			helpMarker(
				"Reset with current stiffness");
		}
		else {
			if (ImGui::Button("Start Detection", ImVec2(160, 25)))
			{
				control_parameter[ONE_FRAME] = true;
				control_parameter[MOVE_OBJ] = false;
				control_parameter[ONLY_MOVE_CURRENT_POSITION] = false;
			}
		}
		if (ImGui::Button("Initial Simulation", ImVec2(160, 25))) {
			control_parameter[INITIAL_SIMULATION] = true;
		}
		ImGui::SameLine();
		helpMarker(
			"Reset with default stiffness");
		ImGui::TreePop();
	}
	ImGui::SetNextItemOpen(true);
	if (ImGui::TreeNode("Simulation control")) {
	//if (!control_parameter[ONLY_COLLISION_TEST]) {
			ImGui::TextWrapped("Pause simulaition rightly after adding a force and relasing the mouse.");
			if (control_parameter[STOP_AFTER_RELEASE_MOUSE]) {
				if (ImGui::Button("Off", ImVec2(160, 25)))
				{
					control_parameter[STOP_AFTER_RELEASE_MOUSE] = false;
				}
			}
			else {
				if (ImGui::Button("On", ImVec2(160, 25)))
				{
					control_parameter[STOP_AFTER_RELEASE_MOUSE] = true;
				}
			}
		//}
		if (control_parameter[OUTPUT_IMAGE]) {
			ImGui::TextWrapped("Press end to stop recording screen.");
			if (ImGui::Button("End##record_screen", ImVec2(160, 25)))
			{
				control_parameter[OUTPUT_IMAGE] = false;
			}
		}
		else {
			ImGui::TextWrapped("Press start to start recording screen.");
			if (ImGui::Button("Start##recording screen", ImVec2(160, 25)))
			{
				control_parameter[OUTPUT_IMAGE] = true;
				control_parameter[START_SIMULATION] = true;
			}
		}
		if (control_parameter[SAVE_OBJ]) {
			ImGui::TextWrapped("Press end to stop saving obj file.");
			if (ImGui::Button("End##record_obj", ImVec2(160, 25)))
			{
				control_parameter[SAVE_OBJ] = false;
			}
		}
		else {
			ImGui::TextWrapped("Press start to start saving obj file.");
			if (ImGui::Button("Start##record_obj", ImVec2(160, 25)))
			{
				control_parameter[SAVE_OBJ] = true;
			}
		}
		ImGui::TreePop();
	}
		ImGui::SetNextItemOpen(true);
		if (ImGui::TreeNode("Set cursor force")) {
			if (!control_parameter[SET_CURSOR_FORCE]) {
				if (ImGui::Button("Set force", ImVec2(160, 25)))
				{
					control_parameter[SET_CURSOR_FORCE] = true;
				}
			}
			ImGui::TreePop();
		}
		//if (!control_parameter[ONLY_COLLISION_TEST]) {
		//ImGui::SetNextItemOpen(true);
		//if (ImGui::TreeNode("Start Test")) {
		//	if (ImGui::Button("Start test", ImVec2(160, 25)))
		//	{
		//		control_parameter[START_TEST] = true;
		//		control_parameter[INITIAL_TEST] = true;
		//	}
		//	ImGui::TreePop();
		//}
		//if (ImGui::Button("Move obj test", ImVec2(160, 25)))
		//{
		//	control_parameter[MOVE_OBJ_SCRIPT] = true;
		//}
		//}
		if (control_parameter[ONLY_COLLISION_TEST]) {
			if (control_parameter[DRAW_VT]) {
				ImGui::Text("Show VT collision pair");
				if (ImGui::Button("Switch to EE", ImVec2(160, 25)))
				{
					control_parameter[DRAW_VT] = false;
				}
			}
			else {
				ImGui::Text("Show EE collision pair");
				if (ImGui::Button("Switch to VT", ImVec2(160, 25)))
				{
					control_parameter[DRAW_VT] = true;
				}
			}
		}
		ImGui::Text("Move Step 1:");
		if (control_parameter[MOVE_OBJ]) {
			if (ImGui::Button("Stop move object", ImVec2(160, 25)))
			{
				control_parameter[MOVE_OBJ] = false;
			}
		}
		else {
			if (ImGui::Button("Start move object", ImVec2(160, 25)))
			{
				control_parameter[MOVE_OBJ] = true;
				control_parameter[ONLY_MOVE_CURRENT_POSITION] = false;
				control_parameter[ROTATION] = false;
				control_parameter[ONLY_ROTATE_CURRENT] = false;
			}
		}
		if (control_parameter[ROTATION]) {
			if (ImGui::Button("Stop rotate object", ImVec2(160, 25)))
			{
				control_parameter[ROTATION] = false;
			}
		}
		else {
			if (ImGui::Button("Start rotate object", ImVec2(160, 25)))
			{
				control_parameter[ROTATION] = true;
				control_parameter[ONLY_MOVE_CURRENT_POSITION] = false;
				control_parameter[MOVE_OBJ] = false;
				control_parameter[ONLY_ROTATE_CURRENT] = false;
			}
		}
		ImGui::Text("Move Step 2:");
		if (control_parameter[ONLY_MOVE_CURRENT_POSITION]) {
			if (ImGui::Button("Stop move current pos", ImVec2(160, 25)))
			{
				control_parameter[ONLY_MOVE_CURRENT_POSITION] = false;
			}
		}
		else {
			if (ImGui::Button("Start move current pos", ImVec2(160, 25)))
			{
				control_parameter[ONLY_MOVE_CURRENT_POSITION] = true;
				control_parameter[MOVE_OBJ] = false;
			}
		}
		if (control_parameter[ONLY_ROTATE_CURRENT]) {
			if (ImGui::Button("Stop rotate current pos", ImVec2(160, 25)))
			{
				control_parameter[ONLY_ROTATE_CURRENT] = false;
			}
		}
		else {
			if (ImGui::Button("Start rotate current pos", ImVec2(160, 25)))
			{
				control_parameter[ONLY_ROTATE_CURRENT] = true;
				control_parameter[MOVE_OBJ] = false;
				control_parameter[ONLY_MOVE_CURRENT_POSITION] = false;
				control_parameter[ROTATION] = false;
			}
		}
		if (control_parameter[ONLY_COLLISION_TEST]) {
			if (control_parameter[DRAW_SPATIAL_HASHING]) {
				if (ImGui::Button("Hide Spatial Hashing Cell", ImVec2(160, 25)))
				{
					control_parameter[DRAW_SPATIAL_HASHING] = false;
				}
			}
			else {
				if (ImGui::Button("Show Spatial Hashing Cell", ImVec2(160, 25)))
				{
					control_parameter[DRAW_SPATIAL_HASHING] = true;
					control_parameter[SPATIAL_HASHING_UPDATE] = true;
				}
			}
			if (!control_parameter[SHORTCUT_INSTRUCTION]) {
				if (ImGui::Button("Shortcut instuct.", ImVec2(160, 25)))
				{
					control_parameter[SHORTCUT_INSTRUCTION] = true;
				}
			}
			if (control_parameter[DRAW_SPATIAL_HASHING]) {
				ImGui::Text("Draw Collision in a cell");
				if (ImGui::Button("<<front.", ImVec2(120, 25))) {
					control_parameter[SEARCH_LEFT_SH_CELL] = true;
				}
				ImGui::SameLine();
				if (ImGui::Button(">>back.", ImVec2(120, 25))) {
					control_parameter[SEARCH_RIGHT_SH_CELL] = true;
				}
				if (control_parameter[DRAW_ALL_PAIRS_IN_A_CELL]) {
					if (ImGui::Button("Only show collision", ImVec2(120, 25))) {
						control_parameter[DRAW_ALL_PAIRS_IN_A_CELL] = false;
					}
				}
				else {
					if (ImGui::Button("Show all primitive", ImVec2(120, 25))) {
						control_parameter[DRAW_ALL_PAIRS_IN_A_CELL] = true;
					}
				}
			}

		}
	ImGui::End();
		if (control_parameter[SET_CURSOR_FORCE]) {
			ImGui::SetNextWindowSize(ImVec2(330, 240));
			ImGui::Begin("Set Cursor Force", &control_parameter[SET_CURSOR_FORCE]);
			ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.9f);
			ImGui::SetNextItemOpen(true);
			if (ImGui::TreeNode("Set force coefficient")) {
				ImGui::SliderFloat("##set force coefficient", force_coe, 0.01, 1.0, "force coefficient = %.2f");
				ImGui::TreePop();
			}
			ImGui::TextWrapped("Force depends on force coefficient and cursor moving speed.");
			ImGui::End();
		}



	if (control_parameter[SHORTCUT_INSTRUCTION]) {
		ImGui::SetNextWindowSize(ImVec2(330, 400));
		ImGui::Begin("Shortcut Instruction", &control_parameter[SHORTCUT_INSTRUCTION]);
		ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.9f);
		ImGui::SetNextItemOpen(true);
		ImGui::TextWrapped("Press \"S\" to start detection");
		ImGui::TextWrapped("Press \"M\" to move camera in rotation and translation mode");
		ImGui::TextWrapped("Press \"R\" to start rotation mode");
		ImGui::TextWrapped("Press \"T\" to start translation mode");
		ImGui::TextWrapped("Press \"Ctrl + T\" to start translation mode (only move current position)");
		ImGui::TextWrapped("Press \"Ctrl + R\" to start rotation mode (only move current position)");
		ImGui::End();
	}
}


void ImGuiWindows::visualizationControlPanel(bool& reset_camera, std::vector<std::vector<bool>>& show_element,
	bool only_collision_test, bool* control_parameter)
{
	ImGui::SetNextWindowPos(ImVec2(0, 710));
	ImGui::SetNextWindowSize(ImVec2(240, 370));
	ImGui::Begin("Visualization");
	ImGui::SetNextItemOpen(true);
	//	if (ImGui::TreeNode("Visualization")) {
	int cloth_no = -1;

	if (control_parameter[SHARP_EDGE_SHADING]) {
		if (ImGui::Button("Swith to soft edge", ImVec2(160, 25)))
		{
			control_parameter[SHARP_EDGE_SHADING] = false;
		}
	}
	else {
		if (ImGui::Button("Swith to sharp edge", ImVec2(160, 25)))
		{
			control_parameter[SHARP_EDGE_SHADING] = true;
		}
	}

	if (ImGui::Button("Reset Camera", ImVec2(160, 25)))
	{
		reset_camera = true;
	}
	std::string tempString;
	for (int i = 0; i < show_element[COLLIDER_].size(); ++i) {
		ImGui::SetNextItemOpen(true);
		tempString = "Collider " + std::to_string(i);
		if (ImGui::TreeNode(tempString.c_str())) {
			if (show_element[COLLIDER_][i]) {
				tempString = "Show##Collider" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[COLLIDER_][i] = false;
				}
			}
			else {
				tempString = "Hide##Collider" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[COLLIDER_][i] = true;
				}
			}
			if (show_element[3 + COLLIDER_][i]) {
				tempString = "Hide WireFrame##Collider" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[3 + COLLIDER_][i] = false;
				}
			}
			else {
				tempString = "Show WireFrame##Collider" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[3 + COLLIDER_][i] = true;
				}
			}
			if (only_collision_test) {
				if (show_element[6 + COLLIDER_][i]) {
					tempString = "Hide Collision##Collider" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[6 + COLLIDER_][i] = false;
					}
				}
				else {
					tempString = "Show Collision##Collider" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[6 + COLLIDER_][i] = true;
					}
				}
				if (show_element[9 + COLLIDER_][i]) {
					tempString = "Hide Ori Pos##Collider" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[9 + COLLIDER_][i] = false;
					}
				}
				else {
					tempString = "Show Ori Pos##Collider" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[9 + COLLIDER_][i] = true;
					}
				}
				if (show_element[12 + COLLIDER_][i]) {
					tempString = "Hide Ori Wireframe##Collider" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[12 + COLLIDER_][i] = false;
					}
				}
				else {
					tempString = "Show Ori Wireframe##Collider" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[12 + COLLIDER_][i] = true;
					}
				}
			}
			ImGui::TreePop();
		}
	}
	for (int i = 0; i < show_element[CLOTH_].size(); ++i) {
		ImGui::SetNextItemOpen(true);
		tempString = "cloth " + std::to_string(i);
		if (ImGui::TreeNode(tempString.c_str())) {
			if (show_element[CLOTH_][i]) {
				tempString = "Show##cloth" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[CLOTH_][i] = false;
				}
			}
			else {
				tempString = "Hide##cloth" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[CLOTH_][i] = true;
				}
			}
			if (show_element[3 + CLOTH_][i]) {
				tempString = "Hide WireFrame##cloth" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[3 + CLOTH_][i] = false;
				}
			}
			else {
				tempString = "Show WireFrame##cloth" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[3 + CLOTH_][i] = true;
				}
			}
			if (only_collision_test) {
				if (show_element[6 + CLOTH_][i]) {
					tempString = "Hide Collision##cloth" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[6 + CLOTH_][i] = false;
					}
				}
				else {
					tempString = "Show Collision##cloth" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[6 + CLOTH_][i] = true;
					}
				}
				if (show_element[9 + CLOTH_][i]) {
					tempString = "Hide Ori Pos##cloth" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[9 + CLOTH_][i] = false;
					}
				}
				else {
					tempString = "Show Ori Pos##cloth" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[9 + CLOTH_][i] = true;
					}
				}
				if (show_element[12 + CLOTH_][i]) {
					tempString = "Hide Ori Wireframe##cloth" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[12 + CLOTH_][i] = false;
					}
				}
				else {
					tempString = "Show Ori Wireframe##cloth" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[12 + CLOTH_][i] = true;
					}
				}
			}
			ImGui::TreePop();
		}
	}
	for (int i = 0; i < show_element[TETRAHEDRON_].size(); ++i) {
		ImGui::SetNextItemOpen(true);
		tempString = "tet " + std::to_string(i);
		if (ImGui::TreeNode(tempString.c_str())) {
			if (show_element[TETRAHEDRON_][i]) {
				tempString = "Show##tetrahedron" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[TETRAHEDRON_][i] = false;
				}
			}
			else {
				tempString = "Hide##tetrahedron" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[TETRAHEDRON_][i] = true;
				}
			}
			if (show_element[3 + TETRAHEDRON_][i]) {
				tempString = "Hide WireFrame##tetrahedron" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[3 + TETRAHEDRON_][i] = false;
				}
			}
			else {
				tempString = "Show WireFrame##tetrahedron" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					show_element[3 + TETRAHEDRON_][i] = true;
				}
			}
			if (only_collision_test) {
				if (show_element[6 + TETRAHEDRON_][i]) {
					tempString = "Hide Collision##tetrahedron" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[6 + TETRAHEDRON_][i] = false;
					}
				}
				else {
					tempString = "Show Collision##tetrahedron" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[6 + TETRAHEDRON_][i] = true;
					}
				}
				if (show_element[9 + TETRAHEDRON_][i]) {
					tempString = "Hide Ori Pos##tetrahedron" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[9 + TETRAHEDRON_][i] = false;
					}
				}
				else {
					tempString = "Show Ori Pos##tetrahedron" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[9 + TETRAHEDRON_][i] = true;
					}
				}
				if (show_element[12 + TETRAHEDRON_][i]) {
					tempString = "Hide Ori Wireframe##tetrahedron" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[12 + TETRAHEDRON_][i] = false;
					}
				}
				else {
					tempString = "Show Ori Wireframe##tetrahedron" + std::to_string(i);
					if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
					{
						show_element[12 + TETRAHEDRON_][i] = true;
					}
				}
			}
			ImGui::TreePop();
		}
	}
	ImGui::End();
}


bool ImGuiWindows::loadModel(std::vector<std::string>& collider_path, std::vector<std::string>& object_path)
{
	bool finished_loading = false;
	ImGui::SetNextWindowPos(ImVec2(SCR_WIDTH - 270, 0));
	ImGui::SetNextWindowSize(ImVec2(270, 70));
	ImGui::Begin("Load Models");
	if (ImGui::Button("Load Model", ImVec2(160, 25))) {
		open_load_model = true;
	}
	ImGui::End();
	if (open_load_model) {
		ImGui::SetNextWindowSize(ImVec2(700, 160));
		ImGui::Begin("Load Model", &open_load_model);
		ImGui::TextWrapped("Load colliders and objects.");
		if (ImGui::Button("Load Collider", ImVec2(160, 25))) {
			m_file_dialog_info.SetTypeFilters({ ".obj" });
			m_file_dialog_info.SetTitle("Load Collider");
			m_file_dialog_info.Open();
			load_collider_is_open = true;
		}
		if (load_collider_is_open) {
			m_file_dialog_info.Display();
			if (m_file_dialog_info.HasSelected()) {
				collider_path.push_back(m_file_dialog_info.GetSelected().string());
				m_file_dialog_info.ClearSelected();
			}
		}
		if (!collider_path.empty()) {
			ImGui::TextWrapped("Collider path:");
			for (int i = 0; i < collider_path.size(); ++i) {
				ImGui::TextWrapped(collider_path[i].c_str());
			}
		}
		if (ImGui::Button("Load object", ImVec2(160, 25))) {
			m_file_dialog_info.SetTypeFilters({ ".obj",".ele" });
			m_file_dialog_info.SetTitle("Load object");
			m_file_dialog_info.Open();
			load_collider_is_open = false;
		}
		if (!load_collider_is_open) {
			m_file_dialog_info.Display();
			if (m_file_dialog_info.HasSelected()) {
				//if (!object_path.empty() && object_path[object_path.size() - 1] == m_file_dialog_info.GetSelected().string()) {
				//	ImGui::TextWrapped("Error, Need to add different path");
				//}
				//else {
				object_path.push_back(m_file_dialog_info.GetSelected().string());
				m_file_dialog_info.ClearSelected();
				//}		
			}
		}
		if (!object_path.empty()) {
			ImGui::TextWrapped("object path:");
			for (int i = 0; i < object_path.size(); ++i) {
				ImGui::TextWrapped(object_path[i].c_str());
			}
		}
		if (ImGui::Button("Cancel", ImVec2(160, 25))) {
			object_path.clear();
			collider_path.clear();
		}
		ImGui::SameLine();
		if (ImGui::Button("Confirm", ImVec2(160, 25))) {
			finished_loading = true;
		}
		ImGui::End();
	}
	if (finished_loading) {
		return true;
	}
	return false;
}

void ImGuiWindows::infoWindow(std::vector<std::array<int, 3>>& cloth_info, std::vector<double>& cloth_mass,
	std::vector<std::array<int, 3>>& tetrahedron_info, std::vector<double>& tetrahedron_mass,
	double time, int* iteration_num, int* set_itr_num, double* convergence_rate, int time_stamp, bool& start_edit, bool& start_simulation)
{
	ImGui::SetNextWindowPos(ImVec2(0, 0));
	ImGui::SetNextWindowSize(ImVec2(280, 330));
	ImGui::Begin("Basic Information This Time step");
	ImGui::Text("stamp: %i", time_stamp);
	ImGui::Text("Time(ms): %f", time);
	ImGui::SetNextItemOpen(true);
	if (ImGui::TreeNode("Iteration")) {
		if (!start_edit) {
			if (ImGui::Button("Edit Itr", ImVec2(160, 25))) {
				start_edit = true;
				start_simulation = false;
			}
		}
		else {
			if (ImGui::Button("Confirm Itr", ImVec2(160, 25))) {
				start_edit = false;
			}
		}
		if (start_edit) {
			ImGui::Text("per substep iteration:");
			ImGui::InputInt("##per substep iteration", &set_itr_num[LOCAL_GLOBAL]);
			ImGui::Text("substep num:");
			ImGui::InputInt("##substep num", &set_itr_num[OUTER]);
			ImGui::Text("Local-global convergence rate:");
			ImGui::InputDouble("##Local-global convergence rate", &convergence_rate[LOCAL_GLOBAL], 0.0f, 0.0f, "%.2e");
			ImGui::Text("Outter convergence rate:");
			ImGui::InputDouble("##Outter convergence rate:", &convergence_rate[OUTER], 0.0f, 0.0f, "%.2e");
		}
		else {
			ImGui::Text("Local-global iteration: %i", iteration_num[LOCAL_GLOBAL]);
			ImGui::Text("Outter iteration: %i", iteration_num[OUTER]);
			ImGui::TextWrapped("Local-global convergence rate: %.2e", convergence_rate[LOCAL_GLOBAL]);
			ImGui::TextWrapped("Outter convergence rate:\n %.2e", convergence_rate[OUTER]);
		}
		ImGui::TreePop();
	}
	std::string tempString;
	for (int i = 0; i < cloth_info.size(); ++i) {
		ImGui::SetNextItemOpen(true);
		tempString = "Cloth " + std::to_string(i);
		if (ImGui::TreeNode(tempString.c_str()))
		{
			ImGui::Text("Vertex: %i", cloth_info[i][0]);
			ImGui::Text("Edge: %i", cloth_info[i][1]);
			ImGui::Text("Face: %i", cloth_info[i][2]);
			ImGui::Text("Total Mass: %.2f", cloth_mass[i]);
			ImGui::TreePop();
		}
	}
	for (int i = 0; i < tetrahedron_info.size(); ++i) {
		ImGui::SetNextItemOpen(true);
		tempString = "Tetrahedron " + std::to_string(i);
		if (ImGui::TreeNode(tempString.c_str()))
		{
			ImGui::Text("Vertex: %i", tetrahedron_info[i][0]);
			ImGui::Text("Tetrahedron: %i", tetrahedron_info[i][1]);
			ImGui::Text("Surface Triangle: %i", tetrahedron_info[i][2]);
			ImGui::Text("Total Mass: %.2f", tetrahedron_mass[i]);
			ImGui::TreePop();
		}
	}
	ImGui::End();
}

void ImGuiWindows::helpMarker(const char* desc)
{
	ImGui::TextDisabled("(?)");
	if (ImGui::IsItemHovered())
	{
		ImGui::BeginTooltip();
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
		ImGui::TextUnformatted(desc);
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}
}



void ImGuiWindows::operationWindow(std::vector<std::array<double, 3>>& cloth_stiffness, std::vector<std::array<double, 3>>& tet_stiffness, double* simulation_parameters, std::vector<std::array<double, 4>>& cloth_collision_stiffness,
	std::vector<std::array<double, 4>>& tet_collision_stiffness,
	bool* set_stiffness, double* temp_stiffness, UpdateObjStiffness& update_obj_stiffness, bool* set_anchor_point, bool tetrahedron_exist,
	double* damp_stiffness, double* rayleigh_damp_stiffness)
{
	ImGui::SetNextWindowPos(ImVec2(0, 330));
	ImGui::SetNextWindowSize(ImVec2(240, 140));
	ImGui::Begin("Simulation");
	ImGui::Text("Time Step: %.3f", simulation_parameters[0]);
	ImGui::Text("Gravity: %f", simulation_parameters[1]);
	//ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.35f, 0.75f, 0.35f, 1.0f));
	if (ImGui::Button("Set Stiffness", ImVec2(160, 25))) {
		set_stiffness[START_SETTING] = true;
	}
	if (tetrahedron_exist) {
		if (!set_anchor_point[0]) {
			if (ImGui::Button("Anchor Vertex", ImVec2(160, 25))) {
				set_anchor_point[0] = true;
			}
		}
		else {
			if (ImGui::Button("Save Anchor", ImVec2(160, 25))) {
				set_anchor_point[1] = true;
				set_anchor_point[0] = false;
			}
		}
	}
	//ImGui::PopStyleColor();
	ImGui::End();
	bool test = true;
	if (set_stiffness[START_SETTING]) {
		ImGui::SetNextWindowSize(ImVec2(500, 330));
		ImGui::Begin("Set Stiffness", set_stiffness + START_SETTING);
		ImGui::Text("Press to set constraint stiffness.");
		if (ImGui::Button("set stiffness", ImVec2(160, 25))) {
			set_stiffness[EDIT] = true;
		}
		ImGui::Text("Default stiffness:");
		std::string tempString;
		for (int i = 0; i < cloth_stiffness.size(); ++i) {
			ImGui::SetNextItemOpen(true);
			tempString = "Cloth " + std::to_string(i) + "##stiffness";
			if (ImGui::TreeNode(tempString.c_str())) {
				ImGui::Text("Stretch: %.2e", cloth_stiffness[i][0]);
				ImGui::SameLine();
				ImGui::Text("Bending: %.2e", cloth_stiffness[i][1]);
				ImGui::SameLine();
				ImGui::Text("Position: %.2e", cloth_stiffness[i][2]);
				ImGui::Text("Collision stiffness:(P point, E edge, T triangle)");
				ImGui::Text("Self-collision:");
				ImGui::Text("PT: %.2e", cloth_collision_stiffness[i][SELF_POINT_TRIANGLE]);
				ImGui::SameLine();
				ImGui::Text("EE: %.2e", cloth_collision_stiffness[i][SELF_EDGE_EDGE]);
				ImGui::SameLine();
				ImGui::Text("PP: %.2e", cloth_collision_stiffness[i][SELF_POINT_POINT]);
				ImGui::Text("body-cloth collision:");
				ImGui::Text("PT: %.2e", cloth_collision_stiffness[i][BODY_POINT_TRIANGLE]);
				ImGui::TreePop();
			}
		}
		for (int i = 0; i < tet_stiffness.size(); ++i) {
			ImGui::SetNextItemOpen(true);
			tempString = "Tet " + std::to_string(i) + "##stiffness";
			if (ImGui::TreeNode(tempString.c_str())) {
				ImGui::Text("ARAP: %.2e", tet_stiffness[i][0]);
				ImGui::Text("Position: %.2e", tet_stiffness[i][1]);
				ImGui::Text("Edge Length: %.2e", tet_stiffness[i][2]);
				ImGui::Text("Collision stiffness:(P point, E edge, T triangle)");
				ImGui::Text("Self-collision:");
				ImGui::Text("PT: %.2e", tet_collision_stiffness[i][SELF_POINT_TRIANGLE]);
				ImGui::SameLine();
				ImGui::Text("EE: %.2e", tet_collision_stiffness[i][SELF_EDGE_EDGE]);
				ImGui::SameLine();
				ImGui::Text("PP: %.2e", tet_collision_stiffness[i][SELF_POINT_POINT]);
				ImGui::Text("body-cloth collision:");
				ImGui::Text("PT: %.2e", tet_collision_stiffness[i][BODY_POINT_TRIANGLE]);
				ImGui::TreePop();
			}
		}

		ImGui::End();
	}
	if (set_stiffness[EDIT]) {
		ImGui::SetNextWindowSize(ImVec2(300, 330));
		ImGui::Begin("Set Stiffness##1", set_stiffness + EDIT);
		ImGui::PushItemWidth(ImGui::GetFontSize() * 7.0f);
		ImGui::Text("Damp stiffness: ");
		if (ImGui::InputDouble("##dampStiff", &temp_stiffness[DAMP_STIFFNESS], 0.0f, 0.0f, "%.4f")) {
			set_stiffness[EDIT_DAMP_STIFFNESS] = false;
		}
		ImGui::SameLine();
		if (!set_stiffness[EDIT_DAMP_STIFFNESS]) {
			if (ImGui::Button("Save##dampStiff", ImVec2(80, 25))) {
				set_stiffness[EDIT_DAMP_STIFFNESS] = true;
				*damp_stiffness = temp_stiffness[DAMP_STIFFNESS];
			}
		}
		else {
			ImGui::Text("Saved");
		}
		ImGui::Text("Rayleigh Damp stiffness: ");
		ImGui::Text("alpha: ");
		ImGui::SameLine();
		if (ImGui::InputDouble("##rayleighDampStiffAlpha", &temp_stiffness[RAYLEIGH_DAMP_STIFFNESS_ALPHA], 0.0f, 0.0f, "%.4f")) {
			set_stiffness[EDIT_DAMP_STIFFNESS] = false;
		}
		ImGui::Text("beta: ");
		ImGui::SameLine();
		if (ImGui::InputDouble("##rayleighDampStiffBeta", &temp_stiffness[RAYLEIGH_DAMP_STIFFNESS_BETA], 0.0f, 0.0f, "%.4f")) {
			set_stiffness[EDIT_DAMP_STIFFNESS] = false;
		}
		if (!set_stiffness[EDIT_DAMP_STIFFNESS]) {
			if (ImGui::Button("Save##RayleidampStiff", ImVec2(80, 25))) {
				set_stiffness[EDIT_DAMP_STIFFNESS] = true;
				rayleigh_damp_stiffness[0] = temp_stiffness[RAYLEIGH_DAMP_STIFFNESS_ALPHA];
				rayleigh_damp_stiffness[1] = temp_stiffness[RAYLEIGH_DAMP_STIFFNESS_BETA];
			}
		}
		else {
			ImGui::Text("Saved");
		}
		ImGui::Text("Length: ");
		//ImGui::SameLine();
		if (ImGui::InputDouble("##Length", &temp_stiffness[LENGTH], 0.0f, 0.0f, "%.2f")) {
			set_stiffness[EDIT_LENGTH] = false;
		}
		ImGui::SameLine();
		if (!set_stiffness[EDIT_LENGTH]) {
			if (ImGui::Button("Save##length", ImVec2(80, 25))) {
				update_obj_stiffness.update_length = true;
				update_obj_stiffness.length_stiffness = temp_stiffness[LENGTH];
				set_stiffness[EDIT_LENGTH] = true;
			}
		}
		else {
			ImGui::Text("Saved");
		}
		ImGui::Text("Bending: ");
		if (ImGui::InputDouble("##Bending", &temp_stiffness[BENDING], 0.0f, 0.0f, "%.2e")) {
			set_stiffness[EDIT_BENDING] = false;
		}
		ImGui::SameLine();
		helpMarker(
			"You can input bending using the scientific notation");
		ImGui::SameLine();
		if (!set_stiffness[EDIT_BENDING]) {
			if (ImGui::Button("Save##bend", ImVec2(80, 25))) {
				set_stiffness[EDIT_BENDING] = true;
				update_obj_stiffness.update_bend = true;
				update_obj_stiffness.bend_stiffness = temp_stiffness[BENDING];
			}
		}
		else {
			ImGui::Text("Saved");
		}
		ImGui::Text("ARAP: ");
		if (ImGui::InputDouble("##ARAP", &temp_stiffness[ARAP], 0.0f, 0.0f, "%.2e")) {
			set_stiffness[EDIT_ARAP] = false;
		}
		ImGui::SameLine();
		helpMarker(
			"You can input ARAP using the scientific notation");
		ImGui::SameLine();
		if (!set_stiffness[EDIT_ARAP]) {
			if (ImGui::Button("Save##ARAP", ImVec2(80, 25))) {
				set_stiffness[EDIT_ARAP] = true;
				update_obj_stiffness.update_ARAP = true;
				update_obj_stiffness.ARAP_stiffness = temp_stiffness[ARAP];
			}
		}
		else {
			ImGui::Text("Saved");
		}
		ImGui::Text("Tet Edge Length: ");
		if (ImGui::InputDouble("##TetEdgeLength", &temp_stiffness[TET_EDGE_LENGTH], 0.0f, 0.0f, "%.2e")) {
			set_stiffness[EDIT_TET_EDGE_LENGTH] = false;
		}
		ImGui::SameLine();
		if (!set_stiffness[EDIT_TET_EDGE_LENGTH]) {
			if (ImGui::Button("Save##TetEdgeLength", ImVec2(80, 25))) {
				set_stiffness[EDIT_TET_EDGE_LENGTH] = true;
				update_obj_stiffness.update_tet_edge_length = true;
				update_obj_stiffness.tet_edge_length_stiffness = temp_stiffness[TET_EDGE_LENGTH];
			}
		}
		else {
			ImGui::Text("Saved");
		}
		if (ImGui::Button("Set Collision", ImVec2(160, 25))) {
			set_stiffness[EDIT_COLLISION] = true;
		}
		if (ImGui::Button("Confirm Stiffness", ImVec2(160, 25))) {
			set_stiffness[STIFFNESS_CONFIRMED] = true;
			set_stiffness[EDIT] = false;
		}
		ImGui::TextWrapped("Press the button above to save all change.");
		ImGui::End();
	}
	else {
		memset(set_stiffness + EDIT_LENGTH, 0, 2);
	}
	if (set_stiffness[EDIT_COLLISION]) {
		ImGui::SetNextWindowSize(ImVec2(270, 230));
		ImGui::Begin("Set Collision", set_stiffness + EDIT_COLLISION);
		ImGui::PushItemWidth(ImGui::GetFontSize() * 7.0f);
		ImGui::Text("Cloth point triangle collision: ");
		if (ImGui::InputDouble("##CPTcollision", &temp_stiffness[SELF_POINT_TRIANGLE], 0.0f, 0.0f, "%.2f")) {
			set_stiffness[EDIT_SELF_POINT_TRIANGLE] = false;
		}
		ImGui::SameLine();
		if (!set_stiffness[EDIT_SELF_POINT_TRIANGLE]) {
			if (ImGui::Button("Save##CPTcollision", ImVec2(80, 25))) {
				set_stiffness[EDIT_SELF_POINT_TRIANGLE] = true;
				update_obj_stiffness.update_collision[SELF_POINT_TRIANGLE] = true;
				update_obj_stiffness.collision_stiffness[SELF_POINT_TRIANGLE] = temp_stiffness[SELF_POINT_TRIANGLE];
			}
		}
		else {
			ImGui::Text("Saved");
		}
		ImGui::Text("Cloth Edge edge collision: ");
		if (ImGui::InputDouble("##CEEcollision", &temp_stiffness[SELF_EDGE_EDGE], 0.0f, 0.0f, "%.2f")) {
			set_stiffness[EDIT_SELF_EDGE_EDGE] = false;
		}
		ImGui::SameLine();
		if (!set_stiffness[EDIT_SELF_EDGE_EDGE]) {
			if (ImGui::Button("Save##CEEcollision", ImVec2(80, 25))) {
				set_stiffness[EDIT_SELF_EDGE_EDGE] = true;
				update_obj_stiffness.update_collision[SELF_EDGE_EDGE] = true;
				update_obj_stiffness.collision_stiffness[SELF_EDGE_EDGE] = temp_stiffness[SELF_EDGE_EDGE];
			}
		}
		else {
			ImGui::Text("Saved");
		}
		ImGui::Text("Body point triangle collision: ");
		if (ImGui::InputDouble("##BPTcollision", &temp_stiffness[BODY_POINT_TRIANGLE], 0.0f, 0.0f, "%.2f")) {
			set_stiffness[EDIT_BODY_POINT_TRIANGLE] = false;
		}
		ImGui::SameLine();
		if (!set_stiffness[EDIT_BODY_POINT_TRIANGLE]) {
			if (ImGui::Button("Save##BPTcollision", ImVec2(80, 25))) {
				set_stiffness[EDIT_BODY_POINT_TRIANGLE] = true;
				update_obj_stiffness.update_collision[BODY_POINT_TRIANGLE] = true;
				update_obj_stiffness.collision_stiffness[BODY_POINT_TRIANGLE] = temp_stiffness[BODY_POINT_TRIANGLE];
			}
		}
		else {
			ImGui::Text("Saved");
		}
		ImGui::End();
	}
	else {
		memset(set_stiffness + EDIT_BODY_POINT_TRIANGLE, 0, 4);
	}
}

void ImGuiWindows::floorInfo(bool& exist, bool& show, bool& normal_direction, unsigned int& dimension, double* value, bool& eidit)
{
	ImGui::SetNextWindowPos(ImVec2(0, 470));
	ImGui::SetNextWindowSize(ImVec2(240, 240));
	ImGui::Begin("Floor");
	if (exist) {
		if (ImGui::Button("Delete Floor", ImVec2(160, 25))) {
			exist = false;
			show = false;
		}
	}
	else {
		if (ImGui::Button("Create Floor", ImVec2(160, 25))) {
			exist = true;
			show = true;
		}
	}
	if (exist) {
		if (show) {
			if (ImGui::Button("Hide Floor", ImVec2(160, 25))) {
				show = false;
			}
		}
		else {
			if (ImGui::Button("Show Floor", ImVec2(160, 25))) {
				show = true;
			}
		}
		ImGui::Text("Select floor direction");
		std::string current_select =std::to_string(dimension);
		if (ImGui::BeginCombo("##floor_dir", current_select.c_str()))
		{
			std::string name;
			for (int n = 0; n < 3; n++)
			{
				name = std::to_string(n);
				bool is_selected = (dimension == n);
				if (ImGui::Selectable(name.c_str(), is_selected))
				{
					dimension = n;
				}
				if (is_selected)
				{
					ImGui::SetItemDefaultFocus();
				}
			}
			ImGui::EndCombo();
		}
		ImGui::Text("Select floor normal direction");

		std::string  direction;
		if (normal_direction) {
			direction = "positive";
		}
		else {
			direction = "negative";
		}
		if (ImGui::BeginCombo("##floor_normal_dir", direction.c_str()))
		{
			std::string name;
			for (int n = 0; n < 2; n++)
			{
				if (n == 0) {
					name = "negative";
				}
				else {
					name = "positive";
				}
				bool is_selected = ((int)normal_direction == n);
				if (ImGui::Selectable(name.c_str(), is_selected))
				{
					normal_direction = n;
				}
				if (is_selected)
				{
					ImGui::SetItemDefaultFocus();
				}
			}
			ImGui::EndCombo();
		}
		if (ImGui::InputDouble("value##floor", value, 0.0f, 0.0f, "%.3f")) {
			eidit = true;
		}
	}


	ImGui::End();
}

void ImGuiWindows::iterationSolverInfoWindow(double& solver_iteration_num, int& itr_sovler_method,
	char** itr_solver_items, char*& itr_solver_item, int item_num, double* conv_rate, bool& record_matrix)
{
	ImGui::SetNextWindowPos(ImVec2(SCR_WIDTH - 540, 0));
	ImGui::SetNextWindowSize(ImVec2(270, 540));
	ImGui::Begin("Iteration Solver Info");
	ImGui::Text("Convergence_rate:");
	ImGui::InputDouble("##Conv_rate", conv_rate, 0.0, 0.0, "%.2e");
	std::string tempString;
	ImGui::TextWrapped("Ave itr num:");
	ImGui::SameLine();
	ImGui::Text("%.2f", solver_iteration_num);
	ImGui::Text("Solver:");
	ImGui::SetNextItemWidth(260);
	if (ImGui::BeginCombo("##Solver", itr_solver_item))
	{
		for (int n = 0; n < item_num; n++)
		{
			bool is_selected = (itr_solver_item == itr_solver_items[n]);
			if (ImGui::Selectable(itr_solver_items[n], is_selected))
			{
				itr_solver_item = itr_solver_items[n];
				itr_sovler_method = n;
			}
			if (is_selected)
			{
				ImGui::SetItemDefaultFocus();
			}
		}
		ImGui::EndCombo();
	}
	if (ImGui::Button("Save Matrix", ImVec2(160, 25))) {
		record_matrix = true;
	}
	ImGui::End();
}
