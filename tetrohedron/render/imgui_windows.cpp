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
		if (ImGui::Button("1 Frame", ImVec2(160, 25)))
		{
			control_parameter[ONE_FRAME] = true;
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
			}
		}
		if (ImGui::Button("Reset Simulation", ImVec2(160, 25))) {
			control_parameter[RESET_SIMULATION] = true;
		}
		ImGui::SameLine();
		helpMarker(
			"Reset with current stiffness");
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
	ImGui::SetNextItemOpen(true);
	if (ImGui::TreeNode("Start Test")) {
		if (ImGui::Button("Start test", ImVec2(160, 25)))
		{
			control_parameter[START_TEST] = true;
			control_parameter[INITIAL_TEST] = true;
		}
		ImGui::TreePop();
	}

	ImGui::End();

	if (control_parameter[SET_CURSOR_FORCE]) {
		ImGui::SetNextWindowSize(ImVec2(330, 240));
		ImGui::Begin("Set Cursor Force", &control_parameter[SET_CURSOR_FORCE]);
		ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.9f);
		ImGui::SetNextItemOpen(true);
		if (ImGui::TreeNode("Set force coefficient")) {
			ImGui::SliderFloat("##set force coefficient", force_coe, 0.5, 1.0f, "force coefficient = %.2f");
			ImGui::TreePop();
		}
		ImGui::TextWrapped("Force depends on force coefficient and cursor moving speed.");
		ImGui::End();
	}
}


void ImGuiWindows::visualizationControlPanel(bool& reset_camera, std::vector<std::vector<bool>>& wireframe, std::vector<std::vector<bool>>& hide)
{
	ImGui::SetNextWindowPos(ImVec2(0, 630));
	ImGui::SetNextWindowSize(ImVec2(240, 370));
	ImGui::Begin("Visualization");
	ImGui::SetNextItemOpen(true);
	//	if (ImGui::TreeNode("Visualization")) {
	int cloth_no = -1;;
	if (ImGui::Button("Reset Camera", ImVec2(160, 25)))
	{
		reset_camera = true;
	}
	std::string tempString;
	for (int i = 0; i < wireframe[COLLIDER].size(); ++i) {
		ImGui::SetNextItemOpen(true);
		tempString = "Collider " + std::to_string(i);
		if (ImGui::TreeNode(tempString.c_str())) {
			if (hide[COLLIDER][i]) {
				tempString = "Show##Collider" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					hide[COLLIDER][i] = false;
				}
			}
			else {
				tempString = "Hide##Collider" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					hide[COLLIDER][i] = true;
				}
			}
			if (wireframe[COLLIDER][i]) {
				tempString = "Hide WireFrame##Collider" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					wireframe[COLLIDER][i] = false;
				}
			}
			else {
				tempString = "Show WireFrame##Collider" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					wireframe[COLLIDER][i] = true;
				}
			}
			ImGui::TreePop();
		}
	}
	for (int i = 0; i < wireframe[TETROHEDRON].size(); ++i) {
		ImGui::SetNextItemOpen(true);
		tempString = "Tetrohedron " + std::to_string(i);
		if (ImGui::TreeNode(tempString.c_str())) {
			if (hide[TETROHEDRON][i]) {
				tempString = "Show##Tetrohedron" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					hide[TETROHEDRON][i] = false;
				}
			}
			else {
				tempString = "Hide##Tetrohedron" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					hide[TETROHEDRON][i] = true;
				}
			}
			if (wireframe[TETROHEDRON][i]) {
				tempString = "Hide WireFrame##Tetrohedron" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					wireframe[TETROHEDRON][i] = false;
				}
			}
			else {
				tempString = "Show WireFrame##Tetrohedron" + std::to_string(i);
				if (ImGui::Button(tempString.c_str(), ImVec2(160, 25)))
				{
					wireframe[TETROHEDRON][i] = true;
				}
			}
			ImGui::TreePop();
		}
	}
	ImGui::End();
}


bool ImGuiWindows::loadModel(std::vector<std::string>& collider_path, std::vector<std::string>& tetrohedron_path)
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
		ImGui::TextWrapped("Load colliders and tetrohedrons.");
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
		if (ImGui::Button("Load Tetrohedron", ImVec2(160, 25))) {
			m_file_dialog_info.SetTypeFilters({ ".ele" });
			m_file_dialog_info.SetTitle("Load Tetrohedron");
			m_file_dialog_info.Open();
			load_collider_is_open = false;
		}
		if (!load_collider_is_open) {
			m_file_dialog_info.Display();
			if (m_file_dialog_info.HasSelected()) {
				if (!tetrohedron_path.empty() && tetrohedron_path[tetrohedron_path.size() - 1] == m_file_dialog_info.GetSelected().string()) {
					ImGui::TextWrapped("Error, Need to add different path");
				}
				else {
					tetrohedron_path.push_back(m_file_dialog_info.GetSelected().string());
					m_file_dialog_info.ClearSelected();
				}		
			}
		}
		if (!tetrohedron_path.empty()) {
			ImGui::TextWrapped("Tetrohedron path:");
			for (int i = 0; i < tetrohedron_path.size(); ++i) {
				ImGui::TextWrapped(tetrohedron_path[i].c_str());
			}
		}
		if (ImGui::Button("Load Tetrohedron", ImVec2(160, 25))) {
			finished_loading = true;
		}
		ImGui::End();
	}
	if (finished_loading) {
		return true;
	}
	return false;	 
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