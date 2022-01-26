#pragma once
#include"basic/camera.h"
#include"load_mesh/preprocessing.h"
#include"render/shadow.h"
#include"load_mesh/writeObj.h"
#include"basic/intersection.h"
#include"basic/cursor.h"
#include"project_dynamic.h"
#include"basic/pick_triangle.h"
#include"basic/set_tetrahedron_anchor.h"

class Scene
{
public:
	Scene();
	size_t time_stamp;
	std::vector<Cloth>cloth;
	std::vector<Collider>collider;
	std::vector<Tetrahedron>tetrahedron;
	Intersection intersection;
	Shadow shadow;
	double camera_center[3];
	bool* control_parameter;
	void initial();
	void reset();
	void loadMesh(std::vector<std::string>& collider_path, std::vector<std::string>& object_path);
	void drawScene(Camera* camera, std::vector<std::vector<bool>>& wireframe, std::vector<std::vector<bool>>& hide, bool start_save_obj);
	void getClothInfo(std::vector<std::array<int, 3>>& mesh_info, std::vector<double>& mass, std::vector<std::array<double, 3>>& mesh_stiffness, double* simulation_parameter, std::vector<std::array<double, 4>>& collision_stiffness);
	void updateCloth(Camera* camera, double* cursor_screen, bool* control_parameter, float force_coe, bool& record_matrix,
		double& ave_iteration);
	void initialIntersection();
	void obtainConvergenceInfo(double* convergence_rate, int* iteration_num);
	void updateConvRate(double* convergence_rate);
	void setTolerance(double* tolerance_ratio);
	void testBVH();
	void obtainCursorIntersection(double* pos, Camera* camera, std::vector<std::vector<bool>>& hide);
	void getTetrahedronInfo(std::vector<std::array<int, 3>>& mesh_info, std::vector<double>& mass, std::vector<std::array<double, 2>>& mesh_stiffness, double* simulation_parameter, std::vector<std::array<double, 4>>& collision_stiffness);
	void selectAnchor(bool* control_parameter, bool* select_anchor, double* screen_pos, bool press_state, bool pre_press_state, Camera* camera, std::vector<bool>& hide);
	void drawSelectRange(bool* select_anchor, bool press_state, bool pre_press_state);
	void updateIterateSolverParameter(double rate, int itr_solver_method);
private:
	Light light;
	int cloth_num, collider_num,tetrahedron_num;
	
	void initialSceneSetting(Preprocessing& preprocessing);
	double time_step;
	size_t last_output_obj_stamp;
	Thread thread;
	WriteObj save_obj;
	void setWireframwColor();

	Cursor cursor;
	void saveObj();
	void updateBuffer();
	ProjectDynamic project_dynamic;


	double ave_edge_length;
	
	void setAveEdgeLength();
	PickTriangle pick_triangle;
	
	void getCursorPos(double* cursor_pos, std::vector<std::array<double, 3>>& vertex, int* vertex_index);
	void setCursorForce(Camera* camera, double* cursor_screen, float force_coe);
	void cursorMovement(Camera* camera, double* cursor_screen, double* force_direction, float force_coe, double* object_position);
	double max_force_magnitude;
	SetTetrahedronAnchor set_tetrahedron_anchor;

	Shader* wireframe_shader;
	Shader* object_shader_front;
	void genShader();
};
