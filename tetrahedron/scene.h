#pragma once
#include<chrono>
#include"basic/camera.h"
#include"load_mesh/preprocessing.h"
#include"render/shadow.h"
#include"load_mesh/writeObj.h"
#include"basic/intersection.h"
#include"basic/cursor.h"
#include"project_dynamic.h"
#include"basic/pick_triangle.h"
#include"basic/set_tetrahedron_anchor.h"
#include"basic/object_chosen_indicator.h"
#include"XPBD/XPBD.h"
#include"./basic/floor.h"
#include"./basic/move_object.h"
#include"collision/test_draw_collision.h"
#include"basic/opengl_input.h"
#include"newton_method.h"
#include"basic/move_model.h"
#include"mesh_struct/graph_color.h"
#include"XPBD_large_system.h"
#include"XPBD_IPC.h"
#include"basic/draw_triangle.h"

//#include<windows.h>
//#include"basic/drawEdge.h"

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
	bool loadMesh(std::string& scene_path, std::vector<std::string>& collider_path, std::vector<std::string>& object_path, double* tolerance_ratio,
		bool* control_parameter, double* initial_stiffness, double* friction_coe, unsigned int* sub_step_per_detection, bool* floor_indicator, unsigned int& floor_dimension, double& floor_value);
	void drawScene(Camera* camera, std::vector<std::vector<bool>>& show_element,
		bool* control_parameter);
	void getClothInfo(std::vector<std::array<int, 3>>& mesh_info, std::vector<double>& mass, std::vector<std::array<double, 6>>& mesh_stiffness, double* simulation_parameter, std::vector<std::array<double, 8>>& collision_stiffness);
	void updateCloth(Camera* camera, Input* input, bool* control_parameter, float force_coe, bool& record_matrix,
		double& ave_iteration, int& select_hash_cell_index); // bool mouse_is_pressed_previous_current_frame
	void initialIntersection();
	void obtainConvergenceInfo(double* convergence_rate, int* iteration_num);
	void updateConvRate(double* convergence_rate);
	void setTolerance(double* tolerance_ratio);
	void testBVH();
	void obtainCursorIntersection(double* pos, Camera* camera, std::vector<std::vector<bool>>& hide);
	void getTetrahedronInfo(std::vector<std::array<int, 3>>& mesh_info, std::vector<double>& mass, std::vector<std::array<double, 6>>& mesh_stiffness, double* simulation_parameter, std::vector<std::array<double, 8>>& collision_stiffness);
	void selectAnchor(bool* control_parameter, bool* select_anchor, double* screen_pos, bool press_state, bool pre_press_state, Camera* camera, std::vector<bool>& hide);
	void drawSelectRange(bool* select_anchor, bool press_state, bool pre_press_state);
	void updateIterateSolverParameter(double rate, int itr_solver_method);
	void resetIntersectionState();

	void testForWritetToArray(int thread_No);
	void testForWritetToArraySingle(int total_thread_num);
	void updateStiffness(UpdateObjStiffness& update_obj_stiffness, std::vector<std::array<double, 6>>& cloth_stiffness,
		std::vector<std::array<double, 6>>& tet_stiffness, std::vector<std::array<double, 8>>& cloth_collision_stiffness,
		std::vector<std::array<double, 8>>& tet_collision_stiffness);
	void updateItrInfo(int* iteration_num);

	void setFloorInfo(bool exist, bool show, bool normal_direction, unsigned int dimension, double value, bool& eidit, bool* control_parameter);

	void pickAxes(double* pos, Camera* camera);

	unsigned int select_dimension_index;
	bool start_rotation;

	Input* input;
	void setDampStiffness(double* damp_stiffness, double* rayleigh_damp_stiffness);
	double* time_per_frame;

	std::chrono::microseconds time_accumulation;

	void readScene(std::string& path);
	void saveScene();
	void saveParameter(std::vector<std::string>& path, std::vector<std::string>& collider_path, std::vector<std::array<double, 6>>& cloth_stiffness, std::vector<std::array<double, 6>>& tet_stiffness,
		std::vector<std::array<double, 8>>& cloth_collision_stiffness, std::vector<std::array<double, 8>>& tet_collision_stiffness, double* tolerance_ratio,
		double* friction_coe);

	Camera* camera;

private:
	double cloth_density;
	double tetrahedron_density;
	unsigned int use_method;
	bool only_test_collision;
	Light light;
	int cloth_num, collider_num, tetrahedron_num;

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
	XPBD xpbd;
	XPBD_IPC xpbd_ipc;
	NewtonMethod newton_method;
	SecondOrderLargeSystem second_order_xpbd_large;

	double ave_edge_length;

	void setAveEdgeLength();
	PickTriangle pick_triangle;

	void getCursorPos(double* cursor_pos, std::vector<std::array<double, 3>>& vertex, int* vertex_index);
	void setCursorForce(Camera* camera, double* cursor_screen, float force_coe);
	void cursorMovement(Camera* camera, double* cursor_screen, double* force_direction, float force_coe,
		double* object_position, double* cursor_pos_in_space);
	double max_force_magnitude;
	SetTetrahedronAnchor set_tetrahedron_anchor;

	Shader* wireframe_shader;
	Shader* object_shader_front_soft_edge;
	Shader* object_shader_front;
	Shader* object_shader_front_sharp_edge;
	Shader* object_shader_texture;
	void genShader();

	unsigned int total_obj_num;
	unsigned int obj_num_except_collider;
	ObjectChosenIndicator object_chosen_indicator;
	void setChosenIndicator();




	unsigned int** test_pair;
	unsigned int test_array_size;
	unsigned int** test_pair_single;

	void compareArray();

	void voidForWritetToArraySingle(int total_thread_num);
	unsigned int iteration_num_for_test_array;
	//unsigned int thread_test[80001];

	void setObjMoveInfo(Camera* camera, double* cursor_screen);
	//DrawCulling draw_culling;
	void getAABB();
	Floor floor;
	MoveObject move_object;
	void moveObj(Camera* camera, double* cursor_screen, bool only_move_vertex_pos);
	void updateObjSimulation(Camera* camera, double* cursor_screen, bool* control_parameter, float force_coe, bool& record_matrix,
		double& ave_iteration);
	void updateSceneCollisionTest(bool* control_parameter);
	TestDrawCollision test_draw_collision;

	//void getCurrentAABB();

	void updateBufferOriPos();

	bool intersect_when_rotation;

	void rotate(Camera* camera, float* angle, bool only_move_vertex_pos);

	bool sameDirection(Camera* camera, unsigned int obj_index);

	void setMove(Camera* camera, Input* input, bool* control_parameter);
	void drawSpatialHashing(Camera* camera, Input* input, int& select_hash_cell);

	MoveModel move_model;

	unsigned int time_record_interval;

	unsigned int time_indicate_for_simu;

	bool use_single_thread = true;
	GraphColor graph_color;
	void setGroup();

	std::vector<std::vector<int>*>anchor_vertex;
	std::vector<MeshStruct*> mesh_struct;

	void reorganizeData();


	void readScene();
	void updateAnchorTet();
	void reflectModel();
	void checkVolume();
	double force_direction[3];
	bool read_force;
	void addExternalForce();
	//Draw_Edge draw_edge_;

	DrawVertex draw_vertex;
	DrawTriangle draw_triangle;
};
