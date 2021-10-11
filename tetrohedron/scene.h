#pragma once
#include"basic/camera.h"
#include"load_mesh/preprocessing.h"
#include"render/shadow.h"
#include"load_mesh/writeObj.h"
#include"basic/intersection.h"
#include"basic/cursor.h"

class Scene
{
public:
	Scene();
	size_t time_stamp;
	std::vector<Cloth>cloth;
	std::vector<Collider>collider;
	std::vector<TetrohedronObject>tetrohedron;
	Intersection intersection;
	Shadow shadow;
	double camera_center[3];
	bool* control_parameter;
	void initialCloth();
	void resetCloth();
	void loadMesh(std::vector<std::string>& collider_path, std::vector<std::string>& object_path);
	void drawScene(Camera* camera, std::vector<std::vector<bool>>& wireframe, std::vector<std::vector<bool>>& hide, bool start_save_obj);
	void getClothInfo(std::vector<std::array<int, 3>>& mesh_info, std::vector<double>& mass, std::vector<std::array<double, 3>>& mesh_stiffness, double* simulation_parameter, std::vector<std::array<double, 4>>& collision_stiffness);
	void updateCloth(Camera* camera, double* cursor_screen, bool* control_parameter, float force_coe);
	void initialIntersection();
private:
	Light light;
	int cloth_num, collider_num,tetrohedron_num;
	
	void initialSceneSetting(Preprocessing& preprocessing);
	double time_step;
	size_t last_output_obj_stamp;
	Thread thread;
	WriteObj save_obj;
	void setWireframwColor();
	std::vector<int>cloth_index_in_object;
	std::vector<int>tetrohedron_index_in_object;

	Cursor cursor;
	void saveObj();
	void updateBuffer();
};