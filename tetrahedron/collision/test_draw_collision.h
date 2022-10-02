#pragma once
#include"collision.h"
#include"draw_collision.h"
#include"draw_spatial_hashing.h"

class TestDrawCollision
{
public:
	
	Collision* collision;

	void initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider,
		std::vector<Tetrahedron>* tetrahedron, Thread* thread, Floor* floor, double* tolerance_ratio,
		Collision* pd_collision, Collision* pbd_collision, Collision* newton_collision, Collision* xpbd_ipc_collision, unsigned int simulation_method);
	void setCollisionData();
	void drawCollision(bool draw_VT, Light& light,  Camera* camera, Shader* object_shader_front, 
		std::vector<std::vector<bool>>& drawCollision, Shadow* shadow,  Shader* wireframe_shader, bool draw_all_collision_pair, bool draw_all_element);

	DrawSpatialHashing draw_spatial_hashing;
	void setForOriSpatialHashing();
	void setForSelectCell();
	void setForSelectCellOne(int index);
	void obtianSpatialHashingCell(Camera* camera, double* cursor_screen);

	void obtainElementsInOneCell(int& index_chosen);
	DrawCollision draw_collision;
private:
	
	std::vector<Cloth>* cloth;
	std::vector<Collider>* collider;
	std::vector<Tetrahedron>* tetrahedron;
	Thread* thread;
	std::vector<unsigned int> select_hash_index;
	std::vector<unsigned int> select_ori_hash_index;

	std::vector<std::vector<unsigned int>> vertex_index_in_one_cell;
	std::vector<std::vector<unsigned int>> triangle_index_in_one_cell;
	std::vector<std::vector<unsigned int>> edge_index_in_one_cell;

	unsigned int total_obj_num;
	bool isIndexNotEmpty();



};
