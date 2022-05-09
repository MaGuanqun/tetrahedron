#pragma once
#include"collision.h"
#include"draw_collision.h"

class TestDrawCollision
{
public:
	
	Collision collision;

	void initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider,
		std::vector<Tetrahedron>* tetrahedron, Thread* thread, Floor* floor, double* tolerance_ratio);
	void setCollisionData();
	void drawCollision(bool draw_VT, Light& light,  Camera* camera, Shader* object_shader_front, std::vector<std::vector<bool>>& show_collision_element);
private:
	DrawCollision draw_collision;
};
