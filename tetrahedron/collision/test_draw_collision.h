#pragma once
#include"collision.h"
#include"draw_collision.h"
#include"draw_spatial_hashing.h"

class TestDrawCollision
{
public:
	
	Collision collision;

	void initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider,
		std::vector<Tetrahedron>* tetrahedron, Thread* thread, Floor* floor, double* tolerance_ratio);
	void setCollisionData();
	void drawCollision(bool draw_VT, Light& light,  Camera* camera, Shader* object_shader_front, 
		std::vector<std::vector<bool>>& drawCollision, Shadow* shadow,  Shader* wireframe_shader);

	DrawSpatialHashing draw_spatial_hashing;
	void setForOriSpatialHashing();
private:
	DrawCollision draw_collision;
	std::vector<Cloth>* cloth;
	std::vector<Collider>* collider;
	std::vector<Tetrahedron>* tetrahedron;
	Thread* thread;

};
