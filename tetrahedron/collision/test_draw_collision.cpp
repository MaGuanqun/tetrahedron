#include"test_draw_collision.h"


void TestDrawCollision::initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider,
	std::vector<Tetrahedron>* tetrahedron, Thread* thread, Floor* floor, double* tolerance_ratio)
{
	collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio);
	draw_collision.initial(cloth, collider, tetrahedron, thread);
	draw_collision.setInPairInfo(collision.point_triangle_target_pos_index.data(), collision.point_triangle_collider_target_pos_index.data(), collision.edge_edge_target_pos_index.data());
}



void TestDrawCollision::setCollisionData()
{
	collision.collisionCulling();
	collision.getCollisionPair();
	draw_collision.setElementIndices();
}



void TestDrawCollision::drawCollision(bool draw_VT, Light& light,Camera* camera, Shader* object_shader_front, std::vector<std::vector<bool>>& show_collision_element)
{
	draw_collision.drawCollision(draw_VT, light, camera, object_shader_front, show_collision_element);
}

