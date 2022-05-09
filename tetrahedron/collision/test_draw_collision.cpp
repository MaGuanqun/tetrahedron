#include"test_draw_collision.h"


void TestDrawCollision::initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider,
	std::vector<Tetrahedron>* tetrahedron, Thread* thread, Floor* floor, double* tolerance_ratio)
{
	collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio);
	draw_collision.initial(cloth, collider, tetrahedron, thread);
	draw_collision.setInPairInfo(collision.point_triangle_target_pos_index.data(), collision.point_triangle_collider_target_pos_index.data(), collision.edge_edge_target_pos_index.data());
	this->cloth = cloth;
	this->collider = collider;
	this->tetrahedron = tetrahedron;
}



void TestDrawCollision::setCollisionData()
{
	collision.collisionCulling();
	collision.getCollisionPair();
	draw_collision.setElementIndices();
}



void TestDrawCollision::drawCollision(bool draw_VT, Light& light,Camera* camera, Shader* object_shader_front, 
	std::vector<std::vector<bool>>& show_element, Shadow* shadow, Shader* wireframe_shader)
{
	draw_collision.drawCollision(draw_VT, light, camera, object_shader_front, show_element);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_CUBE_MAP, shadow->depth_map);
	glEnable(GL_CULL_FACE);
	for (int j = 0; j < cloth->size(); ++j) {
		if (show_element[9+CLOTH_][j]) {
			cloth->data()[j].setSceneShader(light, camera, shadow->far_plane, object_shader_front);
			cloth->data()[j].drawOriPos(camera, object_shader_front);
		}
	}
	glCullFace(GL_BACK);
	for (int j = 0; j < tetrahedron->size(); ++j) {
		if (show_element[9 + TETRAHEDRON_][j]) {
			tetrahedron->data()[j].setSceneShader(light, camera, shadow->far_plane, object_shader_front);
			tetrahedron->data()[j].drawOriPos(camera, object_shader_front);
		}
	}
	for (int j = 0; j < collider->size(); ++j) {
		if (show_element[9 + COLLIDER_][j]) {
			collider->data()[j].setSceneShader(light, camera, shadow->far_plane, object_shader_front);
			collider->data()[j].draw(camera, object_shader_front);
		}
	}
	glDisable(GL_CULL_FACE);
	for (int i = 0; i < collider->size(); ++i) {
		if (show_element[12 + COLLIDER_][i]) {
			collider->data()[i].drawWireframeOriPos(camera, wireframe_shader);
		}
	}
	for (int j = 0; j < cloth->size(); ++j) {
		if (show_element[12 + CLOTH_][j]) {
			cloth->data()[j].drawWireframeOriPos(camera, wireframe_shader);
		}
	}
	for (int j = 0; j < tetrahedron->size(); ++j) {
		if (show_element[12 + TETRAHEDRON_][j]) {
			tetrahedron->data()[j].drawWireframeOriPos(camera, wireframe_shader);
		}
	}
	
}

