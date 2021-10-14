#pragma once
#include"../object/cloth.h"
#include"../object/collider.h"
#include"../object/tetrohedron.h"

class Shadow
{
public:
	Shadow();
	unsigned int depth_map;
	float far_plane;
	float camera_from_origin;
	void setBasic();
	void drawShadow(Camera* camera, std::vector<std::vector<bool>>& hide, std::vector<Cloth>& cloth, std::vector<Collider>& collider,
		std::vector<Tetrohedron>& tetrohedron, std::vector<int>& cloth_index_in_object, std::vector<int>& tetrohedron_index_in_object);
private:
	Shader* shadow_shader;

	const unsigned int SHADOW_WIDTH = 1024, SHADOW_HEIGHT = 1024;
	unsigned int depth_map_FBO;
	float near_plane = 0.1f;

	std::vector<glm::mat4>shadow_view;//This is the vector of lightview for point light
	std::vector<glm::mat4> shadow_transforms;

	glm::mat4 light_projection_matrix, light_view_matrix;
	void FBOdepth();
	void lightSpace(glm::vec3& light_pos);
	void setIn(Camera* camera);
};
