#pragma once
#include<array>
#include<vector>
#include"../external/glm/glm.hpp"
#include"../external/glm/gtc/matrix_transform.hpp"
#include"../external/shader.h"
#include"../object/tetrohedron.h"
#include"camera.h"
#include"draw_vertex.h"

class SetTetrohedronAnchor
{
public:
	SetTetrohedronAnchor();
	void setCorner(double* screen_pos, bool pre_press_state, std::vector<Tetrohedron>& tetrohedron, Camera* camera,
		std::vector<bool>& hide);
	void draw();
private:
	glm::vec3 camera_pos;
	std::vector<std::array<float, 3>>position;
	float draw_first_corner[3];
	float draw_last_corner[3];
	void transferTo3D(double* screen_pos, float* point);

	void setPosition();
	unsigned int VAO, VBO;
	glm::mat4 view;
	glm::mat4 projection;
	Shader* shader;
	void genBuffer();
	void setBufferData();
	void findAllSelectVertex(std::vector<std::array<double, 3>>& vertex_for_render, Camera* camera, std::vector<int>&vertex_index);
	float max_corner[2];
	float min_corner[2];
	void findSurfaceVertex(std::vector<std::array<double, 3>>& vertex_for_render, Camera* camera, std::vector<int>& vertex_index, std::vector<bool>& is_surface);
	DrawVertex draw_vertex;
	glm::vec3 vertex_color;

	float first_corner_in_clip_space[2];
	void setClipSpaceRange(double* screen_pos);
};


