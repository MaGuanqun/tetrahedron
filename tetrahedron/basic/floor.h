#pragma once
#include"camera.h"
#include"../external/shader.h"
#include"../global_struct.h"
#include"../external/stb-master/stb_image.h"
#include"../render/shadow.h"

class Floor
{
public:
	Floor();
	void draw(Camera* camera, Shader* object_shader_front, Shadow* shadow, Light& light, float& far_plane);

	void setFloor(unsigned int dimension, double value, bool normal_direction);

	bool exist;

	bool normal_direction;
	unsigned int dimension;
	double value;
	 
private:
	double half_floor_length;
	std::vector<double> vertex_position;
	std::vector<unsigned int> vertex_index;
	std::vector<double> vertex_texture_coordinate;
	std::vector<double> vertex_normal_;
	MeshMaterial mesh_material;
	void genBuffer();

	unsigned int VAO;
	unsigned int VBO[3];
	unsigned int EBO;
	unsigned int texture1;

	void setBuffer();
	void genTexture();
};
