#include"floor.h"


Floor::Floor()
{
	half_floor_length = 10.0;

	vertex_position.resize(12);
	unsigned int vertex_indi[6] = { 0,1,2,0,2,3 };
	vertex_index.resize(6);
	memcpy(vertex_index.data(), vertex_indi, 24);
	double vertex_texture_coordi[8] = { 0.0,1.0,0.0,0.0,1.0,0.0,1.0,1.0 };
	vertex_texture_coordinate.resize(8);
	memcpy(vertex_texture_coordinate.data(), vertex_texture_coordi, 64);
	vertex_normal_.resize(12);
	genBuffer();
	genTexture();

	setFloor(1, -0.35,true);

	float Kd[3] = { 0.5,0.5,0.5 };
	float Ka[3] = { 0.4,0.4,0.4 };
	float Ks[3] = { 0.3,0.3,0.3 };
	mesh_material.setMaterial(Kd, Ka, Ks);
}


void Floor::setFloor(unsigned int dimension, double value, bool normal_direction)
{
	this->dimension = dimension;
	double n[3] = { 0,0,0 };
	this->value = value;
	this->normal_direction = normal_direction;
	switch (dimension)
	{
	case 0: {
		double a[12] = { value, half_floor_length,half_floor_length, value,-half_floor_length,half_floor_length,
		value,-half_floor_length,-half_floor_length,value,half_floor_length,-half_floor_length };
		memcpy(vertex_position.data(), a, 96);
		if (normal_direction) {
			n[0] = 1.0;
		}
		else {
			n[0] = -1.0;
		}
	}
		break;
	case 1: {
		double a[12] = { -half_floor_length, value, -half_floor_length, half_floor_length, value,-half_floor_length,half_floor_length,
		value,half_floor_length,-half_floor_length,value,half_floor_length };
		memcpy(vertex_position.data(), a, 96);
		if (normal_direction) {
			n[1] = 1.0;
		}
		else {
			n[1] = -1.0;
		}
	}
		break;
	case 2: {
		double a[12] = { half_floor_length, half_floor_length, value, half_floor_length, -half_floor_length, value,-half_floor_length,-half_floor_length,
		value,-half_floor_length,half_floor_length,value };
		memcpy(vertex_position.data(), a, 96);
		if (normal_direction) {
			n[2] = 1.0;
		}
		else {
			n[2] = -1.0;
		}

	}
		break;
	}
	for (unsigned int i = 0; i < 4; ++i) {
		memcpy(vertex_normal_.data() + 3 * i, n, 24);
	}

	setBuffer();
}


void Floor::draw(Camera* camera, Shader* object_shader_front, Shadow* shadow, Light& light, float& far_plane)
{
	object_shader_front->use();
	object_shader_front->setInt("depthMap", 0);
	object_shader_front->setInt("texture1", 1);
	object_shader_front->setFloat("far_plane", far_plane);
	object_shader_front->setBool("lightIsChosen", true);
	object_shader_front->setVec3("lightPos", camera->position);
	object_shader_front->setVec3("light.ambient", light.ambient);
	object_shader_front->setVec3("light.diffuse", light.diffuse);
	object_shader_front->setVec3("light.specular", light.specular);
	object_shader_front->setVec3("viewPos", camera->position);
	object_shader_front->setBool("lightShadowOn", true);
	object_shader_front->setMat4("projection", camera->GetProjectMatrix());
	object_shader_front->setMat4("view", camera->GetViewMatrix());
	object_shader_front->setMat4("model", glm::mat4(1.0));
	object_shader_front->setFloat("transparence", 1.0);
	object_shader_front->setVec3("material.Kd", glm::vec3(mesh_material.Kd[0], mesh_material.Kd[1], mesh_material.Kd[2]));
	object_shader_front->setVec3("material.Ka", glm::vec3(mesh_material.Ka[0], mesh_material.Ka[1], mesh_material.Ka[2]));
	object_shader_front->setVec3("material.Ks", glm::vec3(mesh_material.Ks[0], mesh_material.Ks[1], mesh_material.Ks[2]));
	object_shader_front->setFloat("material.Ns", mesh_material.Ns);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_CUBE_MAP, shadow->depth_map);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, texture1);

	glBindVertexArray(VAO);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawElements(GL_TRIANGLES, vertex_index.size(), GL_UNSIGNED_INT, 0);
}

void Floor::genBuffer()
{
	glGenVertexArrays(1, &VAO);
	glGenBuffers(3, VBO);
	glGenBuffers(1, &EBO);
}


void Floor::setBuffer()
{

	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
	glBufferData(GL_ARRAY_BUFFER, vertex_position.size() * sizeof(double), vertex_position.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, vertex_index.size() * sizeof(unsigned int), vertex_index.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
	glBufferData(GL_ARRAY_BUFFER, vertex_normal_.size() * sizeof(double), vertex_normal_.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[2]);
	glBufferData(GL_ARRAY_BUFFER, vertex_texture_coordinate.size() * sizeof(double), vertex_texture_coordinate.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), (void*)0);
	glBindVertexArray(0);

}


void Floor::genTexture()
{
	glGenTextures(1, &texture1);

	std::string path = "./basic/floor.jpg";
	int width, height, nrChannels;
	unsigned char* data = stbi_load(path.c_str(), &width, &height, &nrChannels, 0);
	if (data)
	{
		GLenum format;
		if (nrChannels == 1)
			format = GL_RED;
		else if (nrChannels == 3)
			format = GL_RGB;
		else if (nrChannels == 4)
			format = GL_RGBA;

		glBindTexture(GL_TEXTURE_2D, texture1);
		glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
		glGenerateMipmap(GL_TEXTURE_2D);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		stbi_image_free(data);
	}
}