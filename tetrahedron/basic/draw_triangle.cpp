#include"draw_triangle.h"

DrawTriangle::DrawTriangle()
{
	genBuffer();


}


void DrawTriangle::genBuffer()
{
	glGenVertexArrays(1, &VAO);
	glGenBuffers(2, VBO);
	glGenBuffers(1, &EBO);
}

void DrawTriangle::drawTriangle(Camera* camera, Shader* object_shader_front, std::vector<std::array<double, 3>>& vertex_for_render,
	std::vector<std::array<int, 3>>& triangle_vertex_index, std::vector<std::array<double, 3>>& norm_for_render, std::vector<unsigned int>& triangle_index, glm::vec3 color)
{
	if (triangle_index.empty()) {
		return;
	}

	std::vector<int> actual_index;
	actual_index.reserve(3 * triangle_index.size());
	for (unsigned int i = 0; i < triangle_index.size(); ++i) {
		actual_index.emplace_back(triangle_vertex_index[triangle_index[i]][0]);
		actual_index.emplace_back(triangle_vertex_index[triangle_index[i]][1]);
		actual_index.emplace_back(triangle_vertex_index[triangle_index[i]][2]);
	}

	setBuffer(actual_index, vertex_for_render, norm_for_render);


	glm::vec3 color_a = 0.1f * color;
	glm::vec3 color_d = 0.3f * color;
	glm::vec3 color_s = 0.6f * color;

	object_shader_front->use();
	object_shader_front->setVec3("viewPos", camera->position);
	object_shader_front->setBool("lightShadowOn", true);
	object_shader_front->setMat4("projection", camera->GetProjectMatrix());
	object_shader_front->setMat4("view", camera->GetViewMatrix());
	object_shader_front->setMat4("model", glm::mat4(1.0));
	object_shader_front->setFloat("transparence", 1.0);
	object_shader_front->setVec3("material.Kd",color_d);
	object_shader_front->setVec3("material.Ka", color_a);
	object_shader_front->setVec3("material.Ks", color_s);
	object_shader_front->setFloat("material.Ns", 32);

	glBindVertexArray(VAO);
	glPolygonMode(GL_FRONT, GL_FILL);
	glCullFace(GL_BACK);
	glDrawElements(GL_TRIANGLES, actual_index.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

}

void DrawTriangle::setBuffer(std::vector<int>& triangle_vertex_index, std::vector<std::array<double, 3>>& vertex_for_render,
	std::vector<std::array<double, 3>>& norm_for_render)
{
	//std::cout << mesh_struct.vertex_position[0][0] << " " << mesh_struct.vertex_position[0][1] << " " << mesh_struct.vertex_position[0][2]<<std::endl;
	//std::cout << mesh_struct.vertex_position[1][0] << " " << mesh_struct.vertex_position[1][1] << " " << mesh_struct.vertex_position[1][2]<<std::endl;
	//std::cout << mesh_struct.edge_vertices[0] << " " << mesh_struct.edge_vertices[1] << std::endl;
		glBindVertexArray(VAO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
		glBufferData(GL_ARRAY_BUFFER, vertex_for_render.size() * sizeof(std::array<double, 3>), vertex_for_render[0].data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangle_vertex_index.size() * sizeof(int), triangle_vertex_index.data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
		glBufferData(GL_ARRAY_BUFFER, norm_for_render.size() * sizeof(std::array<double, 3>), norm_for_render[0].data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindVertexArray(0);
	
	
}