#include"drawEdge.h"

Draw_Edge::Draw_Edge()
{
	genBuffer();
}

void Draw_Edge::genBuffer()
{
	glGenVertexArrays(1, &EE_VAO);
	glGenBuffers(1, &EE_VBO);
	glGenBuffers(1, &EE_EBO);
}


void Draw_Edge::drawEdge(Camera* camera, Shader* wireframe_shader, std::vector<std::array<double, 3>>& vertex_for_render,
	std::vector<unsigned int>& edge_index, glm::vec3 color, std::vector<unsigned int>& edge_vertex_index)
{
	std::vector<unsigned int> edge_v_index(2 * edge_index.size());
	for (unsigned int i = 0; i < edge_index.size(); ++i) {
		memcpy(edge_v_index.data() + (i << 1), edge_vertex_index.data() + (edge_index[i] << 1), 8);
	}
	drawEdge(camera, wireframe_shader, vertex_for_render, edge_v_index, color);
}

void Draw_Edge::drawEdge(Camera* camera, Shader* wireframe_shader, std::vector<std::array<double, 3>>& vertex_for_render,
	std::vector<unsigned int>& edge_vertex_index, glm::vec3 color)
{
	setBuffer(vertex_for_render, edge_vertex_index);
	glLineWidth(4.0);
	wireframe_shader->use();
	wireframe_shader->setMat4("projection", camera->GetProjectMatrix());
	wireframe_shader->setMat4("view", camera->GetViewMatrix());
	wireframe_shader->setMat4("model", glm::mat4(1.0));
	wireframe_shader->setVec3("color", color);
	wireframe_shader->setFloat("transparent", 1.0f);
	glBindVertexArray(EE_VAO);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawElements(GL_LINES, edge_vertex_index.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	glLineWidth(1.0);
}


void Draw_Edge::setBuffer(std::vector<std::array<double, 3>>& vertex_for_render,
	std::vector<unsigned int>& edge_vertex_index)
{
	glBindVertexArray(EE_VAO);
	glBindBuffer(GL_ARRAY_BUFFER, EE_VBO);
	glBufferData(GL_ARRAY_BUFFER, vertex_for_render.size() * sizeof(std::array<double, 3>), vertex_for_render.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EE_EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, edge_vertex_index.size() * sizeof(unsigned int), edge_vertex_index.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);
}