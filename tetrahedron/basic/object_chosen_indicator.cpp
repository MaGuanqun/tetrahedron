#include"object_chosen_indicator.h"

ObjectChosenIndicator::ObjectChosenIndicator()
{
	vertex_num = 50;
	circle_vertices.resize(3);
	for (unsigned int i = 0; i < circle_vertices.size(); ++i) {
		circle_vertices[i].resize(3 * vertex_num);
	}
	genBuffer();
	circle_color.resize(3);
	circle_color[0] = glm::vec3(1.0, 0.1, 0.1);
	circle_color[1] = glm::vec3(0.1, 1.0, 0.1);
	circle_color[2] = glm::vec3(0.1, 0.1, 1.0);

}


void ObjectChosenIndicator::updatePosition(double* AABB)
{
	double center[3] = { 0.5 * (AABB[0] + AABB[3]),0.5 * (AABB[1] + AABB[4]),0.5 * (AABB[2] + AABB[5]) };
	double radius = 0.5 * sqrt((AABB[0] - AABB[3]) * (AABB[0] - AABB[3]) + (AABB[1] - AABB[4]) * (AABB[1] - AABB[4])
		+ (AABB[2] - AABB[5]) * (AABB[2] - AABB[5]));

	for (int i = 0; i < vertex_num; i++)
	{
		circle_vertices[0][3 * i] = center[0]+ radius * sin(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));
		circle_vertices[0][3 * i + 1] = center[1] + radius * cos(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));
		circle_vertices[0][3 * i + 2] = center[2];

		circle_vertices[1][3 * i] = center[0];
		circle_vertices[1][3 * i + 1] = center[1] + radius * sin(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));
		circle_vertices[1][3 * i + 2] = center[2] + radius * cos(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));

		circle_vertices[2][3 * i] = center[0] + radius * cos(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));
		circle_vertices[2][3 * i + 1] = center[1];
		circle_vertices[2][3 * i + 2] = center[2] + radius * sin(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));
	}	
	setBuffer();
}

void ObjectChosenIndicator::setBuffer()
{
	for (unsigned int i = 0; i < 3; ++i) {
		glBindVertexArray(VAO[i]);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[i]);
		glBufferData(GL_ARRAY_BUFFER, circle_vertices[i].size()*sizeof(double), circle_vertices[i].data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}
}


void ObjectChosenIndicator::draw(Shader* shader, Camera* camera)
{
	shader->use();
	shader->setMat4("projection", camera->GetProjectMatrix());
	shader->setMat4("view", camera->GetViewMatrix());
	shader->setMat4("model", glm::mat4(1.0));
	for (int i = 0; i < 3; i++)
	{
		shader->setVec3("color", circle_color[i]);
		glBindVertexArray(VAO[i]);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glLineWidth(3.0f);
		glDrawArrays(GL_LINE_LOOP, 0, vertex_num);
		glBindVertexArray(0);
	}

}

void ObjectChosenIndicator::genBuffer()
{
	for (unsigned int i = 0; i < 3; ++i) {
		glGenVertexArrays(1, &VAO[i]);
		glGenBuffers(1, &VBO[i]);
	}
	
}