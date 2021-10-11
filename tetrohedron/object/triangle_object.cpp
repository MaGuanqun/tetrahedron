#include"triangle_object.h"


void TriangleObject::drawWireframe(Camera* camera)
{
	//setBuffer(cloth_index);
	if (!mesh_struct.triangle_indices.empty()) {
		wireframe_shader->use();
		wireframe_shader->setMat4("projection", camera->GetProjectMatrix());
		wireframe_shader->setMat4("view", camera->GetViewMatrix());
		wireframe_shader->setMat4("model", glm::mat4(1.0));
		wireframe_shader->setVec3("color", wireframe_color);
		glBindVertexArray(VAO);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_TRIANGLES, mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}

void TriangleObject::setBuffer()
{
	if (!mesh_struct.triangle_indices.empty())
	{
		glBindVertexArray(VAO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
		glBufferData(GL_ARRAY_BUFFER, mesh_struct.vertex_for_render.size() * sizeof(std::array<double, 3>), mesh_struct.vertex_for_render[0].data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_struct.triangle_indices.size() * sizeof(int), &mesh_struct.triangle_indices[0], GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
		glBufferData(GL_ARRAY_BUFFER, mesh_struct.vertex_norm_for_render.size() * sizeof(std::array<double, 3>), mesh_struct.vertex_norm_for_render[0].data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindVertexArray(0);
	}
}


void TriangleObject::drawShadow(Camera* camera, Shader* shader)
{
	//setBuffer(cloth_index);
	if (!mesh_struct.triangle_indices.empty()) {
		shader->use();
		shader->setMat4("model", glm::mat4(1.0));
		glBindVertexArray(VAO);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}

void TriangleObject::simpDraw(Camera* camera, Shader* shader)
{
	if (!mesh_struct.triangle_indices.empty()) {
		shader->use();
		shader->setMat4("model", glm::mat4(1.0));
		shader->setMat4("projection", camera->GetProjectMatrix());
		shader->setMat4("view", camera->GetViewMatrix());
		glBindVertexArray(VAO);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}


void TriangleObject::setMaterial(OriMesh& ori_mesh)
{
	material.back_material = ori_mesh.back_material;
	material.front_material = ori_mesh.front_material;
}