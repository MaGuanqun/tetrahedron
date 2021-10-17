#include"tetrohedron.h"


void Tetrohedron::loadMesh(OriMesh& ori_mesh, double density, Thread* thread)
{
	total_thread_num= std::thread::hardware_concurrency();
	this->thread = thread;
	mesh_struct.thread = thread;
	setMeshStruct(density, ori_mesh);
	mesh_struct.setVertex();
	mesh_struct.initialNormalSize();
	mesh_struct.setThreadIndex(total_thread_num);
	mesh_struct.vertex_for_render = mesh_struct.vertex_position;
	mesh_struct.getRenderNormal();
	genBuffer();
	genShader();
	setBuffer();

}



void Tetrohedron::genShader()
{
	object_shader_front = new Shader("./shader/object.vs", "./shader/object.fs");
	*object_shader_front = Shader("./shader/object.vs", "./shader/object.fs");
	wireframe_shader = new Shader("./shader/wireframe.vs", "./shader/wireframe.fs");
	*wireframe_shader = Shader("./shader/wireframe.vs", "./shader/wireframe.fs");
}

void Tetrohedron::drawWireframe(Camera* camera)
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
		glDrawElements(GL_TRIANGLES, 3*mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}



void Tetrohedron::setBuffer()
{
	if (!mesh_struct.triangle_indices.empty())
	{
		glBindVertexArray(VAO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
		glBufferData(GL_ARRAY_BUFFER, mesh_struct.vertex_for_render.size() * sizeof(std::array<double, 3>), mesh_struct.vertex_for_render[0].data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_struct.triangle_indices.size() * sizeof(std::array<int, 3>), mesh_struct.triangle_indices[0].data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
		glBufferData(GL_ARRAY_BUFFER, mesh_struct.vertex_norm_for_render.size() * sizeof(std::array<double, 3>), mesh_struct.vertex_norm_for_render[0].data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindVertexArray(0);
	}
}


void Tetrohedron::setSceneShader(Light& light, Camera* camera, float& far_plane)
{
	object_shader_front->use();
	object_shader_front->setInt("depthMap", 0);
	object_shader_front->setFloat("far_plane", far_plane);
	object_shader_front->setBool("lightIsChosen", true);
	object_shader_front->setVec3("lightPos", camera->position);
	object_shader_front->setVec3("light.ambient", light.ambient);
	object_shader_front->setVec3("light.diffuse", light.diffuse);
	object_shader_front->setVec3("light.specular", light.specular);
}


void Tetrohedron::drawShadow(Camera* camera, Shader* shader)
{
	//setBuffer(cloth_index);
	if (!mesh_struct.triangle_indices.empty()) {
		shader->use();
		shader->setMat4("model", glm::mat4(1.0));
		glBindVertexArray(VAO);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, 3*mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}

void Tetrohedron::simpDraw(Camera* camera, Shader* shader)
{
	if (!mesh_struct.triangle_indices.empty()) {
		shader->use();
		shader->setMat4("model", glm::mat4(1.0));
		shader->setMat4("projection", camera->GetProjectMatrix());
		shader->setMat4("view", camera->GetViewMatrix());
		glBindVertexArray(VAO);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, 3*mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}




void Tetrohedron::setMeshStruct(double density, OriMesh& ori_mesh)
{
	tetrohedron_num = ori_mesh.indices.size() / 4;
	material = ori_mesh.front_material;
	mesh_struct.vertex_position = ori_mesh.vertices;
	this->density = density;
	mesh_struct.indices.resize(ori_mesh.indices.size() / 4);
	memcpy(mesh_struct.indices[0].data(), ori_mesh.indices.data(), 4 * ori_mesh.indices.size());
	mesh_struct.findSurface();
}

void Tetrohedron::draw(Camera* camera)
{
	if (!mesh_struct.triangle_indices.empty()) {
		object_shader_front->use();
		object_shader_front->setVec3("viewPos", camera->position);
		object_shader_front->setBool("lightShadowOn", true);
		object_shader_front->setMat4("projection", camera->GetProjectMatrix());
		object_shader_front->setMat4("view", camera->GetViewMatrix());
		object_shader_front->setMat4("model", glm::mat4(1.0));
		object_shader_front->setFloat("transparence", 1.0);
		object_shader_front->setVec3("material.Kd", glm::vec3(material.Kd[0], material.Kd[1], material.Kd[2]));
		object_shader_front->setVec3("material.Ka", glm::vec3(material.Ka[0], material.Ka[1], material.Ka[2]));
		object_shader_front->setVec3("material.Ks", glm::vec3(material.Ks[0], material.Ks[1], material.Ks[2]));
		object_shader_front->setFloat("material.Ns", material.Ns);
		glBindVertexArray(VAO);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, 3 * mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}