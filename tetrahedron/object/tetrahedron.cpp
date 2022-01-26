#include"tetrahedron.h"


void Tetrahedron::loadMesh(OriMesh& ori_mesh, double density, Thread* thread)
{
	total_thread_num= std::thread::hardware_concurrency();
	this->thread = thread;
	mesh_struct.thread = thread;
	setMeshStruct(density, ori_mesh);
	mesh_struct.setVertex();
	mesh_struct.initialNormalSize();
	mesh_struct.setFace();
	mesh_struct.setEdge();
	mesh_struct.addArounVertex();
	mesh_struct.setThreadIndex(total_thread_num);
	mesh_struct.vertex_for_render = mesh_struct.vertex_position;
	mesh_struct.getRenderNormal();
	mesh_struct.getNormal();
	mesh_struct.prepareForDeformationGradient();
	mass = mesh_struct.setVolumeMass(density);
	genBuffer();
	setBuffer();
	initialHashAABB();
}


//void Tetrahedron::genShader()
//{
//	object_shader_front = new Shader("./shader/object_tetrahedron.vs", "./shader/object_tetrahedron.fs", "./shader/object_tetrahedron.gs");
//	*object_shader_front = Shader("./shader/object_tetrahedron.vs", "./shader/object_tetrahedron.fs", "./shader/object_tetrahedron.gs");
//}

void Tetrahedron::drawWireframe(Camera* camera, Shader* wireframe_shader)
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


void Tetrahedron::initialHashAABB()
{
	PC_radius.resize(4);
	hash_index_for_vertex.resize(mesh_struct.vertices.size());
	for (int i = 0; i < hash_index_for_vertex.size(); ++i) {
		hash_index_for_vertex[i].reserve(32);
	}
	hash_index_for_edge.resize(mesh_struct.edges.size());
	for (int i = 0; i < hash_index_for_edge.size(); ++i) {
		hash_index_for_edge[i].reserve(32);
	}
	triangle_AABB.resize(mesh_struct.faces.size());
	edge_AABB.resize(mesh_struct.edges.size());
	vertex_AABB.resize(mesh_struct.vertices.size());
}

void Tetrahedron::setBuffer()
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
		glBufferData(GL_ARRAY_BUFFER, mesh_struct.vertex_normal_for_render.size() * sizeof(std::array<double, 3>), mesh_struct.vertex_normal_for_render[0].data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindVertexArray(0);
	}
}


void Tetrahedron::setSceneShader(Light& light, Camera* camera, float& far_plane, Shader* object_shader_front)
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


void Tetrahedron::drawShadow(Camera* camera, Shader* shader)
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

void Tetrahedron::simpDraw(Camera* camera, Shader* shader)
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




void Tetrahedron::setMeshStruct(double density, OriMesh& ori_mesh)
{
	tetrahedron_num = ori_mesh.indices.size() / 4;
	material = ori_mesh.front_material;
	mesh_struct.vertex_position = ori_mesh.vertices;
	this->density = density;
	mesh_struct.indices.resize(ori_mesh.indices.size() / 4);
	memcpy(mesh_struct.indices[0].data(), ori_mesh.indices.data(), 4 * ori_mesh.indices.size());
	mesh_struct.findSurface();
}

void Tetrahedron::draw(Camera* camera, Shader* object_shader_front)
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

void Tetrahedron::recordInitialMesh(SingleTetrahedronInfo& single_tetrahedron_info_ref)
{
	ori_vertices = mesh_struct.vertex_position;
	this->single_tetrahedron_info_ref = single_tetrahedron_info_ref;
	memcpy(collision_stiffness, single_tetrahedron_info_ref.collision_stiffness, 32);
	position_stiffness = single_tetrahedron_info_ref.position_stiffness;
	volume_preserve_stiffness = single_tetrahedron_info_ref.volume_preserve_stiffness;
	ARAP_stiffness = single_tetrahedron_info_ref.ARAP_stiffness;
}


void Tetrahedron::initial()
{
	mesh_struct.vertex_position = ori_vertices;
	mesh_struct.vertex_for_render = ori_vertices;
	mesh_struct.anchor_vertex.clear();
	mesh_struct.anchor_position.clear();
	mesh_struct.getRenderNormal();
	mesh_struct.getNormal();
	ARAP_stiffness = single_tetrahedron_info_ref.ARAP_stiffness;
	volume_preserve_stiffness = single_tetrahedron_info_ref.volume_preserve_stiffness;
	position_stiffness = single_tetrahedron_info_ref.position_stiffness;
	std::array<double, 4> collision_stiff_indicator;
	memcpy(collision_stiffness, single_tetrahedron_info_ref.collision_stiffness, 16);

}

void Tetrahedron::reset()
{
	mesh_struct.vertex_position = ori_vertices;
	mesh_struct.vertex_for_render = ori_vertices;
	mesh_struct.getRenderNormal();
	mesh_struct.getNormal();
	for (int i = 0; i < mesh_struct.anchor_vertex.size(); ++i) {
		mesh_struct.anchor_position[i] = mesh_struct.vertex_position[mesh_struct.anchor_vertex[i]];
	}
}
