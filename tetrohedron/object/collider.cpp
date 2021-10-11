#include"collider.h"

void Collider::draw(Camera* camera)
{
	object_shader_front->use();
	object_shader_front->setVec3("viewPos", camera->position);
	object_shader_front->setBool("lightShadowOn", true);
	object_shader_front->setMat4("projection", camera->GetProjectMatrix());
	object_shader_front->setMat4("view", camera->GetViewMatrix());
	object_shader_front->setMat4("model", glm::mat4(1.0));
	object_shader_front->setFloat("transparence", 1.0);
	object_shader_front->setVec3("material.Kd", glm::vec3(material.front_material.Kd[0], material.front_material.Kd[1], material.front_material.Kd[2]));
	object_shader_front->setVec3("material.Ka", glm::vec3(material.front_material.Ka[0], material.front_material.Ka[1], material.front_material.Ka[2]));
	object_shader_front->setVec3("material.Ks", glm::vec3(material.front_material.Ks[0], material.front_material.Ks[1], material.front_material.Ks[2]));
	object_shader_front->setFloat("material.Ns", material.front_material.Ns);

	glBindVertexArray(VAO);
	glPolygonMode(GL_FRONT, GL_FILL);
	glCullFace(GL_BACK);
	glDrawElements(GL_TRIANGLES, mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

}

void Collider::setSceneShader(Light& light, Camera* camera, float& far_plane)
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


void Collider::loadMesh(OriMesh& ori_mesh, Thread* thread)
{
	total_thread_num = std::thread::hardware_concurrency();
	setMeshStruct(ori_mesh);
	mesh_struct.thread = thread;
	mesh_struct.initialNormalSize();
	mesh_struct.setVertex();
	mesh_struct.setFace();
	mesh_struct.setThreadIndex(total_thread_num);
	mesh_struct.vertex_for_render = mesh_struct.vertex_position;
	mesh_struct.getRenderNormal();
	mesh_struct.initialInfo();
	genBuffer();
	setBuffer();
	genShader();
	aabb.resize(mesh_struct.vertices.size());
	ori_vertices = mesh_struct.vertex_position;
}

void Collider::setMeshStruct(OriMesh& ori_mesh)
{
	setMaterial(ori_mesh);
	mesh_struct.vertex_position = ori_mesh.vertices;
	if (!ori_mesh.indices.empty()) {
		mesh_struct.triangle_indices = ori_mesh.indices;
	}
	this->density = density;
}
