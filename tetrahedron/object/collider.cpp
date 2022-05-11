#include"collider.h"

void Collider::draw(Camera* camera, Shader* object_shader_front)
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
	glDrawElements(GL_TRIANGLES, 3 * mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

}

void Collider::drawOriPos(Camera* camera, Shader* object_shader_front)
{
	object_shader_front->use();
	object_shader_front->setVec3("viewPos", camera->position);
	object_shader_front->setBool("lightShadowOn", true);
	object_shader_front->setMat4("projection", camera->GetProjectMatrix());
	object_shader_front->setMat4("view", camera->GetViewMatrix());
	object_shader_front->setMat4("model", glm::mat4(1.0));
	object_shader_front->setFloat("transparence", 0.4);
	object_shader_front->setVec3("material.Kd", glm::vec3(material.front_material.Kd[0], material.front_material.Kd[1], material.front_material.Kd[2]));
	object_shader_front->setVec3("material.Ka", glm::vec3(material.front_material.Ka[0], material.front_material.Ka[1], material.front_material.Ka[2]));
	object_shader_front->setVec3("material.Ks", glm::vec3(material.front_material.Ks[0], material.front_material.Ks[1], material.front_material.Ks[2]));
	object_shader_front->setFloat("material.Ns", material.front_material.Ns);

	glBindVertexArray(VAO_ori);
	glPolygonMode(GL_FRONT, GL_FILL);
	glCullFace(GL_BACK);
	glDrawElements(GL_TRIANGLES, 3 * mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

}


void Collider::setSceneShader(Light& light, Camera* camera, float& far_plane, Shader* object_shader_front)
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
	this->thread = thread;
	total_thread_num = std::thread::hardware_concurrency();
	obj_aabb_per_thread.resize(total_thread_num);
	current_obj_pos_aabb_per_thread.resize(total_thread_num);
	setMeshStruct(ori_mesh);
	mesh_struct.thread = thread;
	mesh_struct.initialNormalSize();
	mesh_struct.setVertex();
	mesh_struct.setFace();
	mesh_struct.setEdge();
	mesh_struct.setThreadIndex(total_thread_num);
	mesh_struct.vertex_for_render = mesh_struct.vertex_position;
	mesh_struct.getRenderNormal();
	mesh_struct.getNormal();
	mesh_struct.initialInfo();
	genBuffer();
	setBuffer();
	ori_vertices = mesh_struct.vertex_position;
	initialHashAABB();
	setRepresentativePrimitve();
	obtainAABBMoveRadius();
}

void Collider::setMeshStruct(OriMesh& ori_mesh)
{
	setMaterial(ori_mesh);
	mesh_struct.vertex_position = ori_mesh.vertices;
	if (!ori_mesh.indices.empty()) {
		mesh_struct.triangle_indices.resize(ori_mesh.indices.size() / 3);
		memcpy(mesh_struct.triangle_indices[0].data(), ori_mesh.indices.data(), 12 * mesh_struct.triangle_indices.size());
	}
	this->density = density;
}

void Collider::obtainAABB(bool has_tolerace)
{
	if (has_tolerace) {
		thread->assignTask(this, VERTEX_AABB);
	}
	else {
		thread->assignTask(this, VERTEX_AABB_WITHOUT_TOLERANCE);
	}
	thread->assignTask(this, EDGE_TRIANGLE_AABB);

	combineObjAABB();
}


void Collider::reset()
{
	memset(rotation_matrix, 0, 72);
	rotation_matrix[0] = rotation_matrix[4] = rotation_matrix[8] = 1.0;
	mesh_struct.vertex_position = ori_vertices;
	mesh_struct.vertex_for_render = ori_vertices;
	mesh_struct.getRenderNormal();
	mesh_struct.getNormal();
	for (int i = 0; i < mesh_struct.anchor_vertex.size(); ++i) {
		mesh_struct.anchor_position[i] = mesh_struct.vertex_position[mesh_struct.anchor_vertex[i]];
	}
	obtainAABBMoveRadius();
}


void Collider::obtainAABBMoveRadius()
{
	thread->assignTask(this, CURRENT_AABB);
	combineCurrentAABBMoveRadius();
}


////EDGE_TRIANGLE_AABB
//void Collider::getTriangleAABBPerThread(int thread_No)
//{
//	std::array<int,3>* face = mesh_struct.triangle_indices.data();
//	int* vertex_index;
//	unsigned int end = mesh_struct.face_index_begin_per_thread[thread_No + 1];
//	for (unsigned int i = mesh_struct.face_index_begin_per_thread[thread_No]; i < end; ++i) {
//		vertex_index = face[i].data();
//		getAABB(triangle_AABB[i].data(), vertex_AABB[vertex_index[0]].data(), 
//			vertex_AABB[vertex_index[1]].data(), vertex_AABB[vertex_index[2]].data());
//	}
//}

void Collider::setTolerance(double* tolerance_ratio, double ave_edge_length)
{
	tolerance = tolerance_ratio[AABB_SELF_PT] * ave_edge_length;
	//tolerance = 0.5 * ave_edge_length;
}

void Collider::initialNeighborPrimitiveRecording(int cloth_num, int tetrahedron_num)
{
	triangle_neighbor_obj_triangle.resize(mesh_struct.triangle_indices.size());
	triangle_neighbor_obj_vertex.resize(mesh_struct.triangle_indices.size());
	collider_triangle_obj_vertex.resize(mesh_struct.triangle_indices.size());
	int obj_num = cloth_num + tetrahedron_num;
	for (int i = 0; i < mesh_struct.triangle_indices.size(); ++i) {
		triangle_neighbor_obj_triangle[i].resize(obj_num);
		triangle_neighbor_obj_vertex[i].resize(obj_num);
		collider_triangle_obj_vertex[i].resize(obj_num);
		for (int j = 0; j < obj_num; ++j) {
			triangle_neighbor_obj_triangle[i][j].reserve(10);
			triangle_neighbor_obj_vertex[i][j].reserve(10);
			collider_triangle_obj_vertex[i][j].reserve(10);
		}
	}
}