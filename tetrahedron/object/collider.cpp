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
	glDrawElements(GL_TRIANGLES, 3* mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
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
	setMeshStruct(ori_mesh);
	mesh_struct.thread = thread;
	mesh_struct.initialNormalSize();
	mesh_struct.setVertex();
	mesh_struct.setFace();
	mesh_struct.setThreadIndex(total_thread_num);
	mesh_struct.vertex_for_render = mesh_struct.vertex_position;
	mesh_struct.getRenderNormal();
	mesh_struct.getNormal();
	mesh_struct.initialInfo();
	genBuffer();
	setBuffer();
	ori_vertices = mesh_struct.vertex_position;
	triangle_AABB.resize(mesh_struct.faces.size());
	vertex_AABB.resize(mesh_struct.vertex_position.size());
}

void Collider::setMeshStruct(OriMesh& ori_mesh)
{
	setMaterial(ori_mesh);
	mesh_struct.vertex_position = ori_mesh.vertices;
	if (!ori_mesh.indices.empty()) {
		mesh_struct.triangle_indices.resize(ori_mesh.indices.size() / 3);
		memcpy(mesh_struct.triangle_indices[0].data(), ori_mesh.indices.data(),12* mesh_struct.triangle_indices.size());
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
}


//VERTEX_AABB
//VERTEX_AABB_WITHOUT_TOLERANCE
void Collider::getVertexAABBPerThread(int thread_No, bool has_tolerance)
{
	std::vector<std::array<double, 3>>* vertex_render = &mesh_struct.vertex_for_render;
	std::vector<std::array<double, 3>>* vertex = &mesh_struct.vertex_position;
	unsigned int end = mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	if (has_tolerance) {
		for (unsigned int i = mesh_struct.vertex_index_begin_per_thread[thread_No]; i < end; ++i) {
			AABB::obtainAABB(vertex_AABB[i].data(),(*vertex_render)[i].data(), (*vertex)[i].data(), tolerance);// 
		}
	}
	else {
		for (unsigned int i = mesh_struct.vertex_index_begin_per_thread[thread_No]; i < end; ++i) {
			AABB::obtainAABB(vertex_AABB[i].data(),(*vertex_render)[i].data(), (*vertex)[i].data());// 
		}
	}
}

//EDGE_TRIANGLE_AABB
void Collider::getTriangleAABBPerThread(int thread_No)
{
	std::array<int,3>* face = mesh_struct.triangle_indices.data();
	int* vertex_index;
	unsigned int end = mesh_struct.face_index_begin_per_thread[thread_No + 1];
	for (unsigned int i = mesh_struct.face_index_begin_per_thread[thread_No]; i < end; ++i) {
		vertex_index = face[i].data();
		getAABB(triangle_AABB[i].data(), vertex_AABB[vertex_index[0]].data(), 
			vertex_AABB[vertex_index[1]].data(), vertex_AABB[vertex_index[2]].data());
	}
}

void Collider::setTolerance(double* tolerance_ratio, double ave_edge_length)
{
	tolerance = tolerance_ratio[BODY_POINT_TRIANGLE] * ave_edge_length;
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