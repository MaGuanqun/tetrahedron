#include"cloth.h"
#include"../mesh_struct/triangle_mesh_struct.h"
#include<thread>

void Cloth::draw(Camera* camera)
{
	if (!mesh_struct.triangle_indices.empty()) {
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
		object_shader_back->use();
		object_shader_back->setVec3("viewPos", camera->position);
		object_shader_back->setBool("lightShadowOn", true);
		object_shader_back->setMat4("projection", camera->GetProjectMatrix());
		object_shader_back->setMat4("view", camera->GetViewMatrix());
		object_shader_back->setMat4("model", glm::mat4(1.0));
		object_shader_front->setFloat("transparence", 1.0);
		object_shader_back->setVec3("material.Kd", glm::vec3(material.back_material.Kd[0], material.back_material.Kd[1], material.back_material.Kd[2]));
		object_shader_back->setVec3("material.Ka", glm::vec3(material.back_material.Ka[0], material.back_material.Ka[1], material.back_material.Ka[2]));
		object_shader_back->setVec3("material.Ks", glm::vec3(material.back_material.Ks[0], material.back_material.Ks[1], material.back_material.Ks[2]));
		object_shader_back->setFloat("material.Ns", material.back_material.Ns);
		glBindVertexArray(VAO);
		glPolygonMode(GL_BACK, GL_FILL);
		glCullFace(GL_FRONT);
		glDrawElements(GL_TRIANGLES, mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}


void Cloth::setSceneShader(Light& light, Camera* camera, float& far_plane)
{
	object_shader_front->use();
	object_shader_front->setInt("depthMap", 0);
	object_shader_front->setFloat("far_plane", far_plane);
	object_shader_front->setBool("lightIsChosen", true);
	object_shader_front->setVec3("lightPos", camera->position);
	object_shader_front->setVec3("light.ambient", light.ambient);
	object_shader_front->setVec3("light.diffuse", light.diffuse);
	object_shader_front->setVec3("light.specular", light.specular);
	object_shader_back->use();
	object_shader_back->setInt("depthMap", 0);
	object_shader_back->setFloat("far_plane", far_plane);
	object_shader_back->setBool("lightIsChosen", true);
	object_shader_back->setVec3("lightPos", camera->position);
	object_shader_back->setVec3("light.ambient", light.ambient);
	object_shader_back->setVec3("light.diffuse", light.diffuse);
	object_shader_back->setVec3("light.specular", light.specular);
}

void Cloth::loadMesh(OriMesh& ori_mesh, double density, Thread* thread)
{
	total_thread_num = std::thread::hardware_concurrency();
	TriangleMeshStruct triangle_mesh_struct;
	mesh_struct = triangle_mesh_struct;

	setMeshStruct(density, ori_mesh);

	mesh_struct.thread = thread;
	mesh_struct.initialNormalSize();
	mesh_struct.setVertex();
	mesh_struct.setFace();
	mesh_struct.setEdge();
	mesh_struct.setThreadIndex(total_thread_num);
	mesh_struct.vertex_for_render = mesh_struct.vertex_position;
	mesh_struct.getRenderNormal();
	mesh_struct.initialInfo();
	genBuffer();
	genShader();
	setBuffer();
	setArea();
	is_vertex_used.resize(mesh_struct.vertices.size(), false);

	//face_need_recal_radius = new bool[mesh_struct.faces.size()];
	//edge_need_recal_radius = new bool[mesh_struct.edges.size()];
	PC_radius.resize(4);
	PC_radius_coe = 0.0001;
	setMass(density);
	mesh_struct.setAnchorPosition();



	aabb.resize(4);
	for (int i = 0; i < 4; ++i) {
		aabb[i].resize(mesh_struct.vertices.size());
	}
	hash_index_for_vertex.resize(mesh_struct.vertices.size());
	for (int i = 0; i < hash_index_for_vertex.size(); ++i) {
		hash_index_for_vertex[i].reserve(32);
	}
	hash_index_for_edge.resize(mesh_struct.edges.size());
	for (int i = 0; i < hash_index_for_edge.size(); ++i) {
		hash_index_for_edge[i].reserve(32);
	}
	update_stiffness_iteration_number.resize(mesh_struct.vertices.size());
}


void Cloth::setMass(double density)
{
	mass = 0.0;
	if (!mesh_struct.faces.empty()) {
		for (int i = 0; i < mesh_struct.faces.size(); ++i) {
			mass += mesh_struct.faces[i].area * density;
		}
	}
	else {
		mass = 1.25;
	}
}

void Cloth::setAnchor(std::vector<int>& anchor_vertex)
{
	mesh_struct.anchor_vertex = anchor_vertex;
	for (int i = 0; i < mesh_struct.anchor_vertex.size(); ++i) {
		mesh_struct.anchor_position.push_back(mesh_struct.vertex_position[mesh_struct.anchor_vertex[i]]);
	}
}

void Cloth::setMeshStruct(double density, OriMesh& ori_mesh)
{
	setMaterial(ori_mesh);
	mesh_struct.vertex_position = ori_mesh.vertices;
	if (!ori_mesh.indices.empty()) {
		mesh_struct.triangle_indices=ori_mesh.indices;
	}
	this->density = density;
}


void Cloth::setArea()
{
	double v01[3], v21[3], v_cross[3];
	for (int i = 0; i < mesh_struct.faces.size(); ++i) {
		SUB(v01, mesh_struct.vertex_position[mesh_struct.triangle_indices[3*i]], mesh_struct.vertex_position[mesh_struct.triangle_indices[3 * i+1]]);
		SUB(v21, mesh_struct.vertex_position[mesh_struct.triangle_indices[3 * i+2]], mesh_struct.vertex_position[mesh_struct.triangle_indices[3 * i+1]]);
		CROSS(v_cross, v21, v01);
		mesh_struct.faces[i].area = 0.5 * sqrt(DOT(v_cross, v_cross));
	}
}


void Cloth::recordInitialMesh(SingleClothInfo& single_cloth_info_ref)
{
	ori_vertices = mesh_struct.vertex_position;
	this->single_cloth_info_ref = single_cloth_info_ref;
	length_stiffness.resize(mesh_struct.edges.size());	
	collision_stiffness.resize(mesh_struct.vertices.size());
	bend_stiffness = single_cloth_info_ref.bending_stiffness;
	position_stiffness = single_cloth_info_ref.position_stiffness;
	std::array<double, 4> collision_stiff;
	std::fill(length_stiffness.begin(), length_stiffness.end(), single_cloth_info_ref.length_stiffness);
	memcpy(collision_stiff.data(), single_cloth_info_ref.collision_stiffness, 32);
	std::fill(collision_stiffness.begin(), collision_stiffness.end(), collision_stiff);
	collision_stiffness_time_step_starts = collision_stiffness;
	std::array<double, 4> collision_stiff_indicator;
	for (int i = 0; i < 4; ++i) {
		collision_stiff_indicator[i] = collision_stiffness_update_indicator * collision_stiff[i];
	}
	collision_stiffness_time_step_starts_indicator.resize(mesh_struct.vertices.size(), collision_stiff_indicator);
}