#include"tetrahedron.h"
#include<vector>

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
	//mesh_struct.addArounVertex();
	mesh_struct.setThreadIndex(total_thread_num);
	mesh_struct.vertex_for_render = mesh_struct.vertex_position;
	mesh_struct.getRenderNormal();
	mesh_struct.getNormal();
	mesh_struct.prepareForDeformationGradient();
	mass = mesh_struct.setVolumeMass(density);
	mesh_struct.recordTetIndexForSurfaceIndex();
	genBuffer();
	setBuffer();	
	setRepresentativePrimitve();
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


void Tetrahedron::setRepresentativePrimitve()
{
	representative_vertex_num.resize(mesh_struct.faces.size(), 0);
	representative_edge_num.resize(mesh_struct.faces.size(), 0);
	setRepresentativeVertex(mesh_struct.surface_triangle_index_in_order, mesh_struct.vertices);
	mesh_struct.setVertexIndexOnSurfaceEdgeTriangle();
	setRepresentativeEdge(mesh_struct.faces, mesh_struct.edges);
	surface_vertex_from_rep_triangle_index.resize(mesh_struct.vertex_index_on_sureface.size(), -1);
	edge_from_rep_triangle_index.resize(mesh_struct.edges.size(), -1);
	
	for (int i = 0; i < mesh_struct.faces.size(); ++i) {
		for (int j = 0; j < representative_vertex_num[i]; ++j) {
			surface_vertex_from_rep_triangle_index[mesh_struct.surface_triangle_index_in_order[i][j]] = i;
		}
		for (int j = 0; j < representative_edge_num[i]; ++j) {
			edge_from_rep_triangle_index[mesh_struct.faces[i].edge[j]] = i;
		}
	}
}

void Tetrahedron::obtainAABB()
{
	thread->assignTask(this, VERTEX_AABB);
	//thread->assignTask(this, EDGE_AABB);
	thread->assignTask(this, EDGE_TRIANGLE_AABB);
}


//VERTEX_AABB
void Tetrahedron::getVertexAABBPerThread(int thread_No)
{
	std::vector<std::array<double, 3>>* vertex_render = &mesh_struct.vertex_for_render;
	std::vector<std::array<double, 3>>* vertex = &mesh_struct.vertex_position;
	int index_end = mesh_struct.vertex_index_on_surface_begin_per_thread[thread_No + 1];
	int* vertex_index_on_surface = mesh_struct.vertex_index_on_sureface.data();
	for (int i = mesh_struct.vertex_index_on_surface_begin_per_thread[thread_No]; i < index_end; ++i) {
		vertex_AABB[i].obtainAABB((*vertex_render)[vertex_index_on_surface[i]].data(), (*vertex)[vertex_index_on_surface[i]].data(), tolerance);	// 
	}
}

//EDGE_TRIANGLE_AABB
void Tetrahedron::getEdgeTriangleAABBPerThread(int thread_No)
{
	std::array<int,2>* edge = mesh_struct.edge_vertex_index_on_surface.data();
	int* vertex_index;
	int index_end = mesh_struct.edge_index_begin_per_thread[thread_No + 1];
	for (int i = mesh_struct.edge_index_begin_per_thread[thread_No]; i < index_end; ++i) {
		vertex_index = edge[i].data();
		getAABB(edge_AABB[i], vertex_AABB[vertex_index[0]], vertex_AABB[vertex_index[1]]);
	}
	std::array<int, 3>* face = mesh_struct.surface_triangle_index_in_order.data();
	index_end = mesh_struct.face_index_begin_per_thread[thread_No + 1];
	for (int i = mesh_struct.face_index_begin_per_thread[thread_No]; i < index_end; ++i) {
		vertex_index = face[i].data();
		getAABB(triangle_AABB[i], vertex_AABB[vertex_index[0]], vertex_AABB[vertex_index[1]], vertex_AABB[vertex_index[2]]);
	}
}



void Tetrahedron::initialHashAABB()
{
	PC_radius.resize(4);
	hash_index_for_vertex.resize(mesh_struct.vertex_index_on_sureface.size());
	for (int i = 0; i < hash_index_for_vertex.size(); ++i) {
		hash_index_for_vertex[i].reserve(32);
	}
	hash_index_for_edge.resize(mesh_struct.edges.size());
	for (int i = 0; i < hash_index_for_edge.size(); ++i) {
		hash_index_for_edge[i].reserve(32);
	}
	triangle_AABB.resize(mesh_struct.faces.size());
	edge_AABB.resize(mesh_struct.edges.size());
	vertex_AABB.resize(mesh_struct.vertex_index_on_sureface.size());
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
	memcpy(sigma_limit, single_tetrahedron_info_ref.sigma_limit, 16);
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


void Tetrahedron::findAllNeighborVertex(int face_index, double cursor_pos[3], double average_edge_length)
{
	std::vector<bool>is_vertex_used(mesh_struct.vertices.size(), false);
	neighbor_vertex.clear();
	neighbor_vertex.push_back(mesh_struct.triangle_indices[face_index][0]);
	is_vertex_used[mesh_struct.triangle_indices[face_index][0]] = true;
	findNeighborVertex(mesh_struct.triangle_indices[face_index][0], 0, is_vertex_used);
	//findInnerVertex(is_vertex_used);
	coe_neighbor_vertex_force.clear();
	coe_neighbor_vertex_force.resize(neighbor_vertex.size());
	double vec3[3];
	for (int i = 0; i < neighbor_vertex.size(); i++) {
		SUB(vec3, cursor_pos, mesh_struct.vertex_for_render[neighbor_vertex[i]]);
		coe_neighbor_vertex_force[i] = gaussian(sqrt(DOT(vec3, vec3)) / average_edge_length, 5.0);
	}
}

void Tetrahedron::findNeighborVertex(int vertex_index, int recursion_deepth, std::vector<bool>& is_vertex_used)
{
	if (recursion_deepth > 2)
		return;
	for (int i = 0; i < mesh_struct.vertices[vertex_index].neighbor_vertex.size(); ++i) {
		if (!is_vertex_used[mesh_struct.vertices[vertex_index].neighbor_vertex[i]]) {
			neighbor_vertex.push_back(mesh_struct.vertices[vertex_index].neighbor_vertex[i]);
			is_vertex_used[mesh_struct.vertices[vertex_index].neighbor_vertex[i]] = true;
			findNeighborVertex(mesh_struct.vertices[vertex_index].neighbor_vertex[i], recursion_deepth + 1, is_vertex_used);
		}
	}
}



void Tetrahedron::setTolerance(double* tolerance_ratio, double ave_edge_length)
{
	for (int i = 0; i < PC_radius.size(); ++i) {
		PC_radius[i] = tolerance_ratio[i] * ave_edge_length;
	}
	tolerance = 0.5 * ave_edge_length;
	//tolerance = tolerance_ratio[SELF_POINT_TRIANGLE] * ave_edge_length;
}

void Tetrahedron::findInnerVertex(std::vector<bool>& is_vertex_used)
{
	int size = neighbor_vertex.size();
	std::vector<int>* vertex_index_;
	for (int i = 0; i < size; ++i) {
		vertex_index_ = &mesh_struct.vertices[neighbor_vertex[i]].tetrahedron;
		for (int j = 0; j < vertex_index_->size(); ++j) {
			if (!is_vertex_used[vertex_index_->data()[j]]) {
				neighbor_vertex.push_back(vertex_index_->data()[j]);
				is_vertex_used[vertex_index_->data()[j]] = true;
			}
		}	
	}
}

void Tetrahedron::initialNeighborPrimitiveRecording(int cloth_num, int tetrahedron_num, int collider_num, bool use_BVH)
{
	int obj_num = cloth_num + tetrahedron_num;
	triangle_neighbor_obj_triangle.resize(mesh_struct.triangle_indices.size());
	for (int i = 0; i < mesh_struct.triangle_indices.size(); ++i) {
		triangle_neighbor_obj_triangle[i].resize(obj_num);
		for (int j = 0; j < obj_num; ++j) {
			triangle_neighbor_obj_triangle[i][j].reserve(10);
		}
	}

	surface_vertex_neighbor_obj_triangle.resize(mesh_struct.vertex_index_on_sureface.size());
	collide_vertex_obj_triangle.resize(mesh_struct.vertex_index_on_sureface.size());
	for (int i = 0; i < surface_vertex_neighbor_obj_triangle.size(); ++i) {
		surface_vertex_neighbor_obj_triangle[i].resize(obj_num);
		collide_vertex_obj_triangle[i].resize(obj_num);
		for (int j = 0; j < obj_num; ++j) {
			surface_vertex_neighbor_obj_triangle[i][j].reserve(10);
			collide_vertex_obj_triangle[i][j].reserve(10);
		}
	}

	edge_neighbor_obj_edge.resize(mesh_struct.edges.size());
	collide_edge_obj_edge.resize(mesh_struct.edges.size());
	for (int i = 0; i < edge_neighbor_obj_edge.size(); ++i) {
		edge_neighbor_obj_edge[i].resize(obj_num);
		collide_edge_obj_edge[i].resize(obj_num);
		for (int j = 0; j < obj_num; ++j) {
			edge_neighbor_obj_edge[i][j].reserve(10);
			collide_edge_obj_edge[i][j].reserve(10);
		}
	}

	if (use_BVH) {
		triangle_neighbor_collider_triangle.resize(mesh_struct.triangle_indices.size());
		for (int i = 0; i < mesh_struct.triangle_indices.size(); ++i) {
			triangle_neighbor_collider_triangle[i].resize(collider_num);
			for (int j = 0; j < collider_num; ++j) {
				triangle_neighbor_collider_triangle[i][j].reserve(10);
			}
		}
		surface_vertex_neighbor_collider_triangle.resize(mesh_struct.vertex_index_on_sureface.size());
		collide_vertex_collider_triangle.resize(mesh_struct.vertex_index_on_sureface.size());
		for (int i = 0; i < surface_vertex_neighbor_collider_triangle.size(); ++i) {
			surface_vertex_neighbor_collider_triangle[i].resize(collider_num);
			collide_vertex_collider_triangle[i].resize(collider_num);
			for (int j = 0; j < collider_num; ++j) {
				surface_vertex_neighbor_collider_triangle[i][j].reserve(10);
				collide_vertex_collider_triangle[i][j].reserve(10);
			}
		}
	}
}