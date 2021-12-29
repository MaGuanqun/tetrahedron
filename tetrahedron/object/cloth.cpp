#include"cloth.h"
#include"../mesh_struct/triangle_mesh_struct.h"
#include<thread>

void Cloth::draw(Camera* camera, Shader* object_shader_front)
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
		glDrawElements(GL_TRIANGLES, 3*mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		//glBindVertexArray(0);
		//object_shader_back->use();
		//object_shader_back->setVec3("viewPos", camera->position);
		//object_shader_back->setBool("lightShadowOn", true);
		//object_shader_back->setMat4("projection", camera->GetProjectMatrix());
		//object_shader_back->setMat4("view", camera->GetViewMatrix());
		//object_shader_back->setMat4("model", glm::mat4(1.0));
		//object_shader_front->setFloat("transparence", 1.0);
		object_shader_front->setVec3("material.Kd", glm::vec3(material.back_material.Kd[0], material.back_material.Kd[1], material.back_material.Kd[2]));
		object_shader_front->setVec3("material.Ka", glm::vec3(material.back_material.Ka[0], material.back_material.Ka[1], material.back_material.Ka[2]));
		object_shader_front->setVec3("material.Ks", glm::vec3(material.back_material.Ks[0], material.back_material.Ks[1], material.back_material.Ks[2]));
		object_shader_front->setFloat("material.Ns", material.back_material.Ns);
		glBindVertexArray(VAO);
		glPolygonMode(GL_BACK, GL_FILL);
		glCullFace(GL_FRONT);
		glDrawElements(GL_TRIANGLES, 3*mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}


void Cloth::setSceneShader(Light& light, Camera* camera, float& far_plane, Shader* object_shader_front)
{
	object_shader_front->use();
	object_shader_front->setInt("depthMap", 0);
	object_shader_front->setFloat("far_plane", far_plane);
	object_shader_front->setBool("lightIsChosen", true);
	object_shader_front->setVec3("lightPos", camera->position);
	object_shader_front->setVec3("light.ambient", light.ambient);
	object_shader_front->setVec3("light.diffuse", light.diffuse);
	object_shader_front->setVec3("light.specular", light.specular);
	//object_shader_back->use();
	//object_shader_back->setInt("depthMap", 0);
	//object_shader_back->setFloat("far_plane", far_plane);
	//object_shader_back->setBool("lightIsChosen", true);
	//object_shader_back->setVec3("lightPos", camera->position);
	//object_shader_back->setVec3("light.ambient", light.ambient);
	//object_shader_back->setVec3("light.diffuse", light.diffuse);
	//object_shader_back->setVec3("light.specular", light.specular);
}

void Cloth::loadMesh(OriMesh& ori_mesh, double density, Thread* thread)
{
	total_thread_num = std::thread::hardware_concurrency();

	setMeshStruct(density, ori_mesh);
	this->thread = thread;
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
	setArea();

	PC_radius.resize(4);
	setMass(density);
	mesh_struct.setAnchorPosition();

	hash_index_for_vertex.resize(mesh_struct.vertices.size());
	for (int i = 0; i < hash_index_for_vertex.size(); ++i) {
		hash_index_for_vertex[i].reserve(32);
	}
	hash_index_for_edge.resize(mesh_struct.edges.size());
	for (int i = 0; i < hash_index_for_edge.size(); ++i) {
		hash_index_for_edge[i].reserve(32);
	}
	update_stiffness_iteration_number.resize(mesh_struct.vertices.size());

	triangle_AABB.resize(mesh_struct.faces.size());
	edge_AABB.resize(mesh_struct.edges.size());
	vertex_AABB.resize(mesh_struct.vertices.size());

	setRepresentativePrimitve();
}

void Cloth::setMass(double density)
{
	mass = 0.0;
	if (!mesh_struct.faces.empty()) {
		for (int i = 0; i < mesh_struct.faces.size(); ++i) {
			mass += mesh_struct.faces[i].area * density;
		}
		std::cout << mesh_struct.faces[0].area + mesh_struct.faces[1].area << std::endl;
	}
	else {
		mass = 1.25;
	}
}

void Cloth::setMeshStruct(double density, OriMesh& ori_mesh)
{
	setMaterial(ori_mesh);
	mesh_struct.vertex_position = ori_mesh.vertices;	
	if (!ori_mesh.indices.empty()) {
		mesh_struct.triangle_indices.resize(ori_mesh.indices.size() / 3);
		memcpy(mesh_struct.triangle_indices[0].data(), ori_mesh.indices.data(),12* mesh_struct.triangle_indices.size());
	}
	this->density = density;
	setAnchor();
}


void Cloth::setArea()
{
	double v01[3], v21[3], v_cross[3];
	int* triangle_indices;
	for (int i = 0; i < mesh_struct.triangle_indices.size(); ++i) {
		triangle_indices = mesh_struct.triangle_indices[i].data();
		SUB(v01, mesh_struct.vertex_position[triangle_indices[0]], mesh_struct.vertex_position[triangle_indices[1]]);
		SUB(v21, mesh_struct.vertex_position[triangle_indices[2]], mesh_struct.vertex_position[triangle_indices[1]]);
		CROSS(v_cross, v21, v01);
		mesh_struct.faces[i].area = 0.5 * sqrt(DOT(v_cross, v_cross));
	}
}


void Cloth::recordInitialMesh(SingleClothInfo& single_cloth_info_ref)
{
	ori_vertices = mesh_struct.vertex_position;
	this->single_cloth_info_ref = single_cloth_info_ref;
	length_stiffness.resize(mesh_struct.edges.size());
	collision_stiffness.resize(4);
	for (int i = 0; i < 4; ++i) {
		collision_stiffness[i].resize(mesh_struct.vertices.size(), single_cloth_info_ref.collision_stiffness[i]);
	}	
	memcpy(collision_stiffness_initial, single_cloth_info_ref.collision_stiffness, 32);
	bend_stiffness = single_cloth_info_ref.bending_stiffness;
	position_stiffness = single_cloth_info_ref.position_stiffness;
	std::fill(length_stiffness.begin(), length_stiffness.end(), single_cloth_info_ref.length_stiffness);	
	collision_stiffness_time_step_starts = collision_stiffness;
	std::array<double, 4> collision_stiff_indicator;
	for (int i = 0; i < 4; ++i) {
		collision_stiff_indicator[i] = collision_stiffness_update_indicator * single_cloth_info_ref.collision_stiffness[i];
	}
	collision_stiffness_time_step_starts_indicator.resize(4);
	for (int i = 0; i < 4; ++i) {
		collision_stiffness_time_step_starts_indicator[i].resize(mesh_struct.vertices.size(), collision_stiff_indicator[i]);
	}
	
}

void Cloth::initial()
{
	mesh_struct.vertex_position = ori_vertices;
	mesh_struct.vertex_for_render = ori_vertices;
	mesh_struct.getRenderNormal();
	mesh_struct.getNormal();
	for (int i = 0; i < mesh_struct.anchor_vertex.size(); ++i) {
		mesh_struct.anchor_position[i] = mesh_struct.vertex_position[mesh_struct.anchor_vertex[i]];
	}
	std::fill(length_stiffness.begin(), length_stiffness.end(), single_cloth_info_ref.length_stiffness);
	bend_stiffness = single_cloth_info_ref.bending_stiffness;	
	std::array<double, 4> collision_stiff_indicator;
	for (int i = 0; i < 4; ++i) {
		collision_stiff_indicator[i] = collision_stiffness_update_indicator * single_cloth_info_ref.collision_stiffness[i];
	}
	for (int i = 0; i < 4; ++i) {
		std::fill(collision_stiffness_time_step_starts[i].begin(), collision_stiffness_time_step_starts[i].end(), single_cloth_info_ref.collision_stiffness[i]);
		std::fill(collision_stiffness_time_step_starts_indicator[i].begin(), collision_stiffness_time_step_starts_indicator[i].end(), collision_stiff_indicator[i]);
	}
}

void Cloth::initialMouseChosenVertex()
{
	coe_neighbor_vertex_force.clear();
	neighbor_vertex.clear();
}

void Cloth::setAnchor()
{
	//mesh_struct.anchor_vertex.push_back(100);
	//mesh_struct.anchor_vertex.push_back(0);
	//mesh_struct.anchor_vertex.push_back(15 * 30 - 1);
	//mesh_struct.anchor_vertex.push_back(15 * 29);
	mesh_struct.anchor_vertex.push_back(121 * 121 - 1);
	mesh_struct.anchor_vertex.push_back(120 * 121);
}

//VERTEX_AABB
void Cloth::getVertexAABBPerThread(int thread_No)
{
	std::vector<std::array<double, 3>>* vertex_render=&mesh_struct.vertex_for_render;
	std::vector<std::array<double, 3>>* vertex=&mesh_struct.vertex_position;
	int index_end= mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	for (int i = mesh_struct.vertex_index_begin_per_thread[thread_No]; i < index_end; ++i) {
		vertex_AABB[i].obtainAABB((*vertex_render)[i].data(), (*vertex)[i].data(), tolerance);	// 
	}
}

//EDGE_TRIANGLE_AABB
void Cloth::getEdgeTriangleAABBPerThread(int thread_No)
{
	MeshStruct::Edge* edge = mesh_struct.edges.data();
	int* vertex_index;
	int index_end = mesh_struct.edge_index_begin_per_thread[thread_No + 1];
	for (int i = mesh_struct.edge_index_begin_per_thread[thread_No]; i < index_end; ++i) {
		vertex_index = edge[i].vertex;
		getAABB(edge_AABB[i], vertex_AABB[vertex_index[0]], vertex_AABB[vertex_index[1]]);
	}
	std::array<int, 3>* face = mesh_struct.triangle_indices.data();
	index_end = mesh_struct.face_index_begin_per_thread[thread_No + 1];
	for (int i = mesh_struct.face_index_begin_per_thread[thread_No]; i < index_end; ++i) {
		vertex_index = face[i].data();
		getAABB(triangle_AABB[i], vertex_AABB[vertex_index[0]], vertex_AABB[vertex_index[1]], vertex_AABB[vertex_index[2]]);
	}
}

//TRIANGLE_AABB
//void Cloth::getTriangleAABBPerThread(int thread_No)
//{
//	std::vector<std::array<int,3>>* face = &mesh_struct.triangle_indices;
//	int* vertex_index;
//	int index_end = mesh_struct.face_index_begin_per_thread[thread_No + 1];
//	for (int i = mesh_struct.face_index_begin_per_thread[thread_No]; i < index_end; ++i) {
//		vertex_index = (*face)[i].data();
//		getAABB(triangle_AABB[i], vertex_AABB[vertex_index[0]], vertex_AABB[vertex_index[1]], vertex_AABB[vertex_index[2]]);
//	}
//}


void Cloth::obtainAABB()
{
	thread->assignTask(this, VERTEX_AABB);
	//thread->assignTask(this, EDGE_AABB);
	thread->assignTask(this, EDGE_TRIANGLE_AABB);
}


void Cloth::setRepresentativePrimitve()
{
	representative_vertex_num.resize(mesh_struct.faces.size(), 0);
	representative_edge_num.resize(mesh_struct.faces.size(), 0);
	setRepresentativeVertex(mesh_struct.faces, mesh_struct.vertices);
	setRepresentativeEdge(mesh_struct.faces, mesh_struct.edges);
	vertex_from_rep_triangle_index.resize(mesh_struct.vertex_position.size(),-1);
	edge_from_rep_triangle_index.resize(mesh_struct.edges.size(),-1);
	for (int i = 0; i < mesh_struct.faces.size(); ++i) {
		for (int j = 0; j < representative_vertex_num[i]; ++j) {
			vertex_from_rep_triangle_index[mesh_struct.faces[i].vertex[j]] = i;
		}
		for (int j = 0; j < representative_edge_num[i]; ++j) {
			edge_from_rep_triangle_index[mesh_struct.faces[i].edge[j]] = i;
		}
	}
}

void Cloth::setRepresentativeVertex(std::vector<MeshStruct::Face>& face, std::vector<MeshStruct::Vertex>& vertex)
{
	int count;
	bool in_this_triangle[3];
	std::vector<bool> is_used(vertex.size(), false);
	for (int i = 0; i < face.size(); ++i) {
		count = 0;
		memset(in_this_triangle, 0, 3);
		for (int j = 0; j < 3; ++j) {
			if (!is_used[face[i].vertex[j]]) {
				count++;
				is_used[face[i].vertex[j]] = true;
				in_this_triangle[j] = true;
			}		
		}
		representative_vertex_num[i] = count;
		setOrder(in_this_triangle, count, face[i].vertex);
	}
}

void Cloth::setRepresentativeEdge(std::vector<MeshStruct::Face>& face, std::vector<MeshStruct::Edge>& edge)
{
	int count;
	bool in_this_triangle[3];
	std::vector<bool> is_used(edge.size(), false);
	for (int i = 0; i < face.size(); ++i) {
		count = 0;
		memset(in_this_triangle, 0, 3);
		for (int j = 0; j < 3; ++j) {
			if (!is_used[face[i].edge[j]]) {
				count++;
				is_used[face[i].edge[j]] = true;
				in_this_triangle[j] = true;
			}
		}
		representative_edge_num[i] = count;
		setOrderEdge(in_this_triangle, count, face[i].edge.data());
	}
}


void Cloth::setOrder(bool* in_this_triangle, int count, int* index)
{
	if (count == 1) {
		if (in_this_triangle[1]) {
			int temp = index[0];
			index[0] = index[1];
			index[1] = index[2];
			index[2] = temp;
		}
		else if (in_this_triangle[2]) {
			int temp = index[0];
			index[0] = index[2];
			index[2] = index[1];
			index[1] = temp;
		}
	}
	else if (count == 2) {
		if (!in_this_triangle[0]) {
			int temp = index[0];
			index[0] = index[1];
			index[1] = index[2];
			index[2] = temp;
		}
		else if (!in_this_triangle[1]) {
			int temp = index[0];
			index[0] = index[2];
			index[2] = index[1];
			index[1] = temp;
		}
	}
}

void Cloth::setOrderEdge(bool* in_this_triangle, int count, int* index)
{
	if (count == 1) {
		if (in_this_triangle[1]) {
			int temp = index[0];
			index[0] = index[1];
			index[1] = temp;
		}
		else if (in_this_triangle[2]) {
			int temp = index[0];
			index[0] = index[2];
			index[2] = temp;
		}
	}
	else if (count == 2) {
		if (!in_this_triangle[0]) {
			int temp = index[0];
			index[0] = index[2];
			index[2] = temp;
		}
		else if (!in_this_triangle[1]) {
			int temp = index[2];
			index[2] = index[1];
			index[1] = temp;
		}
	}
}

void Cloth::setTolerance(double* tolerance_ratio, double ave_edge_length)
{
	for (int i = 0; i < PC_radius.size(); ++i) {
		PC_radius[i] = tolerance_ratio[i] * ave_edge_length;
	}
	tolerance = 0.5*ave_edge_length;
	//tolerance = tolerance_ratio[SELF_POINT_TRIANGLE] * ave_edge_length;
}


void Cloth::findAllNeighborVertex(int face_index, double cursor_pos[3], double average_edge_length)
{
	std::vector<bool>is_vertex_used(mesh_struct.vertices.size(), false);
	neighbor_vertex.clear();
	neighbor_vertex.push_back(mesh_struct.triangle_indices[face_index][0]);
	is_vertex_used[mesh_struct.triangle_indices[face_index][0]] = true;
	findNeighborVertex(mesh_struct.triangle_indices[face_index][0], 0, is_vertex_used);
	coe_neighbor_vertex_force.clear();
	coe_neighbor_vertex_force.resize(neighbor_vertex.size());
	double vec3[3];
	for (int i = 0; i < neighbor_vertex.size(); i++) {
		SUB(vec3, cursor_pos, mesh_struct.vertex_for_render[neighbor_vertex[i]]);
		coe_neighbor_vertex_force[i] = gaussian(sqrt(DOT(vec3, vec3)) / average_edge_length, 5.0);
	}
	//std::cout <<"add "<< neighbor_vertex.size() << std::endl;
}

void Cloth::findNeighborVertex(int vertex_index, int recursion_deepth, std::vector<bool>&is_vertex_used)
{
	if (recursion_deepth > 1)
		return;
	for (int i = 0; i < mesh_struct.vertices[vertex_index].neighbor_vertex.size(); ++i) {
		if (!is_vertex_used[mesh_struct.vertices[vertex_index].neighbor_vertex[i]]) {
			neighbor_vertex.push_back(mesh_struct.vertices[vertex_index].neighbor_vertex[i]);
			is_vertex_used[mesh_struct.vertices[vertex_index].neighbor_vertex[i]] = true;
			findNeighborVertex(mesh_struct.vertices[vertex_index].neighbor_vertex[i], recursion_deepth + 1, is_vertex_used);
		}
	}
}


void Cloth::initialNeighborPrimitiveRecording(int cloth_num, int tetrahedron_num, int collider_num, bool use_BVH)
{
	triangle_neighbor_cloth_triangle.resize(mesh_struct.triangle_indices.size());
	for (int i = 0; i < mesh_struct.triangle_indices.size(); ++i) {
		triangle_neighbor_cloth_triangle[i].resize(cloth_num);
		for (int j = 0; j < cloth_num; ++j) {
			triangle_neighbor_cloth_triangle[i][j].reserve(10);
		}
	}

	vertex_neighbor_cloth_traingle.resize(mesh_struct.vertex_position.size());
	collide_vertex_cloth_triangle.resize(mesh_struct.vertex_position.size());
	for (int i = 0; i < vertex_neighbor_cloth_traingle.size(); ++i) {
		vertex_neighbor_cloth_traingle[i].resize(cloth_num);
		collide_vertex_cloth_triangle[i].resize(cloth_num);
		for (int j = 0; j < cloth_num; ++j) {
			vertex_neighbor_cloth_traingle[i][j].reserve(10);
			collide_vertex_cloth_triangle[i][j].reserve(10);
		}
	}

	edge_neighbor_cloth_edge.resize(mesh_struct.edges.size());
	collide_edge_cloth_edge.resize(mesh_struct.edges.size());
	for (int i = 0; i < edge_neighbor_cloth_edge.size(); ++i) {
		edge_neighbor_cloth_edge[i].resize(cloth_num);
		collide_edge_cloth_edge[i].resize(cloth_num);
		for (int j = 0; j < cloth_num; ++j) {
			edge_neighbor_cloth_edge[i][j].reserve(10);
			collide_edge_cloth_edge[i][j].reserve(10);
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
		vertex_neighbor_collider_triangle.resize(mesh_struct.vertex_position.size());
		collide_vertex_collider_triangle.resize(mesh_struct.vertex_position.size());
		for (int i = 0; i < vertex_neighbor_cloth_traingle.size(); ++i) {
			vertex_neighbor_collider_triangle[i].resize(collider_num);
			collide_vertex_collider_triangle[i].resize(collider_num);
			for (int j = 0; j < collider_num; ++j) {
				vertex_neighbor_collider_triangle[i][j].reserve(10);
				collide_vertex_collider_triangle[i][j].reserve(10);
			}
		}
	}	
}

//void Cloth::test() //ttest representative triangle
//{
//	std::vector<bool> is_used(mesh_struct.vertices.size(), false);
//	for (int i = 0; i < mesh_struct.faces.size(); ++i) {
//		for (int j = 0; j < representative_vertex_num[i]; ++j) {
//			if (!is_used[mesh_struct.faces[i].vertex[j]]) {
//				is_used[mesh_struct.faces[i].vertex[j]] = true;
//			}
//			else {
//				std::cout << "false " << i << std::endl;
//			}
//		}
//	}
//	for (int i = 0; i < mesh_struct.vertices.size(); ++i) {
//		if (!is_used[i]) {
//			std::cout << "false_" << std::endl;
//		}
//	}
//	std::vector<bool> is_used_(mesh_struct.edges.size(), false);
//	for (int i = 0; i < mesh_struct.faces.size(); ++i) {
//		for (int j = 0; j < representative_edge_num[i]; ++j) {
//			if (!is_used_[mesh_struct.faces[i].edge[j]]) {
//				is_used_[mesh_struct.faces[i].edge[j]] = true;
//			}
//			else {
//				std::cout << "false " << i << std::endl;
//			}
//		}
//	}
//	for (int i = 0; i < mesh_struct.edges.size(); ++i) {
//		if (!is_used_[i]) {
//			std::cout << "false_edge" << std::endl;
//		}
//	}
//	std::cout << "test " << std::endl;
//}