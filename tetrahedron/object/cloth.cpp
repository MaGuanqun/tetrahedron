#include"cloth.h"
#include"../mesh_struct/triangle_mesh_struct.h"
#include<thread>


void Cloth::drawOriPos(Camera* camera, Shader* object_shader_front)
{
	if (!mesh_struct.triangle_indices.empty()) {
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

		object_shader_front->setVec3("material.Kd", glm::vec3(material.back_material.Kd[0], material.back_material.Kd[1], material.back_material.Kd[2]));
		object_shader_front->setVec3("material.Ka", glm::vec3(material.back_material.Ka[0], material.back_material.Ka[1], material.back_material.Ka[2]));
		object_shader_front->setVec3("material.Ks", glm::vec3(material.back_material.Ks[0], material.back_material.Ks[1], material.back_material.Ks[2]));
		object_shader_front->setFloat("material.Ns", material.back_material.Ns);
		glBindVertexArray(VAO_ori);
		glPolygonMode(GL_BACK, GL_FILL);
		glCullFace(GL_FRONT);
		glDrawElements(GL_TRIANGLES, 3 * mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}


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
		glDrawElements(GL_TRIANGLES, 3 * mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);

		object_shader_front->setVec3("material.Kd", glm::vec3(material.back_material.Kd[0], material.back_material.Kd[1], material.back_material.Kd[2]));
		object_shader_front->setVec3("material.Ka", glm::vec3(material.back_material.Ka[0], material.back_material.Ka[1], material.back_material.Ka[2]));
		object_shader_front->setVec3("material.Ks", glm::vec3(material.back_material.Ks[0], material.back_material.Ks[1], material.back_material.Ks[2]));
		object_shader_front->setFloat("material.Ns", material.back_material.Ns);
		glBindVertexArray(VAO);
		glPolygonMode(GL_BACK, GL_FILL);
		glCullFace(GL_FRONT);
		glDrawElements(GL_TRIANGLES, 3 * mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
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
	total_thread_num = thread->thread_num;
	obj_aabb_per_thread.resize(total_thread_num);
	current_obj_pos_aabb_per_thread.resize(total_thread_num);
	setMeshStruct(density, ori_mesh);
	this->thread = thread;
	mesh_struct.thread = thread;
	mesh_struct.initialNormalSize();
	mesh_struct.setVertex();
	mesh_struct.setFace();
	mesh_struct.setEdgeForSpring();
	mesh_struct.setEdge();
	mesh_struct.addArounVertex();
	mesh_struct.setThreadIndex(total_thread_num);
	mesh_struct.vertex_for_render = mesh_struct.vertex_position;
	if (!mesh_struct.triangle_indices.empty()) {
		mesh_struct.getRenderNormal();
		mesh_struct.getNormal();
	}




	genBuffer();
	setBuffer();
	setArea();	
	setAnchor();
	mass = mesh_struct.setMass(density);
	mesh_struct.setAnchorPosition();

	//for (unsigned int i = 0; i < mesh_struct.mass.size(); i++) {
	//	std::cout << mesh_struct.mass[i] << std::endl;
	//}
	//for (unsigned int i = 0; i < mesh_struct.edge_vertices.size(); i += 2) {
	//	std::cout << mesh_struct.edge_vertices[i] << " " << mesh_struct.edge_vertices[i + 1] << std::endl;
	//}

	update_stiffness_iteration_number.resize(mesh_struct.vertices.size());
	initialHashAABB();
	setRepresentativePrimitve();
	obtainAABBMoveRadius();
}





void Cloth::setMeshStruct(double density, OriMesh& ori_mesh)
{
	setMaterial(ori_mesh);
	mesh_struct.vertex_position = ori_mesh.vertices;
	if (!ori_mesh.indices.empty()) {
		mesh_struct.triangle_indices.resize(ori_mesh.indices.size() / 3);
		memcpy(mesh_struct.triangle_indices[0].data(), ori_mesh.indices.data(), 12 * mesh_struct.triangle_indices.size());
	}
	this->density = density;

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
	//length_stiffness.resize(mesh_struct.edges.size());
	memcpy(collision_stiffness, single_cloth_info_ref.collision_stiffness, 32);
	memcpy(damp_collision_stiffness, single_cloth_info_ref.collision_stiffness+4, 32);
	//memcpy(collision_stiffness_initial, single_cloth_info_ref.collision_stiffness, 32);
	bend_stiffness = single_cloth_info_ref.bending_stiffness[0];
	damp_bend_stiffness = single_cloth_info_ref.bending_stiffness[1];
	position_stiffness = single_cloth_info_ref.position_stiffness;
	length_stiffness = single_cloth_info_ref.length_stiffness[0];
	damp_length_stiffness = single_cloth_info_ref.length_stiffness[1];
	//std::fill(length_stiffness.begin(), length_stiffness.end(), single_cloth_info_ref.length_stiffness);
	std::array<double, 4> collision_stiff_indicator;
	for (int i = 0; i < 4; ++i) {
		collision_stiff_indicator[i] = collision_stiffness_update_indicator * single_cloth_info_ref.collision_stiffness[i];
	}

}

void Cloth::reset()
{
	memset(rotation_matrix, 0, 72);
	rotation_matrix[0] = 1.0;
	rotation_matrix[4] = 1.0;
	rotation_matrix[8] = 1.0;

	mesh_struct.vertex_position = ori_vertices;
	mesh_struct.vertex_for_render = ori_vertices;
	mesh_struct.getRenderNormal();
	mesh_struct.getNormal();
	for (int i = 0; i < mesh_struct.anchor_vertex.size(); ++i) {
		mesh_struct.anchor_position[i] = mesh_struct.vertex_position[mesh_struct.anchor_vertex[i]];
	}
	obtainAABBMoveRadius();
}


void Cloth::initial()
{
	memset(rotation_matrix, 0, 72);
	rotation_matrix[0] = 1.0;
	rotation_matrix[4] = 1.0;
	rotation_matrix[8] = 1.0;

	mesh_struct.vertex_position = ori_vertices;
	mesh_struct.vertex_for_render = ori_vertices;
	mesh_struct.getRenderNormal();
	mesh_struct.getNormal();
	for (int i = 0; i < mesh_struct.anchor_vertex.size(); ++i) {
		mesh_struct.anchor_position[i] = mesh_struct.vertex_position[mesh_struct.anchor_vertex[i]];
	}
	length_stiffness = single_cloth_info_ref.length_stiffness[0];
	damp_length_stiffness = single_cloth_info_ref.length_stiffness[1];
	//std::fill(length_stiffness.begin(), length_stiffness.end(), single_cloth_info_ref.length_stiffness);
	bend_stiffness = single_cloth_info_ref.bending_stiffness[0];
	damp_bend_stiffness = single_cloth_info_ref.bending_stiffness[1];
	//std::array<double, 4> collision_stiff_indicator;
	//for (int i = 0; i < 4; ++i) {
	//	collision_stiff_indicator[i] = collision_stiffness_update_indicator * single_cloth_info_ref.collision_stiffness[i];
	//}
	memcpy(collision_stiffness, single_cloth_info_ref.collision_stiffness, 32);
	memcpy(damp_collision_stiffness, single_cloth_info_ref.collision_stiffness+4, 32);
	obtainAABBMoveRadius();
}


void Cloth::obtainAABBMoveRadius()
{
	thread->assignTask(this, CURRENT_AABB);
	combineCurrentAABBMoveRadius();
}



void Cloth::initialMouseChosenVertex()
{
	coe_neighbor_vertex_force.clear();
	neighbor_vertex.clear();
}

void Cloth::setAnchor()
{
	//mesh_struct.anchor_vertex.push_back(0);
	//mesh_struct.anchor_vertex.push_back(114);
	//mesh_struct.anchor_vertex.push_back(0);
	//mesh_struct.anchor_vertex.push_back(4);
	mesh_struct.anchor_vertex.push_back(10*10-1);
	mesh_struct.anchor_vertex.push_back(10*9);
	//mesh_struct.anchor_vertex.push_back(50 * 50-1);
	//mesh_struct.anchor_vertex.push_back(49 * 50);
	//mesh_struct.anchor_vertex.push_back(0);

	//for (unsigned int i = 0; i < 20; ++i) {
	//	mesh_struct.anchor_vertex.push_back(mesh_struct.vertex_position.size()-1-i);
	//}


	//mesh_struct.anchor_vertex.push_back(0);
	//mesh_struct.anchor_vertex.push_back(40);
	//mesh_struct.anchor_vertex.push_back(410);
	//mesh_struct.anchor_vertex.push_back(450);

	//mesh_struct.anchor_vertex.push_back(50*50-1);
	//mesh_struct.anchor_vertex.push_back(50 * 49);
	//mesh_struct.anchor_vertex.push_back(0);
	//mesh_struct.anchor_vertex.push_back(49);
	// 
//	mesh_struct.anchor_vertex.push_back(0);
	//mesh_struct.mass_inv[31*31-1] = 0.0;
	//mesh_struct.mass_inv[31*30] = 0.0;
	//mesh_struct.anchor_vertex.push_back(30*30-1);
	//mesh_struct.anchor_vertex.push_back(30 * 29);

	//mesh_struct.anchor_vertex.push_back(15 * 15 - 1);
	//mesh_struct.anchor_vertex.push_back(15 * 14);

	//mesh_struct.anchor_vertex.push_back(100);
	//mesh_struct.anchor_vertex.push_back(0);
	//mesh_struct.anchor_vertex.push_back(15 * 30 - 1);
	//mesh_struct.anchor_vertex.push_back(15 * 29);
	//mesh_struct.anchor_vertex.push_back(101 * 101 - 1);
	//mesh_struct.anchor_vertex.push_back(101 * 100);

	//mesh_struct.anchor_vertex.push_back(101 * 101 - 1);
	//mesh_struct.anchor_vertex.push_back(101 * 100);//131 * 100
}

////VERTEX_AABB
////VERTEX_AABB_WITHOUT_TOLERANCE
//void Cloth::getVertexAABBPerThread(int thread_No, bool has_tolerance)
//{
//	double* aabb = obj_aabb_per_thread[thread_No].data();
//	memset(aabb + 3, 0xFE, 24); //set double to -5.31401e+303
//	memset(aabb, 0x7F, 24); //set double to 1.38242e+306
//
//	std::vector<std::array<double, 3>>* vertex_render=&mesh_struct.vertex_for_render;
//	std::vector<std::array<double, 3>>* vertex=&mesh_struct.vertex_position;
//	unsigned int index_end= mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
//	if (has_tolerance) {
//		for (unsigned int i = mesh_struct.vertex_index_begin_per_thread[thread_No]; i < index_end; ++i) {
//			AABB::obtainAABB(vertex_AABB[i].data(), (*vertex_render)[i].data(), (*vertex)[i].data(), tolerance);	//
//			for (unsigned int j = 0; j < 3; ++j) {
//				if (aabb[j] > vertex_AABB[i][j]) {
//					aabb[j] = vertex_AABB[i][j];
//				}
//			}
//			for (unsigned int j = 3; j < 6; ++j) {
//				if (aabb[j] < vertex_AABB[i][j]) {
//					aabb[j] = vertex_AABB[i][j];
//				}
//			}
//		}
//	}
//	else {
//		for (unsigned int i = mesh_struct.vertex_index_begin_per_thread[thread_No]; i < index_end; ++i) {
//			AABB::obtainAABB(vertex_AABB[i].data(), (*vertex_render)[i].data(), (*vertex)[i].data());	// 
//			for (unsigned int j = 0; j < 3; ++j) {
//				if (aabb[j] > vertex_AABB[i][j]) {
//					aabb[j] = vertex_AABB[i][j];
//				}
//			}
//			for (unsigned int j = 3; j < 6; ++j) {
//				if (aabb[j] < vertex_AABB[i][j]) {
//					aabb[j] = vertex_AABB[i][j];
//				}
//			}
//		}
//	}
//}



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






void Cloth::obtainAABB(bool has_tolerace)
{
	if (has_tolerace) {
		thread->assignTask(this, VERTEX_AABB);
	}
	else {
		thread->assignTask(this, VERTEX_AABB_WITHOUT_TOLERANCE);
	}

	//thread->assignTask(this, EDGE_AABB);
	thread->assignTask(this, EDGE_TRIANGLE_AABB);
	combineObjAABB();
}




void Cloth::setTolerance(double* tolerance_ratio, double ave_edge_length)
{
	for (int i = 0; i < PC_radius.size(); ++i) {
		PC_radius[i] = tolerance_ratio[i] * ave_edge_length;
	}
	//tolerance = 0.5*ave_edge_length;
	tolerance = tolerance_ratio[AABB_SELF_PT] * ave_edge_length;
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

void Cloth::findNeighborVertex(int vertex_index, int recursion_deepth, std::vector<bool>& is_vertex_used)
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


void Cloth::initialNeighborPrimitiveRecording(int cloth_num, int tetrahedron_num, int collider_num, bool use_BVH)
{
	int obj_num = cloth_num + tetrahedron_num;
	triangle_neighbor_obj_triangle.resize(mesh_struct.triangle_indices.size());
	for (int i = 0; i < mesh_struct.triangle_indices.size(); ++i) {
		triangle_neighbor_obj_triangle[i].resize(obj_num);
		for (int j = 0; j < obj_num; ++j) {
			triangle_neighbor_obj_triangle[i][j].reserve(10);
		}
	}

	vertex_neighbor_obj_triangle.resize(mesh_struct.vertex_position.size());
	collide_vertex_obj_triangle.resize(mesh_struct.vertex_position.size());
	for (int i = 0; i < vertex_neighbor_obj_triangle.size(); ++i) {
		vertex_neighbor_obj_triangle[i].resize(obj_num);
		collide_vertex_obj_triangle[i].resize(obj_num);
		for (int j = 0; j < obj_num; ++j) {
			vertex_neighbor_obj_triangle[i][j].reserve(10);
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

	//if (use_BVH) {
	triangle_neighbor_collider_triangle.resize(mesh_struct.triangle_indices.size());
	for (int i = 0; i < mesh_struct.triangle_indices.size(); ++i) {
		triangle_neighbor_collider_triangle[i].resize(collider_num);
		for (int j = 0; j < collider_num; ++j) {
			triangle_neighbor_collider_triangle[i][j].reserve(10);
		}
	}
	//	std::cout<< triangle_neighbor_collider_triangle[0]
	vertex_neighbor_collider_triangle.resize(mesh_struct.vertex_position.size());
	collide_vertex_collider_triangle.resize(mesh_struct.vertex_position.size());
	for (int i = 0; i < vertex_neighbor_collider_triangle.size(); ++i) {
		vertex_neighbor_collider_triangle[i].resize(collider_num);
		collide_vertex_collider_triangle[i].resize(collider_num);
		for (int j = 0; j < collider_num; ++j) {
			vertex_neighbor_collider_triangle[i][j].reserve(10);
			collide_vertex_collider_triangle[i][j].reserve(10);
		}
	}
	//}	
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