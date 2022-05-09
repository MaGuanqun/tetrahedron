#include"triangle_object.h"




void TriangleObject::drawWireframe(Camera* camera, Shader* wireframe_shader)
{
	//setBuffer(cloth_index);
	if (!mesh_struct.triangle_indices.empty()) {
		wireframe_shader->use();
		wireframe_shader->setMat4("projection", camera->GetProjectMatrix());
		wireframe_shader->setMat4("view", camera->GetViewMatrix());
		wireframe_shader->setMat4("model", glm::mat4(1.0));
		wireframe_shader->setVec3("color", wireframe_color);
		wireframe_shader->setFloat("transparent", 1.0f);
		glBindVertexArray(VAO);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_TRIANGLES, 3 * mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}

void TriangleObject::drawWireframeOriPos(Camera* camera, Shader* wireframe_shader)
{
	//setBuffer(cloth_index);
	if (!mesh_struct.triangle_indices.empty()) {
		wireframe_shader->use();
		wireframe_shader->setMat4("projection", camera->GetProjectMatrix());
		wireframe_shader->setMat4("view", camera->GetViewMatrix());
		wireframe_shader->setMat4("model", glm::mat4(1.0));
		wireframe_shader->setVec3("color", wireframe_color);
		wireframe_shader->setFloat("transparent", 0.4f);
		glBindVertexArray(VAO_ori);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_TRIANGLES, 3 * mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}

void TriangleObject::setBufferOriPos()
{
	if (!mesh_struct.triangle_indices.empty())
	{
		glBindVertexArray(VAO_ori);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_ori[0]);
		glBufferData(GL_ARRAY_BUFFER, mesh_struct.vertex_for_render.size() * sizeof(std::array<double, 3>), mesh_struct.vertex_for_render[0].data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_ori);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_struct.triangle_indices.size() * sizeof(std::array<int, 3>), mesh_struct.triangle_indices[0].data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_ori[1]);
		glBufferData(GL_ARRAY_BUFFER, mesh_struct.vertex_normal_for_render.size() * sizeof(std::array<double, 3>), mesh_struct.vertex_normal_for_render[0].data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindVertexArray(0);
	}
}



void TriangleObject::setBuffer()
{
	if (!mesh_struct.triangle_indices.empty())
	{
		glBindVertexArray(VAO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
		glBufferData(GL_ARRAY_BUFFER, mesh_struct.vertex_position.size() * sizeof(std::array<double, 3>), mesh_struct.vertex_position[0].data(), GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_struct.triangle_indices.size() * sizeof(std::array<int, 3>), mesh_struct.triangle_indices[0].data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
		glBufferData(GL_ARRAY_BUFFER, mesh_struct.vertex_normal.size() * sizeof(std::array<double, 3>), mesh_struct.vertex_normal[0].data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindVertexArray(0);
	}

}

void TriangleObject::setRepresentativePrimitve()
{
	representative_vertex_num.resize(mesh_struct.faces.size(), 0);
	representative_edge_num.resize(mesh_struct.faces.size(), 0);
	setRepresentativeVertex(mesh_struct.surface_triangle_index_in_order, mesh_struct.vertices);
	setRepresentativeEdge(mesh_struct.face_edges.data(), mesh_struct.face_edges.size() / 3, mesh_struct.edges.size());
	vertex_from_rep_triangle_index.resize(mesh_struct.vertex_position.size(), -1);
	edge_from_rep_triangle_index.resize(mesh_struct.edges.size(), -1);
	for (int i = 0; i < mesh_struct.faces.size(); ++i) {
		for (int j = 0; j < representative_vertex_num[i]; ++j) {
			vertex_from_rep_triangle_index[mesh_struct.surface_triangle_index_in_order[i][j]] = i;
		}
		for (int j = 0; j < representative_edge_num[i]; ++j) {
			edge_from_rep_triangle_index[mesh_struct.face_edges[3 * i + j]] = i;
		}
	}
}

//CURRENT_AABB
void TriangleObject::getCurrentPosAABB(int thread_No)
{
	double aabb[6];
	memset(aabb + 3, 0xFE, 24); //set double to -5.31401e+303
	memset(aabb, 0x7F, 24); //set double to 1.38242e+306

	std::array<double, 3>* vertex = mesh_struct.vertex_position.data();
	unsigned int end = mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	for (unsigned int i = mesh_struct.vertex_index_begin_per_thread[thread_No]; i < end; ++i) {
		for (unsigned int j = 0; j < 3; ++j) {
			if (aabb[j] > vertex[i].data()[j]) {
				aabb[j] = vertex[i].data()[j];
			}
		}
		for (unsigned int j = 3; j < 6; ++j) {
			if (aabb[j] < vertex[i].data()[j-3]) {
				aabb[j] = vertex[i].data()[j-3];
			}
		}
	}
	memcpy(current_obj_pos_aabb_per_thread[thread_No].data(), aabb, 48);
}


//VERTEX_AABB
//VERTEX_AABB_WITHOUT_TOLERANCE
void TriangleObject::getVertexAABBPerThread(int thread_No, bool has_tolerance)
{
	double* aabb = obj_aabb_per_thread[thread_No].data();
	memset(aabb + 3, 0xFE, 24); //set double to -5.31401e+303
	memset(aabb, 0x7F, 24); //set double to 1.38242e+306

	std::vector<std::array<double, 3>>* vertex_render = &mesh_struct.vertex_for_render;
	std::vector<std::array<double, 3>>* vertex = &mesh_struct.vertex_position;
	unsigned int end = mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	if (has_tolerance) {
		for (unsigned int i = mesh_struct.vertex_index_begin_per_thread[thread_No]; i < end; ++i) {
			AABB::obtainAABB(vertex_AABB[i].data(), (*vertex_render)[i].data(), (*vertex)[i].data(), tolerance);// 
			for (unsigned int j = 0; j < 3; ++j) {
				if (aabb[j] > vertex_AABB[i][j]) {
					aabb[j] = vertex_AABB[i][j];
				}
			}
			for (unsigned int j = 3; j < 6; ++j) {
				if (aabb[j] < vertex_AABB[i][j]) {
					aabb[j] = vertex_AABB[i][j];
				}
			}
		}
	}
	else {
		for (unsigned int i = mesh_struct.vertex_index_begin_per_thread[thread_No]; i < end; ++i) {
			AABB::obtainAABB(vertex_AABB[i].data(), (*vertex_render)[i].data(), (*vertex)[i].data());// 
			for (unsigned int j = 0; j < 3; ++j) {
				if (aabb[j] > vertex_AABB[i][j]) {
					aabb[j] = vertex_AABB[i][j];
				}
			}
			for (unsigned int j = 3; j < 6; ++j) {
				if (aabb[j] < vertex_AABB[i][j]) {
					aabb[j] = vertex_AABB[i][j];
				}
			}
		}
	}
}


void TriangleObject::initialHashAABB()
{
	PC_radius.resize(4);
	hash_index_for_vertex.resize(mesh_struct.vertex_position.size());
	for (int i = 0; i < hash_index_for_vertex.size(); ++i) {
		hash_index_for_vertex[i].reserve(32);
	}
	hash_index_for_edge.resize(mesh_struct.edges.size());
	for (int i = 0; i < hash_index_for_edge.size(); ++i) {
		hash_index_for_edge[i].reserve(32);
	}
	triangle_AABB.resize(mesh_struct.faces.size());
	edge_AABB.resize(mesh_struct.edges.size());
	vertex_AABB.resize(mesh_struct.vertex_position.size());
}


//EDGE_TRIANGLE_AABB
void TriangleObject::getEdgeTriangleAABBPerThread(int thread_No)
{
	unsigned int* edge = mesh_struct.edge_vertices.data();
	unsigned int* vertex_index;
	unsigned int index_end = mesh_struct.edge_index_begin_per_thread[thread_No + 1];
	for (unsigned int i = mesh_struct.edge_index_begin_per_thread[thread_No]; i < index_end; ++i) {
		vertex_index = edge + (i << 1);
		getAABB(edge_AABB[i].data(), vertex_AABB[vertex_index[0]].data(), vertex_AABB[vertex_index[1]].data());
	}
	std::array<int, 3>* face = mesh_struct.triangle_indices.data();
	index_end = mesh_struct.face_index_begin_per_thread[thread_No + 1];
	int* vertex_index_;
	for (unsigned int i = mesh_struct.face_index_begin_per_thread[thread_No]; i < index_end; ++i) {
		vertex_index_ = face[i].data();
		getAABB(triangle_AABB[i].data(), vertex_AABB[vertex_index_[0]].data(), vertex_AABB[vertex_index_[1]].data(), vertex_AABB[vertex_index_[2]].data());
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
		glDrawElements(GL_TRIANGLES, 3 * mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
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
		glDrawElements(GL_TRIANGLES, 3 * mesh_struct.triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}


void TriangleObject::setMaterial(OriMesh& ori_mesh)
{
	material.back_material = ori_mesh.back_material;
	material.front_material = ori_mesh.front_material;
}







//void TriangleObject::genShader()
//{
//	object_shader_front = new Shader("./shader/object_triangle.vs", "./shader/object_triangle.fs");
//	*object_shader_front = Shader("./shader/object_triangle.vs", "./shader/object_triangle.fs");
//}


