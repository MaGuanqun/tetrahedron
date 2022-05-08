#include"draw_collision.h"

void DrawCollision::initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
	Thread* thread)//std::vector<std::vector<std::vector<unsigned int>>>* triangle_patch, std::vector<std::vector<std::vector<unsigned int>>>* patch_vertex
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	this->thread = thread;
	total_obj_num = cloth->size() + tetrahedron->size() + collider->size();
	tetrahedron_end_index = cloth->size() + tetrahedron->size();
	thread_num = thread->thread_num;

	light.ambient = glm::vec3(0.2, 0.2, 0.2);
	light.diffuse = glm::vec3(0.7, 0.7, 0.7);
	light.specular = glm::vec3(0.8, 0.8, 0.8);

	triangle_vertex_index.resize(tetrahedron_end_index);
	collider_triangle_vertex_index.resize(collider->size());
	edge_vertex_index.resize(tetrahedron_end_index);

	vertex_index.resize(tetrahedron_end_index);

	genBuffer();
	reorganzieDataOfObjects();
	initialBoolean();
	draw_vertex.setInObjNum(tetrahedron_end_index);

}

void DrawCollision::initialBoolean()
{
	unsigned int max = tetrahedron_end_index;
	if (max < collider->size()) {
		max = collider->size();
	}
	obj_is_used.resize(max);
	unsigned int count = 0;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		count = cloth->data()[i].mesh_struct.vertex_position.size();
		if (count < cloth->data()[i].mesh_struct.edge_vertices.size() >> 1) {
			count = cloth->data()[i].mesh_struct.edge_vertices.size() >> 1;
		}
		if (count < cloth->data()[i].mesh_struct.triangle_indices.size()) {
			count = cloth->data()[i].mesh_struct.triangle_indices.size();
		}
		if (i < collider->size()) {
			if (count < collider->data()[i].mesh_struct.edge_vertices.size() >> 1) {
				count = collider->data()[i].mesh_struct.edge_vertices.size() >> 1;
			}
			if (count < collider->data()[i].mesh_struct.triangle_indices.size()) {
				count = collider->data()[i].mesh_struct.triangle_indices.size();
			}
			if (count < collider->data()[i].mesh_struct.vertex_position.size()) {
				count = collider->data()[i].mesh_struct.vertex_position.size();
			}
		}
		obj_is_used[i].resize(count,false);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		count = tetrahedron->data()[i].mesh_struct.vertex_position.size();
		if (count < tetrahedron->data()[i].mesh_struct.edge_vertices.size() >> 1) {
			count = tetrahedron->data()[i].mesh_struct.edge_vertices.size() >> 1;
		}
		if (count < tetrahedron->data()[i].mesh_struct.triangle_indices.size()) {
			count = tetrahedron->data()[i].mesh_struct.triangle_indices.size();
		}
		if (i < collider->size()) {
			if (count < collider->data()[i].mesh_struct.edge_vertices.size()>>1) {
				count = collider->data()[i].mesh_struct.edge_vertices.size()>>1;
			}
			if (count < collider->data()[i].mesh_struct.triangle_indices.size()) {
				count = collider->data()[i].mesh_struct.triangle_indices.size();
			}
			if (count < collider->data()[i].mesh_struct.vertex_position.size()) {
				count = collider->data()[i].mesh_struct.vertex_position.size();
			}
		}
		obj_is_used[i+cloth->size()].resize(count, false);
	}
	if (tetrahedron_end_index < collider->size()) {
		for (unsigned int i = tetrahedron_end_index; i < collider->size(); ++i) {
			count = collider->data()[i].mesh_struct.vertex_position.size();
			if (count < collider->data()[i].mesh_struct.edge_vertices.size() >> 1) {
				count = collider->data()[i].mesh_struct.edge_vertices.size() >> 1;
			}
			if (count < collider->data()[i].mesh_struct.triangle_indices.size()) {
				count = collider->data()[i].mesh_struct.triangle_indices.size();
			}
			obj_is_used[i].resize(count, false);
		}	
	}
}

void DrawCollision::reorganzieDataOfObjects()
{
	triangle_indices.resize(tetrahedron_end_index);
	collider_triangle_indices.resize(tetrahedron_end_index);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		triangle_indices[i] = cloth->data()[i].mesh_struct.triangle_indices.data();
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		triangle_indices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.triangle_indices.data();
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		collider_triangle_indices[i ] = collider->data()[i].mesh_struct.triangle_indices.data();
	}

	edge_indices.resize(tetrahedron_end_index);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		edge_indices[i] = cloth->data()[i].mesh_struct.edge_vertices.data();
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		edge_indices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_vertices.data();
	}

	vertex_for_render.resize(total_obj_num);
	vertex_normal_for_render.resize(total_obj_num);
	vertex_number.resize(total_obj_num);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_number[i] = cloth->data()[i].mesh_struct.vertex_for_render.size();
		vertex_for_render[i] = cloth->data()[i].mesh_struct.vertex_for_render.data();
		vertex_normal_for_render[i] = cloth->data()[i].mesh_struct.vertex_normal_for_render.data();
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_number[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.size();
		vertex_for_render[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.data();
		vertex_normal_for_render[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_normal_for_render.data();
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		vertex_number[i + tetrahedron_end_index] = collider->data()[i].mesh_struct.vertex_for_render.size();
		vertex_for_render[i + tetrahedron_end_index] = collider->data()[i].mesh_struct.vertex_for_render.data();
		vertex_normal_for_render[i + tetrahedron_end_index] =collider->data()[i].mesh_struct.vertex_normal_for_render.data();
	}

}

void DrawCollision::setInPairInfo(std::vector<unsigned int>* point_triangle, std::vector<unsigned int>*point_triangle_collider, std::vector<unsigned int>* edge_edge)
{
	point_triangle_target_pos_index = point_triangle;
	point_triangle_collider_target_pos_index = point_triangle_collider;
	edge_edge_target_pos_index = edge_edge;
}



void DrawCollision::resetBooleanVector()
{
	for (unsigned int i = 0; i < obj_is_used.size(); ++i) {
		std::fill(obj_is_used[i].begin(), obj_is_used[i].end(), false);
	}
}

void DrawCollision::setElementIndices()
{
	setIndicesSize();
	resetBooleanVector();
	setTriangleIndices();
	if (!collider->empty()) {
		resetBooleanVector();
		setColliderTriangleIndices();
	}
	resetBooleanVector();
	setEdgeIndices();
	resetBooleanVector();
	setVertexIndices();

	setBuffer();
	draw_vertex.setCollisionVertexData(vertex_for_render, vertex_index);

}



void DrawCollision::setVertexIndices()
{
	unsigned int* pair_index;
	unsigned int size;

	for (unsigned int i = 0; i < thread_num; ++i) {
		pair_index = point_triangle_target_pos_index[i].data() + 1;
		size = point_triangle_target_pos_index[i][0];
		for (unsigned int j = 0; j < size; j += 4) {		
			if (!obj_is_used[*(pair_index + 1)][*pair_index]) {
				vertex_index[*(pair_index + 1)].push_back(*pair_index );
				obj_is_used[*(pair_index + 1)][*pair_index] = true;
			}
			pair_index += 4;
		}
	}

	if (!collider->empty()) {
		for (unsigned int i = 0; i < thread_num; ++i) {
			pair_index = point_triangle_collider_target_pos_index[i].data() + 1;
			size = point_triangle_collider_target_pos_index[i][0];
			for (unsigned int j = 0; j < size; j += 4) {
				if (!obj_is_used[*(pair_index + 1)][*pair_index]) {
					vertex_index[*(pair_index + 1)].push_back(*pair_index);
					obj_is_used[*(pair_index + 1)][*pair_index] = true;
				}
				pair_index += 4;
			}
		}
	}
}


void DrawCollision::setTriangleIndices()
{
	unsigned int* pair_index;
	unsigned int size;

	for (unsigned int i = 0; i < thread_num; ++i) {
		pair_index = point_triangle_target_pos_index[i].data()+1;
		size = point_triangle_target_pos_index[i][0];
		for(unsigned int j=0;j<size;j+=4){		
			if (!obj_is_used[*(pair_index + 3)][*(pair_index + 2)]) {
				triangle_vertex_index[*(pair_index + 3)].resize(triangle_vertex_index[*(pair_index + 3)].size() + 3);
				memcpy(triangle_vertex_index[*(pair_index + 3)].data() + triangle_vertex_index[*(pair_index + 3)].size() - 3,
					triangle_indices[*(pair_index + 3)][*(pair_index + 2)].data(), 12);				
				obj_is_used[*(pair_index + 3)][*(pair_index + 2)] = true;
			}
			pair_index += 4;
		}
	}
}

void DrawCollision::setEdgeIndices()
{
	unsigned int* pair_index;
	unsigned int size;

	for (unsigned int i = 0; i < thread_num; ++i) {
		pair_index = edge_edge_target_pos_index[i].data() + 1;
		size = edge_edge_target_pos_index[i][0];
		for (unsigned int j = 0; j < size; j += 4) {		
			if (!obj_is_used[*(pair_index + 1)][*(pair_index)]) {
				edge_vertex_index[*(pair_index + 1)].emplace_back(edge_indices[*(pair_index + 1)][(*pair_index) << 1]);
				edge_vertex_index[*(pair_index + 1)].emplace_back(edge_indices[*(pair_index + 1)][((*pair_index) << 1) + 1]);
				obj_is_used[*(pair_index + 1)][*(pair_index)] = true;
			}
			if (!obj_is_used[*(pair_index + 3)][*(pair_index + 2)]) {
				edge_vertex_index[*(pair_index + 3)].emplace_back(edge_indices[*(pair_index + 3)][(*(pair_index + 2)) << 1]);
				edge_vertex_index[*(pair_index + 3)].emplace_back(edge_indices[*(pair_index + 3)][((*(pair_index + 2)) << 1) + 1]);
				obj_is_used[*(pair_index + 3)][*(pair_index + 2)] = true;
			}
			pair_index += 4;
		}
	}
}


void DrawCollision::drawVertex(Camera* camera)
{
	draw_vertex.setShaderData(camera);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		draw_vertex.drawCollisionVertex(i, glm::vec3(cloth->data()[i].material.front_material.Kd[2], cloth->data()[i].material.front_material.Kd[1], cloth->data()[i].material.front_material.Kd[0]),
			1.0);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		draw_vertex.drawCollisionVertex(i+cloth->size(), glm::vec3(tetrahedron->data()[i].material.Kd[2], tetrahedron->data()[i].material.Kd[1], tetrahedron->data()[i].material.Kd[0]),
			1.0);
	}

}

void DrawCollision::drawCollision(bool draw_VT, Light& light, float& far_plane, Camera* camera, Shader* object_shader_front)
{
	if (draw_VT) {
		drawVertex(camera);
		drawVT_triangle(light, far_plane, camera, object_shader_front);
	}
	else {

	}
}

void DrawCollision::drawVT_triangle(Light& light, float& far_plane, Camera* camera, Shader* object_shader_front)
{
	object_shader_front->use();
	object_shader_front->setInt("depthMap", 0);
	object_shader_front->setFloat("far_plane", far_plane);
	object_shader_front->setBool("lightIsChosen", true);
	object_shader_front->setVec3("lightPos", camera->position);
	object_shader_front->setVec3("light.ambient", light.ambient);
	object_shader_front->setVec3("light.diffuse", light.diffuse);
	object_shader_front->setVec3("light.specular", light.specular);
	object_shader_front->setVec3("viewPos", camera->position);
	object_shader_front->setBool("lightShadowOn", true);
	object_shader_front->setMat4("projection", camera->GetProjectMatrix());
	object_shader_front->setMat4("view", camera->GetViewMatrix());
	object_shader_front->setMat4("model", glm::mat4(1.0));
	object_shader_front->setFloat("transparence", 1.0);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		object_shader_front->setFloat("material.Ns", cloth->data()[i].material.front_material.Ns);
		object_shader_front->setVec3("material.Kd", glm::vec3(cloth->data()[i].material.front_material.Kd[1], cloth->data()[i].material.front_material.Kd[0], cloth->data()[i].material.front_material.Kd[2]));
		object_shader_front->setVec3("material.Ka", glm::vec3(cloth->data()[i].material.front_material.Ka[1], cloth->data()[i].material.front_material.Ka[0], cloth->data()[i].material.front_material.Ka[2]));
		object_shader_front->setVec3("material.Ks", glm::vec3(cloth->data()[i].material.front_material.Ks[1], cloth->data()[i].material.front_material.Ks[0], cloth->data()[i].material.front_material.Ks[2]));
		glBindVertexArray(VT_VAO[i]);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES,  triangle_vertex_index[i].size(), GL_UNSIGNED_INT, 0);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		object_shader_front->setFloat("material.Ns", tetrahedron->data()[i].material.Ns);
		object_shader_front->setVec3("material.Kd", glm::vec3(tetrahedron->data()[i].material.Kd[1], tetrahedron->data()[i].material.Kd[0], tetrahedron->data()[i].material.Kd[2]));
		object_shader_front->setVec3("material.Ka", glm::vec3(tetrahedron->data()[i].material.Ka[1], tetrahedron->data()[i].material.Ka[0], tetrahedron->data()[i].material.Ka[2]));
		object_shader_front->setVec3("material.Ks", glm::vec3(tetrahedron->data()[i].material.Ks[1], tetrahedron->data()[i].material.Ks[0], tetrahedron->data()[i].material.Ks[2]));
		glBindVertexArray(VT_VAO[i+cloth->size()]);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, triangle_vertex_index[i + cloth->size()].size(), GL_UNSIGNED_INT, 0);
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		object_shader_front->setFloat("material.Ns", collider->data()[i].material.front_material.Ns);
		object_shader_front->setVec3("material.Kd", glm::vec3(collider->data()[i].material.front_material.Kd[1], collider->data()[i].material.front_material.Kd[2], collider->data()[i].material.front_material.Kd[0]));
		object_shader_front->setVec3("material.Ka", glm::vec3(collider->data()[i].material.front_material.Ka[1], collider->data()[i].material.front_material.Ka[2], collider->data()[i].material.front_material.Ka[0]));
		object_shader_front->setVec3("material.Ks", glm::vec3(collider->data()[i].material.front_material.Ks[1], collider->data()[i].material.front_material.Ks[2], collider->data()[i].material.front_material.Ks[0]));
		glBindVertexArray(VT_VAO[i+ tetrahedron_end_index]);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, collider_triangle_vertex_index[i].size(), GL_UNSIGNED_INT, 0);
	}
	glBindVertexArray(0);
}



void DrawCollision::setBuffer()
{
	for (unsigned int i = 0; i < tetrahedron_end_index; ++i) {
		setVertexTriangleBuffer(i);
	}
	if (!collider->empty()) {
		for (unsigned int i = 0; i < collider->size(); ++i) {
			setColliderBuffer(i + tetrahedron_end_index);
		}	
	}
	for (unsigned int i = 0; i < tetrahedron_end_index; ++i) {
		setEdgeEdgeBuffer(i);
	}
}


void DrawCollision::setColliderBuffer(unsigned int obj_index)
{
	glBindVertexArray(VT_VAO[obj_index]);
	glBindBuffer(GL_ARRAY_BUFFER, VT_VBO[3 * obj_index]);
	glBufferData(GL_ARRAY_BUFFER, vertex_number[obj_index]* sizeof(std::array<double, 3>), vertex_for_render[obj_index][0].data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, VT_EBO[obj_index]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, collider_triangle_vertex_index[obj_index- tetrahedron_end_index].size() * sizeof(int), collider_triangle_vertex_index[obj_index- tetrahedron_end_index].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, VT_VBO[3 * obj_index + 1]);
	glBufferData(GL_ARRAY_BUFFER, vertex_number[obj_index] * sizeof(std::array<double, 3>), vertex_normal_for_render[obj_index][0].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);
}

void DrawCollision::setEdgeEdgeBuffer(unsigned int obj_index)
{
	glBindVertexArray(EE_VAO[obj_index]);
	glBindBuffer(GL_ARRAY_BUFFER, EE_VBO[3 * obj_index]);
	glBufferData(GL_ARRAY_BUFFER, vertex_number[obj_index]* sizeof(std::array<double, 3>), vertex_for_render[obj_index][0].data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EE_EBO[obj_index]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, edge_vertex_index[obj_index].size() * sizeof(unsigned int), edge_vertex_index[obj_index].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, EE_VBO[3 * obj_index + 1]);
	glBufferData(GL_ARRAY_BUFFER, vertex_number[obj_index]* sizeof(std::array<double, 3>), vertex_normal_for_render[obj_index][0].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);
}



void DrawCollision::setVertexTriangleBuffer(unsigned int obj_index)
{
	glBindVertexArray(VT_VAO[obj_index]);
	glBindBuffer(GL_ARRAY_BUFFER, VT_VBO[3*obj_index]);
	glBufferData(GL_ARRAY_BUFFER, vertex_number[obj_index]* sizeof(std::array<double, 3>), vertex_for_render[obj_index][0].data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, VT_EBO[obj_index]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,triangle_vertex_index[obj_index].size() * sizeof(int), triangle_vertex_index[obj_index].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, VT_VBO[3*obj_index+1]);
	glBufferData(GL_ARRAY_BUFFER, vertex_number[obj_index] * sizeof(std::array<double, 3>), vertex_normal_for_render[obj_index][0].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);
}


void DrawCollision::setColliderTriangleIndices()
{
	unsigned int* pair_index;
	unsigned int size;

	for (unsigned int i = 0; i < thread_num; ++i) {
		pair_index = point_triangle_collider_target_pos_index[i].data() + 1;
		size = point_triangle_collider_target_pos_index[i][0];
		for (unsigned int j = 0; j < size; j += 4) {		
			if (!obj_is_used[*(pair_index + 3)][*(pair_index + 2)]) {
				collider_triangle_vertex_index[*(pair_index + 3)].resize(collider_triangle_vertex_index[*(pair_index + 3)].size() + 3);
				memcpy(collider_triangle_vertex_index[*(pair_index + 3)].data() + collider_triangle_vertex_index[*(pair_index + 3)].size() - 3,
					collider_triangle_indices[*(pair_index + 3)][*(pair_index + 2)].data(), 12);
				obj_is_used[*(pair_index + 3)][*(pair_index + 2)] = true;
			}			
			pair_index += 4;
		}
	}
}

void DrawCollision::setIndicesSize()
{
	for (unsigned int i = 0; i < tetrahedron_end_index; ++i) {
		triangle_vertex_index[i].clear();
		edge_vertex_index[i].clear();
		vertex_index[i].clear();
	}

	for (unsigned int i = 0; i < cloth->size(); ++i) {
		triangle_vertex_index[i].reserve(cloth->data()[i].mesh_struct.triangle_indices.size()/2);
		edge_vertex_index[i].reserve(cloth->data()[i].mesh_struct.edge_vertices.size()/6);
		vertex_index[i].reserve(cloth->data()[i].mesh_struct.vertex_position.size()/6);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		triangle_vertex_index[i+cloth->size()].reserve(tetrahedron->data()[i].mesh_struct.triangle_indices.size() / 2);
		edge_vertex_index[i+cloth->size()].reserve(tetrahedron->data()[i].mesh_struct.edge_vertices.size() / 2);
		vertex_index[i+cloth->size()].reserve(tetrahedron->data()[i].mesh_struct.vertex_position.size() / 6);
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		collider_triangle_vertex_index[i].clear();
		collider_triangle_vertex_index[i].reserve(collider->data()[i].mesh_struct.triangle_indices.size()/2);
	}


}



void DrawCollision::genBuffer()
{
	VT_VAO.resize(total_obj_num);
	VT_VBO.resize(3*total_obj_num);
	VT_EBO.resize(total_obj_num);

	EE_VAO.resize(total_obj_num);
	EE_VBO.resize(3 * total_obj_num);
	EE_EBO.resize(total_obj_num);

	glGenVertexArrays(total_obj_num, VT_VAO.data());
	glGenBuffers(3* total_obj_num, VT_VBO.data());
	glGenBuffers(total_obj_num, VT_EBO.data());
	glGenVertexArrays(total_obj_num,EE_VAO.data());
	glGenBuffers(3 * total_obj_num, EE_VBO.data());
	glGenBuffers(total_obj_num, EE_EBO.data());


}

