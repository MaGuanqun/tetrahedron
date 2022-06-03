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

	triangle_vertex_index_in_a_cell.resize(tetrahedron_end_index);
	collider_triangle_vertex_index_in_a_cell.resize(collider->size());
	edge_vertex_index_in_a_cell.resize(tetrahedron_end_index);
	vertex_index_in_a_cell.resize(tetrahedron_end_index);

	triangle_vertex_index_all_cell.resize(tetrahedron_end_index);
	collider_triangle_vertex_index_all_cell.resize(collider->size());
	edge_vertex_index_all_cell.resize(tetrahedron_end_index);
	vertex_index_all_cell.resize(tetrahedron_end_index);



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
	obj_is_used2.resize(max);
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
		obj_is_used2[i].resize(count,false);
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
		obj_is_used2[i+cloth->size()].resize(count, false);
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
			obj_is_used2[i].resize(count, false);
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
		vertex_number[i] = cloth->data()[i].mesh_struct.vertex_position.size();
		vertex_for_render[i] = cloth->data()[i].mesh_struct.vertex_position.data();
		vertex_normal_for_render[i] = cloth->data()[i].mesh_struct.vertex_normal_for_render.data();
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_number[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.size();
		vertex_for_render[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		vertex_normal_for_render[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_normal_for_render.data();
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		vertex_number[i + tetrahedron_end_index] = collider->data()[i].mesh_struct.vertex_position.size();
		vertex_for_render[i + tetrahedron_end_index] = collider->data()[i].mesh_struct.vertex_position.data();
		vertex_normal_for_render[i + tetrahedron_end_index] =collider->data()[i].mesh_struct.vertex_normal_for_render.data();
	}

}

void DrawCollision::setInPairInfo(std::vector < std::vector<unsigned int>>* point_triangle, std::vector < std::vector<unsigned int>>*point_triangle_collider, std::vector < std::vector<unsigned int>>* edge_edge)
{
	//std::cout << "initial draw c " << point_triangle << std::endl;
	point_triangle_target_pos_index = point_triangle;
	point_triangle_collider_target_pos_index = point_triangle_collider;
	edge_edge_target_pos_index = edge_edge;
}


void DrawCollision::resetBooleanVector(std::vector<std::vector<bool>>& obj_is_used)
{
	for (unsigned int i = 0; i < obj_is_used.size(); ++i) {
		std::fill(obj_is_used[i].begin(), obj_is_used[i].end(), false);
	}
}

void DrawCollision::setElementIndices()
{
	setIndicesSize(triangle_vertex_index.data(),edge_vertex_index.data(),vertex_index.data(),collider_triangle_vertex_index.data());
	resetBooleanVector(obj_is_used);


	//std::cout << "compare " << point_triangle_target_pos_index->data() << std::endl;

	setTriangleIndices();
	if (!collider->empty()) {
		resetBooleanVector(obj_is_used);
		setColliderTriangleIndices();
	}
	resetBooleanVector(obj_is_used);
	setEdgeIndices();
	resetBooleanVector(obj_is_used);
	setVertexIndices();

	setBuffer();

	//std::cout << vertex_index[0].size() << std::endl;

	draw_vertex.setCollisionVertexData(vertex_for_render, vertex_index);

}



void DrawCollision::setElementInAllCell(std::vector<std::vector<unsigned int>>& vertex_index, std::vector<std::vector<unsigned int>>& triangle_index,
	std::vector<std::vector<unsigned int>>& edge_index)
{
	setIndicesSize(triangle_vertex_index_all_cell.data(), edge_vertex_index_all_cell.data(), vertex_index_all_cell.data(), collider_triangle_vertex_index_all_cell.data());
	resetBooleanVector(obj_is_used);
	setTriangleIndicesInAllCell(triangle_index);
	resetBooleanVector(obj_is_used);
	setVertexIndicesInAllCell(vertex_index);
	resetBooleanVector(obj_is_used);
	setEdgeIndicesInAllCell(edge_index);

	setBufferAllCell();
	draw_vertex.setCollisionVertexDataAllCell(vertex_for_render, vertex_index_all_cell);
}

void DrawCollision::setElementInOneCell(std::vector<std::vector<unsigned int>>& triangle_index,
	std::vector<std::vector<unsigned int>>& edge_index)
{
	setIndicesSize(triangle_vertex_index_in_a_cell.data(), edge_vertex_index_in_a_cell.data(), vertex_index_in_a_cell.data(), collider_triangle_vertex_index_in_a_cell.data());
	resetBooleanVector(obj_is_used);
	resetBooleanVector(obj_is_used2);

	setBooleanTrue(obj_is_used.data(), obj_is_used2.data(), true);
	setTriangleIndicesInOneCell(triangle_index);
	setVertexIndicesInOneCell(vertex_index_all_cell);

	resetBooleanVector(obj_is_used);
	setBooleanTrue(obj_is_used.data(), obj_is_used2.data(), false);
	setEdgeIndicesInOneCell(edge_index);

	setBufferInACell();
	draw_vertex.setCollisionVertexDataInACell(vertex_for_render, vertex_index_in_a_cell);
}



void DrawCollision::setVertexIndicesInAllCell(std::vector<std::vector<unsigned int>>& vertex_index_)
{
	for (unsigned int i = 0; i < vertex_index_.size(); ++i) {
		for (unsigned int j = 0; j < vertex_index_[i].size(); ++j) {
			if (!obj_is_used[i][vertex_index_[i][j]]) {
				vertex_index_all_cell[i].emplace_back(vertex_index_[i][j]);
				obj_is_used[i][vertex_index_[i][j]] = true;
			}
		}
	}
}

void DrawCollision::setVertexIndicesInOneCell(std::vector<std::vector<unsigned int>>& vertex_index_)
{
	for (unsigned int i = 0; i < vertex_index_.size(); ++i) {
		for (unsigned int j = 0; j < vertex_index_[i].size(); ++j) {
			if (obj_is_used[i][vertex_index_[i][j]]) {
				vertex_index_in_a_cell[i].emplace_back(vertex_index_[i][j]);
			}
		}
	}
}

void DrawCollision::setVertexIndices()
{
	unsigned int* pair_index;
	unsigned int size;

	for (unsigned int i = 0; i < thread_num; ++i) {
		pair_index = point_triangle_target_pos_index->data()[i].data() + 1;
		size = point_triangle_target_pos_index->data()[i][0];
		for (unsigned int j = 0; j < size; j += 4) {		
			if (!obj_is_used[*(pair_index + 1)][*pair_index]) {
				vertex_index[*(pair_index + 1)].emplace_back(*pair_index );
				obj_is_used[*(pair_index + 1)][*pair_index] = true;
			}
			pair_index += 4;
		}
	}

	if (!collider->empty()) {
		for (unsigned int i = 0; i < thread_num; ++i) {
			pair_index = point_triangle_collider_target_pos_index->data()[i].data() + 1;
			size = point_triangle_collider_target_pos_index->data()[i][0];
			for (unsigned int j = 0; j < size; j += 4) {
				if (!obj_is_used[*(pair_index + 1)][*pair_index]) {
					vertex_index[*(pair_index + 1)].emplace_back(*pair_index);
					obj_is_used[*(pair_index + 1)][*pair_index] = true;
				}
				pair_index += 4;
			}
		}
	}
}


void DrawCollision::setBooleanTrue(std::vector<bool>* obj_is_used, std::vector<bool>* obj_is_used2, bool for_vt)
{
	unsigned int* pair_index;
	unsigned int size;
	if (for_vt) {
		for (unsigned int i = 0; i < thread_num; ++i) {
			pair_index = point_triangle_target_pos_index->data()[i].data() + 1;
			size = point_triangle_target_pos_index->data()[i][0];
			for (unsigned int j = 0; j < size; j += 4) {
				obj_is_used[*(pair_index + 1)][*pair_index] = true;
				obj_is_used2[*(pair_index + 3)][*(pair_index + 2)] = true;
				pair_index += 4;
			}
		}
	}
	else {
		for (unsigned int i = 0; i < thread_num; ++i) {
			pair_index = edge_edge_target_pos_index->data()[i].data() + 1;
			size = edge_edge_target_pos_index->data()[i][0];
			for (unsigned int j = 0; j < size; j += 4) {
				obj_is_used[*(pair_index + 1)][*pair_index] = true;
				obj_is_used[*(pair_index + 3)][*(pair_index + 2)] = true;
				pair_index += 4;
			}
		}
	}
}


void DrawCollision::setTriangleIndicesInAllCell(std::vector<std::vector<unsigned int>>& triangle_index)
{
	for (unsigned int i = 0; i < triangle_index.size(); ++i) {
		for (unsigned int j = 0; j < triangle_index[i].size(); ++j) {
			if (!obj_is_used[i][triangle_index[i][j]]) {
				triangle_vertex_index_all_cell[i].resize(triangle_vertex_index_all_cell[i].size() + 3);
				memcpy(triangle_vertex_index_all_cell[i].data() + triangle_vertex_index_all_cell[i].size() - 3,
					triangle_indices[i][triangle_index[i][j]].data(), 12);
				obj_is_used[i][triangle_index[i][j]] = true;
			}
		}
	}
}





void DrawCollision::setTriangleIndicesInOneCell(std::vector<std::vector<unsigned int>>& triangle_index)
{
	for (unsigned int i = 0; i < triangle_index.size(); ++i) {
		for (unsigned int j = 0; j < triangle_index[i].size(); ++j) {
			if (obj_is_used2[i][triangle_index[i][j]]) {
				triangle_vertex_index_in_a_cell[i].resize(triangle_vertex_index_in_a_cell[i].size() + 3);
				memcpy(triangle_vertex_index_in_a_cell[i].data() + triangle_vertex_index_in_a_cell[i].size() - 3,
					triangle_indices[i][triangle_index[i][j]].data(), 12);
			}
		}
	}
}


void DrawCollision::setTriangleIndices()
{
	unsigned int* pair_index;
	unsigned int size;

	for (unsigned int i = 0; i < thread_num; ++i) {
		pair_index = point_triangle_target_pos_index->data()[i].data()+1;
		size = point_triangle_target_pos_index->data()[i][0];
		for(unsigned int j=0;j<size;j+=4){		
			if (!obj_is_used[*(pair_index + 3)][*(pair_index + 2)]) {
				triangle_vertex_index[*(pair_index + 3)].resize(triangle_vertex_index[*(pair_index + 3)].size() + 3);
				memcpy(triangle_vertex_index[*(pair_index + 3)].data() + triangle_vertex_index[*(pair_index + 3)].size() - 3,
					triangle_indices[*(pair_index + 3)][*(pair_index + 2)].data(), 12);				
				obj_is_used[*(pair_index + 3)][*(pair_index + 2)] = true;
			}
			pair_index += 4;
		}
		//std::cout<<"draw " << size << std::endl;
	}
}


void DrawCollision::setEdgeIndicesInAllCell(std::vector<std::vector<unsigned int>>& edge_index)
{
	for (unsigned int i = 0; i < edge_index.size(); ++i) {
		for (unsigned int j = 0; j < edge_index[i].size(); ++j) {
			if (!obj_is_used[i][edge_index[i][j]]) {
				edge_vertex_index_all_cell[i].emplace_back(edge_indices[i][edge_index[i][j]<<1]);
				edge_vertex_index_all_cell[i].emplace_back(edge_indices[i][(edge_index[i][j] << 1) + 1]);
				obj_is_used[i][edge_index[i][j]] = true;
			}
		}
	}
}


void DrawCollision::setEdgeIndicesInOneCell(std::vector<std::vector<unsigned int>>& edge_index)
{
	for (unsigned int i = 0; i < edge_index.size(); ++i) {
		for (unsigned int j = 0; j < edge_index[i].size(); ++j) {
			if (obj_is_used[i][edge_index[i][j]]) {
				edge_vertex_index_in_a_cell[i].emplace_back(edge_indices[i][edge_index[i][j] << 1]);
				edge_vertex_index_in_a_cell[i].emplace_back(edge_indices[i][(edge_index[i][j] << 1) + 1]);
			}
		}
	}
}



void DrawCollision::setEdgeIndices()
{
	unsigned int* pair_index;
	unsigned int size;

	for (unsigned int i = 0; i < thread_num; ++i) {
		pair_index = edge_edge_target_pos_index->data()[i].data() + 1;
		size = edge_edge_target_pos_index->data()[i][0];
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


void DrawCollision::drawVertex(Camera* camera, std::vector<std::vector<bool>>& show_collision_element)
{
	draw_vertex.setShaderData(camera);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		if (show_collision_element[6+CLOTH_][i]) {
			draw_vertex.drawAllCollisionVertex(i, glm::vec3(cloth->data()[i].material.front_material.Kd[2], cloth->data()[i].material.front_material.Kd[1], cloth->data()[i].material.front_material.Kd[0]),
				1.0);
		}
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		if (show_collision_element[6+TETRAHEDRON_][i]) {
			draw_vertex.drawAllCollisionVertex(i + cloth->size(), glm::vec3(tetrahedron->data()[i].material.Kd[2], tetrahedron->data()[i].material.Kd[1], tetrahedron->data()[i].material.Kd[0]),
				1.0);
		}
	}
}

void DrawCollision::drawVertexCell(Camera* camera, std::vector<std::vector<bool>>& show_collision_element, bool show_all_element)
{
	draw_vertex.setShaderData(camera);
	if (show_all_element) {
		for (unsigned int i = 0; i < cloth->size(); ++i) {
			if (show_collision_element[6 + CLOTH_][i]) {
				draw_vertex.drawAllCellVertex(i, glm::vec3(cloth->data()[i].material.front_material.Kd[2], cloth->data()[i].material.front_material.Kd[1], cloth->data()[i].material.front_material.Kd[0]),
					1.0);
			}
		}
		for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
			if (show_collision_element[6 + TETRAHEDRON_][i]) {
				draw_vertex.drawAllCellVertex(i + cloth->size(), glm::vec3(tetrahedron->data()[i].material.Kd[2], tetrahedron->data()[i].material.Kd[1], tetrahedron->data()[i].material.Kd[0]),
					1.0);
			}
		}
	}
	else {
		for (unsigned int i = 0; i < cloth->size(); ++i) {
			if (show_collision_element[6 + CLOTH_][i]) {
				draw_vertex.drawOneCellVertex(i, glm::vec3(cloth->data()[i].material.front_material.Kd[2], cloth->data()[i].material.front_material.Kd[1], cloth->data()[i].material.front_material.Kd[0]),
					1.0);
			}
		}
		for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
			if (show_collision_element[6 + TETRAHEDRON_][i]) {
				draw_vertex.drawOneCellVertex(i + cloth->size(), glm::vec3(tetrahedron->data()[i].material.Kd[2], tetrahedron->data()[i].material.Kd[1], tetrahedron->data()[i].material.Kd[0]),
					1.0);
			}
		}
	}
}

void DrawCollision::drawCollision(bool draw_VT, Light& light, Camera* camera, Shader* object_shader_front, Shader* wireframe_shader, std::vector<std::vector<bool>>& show_collision_element)
{
	if (draw_VT) {
		drawVertex(camera, show_collision_element);
		drawVT_triangle(light, camera, object_shader_front, show_collision_element, VT_VAO.data(), triangle_vertex_index.data(),
			collider_triangle_vertex_index.data(), wireframe_shader);
	}
	else {
		drawEdge(light, camera, object_shader_front, show_collision_element,EE_VAO.data(),edge_vertex_index.data());
	}
}

void DrawCollision::drawCollisionCell(bool draw_VT, Light& light, Camera* camera, Shader* object_shader_front, Shader* wireframe_shader,
	std::vector<std::vector<bool>>& show_collision_element,
	bool show_all_element)
{
	if (draw_VT) {
		drawVertexCell(camera, show_collision_element, show_all_element);
		if (show_all_element) {
			drawVT_triangle(light, camera, object_shader_front, show_collision_element, VT_VAO_all_cell.data(), triangle_vertex_index_all_cell.data(), collider_triangle_vertex_index_all_cell.data(),
				wireframe_shader);
		}
		else {
			drawVT_triangle(light, camera, object_shader_front, show_collision_element, VT_VAO_collide_in_a_cell.data(), triangle_vertex_index_in_a_cell.data(), collider_triangle_vertex_index_in_a_cell.data(),
				wireframe_shader);
		}
	}
	else {
		if (show_all_element) {
			drawEdge(light, camera, object_shader_front, show_collision_element, EE_VAO_all_cell.data(), edge_vertex_index_all_cell.data());
		}
		else {
			drawEdge(light, camera, object_shader_front, show_collision_element, EE_VAO_collide_in_a_cell.data(), edge_vertex_index_in_a_cell.data());
		}
	}
}



void DrawCollision::drawEdge(Light& light, Camera* camera, Shader* object_shader_front, std::vector<std::vector<bool>>& show_collision_element,
	unsigned int* EE_VAO, std::vector<int>* edge_vertex_index)
{

	setInShader(light, camera, object_shader_front,1.0);
	object_shader_front->use();
	glLineWidth(4.0);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		if (show_collision_element[6+CLOTH_][i]) {
			object_shader_front->setFloat("material.Ns", cloth->data()[i].material.front_material.Ns);
			object_shader_front->setVec3("material.Kd", glm::vec3(4.0 * cloth->data()[i].material.front_material.Kd[0], cloth->data()[i].material.front_material.Kd[2], cloth->data()[i].material.front_material.Kd[1]));
			object_shader_front->setVec3("material.Ka", glm::vec3(4.0 * cloth->data()[i].material.front_material.Ka[0], 0.3*cloth->data()[i].material.front_material.Ka[2], 0.3*cloth->data()[i].material.front_material.Ka[1]));
			object_shader_front->setVec3("material.Ks", glm::vec3(4.0 * cloth->data()[i].material.front_material.Ks[0], cloth->data()[i].material.front_material.Ks[2], cloth->data()[i].material.front_material.Ks[1]));
			glBindVertexArray(EE_VAO[i]);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glDrawElements(GL_LINES, edge_vertex_index[i].size(), GL_UNSIGNED_INT, 0);
		}
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		if (show_collision_element[6+TETRAHEDRON_][i]) {
			object_shader_front->setFloat("material.Ns", tetrahedron->data()[i].material.Ns);
			object_shader_front->setVec3("material.Kd", glm::vec3(4.0 * tetrahedron->data()[i].material.Kd[0], tetrahedron->data()[i].material.Kd[2], tetrahedron->data()[i].material.Kd[1]));
			object_shader_front->setVec3("material.Ka", glm::vec3(4.0 * tetrahedron->data()[i].material.Ka[0],0.3* tetrahedron->data()[i].material.Ka[2],0.3* tetrahedron->data()[i].material.Ka[1]));
			object_shader_front->setVec3("material.Ks", glm::vec3(4.0 * tetrahedron->data()[i].material.Ks[0], tetrahedron->data()[i].material.Ks[2], tetrahedron->data()[i].material.Ks[1]));
			glBindVertexArray(EE_VAO[i + cloth->size()]);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glDrawElements(GL_LINES, edge_vertex_index[i + cloth->size()].size(), GL_UNSIGNED_INT, 0);
		}
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		if (show_collision_element[6+COLLIDER_][i]) {
			object_shader_front->setFloat("material.Ns", collider->data()[i].material.front_material.Ns);
			object_shader_front->setVec3("material.Kd", glm::vec3(4.0 * collider->data()[i].material.front_material.Kd[1], collider->data()[i].material.front_material.Kd[2], collider->data()[i].material.front_material.Kd[0]));
			object_shader_front->setVec3("material.Ka", glm::vec3(4.0 * collider->data()[i].material.front_material.Ka[1], collider->data()[i].material.front_material.Ka[2], collider->data()[i].material.front_material.Ka[0]));
			object_shader_front->setVec3("material.Ks", glm::vec3(4.0 * collider->data()[i].material.front_material.Ks[1], collider->data()[i].material.front_material.Ks[2], collider->data()[i].material.front_material.Ks[0]));
			glBindVertexArray(EE_VAO[i + tetrahedron_end_index]);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glDrawElements(GL_LINES, edge_vertex_index[i].size(), GL_UNSIGNED_INT, 0);
		}
	}
	glBindVertexArray(0);
	glLineWidth(1.0);
}



void DrawCollision::setInShader(Light& light, Camera* camera, Shader* object_shader_front, double transparent)
{
	object_shader_front->use();
	object_shader_front->setInt("depthMap", 0);
	object_shader_front->setFloat("far_plane", camera->far_plane);
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
	object_shader_front->setFloat("transparence", transparent);
}


void DrawCollision::drawVT_triangleWireframe(Camera* camera, Shader* wireframe_shader, std::vector<std::vector<bool>>& show_collision_element,
	unsigned int* VT_VAO, std::vector<int>* triangle_vertex_index, std::vector<int>* collider_triangle_vertex_index)
{

	wireframe_shader->use();
	wireframe_shader->setMat4("projection", camera->GetProjectMatrix());
	wireframe_shader->setMat4("view", camera->GetViewMatrix());
	wireframe_shader->setMat4("model", glm::mat4(1.0));
	wireframe_shader->setFloat("transparent", 1.0f);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		if (show_collision_element[6 + CLOTH_][i]) {
			wireframe_shader->setVec3("color", cloth->data()[i].wireframe_color);
			glBindVertexArray(VT_VAO[i]);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glDrawElements(GL_TRIANGLES, triangle_vertex_index[i].size(), GL_UNSIGNED_INT, 0);
		}
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		if (show_collision_element[6 + TETRAHEDRON_][i]) {
			wireframe_shader->setVec3("color", tetrahedron->data()[i].wireframe_color);
			glBindVertexArray(VT_VAO[i + cloth->size()]);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glDrawElements(GL_TRIANGLES, triangle_vertex_index[i].size(), GL_UNSIGNED_INT, 0);
		}
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		if (show_collision_element[6 + COLLIDER_][i]) {
			wireframe_shader->setVec3("color", collider->data()[i].wireframe_color);
			glBindVertexArray(VT_VAO[i + tetrahedron_end_index]);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glDrawElements(GL_TRIANGLES, triangle_vertex_index[i].size(), GL_UNSIGNED_INT, 0);
		}
	}
	glBindVertexArray(0);
}


void DrawCollision::drawVT_triangle(Light& light,  Camera* camera, Shader* object_shader_front, std::vector<std::vector<bool>>& show_collision_element,
	unsigned int* VT_VAO, std::vector<int>* triangle_vertex_index, std::vector<int>* collider_triangle_vertex_index, Shader* wireframe_shader)
{
	setInShader(light, camera, object_shader_front,1.0);
	object_shader_front->use();
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		if (show_collision_element[6+CLOTH_][i]) {
			object_shader_front->setFloat("material.Ns", cloth->data()[i].material.front_material.Ns);
			object_shader_front->setVec3("material.Kd", glm::vec3(0.5*cloth->data()[i].material.front_material.Kd[0], cloth->data()[i].material.front_material.Kd[2], cloth->data()[i].material.front_material.Kd[1]));
			object_shader_front->setVec3("material.Ka", glm::vec3(0.5 * cloth->data()[i].material.front_material.Ka[0], cloth->data()[i].material.front_material.Ka[2], cloth->data()[i].material.front_material.Ka[1]));
			object_shader_front->setVec3("material.Ks", glm::vec3(0.5 * cloth->data()[i].material.front_material.Ks[0], cloth->data()[i].material.front_material.Ks[2], cloth->data()[i].material.front_material.Ks[1]));
			glBindVertexArray(VT_VAO[i]);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glDrawElements(GL_TRIANGLES, triangle_vertex_index[i].size(), GL_UNSIGNED_INT, 0);
		}
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		if (show_collision_element[6 + TETRAHEDRON_][i]) {
			object_shader_front->setFloat("material.Ns", tetrahedron->data()[i].material.Ns);
			object_shader_front->setVec3("material.Kd", glm::vec3(0.5 * tetrahedron->data()[i].material.Kd[0], tetrahedron->data()[i].material.Kd[2], tetrahedron->data()[i].material.Kd[1]));
			object_shader_front->setVec3("material.Ka", glm::vec3(0.5 * tetrahedron->data()[i].material.Ka[0], tetrahedron->data()[i].material.Ka[2], tetrahedron->data()[i].material.Ka[1]));
			object_shader_front->setVec3("material.Ks", glm::vec3(0.5 * tetrahedron->data()[i].material.Ks[0], tetrahedron->data()[i].material.Ks[2], tetrahedron->data()[i].material.Ks[1]));
			glBindVertexArray(VT_VAO[i + cloth->size()]);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glDrawElements(GL_TRIANGLES, triangle_vertex_index[i + cloth->size()].size(), GL_UNSIGNED_INT, 0);
		}
		
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		if (show_collision_element[6 + COLLIDER_][i]) {
			object_shader_front->setFloat("material.Ns", collider->data()[i].material.front_material.Ns);
			object_shader_front->setVec3("material.Kd", glm::vec3(4.0 * collider->data()[i].material.front_material.Kd[1], collider->data()[i].material.front_material.Kd[2], collider->data()[i].material.front_material.Kd[0]));
			object_shader_front->setVec3("material.Ka", glm::vec3(4.0 * collider->data()[i].material.front_material.Ka[1], collider->data()[i].material.front_material.Ka[2], collider->data()[i].material.front_material.Ka[0]));
			object_shader_front->setVec3("material.Ks", glm::vec3(4.0 * collider->data()[i].material.front_material.Ks[1], collider->data()[i].material.front_material.Ks[2], collider->data()[i].material.front_material.Ks[0]));
			glBindVertexArray(VT_VAO[i + tetrahedron_end_index]);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glDrawElements(GL_TRIANGLES, collider_triangle_vertex_index[i].size(), GL_UNSIGNED_INT, 0);
		}
	}
	glBindVertexArray(0);

	drawVT_triangleWireframe(camera, wireframe_shader, show_collision_element, VT_VAO, triangle_vertex_index, collider_triangle_vertex_index);

}



void DrawCollision::setBuffer()
{
	for (unsigned int i = 0; i < tetrahedron_end_index; ++i) {
		setBuffer(i,VT_VAO.data(),VT_VBO.data(),VT_EBO.data(),triangle_vertex_index.data());
	}
	if (!collider->empty()) {
		for (unsigned int i = 0; i < collider->size(); ++i) {
			setBuffer(i + tetrahedron_end_index, VT_VAO.data(), VT_VBO.data(), VT_EBO.data(), collider_triangle_vertex_index.data());
		}	
	}
	for (unsigned int i = 0; i < tetrahedron_end_index; ++i) {
		setBuffer(i,EE_VAO.data(),EE_VBO.data(),EE_EBO.data(), edge_vertex_index.data());
	}
}


void DrawCollision::setBufferAllCell()
{
	for (unsigned int i = 0; i < tetrahedron_end_index; ++i) {
		setBuffer(i, VT_VAO_all_cell.data(), VT_VBO_all_cell.data(), VT_EBO_all_cell.data(), triangle_vertex_index_all_cell.data());
	}
	if (!collider->empty()) {
		for (unsigned int i = 0; i < collider->size(); ++i) {
			setBuffer(i + tetrahedron_end_index, VT_VAO_all_cell.data(), VT_VBO_all_cell.data(), VT_EBO_all_cell.data(), collider_triangle_vertex_index_all_cell.data());
		}
	}
	for (unsigned int i = 0; i < tetrahedron_end_index; ++i) {
		setBuffer(i, EE_VAO_all_cell.data(), EE_VBO_all_cell.data(), EE_EBO_all_cell.data(), edge_vertex_index_all_cell.data());
	}
}


void DrawCollision::setBufferInACell()
{
	for (unsigned int i = 0; i < tetrahedron_end_index; ++i) {
		setBuffer(i, VT_VAO_collide_in_a_cell.data(), VT_VBO_collide_in_a_cell.data(), VT_EBO_collide_in_a_cell.data(), triangle_vertex_index_in_a_cell.data());
	}
	if (!collider->empty()) {
		for (unsigned int i = 0; i < collider->size(); ++i) {
			setBuffer(i + tetrahedron_end_index, VT_VAO_collide_in_a_cell.data(), VT_VBO_collide_in_a_cell.data(), VT_EBO_collide_in_a_cell.data(), collider_triangle_vertex_index_in_a_cell.data());
		}
	}
	for (unsigned int i = 0; i < tetrahedron_end_index; ++i) {
		setBuffer(i, EE_VAO_collide_in_a_cell.data(), EE_VBO_collide_in_a_cell.data(), EE_EBO_collide_in_a_cell.data(), edge_vertex_index_in_a_cell.data());
	}
}


void DrawCollision::setBuffer(unsigned int obj_index, unsigned int* VT_VAO, unsigned int* VT_VBO, unsigned int* VT_EBO,
	std::vector<int>* triangle_vertex_index)
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
		pair_index = point_triangle_collider_target_pos_index->data()[i].data() + 1;
		size = point_triangle_collider_target_pos_index->data()[i][0];
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

void DrawCollision::setIndicesSize(std::vector<int>* triangle_vertex_index, std::vector<int>* edge_vertex_index, std::vector<unsigned int>* vertex_index,
	std::vector<int>* collider_triangle_vertex_index)
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



	VT_VAO_collide_in_a_cell.resize(total_obj_num);
	VT_VBO_collide_in_a_cell.resize(3 * total_obj_num);
	VT_EBO_collide_in_a_cell.resize(total_obj_num);

	EE_VAO_collide_in_a_cell.resize(total_obj_num);
	EE_VBO_collide_in_a_cell.resize(3 * total_obj_num);
	EE_EBO_collide_in_a_cell.resize(total_obj_num);

	glGenVertexArrays(total_obj_num, VT_VAO_collide_in_a_cell.data());
	glGenBuffers(3 * total_obj_num, VT_VBO_collide_in_a_cell.data());
	glGenBuffers(total_obj_num, VT_EBO_collide_in_a_cell.data());
	glGenVertexArrays(total_obj_num, EE_VAO_collide_in_a_cell.data());
	glGenBuffers(3 * total_obj_num, EE_VBO_collide_in_a_cell.data());
	glGenBuffers(total_obj_num, EE_EBO_collide_in_a_cell.data());

	VT_VAO_all_cell.resize(total_obj_num);
	VT_VBO_all_cell.resize(3 * total_obj_num);
	VT_EBO_all_cell.resize(total_obj_num);

	EE_VAO_all_cell.resize(total_obj_num);
	EE_VBO_all_cell.resize(3 * total_obj_num);
	EE_EBO_all_cell.resize(total_obj_num);

	glGenVertexArrays(total_obj_num, VT_VAO_all_cell.data());
	glGenBuffers(3 * total_obj_num, VT_VBO_all_cell.data());
	glGenBuffers(total_obj_num, VT_EBO_all_cell.data());
	glGenVertexArrays(total_obj_num, EE_VAO_all_cell.data());
	glGenBuffers(3 * total_obj_num, EE_VBO_all_cell.data());
	glGenBuffers(total_obj_num, EE_EBO_all_cell.data());
}

