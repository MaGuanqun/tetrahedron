#include"drawCulling.h"

void DrawCulling::initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
	Thread* thread)//std::vector<std::vector<std::vector<unsigned int>>>* triangle_patch, std::vector<std::vector<std::vector<unsigned int>>>* patch_vertex
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	this->thread = thread;
	total_obj_num = cloth->size() + tetrahedron->size() + collider->size();
	tetrahedron_end_index = cloth->size() + tetrahedron->size();
	thread_num = thread->thread_num;
	initialUsedIndicator();
	setPalette();

	light.ambient = glm::vec3(0.2, 0.2, 0.2);
	light.diffuse = glm::vec3(0.7, 0.7, 0.7);
	light.specular = glm::vec3(0.8, 0.8, 0.8);

	genBuffer();

	total_displacement.resize(total_obj_num);
	memset(total_displacement[0].data(), 0, 24 * total_obj_num);

	setCellIndexBasic();
}


void DrawCulling::setPalette()
{
	palette.resize(8);
	palette[0] = { 248.0 / 255.0,53.0 / 255.0,146.0 / 255.0 };
	palette[1] = { 252.0 / 255.0,228.0 / 255.0,74.0 / 255.0 };
	palette[2] = { 61.0 / 255.0,147.0 / 255.0,23.0 / 255.0 };
	palette[3] = { 155.0 / 255.0,213.0 / 255.0,241.0 / 255.0 };
	palette[4] = { 47.0 / 255.0,56.0 / 255.0,151.0 / 255.0 };
	palette[5] = { 221.0 / 255.0,206.0 / 255.0,236.0 / 255.0 };
	palette[6] = { 117.0 / 255.0,51.0 / 255.0,7.0 / 255.0 };
	palette[7] = { 56.0 / 255.0,29.0 / 255.0,23.0 / 255.0 };
}

void DrawCulling::initialUsedIndicator()
{
	int total_triangle_num = 0;

	obj_is_used0 = new bool** [thread_num];
	obj_is_used1 = new bool** [thread_num];
	for (int i = 0; i < thread_num; ++i) {
		obj_is_used0[i] = new bool* [tetrahedron_end_index];
		obj_is_used1[i] = new bool* [tetrahedron_end_index];
		for (int j = 0; j < cloth->size(); ++j) {
			obj_is_used0[i][j] = new bool[(*cloth)[j].mesh_struct.triangle_indices.size()];
			obj_is_used1[i][j] = new bool[(*cloth)[j].mesh_struct.triangle_indices.size()];
			memset(obj_is_used0[i][j], 0, (*cloth)[j].mesh_struct.triangle_indices.size());
			memset(obj_is_used1[i][j], 0, (*cloth)[j].mesh_struct.triangle_indices.size());
		}
		for (int j = 0; j < tetrahedron->size(); ++j) {
			obj_is_used0[i][j + cloth->size()] = new bool[(*tetrahedron)[j].mesh_struct.triangle_indices.size()];
			obj_is_used1[i][j + cloth->size()] = new bool[(*tetrahedron)[j].mesh_struct.triangle_indices.size()];
			memset(obj_is_used0[i][j + cloth->size()], 0, (*tetrahedron)[j].mesh_struct.triangle_indices.size());
			memset(obj_is_used1[i][j + cloth->size()], 0, (*tetrahedron)[j].mesh_struct.triangle_indices.size());
		}
	}
	collider_is_used0 = new bool** [thread_num];
	if (!collider->empty()) {
		for (int i = 0; i < thread_num; ++i) {
			collider_is_used0[i] = new bool* [collider->size()];
			for (int j = 0; j < collider->size(); ++j) {
				collider_is_used0[i][j] = new bool[(*collider)[j].mesh_struct.triangle_indices.size()];
				memset(collider_is_used0[i][j], 0, (*collider)[j].mesh_struct.triangle_indices.size());
			}
		}
	}

	for (int j = 0; j < cloth->size(); ++j) {
		total_triangle_num += (*cloth)[j].mesh_struct.triangle_indices.size();
	}
	for (int j = 0; j < tetrahedron->size(); ++j) {
		total_triangle_num += (*tetrahedron)[j].mesh_struct.triangle_indices.size();
	}
	for (int j = 0; j < collider->size(); ++j) {
		total_triangle_num += (*collider)[j].mesh_struct.triangle_indices.size();
	}

	pos_thread.resize(thread_num);
	color_thread.resize(thread_num);
	normal_thread.resize(thread_num);
	for (unsigned int i = 0; i < thread_num; ++i) {
		pos_thread[i].reserve(total_triangle_num);
		color_thread[i].reserve(total_triangle_num);
		normal_thread[i].reserve(total_triangle_num);
	}

	vertex_position.reserve(tetrahedron_end_index);
	triangle_indices.reserve(tetrahedron_end_index);
	vertex_normal_render.reserve(tetrahedron_end_index);
	

	for (int j = 0; j < cloth->size(); ++j) {
		vertex_position.push_back((*cloth)[j].mesh_struct.vertex_position.data());
		triangle_indices.push_back((*cloth)[j].mesh_struct.triangle_indices.data());
		vertex_normal_render.push_back((*cloth)[j].mesh_struct.vertex_normal_for_render.data());
	}
	for (int j = 0; j < tetrahedron->size(); ++j) {
		vertex_position.push_back((*tetrahedron)[j].mesh_struct.vertex_position.data());
		triangle_indices.push_back((*tetrahedron)[j].mesh_struct.triangle_indices.data());
		vertex_normal_render.push_back((*tetrahedron)[j].mesh_struct.vertex_normal_for_render.data());
	}

	for (int j = 0; j < collider->size(); ++j) {
		collider_vertex_position.push_back((*collider)[j].mesh_struct.vertex_position.data());
		collider_triangle_indices.push_back((*collider)[j].mesh_struct.triangle_indices.data());
		collider_vertex_normal_render.push_back((*collider)[j].mesh_struct.vertex_normal_for_render.data());
	}


	pos_start_thread.resize(thread_num);
	pos_start_thread[0] = 0;
}


void DrawCulling::move(unsigned int obj_No, double* displacement)
{
	this->displacement = displacement;
	thread->assignTask(this, MOVE_OBJECT,obj_No);

	SUM_(total_displacement[obj_No], displacement);
	std::cout << "======" << std::endl;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		std::cout << total_displacement[i][0] << ", " << total_displacement[i][1] << ", " << total_displacement[i][2] << std::endl;
	}
}

void DrawCulling::moveScript(unsigned int type)
{
	displacement_test.resize(2);
	switch (time_step)
	{
	case 0:
		if (type == 1) {
			displacement_test[0] = { -0.308742, 0.640283, 0.0271892 };
			displacement_test[1] = { -0.52193, 0.580719, 0.0493092 };
		}
		else if(type==2){
			displacement_test[0] = { 0, 0, 0 };
			displacement_test[1] = { -0.239594, -0.0259783, -0.33549 };
		}
		else if (type == 3) {
			displacement_test[0] = { -0.308742, 0.640283, 0.0271892 };
			displacement_test[1] = { -0.52193, 0.580719, 0.0323092 };
		}
		break;
	case 1:
		if (type == 1) {		
			double a[3] = { -0.510271, 0.631315, 0.0321853 };
			double b[3]= { -0.52193, 0.580719, 0.0493092 };
			SUB(displacement_test[0], a, displacement_test[0]);
			SUB(displacement_test[1], b, displacement_test[1]);
		}
		else if (type == 2) {
			double a[3] = { 0, 0, 0 };
			double b[3] = { -0.0776213, -0.0571768, -0.140402 };
			SUB(displacement_test[0], a, displacement_test[0]);
			SUB(displacement_test[1], b, displacement_test[1]);
		}
		else if (type == 3) {
			double a[3] = { -0.510271, 0.631315, 0.0321853 };
			double b[3] = { -0.52193, 0.580719, 0.0493092 };
			SUB(displacement_test[0], a, displacement_test[0]);
			SUB(displacement_test[1], b, displacement_test[1]);
		}
		break;
	}
	for (unsigned int i = 0; i < displacement_test.size(); ++i) {
		move(i, displacement_test[i].data());
	}
	time_step++;
}




//MOVE_OBJECT
void DrawCulling::move(int thread_No, unsigned int obj_No)
{
	std::array<double, 3>* position;
	unsigned int vertex_start, vertex_end;

	std::array<double, 3>* render_position;
	if (obj_No < cloth->size()) {
		position = cloth->data()[obj_No].mesh_struct.vertex_for_render.data();
		render_position = cloth->data()[obj_No].mesh_struct.vertex_position.data();
		vertex_start = cloth->data()[obj_No].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = cloth->data()[obj_No].mesh_struct.vertex_index_begin_per_thread[thread_No+1];
	}
	else if (obj_No < tetrahedron_end_index) {
		position = tetrahedron->data()[obj_No- cloth->size()].mesh_struct.vertex_for_render.data();
		render_position = tetrahedron->data()[obj_No- cloth->size()].mesh_struct.vertex_position.data();
		vertex_start = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	}
	else {
		position = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_for_render.data();
		render_position = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_position.data();
		vertex_start = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	}
	for (unsigned int i = vertex_start; i < vertex_end; ++i) {
		SUM_(position[i], displacement);
	}
	memcpy(render_position[vertex_start].data(), position[vertex_start].data(), 24 * (vertex_end - vertex_start));
}


void DrawCulling::setInSpatialHashingValue(unsigned int** spatial_hashing_value,
	unsigned int** spatial_hashing_triangle_index, unsigned int** spatial_hashing_value_collider,
	unsigned int** spatial_hashing_triangle_index_collider, std::vector<std::vector<unsigned int>>* prefix_sum,
	std::vector<std::vector<unsigned int>>* prefix_sum_collider, std::vector<unsigned int>* actual_exist_cell_begin_per_thread)
{
	this->spatial_hashing_value = spatial_hashing_value;
	//this->spatial_hashing_obj_index = spatial_hashing_obj_index;
	this->spatial_hashing_triangle_index = spatial_hashing_triangle_index;

	this->spatial_hashing_value_collider = spatial_hashing_value_collider;
	//this->spatial_hashing_obj_index_collider = spatial_hashing_obj_index_collider;
	this->spatial_hashing_triangle_index_collider = spatial_hashing_triangle_index_collider;

	this->prefix_sum = prefix_sum;
	this->prefix_sum_collider = prefix_sum_collider;

	this->actual_exist_cell_begin_per_thread = actual_exist_cell_begin_per_thread;
}

void DrawCulling::setInSpatialHashingValue(std::vector<unsigned int>** spatial_hashing_cell, std::vector<unsigned int>** spatial_hashing_cell_collider,
	unsigned int hash_cell_count)
{
	this->spatial_hashing_cell = spatial_hashing_cell;
	this->spatial_hashing_cell_collider = spatial_hashing_cell_collider;
	this->hash_cell_count = hash_cell_count;

}








void DrawCulling::drawCell(Camera* camera, Shader* shader)
{
	shader->use();
	shader->setVec3("color", glm::vec3(1.0f,1.0f,0.0f));
	shader->setMat4("model", glm::mat4(1.0));
	shader->setMat4("projection", camera->GetProjectMatrix());
	shader->setMat4("view", camera->GetViewMatrix());
	glBindVertexArray(VAO1);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawElements(GL_LINES, edge_index.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}



void DrawCulling::draw(Camera* camera, float& far_plane)
{
	shader->use();
	shader->setInt("depthMap", 0);
	shader->setFloat("far_plane", far_plane);
	shader->setVec3("lightPos", camera->position);
	shader->setVec3("light.ambient", light.ambient);
	shader->setVec3("light.diffuse", light.diffuse);
	shader->setVec3("light.specular", light.specular);
	shader->setVec3("viewPos", camera->position);
	shader->setFloat("transparence", 1.0);
	shader->setMat4("model", glm::mat4(1.0));
	shader->setMat4("projection", camera->GetProjectMatrix());
	shader->setMat4("view", camera->GetViewMatrix());
	glBindVertexArray(VAO);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawArrays(GL_TRIANGLES, 0, pos.size());
	glBindVertexArray(0);

}


void DrawCulling::setThreadDataTogether()
{
	unsigned int size = pos_thread[0].size();
	for (unsigned int i = 1; i < thread_num; ++i) {
		size += pos_thread[i].size();
		pos_start_thread[i] = pos_start_thread[i - 1] + pos_thread[i-1].size();
	}
	pos.resize(size);
	color.resize(size);
	normal.resize(size);
	thread->assignTask(this, SET_DATA_TOGETHER, 0);
}


//SET_DATA_TOGETHER
void DrawCulling::setThreadDataTogether(int thread_No)
{
	memcpy(pos[pos_start_thread[thread_No]].data(), pos_thread[thread_No][0].data(), 24 * pos_thread[thread_No].size());
	memcpy(color[pos_start_thread[thread_No]].data(), color_thread[thread_No][0].data(), 24 * pos_thread[thread_No].size());
	memcpy(normal[pos_start_thread[thread_No]].data(), normal_thread[thread_No][0].data(), 24 * pos_thread[thread_No].size());
}

bool DrawCulling::checkIfExistTwoObject(unsigned int cell_index)
{
	unsigned int current_thread_triangle_num;
	int index = -1;

	if (!collider->empty()) {
		for (unsigned int t = 0; t < thread_num; ++t) {
			current_thread_triangle_num = prefix_sum_collider->data()[t][(cell_index << 1) + 1] - prefix_sum_collider->data()[t][cell_index << 1];
			if (current_thread_triangle_num > 0) {
				return true;
			}		
		}
	}
	for (unsigned int t = 0; t < thread_num; ++t) {
		current_thread_triangle_num = prefix_sum->data()[t][cell_index + 1] - prefix_sum->data()[t][cell_index];
		if (current_thread_triangle_num > 0) {
			for (unsigned int j = prefix_sum->data()[t][cell_index] + 1; j < prefix_sum->data()[t][cell_index + 1]; j+=2) {
				if (index == -1) {
					index = spatial_hashing_triangle_index[t][j];
				}
				else {
					if (spatial_hashing_triangle_index[t][j] != index) {
						return true;
					}
				}
			}
		}
	}
	return false;
}


bool DrawCulling::checkIfExistTwoObjectSP(unsigned int cell_index)
{
	unsigned int current_thread_triangle_num;
	int index = -1;

	if (!collider->empty()) {
		for (unsigned int t = 0; t < thread_num; ++t) {
			if(!spatial_hashing_cell_collider[t][cell_index].empty())
			return true;
		}
	}
	for (unsigned int t = 0; t < thread_num; ++t) {
		if (!spatial_hashing_cell[t][cell_index].empty()) {
			for (unsigned int j = 1; j < spatial_hashing_cell[t][cell_index].size(); j+=2) {
				if (index == -1) {
					index = spatial_hashing_cell[t][cell_index][j];
				}
				else {
					if (spatial_hashing_cell[t][cell_index][j] != index) {
						return true;
					}
				}
			}
		}
	}
	return false;
}


void DrawCulling::drawAABBIntersectBetweenObjects()
{
	thread->assignTask(this, SET_POSITION_COLOR, 0);
	setThreadDataTogether();
//	findTheMaxNumberOfTriangleInCell();
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
	glBufferData(GL_ARRAY_BUFFER, pos.size() * sizeof(std::array<double, 3>), pos[0].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
	glBufferData(GL_ARRAY_BUFFER, color.size() * sizeof(std::array<double, 3>), color[0].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[2]);
	glBufferData(GL_ARRAY_BUFFER, normal.size() * sizeof(std::array<double, 3>), normal[0].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);

	
}




void DrawCulling::setTetrahedronVertex()
{
	std::array<double, 3>* pos_;
	int* indices;
	//tet_triangle_indices.reserve(tetrahedron->data()[0].mesh_struct.indices.size());
	tet_triangle_indices.clear();
	for (unsigned int i = 0; i < thread_num; ++i) {
		for (unsigned int j = 2; j < vertex_tet_pair[i].size(); j+=4) {			
			indices = tetrahedron->data()[0].mesh_struct.indices[vertex_tet_pair[i][j]].data();			
			tet_triangle_indices.emplace_back(indices[0]);
			tet_triangle_indices.emplace_back(indices[2]);
			tet_triangle_indices.emplace_back(indices[1]);

			tet_triangle_indices.emplace_back(indices[0]);
			tet_triangle_indices.emplace_back(indices[3]);
			tet_triangle_indices.emplace_back(indices[2]);

			tet_triangle_indices.emplace_back(indices[0]);
			tet_triangle_indices.emplace_back(indices[1]);
			tet_triangle_indices.emplace_back(indices[3]);

			tet_triangle_indices.emplace_back(indices[1]);
			tet_triangle_indices.emplace_back(indices[2]);
			tet_triangle_indices.emplace_back(indices[3]);
		}
	}
	//std::cout << vertex_tet_pair[0].size() << std::endl;
	//std::cout << tet_triangle_indices.size() << std::endl;
	//for (unsigned int i = 0; i < tet_triangle_indices.size(); ++i) {
	//	std::cout << tet_triangle_indices[i] << std::endl;
	//}
	glBindVertexArray(VAO2);
	glBindBuffer(GL_ARRAY_BUFFER, VBO2);
	glBufferData(GL_ARRAY_BUFFER, tetrahedron->data()[0].mesh_struct.vertex_for_render.size() * sizeof(std::array<double,3>), tetrahedron->data()[0].mesh_struct.vertex_for_render[0].data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO2);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, tet_triangle_indices.size() * sizeof(unsigned int), tet_triangle_indices.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);
}

void DrawCulling::findTheMaxNumberOfTriangleInCell()
{
	unsigned int num1 = 0;
	unsigned int max_num = 0;
	unsigned int num = 0;
	unsigned int hash_start_index;
	for (unsigned int i = 1; i < prefix_sum->data()[0].size(); ++i) {
		num1 = 0;
		for (unsigned int j = 0; j < thread_num; ++j) {
			num += (prefix_sum->data()[j][i] - prefix_sum->data()[j][i - 1]) >> 1;
			num1 += (prefix_sum->data()[j][i] - prefix_sum->data()[j][i - 1]) >> 1;
		}
		if (num1 > 500) {
			hash_start_index = i-1;
			std::cout << "record triangle num " << num1 << std::endl;
			break;
		}
	}
	for (unsigned int j = 0; j < thread_num; ++j) {
		if ((prefix_sum->data()[j][hash_start_index + 1] - prefix_sum->data()[j][hash_start_index]) > 0) {
			hash_index = spatial_hashing_value[j][prefix_sum->data()[j][hash_start_index] >> 1];
			break;
		}
	}
	std::cout << "hash index " << hash_index << " " << hash_start_index << std::endl;
	
	pos.clear();
	color.clear();
	normal.clear();
	unsigned int obj_index, triangle_index;

	bool** obj_is_used0_ = obj_is_used0[0];
	bool** obj_is_used1_ = obj_is_used1[0];

	for (unsigned int i = 0; i < tetrahedron_end_index; ++i) {
		if (i < cloth->size()) {
			memset(obj_is_used0_[i], 0, cloth->data()[i].mesh_struct.triangle_indices.size());
			memset(obj_is_used1_[i], 0, cloth->data()[i].mesh_struct.triangle_indices.size());
		}
		else {
			memset(obj_is_used0_[i], 0, tetrahedron->data()[i - cloth->size()].mesh_struct.triangle_indices.size());
			memset(obj_is_used1_[i], 0, tetrahedron->data()[i - cloth->size()].mesh_struct.triangle_indices.size());
		}
	}

	for (unsigned int t = 0; t < thread_num; ++t) {
		for (unsigned int j = prefix_sum->data()[t][hash_start_index]; j < prefix_sum->data()[t][hash_start_index+1]; j += 2) {
			obj_index = spatial_hashing_triangle_index[t][j + 1];
			triangle_index = spatial_hashing_triangle_index[t][j];
			if (!obj_is_used0_[obj_index][triangle_index]) {
				obj_is_used0_[obj_index][triangle_index] = true;
				for (unsigned int k = 0; k < 3; ++k) {
					pos.emplace_back(vertex_position[obj_index][triangle_indices[obj_index][triangle_index][k]]);
					color.emplace_back(palette[hash_index % palette.size()]);
					normal.emplace_back(vertex_normal_render[obj_index][triangle_indices[obj_index][triangle_index][k]]);
				}
			}
		}
		
	}
	std::cout << "actual Triangle Num " << pos.size() / 3 << std::endl;
}

//SET_POSITION_COLOR
void DrawCulling::setAllTriangle(int thread_No)
{

	bool** obj_is_used0_ = obj_is_used0[thread_No];
	bool** obj_is_used1_ = obj_is_used1[thread_No];
	bool** collider_is_used0_ = collider_is_used0[thread_No];
	for (unsigned int i = 0; i < tetrahedron_end_index; ++i) {
		if (i < cloth->size()) {
			memset(obj_is_used0_[i], 0, cloth->data()[i].mesh_struct.triangle_indices.size());
			memset(obj_is_used1_[i], 0, cloth->data()[i].mesh_struct.triangle_indices.size());
		}
		else {
			memset(obj_is_used0_[i], 0, tetrahedron->data()[i - cloth->size()].mesh_struct.triangle_indices.size());
			memset(obj_is_used1_[i], 0, tetrahedron->data()[i - cloth->size()].mesh_struct.triangle_indices.size());
		}
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		memset(collider_is_used0_[i], 0, collider->data()[i].mesh_struct.triangle_indices.size());
	}
	unsigned int last_cell_index = actual_exist_cell_begin_per_thread->data()[thread_No + 1];

	std::vector<std::array<double, 3>>* pos_=&pos_thread[thread_No];
	std::vector<std::array<double, 3>>* color_ = &color_thread[thread_No];
	std::vector<std::array<double, 3>>* normal_ = &normal_thread[thread_No];

	pos_->clear();
	color_->clear();
	normal_->clear();

	unsigned int current_thread_triangle_num;

	unsigned int obj_index;
	unsigned int triangle_index;
	unsigned int* tri_vertex_index;
	std::array<double, 3>* pos;

	for (unsigned int cell_index = actual_exist_cell_begin_per_thread->data()[thread_No];
		cell_index < last_cell_index; ++cell_index) {
		if (checkIfExistTwoObjectSP(cell_index)) {
			for (unsigned int t = 0; t < thread_num; ++t) {
				if (!spatial_hashing_cell[t][cell_index].empty()) {
					for (unsigned int j = 0; j < spatial_hashing_cell[t][cell_index].size(); j+=2) {
						obj_index = spatial_hashing_cell[t][cell_index][j + 1];
						triangle_index= spatial_hashing_cell[t][cell_index][j];
						if (!obj_is_used0_[obj_index][triangle_index]) {
							obj_is_used0_[obj_index][triangle_index] = true;
							for (unsigned int k = 0; k < 3; ++k) {
								pos_->emplace_back(vertex_position[obj_index][triangle_indices[obj_index][triangle_index][k]]);
								color_->emplace_back(palette[cell_index % palette.size()]);
								normal_->emplace_back(vertex_normal_render[obj_index][triangle_indices[obj_index][triangle_index][k]]);
							}
						}
					}
				}
			}
			if (!collider->empty()) {
				for (unsigned int t = 0; t < thread_num; ++t) {
					if (!spatial_hashing_cell_collider[t][cell_index].empty()) {
						for (unsigned int j = 0; j < spatial_hashing_cell_collider[t][cell_index].size(); j+=2) {
							obj_index = spatial_hashing_cell_collider[t][cell_index][j + 1];
							triangle_index = spatial_hashing_cell_collider[t][cell_index][j];
							if (!collider_is_used0_[obj_index][triangle_index]) {
								collider_is_used0_[obj_index][triangle_index] = true;
								for (unsigned int k = 0; k < 3; ++k) {
									pos_->emplace_back(collider_vertex_position[obj_index][collider_triangle_indices[obj_index][triangle_index][k]]);
									color_->emplace_back(palette[cell_index % palette.size()]);
									normal_->emplace_back(collider_vertex_normal_render[obj_index][collider_triangle_indices[obj_index][triangle_index][k]]);
								}
							}
						}
					}
				}
			}
		}
	}

	//for (unsigned int cell_index = actual_exist_cell_begin_per_thread->data()[thread_No];
	//	cell_index < last_cell_index; ++cell_index) {

	//	if (checkIfExistTwoObject(cell_index)) {
	//		//if (checkIfExistTwoObjectSP(cell_index)) {
	//		for (unsigned int t = 0; t < thread_num; ++t) {
	//			current_thread_triangle_num = prefix_sum->data()[t][cell_index + 1] - prefix_sum->data()[t][cell_index];
	//			if (current_thread_triangle_num > 0) {
	//				for (unsigned int j = prefix_sum->data()[t][cell_index]; j < prefix_sum->data()[t][cell_index + 1]; ++j) {
	//					obj_index = spatial_hashing_obj_index[t][j];
	//					triangle_index = spatial_hashing_triangle_index[t][j];
	//					if (!obj_is_used0_[obj_index][triangle_index]) {
	//						obj_is_used0_[obj_index][triangle_index] = true;
	//						for (unsigned int k = 0; k < 3; ++k) {
	//							pos_->emplace_back(vertex_position[obj_index][triangle_indices[obj_index][triangle_index][k]]);
	//							color_->emplace_back(palette[cell_index % palette.size()]);
	//							normal_->emplace_back(vertex_normal_render[obj_index][triangle_indices[obj_index][triangle_index][k]]);
	//						}
	//					}
	//				}
	//			}
	//		}
	//		if (!collider->empty()) {
	//			for (unsigned int t = 0; t < thread_num; ++t) {
	//				current_thread_triangle_num = prefix_sum_collider->data()[t][(cell_index << 1) + 1] - prefix_sum_collider->data()[t][cell_index << 1];
	//				if (current_thread_triangle_num > 0) {
	//					for (unsigned int j = prefix_sum_collider->data()[t][cell_index << 1]; j < prefix_sum_collider->data()[t][(cell_index << 1) + 1]; ++j) {
	//						obj_index = spatial_hashing_obj_index_collider[t][j];
	//						triangle_index = spatial_hashing_triangle_index_collider[t][j];
	//						if (!collider_is_used0_[obj_index][triangle_index]) {
	//							collider_is_used0_[obj_index][triangle_index] = true;
	//							for (unsigned int k = 0; k < 3; ++k) {
	//								pos_->emplace_back(collider_vertex_position[obj_index][collider_triangle_indices[obj_index][triangle_index][k]]);
	//								color_->emplace_back(palette[cell_index % palette.size()]);
	//								normal_->emplace_back(collider_vertex_normal_render[obj_index][collider_triangle_indices[obj_index][triangle_index][k]]);
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
}



void DrawCulling::genBuffer()
{
	glGenVertexArrays(1, &VAO);
	glGenBuffers(3, VBO);
	shader = new Shader("./shader/mesh_patch.vs", "./shader/mesh_patch.fs");

	glGenVertexArrays(1, &VAO1);
	glGenBuffers(1, &VBO1);
	glGenBuffers(1, &EBO1);

	glGenVertexArrays(1, &VAO2);
	glGenBuffers(1, &VBO2);
	glGenBuffers(1, &EBO2);
}


void DrawCulling::setSingleCellData(double cell_length, unsigned int* cell_number,
	double* initial_aabb)
{
	cell_num0_cell_num1 = cell_number[0] * cell_number[1];
	point_value.resize(24);
	edge_index.resize(24);

	unsigned int* edge_index_start = edge_index.data();	
	reverseHash(hash_index, cell_length, cell_number, initial_aabb, point_value.data());
	for (unsigned int j = 0; j < 24; ++j) {
		(*edge_index_start) = cell_index_basic[j];
		edge_index_start++;
	}
	

	glBindVertexArray(VAO1);
	glBindBuffer(GL_ARRAY_BUFFER, VBO1);
	glBufferData(GL_ARRAY_BUFFER, point_value.size() * sizeof(double), point_value.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO1);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, edge_index.size() * sizeof(unsigned int), edge_index.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);
}


void DrawCulling::setCellData(std::vector<unsigned int>* hash, double cell_length, unsigned int* cell_number,
	double* initial_aabb)
{
	cell_num0_cell_num1 = cell_number[0] * cell_number[1];
	point_value.resize(24 * hash->size());
	edge_index.resize(24 * hash->size());

	unsigned int* edge_index_start = edge_index.data();
	for (unsigned int i = 0; i < hash->size(); ++i) {
		reverseHash(hash->data()[i], cell_length, cell_number, initial_aabb, point_value.data() + 24 * i);
		for (unsigned int j = 0; j < 24; ++j) {
			(*edge_index_start) = cell_index_basic[j] + 8 * i;
			edge_index_start++;
		}
	}

	glBindVertexArray(VAO1);
	glBindBuffer(GL_ARRAY_BUFFER, VBO1);
	glBufferData(GL_ARRAY_BUFFER, point_value.size() * sizeof(double), point_value.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO1);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, edge_index.size() * sizeof(unsigned int), edge_index.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);
}


void DrawCulling::drawTetTriangle(Camera* camera, Shader* object_shader_front, Light& light, float& far_plane)
{
	glm::vec3 color0 = { 0.01,0.1,0.01 };
	glm::vec3 color1 = { 0.1,0.4,0.1 };

	if (!tet_triangle_indices.empty()) {
		object_shader_front->use();	
		//object_shader_front->setInt("depthMap", 0);
		//object_shader_front->setFloat("far_plane", far_plane);
		//object_shader_front->setBool("lightIsChosen", true);
		//object_shader_front->setVec3("lightPos", camera->position);
		//object_shader_front->setVec3("light.ambient", light.ambient);
		//object_shader_front->setVec3("light.diffuse", light.diffuse);
		//object_shader_front->setVec3("light.specular", light.specular);
		//object_shader_front->setVec3("viewPos", camera->position);
		//object_shader_front->setBool("lightShadowOn", true);
		object_shader_front->setMat4("projection", camera->GetProjectMatrix());
		object_shader_front->setMat4("view", camera->GetViewMatrix());
		object_shader_front->setMat4("model", glm::mat4(1.0));
		object_shader_front->setVec3("color", color1);
		//std::cout << color1.x << std::endl;
		//object_shader_front->setFloat("transparence", 1.0);
		//object_shader_front->setVec3("material.Kd", color0);
		//object_shader_front->setVec3("material.Ka", color0);
		//object_shader_front->setVec3("material.Ks", color1);
		//object_shader_front->setFloat("material.Ns", 32);
		glBindVertexArray(VAO2);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, tet_triangle_indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}




void DrawCulling::reverseHash(unsigned int hash_index, double cell_length, unsigned int* cell_number,
	double* initial_aabb, double* point)
{
	unsigned int index[8][3];
	index[0][2] = hash_index / cell_num0_cell_num1;
	index[0][1] = (hash_index % cell_num0_cell_num1) / cell_number[0];
	index[0][0] = (hash_index % cell_num0_cell_num1) % cell_number[0];

	index[1][2] = index[0][2] + 1;
	index[1][1] = index[0][1] + 1;
	index[1][0] = index[0][0] + 1;

	unsigned int count = 0;
	for (unsigned int j0 = 0; j0 < 2; ++j0) {
		for (unsigned int j1 = 0; j1 < 2; ++j1) {
			for (unsigned int j2 = 0; j2 < 2; ++j2) {
				point[3 * count] = cell_length * index[j0][0] + initial_aabb[0];
				point[3 * count + 1] = cell_length * index[j1][1] + initial_aabb[1];
				point[3 * count + 2] = cell_length * index[j2][2] + initial_aabb[2];
				count++;
			}
		}
	}
}

void DrawCulling::setCellIndexBasic()
{
	//----6--7
	//----4--5 
	//2--3
	//0--1

	cell_index_basic[0] = 0;
	cell_index_basic[1] = 1;
	cell_index_basic[2] = 0;
	cell_index_basic[3] = 2;
	cell_index_basic[4] = 1;
	cell_index_basic[5] = 3;
	cell_index_basic[6] = 2;
	cell_index_basic[7] = 3;

	cell_index_basic[8] = 4;
	cell_index_basic[9] = 6;
	cell_index_basic[10] = 4;
	cell_index_basic[11] = 5;

	cell_index_basic[12] = 6;
	cell_index_basic[13] = 7;
	cell_index_basic[14] = 7;
	cell_index_basic[15] = 5;

	cell_index_basic[16] = 2;
	cell_index_basic[17] = 6;
	cell_index_basic[18] = 3;
	cell_index_basic[19] = 7;

	cell_index_basic[20] = 0;
	cell_index_basic[21] = 4;
	cell_index_basic[22] = 1;
	cell_index_basic[23] = 5;
}