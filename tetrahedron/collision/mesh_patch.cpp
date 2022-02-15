#include"mesh_patch.h"

MeshPatch::~MeshPatch()
{
	delete[] patch_is_intersect;
	delete shader;
}

void MeshPatch::initialPatch(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
	Thread* thread)//std::vector<std::vector<std::vector<unsigned int>>>* triangle_patch, std::vector<std::vector<std::vector<unsigned int>>>* patch_vertex
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	this->thread = thread;
	total_obj_num = cloth->size() + tetrahedron->size() + collider->size();
	tetrahedron_end_index = cloth->size() + tetrahedron->size();
	total_thread_num = thread->thread_num;
	setInObjectData();

	std::cout << "patch " << triangle_patch[0].size() << std::endl;
}



void MeshPatch::setInObjectData()
{
	SpatialHashing spatial_hashing;
	double tolerance[4] = { 1.0,1.0,1.0,1.0 }; //here, spatial hashing cell length =tolerance * max_length
	spatial_hashing.setInObject(cloth, collider, tetrahedron, thread, tolerance,8,true);	
	spatial_hashing.buildSpatialHashing();
	spatial_hashing.findPatch(&triangle_patch);
	findVertex();
	//std::cout << triangle_patch.data()[0].size()<<" "<<
	//	(double)tetrahedron->data()[0].mesh_struct.triangle_indices.size()/(double)triangle_patch.data()[0].size() << std::endl;
}


void MeshPatch::findVertex()
{
	patch_index_start_per_thread.resize(triangle_patch.size());
	patch_vertex.resize(triangle_patch.size());
	for (unsigned int i = 0; i < patch_index_start_per_thread.size(); ++i) {
		patch_index_start_per_thread[i].resize(total_thread_num + 1);
		patch_vertex.data()[i].resize(triangle_patch.data()[i].size());
		for (unsigned int j = 0; j < patch_vertex.data()[i].size(); ++j) {
			patch_vertex.data()[i][j].reserve(triangle_patch.data()[i][j].size());
		}
	}
	
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		arrangeIndex(total_thread_num, triangle_patch.data()[i].size(), patch_index_start_per_thread[i].data());
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		arrangeIndex(total_thread_num, triangle_patch.data()[i+cloth->size()].size(), patch_index_start_per_thread[i+cloth->size()].data());
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		arrangeIndex(total_thread_num, triangle_patch.data()[i + cloth->size()+tetrahedron->size()].size(), 
			patch_index_start_per_thread[i + cloth->size() + tetrahedron->size()].data());
	}
	thread->assignTask(this, FIND_VERTEX);
	patch_AABB.resize(total_obj_num);
	patch_is_intersect=new bool*[total_obj_num];
	for (int i = 0; i < total_obj_num; ++i) {
		patch_AABB[i].resize(patch_vertex[i].size());
		patch_is_intersect[i]=new bool[patch_vertex[i].size()];
	}
	thread->assignTask(this, PATCH_AABB);
	
	triangle_index_noted_as_true.resize(total_obj_num);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		triangle_index_noted_as_true[i].reserve(cloth->data()[i].mesh_struct.triangle_indices.size());
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		triangle_index_noted_as_true[i + cloth->size()].reserve(tetrahedron->data()[i].mesh_struct.triangle_indices.size());
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
		triangle_index_noted_as_true[i + cloth->size() + tetrahedron->size()].reserve(collider->data()[i].mesh_struct.triangle_indices.size());
	}
	triangle_size_start_per_thread.resize(total_obj_num); 
	triangle_index_noted_begin_per_thread.resize(total_obj_num);
	//triangle_patch_noted_true.resize(total_obj_num);
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		//triangle_patch_noted_true[i].resize(total_thread_num);
		
		triangle_size_start_per_thread[i].resize(total_thread_num+1,0);
		triangle_index_noted_begin_per_thread[i].resize(total_thread_num+1,0);
		//for (int j = 0; j < total_thread_num; ++j) {			
		//		triangle_patch_noted_true[i][j].reserve(triangle_patch[i].size());
		//}
	}
	
}

//FIND_VERTEX
void MeshPatch::findVertex(int thread_No)
{	
	int vertex_num;
	std::array<int, 3>* triangle_indices;
	std::vector<bool>is_vertex_used;
	for (int i = 0; i < total_obj_num; ++i) {
		if (i < cloth->size()) {
			vertex_num = cloth->data()[i].mesh_struct.vertex_position.size();
			triangle_indices = cloth->data()[i].mesh_struct.triangle_indices.data();
			is_vertex_used.resize(vertex_num, false);
			findVertex(is_vertex_used, triangle_indices, triangle_patch.data()[i].data(),
				patch_index_start_per_thread[i][thread_No], patch_index_start_per_thread[i][thread_No + 1], patch_vertex.data()[i].data());
		}
		else if (i < tetrahedron_end_index) {
			vertex_num = tetrahedron->data()[i-cloth->size()].mesh_struct.vertex_index_on_sureface.size();
			triangle_indices = tetrahedron->data()[i - cloth->size()].mesh_struct.triangle_indices.data();
			is_vertex_used.resize(vertex_num, false);
			findTetVertex(is_vertex_used, triangle_indices, triangle_patch.data()[i].data(),
				patch_index_start_per_thread[i][thread_No], patch_index_start_per_thread[i][thread_No + 1], 
				patch_vertex.data()[i].data(), tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_surface_index);
		}
		else {
			vertex_num = collider->data()[i - tetrahedron_end_index].mesh_struct.vertex_position.size();
			triangle_indices = collider->data()[i - tetrahedron_end_index].mesh_struct.triangle_indices.data();
			is_vertex_used.resize(vertex_num, false);
			findVertex(is_vertex_used, triangle_indices, triangle_patch.data()[i].data(),
				patch_index_start_per_thread[i][thread_No], patch_index_start_per_thread[i][thread_No + 1], patch_vertex.data()[i].data());
		}
	
	}
}


//PATCH_AABB
void MeshPatch::obtainAABB(int thread_No)
{
	std::array<double, 6>* aabb;
	unsigned int start, end;
	std::vector<unsigned int>* vertex_patch;
	std::array<double, 6>* vertex_aabb;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		if (i < cloth->size()) {
			vertex_aabb = cloth->data()[i].vertex_AABB.data();
		}
		else if (i < tetrahedron_end_index) {
			vertex_aabb = tetrahedron->data()[i-cloth->size()].vertex_AABB.data();
		}
		else {
			vertex_aabb = collider->data()[i - tetrahedron_end_index].vertex_AABB.data();
		}
		aabb = patch_AABB[i].data();
		start = patch_index_start_per_thread[i][thread_No];
		end = patch_index_start_per_thread[i][thread_No + 1];
		vertex_patch = patch_vertex.data()[i].data();
		obtainAABB(aabb, start, end, vertex_patch, vertex_aabb);
		memset(patch_is_intersect[i] + start, 0, (end - start));
	}

}

void MeshPatch::obtainAABB(std::array<double, 6>* aabb, unsigned int start, unsigned int end, std::vector<unsigned int>* vertex_patch,
	std::array<double, 6>* vertex_aabb)
{
	for (unsigned int i = start; i < end; ++i) {
		aabb[i] = vertex_aabb[vertex_patch[i][0]];
		for (unsigned int j = 1; j < vertex_patch[i].size(); ++j) {
			getAABB(aabb[i].data(), vertex_aabb[vertex_patch[i][j]].data());
		}
	}
}

void MeshPatch::findVertex(std::vector<bool>& is_vertex_used, std::array<int,3>* triangle_indices,
	std::vector<unsigned int>* patch, unsigned int start, unsigned int end, std::vector<unsigned int>* patch_vertex)
{
	unsigned int triangle_index;
	for (unsigned int i = start; i < end; ++i) {
		for (unsigned int j = 0; j < patch[i].size(); ++j) {
			triangle_index = patch[i][j];
			for (unsigned int k = 0; k < 3; ++k) {
				if (!is_vertex_used[triangle_indices[triangle_index][k]]) {
					patch_vertex[i].push_back(triangle_indices[triangle_index][k]);
					is_vertex_used[triangle_indices[triangle_index][k]] = true;
				}
			}
		}
		for (unsigned int j = 0; j < patch_vertex[i].size(); ++j) {
			is_vertex_used[patch_vertex[i][j]] = false;
		}
		patch_vertex[i].shrink_to_fit();
	}
}

void MeshPatch::findTetVertex(std::vector<bool>& is_vertex_used, std::array<int, 3>* triangle_indices,
	std::vector<unsigned int>* patch, unsigned int start, unsigned int end, std::vector<unsigned int>* patch_vertex, 
	std::vector<int>& surface_vertex_index)
{
	unsigned int triangle_index;
	for (unsigned int i = start; i < end; ++i) {
		for (unsigned int j = 0; j < patch[i].size(); ++j) {
			triangle_index = patch[i][j];
			for (unsigned int k = 0; k < 3; ++k) {
				if (!is_vertex_used[surface_vertex_index[triangle_indices[triangle_index][k]]]) {
					patch_vertex[i].push_back(surface_vertex_index[triangle_indices[triangle_index][k]]);
					is_vertex_used[surface_vertex_index[triangle_indices[triangle_index][k]]] = true;
				}
			}
		}
		for (unsigned int j = 0; j < patch_vertex[i].size(); ++j) {
			is_vertex_used[patch_vertex[i][j]] = false;
		}
		patch_vertex[i].shrink_to_fit();
	}
}

void MeshPatch::genBuffer()
{
	delete shader;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(2, VBO);
	shader = new Shader("./shader/mesh_patch.vs", "./shader/mesh_patch.fs");
}


void MeshPatch::setBuffer(unsigned int obj_index, unsigned int tetrahedron_start_index)
{
	genBuffer();
	std::vector<std::array<double, 3>> pos;
	std::vector<std::array<double, 3>> color;
	unsigned int size;
	if (obj_index < tetrahedron_start_index) {
		size = cloth->data()[obj_index].mesh_struct.triangle_indices.size();
	}
	else {
		size = tetrahedron->data()[obj_index- tetrahedron_start_index].mesh_struct.triangle_indices.size();
	}	
	pos.resize(3 * size);
	color.resize(3 * size);
	std::array<int, 3>* indices; std::array<double, 3>* position;
	if (obj_index < tetrahedron_start_index) {
		indices = cloth->data()[obj_index].mesh_struct.triangle_indices.data();
		position = cloth->data()[obj_index].mesh_struct.vertex_position.data();
	}
	else {
		indices = tetrahedron->data()[obj_index - tetrahedron_start_index].mesh_struct.triangle_indices.data();
		position = tetrahedron->data()[obj_index - tetrahedron_start_index].mesh_struct.vertex_position.data();
	}
	std::vector<unsigned int>* patch = triangle_patch.data()[obj_index].data();
	for (int i = 0; i < size;++i) {
		for (int j = 0; j < 3; ++j) {
			pos[3 * i + j] = { position[indices[i][j]][0], position[indices[i][j]][1], position[indices[i][j]][2] };
		}		
	}
	unsigned int tri_index;
	std::array<double, 3> color_;

	int resi;
	unsigned int k = 256 * 256 * 256 / triangle_patch.data()[obj_index].size();
	unsigned int color_index;
	for (unsigned int i = 0; i < triangle_patch.data()[obj_index].size(); ++i) {
		color_index = i * k;
		resi = color_index % 65536;
		color_ = { 1.0 - (double)(color_index / 65536) / 255.0, (double)(resi / 256) / 255.0, (double)(resi % 256) / 255.0 };	
		for (int j = 0; j < patch[i].size(); ++j) {
			tri_index = patch[i][j];
			for (int k = 0; k < 3; ++k) {
				color[3 * tri_index + k] = color_;
			}
		}
		//std::cout << color_[0] << " " << color_[1] << " " << color_[2] << std::endl;
	}

	
	//test();
	draw_element_num = pos.size();
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
	glBufferData(GL_ARRAY_BUFFER, pos.size() * sizeof(std::array<double, 3>), pos[0].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
	glBufferData(GL_ARRAY_BUFFER, color.size() * sizeof(std::array<double, 3>), color[0].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);	

	//test(obj_index, tetrahedron_start_index);

}

void MeshPatch::test(unsigned int obj_index, unsigned int tetrahedron_start_index)
{
	int triangle_size;
	if (obj_index < tetrahedron_start_index) {
		triangle_size = cloth->data()[obj_index].mesh_struct.triangle_indices.size();
	}
	else {
		triangle_size = tetrahedron->data()[obj_index - tetrahedron_start_index].mesh_struct.triangle_indices.size();
	}
	std::vector<bool> is_used(triangle_size, false);
	for (int i = 0; i < triangle_patch.data()[obj_index].size(); ++i) {
		for (int j = 0; j < triangle_patch.data()[obj_index][i].size(); ++j) {
			if (is_used[triangle_patch.data()[obj_index][i][j]]) {
				std::cout << "one triangle occurs in two patches" << std::endl;
			}
			else {
				is_used[triangle_patch.data()[obj_index][i][j]] = true;
			}
		}
	}
	for (int i = 0; i < is_used.size(); ++i) {
		if (!is_used[i]) {
			std::cout << "patches lost index" << std::endl;
		}
	}
}

//void MeshPatch::test()
//{
//	for (int i = 0; i < triangle_patch.data()[0].size(); ++i) {
//		for (int j = 0; j < triangle_patch.data()[0][i].size(); ++j) {
//			std::cout << triangle_patch.data()[0][i][j] << " ";
//		}
//		std::cout << std::endl;
//	}
//}

void MeshPatch::draw(Camera* camera)
{	
	shader->use();
	shader->setMat4("model", glm::mat4(1.0));
	shader->setMat4("projection", camera->GetProjectMatrix());
	shader->setMat4("view", camera->GetViewMatrix());
	glBindVertexArray(VAO);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawArrays(GL_TRIANGLES, 0, draw_element_num);
	glBindVertexArray(0);
	
}

void MeshPatch::getAABB(double* target, double* aabb0)
{
	for (unsigned int i = 0; i < 3; ++i) {
		if (target[i] > aabb0[i]) {
			target[i] = aabb0[i];
		}
	}
	for (unsigned int i = 3; i < 6; ++i) {	
		if (target[i] < aabb0[i]) {
			target[i] = aabb0[i];
		}
	}
}

void MeshPatch::findAllNotedTrueTriangle()
{
	thread->assignTask(this, DECIDE_TRIANGLE_INDEX_SIZE);
	for (unsigned int j = 0; j < total_obj_num; ++j) {
		for (unsigned int i = 2; i < total_thread_num + 1; ++i) {
			triangle_size_start_per_thread[j][i] += triangle_size_start_per_thread[j][i - 1];
		}
		triangle_index_noted_as_true[j].resize(triangle_size_start_per_thread[j][total_thread_num]);
	}

	thread->assignTask(this, FIND_TRIANGLE_INDEX);

	for (unsigned int j = 0; j < total_obj_num; ++j) {
		arrangeIndex(total_thread_num, triangle_index_noted_as_true[j].size(), triangle_index_noted_begin_per_thread[j].data());
	}
}

//decide noted true triangle size
//DECIDE_TRIANGLE_INDEX_SIZE
void MeshPatch::decideTriangleIndexSize(int thread_No)
{
	unsigned int end;
	bool* is_intersect;
	unsigned int* tri_size;
	//std::vector<unsigned int>* patch;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		tri_size = &triangle_size_start_per_thread[i][thread_No+1];
		(*tri_size) = 0;
		is_intersect = patch_is_intersect[i];
		end = patch_index_start_per_thread[i][thread_No + 1];
		//patch = &triangle_patch_noted_true[i][thread_No];
		//patch->clear();
		for (int j = patch_index_start_per_thread[i][thread_No]; j < end; ++j) {
			if (is_intersect[j]) {
				//patch->push_back(j);
				(*tri_size) += triangle_patch[i][j].size();
			}
		}
	}
}

//FIND_TRIANGLE_INDEX
void MeshPatch::findNotedTrueTriangleIndex(int thread_No)
{
	unsigned int end;
	bool* is_intersect;
	unsigned int* tri_start;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		is_intersect = patch_is_intersect[i];
		end = patch_index_start_per_thread[i][thread_No + 1];
		tri_start = triangle_index_noted_as_true[i].data() + triangle_size_start_per_thread[i][thread_No];
		for (int j = patch_index_start_per_thread[i][thread_No]; j < end; ++j) {
			if (is_intersect[j]) {
				memcpy(tri_start, triangle_patch[i][j].data(), 4 * triangle_patch[i][j].size());
				tri_start += triangle_patch[i][j].size();
			}
		}
	}
}


void MeshPatch::findAllNotedTrueTriangleSingleThread()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		triangle_index_noted_as_true[i].clear();
		for (int j = 0; j < triangle_patch[i].size(); ++j) {
			if (patch_is_intersect[i][j]) {
				triangle_index_noted_as_true[i].insert(triangle_index_noted_as_true[i].end(), triangle_patch[i][j].begin(),
					triangle_patch[i][j].end());
			}
		}
	}
}

void MeshPatch::test()
{
	time_t t = clock();
	for (unsigned int i = 0; i < 1000; ++i) {
		findAllNotedTrueTriangleSingleThread();
	}
	std::cout << "single thread time " << clock() - t << std::endl;
	std::vector<std::vector<unsigned int>> triangle_index_noted_as_true_;
	triangle_index_noted_as_true_ = triangle_index_noted_as_true;
	t = clock();
	for (unsigned int i = 0; i < 1000; ++i) {
		findAllNotedTrueTriangle();
	}
	std::cout << "multi thread time " << clock() - t << std::endl;
	for (int i = 0; i < total_obj_num; ++i) {
		if (triangle_index_noted_as_true_[i].size() != triangle_index_noted_as_true[i].size()) {
			std::cout << "error size"<< triangle_index_noted_as_true_[i].size()<<" "<< triangle_index_noted_as_true[i].size() << std::endl;
		}
		for (int j = 0; j < triangle_index_noted_as_true_[i].size(); ++j) {
			if (triangle_index_noted_as_true_[i][j] != triangle_index_noted_as_true[i][j]) {
				std::cout << " error index " << triangle_index_noted_as_true_[i][j] << " " << triangle_index_noted_as_true[i][j] << std::endl;
			}
		}
	}
	std::cout<<"finished " << triangle_index_noted_as_true_[0].size() << std::endl;
}