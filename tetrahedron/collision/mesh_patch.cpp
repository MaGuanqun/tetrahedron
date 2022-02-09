#include"mesh_patch.h"



void MeshPatch::initialPatch(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
	Thread* thread, std::vector<std::vector<std::vector<unsigned int>>>* triangle_patch)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	this->thread = thread;
	this->triangle_patch = triangle_patch;
	setInObjectData();
}



void MeshPatch::setInObjectData()
{
	SpatialHashing spatial_hashing;
	double tolerance[4] = { 2.0,2.0,2.0,2.0 }; //here, spatial hashing cell length =tolerance * max_length
	spatial_hashing.setInObject(cloth, collider, tetrahedron, thread, tolerance,8,true);	
	spatial_hashing.buildSpatialHashing();
	spatial_hashing.findPatch(triangle_patch);
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
	std::vector<unsigned int>* patch = triangle_patch->data()[obj_index].data();
	for (int i = 0; i < size;++i) {
		for (int j = 0; j < 3; ++j) {
			pos[3 * i + j] = { position[indices[i][j]][0], position[indices[i][j]][1], position[indices[i][j]][2] };
		}		
	}
	unsigned int tri_index;
	std::array<double, 3> color_;

	int resi;
	unsigned int k = 256 * 256 * 256 / triangle_patch->data()[obj_index].size();
	unsigned int color_index;
	for (unsigned int i = 0; i < triangle_patch->data()[obj_index].size(); ++i) {
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

}

void MeshPatch::test()
{
	for (int i = 0; i < triangle_patch->data()[0].size(); ++i) {
		for (int j = 0; j < triangle_patch->data()[0][i].size(); ++j) {
			std::cout << triangle_patch->data()[0][i][j] << " ";
		}
		std::cout << std::endl;
	}
}

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