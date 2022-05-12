#include"draw_spatial_hashing.h"

DrawSpatialHashing::DrawSpatialHashing()
{
	setCellIndexBasic();
	genBuffer();
}


void DrawSpatialHashing::setCellData(std::vector<std::vector<unsigned int>>* hash, double cell_length, unsigned int* cell_number,
	double* initial_aabb)
{
	obtianHashValue(hash);

	cell_num0_cell_num1 = cell_number[0] * cell_number[1];
	point_value.resize(24 * hash_index.size());
	edge_index.resize(24 * hash_index.size());



	unsigned int* edge_index_start = edge_index.data();



	unsigned int count = 0;
	for(std::unordered_set<unsigned int>::iterator i = hash_index.begin(); i != hash_index.end();++i){
		reverseHash(*i, cell_length, cell_number, initial_aabb, point_value.data() + 24 * count);
		for (unsigned int j = 0; j < 24; ++j) {
			(*edge_index_start) = cell_index_basic[j] + 8 * count;
			edge_index_start++;	
		}
		count++;
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

void DrawSpatialHashing::obtianHashValue(std::vector<std::vector<unsigned int>>* hash)
{
	hash_index.clear();
	hash_index.reserve(hash->data()[0].size() / 2);	
	for (unsigned int i = 0; i < hash->size(); ++i) {
		hash_index.insert(hash->data()[i].begin(), hash->data()[i].end());
	}
}




void DrawSpatialHashing::reverseHash(unsigned int hash_index, double cell_length, unsigned int* cell_number,
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
				//std::cout << point[3 * count] << " " << point[3 * count + 1] << " " << point[3 * count + 2] << std::endl;
				count++;
			}
		}
	}

}

void DrawSpatialHashing::setCellIndexBasic()
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

void DrawSpatialHashing::genBuffer()
{
	glGenVertexArrays(1, &VAO1);
	glGenBuffers(1, &VBO1);
	glGenBuffers(1, &EBO1);
}

void DrawSpatialHashing::drawCell(Camera* camera, Shader* shader)
{
	shader->use();

	shader->setVec3("color", glm::vec3(0.98f, 0.54f, 0.27f));
	shader->setMat4("model", glm::mat4(1.0));
	shader->setMat4("projection", camera->GetProjectMatrix());
	shader->setMat4("view", camera->GetViewMatrix());
	shader->setFloat("transparent", 1.0f);
	glBindVertexArray(VAO1);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawElements(GL_LINES, edge_index.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}
