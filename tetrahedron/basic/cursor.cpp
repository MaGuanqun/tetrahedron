#include"cursor.h"
#include"global.h"

Cursor::Cursor()
{
	latitude_num = 20;
	longitude_num = 40;
	vertice_num= (longitude_num + 1) * (latitude_num - 2) + 2;
	ori_vertices_pos.reserve(vertice_num);
	shader = new Shader("./shader/basic.vs", "./shader/basic.fs");
	*shader = Shader("./shader/basic.vs", "./shader/basic.fs");
	genBuffer();
	light.ambient = glm::vec3(0.3, 0.3, 0.3);
	light.diffuse = glm::vec3(0.7, 0.7, 0.7);
	light.specular = glm::vec3(0.95, 0.95, 0.95);
	lines.resize(2);
}

void Cursor::draw(Camera* camera) {
	setBufferData();
	glViewport(WidthOffset, HeightOffset, SCR_WIDTH, SCR_HEIGHT);
	shader->use();
	shader->setMat4("projection", camera->GetProjectMatrix());
	shader->setMat4("view", camera->GetViewMatrix());
	shader->setMat4("model", glm::mat4(1.0));
	shader->setVec3("color", glm::vec3(0.0, 0.0, 1.0));
	shader->setVec3("lightPos", camera->position);
	shader->setVec3("light.ambient", light.ambient);
	shader->setVec3("light.diffuse", light.diffuse);
	shader->setVec3("light.specular", light.specular);
	shader->setVec3("viewPos", camera->position);
	shader->setFloat("transparence", 0.5);
	//glCullFace(GL_BACK);
	glBindVertexArray(VAO1);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

	glBindVertexArray(VAO2);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

	shader->use();
	shader->setVec3("color", glm::vec3(0.0, 0.0, 1.0));
	shader->setFloat("transparence", 1.0);
	glBindVertexArray(VAO3);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawArrays(GL_LINE_STRIP, 0, 2);
	glBindVertexArray(0);
}


void Cursor::translate(double u[3], double cursor_pos[3])
{
	double trans[3];
	SUB(trans, u, center);
	for (int i = 0; i < vertice_num; ++i) {
		SUM(vertices_pos[i], ori_vertices_pos[i], trans);
	}
	SUB(trans, cursor_pos, center);
	for (int i = 0; i < vertice_num; ++i) {
		SUM(cursor_vertices_pos[i], ori_vertices_pos[i], trans);
	}
	memcpy(lines[0].data(), u, 24);
	memcpy(lines[1].data(), cursor_pos, 24);
}

void Cursor::createVertices(double radius, double camera_center[3])
{
	memcpy(center, camera_center, 24);
	//R = radius;
	R = radius;
	int longitude_num_1 = longitude_num + 1;
	int latitude_num_1 = latitude_num - 1;
	for (int j = 1; j < latitude_num_1; ++j) {
		for (int i = 0; i < longitude_num_1; ++i) {
			ori_vertices_pos.push_back(std::array<double, 3>{ R * sin(M_PI * (double)j / (double)latitude_num_1) * cos(2.0*M_PI  * (double)i / (double)longitude_num) + center[0],
			R * sin(M_PI * (double)j / (double)latitude_num_1) * sin(2.0*M_PI * (double)i / (double)longitude_num) + center[1],
			R * cos(M_PI * (double)j / (double)latitude_num_1) + center[2] });
		}
	}
	ori_vertices_pos.push_back(std::array<double, 3>{center[0], center[1], R + center[2]});
	ori_vertices_pos.push_back(std::array<double, 3>{center[0], center[1], center[2] - R});
	for (int j = 1; j < latitude_num - 2; ++j)
	{
		for (int i = 0; i < longitude_num; ++i)
		{
			indices.push_back(longitude_num_1 * (j - 1) + i);
			indices.push_back(longitude_num_1 * j + i);
			indices.push_back(longitude_num_1 * (j - 1) + i + 1);
			indices.push_back(longitude_num_1 * j + i);
			indices.push_back(longitude_num_1 * j + i + 1);
			indices.push_back(longitude_num_1 * (j - 1) + i + 1);
		}
	}
	for (int i = 0; i < longitude_num; ++i) {
		indices.push_back(vertice_num - 2);
		indices.push_back(i);
		indices.push_back(i + 1);
		indices.push_back(longitude_num_1 * (latitude_num - 3) + i);
		indices.push_back(vertice_num - 1);
		indices.push_back(longitude_num_1 * (latitude_num - 3) + i + 1);
	}

	vertices_pos = ori_vertices_pos;
	cursor_vertices_pos = ori_vertices_pos;
	normal.resize(vertice_num);
	getNormal();
}


void Cursor::genBuffer()
{
	glGenVertexArrays(1, &VAO1);
	glGenBuffers(2, VBO1);
	glGenBuffers(1, &EBO1);
	glGenVertexArrays(1, &VAO2);
	glGenBuffers(2, VBO2);
	glGenBuffers(1, &EBO2);
	glGenVertexArrays(1, &VAO3);
	glGenBuffers(1, &VBO3);
}
void Cursor::setBufferData()
{
	glBindVertexArray(VAO1);
	glBindBuffer(GL_ARRAY_BUFFER, VBO1[0]);
	glBufferData(GL_ARRAY_BUFFER, vertices_pos.size() * sizeof(std::array<double, 3>), vertices_pos[0].data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO1);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(int), &indices[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, VBO1[1]);
	glBufferData(GL_ARRAY_BUFFER, normal.size() * sizeof(std::array<double, 3>), normal[0].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);
	glBindVertexArray(VAO2);
	glBindBuffer(GL_ARRAY_BUFFER, VBO2[0]);
	glBufferData(GL_ARRAY_BUFFER, cursor_vertices_pos.size() * sizeof(std::array<double, 3>), cursor_vertices_pos[0].data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO2);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(int), &indices[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, VBO2[1]);
	glBufferData(GL_ARRAY_BUFFER, normal.size() * sizeof(std::array<double, 3>), normal[0].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);
	glBindVertexArray(VAO3);
	glBindBuffer(GL_ARRAY_BUFFER, VBO3);
	glBufferData(GL_ARRAY_BUFFER, lines.size() * sizeof(std::array<double, 3>), &lines[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);

}


void Cursor::getNormal()
{
	double norm[3];
	for (int i = 0; i < vertice_num; ++i) {
		SUB(normal[i], vertices_pos[i], center);
		normalize(normal[i].data());
	}
	

	//double norm[3], e0[3], e2[3];
	//memset(normal[0].data(), 0, 24 * vertice_num);
	//int face_num;
	//face_num = indices.size() / 3;
	//int index;
	//for (int j = 0; j < face_num; ++j) {
	//	index = 3 * j;
	//	SUB(e2, vertices_pos[indices[index+2]], vertices_pos[indices[index+1]]);
	//	SUB(e0, vertices_pos[indices[index]], vertices_pos[indices[index+1]]);
	//	CROSS(norm, e2, e0);
	//	for (int k = 0; k < 3; ++k) {
	//		SUM(normal[indices[index + k]],
	//			normal[indices[index + k]], norm);
	//	}
	//}
	//for (int i = 0; i < vertice_num; ++i) {
	//	normalize(normal[i].data());
	//}

}