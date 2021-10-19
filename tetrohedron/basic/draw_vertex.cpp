#include"draw_vertex.h"

DrawVertex::DrawVertex()
{
	shader = new Shader("./shader/basic.vs", "./shader/basic.fs");
	*shader = Shader("./shader/basic.vs", "./shader/basic.fs");
	setBasicSphere();
	genBuffer();
	light.ambient = glm::vec3(0.3, 0.3, 0.3);
	light.diffuse = glm::vec3(0.7, 0.7, 0.7);
	light.specular = glm::vec3(0.95, 0.95, 0.95);

}

void DrawVertex::setVertex(std::vector<std::array<double, 3>>& v_list, std::vector<bool>& index)
{
	std::array<double, 3>v;
	draw_vertex.clear();
	draw_vertex.reserve(sphere_vertex_num * index.size());
	indices.clear();
	indices.reserve(index.size() * basic_index.size());
	draw_normal.clear();
	draw_normal.reserve(sphere_vertex_num * index.size());
	int num = 0;
	for (int i = 0; i < index.size(); i++) {
		if (index[i]) {
			for (int j = 0; j < sphere_vertex_num; j++) {
				SUM(v.data(), v_list[i].data(), sphere[j].data());
				draw_vertex.push_back(v);
			}
			draw_normal.insert(draw_normal.end(), sphere_normal.begin(), sphere_normal.end());
			int ind = num * sphere_vertex_num;
			for (int j = 0; j < basic_index.size(); ++j) {
				indices.push_back(basic_index[j] + ind);
			}
			num++;
		}
	}
}

void DrawVertex::initialVertex(int vertex_num)
{
	draw_vertex.clear();
	draw_vertex.reserve(sphere_vertex_num * vertex_num);
	draw_normal.clear();
	draw_normal.reserve(sphere_vertex_num * vertex_num);
	indices.clear();
	indices.reserve(vertex_num * basic_index.size());
}
void DrawVertex::setVertexAccumulate(std::vector<std::array<double, 3>>& v_list, std::vector<int>& index)
{

	std::array<double, 3>v;
	int index_size_already = indices.size();
	for (int i = 0; i < index.size(); i++) {
		for (int j = 0; j < sphere_vertex_num; j++) {
			SUM(v.data(), v_list[index[i]].data(), sphere[j].data());
			draw_vertex.push_back(v);
		}
		draw_normal.insert(draw_normal.end(), sphere_normal.begin(), sphere_normal.end());
		int ind = i * sphere_vertex_num;
		for (int j = 0; j < basic_index.size(); ++j) {
			indices.push_back(basic_index[j] + ind + index_size_already);
		}
	}
}


void DrawVertex::setVertex(std::vector<std::array<double, 3>>& v_list, std::vector<int>& index)
{
	std::array<double, 3>v;
	draw_vertex.clear();
	draw_vertex.reserve(sphere_vertex_num * index.size());
	draw_normal.clear();
	draw_normal.reserve(sphere_vertex_num * index.size());
	indices.clear();
	indices.reserve(index.size() * basic_index.size());
	for (int i = 0; i < index.size(); i++) {
		for (int j = 0; j < sphere_vertex_num; j++) {
			SUM(v.data(), v_list[index[i]].data(), sphere[j].data());
			draw_vertex.push_back(v);
		}
		draw_normal.insert(draw_normal.end(), sphere_normal.begin(), sphere_normal.end());
		int ind = i * sphere_vertex_num;
		for (int j = 0; j < basic_index.size(); ++j) {
			indices.push_back(basic_index[j] + ind);
		}
	}
}


void DrawVertex::setVertex(double* vertex, double tolerance)
{
	std::array<double, 3>v;
	draw_vertex.clear();
	draw_vertex.reserve(sphere_vertex_num);
	draw_normal.clear();
	draw_normal.reserve(sphere_vertex_num);
	indices.clear();
	indices.reserve(basic_index.size());
	setSphereRadius(tolerance);


	for (int j = 0; j < sphere_vertex_num; j++) {
		SUM(v, vertex, sphere[j].data());
		draw_vertex.push_back(v);
	}
	draw_normal.insert(draw_normal.end(), sphere_normal.begin(), sphere_normal.end());
	for (int j = 0; j < basic_index.size(); ++j) {
		indices.push_back(basic_index[j]);
	}
}


void DrawVertex::setVertex(std::vector<std::array<double, 3>>& v_list, double tolerance)
{
	std::array<double, 3>v;
	draw_vertex.clear();
	draw_vertex.reserve(sphere_vertex_num * v_list.size());
	draw_normal.clear();
	draw_normal.reserve(sphere_vertex_num * v_list.size());
	indices.clear();
	indices.reserve(v_list.size() * basic_index.size());
	int num = 0;
	setSphereRadius(tolerance);
	for (int i = 0; i < v_list.size(); i++) {		
		for (int j = 0; j < sphere_vertex_num; j++) {
			SUM(v.data(), v_list[i].data(), sphere[j].data());
			draw_vertex.push_back(v);
		}
		draw_normal.insert(draw_normal.end(), sphere_normal.begin(), sphere_normal.end());
		int ind = num * sphere_vertex_num;
		for (int j = 0; j < basic_index.size(); ++j) {
			indices.push_back(basic_index[j] + ind);
		}
		num++;
		
	}
}

void DrawVertex::setSphereRadius(double R)
{
	int longitude_num_1 = longitude_num + 1;
	int latitude_num_1 = latitude_num - 1;
	sphere.clear();
	for (int j = 1; j < latitude_num_1; ++j) {
		for (int i = 0; i < longitude_num_1; ++i) {
			sphere.push_back(std::array{ R * sin(M_PI * (double)j / (double)latitude_num_1) * cos(2.0 * M_PI * (double)i / (double)longitude_num),
			R * sin(M_PI * (double)j / (double)latitude_num_1) * sin(2.0 * M_PI * (double)i / (double)longitude_num),
			R * cos(M_PI * (double)j / (double)latitude_num_1) });
		}
	}
	sphere.push_back(std::array{ 0.0,0.0,R });
	sphere.push_back(std::array{ 0.0,0.0,-R });
}


void DrawVertex::setBasicSphere()
{
	int longitude_num_1 = longitude_num + 1;
	int latitude_num_1 = latitude_num - 1;
	double R = 0.005;
	sphere.reserve(longitude_num_1 * latitude_num_1);
	basic_index.reserve(6*longitude_num * latitude_num);
	for (int j = 1; j < latitude_num_1; ++j) {
		for (int i = 0; i < longitude_num_1; ++i) {
			sphere.push_back(std::array{ R * sin(M_PI * (double)j / (double)latitude_num_1) * cos(2.0 * M_PI * (double)i / (double)longitude_num),
			R * sin(M_PI * (double)j / (double)latitude_num_1) * sin(2.0 * M_PI * (double)i / (double)longitude_num),
			R * cos(M_PI * (double)j / (double)latitude_num_1) });
		}
	}
	sphere.push_back(std::array{ 0.0,0.0,R });
	sphere.push_back(std::array{ 0.0,0.0,-R });
	sphere_vertex_num = sphere.size();
	for (int j = 1; j < latitude_num - 2; ++j)
	{
		for (int i = 0; i < longitude_num; ++i)
		{
			basic_index.push_back(longitude_num_1 * (j - 1) + i);
			basic_index.push_back(longitude_num_1 * j + i);
			basic_index.push_back(longitude_num_1 * (j - 1) + i + 1);
			basic_index.push_back(longitude_num_1 * j + i);
			basic_index.push_back(longitude_num_1 * j + i + 1);
			basic_index.push_back(longitude_num_1 * (j - 1) + i + 1);
		}
	}
	for (int i = 0; i < longitude_num; ++i) {
		basic_index.push_back(sphere_vertex_num - 2);
		basic_index.push_back(i);
		basic_index.push_back(i + 1);
		basic_index.push_back(longitude_num_1 * (latitude_num - 3) + i);
		basic_index.push_back(sphere_vertex_num - 1);
		basic_index.push_back(longitude_num_1 * (latitude_num - 3) + i + 1);
	}
	sphere_normal=sphere;
	for (int i = 0; i < sphere.size(); ++i) {
		normalize(sphere_normal[i].data());
	}
}

void DrawVertex::genBuffer()
{
	glGenVertexArrays(1, &VAO);
	glGenBuffers(2, VBO);
	glGenBuffers(1, &EBO);
}

void DrawVertex::setBuffer()
{
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
	glBufferData(GL_ARRAY_BUFFER, draw_vertex.size() * sizeof(std::array<double, 3>), draw_vertex[0].data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(int), &indices[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
	glBufferData(GL_ARRAY_BUFFER, draw_normal.size() * sizeof(std::array<double, 3>), draw_normal[0].data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
	glBindVertexArray(0);
}
void DrawVertex::draw(Camera* camera,glm::vec3& color,float transparence)
{
	if (!indices.empty()) {
		setBuffer();
		shader->use();
		shader->setMat4("projection", camera->GetProjectMatrix());
		shader->setMat4("view", camera->GetViewMatrix());
		shader->setMat4("model", glm::mat4(1.0));
		shader->setVec3("color", color);
		shader->setVec3("lightPos", camera->position);
		shader->setVec3("light.ambient", light.ambient);
		shader->setVec3("light.diffuse", light.diffuse);
		shader->setVec3("light.specular", light.specular);
		shader->setVec3("viewPos", camera->position);
		shader->setFloat("transparence", transparence);
		glBindVertexArray(VAO);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
}