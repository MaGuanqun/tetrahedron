#include"pick_triangle.h"
PickTriangle::PickTriangle()
{
	initalFBO();

	shader = new Shader("./basic/picking.vs", "./basic/picking.fs");
	*shader = Shader("./basic/picking.vs", "./basic/picking.fs");
	collider_shader = new Shader("./shaders/wireframe.vs", "./shaders/wireframe.fs");
	*collider_shader = Shader("./shaders/wireframe.vs", "./shaders/wireframe.fs");
}

void PickTriangle::initalFBO()
{
	glGenFramebuffers(1, &picking_FBO);
	glBindFramebuffer(GL_FRAMEBUFFER, picking_FBO);

	glGenRenderbuffers(1, &picking_RBO_color);
	glBindRenderbuffer(GL_RENDERBUFFER, picking_RBO_color);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_RGB, SCR_WIDTH, SCR_HEIGHT);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);

	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, picking_RBO_color);

	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void PickTriangle::pickTriangle(std::vector<Cloth>* cloth, std::vector<Collider>* collider, Camera* camera, std::vector<std::vector<bool>>& hide, int* triangle_index, int* pos)
{
	glDisable(GL_BLEND);
	glDisable(GL_MULTISAMPLE);
	glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
	int total_face_num = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		total_face_num += (*cloth)[i].mesh_struct.faces.size();
	}
	int face_index;
	writingFBO(cloth, collider, camera, hide, shader);
	decideTriangle(face_index, total_face_num, &picking_FBO, pos);
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	decideFinalIndi(cloth, face_index, triangle_index);
}

void PickTriangle::writingFBO(std::vector<Cloth>* cloth, std::vector<Collider>* collider, Camera* camera, std::vector<std::vector<bool>>& hide, Shader* shader)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBindFramebuffer(GL_FRAMEBUFFER, picking_FBO);
	picking_base_ID = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		if (!hide[1][i]) {
			shader->use();
			shader->setInt("picking_base_ID", picking_base_ID);
			(*cloth)[i].simpDraw(camera, shader);
		}
		picking_base_ID = (*cloth)[i].mesh_struct.faces.size();
	}
	for (int i = 0; i < collider->size(); ++i) {
		if (!hide[0][i]) {
			collider_shader->use();
			collider_shader->setVec3("color", glm::vec3(1.0f, 1.0f, 1.0f));
			(*collider)[i].simpDraw(camera, collider_shader);
		}
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}


void PickTriangle::decideTriangle(int& triangle_index, int total_cloth_triangle_num, unsigned int* FBO, int* pos)
{
	std::vector<unsigned char> pixel_value(3);
	readPixel(&pixel_value, FBO, pos);
	int num_component[3];
	int triangle;
	num_component[0] = (int)(pixel_value[0]);
	num_component[1] = ((int)(pixel_value[1])) * 256;
	num_component[2] = ((int)(pixel_value[2])) * 256 * 256;
	triangle = num_component[0] + num_component[1] + num_component[2];
	if (triangle > -1 && triangle < total_cloth_triangle_num) {
		triangle_index=triangle;
	}
	else {
		triangle_index = -1;
	}

}

void PickTriangle::readPixel(std::vector<unsigned char>* pixel_value, unsigned int* FBO, int* pos)
{
	glBindFramebuffer(GL_FRAMEBUFFER, *FBO);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(pos[0], pos[1], 1, 1, GL_RGB, GL_UNSIGNED_BYTE, &(*pixel_value)[0]);
	glDeleteFramebuffers(1, FBO);
}

void PickTriangle::decideFinalIndi(std::vector<Cloth>* cloth, int sum_triangle_index, int* triangle_index)
{
	std::vector<int>start_indi(cloth->size() + 1);
	int id = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		start_indi[i] = id;
		id += (*cloth)[i].mesh_struct.faces.size();
	}
	start_indi[cloth->size()] = id;
	for (int j = 0; j < cloth->size(); ++j) {
		if (sum_triangle_index >= start_indi[j] && sum_triangle_index < start_indi[j + 1]) {
			triangle_index[1] = j;
			triangle_index[0] = sum_triangle_index - start_indi[j];
			break;
		}
	}
	if (sum_triangle_index == -1) {
		triangle_index[1] = -1;
		triangle_index[0] = -1;
	}
}