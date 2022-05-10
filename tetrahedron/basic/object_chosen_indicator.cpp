#include"object_chosen_indicator.h"

ObjectChosenIndicator::ObjectChosenIndicator()
{
	initalFBO();
	vertex_num = 50;
	circle_vertices.resize(3);
	for (unsigned int i = 0; i < circle_vertices.size(); ++i) {
		circle_vertices[i].resize(3 * vertex_num);
	}
	genBuffer();
	circle_color.resize(3);
	circle_color[0] = glm::vec3(1.0, 0.1, 0.1);
	circle_color[1] = glm::vec3(0.1, 1.0, 0.1);
	circle_color[2] = glm::vec3(0.1, 0.1, 1.0);

}


void ObjectChosenIndicator::updatePosition(double* AABB)
{
	double center[3] = { 0.5 * (AABB[0] + AABB[3]),0.5 * (AABB[1] + AABB[4]),0.5 * (AABB[2] + AABB[5]) };
	double radius = 0.5 * sqrt((AABB[0] - AABB[3]) * (AABB[0] - AABB[3]) + (AABB[1] - AABB[4]) * (AABB[1] - AABB[4])
		+ (AABB[2] - AABB[5]) * (AABB[2] - AABB[5]));

	for (int i = 0; i < vertex_num; i++)
	{
		circle_vertices[0][3 * i] = center[0]+ radius * sin(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));
		circle_vertices[0][3 * i + 1] = center[1] + radius * cos(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));
		circle_vertices[0][3 * i + 2] = center[2];

		circle_vertices[1][3 * i] = center[0];
		circle_vertices[1][3 * i + 1] = center[1] + radius * sin(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));
		circle_vertices[1][3 * i + 2] = center[2] + radius * cos(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));

		circle_vertices[2][3 * i] = center[0] + radius * cos(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));
		circle_vertices[2][3 * i + 1] = center[1];
		circle_vertices[2][3 * i + 2] = center[2] + radius * sin(DEG_RADIANS(360.0 * (double)i / (double)vertex_num));
	}	
	setBuffer();
}

void ObjectChosenIndicator::setBuffer()
{
	for (unsigned int i = 0; i < 3; ++i) {
		glBindVertexArray(VAO[i]);
		glBindBuffer(GL_ARRAY_BUFFER, VBO[i]);
		glBufferData(GL_ARRAY_BUFFER, circle_vertices[i].size()*sizeof(double), circle_vertices[i].data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}
}


void ObjectChosenIndicator::initalFBO()
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


bool ObjectChosenIndicator::pickAxes(Shader* shader, Camera* camera, unsigned int& dimension, int* pos)
{
	glDisable(GL_BLEND);
	glDisable(GL_MULTISAMPLE);
	glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);

	writingFBO(camera, shader);
	decideAxe(dimension, &picking_FBO, pos);
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	if (dimension == 4) {
		return false;
	}
	return true;
}


void ObjectChosenIndicator::decideAxe(unsigned int& dimension, unsigned int* FBO, int* pos)
{
	dimension = 4;
	std::vector<unsigned char> pixel_value(3);
	readPixel(&pixel_value, FBO, pos);
	for (unsigned int i = 0; i < 3; ++i) {
		if ((int)(pixel_value[i]) > 250) {
			dimension = i;
			return;
		}
	}
}

void ObjectChosenIndicator::readPixel(std::vector<unsigned char>* pixel_value, unsigned int* FBO, int* pos)
{
	glBindFramebuffer(GL_FRAMEBUFFER, *FBO);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(pos[0], pos[1], 1, 1, GL_RGB, GL_UNSIGNED_BYTE, &(*pixel_value)[0]);
	glDeleteFramebuffers(1, FBO);
}


void ObjectChosenIndicator::writingFBO(Camera* camera, Shader* shader)
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBindFramebuffer(GL_FRAMEBUFFER, picking_FBO);	
	draw(shader, camera,4,12.0f);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}


void ObjectChosenIndicator::draw(Shader* shader, Camera* camera, unsigned int dimension, float line_width)
{
	shader->use();
	shader->setMat4("projection", camera->GetProjectMatrix());
	shader->setMat4("view", camera->GetViewMatrix());
	shader->setMat4("model", glm::mat4(1.0));
	shader->setFloat("transparent", 1.0f);

	//std::cout <<"kk "<< dimension << " " << line_width << std::endl;

	for (unsigned int i = 0; i < 3; i++)
	{
		shader->setVec3("color", circle_color[i]);
		glBindVertexArray(VAO[i]);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		if (dimension == i) {
			glLineWidth(2.0f*line_width);
		}
		else {
			glLineWidth(line_width);
		}		
		glDrawArrays(GL_LINE_LOOP, 0, vertex_num);
		glBindVertexArray(0);
	}
	
	glLineWidth(1.0f);
}

void ObjectChosenIndicator::genBuffer()
{
	for (unsigned int i = 0; i < 3; ++i) {
		glGenVertexArrays(1, &VAO[i]);
		glGenBuffers(1, &VBO[i]);
	}
	
}