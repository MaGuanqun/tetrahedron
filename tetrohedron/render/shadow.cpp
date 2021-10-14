#include"shadow.h"

Shadow::Shadow()
{
	shadow_shader = new Shader("./shader/pointShadow.vs", "./shader/pointShadow.fs", "./shader/pointShadow.gs");
	*shadow_shader = Shader("./shader/pointShadow.vs", "./shader/pointShadow.fs", "./shader/pointShadow.gs");
}



void Shadow::setBasic()
{
	light_projection_matrix = glm::perspective(glm::radians(90.0f), (float)SHADOW_WIDTH / (float)SHADOW_HEIGHT, near_plane, far_plane);
	for (int i = 0; i < 6; ++i)
	{
		shadow_view.push_back(glm::mat4(0.0f));
		shadow_transforms.push_back(glm::mat4(0.0f));
	}
	FBOdepth();
}

void Shadow::FBOdepth()
{
	glGenFramebuffers(1, &depth_map_FBO);
	glGenTextures(1, &depth_map);


	glBindTexture(GL_TEXTURE_CUBE_MAP, depth_map);
	for (int i = 0; i < 6; ++i)
	{
		glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_DEPTH_COMPONENT, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
	}
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	glBindFramebuffer(GL_FRAMEBUFFER, depth_map_FBO);
	glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, depth_map, 0);

	// attach depth texture as FBO's depth buffer

	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void Shadow::drawShadow(Camera* camera, std::vector<std::vector<bool>>& hide, std::vector<Cloth>& cloth, std::vector<Collider>& collider, 
	std::vector<Tetrohedron>& tetrohedron, std::vector<int>& cloth_index_in_object, std::vector<int>& tetrohedron_index_in_object)
{
	lightSpace(camera->position);
	glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
	glBindFramebuffer(GL_FRAMEBUFFER, depth_map_FBO);
	glClear(GL_DEPTH_BUFFER_BIT);
	setIn(camera);
	shadow_shader->use();
	for (int i = 0; i < cloth.size(); ++i) {
		if (!hide[1][cloth_index_in_object[i]]) {
			cloth[i].drawShadow(camera, shadow_shader);
		}
	}
	for (int i = 0; i < collider.size(); ++i) {
		if (!hide[0][i]) {
			collider[i].drawShadow(camera, shadow_shader);
		}
	}
	for (int i = 0; i < tetrohedron.size(); ++i) {
		if (!hide[1][tetrohedron_index_in_object[i]]) {
			tetrohedron[i].drawShadow(camera, shadow_shader);
		}
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void Shadow::lightSpace(glm::vec3& light_pos)
{

	shadow_view[0] = glm::lookAt(light_pos, light_pos + camera_from_origin * glm::vec3(1.0f, 0.0f, 0.0f), camera_from_origin * glm::vec3(0.0f, -1.0f, 0.0f));
	shadow_view[1] = glm::lookAt(light_pos, light_pos + camera_from_origin * glm::vec3(-1.0f, 0.0f, 0.0f), camera_from_origin * glm::vec3(0.0f, -1.0f, 0.0f));
	shadow_view[2] = glm::lookAt(light_pos, light_pos + camera_from_origin * glm::vec3(0.0f, 1.0f, 0.0f), camera_from_origin * glm::vec3(0.0f, 0.0f, 1.0f));
	shadow_view[3] = glm::lookAt(light_pos, light_pos + camera_from_origin * glm::vec3(0.0f, -1.0f, 0.0f), camera_from_origin * glm::vec3(0.0f, 0.0f, -1.0f));
	shadow_view[4] = glm::lookAt(light_pos, light_pos + camera_from_origin * glm::vec3(0.0f, 0.0f, 1.0f), camera_from_origin * glm::vec3(0.0f, -1.0f, 0.0f));
	shadow_view[5] = glm::lookAt(light_pos, light_pos + camera_from_origin * glm::vec3(0.0f, 0.0f, -1.0f), camera_from_origin * glm::vec3(0.0f, -1.0f, 0.0f));
	for (int i = 0; i < 6; ++i)
	{
		shadow_transforms[i] = light_projection_matrix * shadow_view[i];
	}
}

void Shadow::setIn(Camera* camera)
{
	shadow_shader->use();
	for (int i = 0; i < 6; ++i)
	{
		shadow_shader->setMat4("shadowMatrices[" + std::to_string(i) + "]", shadow_transforms[i]);
	}
	shadow_shader->setFloat("far_plane", far_plane);
	shadow_shader->setVec3("lightPos", camera->position);
	shadow_shader->setMat4("model", glm::mat4(1.0));

}