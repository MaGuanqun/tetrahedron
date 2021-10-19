#include"set_tetrohedron_anchor.h"
#include"global.h"

SetTetrohedronAnchor::SetTetrohedronAnchor()
{
	position.resize(4);
	shader = new Shader("./shader/wireframe.vs", "./shader/wireframe.fs");
	*shader = Shader("./shader/wireframe.vs", "./shader/wireframe.fs");
	camera_pos = glm::vec3(1.0, 0.0, 0.0);
	view = glm::lookAt(camera_pos, glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 1.0, 0.0));
	projection = glm::ortho(-(float)SCR_WIDTH / (float)SCR_HEIGHT, (float)SCR_WIDTH / (float)SCR_HEIGHT, -1.0f, 1.0f, 0.01f, 3.0f);
	genBuffer();
	vertex_color = glm::vec3(1.0, 0.0, 0.0);
}

void SetTetrohedronAnchor::setCorner(double* screen_pos, bool pre_press_state, std::vector<Tetrohedron>& tetrohedron, Camera* camera, 
	std::vector<bool>& hide)
{
	if (!pre_press_state) {
		transferTo3D(screen_pos, draw_first_corner);
		first_corner_in_clip_space[0] = 2.0f * (screen_pos[0] / (float)SCR_WIDTH -0.5 );
		first_corner_in_clip_space[1] = 2.0f * (0.5f - screen_pos[1] / (float)SCR_HEIGHT);
	}
	else {

		setClipSpaceRange(screen_pos);
		transferTo3D(screen_pos, draw_last_corner);
		setPosition();
		for (int i = 0; i < tetrohedron.size(); ++i) {
			if (!hide[i]) {
				findAllSelectVertex(tetrohedron[i].mesh_struct.vertex_for_render, camera, tetrohedron[i].mesh_struct.anchor_vertex);			
			}
		}	
		int all_select_num = 0;
		for (int i = 0; i < tetrohedron.size(); ++i) {
			if (!hide[i]) {
				all_select_num += tetrohedron[i].mesh_struct.anchor_vertex.size();
			}
		}
		draw_vertex.initialVertex(all_select_num / 5);
		for (int i = 0; i < tetrohedron.size(); ++i) {
			if (!hide[i]) {

				findSurfaceVertex(tetrohedron[i].mesh_struct.vertex_for_render, camera, tetrohedron[i].mesh_struct.anchor_vertex, tetrohedron[i].mesh_struct.vertex_on_surface);
			}
		}
		draw_vertex.draw(camera, vertex_color);
	}	
}

void SetTetrohedronAnchor::setClipSpaceRange(double* screen_pos)
{
	float second_corner[2];
	second_corner[0] = 2.0f * (screen_pos[0] / (float)SCR_WIDTH - 0.5);
	second_corner[1] = 2.0f * (0.5f - screen_pos[1] / (float)SCR_HEIGHT);
	max_corner[0] = myMax(first_corner_in_clip_space[0], second_corner[0]);
	max_corner[1] = myMax(first_corner_in_clip_space[1], second_corner[1]);
	min_corner[0] = myMin(first_corner_in_clip_space[0], second_corner[0]);
	min_corner[1] = myMin(first_corner_in_clip_space[1], second_corner[1]);
}

void SetTetrohedronAnchor::transferTo3D(double* screen_pos, float* point)
{
	float height_coe = 2.0f * (0.5f - screen_pos[1] / (float)SCR_HEIGHT);
	double width_coe = 2.0f * (0.5f - screen_pos[0] / (float)SCR_WIDTH) * ((float)SCR_WIDTH / (float)SCR_HEIGHT);
	point[0] = 0.0;
	point[1] = height_coe;
	point[2] = width_coe;
}

void SetTetrohedronAnchor::setPosition()
{
	float y[2], z[2];
	y[0] = myMin(draw_first_corner[1], draw_last_corner[1]);
	y[1] = myMax(draw_first_corner[1], draw_last_corner[1]);
	z[0] = myMin(draw_first_corner[2], draw_last_corner[2]);
	z[1] = myMax(draw_first_corner[2], draw_last_corner[2]);
	for (int i = 0; i < 2; ++i) {
		position[i].data()[0] = 0.0;
		position[i + 2].data()[0] = 0.0;
		position[i].data()[1] = y[1 - i];
		position[i + 2][1] = y[i];
		position[2 * i][2] = z[i];
		position[2 * i + 1][2] = z[i];
	}	


}

void SetTetrohedronAnchor::genBuffer()
{
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	setBufferData();
}


void SetTetrohedronAnchor::setBufferData()
{
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, position.size() * sizeof(std::array<float,3>), &position[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glBindVertexArray(0);
}

void SetTetrohedronAnchor::draw()
//we use glViewport() to change the position, pay attention to it in main function
{
	setBufferData();
	glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
	glClear(GL_DEPTH_BUFFER_BIT);
	shader->use();
	shader->setMat4("projection", projection);
	shader->setMat4("view", view);
	shader->setMat4("model", glm::mat4(1.0));
	shader->setVec3("color", glm::vec3(0.0f, 1.0f, 0.0f));
	glBindVertexArray(VAO);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawArrays(GL_LINE_LOOP, 0, 4);
	glBindVertexArray(0);

}


void SetTetrohedronAnchor::findAllSelectVertex(std::vector<std::array<double,3>>& vertex_for_render, Camera* camera, std::vector<int>& vertex_index)
{
	double vertex_pos_view[3];
	glm::vec4 vertex_pos;
	glm::vec4 pos_after_transfer;
	vertex_pos[3] = 1.0;
	glm::mat4 transfer = camera->GetProjectMatrix() * camera->GetViewMatrix();
	vertex_index.clear();
	for (int i = 0; i < vertex_for_render.size(); ++i) {
		vertex_pos[0]= vertex_for_render[i][0];
		vertex_pos[1]= vertex_for_render[i][1];
		vertex_pos[2]= vertex_for_render[i][2];
		pos_after_transfer = transfer * vertex_pos;
		pos_after_transfer[0] /= pos_after_transfer[3];
		pos_after_transfer[1] /= pos_after_transfer[3];

		if (pos_after_transfer[0] > min_corner[0] && pos_after_transfer[0]<max_corner[0]
			&& pos_after_transfer[1]>min_corner[1] && pos_after_transfer[1] < max_corner[1]) {
			vertex_index.push_back(i);
		}
	}	
}

void SetTetrohedronAnchor::findSurfaceVertex(std::vector<std::array<double, 3>>& vertex_for_render, Camera* camera, std::vector<int>& vertex_index, std::vector<bool>& is_surface)
{
	std::vector<int> surface_index;
	surface_index.reserve(vertex_index.size() / 5);
	for (int i = 0; i < vertex_index.size(); ++i) {
		if (is_surface[vertex_index[i]]) {
			surface_index.push_back(vertex_index[i]);
		}
	}
	draw_vertex.setVertexAccumulate(vertex_for_render, surface_index);
}