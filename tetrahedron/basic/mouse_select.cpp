#include"mouse_select.h"
MouseSelect::MouseSelect()
{
	camera_pos = glm::vec3(1.0, 0.0, 0.0);
	shader = new Shader("./shader/wireframe.vs", "./shader/wireframe.fs");
	*shader = Shader("./shader/wireframe.vs", "./shader/wireframe.fs");

    view= glm::lookAt(camera_pos, glm::vec3(0.0,0.0,0.0), glm::vec3(0.0,1.0,0.0));
    projection= glm::ortho(-(float)SCR_WIDTH/(float)SCR_HEIGHT, (float)SCR_WIDTH / (float)SCR_HEIGHT, -1.0f, 1.0f, 0.01f, 3.0f);
    position.resize(12, -1.0f);
    for (int i = 0; i < 4; ++i) {
        position[3 * i] = 0.0;
    }
    genBuffer(); 
}

void MouseSelect::setFirstCorner(float x, float y)
{
    transferTo3D(x, y, first_corner);
}

void MouseSelect::setLastCorner(float x, float y)
{
    transferTo3D(x, y, last_corner);
    setPosition();
}

void MouseSelect::transferTo3D(float x, float y, float* point)
{
    float height_coe = 2.0f * (0.5f - y / (float)SCR_HEIGHT);
    double width_coe = 2.0f * (0.5f-x / (float)SCR_WIDTH) * ((float)SCR_WIDTH / (float)SCR_HEIGHT);
    point[0] = 0.0;
    point[1] = height_coe;
    point[2] = width_coe;
}

void MouseSelect::setPosition()
{
    float y[2], z[2];
    y[0] = myMin(first_corner[1], last_corner[1]);
    y[1] = myMax(first_corner[1], last_corner[1]);
    z[0] = myMin(first_corner[2], last_corner[2]);
    z[1] = myMax(first_corner[2], last_corner[2]);
    for (int i = 0; i < 2; ++i) {
        position[3 * i + 1] = y[1 - i];
        position[7 + 3 * i] = y[i];
        position[6 * i + 2] = z[i];
        position[6 * i + 5] = z[i];
    }

}

void MouseSelect::genBuffer()
{
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    setBufferData();
}


void MouseSelect::setBufferData()
{
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, position.size() * sizeof(float), &position[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindVertexArray(0);
}
void MouseSelect::draw()
//we use glViewport() to change the position, pay attention to it in main function
{
    setBufferData();
    glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
    glClear(GL_DEPTH_BUFFER_BIT);
    shader->use();
    shader->setMat4("projection", projection);
    shader->setMat4("view", view);
    shader->setMat4("model", glm::mat4(1.0));
    shader->setVec3("color", glm::vec3(0.0f,1.0f,0.0f));
    glBindVertexArray(VAO);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDrawArrays(GL_LINE_LOOP, 0, 4);
    glBindVertexArray(0);
   
}
