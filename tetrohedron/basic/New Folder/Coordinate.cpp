#include"Coordinate.h"

CoordinateSystem::CoordinateSystem()
{
    Arrow arrowX(0);
    vertex_indices0 = arrowX.vertex_indices;
    vertices0 = arrowX.vertices;
    normals0 = arrowX.normals;
    Arrow arrowY(1);
    vertex_indices1 = arrowY.vertex_indices;
    vertices1 = arrowY.vertices;
    normals1 = arrowY.normals;
    Arrow arrowZ(2);
    vertex_indices2 = arrowZ.vertex_indices;
    vertices2 = arrowZ.vertices;
    normals2 = arrowZ.normals;   
    shader=new Shader("./shader/basic.vs", "./shader/basic.fs");
    *shader = Shader("./shader/basic.vs", "./shader/basic.fs");
    genBuffer();
    light.ambient = glm::vec3(0.4, 0.4, 0.4);
    light.diffuse=glm::vec3(0.8, 0.8, 0.8);
    light.specular=glm::vec3(0.95, 0.95, 0.95);
}

void CoordinateSystem::draw(Camera* camera, glm::vec3 lightPos)
//we use glViewport() to change the position, pay attention to it in main function
{    
    glViewport(-SCR_WIDTH/3, -SCR_HEIGHT/5*2, SCR_WIDTH, SCR_HEIGHT);
    glClear(GL_DEPTH_BUFFER_BIT);
    shader->use();
    shader->setMat4("projection", camera->GetProjectMatrix());
    shader->setMat4("view", glm::lookAt(2.0f * glm::normalize(camera->position-camera->center), glm::vec3(0.0f, 0.0f, 0.0f), camera->up));
    shader->setMat4("model", glm::mat4(1.0));
    shader->setVec3("viewPos", 2.0f*glm::normalize(camera->position-camera->center));
    shader->setVec3("lightPos", lightPos);
    shader->setFloat("transparence", 1.0);
    shader->setVec3("light.specular", light.specular);
    shader->setVec3("light.diffuse", light.diffuse);
    shader->setVec3("light.ambient", light.ambient);
    shader->setVec3("color", glm::vec3(1.0, 0.1, 0.1));
    glBindVertexArray(VAO0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDrawElements(GL_TRIANGLES, vertex_indices0.size(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);    
    //shader->setVec3("color", glm::vec3(0.1, 1.0, 0.1));
    //glBindVertexArray(VAO1);
    //glDrawElements(GL_TRIANGLES, vertex_indices1.size(), GL_UNSIGNED_INT, 0);
    //glBindVertexArray(0);
    //shader->setVec3("color", glm::vec3(0.1, 1.0, 0.1));
    //glBindVertexArray(VAO2);
    //glDrawElements(GL_TRIANGLES, vertex_indices2.size(), GL_UNSIGNED_INT, 0);
    //glBindVertexArray(0);
}

void CoordinateSystem::genBuffer()
{   
    glGenVertexArrays(1, &VAO0);
    glGenBuffers(3, VBO0);
    glGenBuffers(1, &EBO0);
    glGenVertexArrays(1, &VAO1);
    glGenBuffers(3, VBO1);
    glGenBuffers(1, &EBO1);
    glGenVertexArrays(1, &VAO2);
    glGenBuffers(3, VBO2);
    glGenBuffers(1, &EBO2);
    setBufferData();
}

void CoordinateSystem::setBufferData()
{
    glBindVertexArray(VAO0);
    glBindBuffer(GL_ARRAY_BUFFER, VBO0[0]);
    glBufferData(GL_ARRAY_BUFFER, vertices0.size() * sizeof(float), &vertices0[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO0);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, vertex_indices0.size() * sizeof(int), &vertex_indices0[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, VBO0[1]);
    glBufferData(GL_ARRAY_BUFFER, normals0.size() * sizeof(float), &normals0[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindVertexArray(0);

    glBindVertexArray(VAO1);
    glBindBuffer(GL_ARRAY_BUFFER, VBO1[0]);
    glBufferData(GL_ARRAY_BUFFER, vertices1.size() * sizeof(float), &vertices1[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO1);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, vertex_indices1.size() * sizeof(int), &vertex_indices1[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, VBO1[1]);
    glBufferData(GL_ARRAY_BUFFER, normals1.size() * sizeof(float), &normals1[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindVertexArray(0);

    glBindVertexArray(VAO2);
    glBindBuffer(GL_ARRAY_BUFFER, VBO2[0]);
    glBufferData(GL_ARRAY_BUFFER, vertices2.size() * sizeof(float), &vertices2[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO2);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, vertex_indices2.size() * sizeof(int), &vertex_indices2[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, VBO2[1]);
    glBufferData(GL_ARRAY_BUFFER, normals2.size() * sizeof(float), &normals2[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindVertexArray(0);

}

Arrow::Arrow(int type)//type=0, x, type=1, y, type=2, z
{
    setArrow(type);
    setCylinder(type);
    normals = getNorm(vertices, vertex_indices);
}

void Arrow::setArrow(int type)
{
    int i, j;
    int angle_num = 16;//This is the number of vertices of the cone bottom
    float vertice[3] = { 0.0f, 0.0f, 0.0f };
    int arrow_vertex_indices[6] = { 0,0,0,1,1,1 };
    vertice[type] =  0.125f;        
    for (i = 0; i < 3; ++i)
    {
        vertices.push_back(vertice[i]);
    }
    vertice[type] = 0.1f;
    for (i = 0; i < 3; ++i)
    {
        vertices.push_back(vertice[i]);
    }    
    if (type == 2) {
        for (i = 0; i < angle_num; ++i)
        {
            vertice[0] = 0.0125 * cos(2.0 * M_PI / (float)angle_num * i);
            vertice[1] = 0.0125 * sin(2.0 * M_PI / (float)angle_num * i);
            for (j = 0; j < 3; ++j)
            {
                vertices.push_back(vertice[j]);
            }
        }
    }
    else if (type == 0) {
        for (i = 0; i < angle_num; ++i)
        {
            vertice[1] = 0.0125 * cos(2.0 * M_PI / (float)angle_num * i);
            vertice[2] = 0.0125 * sin(2.0 * M_PI / (float)angle_num * i);
            for (j = 0; j < 3; ++j)
            {
                vertices.push_back(vertice[j]);
            }
        }
    }
    else {
        for (i = 0; i < angle_num; ++i)
        {
            vertice[0] = 0.0125 * cos(2.0 * M_PI / (float)angle_num * i);
            vertice[2] = -0.0125 * sin(2.0 * M_PI / (float)angle_num * i);
            for (j = 0; j < 3; ++j)
            {
                vertices.push_back(vertice[j]);
            }
        }
    }
    for (i = 0; i < angle_num - 1; ++i)
    {
        arrow_vertex_indices[1] = 2 + i;
        arrow_vertex_indices[2] = 3 + i;
        arrow_vertex_indices[4] = 3 + i;
        arrow_vertex_indices[5] = 2 + i;
        for (j = 0; j < 6; ++j)
        {
            vertex_indices.push_back(arrow_vertex_indices[j]);
        }
    }

    arrow_vertex_indices[1] = angle_num + 1;
    arrow_vertex_indices[2] = 2;
    arrow_vertex_indices[4] = angle_num + 1;
    arrow_vertex_indices[5] = 2;
    for (j = 0; j < 6; ++j)
    {
        vertex_indices.push_back(arrow_vertex_indices[j]);
    }    
}

void Arrow::setCylinder(int type) {
    int ang_num = 16;
    int basIndex = vertices.size()/3;
    float vertice[3] = { 0.0f, 0.0f, 0.0f};
    vertice[type] = 0.1f;
    for (int i = 0; i < 3; ++i)
    {
        vertices.push_back(vertice[i]);
    }
    vertice[type] = 0.0f;
    for (int i = 0; i < 3; ++i)
    {
        vertices.push_back(vertice[i]);
    }
    if (type == 2) {
        for (int i = 0; i < ang_num; ++i)
        {
            vertice[0] = 0.005 * cos(2.0 * M_PI / (float)ang_num * i);
            vertice[1] = 0.005 * sin(2.0 * M_PI / (float)ang_num * i);
            vertice[2] = 0.1f;
            for (int j = 0; j < 3; ++j)
            {
                vertices.push_back(vertice[j]);
            }
            vertice[2] = 0.0f;
            for (int j = 0; j < 3; ++j)
            {
                vertices.push_back(vertice[j]);
            }
        }
    }
    else if(type==0)
    {
        for (int i = 0; i < ang_num; ++i)
        {
            vertice[1] = 0.005 * cos(2.0 * M_PI / (float)ang_num * i);
            vertice[2] = 0.005 * sin(2.0 * M_PI / (float)ang_num * i);
            vertice[0] = 0.1f;
            for (int j = 0; j < 3; ++j)
            {
                vertices.push_back(vertice[j]);
            }
            vertice[0] = 0.0f;
            for (int j = 0; j < 3; ++j)
            {
                vertices.push_back(vertice[j]);
            }
        }
    }
    else {
        for (int i = 0; i < ang_num; ++i)
        {
            vertice[0] = 0.005 * cos(2.0 * M_PI / (float)ang_num * i);
            vertice[2] = -0.005 * sin(2.0 * M_PI / (float)ang_num * i);
            vertice[1] = 0.1f;
            for (int j = 0; j < 3; ++j)
            {
                vertices.push_back(vertice[j]);
            }
            vertice[1] = 0.0f;
            for (int j = 0; j < 3; ++j)
            {
                vertices.push_back(vertice[j]);
            }
        }
    }
    int cylinder_vertex_indices[6] = { basIndex,basIndex,basIndex,basIndex+1,basIndex+1,basIndex+1 };
    for (int i = 0; i < ang_num - 1; ++i)
    {
        cylinder_vertex_indices[1] = basIndex + 2 + 2 * i;
        cylinder_vertex_indices[2] = basIndex + 3 + 2 * i;
        cylinder_vertex_indices[3] = basIndex + 2 + 2 * (i+1);
        cylinder_vertex_indices[4] = basIndex + 3 + 2 * (i+1);       
        //   0
        // 1   3
        // 2   4
        //   5
        vertex_indices.push_back(cylinder_vertex_indices[0]);
        vertex_indices.push_back(cylinder_vertex_indices[1]);
        vertex_indices.push_back(cylinder_vertex_indices[3]);
        vertex_indices.push_back(cylinder_vertex_indices[1]);
        vertex_indices.push_back(cylinder_vertex_indices[2]);
        vertex_indices.push_back(cylinder_vertex_indices[3]);
        vertex_indices.push_back(cylinder_vertex_indices[3]);
        vertex_indices.push_back(cylinder_vertex_indices[2]);
        vertex_indices.push_back(cylinder_vertex_indices[4]); 
        vertex_indices.push_back(cylinder_vertex_indices[5]);
        vertex_indices.push_back(cylinder_vertex_indices[4]);
        vertex_indices.push_back(cylinder_vertex_indices[2]);        
    }
    cylinder_vertex_indices[1] = basIndex + 2 * ang_num;
    cylinder_vertex_indices[2] = basIndex + 2 * ang_num + 1;
    cylinder_vertex_indices[3] = basIndex + 2;
    cylinder_vertex_indices[4] = basIndex + 3;
    vertex_indices.push_back(cylinder_vertex_indices[0]);
    vertex_indices.push_back(cylinder_vertex_indices[1]);
    vertex_indices.push_back(cylinder_vertex_indices[3]);
    vertex_indices.push_back(cylinder_vertex_indices[1]);
    vertex_indices.push_back(cylinder_vertex_indices[2]);
    vertex_indices.push_back(cylinder_vertex_indices[3]);
    vertex_indices.push_back(cylinder_vertex_indices[3]);
    vertex_indices.push_back(cylinder_vertex_indices[2]);
    vertex_indices.push_back(cylinder_vertex_indices[4]);
    vertex_indices.push_back(cylinder_vertex_indices[5]);
    vertex_indices.push_back(cylinder_vertex_indices[4]);
    vertex_indices.push_back(cylinder_vertex_indices[2]);
}


std::vector<float> Arrow::getNorm(std::vector<float> Position, std::vector<int> vertex_indices)
{
    glm::vec3 a1, a2, a3;
    glm::vec3 temp, normal;
    std::vector<glm::vec3> vertices;
    std::vector<float>normals;
    std::vector<glm::vec3> faceNormal;
    for (int i = 0; i < Position.size() / 3; ++i)
    {
        temp.x = Position[3 * i];
        temp.y = Position[3 * i + 1];
        temp.z = Position[3 * i + 2];
        vertices.push_back(temp);
    }
    for (int i = 0; i < vertex_indices.size()-2; i += 3) {
        a1 = vertices[vertex_indices[i + 1]] - vertices[vertex_indices[i]];
        a2 = vertices[vertex_indices[i + 2]] - vertices[vertex_indices[i]];
        faceNormal.push_back(glm::cross(a1, a2));
    }
    for (int i = 0; i < Position.size() / 3; ++i)
    {
        normal = glm::vec3(0.0f, 0.0f, 0.0f);
        for (int j = 0; j < vertex_indices.size(); ++j)
        {
            if (vertex_indices[j] == i)
            {
                int index = j / 3;
                normal += faceNormal[index];
            }
        }
        glm::normalize(normal);
        normals.push_back(normal.x);
        normals.push_back(normal.y);
        normals.push_back(normal.z);
    }
    return normals;
}



