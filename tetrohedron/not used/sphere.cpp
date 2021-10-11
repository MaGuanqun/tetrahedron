#include"sphere.h"
//#include<learnopengl/stb_image.h>

Sphere::Sphere()
{
    latitude_num = 20;
    longitude_num = 40;
    vertex_num = (longitude_num + 1UL) * (latitude_num - 2UL) + 2UL;
    R = sphere_radius;
    creatGlobe();
    creatNormal();
    setVertices(); 
    //genBuffer();
    //setBufferData();
    ////loadTexture();
    //shader = new Shader("./basic/sphere.vs", "./basic/sphere.fs");
    //*shader = Shader("./basic/sphere.vs", "./basic/sphere.fs");    
}



void Sphere::draw(Camera* camera) {
    shader->use();
    shader->setVec3("viewPos", camera->Position);
    shader->setBool("lightShadowOn", true);
    shader->setMat4("projection", camera->GetProjectMatrix());
    shader->setMat4("view", camera->GetViewMatrix());
    shader->setMat4("model", glm::mat4(1.0));
    glBindVertexArray(VAO0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}
void Sphere::drawShadow(Camera* camera, Shader* shader) {
    shader->use();
    shader->setMat4("projection", camera->GetProjectMatrix());
    shader->setMat4("view", camera->GetViewMatrix());
    shader->setMat4("model", glm::mat4(1.0));
    glBindVertexArray(VAO0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}

void Sphere::creatGlobe()
{
    int longitude_num_1 = longitude_num + 1UL;
    int latitude_num_1 = latitude_num - 1UL;
    vertices.reserve(longitude_num_1 * latitude_num_1);
    indices.reserve(longitude_num * latitude_num);
    for (int j = 1; j < latitude_num_1; ++j)
    {
        for (int i = 0; i < longitude_num_1; ++i)
        {
            vertices.push_back(SpherePoint(Vector3d(R *sin(M_PI * (double)j / (double)latitude_num_1) * cos(2.0*M_PI * (double)i / (double)longitude_num),
                R *sin(M_PI * (double)j / (double)latitude_num_1) * sin(2.0*M_PI * (double)i / (double)longitude_num),
                R * cos(M_PI * (double)j / (double)latitude_num_1)),
                Vector2d((double)i / (double)longitude_num, (double)j / (double)latitude_num_1)));
        }
    }
    vertices.push_back(SpherePoint(Vector3d(0.0, 0.0, R), Vector2d(0.5, 0.0)));
    vertices.push_back(SpherePoint(Vector3d(0.0, 0.0, -R), Vector2d(0.5, 1.0)));
    for (int j = 1; j < latitude_num - 2; ++j)
    {
        for (int i = 0; i < longitude_num; ++i)
        {
            indices.push_back(longitude_num_1 * (j - 1UL) + i);
            indices.push_back(longitude_num_1 * j + i);
            indices.push_back(longitude_num_1 * (j - 1UL) + i + 1UL);
            indices.push_back(longitude_num_1 * j + i);
            indices.push_back(longitude_num_1 * j + i + 1UL);
            indices.push_back(longitude_num_1 * (j - 1UL) + i + 1UL);           
            int indi = indices.size() / 3 - 2;
            vertices[longitude_num_1 * (j - 1UL) + i].face_index.push_back(indi);
            vertices[longitude_num_1 * j + i].face_index.push_back(indi);
            vertices[longitude_num_1 * j + i].face_index.push_back(indi+ 1);
            vertices[longitude_num_1 * (j - 1UL) + i + 1UL].face_index.push_back(indi);
            vertices[longitude_num_1 * (j - 1UL) + i + 1UL].face_index.push_back(indi+ 1);
            vertices[longitude_num_1 * j + i + 1UL].face_index.push_back(indi+ 1);
        }
    }
    for (int i = 0; i < longitude_num; ++i)
    {
        indices.push_back(vertex_num - 2UL);
        indices.push_back(i);        
        indices.push_back(i + 1UL);    
        indices.push_back(longitude_num_1 * (latitude_num - 3UL) + i);
        indices.push_back(vertex_num - 1UL);
        indices.push_back(longitude_num_1 * (latitude_num - 3UL) + i + 1UL);
        
        int indi = indices.size() / 3 - 2;
        vertices[i].face_index.push_back(indi);
        vertices[i+1UL].face_index.push_back(indi);
        vertices[vertex_num - 2UL].face_index.push_back(indi);
        vertices[longitude_num_1 * (latitude_num - 3UL) + i].face_index.push_back(indi + 1);
        vertices[longitude_num_1 * (latitude_num - 3UL) + i + 1UL].face_index.push_back(indi + 1);
        vertices[vertex_num - 1UL].face_index.push_back(indi + 1);
    }
}

void Sphere::creatNormal()
{
    int face_num = indices.size() / 3;
    std::vector<Vector3d>faceNorm;
    faceNorm.reserve(face_num);
    for (int i = 0; i < face_num; ++i) {
        faceNorm.push_back(normalized(crossProduct(vertices[indices[3 * i]].position - vertices[indices[3 * i + 2]].position, vertices[indices[3 * i + 1]].position - vertices[indices[3 * i + 2]].position)));
    }
    Vector3d tempVec3;
    for (int i = 0; i < vertex_num; ++i) {
        tempVec3.setZero();
        for (int j = 0; j < vertices[i].face_index.size(); ++j) {
            tempVec3 += faceNorm[vertices[i].face_index[j]];
        }
        vertices[i].normal = normalized(tempVec3);
    }
}

void Sphere::setVertices()
{
    vertex_position.reserve(3 * vertex_num);
    vertex_texture_coordinate.reserve(2 * vertex_num);
    vertex_normal.reserve(3 * vertex_num);
    for (int i = 0; i < vertex_num; ++i) {
        vertex_position.push_back(vertices[i].position.data()[0]);
        vertex_position.push_back(vertices[i].position.data()[1]);
        vertex_position.push_back(vertices[i].position.data()[2]);
        vertex_texture_coordinate.push_back(vertices[i].texture_coord.data()[0]);
        vertex_texture_coordinate.push_back(vertices[i].texture_coord.data()[1]);
        vertex_normal.push_back(vertices[i].normal.data()[0]);
        vertex_normal.push_back(vertices[i].normal.data()[1]);
        vertex_normal.push_back(vertices[i].normal.data()[2]);
    }
}

void Sphere::setMesh(OriMesh& mesh)
{
    mesh.vertices.resize(vertex_num);
    mesh.face_groups.resize(1);
    mesh.face_groups[0].indices.resize(indices.size());
    mesh.face_groups[0].indices = indices;
    for (int i = 0; i < vertex_num; ++i) {
        mesh.vertices[i][0] = vertex_position[3 * i];
        mesh.vertices[i][1] = vertex_position[3 * i+1];
        mesh.vertices[i][2] = vertex_position[3 * i+2];
    }
    mesh.face_groups[0].front_material.Kd[0] = 0.5;
    mesh.face_groups[0].front_material.Kd[1] = 0.5;
    mesh.face_groups[0].front_material.Kd[2] = 0.5;
    memset(mesh.face_groups[0].front_material.Ka, 0, 24);
    memset(mesh.face_groups[0].front_material.Ks, 0, 24);
    mesh.face_groups[0].front_material.Tf[0] = 1.0;
    mesh.face_groups[0].front_material.Tf[1] = 1.0;
    mesh.face_groups[0].front_material.Tf[2] = 1.0;
    mesh.face_groups[0].front_material.Ni = 1.0;
}

void Sphere::genBuffer()
{
    glGenVertexArrays(1, &VAO0);
    glGenBuffers(3, VBO0);
    glGenBuffers(1, &EBO0);
    //glGenTextures(1, &texture1);
}

void Sphere::setBufferData()
{
    glBindVertexArray(VAO0);
    glBindBuffer(GL_ARRAY_BUFFER, VBO0[0]);
    glBufferData(GL_ARRAY_BUFFER, vertex_position.size() * sizeof(double), &vertex_position[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO0);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(int), &indices[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, VBO0[1]);
    glBufferData(GL_ARRAY_BUFFER, vertex_normal.size() * sizeof(double), &vertex_normal[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, VBO0[2]);
    glBufferData(GL_ARRAY_BUFFER, vertex_texture_coordinate.size() * sizeof(double), &vertex_texture_coordinate[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), (void*)0);
    glBindVertexArray(0);
}

//void Sphere::loadTexture()
//{
//    int width, height, nrChannels;
//    std::string path = "projectImage/worldmap.jpg";
//    unsigned char* data = stbi_load(path.c_str(), &width, &height, &nrChannels, 0);
//    if (data)
//    {
//        GLenum format;
//        if (nrChannels == 1)
//            format = GL_RED;
//        else if (nrChannels == 3)
//            format = GL_RGB;
//        else if (nrChannels == 4)
//            format = GL_RGBA;
//        glBindTexture(GL_TEXTURE_2D, texture1);
//        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
//        glGenerateMipmap(GL_TEXTURE_2D);
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
//        stbi_image_free(data);
//    }
//}

