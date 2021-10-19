#ifndef COORDINATE_H
#define COORDINATE_H

#include<vector>
#include"../external/glm/glm.hpp"
#include"../external/shader.h"
#include"camera.h"

class CoordinateSystem
{
public:
    CoordinateSystem();
    void draw(Camera* camera, glm::vec3& lightPos);

private:
    Shader* shader;
    unsigned int VAO0, VBO0[2], EBO0;
    unsigned int VAO1, VBO1[2], EBO1;
    unsigned int VAO2, VBO2[2], EBO2;
    void genBuffer();
    void setBufferData();
    Light light;
    std::vector<float>vertices0;
    std::vector<float>normals0;
    std::vector<int>vertex_indices0;
    std::vector<float>vertices1;
    std::vector<float>normals1;
    std::vector<int>vertex_indices1;
    std::vector<float>vertices2;
    std::vector<float>normals2;
    std::vector<int>vertex_indices2;
};

class Arrow
{
public:
    Arrow(int type);
    std::vector<float>vertices;
    std::vector<float>normals;
    std::vector<int>vertex_indices;


private:
    void setArrow(int type);
    void setCylinder(int type);
    std::vector<float> getNorm(std::vector<float> Position, std::vector<int> vertex_indices);
};



#endif // !COORDINATE_H

#pragma once