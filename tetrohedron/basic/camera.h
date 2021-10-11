#ifndef CAMERA_H
#define CAMERA_H

#include"../external/glm/glm.hpp"
#include"../external/glm/gtc/matrix_transform.hpp"
#include<iostream>
#include"global.h"

struct Light {
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
};

class Camera
{
private:

    float right;
    glm::vec3 ori_position;
    glm::vec3 ori_up;
    glm::vec3 ori_center;
    glm::vec3 rotationxyz(float angle, glm::vec3& rotate_axe, glm::vec3& vec);
public:
	glm::vec3 position;
	glm::vec3 up;
    glm::vec3 center;
    glm::vec3 y_vec;
 
    float fov;   
    float maxHightValue;//We want to know in our view, the max height value that can show when the plane is cross the origianl point.   
    float far_plane;
    float near_plane;
    Camera(glm::vec3 position = glm::vec3(5.0f, 0.0f, 0.0f), glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3 center = glm::vec3(0.0f, 0.0f, 0.0f), float fov = 45.0f);

    void updateCamera(glm::vec3 position, glm::vec3 up, glm::vec3 center);

    void resetCam();
    void zoomInOut(float zoomValue);
    void getCursorPosInSpace(double* cursor_pos_in_space, double* cursor_screen, double* object_position);

    glm::mat4 GetViewMatrix() {return glm::lookAt(position, center, up);}
    glm::mat4 GetProjectMatrix() { return glm::perspective(glm::radians(fov), (float)SCR_WIDTH / (float)SCR_HEIGHT, near_plane, far_plane); }

    glm::mat4 GetOrthogonalMatrix(){return glm::ortho(-right, right, -maxHightValue, maxHightValue, near_plane, far_plane);}

    void move(float y_pos, float z_pos);
    void rotation(float angleY, float angleZ, int type);  

};
#endif