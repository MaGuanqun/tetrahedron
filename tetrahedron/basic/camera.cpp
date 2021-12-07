#include"camera.h"

Camera::Camera(glm::vec3 position, glm::vec3 up, glm::vec3 center, float fov):
    position(position), center(center), ori_position(position), ori_center(center), fov(fov)
{
    far_plane = 3.0 * sqrt(dot(position - center, position - center));
    near_plane = 0.01f;
    this->up = glm::normalize(up);
    ori_up = this->up;
    maxHightValue = sqrt(dot(position - center, position - center)) * tan(glm::radians(fov / 2.0f));
    right = maxHightValue * (float)SCR_WIDTH / (float)SCR_HEIGHT;
    y_vec = glm::normalize(glm::cross(up, position - center));
}

void Camera::resetCam()
{
    position = ori_position;
    center = ori_center;
    up = ori_up;
    maxHightValue = sqrt(dot(position - center, position - center)) * tan(glm::radians(fov / 2.0f));
    right = maxHightValue * (float)SCR_WIDTH / (float)SCR_HEIGHT;
    y_vec = glm::normalize(glm::cross(up, position - center));
}

void Camera::updateCamera(glm::vec3 position, glm::vec3 up, glm::vec3 center)
{
    this->position = position;
    this->up = up;
    this->center = center;
    ori_position = position;
    ori_center = center;
    ori_up = up;
    far_plane = 2.0 * sqrt(dot(position - center, position - center));
    near_plane = 0.01f;
    maxHightValue = sqrt(dot(position - center, position - center)) * tan(glm::radians(fov / 2.0f));
    right = maxHightValue * (float)SCR_WIDTH / (float)SCR_HEIGHT;
    y_vec = glm::normalize(glm::cross(up, position - center));
}

void Camera::zoomInOut(float zoomValue)
{
    // glm::vec3 temVec3 =Center-Position;       
    glm::vec3 position_center = position - center;
    glm::vec3 temVec32 = glm::normalize(position_center) * zoomValue;
    float valueChange = sqrt(glm::dot(temVec32, temVec32)) / sqrt(dot(position_center, position_center));
    position = temVec32 + center;
    //Center = Position + valueChange * temVec3;
    maxHightValue = sqrt(dot(position_center, position_center)) * tan(glm::radians(fov / 2.0f));
    right = maxHightValue * (float)SCR_WIDTH / (float)SCR_HEIGHT;
}

void Camera::getCursorPosInSpace(double* cursor_pos_in_space, double* cursor_screen, double* object_position)
{
    glm::vec3 object_from_camera;
    object_from_camera.x = object_position[0] - position.x;
    object_from_camera.y = object_position[1] - position.y;
    object_from_camera.z = object_position[2] - position.z;
    glm::vec3 position2center=center-position;
    position2center=normalize(position2center);
    float dis = dot(object_from_camera, position2center);
    glm::vec3 object_plane_center = position2center * dis;
    double max_hight = 2.0 * dis * tan(glm::radians(fov / 2.0f));
    double max_width = max_hight * (float)SCR_WIDTH / (float)SCR_HEIGHT;
    object_plane_center += position;
    double cursor_pos[2];
    cursor_pos[0] = (cursor_screen[0] - (double)(SCR_WIDTH / 2)) / (double)(SCR_WIDTH)*max_width;
    cursor_pos[1] = ((double)(SCR_HEIGHT / 2) - cursor_screen[1]) / (double)(SCR_HEIGHT)*max_hight;
    glm::vec3 x_ = normalize(glm::cross(position2center, up));
    cursor_pos_in_space[0] = object_plane_center[0] + cursor_pos[0] * x_[0] + cursor_pos[1] * up[0];
    cursor_pos_in_space[1] = object_plane_center[1] + cursor_pos[0] * x_[1] + cursor_pos[1] * up[1];
    cursor_pos_in_space[2] = object_plane_center[2] + cursor_pos[0] * x_[2] + cursor_pos[1] * up[2];
}


glm::vec3 Camera::rotationxyz(float angle, glm::vec3& rotate_axe, glm::vec3& vec)
{
    rotate_axe = glm::normalize(rotate_axe);
    glm::mat3x3 ux = glm::mat3x3(0, -rotate_axe.z, rotate_axe.y, rotate_axe.z, 0, -rotate_axe.x, -rotate_axe.y, rotate_axe.x, 0);

    glm::mat3x3 uxu = glm::mat3x3(rotate_axe.x * rotate_axe.x, rotate_axe.x * rotate_axe.y, rotate_axe.x * rotate_axe.z,
        rotate_axe.x * rotate_axe.y, rotate_axe.y * rotate_axe.y, rotate_axe.y * rotate_axe.z,
        rotate_axe.x * rotate_axe.z, rotate_axe.y * rotate_axe.z, rotate_axe.z * rotate_axe.z);

    glm::mat3x3 Rot = cos(angle) * (glm::mat3(1.0f)-uxu) + sin(angle) * ux + uxu;
    return Rot * vec;
}

void Camera::move(float y_pos, float z_pos)
{
    float dis =2.0f * tan(glm::radians(fov / 2.0f)) * sqrt(glm::dot(position - center, position - center));
    glm::vec3 move = y_pos * y_vec / (float)SCR_HEIGHT *dis
        + z_pos / (float)SCR_HEIGHT * up * dis;
    position += move;
    center += move;
}


void Camera::rotation(float angleY, float angleZ, int type) //type to note if the cursor is in the circle
{
    glm::vec3 rotate_axe;
    glm::vec3 position_center = position - center;
    if (type == 1)
    {
        glm::vec3 tem = up - glm::normalize(position_center);

        if (fabs(tem.x) + fabs(tem.y) + fabs(tem.z) > 2e-5f) {
            glm::vec3 y = glm::cross(position_center, up);
            position_center = rotationxyz(angleZ, y, position_center);
            up = rotationxyz(angleZ, y, up);
            rotate_axe = glm::cross(y, position_center);
            position_center = rotationxyz(angleY, rotate_axe, position_center);
            up = rotationxyz(angleY, rotate_axe, up);
        }
    }
    else if (type == 2)
    {
        rotate_axe = position_center;
        //          
        up = rotationxyz(angleY, rotate_axe, up);
        //   Center = rotationxyz(angleY, rotate_axe, Center);
        up = rotationxyz(angleZ, rotate_axe, up);
        //   Center = rotationxyz(angleZ, rotate_axe, Center);
    }
    position = position_center + center;
    up = glm::normalize(up);
    y_vec= glm::normalize(glm::cross(up, position - center));
}