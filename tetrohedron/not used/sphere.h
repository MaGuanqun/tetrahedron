#ifndef EARTH_H
#define EARTH_H
#include<vector>
#include"../basic/global.h"
#include"../basic/camera.h"
#include"../external/shader.h"
#include"../external/Eigen/Dense"
#include"../basic/eigenDenseOperation.h"
#include"../meshStruct.h"
//#define STB_IMAGE_IMPLEMENTATION
using namespace Eigen;
using namespace denseOperation;

class Sphere
{
public:
	double R;
	Sphere();
	void setMesh(OriMesh& mesh);
	void draw(Camera* camera);
	void drawShadow(Camera* camera, Shader* shader);
	Shader* shader;
private:
	struct SpherePoint
	{
		Vector3d position;
		Vector2d texture_coord;
		Vector3d normal;
		std::vector<int>face_index;
		SpherePoint(Vector3d& ori_position, Vector2d& ori_texture_coord) {
			memcpy(position.data(), ori_position.data(), 24);
			memcpy(texture_coord.data(), ori_texture_coord.data(), 16);
		};
	};
	
	unsigned int VAO0, VBO0[3], EBO0;
	std::vector<SpherePoint>vertices;
	std::vector<double>vertex_position;
	std::vector<double>vertex_texture_coordinate;
	std::vector<double>vertex_normal;
	std::vector<int>indices;
	int latitude_num;
	int longitude_num;
	int vertex_num;
	
	//void loadTexture();
	void creatGlobe();
	void creatNormal();
	void genBuffer();
	void setBufferData();
	void setVertices();
	//nsigned int texture1;
};


#endif // !EARTH_H
#pragma once