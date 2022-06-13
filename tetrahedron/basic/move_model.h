#pragma once

#include"../mesh_struct/mesh_struct.h"
#include"../external/Eigen/Dense"

using namespace Eigen;

class MoveModel
{
public:
	MoveModel();
	void rotateCapsule(int t, double* body_center, std::vector<std::array<double, 3>>& vertices);
	void moveBand(int t, MeshStruct* mesh_struct, bool use_PD);
private:
	Matrix3d band_rotate;
	Matrix3d band_rotate_reverse;
	Matrix3d capsule_rotate;
	Matrix3d capsule_rotate_reverse;
	void setBandRotate();
};
