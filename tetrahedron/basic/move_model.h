#pragma once

#include"../object/collider.h"
#include"../mesh_struct/mesh_struct.h"
#include"../external/Eigen/Dense"

using namespace Eigen;

class MoveModel
{
public:
	void sceneRotateCapsule(int t, std::vector<std::array<double, 3>>& ori_capsule_vertices, std::vector<std::array<double, 3>>& capsule_vertices, 
		MeshStruct* band_mesh_struct, bool use_PD, double sub_step_size);
	void updateColliderPosition(std::vector<Collider>& collider);
	void moveSkirt(int t, std::vector<MeshStruct*>& mesh_struct, bool use_PD, double sub_step_size);
	void moveSphere(int t, std::vector<std::array<double, 3>>& ori_capsule_vertices, std::vector<std::array<double, 3>>& capsule_vertices,
		double sub_step_size);

private:
	void setRotate(double step_size,  Matrix3d& capsule_rotate);
	void rotateCapsule(int t, double* body_center, std::vector<std::array<double, 3>>& ori_vertices, std::vector<std::array<double, 3>>& vertices,
		Matrix3d& capsule_rotate, Matrix3d& capsule_rotate_reverse);
	void moveBand(int t, MeshStruct* mesh_struct, bool use_PD, double move_dis);
	void computerModelCenter(double* body_center, std::vector<std::array<double, 3>>& vertices);
	void moveCapsule(int t, std::vector<std::array<double, 3>>& ori_vertices, std::vector<std::array<double, 3>>& vertices, double move_dis);
};
