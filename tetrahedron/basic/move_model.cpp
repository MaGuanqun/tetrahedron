#include"move_model.h"


MoveModel::MoveModel()
{
	setBandRotate();
}

void MoveModel::updateColliderPosition(std::vector<Collider>& collider)
{
	for (unsigned int i = 0; i < collider.size(); ++i) {
		collider[i].mesh_struct.vertex_for_render = collider[i].mesh_struct.vertex_position;
	}
}


void MoveModel::sceneRotateCapsule(int t, std::vector<std::array<double, 3>>& capsule_vertices, MeshStruct* band_mesh_struct, bool use_PD)
{
	double body_center[3];
	moveCapsule(t, capsule_vertices);
	if (t > 100) {
		computerModelCenter(body_center, capsule_vertices);
		rotateCapsule(t, body_center, capsule_vertices);
	}
	moveBand(t, band_mesh_struct, use_PD);
}


void MoveModel::computerModelCenter(double* body_center, std::vector<std::array<double, 3>>& vertices)
{
	memset(body_center, 0, 24);
	double* a = vertices[0].data();
	for (unsigned int i = 0; i < vertices.size(); ++i) 
	{
		body_center[0] += *(a++);
		body_center[1] += *(a++);
		body_center[2] += *(a++);
	}
	DEV_(body_center, (double)vertices.size());
}

void MoveModel::rotateCapsule(int t, double* body_center, std::vector<std::array<double,3>>& vertices)
{
	double move[3];
	memcpy(move, body_center, 24);
	Vector3d body_pos;
	if (t % 1210 < 430) {
		for (int i = 0; i < vertices.size(); ++i) {
			body_pos[0] = vertices[i][0] - move[0];
			body_pos[1] = vertices[i][1] - move[1];
			body_pos[2] = vertices[i][2] - move[2];// -body_center_y
			body_pos = capsule_rotate * body_pos;
			vertices[i][0] = body_pos[0] + move[0];
			vertices[i][1] = body_pos[1] + move[1];// +body_center_y;
			vertices[i][2] = body_pos[2] + move[2];
		}
	}
	else if (t % 1210 > 900) {
		for (int i = 0; i < vertices.size(); ++i) {
			body_pos[0] = vertices[i][0] - move[0];
			body_pos[1] = vertices[i][1] - move[1];
			body_pos[2] = vertices[i][2] - move[2];// -body_center_y
			body_pos = capsule_rotate_reverse * body_pos;
			vertices[i][0] = body_pos[0] + move[0];
			vertices[i][1] = body_pos[1] + move[1];// +body_center_y;
			vertices[i][2] = body_pos[2] + move[2];
		}
	}
}


void MoveModel::moveCapsule(int t, std::vector<std::array<double, 3>>& vertices)
{
	double move_dis = 0.014;
	if (t < 100) {
		for (unsigned int i = 0; i < vertices.size(); ++i) {
			vertices[i][1] -= move_dis;
		}
	}
}


void MoveModel::moveBand(int t, MeshStruct* mesh_struct, bool use_PD)
{
	//double move_dis = 8e-3;
	double move_dis = 0.014;
	//double move_dis = 0.014;
	if (use_PD) {
		if (t < 100) {
			for (int i = 0; i < 4; ++i) {
				if (i % 2 == 0) {
					mesh_struct->anchor_position[i][2] += move_dis;
				}
				else {
					mesh_struct->anchor_position[i][2] -= move_dis;
				}
			}
		}
	}
	else {
		if (t < 100) {
			for (int i = 0; i < 4; ++i) {
				if (i % 2 == 0) {
					mesh_struct->vertex_position[mesh_struct->anchor_vertex[i]][2] += move_dis;
				}
				else {
					mesh_struct->vertex_position[mesh_struct->anchor_vertex[i]][2] -= move_dis;
				}
			}
		}
	}

}



void MoveModel::setBandRotate()
{
	Vector3d rotate_axe(0.0, 1.0, 0.0);
	Matrix3d ux;
	ux << 0, -rotate_axe[2], rotate_axe[1], 
		rotate_axe[2], 0, -rotate_axe[0], 
		-rotate_axe[1], rotate_axe[0], 0;
	Matrix3d uxu;
	uxu << rotate_axe[0] * rotate_axe[0], rotate_axe[0] * rotate_axe[1], rotate_axe[0] * rotate_axe[2],
		rotate_axe[0] * rotate_axe[1], rotate_axe[1] * rotate_axe[1], rotate_axe[1] * rotate_axe[2],
		rotate_axe[0] * rotate_axe[2], rotate_axe[1] * rotate_axe[2], rotate_axe[2] * rotate_axe[2];
	double angle = 0.01 * M_PI;
	band_rotate = cos(angle) *Matrix3d::Identity() + sin(angle) * ux + (1 - cos(angle)) * uxu;
	angle = -angle;
	band_rotate_reverse = cos(angle) * Matrix3d::Identity() + sin(angle) * ux + (1 - cos(angle)) * uxu;


	rotate_axe = Vector3d (0.0, 1.0, 0.0);
	ux << 0, -rotate_axe[2], rotate_axe[1], 
		rotate_axe[2], 0, -rotate_axe[0], 
		-rotate_axe[1], rotate_axe[0], 0;
	uxu << rotate_axe[0] * rotate_axe[0], rotate_axe[0] * rotate_axe[1], rotate_axe[0] * rotate_axe[2],
		rotate_axe[0] * rotate_axe[1], rotate_axe[1] * rotate_axe[1], rotate_axe[1] * rotate_axe[2],
		rotate_axe[0] * rotate_axe[2], rotate_axe[1] * rotate_axe[2], rotate_axe[2] * rotate_axe[2];
	angle = 0.01 * M_PI;
	capsule_rotate = cos(angle) * Matrix3d::Identity() + sin(angle) * ux + (1 - cos(angle)) * uxu;
	angle = -angle;
	capsule_rotate_reverse = cos(angle) * Matrix3d::Identity() + sin(angle) * ux + (1 - cos(angle)) * uxu;

}