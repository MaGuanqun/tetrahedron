#include"move_model.h"


MoveModel::MoveModel()
{
	setBandRotate();
}


//void MoveModel::


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


void MoveModel::moveBand(int t, MeshStruct* mesh_struct, bool use_PD)
{
	//double move_dis = 8e-3;
	double move_dis = 0.013;
	//double move_dis = 0.014;
	if (use_PD) {
		if (t < 100) {
			for (int i = 0; i < 2; ++i) {
				mesh_struct->anchor_position[i][2] -= move_dis;
			}
			for (int i = 2; i < 4; ++i) {
				mesh_struct->anchor_position[i][2] += move_dis;
			}
		}
	}
	else {
		if (t < 100) {
			for (int i = 0; i < 2; ++i) {
				mesh_struct->vertex_position[mesh_struct->anchor_vertex[i]][2] -= move_dis;
			}
			for (int i = 2; i < 4; ++i) {
				mesh_struct->vertex_position[mesh_struct->anchor_vertex[i]][2] += move_dis;
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
	double angle = 0.02 * M_PI;
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
	angle = 0.1 * M_PI;
	capsule_rotate = cos(angle) * Matrix3d::Identity() + sin(angle) * ux + (1 - cos(angle)) * uxu;
	angle = -angle;
	capsule_rotate_reverse = cos(angle) * Matrix3d::Identity() + sin(angle) * ux + (1 - cos(angle)) * uxu;

}