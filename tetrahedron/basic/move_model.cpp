#include"move_model.h"



void MoveModel::updateColliderPosition(std::vector<Collider>& collider)
{
	for (unsigned int i = 0; i < collider.size(); ++i) {
		memcpy(collider[i].mesh_struct.vertex_for_render[0].data(), collider[i].mesh_struct.vertex_position[0].data(), 24 * collider[i].mesh_struct.vertex_position.size());	
	}
}


void MoveModel::moveSkirt(int t, std::vector<MeshStruct*>& mesh_struct, bool use_PD, double sub_step_size)
{
	double move_dis = 0.025;
	double move_dis2 = 0.003;

	double move_dis3 = move_dis/sub_step_size;
	double move_dis4 = move_dis2/sub_step_size;

	if (t < 500) {
		if (use_PD) {
			for (unsigned int j = 0; j < mesh_struct.size(); ++j) {
				for (int i = 0; i < mesh_struct[j]->anchor_position.size(); ++i) {
					if (t % 100 < 50) {
						mesh_struct[j]->anchor_position[i][0] += move_dis3;
					}
					else {
						mesh_struct[j]->anchor_position[i][0] -= move_dis3;
					}
					if (t % 80 < 40) {
						mesh_struct[j]->anchor_position[i][1] += move_dis4;
					}
					else {
						mesh_struct[j]->anchor_position[i][1] -= move_dis4;
					}
				}
			}
		}
		else {
			for (unsigned int j = 0; j < mesh_struct.size(); ++j) {
				for (int i = 0; i < mesh_struct[j]->anchor_vertex.size(); ++i) {
					if (t % 100 < 50) {
						mesh_struct[j]->vertex_position[mesh_struct[j]->anchor_vertex[i]][0] += move_dis3;
					}
					else {
						mesh_struct[j]->vertex_position[mesh_struct[j]->anchor_vertex[i]][0] -= move_dis3;
					}
					if (t % 80 < 40) {
						mesh_struct[j]->vertex_position[mesh_struct[j]->anchor_vertex[i]][1] += move_dis4;
					}
					else {
						mesh_struct[j]->vertex_position[mesh_struct[j]->anchor_vertex[i]][1] -= move_dis4;
					}
				}
			}
		}


	}
	
}


void MoveModel::moveSphere(int t, std::vector<std::array<double, 3>>& ori_capsule_vertices, std::vector<std::array<double, 3>>& capsule_vertices, 
	double sub_step_size)
{
	if (t<200) {
		Matrix3d rotate;
		setRotate(0.01 / sub_step_size, rotate);
		Vector3d body_pos;	
		double body_center[3];
		computerModelCenter(body_center, ori_capsule_vertices);
		for (int i = 0; i < ori_capsule_vertices.size(); ++i) {
			SUB(body_pos, ori_capsule_vertices[i], body_center);
			body_pos = rotate * body_pos;
			SUM(capsule_vertices[i], body_pos, body_center);
			
		}
	}
	else if (t < 350) {
		double move_dis = 0.003 / sub_step_size;
		for (int i = 0; i < ori_capsule_vertices.size(); ++i) {
			capsule_vertices[i][1] = ori_capsule_vertices[i][1] + move_dis;
		}
	}
}



void MoveModel::sceneRotateCapsule(int t, std::vector<std::array<double, 3>>& ori_capsule_vertices, 
	std::vector<std::array<double, 3>>& capsule_vertices, MeshStruct* band_mesh_struct, bool use_PD, double sub_step_size)
{
	Matrix3d capsule_rotate;
	Matrix3d capsule_rotate_reverse;
	setRotate(0.01 / sub_step_size, capsule_rotate);
	setRotate(-0.01 / sub_step_size, capsule_rotate_reverse);
	double body_center[3];
	moveCapsule(t, ori_capsule_vertices, capsule_vertices,0.014 / sub_step_size);
	if (t > 100) {
		computerModelCenter(body_center, ori_capsule_vertices);
		rotateCapsule(t, body_center, ori_capsule_vertices,  capsule_vertices, capsule_rotate, capsule_rotate_reverse);
	}
	moveBand(t, band_mesh_struct, use_PD,0.014 / sub_step_size);
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

void MoveModel::rotateCapsule(int t, double* body_center, std::vector<std::array<double,3>>& ori_vertices, 
	std::vector<std::array<double, 3>>& vertices, Matrix3d& capsule_rotate, Matrix3d& capsule_rotate_reverse)
{
	double move[3];
	memcpy(move, body_center, 24);
	Vector3d body_pos;
	if (t % 1210 < 430) {
		for (int i = 0; i < vertices.size(); ++i) {
			body_pos[0] = ori_vertices[i][0] - move[0];
			body_pos[1] = ori_vertices[i][1] - move[1];
			body_pos[2] = ori_vertices[i][2] - move[2];// -body_center_y
			body_pos = capsule_rotate * body_pos;
			vertices[i][0] = body_pos[0] + move[0];
			vertices[i][1] = body_pos[1] + move[1];// +body_center_y;
			vertices[i][2] = body_pos[2] + move[2];
		}
	}
	else if (t % 1210 > 600 && t%1210<1030) {
		//std::cout << "test " << std::endl;
		for (int i = 0; i < vertices.size(); ++i) {
			body_pos[0] = ori_vertices[i][0] - move[0];
			body_pos[1] = ori_vertices[i][1] - move[1];
			body_pos[2] = ori_vertices[i][2] - move[2];// -body_center_y
			body_pos = capsule_rotate_reverse * body_pos;
			vertices[i][0] = body_pos[0] + move[0];
			vertices[i][1] = body_pos[1] + move[1];// +body_center_y;
			vertices[i][2] = body_pos[2] + move[2];
		}
	}
}


void MoveModel::moveCapsule(int t, std::vector<std::array<double, 3>>& ori_vertices,  std::vector<std::array<double, 3>>& vertices, double move_dis)
{
	if (t < 100) {
		for (unsigned int i = 0; i < vertices.size(); ++i) {
			vertices[i][1] = ori_vertices[i][1] - move_dis;
		}
	}
}


void MoveModel::moveBand(int t, MeshStruct* mesh_struct, bool use_PD, double move_dis)
{
	//double move_dis = 8e-3;
	//double move_dis = 0.014;
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
					mesh_struct->vertex_position[mesh_struct->anchor_vertex[i]][2] =mesh_struct->vertex_for_render[mesh_struct->anchor_vertex[i]][2] + move_dis;
				}
				else {
					mesh_struct->vertex_position[mesh_struct->anchor_vertex[i]][2] = mesh_struct->vertex_for_render[mesh_struct->anchor_vertex[i]][2]  - move_dis;
				}
			}
		}
	}

}



void MoveModel::setRotate(double step_size, Matrix3d& capsule_rotate)
{
	Vector3d rotate_axe(0.0, 1.0, 0.0);
	Matrix3d ux;
	Matrix3d uxu;
	double angle;
	ux << 0, -rotate_axe[2], rotate_axe[1], 
		rotate_axe[2], 0, -rotate_axe[0], 
		-rotate_axe[1], rotate_axe[0], 0;
	uxu << rotate_axe[0] * rotate_axe[0], rotate_axe[0] * rotate_axe[1], rotate_axe[0] * rotate_axe[2],
		rotate_axe[0] * rotate_axe[1], rotate_axe[1] * rotate_axe[1], rotate_axe[1] * rotate_axe[2],
		rotate_axe[0] * rotate_axe[2], rotate_axe[1] * rotate_axe[2], rotate_axe[2] * rotate_axe[2];
	angle = step_size * M_PI;
	capsule_rotate = cos(angle) * Matrix3d::Identity() + sin(angle) * ux + (1 - cos(angle)) * uxu;
	//angle = -angle;
	//capsule_rotate_reverse = cos(angle) * Matrix3d::Identity() + sin(angle) * ux + (1 - cos(angle)) * uxu;
}