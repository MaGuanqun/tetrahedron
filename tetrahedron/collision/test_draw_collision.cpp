#include"test_draw_collision.h"


void TestDrawCollision::initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider,
	std::vector<Tetrahedron>* tetrahedron, Thread* thread, Floor* floor, double* tolerance_ratio,
	Collision* pd_collision, Collision* pbd_collision, Collision* newton_collision, Collision* xpbd_ipc_collision, unsigned int simulation_method)
{
	this->cloth = cloth;
	this->collider = collider;
	this->tetrahedron = tetrahedron;
	this->thread = thread;
	total_obj_num = cloth->size() + tetrahedron->size();
	vertex_index_in_one_cell.resize(total_obj_num);
	triangle_index_in_one_cell.resize(total_obj_num);
	edge_index_in_one_cell.resize(total_obj_num);


	//std::cout << "iniital test " << &pd_collision->point_triangle_target_pos_index << std::endl;

	switch (simulation_method)
	{
	case PD_:
		collision = pd_collision;
			break;
	case XPBD_:
		collision = pbd_collision;
		break;
	case NEWTON_:
		collision = newton_collision;
		break;
	case XPBD_IPC_:
		collision = xpbd_ipc_collision;
		break;
	default:
		collision = new Collision();
		collision->initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio,simulation_method,false);
		collision->spatial_hashing.initialOriHashValue();
		break;
	}
	draw_collision.initial(cloth, collider, tetrahedron, thread);

	//std::cout << &collision->point_triangle_target_pos_index << std::endl;
	draw_collision.setInPairInfo(&collision->point_triangle_target_pos_index, &collision->point_triangle_collider_target_pos_index, &collision->edge_edge_target_pos_index);
}


void TestDrawCollision::setCollisionData()
{
	collision->collisionCulling();
	collision->getCollisionPair();
	draw_collision.setElementIndices();
}


void TestDrawCollision::obtianSpatialHashingCell(Camera* camera, double* cursor_screen)
{
	double start_pos[3] = { camera->position.x, camera->position.y, camera->position.z };
	double dir[3];
	camera->getCursorPosCameraCenterPlane(dir, cursor_screen);
	SUB_(dir, start_pos);
	normalize(dir);
	collision->spatial_hashing.selectCell(start_pos, dir, select_hash_index, select_ori_hash_index);
	//std::cout << "TestDrawCollision::obtianSpatialHashingCell " << select_hash_index.size() << std::endl;

}


bool TestDrawCollision::isIndexNotEmpty()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		if (!triangle_index_in_one_cell[i].empty()) {
			return true;
		}
		if (!vertex_index_in_one_cell[i].empty()) {
			return true;
		}
		if (!edge_index_in_one_cell[i].empty()) {
			return true;
		}
	}
	return false;
}


void TestDrawCollision::obtainElementsInOneCell(int& index_chosen)
{

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_index_in_one_cell[i].clear();
		triangle_index_in_one_cell[i].clear();
		edge_index_in_one_cell[i].clear();
	}
	if (!select_hash_index.empty()) {
		index_chosen = index_chosen % select_hash_index.size();
		if (index_chosen < 0) {
			index_chosen += select_hash_index.size();
		}
		while (true)
		{
			collision->spatial_hashing.findAllElementsInOneCell(select_hash_index[index_chosen], select_ori_hash_index[index_chosen], vertex_index_in_one_cell,
				triangle_index_in_one_cell, edge_index_in_one_cell);
			if (isIndexNotEmpty()) {
				break;
			}
			else {
				index_chosen++;
			}
			if (index_chosen >= select_hash_index.size()) {
				break;
			}
		}
		//std::cout << index_chosen << std::endl;

		draw_collision.setElementInAllCell(vertex_index_in_one_cell, triangle_index_in_one_cell, edge_index_in_one_cell);
		draw_collision.setElementInOneCell(triangle_index_in_one_cell, edge_index_in_one_cell);
	}
}

void TestDrawCollision::setForOriSpatialHashing()
{
	collision->buildSpatialHashingForOri();	
	draw_spatial_hashing.setCellData(&collision->spatial_hashing.ori_hash_value, collision->spatial_hashing.cell_length, collision->spatial_hashing.cell_number, collision->spatial_hashing.scene_aabb);
}

void TestDrawCollision::setForSelectCell()
{
	draw_spatial_hashing.setCellData(select_ori_hash_index, collision->spatial_hashing.cell_length, collision->spatial_hashing.cell_number, collision->spatial_hashing.scene_aabb);

}
void TestDrawCollision::setForSelectCellOne(int index)
{
	if (!select_hash_index.empty()) {
		draw_spatial_hashing.setCellData(select_ori_hash_index[index], collision->spatial_hashing.cell_length, collision->spatial_hashing.cell_number, collision->spatial_hashing.scene_aabb);
	}
}

void TestDrawCollision::drawCollision(bool draw_VT, Light& light, Camera* camera, Shader* object_shader_front, 
	std::vector<std::vector<bool>>& show_element, Shadow* shadow, Shader* wireframe_shader, bool draw_all_collision_pair, bool draw_all_element)
{
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_CUBE_MAP, shadow->depth_map);
	glEnable(GL_CULL_FACE);
	for (int j = 0; j < cloth->size(); ++j) {
		if (show_element[9 + CLOTH_][j]) {
			cloth->data()[j].setSceneShader(light, camera, shadow->far_plane, object_shader_front);
			cloth->data()[j].drawOriPos(camera, object_shader_front);
		}
	}
	glCullFace(GL_BACK);
	for (int j = 0; j < tetrahedron->size(); ++j) {
		if (show_element[9 + TETRAHEDRON_][j]) {
			tetrahedron->data()[j].setSceneShader(light, camera, shadow->far_plane, object_shader_front);
			tetrahedron->data()[j].drawOriPos(camera, object_shader_front);
		}
	}
	for (int j = 0; j < collider->size(); ++j) {
		if (show_element[9 + COLLIDER_][j]) {
			collider->data()[j].setSceneShader(light, camera, shadow->far_plane, object_shader_front);
			collider->data()[j].draw(camera, object_shader_front);
		}
	}
	glDisable(GL_CULL_FACE);
	for (int i = 0; i < collider->size(); ++i) {
		if (show_element[12 + COLLIDER_][i]) {
			collider->data()[i].drawWireframeOriPos(camera, wireframe_shader);
		}
	}
	for (int j = 0; j < cloth->size(); ++j) {
		if (show_element[12 + CLOTH_][j]) {
			cloth->data()[j].drawWireframeOriPos(camera, wireframe_shader);
		}
	}
	for (int j = 0; j < tetrahedron->size(); ++j) {
		if (show_element[12 + TETRAHEDRON_][j]) {
			tetrahedron->data()[j].drawWireframeOriPos(camera, wireframe_shader);
		}
	}
	if (draw_all_collision_pair) {
		draw_collision.drawCollision(draw_VT, light, camera, object_shader_front, wireframe_shader, show_element);
	}
	else {
		draw_collision.drawCollisionCell(draw_VT, light, camera, object_shader_front, wireframe_shader, show_element,draw_all_element);
	}

}



