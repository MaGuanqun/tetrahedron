#include"move_object.h"

void MoveObject::initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
	Thread* thread)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	this->thread = thread;
	total_obj_num = cloth->size() + tetrahedron->size() + collider->size();
	tetrahedron_end_index = cloth->size() + tetrahedron->size();
	thread_num = thread->thread_num;


	total_displacement.resize(total_obj_num);
	memset(total_displacement[0].data(), 0, 24 * total_obj_num);

	
}

void MoveObject::moveScript(unsigned int type)
{
	displacement_test.resize(2);
	displacement_test_render.resize(2);

	if (type == 1) {
		double a[3] = { -0.510271, 0.631315, 0.0321853 };
		double b[3] = { -0.52193, 0.580719, 0.0493092 };
		SUB(displacement_test[0], a, displacement_test[0]);
		SUB(displacement_test[1], b, displacement_test[1]);
	}
	else if (type == 2) {
		displacement_test_render[0] = { 0, 0, 0 };
		displacement_test_render[1] = { -0.190714, 0.17638, -0.324351 };
		//displacement_test_render[1] = { 0, 0.0, 0.0 };

		displacement_test[0] = { 0, 0, 0 };
		//displacement_test[1] = { 0, -0.1, 0 };
		displacement_test[1] = { -0.191217, 0.176467, -0.316329 };
	}
	else if (type == 3) {
		double a[3] = { -0.510271, 0.631315, 0.0321853 };
		double b[3] = { -0.52193, 0.580719, 0.0493092 };
		SUB(displacement_test[0], a, displacement_test[0]);
		SUB(displacement_test[1], b, displacement_test[1]);
	}
	for (unsigned int i = 0; i < displacement_test.size(); ++i) {
		thread->assignTask(this, MOVE_OBJECT2, i);
		if (i < cloth->size()) {
			SUM_(cloth->data()[i].center, displacement_test[i]);
		}
		else if (i < tetrahedron_end_index) {
			SUM_(tetrahedron->data()[i - cloth->size()].center, displacement_test[i]);
		}
		else {
			SUM_(collider->data()[i - tetrahedron_end_index].center, displacement_test[i]);
		}
	}
}

//update both vertex_position & vertex_position_render
void MoveObject::move(unsigned int obj_No, double* displacement, bool only_move_vertex_pos)
{
	this->displacement = displacement;
	thread->assignTask(this, MOVE_OBJECT, obj_No);

	if (obj_No < cloth->size()) {
		SUM_(cloth->data()[obj_No].center, displacement);
	}
	else if (obj_No < tetrahedron_end_index) {
		SUM_(tetrahedron->data()[obj_No - cloth->size()].center, displacement);
	}
	else {
		SUM_(collider->data()[obj_No - tetrahedron_end_index].center, displacement);
	}

	if (!only_move_vertex_pos) {
		if (obj_No < cloth->size()) {
			memcpy(cloth->data()[obj_No].mesh_struct.vertex_for_render[0].data(), cloth->data()[obj_No].mesh_struct.vertex_position[0].data(), 24 * cloth->data()[obj_No].mesh_struct.vertex_position.size());
		}
		else if (obj_No < tetrahedron_end_index) {
			memcpy(tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_for_render[0].data(), tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_position[0].data(),
				24* tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_for_render.size());
		}
		else {
			memcpy(collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_for_render[0].data(), collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_position[0].data(),
				24* collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_for_render.size());
		}
	}
	//SUM_(total_displacement[obj_No], displacement);
	//std::cout << "======" << std::endl;
	//for (unsigned int i = 0; i < total_obj_num; ++i) {
	//	std::cout << total_displacement[i][0] << ", " << total_displacement[i][1] << ", " << total_displacement[i][2] << std::endl;
	//}
}

//MOVE_OBJECT2
void MoveObject::moveDiffInitialCurrent(int thread_No, unsigned int obj_No)
{
	std::array<double, 3>* position;
	unsigned int vertex_start, vertex_end;
	double displacement_[3];
	double displacement_render_[3];
	memcpy(displacement_, displacement_test[obj_No].data(), 24);
	memcpy(displacement_render_, displacement_test_render[obj_No].data(), 24);



	std::array<double, 3>* render_position;
	if (obj_No < cloth->size()) {
		render_position = cloth->data()[obj_No].mesh_struct.vertex_for_render.data();
		position = cloth->data()[obj_No].mesh_struct.vertex_position.data();
		vertex_start = cloth->data()[obj_No].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = cloth->data()[obj_No].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	}
	else if (obj_No < tetrahedron_end_index) {
		render_position = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_for_render.data();
		position = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_position.data();
		vertex_start = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	}
	else {
		render_position = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_for_render.data();
		position = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_position.data();
		vertex_start = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	}
	for (unsigned int i = vertex_start; i < vertex_end; ++i) {
		SUM_(position[i], displacement_);
		SUM_(render_position[i], displacement_render_);

	}

}
//MOVE_OBJECT
void MoveObject::move(int thread_No, unsigned int obj_No)
{
	std::array<double, 3>* position;
	unsigned int vertex_start, vertex_end;
	double displacement_[3];
	memcpy(displacement_, displacement, 24);

	if (obj_No < cloth->size()) {
		position = cloth->data()[obj_No].mesh_struct.vertex_position.data();
		vertex_start = cloth->data()[obj_No].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = cloth->data()[obj_No].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	}
	else if (obj_No < tetrahedron_end_index) {
		position = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_position.data();
		vertex_start = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	}
	else {
		position = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_position.data();
		vertex_start = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	}
	for (unsigned int i = vertex_start; i < vertex_end; ++i) {
		SUM_(position[i], displacement_);
	}

}

void MoveObject::updateRotationMatrix(double angle_move, unsigned int dimension, double* rotation_matrix)
{
	double v[9];
	memset(v, 0, 72);
	v[0] = v[4] = v[8] = 1.0;

	double temp_rotate_matrix[9];//change to row major
	temp_rotate_matrix[0] = rotation_matrix[0];
	temp_rotate_matrix[1] = rotation_matrix[3];
	temp_rotate_matrix[2] = rotation_matrix[6];
	temp_rotate_matrix[3] = rotation_matrix[1];
	temp_rotate_matrix[4] = rotation_matrix[4];
	temp_rotate_matrix[5] = rotation_matrix[7];
	temp_rotate_matrix[6] = rotation_matrix[2];
	temp_rotate_matrix[7] = rotation_matrix[5];
	temp_rotate_matrix[8] = rotation_matrix[8];


	double axis[3];
	switch (dimension)
	{
	case 0:		
		memcpy(axis, temp_rotate_matrix, 24);
		normalize(axis);
		rotateAroundVector(v, axis, angle_move);
		break;
	case 1:		
		memcpy(axis, (temp_rotate_matrix+3), 24);
		normalize(axis);
	//	std::cout << axis[0] << " " << axis[1] << " " << axis[2] << std::endl;
		rotateAroundVector(v, axis, angle_move);
		break;
	case 2:
		memcpy(axis, (temp_rotate_matrix + 6), 24);
		normalize(axis);
		rotateAroundVector(v, axis, angle_move);
		break;
	}

	for (unsigned int i = 0; i < 3; ++i) {
		rotation_matrix[3*i] = DOT(temp_rotate_matrix, (v+3*i));
		rotation_matrix[3*i+1] = DOT((temp_rotate_matrix + 3), (v+3*i));
		rotation_matrix[3*i+2] = DOT((temp_rotate_matrix + 6), (v+3*i));
	}
}


void MoveObject::rotation(double angle_move, unsigned int dimension, bool only_move_vertex_pos)
{
	double* ori_rotation_matrix;

	if (select_object_index < cloth->size()) {
		memcpy(center, cloth->data()[select_object_index].center, 24);
		ori_rotation_matrix = cloth->data()[select_object_index].rotation_matrix;
		
	}
	else if (select_object_index < tetrahedron_end_index) {
		memcpy(center, tetrahedron->data()[select_object_index - cloth->size()].center, 24);
		ori_rotation_matrix = tetrahedron->data()[select_object_index - cloth->size()].rotation_matrix;
	}
	else {
		memcpy(center, collider->data()[select_object_index - tetrahedron_end_index].center, 24);
		ori_rotation_matrix = collider->data()[select_object_index - tetrahedron_end_index].rotation_matrix;
	}
	
	memset(rotation_matrix, 0, 72);
	rotation_matrix[0] = 1.0;
	rotation_matrix[4] = 1.0;
	rotation_matrix[8] = 1.0;

	double axis[3];
	switch (dimension)
	{
	case 0:
		axis[0] = ori_rotation_matrix[0]; axis[1] = ori_rotation_matrix[3]; axis[2] = ori_rotation_matrix[6];
		normalize(axis);
		rotateAroundVector(rotation_matrix, axis, angle_move);
		break;
	case 1:
		axis[0] = ori_rotation_matrix[1]; axis[1] = ori_rotation_matrix[4]; axis[2] = ori_rotation_matrix[7];
		normalize(axis);
		rotateAroundVector(rotation_matrix, axis, angle_move);
		break;
	case 2:
		axis[0] = ori_rotation_matrix[2]; axis[1] = ori_rotation_matrix[5]; axis[2] = ori_rotation_matrix[8];
		normalize(axis);
		rotateAroundVector(rotation_matrix, axis, angle_move);
		break;
	}

	//std::cout << "rotation_matrix " << std::endl;
	//for (unsigned int i = 0; i < 9; ++i) {
	//	std::cout << rotation_matrix[i] << " ";
	//}
	//std::cout << std::endl;

	updateRotationMatrix(angle_move, dimension, ori_rotation_matrix);



	select_dimension = dimension;

	thread->assignTask(this, ROTATE_AROUND_AXIS, select_object_index);
	if (!only_move_vertex_pos) {
		if (select_object_index < cloth->size()) {
			cloth->data()[select_object_index].mesh_struct.vertex_for_render = cloth->data()[select_object_index].mesh_struct.vertex_position;
		}
		else if (select_object_index < tetrahedron_end_index) {
			tetrahedron->data()[select_object_index - cloth->size()].mesh_struct.vertex_for_render = tetrahedron->data()[select_object_index - cloth->size()].mesh_struct.vertex_position;
		}
		else {
			collider->data()[select_object_index - tetrahedron_end_index].mesh_struct.vertex_for_render = collider->data()[select_object_index - tetrahedron_end_index].mesh_struct.vertex_position;
		}
	}

}

//ROTATE_AROUND_AXIS
void MoveObject::rotateAroundAxis(int thread_No, unsigned int obj_No)
{
	double center[3];
	memcpy(center, this->center, 24);


	std::array<double, 3>* position;
	unsigned int vertex_start, vertex_end;

	double rotation_matrix[9];
	memcpy(rotation_matrix, this->rotation_matrix, 72);
	
	if (obj_No < cloth->size()) {
		position = cloth->data()[obj_No].mesh_struct.vertex_position.data();
		vertex_start = cloth->data()[obj_No].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = cloth->data()[obj_No].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	}
	else if (obj_No < tetrahedron_end_index) {
		position = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_position.data();
		vertex_start = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	}
	else {
		position = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_position.data();
		vertex_start = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_index_begin_per_thread[thread_No];
		vertex_end = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_index_begin_per_thread[thread_No + 1];
	}

	double temp_position0, temp_position1, temp_position2;

	for (unsigned int i = vertex_start; i < vertex_end; ++i) {
		position[i][0] -= center[0];
		position[i][1] -= center[1];
		position[i][2] -= center[2];
		temp_position0 = DOT(position[i], rotation_matrix);
		temp_position1 = DOT(position[i], (rotation_matrix+3));
		temp_position2 = DOT(position[i], (rotation_matrix+6));
		position[i].data()[0] = temp_position0 + center[0];
		position[i].data()[1] = temp_position1 + center[1];
		position[i].data()[2] = temp_position2 + center[2];
	}
}


//void MoveObject::moveAroundAxes(int thread_No, int obj_index)
//{
//
//}
//
//
//void MoveObject::moveAroundAxe(int thread_No, int obj_index, unsigned int dimension)
//{
//
//}

//void MoveObject::moveScript(unsigned int type)
//{
//	displacement_test.resize(2);
//	displacement_test_render.resize(2);
//	switch (time_step)
//	{
//	case 0:
//		if (type == 1) {
//			displacement_test[0] = { -0.308742, 0.640283, 0.0271892 };
//			displacement_test[1] = { -0.52193, 0.580719, 0.0493092 };
//		}
//		else if(type==2){
//			displacement_test[0] = { 0, 0, 0 };
//			displacement_test[1] = { -0.239594, -0.0259783, -0.33549 };
//		}
//		else if (type == 3) {
//			displacement_test[0] = { -0.308742, 0.640283, 0.0271892 };
//			displacement_test[1] = { -0.52193, 0.580719, 0.0323092 };
//		}
//		break;
//	case 1:
//		if (type == 1) {		
//			double a[3] = { -0.510271, 0.631315, 0.0321853 };
//			double b[3]= { -0.52193, 0.580719, 0.0493092 };
//			SUB(displacement_test[0], a, displacement_test[0]);
//			SUB(displacement_test[1], b, displacement_test[1]);
//		}
//		else if (type == 2) {
//			double a[3] = { 0, 0, 0 };
//			double b[3] = { -0.0776213, -0.0571768, -0.140402 };
//			SUB(displacement_test[0], a, displacement_test[0]);
//			SUB(displacement_test[1], b, displacement_test[1]);
//		}
//		else if (type == 3) {
//			double a[3] = { -0.510271, 0.631315, 0.0321853 };
//			double b[3] = { -0.52193, 0.580719, 0.0493092 };
//			SUB(displacement_test[0], a, displacement_test[0]);
//			SUB(displacement_test[1], b, displacement_test[1]);
//		}
//		break;
//	}
//	for (unsigned int i = 0; i < displacement_test.size(); ++i) {
//		move(i, displacement_test[i].data());
//	}
//	time_step++;
//}