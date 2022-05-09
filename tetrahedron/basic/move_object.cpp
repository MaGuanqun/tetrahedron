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
	}
}

//update both vertex_position & vertex_position_render
void MoveObject::move(unsigned int obj_No, double* displacement, bool only_move_vertex_pos)
{
	this->displacement = displacement;
	thread->assignTask(this, MOVE_OBJECT, obj_No);

	if (!only_move_vertex_pos) {
		if (obj_No < cloth->size()) {
			cloth->data()[obj_No].mesh_struct.vertex_for_render = cloth->data()[obj_No].mesh_struct.vertex_position;
		}
		else if (obj_No < tetrahedron_end_index) {
			tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_for_render = tetrahedron->data()[obj_No - cloth->size()].mesh_struct.vertex_position;
		}
		else {
			collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_for_render = collider->data()[obj_No - tetrahedron_end_index].mesh_struct.vertex_position;
		}
	}
	SUM_(total_displacement[obj_No], displacement);
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