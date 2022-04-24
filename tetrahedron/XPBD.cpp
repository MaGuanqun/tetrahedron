#include"XPBD.h"



XPBD::XPBD()
{
	gravity_ = 9.8;
	total_thread_num = std::thread::hardware_concurrency();
	sub_step_num = 1;

}


void XPBD::setForPD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider,
	Thread* thread, double* tolerance_ratio)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	sub_time_step = time_step / (double)sub_step_num;
	this->thread = thread;
	collision.initial(cloth, collider, tetrahedron, thread, tolerance_ratio);
	total_obj_num = cloth->size() + tetrahedron->size();
	reorganzieDataOfObjects();
	initialVariable();
}


void XPBD::initialVariable()
{
	f_ext.resize(total_obj_num);
	velocity.resize(total_obj_num);
	total_gravity.resize(total_obj_num);
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		f_ext[i].resize(mesh_struct[i]->vertex_position.size());
		velocity[i].resize(mesh_struct[i]->vertex_position.size());
		total_gravity[i].resize(mesh_struct[i]->vertex_position.size());
	}
}

void XPBD::reorganzieDataOfObjects()
{
	vertex_position.resize(total_obj_num);
	mesh_struct.resize(total_obj_num);
	vertex_index_begin_per_thread.resize(total_obj_num);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_position[i]= cloth->data()[i].mesh_struct.vertex_position.data();
		mesh_struct[i]= &cloth->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i]= cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		mesh_struct[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
	}
}

void XPBD::PBDsolve()
{

}

//SET_POS_PREDICT
void XPBD::setPosPredict(int thread_No)
{
	std::array<double, 3>*vertex_pos;
	unsigned int vertex_end=0;
	double delta_t_2 = sub_time_step * sub_time_step;
	double* mass_inv;
	std::array<double, 3>* f_ext_;
	std::array<double, 3>* velocity_;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		vertex_end = vertex_index_begin_per_thread[i][thread_No + 1];
		mass_inv = mesh_struct[i]->mass.data();
		f_ext_ = f_ext[i].data();
		velocity_ = velocity[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < vertex_end; ++j) {
			for (unsigned int k = 0; k < 3; ++k) {
				vertex_pos[j][k] += sub_time_step * velocity_[j][k] + delta_t_2 * mass_inv[j] * f_ext_[j][k];
			}
		}
	}
}