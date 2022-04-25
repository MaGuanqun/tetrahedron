#include"XPBD.h"



XPBD::XPBD()
{
	gravity_ = 9.8;
	total_thread_num = std::thread::hardware_concurrency();
	sub_step_num = 1;
	iteration_number = 50;
}


void XPBD::setForXPBD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider,
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
	initialClothBending();
	setConstraintIndex();
}


void XPBD::initialClothBending()
{
	lbo_weight.resize(cloth->size());
	vertex_lbo.resize(cloth->size());
	rest_mean_curvature_norm.resize(cloth->size());
	for (unsigned int i = 0; i < cloth->size(); ++i){
		XPBD_constraint.initial_LBO_EdgeCotWeight(cloth->data()[i].mesh_struct, lbo_weight[i], vertex_lbo[i], rest_mean_curvature_norm[i]);
	}

}

void XPBD::setConstraintIndex()
{
	constraint_index_start.resize(5); //bending, edge_length, ARAP, collision
	constraint_index_start[0] = 0;
	unsigned int constraint_number = 0;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		constraint_number += mesh_struct[i]->vertex_position.size();
	}
	constraint_index_start[1] = constraint_number;
	constraint_number = 0;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		constraint_number += mesh_struct[i]->edges.size();
	}
	constraint_index_start[2] = constraint_number + constraint_index_start[1];
	constraint_number = 0;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		constraint_number += tetrahedron->data()[i].mesh_struct.indices.size();
	}
	constraint_index_start[3] = constraint_number + constraint_index_start[2];
	lambda.resize(constraint_index_start[3]);

	constraint_number = 0;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		constraint_number += tetrahedron->data()[i].mesh_struct.vertex_position.size();
	}

	lambda_collision.reserve(constraint_index_start[1] + constraint_number);

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
	thread->assignTask(this, SET_POS_PREDICT);
	for (unsigned int i = 0; i < iteration_number; ++i) {

	}
	
}

void XPBD::solveConstraint()
{
	solveBendingConstraint();
	solveEdgeLengthConstraint();

}


void XPBD::solveEdgeLengthConstraint()
{
	unsigned int size;
	MeshStruct* mesh_struct_;
	std::array<double, 3>* vertex_pos;
	double stiffness;
	unsigned int* edge_vertex_index;
	double* mass_inv;
	double* lambda_ = lambda.data()+ constraint_index_start[1];
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		size = mesh_struct_->edges.size();
		mesh_struct_ = mesh_struct[i];
		vertex_pos = vertex_position[i];
		stiffness = cloth->data()[i].length_stiffness[0];
		edge_vertex_index = mesh_struct_->edge_vertices.data();
		mass_inv = mesh_struct_->mass_inv.data();
		for (unsigned int j = 0; j < size; ++j) {

			XPBD_constraint.solveEdgeLengthConstraint(vertex_pos[edge_vertex_index[j << 1]].data(),
				vertex_pos[edge_vertex_index[(j << 1) + 1]].data(), mesh_struct_->edge_length[j], stiffness, sub_time_step, mass_inv[edge_vertex_index[j << 1]],
				mass_inv[edge_vertex_index[(j << 1) + 1]], *lambda_);
			lambda_++;
		}
	}
}


void XPBD::solveBendingConstraint()
{
	unsigned int size;
	double* mass_inv;
	double* rest_mean_curvature_norm_;
	double* lbo_weight_;
	VectorXd* vertex_lbo_;
	MeshStruct* mesh_struct_;
	std::array<double, 3>*vertex_pos;
	double stiffness;
	double* lambda_ = lambda.data();
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		stiffness = cloth->data()[i].bend_stiffness;
		mesh_struct_ = mesh_struct[i];
		mass_inv = mesh_struct_->mass_inv.data();
		size = mesh_struct_->vertex_position.size();
		rest_mean_curvature_norm_ = rest_mean_curvature_norm[i].data();
		lbo_weight_ = lbo_weight[i].data();
		vertex_lbo_ = vertex_lbo[i].data();
		vertex_pos = vertex_position[i];
		for (unsigned int j = 0; j < size; ++j) {
			XPBD_constraint.solveBendingConstraint(vertex_pos[j].data(), mass_inv[j], vertex_pos, mesh_struct_->vertices[j].neighbor_vertex,
				rest_mean_curvature_norm_[j], lbo_weight_[j], vertex_lbo_[j], stiffness, sub_time_step, mass_inv, *lambda_);
			lambda_++;
		}		
	}
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