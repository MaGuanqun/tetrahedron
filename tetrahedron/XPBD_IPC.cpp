#include"XPBD_IPC.h"

XPBD_IPC::XPBD_IPC()
{
	gravity_ = 9.8;
	sub_step_num = 1;
	iteration_number = 300;

	damping_coe = 0.0;

	perform_collision = false;
	max_iteration_number = 100;
	outer_max_iteration_number = 4;
	XPBD_constraint.epsilon_for_bending = 1e-10;

	velocity_damp = 0.995;
	energy_converge_ratio = 1e-3;
}
void XPBD_IPC::initial()
{
	resetExternalForce();
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
	}
}

void XPBD_IPC::reset()
{
	resetExternalForce();
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
	}
}

void XPBD_IPC::updateItrInfo(int* iteration_num)
{
	iteration_num[LOCAL_GLOBAL] = iteration_number;
	//outer_iteration_number = iteration_num[OUTER];
	//sub_step_num = iteration_num[OUTER];
	sub_time_step = time_step / (double)sub_step_num;
}


void XPBD_IPC::setForXPBD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
	Thread* thread, double* tolerance_ratio)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	sub_time_step = time_step / (double)sub_step_num;
	this->thread = thread;
	total_thread_num = thread->thread_num;

	total_obj_num = cloth->size() + tetrahedron->size();
	reorganzieDataOfObjects();
	initialVariable();
	initialClothBending();
	setConstraintIndex();
	//energy_per_thread.resize(thread->thread_num,0.0);
	if (perform_collision) {
		//collision.energy = energy_per_thread.data();
		collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio, XPBD_);
		//collision.setParameter(&lambda_collision,lambda.data()+ constraint_index_start[3], collision_constraint_index_start.data(), damping_coe, sub_time_step);
	}

	setConvergeCondition();

}

void XPBD_IPC::setConvergeCondition()
{
	converge_condition_ratio = 1e-3;
	double edge_length = calEdgeLength();
	max_move_standard = converge_condition_ratio * edge_length;
	outer_max_move_standard = 5.0 * converge_condition_ratio * edge_length;
}

void XPBD_IPC::initialClothBending()
{
	lbo_weight.resize(cloth->size());
	vertex_lbo.resize(cloth->size());
	rest_mean_curvature_norm.resize(cloth->size());
	//rest_Aq.resize(cloth->size());
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		XPBD_constraint.initial_LBO_EdgeCotWeight(cloth->data()[i].mesh_struct, lbo_weight[i], vertex_lbo[i], rest_mean_curvature_norm[i]);
	}
}

void XPBD_IPC::setConstraintIndex()
{
	constraint_index_start.resize(5); //bending, edge_length, ARAP, floor collision
	constraint_index_start[0] = 0;
	unsigned int constraint_number = 0;
	if (use_bending_based_on_vertex) {
		for (unsigned int i = 0; i < cloth->size(); ++i) {
			constraint_number += mesh_struct[i]->vertex_position.size();
		}
	}
	else {
		for (unsigned int i = 0; i < cloth->size(); ++i) {
			constraint_number += mesh_struct[i]->edge_vertices.size() >> 1;
		}
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
	constraint_number = 0;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		constraint_number += mesh_struct[i]->vertex_position.size();
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		constraint_number += tetrahedron->data()[i].mesh_struct.vertex_index_on_sureface.size();
	}
	constraint_index_start[4] = constraint_number + constraint_index_start[3];
	lambda.resize(constraint_index_start[4]);

	constraint_number = 0;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		constraint_number += tetrahedron->data()[i].mesh_struct.vertex_position.size();
	}

	lambda_collision.reserve(constraint_index_start[1] + constraint_number);
	collision_constraint_index_start.resize(3);
	for (unsigned int i = 0; i < collision_constraint_index_start.size(); ++i) {
		collision_constraint_index_start[i].resize(total_thread_num + 1, 0);
	}
}


void XPBD_IPC::initialVariable()
{
	f_ext.resize(total_obj_num);
	velocity.resize(total_obj_num);
	gravity[0] = 0;
	gravity[1] = -gravity_;
	gravity[2] = 0;

	sn.resize(total_obj_num);


	residual.resize(total_obj_num);

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		f_ext[i].resize(mesh_struct[i]->vertex_position.size());
		velocity[i].resize(mesh_struct[i]->vertex_position.size());
		sn[i].resize(mesh_struct[i]->vertex_position.size());
		residual[i].resize(mesh_struct[i]->vertex_position.size());
	}
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
		memset(f_ext[i][0].data(), 0, 24 * f_ext[i].size());
		memset(sn[i][0].data(), 0, 24 * sn[i].size());
	}

}

void XPBD_IPC::reorganzieDataOfObjects()
{
	vertex_position.resize(total_obj_num);
	initial_vertex_position.resize(total_obj_num);
	mesh_struct.resize(total_obj_num);
	vertex_index_begin_per_thread.resize(total_obj_num);
	record_vertex_position.resize(total_obj_num);
	//record_outer_vertex_position.resize(total_obj_num);
	unfixed_vertex.resize(total_obj_num);

	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position.data();
		initial_vertex_position[i] = cloth->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i] = &cloth->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		//record_outer_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		unfixed_vertex[i] = &cloth->data()[i].mesh_struct.unfixed_point_index;
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		initial_vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		//record_outer_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		unfixed_vertex[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.unfixed_point_index;
	}
	recordVertexPosition();

	collider_mesh_struct.resize(collider->size());
	for (unsigned int i = 0; i < collider->size(); ++i) {
		collider_mesh_struct[i] = &collider->data()[i].mesh_struct;
	}

}
