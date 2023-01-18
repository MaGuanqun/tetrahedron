#include"XPBD_IPC.h"
#include<algorithm>
#include"XPBD/FEM_relate.h"
#include"tet_inversion.h"

XPBD_IPC::XPBD_IPC()
{
	gravity_ =9.8;
	sub_step_num = 1;
	iteration_number = 300;

	damping_coe = 0.0;

	perform_collision = true;

	XPBD_constraint.epsilon_for_bending = 1e-10;

	velocity_damp = 0.995;
	energy_converge_ratio = 1e-3;

	min_inner_iteration = 4;
	min_outer_iteration = 2;

	max_move_standard = 1e-3;
	max_move_standard_inner_itr = 1e-4;
	max_iteration_number = 100;
	outer_max_iteration_number = 40;
	energy_converge_standard = 1e-10;

	second_order_constraint.solve_exact_ARAP_hessian = false;
	min_collision_time = 1e-3;
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
	iteration_num[OUTER] = outer_itr_num;
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
	has_collider = !collider->empty();
	this->floor = floor;

	reorganzieDataOfObjects();
	initialVariable();
	initialClothBending();
	setConstraintIndex();
	//energy_per_thread.resize(thread->thread_num,0.0);
	if (perform_collision) {
		collision.inner_iteration_number = &inner_iteration_number;
		//collision.energy = energy_per_thread.data();
		collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio, XPBD_IPC_,false);
		collision.setCollisionFreeVertex(&record_collision_free_vertex_position_address, &record_vertex_position);
		collision.inner_itr_num_standard = &inner_itr_num_standard;
		collision.min_collision_time = &min_collision_time;


		//collision_compare.inner_iteration_number = &inner_iteration_number;
		//collision_compare.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio, XPBD_IPC_, false);
		//collision_compare.setCollisionFreeVertex(&record_collision_free_vertex_position_address, &record_vertex_position);
		//collision_compare.inner_itr_num_standard = &inner_itr_num_standard;
		//collision_compare.min_collision_time = &min_collision_time;
		//collision.outer_itr_num = &outer_itr_num;
		//collision_compare.outer_itr_num = &outer_itr_num;

		//collision.setParameter(&lambda_collision,lambda.data()+ constraint_index_start[3], collision_constraint_index_start.data(), damping_coe, sub_time_step);
	}

	initalARAPHessianStorages();
	setColorNum();



	initialHessianMap();

	if (!tetrahedron->empty()) {
		color_group_num = tet_color_groups[0]->size();
	}
	else {
		color_group_num = 3;
	}
	

	//thread->assignTask(this, TEST_MULTI);
	//std::cout << "test count " << counter << std::endl;
	//test();
}


void XPBD_IPC::setColorNum()
{
	int color_num = 1;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		if (color_num < tetrahedron->data()[i].mesh_struct.tet_color_group[0].size()) {
			color_num = tetrahedron->data()[i].mesh_struct.tet_color_group[0].size();
		}
	}
	max_tet_color_num = color_num;
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
	grad_max_store.resize(total_thread_num);
	max_dis_record.resize(total_thread_num);
}

void XPBD_IPC::reorganzieDataOfObjects()
{
	energy_per_thread.resize(total_thread_num);

	vertex_position.resize(total_obj_num);
	vertex_position_render.resize(total_obj_num);
	mesh_struct.resize(total_obj_num);
	vertex_index_begin_per_thread.resize(total_obj_num);
	tet_index_begin_per_thread.resize(total_obj_num);
	record_vertex_position.resize(total_obj_num);
	record_vertex_position_address.resize(total_obj_num);
	record_collision_free_vertex_position_address.resize(total_obj_num);
	record_collision_free_vertex_position.resize(total_obj_num);

	//record_outer_vertex_position.resize(total_obj_num);
	//unfixed_vertex.resize(total_obj_num);
	triangle_indices.resize(total_obj_num);
	tet_indices.resize(total_obj_num);
	tet_A.resize(total_obj_num);
	tet_volume.resize(total_obj_num);
	edge_vertices.resize(total_obj_num);
	is_vertex_fixed.resize(total_obj_num);

	vertex_index_surface.resize(total_obj_num);
	mass.resize(total_obj_num);


	triangle_around_triangle.resize(total_obj_num);
	edge_around_triangle.resize(total_obj_num);
	vertices.resize(total_obj_num);
	tet_around_vertex.resize(total_obj_num);
	tet_around_triangle.resize(total_obj_num);
	vertex_surface_to_global.resize(total_obj_num);

	triangle_around_edge.resize(total_obj_num);
	edge_around_edge.resize(total_obj_num);
	tet_around_edge.resize(total_obj_num);


	rest_edge_length.resize(total_obj_num);
	vertex_index_of_a_tet_color_group.resize(total_obj_num);

	vertex_index_of_a_tet_color_per_thread_start_group.resize(total_obj_num);

	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position.data();
		vertex_position_render[i] = cloth->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i] = &cloth->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i] = cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		record_vertex_position_address[i]= record_vertex_position[i].data();
		record_collision_free_vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position;
		record_collision_free_vertex_position_address[i] = record_collision_free_vertex_position[i].data();
		triangle_indices[i] = cloth->data()[i].mesh_struct.triangle_indices.data();
		//record_outer_vertex_position[i].resize(cloth->data()[i].mesh_struct.vertex_position.size());
		//unfixed_vertex[i] = &cloth->data()[i].mesh_struct.unfixed_point_index;
		edge_vertices[i] = cloth->data()[i].mesh_struct.edge_vertices.data();

		is_vertex_fixed[i] = &cloth->data()[i].mesh_struct.is_vertex_fixed;
		vertex_index_surface[i] = cloth->data()[i].mesh_struct.vertex_surface_index.data();
		mass[i] = cloth->data()[i].mesh_struct.mass.data();

		triangle_around_triangle[i] = cloth->data()[i].mesh_struct.face_around_face.data();
		edge_around_triangle[i] = cloth->data()[i].mesh_struct.edge_around_face.data();
		vertices[i] = cloth->data()[i].mesh_struct.vertices.data();

		triangle_around_edge[i] = cloth->data()[i].mesh_struct.face_around_edge.data();
		edge_around_edge[i] = cloth->data()[i].mesh_struct.edge_around_edge.data();

		rest_edge_length[i] = cloth->data()[i].mesh_struct.edge_length.data();

	}

	unfix_tet_index.resize(tetrahedron->size());
	unfixed_tet_vertex_num.resize(tetrahedron->size());

	tet_color_groups.resize(tetrahedron->size());
	//tet_color_groups_label.resize(tetrahedron->size());

	tet_around_tet.resize(tetrahedron->size());

	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		vertex_position_render[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		record_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		record_vertex_position_address[i + cloth->size()] = record_vertex_position[i + cloth->size()].data();
		record_collision_free_vertex_position[i + cloth->size()]=tetrahedron->data()[i].mesh_struct.vertex_position;
		record_collision_free_vertex_position_address[i + cloth->size()] = record_collision_free_vertex_position[i + cloth->size()].data();
		triangle_indices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.triangle_indices.data();
		tet_indices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.indices.data();
		tet_A[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.A.data();
		tet_volume[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.volume.data();
		//record_outer_vertex_position[i + cloth->size()].resize(tetrahedron->data()[i].mesh_struct.vertex_position.size());
		//unfixed_vertex[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.unfixed_point_index;
		edge_vertices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_vertices.data();

		is_vertex_fixed[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.is_vertex_fixed;
		vertex_index_surface[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
		mass[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.mass.data();

		triangle_around_triangle[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.face_around_face.data();
		edge_around_triangle[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_around_face.data();
		vertices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertices.data();
		tet_around_vertex[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_tet_index.data();
		tet_around_triangle[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tet_around_face.data();
		vertex_surface_to_global[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_on_sureface.data();

		triangle_around_edge[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.face_around_edge.data();
		edge_around_edge[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_around_edge.data();
		tet_around_edge[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tet_around_edge.data();
		rest_edge_length[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_length.data();

		
		unfix_tet_index[i] = tetrahedron->data()[i].mesh_struct.unfixied_indices.data();
		unfixed_tet_vertex_num[i] = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();

		tet_color_groups[i] =&tetrahedron->data()[i].mesh_struct.tet_color_group;
		//tet_color_groups_label[i] =tetrahedron->data()[i].mesh_struct.tet_in_collision.data();

		tet_around_tet[i] = tetrahedron->data()[i].mesh_struct.tet_tet_index.data();
		vertex_index_of_a_tet_color_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_of_a_tet_color_group.data();

		tet_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tetrahedron_index_begin_per_thread.data();

		vertex_index_of_a_tet_color_per_thread_start_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_of_a_tet_color_per_thread_start_group.data();
	}
	if (!collider->empty()) {
		triangle_indices_collider.resize(collider->size());
		collider_mesh_struct.resize(collider->size());
		collider_edge_vertices.resize(collider->size());
		vertex_position_collider.resize(collider->size());
		rest_edge_length_collider.resize(collider->size());
		


		for (unsigned int i = 0; i < collider->size(); ++i) {
			collider_mesh_struct[i] = &collider->data()[i].mesh_struct;
			triangle_indices_collider[i] = collider->data()[i].mesh_struct.triangle_indices.data();
			collider_edge_vertices[i] = collider->data()[i].mesh_struct.edge_vertices.data();
			vertex_position_collider[i] = collider->data()[i].mesh_struct.vertex_position.data();
			rest_edge_length_collider[i] = collider->data()[i].mesh_struct.edge_length.data();
		}
	}
	


	int total_vertex_num = 0;
	vertex_num_on_surface_prefix_sum.resize(total_obj_num+1);
	vertex_index_prefix_sum_obj.resize(total_obj_num + 1, 0);

	vertex_num_on_surface_prefix_sum[0] = 0;
	for (int i = 1; i <= total_obj_num; ++i) {
		if (i <= cloth->size()) {
			vertex_num_on_surface_prefix_sum[i] = vertex_num_on_surface_prefix_sum[i - 1] + vertex_index_begin_per_thread[i - 1][total_thread_num];
		}
		else {
			vertex_num_on_surface_prefix_sum[i] = vertex_num_on_surface_prefix_sum[i - 1] + tetrahedron->data()[i-1-cloth->size()].mesh_struct.vertex_index_on_sureface.size();
		}

		vertex_index_prefix_sum_obj[i] = vertex_index_prefix_sum_obj[i-1] + mesh_struct[i-1]->vertex_for_render.size();
	}


	global_vertex_index_start_per_thread.resize(total_thread_num + 1,0);
	arrangeIndex(total_thread_num, vertex_index_prefix_sum_obj[total_obj_num], global_vertex_index_start_per_thread.data());


	int vertex_num;
	record_vertex_position_every_thread.resize(total_obj_num);
	record_vertex_position_num_every_thread.resize(total_obj_num);

	record_vertex_by_thread.resize(total_thread_num);
	record_vertex_update_num_by_thread.resize(total_thread_num);

	for (int i = 0; i < total_obj_num; ++i) {
		vertex_num = mesh_struct[i]->vertex_position.size();
		record_vertex_position_every_thread[i].resize(total_thread_num);
		record_vertex_position_num_every_thread[i].resize(total_thread_num);
		for (int j = 0; j < total_thread_num; ++j) {
			record_vertex_position_every_thread[i][j].resize(vertex_num);
			record_vertex_position_num_every_thread[i][j].resize(vertex_num);
		}
	}


	for (int i = 0; i < total_thread_num; ++i) {
		record_vertex_by_thread[i].resize(total_obj_num);
		record_vertex_update_num_by_thread[i].resize(total_obj_num);
		for (int j = 0; j < total_obj_num; ++j) {
			record_vertex_by_thread[i][j] = record_vertex_position_every_thread[j][i].data();
			record_vertex_update_num_by_thread[i][j] = record_vertex_position_num_every_thread[j][i].data();
		}

	}

	collision.indicate_if_involved_in_last_color = record_vertex_position_num_every_thread.data();
	//collision_compare.indicate_if_involved_in_last_color = record_vertex_position_num_every_thread.data();
	std::cout << "xpbd ipc reorganzieDataOfObjects " << vertex_num_on_surface_prefix_sum[1] << std::endl;
}


void XPBD_IPC::initialRecordPositionForThread()
{
	for (int i = 0; i < total_obj_num; ++i) {
		for (int j = 0; j < total_thread_num; ++j) {
			memset(record_vertex_position_num_every_thread[i][j].data(), 0, record_vertex_position_num_every_thread[i][j].size() << 2);
		}
	}
}



void XPBD_IPC::initialHessianMap()
{
	

	

	common_hessian.reserve(4 * vertex_num_on_surface_prefix_sum.back());
	

	common_grad.resize(total_thread_num);
	for (int i = 0; i < total_thread_num; ++i) {
		common_grad[i].resize(3 * vertex_index_prefix_sum_obj.back());
	}

	if (perform_collision) {
		collision.common_hessian = &common_hessian;
		for (int i = 0; i < total_thread_num; ++i) {
			collision.common_grad[i] = common_grad[i].data();
		}
		collision.floor_hessian = floor_hessian.data();
	}

	floor_hessian.resize(vertex_index_prefix_sum_obj.back());


	//collision_compare.common_hessian = &common_hessian_compare;
	//common_hessian_compare.reserve(4 * vertex_num_on_surface_prefix_sum.back());
	//common_grad_compare.resize(3 * vertex_index_prefix_sum_obj.back());
	//collision_compare.common_grad = common_grad_compare.data();
	//collision_compare.floor_hessian = &compare_floor_hessian;
	//compare_floor_hessian.reserve(vertex_num_on_surface_prefix_sum.back() / 10);


}


void XPBD_IPC::initalARAPHessianStorages()
{
	prefix_sum_of_every_tet_index.resize(tetrahedron->size()+1);
	prefix_sum_of_every_tet_index[0] = 0;
	for (int i = 0; i < tetrahedron->size(); ++i) {
		prefix_sum_of_every_tet_index[i + 1] = prefix_sum_of_every_tet_index[i] + tetrahedron->data()[i].mesh_struct.indices.size();
	}

	tet_hessian.resize(total_obj_num);
	for (int i = 0; i < tetrahedron->size(); ++i) {
		tet_hessian[i + cloth->size()].reserve(tetrahedron->data()[i].mesh_struct.indices.size() * 4);
	}
	//store_tet_arap_hessian.resize(16 * prefix_sum_of_every_tet_index[tetrahedron->size()],0.0);
	//store_tet_arap_grad.resize(12 * prefix_sum_of_every_tet_index[tetrahedron->size()],0.0);

	is_tet_arap_grad_compute = new std::atomic_flag[prefix_sum_of_every_tet_index[tetrahedron->size()]];
	memset(is_tet_arap_grad_compute, 0, prefix_sum_of_every_tet_index[tetrahedron->size()]<<2);

	max_tet_size_of_a_color_group = prefix_sum_of_every_tet_index[tetrahedron->size()];
	//auto t0 = std::chrono::system_clock::now();

	//for (int i = 0; i < total_thread_num; ++i) {
	//	tetHessian(i);
	//}

	//auto t1 = std::chrono::system_clock::now();

	//thread->assignTask(this, UPDATE_TET_HESSIAN);
	tetHessian();
	//auto t2 = std::chrono::system_clock::now();

	//double time_0 = double(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den * 1000.0;
	//double time_1 = double(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den * 1000.0;

	//std::cout << "time " << time_0 << " " << time_1 << std::endl;

}
		
void XPBD_IPC::saveScene(double* force_direction, int obj_No, bool have_force)
{
	std::string file_name;
	save_scene.save_scene_XPBD(*time_stamp, *time_indicate_for_simu, mesh_struct, &velocity,
		collider_mesh_struct, file_name);
	if (have_force) {
		if (obj_No < cloth->size()) {
			save_scene.save_force(file_name, force_direction, cloth->data()[obj_No].coe_neighbor_vertex_force,
				cloth->data()[obj_No].neighbor_vertex, obj_No);
		}
		else {
			save_scene.save_force(file_name, force_direction, tetrahedron->data()[obj_No - cloth->size()].coe_neighbor_vertex_force,
				tetrahedron->data()[obj_No - cloth->size()].neighbor_vertex, obj_No);
		}
	}

}

void XPBD_IPC::readScene(const char* file_name, double* force_direction, int& obj_No)
{
	std::vector<double>force_coe;
	std::vector<int>neighbor_vertex_index;
	save_scene.read_scene_XPBD(file_name, time_stamp, time_indicate_for_simu, mesh_struct,
		&velocity, collider_mesh_struct, *has_force, force_direction, force_coe, neighbor_vertex_index, obj_No);
	for (unsigned int i = 0; i < mesh_struct.size(); ++i) {
		memcpy(mesh_struct[i]->vertex_for_render[0].data(), mesh_struct[i]->vertex_position[0].data(), 24 * mesh_struct[i]->vertex_position.size());
	}
	for (unsigned int i = 0; i < collider_mesh_struct.size(); ++i) {
		memcpy(collider_mesh_struct[i]->vertex_for_render[0].data(), collider_mesh_struct[i]->vertex_position[0].data(), 24 * collider_mesh_struct[i]->vertex_position.size());
	}
	if (*has_force) {
		if (obj_No < cloth->size()) {
			cloth->data()[obj_No].coe_neighbor_vertex_force = force_coe;
			cloth->data()[obj_No].neighbor_vertex = neighbor_vertex_index;
		}
		else {
			tetrahedron->data()[obj_No - cloth->size()].coe_neighbor_vertex_force = force_coe;
			tetrahedron->data()[obj_No - cloth->size()].neighbor_vertex = neighbor_vertex_index;
		}
	}
	updateRenderNormal();
	updateNormal();
	updateRenderVertexNormal();
}

double XPBD_IPC::calEdgeLength()
{
	double max_length = 0;
	double ave_length = 0;
	int edge_num = 0;
	double* edge_length;
	unsigned int edge_size;
	for (int i = 0; i < total_obj_num; ++i) {
		edge_size = mesh_struct[i]->edge_length.size();
		edge_length = mesh_struct[i]->edge_length.data();
		for (int j = 0; j < edge_size; ++j) {
			if (max_length < edge_length[j]) {
				max_length = edge_length[j];
			}
			ave_length += edge_length[j];
		}
		edge_num += edge_size;
	}
	return ave_length / (double)edge_num;
}

void XPBD_IPC::recordInitialPosition()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(vertex_position_render[i][0].data(), vertex_position[i][0].data(), 24 * record_vertex_position[i].size());
	}
}

void XPBD_IPC::updateCollisionFreePosition()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(record_vertex_position[i][0].data(), vertex_position[i][0].data(), 24 * record_vertex_position[i].size());
		memcpy(record_collision_free_vertex_position[i][0].data(), vertex_position[i][0].data(), 24 * record_vertex_position[i].size());
	}
}


//void XPBD_IPC::tempSavePos()
//{
//	temp_save_pos.resize(total_obj_num);
//	for (int i = 0; i < total_obj_num; ++i) {
//		temp_save_pos[i].resize(record_vertex_position[i].size());
//		memcpy(temp_save_pos[i][0].data(), vertex_position[i][0].data(), 24 * record_vertex_position[i].size());
//	}
//}
//
//void XPBD_IPC::tempRestorePos()
//{
//	for (int i = 0; i < total_obj_num; ++i) {
//		memcpy(vertex_position[i][0].data(), temp_save_pos[i][0].data(),  24 * record_vertex_position[i].size());
//	}
//}

void XPBD_IPC::initialCollisionConstriantNum()
{
	lambda_collision.resize(collision.collisionConstraintNumber(collision_constraint_index_start[0].data(), collision_constraint_index_start[1].data(), collision_constraint_index_start[2].data()));
}


void XPBD_IPC::XPBD_IPC_Position_Solve()
{
	updateCollisionFreePosition();
	thread->assignTask(this, SET_POS_PREDICT_);
	updateSn();
	firstNewtonCD();
	iteration_number = 0;
	if (perform_collision) {
		collision.collisionCulling();
	}
	outer_itr_num = 0;

	while (!convergeCondition(outer_itr_num)) {

		if (perform_collision) {
			collision.globalCollisionTime();
			thread->assignTask(this, COLLISION_FREE_POSITION_);
			updateCollisionFreePosition();
			collision.findClosePair();
		}		
		nearly_not_move = false;

		while (!innerConvergeCondition(inner_iteration_number))
		{
			nearly_not_move = true;
			newtonCDTetWithCollision();
			inner_iteration_number++;
		}
		outer_itr_num++;
		//std::cout << inner_iteration_number << std::endl;
		iteration_number += inner_iteration_number;

	}


	collision.globalCollisionTime();
	thread->assignTask(this, COLLISION_FREE_POSITION_);

	thread->assignTask(this, XPBD_IPC_VELOCITY);
	updatePosition();
	updateRenderNormal();

	updateRenderVertexNormal();
}


void XPBD_IPC::initialRecordHessian()
{
	if (floor->exist) {
		memset(floor_hessian.data(), 0, 8 * floor_hessian.size());
	}
	for (int i = 0; i < total_thread_num; ++i) {
		memset(common_grad[i].data(), 0, 8 * common_grad[i].size());
	}


	//common_hessian_compare.clear();
	//compare_floor_hessian.clear();
	//memset(common_grad_compare.data(), 0, 8 * common_grad_compare.size());
	
}



void XPBD_IPC::XPBD_IPC_Block_Solve_Multithread()
{
	dis_record.clear();
	updateCollisionFreePosition();
	thread->assignTask(this, SET_POS_PREDICT_);
	updateSn();
	iteration_number = 0;
	if (perform_collision) {
		collision.initialPairRecordInfo();
	}
	outer_itr_num = 0;
	displacement_satisfied = false;
	inner_itr_num_standard = min_inner_iteration;
	inner_iteration_number = 0;

	if (perform_collision) {
		collision.collisionCulling();
		collision.collisionTimeWithPair();
		thread->assignTask(this, INVERSION_TEST);
		for (int i = 0; i < total_thread_num; ++i) {
			if (collision.collision_time_thread[i] < collision.collision_time) {
				collision.collision_time = collision.collision_time_thread[i];
			}
		}
		thread->assignTask(this, COLLISION_FREE_POSITION_);
	}
	updateCollisionFreePosition();
	if (perform_collision) {
		warmStart();
	}

	//std::cout << "time stamp " << *time_stamp << " outer itr num " << outer_itr_num << " " << inner_iteration_number << std::endl;
	//max_displacement = 0.0;

	while (!convergeCondition(outer_itr_num)) {
		if (outer_itr_num > 0) {
			updateCollisionFreePosition();		}
		//for (int i = 0; i < color_group_num; ++i) {
			solveNewtonCD_tetBlock();
			inner_iteration_number++;
		//}		
		if (perform_collision) {
			collision.collisionCulling();
			collision.collisionTimeWithPair();			
			thread->assignTask(this, INVERSION_TEST);
			for (int i = 0; i < total_thread_num; ++i) {
				if (collision.collision_time_thread[i] < collision.collision_time) {
					collision.collision_time = collision.collision_time_thread[i];
				}
			}
			thread->assignTask(this, COLLISION_FREE_POSITION_);
		}		
		//if (outer_itr_num > (int)min_outer_iteration - 2)
		//{
		//	computeGradient();
		//}
		outer_itr_num++;
	}

	iteration_number += inner_iteration_number;
	thread->assignTask(this, XPBD_IPC_VELOCITY);
	updatePosition();
	updateRenderNormal();

	updateRenderVertexNormal();
	//vertex_trace.push_back(vertex_position[0][3556]);
	//std::cout <<"time stamp "<< * time_stamp << std::endl;
}



void XPBD_IPC::XPBD_IPC_Block_Solve()
{
	//std::cout << "====" << std::endl;

	record_energy.clear();
	recordInitialPosition();
	updateCollisionFreePosition();
	thread->assignTask(this, SET_POS_PREDICT_);
	updateSn();
	iteration_number = 0;
	if (perform_collision) {
		collision.collisionCulling();
	}
	outer_itr_num = 0;
	displacement_satisfied = false;
	while (!convergeCondition(outer_itr_num)) {
		if (perform_collision) {
			collision.collisionCulling();
			collision.globalCollisionTime();
			thread->assignTask(this, COLLISION_FREE_POSITION_);
			updateCollisionFreePosition();
			collision.findClosePair();

			//if (outer_itr_num == 0) {
			//	firstOnlyInertialCollision();
			//}			
		}
		inner_iteration_number = 0;
		nearly_not_move = false;
		previous_energy = energy;
		while (!innerConvergeCondition(inner_iteration_number))
		{
			previous_energy = energy;
			nearly_not_move = true;
			newtonCDBlock();

			if (inner_iteration_number >= min_inner_iteration - 1) {
				computeCurrentEnergy();
			}

			inner_iteration_number++;
			std::cout << "finish one itr " << inner_iteration_number<<" "<< energy << std::endl;

		}

		outer_itr_num++;
		iteration_number += inner_iteration_number;
	}

	if (perform_collision) {
		collision.globalCollisionTime();
		thread->assignTask(this, COLLISION_FREE_POSITION_);
	}
	thread->assignTask(this, XPBD_IPC_VELOCITY);
	updatePosition();
	updateRenderNormal();

	updateRenderVertexNormal();



	//for (int i = 0; i < record_energy.size(); ++i) {
	//	std::cout << record_energy[i] << std::endl;
	//}
}

void XPBD_IPC::XPBD_IPCSolve()
{
	updateCollisionFreePosition();
	thread->assignTask(this, SET_POS_PREDICT_);

	updateSn();
	firstNewtonCD();
	iteration_number = 0;
	if (perform_collision) {
		collision.collisionCulling();
	}
	outer_itr_num = 0;
	computeCurrentEnergy();

	memset(lambda.data(), 0, 8 * lambda.size());
	previous_energy = energy;
	//while (!convergeCondition(outer_itr_num)) {

	//	if (perform_collision) {
	//		collision.globalCollisionTime();
	//		thread->assignTask(this, COLLISION_FREE_POSITION_);
	//		updateCollisionFreePosition();
	//		collision.findClosePair();
	//		collision.saveCollisionPairVolume();
	//		firstOnlyInertialCollision();
	//	}
	//	inner_iteration_number = 0;
	//	nearly_not_move = false;	
		
		while (!innerConvergeCondition(inner_iteration_number))
		{
			previous_energy = energy;
			nearly_not_move = true;
			newtonCDTetWithCollision();
			inner_iteration_number++;
		}		
		outer_itr_num++;
		//std::cout << inner_iteration_number << std::endl;
		iteration_number += inner_iteration_number;
	//}

	if (perform_collision) {
		collision.globalCollisionTime();
		thread->assignTask(this, COLLISION_FREE_POSITION_);
	}
	thread->assignTask(this, XPBD_IPC_VELOCITY);
	updatePosition();
	updateRenderNormal();
	
	updateRenderVertexNormal();
}

bool XPBD_IPC::innerConvergeCondition(unsigned int iteration_num)
{
	if (iteration_num >= color_group_num){    // inner_itr_num_standard) {//max_iteration_number
		return true;
	}
	//if (iteration_num < 1) {
	//	return false;
	//}
	//if (abs(energy - previous_energy) < energy_converge_standard) {
	//	return true;
	//}
	//if (abs(energy - previous_energy) / previous_energy < energy_converge_ratio) {
	//	return true;
	//}

	return false;

}


bool XPBD_IPC::convCondition(unsigned int iteration_num, unsigned int min_itr, double energy, double previous_energy, unsigned int max_itr, double energy_converge_standard,
	double energy_converge_ratio)
{
	if (iteration_num < min_itr) {//max_iteration_number
		return false;
	}

	if (abs(energy - previous_energy) < energy_converge_standard) {
		return true;
	}
	if (abs(energy - previous_energy) / previous_energy < energy_converge_ratio) {
		return true;
	}
	if (iteration_num > max_itr) {
		return true;
	}
	return false;
}

bool XPBD_IPC::checkMaxDisplacement()
{
	unsigned int num;
	std::array<double, 3>* current_pos;
	std::array<double, 3>* previous_pos;
	double* mass_;


	double max = 0.0;
	int vertex_ = -1;
	int obj_ = -1;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		previous_pos = record_collision_free_vertex_position[i].data();
		current_pos = vertex_position[i];
		num = mesh_struct[i]->vertex_position.size();
		mass_ = mesh_struct[i]->mass_inv.data();
		for (unsigned int j = 0; j < num; ++j) {
			if (mass_[j] != 0.0) {
				for (unsigned int k = 0; k < 3; ++k) {

					if (max < abs(previous_pos[j][k] - current_pos[j][k])) {
						max = abs(previous_pos[j][k] - current_pos[j][k]);
						vertex_ = j;
						obj_ = i;
					}

					//max = (std::max)(abs(previous_pos[j][k] - current_pos[j][k]), max);
					//if (abs(previous_pos[j][k] - current_pos[j][k]) > max_move_standard) {
						//std::cout << outer_itr_num<<" "<< abs(previous_pos[j][k] - current_pos[j][k]) << std::endl;
						//return false;
					//}
				}
			}
		}
	}

	//if (record_max_displacement.size() > 4) {
	//	record_max_displacement.clear();
	//	record_max_displace_vertex.clear();
	//}
	//std::cout << "outer itr num " << outer_itr_num<< " max in checkMaxDisplacement "<< vertex_ <<" " << max << std::endl;

	if (max > max_move_standard) {

		//for (int i = 0; i < record_max_displacement.size(); ++i) {
		//	if (abs(record_max_displacement[i]- max)<1e-5 && record_max_displace_vertex[i << 1] == obj_ && record_max_displace_vertex[i + i + 1] == vertex_) {
		//		return true;
		//	}
		//}

		//record_max_displacement.emplace_back(max);
		//record_max_displace_vertex.emplace_back(obj_);
		//record_max_displace_vertex.emplace_back(vertex_);

		return false;
	}

	//record_max_displacement.clear();
	//record_max_displace_vertex.clear();
	return true;
}

//UPDATE_TET_HESSIAN 
void 	XPBD_IPC::tetHessian()
{
	if (tetrahedron->empty()) {
		return;
	}
	unsigned int tet_around_tet_color_group_end;
	unsigned int tet_around_a_group_start;

	MatrixXd Hessian;

	Hessian.resize(4, 4);
	double stiffness = 0.0;
	Matrix<double, 3, 4>* A;
	double* volume;

	unsigned int prefix_sum_start;

	int k = 0;

	std::unordered_map<std::array<int, 2>, double, pair_hash>* hessian;

	std::array<int, 4>*tet_vertex_indices;

	auto address = tet_hessian[cloth->size()].end();
	int* tet_vertex;
	for (int i = 0; i < tetrahedron->size(); ++i) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		prefix_sum_start = prefix_sum_of_every_tet_index[i];

		A = tetrahedron->data()[i].mesh_struct.A.data();

		volume = tet_volume[i + cloth->size()];
		//tet_around_a_group_start = tetrahedron->data()[i].mesh_struct.tetrahedron_index_begin_per_thread[thread_No];
		tet_around_tet_color_group_end = tetrahedron->data()[i].mesh_struct.indices.size();
		hessian = &tet_hessian[i + cloth->size()];

		tet_vertex_indices = tet_indices[i + cloth->size()];

		for (auto j = 0; j < tet_around_tet_color_group_end; ++j) {
			second_order_constraint.setARAPHessian(Hessian, stiffness, A[j],// ,
				volume[j]); //,

			tet_vertex = tet_vertex_indices[j].data();
			for (int k = 0; k < 4; ++k) {
				for (int m = 0; m < 4; ++m) {
					address = hessian->find(std::array{tet_vertex[k], tet_vertex[m] });
					if (address != hessian->end()) {
						address->second += Hessian.data()[(m << 2) + k];
					}
					else {
						hessian->emplace(std::array{ tet_vertex[k], tet_vertex[m] }, Hessian.data()[(m << 2) + k]);
					}
				}
			}
		}
	}
}


void XPBD_IPC::computeInversionForWarmStart()
{
	int obj_No, size;

	std::array<int, 4>* indices;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* ori_vertex_pos;
	double* mass_inv_;
	int* tet_vertex;
	double inversion_time = 1.0;
	double collision_time = 1.0;
	int* record_position_num;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		obj_No = i + cloth->size();
		size = tetrahedron->data()[i].mesh_struct.indices.size();
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		vertex_pos = vertex_position[i + cloth->size()];
		ori_vertex_pos = record_vertex_position[i + cloth->size()].data();
		mass_inv_ = tetrahedron->data()[i].mesh_struct.mass_inv.data();
		record_position_num = collision.indicate_if_involved_in_last_color[i + cloth->size()][0].data();

		for (int j = 0; j < size; ++j) {
			tet_vertex = indices[j].data();
			if (record_position_num[tet_vertex[0]] || record_position_num[tet_vertex[1]] || record_position_num[tet_vertex[2]] || record_position_num[tet_vertex[3]]) {
				//if (mass_inv_[indices[j][0]] != 0.0 || mass_inv_[indices[j][1]] != 0.0 || mass_inv_[indices[j][2]] != 0.0 || mass_inv_[indices[j][3]] != 0.0) {
				if (inversionTest::TetInversionTest(ori_vertex_pos[tet_vertex[0]].data(), ori_vertex_pos[tet_vertex[1]].data(),
					ori_vertex_pos[tet_vertex[2]].data(), ori_vertex_pos[tet_vertex[3]].data(), vertex_pos[tet_vertex[0]].data(),
					vertex_pos[tet_vertex[1]].data(), vertex_pos[tet_vertex[2]].data(), vertex_pos[tet_vertex[3]].data(), &inversion_time)) {
					if (inversion_time < collision_time) {
						collision_time = inversion_time;
					}
				}
			}
		}
	}
	if (collision_time < 1.0) {
		collision_time *= 0.9;
	}


	if (collision_time <  0.0) {
		std::cout << "error for inversion test time warm start()" << std::endl;
	}
	if (collision_time == 0.0) {
		std::cout << "inversion test time warm start() equals 0" << std::endl;
	}

	if (collision_time < collision.collision_time) {
		collision.collision_time = collision_time;
	}
}



//LAST_COLOR_INVERSION_TEST
void XPBD_IPC::computeLastColorInversion(int thread_No)
{
	int obj_No;

	std::array<int, 4>* indices;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* ori_vertex_pos;
	double* mass_inv_;
	int* tet_vertex;
	double inversion_time = 1.0;
	double collision_time = 1.0;

	unsigned int start, end;
	collision.collision_time_thread[thread_No] = 1.0;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		obj_No = i + cloth->size();
		start = tet_index_begin_per_thread[i + cloth->size()][thread_No];
		end = tet_index_begin_per_thread[i + cloth->size()][thread_No + 1];
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		vertex_pos = vertex_position[i + cloth->size()];
		ori_vertex_pos = record_vertex_position[i + cloth->size()].data();
		mass_inv_ = tetrahedron->data()[i].mesh_struct.mass_inv.data();

		for (int j = start; j < end; ++j) {
			if (is_tet_arap_grad_compute[prefix_sum_of_every_tet_index[i] + j].test(std::memory_order_relaxed)) {
				tet_vertex = indices[j].data();
				if (inversionTest::TetInversionTest(ori_vertex_pos[tet_vertex[0]].data(), ori_vertex_pos[tet_vertex[1]].data(),
					ori_vertex_pos[tet_vertex[2]].data(), ori_vertex_pos[tet_vertex[3]].data(), vertex_pos[tet_vertex[0]].data(),
					vertex_pos[tet_vertex[1]].data(), vertex_pos[tet_vertex[2]].data(), vertex_pos[tet_vertex[3]].data(), &inversion_time)) {
					if (inversion_time < collision_time) {
						collision_time = inversion_time;
					}
				}
			}
		}
	}

	if (collision_time < 1.0) {
		collision_time *= 0.9;
	}
	if (collision.collision_time_thread[thread_No] > collision_time) {
		collision.collision_time_thread[thread_No] = collision_time;
	}


	if (collision_time < 0.0) {
		std::cout << "error for inversion test time last color ()" << std::endl;
	}
	if (collision_time == 0.0) {
		std::cout << "inversion test time last color () equas 0" << std::endl;
	}
}


//INVERSION_TEST
void XPBD_IPC::inversionTest(int thread_No)
{
	unsigned int start, end;
	unsigned int* tet_around_a_group;

	MatrixXd grad;
	grad.resize(3, 4);
	std::array<int, 4>* tet_indices_;

	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* ori_vertex_pos;
	double stiffness = 0.0;
	Matrix<double, 3, 4>* A;
	double* volume;

	int prefix_sum_start;

	int k = 0;

	//int color_group_index;
	unsigned int* tet_in_a_group;

	//solve all colors except the last color
	double* grad_address;

	double collision_time = 1.0;
	double inversion_time;
	int* tet_vertex;
	
	collision.collision_time_thread[thread_No] = 1.0;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		vertex_pos = vertex_position[i + cloth->size()];
		ori_vertex_pos = record_collision_free_vertex_position_address[i + cloth->size()];
		tet_indices_ = tet_indices[i + cloth->size()];
		start = tet_index_begin_per_thread[i+cloth->size()][thread_No];
		end = tet_index_begin_per_thread[i+cloth->size()][thread_No+1];
		for (auto j = start; j < end; ++j) {
			tet_vertex = tet_indices_[j].data();
			if(inversionTest::TetInversionTest(ori_vertex_pos[tet_vertex[0]].data(), ori_vertex_pos[tet_vertex[1]].data(),
				ori_vertex_pos[tet_vertex[2]].data(), ori_vertex_pos[tet_vertex[3]].data(), vertex_pos[tet_vertex[0]].data(),
				vertex_pos[tet_vertex[1]].data(), vertex_pos[tet_vertex[2]].data(), vertex_pos[tet_vertex[3]].data(), &inversion_time)) {
				if (inversion_time < collision_time) {
					collision_time = inversion_time;
				}
			}
		}
	}

	if (collision_time < 1.0) {
		collision_time *= 0.9;
	}
	if (collision_time < 0.0) {
		std::cout << "error for inversion test time inversionTest()" << std::endl;
	}
	if (collision_time == 0.0) {
		std::cout << "inversion test time inversionTest() equals 0" << std::endl;
	}

	if(collision.collision_time_thread[thread_No] >collision_time){
		collision.collision_time_thread[thread_No] = collision_time;
	}
}

//TET_GRAD
void XPBD_IPC::tetGrad(int thread_No)
{
	MatrixXd grad;
	grad.resize(3, 4);
	std::array<int, 4>* tet_indices_;

	std::array<double, 3>* vertex_pos;
	double stiffness = 0.0;
	Matrix<double, 3, 4>* A;
	double* volume;

	int prefix_sum_start;
	double* grad_address;

	double* com_grad = common_grad[thread_No].data();
	unsigned int tet_prefix_sum_start;

	int start, end;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		prefix_sum_start = vertex_index_prefix_sum_obj[i + cloth->size()];

		vertex_pos = vertex_position[i + cloth->size()];
		A = tetrahedron->data()[i].mesh_struct.A.data();
		volume = tet_volume[i + cloth->size()];
		tet_indices_ = tet_indices[i + cloth->size()];
		tet_prefix_sum_start = prefix_sum_of_every_tet_index[i];
		start = tet_index_begin_per_thread[i + cloth->size()][thread_No];
		end = tet_index_begin_per_thread[i + cloth->size()][thread_No + 1];
		for (int j = start; j < end; ++j)
		{
			second_order_constraint.setARAPGrad(grad, vertex_pos, stiffness, A[j],// ,
				volume[j], tet_indices_[j].data()); //,
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[j][0]);
			*grad_address += grad.data()[0];
			*(grad_address + 1) += grad.data()[1];
			*(grad_address + 2) += grad.data()[2];
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[j][1]);
			*grad_address += grad.data()[3];
			*(grad_address + 1) += grad.data()[4];
			*(grad_address + 2) += grad.data()[5];
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[j][2]);
			*grad_address += grad.data()[6];
			*(grad_address + 1) += grad.data()[7];
			*(grad_address + 2) += grad.data()[8];
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[j][3]);
			*grad_address += grad.data()[9];
			*(grad_address + 1) += grad.data()[10];
			*(grad_address + 2) += grad.data()[11];
		}
	}
}


//UPDATE_TET_GRAD_SHARED
void XPBD_IPC::tetGradForColor(int thread_No, unsigned int color_No)
{
	MatrixXd grad;
	grad.resize(3, 4);
	std::array<int, 4>* tet_indices_;

	std::array<double, 3>* vertex_pos;
	double stiffness = 0.0;
	Matrix<double, 3, 4>* A;
	double* volume;

	int prefix_sum_start;

	int color_group_index;
	//solve all colors except the last color
	double* grad_address;

	double* com_grad = common_grad[thread_No].data();
	unsigned int tet_prefix_sum_start;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		color_group_index = inner_iteration_number % tet_color_groups[i]->size();
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		prefix_sum_start = vertex_index_prefix_sum_obj[i + cloth->size()];

		vertex_pos = vertex_position[i + cloth->size()];
		A = tetrahedron->data()[i].mesh_struct.A.data();

		volume = tet_volume[i + cloth->size()];
		tet_indices_ = tet_indices[i + cloth->size()];

		tet_prefix_sum_start = prefix_sum_of_every_tet_index[i];

		auto start = tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index][color_No].begin()+tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color_No][thread_No];
		auto end = tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index][color_No].begin() + tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color_No][thread_No+1];
		//start = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color_No][thread_No];

		for (auto j = start; j < end; ++j) {
			second_order_constraint.setARAPGrad(grad, vertex_pos, stiffness, A[*j],// ,
				volume[*j], tet_indices_[*j].data()); //,
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[*j][0]);
			*grad_address += grad.data()[0];
			*(grad_address + 1) += grad.data()[1];
			*(grad_address + 2) += grad.data()[2];
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[*j][1]);
			*grad_address += grad.data()[3];
			*(grad_address + 1) += grad.data()[4];
			*(grad_address + 2) += grad.data()[5];
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[*j][2]);
			*grad_address += grad.data()[6];
			*(grad_address + 1) += grad.data()[7];
			*(grad_address + 2) += grad.data()[8];
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[*j][3]);
			*grad_address += grad.data()[9];
			*(grad_address + 1) += grad.data()[10];
			*(grad_address + 2) += grad.data()[11];
		}


		if (color_No == tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index].size() - 1) {
			for (auto j = start; j < end; ++j) {
				is_tet_arap_grad_compute[tet_prefix_sum_start + *j].test_and_set(std::memory_order_relaxed);
			}
		}

		start = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_groups[color_group_index][color_No].begin() + tetrahedron->data()[i].mesh_struct.tet_around_tet_color_groups_start_per_thread[color_group_index][color_No][thread_No];
		end = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_groups[color_group_index][color_No].begin()+ tetrahedron->data()[i].mesh_struct.tet_around_tet_color_groups_start_per_thread[color_group_index][color_No][thread_No+1];
		for (auto j = start;	j < end; ++j) {
			second_order_constraint.setARAPGrad(grad, vertex_pos, stiffness, A[*j],// ,
				volume[*j], tet_indices_[*j].data());
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[*j][0]);
			*grad_address += grad.data()[0];
			*(grad_address+1) += grad.data()[1];
			*(grad_address+2) += grad.data()[2];
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[*j][1]);
			*grad_address += grad.data()[3];
			*(grad_address + 1) += grad.data()[4];
			*(grad_address + 2) += grad.data()[5];
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[*j][2]);
			*grad_address += grad.data()[6];
			*(grad_address + 1) += grad.data()[7];
			*(grad_address + 2) += grad.data()[8];
			grad_address = com_grad + 3 * (prefix_sum_start + tet_indices_[*j][3]);
			*grad_address += grad.data()[9];
			*(grad_address + 1) += grad.data()[10];
			*(grad_address + 2) += grad.data()[11];
		
		}

		if (color_No == tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index].size() - 1) {
			for (auto j = start; j < end; ++j) {
				is_tet_arap_grad_compute[tet_prefix_sum_start + *j].test_and_set(std::memory_order_relaxed);
			}
		}
	}

}



//void XPBD_IPC::testPrintOut()
//{
//	for (int i = 0; i < tetrahedron->size(); ++i) {
//		for (int j = 0; j < tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups.size(); ++j) {
//			for (int k = 0; k < tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[j].size(); ++k) {
//				std::cout < "print tet_in_a_group_start_per_thread_groups " << std::endl;
//				for (int m = 0; m < total_thread_num; ++m) {
//					std::cout << tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[j][k][m] << " ";
//				}
//				std::cout << std::endl;
//			}
//		}
//
//	}
//}

void XPBD_IPC::setCollisionPairTetNeighborGrad(std::vector<unsigned int>& tet_involved, int start, int end, double* com_grad)
{
	double stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;

	auto start_ = tet_involved.begin() + start;
	auto end_ = tet_involved.begin() + end;
	MatrixXd grad;
	grad.resize(3, 4);
	double* grad_address;

	int* tet_indices_;

	for (auto i = start_; i < end_; i += 2) {
		if (!is_tet_arap_grad_compute[prefix_sum_of_every_tet_index[*i] + *(i+1)].test_and_set(std::memory_order_relaxed)) {
			tet_indices_ = tet_indices[*i][*(i + 1)].data();
			second_order_constraint.setARAPGrad(grad, vertex_position[*i], stiffness, tet_A[*i][*(i+1)],// ,
				tet_volume[*i][*(i + 1)], tet_indices_); //,

			grad_address = com_grad + 3 * (vertex_index_prefix_sum_obj[*i] + tet_indices_[0]);
			*grad_address += grad.data()[0];
			*(grad_address + 1) += grad.data()[1];
			*(grad_address + 2) += grad.data()[2];
			grad_address = com_grad + 3 * (vertex_index_prefix_sum_obj[*i] + tet_indices_[1]);
			*grad_address += grad.data()[3];
			*(grad_address + 1) += grad.data()[4];
			*(grad_address + 2) += grad.data()[5];

			grad_address = com_grad + 3 * (vertex_index_prefix_sum_obj[*i] + tet_indices_[2]);
			*grad_address += grad.data()[6];
			*(grad_address + 1) += grad.data()[7];
			*(grad_address + 2) += grad.data()[8];

			grad_address = com_grad + 3 * (vertex_index_prefix_sum_obj[*i] + tet_indices_[3]);
			*grad_address += grad.data()[9];
			*(grad_address + 1) += grad.data()[10];
			*(grad_address + 2) += grad.data()[11];
		}		
	}
}



//void XPBD_IPC::setCollisionPairTetGrad(int tet_obj_No, unsigned int start, unsigned int end, unsigned int* element, std::vector<unsigned int>* tet_around_an_element)
//{
//	double stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
//	unsigned int prefix_sum_start = prefix_sum_of_every_tet_index[tet_obj_No];
//
//	std::array<double, 3>* vertex_pos = vertex_position[tet_obj_No + cloth->size()];
//	Matrix<double, 3, 4>* A = tetrahedron->data()[tet_obj_No].mesh_struct.A.data();
//
//	double* volume = tet_volume[tet_obj_No + cloth->size()];
//	std::array<int, 4>*  tet_indices_ = tet_indices[tet_obj_No + cloth->size()];
//	MatrixXd grad;
//	grad.resize(3, 4);
//
//	for (auto j = start; j < end; ++j) {
//		for (auto k = tet_around_an_element[element[j]].begin(); k < tet_around_an_element[element[j]].end(); ++k) {
//			if (!is_tet_arap_grad_compute[prefix_sum_start + *k]) {
//				is_tet_arap_grad_compute[prefix_sum_start + *k] = true;
//				second_order_constraint.setARAPGrad(grad, vertex_pos, stiffness, A[*k],// ,
//					volume[*k], tet_indices_[*k].data()); //,
//				memcpy(store_tet_arap_grad.data() + 12 * (prefix_sum_start + *k), grad.data(), 96);//
//			}
//		}
//	}
//}


////UPDATE_TET_GRAD_SHARED_COLLISION
//void XPBD_IPC::tetGradForColorCollision(unsigned int color_No)
//{		// solve tet around all collision pair 
//	int color_group_index;
//	unsigned int start, end;
//	unsigned int* element;
//	std::vector<unsigned int>* tet_around_an_element;
//	
//	int obj_No;
//
//	for (int i = 0; i < tetrahedron->size(); ++i) {
//		color_group_index = inner_iteration_number % tet_color_groups[i]->size();
//		if (color_No != tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index].size() - 1) {
//			continue;
//		}
//
//		obj_No = cloth->size() + i;
//		setCollisionPairTetGrad(i, collision.tet_involved_in_collision_start_per_thread[obj_No][color_group_index][thread_No],
//			collision.tet_involved_in_collision_start_per_thread[obj_No][color_group_index][thread_No + 1],
//			collision.tet_involved_in_collision[obj_No][color_group_index].data(), prefix_sum_of_every_tet_index[i]);
//		//if (has_collider) {
//		//	//tet around vt collider
//		//	end = collision.vertex_index_collide_with_collider_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No + 1];
//		//	start = collision.vertex_index_collide_with_collider_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No];
//		//	element = collision.vertex_index_collide_with_collider[i + cloth->size()].data();
//		//	tet_around_an_element = tet_around_vertex[i + cloth->size()];
//		//	setCollisionPairTetGrad(i, start, end, element, tet_around_an_element);
//		//	//tet around ee collider
//		//	end = collision.edge_index_collide_with_collider_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No + 1];
//		//	start = collision.edge_index_collide_with_collider_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No];
//		//	element = collision.edge_index_collide_with_collider[i + cloth->size()].data();
//		//	tet_around_an_element = tet_around_edge[i + cloth->size()];
//		//	setCollisionPairTetGrad(i, start, end, element, tet_around_an_element);
//		//	//tet around tv collider
//		//	end = collision.triangle_index_collide_with_collider_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No + 1];
//		//	start = collision.triangle_index_collide_with_collider_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No];
//		//	element = collision.triangle_index_collide_with_collider[i + cloth->size()].data();
//		//	tet_around_an_element = tet_around_triangle[i + cloth->size()];
//		//	setCollisionPairTetGrad(i, start, end, element, tet_around_an_element);
//		//}
//		////tet around vt
//		//end = collision.vertex_index_has_VT_pair_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No + 1];
//		//start = collision.vertex_index_has_VT_pair_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No];
//		//element = collision.vertex_index_has_VT_pair[i + cloth->size()].data();
//		//tet_around_an_element = tet_around_vertex[i + cloth->size()];
//		//setCollisionPairTetGrad(i, start, end, element, tet_around_an_element);
//		////tet around ee
//		//end = collision.edge_index_has_EE_pair_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No + 1];
//		//start = collision.edge_index_has_EE_pair_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No];
//		//element = collision.edge_index_has_EE_pair[i + cloth->size()].data();
//		//tet_around_an_element = tet_around_edge[i + cloth->size()];
//		//setCollisionPairTetGrad(i, start, end, element, tet_around_an_element);
//		////tet around tv
//		//end = collision.triangle_index_has_TV_pair_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No + 1];
//		//start = collision.triangle_index_has_TV_pair_start_per_thread[(total_thread_num + 1) * (i + cloth->size()) + thread_No];
//		//element = collision.triangle_index_has_TV_pair[i + cloth->size()].data();
//		//tet_around_an_element = tet_around_triangle[i + cloth->size()];
//		//setCollisionPairTetGrad(i, start, end, element, tet_around_an_element);
//
//	}
//}

//UPDATE_TET_GRAD_SHARED_COLLISION_NEIGHBOR
void XPBD_IPC::tetGradForColorCollisionNeighbor(int thread_No, unsigned int color_No)
{

	int thread_id;
	if (collision.tet_involve_in_collision_start_per_thread[(thread_No + 1) << 1] > collision.tet_involve_in_collision_start_per_thread[thread_No << 1]) {
		thread_id = collision.tet_involve_in_collision_start_per_thread[thread_No << 1];
		setCollisionPairTetNeighborGrad(collision.tet_involve_in_collision[thread_id],
			collision.tet_involve_in_collision_start_per_thread[(thread_No << 1) + 1], collision.tet_involve_in_collision[thread_id].size(), common_grad[thread_No].data());


		for (int i = collision.tet_involve_in_collision_start_per_thread[thread_No << 1] + 1;
			i < collision.tet_involve_in_collision_start_per_thread[(thread_No + 1) << 1]; ++i) {
			setCollisionPairTetNeighborGrad(collision.tet_involve_in_collision[i],
				0, collision.tet_involve_in_collision[i].size(), common_grad[thread_No].data());
		}

		thread_id = collision.tet_involve_in_collision_start_per_thread[(thread_No + 1) << 1];
		setCollisionPairTetNeighborGrad(collision.tet_involve_in_collision[thread_id],
			0, collision.tet_involve_in_collision_start_per_thread[(thread_No << 1) + 3], common_grad[thread_No].data());
	}
	else {
		thread_id = collision.tet_involve_in_collision_start_per_thread[thread_No << 1];
		setCollisionPairTetNeighborGrad(collision.tet_involve_in_collision[thread_id],
			collision.tet_involve_in_collision_start_per_thread[(thread_No << 1) + 1], collision.tet_involve_in_collision_start_per_thread[(thread_No << 1) + 3], common_grad[thread_No].data());
	}
}



//void XPBD_IPC::computeTetHessianInAdvance(int thread_No, int color_No)
//{
//	unsigned int tet_around_tet_color_group_end;
//	unsigned int* tet_around_a_group;
//	MatrixXd Hessian; 
//	MatrixXd grad;
//	Hessian.resize(4, 4);
//	grad.resize(3, 4);
//
//
//	std::array<double, 3>* vertex_pos;
//	double stiffness=0.0;
//	if (!tetrahedron->empty()) {
//		stiffness = tetrahedron->data()[0].collision_stiffness[0];
//	}
//	else {
//		stiffness = cloth->data()[0].collision_stiffness[0];
//	}
//	Matrix<double, 3, 4>* A;
//	double* volume;
//
//	std::array<int, 4>*tet_indices_;
//	unsigned int prefix_sum_start;
//
//
//
//	for (int i = 0; i < tetrahedron->size(); ++i) {
//		if (color_No >= unconnected_tet_index[i]->size()) {
//			continue;
//		}
//		tet_around_tet_color_group_end = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group_start_per_thread[color_No][thread_No + 1];
//		tet_around_a_group = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group[color_No].data();
//
//		vertex_pos = vertex_position[i + cloth->size()];
//		A = tetrahedron->data()[i].mesh_struct.A.data();
//
//		volume = tet_volume[i + cloth->size()];
//		tet_indices_ = tet_indices[i + cloth->size()];
//
//		prefix_sum_start = prefix_sum_of_every_tet_index[i];
//
//		//std::cout << "color " << color_No << " " << tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group_start_per_thread[color_No][thread_No] 
//		//	<< " " << tet_around_tet_color_group_end<<" "<< prefix_sum_start << std::endl;
//
//		for (int j = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group_start_per_thread[color_No][thread_No];
//			j < tet_around_tet_color_group_end; ++j) {
//			if (tet_around_tet_color_group_end > tetrahedron->data()[i].mesh_struct.tet_around_tet_color_group[color_No].size()) {
//				std::cout << "error to large size " << std::endl;
//			}
//
//			if (tet_around_a_group[j]>= tetrahedron->data()[i].mesh_struct.indices.size()) {
//				std::cout << "erro  tet index " << std::endl;
//			}
//
//			second_order_constraint.setARAPGrad(grad, vertex_pos, stiffness, A[tet_around_a_group[j]],
//				volume[tet_around_a_group[j]], tet_indices_[tet_around_a_group[j]].data());
//			//memcpy(store_tet_arap_hessian.data() + 16 * (prefix_sum_start + tet_around_a_group[j]), Hessian.data(),128);
//			if (12 * (prefix_sum_start + tet_around_a_group[j]) >= store_tet_arap_grad.size()) {
//				std::cout<<"error "<<std::endl;
//			}
//			//memcpy(store_tet_arap_grad.data() + 12 * (prefix_sum_start + tet_around_a_group[j]), grad.data(), 96);
//		}
//
//
//		//for (int j = 0; j < tetrahedron->data()[i].mesh_struct.indices.size(); ++j) {
//		//	second_order_constraint.setARAPHessianGrad(Hessian, grad, vertex_pos, stiffness, A[j],
//		//		volume[j], tet_indices_[j].data());
//		//	memcpy(store_tet_arap_hessian.data() + 16 * (prefix_sum_start + j), Hessian.data(), 128);
//		//	memcpy(store_tet_arap_grad.data() + 12 * (prefix_sum_start + j), grad.data(), 96);
//		//}
//	}
//}



bool XPBD_IPC::convergeCondition(unsigned int iteration_num)
{
	if (iteration_num <min_outer_iteration) {//max_iteration_number
		std::cout << iteration_num << " " << max_displacement << std::endl;
		return false;
	}

	if (iteration_num > outer_max_iteration_number - 1) {
		return true;
	}

	if (max_displacement > max_move_standard) {
		std::cout << iteration_num << " " << max_displacement << std::endl;
		if (iteration_num > 100) {
			if (dis_record.size() > 3) {
				dis_record.clear();
			}
			for (int i = 0; i < dis_record.size(); ++i) {
				if (abs(dis_record[i] - max_displacement) < 1e-8) {
					return true;
				}
			}
			dis_record.emplace_back(max_displacement);			
		}		
		return false;
	}
	else {
		std::cout << iteration_num << " " << max_displacement << std::endl;
	}
	//if (grad_max > max_move_standard)
	//{
	//	std::cout << iteration_num << " " << grad_max << std::endl;
	//	return false;
	//}
	return true;

	//return checkMaxDisplacement();

}

void XPBD_IPC::updatePosition()
{

	unsigned int vertex_num;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_num = mesh_struct[i]->vertex_position.size();
		memcpy(vertex_position_render[i][0].data(), vertex_position[i][0].data(), 24 * vertex_num);
	}

	for (unsigned int i = 0; i < collider->size(); ++i) {
		memcpy(collider->data()[i].mesh_struct.vertex_for_render[0].data(), collider->data()[i].mesh_struct.vertex_position[0].data(), 24 * collider->data()[i].mesh_struct.vertex_position.size());
	}
}


double XPBD_IPC::computeWarmStartEnergy()
{
	energy = 1e-15;
	double inertial_energy =computeInertialEnergy();
	energy += inertial_energy;
	if (perform_collision) {
		double barrier_energy = 0.0;
		barrier_energy = computeBarrierEnergy();
		energy += barrier_energy;
	}
	return energy;
}


void XPBD_IPC::computeCurrentEnergy()
{
	energy = 1e-15;
	double inertial_energy = computeInertialEnergy();
	energy += inertial_energy;
	double ARAP_energy = computeCurrentARAPEnergy();
	energy += ARAP_energy;
	if (perform_collision) {
		double barrier_energy = 0.0;
		barrier_energy = computeBarrierEnergy();
		energy += barrier_energy;		
	}
}

void XPBD_IPC::firstNewtonCD()
{
	previous_energy = energy;
	newtonCDTet();
}


//void XPBD_IPC::firstOnlyInertialCollision()
//{
//	unsigned int size;
//	MeshStruct* mesh_struct_;
//	double* volume;
//	std::array<double, 3>* vertex_pos;
//	std::array<double, 3>* last_step_vertex_pos;
//	std::array<double, 3>* initial_vertex_pos;
//	double* mass_inv;
//	std::array<double, 3>* sn_;
//	double* mass;
//	double colliision_stiffness;
//	int* vertex_index;
//	std::vector<bool>* is_surface_vertex;
//	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
//		mesh_struct_ = mesh_struct[i + cloth->size()];
//		size = tetrahedron->data()[i].mesh_struct.vertex_position.size();
//		volume = tetrahedron->data()[i].mesh_struct.volume.data();
//		vertex_pos = vertex_position[i + cloth->size()];
//		last_step_vertex_pos = initial_vertex_position[i + cloth->size()];
//		initial_vertex_pos = record_vertex_position[i + cloth->size()].data();
//		mass_inv = mesh_struct_->mass_inv.data();
//		mass = mesh_struct_->mass.data();
//		sn_ = sn[i + cloth->size()].data();
//		vertex_index = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
//		is_surface_vertex = &tetrahedron->data()[i].mesh_struct.vertex_on_surface;
//		colliision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
//		//std::cout << "=== " << std::endl;
//		for (unsigned int j = 0; j < size; ++j) {
//			if (mass_inv[j] != 0.0) {
//				solveInertialCollision(vertex_pos, initial_vertex_pos[j].data(),
//					last_step_vertex_pos[j].data(), sub_time_step, mass, volume, j, sn_,
//					colliision_stiffness, i + cloth->size(), (*is_surface_vertex)[j], vertex_index[j]);
//			}
//		}
//	}
//}


void XPBD_IPC::newtonCDTetWithCollision()
{
	unsigned int size;
	std::array<int, 4>* indices;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* last_step_vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	double* mass_inv;
	double stiffness;
	Matrix<double, 3, 4>* A;
	std::array<double, 3>* sn_;
	double* mass;
	double colliision_stiffness;
	int* vertex_index;
	std::vector<bool>* is_surface_vertex;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		mesh_struct_ = mesh_struct[i + cloth->size()];
		size = tetrahedron->data()[i].mesh_struct.vertex_position.size();
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i + cloth->size()];
		last_step_vertex_pos = record_collision_free_vertex_position[i + cloth->size()].data();
		initial_vertex_pos = record_vertex_position[i + cloth->size()].data();
		stiffness = tetrahedron->data()[i].ARAP_stiffness;
		A = tetrahedron->data()[i].mesh_struct.A.data();
		mass_inv = mesh_struct_->mass_inv.data();
		mass = mesh_struct_->mass.data();
		sn_ = sn[i + cloth->size()].data();
		vertex_index = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
		is_surface_vertex = &tetrahedron->data()[i].mesh_struct.vertex_on_surface;
		colliision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
		//std::cout << "=== " << std::endl;
		for (unsigned int j = 0; j < size; ++j) {
			if (mass_inv[j] != 0.0) {
					solveNewtonCDTetWithCollision(vertex_pos, initial_vertex_pos[j].data(),
					last_step_vertex_pos[j].data(),	stiffness, sub_time_step, A,
						mesh_struct_->vertex_tet_index[j], indices, mass, volume, j, sn_,
						colliision_stiffness, i + cloth->size(), (*is_surface_vertex)[j],vertex_index[j]);
			}
		}
	}
}




//COMPUTE_RPREVIOUS_COLOR_BARRIER_ENERGY
void XPBD_IPC::computePreviousColorCollisionEnergy(int thread_No)
{
	double energy = 0.0;
	double stiffness = 0.0;
	if (!tetrahedron->empty()) {
		stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		stiffness = cloth->data()[0].collision_stiffness[0];
	}
	//vt
	energy += computeVTEnergy(&collision.record_vt_pair_sum_all_thread, triangle_indices.data(), vertex_position.data(), vertex_position.data(),
		stiffness, collision.record_vt_pair_d_hat.data(),
		collision.vertex_belong_to_color_group, 1, collision.vt_pair_sum_start_per_thread[thread_No],
		collision.vt_pair_sum_start_per_thread[thread_No + 1]);
	//ee
	energy += computeEEEnergy(&collision.record_ee_pair_sum_all_thread, edge_vertices.data(), edge_vertices.data(),
		vertex_position.data(), vertex_position.data(), stiffness, collision.record_ee_pair_d_hat.data(),
		collision.vertex_belong_to_color_group, 1, collision.ee_pair_sum_start_per_thread[thread_No],
		collision.ee_pair_sum_start_per_thread[thread_No+1] );

	if (has_collider) {
		//vt c
		energy += computeVTEnergy(&collision.record_vt_collider_pair_sum_all_thread, triangle_indices_collider.data(), vertex_position.data(), vertex_position_collider.data(),
			stiffness, collision.record_vt_collider_pair_d_hat.data(),
			collision.vertex_belong_to_color_group, 2, collision.vt_collider_pair_sum_start_per_thread[thread_No],
			collision.vt_collider_pair_sum_start_per_thread[thread_No + 1]);

		//tv c
		energy += computeVTEnergy(&collision.record_tv_collider_pair_sum_all_thread, triangle_indices.data(), vertex_position_collider.data(), vertex_position.data(),
			stiffness, collision.record_tv_collider_pair_d_hat.data(),
			collision.vertex_belong_to_color_group, 3, collision.tv_collider_pair_sum_start_per_thread[thread_No],
			collision.tv_collider_pair_sum_start_per_thread[thread_No + 1]);
		//ee c
		energy += computeEEEnergy(&collision.record_ee_collider_pair_sum_all_thread, edge_vertices.data(), collider_edge_vertices.data(),
			vertex_position.data(), vertex_position_collider.data(), stiffness, collision.record_ee_collider_pair_d_hat.data(),
			collision.vertex_belong_to_color_group, 2, collision.ee_collider_pair_sum_start_per_thread[thread_No],
			collision.ee_collider_pair_sum_start_per_thread[thread_No + 1]);
	}
	if (floor->exist) {
		 energy += computeFloorEnergy(1, stiffness, &collision.record_vertex_collide_with_floor_sum_all_thread, collision.floor_pair_sum_start_per_thread[thread_No],
			collision.floor_pair_sum_start_per_thread[thread_No+1]);
	}

	energy_per_thread[thread_No] = energy;
}


double XPBD_IPC::computePreviousColorCollisionEnergy()
{
	double energy = 0.0;
	double stiffness = 0.0;
	thread->assignTask(this, COMPUTE_RPREVIOUS_COLOR_BARRIER_ENERGY);
	for (int i = 0; i < total_thread_num; ++i)	{
		energy += energy_per_thread[i];
	}
	return energy;
}


double XPBD_IPC::computeBarrierEnergy()
{
	double energy = 0.0;
	thread->assignTask(this, COMPUTE_BARRIER_ENERGY);
	for (int i = 0; i < total_thread_num; ++i) {
		energy += energy_per_thread[i];
	}
	return energy;
}




double XPBD_IPC::computeEECollisionEnergy(unsigned int**edge_edge_pair_by_vertex_, unsigned int** edge_edge_pair_num_record_,
	std::array<double, 3>** edge_0_position, std::array<double, 3>** edge_1_position, unsigned int close_pair_num,
	unsigned int** edge_0_vertex, unsigned int** edge_1_vertex, bool is_self)
{
	double energy = 0.0;
	double	collision_stiffness;
	if (!tetrahedron->empty()) {
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	int size;
	unsigned int* edge_edge_pair_by_vertex;
	unsigned int* edge_edge_pair_num_record;


	unsigned int* pair_index;
	unsigned int pair_num;

	unsigned int* edge_vertex_index;
	unsigned int* edge_vertex_1_index;
	
	if (is_self) {
		for (int i = 0; i < total_obj_num; ++i) {
			edge_edge_pair_by_vertex = edge_edge_pair_by_vertex_[i];
			edge_edge_pair_num_record = edge_edge_pair_num_record_[i];
			size = mesh_struct[i]->triangle_indices.size();

			for (int j = 0; j < size; ++j) {
				pair_index = edge_edge_pair_by_vertex + close_pair_num * j;
				pair_num = edge_edge_pair_num_record[j];
				edge_vertex_index = edge_0_vertex[i] + (j << 1);
				for (int k = 0; k < pair_num; k += 3) {
					if (i > pair_index[k]) {
						break;
					}
					if (i == pair_index[k] && j > pair_index[k + 1]) {
						break;
					}
					edge_vertex_1_index = edge_1_vertex[pair_index[k]] + (pair_index[k + 1] << 1);
					energy += compute_energy.computeBarrierEnergy(edge_0_position[i][edge_vertex_index[0]].data(),
						edge_0_position[i][edge_vertex_index[1]].data(),
						edge_1_position[pair_index[k]][edge_vertex_1_index[0]].data(),
						edge_1_position[pair_index[k]][edge_vertex_1_index[1]].data(), collision_stiffness, collision.d_hat_2, false);
				}
			}
		}
	}
	else{
		for (int i = 0; i < total_obj_num; ++i) {
			edge_edge_pair_by_vertex = edge_edge_pair_by_vertex_[i];
			edge_edge_pair_num_record = edge_edge_pair_num_record_[i];
			size = mesh_struct[i]->triangle_indices.size();

			for (int j = 0; j < size; ++j) {
				pair_index = edge_edge_pair_by_vertex + close_pair_num * j;
				pair_num = edge_edge_pair_num_record[j];
				edge_vertex_index = edge_0_vertex[i] + (j << 1);
				for (int k = 0; k < pair_num; k += 2) {
					edge_vertex_1_index = edge_1_vertex[pair_index[k]] + (pair_index[k + 1] << 1);
					energy += compute_energy.computeBarrierEnergy(edge_0_position[i][edge_vertex_index[0]].data(),
						edge_0_position[i][edge_vertex_index[1]].data(),
						edge_1_position[pair_index[k]][edge_vertex_1_index[0]].data(),
						edge_1_position[pair_index[k]][edge_vertex_1_index[1]].data(), collision_stiffness, collision.d_hat_2, false);
				}
			}
		}
	}
		return energy;
}



double XPBD_IPC::computeVTCollisionEnergy(unsigned int** vertex_triangle_pair_by_vertex_, unsigned int** vertex_triangle_pair_num_record_,
	std::array<double,3>** triangle_position, std::array<double,3>** vertex_position, unsigned int close_pair_num, bool is_TV,
	std::array<int, 3>** triangle_vertex)
{
	double energy=0.0;
	double	collision_stiffness;
	if (!tetrahedron->empty()) {
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	int size;
	unsigned int* vertex_triangle_pair_by_vertex;
	unsigned int* vertex_triangle_pair_num_record;


	unsigned int* pair_index;
	unsigned int pair_num;

	int* triangle_vertex_index;
	if (!is_TV) {
		for (int i = 0; i < total_obj_num; ++i) {
			vertex_triangle_pair_by_vertex = vertex_triangle_pair_by_vertex_[i];
			vertex_triangle_pair_num_record = vertex_triangle_pair_num_record_[i];
			if (i < cloth->size()) {
				size = mesh_struct[i]->vertex_position.size();
				for (int j = 0; j < size; ++j) {
					pair_index = vertex_triangle_pair_by_vertex + close_pair_num * j;
					pair_num = vertex_triangle_pair_num_record[j];
					for (int k = 0; k < pair_num; k += 2) {
						triangle_vertex_index = triangle_vertex[pair_index[k]][pair_index[k + 1]].data();
						energy+= compute_energy.computeBarrierEnergy(vertex_position[i][j].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[0]].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[1]].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[2]].data(), collision_stiffness, collision.d_hat_2, true);
					}
				}
			}
			else {
				vertex_triangle_pair_by_vertex = vertex_triangle_pair_by_vertex_[i];
				vertex_triangle_pair_num_record = vertex_triangle_pair_num_record_[i];
				int j;
				size = tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_index_on_sureface.size();
				for (int m = 0; m < size; ++m) {
					j = vertex_surface_to_global[i][m];
					pair_index = vertex_triangle_pair_by_vertex + close_pair_num * m;
					pair_num = vertex_triangle_pair_num_record[m];
					for (int k = 0; k < pair_num; k += 2) {
						triangle_vertex_index = triangle_vertex[pair_index[k]][pair_index[k + 1]].data();
						energy += compute_energy.computeBarrierEnergy(vertex_position[i][j].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[0]].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[1]].data(),
							triangle_position[pair_index[k]][triangle_vertex_index[2]].data(), collision_stiffness, collision.d_hat_2, true);
					}
				}
			}
		}
	}
	else {
		for (int i = 0; i < total_obj_num; ++i) {
			vertex_triangle_pair_by_vertex = vertex_triangle_pair_by_vertex_[i];
			vertex_triangle_pair_num_record = vertex_triangle_pair_num_record_[i];
			size = mesh_struct[i]->triangle_indices.size();

			for (int j = 0; j < size; ++j) {
				pair_index = vertex_triangle_pair_by_vertex + close_pair_num * j;
				pair_num = vertex_triangle_pair_num_record[j];
				triangle_vertex_index = triangle_vertex[i][j].data();
				for (int k = 0; k < pair_num; k += 2) {
					energy += compute_energy.computeBarrierEnergy(vertex_position[pair_index[k]][pair_index[k+1]].data(),
						triangle_position[i][triangle_vertex_index[0]].data(),
						triangle_position[i][triangle_vertex_index[1]].data(),
						triangle_position[i][triangle_vertex_index[2]].data(), collision_stiffness, collision.d_hat_2, true);
				}
			}

		}
	}
	return energy;
}




void XPBD_IPC::newtonVTCollisionBlock()
{
	int size;
	int pair_num;
	unsigned int* pair_index;
	unsigned int* vertex_triangle_pair_by_vertex;
	std::atomic_uint* vertex_triangle_pair_num_record;
	double stiffness = 0;
	double collision_stiffness =0.0;

	MeshStruct::Vertex* vertex;

	if (!tetrahedron->empty()) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}	
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	std::vector<unsigned int>tet_around;

	std::vector<unsigned int>* tet_around_vertex_;
	std::vector<unsigned int>* tet_around_triangle_;

	for (int i = 0; i <total_obj_num ; ++i) {
		vertex_triangle_pair_by_vertex = collision.vertex_triangle_pair_by_vertex[i];
		vertex_triangle_pair_num_record = collision.vertex_triangle_pair_num_record[i];
		vertex = mesh_struct[i]->vertices.data();
		if (i < cloth->size()) {
			size = mesh_struct[i]->vertex_position.size();
			tet_around_vertex_ = &tet_around;
			for (int j = 0; j < size; ++j) {
				pair_index = vertex_triangle_pair_by_vertex + collision.close_vt_pair_num * j;
				pair_num = vertex_triangle_pair_num_record[j];
				for (int k = 0; k < pair_num; k += 2) {
					if(pair_index[k]<cloth->size()){
						tet_around_triangle_ = &tet_around;;
					}
					else {
						tet_around_triangle_ = &tet_around_triangle[pair_index[k]][pair_index[k + 1]];
					}
					solveVT_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &vertex[j].face, &triangle_around_triangle[pair_index[k]][pair_index[k + 1]],
						&vertex[j].edge, &edge_around_triangle[pair_index[k]][pair_index[k + 1]], tet_around_vertex_, tet_around_triangle_, collision.d_hat_2,
						false,false);
				}				
			}
		}
		else {
			int j;
			size = tetrahedron->data()[i-cloth->size()].mesh_struct.vertex_index_on_sureface.size();
			for (int m = 0; m < size; ++m) {
				j = vertex_surface_to_global[i][m];
				tet_around_vertex_ = &tet_around_vertex[i][j];
				pair_index = vertex_triangle_pair_by_vertex + collision.close_vt_pair_num * m;
				pair_num = vertex_triangle_pair_num_record[m];
				for (int k = 0; k < pair_num; k += 2) {
					if (pair_index[k] < cloth->size()) {
						tet_around_triangle_ = &tet_around;
					}
					else {
						tet_around_triangle_ = &tet_around_triangle[pair_index[k]][pair_index[k + 1]];
					}
					solveVT_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &vertex[j].face, &triangle_around_triangle[pair_index[k]][pair_index[k + 1]],
						&vertex[j].edge, &edge_around_triangle[pair_index[k]][pair_index[k + 1]], tet_around_vertex_, tet_around_triangle_, collision.d_hat_2,
						false,false);
				}
			}		
		}		
	}
}


void XPBD_IPC::newtonTVColliderCollisionBlock()
{
	int size;
	int pair_num;
	unsigned int* pair_index;
	unsigned int* vertex_triangle_pair_by_triangle;
	std::atomic_uint* vertex_triangle_pair_num_record;
	double stiffness = 0;
	double collision_stiffness = 0.0;

	MeshStruct::Vertex* vertex;

	if (!tetrahedron->empty()) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	std::vector<unsigned int>tet_around;

	std::vector<unsigned int>* tet_around_v_;
	std::vector<unsigned int>* tet_around_t_;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_triangle_pair_by_triangle = collision.triangle_vertex_collider_pair_by_triangle[i];
		vertex_triangle_pair_num_record = collision.triangle_vertex_collider_pair_num_record[i];
		size = mesh_struct[i]->triangle_indices.size();
		if (i < cloth->size()) {
			tet_around_t_ = &tet_around;
			tet_around_v_ = &tet_around;
			for (unsigned int j = 0; j < size; ++j) {
				pair_index = vertex_triangle_pair_by_triangle + collision.close_tv_collider_pair_num * j;
				pair_num = vertex_triangle_pair_num_record[j];
				for (unsigned int k = 0; k < pair_num; k += 2) {
					solveVT_collisionBlock(pair_index[k], pair_index[k + 1], i, j,
						stiffness, sub_time_step, collision_stiffness, &tet_around , &triangle_around_triangle[i][j],
						&tet_around, &edge_around_triangle[i][j], tet_around_v_, tet_around_t_, collision.d_hat_2,
						true, false);
				}
			}
		}
		else {
			tet_around_v_ = &tet_around;
			for (unsigned int j = 0; j < size; ++j) {
				pair_index = vertex_triangle_pair_by_triangle + collision.close_tv_collider_pair_num * j;
				pair_num = vertex_triangle_pair_num_record[j];
				tet_around_t_ = &tet_around_triangle[i][j];
				for (unsigned int k = 0; k < pair_num; k += 2) {
					solveVT_collisionBlock(pair_index[k], pair_index[k + 1], i, j,
						stiffness, sub_time_step, collision_stiffness, &tet_around, &triangle_around_triangle[i][j],
						&tet_around, &edge_around_triangle[i][j], tet_around_v_, tet_around_t_, collision.d_hat_2,
						true, false);
				}
			}
			
		}
	}
}


void XPBD_IPC::newtonVTColliderCollisionBlock()
{
	int size;
	int pair_num;
	unsigned int* pair_index;
	unsigned int* vertex_triangle_pair_by_vertex;
	std::atomic_uint* vertex_triangle_pair_num_record;
	double stiffness = 0;
	double collision_stiffness = 0.0;

	MeshStruct::Vertex* vertex;

	if (!tetrahedron->empty()) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	std::vector<unsigned int>tet_around;

	std::vector<unsigned int>* tet_around_vertex_;
	std::vector<unsigned int>* tet_around_triangle_;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_triangle_pair_by_vertex = collision.vertex_obj_triangle_collider_pair_by_vertex[i];
		vertex_triangle_pair_num_record = collision.vertex_obj_triangle_collider_num_record[i];
		vertex = mesh_struct[i]->vertices.data();
		if (i < cloth->size()) {
			size = mesh_struct[i]->vertex_position.size();
			tet_around_vertex_ = &tet_around;
			tet_around_triangle_ = &tet_around;
			for (unsigned int j = 0; j < size; ++j) {
				pair_index = vertex_triangle_pair_by_vertex + collision.close_vt_collider_pair_num * j;
				pair_num = vertex_triangle_pair_num_record[j];
				for (unsigned int k = 0; k < pair_num; k += 2) {
					solveVT_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &vertex[j].face, &triangle_around_triangle[pair_index[k]][pair_index[k + 1]],
						&vertex[j].edge, &edge_around_triangle[pair_index[k]][pair_index[k + 1]], tet_around_vertex_, tet_around_triangle_, collision.d_hat_2,
						false,true);
				}
			}
		}
		else {
			int j;
			tet_around_triangle_ = &tet_around;
			size = tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_index_on_sureface.size();
			for (unsigned int m = 0; m < size; ++m) {
				j = vertex_surface_to_global[i][m];
				tet_around_vertex_ = &tet_around_vertex[i][j];
				pair_index = vertex_triangle_pair_by_vertex + collision.close_vt_collider_pair_num * m;
				pair_num = vertex_triangle_pair_num_record[m];
				for (unsigned int k = 0; k < pair_num; k += 2) {

					solveVT_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &vertex[j].face, &tet_around,
						&vertex[j].edge, &tet_around, tet_around_vertex_, tet_around_triangle_,
						collision.d_hat_2,false ,true);
				}
			}
		}
	}
}



void XPBD_IPC::newtonEEColliderCollisionBlock()
{
	int size;
	int pair_num;
	unsigned int* pair_index;
	unsigned int* edge_edge_pair_by_edge;
	std::atomic_uint* edge_edge_pair_num_record;
	double stiffness = 0;
	double collision_stiffness = 0.0;



	if (!tetrahedron->empty()) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	std::vector<unsigned int>* tet_around_edge_0;

	std::vector<unsigned int>* tet_around_edge_1;

	std::vector<unsigned int>* triangle_around_edge_;
	std::vector<unsigned int>* edge_around_edge_;

	std::vector<unsigned int>tet_around;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		edge_edge_pair_by_edge = collision.edge_edge_collider_pair_by_edge[i];
		edge_edge_pair_num_record = collision.edge_edge_collider_pair_num_record[i];
		triangle_around_edge_ = triangle_around_edge[i];
		edge_around_edge_ = edge_around_edge[i];
		size = mesh_struct[i]->edge_length.size();
		if (i < cloth->size()) {
			tet_around_edge_0 = &tet_around;
			tet_around_edge_1 = &tet_around;
			for (unsigned int j = 0; j < size; ++j) {
				pair_index = edge_edge_pair_by_edge + collision.close_ee_collider_pair_num * j;
				pair_num = edge_edge_pair_num_record[j];
				for (int k = 0; k < pair_num; k += 2) {
					solveEE_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &triangle_around_edge_[j], &tet_around,
						&edge_around_edge_[j], &tet_around,
						tet_around_edge_0, tet_around_edge_1, collision.d_hat_2, false, true);
				}
			}
		}
		else {
			tet_around_edge_1 = &tet_around;
			for (unsigned int j = 0; j < size; ++j) {

				tet_around_edge_0 = &tet_around_edge[i][j];
				pair_index = edge_edge_pair_by_edge + collision.close_ee_collider_pair_num * j;
				pair_num = edge_edge_pair_num_record[j];
				for (int k = 0; k < pair_num; k += 2) {
					solveEE_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &triangle_around_edge_[j], &tet_around,
						&edge_around_edge_[j], &tet_around,
						tet_around_edge_0, tet_around_edge_1, collision.d_hat_2, false, true);
				}
			}
		}
	}

}


void XPBD_IPC::newtonEECollisionBlock()
{
	int size;
	int pair_num;
	unsigned int* pair_index;
	unsigned int* edge_edge_pair_by_edge;
	std::atomic_uint* edge_edge_pair_num_record;
	double stiffness = 0;
	double collision_stiffness = 0.0;



	if (!tetrahedron->empty()) {
		stiffness = 0.5 * tetrahedron->data()[0].ARAP_stiffness;
		collision_stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness[0];
	}
	std::vector<unsigned int>* tet_around_edge_0;

	std::vector<unsigned int>* tet_around_edge_1;

	std::vector<unsigned int>* triangle_around_edge_;
	std::vector<unsigned int>* edge_around_edge_;

	std::vector<unsigned int>tet_around;

	for (int i = 0; i < total_obj_num; ++i) {
		edge_edge_pair_by_edge = collision.edge_edge_pair_by_edge[i];
		edge_edge_pair_num_record = collision.edge_edge_pair_number_record[i];
		triangle_around_edge_ = triangle_around_edge[i];
		edge_around_edge_ = edge_around_edge[i];
		size = mesh_struct[i]->edge_length.size();
		if (i < cloth->size()) {
			tet_around_edge_0 = &tet_around;
			for (int j = 0; j < size; ++j) {
				pair_index = edge_edge_pair_by_edge + collision.close_ee_pair_num * j;
				pair_num = edge_edge_pair_num_record[j];
				for (int k = 0; k < pair_num; k += 3) {
					if (i > pair_index[k]) {
						break;
					}
					else if (i == pair_index[k]) {
						if (j > pair_index[k + 1]) {
							break;
						}
					}
					if (pair_index[k] < cloth->size()) {
						tet_around_edge_1 = &tet_around;
					}
					else {
						tet_around_edge_1 = &tet_around_edge[pair_index[k]][pair_index[k + 1]];
					}
					solveEE_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &triangle_around_edge_[j], &triangle_around_edge[pair_index[k]][pair_index[k + 1]],
						&edge_around_edge_[j], &edge_around_edge[pair_index[k]][pair_index[k + 1]],
						tet_around_edge_0, tet_around_edge_1, collision.d_hat_2,false,false);
				}
			}
		}
		else {
			for (int j = 0; j < size; ++j) {
				tet_around_edge_0 = &tet_around_edge[i][j];
				pair_index = edge_edge_pair_by_edge + collision.close_ee_pair_num * j;
				pair_num = edge_edge_pair_num_record[j];
				for (int k = 0; k < pair_num; k += 3) {
					if (i > pair_index[k]) {
						break;
					}
					else if (i == pair_index[k]) {
						if (j > pair_index[k + 1]) {
							break;
						}
					}
					if (pair_index[k] < cloth->size()) {
						tet_around_edge_1 = &tet_around;
					}
					else {
						tet_around_edge_1 = &tet_around_edge[pair_index[k]][pair_index[k + 1]];
					}
					solveEE_collisionBlock(i, j, pair_index[k], pair_index[k + 1],
						stiffness, sub_time_step, collision_stiffness, &triangle_around_edge_[j], &triangle_around_edge[pair_index[k]][pair_index[k + 1]],
						&edge_around_edge_[j], &edge_around_edge[pair_index[k]][pair_index[k + 1]],
						tet_around_edge_0, tet_around_edge_1, collision.d_hat_2, false, false);
				}
			}
		}
	}
}


void XPBD_IPC::test()
{
	//std::atomic_uint* k;
	//k = new std::atomic_uint[3];
	//std::cout << "size " << sizeof(unsigned int) << " " << sizeof(std::atomic_uint) << std::endl;

	//for (int i = 0; i < 3; ++i) {		
	//	k[i].store(1);
	//}
	//for (int i = 0; i < 3; ++i) {
	//	std::cout << k[i] << std::endl;
	//}
	//memset(k+1, 0, 8);
	//for (int i = 0; i < 3; ++i) {
	//	std::cout <<k[i] << std::endl;
	//}
}


//TEST_MULTI
//void XPBD_IPC::testMulti(int thread_No)
//{
//	std::unique_lock<std::mutex> l(thread->m, std::defer_lock);
//
//	for (int i = 0; i < 100; i++) {
//		l.lock();
//		std::this_thread::sleep_for(std::chrono::milliseconds(1));
//		counter++;
//		l.unlock();
//	}
//
//	
//}



void XPBD_IPC::newtonCDBlock()
{
	//newtonEECollisionBlock();
	//newtonVTCollisionBlock();
	//if (has_collider) {
	//	newtonEEColliderCollisionBlock();
	//	newtonVTColliderCollisionBlock();
	//	newtonTVColliderCollisionBlock();
	//}
	newtonCDTetBlock();
}



void XPBD_IPC::computeGradient()
{
	for (int i = 0; i < total_thread_num; ++i) {
		memset(common_grad[i].data(), 0, 8 * common_grad[i].size());
	}
	thread->assignTask(this, TET_GRAD);
	if (perform_collision) {
		thread->assignTask(&collision, COLLISION_GRAD);
	}
	thread->assignTask(this, SUM_ALL_GRAD);

	thread->assignTask(this, GRAD_NORM);

	grad_max = grad_max_store[0];
	for (int i = 1; i < total_thread_num; ++i) {
		if (grad_max < grad_max_store[i]) {
			grad_max = grad_max_store[i];
		}
	}
}



void XPBD_IPC::solveNewtonCD_tetBlock()
{
	double ori_energy = 0.0;
	max_displacement = 0.0;
	memset(max_dis_record.data(), 0, 8 * max_dis_record.size());
	for (unsigned int i = 0; i < max_tet_color_num-1; ++i) {		//
		initialRecordHessian();
		//update shared tet, collision hessian 
		//tetGradForColor(i);
		if (perform_collision) {
			collision.computeHessian(i);
			//collision_compare.computeHessian(i);
			//compareIfRecordHessianIsRight(i);
		}
		thread->assignTask(this, UPDATE_TET_GRAD_SHARED, i);
		thread->assignTask(this, SUM_ALL_GRAD);		
		ori_energy = computePreviousColorEnergy(i);
		thread->assignTask(this, SOLVE_TET_BLOCK, i);	

		thread->assignTask(this, MAX_DISPLACEMENT_COLOR, i);	
		lineSearchFirstColor(i, ori_energy);



	}

	//testPrintOut();
	initialRecordHessian();
	memset(is_tet_arap_grad_compute, 0, prefix_sum_of_every_tet_index[tetrahedron->size()]<<2);
	
	//tetGradForColor(max_tet_color_num-1);
	if (perform_collision) {
		//thread->assignTask(this, UPDATE_TET_GRAD_SHARED_COLLISION, max_tet_color_num - 1);
		//thread->assignTask(this, UPDATE_TET_GRAD_SHARED_COLLISION_NEIGHBOR, max_tet_color_num - 1);
		//tetGradForColorCollisionNeighbor(max_tet_color_num - 1);
		collision.computeHessian(max_tet_color_num - 1);
		//collision_compare.computeHessian(max_tet_color_num - 1);
		compareIfRecordHessianIsRight(-1);
	}
	thread->assignTask(this, UPDATE_TET_GRAD_SHARED, max_tet_color_num - 1);
	thread->assignTask(this, UPDATE_TET_GRAD_SHARED_COLLISION_NEIGHBOR, max_tet_color_num - 1);
	thread->assignTask(this, SUM_ALL_GRAD);

	initialRecordPositionForThread();

	//newtonCDTetBlockAGroupCollision(max_tet_color_num - 1);
	thread->assignTask(this, SOLVE_TET_BLOCK_COLLISION, max_tet_color_num-1);

	for (int i = 0; i < total_thread_num; ++i) {
		if (max_displacement < max_dis_record[i]) {
			max_displacement = max_dis_record[i];
		}
	}

	lineSearchLastColor();

	// 
	//for (int i = 0; i < time.size(); ++i) {
	//	std::cout << time[i] << " ";
	//}
	//std::cout << std::endl;

	//system("pause");
	//newtonCDTetBlockAGroupTest();

}


void XPBD_IPC::lineSearchFirstColor(unsigned int color_No, double ori_energy)
{
	//compute energy line search
	previousColorCollisionInversionTime(color_No);

	//std::cout << "color collision time " << color_No << " " << collision.collision_time << std::endl;

	double current_energy = 0.0;
	current_energy = computePreviousColorEnergy(color_No);
	if (current_energy > energy_converge_standard) {

		double record_collision_time = collision.collision_time;
		collision.collision_time = 0.5;
		//ori_energy += 1e-15;
		while (current_energy > ori_energy)
		{
			record_collision_time *= 0.5;
			thread->assignTask(&collision, UPDATE_COLOR_POSITION, color_No);
			if (record_collision_time < 1e-5) {
				std::cout << "previous color very close " << color_No << std::endl;

				//collision.collision_time = 0.0;
				//thread->assignTask(&collision, UPDATE_COLOR_POSITION, color_No);


				//for (int k = 0; k < mesh_struct[0]->vertex_for_render.size(); ++k) {
				//	for (int kk = 0; kk < 3; ++kk) {
				//		if (abs(record_vertex_position[0][k][kk] - vertex_position[0][k][kk])>1e-6) {
				//			std::cout << abs(record_vertex_position[0][k][kk] - vertex_position[0][k][kk])<<" "<<k << std::endl;
				//		}		
				//	}
				//}
				break;
			}
			current_energy = computePreviousColorEnergy(color_No);

		}
		//std::cout << "actual collision time " << color_No << " " << record_collision_time << std::endl;
	}
	//else {
	//	std::cout << "actual collision time " << color_No << " " << collision.collision_time << std::endl;
	//}
	
	thread->assignTask(&collision, UPDATE_RECORD_POSITION_PREVIOUS_COLOR, color_No);
	//collision.upodateRecordPositionFirstColor(color_No);
}

void XPBD_IPC::lineSearchLastColor()
{
	double ori_energy = 0.0, current_energy = 0.0;

	ori_energy = 0.0, current_energy = 0.0;
	thread->assignTask(this, UPDATE_LAST_COLOR_VERTEX_BELONG);
	if (perform_collision) {
		ori_energy = computeLastColorEnergy();
	}
	thread->assignTask(this, UPDATE_POSITION_AVERAGE);



	if (perform_collision) {
		//if (collision.collision_time > 0.0) {
		allPairCollisionInversionTime();

		//std::cout << "last color collision time " << collision.collision_time << std::endl;

		current_energy = computeLastColorEnergy();
		if (current_energy > energy_converge_standard) {
			double record_collision_time = collision.collision_time;
			collision.collision_time = 0.5;
			//line search, make sure energy is decrease
			//ori_energy += 1e-15;
			while (current_energy > ori_energy)
			{
				record_collision_time *= 0.5;
				thread->assignTask(&collision, COLLISION_FREE_POSITION_LAST_COLOR);
				if (record_collision_time < 1e-5) {
					std::cout << "last color record_collision_time is too small " << ori_energy << " " << current_energy << std::endl;

					//collision.collision_time = 0.0;
					//thread->assignTask(&collision, COLLISION_FREE_POSITION_LAST_COLOR);

					break;
				}
				current_energy = computeLastColorEnergy();
			}
			//std::cout << "final color actual collision time " << record_collision_time << std::endl;
		}		
		//else {
		//	std::cout << "final color actual collision time " << collision.collision_time << std::endl;
		//}
		//update record position
		thread->assignTask(&collision, UPDATE_RECORD_VERTEX_POSITION);
		//}	
	}
}

double XPBD_IPC::computeInertialEnergyWarmStart()
{
	double energy = 0.0;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* sn_;
	unsigned int vertex_end;
	double* mass;
	int* record_position_num;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		sn_ = sn[i].data();
		vertex_end = vertex_index_begin_per_thread[i][total_thread_num];
		mass = mesh_struct[i]->mass.data();
		record_position_num = record_vertex_position_num_every_thread[i][0].data();
		for (unsigned int j = 0; j < vertex_end; ++j) {
			if (record_position_num[j]) {
				energy += mass[j] * (EDGE_LENGTH(vertex_pos[j], sn_[j]));
			}		
		}
	}
	return 0.5 * energy / (sub_time_step * sub_time_step);
}


//UPDATE_LAST_COLOR_VERTEX_BELONG
void XPBD_IPC::lastColorVertexBelongToGroup(int thread_No)
{
	int start, end;
	std::vector<std::array<double, 3>>* record_position;
	std::vector<int>* record_position_num;

	int vertex_num;
	double pos[3];
	bool* belong;
	for (int i = 0; i < total_obj_num; ++i) {
		start = vertex_index_begin_per_thread[i][thread_No];
		end = vertex_index_begin_per_thread[i][thread_No + 1];

		record_position_num = record_vertex_position_num_every_thread[i].data();
		belong = collision.vertex_belong_to_color_group[i];
		for (int j = start; j < end; ++j) {
			memset(pos, 0, 24);
			vertex_num = 0;
			for (int k = 0; k < total_thread_num; ++k) {
				if (record_position_num[k][j]) {
					belong[k] = true;
					break;
				}
			}
		}
	}
}

//UPDATE_POSITION_AVERAGE
void XPBD_IPC::updatePositionAverage(int thread_No)
{
	int start, end;
	std::vector<std::array<double, 3>>* record_position;
	std::vector<int>* record_position_num;

	int vertex_num;
	double pos[3];

	std::array<double, 3>* vertex_pos;

	for (int i = 0; i < total_obj_num; ++i) {
		start = vertex_index_begin_per_thread[i][thread_No];
		end = vertex_index_begin_per_thread[i][thread_No + 1];

		record_position = record_vertex_position_every_thread[i].data();
		record_position_num = record_vertex_position_num_every_thread[i].data();

		vertex_pos = vertex_position[i];

		for (int j = start; j < end; ++j) {
			memset(pos, 0, 24);
			vertex_num = 0;
			for (int k = 0; k < total_thread_num; ++k) {
				if (record_position_num[k][j]) {
					vertex_num += record_position_num[k][j];
					SUM_(pos, record_position[k][j]);
				}
			}

			if (vertex_num) {
				DEV(vertex_pos[j], pos, (double)vertex_num);
				record_position_num[0][j] = vertex_num;
			}
		}
	}
}


void XPBD_IPC::checkExceedFloor()
{
	for (int i = 0; i < total_obj_num; ++i) {
		for (int j = 0; j < vertex_index_begin_per_thread[i][total_thread_num]; ++j) {
			if(vertex_position[i][j][floor->dimension]<floor->value+collision.tolerance){
				std::cout << "error vertex" << j << std::endl;
			}
		}
	}

}


//WARM_START
void XPBD_IPC::solveBlockForWarmStart(int thread_No)
{
	std::array<double, 3>** record_vertex_position = record_vertex_by_thread[thread_No].data();
	int** record_vertex_num = record_vertex_update_num_by_thread[thread_No].data();
	solveVT_BlockPerThread(record_vertex_position, record_vertex_num, collision.record_vt_pair_sum_all_thread.data(), collision.vt_pair_sum_start_per_thread[thread_No],
		collision.vt_pair_sum_start_per_thread[thread_No + 1],
		collision.vt_hessian_index_exist.data(), true, thread_No);
	solveEE_BlockPerThread(record_vertex_position, record_vertex_num, collision.record_ee_pair_sum_all_thread.data(), collision.ee_pair_sum_start_per_thread[thread_No], 
		collision.ee_pair_sum_start_per_thread[thread_No+1], collision.ee_hessian_index_exist.data(), true, thread_No);

	if (has_collider || floor->exist) {
		solvecollider_BlockPerThread(record_vertex_position, record_vertex_num, &collision.vertex_index_collide_collider_sum, true, 0,
			collision.vertex_c_start_per_thread[thread_No], collision.vertex_c_start_per_thread[thread_No + 1], thread_No);
	}
	if (has_collider){
		solvecollider_BlockPerThread(record_vertex_position, record_vertex_num, &collision.edge_index_collide_collider_sum, true, 1,
			collision.edge_c_start_per_thread[thread_No], collision.edge_c_start_per_thread[thread_No + 1], thread_No);
		solvecollider_BlockPerThread(record_vertex_position, record_vertex_num, &collision.triangle_index_collide_collider_sum, true, 2,
			collision.triangle_c_start_per_thread[thread_No], collision.triangle_c_start_per_thread[thread_No + 1], thread_No);
	}



}


//WARM START
void XPBD_IPC::warmStart()
{
	//collision.initialVertexBelongColorGroup();
	unsigned itr_num = 0;
	double initial_energy = 0.0;

	previous_energy = 1e-15;
	energy = previous_energy;

	while (!convCondition(itr_num,4,energy,previous_energy,100, energy_converge_standard, energy_converge_ratio))
	{
		initialRecordHessian();
		collision.computeHessian(max_tet_color_num - 1);
		thread->assignTask(this, SUM_ALL_GRAD);
		//collision_compare.computeHessian(max_tet_color_num - 1);

		//compareIfRecordHessianIsRight(-2);

		initialRecordPositionForThread();

		
		//all self collision pair
		thread->assignTask(this, WARM_START);

		for (unsigned int i = 0; i < total_obj_num; ++i) {
			memcpy(vertex_position[i][0].data(), sn[i][0].data(), sn[i].size() * 24);
		}

		thread->assignTask(this, UPDATE_POSITION_AVERAGE);
		allPairCollisionTimeWarmStart();
		for (unsigned int i = 0; i < total_obj_num; ++i) {
			memcpy(this->record_vertex_position[i][0].data(), vertex_position[i][0].data(), this->record_vertex_position[i].size() * 24);
		}
		previous_energy=energy;
		energy = computeWarmStartEnergy();

		itr_num++;
	}

	//std::cout << "warm start itr " << itr_num << std::endl;

	collision.collision_time = 1.0;
	thread->assignTask(this, INVERSION_TEST);
	//inversionTest();
	collision.obtainCollisionTimeFromEveryThread();
	std::cout << "inversion test time " << collision.collision_time << std::endl;
	thread->assignTask(this, COLLISION_FREE_POSITION_);
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(this->record_vertex_position[i][0].data(), vertex_position[i][0].data(), this->record_vertex_position[i].size() * 24);
	}
}


void XPBD_IPC::previousColorCollisionInversionTime(unsigned int color_No)
{
	
	collision.collisionTimeColor(color_No);
	thread->assignTask(this, PREVIOUS_COLOR_INVERSION_TEST,color_No);
	//computePreviousColorInversion(color_No);
	for (int i = 0; i < total_thread_num; ++i) {
		if (collision.collision_time_thread[i] < collision.collision_time) {
			collision.collision_time = collision.collision_time_thread[i];
		}
	}

	if (collision.collision_time < 1.0) {
		thread->assignTask(&collision, UPDATE_COLOR_POSITION, color_No);
	}
}




void XPBD_IPC::allPairCollisionInversionTime()
{
	collision.allPairCollisionTime();

	thread->assignTask(this, LAST_COLOR_INVERSION_TEST);
	for (int i = 0; i < total_thread_num; ++i) {
		if (collision.collision_time_thread[i] < collision.collision_time) {
			collision.collision_time = collision.collision_time_thread[i];
		}
	}
	//computeLastColorInversion();
	if (collision.collision_time < 1.0) {
		thread->assignTask(&collision, COLLISION_FREE_POSITION_LAST_COLOR);
	}
}

void XPBD_IPC::allPairCollisionTimeWarmStart()
{
	collision.allPairCollisionTime();
	std::cout << "warm start collision time " << collision.collision_time << std::endl;
	if (collision.collision_time < 1.0) {
		thread->assignTask(this, COLLISION_FREE_POSITION_FROM_RECORD);
	}
}




//SOLVE_TET_BLOCK_COLLISION
void XPBD_IPC::newtonCDTetBlockAGroupCollision(int thread_No, int color)
{
	unsigned int* tet_group;
	int start, end;

	std::array<int, 4>* indices;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;

	double stiffness;
	double collision_stiffness;;
	Matrix<double, 3, 4>* A;
	std::array<double, 3>* sn_;
	double* mass;
	unsigned int* unfixed_vertex_num;
	std::vector<unsigned int>* tet_neightbor_tet;
	std::vector<unsigned int>* tet_neightbor_tet_common_vertex;
	std::array<int, 4>* unfixed_vertex_index;
	std::array<int, 4>* unfixed_actual_vertex_index;


	std::vector<unsigned int>* triangle_of_a_tet;
	std::vector<unsigned int>* edge_of_a_tet;
	unsigned int obj_No;
	int* vertex_index_on_surface;
	int j;


	int color_group_index;

	std::array<double, 3>* vertex_pos_record;
	std::array<double, 3>* free_pos;
	int* vertex_pos_num_record;

	unsigned int prefix_vertex;

	double max_dis = 0.0;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		color_group_index = inner_iteration_number % tet_color_groups[i]->size();
		if (color != tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index].size() - 1) {
			continue;
		}

		collision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
		obj_No = i + cloth->size();

		mesh_struct_ = mesh_struct[i + cloth->size()];
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i + cloth->size()];
		stiffness = 0.5 * tetrahedron->data()[i].ARAP_stiffness;

		A = tetrahedron->data()[i].mesh_struct.A.data();
		mass = mesh_struct_->mass.data();
		sn_ = sn[i + cloth->size()].data();
		unfixed_vertex_num = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();
		tet_neightbor_tet = tetrahedron->data()[i].mesh_struct.tet_tet_index.data();
		tet_neightbor_tet_common_vertex = tetrahedron->data()[i].mesh_struct.tet_neighbor_tet_vertex_order.data();
		unfixed_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_indices.data();
		unfixed_actual_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_actual_indices.data();
		triangle_of_a_tet = tetrahedron->data()[i].mesh_struct.triangle_index_of_a_tet.data();
		edge_of_a_tet = tetrahedron->data()[i].mesh_struct.edge_index_of_a_tet.data();

		vertex_index_on_surface = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();

		vertex_pos_record = record_vertex_position_every_thread[obj_No][thread_No].data();
		vertex_pos_num_record = record_vertex_position_num_every_thread[obj_No][thread_No].data();

		free_pos = record_collision_free_vertex_position[obj_No].data();

		//the original tet in last color
		end = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color][thread_No + 1];
		start = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color][thread_No];
		tet_group = tet_color_groups[i]->data()[color_group_index][color].data();

		prefix_vertex = vertex_index_prefix_sum_obj[obj_No];

		for (int k = start; k < end; ++k) {
			j = tet_group[k];
			solveTetBlockCollision(vertex_pos, stiffness, sub_time_step, tet_indices[i], mass, A, tet_neightbor_tet[j],
				volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
				unfixed_vertex_index[j].data(), unfixed_vertex_num[j], &(triangle_of_a_tet[j]), &(edge_of_a_tet[j]), collision_stiffness,
				obj_No, unfixed_actual_vertex_index[j].data(), vertex_index_on_surface, tet_hessian[obj_No], common_hessian,
				common_grad[0].data(), vertex_pos_record, vertex_pos_num_record, prefix_vertex, floor_hessian.data(),
				max_dis, free_pos);
		}
	}

	if (max_dis > max_dis_record[thread_No]) {
		max_dis_record[thread_No] = max_dis;
	}

	if (perform_collision) {
		int obj_No_0, obj_No_1, index_0, index_1;
		std::array<double, 3>** record_vertex_position = record_vertex_by_thread[thread_No].data();
		int** record_vertex_num = record_vertex_update_num_by_thread[thread_No].data();
		//all self collision pair
		solveVT_BlockPerThread(record_vertex_position, record_vertex_num, collision.record_vt_pair_sum_all_thread.data(), collision.vt_pair_sum_start_per_thread[thread_No],
			collision.vt_pair_sum_start_per_thread[thread_No + 1],
			collision.vt_hessian_index_exist.data(),false, thread_No);
		solveEE_BlockPerThread(record_vertex_position, record_vertex_num, collision.record_ee_pair_sum_all_thread.data(), collision.ee_pair_sum_start_per_thread[thread_No],
			collision.ee_pair_sum_start_per_thread[thread_No + 1],
			collision.ee_hessian_index_exist.data(),false, thread_No);
		if (has_collider || floor->exist) {
			solvecollider_BlockPerThread(record_vertex_position, record_vertex_num, &collision.vertex_index_collide_collider_sum, false, 0,
				collision.vertex_c_start_per_thread[thread_No], collision.vertex_c_start_per_thread[thread_No + 1],thread_No);
		}
		if (has_collider) {
			solvecollider_BlockPerThread(record_vertex_position, record_vertex_num, &collision.edge_index_collide_collider_sum, false, 1,
				collision.edge_c_start_per_thread[thread_No], collision.edge_c_start_per_thread[thread_No + 1],thread_No);
			solvecollider_BlockPerThread(record_vertex_position, record_vertex_num, &collision.triangle_index_collide_collider_sum, false, 2,
				collision.triangle_c_start_per_thread[thread_No], collision.triangle_c_start_per_thread[thread_No + 1],thread_No);
		}
	}
}




// 0v, 1e, 2t
void XPBD_IPC::solvecollider_BlockPerThread(std::array<double, 3>** record_vertex_position, int** record_vertex_num, std::vector<unsigned int>* pair, 
	bool only_solve_collision_pair, int type, int start, int end, int thread_No)
{
	auto start_ = pair->begin() + start;
	auto end_ = pair->begin() + end;
	for (auto j = start_; j < end_; j += 2){

			solveCollisionWithColliderBlock(sub_time_step, *j, *(j+1), type, record_vertex_position, record_vertex_num, only_solve_collision_pair,
				thread_No);
		
	}
}





void XPBD_IPC::solveVT_BlockPerThread(std::array<double, 3>** record_vertex_position, int** record_vertex_num, unsigned int* pair, unsigned int start, unsigned int end,
	char* vt_hessian_index_exist, bool only_solve_collision_pair, int thread_No)
{
	int obj_No_0, obj_No_1, index_0, index_1;
	char* record_index;
	for (int k = start; k < end; k += 4) {
		obj_No_0 = pair[k];
		index_0 = pair[k + 1];
		obj_No_1 = pair[k + 2];
		index_1 = pair[k + 3];


		record_index = vt_hessian_index_exist + 5 * (k >> 2);
		if (!(*record_index)) {
			continue;
		}
		solveVT_Block(obj_No_0, index_0,
			obj_No_1, index_1, sub_time_step,
			record_vertex_position, record_vertex_num, only_solve_collision_pair, thread_No);
	}
}


void XPBD_IPC::solveEE_BlockPerThread(std::array<double, 3>** record_vertex_position, int** record_vertex_num, unsigned int* pair, unsigned int start, unsigned int end,
	char* ee_hessian_record_index_exist, bool only_solve_collision_pair, int thread_No)
{
	int obj_No_0, obj_No_1, index_0, index_1;
	char* record_index;
	for (int k = start; k < end; k += 4) {
		obj_No_0 = pair[k];
		index_0 = pair[k + 1];
		obj_No_1 = pair[k + 2];
		index_1 = pair[k + 3];
		record_index = ee_hessian_record_index_exist + 5 * (k >> 2);
		if (!(* record_index)) {
			continue;
		}
		solveEE_Block(obj_No_0, index_0,
			obj_No_1, index_1, sub_time_step,
			record_vertex_position, record_vertex_num, only_solve_collision_pair, thread_No);
	}
}


//SOLVE_TET_BLOCK
void XPBD_IPC::newtonCDTetBlockAGroup(int thread_No, int color)
{

	unsigned int* tet_group;
	int start, end;

	std::array<int, 4>* indices;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* record_vertex_pos;

	double stiffness;
	double collision_stiffness;;
	Matrix<double, 3, 4>* A;
	std::array<double, 3>* sn_;
	double* mass;
	double* mass_inv;
	unsigned int* unfixed_vertex_num;
	std::vector<unsigned int>* tet_neightbor_tet;
	std::vector<unsigned int>* tet_neightbor_tet_common_vertex;
	std::array<int, 4>* unfixed_vertex_index;
	std::array<int, 4>* unfixed_actual_vertex_index;
	unsigned int obj_No;

	int j;

	int color_group_index;


	//char* is_tet_involved_in_collision;

	std::vector<unsigned int>* triangle_of_a_tet;
	std::vector<unsigned int>* edge_of_a_tet;

	int* vertex_index_on_surface;

	unsigned int prefix_vertex;
	char* indicate_collide_with_floor_;
	for (int i = 0; i < tetrahedron->size(); ++i) {
		color_group_index = inner_iteration_number % tet_color_groups[i]->size();
		if (color >= tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index].size()-1) {
			continue;
		}

		collision_stiffness = tetrahedron->data()[i].collision_stiffness[0];

		tet_group = tet_color_groups[i]->data()[color_group_index][color].data();
		//is_tet_involved_in_collision = tet_color_groups_label[i][color_group_index][color].data();

		end = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color][thread_No + 1];
		start = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color][thread_No];

		mesh_struct_ = mesh_struct[i + cloth->size()];
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i + cloth->size()];
		stiffness = 0.5 * tetrahedron->data()[i].ARAP_stiffness;

		A = tetrahedron->data()[i].mesh_struct.A.data();
		mass = mesh_struct_->mass.data();
		mass_inv = mesh_struct_->mass_inv.data();
		sn_ = sn[i + cloth->size()].data();
		unfixed_vertex_num = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();
		tet_neightbor_tet = tetrahedron->data()[i].mesh_struct.tet_tet_index.data();
		tet_neightbor_tet_common_vertex = tetrahedron->data()[i].mesh_struct.tet_neighbor_tet_vertex_order.data();
		unfixed_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_indices.data();
		unfixed_actual_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_actual_indices.data();
		obj_No = i + cloth->size();

		triangle_of_a_tet = tetrahedron->data()[i].mesh_struct.triangle_index_of_a_tet.data();
		edge_of_a_tet = tetrahedron->data()[i].mesh_struct.edge_index_of_a_tet.data();

		vertex_index_on_surface = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();

		prefix_vertex = vertex_index_prefix_sum_obj[obj_No];
		record_vertex_pos = record_vertex_position[obj_No].data();

		indicate_collide_with_floor_ = collision.indicate_vertex_collide_with_floor[obj_No].data();

		for (int k = start; k < end; ++k) {
			j = tet_group[k];
			solveTetBlock(vertex_pos, stiffness, sub_time_step, mass, mass_inv, A, tet_neightbor_tet[j],
				volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
				unfixed_vertex_index[j].data(), unfixed_vertex_num[j],  collision_stiffness,
				obj_No, unfixed_actual_vertex_index[j].data(), tet_hessian[obj_No], common_hessian,
				common_grad[0].data(), &(triangle_of_a_tet[j]), &(edge_of_a_tet[j]), vertex_index_on_surface, prefix_vertex, floor_hessian.data(),
				record_vertex_pos, indicate_collide_with_floor_,color);
		}
	}
}


//void XPBD_IPC::newtonCDTetBlockAGroupTest(int  thread_No, int color)
//{
//	unsigned int* tet_group;
//	int start, end;
//
//	std::array<int, 4>* indices;
//	MeshStruct* mesh_struct_;
//	double* volume;
//	std::array<double, 3>* vertex_pos;
//
//	double stiffness;
//	double collision_stiffness;;
//	Matrix<double, 3, 4>* A;
//	std::array<double, 3>* sn_;
//	double* mass;
//	unsigned int* unfixed_vertex_num;
//	std::vector<unsigned int>* tet_neightbor_tet;
//	std::vector<unsigned int>* tet_neightbor_tet_common_vertex;
//	std::array<int, 4>* unfixed_vertex_index;
//	std::array<int, 4>* unfixed_actual_vertex_index;
//
//
//	std::vector<unsigned int>* triangle_of_a_tet;
//	std::vector<unsigned int>* edge_of_a_tet;
//	unsigned int obj_No;
//	int* vertex_index_on_surface;
//
//	std::array<double, 3>* ori_vertex_pos;
//
//	int end_index;
//
//	for (int i = 0; i < tetrahedron->size(); ++i) {
//
//		collision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
//
//		mesh_struct_ = mesh_struct[i + cloth->size()];
//		indices = tetrahedron->data()[i].mesh_struct.indices.data();
//		volume = tetrahedron->data()[i].mesh_struct.volume.data();
//		vertex_pos = vertex_position[i + cloth->size()];
//		stiffness = 0.5 * tetrahedron->data()[i].ARAP_stiffness;
//
//		A = tetrahedron->data()[i].mesh_struct.A.data();
//		mass = mesh_struct_->mass.data();
//		sn_ = sn[i + cloth->size()].data();
//		unfixed_vertex_num = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();
//		tet_neightbor_tet = tetrahedron->data()[i].mesh_struct.tet_tet_index.data();
//		tet_neightbor_tet_common_vertex = tetrahedron->data()[i].mesh_struct.tet_neighbor_tet_vertex_order.data();
//		unfixed_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_indices.data();
//		unfixed_actual_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_actual_indices.data();
//		triangle_of_a_tet = tetrahedron->data()[i].mesh_struct.triangle_index_of_a_tet.data();
//		edge_of_a_tet = tetrahedron->data()[i].mesh_struct.edge_index_of_a_tet.data();
//		obj_No = i + cloth->size();
//		vertex_index_on_surface = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
//		ori_vertex_pos = initial_vertex_position[i + cloth->size()];
//		int j;
//		end_index = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread[color][thread_No + 1];
//		for (int m = 0; m < end_index; ++m) {
//			j = tetrahedron->data()[i].mesh_struct.unconnected_tet_index[color][m];
//			solveTetBlock(vertex_pos, stiffness, sub_time_step, mass, A, tet_neightbor_tet[j],
//				volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
//				unfixed_vertex_index[j].data(), unfixed_vertex_num[j], &(triangle_of_a_tet[j]), &(edge_of_a_tet[j]), collision_stiffness,
//				obj_No, unfixed_actual_vertex_index[j].data(), vertex_index_on_surface, store_tet_arap_hessian.data() + prefix_sum_of_every_tet_index[i] * 16, 
//				store_tet_arap_grad.data() + prefix_sum_of_every_tet_index[i] * 12);//
//			//
//		}
//
//
//
//
//		//for (unsigned int j = 0; j < tetrahedron->data()[i].mesh_struct.indices.size(); ++j) {
//		//	solveNewtonCD_tetBlock(vertex_pos, stiffness, sub_time_step, mass, A, tet_neightbor_tet[j],
//		//		indices, volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
//		//		unfixed_vertex_index[j].data(), unfixed_vertex_num[j], &(triangle_of_a_tet[j]), &(edge_of_a_tet[j]), collision_stiffness,
//		//		obj_No, unfixed_actual_vertex_index[j].data(), vertex_index_on_surface, ori_vertex_pos);
//		//}
//	}
//}

//void XPBD_IPC::newtonCDTetBlockTest(int color_No)
//{
//	unsigned int size;
//	std::array<int, 4>* indices;
//	MeshStruct* mesh_struct_;
//	double* volume;
//	std::array<double, 3>* vertex_pos;
//	std::array<double, 3>* ori_vertex_pos;
//	double* mass_inv;
//	double stiffness;
//	double collision_stiffness;;
//	Matrix<double, 3, 4>* A;
//	std::array<double, 3>* sn_;
//	double* mass;
//	unsigned int* unfixed_vertex_num;
//	std::vector<unsigned int>* tet_neightbor_tet;
//	std::vector<unsigned int>* tet_neightbor_tet_common_vertex;
//	std::array<int, 4>* unfixed_vertex_index;
//	std::array<int, 4>* unfixed_actual_vertex_index;
//
//
//	std::vector<unsigned int>* triangle_of_a_tet;
//	std::vector<unsigned int>* edge_of_a_tet;
//	unsigned int obj_No;
//	int* vertex_index_on_surface;
//
//
//	int j;
//
//	unsigned int* color;
//	int color_group_index;
//	for (int i = 0; i < tetrahedron->size(); ++i) {
//
//		color_group_index = inner_iteration_number % tet_color_groups[i]->size();
//
//		collision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
//		mesh_struct_ = mesh_struct[i + cloth->size()];
//		size = tetrahedron->data()[i].mesh_struct.indices.size();
//		indices = tetrahedron->data()[i].mesh_struct.indices.data();
//		volume = tetrahedron->data()[i].mesh_struct.volume.data();
//		vertex_pos = vertex_position[i + cloth->size()];
//		stiffness = 0.5 * tetrahedron->data()[i].ARAP_stiffness;
//
//		A = tetrahedron->data()[i].mesh_struct.A.data();
//		mass_inv = mesh_struct_->mass_inv.data();
//		mass = mesh_struct_->mass.data();
//		sn_ = sn[i + cloth->size()].data();
//		unfixed_vertex_num = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();
//		tet_neightbor_tet = tetrahedron->data()[i].mesh_struct.tet_tet_index.data();
//		tet_neightbor_tet_common_vertex = tetrahedron->data()[i].mesh_struct.tet_neighbor_tet_vertex_order.data();
//		unfixed_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_indices.data();
//		unfixed_actual_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_actual_indices.data();
//		triangle_of_a_tet = tetrahedron->data()[i].mesh_struct.triangle_index_of_a_tet.data();
//		edge_of_a_tet = tetrahedron->data()[i].mesh_struct.edge_index_of_a_tet.data();
//		obj_No = i + cloth->size();
//		vertex_index_on_surface = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
//		ori_vertex_pos = record_vertex_position[i + cloth->size()].data();
//
//		color = tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index][color_No].data();
//
//
//		if (color_No >= tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index].size() - 1) {
//			continue;
//		}
//
//		for (unsigned int k = 0; k < tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index][color_No].size(); ++k) {
//			j = color[k];
//			solveNewtonCD_tetBlock(vertex_pos, stiffness, sub_time_step, mass, A, tet_neightbor_tet[j],
//				indices, volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
//				unfixed_vertex_index[j].data(), unfixed_vertex_num[j], &(triangle_of_a_tet[j]), &(edge_of_a_tet[j]), collision_stiffness,
//				obj_No, unfixed_actual_vertex_index[j].data(), vertex_index_on_surface, ori_vertex_pos,
//				store_tet_arap_hessian.data()+16*prefix_sum_of_every_tet_index[i]);
//		}
//	}
//}
//



void XPBD_IPC::newtonCDTetBlock()
{
	//unsigned int size;
	//std::array<int, 4>* indices;
	//MeshStruct* mesh_struct_;
	//double* volume;
	//std::array<double, 3>* vertex_pos;
	//std::array<double, 3>* ori_vertex_pos;
	//double* mass_inv;
	//double stiffness;
	//double collision_stiffness;;
	//Matrix<double, 3, 4>* A;
	//std::array<double, 3>* sn_;
	//double* mass;
	//unsigned int* unfixed_vertex_num;
	//std::vector<unsigned int>* tet_neightbor_tet;
	//std::vector<unsigned int>* tet_neightbor_tet_common_vertex;
	//std::array<int,4>* unfixed_vertex_index;
	//std::array<int,4>* unfixed_actual_vertex_index;


	//std::vector<unsigned int>* triangle_of_a_tet;
	//std::vector<unsigned int>* edge_of_a_tet;
	//unsigned int obj_No;
	//int* vertex_index_on_surface;

	//for (int i = 0; i < tetrahedron->size(); ++i) {
	//	collision_stiffness = tetrahedron->data()[i].collision_stiffness[0];
	//	mesh_struct_ = mesh_struct[i + cloth->size()];
	//	size = tetrahedron->data()[i].mesh_struct.indices.size();
	//	indices = tetrahedron->data()[i].mesh_struct.indices.data();
	//	volume = tetrahedron->data()[i].mesh_struct.volume.data();
	//	vertex_pos = vertex_position[i + cloth->size()];
	//	stiffness =0.5* tetrahedron->data()[i].ARAP_stiffness;
	//	
	//	A = tetrahedron->data()[i].mesh_struct.A.data();
	//	mass_inv = mesh_struct_->mass_inv.data();
	//	mass = mesh_struct_->mass.data();
	//	sn_ = sn[i + cloth->size()].data();
	//	unfixed_vertex_num = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();
	//	tet_neightbor_tet = tetrahedron->data()[i].mesh_struct.tet_tet_index.data();
	//	tet_neightbor_tet_common_vertex= tetrahedron->data()[i].mesh_struct.tet_neighbor_tet_vertex_order.data();
	//	unfixed_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_indices.data();
	//	unfixed_actual_vertex_index = tetrahedron->data()[i].mesh_struct.unfixied_actual_indices.data();
	//	triangle_of_a_tet = tetrahedron->data()[i].mesh_struct.triangle_index_of_a_tet.data();
	//	edge_of_a_tet = tetrahedron->data()[i].mesh_struct.edge_index_of_a_tet.data();
	//	obj_No = i + cloth->size();
	//	vertex_index_on_surface = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();
	//	ori_vertex_pos = record_vertex_position[i + cloth->size()].data();
	//	for (unsigned int j = 0; j < size; ++j) {
	//		solveNewtonCD_tetBlock(vertex_pos, stiffness, sub_time_step, mass, A, tet_neightbor_tet[j],
	//			indices, volume, j, sn_, tet_neightbor_tet_common_vertex[j].data(), indices[j].data(),
	//			unfixed_vertex_index[j].data(), unfixed_vertex_num[j],&(triangle_of_a_tet[j]),&(edge_of_a_tet[j]), collision_stiffness,
	//			obj_No, unfixed_actual_vertex_index[j].data(), vertex_index_on_surface, ori_vertex_pos,
	//			store_tet_arap_hessian.data()+16* prefix_sum_of_every_tet_index[i]);
	//	}
	//}
}



void XPBD_IPC::newtonCDTet()
{
	unsigned int size;
	std::array<int, 4>* indices;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;
	double* mass_inv;
	double stiffness;
	Matrix<double, 3, 4>* A;
	std::array<double, 3>* sn_;
	double* mass;
	double* lambda_= lambda.data() + constraint_index_start[2];
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		mesh_struct_ = mesh_struct[i + cloth->size()];
		size = tetrahedron->data()[i].mesh_struct.vertex_position.size();
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i + cloth->size()];
		stiffness = tetrahedron->data()[i].ARAP_stiffness;
		A = tetrahedron->data()[i].mesh_struct.A.data();
		mass_inv = mesh_struct_->mass_inv.data();
		mass = mesh_struct_->mass.data();
		sn_ = sn[i + cloth->size()].data();
		for (unsigned int j = 0; j < size; ++j) {
			if (mass_inv[j] != 0.0) {
				solveNewtonCD_tet(vertex_pos, stiffness, sub_time_step, A,
					mesh_struct_->vertex_tet_index[j], indices, mass, volume, j, sn_,lambda_);
			}
		}
		lambda_ += tetrahedron->data()[i].mesh_struct.indices.size();
	}

}


bool XPBD_IPC::getFloorHessian(double& Hessian, double& grad, double* vertex_position, double floor_value,
	double* last_step_position, unsigned int dimension, double collision_stiffness, bool direction, double tolerance)
{
	if(direction){
		if (vertex_position[dimension] > floor_value+ tolerance) {
			return false;
		}
	}
	else {
		if (vertex_position[dimension] < floor_value- tolerance) {
			return false;
		}
	}
	double d_current = vertex_position[dimension] - floor_value;
	double d_ori = abs(last_step_position[dimension] - floor_value);// 
	if (d_ori < tolerance) {
		d_ori = tolerance;
	}
	double ln_;

	if (d_current < 0) {
		d_current = -d_current;
		ln_ = log(d_current / d_ori);
		grad =-collision_stiffness * (d_ori - d_current) * (2.0 * ln_ - d_ori / d_current + 1.0);
		Hessian =- collision_stiffness *(-2.0 * ln_ + (d_ori - d_current) * (d_ori + 3.0 * d_current) / (d_current * d_current));
	}
	else {
		ln_ = log(d_current / d_ori);
		grad = collision_stiffness * (d_ori - d_current) * (2.0 * ln_ - d_ori / d_current + 1.0);
		Hessian = collision_stiffness * (-2.0 * ln_ + (d_ori - d_current) * (d_ori + 3.0 * d_current) / (d_current * d_current));
	}
	return true;
}


double XPBD_IPC::getCollisionTime(std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet,
	int* tet_actual_unfixed_vertex_indices,
	int unfixed_tet_vertex_num,
	int** vertex_index_on_surface,
	std::array<double, 3>** current_vertex_position,
	std::array<double, 3>** initial_vertex_position)
{
	double collision_time = 1.0;
	int* triangle_;		unsigned int* edge_;

	for (int i = 0; i < unfixed_tet_vertex_num; i+=2) {
		if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i + 1]] == -1) {
			continue;
		}
		collision.VTCollisionTimeOneVertex(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data(),
			current_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data(), collision_time,
			collision.vertex_triangle_pair_num_record[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i + 1]],
			collision.vertex_triangle_pair_by_vertex[tet_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * tet_actual_unfixed_vertex_indices[i + 1],
			initial_vertex_position, vertex_position.data(), triangle_indices.data(),false, tet_actual_unfixed_vertex_indices[i + 1]);

		if (has_collider) {
			collision.VTCollisionTimeOneVertex(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data(),
				current_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data(), collision_time,
				collision.vertex_obj_triangle_collider_num_record[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_obj_triangle_collider_pair_by_vertex[tet_actual_unfixed_vertex_indices[i]] + 
				collision.close_vt_collider_pair_num * tet_actual_unfixed_vertex_indices[i + 1],
				vertex_position_collider.data(), vertex_position_collider.data(), triangle_indices_collider.data(),true, tet_actual_unfixed_vertex_indices[i + 1]);
		}
		if (floor->exist) {
			collision.floorCollisionTime(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data()[floor->dimension],
				current_vertex_position[tet_actual_unfixed_vertex_indices[i]][tet_actual_unfixed_vertex_indices[i+1]].data()[floor->dimension],
				floor->normal_direction, floor->value, collision_time, collision.tolerance);
		}
	}

	for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); i+=2) {
		triangle_ = triangle_indices[*i][*(i+1)].data();
		collision.TVCollisionTimeOneTriangle(initial_vertex_position[*i][triangle_[0]].data(), initial_vertex_position[*i][triangle_[1]].data(),
			initial_vertex_position[*i][triangle_[2]].data(),
			current_vertex_position[*i][triangle_[0]].data(), current_vertex_position[*i][triangle_[1]].data(),
			current_vertex_position[*i][triangle_[2]].data(), collision_time,
			collision.triangle_vertex_pair_num_record[*i][*(i+1)],
			collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i+1)),
			initial_vertex_position, vertex_position.data());
		if (has_collider) {
			collision.TVCollisionTimeOneTriangle(initial_vertex_position[*i][triangle_[0]].data(), initial_vertex_position[*i][triangle_[1]].data(),
				initial_vertex_position[*i][triangle_[2]].data(),
				current_vertex_position[*i][triangle_[0]].data(), current_vertex_position[*i][triangle_[1]].data(),
				current_vertex_position[*i][triangle_[2]].data(), collision_time,
				collision.triangle_vertex_collider_pair_num_record[*i][*(i + 1)],
				collision.triangle_vertex_collider_pair_by_triangle[*i] + collision.close_tv_collider_pair_num * (*(i + 1)),
				vertex_position_collider.data(), vertex_position_collider.data());
		}
	}

	for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); i+=2) {
		edge_ = edge_vertices[*i] + ((*(i+1)) << 1);
		collision.EECollisionTimeOneEdgeAll(initial_vertex_position[*i][*edge_].data(),
			initial_vertex_position[*i][*(edge_ + 1)].data(),
			current_vertex_position[*i][*edge_].data(),
			current_vertex_position[*i][*(edge_ + 1)].data(), collision_time,
			collision.edge_edge_pair_number_record[*i][*(i+1)],
			collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i+1)),
			initial_vertex_position, vertex_position.data(), edge_vertices.data());
		if (has_collider) {
			collision.EECollisionTimeOneEdgeAll(initial_vertex_position[*i][*edge_].data(),
				initial_vertex_position[*i][*(edge_ + 1)].data(),
				current_vertex_position[*i][*edge_].data(),
				current_vertex_position[*i][*(edge_ + 1)].data(), collision_time,
				collision.edge_edge_collider_pair_num_record[*i][*(i + 1)],
				collision.edge_edge_collider_pair_by_edge[*i] + collision.close_ee_collider_pair_num * (*(i + 1)),
				vertex_position_collider.data(), vertex_position_collider.data(), collider_edge_vertices.data());
		}
	}
	return collision_time;
}


double XPBD_IPC::getInversionTime(unsigned int tet_index, std::vector<unsigned int>* neighbor_tet_index,
	std::array<int,4>* tet_vertex, std::array<double, 3>* current_vertex_position,
	std::array<double, 3>* initial_vertex_position)
{
	double collision_time = 1.0;
	int* vertex_index;
	double inversion_time=1.0;

	vertex_index = tet_vertex[tet_index].data();
	if (inversionTest::TetInversionTest(initial_vertex_position[vertex_index[0]].data(), initial_vertex_position[vertex_index[1]].data(),
		initial_vertex_position[vertex_index[2]].data(), initial_vertex_position[vertex_index[3]].data(),
		current_vertex_position[vertex_index[0]].data(), current_vertex_position[vertex_index[1]].data(),
		current_vertex_position[vertex_index[2]].data(), current_vertex_position[vertex_index[3]].data(), &inversion_time)) {
		if (inversion_time < collision_time) {
			collision_time = inversion_time;
		}
	}
	for (auto i = neighbor_tet_index->begin(); i < neighbor_tet_index->end(); ++i) {
		vertex_index = tet_vertex[*i].data();
		
		if (inversionTest::TetInversionTest(initial_vertex_position[vertex_index[0]].data(), initial_vertex_position[vertex_index[1]].data(),
			initial_vertex_position[vertex_index[2]].data(), initial_vertex_position[vertex_index[3]].data(),
			current_vertex_position[vertex_index[0]].data(), current_vertex_position[vertex_index[1]].data(),
			current_vertex_position[vertex_index[2]].data(), current_vertex_position[vertex_index[3]].data(), &inversion_time)) {
			if (inversion_time < collision_time) {
				collision_time= inversion_time;
			}
		}
	}

	return collision_time;
}


double XPBD_IPC::getCollisionTime(std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet,
	unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	int unfixed_tet_vertex_num,
	int* vertex_index_on_surface, std::array<double, 3>* current_vertex_position,
	std::array<double, 3>* initial_vertex_position, char* indicate_collide_with_floor)
{
	double collision_time = 1.0;
	int* triangle_;		unsigned int* edge_;

	for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
		if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
			continue;
		}
		collision.VTCollisionTimeOneVertex(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]].data(),
			current_vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), collision_time,
			collision.vertex_triangle_pair_num_record[obj_No][tet_actual_unfixed_vertex_indices[i]],
			collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * tet_actual_unfixed_vertex_indices[i],
			this->record_vertex_position_address.data(), vertex_position.data(), triangle_indices.data(),false, tet_actual_unfixed_vertex_indices[i]);
		if (has_collider) {
			collision.VTCollisionTimeOneVertex(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]].data(),
				current_vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), collision_time,
				collision.vertex_obj_triangle_collider_num_record[obj_No][tet_actual_unfixed_vertex_indices[i]],
				collision.vertex_obj_triangle_collider_pair_by_vertex[obj_No] + collision.close_vt_collider_pair_num * tet_actual_unfixed_vertex_indices[i],
				vertex_position_collider.data(), vertex_position_collider.data(), triangle_indices_collider.data(),true, tet_actual_unfixed_vertex_indices[i]);
		}
		if (floor->exist) {
			if (indicate_collide_with_floor[tet_actual_unfixed_vertex_indices[i]]) {
				collision.floorCollisionTime(initial_vertex_position[tet_actual_unfixed_vertex_indices[i]].data()[floor->dimension],
					current_vertex_position[tet_actual_unfixed_vertex_indices[i]].data()[floor->dimension],
					floor->normal_direction, floor->value, collision_time, collision.tolerance);
			}
		}
	}

	for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
		triangle_ = triangle_indices[obj_No][*i].data();
		collision.TVCollisionTimeOneTriangle(initial_vertex_position[triangle_[0]].data(), initial_vertex_position[triangle_[1]].data(),
			initial_vertex_position[triangle_[2]].data(),
			current_vertex_position[triangle_[0]].data(), current_vertex_position[triangle_[1]].data(),
			current_vertex_position[triangle_[2]].data(), collision_time,
			collision.triangle_vertex_pair_num_record[obj_No][*i],
			collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * (*i),
			this->record_vertex_position_address.data(), vertex_position.data());
		if (has_collider) {
			collision.TVCollisionTimeOneTriangle(initial_vertex_position[triangle_[0]].data(), initial_vertex_position[triangle_[1]].data(),
				initial_vertex_position[triangle_[2]].data(),
				current_vertex_position[triangle_[0]].data(), current_vertex_position[triangle_[1]].data(),
				current_vertex_position[triangle_[2]].data(), collision_time,
				collision.triangle_vertex_collider_pair_num_record[obj_No][*i],
				collision.triangle_vertex_collider_pair_by_triangle[obj_No] + collision.close_tv_collider_pair_num * (*i),
				vertex_position_collider.data(), vertex_position_collider.data());
		}
	}

	for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
		edge_ = edge_vertices[obj_No] + ((*i) << 1);
		collision.EECollisionTimeOneEdgeAll(initial_vertex_position[*edge_].data(),
			initial_vertex_position[*(edge_ + 1)].data(),
			current_vertex_position[*edge_].data(),
			current_vertex_position[*(edge_ + 1)].data(), collision_time,
			collision.edge_edge_pair_number_record[obj_No][*i],
			collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * (*i),
			this->record_vertex_position_address.data(), vertex_position.data(), edge_vertices.data());
		if (has_collider) {
			collision.EECollisionTimeOneEdgeAll(initial_vertex_position[*edge_].data(),
				initial_vertex_position[*(edge_ + 1)].data(),
				current_vertex_position[*edge_].data(),
				current_vertex_position[*(edge_ + 1)].data(), collision_time,
				collision.edge_edge_collider_pair_num_record[obj_No][*i],
				collision.edge_edge_collider_pair_by_edge[obj_No] + collision.close_ee_collider_pair_num * (*i),
				vertex_position_collider.data(), vertex_position_collider.data(), collider_edge_vertices.data());
		}
	}

	return collision_time;
}


void XPBD_IPC::getCollisionBlockTetHessian(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* tets, double stiffness, int* pair_actual_unfixed_vertex_indices, //pair_actual_unfixed_vertex_indices first obj, second primitive index
	int unfixed_vertex_num)//unfixed_vertex_num size double )
{
	for (auto i = tets->begin(); i < tets->end(); i += 2) {

		getARAPCollisionHessianForPair(Hessian, grad, stiffness, *i, *(i + 1), tet_indices[*i][*(i + 1)].data(),
			pair_actual_unfixed_vertex_indices, unfixed_vertex_num, vertex_position[*i], tet_A[*i], tet_volume[*i]);
	}
}

void XPBD_IPC::getCollisionBlockCollisionHessianTest(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangles,
	std::vector<unsigned int>* edges,
	double collision_stiffness, int* pair_actual_unfixed_vertex_indices, //pair_actual_unfixed_vertex_indices first obj, second primitive index
	int unfixed_vertex_num, double d_hat_2, int** vertex_index_on_surface, unsigned int vertex_obj_No, unsigned int tri_obj_No)//unfixed_vertex_num size double
{
	if (has_collider) {
		for (int i = 0; i < unfixed_vertex_num; i += 2) {
			getVTCollisionHessainForPairTest(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]].data(),
				collision_stiffness, collision.vertex_triangle_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_triangle_pair_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]],
				i >> 1, vertex_obj_No, tri_obj_No, pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				collision.vertex_obj_triangle_collider_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_collider_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_obj_triangle_collider_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (int i = 0; i < unfixed_vertex_num; i += 2) {
			getVTCollisionHessainForPairTest(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]].data(),
				collision_stiffness, collision.vertex_triangle_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_triangle_pair_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]],
				i >> 1, vertex_obj_No, tri_obj_No, pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				emp, 0);
		}
	}
	int* triangle_;
	if (has_collider) {
		for (auto i = triangles->begin(); i < triangles->end(); i += 2) {
			triangle_ = triangle_indices[*i][*(i + 1)].data();
			getTVCollisionHessainForPairTest(Hessian, grad, vertex_position[*i][triangle_[0]].data(), vertex_position[*i][triangle_[1]].data(),
				vertex_position[*i][triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i + 1)),
				collision.triangle_vertex_pair_num_record[*i][*(i + 1)], *i,
				vertex_obj_No, tri_obj_No, triangle_, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
				d_hat_2,
				collision.triangle_vertex_collider_pair_by_triangle[*i] + collision.close_tv_collider_pair_num * (*(i + 1)),
				collision.triangle_vertex_collider_pair_num_record[*i][*(i + 1)]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = triangles->begin(); i < triangles->end(); i += 2) {
			triangle_ = triangle_indices[*i][*(i + 1)].data();
			getTVCollisionHessainForPairTest(Hessian, grad, vertex_position[*i][triangle_[0]].data(), vertex_position[*i][triangle_[1]].data(),
				vertex_position[*i][triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i + 1)),
				collision.triangle_vertex_pair_num_record[*i][*(i + 1)], *i,
				vertex_obj_No, tri_obj_No, triangle_, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
				d_hat_2,
				emp, 0);
		}
	}
	unsigned int* edge_;
	if (has_collider) {
		for (auto i = edges->begin(); i < edges->end(); i += 2) {
			edge_ = edge_vertices[*i] + ((*(i + 1)) << 1);
			getEECollisionHessainForPairTest(Hessian, grad, *i, edge_,
				collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i + 1)),
				collision.edge_edge_pair_number_record[*i][*(i + 1)],
				pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				vertex_obj_No, tri_obj_No,
				i - edges->begin(), edges, vertex_position[*i][*edge_].data(),
				vertex_position[*i][*(edge_ + 1)].data(), collision_stiffness, collision.edge_edge_collider_pair_by_edge[*i] + collision.close_ee_collider_pair_num * (*(i + 1)),
				collision.edge_edge_collider_pair_num_record[*i][*(i + 1)], rest_edge_length[*i][*(i + 1)]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = edges->begin(); i < edges->end(); i += 2) {
			edge_ = edge_vertices[*i] + ((*(i + 1)) << 1);
			getEECollisionHessainForPairTest(Hessian, grad, *i, edge_,
				collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i + 1)),
				collision.edge_edge_pair_number_record[*i][*(i + 1)],
				pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				vertex_obj_No, tri_obj_No,
				i - edges->begin(), edges, vertex_position[*i][*edge_].data(),
				vertex_position[*i][*(edge_ + 1)].data(), collision_stiffness, emp, 0, rest_edge_length[*i][*(i + 1)]);
		}
	}

	if (floor->exist) {
		for (int i = 0; i < unfixed_vertex_num; i += 2) {
			getFloorHessianForTet(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]].data(),
				floor->value,
				floor->dimension, collision_stiffness, floor->normal_direction, collision.d_hat_2, i >> 1, unfixed_vertex_num >> 1);
		}
	}
}


void XPBD_IPC::getCollisionBlockCollisionHessian(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangles,
	std::vector<unsigned int>* edges,
	double collision_stiffness, int* pair_actual_unfixed_vertex_indices, //pair_actual_unfixed_vertex_indices first obj, second primitive index
	int unfixed_vertex_num, double d_hat_2, int** vertex_index_on_surface, unsigned int vertex_obj_No, unsigned int tri_obj_No)//unfixed_vertex_num size double
{
	if (has_collider) {
		for (int i = 0; i < unfixed_vertex_num; i += 2) {
			getVTCollisionHessainForPair(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]].data(),
				collision_stiffness, collision.vertex_triangle_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_triangle_pair_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]],
				i >> 1, vertex_obj_No, tri_obj_No, pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				collision.vertex_obj_triangle_collider_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_collider_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_obj_triangle_collider_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (int i = 0; i < unfixed_vertex_num; i += 2) {
			getVTCollisionHessainForPair(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]].data(),
				collision_stiffness, collision.vertex_triangle_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]],
				collision.vertex_triangle_pair_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_index_on_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]]],
				i >> 1, vertex_obj_No, tri_obj_No, pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				emp, 0);
		}
	}
	int* triangle_;
	if (has_collider) {
		for (auto i = triangles->begin(); i < triangles->end(); i += 2) {
			triangle_ = triangle_indices[*i][*(i + 1)].data();
			getTVCollisionHessainForPair(Hessian, grad, vertex_position[*i][triangle_[0]].data(), vertex_position[*i][triangle_[1]].data(),
				vertex_position[*i][triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i + 1)),
				collision.triangle_vertex_pair_num_record[*i][*(i + 1)], *i,
				vertex_obj_No, tri_obj_No, triangle_, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
				d_hat_2,
				collision.triangle_vertex_collider_pair_by_triangle[*i] + collision.close_tv_collider_pair_num * (*(i + 1)),
				collision.triangle_vertex_collider_pair_num_record[*i][*(i + 1)]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = triangles->begin(); i < triangles->end(); i += 2) {
			triangle_ = triangle_indices[*i][*(i + 1)].data();
			getTVCollisionHessainForPair(Hessian, grad, vertex_position[*i][triangle_[0]].data(), vertex_position[*i][triangle_[1]].data(),
				vertex_position[*i][triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i + 1)),
				collision.triangle_vertex_pair_num_record[*i][*(i + 1)], *i,
				vertex_obj_No, tri_obj_No, triangle_, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
				d_hat_2,
				emp,0);

		}
	}
	unsigned int* edge_;
	if (has_collider) {
		for (auto i = edges->begin(); i < edges->end(); i += 2) {
			edge_ = edge_vertices[*i] + ((*(i + 1)) << 1);
			getEECollisionHessainForPair(Hessian, grad, *i, edge_,
				collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i + 1)),
				collision.edge_edge_pair_number_record[*i][*(i + 1)],
				pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				vertex_obj_No, tri_obj_No,
				i - edges->begin(), edges, vertex_position[*i][*edge_].data(),
				vertex_position[*i][*(edge_ + 1)].data(), collision_stiffness, collision.edge_edge_collider_pair_by_edge[*i] + collision.close_ee_collider_pair_num * (*(i + 1)),
				collision.edge_edge_collider_pair_num_record[*i][*(i + 1)], rest_edge_length[*i][*(i+1)]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = edges->begin(); i < edges->end(); i += 2) {
			edge_ = edge_vertices[*i] + ((*(i + 1)) << 1);
			getEECollisionHessainForPair(Hessian, grad, *i, edge_,
				collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i + 1)),
				collision.edge_edge_pair_number_record[*i][*(i + 1)],
				pair_actual_unfixed_vertex_indices, unfixed_vertex_num, d_hat_2,
				vertex_obj_No, tri_obj_No,
				i - edges->begin(), edges, vertex_position[*i][*edge_].data(),
				vertex_position[*i][*(edge_ + 1)].data(), collision_stiffness,emp,	0, rest_edge_length[*i][*(i + 1)]);
		}
	}
	if (floor->exist) {
		for (int i = 0; i < unfixed_vertex_num; i+=2) {
			getFloorHessianForTet(Hessian, grad, vertex_position[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i+1]].data(), 
				floor->value,
				floor->dimension, collision_stiffness, floor->normal_direction, collision.d_hat_2, i>>1, unfixed_vertex_num>>1);
		}
	}
}



void XPBD_IPC::getCollisionHessianFromRecord(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet,
	double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	int unfixed_tet_vertex_num, 
	int* vertex_index_on_surface, unsigned int* vt_prefix_sum, unsigned int* ee_prefix_sum, unsigned int* tv_collider_prefix_sum, unsigned int* ee_collider_prefix_sum, 
	std::array<double, 3>* vertex_position)
{
	//int surface_index;
	//if (has_collider) {
	//	for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
	//		if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
	//			continue;
	//		}
	//		surface_index = vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]];
	//		getVTCollisionHessainForTetFromRecord(Hessian, grad, 
	//			collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * surface_index,
	//			collision.vertex_triangle_pair_num_record[obj_No][surface_index],
	//			i, obj_No, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num,
	//			collision.vertex_obj_triangle_collider_num_record[obj_No][surface_index],
	//			collision.vt_hessian_record.data() + 144 * vt_prefix_sum[surface_index],
	//			collision.vt_grad_record.data() + 12 * vt_prefix_sum[surface_index],
	//			collision.vt_hessian_record_index.data() + 5 * vt_prefix_sum[surface_index],
	//			collision.vt_colldier_hessian_record.data() + (collision.vertex_num_on_surface_prefix_sum[obj_No] + surface_index) * 9,
	//			collision.vt_colldier_grad_record.data() + (collision.vertex_num_on_surface_prefix_sum[obj_No] + surface_index) * 3, 
	//			collision.vertex_obj_triangle_collider_pair_by_vertex[obj_No] + collision.close_vt_collider_pair_num * surface_index);
	//	}
	//}
	//else {
	//	unsigned int emp[1] = { 0 }; 
	//	double temp_double[1] = { 0.0 };
	//	for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
	//		if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
	//			continue;
	//		}
	//		surface_index = vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]];
	//		getVTCollisionHessainForTetFromRecord(Hessian, grad, 
	//			collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * surface_index,
	//			collision.vertex_triangle_pair_num_record[obj_No][surface_index],
	//			i, obj_No, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num,
	//			0,
	//			collision.vt_hessian_record.data() + 144 * vt_prefix_sum[surface_index],
	//			collision.vt_grad_record.data() + 12 * vt_prefix_sum[surface_index],
	//			collision.vt_hessian_record_index.data() + 5 * vt_prefix_sum[surface_index],
	//			temp_double, temp_double,emp);
	//	}
	//}

	//int* triangle_;
	//if (has_collider) {
	//	for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
	//		triangle_ = triangle_indices[obj_No][*i].data();
	//		getTVCollisionHessainForTetFromRecord(Hessian, grad,
	//			collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * (*i),
	//			collision.triangle_vertex_pair_num_record[obj_No][*i], obj_No, triangle_, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num,
	//			collision.triangle_vertex_collider_pair_by_triangle[obj_No] + collision.close_tv_collider_pair_num * (*i),
	//			collision.triangle_vertex_collider_pair_num_record[obj_No][*i],
	//			collision.tv_colldier_hessian_record_index.data() + 4 * tv_collider_prefix_sum[*i],
	//			collision.tv_colldier_hessian_record.data() + 81 * tv_collider_prefix_sum[*i],
	//			collision.tv_colldier_grad_record.data() + 9 * tv_collider_prefix_sum[*i]);
	//	}
	//}
	//else {
	//	unsigned int emp[1] = { 0 };
	//	double temp_double[1] = { 0.0 }; int temp_int[1] = { 0 };
	//	for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
	//		triangle_ = triangle_indices[obj_No][*i].data();
	//		getTVCollisionHessainForTetFromRecord(Hessian, grad,
	//			collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * (*i),
	//			collision.triangle_vertex_pair_num_record[obj_No][*i], obj_No, triangle_, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num,
	//			emp,0,
	//			temp_int,temp_double, temp_double);
	//	}
	//}


	//unsigned int* edge_;
	//if (has_collider) {
	//	for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
	//		edge_ = edge_vertices[obj_No] + ((*i) << 1);
	//		getEECollisionHessainForTetFromRecord(Hessian, grad, 
	//			obj_No, *i, edge_,
	//			collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * (*i),
	//			collision.edge_edge_pair_number_record[obj_No][*i],
	//			tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, 
	//			i - edge_of_a_tet->begin(), edge_of_a_tet, collision.edge_edge_collider_pair_by_edge[obj_No] + collision.close_ee_collider_pair_num * (*i),
	//			collision.edge_edge_collider_pair_num_record[obj_No][*i],
	//			collision.ee_hessian_record.data()+144*ee_prefix_sum[*i], collision.ee_grad_record.data()+12*ee_prefix_sum[*i],
	//			collision.ee_hessian_record_index.data()+5*ee_prefix_sum[*i],
	//			collision.ee_collider_hessian_record_index.data()+3*ee_collider_prefix_sum[*i],
	//			collision.ee_collider_hessian_record.data()+36*ee_collider_prefix_sum[*i],
	//			collision.ee_collider_grad_record.data()+6*ee_collider_prefix_sum[*i]);
	//	}
	//}
	//else {
	//	unsigned int emp[1] = { 0 };
	//	double temp_double[1] = { 0.0 }; int temp_int[1] = { 0 };
	//	for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
	//		edge_ = edge_vertices[obj_No] + ((*i) << 1);
	//		getEECollisionHessainForTetFromRecord(Hessian, grad, 
	//			obj_No, *i, edge_,
	//			collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * (*i),
	//			collision.edge_edge_pair_number_record[obj_No][*i],
	//			tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num,
	//			i - edge_of_a_tet->begin(), edge_of_a_tet, emp, 0,
	//			collision.ee_hessian_record.data() + 144 * ee_prefix_sum[*i], collision.ee_grad_record.data() + 12 * ee_prefix_sum[*i],
	//			collision.ee_hessian_record_index.data() + 5 * ee_prefix_sum[*i],
	//			temp_int, temp_double, temp_double);
	//	}
	//}

	//if (floor->exist) {
	//	int locate;
	//	for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
	//		if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
	//			continue;
	//		}
	//		locate = vertex_num_on_surface_prefix_sum[obj_No] + vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]];
	//		grad[3 * i + floor->dimension] += collision.floor_grad_record[locate];
	//		Hessian.data()[(3 * i + floor->dimension) * (3 * unfixed_tet_vertex_num + 1)] += collision.floor_hessian_record[locate];
	//	}
	//}
}




void XPBD_IPC::getCollisionHessian(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet,
	double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	int unfixed_tet_vertex_num, double d_hat_2,
	int* vertex_index_on_surface, std::array<double, 3>* vertex_position)
{
	if (has_collider) {
		for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
			if (tet_actual_unfixed_vertex_indices[i] == -1) {
				continue;
			}
			getVTCollisionHessainForTet(Hessian, grad, vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), collision_stiffness,
				collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * tet_actual_unfixed_vertex_indices[i],
				collision.vertex_triangle_pair_num_record[obj_No][tet_actual_unfixed_vertex_indices[i]],
				i, obj_No, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				collision.vertex_obj_triangle_collider_pair_by_vertex[obj_No] + collision.close_vt_collider_pair_num * tet_actual_unfixed_vertex_indices[i],
				collision.vertex_obj_triangle_collider_num_record[obj_No][tet_actual_unfixed_vertex_indices[i]]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
			if (tet_actual_unfixed_vertex_indices[i] == -1) {
				continue;
			}
			getVTCollisionHessainForTet(Hessian, grad, vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), collision_stiffness,
				collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * tet_actual_unfixed_vertex_indices[i],
				collision.vertex_triangle_pair_num_record[obj_No][tet_actual_unfixed_vertex_indices[i]],
				i, obj_No, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				emp, 0);
		}
	}
	
	int* triangle_;
	if (has_collider) {
		for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
			triangle_ = triangle_indices[obj_No][*i].data();
			getTVCollisionHessainForTet(Hessian, grad, vertex_position[triangle_[0]].data(), vertex_position[triangle_[1]].data(),
				vertex_position[triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * (*i),
				collision.triangle_vertex_pair_num_record[obj_No][*i], obj_No,
				triangle_, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				collision.triangle_vertex_collider_pair_by_triangle[obj_No] + collision.close_tv_collider_pair_num * (*i),
				collision.triangle_vertex_collider_pair_num_record[obj_No][*i]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
			triangle_ = triangle_indices[obj_No][*i].data();
			getTVCollisionHessainForTet(Hessian, grad, vertex_position[triangle_[0]].data(), vertex_position[triangle_[1]].data(),
				vertex_position[triangle_[2]].data(), collision_stiffness,
				collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * (*i),
				collision.triangle_vertex_pair_num_record[obj_No][*i], obj_No,
				triangle_, tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				emp, 0);
		}
	}




	unsigned int* edge_;
	if (has_collider) {
		for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
			edge_ = edge_vertices[obj_No] + ((*i) << 1);
			getEECollisionHessainForTet(Hessian, grad, obj_No, edge_,
				collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * (*i),
				collision.edge_edge_pair_number_record[obj_No][*i],
				tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				i - edge_of_a_tet->begin(), edge_of_a_tet, vertex_position[*edge_].data(),
				vertex_position[*(edge_ + 1)].data(), collision_stiffness, collision.edge_edge_collider_pair_by_edge[obj_No] + collision.close_ee_collider_pair_num * (*i), 
				collision.edge_edge_collider_pair_num_record[obj_No][*i],rest_edge_length[obj_No][*i]);
		}
	}
	else {
		unsigned int emp[1] = { 0 };
		for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
			edge_ = edge_vertices[obj_No] + ((*i) << 1);
			getEECollisionHessainForTet(Hessian, grad, obj_No, edge_,
				collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * (*i),
				collision.edge_edge_pair_number_record[obj_No][*i],
				tet_actual_unfixed_vertex_indices, unfixed_tet_vertex_num, d_hat_2,
				i - edge_of_a_tet->begin(), edge_of_a_tet, vertex_position[*edge_].data(),
				vertex_position[*(edge_ + 1)].data(), collision_stiffness, emp,
				0,rest_edge_length[obj_No][*i]);
		}
	}

	//if (Hessian.norm() != 0.0) {
	//	std::cout << "compare correct " << std::endl;
	//	std::cout << grad.transpose() << std::endl;
	//	std::cout << Hessian << std::endl;
	//}


	if (floor->exist) {
		for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
			if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
				continue;
			}
			getFloorHessianForTet(Hessian, grad, vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), floor->value,
				floor->dimension, collision_stiffness, floor->normal_direction, collision.d_hat_2,i, unfixed_tet_vertex_num);
		}		
	}

}


void XPBD_IPC::getFloorHessian(MatrixXd& Hessian, VectorXd& grad, int* tet_actual_unfixed_vertex_indices,
	int unfixed_tet_vertex_num, double d_hat_2,
	int* vertex_index_on_surface, std::array<double, 3>* vertex_position, double collision_stiffness)
{
	if (floor->exist) {
		for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
			if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
				continue;
			}
			getFloorHessianForTet(Hessian, grad, vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), floor->value,
				floor->dimension, collision_stiffness, floor->normal_direction, collision.d_hat_2, i, unfixed_tet_vertex_num);
		}
	}
}



void XPBD_IPC::getCollisionHessian(Matrix3d& Hessian, Vector3d& grad, std::array<double, 3>* vertex_position, 
	double* last_step_vertex_position,
	double collision_stiffness, unsigned int obj_No,
	unsigned int vertex_index,	unsigned int vertex_index_on_surface)
{
	Matrix3d a = Hessian;
	getVTCollisionHessain(Hessian, grad, vertex_position[vertex_index].data(), collision_stiffness,
		collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * vertex_index_on_surface, collision.vertex_triangle_pair_num_record[obj_No][vertex_index_on_surface],
		collision.VT_volume[obj_No].data() + collision.VT_start_index[obj_No][vertex_index_on_surface],obj_No,vertex_index);
	int* triangle_;
	std::vector<unsigned int>* triangle = &mesh_struct[obj_No]->vertices[vertex_index].face;
	for (unsigned int i = 0; i < triangle->size(); ++i) {
		triangle_ = triangle_indices[obj_No][triangle->data()[i]].data();
		getTVCollisionHessain(Hessian, grad, vertex_position[triangle_[0]].data(), vertex_position[triangle_[1]].data(),
			vertex_position[triangle_[2]].data(),
			findVertexNo(vertex_index, triangle_, 3),
			collision_stiffness, collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * triangle->data()[i],
			collision.triangle_vertex_pair_num_record[obj_No][triangle->data()[i]],
			collision.TV_volume[obj_No].data() + collision.TV_start_index[obj_No][triangle->data()[i]]);
	}
	std::vector<unsigned int>* edge = &mesh_struct[obj_No]->vertices[vertex_index].edge;
	unsigned int* edge_;
	for (unsigned int i = 0; i < edge->size(); ++i) {
		edge_ = edge_vertices[obj_No] + (edge->data()[i] << 1);
		getEECollisionHessian(Hessian, grad, vertex_position[edge_[0]].data(), vertex_position[edge_[1]].data(),
			collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * edge->data()[i],
			collision.edge_edge_pair_number_record[obj_No][edge->data()[i]], 
			collision.EE_volume[obj_No].data() + collision.EE_start_index[obj_No][edge->data()[i]], 
			collision_stiffness, obj_No, edge->data()[i],
			findVertexNo(vertex_index, edge_, 2));
	}
	if (!collider->empty()) {
		getVT_ColiderCollisionHessain(Hessian, grad, vertex_position[vertex_index].data(), collision_stiffness,
			collision.vertex_obj_triangle_collider_pair_by_vertex[obj_No] + collision.close_vt_collider_pair_num * vertex_index_on_surface,
			collision.vertex_obj_triangle_collider_num_record[obj_No][vertex_index_on_surface],
			collision.VT_collider_volume[obj_No].data() + collision.VT_collider_start_index[obj_No][vertex_index_on_surface]);
	}

	//floor
	if (floor->exist) {
		double Hessian_, grad_;
		if (getFloorHessian(Hessian_, grad_, vertex_position[vertex_index].data(), floor->value,
			last_step_vertex_position, floor->dimension, collision_stiffness, floor->normal_direction, collision.d_hat)) {
			grad.data()[floor->dimension] -= grad_;
			Hessian.data()[floor->dimension << 2] += Hessian_;
		}	
	}

}



void XPBD_IPC::solveInertialCollision(std::array<double, 3>* vertex_position,
	double* record_vertex_position_,
	double* last_step_vertex_position,
	double dt, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn, double collision_stiffness, unsigned int obj_No,
	bool vertex_on_surface, unsigned int vertex_index_on_surface)
{
	Matrix3d Hessian;
	Vector3d grad;
	Hessian.setZero();
	grad.setZero();
	//if (vertex_on_surface) {
	//	getCollisionHessian(Hessian, grad, vertex_position, last_step_vertex_position, collision_stiffness, obj_No, vertex_index, vertex_index_on_surface);
	//}
	double mass_dt_2 = mass[vertex_index] / (dt * dt);

	Hessian.data()[0] += mass_dt_2;
	Hessian.data()[4] += mass_dt_2;
	Hessian.data()[8] += mass_dt_2;
	grad.data()[0] -= mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	grad.data()[1] -= mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	grad.data()[2] -= mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <Matrix3d> linear(Hessian);
	Vector3d result = linear.solve(grad);

	SUM_(vertex_position[vertex_index], result);
	if (vertex_on_surface) {
		collision.collisionFreeOneVertex(obj_No, vertex_index, vertex_index_on_surface,
			record_vertex_position_, vertex_position[vertex_index].data(),
			record_vertex_position[obj_No].data(), vertex_position, this->vertex_position.data());
	}
	memcpy(record_vertex_position_, vertex_position[vertex_index].data(), 24);
}




void XPBD_IPC::solveNewtonCDTetWithCollision(std::array<double, 3>* vertex_position, 
	double* record_vertex_position_,
	double* last_step_vertex_position,
	double ARAP_stiffness, double dt,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* tet_vertex_indices, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn, double collision_stiffness, unsigned int obj_No,
	bool vertex_on_surface, unsigned int vertex_index_on_surface)
{
	Matrix3d Hessian;
	Vector3d grad;
	Hessian.setZero();
	grad.setZero();

	getARAPHessian(Hessian, grad, vertex_position, ARAP_stiffness, A, tet_indices, tet_vertex_indices, volume, vertex_index, obj_No);

	if (perform_collision) {
		if (vertex_on_surface) {
			getCollisionHessian(Hessian, grad, vertex_position, last_step_vertex_position, collision_stiffness, obj_No, vertex_index, vertex_index_on_surface);
		}
	}


	double mass_dt_2 = mass[vertex_index] / (dt * dt);

	Hessian.data()[0] += mass_dt_2;
	Hessian.data()[4] += mass_dt_2;
	Hessian.data()[8] += mass_dt_2;
	grad.data()[0] -= mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	grad.data()[1] -= mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	grad.data()[2] -= mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <Matrix3d> linear(Hessian);
	Vector3d result = linear.solve(grad);

	SUM_(vertex_position[vertex_index], result);
	
//	if (obj_No == 0) {
//		if (inner_iteration_number == 1) {
//			if (abs(result[1]) < 1e-4) {
//				std::cout << "index " << vertex_index << std::endl;
//			}
//		}
//		
///*		if (vertex_index==11 || vertex_index == 8) {
//			std::cout <<vertex_index<<" "<< sn[vertex_index][1] << " " << vertex_position[vertex_index][1] << std::endl;
//			std::cout << result[1] << " " << vertex_index<<" mdt2 "<<mass_dt_2 << std::endl;
//		}	*/	
//	}

	if (nearly_not_move) {
		for (unsigned int i = 0; i < 3; ++i) {
			if (abs(result[i]) > max_move_standard_inner_itr) {
				nearly_not_move = false;
				break;
			}
		}
	}	

	if (perform_collision) {
		if (vertex_on_surface) {
			collision.collisionFreeOneVertex(obj_No, vertex_index, vertex_index_on_surface,
				record_vertex_position_, vertex_position[vertex_index].data(),
				record_vertex_position[obj_No].data(), vertex_position, this->vertex_position.data());
		}
		memcpy(record_vertex_position_, vertex_position[vertex_index].data(), 24);
	}



}





void XPBD_IPC::getEECollisionHessian(Matrix3d& Hessian, Vector3d& grad, double* pos0, double* pos1, unsigned int* EE, unsigned int num,
	double* ori_volume, double stiffness, unsigned int obj_index, unsigned int edge_index, unsigned int vertex_no)
{
	Matrix3d Hessian_single;
	Vector3d grad_single;
	unsigned int* edge_vertex;
	double volume;
	for (unsigned int i = 0; i < num; i += 3) {
		volume = (ori_volume[i >> 1] > collision.volume_boundary ? ori_volume[i >> 1] : collision.volume_boundary);
		edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
		if (obj_index < EE[i] || (obj_index == EE[i] && edge_index < EE[i + 1])) {
			if (second_order_constraint.getCollisionPairHessian(pos0, pos1, vertex_position[EE[i]][edge_vertex[0]].data(),
				vertex_position[EE[i]][edge_vertex[1]].data(),
				volume, Hessian_single, grad_single, vertex_no)) {
				Hessian_single *= stiffness;
				grad_single *= stiffness;
				grad -= grad_single;
				Hessian += Hessian_single;
			}
		}
		else {
		if (second_order_constraint.getCollisionPairHessian(vertex_position[EE[i]][edge_vertex[0]].data(),
			vertex_position[EE[i]][edge_vertex[1]].data(), pos0, pos1,
			volume, Hessian_single, grad_single, vertex_no + 2)) {
			Hessian_single *= stiffness;
			grad_single *= stiffness;
			grad -= grad_single;
			Hessian += Hessian_single;
		}
		}
	}
}

void XPBD_IPC::getTVCollisionHessain(Matrix3d& Hessian, Vector3d& grad,
	double* pos_0, double* pos_1, double* pos_2,
	unsigned int vertex_no, double stiffness, unsigned int* TV, unsigned int num, double* ori_volume)
{
	Matrix3d Hessian_single;
	Vector3d grad_single;
	double volume;
	for (unsigned int i = 0; i < num; i += 2) {
		volume = (ori_volume[i >> 1] > collision.volume_boundary ? ori_volume[i >> 1] : collision.volume_boundary);
		if (second_order_constraint.getCollisionPairHessian(pos_0, pos_1, pos_2,
			vertex_position[TV[i]][TV[i + 1]].data(),
			volume, Hessian_single, grad_single, vertex_no)) {
			Hessian_single *= stiffness;
			grad_single *= stiffness;
			grad -= grad_single;
			Hessian += Hessian_single;
		}
	}
}



void XPBD_IPC::getVT_ColiderCollisionHessain(Matrix3d& Hessian, Vector3d& grad, double* vertex_position_, double stiffness,
	unsigned int* VT, unsigned int num, double* ori_volume)
{
	Matrix3d Hessian_single;
	Vector3d grad_single;
	int* triangle_vertex;
	double volume;
	for (unsigned int i = 0; i < num; i += 2) {
		triangle_vertex = triangle_indices_collider[VT[i]][VT[i + 1]].data();
		volume = (ori_volume[i >> 1] > collision.volume_boundary ? ori_volume[i >> 1] : collision.volume_boundary);
		if (second_order_constraint.getCollisionPairHessian(vertex_position_collider[VT[i]][triangle_vertex[0]].data(),
			vertex_position_collider[VT[i]][triangle_vertex[1]].data(), vertex_position_collider[VT[i]][triangle_vertex[2]].data(), vertex_position_,
			volume, Hessian_single, grad_single, 3)) {
			Hessian_single *= stiffness;
			grad_single *= stiffness;
			grad -= grad_single;
			Hessian += Hessian_single;
		}
	}
}


void XPBD_IPC::getEECollisionHessainForPairTest(MatrixXd& Hessian, VectorXd& grad, unsigned int ee_obj_No, unsigned int* edge_vertex_index,
	unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2,
	unsigned int obj_No_0, unsigned int obj_No_1,
	int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, double* ea0, double* ea1, double stiffness,
	unsigned int* EE_collider, int num_collider, double edge_length_0)//edge_order_in_tet should x2
{
	unsigned int* edge_vertex;
	int vertex_order_in_tet[4];
	memset(vertex_order_in_tet, 0xff, 8);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex_index, ee_obj_No, vertex_order_in_tet);
	for (int i = 0; i < num; i += 3) {
		double barrier_=0.0;
		edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
		memset(vertex_order_in_tet + 2, 0xff, 8);
		if (EE[i] == obj_No_0 || EE[i] == obj_No_1) {
			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, EE[i + 1], EE[i])) {
				continue;
			}
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex, EE[i], vertex_order_in_tet + 2);
		}
		second_order_constraint.computeEEBarrierGradientHessianTest(ea0, ea1, vertex_position[EE[i]][*edge_vertex].data(),
			vertex_position[EE[i]][edge_vertex[1]].data(), Hessian, grad,
			vertex_order_in_tet, stiffness, d_hat_2, edge_length_0, rest_edge_length[EE[i]][EE[i + 1]],
			barrier_);

		if (barrier_ != 0.0) {
			std::cout << "EE " << barrier_ << std::endl;
		}

	}

	if (has_collider) {
		memset(vertex_order_in_tet + 2, 0xff, 8);
		for (int i = 0; i < num_collider; i += 2) {
			edge_vertex = collider_edge_vertices[EE_collider[i]] + (EE_collider[i + 1] << 1);
			second_order_constraint.computeEEBarrierGradientHessian(ea0, ea1, vertex_position_collider[EE_collider[i]][*edge_vertex].data(),
				vertex_position_collider[EE_collider[i]][edge_vertex[1]].data(), Hessian, grad,
				vertex_order_in_tet, stiffness, d_hat_2, edge_length_0,
				rest_edge_length_collider[EE_collider[i]][EE_collider[i + 1]]);
		}
	}
}





void XPBD_IPC::getEECollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, unsigned int ee_obj_No, unsigned int* edge_vertex_index,
	unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2,
	unsigned int obj_No_0, unsigned int obj_No_1,
	int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, double* ea0, double* ea1, double stiffness,
	unsigned int* EE_collider, int num_collider, double edge_length_0)//edge_order_in_tet should x2
{
	unsigned int* edge_vertex;
	int vertex_order_in_tet[4];
	memset(vertex_order_in_tet, 0xff, 8);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex_index,ee_obj_No, vertex_order_in_tet);
	for (int i = 0; i < num; i += 3) {
		edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
		memset(vertex_order_in_tet + 2, 0xff, 8);
		if (EE[i] == obj_No_0 || EE[i]==obj_No_1) {
			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, EE[i + 1],EE[i])) {
				continue;
			}
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex, EE[i], vertex_order_in_tet + 2);
		}
		second_order_constraint.computeEEBarrierGradientHessian(ea0, ea1, vertex_position[EE[i]][*edge_vertex].data(),
			vertex_position[EE[i]][edge_vertex[1]].data(), Hessian, grad,
			vertex_order_in_tet, stiffness, d_hat_2, edge_length_0,rest_edge_length[EE[i]][EE[i+1]]);

	}

	if (has_collider) {
		memset(vertex_order_in_tet + 2, 0xff, 8);
		for (int i = 0; i < num_collider; i += 2) {
			edge_vertex = collider_edge_vertices[EE_collider[i]] + (EE_collider[i + 1] << 1);
			second_order_constraint.computeEEBarrierGradientHessian(ea0, ea1, vertex_position_collider[EE_collider[i]][*edge_vertex].data(),
				vertex_position_collider[EE_collider[i]][edge_vertex[1]].data(), Hessian, grad,
				vertex_order_in_tet, stiffness, d_hat_2, edge_length_0,
				rest_edge_length_collider[EE_collider[i]][EE_collider[i+1]]);
		}
	}
}

void XPBD_IPC::getARAPCollisionHessianForPair(MatrixXd& Hessian, VectorXd& grad, double stiffness, int tet_obj, int tet_index, int* tet_vertex, int* tet_unfixed_vertex_indices,
	int unfixed_tet_vertex_num, std::array<double,3>* vertex_position, Matrix<double, 3, 4>* A, double* volume)
{
	int triangle_vertex_order_in_tet[4];
	memset(triangle_vertex_order_in_tet, 0xff, 16);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, tet_vertex, tet_obj, triangle_vertex_order_in_tet, 4);
	second_order_constraint.solveTetCertainVertices(vertex_position, stiffness, A[tet_index], triangle_vertex_order_in_tet, tet_vertex, Hessian,
		volume[tet_index], grad);
}


void XPBD_IPC::getTVCollisionHessainForPairTest(MatrixXd& Hessian, VectorXd& grad, double* t0, double* t1, double* t2, double stiffness,
	unsigned int* TV, int num, unsigned int tri_obj, unsigned int obj_No_0, unsigned int obj_No_1, int* triangle_indices,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* TV_collider, int collider_num)
{
	int triangle_vertex_order_in_tet[4];
	memset(triangle_vertex_order_in_tet, 0xff, 16);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_indices, tri_obj, triangle_vertex_order_in_tet + 1, 3);

	for (int i = 0; i < num; i += 2) {
		if (TV[i] == obj_No_0 || TV[i] == obj_No_1) {
			if (vertexInPair(unfixed_tet_vertex_num, TV[i + 1], tet_unfixed_vertex_indices, TV[i])) {
				continue;
			}
		}

		double barrier = 0.0;

		second_order_constraint.computeVTBarrierGradientHessianTest(Hessian, grad, vertex_position[TV[i]][TV[i + 1]].data(),
			t0, t1, t2,
			d_hat_2, triangle_vertex_order_in_tet, stiffness, barrier);
		if (barrier != 0.0) {
			std::cout << "TV " << barrier << std::endl;
		}
	}

	if (has_collider) {
		for (int i = 0; i < collider_num; i += 2) {
			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_collider[TV_collider[i]][TV_collider[i + 1]].data(),
				t0, t1, t2,
				d_hat_2, triangle_vertex_order_in_tet, stiffness);
		}
	}

}




void XPBD_IPC::getTVCollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, double* t0, double* t1, double* t2, double stiffness,
	unsigned int* TV, int num, unsigned int tri_obj, unsigned int obj_No_0, unsigned int obj_No_1, int* triangle_indices,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* TV_collider, int collider_num)
{
	int triangle_vertex_order_in_tet[4];
	memset(triangle_vertex_order_in_tet, 0xff, 16);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_indices, tri_obj, triangle_vertex_order_in_tet+1,3);

	for (int i = 0; i < num; i += 3) {
		if (TV[i] == obj_No_0 || TV[i] == obj_No_1) {
			if (vertexInPair(unfixed_tet_vertex_num, TV[i + 1], tet_unfixed_vertex_indices,TV[i])) {
				continue;
			}
		}
		second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position[TV[i]][TV[i + 1]].data(),
			t0, t1, t2,
			d_hat_2, triangle_vertex_order_in_tet, stiffness);
	}

	if (has_collider) {
		for (int i = 0; i < collider_num; i += 2) {
			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_collider[TV_collider[i]][TV_collider[i + 1]].data(),
				t0, t1, t2,
				d_hat_2, triangle_vertex_order_in_tet, stiffness);
		}
	}

}

void XPBD_IPC::getVTCollisionHessainForPairTest(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int vertex_obj_No, unsigned int tri_obj_No,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* VT_collider, int num_collider)
{
	int* triangle_vertex;
	int triangle_vertex_order_in_tet[4];
	triangle_vertex_order_in_tet[0] = vertex_order_in_matrix;
	
	for (int i = 0; i < num; i += 2) {
		double barrier = 0.0;
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		if (VT[i] == tri_obj_No || VT[i] == vertex_obj_No) {
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, VT[i], triangle_vertex_order_in_tet + 1, 3);
		}
		second_order_constraint.computeVTBarrierGradientHessianTest(Hessian, grad, vertex_position_, vertex_position[VT[i]][triangle_vertex[0]].data(),
			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(),
			d_hat_2, triangle_vertex_order_in_tet, stiffness, barrier);

		if (barrier != 0.0) {
			std::cout << barrier << std::endl;
		}

	}

	if (has_collider) {
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		for (int i = 0; i < num_collider; i += 2) {
			triangle_vertex = triangle_indices_collider[VT_collider[i]][VT_collider[i + 1]].data();
			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position_collider[VT_collider[i]][triangle_vertex[0]].data(),
				vertex_position_collider[VT_collider[i]][triangle_vertex[1]].data(), vertex_position_collider[VT_collider[i]][triangle_vertex[2]].data(),
				d_hat_2, triangle_vertex_order_in_tet, stiffness);
		}
	}

}



void XPBD_IPC::getVTCollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int vertex_obj_No, unsigned int tri_obj_No,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* VT_collider, int num_collider)
{
	int* triangle_vertex;
	int triangle_vertex_order_in_tet[4];
	triangle_vertex_order_in_tet[0] = vertex_order_in_matrix;

	for (int i = 0; i < num; i += 2) {

		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		if (VT[i] == tri_obj_No || VT[i] == vertex_obj_No) {
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, VT[i], triangle_vertex_order_in_tet+1,3);
		}
		second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position[VT[i]][triangle_vertex[0]].data(),
			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(),
			d_hat_2, triangle_vertex_order_in_tet, stiffness);
	}

	if (has_collider) {
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		for (int i = 0; i < num_collider; i += 2) {
			triangle_vertex = triangle_indices_collider[VT_collider[i]][VT_collider[i + 1]].data();
			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position_collider[VT_collider[i]][triangle_vertex[0]].data(),
				vertex_position_collider[VT_collider[i]][triangle_vertex[1]].data(), vertex_position_collider[VT_collider[i]][triangle_vertex[2]].data(),
				d_hat_2, triangle_vertex_order_in_tet, stiffness);
		}
	}
}


void XPBD_IPC::getTVCollisionHessainForTetFromRecord(MatrixXd& Hessian, VectorXd& grad,
	unsigned int* TV, int num, unsigned int obj_No, int* triangle_indices,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, unsigned int* TV_collider, int collider_num,
	int* tv_collider_hessian_record_index, double* tv_collider_hessian_record, double* tv_collider_grad_record)
{
	//int triangle_vertex_order_in_tet[4];
	//memset(triangle_vertex_order_in_tet, 0xff, 16);
	//checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_indices,
	//	triangle_vertex_order_in_tet + 1, 3);
	//int pair_index;
	//int* vt_record_index;

	//for (int i = 0; i < num; i += 3) {
	//	if (TV[i] == obj_No) {
	//		if (vertexInTet(unfixed_tet_vertex_num, TV[i + 1], tet_unfixed_vertex_indices)) {
	//			continue;
	//		}
	//	}
	//	if (TV[i] < cloth->size()) {
	//		pair_index = collision.vertex_triangle_pair_num_record_prefix_sum[TV[i]][TV[i + 1]] + TV[i + 2];
	//	}
	//	else {
	//		pair_index = collision.vertex_triangle_pair_num_record_prefix_sum[TV[i]][vertex_index_surface[TV[i]][TV[i + 1]]]+TV[i + 2];
	//	}
	//	vt_record_index = collision.vt_hessian_record_index.data() + 5 * pair_index;

	//	if (*vt_record_index != 0) {
	//		setTetHessianFromBarrierHessian(Hessian, grad.data(), collision.vt_hessian_record.data() + 144 * pair_index,
	//			collision.vt_grad_record.data() + 12 * pair_index, triangle_vertex_order_in_tet,
	//			vt_record_index + 1, *vt_record_index);

	//	}
	//}
	//if (has_collider) {
	//	for (int i = 0; i < collider_num; i += 2) {
	//		if ((*tv_collider_hessian_record_index) != 0) {

	//			setTetHessianFromBarrierHessian(Hessian, grad.data(), tv_collider_hessian_record, tv_collider_grad_record,
	//				triangle_vertex_order_in_tet,
	//				tv_collider_hessian_record_index + 1, *tv_collider_hessian_record_index);
	//		}
	//		tv_collider_hessian_record_index += 4;
	//		tv_collider_hessian_record += 81;
	//		tv_collider_grad_record += 9;
	//	}
	//}

}


void XPBD_IPC::getVTCollisionHessainForTetFromRecord(MatrixXd& Hessian, VectorXd& grad,
	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int obj_No, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num,
	int num_collider, double* vt_hessian_record, double* vt_grad_record, int* vt_hessian_record_index,
	double* vt_collider_hessian_record, double* vt_collider_grad_record, unsigned int* vt_collider)
{
	int* triangle_vertex;
	int triangle_vertex_order_in_tet[4];
	triangle_vertex_order_in_tet[0] = vertex_order_in_matrix;
	for (unsigned int i = 0; i < num; i += 2) {
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		if (VT[i] == obj_No) {
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, triangle_vertex_order_in_tet + 1, 3);
		}

		if ((*vt_hessian_record_index) != 0) {
			setTetHessianFromBarrierHessian(Hessian, grad.data(), vt_hessian_record, vt_grad_record, triangle_vertex_order_in_tet,
				vt_hessian_record_index + 1, *vt_hessian_record_index);		
		}

		vt_hessian_record_index += 5;
		vt_grad_record += 12;
		vt_hessian_record += 144;
	}

	if (has_collider) {
		int vertex_pos = 3 * vertex_order_in_matrix * (Hessian.cols() + 1);
		if (num_collider > 0) {
			for (unsigned int i = 0; i < 3; ++i) {
				grad.data()[3 * vertex_order_in_matrix + i] += vt_collider_grad_record[i];
				Hessian.data()[vertex_pos + i] += vt_collider_hessian_record[i];
				Hessian.data()[vertex_pos + Hessian.cols() + i] += vt_collider_hessian_record[i + 3];
				Hessian.data()[vertex_pos + Hessian.cols() + Hessian.cols() + i] += vt_collider_hessian_record[i + 6];
			}
		}
	}
}


void XPBD_IPC::setHessianFromTetHessian(MatrixXd& Hessian_system, double* grad_system, double* Hessian_, double* grad_, int* vertex_in_sys)
{
	double* sys_start;
	double* record_start;
	double value;
	for (int i = 0; i < 4; ++i) {
		if (vertex_in_sys[i] != -1) {
			sys_start = grad_system + 3 * vertex_in_sys[i];
			record_start = grad_ + 3 * i;
			*sys_start += *record_start;
			*(sys_start + 1) += *(record_start + 1);
			*(sys_start + 2) += *(record_start + 2);

			for (int j = 0; j < 4; ++j) {
				if (vertex_in_sys[j] != -1) {
					value = Hessian_[i * 4 + j];
					sys_start = &Hessian_system(3 * vertex_in_sys[j], 3 * vertex_in_sys[i]);
					*sys_start += value;
					*(sys_start + Hessian_system.cols() + 1) += value;
					*(sys_start + ((Hessian_system.cols() + 1)<<1)) += value;
				}
			}
			
		}

	}
}


void XPBD_IPC::setTetHessianFromBarrierHessian(MatrixXd& Hessian_system, double* grad_system, double* Hessian_, double* grad_,
	int* triangle_vertex_order_in_system, int* vertex_in_pair, int vertex_in_use)
{
	double* sys_start;
	double* record_start;

	int size_sys = Hessian_system.cols();
	int size_record = vertex_in_use * 3;

	for (int i = 0; i < vertex_in_use; ++i) {
		if (triangle_vertex_order_in_system[vertex_in_pair[i]] != -1) {
			for (int j = 0; j < vertex_in_use; ++j) {
				if (triangle_vertex_order_in_system[vertex_in_pair[j]] != -1) {
					sys_start = &Hessian_system(3 * triangle_vertex_order_in_system[vertex_in_pair[i]], 3 * triangle_vertex_order_in_system[vertex_in_pair[j]]);
					record_start = Hessian_ + 9 * j * vertex_in_use + 3 * i;

					*sys_start += *record_start;
					*(sys_start + 1) += *(record_start + 1);
					*(sys_start + 2) += *(record_start + 2);

					*(sys_start + size_sys) += *(record_start + size_record);
					*(sys_start + size_sys + 1) += *(record_start + size_record + 1);
					*(sys_start + size_sys + 2) += *(record_start + size_record + 2);

					*(sys_start + (size_sys<<1)) += *(record_start + (size_record<<1));
					*(sys_start + (size_sys << 1) + 1) += *(record_start + (size_record << 1) + 1);
					*(sys_start + (size_sys << 1) + 2) += *(record_start + (size_record << 1) + 2);
				}

			}
			sys_start = grad_system + 3 * triangle_vertex_order_in_system[vertex_in_pair[i]];
			record_start = grad_ + 3 * i;

			*sys_start += *record_start;
			*(sys_start + 1) += *(record_start + 1);
			*(sys_start + 2) += *(record_start + 2);
		}
	}
}


void XPBD_IPC::getVTCollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int obj_No,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* VT_collider, int num_collider)
{
	int* triangle_vertex;
	int triangle_vertex_order_in_tet[4];
	triangle_vertex_order_in_tet[0] = vertex_order_in_matrix;
	for (unsigned int i = 0; i < num; i += 2) {
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		memset(triangle_vertex_order_in_tet +1, 0xff, 12);
		if (VT[i] == obj_No) {
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, triangle_vertex_order_in_tet+1,3);
		}
		second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position[VT[i]][triangle_vertex[0]].data(),
			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(),
			d_hat_2, triangle_vertex_order_in_tet, stiffness);
	}
	if (has_collider) {
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		for (unsigned int i = 0; i < num_collider; i += 2) {
			triangle_vertex = triangle_indices_collider[VT_collider[i]][VT_collider[i + 1]].data();
			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position_collider[VT_collider[i]][triangle_vertex[0]].data(),
				vertex_position_collider[VT_collider[i]][triangle_vertex[1]].data(), vertex_position_collider[VT_collider[i]][triangle_vertex[2]].data(),
				d_hat_2, triangle_vertex_order_in_tet, stiffness);
		}
	}
}


//void XPBD_IPC::getVTCollisionHessainForPair(MatrixXd& Hessian, VectorXd& grad, double* vertex_position_, double stiffness,
//	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, 
//	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2)
//{
//	int* triangle_vertex;
//	int triangle_vertex_order_in_tet[4];
//	triangle_vertex_order_in_tet[0] = vertex_order_in_matrix;
//	for (unsigned int i = 0; i < num; i += 2) {
//		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
//		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
//		checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, VT[i], triangle_vertex_order_in_tet);		
//		second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_, vertex_position[VT[i]][triangle_vertex[0]].data(),
//			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(),
//			d_hat_2, triangle_vertex_order_in_tet, stiffness);
//	}
//}




void XPBD_IPC::getFloorHessianForTet(MatrixXd& Hessian, VectorXd& grad, double* vertex_position, double floor_value,
	unsigned int dimension, double collision_stiffness, bool direction, double d_hat_2, unsigned int vertex_order_in_matrix, unsigned int unfixed_vertex_num)
{
	double distance = (vertex_position[dimension] - floor_value) * (vertex_position[dimension] - floor_value);

	if (distance > d_hat_2) {
		return;
	}
	
	double h, g;
	barrierGradHessian(distance, d_hat_2, g, h);

	//if (vertex_order_in_matrix == 1) {
		//std::cout <<"vertex_order_in_matrix "<< vertex_order_in_matrix<<" "<< g << " " << h << std::endl;
	//}

	double grad_d, h_d;
	grad_d = 2 * (vertex_position[dimension] - floor_value);
	h_d = 2;
	grad[3 * vertex_order_in_matrix + dimension] += collision_stiffness * g * grad_d;
	Hessian.data()[(3 * vertex_order_in_matrix + dimension) * (3 * unfixed_vertex_num + 1)] += collision_stiffness * (h * grad_d * grad_d + g * h_d);

	//std::cout << d_hat_2<<" "<< collision_stiffness * (h * grad_d * grad_d + g * h_d)<<" "<<collision_stiffness<<" "<<floor_value << std::endl;

}




void XPBD_IPC::getTVCollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, double* t0, double* t1, double* t2, double stiffness,
	unsigned int* TV, int num, unsigned int obj_No, int* triangle_indices,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2, unsigned int* TV_collider, int collider_num)
{
	int triangle_vertex_order_in_tet[4];
	memset(triangle_vertex_order_in_tet, 0xff, 16);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_indices, triangle_vertex_order_in_tet+1,3);

	for (int i = 0; i < num; i += 3) {
		if (TV[i] == obj_No) {
			if (vertexInTet(unfixed_tet_vertex_num, TV[i + 1], tet_unfixed_vertex_indices)) {
				continue;
			}
		}
		second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position[TV[i]][TV[i + 1]].data(),
			t0, t1, t2,
			d_hat_2, triangle_vertex_order_in_tet, stiffness);
	}
	
	if (has_collider) {
		for (int i = 0; i < collider_num; i += 2) {

			//MatrixXd hessian_test = Hessian;
			//hessian_test.setZero();

			second_order_constraint.computeVTBarrierGradientHessian(Hessian, grad, vertex_position_collider[TV_collider[i]][TV_collider[i + 1]].data(),
				t0, t1, t2,
				d_hat_2, triangle_vertex_order_in_tet, stiffness);

			//Hessian += hessian_test;
			//std::cout << "hessian for test " << TV_collider[i + 1] << " " <<  std::endl;
			//std::cout << t0[0] << " " << t0[1] << " " << t0[2] << std::endl;
			//std::cout << t1[0] << " " << t1[1] << " " << t1[2] << std::endl;
			//std::cout << t2[0] << " " << t2[1] << " " << t2[2] << std::endl;
			//std::cout << "stiffness " << stiffness << std::endl;
			//std::cout << hessian_test << std::endl;


		}
	}

}



//void XPBD_IPC::checkIfRight(int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, unsigned int obj_No, int* triangle_vertex_order_in_tet)
//{
//
//}



void XPBD_IPC ::getEECollisionHessainForTetFromRecord(MatrixXd& Hessian, VectorXd& grad, 
	unsigned int obj_No, unsigned int edge_index, 
	unsigned int* edge_vertex_index,
	unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, 
	int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, 
	unsigned int* EE_collider, int num_collider, 
	double* ee_hessian_record_, double* ee_grad_record_, int* ee_hessian_record_index_,
	int* ee_collider_hessian_record_index, double* ee_collider_hessian_record, double* ee_collider_grad_record)
{
	//unsigned int* edge_vertex;
	//int vertex_order_in_tet[4];
	//memset(vertex_order_in_tet, 0xff, 8);
	//checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex_index, vertex_order_in_tet);
	//int back_up_order[2];
	//memcpy(back_up_order, vertex_order_in_tet, 8);
	//double* ee_hessian_record; double* ee_grad_record; int* ee_hessian_record_index;
	//int pair_index;
	//for (int i = 0; i < num; i += 3) {
	//	edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
	//	if (EE[i] > obj_No || (EE[i] == obj_No && EE[i + 1] > edge_index)) {


	//		memset(vertex_order_in_tet + 2, 0xff, 8);
	//		if (EE[i] == obj_No) {
	//			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, EE[i + 1])) {
	//				continue;
	//			}
	//			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex, vertex_order_in_tet + 2);
	//		}

	//		ee_hessian_record_index = ee_hessian_record_index_ + i / 3 * 5;
	//		ee_grad_record = ee_grad_record_ + i / 3 * 12;
	//		ee_hessian_record = ee_hessian_record_ + i / 3 * 144;

	//	}
	//	else {

	//		memcpy(vertex_order_in_tet + 2, back_up_order, 8);
	//		memset(vertex_order_in_tet, 0xff, 8);

	//		if (EE[i] == obj_No) {
	//			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, EE[i + 1])) {
	//				continue;
	//			}
	//			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex, vertex_order_in_tet);
	//		}

	//		pair_index = collision.edge_edge_pair_num_record_prefix_sum[EE[i]][EE[i + 1]] + EE[i + 2] / 3;
	//		ee_hessian_record_index = collision.ee_hessian_record_index.data() + 5 * pair_index;
	//		ee_grad_record = collision.ee_grad_record.data() + 12 * pair_index;
	//		ee_hessian_record = collision.ee_hessian_record.data() + 144 * pair_index;
	//	}
	//	if ((*ee_hessian_record_index) != 0) {

	//		setTetHessianFromBarrierHessian(Hessian, grad.data(), ee_hessian_record, ee_grad_record, vertex_order_in_tet,
	//			ee_hessian_record_index + 1, *ee_hessian_record_index);
	//	}
	//}
	//if (has_collider) {
	//	memcpy(vertex_order_in_tet, back_up_order, 8);
	//	memset(vertex_order_in_tet + 2, 0xff, 8);
	//	for (int i = 0; i < num_collider; i += 2) {
	//		if ((*ee_collider_hessian_record_index) != 0) {
	//			setTetHessianFromBarrierHessian(Hessian, grad.data(), ee_collider_hessian_record, ee_collider_grad_record, vertex_order_in_tet,
	//				ee_collider_hessian_record_index + 1, *ee_collider_hessian_record_index);
	//		}
	//		ee_collider_hessian_record_index += 3;
	//		ee_collider_grad_record += 6;
	//		ee_collider_hessian_record += 36;
	//	}
	//}

}


void XPBD_IPC::getEECollisionHessainForTet(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_No, unsigned int* edge_vertex_index,
	unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, double d_hat_2,
	int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, double*ea0, double* ea1, double stiffness,
	unsigned int* EE_collider, int num_collider, double edge_length_0)
{
	unsigned int* edge_vertex;
	int vertex_order_in_tet[4];
	memset(vertex_order_in_tet, 0xff, 8);	

	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex_index, vertex_order_in_tet);

	for (int i = 0; i < num; i += 3) {
		edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
		memset(vertex_order_in_tet + 2, 0xff, 8);
		if (EE[i] == obj_No) {
			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, EE[i + 1])) {
				continue;
			}
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex, vertex_order_in_tet+2);
		}
		second_order_constraint.computeEEBarrierGradientHessian(ea0, ea1,  vertex_position[EE[i]][*edge_vertex].data(), 
			vertex_position[EE[i]][edge_vertex[1]].data(), Hessian, grad,
			vertex_order_in_tet, stiffness, d_hat_2, edge_length_0,rest_edge_length[EE[i]][EE[i+1]]);

		//std::cout << Hessian << std::endl;
		//std::cout << "++++" << std::endl;

	}

	if (has_collider) {
		
		memset(vertex_order_in_tet + 2, 0xff, 8);
		for (int i = 0; i < num_collider; i += 2) {
			edge_vertex = collider_edge_vertices[EE_collider[i]] + (EE_collider[i + 1] << 1);
			second_order_constraint.computeEEBarrierGradientHessian(ea0, ea1, vertex_position_collider[EE_collider[i]][*edge_vertex].data(),
				vertex_position_collider[EE_collider[i]][edge_vertex[1]].data(), Hessian, grad,
				vertex_order_in_tet, stiffness, d_hat_2, edge_length_0,rest_edge_length_collider[EE[i]][EE[i+1]]);
			//std::cout << "===========------" << std::endl;
			//std::cout << Hessian << std::endl;

		}
	}
}


bool XPBD_IPC::edgeInSameTetDuplicate(int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet,
	unsigned int compare_edge_index)
{
	for (auto i = edge_of_a_tet->begin() + edge_order_in_tet + 1; i < edge_of_a_tet->end(); ++i) {
		if (*i == compare_edge_index) {
			return true;
		}
	}
	return false;
}


bool XPBD_IPC::edgeInSameTetDuplicate(int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet,
	unsigned int compare_edge_index, unsigned int compare_edge_obj)
{
	for (auto i = edge_of_a_tet->begin() + edge_order_in_tet + 2; i < edge_of_a_tet->end(); i+=2) {
		if (*i == compare_edge_obj &&  *(i+1) == compare_edge_index) {
			return true;
		}
	}
	return false;
}

bool XPBD_IPC::vertexInTet(int unfixed_tet_vertex_num, int vertex_No, int* tet_unfixed_vertex_indices)
{
	for (int i = 0; i < unfixed_tet_vertex_num; ++i) {
		if (vertex_No == tet_unfixed_vertex_indices[i]) {
			return true;
		}
	}
	return false;
}

bool XPBD_IPC::vertexInPair(int unfixed_tet_vertex_num, int vertex_No, int* tet_unfixed_vertex_indices, int obj_No)
{
	for (int i = 0; i < unfixed_tet_vertex_num; i+=2) {
		if (vertex_No == tet_unfixed_vertex_indices[i+1] && obj_No == tet_unfixed_vertex_indices[i]) {
			return true;
		}
	}
	return false;
}


void XPBD_IPC::checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, int* element_indices, int obj_No,
	int* triangle_vertex_order_in_system, int size_num)
{
	for (int i = 0; i < size_num; ++i) {
		for (int j = 0; j < unfixed_tet_vertex_num; j+=2) {
			if (tet_unfixed_vertex_indices[j+1] == element_indices[i] && tet_unfixed_vertex_indices[j]==obj_No) {
				triangle_vertex_order_in_system[i] = (j>>1);
				break;
			}
		}
	}

}


void XPBD_IPC::checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, int* element_indices,
	int* triangle_vertex_order_in_system, int size_num)
{
	for (int i = 0; i < size_num; ++i) {
		for (int j = 0; j < unfixed_tet_vertex_num; ++j) {
			if (tet_unfixed_vertex_indices[j] == element_indices[i]) {
				triangle_vertex_order_in_system[i] = j;
				break;
			}
		}
	}

}

void XPBD_IPC::checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, unsigned int* element_indices,
	int* triangle_vertex_order_in_system)
{
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < unfixed_tet_vertex_num; ++j) {
			if (tet_unfixed_vertex_indices[j] == element_indices[i]) {
				triangle_vertex_order_in_system[i] = j;
				break;
			}
		}
	}

}

void XPBD_IPC::checkPairIndexInSys(int unfixed_tet_vertex_num, int* tet_unfixed_vertex_indices, unsigned int* element_indices, int obj_No,
	int* triangle_vertex_order_in_system)
{
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < unfixed_tet_vertex_num; j+=2) {
			if (tet_unfixed_vertex_indices[j+1] == element_indices[i] && tet_unfixed_vertex_indices[j]==obj_No) {
				triangle_vertex_order_in_system[i] = (j>>1);
				break;
			}
		}
	}

}

void XPBD_IPC::getVTCollisionHessain(Matrix3d& Hessian, Vector3d& grad, double* vertex_position_, double stiffness, 
	unsigned int* VT, unsigned int num, double* ori_volume, unsigned int obj_No, unsigned int vertex_index)
{
	Matrix3d Hessian_single;
	Vector3d grad_single;
	int* triangle_vertex;
	double volume;
	for (unsigned int i = 0; i < num; i += 2) {
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		volume =( ori_volume[i >> 1]>collision.volume_boundary?ori_volume[i >> 1] :collision.volume_boundary);
		if (second_order_constraint.getCollisionPairHessian( vertex_position[VT[i]][triangle_vertex[0]].data(),
			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(), vertex_position_,
			volume, Hessian_single, grad_single, 3)) {
			Hessian_single *= stiffness;
			grad_single *= stiffness;
			grad -= grad_single;
			Hessian += Hessian_single;
		}
	}


}


void XPBD_IPC::getARAPHessian(Matrix3d& Hessian, Vector3d& grad, std::array<double, 3>* vertex_position, double stiffness,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, 
	double* volume, unsigned int vertex_index, unsigned int obj_No)
{
	unsigned int tet_index, vertex_no;
	Matrix3d Hessian_single;
	Vector3d grad_single;
	double C;
	for (unsigned int i = 0; i < tet_indices.size(); ++i) {
		tet_index = tet_indices[i];
		vertex_no = findVertexNo(vertex_index, indices[tet_index].data(),4);

		if (second_order_constraint.getARAPGradHessianNewton(vertex_position[indices[tet_index][0]].data(), vertex_position[indices[tet_index][1]].data(),
			vertex_position[indices[tet_index][2]].data(), vertex_position[indices[tet_index][3]].data(),
			A[tet_index], Hessian_single, grad_single, C, vertex_no)) {
			Hessian += (0.5 * stiffness * volume[tet_index]) * Hessian_single;
			grad -= (0.5 * stiffness * volume[tet_index]) * grad_single;
		}
		//if (obj_No == 0 && vertex_index == 11) {
		//	std::cout << C << " " << grad_single.transpose() << std::endl;
		//	std::cout << Hessian_single << std::endl;
		//}

	}
}


void XPBD_IPC::solveEE_collisionBlock(unsigned int obj_No_0, unsigned int primitive_0_index, unsigned int obj_No_1, unsigned int primitive_1_index,
	double stiffness, double dt, double collision_stiffne, std::vector<unsigned int>* triangle_around_0, std::vector<unsigned int>* triangle_around_1,
	std::vector<unsigned int>* edge_around_0, std::vector<unsigned int>* edge_around_1,
	std::vector<unsigned int>* tet_around_0, std::vector<unsigned int>* tet_around_1,
	double d_hat_2, bool edge_0_collider, bool edge_1_collider)
{
	MatrixXd Hessian;
	VectorXd grad;
	int unfixed_pair_vertex_index[8];
	memset(unfixed_pair_vertex_index, 0xff, 32);
	int unfixed_num = 0;

	if (!edge_0_collider) {
		unsigned int* edge_vertex_0 = edge_vertices[obj_No_0] + (primitive_0_index << 1);
		if (!(*is_vertex_fixed[obj_No_0])[*edge_vertex_0]) {
			unfixed_pair_vertex_index[unfixed_num] = obj_No_0;
			unfixed_pair_vertex_index[unfixed_num + 1] = *edge_vertex_0;
			unfixed_num += 2;
		}
		if (!(*is_vertex_fixed[obj_No_0])[edge_vertex_0[1]]) {
			unfixed_pair_vertex_index[unfixed_num] = obj_No_0;
			unfixed_pair_vertex_index[unfixed_num + 1] = edge_vertex_0[1];
			unfixed_num += 2;
		}
	}

	if (!edge_1_collider) {
		unsigned int* edge_vertex_1 = edge_vertices[obj_No_1] + (primitive_1_index << 1);
		if (!(*is_vertex_fixed[obj_No_1])[*edge_vertex_1]) {
			unfixed_pair_vertex_index[unfixed_num] = obj_No_1;
			unfixed_pair_vertex_index[unfixed_num + 1] = *edge_vertex_1;
			unfixed_num += 2;
		}
		if (!(*is_vertex_fixed[obj_No_1])[edge_vertex_1[1]]) {
			unfixed_pair_vertex_index[unfixed_num] = obj_No_1;
			unfixed_pair_vertex_index[unfixed_num + 1] = edge_vertex_1[1];
			unfixed_num += 2;
		}
	}

	if (unfixed_num == 0) {
		return;
	}

	grad.resize(3 * (unfixed_num >> 1));
	grad.setZero();
	Hessian.resize(3 * (unfixed_num >> 1), 3 * (unfixed_num >> 1));
	Hessian.setZero();

	std::vector<unsigned int> around_triangle, around_edge, around_tet;

	getCollisionPairHessian(Hessian, grad, obj_No_0, obj_No_1, collision_stiffne, triangle_around_0, triangle_around_1,
		edge_around_0, edge_around_1,
		tet_around_0, tet_around_1,
		d_hat_2, unfixed_pair_vertex_index, unfixed_num, &around_triangle, &around_edge, &around_tet, stiffness,
		edge_0_collider, edge_1_collider);

	double mass_dt_2;
	int obj_index_, vertex_index_;

	int i;

	for (int j = 0; j < unfixed_num; j += 2) {
		i = (j >> 1);
		vertex_index_ = unfixed_pair_vertex_index[j + 1];
		obj_index_ = unfixed_pair_vertex_index[j];
		mass_dt_2 = mass[obj_index_][vertex_index_] / (dt * dt);
		Hessian(3 * i, 3 * i) += mass_dt_2;
		Hessian(3 * i + 1, 3 * i + 1) += mass_dt_2;
		Hessian(3 * i + 2, 3 * i + 2) += mass_dt_2;
		grad.data()[3 * i] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][0] - sn[obj_index_][vertex_index_][0]);
		grad.data()[3 * i + 1] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][1] - sn[obj_index_][vertex_index_][1]);
		grad.data()[3 * i + 2] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][2] - sn[obj_index_][vertex_index_][2]);
	}

	//HouseholderQR <MatrixXd> linear(Hessian);
	LLT <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(grad);

	double tt = result.norm();


	//if (edge_1_collider == false && primitive_0_index == 3533 && primitive_1_index == 6302) {
	//	//std::cout<<ed
	//	//double e0[3]; double e1[3];
	//	//std::vector<unsigned int> around_triangle_, around_edge_, around_tet_;
	//	//MatrixXd hessen_record = Hessian; hessen_record.setZero();
	//	//VectorXd grad_record = grad; grad_record.setZero();
	//	//getCollisionPairHessianTest(hessen_record, grad_record, obj_No_0, obj_No_1, collision_stiffne, triangle_around_0, triangle_around_1,
	//	//	edge_around_0, edge_around_1,
	//	//	tet_around_0, tet_around_1,
	//	//	d_hat_2, unfixed_pair_vertex_index, unfixed_num, &around_triangle_, &around_edge_, &around_tet_, stiffness,
	//	//	edge_0_collider, edge_1_collider);
	//}

	//if (tt > 1e-2) {
	//	std::cout <<"displacement ee"<< edge_1_collider<<" "<< primitive_0_index << " " << primitive_1_index << " " << tt << std::endl;
	//}
	

	for (int i = 0; i < unfixed_num; i += 2) {
		vertex_index_ = unfixed_pair_vertex_index[i + 1];
		obj_index_ = unfixed_pair_vertex_index[i];
		SUB_(vertex_position[obj_index_][vertex_index_], (result.data() + 3 * (i >> 1)));
	}

	double t = getCollisionTime(&around_triangle, &around_edge, unfixed_pair_vertex_index, unfixed_num, vertex_index_surface.data(),
		vertex_position.data(), record_vertex_position_address.data());

	if (t < 1.0) {
		t *= 0.9;
		double* p_c; double* p_i;
		for (int i = 0; i < unfixed_num; i += 2) {
			p_i = record_vertex_position_address[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data();
			p_c = vertex_position[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data();
			COLLISION_POS(p_c, t, p_i, p_c);
			memcpy(p_i, p_c, 24);
		}
	}
	else {
		for (int i = 0; i < unfixed_num; i += 2) {
			memcpy(record_vertex_position_address[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data(),
				vertex_position[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data(), 24);
		}
	}

}





void XPBD_IPC::solveVT_collisionBlock(unsigned int vertex_obj_no, unsigned int vertex_index, unsigned int triangle_obj_No, unsigned int triangle_index,
	double stiffness, double dt, double collision_stiffne, std::vector<unsigned int>* triangle_around_vertex, std::vector<unsigned int>* triangle_around_triangle,
	std::vector<unsigned int>* edge_around_vertex, std::vector<unsigned int>* edge_around_triangle, 
	std::vector<unsigned int>* tet_around_vertex, std::vector<unsigned int>*tet_around_triangle,
	double d_hat_2, bool vertex_collider, bool triangle_collider)
{
	MatrixXd Hessian;
	VectorXd grad;
	int unfixed_pair_vertex_index[8];
	memset(unfixed_pair_vertex_index, 0xff, 32);
	int unfixed_num = 0;
	int* tri_indices = triangle_indices[triangle_obj_No][triangle_index].data();
	if (!vertex_collider) {
		if (!(*is_vertex_fixed[vertex_obj_no])[vertex_index]) {
			unfixed_pair_vertex_index[unfixed_num] = vertex_obj_no;
			unfixed_pair_vertex_index[unfixed_num + 1] = vertex_index;
			unfixed_num += 2;
		}
	}
	if (!triangle_collider) {
		for (unsigned int i = 0; i < 3; ++i) {
			if (!(*is_vertex_fixed[triangle_obj_No])[tri_indices[i]]) {
				unfixed_pair_vertex_index[unfixed_num] = triangle_obj_No;
				unfixed_pair_vertex_index[unfixed_num + 1] = tri_indices[i];
				unfixed_num += 2;
			}
		}
	}
	if (unfixed_num == 0) {
		return;
	}

	grad.resize(3*(unfixed_num >> 1));
	grad.setZero();
	Hessian.resize(3 * (unfixed_num >> 1), 3 * (unfixed_num >> 1));
	Hessian.setZero();
	
	std::vector<unsigned int> around_triangle, around_edge, around_tet;




	getCollisionPairHessian(Hessian, grad, vertex_obj_no, triangle_obj_No, collision_stiffne, triangle_around_vertex, triangle_around_triangle,
		edge_around_vertex, edge_around_triangle,
		tet_around_vertex, tet_around_triangle,
		d_hat_2,unfixed_pair_vertex_index, unfixed_num,&around_triangle,&around_edge,&around_tet, stiffness,vertex_collider, triangle_collider);


	double mass_dt_2;
	int obj_index_,vertex_index_;

	int i;
	for (int j = 0; j < unfixed_num; j+=2) {
		i =( j >> 1);
		vertex_index_ = unfixed_pair_vertex_index[j+1];
		obj_index_= unfixed_pair_vertex_index[j];
		mass_dt_2 = mass[obj_index_][vertex_index_] / (dt * dt);
		Hessian(3 * i, 3 * i) += mass_dt_2;
		Hessian(3 * i + 1, 3 * i + 1) += mass_dt_2;
		Hessian(3 * i + 2, 3 * i + 2) += mass_dt_2;
		grad.data()[3 * i] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][0] - sn[obj_index_][vertex_index_][0]);
		grad.data()[3 * i + 1] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][1] - sn[obj_index_][vertex_index_][1]);
		grad.data()[3 * i + 2] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][2] - sn[obj_index_][vertex_index_][2]);
	}

	//HouseholderQR <MatrixXd> linear(Hessian);
	LLT <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(grad);
	//std::cout << Hessian << std::endl;
	//std::cout << "++" << std::endl;
	//std::cout << grad.transpose() << std::endl;

	double tt = result.norm();

	//if (tt > 1e-2) {
	//	std::cout << "displacement vt" << triangle_collider <<" "<< vertex_collider << " " << vertex_index << " " << triangle_index << " " << tt << std::endl;
	//}

	for (int i = 0; i < unfixed_num; i+=2) {
		vertex_index_ = unfixed_pair_vertex_index[i + 1];
		obj_index_ = unfixed_pair_vertex_index[i];
		SUB_(vertex_position[obj_index_][vertex_index_], (result.data() + 3 * (i>>1)));
	}

	double t = getCollisionTime(&around_triangle, &around_edge, unfixed_pair_vertex_index, unfixed_num, vertex_index_surface.data(),
		vertex_position.data(), record_vertex_position_address.data());

	if (t < 1.0) {
		t *= 0.9;
		double* p_c; double* p_i;
		for (int i = 0; i < unfixed_num; i += 2) {
			p_i = record_vertex_position_address[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data();
			p_c = vertex_position[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i + 1]].data();
			COLLISION_POS(p_c, t, p_i, p_c);
			memcpy(p_i, p_c, 24);
		}
	}
	else {
		for (int i = 0; i < unfixed_num; i += 2) {
			memcpy(record_vertex_position_address[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i+1]].data(),
				vertex_position[unfixed_pair_vertex_index[i]][unfixed_pair_vertex_index[i+1]].data(), 24);
		}
	}
}




void XPBD_IPC::getCollisionPairHessianTest(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_0, unsigned int obj_1,
	double collision_stiffness, std::vector<unsigned int>* triangle_around_0, std::vector<unsigned int>* triangle_around_1,
	std::vector<unsigned int>* edge_around_0, std::vector<unsigned int>* edge_around_1,
	std::vector<unsigned int>* tet_around_0, std::vector<unsigned int>* tet_around_1,
	double d_hat_2, int* unfixed_pair_vertex_index, int unfixed_num, std::vector<unsigned int>* around_triangle, std::vector<unsigned int>* around_edge,
	std::vector<unsigned int>* around_tet, double arap_stiffness, bool obj_0_collider, bool obj_1_collider)
{
	if (obj_0_collider) {
		for (auto i = triangle_around_1->begin(); i < triangle_around_1->end(); ++i) {
			around_triangle->emplace_back(obj_1);
			around_triangle->emplace_back(*i);
		}
		for (auto i = edge_around_1->begin(); i < edge_around_1->end(); ++i) {
			around_edge->emplace_back(obj_1);
			around_edge->emplace_back(*i);
		}
		for (auto i = tet_around_1->begin(); i < tet_around_1->end(); ++i) {
			around_tet->emplace_back(obj_1);
			around_tet->emplace_back(*i);
		}
	}
	else if (obj_1_collider) {
		for (auto i = triangle_around_0->begin(); i < triangle_around_0->end(); ++i) {
			around_triangle->emplace_back(obj_0);
			around_triangle->emplace_back(*i);
		}
		for (auto i = edge_around_0->begin(); i < edge_around_0->end(); ++i) {
			around_edge->emplace_back(obj_0);
			around_edge->emplace_back(*i);
		}
		for (auto i = tet_around_0->begin(); i < tet_around_0->end(); ++i) {
			around_tet->emplace_back(obj_0);
			around_tet->emplace_back(*i);
		}
	}
	else {
		comparePrimitiveAroundPrimitveTogether(triangle_around_0, triangle_around_1, obj_0, obj_1,
			around_triangle);

		comparePrimitiveAroundPrimitveTogether(edge_around_0, edge_around_1, obj_0, obj_1,
			around_edge);

		comparePrimitiveAroundPrimitveTogether(tet_around_0, tet_around_1, obj_0, obj_1,
			around_tet);
	}
	getCollisionBlockCollisionHessianTest(Hessian, grad, around_triangle, around_edge, collision_stiffness,
		unfixed_pair_vertex_index, unfixed_num, d_hat_2, vertex_index_surface.data(), obj_0, obj_1);

	getCollisionBlockTetHessian(Hessian, grad, around_tet, arap_stiffness, unfixed_pair_vertex_index, unfixed_num);
}


void XPBD_IPC::getCollisionPairHessian(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_0, unsigned int obj_1,
	double collision_stiffness, std::vector<unsigned int>* triangle_around_0, std::vector<unsigned int>* triangle_around_1,
	std::vector<unsigned int>* edge_around_0, std::vector<unsigned int>* edge_around_1, 
	std::vector<unsigned int>* tet_around_0, std::vector<unsigned int>* tet_around_1,
	double d_hat_2, int* unfixed_pair_vertex_index, int unfixed_num, std::vector<unsigned int>* around_triangle, std::vector<unsigned int>* around_edge,
	std::vector<unsigned int>* around_tet, double arap_stiffness, bool obj_0_collider, bool obj_1_collider)
{
	if (obj_0_collider) {
		for (auto i = triangle_around_1->begin(); i < triangle_around_1->end(); ++i) {
			around_triangle->emplace_back(obj_1);
			around_triangle->emplace_back(*i);
		}
		for (auto i = edge_around_1->begin(); i < edge_around_1->end(); ++i) {
			around_edge->emplace_back(obj_1);
			around_edge->emplace_back(*i);
		}
		for (auto i = tet_around_1->begin(); i < tet_around_1->end(); ++i) {
			around_tet->emplace_back(obj_1);
			around_tet->emplace_back(*i);
		}
	}
	else if(obj_1_collider) {
		for (auto i = triangle_around_0->begin(); i < triangle_around_0->end(); ++i) {
			around_triangle->emplace_back(obj_0);
			around_triangle->emplace_back(*i);
		}
		for (auto i = edge_around_0->begin(); i < edge_around_0->end(); ++i) {
			around_edge->emplace_back(obj_0);
			around_edge->emplace_back(*i);
		}
		for (auto i = tet_around_0->begin(); i < tet_around_0->end(); ++i) {
			around_tet->emplace_back(obj_0);
			around_tet->emplace_back(*i);
		}
	}
	else {
		comparePrimitiveAroundPrimitveTogether(triangle_around_0, triangle_around_1, obj_0, obj_1,
			around_triangle);
		comparePrimitiveAroundPrimitveTogether(edge_around_0, edge_around_1, obj_0, obj_1,
			around_edge);
		comparePrimitiveAroundPrimitveTogether(tet_around_0, tet_around_1, obj_0, obj_1,
			around_tet);
	}

	getCollisionBlockCollisionHessian(Hessian, grad, around_triangle, around_edge, collision_stiffness,
		unfixed_pair_vertex_index, unfixed_num, d_hat_2, vertex_index_surface.data(), obj_0, obj_1);

	//FEM::SPDprojection(Hessian);

	getCollisionBlockTetHessian(Hessian, grad, around_tet, arap_stiffness, unfixed_pair_vertex_index, unfixed_num);
}





void XPBD_IPC::comparePrimitiveAroundPrimitveTogether(std::vector<unsigned int>* primitive_around_1, std::vector<unsigned int>* primitive_around_2,
	unsigned int obj_1, unsigned int obj_2, std::vector<unsigned int>* primitive_together)
{
	primitive_together->reserve((primitive_around_1->size() + primitive_around_2->size()) << 1);
	if (obj_1 != obj_2) {
		for (auto i = primitive_around_1->begin(); i < primitive_around_1->end(); ++i) {
			primitive_together->emplace_back(obj_1);
			primitive_together->emplace_back(*i);
		}
		for (auto i = primitive_around_2->begin(); i < primitive_around_2->end(); ++i) {
			primitive_together->emplace_back(obj_2);
			primitive_together->emplace_back(*i);
		}
	}
	else {
		std::set_union(primitive_around_1->begin(), primitive_around_1->end(),
			primitive_around_2->begin(), primitive_around_2->end(), std::back_inserter(*primitive_together));
		primitive_together->resize(primitive_together->size() << 1);
		for (int i = primitive_together->size() - 1; i > -1; i -= 2) {
			primitive_together->data()[i] = primitive_together->data()[i >> 1];
			primitive_together->data()[i - 1] = obj_1;
		}
	}
}



void XPBD_IPC::solveBlockWithPair(double dt,
	std::array<double, 3>** record_vertex_position,
	int** record_vertex_num, int* unfixed_pair_vertex_index, int unfixed_num,double* common_grad,
	std::unordered_map<std::array<unsigned int, 2>, StoreHessianWithOrderInConstraint, pair_hash>& collision_hessian,
	double* floor_map, bool only_solve_collision_pair, double& max_dis)
{
	MatrixXd Hessian;
	VectorXd grad;

	int real_unfiex_num = unfixed_num >> 1;

	grad.resize(3 * real_unfiex_num);
	grad.setZero();
	Hessian.resize(3 * real_unfiex_num, 3 * real_unfiex_num);
	Hessian.setZero();

	if (!only_solve_collision_pair) {
		auto k = tet_hessian[cloth->size()].end();
		for (int i = 0; i < unfixed_num; i += 2) {
			for (int j = 0; j < unfixed_num; j += 2) {
				if (unfixed_pair_vertex_index[i] != unfixed_pair_vertex_index[j]) {
					continue;
				}
				k = tet_hessian[unfixed_pair_vertex_index[i]].find(std::array{ unfixed_pair_vertex_index[i + 1],unfixed_pair_vertex_index[j + 1] });
				if (k != tet_hessian[unfixed_pair_vertex_index[i]].end()) {
					Hessian.data()[3 * ((j * Hessian.cols() + i) >> 1)] = k->second;
					Hessian.data()[3 * ((j * Hessian.cols() + i) >> 1) + Hessian.cols() + 1] = k->second;
					Hessian.data()[3 * ((j * Hessian.cols() + i) >> 1) + Hessian.cols() + Hessian.cols() + 2] = k->second;
				}				
			}			
		}
	}

	for (int i = 0; i < unfixed_num; i += 2) {
		memcpy(grad.data() + 3 * (i >> 1), common_grad + 3 * (vertex_index_prefix_sum_obj[unfixed_pair_vertex_index[i]] + unfixed_pair_vertex_index[i + 1]), 24);
	}


	if (perform_collision) {
		auto k_ = collision_hessian.end();
		double* add_of_hessian;
		double* hessian_in_collision;

		unsigned int index_0, index_1;
		for (int i = 0; i < unfixed_num; i+=2) {
			index_0 = vertex_index_prefix_sum_obj[unfixed_pair_vertex_index[i]] + unfixed_pair_vertex_index[i + 1];
			for (int j = 0; j < unfixed_num; j+=2) {
				index_1 = vertex_index_prefix_sum_obj[unfixed_pair_vertex_index[j]] + unfixed_pair_vertex_index[j + 1];
				k_ = collision_hessian.find(std::array{index_0,index_1	});
				if (k_ != collision_hessian.end()) {
					add_of_hessian = Hessian.data() + 3 *( (j * Hessian.cols() + i)>>1);
					hessian_in_collision = k_->second.hessian.data();
					for (int m = 0; m < 3; ++m) {
						*add_of_hessian += *hessian_in_collision;
						*(add_of_hessian + 1) += *(hessian_in_collision + 1);
						*(add_of_hessian + 2) += *(hessian_in_collision + 2);
						hessian_in_collision += 3;
						add_of_hessian += Hessian.cols();
					}
					if (i != j) {
						add_of_hessian = Hessian.data() + 3 * ((i * Hessian.cols() + j) >> 1);
						hessian_in_collision = k_->second.hessian.data();
						for (int m = 0; m < 3; ++m) {
							*add_of_hessian += *hessian_in_collision;
							*(add_of_hessian + 1) += *(hessian_in_collision + 3);
							*(add_of_hessian + 2) += *(hessian_in_collision + 6);
							hessian_in_collision ++;
							add_of_hessian += Hessian.cols();
						}
					}
				}
			}			
		}

		if (floor->exist) {
			for (int i = 0; i < unfixed_num; i += 2) {
				Hessian.data()[(3 * (i >> 1) + floor->dimension) * (Hessian.cols() + 1)] += floor_map[vertex_index_prefix_sum_obj[unfixed_pair_vertex_index[i]] + unfixed_pair_vertex_index[i + 1]];		
			}
		}
	}

	//if (unfixed_num == 4) {
	//	if (unfixed_pair_vertex_index[1] == 1973 || unfixed_pair_vertex_index[3] == 1973) {
	//		std::cout << "==" << unfixed_pair_vertex_index[1] << std::endl;
	//		std::cout << grad.transpose() << std::endl;
	//		std::cout << Hessian << std::endl;
	//	}
	//}


	double mass_dt_2;
	int obj_index_, vertex_index_;

	int i;
	for (int j = 0; j < unfixed_num; j += 2) {
		i = (j >> 1);
		vertex_index_ = unfixed_pair_vertex_index[j + 1];
		obj_index_ = unfixed_pair_vertex_index[j];
		mass_dt_2 = mass[obj_index_][vertex_index_] / (dt * dt);
		Hessian(3 * i, 3 * i) += mass_dt_2;
		Hessian(3 * i + 1, 3 * i + 1) += mass_dt_2;
		Hessian(3 * i + 2, 3 * i + 2) += mass_dt_2;
		grad.data()[3 * i] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][0] - sn[obj_index_][vertex_index_][0]);
		grad.data()[3 * i + 1] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][1] - sn[obj_index_][vertex_index_][1]);
		grad.data()[3 * i + 2] += mass_dt_2 * (vertex_position[obj_index_][vertex_index_][2] - sn[obj_index_][vertex_index_][2]);
	}

	LLT <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(grad);




	//if (unfixed_num == 4) {
	//	if (unfixed_pair_vertex_index[1] == 1973 || unfixed_pair_vertex_index[3] == 1973) {
	//		std::cout << std::setprecision(12) << mass_dt_2<<" "<< sn[obj_index_][ unfixed_pair_vertex_index[1]][0]<<" "<< sn[obj_index_][ unfixed_pair_vertex_index[1]][1]<<" "<< sn[obj_index_][ unfixed_pair_vertex_index[1]][2] << std::endl;
	//		std::cout << std::setprecision(12) << mass_dt_2<<" "<< sn[obj_index_][ unfixed_pair_vertex_index[3]][0]<<" "<< sn[obj_index_][ unfixed_pair_vertex_index[3]][1]<<" "<< sn[obj_index_][ unfixed_pair_vertex_index[3]][2] << std::endl;
	//		//std::cout << "==" << unfixed_pair_vertex_index[1] << std::endl;
	//		//std::cout << grad.transpose() << std::endl;
	//		//std::cout << Hessian << std::endl;
	//		std::cout << result.transpose() << std::endl;
	//	}
	//}
	//
		//std::cout << unfixed_pair_vertex_index[1] << std::endl;
		//std::cout << vertex_position[obj_index_][vertex_index_][1] - *(result.data() + 3 * (i >> 1) + 1) - sn[obj_index_][vertex_index_][1] << std::endl;
	


	double dis;
	for (int i = 0; i < unfixed_num; i += 2) {
		vertex_index_ = unfixed_pair_vertex_index[i + 1];
		obj_index_ = unfixed_pair_vertex_index[i];
		if (record_vertex_num[obj_index_][vertex_index_]) {
			record_vertex_position[obj_index_][vertex_index_][0] += vertex_position[obj_index_][vertex_index_][0] - *(result.data() + 3 * (i >> 1));
			record_vertex_position[obj_index_][vertex_index_][1] += vertex_position[obj_index_][vertex_index_][1] - *(result.data() + 3 * (i >> 1) + 1);
			record_vertex_position[obj_index_][vertex_index_][2] += vertex_position[obj_index_][vertex_index_][2] - *(result.data() + 3 * (i >> 1) + 2);
		}
		else {
			SUB(record_vertex_position[obj_index_][vertex_index_], vertex_position[obj_index_][vertex_index_], (result.data() + 3 * (i>>1)));
		}
		record_vertex_num[obj_index_][vertex_index_] ++;

		dis = abs(vertex_position[obj_index_][vertex_index_][0] - *(result.data() + 3 * (i >> 1)) - record_collision_free_vertex_position[obj_index_][vertex_index_][0]);
		if (dis > max_dis) {
			max_dis = dis;
		}
		dis = abs(vertex_position[obj_index_][vertex_index_][1] - *(result.data() + 3 * (i >> 1) + 1) - record_collision_free_vertex_position[obj_index_][vertex_index_][1]);
		if (dis > max_dis) {
			max_dis = dis;
		}
		dis = abs(vertex_position[obj_index_][vertex_index_][2] - *(result.data() + 3 * (i >> 1) + 2) - record_collision_free_vertex_position[obj_index_][vertex_index_][2]);
		if (dis > max_dis) {
			max_dis = dis;
		}
	}
}


void XPBD_IPC::getVTUnfixedPairIndex(int* unfixed_pair_vertex_index, int& unfixed_num, unsigned int vertex_obj_no,
	unsigned int vertex_index, unsigned int triangle_obj_No, unsigned int triangle_index)
{
	memset(unfixed_pair_vertex_index, 0xff, 32);
	
	if (!(*is_vertex_fixed[vertex_obj_no])[vertex_index]) {
		unfixed_pair_vertex_index[unfixed_num] = vertex_obj_no;
		unfixed_pair_vertex_index[unfixed_num + 1] = vertex_index;
		unfixed_num += 2;
	}

	int* tri_indices = triangle_indices[triangle_obj_No][triangle_index].data();
	for (unsigned int i = 0; i < 3; ++i) {
		if (!(*is_vertex_fixed[triangle_obj_No])[tri_indices[i]]) {
			unfixed_pair_vertex_index[unfixed_num] = triangle_obj_No;
			unfixed_pair_vertex_index[unfixed_num + 1] = tri_indices[i];
			unfixed_num += 2;
		}
	}
}

//type 0v, 1e, 2t
void XPBD_IPC::solveCollisionWithColliderBlock(double dt, unsigned int obj_no, unsigned int index, int type, std::array<double, 3>** record_vertex_position,
	int** record_vertex_num, bool only_solve_collision_pair, int thread_No)
{
	int unfixed_pair_vertex_index[6];
	memset(unfixed_pair_vertex_index, 0, 24);
	int unfixed_num = 0;
	switch (type)
	{
		//v
	case 0:
		if (!(*is_vertex_fixed[obj_no])[index]) {
			unfixed_pair_vertex_index[0] = obj_no;
			unfixed_pair_vertex_index[1] = index;
			unfixed_num = 2;
		}
		break;
		//e
	case 1: {
		unsigned int* edge_vertex_0 = edge_vertices[obj_no] + (index << 1);
		if ((!(*is_vertex_fixed[obj_no])[*edge_vertex_0])) {
			unfixed_pair_vertex_index[unfixed_num] = obj_no;
			unfixed_pair_vertex_index[unfixed_num + 1] = *edge_vertex_0;
			unfixed_num += 2;
		}
		if ((!(*is_vertex_fixed[obj_no])[edge_vertex_0[1]])) {
			unfixed_pair_vertex_index[unfixed_num] = obj_no;
			unfixed_pair_vertex_index[unfixed_num + 1] = edge_vertex_0[1];
			unfixed_num += 2;
		}
	}
		break;
		//t
	case 2:
		int* tri_indices = triangle_indices[obj_no][index].data();
		for (unsigned int i = 0; i < 3; ++i) {
			if (!(*is_vertex_fixed[obj_no])[tri_indices[i]]) {
				unfixed_pair_vertex_index[unfixed_num] = obj_no;
				unfixed_pair_vertex_index[unfixed_num + 1] = tri_indices[i];
				unfixed_num += 2;
			}
		}
		break;
	}

	if (unfixed_num == 0) {
		return;
	}
	solveBlockWithPair(
		dt, record_vertex_position, record_vertex_num, unfixed_pair_vertex_index, unfixed_num, common_grad[0].data(),
		common_hessian, floor_hessian.data(), only_solve_collision_pair, max_dis_record[thread_No]);
}

void XPBD_IPC::solveVT_Block(unsigned int vertex_obj_no, unsigned int vertex_index, unsigned int triangle_obj_No, unsigned int triangle_index,
	double dt, 
	std::array<double, 3>** record_vertex_position,
	int** record_vertex_num, bool only_solve_collision_pair, int thread_No)
{

	int unfixed_pair_vertex_index[8];
	int unfixed_num = 0;

	getVTUnfixedPairIndex(unfixed_pair_vertex_index, unfixed_num, vertex_obj_no, vertex_index, triangle_obj_No, triangle_index);

	if (unfixed_num == 0) {
		return;
	}
	solveBlockWithPair(
		dt, record_vertex_position, record_vertex_num, unfixed_pair_vertex_index, unfixed_num,common_grad[0].data(),
		common_hessian,floor_hessian.data(), only_solve_collision_pair, max_dis_record[thread_No]);
}


void XPBD_IPC::getEEUnfixedPairIndex(int* unfixed_pair_vertex_index,int& unfixed_num, unsigned int obj_No_0, 
	unsigned int primitive_0_index, unsigned int obj_No_1, unsigned int primitive_1_index)
{
	memset(unfixed_pair_vertex_index, 0xff, 32);
	unsigned int* edge_vertex_0 = edge_vertices[obj_No_0] + (primitive_0_index << 1);

	//bool has[4];
	//memset(has, 0, 4);
	//for (int i = 0; i < *hessian_record_index; ++i) {
	//	has[hessian_record_index[i + 1]] = true;
	//}
	if ((!(*is_vertex_fixed[obj_No_0])[*edge_vertex_0])) {
		unfixed_pair_vertex_index[unfixed_num] = obj_No_0;
		unfixed_pair_vertex_index[unfixed_num + 1] = *edge_vertex_0;
		unfixed_num += 2;
	}
	if ((!(*is_vertex_fixed[obj_No_0])[edge_vertex_0[1]])) {
		unfixed_pair_vertex_index[unfixed_num] = obj_No_0;
		unfixed_pair_vertex_index[unfixed_num + 1] = edge_vertex_0[1];
		unfixed_num += 2;
	}


	unsigned int* edge_vertex_1 = edge_vertices[obj_No_1] + (primitive_1_index << 1);

	if ((!(*is_vertex_fixed[obj_No_1])[*edge_vertex_1])) {
		unfixed_pair_vertex_index[unfixed_num] = obj_No_1;
		unfixed_pair_vertex_index[unfixed_num + 1] = *edge_vertex_1;
		unfixed_num += 2;
	}
	if ((!(*is_vertex_fixed[obj_No_1])[edge_vertex_1[1]])) {
		unfixed_pair_vertex_index[unfixed_num] = obj_No_1;
		unfixed_pair_vertex_index[unfixed_num + 1] = edge_vertex_1[1];
		unfixed_num += 2;
	}

}




void XPBD_IPC::solveEE_Block(unsigned int obj_No_0, unsigned int primitive_0_index, unsigned int obj_No_1, unsigned int primitive_1_index,
	double dt,
	std::array<double, 3>** record_vertex_position, 
	int** record_vertex_num, bool only_solve_collision_pair, int thread_No)
{
	MatrixXd Hessian;
	VectorXd grad;
	int unfixed_pair_vertex_index[8];
	int unfixed_num = 0;

	
	getEEUnfixedPairIndex(unfixed_pair_vertex_index, unfixed_num, obj_No_0, primitive_0_index, obj_No_1, primitive_1_index);

	if (unfixed_num == 0) {
		return;
	}

	solveBlockWithPair(
		dt, record_vertex_position, record_vertex_num, unfixed_pair_vertex_index, unfixed_num,common_grad[0].data(),
		common_hessian, floor_hessian.data(), only_solve_collision_pair, max_dis_record[thread_No]);
}





void XPBD_IPC::getHessianForCollisionBlock(MatrixXd& Hessian, VectorXd& grad, unsigned int obj_0, unsigned int obj_1,
	std::vector<unsigned int>* triangle_around_0, std::vector<unsigned int>* triangle_around_1,
	std::vector<unsigned int>* edge_around_0, std::vector<unsigned int>* edge_around_1,
	std::vector<unsigned int>* tet_around_0, std::vector<unsigned int>* tet_around_1,
	int* unfixed_pair_vertex_index, int unfixed_num, std::vector<unsigned int>* around_triangle, std::vector<unsigned int>* around_edge,
	std::vector<unsigned int>* around_tet)
{
	comparePrimitiveAroundPrimitveTogether(triangle_around_0, triangle_around_1, obj_0, obj_1,
		around_triangle);

	comparePrimitiveAroundPrimitveTogether(edge_around_0, edge_around_1, obj_0, obj_1,
		around_edge);

	comparePrimitiveAroundPrimitveTogether(tet_around_0, tet_around_1, obj_0, obj_1,
		around_tet);

	getCollisionBlockCollisionHessianFromRecord(Hessian, grad, around_triangle, around_edge, obj_0, obj_1,
		unfixed_pair_vertex_index, unfixed_num, Hessian, grad);


	//for (auto i = around_tet->begin(); i < around_tet->end(); i += 2) {
	//	getARAPCollisionHessianForPairFromRecord(Hessian, grad, *i, tet_indices[*i][*(i + 1)].data(),
	//		unfixed_pair_vertex_index, unfixed_num,  
	//		store_tet_arap_hessian.data()+16* (prefix_sum_of_every_tet_index[*i - cloth->size()]+*(i+1)),
	//		store_tet_arap_grad.data() + 12 * (prefix_sum_of_every_tet_index[*i - cloth->size()]+ *(i + 1)));
	//}


}


void XPBD_IPC::getARAPCollisionHessianForPairFromRecord(MatrixXd& Hessian, VectorXd& grad, int tet_obj, int* tet_vertex, int* tet_unfixed_vertex_indices,
	int unfixed_tet_vertex_num, double* Hessian_, double* grad_)
{
	int triangle_vertex_order_in_tet[4];
	memset(triangle_vertex_order_in_tet, 0xff, 16);
	checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, tet_vertex, tet_obj, triangle_vertex_order_in_tet, 4);
	setHessianFromTetHessian(Hessian, grad.data(), Hessian_, grad_, triangle_vertex_order_in_tet);
}


void XPBD_IPC::getCollisionBlockCollisionHessianFromRecord(MatrixXd& Hessian, VectorXd& grad, std::vector<unsigned int>* triangles,
	std::vector<unsigned int>* edges,
	unsigned int obj_No_0, unsigned int obj_No_1, int* pair_actual_unfixed_vertex_indices,
	int unfixed_vertex_num,MatrixXd& test_m, VectorXd& grad_m)
{
	////vt
	//int vertex_surface_index;
	//if (has_collider) {
	//	for (int i = 0; i < unfixed_vertex_num; i += 2) {
	//		vertex_surface_index =vertex_index_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]];
	//		if (vertex_surface_index == -1) {
	//			continue;
	//		}
	//		getVTCollisionHessainForPairFromRecord(Hessian, grad,
	//			collision.vertex_triangle_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * vertex_surface_index,
	//			collision.vertex_triangle_pair_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_surface_index],
	//			i >> 1, obj_No_0, obj_No_1, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
	//			collision.vt_hessian_record.data()+ 144* collision.vertex_triangle_pair_num_record_prefix_sum[pair_actual_unfixed_vertex_indices[i]][vertex_surface_index],
	//			collision.vt_grad_record.data()+ 12* collision.vertex_triangle_pair_num_record_prefix_sum[pair_actual_unfixed_vertex_indices[i]][vertex_surface_index],
	//			collision.vt_hessian_record_index.data()+5* collision.vertex_triangle_pair_num_record_prefix_sum[pair_actual_unfixed_vertex_indices[i]][vertex_surface_index],
	//			collision.vt_colldier_hessian_record.data() + 9*(vertex_num_on_surface_prefix_sum[pair_actual_unfixed_vertex_indices[i]] + vertex_surface_index),
	//			collision.vt_colldier_grad_record.data() + 3*(vertex_num_on_surface_prefix_sum[pair_actual_unfixed_vertex_indices[i]] + vertex_surface_index));
	//	}
	//}
	//else {
	//	double emp[1] = { 0 };
	//	for (int i = 0; i < unfixed_vertex_num; i += 2) {
	//		vertex_surface_index = vertex_index_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]];
	//		if (vertex_surface_index == -1) {
	//			continue;
	//		}
	//		getVTCollisionHessainForPairFromRecord(Hessian, grad,
	//			collision.vertex_triangle_pair_by_vertex[pair_actual_unfixed_vertex_indices[i]] + collision.close_vt_pair_num * vertex_surface_index,
	//			collision.vertex_triangle_pair_num_record[pair_actual_unfixed_vertex_indices[i]][vertex_surface_index],
	//			i >> 1, obj_No_0, obj_No_1, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
	//			collision.vt_hessian_record.data() + 144 * collision.vertex_triangle_pair_num_record_prefix_sum[pair_actual_unfixed_vertex_indices[i]][vertex_surface_index],
	//			collision.vt_grad_record.data() + 12 * collision.vertex_triangle_pair_num_record_prefix_sum[pair_actual_unfixed_vertex_indices[i]][vertex_surface_index],
	//			collision.vt_hessian_record_index.data() + 5 * collision.vertex_triangle_pair_num_record_prefix_sum[pair_actual_unfixed_vertex_indices[i]][vertex_surface_index],
	//			emp, emp);
	//	}
	//}

	////tv
	//int* triangle_;
	//if (has_collider) {
	//	for (auto i = triangles->begin(); i < triangles->end(); i += 2) {
	//		triangle_ = triangle_indices[*i][*(i + 1)].data();
	//		getTVCollisionHessainForPairFromRecord(test_m, grad_m,
	//			collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i + 1)),
	//			collision.triangle_vertex_pair_num_record[*i][*(i + 1)], *i,
	//			obj_No_0, obj_No_1, triangle_, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
	//			collision.triangle_vertex_collider_pair_num_record[*i][*(i + 1)],
	//			collision.tv_colldier_hessian_record_index.data()+ 4 * collision.triangle_vertex_collider_num_record_prefix_sum[*i][*(i+1)],
	//			collision.tv_colldier_hessian_record.data()+ 81 * collision.triangle_vertex_collider_num_record_prefix_sum[*i][*(i+1)],
	//			collision.tv_colldier_grad_record.data()+ 9 * collision.triangle_vertex_collider_num_record_prefix_sum[*i][*(i+1)]);
	//	}
	//}
	//else {
	//	double emp[1] = { 0 };
	//	int emp_[1] = { 0 };
	//	for (auto i = triangles->begin(); i < triangles->end(); i += 2) {
	//		triangle_ = triangle_indices[*i][*(i + 1)].data();
	//		getTVCollisionHessainForPairFromRecord(test_m, grad_m,
	//			collision.triangle_vertex_pair_by_triangle[*i] + collision.close_tv_pair_num * (*(i + 1)),
	//			collision.triangle_vertex_pair_num_record[*i][*(i + 1)], *i,
	//			obj_No_0, obj_No_1, triangle_, pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
	//			0,emp_,emp,emp);
	//	}
	//}

	////ee
	//unsigned int* edge_;
	//if (has_collider) {
	//	for (auto i = edges->begin(); i < edges->end(); i += 2) {
	//		edge_ = edge_vertices[*i] + ((*(i + 1)) << 1);
	//		getEECollisionHessainForPairFromRecord(Hessian, grad, *i, *(i+1), edge_,
	//			collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i + 1)),
	//			collision.edge_edge_pair_number_record[*i][*(i + 1)],
	//			pair_actual_unfixed_vertex_indices, unfixed_vertex_num, 
	//			obj_No_0, obj_No_1,
	//			i - edges->begin(), edges, collision.edge_edge_collider_pair_by_edge[*i] + collision.close_ee_collider_pair_num * (*(i + 1)),
	//			collision.edge_edge_collider_pair_num_record[*i][*(i + 1)], 
	//			collision.ee_hessian_record.data()+ 144 * collision.edge_edge_pair_num_record_prefix_sum[*i][*(i+1)],
	//			collision.ee_grad_record.data()+ 12 * collision.edge_edge_pair_num_record_prefix_sum[*i][*(i+1)],
	//			collision.ee_hessian_record_index.data()+ 5 * collision.edge_edge_pair_num_record_prefix_sum[*i][*(i+1)],
	//			collision.ee_collider_hessian_record_index.data() + 3 * collision.edge_edge_collider_num_record_prefix_sum[*i][*(i + 1)],
	//			collision.ee_collider_hessian_record.data() + 36 * collision.edge_edge_collider_num_record_prefix_sum[*i][*(i + 1)],
	//			collision.ee_collider_grad_record.data() + 6 * collision.edge_edge_collider_num_record_prefix_sum[*i][*(i + 1)]	);
	//	}
	//}
	//else {
	//	 int emp_[1] = { 0 };
	//	double emp[1] = { 0 };
	//	for (auto i = edges->begin(); i < edges->end(); i += 2) {
	//		edge_ = edge_vertices[*i] + ((*(i + 1)) << 1);
	//		getEECollisionHessainForPairFromRecord(Hessian, grad, *i, *(i + 1), edge_,
	//			collision.edge_edge_pair_by_edge[*i] + collision.close_ee_pair_num * (*(i + 1)),
	//			collision.edge_edge_pair_number_record[*i][*(i + 1)],
	//			pair_actual_unfixed_vertex_indices, unfixed_vertex_num,
	//			obj_No_0, obj_No_1,
	//			i - edges->begin(), edges, 
	//			collision.edge_edge_collider_pair_by_edge[*i] + collision.close_ee_collider_pair_num * (*(i + 1)),
	//			collision.edge_edge_collider_pair_num_record[*i][*(i + 1)], 
	//			collision.ee_hessian_record.data() + 144 * collision.edge_edge_pair_num_record_prefix_sum[*i][*(i + 1)],
	//			collision.ee_grad_record.data() + 12 * collision.edge_edge_pair_num_record_prefix_sum[*i][*(i + 1)],
	//			collision.ee_hessian_record_index.data() + 5 * collision.edge_edge_pair_num_record_prefix_sum[*i][*(i + 1)],
	//			emp_, emp, emp);
	//	}
	//}


	//if (floor->exist) {
	//	int locate;

	//	for (int i = 0; i < unfixed_vertex_num; i += 2) {
	//		vertex_surface_index = vertex_index_surface[pair_actual_unfixed_vertex_indices[i]][pair_actual_unfixed_vertex_indices[i + 1]];
	//		if (vertex_surface_index == -1) {
	//			continue;
	//		}
	//		locate = vertex_num_on_surface_prefix_sum[pair_actual_unfixed_vertex_indices[i]] + vertex_surface_index;
	//		grad[3 *(i>>1) + floor->dimension] += collision.floor_grad_record[locate];
	//		Hessian.data()[(3 * (i>>1) + floor->dimension) * (3 * (unfixed_vertex_num>>1 )+ 1)] += collision.floor_hessian_record[locate];
	//	}
	//}
}


void XPBD_IPC::getEECollisionHessainForPairFromRecord(MatrixXd& Hessian, VectorXd& grad, unsigned int ee_obj_No, 
	unsigned int edge_index,
	unsigned int* edge_vertex_index,
	unsigned int* EE, int num, int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, 
	unsigned int obj_No_0, unsigned int obj_No_1,
	int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet,	unsigned int* EE_collider, int num_collider,
	double* ee_hessian_record_, double* ee_grad_record_, int* ee_hessian_record_index_,
	int* ee_collider_hessian_record_index, double* ee_collider_hessian_record, double* ee_collider_grad_record)
{
	//unsigned int* edge_vertex;
	//int vertex_order_in_tet[4];
	//memset(vertex_order_in_tet, 0xff, 8);
	//checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex_index, ee_obj_No, vertex_order_in_tet);

	//int back_up_order[2];
	//memcpy(back_up_order, vertex_order_in_tet, 8);

	//double* ee_hessian_record; double* ee_grad_record; int* ee_hessian_record_index;

	//int pair_index;

	//for (int i = 0; i < num; i += 3) {
	//	edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
	//	if (EE[i] > ee_obj_No || (EE[i] == ee_obj_No && EE[i + 1] > edge_index)) {
	//		ee_hessian_record_index = ee_hessian_record_index_ + i / 3 * 5;
	//		ee_grad_record = ee_grad_record_ + i / 3 * 12;
	//		ee_hessian_record = ee_hessian_record_ + i / 3 * 144;
	//		memset(vertex_order_in_tet + 2, 0xff, 8);
	//		if (EE[i] == obj_No_0 || EE[i] == obj_No_1) {
	//			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, EE[i + 1], EE[i])) {
	//				continue;
	//			}
	//			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex, EE[i], vertex_order_in_tet + 2);
	//		}


	//	}
	//	else {
	//		memcpy(vertex_order_in_tet + 2, back_up_order, 8);
	//		memset(vertex_order_in_tet, 0xff, 8);
	//		pair_index = collision.edge_edge_pair_num_record_prefix_sum[EE[i]][EE[i + 1]] + EE[i + 2] / 3;
	//		ee_hessian_record_index = collision.ee_hessian_record_index.data() + 5 * pair_index;
	//		ee_grad_record = collision.ee_grad_record.data() + 12 * pair_index;
	//		ee_hessian_record = collision.ee_hessian_record.data() + 144 * pair_index;
	//		if (EE[i] == obj_No_0 || EE[i] == obj_No_1) {
	//			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, EE[i + 1], EE[i])) {
	//				continue;
	//			}
	//			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, edge_vertex, EE[i], vertex_order_in_tet);
	//		}
	//	}
	//	if ((*ee_hessian_record_index) != 0) {
	//		setTetHessianFromBarrierHessian(Hessian, grad.data(), ee_hessian_record, ee_grad_record, vertex_order_in_tet,
	//			ee_hessian_record_index + 1, *ee_hessian_record_index);
	//	}
	//}

	//if (has_collider) {
	//	memcpy(vertex_order_in_tet, back_up_order, 8);
	//	memset(vertex_order_in_tet + 2, 0xff, 8);
	//	for (int i = 0; i < num_collider; i += 2) {
	//		if ((*ee_collider_hessian_record_index) != 0) {
	//			setTetHessianFromBarrierHessian(Hessian, grad.data(), ee_collider_hessian_record, ee_collider_grad_record, vertex_order_in_tet,
	//				ee_collider_hessian_record_index + 1, *ee_collider_hessian_record_index);
	//		}
	//		ee_collider_hessian_record_index += 3;
	//		ee_collider_grad_record += 6;
	//		ee_collider_hessian_record += 36;
	//	}
	//}
}



void XPBD_IPC::getTVCollisionHessainForPairFromRecord(MatrixXd& Hessian, VectorXd& grad, 
	unsigned int* TV, int num, unsigned int tri_obj, unsigned int obj_No_0, unsigned int obj_No_1, int* triangle_indices,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, int collider_num,
	int* tv_collider_hessian_record_index, double* tv_collider_hessian_record, double* tv_collider_grad_record)
{
	//int triangle_vertex_order_in_tet[4];
	//memset(triangle_vertex_order_in_tet, 0xff, 16);
	//checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_indices, tri_obj, triangle_vertex_order_in_tet + 1, 3);
	//int pair_index;
	//int* vt_record_index;

	//for (int i = 0; i < num; i += 3) {
	//	if (TV[i] == obj_No_0 || TV[i] == obj_No_1) {
	//		if (vertexInPair(unfixed_tet_vertex_num, TV[i + 1], tet_unfixed_vertex_indices, TV[i])) {
	//			continue;
	//		}
	//	}

	//	if (TV[i] < cloth->size()) {
	//		pair_index = collision.vertex_triangle_pair_num_record_prefix_sum[TV[i]][TV[i + 1]] + TV[i + 2];
	//	}
	//	else {
	//		pair_index = collision.vertex_triangle_pair_num_record_prefix_sum[TV[i]][vertex_index_surface[TV[i]][TV[i + 1]]] + TV[i + 2];
	//	}
	//	vt_record_index = collision.vt_hessian_record_index.data() + 5 * pair_index;
	//	if ((*vt_record_index) != 0) {
	//		setTetHessianFromBarrierHessian(Hessian, grad.data(), collision.vt_hessian_record.data() + 144 * pair_index,
	//			collision.vt_grad_record.data() + 12 * pair_index, triangle_vertex_order_in_tet,
	//			vt_record_index + 1, *vt_record_index);
	//	}
	//}
	//if (has_collider) {
	//	for (int i = 0; i < collider_num; i += 2) {
	//		if ((*tv_collider_hessian_record_index) != 0) {
	//			setTetHessianFromBarrierHessian(Hessian, grad.data(), tv_collider_hessian_record, tv_collider_grad_record, triangle_vertex_order_in_tet,
	//				tv_collider_hessian_record_index + 1, *tv_collider_hessian_record_index);
	//		}
	//		tv_collider_hessian_record_index += 4;
	//		tv_collider_hessian_record += 9;
	//		tv_collider_grad_record += 81;
	//	}
	//}
}


void XPBD_IPC::getVTCollisionHessainForPairFromRecord(MatrixXd& Hessian, VectorXd& grad,
	unsigned int* VT, unsigned int num, unsigned int vertex_order_in_matrix, unsigned int vertex_obj_No, unsigned int tri_obj_No,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num,
	double* vt_hessian_record, double* vt_grad_record, int* vt_hessian_record_index,
	double* vt_collider_hessian_record, double* vt_collider_grad_record)
{
	int* triangle_vertex;
	int triangle_vertex_order_in_tet[4];
	triangle_vertex_order_in_tet[0] = vertex_order_in_matrix;
	for (int i = 0; i < num; i += 2) {
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		memset(triangle_vertex_order_in_tet + 1, 0xff, 12);
		if (VT[i] == tri_obj_No || VT[i] == vertex_obj_No) {
			checkPairIndexInSys(unfixed_tet_vertex_num, tet_unfixed_vertex_indices, triangle_vertex, VT[i], triangle_vertex_order_in_tet + 1, 3);
		}
		if ((*vt_hessian_record_index) != 0) {
			setTetHessianFromBarrierHessian(Hessian, grad.data(), vt_hessian_record, vt_grad_record, triangle_vertex_order_in_tet,
				vt_hessian_record_index + 1, *vt_hessian_record_index);
		}
		vt_hessian_record_index += 5;
		vt_grad_record += 12;
		vt_hessian_record += 144;
	}

	if (has_collider) {
		int vertex_pos = 3 * vertex_order_in_matrix * (Hessian.cols() + 1);
		for (unsigned int i = 0; i < 3; ++i) {
			grad.data()[3 * vertex_order_in_matrix + i] += vt_collider_grad_record[i];
			Hessian.data()[vertex_pos + i] += vt_collider_hessian_record[i];
			Hessian.data()[vertex_pos + Hessian.cols() + i] += vt_collider_hessian_record[i + 3];
			Hessian.data()[vertex_pos + Hessian.cols() + Hessian.cols() + i] += vt_collider_hessian_record[i + 6];
		}
	}
}



void XPBD_IPC::solveTetBlockCollision(std::array<double, 3>* vertex_position, double stiffness, double dt, std::array<int, 4>* indices,
	double* mass,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices,
	double* volume, unsigned int tet_index, std::array<double, 3>* sn, unsigned int* common_vertex_in_order,
	int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num, std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet, double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	int* vertex_index_on_surface, std::unordered_map<std::array<int, 2>, double, pair_hash>& tet_hessian,
	std::unordered_map<std::array<unsigned int, 2>, StoreHessianWithOrderInConstraint, pair_hash>& collision_hessian,
	double* common_grad, std::array<double, 3>* record_vertex_position,
	int* record_vertex_num, unsigned int prefix_sum_vetex_obj, double* floor_map, double& max_dis, std::array<double, 3>* collision_free_pos)//
{
	if (unfixed_vertex_num == 0) {
		return;
	}

	MatrixXd Hessian;
	VectorXd grad;
	grad.resize(3 * unfixed_vertex_num);
	Hessian.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);
	Hessian.setZero();
	grad.setZero();


	auto k = tet_hessian.end();
	for (int i = 0; i < unfixed_vertex_num; ++i) {
		for (int j = 0; j < unfixed_vertex_num; ++j) {
			k = tet_hessian.find(std::array{ tet_actual_unfixed_vertex_indices[i],tet_actual_unfixed_vertex_indices[j] });
			if (k == tet_hessian.end()) {
				std::cout << "error, cannot find the hessian of a tet solveTetBlockCollision()" << std::endl;
			}
			Hessian.data()[3 * (j * Hessian.cols() + i)] = k->second;
			Hessian.data()[3 * (j * Hessian.cols() + i) + Hessian.cols() + 1] = k->second;
			Hessian.data()[3 * (j * Hessian.cols() + i) + Hessian.cols() + Hessian.cols() + 2] = k->second;
		}
		memcpy(grad.data() + 3 * i, common_grad + 3 * (prefix_sum_vetex_obj + tet_actual_unfixed_vertex_indices[i]), 24);
	}

	if (perform_collision) {
		auto k_ = collision_hessian.end();
		double* add_of_hessian;
		double* hessian_in_collision;
		for (int i = 0; i < unfixed_vertex_num; ++i) {
			for (int j = 0; j < unfixed_vertex_num; ++j) {
				k_ = collision_hessian.find(std::array{ prefix_sum_vetex_obj + tet_actual_unfixed_vertex_indices[i], prefix_sum_vetex_obj + tet_actual_unfixed_vertex_indices[j] });
				if (k_ != collision_hessian.end()) {
					if (k_->second.is_update) {
						add_of_hessian = Hessian.data() + 3 * (j * Hessian.cols() + i);
						hessian_in_collision = k_->second.hessian.data();
						for (int m = 0; m < 3; ++m) {
							*add_of_hessian += *hessian_in_collision;
							*(add_of_hessian + 1) += *(hessian_in_collision + 1);
							*(add_of_hessian + 2) += *(hessian_in_collision + 2);
							hessian_in_collision += 3;
							add_of_hessian += Hessian.cols();
						}
						if (i != j) {
							add_of_hessian = Hessian.data() + 3 * (i * Hessian.cols() + j);
							hessian_in_collision = k_->second.hessian.data();
							for (int m = 0; m < 3; ++m) {
								*add_of_hessian += *hessian_in_collision;
								*(add_of_hessian + 1) += *(hessian_in_collision + 3);
								*(add_of_hessian + 2) += *(hessian_in_collision + 6);
								hessian_in_collision++;
								add_of_hessian += Hessian.cols();
							}
						}
					}
				}
			}			
		}
		if (floor->exist) {
			for (int i = 0; i < unfixed_vertex_num; ++i) {
					Hessian.data()[(3 * i + floor->dimension) * (Hessian.cols() + 1)] += floor_map[prefix_sum_vetex_obj + tet_actual_unfixed_vertex_indices[i]];
			}
		}
		//MatrixXd Hessian_collision;
		//VectorXd grad_collision;
		//Hessian_collision.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);
		//Hessian_collision.setZero();
		//grad_collision.resize(3 * unfixed_vertex_num);
		//grad_collision.setZero();
		//second_order_constraint.solveCD_ARAP_blockTest(Hessian_collision, grad_collision, vertex_position, stiffness, A[tet_index],
		//	volume[tet_index], tet_actual_unfixed_vertex_indices,
		//	unfixed_tet_vertex_index, unfixed_vertex_num);
		//for (auto i = neighbor_tet_indices.begin(); i < neighbor_tet_indices.end(); ++i) {
		//	second_order_constraint.solveCertainHessianForNeighborTet(vertex_position, stiffness, A[*i], common_vertex_in_order, indices[*i].data(),
		//		Hessian_collision, volume[*i], grad_collision);		}
		//getCollisionHessian(Hessian_collision, grad_collision, triangle_of_a_tet, edge_of_a_tet, collision_stiffness, obj_No, tet_actual_unfixed_vertex_indices,
		//unfixed_vertex_num,
		//collision.d_hat_2, vertex_index_on_surface, vertex_position);
		//if ((Hessian_collision - Hessian_test).norm() > 1e-10) {
		//	std::cout<< (Hessian_collision - Hessian_test).norm() << ",,,,,,,,,,,,,,,,,,,,,,,,,,,,," << std::endl;
		//	std::cout << Hessian_collision << std::endl;
		//	std::cout << "====" << std::endl;
		//	std::cout << Hessian_test << std::endl;
		//	std::cout << "====" << std::endl;
		//}
		//if ((grad_collision - grad).norm() > 1e-10) {
		//	std::cout << grad_collision.transpose() << std::endl;
		//	std::cout << "============" << std::endl;
		//	std::cout << grad.transpose() << std::endl;
		//}
	}
	//if (Hessian_collision_test.norm() != 0.0) {
	//	std::cout << "show" << std::endl;
	//	std::cout << grad_collision_test.transpose() << std::endl;
	//	std::cout << Hessian_collision_test << std::endl;
	//}
	//if (perform_collision)
	//{
	//	MatrixXd Hessian_collision;
	//	VectorXd grad_collision;
	//	Hessian_collision.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);
	//	Hessian_collision.setZero();
	//	grad_collision.resize(3 * unfixed_vertex_num);
	//	grad_collision.setZero();
	//	getCollisionHessian(Hessian_collision, grad_collision, triangle_of_a_tet, edge_of_a_tet, collision_stiffness, obj_No, tet_actual_unfixed_vertex_indices,
	//		unfixed_vertex_num,
	//		collision.d_hat_2, vertex_index_on_surface, vertex_position);
	//	//if (Hessian_collision.norm() != 0.0) {
	//	//	std::cout << "compare correct " << std::endl;
	//	//	std::cout << grad_collision.transpose() << std::endl;
	//	//	std::cout << Hessian_collision << std::endl;
	//	//}
	//	//FEM::SPDprojection(Hessian_collision);
	//	Hessian += Hessian_collision;
	//	grad += grad_collision;
	//}

	double mass_dt_2;
	int vertex_index;
	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_actual_unfixed_vertex_indices[i];
		mass_dt_2 = mass[vertex_index] / (dt * dt);
		Hessian(3 * i, 3 * i) += mass_dt_2;
		Hessian(3 * i + 1, 3 * i + 1) += mass_dt_2;
		Hessian(3 * i + 2, 3 * i + 2) += mass_dt_2;
		grad.data()[3 * i] += mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
		grad.data()[3 * i + 1] += mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
		grad.data()[3 * i + 2] += mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);
	}
	LLT <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(grad);

	
	double dis;
	
	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_actual_unfixed_vertex_indices[i];
		if (record_vertex_num[vertex_index]) {
			record_vertex_position[vertex_index][0] += vertex_position[vertex_index][0] - *(result.data() + 3 * i);
			record_vertex_position[vertex_index][1] += vertex_position[vertex_index][1] - *(result.data() + 3 * i + 1);
			record_vertex_position[vertex_index][2] += vertex_position[vertex_index][2] - *(result.data() + 3 * i + 2);
		}
		else {
			SUB(record_vertex_position[vertex_index], vertex_position[vertex_index], (result.data() + 3 * i));
		}
		record_vertex_num[vertex_index] ++;

		dis = abs(vertex_position[vertex_index][0] - *(result.data() + 3 * i) - collision_free_pos[vertex_index][0]);
		if (dis > max_dis) {
			max_dis = dis;
		}
		dis = abs(vertex_position[vertex_index][1] - *(result.data() + 3 * i + 1) - collision_free_pos[vertex_index][1]);
		if (dis > max_dis) {
			max_dis = dis;
		}
		dis = abs(vertex_position[vertex_index][2] - *(result.data() + 3 * i + 2) - collision_free_pos[vertex_index][2]);
		if (dis > max_dis) {
			max_dis = dis;
		}
	}

	//for (int i = 0; i < unfixed_vertex_num; ++i) {
	//	vertex_index = tet_actual_unfixed_vertex_indices[i];
	//	SUB_(vertex_position[vertex_index], (result.data() + 3 * i));
	//}
	////unsigned int* common_vertex_in_order_ = common_vertex_in_order;
	////if (result.norm() > 1e-1) {
	////	std::cout << result.norm() << std::endl;
	////	second_order_constraint.solveCD_ARAP_blockTest(Hessian, grad, vertex_position, stiffness, A[tet_index],
	////		volume[tet_index], tet_vertex_index,
	////		unfixed_tet_vertex_index, unfixed_vertex_num);
	////	for (auto i = neighbor_tet_indices.begin(); i < neighbor_tet_indices.end(); ++i) {
	////		std::cout << "i " << i- neighbor_tet_indices.begin() << std::endl;
	////		second_order_constraint.solveCertainHessianForNeighborTetTest(vertex_position, stiffness, A[*i], common_vertex_in_order_,
	////			tet_indices[obj_No][*i].data(),
	////			Hessian, volume[*i], grad);
	////	}
	////}


}



void XPBD_IPC::compareVector(std::vector<unsigned int>* a, std::vector<unsigned int>* b, int type)
{
	if (a->size() != b->size()) {
		std::cout << "type " << type << "not right "<<a->size()<<" "<<b->size() << std::endl;
	}
	else {
		for (int i = 0; i <a->size(); i += 4) {
			if (a->data()[i] != b->data()[i]) {
				std::cout << "error record pair " << type << " " << a->data()[i] << " " << b->data()[i] << std::endl;
			}
		}
	}
}

void XPBD_IPC::compareIfRecordHessianIsRight(int color)
{
	//compareVector(&collision.record_vt_pair[0], &collision_compare.record_vt_pair[0], 0);
	//compareVector(&collision.record_ee_pair[0], &collision_compare.record_ee_pair[0], 1);

	//compareVector(&collision.record_tv_collider_pair[0], &collision_compare.record_tv_collider_pair[0], 2);
	//compareVector(&collision.record_ee_collider_pair[0], &collision_compare.record_ee_collider_pair[0], 3);
	//compareVector(&collision.record_vt_collider_pair[0], &collision_compare.record_vt_collider_pair[0], 4);

	//double value_sum;
	//for (int i = 0; i < common_grad_compare.size(); ++i) {
	//	value_sum = 0.0;

	//	for (int k = 0; k < total_thread_num; ++k) {
	//		value_sum += common_grad[k][i];
	//	}
	//	if (abs(common_grad_compare[i] - value_sum) > 1e-8) {
	//		std::cout << "wrong in common grad "<<color<<" "<<i/3 << " " << abs(common_grad_compare[i] - common_grad[0][i]) << std::endl;
	//	}
	//}

	//auto k = common_hessian_compare.end();
	//for (auto i = common_hessian.begin(); i != common_hessian.end(); i++) {
	//	if (i->second.is_update) {
	//		k = common_hessian_compare.find(i->first);
	//		if (k == common_hessian_compare.end()) {
	//			std::cout << "error, cannot find pair in compare " << i->first[0] << " " << i->first[1] << std::endl;
	//			continue;
	//		}
	//		
	//		double value = 0.0;
	//		for (int j = 0; j < 9; ++j) {
	//			value += abs(i->second.hessian[j]-k->second[j]);
	//		}
	//		if (value > 1e-8) {
	//			std::cout << "inner iter " << inner_iteration_number <<"outer "<<outer_itr_num<<  std::endl;


	//			std::cout << "ori hessian show " << std::endl;
	//			for (int m = 0; m < 9; ++m) {
	//				std::cout << i->second.hessian[m] << " ";
	//			}
	//			std::cout << std::endl;
	//			std::cout << "compare hessian show " << std::endl;
	//			for (int m = 0; m < 9; ++m) {
	//				std::cout << k->second[m] << " ";
	//			}
	//			std::cout << std::endl;
	//			std::cout<<"error to sum all hessian in compare "<<color<<" " << value << " " << i->first[0] << " " << i->first[1] << std::endl;
	//		}
	//	}
	//}


	//auto k2 = common_hessian.end();
	//for (auto i = common_hessian_compare.begin(); i != common_hessian_compare.end(); i++) {
	//	if (i->first[0]<= i->first[1]) {
	//		k2 = common_hessian.find(i->first);
	//		if (k2 == common_hessian.end()) {
	//			std::cout << "error, cannot find pair  in ori" << i->first[0] << " " << i->first[1] << std::endl;
	//			continue;
	//		}
	//		if (!k2->second.is_update) {
	//			std::cout << "error, cannot find pair  in ori not update" << i->first[0] << " " << i->first[1] << std::endl;
	//			continue;
	//		}
	//		double value = 0.0;
	//		for (int j = 0; j < 9; ++j) {
	//			value += abs(i->second[j] - k2->second.hessian[j]);
	//		}
	//		if (value > 1e-8) {
	//			std::cout << "error to sum all hessian in ori " << value << " " << i->first[0] << " " << i->first[1] << std::endl;
	//		}
	//	}
	//}

	//
	//for (auto i = compare_floor_hessian.begin(); i != compare_floor_hessian.end(); ++i) {
	//	if (abs(i->second - floor_hessian[i->first]) > 1e-8) {
	//		std::cout << "floor error in ori " << i->first << " " << abs(i->second - floor_hessian[i->first]) << std::endl;
	//	}
	//}
	//auto k3 = compare_floor_hessian.end();
	//for (auto i =0; i < floor_hessian.size(); ++i) {
	//	if (floor_hessian[i] != 0.0) {
	//		k3 = compare_floor_hessian.find(i);
	//		if (k3 == compare_floor_hessian.end()) {
	//			std::cout << "floor error cannot find in compare " <<i << std::endl;
	//			continue;
	//		}
	//		if (abs(k3->second - floor_hessian[i]) > 1e-8) {
	//			std::cout << "floor error in compare " << abs(k3->second - floor_hessian[i]) << " " <<i<< std::endl;
	//		}
	//	}
	//}


}



void XPBD_IPC::solveTetBlock(std::array<double, 3>* vertex_position, double stiffness, double dt,
	double* mass, double* mass_inv,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices,
	double* volume, unsigned int tet_index, std::array<double, 3>* sn, unsigned int* common_vertex_in_order,
	int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num,
	double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	std::unordered_map<std::array<int, 2>, double, pair_hash>& tet_hessian, 
	std::unordered_map<std::array<unsigned int, 2>, StoreHessianWithOrderInConstraint, pair_hash>& collision_hessian,
	double* common_grad, std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet, int* vertex_index_on_surface, unsigned int prefix_sum_vetex_obj,
	double* floor_map, std::array<double, 3>* record_ori_pos, char* indicate_collide_with_floor,
	int color_No)//
{
	if (unfixed_vertex_num == 0) {
		return;
	}

	MatrixXd Hessian;
	VectorXd grad;
	grad.resize(3 * unfixed_vertex_num);
	Hessian.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);
	Hessian.setZero();
	grad.setZero();

	auto k = tet_hessian.end();
	for (int i = 0; i < unfixed_vertex_num; ++i) {
		for (int j = 0; j < unfixed_vertex_num; ++j) {
			k = tet_hessian.find(std::array{ tet_actual_unfixed_vertex_indices[i],tet_actual_unfixed_vertex_indices[j] });
			if(k== tet_hessian.end())	{
				std::cout << "error, cannot find the hessian of a tet solveTetBlock()" << std::endl;
			}
			Hessian.data()[3 * (j * Hessian.cols() + i)] = k->second;
			Hessian.data()[3 * (j * Hessian.cols() + i) + Hessian.cols() +1] = k->second;
			Hessian.data()[3 * (j * Hessian.cols() + i) + Hessian.cols()+ Hessian.cols() +2] = k->second;
		}
		memcpy(grad.data()+3*i, common_grad + 3 * (prefix_sum_vetex_obj + tet_actual_unfixed_vertex_indices[i]), 24);
	}

	if (perform_collision) {
		auto k_ = collision_hessian.end();
		double* add_of_hessian;
		double* hessian_in_collision;
		for (int i = 0; i < unfixed_vertex_num; ++i) {
			for (int j = 0; j < unfixed_vertex_num; ++j) {
				k_ = collision_hessian.find(std::array{ prefix_sum_vetex_obj + tet_actual_unfixed_vertex_indices[i], prefix_sum_vetex_obj + tet_actual_unfixed_vertex_indices[j]});
				if (k_ != collision_hessian.end()) {
					if (k_->second.is_update) {
						add_of_hessian = Hessian.data() + 3 * (j * Hessian.cols() + i);
						hessian_in_collision = k_->second.hessian.data();
						for (int m = 0; m < 3; ++m) {
							*add_of_hessian += *hessian_in_collision;
							*(add_of_hessian + 1) += *(hessian_in_collision + 1);
							*(add_of_hessian + 2) += *(hessian_in_collision + 2);
							hessian_in_collision += 3;
							add_of_hessian += Hessian.cols();
						}
						if (i != j) {
							add_of_hessian = Hessian.data() + 3 * (i * Hessian.cols() + j);
							hessian_in_collision = k_->second.hessian.data();
							for (int m = 0; m < 3; ++m) {
								*add_of_hessian += *hessian_in_collision;
								*(add_of_hessian + 1) += *(hessian_in_collision + 3);
								*(add_of_hessian + 2) += *(hessian_in_collision + 6);
								hessian_in_collision++;
								add_of_hessian += Hessian.cols();
							}
						}
					}					
				}
			}	
		}

		if (floor->exist) {
			for (int i = 0; i < unfixed_vertex_num; ++i) {
					Hessian.data()[(3 * i + floor->dimension) * (Hessian.cols() + 1)] += floor_map[prefix_sum_vetex_obj + tet_actual_unfixed_vertex_indices[i]];
			}
		}
	}

	double mass_dt_2;
	int vertex_index;
	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_actual_unfixed_vertex_indices[i];
		mass_dt_2 = mass[vertex_index] / (dt * dt);
		Hessian(3 * i, 3 * i) += mass_dt_2;
		Hessian(3 * i + 1, 3 * i + 1) += mass_dt_2;
		Hessian(3 * i + 2, 3 * i + 2) += mass_dt_2;
		grad.data()[3 * i] += mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
		grad.data()[3 * i + 1] += mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
		grad.data()[3 * i + 2] += mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);
	}

	LLT <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(grad);
	
	//double energy_initial = 1e-15;

	//if (perform_collision) {
	//	energy_initial = computeBlockCurrentEnergy(vertex_position, stiffness, dt, mass, A, neighbor_tet_indices, volume, tet_index, sn, tet_vertex_index, unfixed_vertex_num,
	//		collision_stiffness, obj_No, tet_actual_unfixed_vertex_indices, triangle_of_a_tet, edge_of_a_tet, vertex_index_on_surface, mass_inv, tet_indices[obj_No],
	//		collision.indicate_vertex_collide_with_floor[obj_No].data(), collision.record_vertex_collide_with_floor_d_hat[obj_No].data());
	//}

	for (int j = 0; j < result.size(); ++j) {
		if (abs(result[j]) > 0.3) {
			std::cout << result.transpose() << std::endl;
			std::cout << "move so large " << inner_iteration_number << " " << tet_index << " " << tet_actual_unfixed_vertex_indices[0] << " " << tet_actual_unfixed_vertex_indices[1] << " " << tet_actual_unfixed_vertex_indices[2] << " " << tet_actual_unfixed_vertex_indices[3] << std::endl;
			for (int i = 0; i < unfixed_vertex_num; ++i) {
				vertex_index = tet_actual_unfixed_vertex_indices[i];
				std::cout << vertex_position[vertex_index][0] - sn[vertex_index][0] << " " << vertex_position[vertex_index][1] - sn[vertex_index][1] << " " <<
					vertex_position[vertex_index][2] - sn[vertex_index][2] << std::endl;


			}
			break;
		}
	}

	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_actual_unfixed_vertex_indices[i];
		SUB_(vertex_position[vertex_index], (result.data() + 3 * i));
	}




	//if (perform_collision) {

		//double t = getCollisionTime(triangle_of_a_tet, edge_of_a_tet, obj_No, tet_actual_unfixed_vertex_indices,
		//	unfixed_vertex_num, vertex_index_on_surface, vertex_position, record_ori_pos, indicate_collide_with_floor);
		//double t1 = getInversionTime(tet_index, &neighbor_tet_indices, tet_indices[obj_No], vertex_position, record_ori_pos);
		//if (t1 < t) {
		//	t = t1;
		//}
		//if (t < 1.0) {
		//	for (int i = 0; i < unfixed_vertex_num; ++i) {
		//		vertex_index = tet_actual_unfixed_vertex_indices[i];
		//		COLLISION_POS(vertex_position[vertex_index], t, record_ori_pos[vertex_index], vertex_position[vertex_index]);
		//	}
		//}
		
		//if (energy_initial > energy_converge_standard) {
		//	double energy_current = computeBlockCurrentEnergy(vertex_position, stiffness, dt, mass, A, neighbor_tet_indices, volume, tet_index, sn, tet_vertex_index, unfixed_vertex_num,
		//		collision_stiffness, obj_No, tet_actual_unfixed_vertex_indices, triangle_of_a_tet, edge_of_a_tet, vertex_index_on_surface, mass_inv, tet_indices[obj_No],
		//		collision.indicate_vertex_collide_with_floor[obj_No].data(), collision.record_vertex_collide_with_floor_d_hat[obj_No].data());
		//	while (energy_current > energy_initial)
		//	{
		//		std::cout << "previous color " <<t<<" "<<  energy_current - energy_initial << std::endl;
		//		t *= 0.5;
		//		for (int i = 0; i < unfixed_vertex_num; ++i) {
		//			vertex_index = tet_actual_unfixed_vertex_indices[i];
		//			COLLISION_POS(vertex_position[vertex_index], 0.5, record_ori_pos[vertex_index], vertex_position[vertex_index]);
		//		}
		//		energy_current = computeBlockCurrentEnergy(vertex_position, stiffness, dt, mass, A, neighbor_tet_indices, volume, tet_index, sn, tet_vertex_index, unfixed_vertex_num,
		//			collision_stiffness, obj_No, tet_actual_unfixed_vertex_indices, triangle_of_a_tet, edge_of_a_tet, vertex_index_on_surface, mass_inv, tet_indices[obj_No],
		//			collision.indicate_vertex_collide_with_floor[obj_No].data(), collision.record_vertex_collide_with_floor_d_hat[obj_No].data());
		//		if (t < 1e-6) {
		//			std::cout << "===get out " << std::endl;
		//			break;
		//		}
		//	}
		//}
		
		//compute energy
		//for (int i = 0; i < unfixed_vertex_num; ++i) {
		//	memcpy(record_ori_pos[tet_actual_unfixed_vertex_indices[i]].data(), vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), 24);
		//}
	//}
	////if (result.norm() > 1e-1) {
	////	std::cout << result.norm() << std::endl;
	////	second_order_constraint.solveCD_ARAP_blockTest(Hessian, grad, vertex_position, stiffness, A[tet_index],
	////		volume[tet_index], tet_vertex_index,
	////		unfixed_tet_vertex_index, unfixed_vertex_num);
	////	for (auto i = neighbor_tet_indices.begin(); i < neighbor_tet_indices.end(); ++i) {
	////		std::cout << "i " << i- neighbor_tet_indices.begin() << std::endl;
	////		second_order_constraint.solveCertainHessianForNeighborTetTest(vertex_position, stiffness, A[*i], common_vertex_in_order_,
	////			tet_indices[obj_No][*i].data(),
	////			Hessian, volume[*i], grad);
	////	}
	////}
}




//PREVIOUS_COLOR_INVERSION_TEST
void XPBD_IPC::computePreviousColorInversion(int thread_No, int color_No)
{
	int obj_No, size;
	unsigned int start, end;
	unsigned int* tet_around_a_group;
	std::array<int, 4>* indices;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* ori_vertex_pos;
	double* mass_inv_;
	int* tet_vertex;
	double inversion_time = 1.0;
	double collision_time = 1.0;
	int color_group_index;
	unsigned int* tet_in_a_group;
	int j;

	collision.collision_time_thread[thread_No] = 1.0;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		obj_No = i + cloth->size();
		size = tetrahedron->data()[i].mesh_struct.indices.size();
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		vertex_pos = vertex_position[i + cloth->size()];
		ori_vertex_pos = record_vertex_position[i + cloth->size()].data();
		mass_inv_ = tetrahedron->data()[i].mesh_struct.mass_inv.data();
		color_group_index = inner_iteration_number % tet_color_groups[i]->size();

		tet_around_a_group = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_groups[color_group_index][color_No].data();
		start = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_groups_start_per_thread[color_group_index][color_No][thread_No];
		end = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_groups_start_per_thread[color_group_index][color_No][thread_No + 1];
		for (auto k = start; k < end; ++k) {
			j = tet_around_a_group[k];
			tet_vertex = indices[j].data();
			if (mass_inv_[tet_vertex[0]] != 0.0 || mass_inv_[tet_vertex[1]] != 0.0 || mass_inv_[tet_vertex[2]] != 0.0 || mass_inv_[tet_vertex[3]] != 0.0) {
				if (inversionTest::TetInversionTest(ori_vertex_pos[tet_vertex[0]].data(), ori_vertex_pos[tet_vertex[1]].data(),
					ori_vertex_pos[tet_vertex[2]].data(), ori_vertex_pos[tet_vertex[3]].data(), vertex_pos[tet_vertex[0]].data(),
					vertex_pos[tet_vertex[1]].data(), vertex_pos[tet_vertex[2]].data(), vertex_pos[tet_vertex[3]].data(), &inversion_time)) {
					if (inversion_time < collision_time) {
						collision_time = inversion_time;

						if (inversion_time == 0) {
							std::cout << "tet around inv test time zero" << j << " " << tet_vertex[0] << " " << tet_vertex[1] << " " << tet_vertex[2] << " " << tet_vertex[3] << std::endl;
						}

					}
				}
			}
		}

		end = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color_No][thread_No+1];
		start = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color_No][thread_No];
		tet_in_a_group = tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index][color_No].data();
		for (auto k = start; k < end; ++k) {
			j = tet_in_a_group[k];
			tet_vertex = indices[j].data();
			if (mass_inv_[tet_vertex[0]] != 0.0 || mass_inv_[tet_vertex[1]] != 0.0 || mass_inv_[tet_vertex[2]] != 0.0 || mass_inv_[tet_vertex[3]] != 0.0) {
				if (inversionTest::TetInversionTest(ori_vertex_pos[tet_vertex[0]].data(), ori_vertex_pos[tet_vertex[1]].data(),
					ori_vertex_pos[tet_vertex[2]].data(), ori_vertex_pos[tet_vertex[3]].data(), vertex_pos[tet_vertex[0]].data(),
					vertex_pos[tet_vertex[1]].data(), vertex_pos[tet_vertex[2]].data(), vertex_pos[tet_vertex[3]].data(), &inversion_time)) {
					if (inversion_time < collision_time) {
						collision_time = inversion_time;

						if (inversion_time == 0) {
							std::cout << "tet in inv test time zero" << j << " " << tet_vertex[0] << " " << tet_vertex[1] << " " << tet_vertex[2] << " " << tet_vertex[3] << std::endl;
						}

					}
				}
			}
		}
	}
	if (collision_time < 1.0) {
		collision_time *= 0.9;
	}

	if (collision_time < 0.0) {
		std::cout << "error for inversion test time previous color()" << std::endl;
	}

	if (collision_time == 0.0) {
		std::cout << "inversion test time previous color() equals 0" << std::endl;
	}

	if (collision_time < collision.collision_time_thread[thread_No]) {
		collision.collision_time_thread[thread_No] = collision_time;
	}
}

//FIRST_COLOR_ARAP_ENERGY
void XPBD_IPC::computePreviousColorARAPEnergy(int thread_No, unsigned int color_No)
{
	unsigned int start, end;
	unsigned int* tet_around_a_group;

	MatrixXd grad;
	grad.resize(3, 4);
	std::array<int, 4>* indices;

	std::array<double, 3>* vertex_pos;
	double stiffness = 0.0;
	Matrix<double, 3, 4>* A;
	double* volume;

	unsigned int prefix_sum_start;

	int j = 0;

	int color_group_index;
	unsigned int* tet_in_a_group;
	double* mass_inv_;

	double energy = 0.0;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		color_group_index = inner_iteration_number % tet_color_groups[i]->size();
		if (color_No >= tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index].size() - 1) {
			continue;
		}

		stiffness = tetrahedron->data()[0].ARAP_stiffness;
		vertex_pos = vertex_position[i + cloth->size()];
		A = tetrahedron->data()[i].mesh_struct.A.data();

		volume = tet_volume[i + cloth->size()];
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		mass_inv_ = tetrahedron->data()[i].mesh_struct.mass_inv.data();

		tet_around_a_group = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_groups[color_group_index][color_No].data();
		start = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_groups_start_per_thread[color_group_index][color_No][thread_No];
		end = tetrahedron->data()[i].mesh_struct.tet_around_tet_color_groups_start_per_thread[color_group_index][color_No][thread_No + 1];


		for (auto k = start; k < end; ++k) {
			j = tet_around_a_group[k];
			if (mass_inv_[indices[j][0]] != 0.0 || mass_inv_[indices[j][1]] != 0.0 || mass_inv_[indices[j][2]] != 0.0 || mass_inv_[indices[j][3]] != 0.0) {
				energy += compute_energy.ARAPEnergy(vertex_pos[indices[j][0]].data(), vertex_pos[indices[j][1]].data(),
					vertex_pos[indices[j][2]].data(), vertex_pos[indices[j][3]].data(), A[j], volume[j], stiffness);
			}
		}

		end = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color_No][thread_No + 1];
		start = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color_No][thread_No];
		tet_in_a_group = tetrahedron->data()[i].mesh_struct.tet_color_group[color_group_index][color_No].data();

		for (auto k = start; k < end; ++k) {
			j = tet_in_a_group[k];
			if (mass_inv_[indices[j][0]] != 0.0 || mass_inv_[indices[j][1]] != 0.0 || mass_inv_[indices[j][2]] != 0.0 || mass_inv_[indices[j][3]] != 0.0) {
				energy += compute_energy.ARAPEnergy(vertex_pos[indices[j][0]].data(), vertex_pos[indices[j][1]].data(),
					vertex_pos[indices[j][2]].data(), vertex_pos[indices[j][3]].data(), A[j], volume[j], stiffness);
			}
		}
	}

	energy_per_thread[thread_No] = 0.5*energy;
}



//LAST_COLOR_ARAP_ENERGY
void XPBD_IPC::computeLastColorARAPEnergy(int thread_No)
{
	int obj_No, start,end;
	double energy = 0.0;

	std::array<int, 4>* indices;
	std::array<double, 3>* vertex_pos;
	Matrix<double, 3, 4>* A;
	double* volume;
	double stiffness;
	double* mass_inv_;
	
	for (int i = 0; i < tetrahedron->size(); ++i) {
		obj_No = i + cloth->size();
		start = tet_index_begin_per_thread[obj_No][thread_No];
		end = tet_index_begin_per_thread[obj_No][thread_No+1];
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		vertex_pos = vertex_position[i + cloth->size()];
		A = tetrahedron->data()[i].mesh_struct.A.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		stiffness = tetrahedron->data()[i].ARAP_stiffness;
		mass_inv_ = tetrahedron->data()[i].mesh_struct.mass_inv.data();

		for (int j = start; j < end; ++j) {
			if (is_tet_arap_grad_compute[prefix_sum_of_every_tet_index[i] + j].test(std::memory_order_relaxed)) {
				if (mass_inv_[indices[j][0]] != 0.0 || mass_inv_[indices[j][1]] != 0.0 || mass_inv_[indices[j][2]] != 0.0 || mass_inv_[indices[j][3]] != 0.0) {
					energy += compute_energy.ARAPEnergy(vertex_pos[indices[j][0]].data(), vertex_pos[indices[j][1]].data(),
						vertex_pos[indices[j][2]].data(), vertex_pos[indices[j][3]].data(), A[j], volume[j], stiffness);
				}
			}
		}
	}

	energy_per_thread[thread_No] = 0.5 * energy;
}

double XPBD_IPC::computeLastColorARAPEnergy()
{
	double energy = 0.0;
	thread->assignTask(this, LAST_COLOR_ARAP_ENERGY);
	for (int i = 0; i < total_thread_num; ++i) {
		energy += energy_per_thread[i];
	}
	return energy;
}


double XPBD_IPC::computeFloorEnergy(int type,double collision_stiffness, std::vector<unsigned int>* record_vertex_collide_with_floor, int start, int end)
{
	auto start_ = record_vertex_collide_with_floor->begin() + start;
	auto end_ = record_vertex_collide_with_floor->begin() + end;

	double energy = 0.0;
	if (type == 0) {
		for (auto i = start_; i < end_; i += 2) {
			energy += compute_energy.computeFloorBarrierEnergy(vertex_position[*i][*(i + 1)][floor->dimension], collision.record_vertex_collide_with_floor_d_hat[*i][*(i + 1)],
				collision_stiffness, floor->value);
		}
	}
	else {
		for (auto i = start_; i < end_; i += 2) {
			if (collision.vertex_belong_to_color_group[*i][*(i + 1)]) {
				energy += compute_energy.computeFloorBarrierEnergy(vertex_position[*i][*(i + 1)][floor->dimension], collision.record_vertex_collide_with_floor_d_hat[*i][*(i + 1)],
					collision_stiffness, floor->value);
			}
		}
	}

	return energy;

}

//0 total, 1 inovled in this color
double XPBD_IPC::computeFloorEnergy(int type)
{
	double stiffness = 0.0;
	if (!tetrahedron->empty()) {
		stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		stiffness = cloth->data()[0].collision_stiffness[0];
	}
	double energy = 0.0;
	if (type == 0) {
		for (auto j = collision.record_vertex_collide_with_floor.begin(); j < collision.record_vertex_collide_with_floor.end(); ++j) {
			for (auto i = j->begin(); i < j->end(); i += 2) {
				energy += compute_energy.computeFloorBarrierEnergy(vertex_position[*i][*(i + 1)][floor->dimension], collision.record_vertex_collide_with_floor_d_hat[*i][*(i + 1)],
					stiffness, floor->value);
			}			
		}
	}
	else {
		for (auto j = collision.record_vertex_collide_with_floor.begin(); j < collision.record_vertex_collide_with_floor.end(); ++j) {
			for (auto i = j->begin(); i < j->end(); i += 2) {
				if (collision.vertex_belong_to_color_group[*i][*(i + 1)]) {
					energy += compute_energy.computeFloorBarrierEnergy(vertex_position[*i][*(i + 1)][floor->dimension], collision.record_vertex_collide_with_floor_d_hat[*i][*(i + 1)],
						stiffness, floor->value);
				}
			}
		}
	}
	return energy;
}





double XPBD_IPC::computePreviousColorEnergy(unsigned int color_no)
{
	double energy = 0.0;
	energy += computePreviousColorInertialEnergy(color_no);

	thread->assignTask(this, FIRST_COLOR_ARAP_ENERGY, color_no);
	for (int i = 0; i < total_thread_num; ++i) {
		energy += energy_per_thread[i];
	}

	if (perform_collision) {
		energy += computePreviousColorCollisionEnergy();
	}

	return energy;
}



double XPBD_IPC::computeLastColorEnergy()
{
	double energy = 0.0;
	energy += computeColorInertialEnergy();
	//ARAP
	energy += computeLastColorARAPEnergy();
	if (perform_collision) {
		energy += computeBarrierEnergy();
	}

	return energy;
}


//COMPUTE_BARRIER_ENERGY
void XPBD_IPC::computeBarrierEnergy(int thread_No)
{
	double energy = 0.0;
	double stiffness = 0.0;
	if (!tetrahedron->empty()) {
		stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		stiffness = cloth->data()[0].collision_stiffness[0];
	}
	//vt
	energy += computeVTEnergy(&collision.record_vt_pair_sum_all_thread, triangle_indices.data(), vertex_position.data(), vertex_position.data(),
		stiffness, collision.record_vt_pair_d_hat.data(),
		collision.vertex_belong_to_color_group, 0, collision.vt_pair_sum_start_per_thread[thread_No],
		collision.vt_pair_sum_start_per_thread[thread_No + 1]);
	//ee
	energy += computeEEEnergy(&collision.record_ee_pair_sum_all_thread, edge_vertices.data(), edge_vertices.data(),
		vertex_position.data(), vertex_position.data(), stiffness, collision.record_ee_pair_d_hat.data(),
		collision.vertex_belong_to_color_group, 0, collision.ee_pair_sum_start_per_thread[thread_No],
		collision.ee_pair_sum_start_per_thread[thread_No + 1]);

	if (has_collider) {
		//vt c
		energy += computeVTEnergy(&collision.record_vt_collider_pair_sum_all_thread, triangle_indices_collider.data(), vertex_position.data(), vertex_position_collider.data(),
			stiffness, collision.record_vt_collider_pair_d_hat.data(),
			collision.vertex_belong_to_color_group, 0, collision.vt_collider_pair_sum_start_per_thread[thread_No],
			collision.vt_collider_pair_sum_start_per_thread[thread_No + 1]);

		//tv c
		energy += computeVTEnergy(&collision.record_tv_collider_pair_sum_all_thread, triangle_indices.data(), vertex_position_collider.data(), vertex_position.data(),
			stiffness, collision.record_tv_collider_pair_d_hat.data(),
			collision.vertex_belong_to_color_group, 0, collision.tv_collider_pair_sum_start_per_thread[thread_No],
			collision.tv_collider_pair_sum_start_per_thread[thread_No + 1]);
		//ee c
		energy += computeEEEnergy(&collision.record_ee_collider_pair_sum_all_thread, edge_vertices.data(), collider_edge_vertices.data(),
			vertex_position.data(), vertex_position_collider.data(), stiffness, collision.record_ee_collider_pair_d_hat.data(),
			collision.vertex_belong_to_color_group, 0, collision.ee_collider_pair_sum_start_per_thread[thread_No],
			collision.ee_collider_pair_sum_start_per_thread[thread_No + 1]);
	}
	if (floor->exist) {
		energy += computeFloorEnergy(0, stiffness, &collision.record_vertex_collide_with_floor_sum_all_thread, collision.floor_pair_sum_start_per_thread[thread_No],
			collision.floor_pair_sum_start_per_thread[thread_No + 1]);
	}
	energy_per_thread[thread_No] = energy;
}




//type 0: compute all pair 1: ee in this color 2: ee_c in this color
double XPBD_IPC::computeEEEnergy(std::vector<unsigned int>* record_pair, unsigned int** edge_v_0, unsigned int** edge_v_1,
	std::array<double, 3>** e0_current_pos, std::array<double, 3>** e1_current_pos, double collision_stiffness, double* d_hat,
	bool** belong_to_color_group, int type, unsigned int start, unsigned int end)
{
	double energy = 0.0;
	unsigned int* edge_0_vertex;
	unsigned int* edge_1_vertex;

	auto pair_end = record_pair->begin() + end;
	auto pair_start = record_pair->begin() + start;
	double* d_hat_ = d_hat + (start >> 2);

	switch (type)
	{
	case 0:
	{
		for (auto i = pair_start; i < pair_end; i += 4) {
			edge_0_vertex = edge_v_0[*i] + ((*(i + 1)) << 1);
			edge_1_vertex = edge_v_1[*(i + 2)] + ((*(i + 3)) << 1);
			energy += compute_energy.computeBarrierEnergy(e0_current_pos[*i][*edge_0_vertex].data(),
				e0_current_pos[*i][*(edge_0_vertex + 1)].data(), e1_current_pos[*(i + 2)][*edge_1_vertex].data(),
				e1_current_pos[*(i + 2)][*(edge_1_vertex + 1)].data(), collision_stiffness, *d_hat_, false);
			d_hat_++;
		}
	}
	break;

	case 1:
	{
		for (auto i = pair_start; i < pair_end; i += 4) {
			edge_0_vertex = edge_v_0[*i] + ((*(i + 1)) << 1);
			edge_1_vertex = edge_v_1[*(i + 2)] + ((*(i + 3)) << 1);
			if (belong_to_color_group[*i][edge_0_vertex[0]] || belong_to_color_group[*i][edge_0_vertex[1]]
				|| belong_to_color_group[*(i + 2)][edge_1_vertex[0]] || belong_to_color_group[*(i + 2)][edge_1_vertex[1]]) {
				energy += compute_energy.computeBarrierEnergy(e0_current_pos[*i][*edge_0_vertex].data(),
					e0_current_pos[*i][*(edge_0_vertex + 1)].data(), e1_current_pos[*(i + 2)][*edge_1_vertex].data(),
					e1_current_pos[*(i + 2)][*(edge_1_vertex + 1)].data(), collision_stiffness, *d_hat_, false);
			}
			d_hat_++;
		}
	}
	break;

	case 2:
	{
		for (auto i = pair_start; i < pair_end; i += 4) {
			edge_0_vertex = edge_v_0[*i] + ((*(i + 1)) << 1);
			if (belong_to_color_group[*i][edge_0_vertex[0]] || belong_to_color_group[*i][edge_0_vertex[1]]) {
				edge_1_vertex = edge_v_1[*(i + 2)] + ((*(i + 3)) << 1);
				energy += compute_energy.computeBarrierEnergy(e0_current_pos[*i][*edge_0_vertex].data(),
					e0_current_pos[*i][*(edge_0_vertex + 1)].data(), e1_current_pos[*(i + 2)][*edge_1_vertex].data(),
					e1_current_pos[*(i + 2)][*(edge_1_vertex + 1)].data(), collision_stiffness, *d_hat_, false);

			}
			d_hat_++;
		}
	}
	break;
	}

	return energy;
}

//type 0: compute all pair 1: ee in this color 2: ee_c in this color
double XPBD_IPC::computeEEEnergy(std::vector<std::vector<unsigned int>>* record_pair, unsigned int** edge_v_0, unsigned int** edge_v_1,
	 std::array<double, 3>** e0_current_pos, std::array<double, 3>** e1_current_pos, double collision_stiffness, std::vector<std::vector<double>>* d_hat, bool** belong_to_color_group,
	int type)
{
	double energy = 0.0;
	unsigned int* edge_0_vertex;
	unsigned int* edge_1_vertex;
	double* d_hat_;
	//if (type == 0) {

	switch (type)
	{
	case 0:
	{
		for (auto j = record_pair->begin(); j < record_pair->end(); ++j) {
			d_hat_ = d_hat->data()[j - record_pair->begin()].data();
			for (auto i = j->begin(); i < j->end(); i += 4) {
				edge_0_vertex = edge_v_0[*i] + ((*(i + 1)) << 1);
				edge_1_vertex = edge_v_1[*(i + 2)] + ((*(i + 3)) << 1);
				energy += compute_energy.computeBarrierEnergy(e0_current_pos[*i][*edge_0_vertex].data(),
					e0_current_pos[*i][*(edge_0_vertex + 1)].data(), e1_current_pos[*(i + 2)][*edge_1_vertex].data(),
					e1_current_pos[*(i + 2)][*(edge_1_vertex + 1)].data(), collision_stiffness, d_hat_[(i - j->begin()) >> 2], false);
			}
		}
	}
		break;

	case 1:
	{
		for (auto j = record_pair->begin(); j < record_pair->end(); ++j) {
			d_hat_ = d_hat->data()[j - record_pair->begin()].data();
			for (auto i = j->begin(); i < j->end(); i += 4) {
				edge_0_vertex = edge_v_0[*i] + ((*(i + 1)) << 1);				
				edge_1_vertex = edge_v_1[*(i + 2)] + ((*(i + 3)) << 1);
				if (belong_to_color_group[*i][edge_0_vertex[0]] || belong_to_color_group[*i][edge_0_vertex[1]]
					|| belong_to_color_group[*(i + 2)][edge_1_vertex[0]] || belong_to_color_group[*(i + 2)][edge_1_vertex[1]]) {
					energy += compute_energy.computeBarrierEnergy(e0_current_pos[*i][*edge_0_vertex].data(),
						e0_current_pos[*i][*(edge_0_vertex + 1)].data(), e1_current_pos[*(i + 2)][*edge_1_vertex].data(),
						e1_current_pos[*(i + 2)][*(edge_1_vertex + 1)].data(), collision_stiffness, d_hat_[(i - j->begin()) >> 2], false);
				}
			}
		}
	}
	break;

	case 2:
	{
		for (auto j = record_pair->begin(); j < record_pair->end(); ++j) {
			d_hat_ = d_hat->data()[j - record_pair->begin()].data();
			for (auto i = j->begin(); i < j->end(); i += 4) {
				edge_0_vertex = edge_v_0[*i] + ((*(i + 1)) << 1);
				if (belong_to_color_group[*i][edge_0_vertex[0]] || belong_to_color_group[*i][edge_0_vertex[1]]) {
					edge_1_vertex = edge_v_1[*(i + 2)] + ((*(i + 3)) << 1);
					energy += compute_energy.computeBarrierEnergy(e0_current_pos[*i][*edge_0_vertex].data(),
						e0_current_pos[*i][*(edge_0_vertex + 1)].data(), e1_current_pos[*(i + 2)][*edge_1_vertex].data(),
						e1_current_pos[*(i + 2)][*(edge_1_vertex + 1)].data(), collision_stiffness, d_hat_[(i - j->begin()) >> 2], false);
				}
			}
		}
	}
	break;
	}

	return energy;
}

//type 0: compute all pair 1: vt in this color 2: vt_c in this color 3: tv_c in this color
double XPBD_IPC::computeVTEnergy(std::vector<unsigned int>* record_vt_pair, std::array<int, 3>** triangle_indices,
	std::array<double, 3>** v_current_pos, std::array<double, 3>** t_current_pos, double collision_stiffness, double* d_hat,
	bool** belong_to_color_group, int type, unsigned int start, unsigned int end)
{
	double energy = 0.0;
	int* indices;
	auto pair_end = record_vt_pair->begin() + end;
	auto pair_start = record_vt_pair->begin() + start;
	double* d_hat_ = d_hat + (start >> 2);

	switch (type)
	{
	case 0:
	{
		for (auto i = pair_start; i < pair_end; i+=4) {
			indices = triangle_indices[*(i + 2)][*(i + 3)].data();
			energy += compute_energy.computeBarrierEnergy(v_current_pos[*i][*(i + 1)].data(),
				t_current_pos[*(i + 2)][indices[0]].data(),
				t_current_pos[*(i + 2)][indices[1]].data(),
				t_current_pos[*(i + 2)][indices[2]].data(), collision_stiffness, *d_hat_, true);
			d_hat_++;
		}
	}
	break;

	case 1:
	{
		for (auto i = pair_start; i < pair_end; i += 4) {
			indices = triangle_indices[*(i + 2)][*(i + 3)].data();
			if (belong_to_color_group[*i][*(i + 1)] || belong_to_color_group[*(i + 2)][indices[0]] ||
				belong_to_color_group[*(i + 2)][indices[1]] || belong_to_color_group[*(i + 2)][indices[2]]) {
				energy += compute_energy.computeBarrierEnergy(v_current_pos[*i][*(i + 1)].data(),
					t_current_pos[*(i + 2)][indices[0]].data(),
					t_current_pos[*(i + 2)][indices[1]].data(),
					t_current_pos[*(i + 2)][indices[2]].data(), collision_stiffness, *d_hat_, true);
			}
			d_hat_++;
		}
	}
	break;

	case 2:
	{
		for (auto i = pair_start; i < pair_end; i += 4) {
			if (belong_to_color_group[*i][*(i + 1)]) {
				indices = triangle_indices[*(i + 2)][*(i + 3)].data();
				energy += compute_energy.computeBarrierEnergy(v_current_pos[*i][*(i + 1)].data(),
					t_current_pos[*(i + 2)][indices[0]].data(),
					t_current_pos[*(i + 2)][indices[1]].data(),
					t_current_pos[*(i + 2)][indices[2]].data(), collision_stiffness, *d_hat_, true);

			}
			d_hat_++;
		}
	}
	break;

	case 3:
		for (auto i = pair_start; i < pair_end; i += 4) {
			indices = triangle_indices[*(i + 2)][*(i + 3)].data();
			if (belong_to_color_group[*(i + 2)][indices[0]] || belong_to_color_group[*(i + 2)][indices[1]] || belong_to_color_group[*(i + 2)][indices[2]]) {
				energy += compute_energy.computeBarrierEnergy(v_current_pos[*i][*(i + 1)].data(),
					t_current_pos[*(i + 2)][indices[0]].data(),
					t_current_pos[*(i + 2)][indices[1]].data(),
					t_current_pos[*(i + 2)][indices[2]].data(), collision_stiffness, *d_hat_, true);
			}
			d_hat_++;
		}
		break;
	}

	return energy;
}


//type 0: compute all pair 1: vt in this color 2: vt_c in this color 3: tv_c in this color
double XPBD_IPC::computeVTEnergy(std::vector<std::vector<unsigned int>>* record_vt_pair, std::array<int, 3>** triangle_indices,
	std::array<double, 3>** v_current_pos,	std::array<double, 3>** t_current_pos, double collision_stiffness, std::vector<std::vector<double>>* d_hat, 
	bool** belong_to_color_group, int type)
{
	double energy = 0.0;
	int* indices;
	double* d_hat_;
	switch (type)
	{
	case 0:
	{
		for (auto j = record_vt_pair->begin(); j < record_vt_pair->end(); ++j) {
			d_hat_ = d_hat->data()[j - record_vt_pair->begin()].data();
			for (auto i = j->begin(); i < j->end(); i += 4) {
				indices = triangle_indices[*(i + 2)][*(i + 3)].data();
				energy += compute_energy.computeBarrierEnergy(v_current_pos[*i][*(i + 1)].data(),
					t_current_pos[*(i + 2)][indices[0]].data(),
					t_current_pos[*(i + 2)][indices[1]].data(),
					t_current_pos[*(i + 2)][indices[2]].data(), collision_stiffness, d_hat_[(i - j->begin()) >> 2], true);
			}
		}
	}
		break;

	case 1:
	{
		for (auto j = record_vt_pair->begin(); j < record_vt_pair->end(); ++j) {
			d_hat_ = d_hat->data()[j - record_vt_pair->begin()].data();
			for (auto i = j->begin(); i < j->end(); i += 4) {
				indices = triangle_indices[*(i + 2)][*(i + 3)].data();
				if (belong_to_color_group[*i][*(i + 1)] || belong_to_color_group[*(i + 2)][indices[0]] || 
					belong_to_color_group[*(i + 2)][indices[1]] || belong_to_color_group[*(i + 2)][indices[2]]) {
					energy += compute_energy.computeBarrierEnergy(v_current_pos[*i][*(i + 1)].data(),
						t_current_pos[*(i + 2)][indices[0]].data(),
						t_current_pos[*(i + 2)][indices[1]].data(),
						t_current_pos[*(i + 2)][indices[2]].data(), collision_stiffness, d_hat_[(i - j->begin()) >> 2], true);
				}				
			}
		}
	}
	break;

	case 2:
	{
		for (auto j = record_vt_pair->begin(); j < record_vt_pair->end(); ++j) {
			d_hat_ = d_hat->data()[j - record_vt_pair->begin()].data();
			for (auto i = j->begin(); i < j->end(); i += 4) {
				if (belong_to_color_group[*i][*(i + 1)]) {
					indices = triangle_indices[*(i + 2)][*(i + 3)].data();
					energy += compute_energy.computeBarrierEnergy(v_current_pos[*i][*(i + 1)].data(),
						t_current_pos[*(i + 2)][indices[0]].data(),
						t_current_pos[*(i + 2)][indices[1]].data(),
						t_current_pos[*(i + 2)][indices[2]].data(), collision_stiffness, d_hat_[(i - j->begin()) >> 2], true);

				}
			}
		}
	}
	break;

	case 3:
		for (auto j = record_vt_pair->begin(); j < record_vt_pair->end(); ++j) {
			d_hat_ = d_hat->data()[j - record_vt_pair->begin()].data();
			for (auto i = j->begin(); i < j->end(); i += 4) {
				indices = triangle_indices[*(i + 2)][*(i + 3)].data();
				if (belong_to_color_group[*(i + 2)][indices[0]] || belong_to_color_group[*(i + 2)][indices[1]] || belong_to_color_group[*(i + 2)][indices[2]]) {
					energy += compute_energy.computeBarrierEnergy(v_current_pos[*i][*(i + 1)].data(),
						t_current_pos[*(i + 2)][indices[0]].data(),
						t_current_pos[*(i + 2)][indices[1]].data(),
						t_current_pos[*(i + 2)][indices[2]].data(), collision_stiffness, d_hat_[(i - j->begin()) >> 2], true);
				}
			}
		}
	break;
	}

	return energy;
}


double XPBD_IPC::computePreviousColorInertialEnergy(unsigned int color_No)
{
	double energy = 0.0;
	thread->assignTask(this, PREVIOUS_COLOR_INERTIAL_ENERGY, color_No);
	for (int i = 0; i < total_thread_num; ++i) {
		energy += energy_per_thread[i];
	}
	return energy;
}


//MAX_DISPLACEMENT_COLOR
void XPBD_IPC::maxDisplacement(int thread_No, unsigned int color_No)
{
	int start, end;
	bool* vertex_belong;
	double max_dis=0.0;
	double dis;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* record_vertex_pos;

	for (int i = 0; i < total_obj_num; ++i) {
		start = vertex_index_begin_per_thread[i][thread_No];
		end = vertex_index_begin_per_thread[i][thread_No + 1];
		vertex_belong = collision.vertex_belong_to_color_group[i];
		vertex_pos = vertex_position[i];
		record_vertex_pos = record_collision_free_vertex_position[i].data();
		for (int j = start; j < end; ++j) {
			if (vertex_belong[j]) {
				dis = abs(vertex_pos[j][0] - record_vertex_pos[j][0]);
				if (dis > max_dis) {
					max_dis = dis;
				}
				dis = abs(vertex_pos[j][1] - record_vertex_pos[j][1]);
				if (dis > max_dis) {
					max_dis = dis;
				}
				dis = abs(vertex_pos[j][2] - record_vertex_pos[j][2]);
				if (dis > max_dis) {
					max_dis = dis;
				}
			}
		}
	}
	if (max_dis_record[thread_No] < max_dis) {
		max_dis_record[thread_No] = max_dis;
	}
}

//PREVIOUS_COLOR_INERTIAL_ENERGY
void XPBD_IPC::computePreviousColorInertialEnergy(int thread_No, unsigned int color_No)
{
	double energy = 0.0;
	std::vector<unsigned int>* index_group;
	int color_group_index;
	int i;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* sn_;
	double* mass;	

	for (int tet_No = 0; tet_No < tetrahedron->size(); ++tet_No) {
		i = tet_No + cloth->size();
		vertex_pos = vertex_position[i];
		sn_ = sn[i].data();
		mass = mesh_struct[i]->mass.data();
		color_group_index = inner_iteration_number % tet_color_groups[tet_No]->size();
		index_group = &vertex_index_of_a_tet_color_group[i][color_group_index][color_No];

		auto begin = index_group->begin() + vertex_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color_No][thread_No];
		auto end = index_group->begin() + vertex_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color_No][thread_No + 1];

		for (auto j = begin; j < end; ++j) {
			energy += mass[*j] * (EDGE_LENGTH(vertex_pos[*j], sn_[*j]));
		}
	}

	energy_per_thread[thread_No]= 0.5 * energy / (sub_time_step * sub_time_step);
}


double XPBD_IPC::computeColorInertialEnergy()
{
	double energy = 0.0;
	thread->assignTask(this, LAST_COLOR_INERTIAL_ENERGY);
	for (int i = 0; i < total_thread_num; ++i) {
		energy += energy_per_thread[i];
	}
	return energy;
}




//LAST_COLOR_INERTIAL_ENERGY
void  XPBD_IPC::computeColorInertialEnergy(int thread_No)
{
	double energy = 0.0;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* sn_;
	unsigned int vertex_start;
	unsigned int vertex_end;
	double* mass;
	bool* belong_color_group;
	bool* belong_;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		sn_ = sn[i].data();
		vertex_end = vertex_index_begin_per_thread[i][thread_No+1];
		vertex_start = vertex_index_begin_per_thread[i][thread_No];
		mass = mesh_struct[i]->mass.data();
		belong_color_group = collision.vertex_belong_to_color_group[i];
		for (unsigned int j = vertex_start; j < vertex_end; ++j) {
			if (belong_color_group[j]) {
				energy += mass[j] * (EDGE_LENGTH(vertex_pos[j], sn_[j]));
			}
		}
	}

	energy_per_thread[thread_No] = 0.5 * energy / (sub_time_step * sub_time_step);
}



double XPBD_IPC::computeBlockCurrentEnergy(std::array<double, 3>* vertex_position, double stiffness, double dt,
	double* mass,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices,
	double* volume, unsigned int tet_index, std::array<double, 3>* sn, 
	int* tet_vertex_index, unsigned int unfixed_vertex_num,
	double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet, int* vertex_index_on_surface, 
	double* mass_inv_, std::array<int,4>*tet_vertex_indices,
	char* indicate_vertex_collide_with_floor, double* record_vertex_collide_with_floor_d_hat)
{
	double energy = 0.0;
	//intertial
	int vertex_index;
	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_actual_unfixed_vertex_indices[i];
		energy += mass[vertex_index] * (EDGE_LENGTH(vertex_position[vertex_index], sn[vertex_index]));
	}
	energy = 0.5 * energy / (dt * dt);

	//ARAP:
	for (auto i = neighbor_tet_indices.begin(); i < neighbor_tet_indices.end(); ++i) {
		if (mass_inv_[tet_vertex_indices[*i][0]] != 0.0 || mass_inv_[tet_vertex_indices[*i][1]] != 0.0 || mass_inv_[tet_vertex_indices[*i][2]] != 0.0 || mass_inv_[tet_vertex_indices[*i][3]] != 0.0) {
			energy += 0.5*compute_energy.ARAPEnergy(vertex_position[tet_vertex_indices[*i][0]].data(), vertex_position[tet_vertex_indices[*i][1]].data(),
				vertex_position[tet_vertex_indices[*i][2]].data(), vertex_position[tet_vertex_indices[*i][3]].data(), A[*i], volume[*i], stiffness);
		}
	}
	energy += 0.5*compute_energy.ARAPEnergy(vertex_position[tet_vertex_index[0]].data(), vertex_position[tet_vertex_index[1]].data(),
		vertex_position[tet_vertex_index[2]].data(), vertex_position[tet_vertex_index[3]].data(), A[tet_index], volume[tet_index], stiffness);

	//BARRIER
	//VT
	for (int i = 0; i < unfixed_vertex_num; ++i) {
		if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
			continue;
		}
		energy += computeVTCollisionEnergyPerElement(vertex_position[tet_actual_unfixed_vertex_indices[i]].data(),
			collision.vertex_triangle_pair_by_vertex[obj_No] + collision.close_vt_pair_num * tet_actual_unfixed_vertex_indices[i],
			this->vertex_position.data(), triangle_indices.data(), collision.vertex_triangle_pair_num_record[obj_No][tet_actual_unfixed_vertex_indices[i]],
			collision.vertex_triangle_pair_d_hat[obj_No] + collision.vt_d_hat_num * tet_actual_unfixed_vertex_indices[i], collision_stiffness);
	}

	//TV
	int* triangle_;
	for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
		triangle_ = triangle_indices[obj_No][*i].data();
		energy += computeTVCollisionEnergyPerElement(vertex_position[triangle_[0]].data(), vertex_position[triangle_[1]].data(),
			vertex_position[triangle_[2]].data(), collision.triangle_vertex_pair_by_triangle[obj_No] + collision.close_tv_pair_num * (*i),
			collision.triangle_vertex_pair_num_record[obj_No][*i], collision.triangle_vertex_pair_d_hat[obj_No] + collision.tv_d_hat_num * (*i), collision_stiffness,
			tet_actual_unfixed_vertex_indices, unfixed_vertex_num, this->vertex_position.data(), obj_No);
	}

	//EE
	unsigned int* edge_;
	for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
		edge_ = edge_vertices[obj_No] + ((*i) << 1);
		energy += computeEECollisionEnergyPerElement(vertex_position[*edge_].data(),
			vertex_position[*(edge_ + 1)].data(), collision.edge_edge_pair_by_edge[obj_No] + collision.close_ee_pair_num * (*i),
			collision.edge_edge_pair_number_record[obj_No][*i], collision.edge_edge_pair_d_hat[obj_No] + collision.ee_d_hat_num * (*i),
			collision_stiffness, i - edge_of_a_tet->begin(), edge_of_a_tet, this->vertex_position.data(), obj_No, edge_vertices.data());
	}

	if (has_collider) {
		//vt collider
		for (int i = 0; i < unfixed_vertex_num; ++i) {
			if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
				continue;
			}
			energy += computeVTCollisionEnergyPerElement(vertex_position[tet_actual_unfixed_vertex_indices[i]].data(),
				collision.vertex_obj_triangle_collider_pair_by_vertex[obj_No] + collision.close_vt_collider_pair_num * tet_actual_unfixed_vertex_indices[i],
				vertex_position_collider.data(), triangle_indices_collider.data(), collision.vertex_obj_triangle_collider_num_record[obj_No][tet_actual_unfixed_vertex_indices[i]],
				collision.vertex_obj_triangle_collider_pair_d_hat[obj_No] + collision.vt_collider_d_hat_num * tet_actual_unfixed_vertex_indices[i], collision_stiffness);
		}
		//TV collider
		for (auto i = triangle_of_a_tet->begin(); i < triangle_of_a_tet->end(); ++i) {
			triangle_ = triangle_indices[obj_No][*i].data();
			energy += computeTVColliderCollisionEnergyPerElement(vertex_position[triangle_[0]].data(), vertex_position[triangle_[1]].data(),
				vertex_position[triangle_[2]].data(), collision.triangle_vertex_collider_pair_by_triangle[obj_No] + collision.close_tv_collider_pair_num * (*i),
				collision.triangle_vertex_collider_pair_num_record[obj_No][*i], collision.triangle_vertex_collider_pair_d_hat[obj_No] + collision.tv_collider_d_hat_num * (*i), collision_stiffness,
			vertex_position_collider.data());
		}
		//ee collider
		for (auto i = edge_of_a_tet->begin(); i < edge_of_a_tet->end(); ++i) {
			edge_ = edge_vertices[obj_No] + ((*i) << 1);
			energy += computeEEColliderCollisionEnergyPerElement(vertex_position[*edge_].data(),
				vertex_position[*(edge_ + 1)].data(), collision.edge_edge_collider_pair_by_edge[obj_No] + collision.close_ee_collider_pair_num * (*i),
				collision.edge_edge_collider_pair_num_record[obj_No][*i], collision.edge_edge_collider_pair_d_hat[obj_No] + collision.ee_collider_d_hat_num * (*i),
				collision_stiffness, vertex_position_collider.data(), collider_edge_vertices.data());
		}
	}

	if (floor->exist) {
		for (int i = 0; i < unfixed_vertex_num; ++i) {
			if (vertex_index_on_surface[tet_actual_unfixed_vertex_indices[i]] == -1) {
				continue;
			}
			if (!indicate_vertex_collide_with_floor[tet_actual_unfixed_vertex_indices[i]]) {
				continue;
			}
			energy += compute_energy.computeFloorBarrierEnergy(vertex_position[tet_actual_unfixed_vertex_indices[i]][floor->dimension], record_vertex_collide_with_floor_d_hat[tet_actual_unfixed_vertex_indices[i]],
				collision_stiffness, floor->value);
		}
	}

	return energy;

}

double XPBD_IPC::computeVTCollisionEnergyPerElement(double* pos0, unsigned int* triangle_index, std::array<double, 3>** pos_t, std::array<int, 3>** triangle_indices,
	unsigned int num, double* d_hat, double collision_stiffness)
{
	double energy = 0.0;
	int* indices;
	for (int i = 0; i < num; i += 2) {
		indices = triangle_indices[triangle_index[i]][triangle_index[i + 1]].data();
		energy += compute_energy.computeBarrierEnergy(pos0,
			pos_t[triangle_index[i]][indices[0]].data(),
			pos_t[triangle_index[i]][indices[1]].data(),
			pos_t[triangle_index[i]][indices[2]].data(), collision_stiffness, d_hat[i>>1], true);
	}
	return energy;
}

double XPBD_IPC::computeTVCollisionEnergyPerElement(double* pos0, double* pos1, double* pos2, unsigned int* vertex_index_,
	unsigned int num, double* d_hat, double collision_stiffness,
	int* tet_unfixed_vertex_indices, int unfixed_tet_vertex_num, std::array<double, 3>** pos_v, unsigned int obj_No)
{
	double energy = 0.0;
	for (int i = 0; i < num; i += 2) {
		if (vertex_index_[i] == obj_No) {
			if (vertexInTet(unfixed_tet_vertex_num, vertex_index_[i + 1], tet_unfixed_vertex_indices)) {
				continue;
			}
		}
		energy += compute_energy.computeBarrierEnergy(pos_v[vertex_index_[i]][vertex_index_[i+1]].data(),
			pos0, pos1, pos2, collision_stiffness, d_hat[i >> 1], true);

	}
	return energy;
}


double XPBD_IPC::computeTVColliderCollisionEnergyPerElement(double* pos0, double* pos1, double* pos2, unsigned int* vertex_index_,
	unsigned int num, double* d_hat, double collision_stiffness, std::array<double, 3>** pos_v)
{
	double energy = 0.0;
	for (int i = 0; i < num; i += 2) {
		energy += compute_energy.computeBarrierEnergy(pos_v[vertex_index_[i]][vertex_index_[i + 1]].data(),
			pos0, pos1, pos2, collision_stiffness, d_hat[i >> 1], true);

	}
	return energy;
}

double XPBD_IPC::computeEEColliderCollisionEnergyPerElement(double* pos0, double* pos1, unsigned int* edge_index,
	unsigned int num, double* d_hat, double collision_stiffness, std::array<double, 3>** pos_e, unsigned int** edge_1_vertex)
{
	double energy = 0.0;
	unsigned int* edge_vertex;
	for (int i = 0; i < num; i += 2) {
		edge_vertex = edge_1_vertex[edge_index[i]] + (edge_index[i + 1] << 1);
		energy += compute_energy.computeBarrierEnergy(pos0, pos1,
			pos_e[edge_index[i]][edge_vertex[0]].data(), pos_e[edge_index[i]][edge_vertex[1]].data(),
			collision_stiffness, d_hat[i >> 1], false);
	}
	return energy;
}


double XPBD_IPC::computeEECollisionEnergyPerElement(double* pos0, double* pos1, unsigned int* edge_index,
	unsigned int num, double* d_hat, double collision_stiffness,
	int edge_order_in_tet, std::vector<unsigned int>* edge_of_a_tet, std::array<double, 3>** pos_e, unsigned int obj_No, unsigned int** edge_1_vertex)
{
	double energy = 0.0;
	unsigned int* edge_vertex;
	for (int i = 0; i < num; i += 2) {
		if (edge_index[i] == obj_No) {
			if (edgeInSameTetDuplicate(edge_order_in_tet, edge_of_a_tet, edge_index[i + 1])) {
				continue;
			}
		}
		edge_vertex = edge_1_vertex[edge_index[i]] + (edge_index[i + 1] << 1);
		energy += compute_energy.computeBarrierEnergy(	pos0, pos1, 
			pos_e[edge_index[i]][edge_vertex[0]].data(), pos_e[edge_index[i]][edge_vertex[1]].data(),
			collision_stiffness, d_hat[i >> 1], false);
	}
	return energy;
}


void XPBD_IPC::solveNewtonCD_tetBlock(std::array<double, 3>* vertex_position, double stiffness, double dt,
	double* mass,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& neighbor_tet_indices, std::array<int, 4>* indices,
	double* volume, unsigned int tet_index, std::array<double, 3>* sn, unsigned int* common_vertex_in_order,
	int* tet_vertex_index, int* unfixed_tet_vertex_index, unsigned int unfixed_vertex_num, std::vector<unsigned int>* triangle_of_a_tet,
	std::vector<unsigned int>* edge_of_a_tet, double collision_stiffness, unsigned int obj_No, int* tet_actual_unfixed_vertex_indices,
	int* vertex_index_on_surface, std::array<double, 3>* record_ori_pos, double* hessian_record, char* indicate_collide_with_floor)
{
	if (unfixed_vertex_num == 0) {
		return;
	}

	MatrixXd Hessian;
	VectorXd grad;
	grad.resize(3 * unfixed_vertex_num);
	Hessian.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);

	second_order_constraint.solveCD_ARAP_block(Hessian, grad, vertex_position, stiffness, A[tet_index],
		volume[tet_index],  tet_vertex_index,
		unfixed_tet_vertex_index, unfixed_vertex_num, hessian_record+16* tet_index);

	for (auto i = neighbor_tet_indices.begin(); i < neighbor_tet_indices.end(); ++i) {
		second_order_constraint.solveCertainHessianForNeighborTet(vertex_position, stiffness, A[*i], common_vertex_in_order, indices[*i].data(),
			Hessian, volume[*i], grad);
	}

	//MatrixXd test = Hessian;

	if (perform_collision)
	{
		MatrixXd Hessian_collision;
		Hessian_collision.resize(3 * unfixed_vertex_num, 3 * unfixed_vertex_num);
		Hessian_collision.setZero();
		getCollisionHessian(Hessian_collision, grad, triangle_of_a_tet, edge_of_a_tet, collision_stiffness, obj_No, tet_actual_unfixed_vertex_indices,
			unfixed_vertex_num,
			collision.d_hat_2, vertex_index_on_surface, vertex_position);
		//FEM::SPDprojection(Hessian_collision);
		Hessian += Hessian_collision;
	}



	//if (tet_index == 11607) {
	//	std::cout << Hessian - test << std::endl;
	//}
	double mass_dt_2;
	int vertex_index;
	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_actual_unfixed_vertex_indices[i];
		mass_dt_2 = mass[vertex_index] / (dt * dt);
		Hessian(3 * i, 3 * i) += mass_dt_2;
		Hessian(3 * i + 1, 3 * i + 1) += mass_dt_2;
		Hessian(3 * i + 2, 3 * i + 2) += mass_dt_2;
		grad.data()[3 * i] += mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
		grad.data()[3 * i + 1] += mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
		grad.data()[3 * i + 2] += mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);
	}

	//ColPivHouseholderQR<MatrixXd> linear(Hessian);
	LLT <MatrixXd> linear(Hessian);
	VectorXd result = linear.solve(grad);


	//if (tet_index == 11607) {
	//	JacobiSVD<MatrixXd > k(Hessian);
	//	std::cout << "svd "<< k.singularValues().transpose() << std::endl;
	//
	//}


	for (int i = 0; i < unfixed_vertex_num; ++i) {
		vertex_index = tet_actual_unfixed_vertex_indices[i];
		SUB_(vertex_position[vertex_index], (result.data() + 3 * i));
	}

	if (perform_collision) {

		double t = getCollisionTime(triangle_of_a_tet, edge_of_a_tet, obj_No, tet_actual_unfixed_vertex_indices,
			unfixed_vertex_num, vertex_index_on_surface, vertex_position, record_ori_pos, indicate_collide_with_floor);

		if (t < 1.0) {
			t *= 0.9;
			for (int i = 0; i < unfixed_vertex_num; ++i) {
				vertex_index = tet_actual_unfixed_vertex_indices[i];
				COLLISION_POS(vertex_position[vertex_index], t, record_ori_pos[vertex_index], vertex_position[vertex_index]);
				memcpy(record_ori_pos[vertex_index].data(), vertex_position[vertex_index].data(), 24);
			}
		}
		else {
			for (int i = 0; i < unfixed_vertex_num; ++i) {
				memcpy(record_ori_pos[tet_actual_unfixed_vertex_indices[i]].data(), vertex_position[tet_actual_unfixed_vertex_indices[i]].data(), 24);
			}
		}
	}
}



void XPBD_IPC::solveNewtonCD_tet(std::array<double, 3>* vertex_position, double stiffness, double dt,
	Matrix<double, 3, 4>* A, std::vector<unsigned int>& tet_indices, std::array<int, 4>* indices, double* mass,
	double* volume, unsigned int vertex_index, std::array<double, 3>* sn, double* lambda)
{
	Matrix3d Hessian;
	Vector3d grad;
	Hessian.setZero();
	grad.setZero();

	getARAPHessian(Hessian, grad, vertex_position, stiffness, A, tet_indices, indices, volume, vertex_index,0);


	double mass_dt_2 = mass[vertex_index] / (dt * dt);

	Hessian.data()[0] += mass_dt_2;
	Hessian(1, 1) += mass_dt_2;
	Hessian(2, 2) += mass_dt_2;
	grad.data()[0] -= mass_dt_2 * (vertex_position[vertex_index][0] - sn[vertex_index][0]);
	grad.data()[1] -= mass_dt_2 * (vertex_position[vertex_index][1] - sn[vertex_index][1]);
	grad.data()[2] -= mass_dt_2 * (vertex_position[vertex_index][2] - sn[vertex_index][2]);

	ColPivHouseholderQR <Matrix3d> linear(Hessian);
	Vector3d result = linear.solve(grad);
	SUM_(vertex_position[vertex_index], result);

}


//XPBD_IPC_VELOCITY
void XPBD_IPC::computeVelocity(int thread_No)
{
	double damp_coe = velocity_damp;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	unsigned int vertex_end = 0;
	std::array<double, 3>* velocity_;
	double delta_t = sub_time_step / damp_coe;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		initial_vertex_pos = vertex_position_render[i];
		vertex_end = vertex_index_begin_per_thread[i][thread_No + 1];
		velocity_ = velocity[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < vertex_end; ++j) {
			for (unsigned int k = 0; k < 3; ++k) {
				velocity_[j][k] = (vertex_pos[j][k] - initial_vertex_pos[j][k]) / delta_t;
			}
		}
	}
}


void XPBD_IPC::updateSn()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memcpy(sn[i][0].data(), vertex_position[i][0].data(), sn[i].size() * 24);
	}
}

//INERTIAL_ENERGY
void XPBD_IPC::inertialEnergyPerThread(int thread_No)
{
	double energy = 0.0;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* sn_;
	unsigned int vertex_end;
	double* mass;
	unsigned int start;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		sn_ = sn[i].data();
		vertex_end = vertex_index_begin_per_thread[i][thread_No+1];
		mass = mesh_struct[i]->mass.data();
		start = vertex_index_begin_per_thread[i][thread_No];
		for (unsigned int j = start; j < vertex_end; ++j) {
			energy += mass[j] * (EDGE_LENGTH(vertex_pos[j], sn_[j]));
		}
	}
	energy_per_thread[thread_No] =0.5* energy / (sub_time_step * sub_time_step);
}



double XPBD_IPC::computeInertialEnergy()
{
	double energy = 0.0;

	thread->assignTask(this, INERTIAL_ENERGY);
	for (unsigned int i = 0; i < total_thread_num; ++i) {
		energy += energy_per_thread[i];
	}
	return energy;
}
//SET_POS_PREDICT_
void XPBD_IPC::setPosPredict(int thread_No)
{
	std::array<double, 3>* vertex_pos;
	unsigned int vertex_end = 0;

	double delta_t = time_step;
	double delta_t_2 = delta_t * delta_t;
	double* mass_inv;
	std::array<double, 3>* f_ext_;
	std::array<double, 3>* velocity_;
	double gravity__[3];
	memcpy(gravity__, gravity, 24);
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		vertex_end = vertex_index_begin_per_thread[i][thread_No + 1];
		mass_inv = mesh_struct[i]->mass_inv.data();
		f_ext_ = f_ext[i].data();
		velocity_ = velocity[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < vertex_end; ++j) {
			if (mass_inv[j] != 0) {
				for (unsigned int k = 0; k < 3; ++k) {
					vertex_pos[j][k] += delta_t * velocity_[j][k] + delta_t_2 * (mass_inv[j] * f_ext_[j][k] + gravity__[k]);
				}
			}
		}
	}
}


void XPBD_IPC::resetExternalForce()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(f_ext[i][0].data(), 0, 24 * f_ext[i].size());
	}
}

void XPBD_IPC::initialDHatTolerance(double ave_edge_length)
{
	if (perform_collision) {
		collision.initialDHatTolerance(ave_edge_length);
	}
}

void XPBD_IPC::updateTetrahedronAnchorVertices()
{
	double* mass_inv;
	int* anchor_vertex;
	unsigned int anchor_vertex_size;
	for (unsigned int i = cloth->size(); i < cloth->size() + tetrahedron->size(); ++i) {
		mass_inv = mesh_struct[i]->mass_inv.data();
		anchor_vertex_size = mesh_struct[i]->anchor_vertex.size();
		anchor_vertex = mesh_struct[i]->anchor_vertex.data();
		mesh_struct[i]->resetMassInv();

	}
}



void XPBD_IPC::addExternalForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int obj_No)
{
	if (!coe.empty()) {
		for (unsigned int i = 0; i < coe.size(); ++i) {
			for (unsigned int j = 0; j < 3; ++j) {
				f_ext[obj_No][neighbor_vertex[i]][j] += coe[i] * neighbor_vertex_force_direction[j];
			}
		}
	}
}


void XPBD_IPC::updateRenderVertexNormal()
{
	for (unsigned int i = 0; i < collider->size(); ++i) {
		thread->assignTask(&(*collider)[i].mesh_struct, VERTEX_NORMAL_FROM_RENDER);
	}
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, VERTEX_NORMAL_FROM_RENDER);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, VERTEX_NORMAL_FROM_RENDER);
	}
}


void XPBD_IPC::updateRenderNormal()
{
	for (unsigned int i = 0; i < collider->size(); ++i) {
		thread->assignTask(&(*collider)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
}

void XPBD_IPC::updateNormal()
{
	for (unsigned int i = 0; i < collider->size(); ++i) {
		thread->assignTask(&(*collider)[i].mesh_struct, FACE_NORMAL);
	}
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL);
	}
}


//ARAP_ENERGY
void XPBD_IPC::computeARAPEnergyPerThread(int thread_No)
{
	double energy = 0.0;
	unsigned int end,start;
	std::array<int, 4>* indices;
	std::array<double, 3>* vertex_pos;
	Matrix<double, 3, 4>* A;
	double* volume;
	double stiffness;
	double* mass_inv_;
	double tet_energy;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		end = tet_index_begin_per_thread[i+cloth->size()][thread_No+1];
		start = tet_index_begin_per_thread[i + cloth->size()][thread_No];
		indices = tetrahedron->data()[i].mesh_struct.indices.data();
		vertex_pos = vertex_position[i + cloth->size()];
		A = tetrahedron->data()[i].mesh_struct.A.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		stiffness = tetrahedron->data()[i].ARAP_stiffness;
		mass_inv_ = tetrahedron->data()[i].mesh_struct.mass_inv.data();
		for (unsigned int j = start; j < end; ++j) {
			if (mass_inv_[indices[j][0]] != 0.0 || mass_inv_[indices[j][1]] != 0.0 || mass_inv_[indices[j][2]] != 0.0 || mass_inv_[indices[j][3]] != 0.0) {
				tet_energy = compute_energy.ARAPEnergy(vertex_pos[indices[j][0]].data(), vertex_pos[indices[j][1]].data(),
					vertex_pos[indices[j][2]].data(), vertex_pos[indices[j][3]].data(), A[j], volume[j], stiffness);
				energy += tet_energy;
			}
		}
	}
	energy_per_thread[thread_No] =0.5* energy;
}

double XPBD_IPC::computeCurrentARAPEnergy()
{
	double energy = 0.0;
	thread->assignTask(this, ARAP_ENERGY);
	for(int i=0;i<total_thread_num;++i){
		energy += energy_per_thread[i];
	}
	return energy;
}





//COLLISION_FREE_POSITION_FROM_RECORD
void XPBD_IPC::computeCollisionFreePositionFromRecord(int thread_No)
{
	double collision_time = collision.collision_time;

	if (collision_time == 1.0) {
		return;
	}
	unsigned int index_end;
	std::array<double, 3>* q_pre;
	std::array<double, 3>* q_end;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		index_end = vertex_index_begin_per_thread[i][thread_No + 1];
		q_end = vertex_position[i];
		q_pre = record_vertex_position[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			q_end[j][0] = q_pre[j][0] + collision_time * (q_end[j][0] - q_pre[j][0]);
			q_end[j][1] = q_pre[j][1] + collision_time * (q_end[j][1] - q_pre[j][1]);
			q_end[j][2] = q_pre[j][2] + collision_time * (q_end[j][2] - q_pre[j][2]);

		}
	}
}



//COLLISION_FREE_POSITION_
void XPBD_IPC::computeCollisionFreePosition(int thread_No)
{
	double collision_time = collision.collision_time;

	if (collision_time == 1.0) {
		return;
	}
	unsigned int index_end;
	std::array<double, 3>* q_pre;
	std::array<double, 3>* q_end;

	for (unsigned int i = 0; i <total_obj_num; ++i) {
		index_end = vertex_index_begin_per_thread[i][thread_No + 1];
		q_end = vertex_position[i];
		q_pre = record_collision_free_vertex_position[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < index_end; ++j) {
			q_end[j][0] = q_pre[j][0] + collision_time * (q_end[j][0] - q_pre[j][0]);
			q_end[j][1] = q_pre[j][1] + collision_time * (q_end[j][1] - q_pre[j][1]);
			q_end[j][2] = q_pre[j][2] + collision_time * (q_end[j][2] - q_pre[j][2]);
		}
	}
}


//GRAD_NORM
void XPBD_IPC::sumWithInertial(int thread_No)
{
	unsigned int index_end;
	unsigned int index_start;

	std::array<double, 3>* pos;
	std::array<double, 3>* sn_;
	unsigned int prefix_sum;
	double max = 0.0;
	double actual_grad;
	double* mass_;
	double mass_dt;
	double dt_2 = sub_time_step * sub_time_step;
	double* alread_grad = common_grad[0].data();
	std::vector<bool> * is_vertex_fixed_;
	for (int i = 0; i < total_obj_num; ++i) {
		index_end = vertex_index_begin_per_thread[i][thread_No + 1];
		index_start = vertex_index_begin_per_thread[i][thread_No];
		pos = vertex_position[i];
		sn_ = sn[i].data();
		prefix_sum = vertex_index_prefix_sum_obj[i];
		mass_ = mass[i];
		is_vertex_fixed_ = is_vertex_fixed[i];
		for (unsigned int j = index_start; j < index_end; ++j) {
			if (!(*is_vertex_fixed_)[j]) {
				mass_dt = mass_[j] / dt_2;
				actual_grad = abs(alread_grad[3 * (prefix_sum + j)] + mass_dt * (pos[j][0] - sn_[j][0]));
				if (actual_grad > max) {
					max = actual_grad;
				}
				actual_grad = abs(alread_grad[3 * (prefix_sum + j) + 1] + mass_dt * (pos[j][1] - sn_[j][1]));
				if (actual_grad > max) {
					max = actual_grad;
				}
				actual_grad = abs(alread_grad[3 * (prefix_sum + j) + 2] + mass_dt * (pos[j][2] - sn_[j][2]));
				if (actual_grad > max) {
					max = actual_grad;
				}
			}		
		}
	}

	grad_max_store[thread_No] = max;
}



//SUM_ALL_GRAD
void XPBD_IPC::sumAllGrad(int thread_No)
{
	int start = 3*global_vertex_index_start_per_thread[thread_No];
	int end =3* global_vertex_index_start_per_thread[thread_No + 1];
	for (int i = start; i < end; i++) {
		for (int j = 1; j < total_thread_num; ++j) {
			common_grad[0][i] += common_grad[j][i];				
		}
	}
}