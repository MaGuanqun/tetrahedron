#include"IPC.h"
#include"./basic/write_txt.h"
#include"XPBD/FEM_relate.h"


IPC::IPC()
{
	gravity_ = 9.8;
	sub_step_num = 1;
	iteration_number = 300;

	damping_coe = 0.0;

	perform_collision = true;

	velocity_damp = 0.995;
	energy_converge_ratio = 1e-3;

	min_inner_iteration = 4;
	min_outer_iteration = 2;

	max_move_standard = 1e-3;

	max_iteration_number = 100;
	outer_max_iteration_number = 40;
	energy_converge_standard = 1e-10;

	second_order_constraint.solve_exact_ARAP_hessian = false;
	time_step = 0.01;
}

void IPC::setForIPC(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider, Floor* floor,
	Thread* thread, double* tolerance_ratio)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	this->thread = thread;
	total_obj_num = cloth->size() + tetrahedron->size();
	total_thread_num = thread->thread_num;
	reorganzieDataOfObjects();
	energy_per_thread.resize(total_thread_num);

	initialVariable();
	setARAPWithInertial();

	//recordEdgeHessian();
	//computeGravity();
	//setHessian();

	if (perform_collision) {
		collision.initial(cloth, collider, tetrahedron, thread, floor, tolerance_ratio, IPC_, false);
	}
}



void IPC::setARAPWithInertial()
{
	std::array<double, 3>* vertex_pos;

	std::vector<Triplet<double>>* hessian_nnz_ = &hessian_nnz;
	double stiffness;
	unsigned int size;
	double* inv_mass_;

	bool is_unfixed[4];
	std::array<int, 4>* indices;
	unsigned int vertex_index_start;

	int vertex_position_in_system[4];


	double use_for_temp[9];

	double* volume;
	Matrix<double, 3, 4>* A;
	std::vector<bool>* this_is_vertex_fixed;

	hessian_nnz.reserve(18*vertex_index_prefix_sum_obj[total_obj_num]);


	for (unsigned int obj_No = 0; obj_No < tetrahedron->size(); ++obj_No) {
		vertex_pos = vertex_position[obj_No + cloth->size()];
		stiffness = 0.5*tetrahedron->data()[obj_No].ARAP_stiffness;
		size = tetrahedron->data()[obj_No].mesh_struct.indices.size();
		indices = tet_indices[obj_No];
		vertex_index_start = vertex_index_prefix_sum_obj[cloth->size() + obj_No];
		volume = tet_volume[cloth->size() + obj_No];
		A = tet_A[cloth->size() + obj_No];
		this_is_vertex_fixed = is_vertex_fixed[cloth->size() + obj_No];
		for (unsigned int i = 0; i < size; ++i) {
			for (unsigned int j = 0; j < 4; ++j) {
				is_unfixed[j] = !(* this_is_vertex_fixed)[indices[i][j]];
				vertex_position_in_system[j] = 3 * (vertex_index_start + indices[i][j]);
			}
			if (is_unfixed[0] || is_unfixed[1] || is_unfixed[2] || is_unfixed[3]) {
				computeARAPHessian(vertex_pos[indices[i][0]].data(), vertex_pos[indices[i][1]].data(), vertex_pos[indices[i][2]].data(), vertex_pos[indices[i][3]].data(),
					hessian_nnz_, vertex_position_in_system,
					stiffness * volume[i], A[i], is_unfixed, volume[i]);
			}
		}
	}

	double* mass_;
	double mass_dt_2;
	int index;
	for (int i = 0; i < total_obj_num; ++i) {
		size = mesh_struct[i]->vertex_position.size();
		mass_ = mass[i];
		vertex_index_start = vertex_index_prefix_sum_obj[i];
		for (int j = 0; j < size; ++j) {			
			mass_dt_2 = mass_[j] / (time_step * time_step);
			index = 3 * (vertex_index_start + j);
			hessian_nnz.emplace_back(Triplet<double>(index, index, mass_dt_2));
			hessian_nnz.emplace_back(Triplet<double>(index + 1, index + 1, mass_dt_2));
			hessian_nnz.emplace_back(Triplet<double>(index + 2, index + 2, mass_dt_2));
		}
	}

	Hessian.setFromTriplets(hessian_nnz.begin(), hessian_nnz.end());
	hessian_nnz.clear();
	
	double* value = Hessian.valuePtr();
	int* outer = Hessian.outerIndexPtr();
	int* inner = Hessian.innerIndexPtr();
	int col_size = Hessian.outerSize();
	int start;
	int num;

	for(int i = 0; i < col_size; ++i) {
		start = outer[i];
		if (i == col_size) {
			num = Hessian.nonZeros();
		}
		else {
			num = outer[i + 1];
		}		
		for (int j = start; j < num; ++j) {
			hessian_nnz.emplace_back(Triplet<double>(inner[j], i, value[j]));
		}
	}
	record_size_of_arap_triplet = hessian_nnz.size();
}


void IPC::computeARAPHessian(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, double* vertex_position_3,
	std::vector<Triplet<double>>* hessian_nnz,
	int* vertex_index, double stiffness, Matrix<double, 3, 4>& A, bool* is_unfixed, double volume)
{
	MatrixXd Hessian;
	Hessian.resize(4, 4);
	second_order_constraint.setARAPHessian(Hessian, stiffness, A,// ,
		volume); //,
	Hessian *= stiffness;
	for (int i = 0; i < 4; ++i) {
		if (is_unfixed[i]) {
			for (int j = 0; j < 4; ++j) {
				if (is_unfixed[j]) {
					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j], vertex_index[i], Hessian.data()[j+i*4]));
					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j]+1, vertex_index[i] + 1, Hessian.data()[j+i*4]));
					hessian_nnz->emplace_back(Triplet<double>(vertex_index[j] + 2, vertex_index[i] + 2, Hessian.data()[j+i*4]));
				}
			}
		}
	}
}

void IPC::initialVariable()
{

	Hessian.resize(3 * vertex_index_prefix_sum_obj[total_obj_num], 3 * vertex_index_prefix_sum_obj[total_obj_num]);
	b.resize(3 * vertex_index_prefix_sum_obj[total_obj_num]);
	sn.resize(3 * vertex_index_prefix_sum_obj[total_obj_num]);
	f_ext.resize(3 * vertex_index_prefix_sum_obj[total_obj_num]);
	velocity.resize(3 * vertex_index_prefix_sum_obj[total_obj_num]);
	velocity.setZero();
	acceleration.resize(3 * vertex_index_prefix_sum_obj[total_obj_num]);
	acceleration.setZero();

	grad_max_store.resize(total_thread_num);
	max_dis_record.resize(total_thread_num);
}

void IPC::reorganzieDataOfObjects()
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
		record_vertex_position_address[i] = record_vertex_position[i].data();
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
		record_collision_free_vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position;
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

		tet_color_groups[i] = &tetrahedron->data()[i].mesh_struct.tet_color_group;
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
	vertex_num_on_surface_prefix_sum.resize(total_obj_num + 1);
	vertex_index_prefix_sum_obj.resize(total_obj_num + 1, 0);

	vertex_num_on_surface_prefix_sum[0] = 0;
	for (int i = 1; i <= total_obj_num; ++i) {
		if (i <= cloth->size()) {
			vertex_num_on_surface_prefix_sum[i] = vertex_num_on_surface_prefix_sum[i - 1] + vertex_index_begin_per_thread[i - 1][total_thread_num];
		}
		else {
			vertex_num_on_surface_prefix_sum[i] = vertex_num_on_surface_prefix_sum[i - 1] + tetrahedron->data()[i - 1 - cloth->size()].mesh_struct.vertex_index_on_sureface.size();
		}

		vertex_index_prefix_sum_obj[i] = vertex_index_prefix_sum_obj[i - 1] + mesh_struct[i - 1]->vertex_for_render.size();
	}


	global_vertex_index_start_per_thread.resize(total_thread_num + 1, 0);
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

