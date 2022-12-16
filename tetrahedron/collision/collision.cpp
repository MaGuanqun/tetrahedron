#include"collision.h"
#include"../basic/basic_function.h"


Collision::Collision()
{
	d_hat = 1e-2;
	d_hat_2 = d_hat * d_hat;
	tolerance = 1e-5 * d_hat;
	collision_pair_around_pair_size_per_pair = 20;
}

void Collision::initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider,
	std::vector<Tetrahedron>* tetrahedron, Thread* thread, Floor* floor,  double* tolerance_ratio, unsigned int use_method, bool record_pair_by_element)
{

	//CCD::test();
	//dcd.test();

	this->use_method = use_method;

	testNearestPoint();
	conservative_rescaling = 0.9;

	collision_time = 1.0;
	has_collider = !collider->empty();
	tetrahedron_begin_obj_index = cloth->size();
	total_obj_num = cloth->size() + tetrahedron->size();
	total_obj_with_collider = total_obj_num + collider->size();
	this->cloth = cloth;
	this->collider = collider;
	this->tetrahedron = tetrahedron;
	this->thread = thread;
	this->floor = floor;
	thread_num = thread->thread_num;
	max_index_number_in_one_cell = 1200;
	max_index_number_in_one_cell_collider = 400;

	estimate_coeff_for_vt_pair_num = 600;
	estimate_coeff_for_vt_collider_pair_num =400;
	estimate_coeff_for_ee_pair_num =1000;
	estimate_coeff_for_tv_pair_num =500;

	estimate_coeff_for_ee_collider_pair_num = 400;
	estimate_coeff_for_tv_collider_pair_num = 300;


	estimate_coeff_for_vt_pair_num_exist = estimate_coeff_for_vt_pair_num/2;
	estimate_coeff_for_vt_collider_pair_num_exist= estimate_coeff_for_vt_collider_pair_num/2;
	estimate_coeff_for_ee_pair_num_exist = estimate_coeff_for_ee_pair_num/2;
	estimate_coeff_for_tv_pair_num_exist = estimate_coeff_for_tv_pair_num /2;

	estimate_coeff_for_ee_collider_pair_num_exist = estimate_coeff_for_ee_collider_pair_num/2;
	estimate_coeff_for_tv_collider_pair_num_exist = estimate_coeff_for_tv_collider_pair_num/2;


	close_vt_pair_num = estimate_coeff_for_vt_pair_num/4;
	close_vt_collider_pair_num =estimate_coeff_for_vt_collider_pair_num/4;
	close_ee_pair_num=estimate_coeff_for_ee_pair_num/4;
	close_tv_pair_num = estimate_coeff_for_tv_pair_num/4;

	close_ee_collider_pair_num= estimate_coeff_for_ee_collider_pair_num / 4;
	close_tv_collider_pair_num = estimate_coeff_for_tv_collider_pair_num / 4;

	this->record_pair_by_element = record_pair_by_element;
	use_BVH = false;
	reorganzieDataOfObjects();
	//findPatchOfObjects();
	if (use_BVH) {
		initialBVH(cloth, collider, tetrahedron, thread);
	}

	if (use_method == PD_) {
		initialTargetPos(cloth, tetrahedron, thread);
	}
	initialSpatialHashing(cloth, collider, tetrahedron, thread, tolerance_ratio);

	this->tolerance_ratio = tolerance_ratio;

	//testPairEven();
	if (record_pair_by_element) {
		initialPairRecord();
		initialHessianRecord();
	}
	else {
		if (use_method == PD_) {
			initialPair();
		}


	}

	if (use_BVH) {
		initialNeighborPrimitive();
	}
	collision_time_thread.resize(thread->thread_num);


	initialFloorPairInfo();
	//collision_constraint.testPT();
	//approx_CCD.test();
	//////std::cout <<"floor coordinate "<< (*collider)[0].ori_vertices[0][1] << std::endl;

	//draw_culling.setInSpatialHashingValue(spatial_hashing.spatial_hashing_value,
	//	spatial_hashing.spatial_hashing_triangle_index, spatial_hashing.spatial_hashing_value_collider,
	//	spatial_hashing.spatial_hashing_triangle_index_collider,
	//	&spatial_hashing.prefix_sum, &spatial_hashing.prefix_sum_collider, &spatial_hashing.cell_begin_per_thread);
	//the above last input variable should be actual_exist_cell_begin_per_thread(sorting) /cell_begin_per_thread(unsorting)
	//draw_culling.setInSpatialHashingValue(spatial_hashing.spatial_hashing_cell, spatial_hashing.spatial_hashing_cell_collider,
	//	spatial_hashing.hash_cell_count);
	//draw_culling->vertex_tet_pair = spatial_hashing.vertex_tet_pair.data();


	edge_edge_count.resize(thread_num, 0);
	vertex_triangle_count.resize(thread_num, 0);

	vertex_triangle_pair_index_start_per_thread.resize(2 * (thread_num + 1),0);
	vertex_obj_triangle_collider_pair_index_start_per_thread.resize(2 * (thread_num + 1),0);
	vertex_collider_triangle_obj_pair_index_start_per_thread.resize(2 * (thread_num + 1),0);
	edge_edge_pair_index_start_per_thread.resize(2 * (thread_num + 1),0);
	edge_edge_pair_collider_index_start_per_thread.resize(2 * (thread_num + 1),0);

	vertex_edge_pair_index_start_per_thread.resize(2 * (thread_num + 1), 0);
	vertex_obj_edge_collider_pair_index_start_per_thread.resize(2 * (thread_num + 1), 0);
	vertex_collider_edge_obj_pair_index_start_per_thread.resize(2 * (thread_num + 1), 0);
	vertex_vertex_pair_index_start_per_thread.resize(2 * (thread_num + 1), 0);
	vertex_vertex_pair_collider_index_start_per_thread.resize(2 * (thread_num + 1), 0);

	//target_position_element_start_per_thread.resize(2 * (thread_num + 1), 0);
	point_triangle_target_position_element_start_per_thread.resize(2 * (thread_num + 1), 0);
	point_collider_triangle_target_position_element_start_per_thread.resize(2 * (thread_num + 1), 0);
	edge_edge_target_position_element_start_per_thread.resize(2 * (thread_num + 1), 0);
	point_triangle_collider_target_position_element_start_per_thread.resize(2 * (thread_num + 1), 0);
	edge_edge_collider_target_position_element_start_per_thread.resize(2 * (thread_num + 1), 0);
	

	//target_position_start_per_thread.resize(2 * (thread_num + 1), 0);



	if (CCD_compare) {
		vertex_for_render_eigen.resize(total_obj_num);
		vertex_position_eigen.resize(total_obj_num);
		for (int i = 0; i < total_obj_num; ++i) {
			if (i < cloth->size()) {
				vertex_for_render_eigen[i].resize(cloth->data()[i].mesh_struct.vertex_for_render.size());
				vertex_position_eigen[i].resize(cloth->data()[i].mesh_struct.vertex_for_render.size());

			}
			else {
				vertex_for_render_eigen[i].resize(tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_for_render.size());
				vertex_position_eigen[i].resize(tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_for_render.size());

			}
		}
		updateEigenPosition();
	}


	if (cloth->empty()) {
		collision_stiffness = tetrahedron->data()[0].collision_stiffness;
	}
	else {
		collision_stiffness = cloth->data()[0].collision_stiffness;
	}

	//record_VT_collision_time.resize(thread_num);
	//record_EE_collision_time.resize(thread_num);
	//record_VTCollider_collision_time.resize(thread_num);

	if (use_method == XPBD_IPC_)
	{
			resizeFloorCollisionHessianRecord();
			initialPairCompress();
			recordPair();
		
	}
}



void Collision::recordPair()
{
	record_vt_pair.resize(thread_num);
	record_ee_pair.resize(thread_num);
	record_vt_collider_pair.resize(thread_num);
	record_tv_collider_pair.resize(thread_num);
	record_ee_collider_pair.resize(thread_num);

	record_vt_pair_d_hat.resize(thread_num);
	record_ee_pair_d_hat.resize(thread_num);
	record_vt_collider_pair_d_hat.resize(thread_num);
	record_tv_collider_pair_d_hat.resize(thread_num);
	record_ee_collider_pair_d_hat.resize(thread_num);

	int total_triangle_num = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		total_triangle_num += cloth->data()[i].mesh_struct.triangle_indices.size();
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		total_triangle_num += tetrahedron->data()[i].mesh_struct.triangle_indices.size();
	}
	//for (int i = 0; i < collider->size(); ++i) {
	//	total_triangle_num += collider->data()[i].mesh_struct.triangle_indices.size();
	//}
	for (unsigned int i = 0; i < thread_num; ++i) {
		record_ee_pair[i].reserve(2 * close_ee_pair_num * total_triangle_num);// 
		record_vt_pair[i].reserve(close_vt_pair_num * total_triangle_num / 2);// 
		record_ee_pair_d_hat[i].reserve( close_ee_pair_num/2 * total_triangle_num);// 
		record_vt_pair_d_hat[i].reserve(close_vt_pair_num * total_triangle_num / 8);// 


		if (has_collider) {
			record_ee_collider_pair[i].reserve(close_ee_collider_pair_num * total_triangle_num);
			record_vt_collider_pair[i].reserve(close_vt_collider_pair_num * total_triangle_num / 2);
			record_tv_collider_pair[i].reserve(close_tv_collider_pair_num * total_triangle_num / 2);

			record_ee_collider_pair_d_hat[i].reserve(close_ee_collider_pair_num * total_triangle_num/4);
			record_vt_collider_pair_d_hat[i].reserve(close_vt_collider_pair_num * total_triangle_num / 8);
			record_tv_collider_pair_d_hat[i].reserve(close_tv_collider_pair_num * total_triangle_num / 8);
		}
		else {
			record_ee_collider_pair[i].reserve(1);
			record_vt_collider_pair[i].reserve(1);
			record_tv_collider_pair[i].reserve(1);
			record_ee_collider_pair_d_hat[i].reserve(1);
			record_vt_collider_pair_d_hat[i].reserve(1);
			record_tv_collider_pair_d_hat[i].reserve(1);
		}
	}

	P1 = 73856093; //2147483647
	P2 = 19349663; //500000003


	pair_hash_table_size= 99991;//999983
	pair_hash_table_cell_size = 15;
	
	vt_hash_size_record.resize(pair_hash_table_size);
	ee_hash_size_record.resize(pair_hash_table_size);
	vt_hash_record.resize(pair_hash_table_size * pair_hash_table_cell_size);
	ee_hash_record.resize(pair_hash_table_size * pair_hash_table_cell_size);


	if (has_collider) {
		vt_collider_hash_size_record.resize(pair_hash_table_size);
		tv_collider_hash_size_record.resize(pair_hash_table_size);
		ee_collider_hash_size_record.resize(pair_hash_table_size);
		vt_collider_hash_record.resize(pair_hash_table_size * pair_hash_table_cell_size);
		tv_collider_hash_record.resize(pair_hash_table_size * pair_hash_table_cell_size);
		ee_collider_hash_record.resize(pair_hash_table_size * pair_hash_table_cell_size);
	}

	vertex_index_prefix_sum_obj.resize(total_obj_num+1,0);
	triangle_index_prefix_sum_obj.resize(total_obj_num+1,0);
	edge_index_prefix_sum_obj.resize(total_obj_num+1,0);
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_index_prefix_sum_obj[i + 1] = vertex_index_prefix_sum_obj[i] + mesh_struct[i]->vertex_for_render.size();
		triangle_index_prefix_sum_obj[i + 1] = triangle_index_prefix_sum_obj[i] + mesh_struct[i]->triangle_indices.size();
		edge_index_prefix_sum_obj[i + 1] = edge_index_prefix_sum_obj[i] + mesh_struct[i]->edge_vertices.size() / 2;
	}

	if (has_collider) {
		vertex_index_prefix_sum_obj_collider.resize(collider->size()+1, 0);
		triangle_index_prefix_sum_obj_collider.resize(collider->size()+1, 0);
		edge_index_prefix_sum_obj_collider.resize(collider->size()+1, 0);

		for (unsigned int i = 0; i < collider->size(); ++i) {
			vertex_index_prefix_sum_obj_collider[i + 1] = vertex_index_prefix_sum_obj_collider[i] + collider->data()[i].mesh_struct.vertex_for_render.size();
			triangle_index_prefix_sum_obj_collider[i + 1] = triangle_index_prefix_sum_obj_collider[i] + collider->data()[i].mesh_struct.triangle_indices.size();
			edge_index_prefix_sum_obj_collider[i + 1] = edge_index_prefix_sum_obj_collider[i] + collider->data()[i].mesh_struct.edge_vertices.size() / 2;
		}
	}



	indicate_vertex_collide_with_floor.resize(total_obj_num);
	for (int i = 0; i < total_obj_num; ++i) {		
		indicate_vertex_collide_with_floor[i].resize(vertex_index_start_per_thread[i][thread_num]);
	}
	record_vertex_collide_with_floor.resize(thread_num);
	record_vertex_collide_with_floor_d_hat.resize(thread_num);
	for (int i = 0; i < thread_num; ++i) {
		record_vertex_collide_with_floor[i].reserve(vertex_index_prefix_sum_obj[total_obj_num]/thread_num/20);
		record_vertex_collide_with_floor_d_hat[i].resize(vertex_index_start_per_thread[i][thread_num]);
	}
}




void Collision::initialPairRecordInfo()
{
	memset(vt_hash_size_record.data(), 0, 4 * pair_hash_table_size);
	memset(ee_hash_size_record.data(), 0, 4 * pair_hash_table_size);
	if (has_collider) {
		memset(vt_collider_hash_size_record.data(), 0, 4 * pair_hash_table_size);
		memset(tv_collider_hash_size_record.data(), 0, 4 * pair_hash_table_size);
		memset(ee_collider_hash_size_record.data(), 0, 4 * pair_hash_table_size);
	}
	for (int i = 0; i < total_obj_num; ++i) {
		memset(record_vertex_collide_with_floor[i].data(), 0, record_vertex_collide_with_floor[i].size());
	}

	for (int i = 0; i < thread_num; ++i) {
		record_vt_pair[i].clear();
		record_ee_pair[i].clear();
		record_vt_pair_d_hat[i].clear();
		record_ee_pair_d_hat[i].clear();
		record_vertex_collide_with_floor[i].clear();
		record_vertex_collide_with_floor_d_hat[i].clear();
		if (has_collider) {
			record_vt_collider_pair[i].clear();
			record_tv_collider_pair[i].clear();
			record_ee_collider_pair[i].clear();
			record_vt_collider_pair_d_hat[i].clear();
			record_tv_collider_pair_d_hat[i].clear();
			record_ee_collider_pair_d_hat[i].clear();
		}
	}
}


//void Collision::setParameter(std::vector<double>* lambda, double* floor_lambda, std::vector<unsigned int>* collision_lambda_index_start, double damp_stiffness, double dt)
//{
//	this->floor_lambda = floor_lambda;
//	this->lambda = lambda;
//	this->collision_lambda_index_start = collision_lambda_index_start;
//	//this->damp_stiffness= damp_stiffness;
//	this->dt = dt;
//}

void Collision::resizeFloorCollisionHessianRecord()
{
	floor_hessian_record.resize(vertex_num_on_surface_prefix_sum[total_obj_num]);
	floor_grad_record.resize(vertex_num_on_surface_prefix_sum[total_obj_num]);

	//is_vertex_collide_with_floor = new bool[vertex_num_on_surface_prefix_sum[total_obj_num]];	
}

void Collision::updateEigenPosition()
{
	for (int i = 0; i < total_obj_num; ++i) {
		if (i < cloth->size()) {
			memcpy(vertex_for_render_eigen[i][0].data(), cloth->data()[i].mesh_struct.vertex_for_render[0].data(), 24 * vertex_for_render_eigen[i].size());
			memcpy(vertex_position_eigen[i][0].data(), cloth->data()[i].mesh_struct.vertex_position[0].data(), 24 * vertex_for_render_eigen[i].size());
		}
		else {
			memcpy(vertex_for_render_eigen[i][0].data(), tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_for_render[0].data(), 24 * vertex_for_render_eigen[i].size());
			memcpy(vertex_position_eigen[i][0].data(), tetrahedron->data()[i - cloth->size()].mesh_struct.vertex_position[0].data(), 24 * vertex_for_render_eigen[i].size());
		}
	}
}


void Collision::initialFloorPairInfo()
{
	floor_collision_vertex.resize(thread_num);
	floor_collision_record.resize(thread_num);
	unsigned int vertex_num=0;
	for (int i = 0; i < total_obj_num; ++i) {
		vertex_num += vertex_index_start_per_thread[i][thread_num];
	}
	vertex_num = vertex_num / thread_num *2;
	for (int i = 0; i < thread_num; ++i) {
		floor_collision_vertex[i].reserve(vertex_num);
		floor_collision_record[i].reserve(4*vertex_num);
	}
	//int total_triangle_num = 0;
	//for (int i = 0; i < cloth->size(); ++i) {
	//	total_triangle_num += cloth->data()[i].mesh_struct.triangle_indices.size();
	//}
	//for (int i = 0; i < tetrahedron->size(); ++i) {
	//	total_triangle_num += tetrahedron->data()[i].mesh_struct.triangle_indices.size();
	//}

	//point_triangle_pair = new unsigned int* [thread_num];
	//point_obj_triangle_collider_pair = new unsigned int* [thread_num];
	//point_collider_triangle_obj_pair = new unsigned int* [thread_num];
	//edge_edge_pair = new unsigned int* [thread_num];
	//edge_edge_collider_pair = new unsigned int* [thread_num];



	//for (int i = 0; i < thread_num; ++i) {
	//	point_triangle_pair[i] = new unsigned int[2 * max_index_number_in_one_cell * total_triangle_num];
	//	edge_edge_pair[i] = new unsigned int[4 * max_index_number_in_one_cell * total_triangle_num];

	//	memset(point_triangle_pair[i], 0, 8 * max_index_number_in_one_cell * total_triangle_num);
	//	memset(edge_edge_pair[i], 0, 16 * max_index_number_in_one_cell * total_triangle_num);

	//	if (has_collider) {
	//		point_obj_triangle_collider_pair[i] = new unsigned int[2 * max_index_number_in_one_cell_collider * total_triangle_num];
	//		edge_edge_collider_pair[i] = new unsigned int[4 * max_index_number_in_one_cell_collider * total_triangle_num];
	//		point_collider_triangle_obj_pair[i] = new unsigned int[2 * max_index_number_in_one_cell_collider * total_triangle_num];

	//		memset(point_obj_triangle_collider_pair[i], 0, 8 * max_index_number_in_one_cell_collider * total_triangle_num);
	//		memset(point_collider_triangle_obj_pair[i], 0, 8 * max_index_number_in_one_cell_collider * total_triangle_num);
	//		memset(edge_edge_collider_pair[i], 0, 16 * max_index_number_in_one_cell_collider * total_triangle_num);
	//	}
	//	else {
	//		point_obj_triangle_collider_pair[i] = new unsigned int[1];
	//		point_collider_triangle_obj_pair[i] = new unsigned int[1];
	//		edge_edge_collider_pair[i] = new unsigned int[1];
	//		memset(point_obj_triangle_collider_pair[i], 0, 4);
	//		memset(point_collider_triangle_obj_pair[i], 0, 4);
	//		memset(edge_edge_collider_pair[i], 0, 4);
	//	}
	//}
}

void Collision::setCollisionFreeVertex(std::vector<std::array<double, 3>*>* record_vertex_position, std::vector<std::vector<std::array<double, 3>>>* record_vertex_for_this_color)
{
	vertex_collision_free.resize(total_obj_num);
	for (int i = 0; i < total_obj_num; ++i) {
		vertex_collision_free[i] = record_vertex_position->data()[i];
	}

	vertex_record_for_this_color.resize(total_obj_num);
	for (int i = 0; i < total_obj_num; ++i) {
		vertex_record_for_this_color[i] = record_vertex_for_this_color->data()[i].data();
	}

}
void Collision::setCollisionFreeVertex(std::vector< std::vector<std::array<double, 3>>>* record_vertex_position)
{
	vertex_collision_free.resize(total_obj_num);
	for (int i = 0; i < total_obj_num; ++i) {
		vertex_collision_free[i] = record_vertex_position->data()[i].data();
	}
}

void Collision::reorganzieDataOfObjects()
{
	obj_tri_aabb.resize(total_obj_num);
	vertex_aabb.resize(total_obj_num);
	edge_aabb.resize(total_obj_num);
	representative_vertex_num.resize(total_obj_num);
	representative_edge_num.resize(total_obj_num);
	triangle_index_in_order.resize(total_obj_num);
	faces.resize(total_obj_num);
	edges.resize(total_obj_num);

	face_edges.resize(total_obj_num);
	edge_vertices.resize(total_obj_num);
	vertex_for_render.resize(total_obj_num);
	vertex_position.resize(total_obj_num);
	triangle_indices.resize(total_obj_num);
	triangle_normal_render.resize(total_obj_num);
	triangle_normal.resize(total_obj_num);

	triangle_normal_not_normalized.resize(total_obj_num);
	triangle_normal_render_not_normalized.resize(total_obj_num);
	cross_for_approx_CCD.resize(total_obj_num);
	mass.resize(total_obj_num);
	mass_inv.resize(total_obj_num);

	triangle_normal_magnitude_reciprocal.resize(total_obj_num);

	vertex_index_start_per_thread.resize(total_obj_num);
	all_vertex_index_start_per_thread.resize(total_obj_num);
	edge_index_start_per_thread.resize(total_obj_num);
	triangle_index_start_per_thread.resize(total_obj_num);


	mesh_struct.resize(total_obj_num);

	damp_collision.resize(total_obj_num);

	vertex_index_on_surface.resize(total_obj_num);
	general_index_to_surface_index.resize(total_obj_num);

	rest_edge_length.resize(total_obj_num);


	triangle_index_of_a_tet_color_group.resize(total_obj_num);
	edge_index_of_a_tet_color_group.resize(total_obj_num);
	surface_vertex_index_of_a_tet_color_group.resize(total_obj_num);
	vertex_index_of_a_tet_color_group.resize(total_obj_num);

	triangle_index_of_a_tet_color_per_thread_start_group.resize(total_obj_num);
	edge_index_of_a_tet_color_per_thread_start_group.resize(total_obj_num);
	surface_vertex_index_of_a_tet_color_per_thread_start_group.resize(total_obj_num);
	vertex_index_of_a_tet_color_per_thread_start_group.resize(total_obj_num);

	total_vertex_num.resize(total_obj_num);

	tet_color_groups.resize(total_obj_num);
	tet_color_groups_label.resize(total_obj_num);


	triangle_around_triangle.resize(total_obj_num);
	edge_around_triangle.resize(total_obj_num);
	tet_around_vertex.resize(total_obj_num);
	tet_around_triangle.resize(total_obj_num);
	triangle_around_edge.resize(total_obj_num);
	edge_around_edge.resize(total_obj_num);
	tet_around_edge.resize(total_obj_num);

	tet_order_in_color_group.resize(total_obj_num);

	for (int i = 0; i < cloth->size(); ++i) {
		obj_tri_aabb[i] = cloth->data()[i].triangle_AABB.data();
		vertex_aabb[i] = cloth->data()[i].vertex_AABB.data();
		edge_aabb[i] = cloth->data()[i].edge_AABB.data();
		representative_vertex_num[i] = cloth->data()[i].representative_vertex_num.data();
		representative_edge_num[i] = cloth->data()[i].representative_edge_num.data();
		triangle_index_in_order[i] = cloth->data()[i].mesh_struct.surface_triangle_index_in_order.data();
		faces[i] = cloth->data()[i].mesh_struct.faces.data();
		edges[i] = cloth->data()[i].mesh_struct.edges.data();
		face_edges[i] = cloth->data()[i].mesh_struct.face_edges.data();
		edge_vertices[i] = cloth->data()[i].mesh_struct.edge_vertices.data();
		vertex_for_render[i] = cloth->data()[i].mesh_struct.vertex_for_render.data();
		vertex_position[i] = cloth->data()[i].mesh_struct.vertex_position.data();
		triangle_normal_render[i] = cloth->data()[i].mesh_struct.face_normal_for_render.data();
		triangle_normal[i] = cloth->data()[i].mesh_struct.face_normal.data();
		triangle_indices[i] = cloth->data()[i].mesh_struct.triangle_indices.data();
		triangle_normal_not_normalized[i] = cloth->data()[i].mesh_struct.ori_face_normal.data();
		triangle_normal_render_not_normalized[i] = cloth->data()[i].mesh_struct.ori_face_normal_for_render.data();
		cross_for_approx_CCD[i] = cloth->data()[i].mesh_struct.cross_for_approx_CCD.data();
		mass[i] = cloth->data()[i].mesh_struct.mass.data();
		mass_inv[i] = cloth->data()[i].mesh_struct.mass_inv.data();
		triangle_normal_magnitude_reciprocal[i] = cloth->data()[i].mesh_struct.triangle_normal_magnitude_reciprocal.data();
		vertex_index_start_per_thread[i] = cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		all_vertex_index_start_per_thread[i] = cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		edge_index_start_per_thread[i] = cloth->data()[i].mesh_struct.edge_index_begin_per_thread.data();
		triangle_index_start_per_thread[i] = cloth->data()[i].mesh_struct.face_index_begin_per_thread.data();
		damp_collision[i] = cloth->data()[i].damp_collision_stiffness;
		mesh_struct[i] = &cloth->data()[i].mesh_struct;

		rest_edge_length[i] = cloth->data()[i].mesh_struct.edge_length.data();

		total_vertex_num[i] = cloth->data()[i].mesh_struct.vertex_position.size();

		triangle_around_triangle[i] = cloth->data()[i].mesh_struct.face_around_face.data();
		edge_around_triangle[i] = cloth->data()[i].mesh_struct.edge_around_face.data();
		triangle_around_edge[i] = cloth->data()[i].mesh_struct.face_around_edge.data();
		edge_around_edge[i] = cloth->data()[i].mesh_struct.edge_around_edge.data();

	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		obj_tri_aabb[i + cloth->size()] = tetrahedron->data()[i].triangle_AABB.data();
		vertex_aabb[i + cloth->size()] = tetrahedron->data()[i].vertex_AABB.data();
		edge_aabb[i + cloth->size()] = tetrahedron->data()[i].edge_AABB.data();
		representative_vertex_num[i + cloth->size()] = tetrahedron->data()[i].representative_vertex_num.data();
		representative_edge_num[i + cloth->size()] = tetrahedron->data()[i].representative_edge_num.data();
		triangle_index_in_order[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.surface_triangle_index_in_order.data();
		faces[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.faces.data();
		edges[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edges.data();
		face_edges[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.face_edges.data();
		edge_vertices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_vertices.data();

		vertex_for_render[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.data();
		vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		triangle_normal_render[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.face_normal_for_render.data();
		triangle_normal[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.face_normal.data();
		triangle_indices[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.triangle_indices.data();

		triangle_normal_not_normalized[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.ori_face_normal.data();
		triangle_normal_render_not_normalized[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.ori_face_normal_for_render.data();
		cross_for_approx_CCD[i+ cloth->size()] = tetrahedron->data()[i].mesh_struct.cross_for_approx_CCD.data();
		mass[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.mass.data();
		mass_inv[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.mass_inv.data();
		triangle_normal_magnitude_reciprocal[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.triangle_normal_magnitude_reciprocal.data();
		vertex_index_start_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_on_surface_begin_per_thread.data();
		all_vertex_index_start_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
		edge_index_start_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_index_begin_per_thread.data();
		triangle_index_start_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.face_index_begin_per_thread.data();
		damp_collision[i + cloth->size()] = tetrahedron->data()[i].damp_collision_stiffness;
		mesh_struct[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct;
		vertex_index_on_surface[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_on_sureface.data();
		general_index_to_surface_index[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_surface_index.data();

		rest_edge_length[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_length.data();

		triangle_index_of_a_tet_color_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.triangle_index_of_a_tet_color_group.data();
		edge_index_of_a_tet_color_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_index_of_a_tet_color_group.data();
		surface_vertex_index_of_a_tet_color_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.surface_vertex_index_of_a_tet_color_group.data();
		vertex_index_of_a_tet_color_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_of_a_tet_color_group.data();

		triangle_index_of_a_tet_color_per_thread_start_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.triangle_index_of_a_tet_color_per_thread_start_group.data();
		edge_index_of_a_tet_color_per_thread_start_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_index_of_a_tet_color_per_thread_start_group.data();
		surface_vertex_index_of_a_tet_color_per_thread_start_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.surface_vertex_index_of_a_tet_color_per_thread_start_group.data();
		vertex_index_of_a_tet_color_per_thread_start_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_of_a_tet_color_per_thread_start_group.data();
		total_vertex_num[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.size();

		tet_color_groups[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct.tet_color_group;
		tet_color_groups_label[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tet_in_collision.data();


		triangle_around_triangle[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.face_around_face.data();
		edge_around_triangle[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_around_face.data();
		tet_around_vertex[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_tet_index.data();
		tet_around_triangle[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tet_around_face.data();

		triangle_around_edge[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.face_around_edge.data();
		edge_around_edge[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.edge_around_edge.data();
		tet_around_edge[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tet_around_edge.data();
		tet_order_in_color_group[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.tet_order_in_color_group.data();

	}

	if (has_collider) {
		obj_tri_aabb_collider.resize(collider->size());
		collider_face_edges.resize(collider->size());
		collider_edge_vertices.resize(collider->size());
		representative_vertex_num_collider.resize(collider->size());
		representative_edge_num_collider.resize(collider->size());
		triangle_index_in_order_collider.resize(collider->size());

		vertex_aabb_collider.resize(collider->size());
		edge_aabb_collider.resize(collider->size());

		vertex_for_render_collider.resize(collider->size());
		vertex_position_collider.resize(collider->size());
		triangle_normal_render_collider.resize(collider->size());
		triangle_normal_collider.resize(collider->size());
		triangle_indices_collider.resize(collider->size());
		triangle_normal_not_normalized_collider.resize(collider->size());
		triangle_normal_render_not_normalized_collider.resize(collider->size());
		cross_for_approx_CCD_collider.resize(collider->size());
		collider_triangle_normal_magnitude_reciprocal.resize(collider->size());

		rest_edge_length_collider.resize(collider->size());

		for (int i = 0; i < collider->size(); ++i) {
			obj_tri_aabb_collider[i] = collider->data()[i].triangle_AABB.data();
			vertex_aabb_collider[i] = collider->data()[i].vertex_AABB.data();
			edge_aabb_collider[i] = collider->data()[i].edge_AABB.data();
			collider_face_edges[i] = collider->data()[i].mesh_struct.face_edges.data();
			collider_edge_vertices[i] = collider->data()[i].mesh_struct.edge_vertices.data();
			representative_vertex_num_collider[i] = collider->data()[i].representative_vertex_num.data();
			representative_edge_num_collider[i] = collider->data()[i].representative_edge_num.data();
			triangle_index_in_order_collider[i] = collider->data()[i].mesh_struct.surface_triangle_index_in_order.data();

			vertex_for_render_collider[i] = collider->data()[i].mesh_struct.vertex_for_render.data();
			vertex_position_collider[i] = collider->data()[i].mesh_struct.vertex_position.data();
			triangle_normal_render_collider[i] = collider->data()[i].mesh_struct.face_normal_for_render.data();
			triangle_normal_collider[i] = collider->data()[i].mesh_struct.face_normal.data();
			triangle_indices_collider[i] = collider->data()[i].mesh_struct.triangle_indices.data();

			triangle_normal_not_normalized_collider[i] = collider->data()[i].mesh_struct.ori_face_normal.data();
			triangle_normal_render_not_normalized_collider[i] = collider->data()[i].mesh_struct.ori_face_normal_for_render.data();
			cross_for_approx_CCD_collider[i] = collider->data()[i].mesh_struct.cross_for_approx_CCD.data();
			collider_triangle_normal_magnitude_reciprocal[i] = collider->data()[i].mesh_struct.triangle_normal_magnitude_reciprocal.data();

			rest_edge_length_collider[i] = collider->data()[i].mesh_struct.edge_length.data();
		}
	}

}

void Collision::initialDHatTolerance(double ave_edge_length)
{
	for (int i = 0; i < 4; ++i) {
		tolerance_radius[i] = tolerance_ratio[i] * ave_edge_length;
	}	

	volume_boundary = 0.5 * d_hat*ave_edge_length;



	eta = 0.01;	
	tolerance_2 = tolerance * tolerance;
	epsilon = 0.001;

}


void Collision::findPatchOfObjects()
{
	//getAABBWithoutTolerance();
	//mesh_patch.initialPatch(cloth, collider, tetrahedron, thread);
	//mesh_patch.setBuffer(0,tetrahedron_begin_obj_index);

}



void Collision::initialNeighborPrimitive()
{
	for (int i = 0; i < cloth->size(); ++i) {
		(*cloth)[i].initialNeighborPrimitiveRecording(cloth->size(), tetrahedron->size(), collider->size(), true);
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		(*tetrahedron)[i].initialNeighborPrimitiveRecording(cloth->size(), tetrahedron->size(), collider->size(), true);
	}
}


void Collision::initialTargetPos(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, Thread* thread)
{
	
	obj_target_pos_per_thread.resize(thread_num);
	for (int i = 0; i < obj_target_pos_per_thread.size(); ++i) {
		obj_target_pos_per_thread[i].initialSet(total_obj_num);
		for (int j = 0; j < total_obj_num; ++j) {
			if (j < tetrahedron_begin_obj_index) {
				obj_target_pos_per_thread[i].initialSet2(j, (*cloth)[j].ori_vertices.size());
			}
			else {
				obj_target_pos_per_thread[i].initialSet2(j, (*tetrahedron)[j - tetrahedron_begin_obj_index].ori_vertices.size());
			}
		}
		obj_target_pos_per_thread[i].initial();
	}
	obj_target_pos.initialSet(total_obj_num);
	for (int j = 0; j < total_obj_num; ++j) {
		if (j < tetrahedron_begin_obj_index) {
			obj_target_pos.initialSet2(j, (*cloth)[j].ori_vertices.size());
		}
		else {
			obj_target_pos.initialSet2(j, (*tetrahedron)[j - tetrahedron_begin_obj_index].ori_vertices.size());
		}
	}
	obj_target_pos.initial();
}

void Collision::initialBVH(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron, Thread* thread)
{
	obj_BVH.resize(total_obj_num);
	collider_BVH.resize(collider->size());
	for (int i = 0; i < total_obj_num; ++i) {
		if (i < tetrahedron_begin_obj_index) {
			obj_BVH[i].init((*cloth)[i].mesh_struct.faces.size(), (*cloth)[i].mesh_struct.face_index_begin_per_thread, thread);
		}
		else {
			obj_BVH[i].init((*tetrahedron)[i - tetrahedron_begin_obj_index].mesh_struct.triangle_indices.size(), (*tetrahedron)[i - tetrahedron_begin_obj_index].mesh_struct.face_index_begin_per_thread, thread);
		}
	}
	for (int i = 0; i < collider->size(); ++i) {
		collider_BVH[i].init((*collider)[i].mesh_struct.faces.size(), (*collider)[i].mesh_struct.face_index_begin_per_thread, thread);
	}
}

void Collision::initialSpatialHashing(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron, Thread* thread,
	double* tolerance_ratio)
{
	spatial_hashing.setInObject(cloth, collider, tetrahedron, thread, tolerance_ratio, 8, 
		max_index_number_in_one_cell, max_index_number_in_one_cell_collider, estimate_coeff_for_vt_pair_num, estimate_coeff_for_vt_collider_pair_num,
		estimate_coeff_for_tv_pair_num, estimate_coeff_for_ee_pair_num, record_pair_by_element, estimate_coeff_for_tv_collider_pair_num,
		estimate_coeff_for_ee_collider_pair_num);
}


void Collision::getAABB()
{
	for (int i = 0; i < cloth->size(); ++i) {
		(*cloth)[i].obtainAABB(true);
	}
	for (int i = 0; i < collider->size(); ++i) {
		(*collider)[i].obtainAABB(true);
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		(*tetrahedron)[i].obtainAABB(true);
	}
	//thread->assignTask(&mesh_patch, PATCH_AABB);
}

void Collision::getAABBWithoutTolerance()
{
	for (int i = 0; i < cloth->size(); ++i) {
		(*cloth)[i].obtainAABB(false);
	}
	for (int i = 0; i < collider->size(); ++i) {
		(*collider)[i].obtainAABB(false);
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		(*tetrahedron)[i].obtainAABB(false);
	}
}

void Collision::buildBVH()
{
	bool a[2];
	double* aabb;
	for (int i = 0; i < total_obj_num; ++i) {
		if (i < tetrahedron_begin_obj_index) {
			aabb = cloth->data()[i].obj_aabb;
		}
		else {
			aabb = tetrahedron->data()[i - tetrahedron_begin_obj_index].obj_aabb;
		}
		obj_BVH[i].buildBVH(obj_tri_aabb[i]);
	}
	//testBVHUpdate();
	for (int i = 0; i < collider->size(); ++i) {
		collider_BVH[i].buildBVH(obj_tri_aabb_collider[i]);
	}
}




void Collision::globalCollision()
{
	getAABB();
	time_t t1 = clock();
	getSceneAABB();
	//buildBVH();	
		//for (int i = 0; i < 100; ++i) {
	spatial_hashing.buildSpatialHashing(scene_aabb);
	//}				
//////std::cout << "build " << clock() - t1 << std::endl;
//t1 = clock();
//for (int i = 0; i < 100; ++i) {
//thread->assignTask(this, FIND_TRIANGLE_PAIRS);
//}
//////std::cout << "find triangle pair " << clock() - t1 << std::endl;

//testCollision();

//t1 = clock();
//for (int i = 0; i < 100; ++i) {
	//thread->assignTask(this, FIND_PRIMITIVE_AROUND);

	//}
	//////std::cout << "find around primitive " << clock() - t1 << std::endl;
	////////std::cout << "search " << clock() - t1 << std::endl;




	thread->assignTask(this, GLOBAL_COLLISION_DETECTION);
	sumTargetPosition();
}


void Collision::totalCount()
{
	//unsigned int vertex_triangle_count_total = 0;
	//unsigned int edge_edge_count_count_total = 0;

	//unsigned int triangle_triangle_count_total = 0;

	unsigned int vertex_triangle_count_total_final = 0;
	unsigned int edge_edge_count_count_total_final = 0;

	for (int i = 0; i < thread_num; ++i) {
		vertex_triangle_count_total_final += spatial_hashing.vertex_triangle_pair[i][0] >> 2;
		edge_edge_count_count_total_final += spatial_hashing.edge_edge_pair[i][0] >> 2;
		//triangle_triangle_count_total += spatial_hashing.triangle_pair[i][0] >> 2;
	}


}

void Collision::testRepeatability()
{
	unsigned int pair_num = 0;
	unsigned int edge_pair_num = 0;
	for (int i = 0; i < thread_num; ++i) {
		pair_num += spatial_hashing.vertex_triangle_pair[i][0] / 4;
		edge_pair_num += spatial_hashing.edge_edge_pair[i][0] / 4;
	}
	std::vector<TriangleElementPair> triangle_pair(pair_num);
	std::vector<TriangleElementPair> edge_pair(edge_pair_num);
	unsigned int index = 0;
	unsigned int index_edge = 0;

	unsigned int* tri_pair_;
	unsigned int* edge_pair_;

	for (int i = 0; i < thread_num; ++i) {
		tri_pair_ = spatial_hashing.vertex_triangle_pair[i] + 1;
		edge_pair_ = spatial_hashing.edge_edge_pair[i] + 1;
		for (int j = 0; j < spatial_hashing.vertex_triangle_pair[i][0]; j += 4) {
			memcpy(triangle_pair[index].index, tri_pair_ + j, 16);
			index++;
		}
		for (int j = 0; j < spatial_hashing.edge_edge_pair[i][0]; j += 4) {

			if (edge_pair_[j + 1] < edge_pair_[j + 3] ||
				(edge_pair_[j + 1] == edge_pair_[j + 3] &&
					edge_pair_[j] < edge_pair_[j + 2])) {
				memcpy(edge_pair[index_edge].index, edge_pair_ + j, 16);
			}
			else {
				edge_pair[index_edge].index[0] = edge_pair_[j + 2];
				edge_pair[index_edge].index[1] = edge_pair_[j + 3];
				edge_pair[index_edge].index[2] = edge_pair_[j];
				edge_pair[index_edge].index[3] = edge_pair_[j + 1];
			}
			index_edge++;
		}
	}
	std::sort(triangle_pair.begin(), triangle_pair.end());
	std::sort(edge_pair.begin(), edge_pair.end());


	TriangleElementPair a;
	std::vector<unsigned int> count;
	std::vector<unsigned int> count_edge;
	a = triangle_pair[0];
	count.push_back(1);
	for (int i = 1; i < triangle_pair.size(); ++i) {
		if (triangle_pair[i] == a) {
			count.back()++;
		}
		else {
			a = triangle_pair[i];
			count.push_back(1);
		}
	}

	a = edge_pair[0];
	count_edge.push_back(1);
	for (int i = 1; i < edge_pair.size(); ++i) {
		if (edge_pair[i] == a) {
			count_edge.back()++;
		}
		else {
			a = edge_pair[i];
			count_edge.push_back(1);
		}
	}

	int count_2 = 0;
	int count_2_collider = 0;
	int count_2_larger = 0;
	int count_2_larger_collider = 0;
	for (int i = 0; i < count.size(); ++i) {
		if (count[i] == 2) {
			count_2++;
		}
		else if (count[i] > 2) {
			count_2_larger++;
		}
	}
	for (int i = 0; i < count_edge.size(); ++i) {
		if (count_edge[i] == 2) {
			count_2_collider++;
		}
		else if (count_edge[i] > 2) {
			count_2_larger_collider++;
		}
	}



	////std::cout << "vertex triangle Repeatability 2: " << (double)count_2 / (double)count.size() << " larger than 2 " << (double)count_2_larger / (double)count.size() << std::endl;
	////std::cout << count_2 << " " << count_2_larger << " " << count.size() << std::endl;
	////std::cout << "total VT pair " << triangle_pair.size() << std::endl;
	////std::cout << "edge edge Repeatability 2: " << (double)count_2_collider / (double)count_edge.size() << " larger than 2 " << (double)count_2_larger_collider / (double)count_edge.size() << std::endl;
	////std::cout << count_2_collider << " " << count_2_larger_collider << " " << count_edge.size() << std::endl;

	////std::cout << "total edge pair " << edge_pair.size() << std::endl;

	std::vector<std::vector<unsigned int>> prfix_sum_vertex(total_obj_num);
	std::vector<std::vector<unsigned int>> prfix_sum_edge(total_obj_num);	
	for (int i = 0; i < total_obj_num; ++i) {
		if (i < cloth->size()) {
			prfix_sum_vertex[i].resize(cloth->data()[i].mesh_struct.vertex_position.size() + 1, 0);
			prfix_sum_edge[i].resize(cloth->data()[i].mesh_struct.edges.size() + 1, 0);
		}
		else {
			prfix_sum_vertex[i].resize(tetrahedron->data()[i- cloth->size()].mesh_struct.vertex_position.size() + 1, 0);
			prfix_sum_edge[i].resize(tetrahedron->data()[i- cloth->size()].mesh_struct.edges.size() + 1, 0);
		}
	}
	for (int i = 0; i < triangle_pair.size(); ++i) {
		prfix_sum_vertex[triangle_pair[i].index[1]][triangle_pair[i].index[0] + 1]++;
	}
	for (int i = 0; i < edge_pair.size(); ++i) {
		prfix_sum_edge[edge_pair[i].index[1]][edge_pair[i].index[0] + 1]++;
	}
	//for (int i = 0; i < total_obj_num; ++i) {
		for (int j = 1; j < prfix_sum_vertex[0].size() + 1; ++j) {
			prfix_sum_vertex[0][j] += prfix_sum_vertex[0][j - 1];
		}
		for (int j = 1; j < prfix_sum_edge[0].size() + 1; ++j) {
			prfix_sum_edge[0][j] += prfix_sum_edge[0][j - 1];
		}
	//}
		for (int i = 1; i < total_obj_num; ++i) {
			prfix_sum_vertex[i][0] = prfix_sum_vertex[i - 1].back();
			prfix_sum_edge[i][0] = prfix_sum_edge[i - 1].back();
			for (int j = 1; j < prfix_sum_vertex[i].size() + 1; ++j) {
				prfix_sum_vertex[i][j] += prfix_sum_vertex[i][j - 1];
			}
			for (int j = 1; j < prfix_sum_edge[i].size() + 1; ++j) {
				prfix_sum_edge[i][j] += prfix_sum_edge[i][j - 1];
			}
		}



	int surface_vertex_index_on_global;
	std::vector<std::vector<std::vector<int>>>* tri_obj_tri;
	std::vector<std::vector<std::vector<int>>>* edge_obj_edge;
	bool need_break;
	unsigned int tri_pair_num_ = 0;
	unsigned int edge_pair_num_ = 0;
	for (int m = cloth->size(); m < total_obj_num; ++m) { //obj 0
		if (m < cloth->size()) {
			tri_obj_tri = &cloth->data()[m].vertex_neighbor_obj_triangle;
		}
		else {
			tri_obj_tri = &tetrahedron->data()[m-cloth->size()].surface_vertex_neighbor_obj_triangle;
		}		
		for (int i = 0; i < tri_obj_tri->size(); ++i) {		//vertex 0			
			surface_vertex_index_on_global = vertex_index_on_surface[m][i];
			for (int j = 0; j < tri_obj_tri->data()[i].size(); ++j) {		//obj 1		
				for (int k = 0; k < tri_obj_tri->data()[i][j].size(); ++k) {  // triangle 1
					tri_pair_num_++;
					//obj_index, i, j, tri_obj->data()[i][j][k]
					need_break = false;
					for (int l = prfix_sum_vertex[m][surface_vertex_index_on_global]; l < prfix_sum_vertex[m][surface_vertex_index_on_global + 1]; ++l) {
						if (triangle_pair[l].index[3] == j && triangle_pair[l].index[2] == tri_obj_tri->data()[i][j][k]) {
							need_break = true;
							break;
						}
					}
					if (!need_break) {
						////std::cout << "does not find VT pair: " << m << " " << surface_vertex_index_on_global << " " << j << " " << tri_obj_tri->data()[i][j][k] << std::endl;
					}
				}
			}
		}
		edge_obj_edge = &tetrahedron->data()[m - cloth->size()].edge_neighbor_obj_edge;
		for (int i = 0; i < edge_obj_edge->size(); ++i) {		//edge 0			
			for (int j = 0; j < edge_obj_edge->data()[i].size(); ++j) {		//obj 1		
				for (int k = 0; k < edge_obj_edge->data()[i][j].size(); ++k) {  // edge 1
					//obj_index, i, j, tri_obj->data()[i][j][k]
					edge_pair_num_++;
					need_break = false;
					for (int l = prfix_sum_edge[m][i]; l < prfix_sum_edge[m][i + 1]; ++l) {
						if (edge_pair[l].index[3] == j && edge_pair[l].index[2] == edge_obj_edge->data()[i][j][k]) {
							need_break = true;
							break;
						}
					}
					if (!need_break) {
						////std::cout << "does not find EE pair: " << m << " " << i << " " << j << " " << edge_obj_edge->data()[i][j][k] << std::endl;
					}
				}
			}
		}
	}
	////std::cout << "edge pair num " << edge_pair_num_ << std::endl;
	////std::cout << "triangle pair num " << tri_pair_num_ << std::endl;
	////std::cout << "SH is correct " << std::endl;

	//for (int i = 0; i < thread_num; ++i) {
	//	for (int j = 0; j < spatial_hashing.triangle_pair[i].size(); j+=4) {
	//		if (spatial_hashing.triangle_pair[i][j] == 1 || spatial_hashing.triangle_pair[i][j+2]==1) {
	//			////std::cout <<i<<" "<< spatial_hashing.triangle_pair[i][j + 1] << " " << spatial_hashing.triangle_pair[i][j] << " "
	//				<< spatial_hashing.triangle_pair[i][j + 3] << " " << spatial_hashing.triangle_pair[i][j + 2] << std::endl;
	//		}
	//	}
	//}
	//////std::cout << cloth->data()[0].triangle_neighbor_obj_triangle[0][0].size() << std::endl;
	//for (int i = 0; i < cloth->data()[0].triangle_neighbor_obj_triangle[0][0].size(); ++i) {
	//	////std::cout << cloth->data()[0].triangle_neighbor_obj_triangle[0][0][i]<<" " << std::endl;
	//}
	//////std::cout <<(int) AABB::AABB_intersection(obj_tri_aabb[0][0].data(), obj_tri_aabb[0][6].data()) << std::endl;

	

	//for (int i = 0; i < triangle_pair.size(); ++i) {
	//	////std::cout << triangle_pair[i].index[1] << " " << triangle_pair[i].index[0] << " " << triangle_pair[i].index[3] << " " << triangle_pair[i].index[2] << std::endl;
	//}
	//TriangleElementPair a;
	//std::vector<unsigned int> count;
	//std::vector<unsigned int> count_collider;
	//a = triangle_pair[0];
	//count.push_back(1);
	//for (int i = 1; i < triangle_pair.size(); ++i) {
	//	if (triangle_pair[i] == a) {
	//		count.back()++;
	//	}
	//	else {
	//		a = triangle_pair[i];
	//		count.push_back(1);
	//	}
	//}
	//if (!triangle_pair_collider.empty()) {
	//	a = triangle_pair_collider[0];
	//	count_collider.push_back(1);
	//	for (int i = 1; i < triangle_pair_collider.size(); ++i) {
	//		if (triangle_pair_collider[i] == a) {
	//			count_collider.back()++;
	//		}
	//		else {
	//			a = triangle_pair_collider[i];
	//			count_collider.push_back(1);
	//		}
	//	}
	//}
	//int count_2 = 0;
	//int count_2_collider = 0;
	//int count_2_larger = 0;
	//int count_2_larger_collider = 0;
	//for (int i = 0; i < count.size(); ++i) {
	//	if (count[i] == 2) {
	//		count_2++;
	//	}
	//	else if (count[i] > 2) {
	//		count_2_larger++;
	//	}
	//}
	//if (!triangle_pair_collider.empty()) {
	//	for (int i = 0; i < count_collider.size(); ++i) {
	//		if (count_collider[i] == 2) {
	//			count_2_collider++;
	//		}
	//		else if (count_collider[i] > 2) {
	//			count_2_larger_collider++;
	//		}
	//	}
	//}
	//////std::cout << "Repeatability 2: " << (double)count_2 / (double)count.size() << " larger than 2 " << (double)count_2_larger / (double)count.size() << std::endl;
	//////std::cout << count_2 << " " << count_2_larger << " " << count.size() << std::endl;
	//////std::cout << "Collider Repeatability 2: " << (double)count_2_collider / (double)count_collider.size() << " " << (double)count_2_larger_collider / (double)count_collider.size() << std::endl;


}

void Collision::findInSP()
{
	unsigned int* vertex_tri_pair = spatial_hashing.vertex_triangle_pair[0]+1;
	unsigned int vertex_num;
	unsigned int num;
	unsigned int real_vertex_index;
	unsigned int* record_triangle;
	for (int i = 0; i < total_obj_num; ++i) {
			vertex_num = vertex_index_start_per_thread[i][thread_num];
		for (int j = 0; j < vertex_num; ++j) {
			num = spatial_hashing.vertex_triangle_pair_num_record[i][j];
			record_triangle = spatial_hashing.vertex_triangle_pair_by_vertex[i] + j * estimate_coeff_for_vt_pair_num;
			if (i < cloth->size()) {
				real_vertex_index = j;
			}
			else {
				real_vertex_index = vertex_index_on_surface[i][j];
			}

			for (int k = 0; k < num; k += 2) {
				if (vertex_tri_pair[0] == real_vertex_index && vertex_tri_pair[1] == i && vertex_tri_pair[2] == record_triangle[k + 1]
					&& vertex_tri_pair[3] == record_triangle[k]) {
					vertex_tri_pair += 4;
				}
				else {
					////std::cout << "vertex triangle not consistent error " << std::endl;
				}
			}

		}
	}
	////std::cout << "finished testing " << std::endl;
}

void Collision::testIfSPRight()
{
	unsigned int* vertex_tri_pair;
	for (int i = 0; i < thread_num; ++i) {
		vertex_tri_pair = spatial_hashing.vertex_triangle_pair[i] + 1;
		for (int j = 0; j < spatial_hashing.vertex_triangle_pair[i][0]; j += 4) {
			if (vertex_tri_pair[j + 1] < cloth->size()) {
				findVertexTriangleInNewSP(vertex_tri_pair[j + 1], vertex_tri_pair[j], vertex_tri_pair[j], vertex_tri_pair[j + 3], vertex_tri_pair[j + 2]);
			}
			else {
				findVertexTriangleInNewSP(vertex_tri_pair[j + 1], tetrahedron->data()[vertex_tri_pair[j + 1]-cloth->size()].mesh_struct.vertex_surface_index[vertex_tri_pair[j]], vertex_tri_pair[j], vertex_tri_pair[j + 3], vertex_tri_pair[j + 2]);
			}
		}
	}
	unsigned int* edge_edge_pair;
	for (int i = 0; i < thread_num; ++i) {
		edge_edge_pair = spatial_hashing.edge_edge_pair[i] + 1;
		for (int j = 0; j < spatial_hashing.edge_edge_pair[i][0]; j += 4) {
			findEdgeEdgeInNewSP(edge_edge_pair[j + 1], edge_edge_pair[j], edge_edge_pair[j + 3], edge_edge_pair[j + 2]);
			//findEdgeEdgeInBVH(edge_edge_pair[j + 1], edge_edge_pair[j], edge_edge_pair[j + 3], edge_edge_pair[j + 2]);
		}
	}
	findInSP();

	//std::vector<unsigned int>* tri_tri_pair = spatial_hashing.triangle_pair;
	//std::vector<unsigned int>* tri_tri_pair_collider = spatial_hashing.triangle_pair_with_collider;
	//for (int i = 0; i < thread_num; ++i) {
	//	for (int j = 0; j < tri_tri_pair[i].size(); j+=4) {
	//		findInBVH(tri_tri_pair[i][j + 1], tri_tri_pair[i][j], tri_tri_pair[i][j + 3], tri_tri_pair[i][j + 2],false);
	//	}
	//	for (int j = 0; j < tri_tri_pair_collider[i].size(); j += 4) {
	//		findInBVH(tri_tri_pair_collider[i][j + 1], tri_tri_pair_collider[i][j], tri_tri_pair_collider[i][j + 3], tri_tri_pair_collider[i][j + 2], true);
	//	}
	//}
	//////std::cout << "find pair in BVH " << std::endl;
	//std::vector<std::vector<std::vector<unsigned int>>>* tri_obj_tri;
	//for (int i = 0; i < total_obj_num; ++i) {
	//	if (i < cloth->size()) {
	//		tri_obj_tri = &cloth->data()[i].triangle_neighbor_obj_triangle;
	//	}
	//	else {
	//		tri_obj_tri = &tetrahedron->data()[i-cloth->size()].triangle_neighbor_obj_triangle;
	//	}
	//	findInSP(tri_obj_tri, i);
	//	//if (i < cloth->size()) {
	//	//	tri_obj_tri = &cloth->data()[i].triangle_neighbor_collider_triangle;
	//	//}
	//	//else {
	//	//	tri_obj_tri = &tetrahedron->data()[i - cloth->size()].triangle_neighbor_collider_triangle;
	//	//}		
	//	//findInSPCollider(tri_obj_tri, i);
	//}
	//for (int i = 0; i < thread_num; ++i) {
	//	for (int j = 0; j < spatial_hashing.triangle_pair[i].size(); j+=4) {
	//		if (spatial_hashing.triangle_pair[i][j] == 1 || spatial_hashing.triangle_pair[i][j+2]==1) {
	//			////std::cout <<i<<" "<< spatial_hashing.triangle_pair[i][j + 1] << " " << spatial_hashing.triangle_pair[i][j] << " "
	//				<< spatial_hashing.triangle_pair[i][j + 3] << " " << spatial_hashing.triangle_pair[i][j + 2] << std::endl;
	//		}
	//	}
	//}
	//////std::cout << cloth->data()[0].triangle_neighbor_obj_triangle[0][0].size() << std::endl;
	//for (int i = 0; i < cloth->data()[0].triangle_neighbor_obj_triangle[0][0].size(); ++i) {
	//	////std::cout << cloth->data()[0].triangle_neighbor_obj_triangle[0][0][i]<<" " << std::endl;
	//}
	////////std::cout <<(int) AABB::AABB_intersection(obj_tri_aabb[0][0].data(), obj_tri_aabb[0][6].data()) << std::endl;
	//////std::cout << "test is right" << std::endl;


}


void Collision::findInSPCollider(std::vector<std::vector<std::vector<unsigned int>>>* tri_obj, unsigned int obj_index)
{
	//std::vector<unsigned int>* triangle_pair = spatial_hashing.triangle_pair_with_collider;
	//bool need_break;
	//bool found_one;
	//for (int i = 0; i < tri_obj->size(); ++i) {
	//	//
	//	for (int j = 0; j < tri_obj->data()[i].size(); ++j) {
	//		found_one = true;
	//		for (int k = 0; k < tri_obj->data()[i][j].size(); ++k) {
	//			//obj_index, i, j, tri_obj->data()[i][j][k]
	//			need_break = false;
	//			for (int l = 0; l < thread_num; ++l) {
	//				for (int m = 0; m < triangle_pair[l].size(); m += 4) {
	//					if (obj_index == triangle_pair[l][m + 1] && i == triangle_pair[l][m] &&
	//						j == triangle_pair[l][m + 3] && tri_obj->data()[i][j][k] == triangle_pair[l][m + 2]) {
	//						need_break = true;
	//						break;
	//					}
	//				}
	//				if (need_break) {
	//					break;
	//				}
	//			}
	//			if (!need_break) {
	//				////std::cout << "does not find the collider pair: " << obj_index << " " << i << " " << j << " " << tri_obj->data()[i][j][k] << std::endl;
	//			}
	//		}
	//	}
	//}
}



void Collision::findInSP(std::vector<std::vector<std::vector<unsigned int>>>* tri_obj, unsigned int obj_index)
{
	//std::vector<unsigned int>* triangle_pair = spatial_hashing.triangle_pair;
	//bool need_break;
	//bool found_one;
	//for (int i = 0; i < tri_obj->size(); ++i) {
	//	//
	//	for (int j = 0; j < tri_obj->data()[i].size(); ++j) {
	//		found_one = true;
	//		for (int k = 0; k < tri_obj->data()[i][j].size(); ++k) {
	//			//obj_index, i, j, tri_obj->data()[i][j][k]
	//			need_break = false;
	//			for (int l = 0; l < thread_num; ++l) {
	//				for (int m = 0; m < triangle_pair[l].size(); m += 4) {
	//					if (obj_index == triangle_pair[l][m + 1] && i == triangle_pair[l][m] &&
	//						j == triangle_pair[l][m + 3] && tri_obj->data()[i][j][k] == triangle_pair[l][m + 2]) {
	//						need_break = true;
	//						break;
	//					}
	//					if (obj_index == triangle_pair[l][m + 3] && i == triangle_pair[l][m + 2] &&
	//						j == triangle_pair[l][m + 1] && tri_obj->data()[i][j][k] == triangle_pair[l][m]) {
	//						need_break = true;
	//						break;
	//					}
	//				}
	//				if (need_break) {
	//					break;
	//				}
	//			}
	//			if (!need_break) {
	//				////std::cout << "does not find the pair: " << obj_index << " " << i << " " << j << " " << tri_obj->data()[i][j][k] << std::endl;
	//			}
	//		}
	//	}
	//}
}

void Collision::findEdgeEdgeInNewSP(unsigned int obj_0, unsigned int edge_index_0, unsigned int obj_1, unsigned int edge_index_1)
{
	unsigned int* neighbor;
	bool find = false;

	unsigned int num = spatial_hashing.edge_edge_pair_num_record[obj_0][edge_index_0];
	neighbor = spatial_hashing.edge_edge_pair_by_edge[obj_0] + edge_index_0 * estimate_coeff_for_ee_pair_num;

	for (int i = 0; i < num; i += 2) {
		if (neighbor[i] == obj_1 && neighbor[i + 1] == edge_index_1) {
			find = true;
		}
	}
	if (!find) {
		////std::cout << "does not find edge edge pair: " << obj_0 << " " << edge_index_0 << " " << obj_1 << " " << edge_index_1 << std::endl;
	}
}

void Collision::findEdgeEdgeInBVH(unsigned int obj_0, unsigned int edge_index_0, unsigned int obj_1, unsigned int edge_index_1)
{
	std::vector<int>* neighbor;
	bool find = false;
	if (obj_0 < cloth->size()) {
		neighbor = &cloth->data()[obj_0].edge_neighbor_obj_edge[edge_index_0][obj_1];
	}
	else {
		neighbor = &tetrahedron->data()[obj_0 - cloth->size()].edge_neighbor_obj_edge[edge_index_0][obj_1];
	}
	for (int i = 0; i < neighbor->size(); ++i) {
		if (neighbor->data()[i] == edge_index_1) {
			find = true;
		}
	}
	if (!find) {
		////std::cout << "does not find edge edge pair: " << obj_0 << " " << edge_index_0 << " " << obj_1 << " " << edge_index_1 << std::endl;
	}
}


void Collision::findVertexTriangleInNewSP(unsigned int obj_0, unsigned int vertex_index_on_surface, unsigned int vertex_index_0, unsigned int obj_1, unsigned int tri_index1)
{
	unsigned int* neighbor;
	bool find = false;
	unsigned int num = spatial_hashing.vertex_triangle_pair_num_record[obj_0][vertex_index_on_surface];
	neighbor = spatial_hashing.vertex_triangle_pair_by_vertex[obj_0] + vertex_index_on_surface * estimate_coeff_for_vt_pair_num;
	
	for (int i = 0; i < num; i+=2) {
		if (neighbor[i]== obj_1 &&  neighbor[i+1] == tri_index1) {
			find = true;
		}
	}
	if (!find) {
		////std::cout << "does not find vertex triangle pair: " << obj_0 << " " << vertex_index_0 << " " << obj_1 << " " << tri_index1 << std::endl;
	}
}


void Collision::findVertexTriangleInBVH(unsigned int obj_0, unsigned int vertex_index_0, unsigned int obj_1, unsigned int tri_index1)
{
	std::vector<int>* neighbor;
	bool find = false;
	if (obj_0 < cloth->size()) {
		neighbor = &cloth->data()[obj_0].vertex_neighbor_obj_triangle[vertex_index_0][obj_1];
	}
	else {
		//////std::cout << vertex_index_0 << " " << tetrahedron->data()[0].mesh_struct.vertex_surface_index.size() << std::endl;
		int vertex_index_in_surface = tetrahedron->data()[obj_0 - cloth->size()].mesh_struct.vertex_surface_index[vertex_index_0];
		//////std::cout << vertex_index_in_surface << std::endl;
		neighbor = &tetrahedron->data()[obj_0 - cloth->size()].surface_vertex_neighbor_obj_triangle[vertex_index_in_surface][obj_1];
	}
	for (int i = 0; i < neighbor->size(); ++i) {
		if (neighbor->data()[i] == tri_index1) {
			find = true;
		}
	}
	if (!find) {
		////std::cout << "does not find vertex triangle pair: " << obj_0 << " " << vertex_index_0 << " " << obj_1 << " " << tri_index1 << std::endl;
	}
}


//void Collision::findVertexTriangleInSH()
//{
//
//}


void Collision::findInBVH(unsigned int obj_0, unsigned int tri_index_0, unsigned int obj_1, unsigned int tri_index1, bool with_collider)
{
	std::vector<unsigned int>* neighbor;
	bool find = false;
	if (!with_collider) {
		if (obj_0 < cloth->size()) {
			neighbor = &cloth->data()[obj_0].triangle_neighbor_obj_triangle[tri_index_0][obj_1];
		}
		else {
			neighbor = &tetrahedron->data()[obj_0 - cloth->size()].triangle_neighbor_obj_triangle[tri_index_0][obj_1];
		}
		for (int i = 0; i < neighbor->size(); ++i) {
			if (neighbor->data()[i] == tri_index1) {
				find = true;
			}
		}
		if (obj_1 < cloth->size()) {
			neighbor = &cloth->data()[obj_1].triangle_neighbor_obj_triangle[tri_index1][obj_0];
		}
		else {
			neighbor = &tetrahedron->data()[obj_1 - cloth->size()].triangle_neighbor_obj_triangle[tri_index1][obj_0];
		}
		for (int i = 0; i < neighbor->size(); ++i) {
			if (neighbor->data()[i] == tri_index_0) {
				find = true;
			}
		}
		if (!find) {
			////std::cout << "does not find the pair: " << obj_0 << " " << tri_index_0 << " " << obj_1 << " " << tri_index1 << std::endl;
		}
	}
	else {
		if (obj_0 < cloth->size()) {
			neighbor = &cloth->data()[obj_0].triangle_neighbor_collider_triangle[tri_index_0][obj_1];
		}
		else {
			neighbor = &tetrahedron->data()[obj_0 - cloth->size()].triangle_neighbor_collider_triangle[tri_index_0][obj_1];
		}
		for (int i = 0; i < neighbor->size(); ++i) {
			if (neighbor->data()[i] == tri_index1) {
				find = true;
			}
		}
		if (!find) {
			////std::cout << "does not find the collider pair: " << obj_0 << " " << tri_index_0 << " " << obj_1 << " " << tri_index1 << std::endl;
		}
	}

}

void Collision::testCulling()
{
	for (int i = 0; i < (*collider)[0].triangle_AABB.size(); ++i) {
		for (int j = 0; j < (*collider)[0].triangle_neighbor_obj_vertex[i][0].size(); ++j) {
			if ((*collider)[0].triangle_neighbor_obj_vertex[i][0][j] == 13) {
				////std::cout << i << std::endl;
			}
		}
	}
}


void Collision::buildSpatialHashingForOri()
{
	getAABB();
	getSceneAABB();
	spatial_hashing.buildSpatialHashingForOri(scene_aabb);
 }


void Collision::collisionCulling()
{
	time_t t = clock();
	time_t t1 = clock();
//	t = clock();
//	for (int i = 0; i < 100; ++i) {
	getAABB();
	//////std::cout << "end here " << std::endl;
	getSceneAABB();
	//////std::cout << "end here 2 " << std::endl;
//	}
	//t1 = clock();
	//////std::cout << "AABB " << t1 - t << std::endl;
	//
	//////std::cout <<"aabb "<< scene_aabb[0] << " " << scene_aabb[1] << " " << scene_aabb[2] << " " << scene_aabb[3] << " " << scene_aabb[4] << " " << scene_aabb[5] << std::endl;

	spatial_hashing.buildSpatialHashing(scene_aabb);

	//unsigned int num = 0;
	//////std::cout << "pair num " << std::endl;
	//for (int j = 0; j < thread_num; ++j) {
	//	//////std::cout << spatial_hashing.vertex_triangle_pair[j][0] / 4  << std::endl;
	//	num += spatial_hashing.vertex_triangle_pair[j][0] / 4;
	//}
	//////std::cout << "total v-t pair " << num << std::endl;
	//for (int j = 0; j < thread_num; ++j) {
	//	for (int i = 1; i < spatial_hashing.vertex_triangle_pair[j][0] + 1; i += 4) {
	//		//if (spatial_hashing.vertex_triangle_pair[j][i + 1] != spatial_hashing.vertex_triangle_pair[j][i + 3]) {
	//			////std::cout << spatial_hashing.vertex_triangle_pair[j][i] << " " << spatial_hashing.vertex_triangle_pair[j][i + 1] << " "
	//				<< spatial_hashing.vertex_triangle_pair[j][i + 2] << " " << spatial_hashing.vertex_triangle_pair[j][i + 3] << std::endl;
	//		//}
	//	}		
	//}


	//num = 0;
	//for (int j = 0; j < thread_num; ++j) {
	//	////std::cout << spatial_hashing.edge_edge_pair[j][0] / 4 << std::endl;
	//	num += spatial_hashing.edge_edge_pair[j][0] / 4;
	//}
	//////std::cout << "total e-e pair " << num << std::endl;

	//for (int j = 0; j < thread_num; ++j) {
	//	spatial_hashing.vertex_triangle_pair[j][0] = 0;
	//	spatial_hashing.edge_edge_pair[j][0] = 0;
	//	spatial_hashing.vertex_obj_triangle_collider_pair[j][0] = 0;
	//}
	//spatial_hashing.testColliderPair();




	if (!record_pair_by_element) {
		setPairIndexEveryThread();
	}
	//if (use_BVH) {
	//	buildBVH();
	//	thread->assignTask(this, FIND_TRIANGLE_PAIRS);
	//	thread->assignTask(this, FIND_PRIMITIVE_AROUND);
		//testIfSPRight();
	//	testRepeatability();
	//}
	


	//t = clock();	
	//for (int i = 0; i < 1000; ++i) {
	//for (int j = 0; j < thread_num; ++j) {
	//		findAllVertexVertexEdgePairs(j);
	//}		
	//}
	//t1 = clock();
	//////std::cout << "find all vertex vertex edge pairs single thread " <<  t1 - t << std::endl;	
	//t = clock();
	//for (int i = 0; i < 1000; ++i) {



	//thread->assignTask(this, FIND_VERTEX_VERTEX_VERTEX_EDGE_PAIRS);


	//}
	//t1 = clock();
	//////std::cout << "find all vertex vertex edge pairs multi thread "<<t1-t << std::endl;
	//num = 0;
	//unsigned int num1 = 0;
	//for (int j = 0; j < thread_num; ++j) {
	//	////std::cout << vertex_vertex_pair[j][0]/4 << " " << vertex_edge_pair[j][0]/4 << std::endl;
	//	num += vertex_vertex_pair[j][0] / 4;
	//	num1 += vertex_edge_pair[j][0] / 4;
	//}
	//////std::cout << "total v-v pair " << num <<" total v-e pair "<<num1 << std::endl;


	//setVertexVertexEdgePairIndexEveryThread();





	//
	//t = clock();
	//for (int i = 0; i < 10; ++i) {
	//
	//}
	//t1 = clock();
	//////std::cout << "find primitive around multi thread" << t1 - t << std::endl;
	//t = clock();
	//for (int i = 0; i < 10; ++i) {
	//	for (int j = 0; j < thread_num; ++j) {
	//		findPointTriangleEdgeEdgePair(j);
	//	}	
	//}
	//t1 = clock();
	//////std::cout << "find primitive around single thread" << t1 - t << std::endl;
}





void Collision::getSceneAABB()
{
	memset(scene_aabb + 3, 0xFE, 24); //set double to -5.31401e+303
	memset(scene_aabb, 0x7F, 24); //set double to 1.38242e+306
	double* aabb_;
	for (int i = 0; i < total_obj_with_collider; ++i) {
		if (i < tetrahedron_begin_obj_index) {
			aabb_ = cloth->data()[i].obj_aabb;
		}
		else if (i < total_obj_num) {
			aabb_ = tetrahedron->data()[i - tetrahedron_begin_obj_index].obj_aabb;
		}
		else {
			aabb_ = collider->data()[i - total_obj_num].obj_aabb;
		}
		for (int j = 0; j < 3; ++j) {
			if (scene_aabb[j] > aabb_[j]) {
				scene_aabb[j] = aabb_[j];
			}
		}
		for (int j = 3; j < 6; ++j) {
			if (scene_aabb[j] < aabb_[j]) {
				scene_aabb[j] = aabb_[j];
			}
		}
	}
	for (int i = 0; i < 3; ++i) {
		scene_aabb[i + 3] += 0.1;
		scene_aabb[i] -= 0.1;
	}
	//////std::cout << scene_aabb[1] << std::endl;
	//for (int i = 0; i < 3; ++i) {
	//	scene_aabb[i + 3] = scene_aabb[i] + spatial_hashing.cell_length * (double)((unsigned int)floor((scene_aabb[i + 3] - scene_aabb[i]) / cell_length) + 1);
	//}
	//////std::cout << scene_aabb[4] << std::endl;
}


void Collision::setPairIndexEveryThread()
{
	setPairIndexEveryThread(spatial_hashing.vertex_triangle_pair, vertex_triangle_pair_index_start_per_thread, ave_pair_num[0]);
	setPairIndexEveryThread(spatial_hashing.edge_edge_pair, edge_edge_pair_index_start_per_thread, ave_pair_num[1]);
	if (has_collider) {
		setPairIndexEveryThread(spatial_hashing.vertex_obj_triangle_collider_pair, vertex_obj_triangle_collider_pair_index_start_per_thread, ave_pair_num[2]);
		setPairIndexEveryThread(spatial_hashing.vertex_collider_triangle_obj_pair, vertex_collider_triangle_obj_pair_index_start_per_thread, ave_pair_num[3]);
		setPairIndexEveryThread(spatial_hashing.edge_edge_pair_collider, edge_edge_pair_collider_index_start_per_thread, ave_pair_num[4]);

		//for (int i = 0; i < thread_num; ++i) {
		//	////std::cout << edge_edge_pair_collider_index_start_per_thread[i] << " ";
		//}
		//////std::cout << std::endl;
		//for (int i = 0; i < thread_num+1; ++i) {
		//	////std::cout << vertex_obj_triangle_collider_pair_index_start_per_thread[i << 1] << " " << vertex_obj_triangle_collider_pair_index_start_per_thread[i * 2 + 1] << std::endl;
		//}


	}

	

}


void Collision::setVertexVertexEdgePairIndexEveryThread()
{
	unsigned int pair_;
	setPairIndexEveryThread(vertex_edge_pair, vertex_edge_pair_index_start_per_thread, pair_);
	setPairIndexEveryThread(vertex_vertex_pair, vertex_vertex_pair_index_start_per_thread, pair_);
	if (has_collider) {
		setPairIndexEveryThread(vertex_obj_edge_collider_pair, vertex_obj_edge_collider_pair_index_start_per_thread, pair_);
		setPairIndexEveryThread(vertex_collider_edge_obj_pair, vertex_collider_edge_obj_pair_index_start_per_thread, pair_);
		setPairIndexEveryThread(vertex_vertex_pair_collider, vertex_vertex_pair_collider_index_start_per_thread, pair_);
	}
}

void Collision::setTargetPositionEven()
{
	//setPairIndexEveryThread(target_position_index.data(), target_position_element_start_per_thread);
	//memset(point_triangle_target_position_element_start_per_thread.data(), 0, 8 * (thread_num + 1));
	setPairIndexEveryThread(point_triangle_target_pos_index.data(), point_triangle_target_position_element_start_per_thread);
	setPairIndexEveryThread(edge_edge_target_pos_index.data(), edge_edge_target_position_element_start_per_thread);
	if (has_collider) {
		//setPairIndexEveryThread(point_collider_triangle_target_pos_index.data(), point_collider_triangle_target_position_element_start_per_thread);
		//setPairIndexEveryThread(edge_edge_collider_target_pos_index.data(), edge_edge_collider_target_position_element_start_per_thread);
		setPairIndexEveryThread(point_triangle_collider_target_pos_index.data(), point_triangle_collider_target_position_element_start_per_thread);
	}
	

}


//make pair num in every thread even
void Collision::setIndexEveryThread(std::vector<unsigned int>* pair, std::vector<unsigned int>& pair_index_start_per_thread)
{
	std::vector<unsigned int>pair_start(thread_num);
	unsigned int pair_num = pair[0][0];
	for (int i = 1; i < thread_num; ++i) {
		pair_num += pair[i][0];
	}
	countInEveryThread(thread_num, pair_num, pair_start.data());

	unsigned int num;
	for (int i = 0; i < thread_num; ++i) {
		num = pair_start[i];
		if (num <= (pair[pair_index_start_per_thread[i << 1]][0] - pair_index_start_per_thread[(i << 1) + 1]))
		{
			pair_index_start_per_thread[(i + 1) << 1] = pair_index_start_per_thread[i << 1];
			pair_index_start_per_thread[(i << 1) + 3] = pair_index_start_per_thread[(i << 1) + 1] + num;
		}
		else {
			num -= (pair[pair_index_start_per_thread[i << 1]][0] - pair_index_start_per_thread[(i << 1) + 1]);
			for (int j = pair_index_start_per_thread[i << 1] + 1; j < thread_num; ++j) {
				if (num <= (pair[j][0])) {
					pair_index_start_per_thread[(i + 1) << 1] = j;
					pair_index_start_per_thread[(i << 1) + 3] = num;
					break;
				}
				else {
					num -= pair[j][0];
				}
			}
		}
	}
}




void Collision::testPairEven()
{
	std::vector<std::vector<unsigned int>> pair(thread_num);
	for (int i = 0; i < thread_num; ++i) {
		pair[i].resize(3,0);
	}
	for (int i = 0; i < 4; ++i) {
		pair[i][0] = i << 2;
	}
	std::vector<unsigned int> pair_index_start_per_thread((thread_num+1) << 1, 0);
	setPairIndexEveryThread(pair.data(), pair_index_start_per_thread);
	for (int i = 0; i < thread_num; ++i) {
		////std::cout << pair[i][0] << " ";
	}
	////std::cout << std::endl;
	for (int i = 0; i < thread_num+1; ++i) {
		////std::cout << pair_index_start_per_thread[i << 1] << " " << pair_index_start_per_thread[(i << 1) + 1]<<std::endl;
	}
	for (int i = 0; i < thread_num; ++i) {
		pair[i][0]=0;
	}
	for (int i = 4; i < thread_num; ++i) {
		pair[i][0] = i << 2;
	}
	setPairIndexEveryThread(pair.data(), pair_index_start_per_thread);
	for (int i = 0; i < thread_num; ++i) {
		////std::cout << pair[i][0] << " ";
	}
	////std::cout << std::endl;
	for (int i = 0; i < thread_num + 1; ++i) {
		////std::cout << pair_index_start_per_thread[i << 1] << " " << pair_index_start_per_thread[(i << 1) + 1] << std::endl;
	}
	for (int i = 0; i < thread_num; ++i) {
		pair[i][0] = 0;
	}
	for (int i = 0; i < thread_num; ++i) {
		pair[i][0] = 4*(5+rand()%5);
	}
	pair[thread_num - 1][0] = 4;
	for (int i = 2; i < 6; ++i) {
		pair[i][0] = 0;
	}

	setPairIndexEveryThread(pair.data(), pair_index_start_per_thread);
	for (int i = 0; i < thread_num; ++i) {
		////std::cout << pair[i][0] << " ";
	}
	////std::cout << std::endl;
	for (int i = 0; i < thread_num + 1; ++i) {
		////std::cout << pair_index_start_per_thread[i << 1] << " " << pair_index_start_per_thread[(i << 1) + 1] << std::endl;
	}
}


//make pair num in every thread even
void Collision::setPairIndexEveryThread(std::vector<unsigned int>* pair, std::vector<unsigned int>& pair_index_start_per_thread)
{
	std::vector<unsigned int>pair_start(thread_num);
	unsigned int pair_num = pair[0][0] >> 2;
	for (int i = 1; i < thread_num; ++i) {
		pair_num += pair[i][0] >> 2;
	}
	countInEveryThread(thread_num, pair_num, pair_start.data());
	//unsigned int total_num;

	unsigned int num;
	for (int i = 0; i < thread_num; ++i) {
		num = pair_start[i];
		//if (pair_num > 0) {
		//	//std::cout<<"pair "<< num << " " << pair[pair_index_start_per_thread[i << 1]][0] - pair_index_start_per_thread[(i << 1) + 1] << std::endl;
		//}

		if (num <= (pair[pair_index_start_per_thread[i << 1]][0] - pair_index_start_per_thread[(i << 1) + 1]) >> 2)
		{
			pair_index_start_per_thread[(i + 1) << 1] = pair_index_start_per_thread[i << 1];
			pair_index_start_per_thread[(i << 1) + 3] = pair_index_start_per_thread[(i << 1) + 1] + (num << 2);
		}
		else {
			num -= (pair[pair_index_start_per_thread[i << 1]][0] - pair_index_start_per_thread[(i << 1) + 1]) >> 2;
			for (int j = pair_index_start_per_thread[i << 1] + 1; j < thread_num; ++j) {
				if (num <= (pair[j][0] >> 2)) {
					pair_index_start_per_thread[(i + 1) << 1] = j;
					pair_index_start_per_thread[(i << 1) + 3] = num << 2;
					break;
				}
				else {
					num -= pair[j][0] >> 2;
				}
			}
		}

	}
}


//make pair num in every thread even
void Collision::setPairIndexEveryThread(unsigned int** pair, std::vector<unsigned int>&pair_index_start_per_thread, unsigned int& ave_pair_num)
{
	std::vector<unsigned int>pair_start(thread_num);
	unsigned int pair_num = pair[0][0]>>2;
	for (int i = 1; i < thread_num; ++i) {
		pair_num += pair[i][0]>>2;
	}
	countInEveryThread(thread_num, pair_num, pair_start.data());
	ave_pair_num = pair_start[0];

	unsigned int num;
	for (int i = 0; i < thread_num; ++i) {
		num = pair_start[i];
		if (num <= (pair[pair_index_start_per_thread[i << 1]][0]- pair_index_start_per_thread[(i << 1) + 1]) >> 2)
		{
			pair_index_start_per_thread[(i + 1) << 1] = pair_index_start_per_thread[i << 1];
			pair_index_start_per_thread[(i << 1) + 3] = pair_index_start_per_thread[(i << 1) + 1] + (num<<2);
		}
		else {
			num -= (pair[pair_index_start_per_thread[i << 1]][0] - pair_index_start_per_thread[(i << 1) + 1]) >> 2;
			for (int j = pair_index_start_per_thread[i << 1] + 1; j < thread_num; ++j) {
				if (num <= (pair[j][0] >> 2)) {
					pair_index_start_per_thread[(i + 1) << 1] = j;
					pair_index_start_per_thread[(i << 1) + 3] = num<<2;
					break;
				}
				else {
					num -= pair[j][0] >> 2;
				}
			}
		}
	}
	//////std::cout << "pair num " << std::endl;
	//for (int i = 0; i < thread_num; ++i) {
	//	////std::cout << pair_start[i] << " ";
	//}
	//////std::cout << std::endl;
	//for (int i = 0; i < thread_num; ++i) {
	//	////std::cout << pair[i][0] / 4<<" ";
	//}
	//////std::cout << std::endl;
	//for (int i = 0; i < thread_num; ++i) {
	//	////std::cout << pair_index_start_per_thread[i * 2 + 2] << " " << pair_index_start_per_thread[i * 2 + 3]/4 << " ";
	//}
	//////std::cout << std::endl;

}

//recrod all pair with collision, compute collision_time
void Collision::collisionTimeWithPair()
{
	thread->assignTask(this, GLOBAL_COLLISION_TIME_ADD_PAIR);
	collision_time = collision_time_thread[0];
	for (int i = 1; i < thread_num; ++i) {
		if (collision_time > collision_time_thread[i]) {
			collision_time = collision_time_thread[i];
		}
	}
	if (collision_time > 1.0) {
		collision_time = 1.0;
	}
	if (collision_time < 1.0) {
		collision_time *= 0.9;
	}


	if (collision_time == 0.0) {
		std::cout << "attention: collision time equals zero" << std::endl;
	}

	if (collision_time != 1.0) {
		std::cout << "--collision time " << collision_time << std::endl;
	}
}


void Collision::globalCollisionTime()
{
	//for (int i = 0; i < cloth->size(); ++i) {
	//	thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	//}
	//for (int i = 0; i < tetrahedron->size(); ++i) {
	//	thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL);
	//}
	//time_t t, t1;
	//updateEigenPosition();	
	//t = clock();
	//for (int j = 0; j < 10; ++j) {
	//	for (int i = 0; i < thread_num; ++i) {		
	//		collisionTime(i);
	//	}		
	//}
	//t1 = clock();
	//////std::cout << "ccd collision detection time single thread " << t1 - t << std::endl;
	//t = clock();
	//for (int j = 0; j < 10; ++j) {
	//	for (int i = 0; i < thread_num; ++i) {
	//		collisionTimeCompare(i);
	//	}
	//}
	//t1 = clock();
	//////std::cout << "ccd collision detection time IPC ori single thread " << t1 - t << std::endl;
	//t = clock();
	//for (int j = 0; j < 10; ++j) {
	//if (!record_pair_by_element) {
	//	if (use_method == PD_) {
	//		initialCollisionTimeRecord();
	//	}
	//}

	thread->assignTask(this, GLOBAL_COLLISION_TIME);
	//}
	//t1 = clock();
	//////std::cout << "ccd collision detection time multi thread " << t1 - t << std::endl;
	collision_time = collision_time_thread[0];
	for (int i = 1; i < thread_num; ++i) {
		if (collision_time > collision_time_thread[i]) {
			collision_time = collision_time_thread[i];
		}
	}
	if (collision_time > 1.0) {
		collision_time = 1.0;
	}
	if (collision_time < 1.0) {
		collision_time *= 0.9;
	}


	if (collision_time == 0.0) {
		std::cout << "attention: collision time equals zero" << std::endl;
	}

	if (collision_time != 1.0) {
		std::cout << "--collision time " << collision_time << std::endl;
	}
}



void Collision::closePairCollisionTime()
{	
	thread->assignTask(this, CLOSE_PAIR_COLLISION_TIME);
	collision_time = collision_time_thread[0];
	for (int i = 1; i < thread_num; ++i) {
		if (collision_time > collision_time_thread[i]) {
			collision_time = collision_time_thread[i];
		}
	}
	if (collision_time > 1.0) {
		collision_time = 1.0;
	}
	if (collision_time < 1.0) {
		collision_time *= 0.9;
	}
	if (collision_time == 0.0) {
		std::cout << "attention: collision time equals zero, color num inner itr "<<*inner_iteration_number << std::endl;
	}
	if (collision_time != 1.0) {
		std::cout << "last color collision time " << collision_time << " " << std::endl;
	}


	thread->assignTask(this, COLLISION_FREE_POSITION_LAST_COLOR);

	//double dist2_cur = CCD::internal::pointTriangleDistanceUnclassified(vertex_position[0][4048].data(), vertex_position[0][triangle_indices[0][7842][0]].data(),
	//	vertex_position[0][triangle_indices[0][7842][1]].data(), vertex_position[0][triangle_indices[0][7842][2]].data());
	//std::cout << "distance of the chosen VT " << dist2_cur << " " << tolerance * tolerance << std::endl;
}


//COLLISION_FREE_POSITION_LAST_COLOR
void Collision::computeCollisionFreePositionForColor(int thread_No)
{
	unsigned int index_end;

	double collision_time = this->collision_time;
	std::array<double, 3>* q_pre;
	std::array<double, 3>* q_end;


	int* record_position_num;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		index_end = all_vertex_index_start_per_thread[i][thread_No + 1];
		q_end = vertex_position[i];
		q_pre = vertex_record_for_this_color[i];
		record_position_num = indicate_if_involved_in_last_color[i][0].data();
		for (unsigned int j = all_vertex_index_start_per_thread[i][thread_No]; j < index_end; ++j) {
			if (record_position_num[j]) {
				q_end[j][0] = q_pre[j][0] + collision_time * (q_end[j][0] - q_pre[j][0]);
				q_end[j][1] = q_pre[j][1] + collision_time * (q_end[j][1] - q_pre[j][1]);
				q_end[j][2] = q_pre[j][2] + collision_time * (q_end[j][2] - q_pre[j][2]);
				memcpy(q_pre[j].data(), q_end[j].data(), 24);
			}
		}
	}
}


void Collision::floorCollisionTime(int color)
{
	if (!floor->exist) {
		return;
	}


	thread->assignTask(this, FLOOR_COLLISION_TIME,color);
	collision_time = collision_time_thread[0];
	for (int i = 1; i < thread_num; ++i) {
		if (collision_time > collision_time_thread[i]) {
			collision_time = collision_time_thread[i];
		}
	}
	if (collision_time > 1.0) {
		collision_time = 1.0;
	}
	if (collision_time < 1.0) {
		collision_time *= 0.9;
	}
	if (collision_time == 0.0) {
		std::cout << "attention: floor collision time equals zero, color num " << std::endl;
	}

	thread->assignTask(this, UPDATE_POSITION_FOR_FLOOR_COLLISION, color);
	
}

//FLOOR_COLLISION_TIME
void Collision::floorCollisionTime(int thread_No, int color)
{
	unsigned int* tet_group;
	int color_group_index;
	int start, end;
	int j;

	char* is_tet_involved_in_collision;
	int* tet_vertex;

	std::array<int, 4>* tet_indices;
	unsigned int* unfixed_vertex_num;

	int num;

	int* global_to_normal;

	std::array<double, 3>* initial_pos;
	std::array<double, 3>* current_pos;

	double collision_time = 2.0;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		color_group_index = *inner_iteration_number % tet_color_groups[i]->size();
		end = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color][thread_No + 1];
		start = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color][thread_No];

		tet_group = tet_color_groups[i]->data()[color_group_index][color].data();

		is_tet_involved_in_collision = tet_color_groups_label[i][color_group_index][color].data();

		tet_indices = tetrahedron->data()[i].mesh_struct.unfixied_actual_indices.data();
		unfixed_vertex_num = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();
		global_to_normal = general_index_to_surface_index[i + cloth->size()];
		initial_pos = vertex_record_for_this_color[i + cloth->size()];
		current_pos = vertex_position[i + cloth->size()];

		for (int k = start; k < end; ++k) {			
			if (!is_tet_involved_in_collision[k]) {
				j = tet_group[k];
				tet_vertex = tet_indices[j].data();
				num = unfixed_vertex_num[j];
				for (int m = 0; m < num; ++m) {
					if (global_to_normal[tet_vertex[m]] != -1) {
						floorCollisionTime(initial_pos[tet_vertex[m]].data(), current_pos[tet_vertex[m]].data(), floor->dimension,
							floor->normal_direction, floor->value, collision_time, tolerance);
					}
				}
			}
		}
	}

	collision_time_thread[thread_No] = collision_time;
}


//UPDATE_POSITION_FOR_FLOOR_COLLISION
void Collision::updatePositionForFloor(int thread_No, int color)
{
	unsigned int* tet_group;
	int color_group_index;
	int start, end;
	int j;

	char* is_tet_involved_in_collision;
	int* tet_vertex;

	std::array<int, 4>* tet_indices;
	unsigned int* unfixed_vertex_num;

	int num;

	int* global_to_normal;

	std::array<double, 3>* initial_pos;
	std::array<double, 3>* current_pos;

	int vertex_index;
	double collision_time = this->collision_time;

	for (int i = 0; i < tetrahedron->size(); ++i) {
		color_group_index = *inner_iteration_number % tet_color_groups[i]->size();
		end = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color][thread_No + 1];
		start = tetrahedron->data()[i].mesh_struct.tet_in_a_group_start_per_thread_groups[color_group_index][color][thread_No];

		tet_group = tet_color_groups[i]->data()[color_group_index][color].data();

		is_tet_involved_in_collision = tet_color_groups_label[i][color_group_index][color].data();

		tet_indices = tetrahedron->data()[i].mesh_struct.unfixied_actual_indices.data();
		unfixed_vertex_num = tetrahedron->data()[i].mesh_struct.tet_unfixed_vertex_num.data();

		initial_pos = vertex_record_for_this_color[i + cloth->size()];
		current_pos = vertex_position[i + cloth->size()];

		for (int k = start; k < end; ++k) {
			if (!is_tet_involved_in_collision[k]) {
				j = tet_group[k];
				tet_vertex = tet_indices[j].data();
				num = unfixed_vertex_num[j];
				for (int m = 0; m < num; ++m) {
					vertex_index = tet_vertex[m];
					current_pos[vertex_index][0] = initial_pos[vertex_index][0] + collision_time * (current_pos[vertex_index][0] - initial_pos[vertex_index][0]);
					current_pos[vertex_index][1] = initial_pos[vertex_index][1] + collision_time * (current_pos[vertex_index][1] - initial_pos[vertex_index][1]);
					current_pos[vertex_index][2] = initial_pos[vertex_index][2] + collision_time * (current_pos[vertex_index][2] - initial_pos[vertex_index][2]);
					
					memcpy(initial_pos[vertex_index].data(), current_pos[vertex_index].data(), 24);

				}
			}
		}
	}
}


void Collision::collisionTimeColor(int color)
{
	thread->assignTask(this, COLOR_COLLISION_TIME,color);
	collision_time = collision_time_thread[0];
	for (int i = 1; i < thread_num; ++i) {
		if (collision_time > collision_time_thread[i]) {
			collision_time = collision_time_thread[i];
		}
	}
	if (collision_time > 1.0) {
		collision_time = 1.0;
	}
	if (collision_time < 1.0) {
		collision_time *= 0.9;
	}
	if (collision_time == 0.0) {
		std::cout << "attention: color collision time equals zero, color num "<<color << std::endl;
	}
	if (collision_time != 1.0) {
		std::cout << "color collision time " << collision_time << " " << color << std::endl;
	}


	thread->assignTask(this, UPDATE_COLOR_POSITION, color);

	//double dist2_cur = CCD::internal::pointTriangleDistanceUnclassified(vertex_position[0][4048].data(), vertex_position[0][triangle_indices[0][7842][0]].data(),
	//	vertex_position[0][triangle_indices[0][7842][1]].data(), vertex_position[0][triangle_indices[0][7842][2]].data());
	//std::cout <<"distance of the chosen VT "<<  dist2_cur << " " << tolerance * tolerance << std::endl;

}


//UPDATE_COLOR_POSITION
void Collision::updatePositionColor(int thread_No, int color)
{
	int i;
	unsigned int end_per_thread;
	unsigned int* index_of_a_tet_color;

	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* ini_vertex_pos;


	unsigned int* edge_vertices;



	double* edge_length;


	unsigned int* edge_vertex;

	double collision_time = this->collision_time;

	bool has_collider = this->has_collider;


	int edge_index;
	int color_group_index;

	for (int tet_obj_no = 0; tet_obj_no < tetrahedron->size(); ++tet_obj_no) {
		i = tet_obj_no + cloth->size();
		color_group_index = *inner_iteration_number % tet_color_groups[i]->size();
		//if (color >= tet_color_groups[i]->data()[color_group_index].size()-1) {
		//	continue;
		//}
		vertex_pos = vertex_position[i];
		ini_vertex_pos = vertex_record_for_this_color[i];
		index_of_a_tet_color = vertex_index_of_a_tet_color_group[i][color_group_index][color].data();
		end_per_thread = vertex_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color][thread_No + 1];
		for (int j = vertex_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color][thread_No]; j < end_per_thread; ++j) {
			edge_index = index_of_a_tet_color[j];
			vertex_pos[edge_index][0] = ini_vertex_pos[edge_index][0] + collision_time * (vertex_pos[edge_index][0] - ini_vertex_pos[edge_index][0]);
			vertex_pos[edge_index][1] = ini_vertex_pos[edge_index][1] + collision_time * (vertex_pos[edge_index][1] - ini_vertex_pos[edge_index][1]);
			vertex_pos[edge_index][2] = ini_vertex_pos[edge_index][2] + collision_time * (vertex_pos[edge_index][2] - ini_vertex_pos[edge_index][2]);
			memcpy(ini_vertex_pos[edge_index].data(), vertex_pos[edge_index].data(), 24);
		}
	}
}






void Collision::saveCollisionPairVolume()
{
	storeVolume();
	thread->assignTask(this, COMPUTE_VOLUME);
}

void Collision::storeVolume()
{
	//store the start index of volume
	for (int i = 0; i < total_obj_num; ++i) {
		Basic::cumulativeSumDoubleRecord(vertex_triangle_pair_num_record[i], 
			vertex_index_start_per_thread[i][thread_num], VT_start_index[i]);
		Basic::cumulativeSumDoubleRecord(triangle_vertex_pair_num_record[i],
			(unsigned int)mesh_struct[i]->triangle_indices.size(), TV_start_index[i]);
		//Basic::cumulativeSumDoubleRecord(edge_edge_pair_num_record[i],
		//	(unsigned int)(mesh_struct[i]->edge_vertices.size()>>1), EE_start_index[i]);
		if (!collider->empty()) {
			Basic::cumulativeSumDoubleRecord(vertex_obj_triangle_collider_num_record[i],
				vertex_index_start_per_thread[i][thread_num], VT_collider_start_index[i]);
		}
	}

	//resize volume
	for (int i = 0; i < total_obj_num; ++i) {
		VT_volume[i].resize(VT_start_index[i][vertex_index_start_per_thread[i][thread_num]]);

		TV_volume[i].resize(TV_start_index[i][mesh_struct[i]->triangle_indices.size()]);
		EE_volume[i].resize(EE_start_index[i][mesh_struct[i]->edge_vertices.size()/2]);
		if (!collider->empty()) {
			VT_collider_volume[i].resize(VT_collider_start_index[i][vertex_index_start_per_thread[i][thread_num]]);
		}
	}	
}

//COMPUTE_VOLUME
void Collision::computeVolume(int thread_No)
{
	computeVTVolume(thread_No);
	computeTVVolume(thread_No);
	//computeEEVolume(thread_No);
	if (!collider->empty()) {
		computeVTColliderVolume(thread_No);
	}

}

//void Collision::computeEEVolume(int thread_No)
//{
//	unsigned int edge_num;
//	std::array<double, 3>* vertex_pos;
//	unsigned int* pair_num;
//	unsigned int* primitive_record;
//	double* volume;
//	unsigned int* start_index;
//	unsigned int* edge_vertex;
//
//	for (int i = 0; i < total_obj_num; ++i) {
//		edge_num =edge_index_start_per_thread[i][thread_No + 1];
//		vertex_pos = vertex_for_render[i];
//		pair_num =edge_edge_pair_num_record[i];
//		primitive_record = edge_edge_pair_by_edge[i];
//		volume = EE_volume[i].data();
//		start_index = EE_start_index[i];
//		edge_vertex = edge_vertices[i];
//		for (int j = edge_index_start_per_thread[i][thread_No]; j < edge_num; ++j) {
//			computeEEVolume(vertex_pos[edge_vertex[j<<1]].data(), vertex_pos[edge_vertex[(j << 1)+1]].data(),
//				pair_num[j], primitive_record + close_ee_pair_num * j, volume + start_index[j],vertex_for_render, edge_vertices);
//		}		
//	}
//}


void Collision::computeTVVolume(int thread_No)
{
	unsigned int triangle_num;
	std::array<double, 3>* vertex_pos;
	unsigned int* pair_num;
	unsigned int* primitive_record;
	double* volume;
	unsigned int* start_index;
	std::array<int,3>* triangle_vertex;
	for (int i = 0; i < total_obj_num; ++i) {
		triangle_num = triangle_index_start_per_thread[i][thread_No + 1];
		vertex_pos = vertex_for_render[i];
		pair_num = triangle_vertex_pair_num_record[i];
		primitive_record = triangle_vertex_pair_by_triangle[i];
		volume =TV_volume[i].data();
		start_index = TV_start_index[i];
		triangle_vertex = triangle_indices[i];
		for (int j = triangle_index_start_per_thread[i][thread_No]; j < triangle_num; ++j) {		
			computeTVVolume(vertex_pos[triangle_vertex[j][0]].data(), vertex_pos[triangle_vertex[j][1]].data(),
				vertex_pos[triangle_vertex[j][2]].data(),
				pair_num[j], primitive_record + close_tv_pair_num * j, volume + start_index[j], vertex_for_render);
		}
	}
}




void Collision::computeVTColliderVolume(int thread_No)
{
	unsigned int vertex_num;
	std::array<double, 3>* vertex_pos;
	unsigned int* pair_num;
	unsigned int* primitive_record;
	double* volume;
	unsigned int* start_index;
	unsigned int* vertex_index_on_surface_;
	for (int i = 0; i < total_obj_num; ++i) {
		vertex_num = vertex_index_start_per_thread[i][thread_No + 1];
		vertex_pos = vertex_for_render[i];
		pair_num = vertex_obj_triangle_collider_num_record[i];
		primitive_record = vertex_obj_triangle_collider_pair_by_vertex[i];
		volume = VT_collider_volume[i].data();
		start_index = VT_collider_start_index[i];
		if (i < cloth->size()) {
			for (int j = vertex_index_start_per_thread[i][thread_No]; j < vertex_num; ++j) {
				computeVTVolume(vertex_pos[j].data(), pair_num[j], primitive_record + close_vt_collider_pair_num * j, 
					volume + start_index[j],
					vertex_for_render_collider);
			}
		}
		else {
			vertex_index_on_surface_ = vertex_index_on_surface[i];
			for (int j = vertex_index_start_per_thread[i][thread_No]; j < vertex_num; ++j) {
				computeVTVolume(vertex_pos[vertex_index_on_surface_[j]].data(), pair_num[j], primitive_record + close_vt_collider_pair_num * j, volume + start_index[j],
					vertex_for_render_collider);
			}
		}
	}
}




void Collision::computeVTVolume(int thread_No)
{
	unsigned int vertex_num;
	std::array<double, 3>* vertex_pos;
	unsigned int* pair_num;
	unsigned int* primitive_record;
	double* volume;
	unsigned int* start_index;
	unsigned int* vertex_index_on_surface_;
	for (int i = 0; i < total_obj_num; ++i) {
		vertex_num = vertex_index_start_per_thread[i][thread_No + 1];
		vertex_pos = vertex_for_render[i];
		pair_num = vertex_triangle_pair_num_record[i];
		primitive_record = vertex_triangle_pair_by_vertex[i];
		volume = VT_volume[i].data();
		start_index = VT_start_index[i];
		if (i < cloth->size()) {
			for (int j = vertex_index_start_per_thread[i][thread_No]; j < vertex_num; ++j) {
				computeVTVolume(vertex_pos[j].data(), pair_num[j], primitive_record + close_vt_pair_num * j, volume + start_index[j],
					vertex_for_render);
			}
		}
		else {
			vertex_index_on_surface_ = vertex_index_on_surface[i];
			for (int j = vertex_index_start_per_thread[i][thread_No]; j < vertex_num; ++j) {
				computeVTVolume(vertex_pos[vertex_index_on_surface_[j]].data(), pair_num[j], primitive_record + close_vt_pair_num * j, volume + start_index[j],
					vertex_for_render);
			}
		}
	}
}



void Collision::computeVTVolume(double* vertex_position_, unsigned int pair_num, unsigned int* triangle_record, double* volume,
 std::vector<std::array<double,3>*>& vertex_position)
{
	int* triangle_index;
	for (int i = 0; i < pair_num; i += 2) {
		triangle_index = triangle_indices[triangle_record[i]][triangle_record[i + 1]].data();
		*(volume++)=abs(getTetCubeVolume(vertex_position_, vertex_position[triangle_record[i]][triangle_index[0]].data(),
			vertex_position[triangle_record[i]][triangle_index[1]].data(), vertex_position[triangle_record[i]][triangle_index[2]].data()));
	}
}

void Collision::computeTVVolume(double* vertex_position_0, double* vertex_position_1, double* vertex_position_2, unsigned int pair_num,
	unsigned int* vertex_record, double* volume, std::vector<std::array<double, 3>*>& vertex_position)
{
	for (int i = 0; i < pair_num; i += 2) {
		*(volume++) =abs( getTetCubeVolume(vertex_position_0, vertex_position_1,
			vertex_position_2, vertex_position[vertex_record[i]][vertex_record[i+1]].data()));
	}
}

//void Collision::computeEEVolume(double* vertex_position_0, double* vertex_position_1,  unsigned int pair_num,
//	unsigned int* edge_record, double* volume, std::vector<std::array<double, 3>*>& vertex_position, std::vector<unsigned int*>& edge_vertices)
//{
//	unsigned int* edge_index;
//	for (int i = 0; i < pair_num; i += 2) {
//		edge_index = edge_vertices[edge_record[i]] + (edge_record[i + 1] << 1);
//		*(volume++) = abs(getTetCubeVolume(vertex_position_0, vertex_position_1,
//			vertex_position[edge_record[i]][*edge_index].data(),vertex_position[edge_record[i]][edge_index[1]].data()));
//	}
//}

void Collision::getCollisionPair()
{
	//for (int i = 0; i < collider->size(); ++i) {
	//	thread->assignTask(&(*collider)[i].mesh_struct, FACE_NORMAL);
	//}
	//for (int i = 0; i < cloth->size(); ++i) {
	//	thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	//}
	//for (int i = 0; i < tetrahedron->size(); ++i) {
	//	thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL);
	//}
	thread->assignTask(this, FIND_COLLISION_PAIR);
	setTargetPositionEven();	
}


unsigned int Collision::collisionConstraintNumber(unsigned int* point_triangle_collider_constraint, unsigned int* point_triangle_constraint, unsigned int* edge_edge_constraint)
{
	unsigned int count;
	
	memset(point_triangle_collider_constraint, 0, 4 * (thread_num + 1));
	memset(point_triangle_constraint, 0, 4 * (thread_num + 1));
	memset(edge_edge_constraint, 0, 4 * (thread_num + 1));

	if (has_collider) {
		for (int thread_No = 0; thread_No < thread_num; ++thread_No) {
			count = 0;
			if (vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]) {
				count += spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]][0]
					- vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1];
				for (int i = vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1] + 1;
					i < vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
					count += spatial_hashing.vertex_obj_triangle_collider_pair[i][0];
				}
				count += vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3];
			}
			else {
				count += vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3] -
					vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1];
			}
			count >>= 2;
			point_triangle_collider_constraint[thread_No + 1] = point_triangle_collider_constraint[thread_No] + count;
		}		
	}
	//for (int i = 0; i < thread_num; ++i) {
	//	////std::cout << point_triangle_collider_target_pos_index[i][0] / 4 << " ";
	//}
	//////std::cout << std::endl;
	//for (int i = 0; i < thread_num+1; ++i) {
	//	////std::cout << point_triangle_collider_constraint[i] << " ";
	//}
	//////std::cout << std::endl;

	point_triangle_constraint[0] = point_triangle_collider_constraint[thread_num];
	for (int thread_No = 0; thread_No < thread_num; ++thread_No) {
		count = 0;
		if (vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_triangle_pair_index_start_per_thread[thread_No << 1]) {
			count += spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]][0] - vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1];
			for (int i = vertex_triangle_pair_index_start_per_thread[thread_No << 1] + 1;
				i < vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
				count += spatial_hashing.vertex_triangle_pair[i][0];
			}
			count += vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3];
		}
		else {
			count += vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3] - vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1];
		}
		count >> 2;
		point_triangle_constraint[thread_No + 1] = point_triangle_constraint[thread_No] + count;
	}

	edge_edge_constraint[0] = point_triangle_constraint[thread_num];
	for (int thread_No = 0; thread_No < thread_num; ++thread_No) {
		count = 0;
		if (edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1] > edge_edge_pair_index_start_per_thread[thread_No << 1]) {
			count += spatial_hashing.edge_edge_pair[edge_edge_pair_index_start_per_thread[thread_No << 1]][0]
				- edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1];
			for (int i = edge_edge_pair_index_start_per_thread[thread_No << 1] + 1;
				i < edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
				count += spatial_hashing.edge_edge_pair[i][0];
			}
			count += edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3];
		}
		else {
			count += edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3] -
				edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1];
		}
		count >> 2;
		edge_edge_constraint[thread_No + 1] = edge_edge_constraint[thread_No] + count;
	}

	//////std::cout << point_triangle_collider_constraint[thread_num] << " " << point_triangle_constraint[thread_num] << " " << edge_edge_constraint[thread_num] << std::endl;

	return edge_edge_constraint[thread_num];
	//point_triangle_collider_constraint = 0;
	//point_triangle_constraint = 0;
	//edge_edge_constraint = 0;
	//for (int i = 0; i < thread_num; ++i) {
	//	point_triangle_constraint += (point_triangle_target_pos_index[i][0] >> 2);
	//	point_triangle_collider_constraint += (point_triangle_collider_target_pos_index[i][0] >> 2);
	//	edge_edge_constraint += (edge_edge_target_pos_index[i][0] >> 2);
	//}
	//return point_triangle_collider_constraint+ point_triangle_constraint+ edge_edge_constraint;
}








void Collision::solveCollisionConstraintDCD()
{
	for (int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL);
	}

	//time_t t = clock();
	//for (int i = 0; i < 10; ++i) {
		thread->assignTask(this, COLLISION_CONSTRAINT);
		sumTargetPosition();
		setTargetPositionEven();
	//}
	//time_t t1 = clock();
	//////std::cout << "DCD " << t1 - t << std::endl;
	//for (int i = 0; i < thread_num; ++i) {
	//	collisionConstraint(i);
	//}


	
}

void Collision::reSolveCollisionConstraintDCD()
{
	for (int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL);
	}
	thread->assignTask(this, RE_COLLISION_CONSTRAINT);

	//for (int i = 0; i < thread_num; ++i) {
	//	re_collisionConstraint(i);
	//}

	resumTargetPosition();
}


void Collision::XPBDsolveCollisionConstraint()
{
	for (int i = 0; i < thread_num; ++i) {
		point_triangle_target_pos_index[i][0] = 0;
		edge_edge_target_pos_index[i][0] = 0;
		if (has_collider) {
			point_triangle_collider_target_pos_index[i][0] = 0;
			//point_collider_triangle_target_pos_index[thread_No][0] = 0;
			//edge_edge_collider_target_pos_index[thread_No][0] = 0;
		}
	}


	if (floor->exist) {
		XPBDfloorCollisionResponse();
	}
	if (has_collider) {
		for (int i = 0; i < thread_num; ++i) {
			solveXPBDpointTriangleColliderResponse(i);
		}
	}
	for (int i = 0; i < thread_num; ++i) {
		solveXPBDpointTriangleResponse(i);
	}
	for (int i = 0; i < thread_num; ++i) {			
		solveXPBDedgeEdgeResponse(i);			
	}

	setTargetPositionEven();

	//unsigned int j = 0;
	//unsigned int k = 0;
	//for (int i = 0; i < thread_num; ++i) {
	//	j += point_triangle_target_pos_index[i][0];
	//	k += edge_edge_target_pos_index[i][0];
	//}
	//////std::cout << "vt " << j / 4 << " ee " << k / 4 << std::endl;
}

void Collision::re_XPBDsolveCollisionConstraint()
{
	if (floor->exist) {
		XPBDfloorCollisionResponse();
	}
	if (has_collider) {
		for (int i = 0; i < thread_num; ++i) {
			re_solveXPBDpointTriangleColliderResponse(i);
		}
	}
	for (int i = 0; i < thread_num; ++i) {
		re_solveXPBDpointTriangleResponse(i);
	}
	for (int i = 0; i < thread_num; ++i) {
		re_solveXPBDedgeEdgeResponse(i);
	}
}



void Collision::XPBDfloorCollisionResponse()
{
	double energy_=0.0;
	double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	int* indices;
	double target_pos_v[3];
	double target_pos_tri[3][3];

	unsigned int* pair;
	//double* lambda_ = lambda->data() + collision_lambda_index_start[0];

	unsigned int dimension = floor->dimension;
    bool normal_direction = floor->normal_direction;

	//////std::cout << dimension << " " << normal_direction << std::endl;
	double friction_coe_floor = friction_coe[2];
	unsigned int size = 0;
	double floor_value = floor->value;
	for (int i = 0; i < cloth->size(); ++i) {
		size = cloth->data()[i].mesh_struct.vertex_for_render.size();
		for (int j = 0; j < size; ++j) {
			dcd.XPBDFloor(vertex_for_render[i][j].data(), vertex_position[i][j].data(), dimension, normal_direction, tolerance,floor_value, friction_coe_floor);
		}	
	}
	unsigned int* index;
	unsigned int obj_index;
	for (int i = 0; i < tetrahedron->size(); ++i) {
		index = vertex_index_on_surface[i+cloth->size()];
		size = vertex_index_start_per_thread[i+cloth->size()][thread_num];
		obj_index = i + cloth->size();
		for (int j = 0; j < size; ++j) {
			dcd.XPBDFloor(vertex_for_render[obj_index][index[j]].data(), vertex_position[obj_index][index[j]].data(), dimension, normal_direction,
				tolerance, floor_value, friction_coe_floor);
		}
	}

}

void Collision::solveCollisionConstraintForIPC()
{
	for (int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
	thread->assignTask(this, COLLISION_CONSTRAINT_IPC);
	sumTargetPosition();
	setTargetPositionEven();
}

void Collision::re_solveCollisionConstraintForIPC()
{
	thread->assignTask(this,	RE_COLLISION_CONSTRAINT_IPC);
	resumTargetPosition();
}


void Collision::solveCollisionConstraint()
{
	for (int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL_RENDER);
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL_RENDER);
	}


	thread->assignTask(this, COLLISION_CONSTRAINT);
	sumTargetPosition();
	setTargetPositionEven();

	//testIfBuildCollisionConstraint()
}


void Collision::testIfBuildCollisionConstraint()
{
	for (int i = 0; i < obj_target_pos.b_sum[0].size(); ++i) {
		if (obj_target_pos.need_update[0][i]) {
			//////std::cout << "build collision constraint "<<i << std::endl;
		}
	}
}

//FIND_COLLISION_PAIR
void Collision::getCollisionPair(int thread_No)
{
	checkTargetPosSize(thread_No);
	point_triangle_target_pos_index[thread_No][0] = 0;
	edge_edge_target_pos_index[thread_No][0] = 0;
	if (has_collider) {
		//point_collider_triangle_target_pos_index[thread_No][0] = 0;
		point_triangle_collider_target_pos_index[thread_No][0] = 0;
		//edge_edge_collider_target_pos_index[thread_No][0] = 0;
	}
	pointTriangleCollisionPair(thread_No);
	edgeEdgeCollisionPair(thread_No);
	if (has_collider) {
		pointTriangleColliderPair(thread_No);
	}
}



//FIND_CLOSE_PAIR
void Collision::findClosePair(int thread_No)
{
	findVT_ClosePair(thread_No);
	findEE_ClosePair(thread_No);
	if (has_collider) {
		findVT_ColliderClosePair(thread_No);
		findEE_ColliderClosePair(thread_No);
		findTV_ColliderClosePair(thread_No);
	}
}


void Collision::findClosePair()
{
	initialPairByElement();


	thread->assignTask(this, FIND_CLOSE_PAIR);

	findEE_ClosePair_ReverseOrder();
	findTV_ClosePair();

	//for (int i = 0; i < total_obj_num; ++i) {
	//	memset(vertex_triangle_pair_num_record[i], 0, 4 * vertex_index_start_per_thread[i][thread_num]);
	//	memset(triangle_vertex_pair_num_record[i], 0, 4 * mesh_struct[i]->triangle_indices.size());
	//	memset(edge_edge_pair_number_record[i], 0, 4 * edge_index_start_per_thread[i][thread_num]);
	//	if (has_collider) {
	//		memset(vertex_obj_triangle_collider_num_record[i], 0, 4 * vertex_index_start_per_thread[i][thread_num]);
	//		memset(triangle_vertex_collider_pair_num_record[i], 0, 4 * triangle_index_start_per_thread[i][thread_num]);
	//		memset(edge_edge_collider_pair_num_record[i], 0, 4 * edge_index_start_per_thread[i][thread_num]);
	//	}
	//}

	thread->assignTask(this, PREFIX_SUM_ALL_PAIRS);
	recordPairCompress();
	resizeHessianRecordIndex();
	setSelfCollisionPairPerThread();
	if (has_collider) {
		thread->assignTask(this, SET_ELEMENT_COLLIDE_WITH_COLLIDER);
	}
	addFlagToColorGroup();
	//testColliderPair();
}


void Collision::initialPairCompress()
{
	temp_save_ee_pair_compress_per_thread.resize(thread_num - 1);
	for (int i = 0; i < thread_num - 1; ++i) {
		temp_save_ee_pair_compress_per_thread[i].reserve(100);
	}

	tet_involved_in_collision.resize(total_obj_num);
	tet_involved_in_collision_start_per_thread.resize(total_obj_num);

	int i;
	for (int tet= 0; tet < tetrahedron->size(); ++tet) {
		i = tet + cloth->size();
		tet_involved_in_collision[i].resize(tet_color_groups[i]->size());
		for (auto j = tet_involved_in_collision[i].begin(); j < tet_involved_in_collision[i].end(); ++j) {
			j->reserve(tetrahedron->data()[tet].mesh_struct.indices.size() / 8);
		}
		tet_involved_in_collision_start_per_thread[i].resize(tet_color_groups[i]->size());

		for (auto j = tet_involved_in_collision_start_per_thread[i].begin(); j < tet_involved_in_collision_start_per_thread[i].end(); ++j) {
			j->resize(thread_num + 1);
		}

	}

	triangle_index_has_TV_pair.resize(total_obj_num);
	edge_index_has_EE_pair.resize(total_obj_num);
	vertex_index_has_VT_pair.resize(total_obj_num);
	temp_save_element_in_collision_pair.resize(thread_num - 1);
	triangle_index_has_TV_pair_start_per_thread.resize((thread_num + 1) * total_obj_num);
	edge_index_has_EE_pair_start_per_thread.resize((thread_num + 1) * total_obj_num);
	vertex_index_has_VT_pair_start_per_thread.resize((thread_num + 1) * total_obj_num);

	triangle_index_collide_with_collider_start_per_thread.resize(total_obj_num * (thread_num + 1));
	edge_index_collide_with_collider_start_per_thread.resize(total_obj_num * (thread_num + 1));
	vertex_index_collide_with_collider_start_per_thread.resize(total_obj_num * (thread_num + 1));

	for (int i = 0; i < total_obj_num; ++i) {
		triangle_index_has_TV_pair[i].reserve(mesh_struct[i]->triangle_indices.size() / 3);
		edge_index_has_EE_pair[i].reserve(mesh_struct[i]->edge_vertices.size() / 6);
		vertex_index_has_VT_pair[i].reserve(mesh_struct[i]->triangle_indices.size() / 4);
	}

	for (int j = 0; j < thread_num - 1; ++j) {
		temp_save_element_in_collision_pair[j].resize(total_obj_num);
		for (int i = 0; i < total_obj_num; ++i) {
			temp_save_element_in_collision_pair[j][i].reserve(mesh_struct[i]->edge_vertices.size() / 6);
		}
	}


}

void Collision::recordPairCompress()
{
	vt_pair_compressed_record.resize(
		5 * vertex_triangle_pair_num_record_prefix_sum[total_obj_num-1][vertex_index_start_per_thread[total_obj_num-1][thread_num]]);

	if (ee_pair_compressed_record.capacity() < 3 * edge_edge_pair_num_record_prefix_sum[total_obj_num-1]
		[vertex_index_start_per_thread[total_obj_num-1][thread_num]]) {
		ee_pair_compressed_record.reserve(
				3 * edge_edge_pair_num_record_prefix_sum[total_obj_num-1][vertex_index_start_per_thread[total_obj_num-1][thread_num]]);
	}

	ee_pair_compressed_record.clear();
	for (int i = 0; i < thread_num - 1; ++i) {
		if (temp_save_ee_pair_compress_per_thread[i].capacity() < ee_pair_compressed_record.capacity() / thread_num) {
			temp_save_ee_pair_compress_per_thread[i].reserve(ee_pair_compressed_record.capacity() / thread_num);
		}
		temp_save_ee_pair_compress_per_thread[i].clear();
	}	
	thread->assignTask(this, RECORD_VT_PAIR_COMPRESS);

	for (int i = 0; i < temp_save_element_in_collision_pair.size(); ++i) {
		for (int j = 0; j < total_obj_num; ++j) {
			vertex_index_has_VT_pair[j].insert(vertex_index_has_VT_pair[j].end(), temp_save_element_in_collision_pair[i][j].begin(), temp_save_element_in_collision_pair[i][j].end());
		}		
	}

	//std::cout << "compress pair " << std::endl;
	//for (int j = 0; j < vt_pair_compressed_record.size(); j+=5) {
	//	std::cout << vt_pair_compressed_record[j] << " " << vt_pair_compressed_record[j + 1] << " " << vt_pair_compressed_record[j + 2] << " " << vt_pair_compressed_record[j + 3] << std::endl;
	//}
	//std::cout << "vertex_index_has_VT_pair pair " << std::endl;
	//for (int j = 0; j < total_obj_num; ++j) {
	//	for (int i = 0; i < vertex_index_has_VT_pair[j].size(); i++) {
	//		std::cout << j << " " << vertex_index_has_VT_pair[j][i] << std::endl;
	//	}
	//}

	thread->assignTask(this, RECORD_EE_PAIR_COMPRESS);
	for (int i = 0; i < temp_save_ee_pair_compress_per_thread.size(); ++i) {
		ee_pair_compressed_record.insert(ee_pair_compressed_record.end(), temp_save_ee_pair_compress_per_thread[i].begin(), temp_save_ee_pair_compress_per_thread[i].end());
	}
	for (int i = 0; i < temp_save_element_in_collision_pair.size(); ++i) {
		for (int j = 0; j < total_obj_num; ++j) {
			edge_index_has_EE_pair[j].insert(edge_index_has_EE_pair[j].end(), temp_save_element_in_collision_pair[i][j].begin(), temp_save_element_in_collision_pair[i][j].end());
		}
	}
	thread->assignTask(this, RECORD_TRIANGLE_HAS_COLLISION_PAIR);
	for (int i = 0; i < temp_save_element_in_collision_pair.size(); ++i) {
		for (int j = 0; j < total_obj_num; ++j) {
			triangle_index_has_TV_pair[j].insert(triangle_index_has_TV_pair[j].end(), temp_save_element_in_collision_pair[i][j].begin(), temp_save_element_in_collision_pair[i][j].end());
		}
	}


	for (int i = 0; i < total_obj_num; ++i) {
		arrangeIndex(thread_num, triangle_index_has_TV_pair[i].size(), triangle_index_has_TV_pair_start_per_thread.data() + i * (thread_num + 1));
		arrangeIndex(thread_num, edge_index_has_EE_pair[i].size(), edge_index_has_EE_pair_start_per_thread.data() + i * (thread_num + 1));
		arrangeIndex(thread_num, vertex_index_has_VT_pair[i].size(), vertex_index_has_VT_pair_start_per_thread.data() + i * (thread_num + 1));
	}

	
	

}


//RECORD_VT_PAIR_COMPRESS
void Collision::recordVTPairCompress(int thread_No)
{
	if (thread_No == 0) {
		recordVTCollisionPairCompress(thread_No, vertex_index_start_per_thread.data(), vertex_triangle_pair_num_record, vertex_triangle_pair_by_vertex,
			vertex_triangle_pair_num_record_prefix_sum, vt_pair_compressed_record.data(), vertex_index_on_surface.data(), vertex_index_has_VT_pair.data());
	}
	else {
		recordVTCollisionPairCompress(thread_No, vertex_index_start_per_thread.data(), vertex_triangle_pair_num_record, vertex_triangle_pair_by_vertex,
			vertex_triangle_pair_num_record_prefix_sum, vt_pair_compressed_record.data(), vertex_index_on_surface.data(), temp_save_element_in_collision_pair[thread_No-1].data());
	}


}


//RECORD_EE_PAIR_COMPRESS
void Collision::recordEEPairCompress(int thread_No)
{
	if (thread_No == 0) {
		recordEECollisionPairCompress(thread_No, edge_index_start_per_thread.data(), edge_edge_pair_number_record, edge_edge_pair_by_edge,
			&ee_pair_compressed_record, edge_index_has_EE_pair.data());
	}
	else {
		recordEECollisionPairCompress(thread_No, edge_index_start_per_thread.data(), edge_edge_pair_number_record, edge_edge_pair_by_edge,
			&temp_save_ee_pair_compress_per_thread[thread_No - 1], temp_save_element_in_collision_pair[thread_No - 1].data());
	}

}

//RECORD_TRIANGLE_HAS_COLLISION_PAIR
void Collision::recordTriangleHasTVPair(int thread_No)
{
	if (thread_No == 0) {
		recordTriangleHasTVPair(thread_No, triangle_index_start_per_thread.data(), triangle_vertex_pair_num_record, 
			triangle_index_has_TV_pair.data());
	}
	else {
		recordTriangleHasTVPair(thread_No, triangle_index_start_per_thread.data(), triangle_vertex_pair_num_record,
			temp_save_element_in_collision_pair[thread_No - 1].data());
	}
}



void Collision::recordTriangleHasTVPair(int thread_No, unsigned int** start_per_thread, unsigned int** pair_num_record, 
	std::vector<unsigned int>* has_pair)
{
	int start, end;
	unsigned int* num_record;
	unsigned int* vertex_surface_to_global_;
	unsigned int* pair_;
	std::vector<unsigned int>* has_pair_;
	for (int i = 0; i < total_obj_num; ++i) {
		num_record = pair_num_record[i];
		start = start_per_thread[i][thread_No];
		end = start_per_thread[i][thread_No + 1];
		has_pair_ = &has_pair[i];
		has_pair_->clear();
		for (int j = start; j < end; ++j) {
			if (num_record[j] != 0) {
				has_pair_->emplace_back(j);
			}
		}
	}
}


void Collision::recordEECollisionPairCompress(int thread_No, unsigned int** start_per_thread, unsigned int** pair_num_record, unsigned int** pair,
	std::vector<unsigned int>* pair_compres_record, std::vector<unsigned int>* has_pair)
{
	int start, end;
	unsigned int* num_record;
	unsigned int* vertex_surface_to_global_;
	unsigned int* pair_;
	std::vector<unsigned int>* has_pair_;
	for (int i = 0; i < total_obj_num; ++i) {
		num_record = pair_num_record[i];
		start = start_per_thread[i][thread_No];
		end = start_per_thread[i][thread_No + 1];
		has_pair_ = &has_pair[i];
		has_pair_->clear();
		for (int j = start; j < end; ++j) {
			if (num_record[j] != 0) {
				pair_ = pair[i] + j * close_ee_pair_num;
				for (int k = 0; k < num_record[j]; k += 3) {
					if (pair_[k] < i || (pair_[k] == i && pair_[k + 1] < j)) {
						break;
					}
					pair_compres_record->emplace_back(i);
					pair_compres_record->emplace_back(j);
					pair_compres_record->emplace_back(pair_[k]);
					pair_compres_record->emplace_back(pair_[k + 1]);
					pair_compres_record->emplace_back(k/3);
				}
				has_pair_->emplace_back(j);
			}
		}		
	}
}



void Collision::recordVTCollisionPairCompress(int thread_No, unsigned int** start_per_thread, unsigned int** pair_num_record, unsigned int** pair,
	unsigned int** prefix_sum, unsigned int* pair_compress_record, unsigned int** vertex_surface_to_global, std::vector<unsigned int>* has_pair)
{
	int start, end;
	unsigned int* num_record;
	unsigned int* prefix_sum_;
	unsigned int* record;
	unsigned int* vertex_surface_to_global_;
	unsigned int* pair_;
	std::vector<unsigned int>* has_pair_;

	for (int i = 0; i < total_obj_num; ++i) {
		num_record = pair_num_record[i];
		prefix_sum_ = prefix_sum[i];
		start = start_per_thread[i][thread_No];
		end = start_per_thread[i][thread_No + 1];
		has_pair_ = &has_pair[i];
		has_pair_->clear();
		if (i < cloth->size()) {
			for (int j = start; j < end; ++j) {
				if (num_record[j] != 0) {
					pair_ = pair[i] + j * close_vt_pair_num;

					for (int k = 0; k < num_record[j]; k += 2) {
						record = pair_compress_record + 5 * (prefix_sum_[j] + (k >> 1));
						*record = i;
						*(record + 1) = j;
						*(record + 2) = pair_[k];
						*(record + 3) = pair_[k + 1];
						*(record + 4) = k >> 1;
					}
					has_pair_->emplace_back(j);

				}
			}
		}
		else {
			vertex_surface_to_global_ = vertex_surface_to_global[i];
			for (int j = start; j < end; ++j) {
				if (num_record[j] != 0) {
					pair_ = pair[i] + j * close_vt_pair_num;
					for (int k = 0; k < num_record[j]; k += 2) {
						record = pair_compress_record + 5 * (prefix_sum_[j] + (k >> 1));
						*record = i;
						*(record + 1) = vertex_surface_to_global_[j];
						*(record + 2) = pair_[k];
						*(record + 3) = pair_[k + 1];
						*(record + 4) = k >> 1;

						//std::cout << "collision pair record " << i << " " << vertex_surface_to_global_[j] << " " << pair_[k] << " " << pair_[k + 1] << std::endl;

					}
					has_pair_->emplace_back(vertex_surface_to_global_[j]);
				}				
			}
		}	
	}
}



void Collision::computeFloorHessian(int color_No)
{
	memset(floor_grad_record.data(), 0, floor_grad_record.size() << 3);
	memset(floor_hessian_record.data(), 0, floor_hessian_record.size() << 3);
}



void Collision::computeHessian(int color_No)
{

	updateVertexBelongColorGroup(color_No);
	//int color_group_index = 0;
	//if (!tetrahedron->empty()) {
	//	color_group_index = *inner_iteration_number % tet_color_groups[cloth->size()]->size();

	//}
	//else {
	//	std::cout << "have not achieve when there is no tet, Collision::computeHessianPerThread(" << std::endl;
	//}

	if (floor->exist) {
		memset(floor_grad_record.data(), 0, floor_grad_record.size() << 3);
		memset(floor_hessian_record.data(), 0, floor_hessian_record.size() << 3);
	}



	thread->assignTask(this, UPDATE_COLLISION_HESSIAN_COLOR, color_No);
	//std::cout << " if_used.size() " << ee_pair_compressed_record.size()/5<<" "<< if_used.size()/5 << std::endl;

	//std::cout << edge_edge_pair_num_record_prefix_sum[total_obj_num - 1][mesh_struct[total_obj_num - 1]->edge_length.size()] << " "
	//	<< ee_hessian_record_index.size() / 5 << std::endl;

	//for (int i = 0; i < thread_num; ++i) {
	//	computeHessianPerThread(i, color_No);
	//}

}



void Collision::updateVertexBelongColorGroup(int color_No)
{
	std::vector<unsigned int>*index_group;
	bool* belong;
	int color_group_index;

	int i;
	for (int tet_No = 0; tet_No < tetrahedron->size(); ++tet_No) {
		i = tet_No + cloth->size();
		color_group_index = *inner_iteration_number % tet_color_groups[i]->size();
		belong = vertex_belong_to_color_group[i];
		memset(belong, 0, total_vertex_num[i]);
		//if (color_No > tet_color_groups[i]->data()[color_group_index].size()-1) {
		//	continue;
		//}
		index_group = &surface_vertex_index_of_a_tet_color_group[i][color_group_index][color_No];
		for (auto j = index_group->begin(); j < index_group->end(); ++j) {
			belong[*j] = true;
		}
	}
}


void Collision::resizeHessianRecordIndex()
{
	int size = vertex_triangle_pair_num_record_prefix_sum[total_obj_num - 1][vertex_index_start_per_thread[total_obj_num - 1][thread_num]];
	vt_hessian_record_index.resize(5 * size);
	vt_hessian_record.resize(144 * size);
	vt_grad_record.resize(12 * size);



	size = ee_pair_compressed_record.size()/5*2;
	ee_hessian_record_index.resize(5 * size);
	ee_hessian_record.resize(144 * size);
	ee_grad_record.resize(12 * size);

	if (has_collider) {
		size= triangle_vertex_collider_num_record_prefix_sum[total_obj_num - 1][triangle_index_start_per_thread[total_obj_num - 1][thread_num]];
		tv_colldier_hessian_record_index.resize(4 * size);
		tv_colldier_hessian_record.resize(81 * size);
		tv_colldier_grad_record.resize(9 * size);

		size = edge_edge_collider_num_record_prefix_sum[total_obj_num - 1][edge_index_start_per_thread[total_obj_num - 1][thread_num]];
		ee_collider_hessian_record_index.resize(3 * size);
		ee_collider_hessian_record.resize(36 * size);
		ee_collider_grad_record.resize(6 * size);
	}
}



//void Collision::computeFloorHessian(int thread_No, int color_No)
//{
//	if (floor->exist) {
//
//		int color_group_index = 0;
//		color_group_index = *inner_iteration_number % tet_color_groups[cloth->size()]->size();
//
//
//		double stiffness;
//		if (!tetrahedron->empty()) {
//			stiffness = tetrahedron->data()[0].collision_stiffness[0];
//		}
//		else {
//			stiffness = cloth->data()[0].collision_stiffness[0];
//		}
//
//		double floor_value = floor->value;
//		int floor_dimension = floor->dimension;
//
//		std::array<double, 3>* vertex_pos;
//		int end_per_thread;
//		int vertex_prefix_sum_this_obj;
//		int i;
//
//		//surface_vertex_index_of_a_tet_color_group
//		unsigned int* surface_vertex;
//
//		int* global_to_surface_index;
//
//		for (int tet_obj_no = 0; tet_obj_no < tetrahedron->size(); ++tet_obj_no) {
//			i = tet_obj_no + cloth->size();
//			vertex_pos = vertex_position[i];
//			vertex_prefix_sum_this_obj = vertex_num_on_surface_prefix_sum[i];
//			end_per_thread = surface_vertex_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color_No][thread_No + 1];
//
//			surface_vertex = surface_vertex_index_of_a_tet_color_group[i][color_group_index][color_No].data();
//			global_to_surface_index = general_index_to_surface_index[i];
//			for (int j = surface_vertex_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color_No][thread_No]; j < end_per_thread; ++j) {
//				computeFloorHessian(d_hat_2, stiffness, floor_value, floor_hessian_record[vertex_prefix_sum_this_obj + global_to_surface_index[surface_vertex[j]]], floor_grad_record[vertex_prefix_sum_this_obj + global_to_surface_index[surface_vertex[j]]],
//					vertex_pos[surface_vertex[j]][floor_dimension]);
//			}
//		}
//	}
//}


void Collision::computeHessianPreviousThread(int thread_No, int color_No)
{
	int i = -1; //obj_no
	int color_group_index = 0;
	double stiffness = 0.0;


	color_group_index = *inner_iteration_number % tet_color_groups[cloth->size()]->size();
	if (!tetrahedron->empty()) {
		stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		stiffness = cloth->data()[0].collision_stiffness[0];
	}



	int j;//vertex_index
	int* vertex_index_on_surface;

	std::array<double, 3>* vertex_pos;



	unsigned int* edge_vertex;
	unsigned int* edge_vertex_record;
	unsigned int* edge_vertex_compare;
	double* edge_length;
	int surface_index;

	unsigned int* prefix_sum;
	unsigned int* prefix_sum_collider;

	bool has_collider = this->has_collider;
	bool has_floor = floor->exist;
	double floor_value = floor->value;
	int floor_dimension = floor->dimension;


	int* triangle_vertex;
	MatrixXd Hessian; VectorXd grad;

	int start, end;


	//vt
	start = vt_per_thread_start_index[thread_No];
	end = vt_per_thread_start_index[thread_No + 1];

	int obj_2;

	if (start < end) {
		i = vt_pair_compressed_record[start];
		vertex_pos = vertex_position[i];
		vertex_index_on_surface = this->general_index_to_surface_index[i];
		prefix_sum = vertex_triangle_pair_num_record_prefix_sum[i];

		for (int index = start; index < end; index += 5) {
			if (i != vt_pair_compressed_record[index]) {
				i = vt_pair_compressed_record[index];
				vertex_pos = vertex_position[i];
				vertex_index_on_surface = this->general_index_to_surface_index[i];
				prefix_sum = vertex_triangle_pair_num_record_prefix_sum[i];
			}
			j = vt_pair_compressed_record[index + 1];
			surface_index = vertex_index_on_surface[j];

			obj_2 = vt_pair_compressed_record[index + 2];
			triangle_vertex = triangle_indices[obj_2][vt_pair_compressed_record[index + 3]].data();

			if (vertex_belong_to_color_group[i][j] || vertex_belong_to_color_group[obj_2][triangle_vertex[0]]
				|| vertex_belong_to_color_group[obj_2][triangle_vertex[1]] || vertex_belong_to_color_group[obj_2][triangle_vertex[2]])
			{
				if (second_order_constraint.computeBarrierVTGradientHessian(Hessian, grad, vertex_pos[j].data(), vertex_position[obj_2][triangle_vertex[0]].data(),
					vertex_position[obj_2][triangle_vertex[1]].data(), vertex_position[obj_2][triangle_vertex[2]].data(), d_hat_2,
					vt_hessian_record_index.data() + 5 * (prefix_sum[surface_index] + vt_pair_compressed_record[index + 4]), stiffness)) {
					memcpy(vt_hessian_record.data() + (prefix_sum[surface_index] + vt_pair_compressed_record[index + 4]) * 144, Hessian.data(), Hessian.size() << 3);
					memcpy(vt_grad_record.data() + (prefix_sum[surface_index] + vt_pair_compressed_record[index + 4]) * 12, grad.data(), grad.size() << 3);
				}
			}
		}
	}

	//ee
	start = ee_per_thread_start_index[thread_No];
	end = ee_per_thread_start_index[thread_No + 1];

	if (start < end) {
		i = ee_pair_compressed_record[start];
		vertex_pos = vertex_position[i];
		prefix_sum = edge_edge_pair_num_record_prefix_sum[i];
		edge_vertex = this->edge_vertices[i];
		edge_length = rest_edge_length[i];

		for (int index = start; index < end; index += 5) {
			if (i != ee_pair_compressed_record[index]) {
				i = ee_pair_compressed_record[index];
				vertex_pos = vertex_position[i];
				prefix_sum = edge_edge_pair_num_record_prefix_sum[i];
				edge_vertex = this->edge_vertices[i];
				edge_length = rest_edge_length[i];
			}

			j = ee_pair_compressed_record[index + 1];
			edge_vertex_record = edge_vertex + (j << 1);

			obj_2 = ee_pair_compressed_record[index + 2];

			edge_vertex_compare = edge_vertices[obj_2] + (ee_pair_compressed_record[index + 3] << 1);

			if (vertex_belong_to_color_group[i][edge_vertex_record[0]] || vertex_belong_to_color_group[i][edge_vertex_record[1]]
				|| vertex_belong_to_color_group[obj_2][edge_vertex_compare[0]] || vertex_belong_to_color_group[obj_2][edge_vertex_compare[1]]) {
				if (second_order_constraint.computeBarrierEEGradientHessian(vertex_pos[edge_vertex_record[0]].data(), vertex_pos[edge_vertex_record[1]].data(),
					vertex_position[obj_2][edge_vertex_compare[0]].data(),
					vertex_position[obj_2][edge_vertex_compare[1]].data(), Hessian, grad,
					ee_hessian_record_index.data() + 5 * (prefix_sum[j] + ee_pair_compressed_record[index + 4]), stiffness, d_hat_2,
					edge_length[j], rest_edge_length[obj_2][ee_pair_compressed_record[index + 3]])) {
					memcpy(ee_hessian_record.data() + 144 * (prefix_sum[j] + ee_pair_compressed_record[index + 4]), Hessian.data(), Hessian.size() << 3);
					memcpy(ee_grad_record.data() + 12 * (prefix_sum[j] + ee_pair_compressed_record[index + 4]), grad.data(), grad.size() << 3);

					//if (ee_pair_compressed_record[index] == 0 && ee_pair_compressed_record[index + 1] == 28) {
					//	if (ee_pair_compressed_record[index + 2] == 1 && (ee_pair_compressed_record[index + 3] == 102 || ee_pair_compressed_record[index + 3] == 106)) {
					//		int* add = ee_hessian_record_index.data() + 5 * (prefix_sum[j] + ee_pair_compressed_record[index + 4]);
					//		std::cout << "index " << ee_pair_compressed_record[index + 3] << std::endl;
					//		std::cout << add[0] << " " << add[1] << " " << add[2] << " " << add[3] << " " << add[4] << std::endl;
					//	}
					//}

				}
			}
		}
	}
	int vertex_prefix_sum_this_obj;

	if (has_collider)
	{
		unsigned int* element_collide_with_collider;
		unsigned int* pair_by_element;
		unsigned int* pair_by_element_num;

		bool* vertex_belong;

		for (int tet_obj_no = 0; tet_obj_no < tetrahedron->size(); ++tet_obj_no) {

			i = tet_obj_no + cloth->size();
			vertex_pos = vertex_position[i];
			vertex_belong = vertex_belong_to_color_group[i];
			//vt collider
			vertex_prefix_sum_this_obj = vertex_num_on_surface_prefix_sum[i];
			start = vertex_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No];
			end = vertex_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No + 1];
			vertex_index_on_surface = this->general_index_to_surface_index[i];
			element_collide_with_collider = vertex_index_collide_with_collider[i].data();
			pair_by_element = vertex_obj_triangle_collider_pair_by_vertex[i];
			pair_by_element_num = vertex_obj_triangle_collider_num_record[i];
			for (int k = start; k < end; k++) {
				j = element_collide_with_collider[k];
				if (vertex_belong[j]) {
					surface_index = vertex_index_on_surface[j];
					computeVTColliderHessian(pair_by_element
						+ close_vt_collider_pair_num * surface_index,
						pair_by_element_num[surface_index], d_hat_2,
						vertex_pos[j].data(), vt_colldier_hessian_record.data() + (vertex_prefix_sum_this_obj + surface_index) * 9,
						stiffness, vt_colldier_grad_record.data() + 3 * (surface_index + vertex_prefix_sum_this_obj));
				}
			}

			//ee collider
			start = edge_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No];
			end = edge_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No + 1];
			pair_by_element = edge_edge_collider_pair_by_edge[i];
			pair_by_element_num = edge_edge_collider_pair_num_record[i];
			element_collide_with_collider = edge_index_collide_with_collider[i].data();
			edge_vertex = this->edge_vertices[i];
			prefix_sum_collider = edge_edge_collider_num_record_prefix_sum[i];
			edge_length = rest_edge_length[i];
			for (int k = start; k < end; k++) {
				j = element_collide_with_collider[k];
				if (vertex_belong[edge_vertex[j << 1]] || vertex_belong[edge_vertex[(j << 1) + 1]]) {
					computeEEColliderHessian(pair_by_element + close_ee_collider_pair_num * j,
						pair_by_element_num[j],
						d_hat_2, vertex_pos[edge_vertex[j << 1]].data(),
						vertex_pos[edge_vertex[(j << 1) + 1]].data(),
						ee_collider_hessian_record.data() + 36 * prefix_sum_collider[j],
						ee_collider_hessian_record_index.data() + 3 * prefix_sum_collider[j], stiffness, ee_collider_grad_record.data() + 6 * prefix_sum_collider[j],
						edge_length[j]);
				}
			}

			//tv collider
			start = triangle_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No];
			end = triangle_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No + 1];
			pair_by_element = triangle_vertex_collider_pair_by_triangle[i];
			pair_by_element_num = triangle_vertex_collider_pair_num_record[i];
			element_collide_with_collider = triangle_index_collide_with_collider[i].data();
			prefix_sum_collider = triangle_vertex_collider_num_record_prefix_sum[i];

			for (int k = start; k < end; k++) {
				j = element_collide_with_collider[k];
				triangle_vertex = triangle_indices[i][j].data();
				if (vertex_belong[triangle_vertex[0]] || vertex_belong[triangle_vertex[1]] || vertex_belong[triangle_vertex[2]]) {
					computeTVColliderHessian(pair_by_element + close_tv_collider_pair_num * j,
						pair_by_element_num[j], d_hat_2, vertex_pos[triangle_vertex[0]].data(),
						vertex_pos[triangle_vertex[1]].data(), vertex_pos[triangle_vertex[2]].data(),
						tv_colldier_hessian_record.data() + 81 * prefix_sum_collider[j], stiffness, tv_colldier_grad_record.data() + 9 * prefix_sum_collider[j],
						tv_colldier_hessian_record_index.data() + 4 * prefix_sum_collider[j]);
				}
			}
		}
	}

	if (has_floor) {
		int end_per_thread;
		unsigned int* surface_index_to_global;
		bool* vertex_belong;
		for (int tet_obj_no = 0; tet_obj_no < tetrahedron->size(); ++tet_obj_no) {
			i = tet_obj_no + cloth->size();
			vertex_belong = vertex_belong_to_color_group[i];
			vertex_pos = vertex_position[i];
			vertex_prefix_sum_this_obj = vertex_num_on_surface_prefix_sum[i];
			end_per_thread = vertex_index_start_per_thread[i][thread_No + 1];
			surface_index_to_global = this->vertex_index_on_surface[i];
			for (int j = vertex_index_start_per_thread[i][thread_No]; j < end_per_thread; ++j) {
				if (vertex_belong[surface_index_to_global[j]]) {
					computeFloorHessian(d_hat_2, stiffness, floor_value, floor_hessian_record[vertex_prefix_sum_this_obj + j], floor_grad_record[vertex_prefix_sum_this_obj + j],
						vertex_pos[surface_index_to_global[j]][floor_dimension]);
				}
			}
		}
	}
}


void Collision::computeLastColorHessianPerThread(int thread_No, int color_No)
{
	int i = -1; //obj_no
	int color_group_index = 0;
	double stiffness = 0.0;


	color_group_index = *inner_iteration_number % tet_color_groups[cloth->size()]->size();
	if (!tetrahedron->empty()) {
		stiffness = tetrahedron->data()[0].collision_stiffness[0];
	}
	else {
		stiffness = cloth->data()[0].collision_stiffness[0];
	}

	int j;//vertex_index
	int* vertex_index_on_surface;

	std::array<double, 3>* vertex_pos;



	unsigned int* edge_vertex;
	unsigned int* edge_vertex_record;
	unsigned int* edge_vertex_compare;
	double* edge_length;
	int surface_index;

	unsigned int* prefix_sum;
	unsigned int* prefix_sum_collider;

	bool has_collider = this->has_collider;
	bool has_floor = floor->exist;
	double floor_value = floor->value;
	int floor_dimension = floor->dimension;


	int* triangle_vertex;
	MatrixXd Hessian; VectorXd grad;

	int start, end;


	//vt
	start = vt_per_thread_start_index[thread_No];
	end = vt_per_thread_start_index[thread_No + 1];

	if (start < end) {
		i = vt_pair_compressed_record[start];
		vertex_pos = vertex_position[i];
		vertex_index_on_surface = this->general_index_to_surface_index[i];
		prefix_sum = vertex_triangle_pair_num_record_prefix_sum[i];

		for (int index = start; index < end; index += 5) {
			if (i != vt_pair_compressed_record[index]) {
				i = vt_pair_compressed_record[index];
				vertex_pos = vertex_position[i];
				vertex_index_on_surface = this->general_index_to_surface_index[i];
				prefix_sum = vertex_triangle_pair_num_record_prefix_sum[i];
			}

			j = vt_pair_compressed_record[index + 1];
			surface_index = vertex_index_on_surface[j];
			triangle_vertex = triangle_indices[vt_pair_compressed_record[index + 2]][vt_pair_compressed_record[index + 3]].data();
			if (second_order_constraint.computeBarrierVTGradientHessian(Hessian, grad, vertex_pos[j].data(), vertex_position[vt_pair_compressed_record[index + 2]][triangle_vertex[0]].data(),
				vertex_position[vt_pair_compressed_record[index + 2]][triangle_vertex[1]].data(), vertex_position[vt_pair_compressed_record[index + 2]][triangle_vertex[2]].data(), d_hat_2,
				vt_hessian_record_index.data() + 5 * (prefix_sum[surface_index] + vt_pair_compressed_record[index + 4]), stiffness)) {
				memcpy(vt_hessian_record.data() + (prefix_sum[surface_index] + vt_pair_compressed_record[index + 4]) * 144, Hessian.data(), Hessian.size() << 3);
				memcpy(vt_grad_record.data() + (prefix_sum[surface_index] + vt_pair_compressed_record[index + 4]) * 12, grad.data(), grad.size() << 3);
			}
		}


	}

	//ee
	start = ee_per_thread_start_index[thread_No];
	end = ee_per_thread_start_index[thread_No + 1];

	if (start < end) {
		i = ee_pair_compressed_record[start];
		vertex_pos = vertex_position[i];
		prefix_sum = edge_edge_pair_num_record_prefix_sum[i];
		edge_vertex = this->edge_vertices[i];
		edge_length = rest_edge_length[i];

		for (int index = start; index < end; index += 5) {
			if (i != ee_pair_compressed_record[index]) {
				i = ee_pair_compressed_record[index];
				vertex_pos = vertex_position[i];
				prefix_sum = edge_edge_pair_num_record_prefix_sum[i];
				edge_vertex = this->edge_vertices[i];
				edge_length = rest_edge_length[i];
			}

			j = ee_pair_compressed_record[index + 1];
			edge_vertex_record = edge_vertex + (j << 1);
			edge_vertex_compare = edge_vertices[ee_pair_compressed_record[index + 2]] + (ee_pair_compressed_record[index + 3] << 1);

			if (second_order_constraint.computeBarrierEEGradientHessian(vertex_pos[edge_vertex_record[0]].data(), vertex_pos[edge_vertex_record[1]].data(),
				vertex_position[ee_pair_compressed_record[index + 2]][edge_vertex_compare[0]].data(),
				vertex_position[ee_pair_compressed_record[index + 2]][edge_vertex_compare[1]].data(), Hessian, grad,
				ee_hessian_record_index.data() + 5 * (prefix_sum[j] + ee_pair_compressed_record[index + 4]), stiffness, d_hat_2,
				edge_length[j], rest_edge_length[ee_pair_compressed_record[index + 2]][ee_pair_compressed_record[index + 3]])) {
				memcpy(ee_hessian_record.data() + 144 * (prefix_sum[j] + ee_pair_compressed_record[index + 4]), Hessian.data(), Hessian.size() << 3);
				memcpy(ee_grad_record.data() + 12 * (prefix_sum[j] + ee_pair_compressed_record[index + 4]), grad.data(), grad.size() << 3);

				//if (ee_pair_compressed_record[index] == 0 && ee_pair_compressed_record[index + 1] == 28) {
				//	if (ee_pair_compressed_record[index + 2] == 1 && (ee_pair_compressed_record[index + 3] == 102 || ee_pair_compressed_record[index + 3] == 106)) {
				//		int* add = ee_hessian_record_index.data() + 5 * (prefix_sum[j] + ee_pair_compressed_record[index + 4]);
				//		std::cout << "index " << ee_pair_compressed_record[index + 3] << std::endl;
				//		std::cout << add[0] << " " << add[1] << " " << add[2] << " " << add[3] << " " << add[4] << std::endl;
				//	}
				//}


			}
		}
	}
	int vertex_prefix_sum_this_obj;

	if (has_collider)
	{
		unsigned int* element_collide_with_collider;
		unsigned int* pair_by_element;
		unsigned int* pair_by_element_num;

		for (int tet_obj_no = 0; tet_obj_no < tetrahedron->size(); ++tet_obj_no) {

			i = tet_obj_no + cloth->size();
			vertex_pos = vertex_position[i];
			//vt collider
			vertex_prefix_sum_this_obj = vertex_num_on_surface_prefix_sum[i];
			start = vertex_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No];
			end = vertex_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No + 1];
			vertex_index_on_surface = this->general_index_to_surface_index[i];
			element_collide_with_collider = vertex_index_collide_with_collider[i].data();
			pair_by_element = vertex_obj_triangle_collider_pair_by_vertex[i];
			pair_by_element_num = vertex_obj_triangle_collider_num_record[i];
			for (int k = start; k < end; k++) {
				j = element_collide_with_collider[k];
				surface_index = vertex_index_on_surface[j];
				computeVTColliderHessian(pair_by_element
					+ close_vt_collider_pair_num * surface_index,
					pair_by_element_num[surface_index], d_hat_2,
					vertex_pos[j].data(), vt_colldier_hessian_record.data() + (vertex_prefix_sum_this_obj + surface_index) * 9,
					stiffness, vt_colldier_grad_record.data() + 3 * (surface_index + vertex_prefix_sum_this_obj));
			}

			//ee collider
			start = edge_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No];
			end = edge_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No + 1];
			pair_by_element = edge_edge_collider_pair_by_edge[i];
			pair_by_element_num = edge_edge_collider_pair_num_record[i];
			element_collide_with_collider = edge_index_collide_with_collider[i].data();
			edge_vertex = this->edge_vertices[i];
			prefix_sum_collider = edge_edge_collider_num_record_prefix_sum[i];
			edge_length = rest_edge_length[i];
			for (int k = start; k < end; k++) {
				j = element_collide_with_collider[k];
				computeEEColliderHessian(pair_by_element + close_ee_collider_pair_num * j,
					pair_by_element_num[j],
					d_hat_2, vertex_pos[edge_vertex[j << 1]].data(),
					vertex_pos[edge_vertex[(j << 1) + 1]].data(),
					ee_collider_hessian_record.data() + 36 * prefix_sum_collider[j],
					ee_collider_hessian_record_index.data() + 3 * prefix_sum_collider[j], stiffness, ee_collider_grad_record.data() + 6 * prefix_sum_collider[j],
					edge_length[j]);

			}

			//tv collider
			start = triangle_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No];
			end = triangle_index_collide_with_collider_start_per_thread[(thread_num + 1) * i + thread_No + 1];
			pair_by_element = triangle_vertex_collider_pair_by_triangle[i];
			pair_by_element_num = triangle_vertex_collider_pair_num_record[i];
			element_collide_with_collider = triangle_index_collide_with_collider[i].data();
			prefix_sum_collider = triangle_vertex_collider_num_record_prefix_sum[i];

			for (int k = start; k < end; k++) {
				j = element_collide_with_collider[k];
				triangle_vertex = triangle_indices[i][j].data();
				computeTVColliderHessian(pair_by_element + close_tv_collider_pair_num * j,
					pair_by_element_num[j], d_hat_2, vertex_pos[triangle_vertex[0]].data(),
					vertex_pos[triangle_vertex[1]].data(), vertex_pos[triangle_vertex[2]].data(),
					tv_colldier_hessian_record.data() + 81 * prefix_sum_collider[j], stiffness, tv_colldier_grad_record.data() + 9 * prefix_sum_collider[j],
					tv_colldier_hessian_record_index.data() + 4 * prefix_sum_collider[j]);
			}
		}
	}

	if (has_floor) {
		int end_per_thread;
		unsigned int* surface_index_to_global;
		for (int tet_obj_no = 0; tet_obj_no < tetrahedron->size(); ++tet_obj_no) {
			i = tet_obj_no + cloth->size();
			vertex_pos = vertex_position[i];
			vertex_prefix_sum_this_obj = vertex_num_on_surface_prefix_sum[i];
			end_per_thread = vertex_index_start_per_thread[i][thread_No + 1];
			surface_index_to_global = this->vertex_index_on_surface[i];
			for (int j = vertex_index_start_per_thread[i][thread_No]; j < end_per_thread; ++j) {
				computeFloorHessian(d_hat_2, stiffness, floor_value, floor_hessian_record[vertex_prefix_sum_this_obj + j], floor_grad_record[vertex_prefix_sum_this_obj + j],
					vertex_pos[surface_index_to_global[j]][floor_dimension]);
			}
		}
	}
}

//UPDATE_COLLISION_HESSIAN_COLOR
void Collision::computeHessianPerThread(int thread_No, int color_No)
{
	int color_group_index = 0;
	color_group_index = *inner_iteration_number % tet_color_groups[cloth->size()]->size();
	if (color_No < tet_color_groups[cloth->size()]->data()[color_group_index].size() - 1) {
		computeHessianPreviousThread(thread_No, color_No);
	}
	else {
		computeLastColorHessianPerThread(thread_No, color_No);
	}

}


void Collision::computeVTHessian(unsigned int* VT, unsigned int num, double d_hat_2, double* vertex_position_,  double* hessian_record, int* hessian_record_index, double stiffness,  double* grad_record)
{
	int* triangle_vertex;
	MatrixXd Hessian; VectorXd grad;
	for (int i = 0; i < num; i+=2) {
		triangle_vertex = triangle_indices[VT[i]][VT[i + 1]].data();
		if (second_order_constraint.computeBarrierVTGradientHessian(Hessian, grad, vertex_position_, vertex_position[VT[i]][triangle_vertex[0]].data(),
			vertex_position[VT[i]][triangle_vertex[1]].data(), vertex_position[VT[i]][triangle_vertex[2]].data(), d_hat_2, hessian_record_index, stiffness)) {
			memcpy(hessian_record, Hessian.data(), Hessian.size() << 3);
			memcpy(grad_record, grad.data(), grad.size() << 3);
		}
		hessian_record += 144;
		grad_record += 12;
		hessian_record_index += 5;
	}	
}


void Collision::computeTVHessian(unsigned int* TV, unsigned int num, double d_hat_2, double* vertex_position_0,
	double* vertex_position_1, double* vertex_position_2, double stiffness)
{
	MatrixXd Hessian; VectorXd grad;
	int surface_index;
	for (int i = 0; i < num; i += 3) {
		//if (!vertex_belong_to_color_group[TV[i]][TV[i + 1]]) {
			if (TV[i] < cloth->size()) {
				surface_index = TV[i + 1];
			}
			else {
				surface_index = general_index_to_surface_index[TV[i]][TV[i + 1]];
			}
			if (second_order_constraint.computeBarrierVTGradientHessian(Hessian, grad, vertex_position[TV[i]][TV[i + 1]].data(),
				vertex_position_0, vertex_position_1, vertex_position_2, d_hat_2, vt_hessian_record_index.data()+5*(vertex_triangle_pair_num_record_prefix_sum[TV[i]][surface_index]
					+TV[i+2]), stiffness)) {
				memcpy(vt_hessian_record.data()+ 144 * (vertex_triangle_pair_num_record_prefix_sum[TV[i]][surface_index]
					+ TV[i + 2]), Hessian.data(), Hessian.size() << 3);
				memcpy(vt_grad_record.data()+ 12 * (vertex_triangle_pair_num_record_prefix_sum[TV[i]][surface_index]
					+ TV[i + 2]), grad.data(), grad.size() << 3);
			}
		//}
	}
}






void Collision::computeTVColliderHessian(unsigned int* TV, unsigned int num, double d_hat_2, double* vertex_position_0,
	double* vertex_position_1, double* vertex_position_2,
	double* hessian_record,
	double stiffness, double* grad_record, int* hessian_record_index)
{
	int vertex_in_pair[5];
	MatrixXd Hessian; VectorXd grad;
	for (int i = 0; i < num; i += 2) {
		if (second_order_constraint.computeBarrierVTGradientHessian(Hessian, grad, vertex_position_collider[TV[i]][TV[i + 1]].data(),
			vertex_position_0, vertex_position_1, vertex_position_2, d_hat_2, vertex_in_pair, stiffness)) {
			setTVColliderHessianFix(Hessian, grad, hessian_record, vertex_in_pair, grad_record, hessian_record_index);
		}
		else {
			hessian_record_index[0] = 0;
		}
		hessian_record += 81;
		grad_record += 9;
		hessian_record_index += 4;
	}
}


void Collision::computeFloorHessian(double d_hat, double stiffness, double floor_value, double& hessian, double& grad,  double position)
{
	double distance = (position- floor_value) * (position - floor_value);
	if (distance > d_hat_2) {
		return;
	}
	double h, g;
	barrierGradHessian(distance, d_hat_2, g, h);
	double grad_d;
	grad_d = 2 * (position - floor_value);

	grad = stiffness * g * grad_d;
	hessian = stiffness * (h * grad_d * grad_d + g+g);
}


void Collision::computeVTColliderHessian(unsigned int* VT, unsigned int num, double d_hat_2, double* vertex_position_, 
	double* hessian_record, 
	double stiffness, double* grad_record)
{
	int* triangle_vertex;
	MatrixXd Hessian; VectorXd grad;
	Matrix3d Hessian_total; Vector3d grad_total;
	Hessian_total.setZero();
	grad_total.setZero();
	int hessian_record_index[5];
	for (int i = 0; i < num; i+=2) {
		triangle_vertex = triangle_indices_collider[VT[i]][VT[i + 1]].data();
		if (second_order_constraint.computeBarrierVTGradientHessian(Hessian, grad, vertex_position_, vertex_position_collider[VT[i]][triangle_vertex[0]].data(),
			vertex_position_collider[VT[i]][triangle_vertex[1]].data(), vertex_position_collider[VT[i]][triangle_vertex[2]].data(), d_hat_2, hessian_record_index, stiffness)) {
			Hessian_total += Hessian.block<3, 3>(0, 0);
			grad_total += grad.segment(0, 3);
		}
	}
	memcpy(hessian_record, Hessian_total.data(), 72);
	memcpy(grad_record, grad_total.data(), 24);
}


void Collision::computeEEHessian(unsigned int* EE, unsigned int num, double d_hat_2, double* ea0, double* ea1,
	double* hessian_record, int* hessian_record_index, double stiffness, double* grad_record, double edge_length_0)
{
	unsigned int* edge_vertex;
	MatrixXd Hessian; VectorXd grad;

	for (int i = 0; i < num; i+=3) {
		edge_vertex = edge_vertices[EE[i]] + (EE[i + 1] << 1);
		if(second_order_constraint.computeBarrierEEGradientHessian(ea0, ea1, vertex_position[EE[i]][edge_vertex[0]].data(),
			vertex_position[EE[i]][edge_vertex[1]].data(),Hessian,grad, hessian_record_index,stiffness,d_hat_2,
			edge_length_0, rest_edge_length[EE[i]][EE[i+1]])) {
			memcpy(hessian_record, Hessian.data(), Hessian.size() << 3);
			memcpy(grad_record, grad.data(), grad.size() << 3);
		}		
		hessian_record += 144;
		grad_record += 12;
		hessian_record_index += 5;
	}
}


void Collision::setTVColliderHessianFix(MatrixXd& Hessian, VectorXd& grad, double* hessian_record,
	int* curent_hessian_record_local, double* grad_record, int* hessian_record_index)
{
	int size = grad.size() - 3;
	for (int i = 0; i < size; ++i) {
		memcpy(hessian_record + i * size, Hessian.data() + grad.size()* (i+3)  + 3, size << 3);
	}
	memcpy(grad_record, grad.data() + 3, size << 3);

	memcpy(hessian_record_index + 1, curent_hessian_record_local + 2, (curent_hessian_record_local[0] - 1) << 2);
	hessian_record_index[0]=curent_hessian_record_local[0] - 1;
}


void Collision::setEEColliderHessianFix(MatrixXd& Hessian, VectorXd& grad, double* hessian_record,
	int* curent_hessian_record_local, double* grad_record, int* hessian_record_index)
{
	if (curent_hessian_record_local[0] < 4) {
		if (curent_hessian_record_local[1] == 0) {
			memcpy(hessian_record, Hessian.data(), 24);
			memcpy(hessian_record +3, Hessian.data()+Hessian.cols(), 24);
			memcpy(hessian_record + 6, Hessian.data() + Hessian.cols() * 2, 24);

			memcpy(grad_record, grad.data(), 24);

			hessian_record_index[0] = 1;
			hessian_record_index[1] = 0;
		}
		else if (curent_hessian_record_local[1] == 1) {
			memcpy(hessian_record, Hessian.data(), 24);
			memcpy(hessian_record + 3, Hessian.data() + Hessian.cols(), 24);
			memcpy(hessian_record + 6, Hessian.data() + Hessian.cols() * 2, 24);

			memcpy(grad_record, grad.data(), 24);

			hessian_record_index[0] = 1;
			hessian_record_index[1] = 1;
		}
		else if (curent_hessian_record_local[1] > 1) {
			for (int i = 0; i < 6; ++i) {
				memcpy(hessian_record + 6 * i, Hessian.data()+9*i +30, 48);
			}

			memcpy(grad_record, grad.data()+3, 48);

			hessian_record_index[0] = 2;
			hessian_record_index[1] = 0;
			hessian_record_index[2] = 1;
		}
	}
	else {
		for (int i = 0; i < 6; ++i) {
			memcpy(hessian_record + 6 * i, Hessian.data() + 12 * i, 48);
		}

		memcpy(grad_record, grad.data(), 48);

		hessian_record_index[0] = 2;
		hessian_record_index[1] = 0;
		hessian_record_index[2] = 1;
	}

}


void Collision::computeEEColliderHessian(unsigned int* EE, unsigned int num, double d_hat_2, double* ea0, double* ea1,
	double* hessian_record, int* hessian_record_index, double stiffness, double* grad_record, 
	double edge_length_0)
{
	unsigned int* edge_vertex;
	MatrixXd Hessian; VectorXd grad;

	int hessian_record_[5];

	for (int i = 0; i < num; i += 2) {
		edge_vertex = collider_edge_vertices[EE[i]] + (EE[i + 1] << 1);
		if (second_order_constraint.computeBarrierEEGradientHessian(ea0, ea1, vertex_position_collider[EE[i]][edge_vertex[0]].data(),
			vertex_position_collider[EE[i]][edge_vertex[1]].data(), Hessian, grad, hessian_record_, stiffness, d_hat_2,
			edge_length_0, rest_edge_length_collider[EE[i]][EE[i + 1]])) {

			setEEColliderHessianFix(Hessian, grad, hessian_record, hessian_record_, 
				grad_record, hessian_record_index);
		}
		else {
			hessian_record_index[0] = 0;
		}
		hessian_record += 36;
		grad_record += 6;
		hessian_record_index += 3;
	}
}



void Collision::findVT_ClosePair(int thread_No)
{
	unsigned int start,end;
	unsigned int* vertex_triangle;
	unsigned int* vertex_triangle_num;
	unsigned int* vertex_index_on_surface_;

	bool* pair_exist;

	for (int i = 0; i < total_obj_num; ++i) {
		start= vertex_index_start_per_thread[i][thread_No];
		end = vertex_index_start_per_thread[i][thread_No+1];
		vertex_triangle = spatial_hashing.vertex_triangle_pair_by_vertex[i];
		vertex_triangle_num = spatial_hashing.vertex_triangle_pair_num_record[i];
		pair_exist = spatial_hashing.is_used_vertex_triangle_pair_by_vertex[i];

		if (i < cloth->size()) {
			for (int j = start; j < end; ++j) {
				findVT_ClosePairSingleVertex(vertex_position[i][j].data(), vertex_triangle + estimate_coeff_for_vt_pair_num * j,
					vertex_triangle_num[j], vertex_position, triangle_indices, vertex_triangle_pair_by_vertex[i] + j * close_vt_pair_num,
					vertex_triangle_pair_num_record[i][j], pair_exist + estimate_coeff_for_vt_pair_num_exist * j);
			}
		}
		else {
			vertex_index_on_surface_ = vertex_index_on_surface[i];
			for (int j = start; j < end; ++j) {				
				findVT_ClosePairSingleVertex(vertex_position[i][vertex_index_on_surface_[j]].data(), vertex_triangle + estimate_coeff_for_vt_pair_num * j,
					vertex_triangle_num[j], vertex_position, triangle_indices, vertex_triangle_pair_by_vertex[i] + j * close_vt_pair_num,
					vertex_triangle_pair_num_record[i][j], pair_exist + estimate_coeff_for_vt_pair_num_exist * j);
			}
		}

	
	}
}


void Collision::initialPairByElement()
{
	for (int i = 0; i < total_obj_num; ++i) {
		memset(vertex_triangle_pair_num_record[i], 0, 4 * vertex_index_start_per_thread[i][thread_num]);
		memset(triangle_vertex_pair_num_record[i], 0, 4 * mesh_struct[i]->triangle_indices.size());
		memset(edge_edge_pair_number_record[i], 0, 4 * edge_index_start_per_thread[i][thread_num]);
		if (has_collider) {
			memset(vertex_obj_triangle_collider_num_record[i], 0, 4 * vertex_index_start_per_thread[i][thread_num]);
			memset(triangle_vertex_collider_pair_num_record[i], 0, 4 * triangle_index_start_per_thread[i][thread_num]);
			memset(edge_edge_collider_pair_num_record[i], 0, 4 * edge_index_start_per_thread[i][thread_num]);
		}		
	}
}

void Collision::findTV_ClosePair()
{

	for (int i = 0; i < total_obj_num; ++i) {
		findTV_ClosePair(i, vertex_triangle_pair_by_vertex[i], vertex_triangle_pair_num_record[i], 
			triangle_vertex_pair_by_triangle,triangle_vertex_pair_num_record, close_vt_pair_num, 
			close_tv_pair_num,i>=cloth->size());
	}
}

void Collision::findTV_ClosePair(int obj_No, unsigned int* vt_pair_initial, unsigned int* vt_pair_num, unsigned int** tv_pair,
	unsigned int** tv_pair_num, unsigned int total_length_every_element_vt, unsigned int total_length_every_element_tv, bool is_tet)
{
	unsigned int vertex_end, i;
	unsigned int total_num;

	unsigned int* temp_address;

	unsigned int* vt_pair;
	vertex_end = vertex_index_start_per_thread[obj_No][thread_num];
	for (int ind = 0; ind < vertex_end; ++ind) {
		if (is_tet) {
			i = vertex_index_on_surface[obj_No][ind];
		}
		else {
			i = ind;
		}
		total_num = vt_pair_num[ind];
		vt_pair = vt_pair_initial + ind * total_length_every_element_vt;
		for (int j = 0; j < total_num; j += 2) {
			temp_address = tv_pair[vt_pair[j]] + vt_pair[j + 1] * total_length_every_element_tv + tv_pair_num[vt_pair[j]][vt_pair[j + 1]];
			*(temp_address++) = obj_No;
			*(temp_address++) = i;
			*temp_address = (j>>1);
			tv_pair_num[vt_pair[j]][vt_pair[j + 1]] += 3;
		}
	}
}






void Collision::findEE_ClosePair(int thread_No)
{
	unsigned int start;
	unsigned int end;
	unsigned int* edge_index;
	unsigned int* edge_num;
	unsigned int* edge_vertex;
	bool* pair_exist;
	for (int i = 0; i < total_obj_num; ++i) {
		start = edge_index_start_per_thread[i][thread_No];
		end = edge_index_start_per_thread[i][thread_No + 1];
		edge_index = spatial_hashing.edge_edge_pair_by_edge[i];
		edge_num = spatial_hashing.edge_edge_pair_num_record[i];
		edge_vertex = edge_vertices[i];
		pair_exist = spatial_hashing.is_used_edge_edge_pair_by_edge[i];
		for (int j = start; j < end; ++j) {
			findEE_ClosePairSingleEdge(vertex_position[i][edge_vertex[j << 1]].data(),
				vertex_position[i][edge_vertex[(j << 1) + 1]].data(),
				edge_index + estimate_coeff_for_ee_pair_num * j,
				edge_num[j], vertex_position, edge_vertices,
				edge_edge_pair_by_edge[i] + j * close_ee_pair_num,
				edge_edge_pair_number_record[i][j], pair_exist + estimate_coeff_for_ee_pair_num_exist * j,3);
		}
	}
}

void Collision::findEE_ClosePair_ReverseOrder()
{
	unsigned int ee_num;
	unsigned int* edge_index;
	unsigned int* edge_num;
	unsigned int* edge_vertex;
	bool* pair_exist;

	for (int i = 0; i < total_obj_num; ++i) {
		ee_num = edge_index_start_per_thread[i][thread_num];
		edge_index = spatial_hashing.edge_edge_pair_by_edge[i];
		edge_num = spatial_hashing.edge_edge_pair_num_record[i];
		edge_vertex = edge_vertices[i];
		pair_exist = spatial_hashing.is_used_edge_edge_pair_by_edge[i];
		for (int j = 0; j < ee_num; ++j) {
			findEE_ClosePairSingleEdgeByReverse(i,j,edge_edge_pair_by_edge[i] + j * close_ee_pair_num,
				edge_edge_pair_number_record[i][j]);
		}
	}
}


void Collision::findEE_ColliderClosePair(int thread_No)
{
	unsigned int start,end;
	unsigned int* edge_index;
	unsigned int* edge_num;
	unsigned int* edge_vertex;
	bool* pair_exist;
	for (int i = 0; i < total_obj_num; ++i) {
		start = edge_index_start_per_thread[i][thread_No];
		end = edge_index_start_per_thread[i][thread_No + 1];
		edge_index = spatial_hashing.edge_obj_edge_collider_pair_by_edge[i];
		edge_num = spatial_hashing.edge_obj_edge_collider_num_record[i];
		edge_vertex = edge_vertices[i];
		pair_exist = spatial_hashing.is_used_edge_obj_edge_collider_pair_by_edge[i];

		for (int j = start; j < end; ++j) {
			findEE_ClosePairSingleEdge(vertex_position[i][edge_vertex[j << 1]].data(),
				vertex_position[i][edge_vertex[(j << 1) + 1]].data(),
				edge_index + estimate_coeff_for_ee_collider_pair_num * j,
				edge_num[j], vertex_position_collider, collider_edge_vertices,
				edge_edge_collider_pair_by_edge[i] + j * close_ee_collider_pair_num,
				edge_edge_collider_pair_num_record[i][j],
				pair_exist + estimate_coeff_for_ee_collider_pair_num_exist * j,2);
		}
	}

	
	//if (edge_edge_collider_pair_num_record[0][0] != 0) {
	//	std::cout << "record close edge edge collider num " << edge_edge_collider_pair_num_record[0][0] << std::endl;
	//	for (int i = 0; i < edge_edge_collider_pair_num_record[0][0]; ++i) {
	//		std::cout << edge_edge_collider_pair_by_edge[0][i] << " ";
	//	}
	//	std::cout << std::endl;
	//}
}


void Collision::testColliderPair()
{
	unsigned int ee_num;
	unsigned int* edge_index;
	unsigned int* edge_num;
	unsigned int* start;
	unsigned int num;

	unsigned int* edge_vertex;
	unsigned int* edge_vertex_1;
	for (int k = 0; k < total_obj_num; ++k) {
		ee_num = edge_index_start_per_thread[k][thread_num];
		edge_index = edge_edge_collider_pair_by_edge[k];
		edge_num = edge_edge_collider_pair_num_record[k];
		for (int j = 0; j < ee_num; ++j) {
			start = edge_index + close_ee_collider_pair_num * j;
			if (edge_num[j] > 0) {
				std::cout << "num " << edge_num[j] << std::endl;
			}
			edge_vertex = edge_vertices[k]+(j << 1);
			for (int i = 0; i < edge_num[j]; i += 2) {
				edge_vertex_1 = collider_edge_vertices[start[i]] + (start[i + 1] << 1);
				double distance;
				distance = CCD::internal::edgeEdgeDistanceUnclassified(vertex_position[k][*edge_vertex].data(), vertex_position[k][*(edge_vertex + 1)].data(),
					vertex_position_collider[start[i]][edge_vertex_1[0]].data(), vertex_position_collider[start[i]][edge_vertex_1[1]].data());
				if (distance < d_hat_2) {

				}
				else {
					std::cout << "ee collider distance too far away from "<<distance<<" "<< d_hat_2 << std::endl;
				}
			}
		}
	}

	int* triangle_vertex;

	for (int k = 0; k < total_obj_num; ++k) {
		ee_num = triangle_index_start_per_thread[k][thread_num];
		edge_index = triangle_vertex_collider_pair_by_triangle[k];
		edge_num = triangle_vertex_collider_pair_num_record[k];
		for (int j = 0; j < ee_num; ++j) {
			start = edge_index + close_tv_collider_pair_num * j;
			if (edge_num[j] > 0) {
				std::cout << "num " << edge_num[j] << std::endl;
			}
			triangle_vertex = triangle_indices[k][j].data();
			for (int i = 0; i < edge_num[j]; i += 2) {
				double distance = CCD::internal::pointTriangleDistanceUnclassified(vertex_position_collider[start[i]][start[i + 1]].data(),
					vertex_position[k][triangle_vertex[0]].data(), vertex_position[k][triangle_vertex[1]].data(),
					vertex_position[k][triangle_vertex[2]].data());
				if (distance < d_hat_2) {

				}
				else {
					std::cout << "tv collider distance too far away from "<< distance<<" "<< d_hat_2 << std::endl;
				}
			}
		}
	}

}


void Collision::findTV_ColliderClosePair(int thread_No)
{
	unsigned int start,end;
	unsigned int* vertex_index;
	unsigned int* vertex_num;
	std::array<int,3>* triangle_vertex;
	bool* pair_exist;
	for (int i = 0; i < total_obj_num; ++i) {
		start = triangle_index_start_per_thread[i][thread_No];
		end = triangle_index_start_per_thread[i][thread_No + 1];
		vertex_index = spatial_hashing.triangle_obj_vertex_collider_pair_by_triangle[i];
		vertex_num = spatial_hashing.triangle_obj_vertex_collider_num_record[i];
		triangle_vertex = triangle_indices[i];
		pair_exist = spatial_hashing.is_used_triangle_obj_vertex_collider_pair_by_triangle[i];
		for (int j = start; j < end; ++j) {
			findTV_ClosePairSingleTriangle(vertex_position[i][triangle_vertex[j][0]].data(),
				vertex_position[i][triangle_vertex[j][1]].data(), vertex_position[i][triangle_vertex[j][2]].data(),
				vertex_index + estimate_coeff_for_tv_collider_pair_num * j,
				vertex_num[j], vertex_position_collider, 
				triangle_vertex_collider_pair_by_triangle[i] + j * close_tv_collider_pair_num,
				triangle_vertex_collider_pair_num_record[i][j],
				pair_exist + estimate_coeff_for_tv_collider_pair_num_exist * j);
		}
	}
}



void Collision::findVT_ColliderClosePair(int thread_No)
{
	unsigned int start,end;
	unsigned int* vertex_triangle;
	unsigned int* vertex_triangle_num;
	unsigned int* vertex_index_on_surface_;
	bool* pair_exist;
	for (int i = 0; i < total_obj_num; ++i) {
		start = vertex_index_start_per_thread[i][thread_No];
		end = vertex_index_start_per_thread[i][thread_No+1];
		vertex_triangle = spatial_hashing.vertex_obj_triangle_collider_pair_by_vertex[i];
		vertex_triangle_num = spatial_hashing.vertex_obj_triangle_collider_num_record[i];
		pair_exist = spatial_hashing.is_used_vertex_obj_triangle_collider_pair_by_vertex[i];
		if (i < cloth->size()) {
			for (int j = start; j < end; ++j) {
				findVT_ClosePairSingleVertex(vertex_position[i][j].data(), vertex_triangle + estimate_coeff_for_vt_collider_pair_num * j,
					vertex_triangle_num[j], vertex_position_collider, triangle_indices_collider,
					vertex_obj_triangle_collider_pair_by_vertex[i] + j * close_vt_collider_pair_num,
					vertex_obj_triangle_collider_num_record[i][j],
					pair_exist + estimate_coeff_for_vt_collider_pair_num_exist * j);
			}
		}
		else {
			vertex_index_on_surface_ = vertex_index_on_surface[i];
			for (int j = start; j < end; ++j) {
				findVT_ClosePairSingleVertex(vertex_position[i][vertex_index_on_surface_[j]].data(), vertex_triangle + estimate_coeff_for_vt_collider_pair_num * j,
					vertex_triangle_num[j], vertex_position_collider, triangle_indices_collider,
					vertex_obj_triangle_collider_pair_by_vertex[i] + j * close_vt_collider_pair_num,
					vertex_obj_triangle_collider_num_record[i][j],
					pair_exist + estimate_coeff_for_vt_collider_pair_num_exist * j);
			}
		}
	}
}


void Collision::findVT_ClosePairSingleVertex(double* current_position, unsigned int* trianlge_index, 
	unsigned int triangle_num, std::vector<std::array<double,3>*>& position, 
	std::vector<std::array<int, 3>*>& triangle_vertex_index,
	unsigned int* close_triangle_index, unsigned int& close_triangle_num, bool* is_pair_exist)
{
	int* triangle_vertex;
	for (int i = 0; i < triangle_num; i+=2) {
		triangle_vertex = triangle_vertex_index[trianlge_index[i]][trianlge_index[i + 1]].data();
		double distance = CCD::internal::pointTriangleDistanceUnclassified(current_position, position[trianlge_index[i]][triangle_vertex[0]].data(),
			position[trianlge_index[i]][triangle_vertex[1]].data(), position[trianlge_index[i]][triangle_vertex[2]].data());

		if (distance < d_hat_2) {
			memcpy(close_triangle_index, trianlge_index + i, 8);
			close_triangle_index += 2;
			close_triangle_num += 2;			

			if (!*is_pair_exist) {
				*is_pair_exist = true;
			}
		}
		is_pair_exist++;
	}
}


void Collision::findEE_ClosePairSingleEdgeByReverse(int obj_index, int edge_index, unsigned int* close_edge_index, unsigned int& close_edge_num)
{
	unsigned int* edge_pair_address;
	
	for (int i = 0; i < close_edge_num; i += 3) {
		if (close_edge_index[i] < obj_index) {
			return;
		}
		if (close_edge_index[i] == obj_index && close_edge_index[i+1] < edge_index) {
			return;
		}
		edge_pair_address = edge_edge_pair_by_edge[close_edge_index[i]] + 
			close_edge_index[i + 1] * close_ee_pair_num + edge_edge_pair_number_record[close_edge_index[i]][close_edge_index[i + 1]];
		edge_edge_pair_number_record[close_edge_index[i]][close_edge_index[i + 1]] += 3;
		*edge_pair_address = obj_index;
		*(edge_pair_address + 1) = edge_index;
		*(edge_pair_address + 2) = i;
	}
}


void Collision::findEE_ClosePairSingleEdge(double* current_position_a0, double* current_position_a1, 
	unsigned int* edge_index, unsigned int edge_num, std::vector<std::array<double, 3>*>& position,
	std::vector<unsigned int*>& edge_vertex_index,
	unsigned int* close_edge_index, unsigned int& close_edge_num, bool* is_pair_exist, int move_size)
{
	unsigned int* edge_vertex;
	for (int i = 0; i < edge_num; i += 2) {
		edge_vertex = edge_vertex_index[edge_index[i]]+ (edge_index[i+1]<<1);
		if (CCD::internal::edgeEdgeDistanceUnclassified(current_position_a0, current_position_a1, 
			position[edge_index[i]][edge_vertex[0]].data(),position[edge_index[i]][edge_vertex[1]].data()) < d_hat_2) {
			memcpy(close_edge_index, edge_index + i, 8);
			close_edge_index += move_size;
			close_edge_num += move_size;
			if (!*is_pair_exist) {
				*is_pair_exist = true;
			}
		}
		is_pair_exist++;
	}
}

void Collision::findTV_ClosePairSingleTriangle(double* current_position_a0, double* current_position_a1, double* current_position_a2,
	unsigned int* vertex_index, unsigned int vertex_num, std::vector<std::array<double, 3>*>& position,
	unsigned int* close_vertex_index, unsigned int& close_vertex_num, bool* is_pair_exist)
{
	for (int i = 0; i < vertex_num; i += 2) {
		if (CCD::internal::pointTriangleDistanceUnclassified(position[vertex_index[i]][vertex_index[i+1]].data(), current_position_a0, current_position_a1,
			current_position_a2) < d_hat_2) {
			memcpy(close_vertex_index, vertex_index + i, 8);
			close_vertex_index += 2;
			close_vertex_num += 2;
			if (!*is_pair_exist) {
				*is_pair_exist = true;
			}
		}
		is_pair_exist++;
	}
}


//COLLISION_CONSTRAINT_IPC
void Collision::collisionConstraintIPC(int thread_No)
{
	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
	target_pos->initial();
	checkTargetPosSize(thread_No);

	point_triangle_target_pos_index[thread_No][0] = 0;
	edge_edge_target_pos_index[thread_No][0] = 0;
	if (has_collider) {
		point_triangle_collider_target_pos_index[thread_No][0] = 0;
		//point_collider_triangle_target_pos_index[thread_No][0] = 0;
		//edge_edge_collider_target_pos_index[thread_No][0] = 0;
	}

	pointTriangleResponseForIPC(thread_No, target_pos);
	edgeEdgeResponseForIPC(thread_No, target_pos);

	if (has_collider) {
		pointTriangleColliderResponseForIPC(thread_No, target_pos);
	}
	if (floor->exist) {
		floorCollisionVertexForIPC(thread_No);
	}
}



//COLLISION_CONSTRAINT
void Collision::collisionConstraint(int thread_No)
{
	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
	target_pos->initial();
	checkTargetPosSize(thread_No);
	//target_position_index[thread_No][0] = 0;
	point_triangle_target_pos_index[thread_No][0] = 0;
	edge_edge_target_pos_index[thread_No][0] = 0;
	if (has_collider) {
		point_triangle_collider_target_pos_index[thread_No][0] = 0;
		//point_collider_triangle_target_pos_index[thread_No][0] = 0;
		//edge_edge_collider_target_pos_index[thread_No][0] = 0;
	}

	pointTriangleResponse(thread_No, target_pos);
	edgeEdgeResponse(thread_No, target_pos);	

	if (has_collider) {
		pointTriangleColliderResponse(thread_No, target_pos);
		//pointColliderTriangleResponse(thread_No, target_pos);	
		//edgeEdgeColliderResponse(thread_No, target_pos);
	}
	if (floor->exist) {
		floorCollisionVertex(thread_No);
	}
}



//RE_COLLISION_CONSTRAINT_IPC
void Collision::re_collisionConstraintIPC(int thread_No)
{
	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
	target_pos->partialInitial();
	re_pointTriangleResponseForIPC(thread_No, target_pos);
	re_edgeEdgeResponseForIPC(thread_No, target_pos);
	if (has_collider) {
		re_pointTriangleColliderResponseForIPC(thread_No, target_pos);
		//re_pointColliderTriangleResponse(thread_No, target_pos);
		//re_edgeEdgeColliderResponse(thread_No, target_pos);
	}
	if (floor->exist) {
		re_FloorCollisionVertexForIPC(thread_No);
	}

}


//RE_COLLISION_CONSTRAINT
void Collision::re_collisionConstraint(int thread_No)
{
	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
	target_pos->partialInitial();
	re_pointTriangleResponse(thread_No, target_pos);
	re_edgeEdgeResponse(thread_No, target_pos);	
	if (has_collider) {
		re_pointTriangleColliderResponse(thread_No, target_pos);
		//re_pointColliderTriangleResponse(thread_No, target_pos);
		//re_edgeEdgeColliderResponse(thread_No, target_pos);
	}
	if (floor->exist) {
		re_FloorCollisionVertex(thread_No);
	}

}



void Collision::collisionEnergy()
{
	////std::cout <<"energy 0 "<< obj_target_pos.collision_energy << std::endl;
	obj_target_pos.collision_energy = 0;
	thread->assignTask(this, COMPUTE_COLLISION_ENERGY);
	for (int i = 0; i < thread_num; ++i) {
		obj_target_pos.collision_energy += obj_target_pos_per_thread[i].collision_energy;
	}
	////std::cout << "energy 1 " << obj_target_pos.collision_energy << std::endl;
}

//COMPUTE_COLLISION_ENERGY
void Collision::collisionEnergy(int thread_No)
{
	//double energy = 0;
	//if (target_position_element_start_per_thread[(thread_No + 1) << 1] > target_position_element_start_per_thread[thread_No << 1]) {
	//	computeEnergy(target_position_element_start_per_thread[thread_No << 1], target_position_element_start_per_thread[(thread_No << 1) + 1],
	//		target_position_index[target_position_element_start_per_thread[thread_No << 1]][0], energy);
	//	for (int i = target_position_element_start_per_thread[thread_No << 1] + 1;
	//		i < target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
	//		computeEnergy(i, 0,
	//			target_position_index[i][0], energy);
	//	}
	//	computeEnergy(target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
	//		target_position_element_start_per_thread[(thread_No << 1) + 3], energy);
	//}
	//else {
	//	computeEnergy(target_position_element_start_per_thread[thread_No << 1], target_position_element_start_per_thread[(thread_No << 1) + 1],
	//		target_position_element_start_per_thread[(thread_No << 1) + 3], energy);
	//}
	//obj_target_pos_per_thread[thread_No].collision_energy = 0.5*energy;
}





void Collision::checkTargetPosSize(int thread_No)
{
	//unsigned int target_pos_size = (ave_pair_num[0]+ ave_pair_num[1]) << 2;
	//if (has_collider) {
	//	target_pos_size += ave_pair_num[2];
	//	target_pos_size += 3 * ave_pair_num[3];
	//	target_pos_size += 2 * ave_pair_num[4];
	//}
	if (point_triangle_target_pos_index[thread_No].size() < ave_pair_num[0]) {
		point_triangle_target_pos_index[thread_No].resize(ave_pair_num[0]+1);
	}
	if (point_triangle_collider_target_pos_index[thread_No].size() < ave_pair_num[2]) {
		point_triangle_collider_target_pos_index[thread_No].resize(ave_pair_num[2]+1);
	}
	//if (point_collider_triangle_target_pos_index[thread_No].size() < ave_pair_num[3]) {
	//	point_collider_triangle_target_pos_index[thread_No].resize(ave_pair_num[3] + 1);
	//}

	if (edge_edge_target_pos_index[thread_No].size() < ave_pair_num[1]) {
		edge_edge_target_pos_index[thread_No].resize(ave_pair_num[1]+1);
	}
	//if (edge_edge_collider_target_pos_index[thread_No].size() < ave_pair_num[4]) {
	//	edge_edge_collider_target_pos_index[thread_No].resize(ave_pair_num[4]+1);
	//}
	//if (target_position_and_stiffness[thread_No].size() < target_pos_size <<2) {
	//	target_position_and_stiffness[thread_No].resize(target_pos_size <<2);
	//}
}



void Collision::solveXPBDpointTriangleResponse(int thread_No)
{
	if (vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_triangle_pair_index_start_per_thread[thread_No << 1]) {
		solveXPBDpointTriangleResponse(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]][0]);
		for (int i = vertex_triangle_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			solveXPBDpointTriangleResponse(thread_No, i, 0,
				spatial_hashing.vertex_triangle_pair[i][0]);
		}
		solveXPBDpointTriangleResponse(thread_No, vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
	else {
		solveXPBDpointTriangleResponse(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
}


void Collision::re_solveXPBDpointTriangleResponse(int thread_No)
{
	if (point_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1] > point_triangle_target_position_element_start_per_thread[thread_No << 1]) {
		re_solveXPBDpointTriangleResponse(point_triangle_target_position_element_start_per_thread[thread_No << 1], point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_target_pos_index[point_triangle_target_position_element_start_per_thread[thread_No << 1]][0]);
		for (int i = point_triangle_target_position_element_start_per_thread[thread_No << 1] + 1;
			i < point_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
			re_solveXPBDpointTriangleResponse(i, 0,
				point_triangle_target_pos_index[i][0]);
		}
		re_solveXPBDpointTriangleResponse(point_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
			point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 3]);
	}
	else {
		re_solveXPBDpointTriangleResponse(point_triangle_target_position_element_start_per_thread[thread_No << 1], point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 3]);
	}
}




void Collision::re_pointColliderTriangleResponse(int thread_No, TargetPosition* target_pos)
{
	//if (point_collider_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1] > point_collider_triangle_target_position_element_start_per_thread[thread_No << 1]) {
	//	re_pointColliderTriangleResponse(point_collider_triangle_target_position_element_start_per_thread[thread_No << 1], point_collider_triangle_target_position_element_start_per_thread[(thread_No << 1) + 1],
	//		point_collider_triangle_target_pos_index[point_collider_triangle_target_position_element_start_per_thread[thread_No << 1]][0], target_pos);
	//	for (int i = point_collider_triangle_target_position_element_start_per_thread[thread_No << 1] + 1;
	//		i < point_collider_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
	//		re_pointColliderTriangleResponse(i, 0,
	//			point_collider_triangle_target_pos_index[i][0], target_pos);
	//	}
	//	re_pointColliderTriangleResponse(point_collider_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
	//		point_collider_triangle_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	//}
	//else {
	//	re_pointColliderTriangleResponse(point_collider_triangle_target_position_element_start_per_thread[thread_No << 1], point_collider_triangle_target_position_element_start_per_thread[(thread_No << 1) + 1],
	//		point_collider_triangle_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	//}
}


void Collision::re_pointTriangleResponseForIPC(int thread_No, TargetPosition* target_pos)
{
	if (point_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1] > point_triangle_target_position_element_start_per_thread[thread_No << 1]) {
		re_pointTriangleResponseForIPC(point_triangle_target_position_element_start_per_thread[thread_No << 1], point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_target_pos_index[point_triangle_target_position_element_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = point_triangle_target_position_element_start_per_thread[thread_No << 1] + 1;
			i < point_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
			re_pointTriangleResponseForIPC(i, 0,
				point_triangle_target_pos_index[i][0], target_pos);
		}
		re_pointTriangleResponseForIPC(point_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
			point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		re_pointTriangleResponseForIPC(point_triangle_target_position_element_start_per_thread[thread_No << 1], point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
}


void Collision::re_pointTriangleResponse(int thread_No, TargetPosition* target_pos)
{
	if (point_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1] > point_triangle_target_position_element_start_per_thread[thread_No << 1]) {
		re_pointTriangleResponse(point_triangle_target_position_element_start_per_thread[thread_No << 1], point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_target_pos_index[point_triangle_target_position_element_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = point_triangle_target_position_element_start_per_thread[thread_No << 1] + 1;
			i < point_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
			re_pointTriangleResponse(i, 0,
				point_triangle_target_pos_index[i][0], target_pos);
		}
		re_pointTriangleResponse(point_triangle_target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
			point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		re_pointTriangleResponse(point_triangle_target_position_element_start_per_thread[thread_No << 1], point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	/*re_pointTriangleResponse(thread_No, 0,
		point_triangle_target_pos_index[thread_No][0], target_pos);*/
}


void Collision::pointTriangleCollisionPair(int thread_No)
{
	if (vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_triangle_pair_index_start_per_thread[thread_No << 1]) {
		pointTriangleCollisionPair(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]][0]);
		for (int i = vertex_triangle_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			pointTriangleCollisionPair(thread_No, i, 0,
				spatial_hashing.vertex_triangle_pair[i][0]);
		}
		pointTriangleCollisionPair(thread_No, vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
	else {
		pointTriangleCollisionPair(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
}


void Collision::pointTriangleResponseForIPC(int thread_No, TargetPosition* target_pos)
{
	if (vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_triangle_pair_index_start_per_thread[thread_No << 1]) {
		pointTriangleResponseForIPC(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = vertex_triangle_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			pointTriangleResponseForIPC(thread_No, i, 0,
				spatial_hashing.vertex_triangle_pair[i][0], target_pos);
		}
		pointTriangleResponseForIPC(thread_No, vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		pointTriangleResponseForIPC(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
}


void Collision::pointTriangleResponse(int thread_No, TargetPosition* target_pos)
{
	if (vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_triangle_pair_index_start_per_thread[thread_No << 1]) {
		pointTriangleResponse(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = vertex_triangle_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			pointTriangleResponse(thread_No, i, 0,
				spatial_hashing.vertex_triangle_pair[i][0], target_pos);
		}
		pointTriangleResponse(thread_No, vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		pointTriangleResponse(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
}

void Collision::re_pointTriangleColliderResponseForIPC(int thread_No, TargetPosition* target_pos)
{
	if (point_triangle_collider_target_position_element_start_per_thread[(thread_No + 1) << 1] > point_triangle_collider_target_position_element_start_per_thread[thread_No << 1]) {
		re_pointTriangleColliderResponseForIPC(point_triangle_collider_target_position_element_start_per_thread[thread_No << 1], point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_collider_target_pos_index[point_triangle_collider_target_position_element_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = point_triangle_collider_target_position_element_start_per_thread[thread_No << 1] + 1;
			i < point_triangle_collider_target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
			re_pointTriangleColliderResponseForIPC(i, 0,
				point_triangle_collider_target_pos_index[i][0], target_pos);
		}
		re_pointTriangleColliderResponseForIPC(point_triangle_collider_target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
			point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		re_pointTriangleColliderResponseForIPC(point_triangle_collider_target_position_element_start_per_thread[thread_No << 1], point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
}

void Collision::re_pointTriangleColliderResponse(int thread_No, TargetPosition* target_pos)
{
	if (point_triangle_collider_target_position_element_start_per_thread[(thread_No + 1) << 1] > point_triangle_collider_target_position_element_start_per_thread[thread_No << 1]) {
		re_pointTriangleColliderResponse(point_triangle_collider_target_position_element_start_per_thread[thread_No << 1], point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_collider_target_pos_index[point_triangle_collider_target_position_element_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = point_triangle_collider_target_position_element_start_per_thread[thread_No << 1] + 1;
			i < point_triangle_collider_target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
			re_pointTriangleColliderResponse(i, 0,
				point_triangle_collider_target_pos_index[i][0], target_pos);
		}
		re_pointTriangleColliderResponse(point_triangle_collider_target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
			point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		re_pointTriangleColliderResponse(point_triangle_collider_target_position_element_start_per_thread[thread_No << 1], point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
}


void Collision::re_solveXPBDpointTriangleColliderResponse(int thread_No)
{
	if (point_triangle_collider_target_position_element_start_per_thread[(thread_No + 1) << 1] > point_triangle_collider_target_position_element_start_per_thread[thread_No << 1]) {
		re_solveXPBDpointTriangleColliderResponse(point_triangle_collider_target_position_element_start_per_thread[thread_No << 1], point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_collider_target_pos_index[point_triangle_collider_target_position_element_start_per_thread[thread_No << 1]][0]);
		for (int i = point_triangle_collider_target_position_element_start_per_thread[thread_No << 1] + 1;
			i < point_triangle_collider_target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
			re_solveXPBDpointTriangleColliderResponse(i, 0,
				point_triangle_collider_target_pos_index[i][0]);
		}
		re_solveXPBDpointTriangleColliderResponse(point_triangle_collider_target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
			point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 3]);
	}
	else {
		re_solveXPBDpointTriangleColliderResponse(point_triangle_collider_target_position_element_start_per_thread[thread_No << 1], point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 1],
			point_triangle_collider_target_position_element_start_per_thread[(thread_No << 1) + 3]);
	}
}

void Collision::solveXPBDpointTriangleColliderResponse(int thread_No)
{
	if (vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]) {
		solveXPBDpointTriangleColliderResponse(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]][0]);
		for (int i = vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			solveXPBDpointTriangleColliderResponse(thread_No, i, 0,
				spatial_hashing.vertex_obj_triangle_collider_pair[i][0]);
		}
		solveXPBDpointTriangleColliderResponse(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
	else {
		solveXPBDpointTriangleColliderResponse(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
}



void Collision::pointTriangleColliderPair(int thread_No)
{
	if (vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]) {
		pointTriangleColliderPair(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]][0]);
		for (int i = vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			pointTriangleColliderPair(thread_No, i, 0,
				spatial_hashing.vertex_obj_triangle_collider_pair[i][0]);
		}
		pointTriangleColliderPair(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
	else {
		pointTriangleColliderPair(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
}


void Collision::pointTriangleColliderResponseForIPC(int thread_No, TargetPosition* target_pos)
{

	if (vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]) {
		pointTriangleColliderResponseForIPC(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			pointTriangleColliderResponseForIPC(thread_No, i, 0,
				spatial_hashing.vertex_obj_triangle_collider_pair[i][0], target_pos);
		}
		pointTriangleColliderResponseForIPC(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		pointTriangleColliderResponseForIPC(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
}


void Collision::pointTriangleColliderResponse(int thread_No, TargetPosition* target_pos)
{
	if (vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]) {
		pointTriangleColliderResponse(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			pointTriangleColliderResponse(thread_No, i, 0,
				spatial_hashing.vertex_obj_triangle_collider_pair[i][0], target_pos);
		}
		pointTriangleColliderResponse(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		pointTriangleColliderResponse(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
}

void Collision::pointColliderTriangleResponse(int thread_No, TargetPosition* target_pos)
{
	//if (vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1]) {
	//	pointColliderTriangleResponse(thread_No, vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1], vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 1],
	//		spatial_hashing.vertex_collider_triangle_obj_pair[vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1]][0], target_pos);
	//	for (int i = vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1] + 1;
	//		i < vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
	//		pointColliderTriangleResponse(thread_No, i, 0,
	//			spatial_hashing.vertex_collider_triangle_obj_pair[i][0], target_pos);
	//	}
	//	pointColliderTriangleResponse(thread_No, vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
	//		vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	//}
	//else {
	//	pointColliderTriangleResponse(thread_No, vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1], vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 1],
	//		vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	//}
}



void Collision::solveXPBDedgeEdgeResponse(int thread_No)
{
	if (edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1] > edge_edge_pair_index_start_per_thread[thread_No << 1]) {
		solveXPBDedgeEdgeResponse(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.edge_edge_pair[edge_edge_pair_index_start_per_thread[thread_No << 1]][0]);
		for (int i = edge_edge_pair_index_start_per_thread[thread_No << 1] + 1;
			i < edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			solveXPBDedgeEdgeResponse(thread_No, i, 0,
				spatial_hashing.edge_edge_pair[i][0]);
		}
		solveXPBDedgeEdgeResponse(thread_No, edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
	else {
		solveXPBDedgeEdgeResponse(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
}



void Collision::re_solveXPBDedgeEdgeResponse(int thread_No)
{
	if (edge_edge_target_position_element_start_per_thread[(thread_No + 1) << 1] > edge_edge_target_position_element_start_per_thread[thread_No << 1]) {
		re_solveXPBDedgeEdgeResponse(edge_edge_target_position_element_start_per_thread[thread_No << 1], edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 1],
			edge_edge_target_pos_index[edge_edge_target_position_element_start_per_thread[thread_No << 1]][0]);
		for (int i = edge_edge_target_position_element_start_per_thread[thread_No << 1] + 1;
			i < edge_edge_target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
			re_solveXPBDedgeEdgeResponse(i, 0,
				edge_edge_target_pos_index[i][0]);
		}
		re_solveXPBDedgeEdgeResponse(edge_edge_target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
			edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 3]);
	}
	else {
		re_solveXPBDedgeEdgeResponse(edge_edge_target_position_element_start_per_thread[thread_No << 1], edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 1],
			edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 3]);
	}
}



void Collision::re_edgeEdgeColliderResponse(int thread_No, TargetPosition* target_pos)
{
	//if (edge_edge_collider_target_position_element_start_per_thread[(thread_No + 1) << 1] > edge_edge_collider_target_position_element_start_per_thread[thread_No << 1]) {
	//	re_edgeEdgeColliderResponse(edge_edge_collider_target_position_element_start_per_thread[thread_No << 1], edge_edge_collider_target_position_element_start_per_thread[(thread_No << 1) + 1],
	//		edge_edge_collider_target_pos_index[edge_edge_collider_target_position_element_start_per_thread[thread_No << 1]][0], target_pos);
	//	for (int i = edge_edge_collider_target_position_element_start_per_thread[thread_No << 1] + 1;
	//		i < edge_edge_collider_target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
	//		re_edgeEdgeColliderResponse(i, 0,
	//			edge_edge_collider_target_pos_index[i][0], target_pos);
	//	}
	//	re_edgeEdgeColliderResponse(edge_edge_collider_target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
	//		edge_edge_collider_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	//}
	//else {
	//	re_edgeEdgeColliderResponse(edge_edge_collider_target_position_element_start_per_thread[thread_No << 1], edge_edge_collider_target_position_element_start_per_thread[(thread_No << 1) + 1],
	//		edge_edge_collider_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	//}
}


void Collision::re_edgeEdgeResponseForIPC(int thread_No, TargetPosition* target_pos)
{
	if (edge_edge_target_position_element_start_per_thread[(thread_No + 1) << 1] > edge_edge_target_position_element_start_per_thread[thread_No << 1]) {
		re_edgeEdgeResponseForIPC(edge_edge_target_position_element_start_per_thread[thread_No << 1], edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 1],
			edge_edge_target_pos_index[edge_edge_target_position_element_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = edge_edge_target_position_element_start_per_thread[thread_No << 1] + 1;
			i < edge_edge_target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
			re_edgeEdgeResponseForIPC(i, 0,
				edge_edge_target_pos_index[i][0], target_pos);
		}
		re_edgeEdgeResponseForIPC(edge_edge_target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
			edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		re_edgeEdgeResponseForIPC(edge_edge_target_position_element_start_per_thread[thread_No << 1], edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 1],
			edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	/*re_pointTriangleResponse(thread_No, 0,
		point_triangle_target_pos_index[thread_No][0], target_pos);*/
}


void Collision::re_edgeEdgeResponse(int thread_No, TargetPosition* target_pos)
{
	if (edge_edge_target_position_element_start_per_thread[(thread_No + 1) << 1] > edge_edge_target_position_element_start_per_thread[thread_No << 1]) {
		re_edgeEdgeResponse(edge_edge_target_position_element_start_per_thread[thread_No << 1], edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 1],
			edge_edge_target_pos_index[edge_edge_target_position_element_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = edge_edge_target_position_element_start_per_thread[thread_No << 1] + 1;
			i < edge_edge_target_position_element_start_per_thread[(thread_No + 1) << 1]; ++i) {
			re_edgeEdgeResponse(i, 0,
				edge_edge_target_pos_index[i][0], target_pos);
		}
		re_edgeEdgeResponse(edge_edge_target_position_element_start_per_thread[(thread_No + 1) << 1], 0,
			edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		re_edgeEdgeResponse(edge_edge_target_position_element_start_per_thread[thread_No << 1], edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 1],
			edge_edge_target_position_element_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	/*re_pointTriangleResponse(thread_No, 0,
		point_triangle_target_pos_index[thread_No][0], target_pos);*/
}


void Collision::edgeEdgeCollisionPair(int thread_No)
{
	if (edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1] > edge_edge_pair_index_start_per_thread[thread_No << 1]) {
		edgeEdgeCollisionPair(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.edge_edge_pair[edge_edge_pair_index_start_per_thread[thread_No << 1]][0]);
		for (int i = edge_edge_pair_index_start_per_thread[thread_No << 1] + 1;
			i < edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			edgeEdgeCollisionPair(thread_No, i, 0,
				spatial_hashing.edge_edge_pair[i][0]);
		}
		edgeEdgeCollisionPair(thread_No, edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
	else {
		edgeEdgeCollisionPair(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3]);
	}
}

void Collision::edgeEdgeResponseForIPC(int thread_No, TargetPosition* target_pos)
{
	if (edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1] > edge_edge_pair_index_start_per_thread[thread_No << 1]) {
		edgeEdgeResponseForIPC(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.edge_edge_pair[edge_edge_pair_index_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = edge_edge_pair_index_start_per_thread[thread_No << 1] + 1;
			i < edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			edgeEdgeResponseForIPC(thread_No, i, 0,
				spatial_hashing.edge_edge_pair[i][0], target_pos);
		}
		edgeEdgeResponseForIPC(thread_No, edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		edgeEdgeResponseForIPC(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
}


void Collision::edgeEdgeResponse(int thread_No, TargetPosition* target_pos)
{
	if (edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1] > edge_edge_pair_index_start_per_thread[thread_No << 1]) {
		edgeEdgeResponse(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.edge_edge_pair[edge_edge_pair_index_start_per_thread[thread_No << 1]][0], target_pos);
		for (int i = edge_edge_pair_index_start_per_thread[thread_No << 1] + 1;
			i < edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			edgeEdgeResponse(thread_No, i, 0,
				spatial_hashing.edge_edge_pair[i][0], target_pos);
		}
		edgeEdgeResponse(thread_No, edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
	else {
		edgeEdgeResponse(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	}
}

void Collision::edgeEdgeColliderResponse(int thread_No, TargetPosition* target_pos)
{
	//if (edge_edge_pair_collider_index_start_per_thread[(thread_No + 1) << 1] > edge_edge_pair_collider_index_start_per_thread[thread_No << 1]) {
	//	edgeEdgeColliderResponse(thread_No, edge_edge_pair_collider_index_start_per_thread[thread_No << 1], edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 1],
	//		spatial_hashing.edge_edge_pair_collider[edge_edge_pair_collider_index_start_per_thread[thread_No << 1]][0], target_pos);
	//	for (int i = edge_edge_pair_collider_index_start_per_thread[thread_No << 1] + 1;
	//		i < edge_edge_pair_collider_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
	//		edgeEdgeColliderResponse(thread_No, i, 0,
	//			spatial_hashing.edge_edge_pair_collider[i][0], target_pos);
	//	}
	//	edgeEdgeColliderResponse(thread_No, edge_edge_pair_collider_index_start_per_thread[(thread_No + 1) << 1], 0,
	//		edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	//}
	//else {
	//	edgeEdgeColliderResponse(thread_No, edge_edge_pair_collider_index_start_per_thread[thread_No << 1], edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 1],
	//		edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 3], target_pos);
	//}
}

void Collision::computeEnergy(unsigned int pair_thread_No, int index_start, int index_end, double& energy)
{
	//unsigned int* target_pos_index_record_ = target_position_index[pair_thread_No].data() + 1;
	//double* target_position_ = target_position_and_stiffness[pair_thread_No].data();
	//unsigned int position_index_ = index_start <<2;
	//index_end <<= 1;
	//double* current_vertex_position;
	//for (int i = index_start<<1; i < index_end; i += 2) {
	//	current_vertex_position = vertex_position[*(target_pos_index_record_ + i + 1)][*(target_pos_index_record_ + i)].data();
	//	energy += *(target_position_ + position_index_+3)* EDGE_LENGTH(current_vertex_position, (target_position_ + position_index_));
	//	position_index_ += 4;
	//	////std::cout << *(target_position_ + position_index_ + 3) << std::endl;
	//}
}

void Collision::re_solveXPBDpointTriangleResponse(unsigned int pair_thread_No, int index_start, int index_end)
{
	unsigned int* target_pos_index_record_ = point_triangle_target_pos_index[pair_thread_No].data() + 1;
	double tolerance = tolerance_radius[SELF_POINT_TRIANGLE];
	double epsilon_ = epsilon;
	int* indices;

	unsigned int* pair;
	//	double* lambda_ = lambda->data() + collision_lambda_index_start[1];triangle_normal_not_normalized
	double energy_ = 0.0;
	double friction_coe_self = friction_coe[0];
	for (int i = index_start; i < index_end; i += 4) {
		pair = target_pos_index_record_ + i;
		indices = triangle_indices[*(pair + 3)][*(pair + 2)].data();
		dcd.XPBDpointSelfTriangle(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
			vertex_for_render[*(pair + 3)][indices[0]].data(), vertex_for_render[*(pair + 3)][indices[1]].data(),
			vertex_for_render[*(pair + 3)][indices[2]].data(),
			vertex_position[*(pair + 3)][indices[0]].data(), vertex_position[*(pair + 3)][indices[1]].data(),
			vertex_position[*(pair + 3)][indices[2]].data(), triangle_normal_render[*(pair + 3)][*(pair + 2)].data(),
			triangle_normal[*(pair + 3)][*(pair + 2)].data(),
			tolerance,
			mass_inv[*(pair + 1)][*pair], mass_inv[*(pair + 3)][indices[0]], mass_inv[*(pair + 3)][indices[1]], mass_inv[*(pair + 3)][indices[2]],
			//1.0, 1.0, 1.0, 1.0,
			triangle_normal_magnitude_reciprocal[*(pair + 3)][*(pair + 2)], friction_coe_self);
	}
}
void Collision::solveXPBDpointTriangleResponse(unsigned int thread_No, unsigned int pair_thread_No, int index_start, int index_end)
{
	unsigned int* pair_ = spatial_hashing.vertex_triangle_pair[pair_thread_No] + 1;
	double tolerance = tolerance_radius[SELF_POINT_TRIANGLE];
	unsigned int* target_pos_index_ = point_triangle_target_pos_index[thread_No].data() + 1 + point_triangle_target_pos_index[thread_No][0];
	int* indices;
	unsigned int* pair;
//	double* lambda_ = lambda->data() + collision_lambda_index_start[1];
	double friction_coe_self = friction_coe[0];
	for (int i = index_start; i < index_end; i += 4) {
		pair = pair_ + i;
		indices = triangle_indices[*(pair + 3)][*(pair + 2)].data();
		if (dcd.XPBDpointSelfTriangle(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
			vertex_for_render[*(pair + 3)][indices[0]].data(), vertex_for_render[*(pair + 3)][indices[1]].data(),
			vertex_for_render[*(pair + 3)][indices[2]].data(),
			vertex_position[*(pair + 3)][indices[0]].data(), vertex_position[*(pair + 3)][indices[1]].data(),
			vertex_position[*(pair + 3)][indices[2]].data(), triangle_normal_render[*(pair + 3)][*(pair + 2)].data(),
			triangle_normal[*(pair + 3)][*(pair + 2)].data(),
			tolerance,
			mass_inv[*(pair + 1)][*pair], mass_inv[*(pair + 3)][indices[0]], mass_inv[*(pair + 3)][indices[1]], mass_inv[*(pair + 3)][indices[2]],
			//1.0, 1.0, 1.0, 1.0,
			triangle_normal_magnitude_reciprocal[*(pair + 3)][*(pair + 2)], friction_coe_self)) {
			memcpy(target_pos_index_, pair, 16);
			target_pos_index_ += 4;

		};
	}
	point_triangle_target_pos_index[thread_No][0] = target_pos_index_ - point_triangle_target_pos_index[thread_No].data() - 1;
}



void Collision::re_pointColliderTriangleResponse(unsigned int pair_thread_No, int index_start, int index_end, TargetPosition* target_pos)
{
	//unsigned int* target_pos_index_record_ = point_collider_triangle_target_pos_index[pair_thread_No].data() + 1;
	//double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	//int* indices;
	//double target_pos_tri[3][3];
	//double stiffness_initial = collision_stiffness[BODY_POINT_TRIANGLE];
	//double stiffness;
	//std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	//unsigned int* pair;
	//stiffness = stiffness_initial;
	//for (int i = index_start; i < index_end; i += 4) {
	//	pair = target_pos_index_record_ + i;
	//	indices = triangle_indices[*(pair + 3)][*(pair + 2)].data();
	//	if (dcd.pointColliderTriangle(vertex_for_render_collider[*(pair + 1)][*pair].data(), vertex_position_collider[*(pair + 1)][*pair].data(),
	//		vertex_for_render[*(pair + 3)][indices[0]].data(), vertex_for_render[*(pair + 3)][indices[1]].data(),
	//		vertex_for_render[*(pair + 3)][indices[2]].data(),
	//		vertex_position[*(pair + 3)][indices[0]].data(), vertex_position[*(pair + 3)][indices[1]].data(),
	//		vertex_position[*(pair + 3)][indices[2]].data(), triangle_normal_render[*(pair + 3)][*(pair + 2)].data(),
	//		triangle_normal[*(pair + 3)][*(pair + 2)].data(),
	//		target_pos_tri[0], target_pos_tri[1], target_pos_tri[2], tolerance,
	//		mass[*(pair + 3)][indices[0]], mass[*(pair + 3)][indices[1]], mass[*(pair + 3)][indices[2]],
	//		//1.0,1.0,1.0,1.0,
	//		triangle_normal_magnitude_reciprocal[*(pair + 3)][*(pair + 2)]))
	//	{
	//		for (int j = 0; j < 3; ++j) {
	//			addTargetPosToSystem(vertex_b_sum[*(pair + 3)][indices[j]].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][indices[j]].data(),
	//				target_pos_tri[j], stiffness);
	//		}
	//	}
	//	else {
	//		for (int j = 0; j < 3; ++j) {
	//			construct_b_sum(vertex_b_sum[*(pair + 3)][indices[j]].data(), vertex_position[*(pair + 3)][indices[j]].data(), stiffness);
	//		}
	//	}
	//}
}


void Collision::re_pointTriangleResponseForIPC(unsigned int pair_thread_No, int index_start, int index_end, TargetPosition* target_pos)
{
	unsigned int* target_pos_index_record_ = point_triangle_target_pos_index[pair_thread_No].data() + 1;
	double* target_pos_record_ = point_triangle_target_pos_record[pair_thread_No].data();
	int* indices;
	double* target_pos_v;
	double stiffness_initial = collision_stiffness[SELF_POINT_TRIANGLE];
	double stiffness;
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	unsigned int* pair;

	for (int i = index_start; i < index_end; i += 4) {
		pair = target_pos_index_record_ + i;
		target_pos_v = target_pos_record_ + (i>>2)*13;
		stiffness = target_pos_v[12];
		indices = triangle_indices[*(pair + 3)][*(pair + 2)].data();
		addTargetPosToSystem(vertex_b_sum[*(pair + 1)][*pair].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*pair].data(),
			target_pos_v, stiffness);
		for (int j = 0; j < 3; ++j) {
			addTargetPosToSystem(vertex_b_sum[*(pair + 3)][indices[j]].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][indices[j]].data(),
				target_pos_v+3*(j+1), stiffness);
		}

	}
}


void Collision::re_pointTriangleResponse(unsigned int pair_thread_No, int index_start, int index_end, TargetPosition* target_pos)
{
	unsigned int* target_pos_index_record_ = point_triangle_target_pos_index[pair_thread_No].data() + 1;
	double tolerance = tolerance_radius[SELF_POINT_TRIANGLE];
	int* indices;
	double target_pos_v[3];
	double target_pos_tri[3][3];
	double stiffness_initial = collision_stiffness[SELF_POINT_TRIANGLE];
	double stiffness;
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	unsigned int* pair;

	for (int i = index_start; i < index_end; i += 4) {
		stiffness = stiffness_initial;
		pair = target_pos_index_record_ + i;
		indices = triangle_indices[*(pair + 3)][*(pair + 2)].data();
		//if (collision_constraint.pointTriangleResponse(vertex_for_render[*(pair+ 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
		//	vertex_for_render[*(pair + 3)][indices[0]].data(), vertex_for_render[pair[i + 3]][indices[1]].data(),
		//	vertex_for_render[*(pair + 3)][indices[2]].data(),
		//	vertex_position[*(pair + 3)][indices[0]].data(), vertex_position[*(pair + 3)][indices[1]].data(),
		//	vertex_position[*(pair + 3)][indices[2]].data(), triangle_normal_render[*(pair + 3)][*(pair + 2)].data(),
		//	target_pos_v, target_pos_tri[0], target_pos_tri[1], target_pos_tri[2], d_hat_, stiffness, epsilon_,
		//	mass[*(pair + 1)][*pair],
		//	mass[*(pair + 3)][indices[0]], mass[*(pair + 3)][indices[1]], mass[*(pair + 3)][indices[2]]))
		if (dcd.pointSelfTriangle(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
			vertex_for_render[*(pair + 3)][indices[0]].data(), vertex_for_render[*(pair + 3)][indices[1]].data(),
			vertex_for_render[*(pair + 3)][indices[2]].data(),
			vertex_position[*(pair + 3)][indices[0]].data(), vertex_position[*(pair + 3)][indices[1]].data(),
			vertex_position[*(pair + 3)][indices[2]].data(), triangle_normal_render[*(pair + 3)][*(pair + 2)].data(),
			triangle_normal[*(pair + 3)][*(pair + 2)].data(),
			target_pos_v, target_pos_tri[0], target_pos_tri[1], target_pos_tri[2], tolerance,
			mass[*(pair + 1)][*pair],mass[*(pair + 3)][indices[0]], mass[*(pair + 3)][indices[1]], mass[*(pair + 3)][indices[2]],
			//1.0,1.0,1.0,1.0,
			triangle_normal_magnitude_reciprocal[*(pair + 3)][*(pair + 2)]))
		{
			addTargetPosToSystem(vertex_b_sum[*(pair + 1)][*pair].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*pair].data(),
				target_pos_v, stiffness);
			for (int j = 0; j < 3; ++j) {
				addTargetPosToSystem(vertex_b_sum[*(pair + 3)][indices[j]].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][indices[j]].data(),
					target_pos_tri[j], stiffness);
			}
		}
		else {
			construct_b_sum(vertex_b_sum[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(), stiffness);
			for (int j = 0; j < 3; ++j) {
				construct_b_sum(vertex_b_sum[*(pair + 3)][indices[j]].data(), vertex_position[*(pair + 3)][indices[j]].data(),stiffness);
			}
		}
	}
}

void Collision::re_solveXPBDedgeEdgeResponse(unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index)
{
	unsigned int* target_pos_index_record_ = edge_edge_target_pos_index[pair_thread_No].data() + 1;
	double tolerance = tolerance_radius[SELF_EDGE_EDGE];
	unsigned int* indices;
	unsigned int* compare_indices;

	unsigned int* pair;

	double friction_coe_self = friction_coe[0];

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		pair = target_pos_index_record_ + i;
		indices = edge_vertices[*(pair + 1)] + ((*(pair)) << 1);
		compare_indices = edge_vertices[*(pair + 3)] + ((*(pair + 2)) << 1);
		dcd.XPBDedgeEdge(vertex_position[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
			vertex_for_render[*(pair + 1)][*indices].data(), vertex_for_render[*(pair + 1)][*(indices + 1)].data(),
			vertex_position[*(pair + 3)][*compare_indices].data(), vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
			vertex_for_render[*(pair + 3)][*compare_indices].data(), vertex_for_render[*(pair + 3)][*(compare_indices + 1)].data(),
			tolerance,
			mass_inv[*(pair + 1)][*indices], mass_inv[*(pair + 1)][*(indices + 1)], mass_inv[*(pair + 3)][*compare_indices], mass_inv[*(pair + 3)][*(compare_indices + 1)], friction_coe_self);
	}
}


void Collision::solveXPBDedgeEdgeResponse(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index)
{
	unsigned int* pair_ = spatial_hashing.edge_edge_pair[pair_thread_No] + 1;
	double tolerance = tolerance_radius[SELF_EDGE_EDGE];
	unsigned int* indices;
	unsigned int* compare_indices;
	unsigned int* target_pos_index_ = edge_edge_target_pos_index[thread_No].data() + 1 + edge_edge_target_pos_index[thread_No][0];

	unsigned int* pair;
	//double* lambda_ = lambda->data() + collision_lambda_index_start[2];
	double friction_coe_self = friction_coe[0];
	for (int i = start_pair_index; i < end_pair_index; i += 4) {

		pair = pair_ + i;
		indices = edge_vertices[*(pair + 1)] + ((*(pair)) << 1);
		compare_indices = edge_vertices[*(pair + 3)] + ((*(pair + 2)) << 1);
		if (dcd.XPBDedgeEdge(vertex_position[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
			vertex_for_render[*(pair + 1)][*indices].data(), vertex_for_render[*(pair + 1)][*(indices + 1)].data(),
			vertex_position[*(pair + 3)][*compare_indices].data(), vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
			vertex_for_render[*(pair + 3)][*compare_indices].data(), vertex_for_render[*(pair + 3)][*(compare_indices + 1)].data(),
			tolerance,
			mass_inv[*(pair + 1)][*indices], mass_inv[*(pair + 1)][*(indices + 1)], mass_inv[*(pair + 3)][*compare_indices], mass_inv[*(pair + 3)][*(compare_indices + 1)], friction_coe_self))
		{
			memcpy(target_pos_index_, pair, 16);
			target_pos_index_ += 4;
		}
	}

	edge_edge_target_pos_index[thread_No][0] = target_pos_index_ - edge_edge_target_pos_index[thread_No].data() - 1;
}



void Collision::re_edgeEdgeColliderResponse(unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, TargetPosition* target_pos)
{
	//unsigned int* target_pos_index_record_ = edge_edge_collider_target_pos_index[pair_thread_No].data() + 1;
	//double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	//unsigned int* indices;
	//unsigned int* compare_indices;
	//double target_pos_edge_0[3];
	//double target_pos_edge_1[3];
	//double target_pos_compare_edge_0[3];
	//double target_pos_compare_edge_1[3];
	//double stiffness;
	//double stiffness_initial = collision_stiffness[BODY_POINT_TRIANGLE];
	//std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	//unsigned int* pair;
	//stiffness = stiffness_initial;
	//for (int i = start_pair_index; i < end_pair_index; i += 4) {	
	//	pair = target_pos_index_record_ + i;
	//	indices = edge_vertices[*(pair + 1)] + ((*pair) << 1);
	//	compare_indices = collider_edge_vertices[*(pair + 3)] + (2*(*(pair + 2)));
	//	if (dcd.edgeEdgeCollider(target_pos_edge_0, target_pos_edge_1,
	//		vertex_position[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
	//		vertex_for_render[*(pair + 1)][*indices].data(), vertex_for_render[*(pair + 1)][*(indices + 1)].data(),
	//		vertex_position_collider[*(pair + 3)][*compare_indices].data(), vertex_position_collider[*(pair + 3)][*(compare_indices + 1)].data(),
	//		vertex_for_render_collider[*(pair + 3)][*compare_indices].data(), vertex_for_render_collider[*(pair + 3)][*(compare_indices + 1)].data(),
	//		tolerance, mass[*(pair + 1)][*indices], mass[*(pair + 1)][*(indices + 1)]))
	//	{
	//		addTargetPosToSystem(vertex_b_sum[*(pair + 1)][*indices].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*indices].data(),
	//			target_pos_edge_0, stiffness);
	//		addTargetPosToSystem(vertex_b_sum[*(pair + 1)][*(indices + 1)].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*(indices + 1)].data(),
	//			target_pos_edge_1, stiffness);
	//	}
	//	else {
	//		construct_b_sum(vertex_b_sum[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*indices].data(), stiffness);
	//		construct_b_sum(vertex_b_sum[*(pair + 1)][*(indices + 1)].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
	//			stiffness);
	//	}
	//}
}


void Collision::re_edgeEdgeResponseForIPC(unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, TargetPosition* target_pos)
{
	unsigned int* target_pos_index_record_ = edge_edge_target_pos_index[pair_thread_No].data() + 1;
	double* target_pos_record_ = edge_edge_target_pos_record[pair_thread_No].data();
	double tolerance = tolerance_radius[SELF_EDGE_EDGE];
	unsigned int* indices;
	unsigned int* compare_indices;
	double* target_pos_edge;
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	double stiffness_initial = collision_stiffness[SELF_EDGE_EDGE];
	double stiffness;
	unsigned int* pair;

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		pair = target_pos_index_record_ + i;
		indices = edge_vertices[*(pair + 1)] + ((*(pair)) << 1);
		compare_indices = edge_vertices[*(pair + 3)] + ((*(pair + 2)) << 1);
		target_pos_edge = target_pos_record_ + (i >> 2) * 13;
		stiffness = target_pos_edge[12];
		addTargetPosToSystem(vertex_b_sum[*(pair + 1)][*indices].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*indices].data(),
			target_pos_edge, stiffness);
		addTargetPosToSystem(vertex_b_sum[*(pair + 1)][*(indices + 1)].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*(indices + 1)].data(),
			target_pos_edge+3, stiffness);
		addTargetPosToSystem(vertex_b_sum[*(pair + 3)][*compare_indices].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][*compare_indices].data(),
			target_pos_edge+6, stiffness);
		addTargetPosToSystem(vertex_b_sum[*(pair + 3)][*(compare_indices + 1)].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
			target_pos_edge+9, stiffness);
	
	}
}



void Collision::re_edgeEdgeResponse(unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, TargetPosition* target_pos)
{
	//unsigned int* target_pos_index_ = target_position_index[thread_No].data() + 1 + (target_position_index[thread_No][0] << 1);
	unsigned int* target_pos_index_record_ = edge_edge_target_pos_index[pair_thread_No].data() + 1;
	double tolerance = tolerance_radius[SELF_EDGE_EDGE];
	unsigned int* indices;
	unsigned int* compare_indices;
	double target_pos_edge_0[3];
	double target_pos_edge_1[3];
	double target_pos_compare_edge_0[3];
	double target_pos_compare_edge_1[3];
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	double stiffness_initial = collision_stiffness[SELF_EDGE_EDGE];
	double stiffness;
	unsigned int* pair;
	stiffness = stiffness_initial;
	for (int i = start_pair_index; i < end_pair_index; i += 4) {	
		pair = target_pos_index_record_ + i;
		indices = edge_vertices[*(pair + 1)] + ((*(pair)) << 1);
		compare_indices = edge_vertices[*(pair + 3)] + ((*(pair + 2)) << 1);
		if (dcd.edgeEdge(target_pos_edge_0, target_pos_edge_1, target_pos_compare_edge_0, target_pos_compare_edge_1,
			vertex_position[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
			vertex_for_render[*(pair + 1)][*indices].data(), vertex_for_render[*(pair + 1)][*(indices + 1)].data(),
			vertex_position[*(pair + 3)][*compare_indices].data(), vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
			vertex_for_render[*(pair + 3)][*compare_indices].data(), vertex_for_render[*(pair + 3)][*(compare_indices + 1)].data(),
			tolerance, 
			//1.0, 1.0, 1.0, 1.0
			mass[*(pair + 1)][*indices], mass[*(pair + 1)][*(indices + 1)],	mass[*(pair + 3)][*compare_indices], mass[*(pair + 3)][*(compare_indices + 1)]
		))
		{
			addTargetPosToSystem(vertex_b_sum[*(pair + 1)][*indices].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*indices].data(),
				target_pos_edge_0, stiffness);
			addTargetPosToSystem(vertex_b_sum[*(pair + 1)][*(indices + 1)].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*(indices + 1)].data(),
				target_pos_edge_1, stiffness);
			addTargetPosToSystem(vertex_b_sum[*(pair + 3)][*compare_indices].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][*compare_indices].data(),
				target_pos_compare_edge_0, stiffness);
			addTargetPosToSystem(vertex_b_sum[*(pair + 3)][*(compare_indices + 1)].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
				target_pos_compare_edge_1, stiffness);
		}
		else{
			construct_b_sum(vertex_b_sum[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*indices].data(),stiffness);
			construct_b_sum(vertex_b_sum[*(pair + 1)][*(indices + 1)].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
				stiffness);
			construct_b_sum(vertex_b_sum[*(pair + 3)][*compare_indices].data(), vertex_position[*(pair + 3)][*compare_indices].data(),
				stiffness);
			construct_b_sum(vertex_b_sum[*(pair + 3)][*(compare_indices + 1)].data(), vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
				stiffness);
		}
	}
}





void Collision::pointTriangleCollisionPair(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index)
{
	//////std::cout << thread_No<<" "<< start_pair_index << " " << end_pair_index << std::endl;
	unsigned int* target_pos_index_ = point_triangle_target_pos_index[thread_No].data() + 1 + point_triangle_target_pos_index[thread_No][0];
	//unsigned int* target_pos_index_ = target_position_index[thread_No].data() + 1 + (target_position_index[thread_No][0] << 1);
	//double* target_pos_ = target_position_and_stiffness[thread_No].data() + (target_position_index[thread_No][0]<<2);

	unsigned int* pair_ = spatial_hashing.vertex_triangle_pair[pair_thread_No] + 1;
	//double d_hat_ = d_hat;
	double tolerance = tolerance_radius[SELF_POINT_TRIANGLE];
	int* indices;
	unsigned int* pair;
	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		pair = pair_ + i;
		indices = triangle_indices[*(pair + 3)][*(pair + 2)].data();
		if (dcd.checkPointTriangle(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
			vertex_for_render[*(pair + 3)][indices[0]].data(), vertex_for_render[*(pair + 3)][indices[1]].data(),
			vertex_for_render[*(pair + 3)][indices[2]].data(),
			vertex_position[*(pair + 3)][indices[0]].data(), vertex_position[*(pair + 3)][indices[1]].data(),
			vertex_position[*(pair + 3)][indices[2]].data(), triangle_normal_render[*(pair + 3)][*(pair + 2)].data(),
			tolerance))
		{
			memcpy(target_pos_index_, pair, 16);
			target_pos_index_ += 4;
		}
	}

	point_triangle_target_pos_index[thread_No][0] = target_pos_index_ - point_triangle_target_pos_index[thread_No].data() - 1;
}




void Collision::floorCollisionVertex(int thread_No)
{
	unsigned int start, end;
	std::array<double, 3>* vertex_pos;

	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	std::vector<double>* record_stiffness = target_pos->stiffness.data();
	bool** vetex_need_update = target_pos->need_update;

	double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	int* indices;
	double target_pos_v[3];
	double stiffness_initial = collision_stiffness[BODY_POINT_TRIANGLE];
	unsigned int dimension = floor->dimension;
	bool normal_direction = floor->normal_direction;
	double floor_value = floor->value;
	std::vector<unsigned int>* floor_collision_vertex_= &floor_collision_vertex[thread_No];
	floor_collision_vertex_->clear();

	for (int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		start = vertex_index_start_per_thread[i][thread_No];
		end = vertex_index_start_per_thread[i][thread_No + 1];
		for (int j = start; j < end; ++j) {
			if (dcd.PDFloor(target_pos_v, vertex_pos[j].data(), dimension, normal_direction, tolerance, floor_value)) {
				floor_collision_vertex_->emplace_back(j);
				floor_collision_vertex_->emplace_back(i);
				addTargetPosToSystemTotal(vertex_b_sum[i][j].data(), target_pos->collision_energy, vertex_position[i][j].data(),
					target_pos_v, stiffness_initial, record_stiffness[i][j], vetex_need_update[i][j]);
			}
		}
	}
}


void Collision::re_FloorCollisionVertexForIPC(int thread_No)
{
	unsigned int start, end;

	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	std::vector<double>* record_stiffness = target_pos->stiffness.data();
	int* indices;
	double* target_pos_v;

	unsigned int total_size = floor_collision_vertex[thread_No].size();
	unsigned int* floor_collide_vertex = floor_collision_vertex[thread_No].data();

	double* collision_record_ = floor_collision_record[thread_No].data();

	for (int j = 0; j < total_size; j += 2) {
		target_pos_v = collision_record_ + (j << 1);
		
		////std::cout << "find record " << target_pos_v[3] << " " << floor_collide_vertex[j] << " " << target_pos_v[1]<<" "<< 
		//	vertex_position[floor_collide_vertex[j + 1]][floor_collide_vertex[j]][1] << std::endl;

		addTargetPosToSystem(vertex_b_sum[floor_collide_vertex[j + 1]][floor_collide_vertex[j]].data(), target_pos->collision_energy,
			vertex_position[floor_collide_vertex[j + 1]][floor_collide_vertex[j]].data(), target_pos_v, target_pos_v[3]);	

	}

}


void Collision::re_FloorCollisionVertex(int thread_No)
{
	unsigned int start, end;

	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	std::vector<double>* record_stiffness = target_pos->stiffness.data();
	bool** vetex_need_update = target_pos->need_update;

	double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	int* indices;
	double target_pos_v[3];
	double stiffness_initial = collision_stiffness[BODY_POINT_TRIANGLE];
	unsigned int dimension = floor->dimension;
	bool normal_direction = floor->normal_direction;
	double floor_value = floor->value;
	std::vector<unsigned int>* floor_collision_vertex_ = &floor_collision_vertex[thread_No];

	unsigned int total_size = floor_collision_vertex[thread_No].size();
	unsigned int* floor_collide_vertex= floor_collision_vertex[thread_No].data();

	for (int j = 0; j < total_size; j += 2) {
		if (dcd.PDFloor(target_pos_v, vertex_position[floor_collide_vertex[j+1]][floor_collide_vertex[j]].data(), 
			dimension, normal_direction, tolerance, floor_value)) {
			addTargetPosToSystem(vertex_b_sum[floor_collide_vertex[j + 1]][floor_collide_vertex[j]].data(), target_pos->collision_energy, 
				vertex_position[floor_collide_vertex[j + 1]][floor_collide_vertex[j]].data(),	target_pos_v, stiffness_initial);
		}
		else {
			construct_b_sum(vertex_b_sum[floor_collide_vertex[j + 1]][floor_collide_vertex[j]].data(),
				vertex_position[floor_collide_vertex[j + 1]][floor_collide_vertex[j]].data(), stiffness_initial);
		}
	}

}

void Collision::pointColliderTriangleResponse(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, TargetPosition* target_pos)
{
	//unsigned int* target_pos_index_ = point_collider_triangle_target_pos_index[thread_No].data() + 1 + point_collider_triangle_target_pos_index[thread_No][0];
	//unsigned int* pair_ = spatial_hashing.vertex_collider_triangle_obj_pair[pair_thread_No] + 1;
	////double d_hat_ = d_hat;
	//double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	//int* indices;
	//double target_pos_tri[3][3];
	//double stiffness_initial = collision_stiffness[BODY_POINT_TRIANGLE];
	//double stiffness;
	//std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	//std::vector<double>* record_stiffness = target_pos->stiffness.data();
	//bool** vetex_need_update = target_pos->need_update;

	//unsigned int* pair;
	//stiffness = stiffness_initial;

	//for (int i = start_pair_index; i < end_pair_index; i += 4) {
	//	pair = pair_ + i;
	//	indices = triangle_indices[*(pair + 3)][*(pair + 2)].data();
	//	if (dcd.pointColliderTriangle(vertex_for_render_collider[*(pair + 1)][*pair].data(), vertex_position_collider[*(pair + 1)][*pair].data(),
	//		vertex_for_render[*(pair + 3)][indices[0]].data(), vertex_for_render[*(pair + 3)][indices[1]].data(),
	//		vertex_for_render[*(pair + 3)][indices[2]].data(),
	//		vertex_position[*(pair + 3)][indices[0]].data(), vertex_position[*(pair + 3)][indices[1]].data(),
	//		vertex_position[*(pair + 3)][indices[2]].data(), triangle_normal_render[*(pair + 3)][*(pair + 2)].data(), triangle_normal[*(pair + 3)][*(pair + 2)].data(),
	//		target_pos_tri[0], target_pos_tri[1], target_pos_tri[2],tolerance,
	//		mass[*(pair + 3)][indices[0]], mass[*(pair + 3)][indices[1]], mass[*(pair + 3)][indices[2]], triangle_normal_magnitude_reciprocal[*(pair + 3)][*(pair + 2)]
	//	))
	//	{
	//		memcpy(target_pos_index_, pair, 16);
	//		target_pos_index_ += 4;
	//		for (int j = 0; j < 3; ++j) {
	//			addTargetPosToSystemTotal(vertex_b_sum[*(pair + 3)][indices[j]].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][indices[j]].data(),
	//				target_pos_tri[j], stiffness, record_stiffness[*(pair + 3)][indices[j]], vetex_need_update[*(pair + 3)][indices[j]]);
	//		}
	//	}
	//}
	//point_collider_triangle_target_pos_index[thread_No][0] = target_pos_index_ - point_collider_triangle_target_pos_index[thread_No].data() - 1;
}
void Collision::floorCollisionVertexForIPC(int thread_No)
{
	unsigned int start, end;
	std::array<double, 3>* vertex_pos;

	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	std::vector<double>* record_stiffness = target_pos->stiffness.data();
	bool** vetex_need_update = target_pos->need_update;

	double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	int* indices;
	double target_pos_v[3];
	double stiffness_initial = collision_stiffness[BODY_POINT_TRIANGLE];
	unsigned int dimension = floor->dimension;
	bool normal_direction = floor->normal_direction;
	double floor_value = floor->value;
	std::vector<unsigned int>* floor_collision_vertex_ = &floor_collision_vertex[thread_No];
	floor_collision_vertex_->clear();

	std::vector<double>* floor_collision_record_ = &floor_collision_record[thread_No];
	floor_collision_record_->clear();
	
	std::array<double, 3>* initial_vertex_pos;
	double stiffness;
	unsigned int* vertex_index_on_surface_;

	std::array<double, 3>* free_vertex_pos;

	for (int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		free_vertex_pos = vertex_collision_free[i];
		initial_vertex_pos = vertex_for_render[i];
		start = vertex_index_start_per_thread[i][thread_No];
		end = vertex_index_start_per_thread[i][thread_No + 1];

		if (i < cloth->size()) {
			for (int j = start; j < end; ++j) {
				stiffness = stiffness_initial;
				if (collision_constraint.floorResponse(target_pos_v, vertex_pos[j].data(), initial_vertex_pos[j].data(),
					dimension, normal_direction, floor_value, d_hat, stiffness, free_vertex_pos[j].data())) {

					addTargetPosToSystemTotal(vertex_b_sum[i][j].data(), target_pos->collision_energy, vertex_position[i][j].data(),
						target_pos_v, stiffness, record_stiffness[i][j], vetex_need_update[i][j]);
					floor_collision_vertex_->emplace_back(j);
					floor_collision_vertex_->emplace_back(i);

					floor_collision_record_->emplace_back(target_pos_v[0]);
					floor_collision_record_->emplace_back(target_pos_v[1]);
					floor_collision_record_->emplace_back(target_pos_v[2]);
					floor_collision_record_->emplace_back(stiffness);

					

				}
			}
		}
		else {
			vertex_index_on_surface_ = vertex_index_on_surface[i];
			int j;
			for (int k = start; k < end; ++k) {
				j = vertex_index_on_surface_[k];
				stiffness = stiffness_initial;
				if (collision_constraint.floorResponse(target_pos_v, vertex_pos[j].data(), initial_vertex_pos[j].data(),
					dimension, normal_direction, floor_value, d_hat, stiffness, free_vertex_pos[j].data())) {

					addTargetPosToSystemTotal(vertex_b_sum[i][j].data(), target_pos->collision_energy, vertex_position[i][j].data(),
						target_pos_v, stiffness, record_stiffness[i][j], vetex_need_update[i][j]);
					floor_collision_vertex_->emplace_back(j);
					floor_collision_vertex_->emplace_back(i);

					floor_collision_record_->emplace_back(target_pos_v[0]);
					floor_collision_record_->emplace_back(target_pos_v[1]);
					floor_collision_record_->emplace_back(target_pos_v[2]);
					floor_collision_record_->emplace_back(stiffness);

					////std::cout << "record stiffness " << stiffness << " "<<j<<" " << target_pos_v[1] << std::endl;

				}
			}
		}
	}
}


void Collision::pointTriangleResponseForIPC(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, TargetPosition* target_pos)
{
	//unsigned int* target_pos_index_ = point_triangle_target_pos_index[thread_No].data() + 1 + point_triangle_target_pos_index[thread_No][0];
	//double* target_pos_record = point_triangle_target_pos_record[thread_No].data()+(point_triangle_target_pos_index[thread_No][0]>>2)*13;

	//double* collision_time = record_VT_collision_time[pair_thread_No].data();

	//unsigned int* pair_ = spatial_hashing.vertex_triangle_pair[pair_thread_No] + 1;
	//double d_hat_ = d_hat;
	//double tolerance = tolerance_radius[SELF_POINT_TRIANGLE];
	//double epsilon_ = epsilon;
	//int* indices;
	//double target_pos_v[12];
	//double stiffness_initial = collision_stiffness[SELF_POINT_TRIANGLE];
	//double stiffness;
	//std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	//std::vector<double>* record_stiffness = target_pos->stiffness.data();
	//bool** vetex_need_update = target_pos->need_update;

	//unsigned int* pair;



	////indices = triangle_indices[1][3].data();
	////////std::cout << indices[0] << " " << indices[1] << " " << indices[2] << std::endl;
	////////std::cout << "////// "<< start_pair_index <<" "<< end_pair_index << std::endl;
	////////std::cout << tetrahedron->data()[0].mesh_struct.vertex_for_render[2][0] << std::endl;
	////////std::cout << vertex_for_render[0][2][0] << " " << vertex_for_render[0][2][1] << " " << vertex_for_render[0][2][2] << std::endl;
	////////std::cout << vertex_for_render[1][indices[0]][0] << " " << vertex_for_render[1][indices[0]][1] << " " << vertex_for_render[1][indices[0]][2] << std::endl;
	////////std::cout << vertex_for_render[1][indices[1]][0] << " " << vertex_for_render[1][indices[1]][1] << " " << vertex_for_render[1][indices[1]][2] << std::endl;
	////////std::cout << vertex_for_render[1][indices[2]][0] << " " << vertex_for_render[1][indices[2]][1] << " " << vertex_for_render[1][indices[2]][2] << std::endl;

	////////std::cout << start_pair_index<<" "<< end_pair_index << std::endl;




	//for (int i = start_pair_index; i < end_pair_index; i += 4) {
	//	stiffness = stiffness_initial;
	//	pair = pair_ + i;
	//	indices = triangle_indices[*(pair + 3)][*(pair + 2)].data();

	//	//if (pair[0] == 1 && pair[1] == 0 && pair[2] == 0 && pair[3] == 1) {
	//	//	//std::cout << "occurs in response section " << std::endl;
	//	//}
	//	//if (pair[1] == 0 && pair[0] == 1 && pair[2] == 0 && pair[3] == 1) {
	//	//	//std::cout << indices[0] << " " << indices[1] << " " << indices[2] << std::endl;
	//	//	std::cout<<"compute position "<<  vertex_for_render[pair[ 1]][pair[0]][0] << " " << vertex_for_render[pair[ 1]][pair[0]][1] << " " << vertex_for_render[pair[ 1]][pair[0]][2] << std::endl;
	//	//	//std::cout << vertex_for_render[pair[ 3]][indices[0]][0] << " " << vertex_for_render[pair[ 3]][indices[0]][1] << " " << vertex_for_render[pair[ 3]][indices[0]][2] << std::endl;
	//	//	//std::cout << vertex_for_render[pair[ 3]][indices[1]][0] << " " << vertex_for_render[pair[ 3]][indices[1]][1] << " " << vertex_for_render[pair[ 3]][indices[1]][2] << std::endl;
	//	//	//std::cout << vertex_for_render[pair[ 3]][indices[2]][0] << " " << vertex_for_render[pair[ 3]][indices[2]][1] << " " << vertex_for_render[pair[ 3]][indices[2]][2] << std::endl;
	//	//}


	//	//if (*(pair + 3) == 1 && *(pair + 2) == 73 && *pair == 18 && *(pair + 1) == 0) {

	//	//		//std::cout << "collision_time here" << collision_time[i >> 2] << std::endl;

	//	//}




	//	if (collision_constraint.pointTriangleResponse(vertex_for_render[*(pair+ 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
	//		vertex_for_render[*(pair + 3)][indices[0]].data(), vertex_for_render[*(pair + 3)][indices[1]].data(),
	//		vertex_for_render[*(pair + 3)][indices[2]].data(),
	//		vertex_position[*(pair + 3)][indices[0]].data(), vertex_position[*(pair + 3)][indices[1]].data(),
	//		vertex_position[*(pair + 3)][indices[2]].data(), 
	//		vertex_collision_free[*(pair + 1)][*pair].data(),
	//		vertex_collision_free[*(pair + 3)][indices[0]].data(), vertex_collision_free[*(pair + 3)][indices[1]].data(),
	//		vertex_collision_free[*(pair + 3)][indices[2]].data(),
	//		triangle_normal_render[*(pair + 3)][*(pair + 2)].data(),
	//		target_pos_v, target_pos_v+3, target_pos_v+6, target_pos_v+9, d_hat_, stiffness, epsilon_,
	//		mass[*(pair + 1)][*pair],
	//		mass[*(pair + 3)][indices[0]], mass[*(pair + 3)][indices[1]], mass[*(pair + 3)][indices[2]], pair,
	//		collision_time[i>>2], this->collision_time))
	//	{



	//		//if (*(pair + 3) == 1 &&(* (pair + 2) == 15|| *(pair + 1) == 0 || *(pair +0) == 10)) {
	//		//if (*(pair + 3) == 0 && (indices[0] == 10 || indices[1] ==10 || indices[2] == 10)) {
	//		//	std::cout << *pair << " " << *(pair + 1) << " " << *(pair + 2) << " " << *(pair + 3) <<" "<< indices[0]
	//		//		<< " " << indices[1] <<" "<<indices[2]<<" "<<collision_time[i>>2] << std::endl;
	//		//	for (unsigned int k = 0; k < 3; ++k) {
	//		//		if (indices[k] == 10) {
	//		//			double* pos = target_pos_v + 3 * k + 3;
	//		//			std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//		//			pos = vertex_collision_free[*(pair + 3)][indices[k]].data();
	//		//			std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//		//			pos = vertex_position[*(pair + 3)][indices[k]].data();
	//		//			std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//		//		}
	//		//	}
	//		//}
	//		//if (*(pair + 1) == 0 && *( pair +0)== 10) {
	//			//std::cout << *pair << " " << *(pair + 1) << " " << *(pair + 2) << " " << *(pair + 3) << " " << indices[0]
	//			//	<< " " << indices[1] << " " << indices[2] << std::endl;
	//			//double* pos = target_pos_v;
	//			//std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//			//pos = vertex_collision_free[*(pair + 1)][*pair].data();
	//			//std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//			//pos = vertex_position[*(pair + 3)][*pair].data();
	//			//std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//		


	//		//if ( *(pair + 1) == 0 && *(pair + 0) == 10) {//* (pair + 3) == 1 && *(pair + 2) == 15 &&
	//		//	//if (*(pair + 3) == 0 && (indices[0] == 10 || indices[1] ==10 || indices[2] == 10)) {
	//		//	//	//std::cout << *pair << " " << *(pair + 1) << " " << *(pair + 2) << " " << *(pair + 3) <<" "<< indices[0]
	//		//	//		<< " " << indices[1] <<" "<<indices[2] << std::endl;
	//		//	//	for (unsigned int k = 0; k < 3; ++k) {
	//		//	//		if (indices[k] == 0) {
	//		//	//			double* pos = target_pos_v + 3 * k + 3;
	//		//	//			//std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//		//	//			pos = vertex_collision_free[*(pair + 3)][indices[k]].data();
	//		//	//			//std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//		//	//			pos = vertex_position[*(pair + 3)][indices[k]].data();
	//		//	//			//std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//		//	//		}
	//		//	//	}
	//		//	//}
	//		//	//if (*(pair + 1) == 0 && *( pair +0)== 10) {
	//		//	std::cout << *pair << " " << *(pair + 1) << " " << *(pair + 2) << " " << *(pair + 3) << " " << indices[0]
	//		//		<< " " << indices[1] << " " << indices[2] << std::endl;
	//		//	double* pos = vertex_collision_free[*(pair + 1)][*pair].data();
	//		//	std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//		//	pos = target_pos_v;
	//		//	std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//		//	pos = vertex_position[*(pair + 1)][*pair].data();
	//		//	std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	//		//}


	//		addTargetPosToSystemTotal(vertex_b_sum[*(pair + 1)][*pair].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*pair].data(),
	//			target_pos_v, stiffness, record_stiffness[*(pair + 1)][*pair], vetex_need_update[*(pair + 1)][*pair]);
	//		for (int j = 0; j < 3; ++j) {
	//			addTargetPosToSystemTotal(vertex_b_sum[*(pair + 3)][indices[j]].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][indices[j]].data(),
	//				target_pos_v+3*(j+1), stiffness, record_stiffness[*(pair + 3)][indices[j]], vetex_need_update[*(pair + 3)][indices[j]]);
	//		}
	//		memcpy(target_pos_index_, pair, 16);
	//		target_pos_index_ += 4;
	//		memcpy(target_pos_record, target_pos_v, 96);
	//		target_pos_record[12] = stiffness;
	//		target_pos_record += 13;
	//		////std::cout << "stiffness " << stiffness << std::endl;


	//	}
	//}


	//point_triangle_target_pos_index[thread_No][0] = target_pos_index_ - point_triangle_target_pos_index[thread_No].data() - 1;
}




void Collision::pointTriangleResponse(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, TargetPosition* target_pos)
{
	//////std::cout << thread_No<<" "<< start_pair_index << " " << end_pair_index << std::endl;
	unsigned int* target_pos_index_ = point_triangle_target_pos_index[thread_No].data() + 1 + point_triangle_target_pos_index[thread_No][0];
	//unsigned int* target_pos_index_ = target_position_index[thread_No].data() + 1 + (target_position_index[thread_No][0] << 1);
	//double* target_pos_ = target_position_and_stiffness[thread_No].data() + (target_position_index[thread_No][0]<<2);



	unsigned int* pair_ = spatial_hashing.vertex_triangle_pair[pair_thread_No] + 1;
	//double d_hat_ = d_hat;
	double tolerance = tolerance_radius[SELF_POINT_TRIANGLE];
	double epsilon_ = epsilon;
	int* indices;
	double target_pos_v[3];
	double target_pos_tri[3][3];
	double stiffness_initial = collision_stiffness[SELF_POINT_TRIANGLE];
	double stiffness;
	std::vector<std::array<double,3>>* vertex_b_sum = target_pos->b_sum.data();
	std::vector<double>* record_stiffness = target_pos->stiffness.data();
	bool** vetex_need_update = target_pos->need_update;
	
	unsigned int* pair;

	for (int i = start_pair_index; i < end_pair_index; i+=4) {
		stiffness = stiffness_initial;
		pair = pair_ + i;
		indices = triangle_indices[*(pair+3)][*(pair+ 2)].data();
		//if (collision_constraint.pointTriangleResponse(vertex_for_render[*(pair+ 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
		//	vertex_for_render[*(pair + 3)][indices[0]].data(), vertex_for_render[pair[i + 3]][indices[1]].data(),
		//	vertex_for_render[*(pair + 3)][indices[2]].data(),
		//	vertex_position[*(pair + 3)][indices[0]].data(), vertex_position[*(pair + 3)][indices[1]].data(),
		//	vertex_position[*(pair + 3)][indices[2]].data(), triangle_normal_render[*(pair + 3)][*(pair + 2)].data(),
		//	target_pos_v, target_pos_tri[0], target_pos_tri[1], target_pos_tri[2], d_hat_, stiffness, epsilon_,
		//	mass[*(pair + 1)][*pair],
		//	mass[*(pair + 3)][indices[0]], mass[*(pair + 3)][indices[1]], mass[*(pair + 3)][indices[2]]))
		if (dcd.pointSelfTriangle(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
			vertex_for_render[*(pair + 3)][indices[0]].data(), vertex_for_render[*(pair + 3)][indices[1]].data(),
			vertex_for_render[*(pair + 3)][indices[2]].data(),
			vertex_position[*(pair + 3)][indices[0]].data(), vertex_position[*(pair + 3)][indices[1]].data(),
			vertex_position[*(pair + 3)][indices[2]].data(), triangle_normal_render[*(pair + 3)][*(pair + 2)].data(),
			triangle_normal[*(pair + 3)][*(pair + 2)].data(),
			target_pos_v, target_pos_tri[0], target_pos_tri[1], target_pos_tri[2], tolerance,
			mass[*(pair + 1)][*pair], mass[*(pair + 3)][indices[0]], mass[*(pair + 3)][indices[1]], mass[*(pair + 3)][indices[2]],
			//1.0,1.0,1.0,1.0,
			triangle_normal_magnitude_reciprocal[*(pair + 3)][*(pair + 2)]))
		{
			//////std::cout << "should not occur point triangle collision" << std::endl;
			/*////std::cout << "===========" << std::endl;
			////std::cout << *(pair + 1) << " "<<*(pair + 3)<<" "<< *(pair)<<" "<< indices[0]<<" "<< indices[1]<<" "<< indices[2] << std::endl;
			////std::cout << vertex_for_render[*(pair + 1)][*pair][0] << " " << vertex_for_render[*(pair + 1)][*pair][1] << " "
				<< vertex_for_render[*(pair + 1)][*pair][2] << std::endl;
			////std::cout << vertex_position[*(pair + 1)][*pair][0] << " " << vertex_position[*(pair + 1)][*pair][1] << " "
				<< vertex_position[*(pair + 1)][*pair][2] << std::endl;
			for (int j = 0; j < 3; ++j) {
				////std::cout << vertex_for_render[*(pair + 3)][indices[j]][0] << " " << vertex_for_render[*(pair + 3)][indices[j]][1] << " "
					<< vertex_for_render[*(pair + 3)][indices[j]][2] << std::endl;
				////std::cout << vertex_position[*(pair + 3)][indices[j]][0] << " " << vertex_position[*(pair + 3)][indices[j]][1] << " "
					<< vertex_position[*(pair + 3)][indices[j]][2] << std::endl;
			}*/
			addTargetPosToSystemTotal(vertex_b_sum[*(pair + 1)][*pair].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*pair].data(),
				target_pos_v, stiffness, record_stiffness[*(pair + 1)][*pair], vetex_need_update[*(pair + 1)][*pair]);

			//*(target_pos_index_++) = *pair;
			//*(target_pos_index_++) = *(pair + 1);

			memcpy(target_pos_index_, pair, 16);
			target_pos_index_ += 4;


			//memcpy(target_pos_, target_pos_v, 24);
			//target_pos_ += 3;
			//*(target_pos_++) = stiffness;
			for (int j = 0; j < 3; ++j) {
				addTargetPosToSystemTotal(vertex_b_sum[*(pair + 3)][indices[j]].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][indices[j]].data(),
					target_pos_tri[j], stiffness, record_stiffness[*(pair + 3)][indices[j]], vetex_need_update[*(pair + 3)][indices[j]]);

				//*(target_pos_index_++) = indices[j];
				//*(target_pos_index_++) = *(pair + 3);

				//memcpy(target_pos_, target_pos_tri[j], 24);
				//target_pos_ += 3;
				//*(target_pos_++) = stiffness;
			}
		}
	}

	point_triangle_target_pos_index[thread_No][0] = target_pos_index_ - point_triangle_target_pos_index[thread_No].data() - 1;
	//target_position_index[thread_No][0] = (target_pos_index_ - target_position_index[thread_No].data() - 1) >> 1;
}



void Collision::re_solveXPBDpointTriangleColliderResponse(unsigned int pair_thread_No, int index_start, int index_end)
{
	unsigned int* target_pos_index_record_ = point_triangle_collider_target_pos_index[pair_thread_No].data() + 1;
	double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	int* indices;
	double target_pos_v[3];
	double target_pos_tri[3][3];
	unsigned int* pair;
	double friction_coe_collider = friction_coe[1];
	//double* lambda_ = lambda->data() + collision_lambda_index_start[0];
	for (int i = index_start; i < index_end; i += 4) {
		pair = target_pos_index_record_ + i;
		indices = triangle_indices_collider[*(pair + 3)][*(pair + 2)].data();
		dcd.XPBDpointTriangleCollider(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
			vertex_for_render_collider[*(pair + 3)][indices[0]].data(), vertex_for_render_collider[*(pair + 3)][indices[1]].data(),
			vertex_for_render_collider[*(pair + 3)][indices[2]].data(),
			vertex_position_collider[*(pair + 3)][indices[0]].data(), vertex_position_collider[*(pair + 3)][indices[1]].data(),
			vertex_position_collider[*(pair + 3)][indices[2]].data(),
			triangle_normal_render_collider[*(pair + 3)][*(pair + 2)].data(),
			triangle_normal_collider[*(pair + 3)][*(pair + 2)].data(), tolerance, friction_coe_collider);
	}
}



void Collision::solveXPBDpointTriangleColliderResponse(unsigned int thread_No, unsigned int pair_thread_No, int index_start, int index_end)
{
	unsigned int* pair_ = spatial_hashing.vertex_obj_triangle_collider_pair[pair_thread_No] + 1;
	unsigned int* target_pos_index_ = point_triangle_collider_target_pos_index[thread_No].data() + 1 + point_triangle_collider_target_pos_index[thread_No][0];
	double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	int* indices;
	double target_pos_v[3];
	double target_pos_tri[3][3];

	unsigned int* pair;
	//double* lambda_ = lambda->data() + collision_lambda_index_start[0];
	double friction_coe_collider = friction_coe[1];
	for (int i = index_start; i < index_end; i += 4) {
		pair = pair_ + i;
		indices = triangle_indices_collider[*(pair + 3)][*(pair + 2)].data();
		if (dcd.XPBDpointTriangleCollider(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
			vertex_for_render_collider[*(pair + 3)][indices[0]].data(), vertex_for_render_collider[*(pair + 3)][indices[1]].data(),
			vertex_for_render_collider[*(pair + 3)][indices[2]].data(),
			vertex_position_collider[*(pair + 3)][indices[0]].data(), vertex_position_collider[*(pair + 3)][indices[1]].data(),
			vertex_position_collider[*(pair + 3)][indices[2]].data(),
			triangle_normal_render_collider[*(pair + 3)][*(pair + 2)].data(),
			triangle_normal_collider[*(pair + 3)][*(pair + 2)].data(), tolerance, friction_coe_collider)) {
			memcpy(target_pos_index_, pair, 16);
			target_pos_index_ += 4;
		}
	}

	point_triangle_collider_target_pos_index[thread_No][0] = target_pos_index_ - point_triangle_collider_target_pos_index[thread_No].data() - 1;
}


void Collision::re_pointTriangleColliderResponseForIPC(unsigned int pair_thread_No, int index_start, int index_end, TargetPosition* target_pos)
{
	unsigned int* target_pos_index_record_ = point_triangle_collider_target_pos_index[pair_thread_No].data() + 1;
	double* target_pos_record_ = point_triangle_collider_target_pos_record[pair_thread_No].data();

	double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	double epsilon_ = epsilon;
	int* indices;
	double* target_pos_v;
	double stiffness;
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	std::vector<double>* record_stiffness = target_pos->stiffness.data();
	unsigned int* pair;

	for (int i = index_start; i < index_end; i += 4) {
		target_pos_v = target_pos_record_ + i;
		stiffness = target_pos_v[3];
		pair = target_pos_index_record_ + i;

		indices = triangle_indices_collider[*(pair + 3)][*(pair + 2)].data();		
		addTargetPosToSystem(vertex_b_sum[*(pair + 1)][*pair].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*pair].data(),
				target_pos_v, stiffness);
	}
}


void Collision::re_pointTriangleColliderResponse(unsigned int pair_thread_No, int index_start, int index_end, TargetPosition* target_pos)
{
	unsigned int* target_pos_index_record_ = point_triangle_collider_target_pos_index[pair_thread_No].data() + 1;
	double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	double epsilon_ = epsilon;
	int* indices;
	double target_pos_v[3];
	double target_pos_tri[3][3];
	double stiffness_initial = collision_stiffness[BODY_POINT_TRIANGLE];
	double stiffness;
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	std::vector<double>* record_stiffness = target_pos->stiffness.data();
	unsigned int* pair;

	for (int i = index_start; i < index_end; i += 4) {
		stiffness = stiffness_initial;
		pair = target_pos_index_record_ + i;
		indices = triangle_indices_collider[*(pair + 3)][*(pair + 2)].data();
		//if (collision_constraint.pointTriangleResponse(vertex_for_render[*(pair+ 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
		//	vertex_for_render[*(pair + 3)][indices[0]].data(), vertex_for_render[pair[i + 3]][indices[1]].data(),
		//	vertex_for_render[*(pair + 3)][indices[2]].data(),
		//	vertex_position[*(pair + 3)][indices[0]].data(), vertex_position[*(pair + 3)][indices[1]].data(),
		//	vertex_position[*(pair + 3)][indices[2]].data(), triangle_normal_render[*(pair + 3)][*(pair + 2)].data(),
		//	target_pos_v, target_pos_tri[0], target_pos_tri[1], target_pos_tri[2], d_hat_, stiffness, epsilon_,
		//	mass[*(pair + 1)][*pair],
		//	mass[*(pair + 3)][indices[0]], mass[*(pair + 3)][indices[1]], mass[*(pair + 3)][indices[2]]))
		if (dcd.pointTriangleCollider(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
			vertex_for_render_collider[*(pair + 3)][indices[0]].data(), vertex_for_render_collider[*(pair+ 3)][indices[1]].data(),
			vertex_for_render_collider[*(pair + 3)][indices[2]].data(),
			vertex_position_collider[*(pair + 3)][indices[0]].data(), vertex_position_collider[*(pair + 3)][indices[1]].data(),
			vertex_position_collider[*(pair + 3)][indices[2]].data(),
			triangle_normal_render_collider[*(pair + 3)][*(pair + 2)].data(),
			triangle_normal_collider[*(pair + 3)][*(pair + 2)].data(),
			target_pos_v, tolerance))
		{
			addTargetPosToSystem(vertex_b_sum[*(pair + 1)][*pair].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*pair].data(),
				target_pos_v, stiffness);
		}
		else {
			construct_b_sum(vertex_b_sum[*(pair + 1)][*pair].data(),
				vertex_position[*(pair + 1)][*pair].data(), stiffness);
		}
	}
}


void Collision::pointTriangleColliderPair(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index)
{
	unsigned int* target_pos_index_ = point_triangle_collider_target_pos_index[thread_No].data() + 1 + point_triangle_collider_target_pos_index[thread_No][0];
	unsigned int* pair_ = spatial_hashing.vertex_obj_triangle_collider_pair[pair_thread_No] + 1;
	//double d_hat_ = d_hat;
	double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	int* indices;

	unsigned int* pair;

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		pair = pair_ + i;
		indices = triangle_indices_collider[*(pair + 3)][*(pair + 2)].data();
		if (dcd.checkPointTriangle(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
			vertex_for_render_collider[*(pair + 3)][indices[0]].data(), vertex_for_render_collider[*(pair + 3)][indices[1]].data(),
			vertex_for_render_collider[*(pair + 3)][indices[2]].data(),
			vertex_position_collider[*(pair + 3)][indices[0]].data(), vertex_position_collider[*(pair + 3)][indices[1]].data(),
			vertex_position_collider[*(pair + 3)][indices[2]].data(),
			triangle_normal_render_collider[*(pair + 3)][*(pair + 2)].data(),tolerance))
		{
			memcpy(target_pos_index_, pair, 16);
			target_pos_index_ += 4;
		}
	}
	point_triangle_collider_target_pos_index[thread_No][0] = target_pos_index_ - point_triangle_collider_target_pos_index[thread_No].data() - 1;
}


void Collision::pointTriangleColliderResponseForIPC(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, TargetPosition* target_pos)
{
	//unsigned int* target_pos_index_ = point_triangle_collider_target_pos_index[thread_No].data() + 1 + point_triangle_collider_target_pos_index[thread_No][0];
	//double* target_pos_record= point_triangle_collider_target_pos_record[thread_No].data() + point_triangle_collider_target_pos_index[thread_No][0];

	//double* collision_time = record_VTCollider_collision_time[pair_thread_No].data();

	//unsigned int* pair_ = spatial_hashing.vertex_obj_triangle_collider_pair[pair_thread_No] + 1;
	//double d_hat_ = d_hat;
	//double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	//double epsilon_ = epsilon;
	//int* indices;
	//double target_pos_v[3];
	//double stiffness;
	//double stiffness_initial = collision_stiffness[BODY_POINT_TRIANGLE];
	//std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	//std::vector<double>* record_stiffness = target_pos->stiffness.data();
	//bool** vetex_need_update = target_pos->need_update;

	//unsigned int* pair;

	//for (int i = start_pair_index; i < end_pair_index; i += 4) {
	//	stiffness = stiffness_initial;
	//	pair = pair_ + i;
	//	indices = triangle_indices_collider[*(pair + 3)][*(pair + 2)].data();
	//	if (collision_constraint.pointTriangleColliderResponse(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
	//		vertex_collision_free[*(pair + 1)][*pair].data(),
	//		vertex_for_render_collider[*(pair + 3)][indices[0]].data(), vertex_for_render_collider[pair[i + 3]][indices[1]].data(),
	//		vertex_for_render_collider[*(pair + 3)][indices[2]].data(),
	//		triangle_normal_render_collider[*(pair + 3)][*(pair + 2)].data(),
	//		target_pos_v, d_hat_, stiffness, epsilon_, *pair, collision_time[i>>2]))
	//	{
	//		////std::cout << *(pair + 2)<<" "<< vertex_for_render[*(pair + 1)][*pair][1] << " " << target_pos_v[1]<<" "<<stiffness << std::endl;

	//		addTargetPosToSystemTotal(vertex_b_sum[*(pair + 1)][*pair].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*pair].data(),
	//			target_pos_v, stiffness, record_stiffness[*(pair + 1)][*pair], vetex_need_update[*(pair + 1)][*pair]);
	//		memcpy(target_pos_index_, pair, 16);
	//		target_pos_index_ += 4;
	//		memcpy(target_pos_record, target_pos_v, 24);
	//		target_pos_record[3] = stiffness;
	//		target_pos_record += 4;

	//		if (*pair == chosen_show_vertex) {
	//			draw_target_position.clear();
	//			draw_target_position.push_back({ target_pos_v[0],target_pos_v[1],target_pos_v[2] });
	//		}

	//	}
	//}
	//point_triangle_collider_target_pos_index[thread_No][0] = target_pos_index_ - point_triangle_collider_target_pos_index[thread_No].data() - 1;
}



void Collision::pointTriangleColliderResponse(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, TargetPosition* target_pos)
{
	unsigned int* target_pos_index_ = point_triangle_collider_target_pos_index[thread_No].data() + 1 + point_triangle_collider_target_pos_index[thread_No][0];
	//unsigned int* target_pos_index_ = target_position_index[thread_No].data() + 1 + (target_position_index[thread_No][0] << 1);
	//double* target_pos_ = target_position_and_stiffness[thread_No].data() + (target_position_index[thread_No][0] <<2);

	unsigned int* pair_ = spatial_hashing.vertex_obj_triangle_collider_pair[pair_thread_No] + 1;
	//double d_hat_ = d_hat;
	double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	double epsilon_ = epsilon;
	int* indices;
	double target_pos_v[3];
	double stiffness;
	double stiffness_initial = collision_stiffness[BODY_POINT_TRIANGLE];
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	std::vector<double>* record_stiffness = target_pos->stiffness.data();
	bool** vetex_need_update = target_pos->need_update;

	unsigned int* pair;

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		stiffness = stiffness_initial;
		pair = pair_ + i;
		indices = triangle_indices_collider[*(pair + 3)][*(pair + 2)].data();
		//if (collision_constraint.pointTriangleColliderResponse(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
		//	vertex_for_render_collider[*(pair + 3)][indices[0]].data(), vertex_for_render_collider[pair[i + 3]][indices[1]].data(),
		//	vertex_for_render_collider[*(pair + 3)][indices[2]].data(),
		//	triangle_normal_render_collider[*(pair + 3)][*(pair + 2)].data(),
		//	target_pos_v, d_hat_, stiffness, epsilon_))
		if (dcd.pointTriangleCollider(vertex_for_render[*(pair + 1)][*pair].data(), vertex_position[*(pair + 1)][*pair].data(),
			vertex_for_render_collider[*(pair + 3)][indices[0]].data(), vertex_for_render_collider[*(pair+ 3)][indices[1]].data(),
			vertex_for_render_collider[*(pair + 3)][indices[2]].data(),
			vertex_position_collider[*(pair + 3)][indices[0]].data(), vertex_position_collider[*(pair + 3)][indices[1]].data(),
			vertex_position_collider[*(pair + 3)][indices[2]].data(),
			triangle_normal_render_collider[*(pair + 3)][*(pair + 2)].data(),
			triangle_normal_collider[*(pair + 3)][*(pair + 2)].data(),
			target_pos_v, tolerance))
		{
			addTargetPosToSystemTotal(vertex_b_sum[*(pair + 1)][*pair].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*pair].data(),
				target_pos_v, stiffness, record_stiffness[*(pair + 1)][*pair], vetex_need_update[*(pair + 1)][*pair]);

			memcpy(target_pos_index_, pair, 16);
			target_pos_index_ += 4;



			//*(target_pos_index_++) = *pair;
			//*(target_pos_index_++) = *(pair + 1);
			//memcpy(target_pos_, target_pos_v, 24);
			//target_pos_ += 3;
			//*(target_pos_++) = stiffness;
		}
	}

	//target_position_index[thread_No][0] = (target_pos_index_ - target_position_index[thread_No].data() - 1) >> 1;
	point_triangle_collider_target_pos_index[thread_No][0] = target_pos_index_ - point_triangle_collider_target_pos_index[thread_No].data() - 1;
}




void Collision::edgeEdgeCollisionPair(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index)
{
	//unsigned int* target_pos_index_ = target_position_index[thread_No].data() + 1 + (target_position_index[thread_No][0] << 1);
	unsigned int* target_pos_index_ = edge_edge_target_pos_index[thread_No].data() + 1 + edge_edge_target_pos_index[thread_No][0];
	//double* target_pos_ = target_position_and_stiffness[thread_No].data() + (target_position_index[thread_No][0] << 2);
	unsigned int* pair_ = spatial_hashing.edge_edge_pair[pair_thread_No] + 1;
	//double d_hat_ = d_hat;
	double tolerance = tolerance_radius[SELF_EDGE_EDGE];
	unsigned int* indices;
	unsigned int* compare_indices;
	unsigned int* pair;
	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		pair = pair_ + i;
		indices = edge_vertices[*(pair + 1)] + ((*(pair)) << 1);
		compare_indices = edge_vertices[*(pair + 3)] + ((*(pair + 2)) << 1);
		//if (collision_constraint.edgeEdgeResponse(target_pos_edge_0, target_pos_edge_1, target_pos_compare_edge_0, target_pos_compare_edge_1,
		//	vertex_position[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
		//	vertex_for_render[*(pair + 1)][*indices].data(), vertex_for_render[*(pair + 1)][*(indices + 1)].data(),
		//	vertex_position[*(pair + 3)][*compare_indices].data(), vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
		//	vertex_for_render[*(pair + 3)][*compare_indices].data(), vertex_for_render[*(pair + 3)][*(compare_indices + 1)].data(),
		//	d_hat_, stiffness, epsilon_,
		//	mass[*(pair + 1)][*indices], mass[*(pair + 1)][*(indices + 1)],
		//	mass[*(pair + 3)][*compare_indices], mass[*(pair + 3)][*(compare_indices + 1)]))
		if (dcd.checkEdgeEdge(
			vertex_position[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
			vertex_for_render[*(pair + 1)][*indices].data(), vertex_for_render[*(pair + 1)][*(indices + 1)].data(),
			vertex_position[*(pair + 3)][*compare_indices].data(), vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
			vertex_for_render[*(pair + 3)][*compare_indices].data(), vertex_for_render[*(pair + 3)][*(compare_indices + 1)].data(),
			tolerance))
		{
			memcpy(target_pos_index_, pair, 16);
			target_pos_index_ += 4;
		}
	}
	edge_edge_target_pos_index[thread_No][0] = target_pos_index_ - edge_edge_target_pos_index[thread_No].data() - 1;
}



void Collision::edgeEdgeResponseForIPC(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, TargetPosition* target_pos)
{
	//unsigned int* target_pos_index_ = edge_edge_target_pos_index[thread_No].data() + 1 + edge_edge_target_pos_index[thread_No][0];
	//double* target_pos_index_record = edge_edge_target_pos_record[thread_No].data() + (edge_edge_target_pos_index[thread_No][0]>>2)*13;

	//unsigned int* pair_ = spatial_hashing.edge_edge_pair[pair_thread_No] + 1;
	//double d_hat_ = d_hat;
	//double epsilon_ = epsilon;
	//unsigned int* indices;
	//unsigned int* compare_indices;
	//double target_pos_edge[12];
	//double stiffness;
	//double stiffness_initial = collision_stiffness[SELF_EDGE_EDGE];
	//std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	//std::vector<double>* record_stiffness = target_pos->stiffness.data();
	//bool** vetex_need_update = target_pos->need_update;
	//unsigned int* pair;
	//double* collision_time = record_EE_collision_time[pair_thread_No].data();
	//for (int i = start_pair_index; i < end_pair_index; i += 4) {
	//	stiffness = stiffness_initial;
	//	pair = pair_ + i;
	//	indices = edge_vertices[*(pair + 1)] + ((*pair) << 1);
	//	compare_indices = edge_vertices[*(pair + 3)] + ((*(pair + 2)) << 1);

	//	if (collision_constraint.edgeEdgeResponse(target_pos_edge, target_pos_edge+3, target_pos_edge+6, target_pos_edge+9,
	//		vertex_position[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
	//		vertex_for_render[*(pair + 1)][*indices].data(), vertex_for_render[*(pair + 1)][*(indices + 1)].data(),
	//		vertex_position[*(pair + 3)][*compare_indices].data(), vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
	//		vertex_for_render[*(pair + 3)][*compare_indices].data(), vertex_for_render[*(pair + 3)][*(compare_indices + 1)].data(),
	//		vertex_collision_free[*(pair + 1)][*indices].data(), vertex_collision_free[*(pair + 1)][*(indices + 1)].data(),
	//		vertex_collision_free[*(pair + 3)][*compare_indices].data(), vertex_collision_free[*(pair + 3)][*(compare_indices + 1)].data(),
	//		d_hat_, stiffness, epsilon_,
	//		mass[*(pair + 1)][*indices], mass[*(pair + 1)][*(indices + 1)],
	//		mass[*(pair + 3)][*compare_indices], mass[*(pair + 3)][*(compare_indices + 1)], collision_time[i>>2], pair))
	//	{
	//		addTargetPosToSystemTotal(vertex_b_sum[*(pair + 1)][*indices].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*indices].data(),
	//			target_pos_edge, stiffness, record_stiffness[*(pair + 1)][*indices], vetex_need_update[*(pair + 1)][*indices]);
	//		addTargetPosToSystemTotal(vertex_b_sum[*(pair + 1)][*(indices + 1)].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*(indices + 1)].data(),
	//			target_pos_edge+3, stiffness, record_stiffness[*(pair + 1)][*(indices + 1)], vetex_need_update[*(pair + 1)][*(indices + 1)]);
	//		addTargetPosToSystemTotal(vertex_b_sum[*(pair + 3)][*compare_indices].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][*compare_indices].data(),
	//			target_pos_edge+6, stiffness, record_stiffness[*(pair + 3)][*compare_indices], vetex_need_update[*(pair + 3)][*compare_indices]);
	//		addTargetPosToSystemTotal(vertex_b_sum[*(pair + 3)][*(compare_indices + 1)].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
	//			target_pos_edge+9, stiffness, record_stiffness[*(pair + 3)][*(compare_indices + 1)], vetex_need_update[*(pair + 3)][*(compare_indices + 1)]);
	//		memcpy(target_pos_index_, pair, 16);
	//		target_pos_index_ += 4;
	//		memcpy(target_pos_index_record, target_pos_edge, 96);
	//		target_pos_index_record[12] = stiffness;
	//		target_pos_index_record += 13;

	//	}
	//}
	//edge_edge_target_pos_index[thread_No][0] = target_pos_index_ - edge_edge_target_pos_index[thread_No].data() - 1;
}


void Collision::edgeEdgeResponse(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, TargetPosition* target_pos)
{
	//unsigned int* target_pos_index_ = target_position_index[thread_No].data() + 1 + (target_position_index[thread_No][0] << 1);
	unsigned int* target_pos_index_ = edge_edge_target_pos_index[thread_No].data() + 1 + edge_edge_target_pos_index[thread_No][0];
	//double* target_pos_ = target_position_and_stiffness[thread_No].data() + (target_position_index[thread_No][0] << 2);
	unsigned int* pair_ = spatial_hashing.edge_edge_pair[pair_thread_No] + 1;
	//double d_hat_ = d_hat;
	double tolerance = tolerance_radius[SELF_EDGE_EDGE];
	double epsilon_ = epsilon;
	unsigned int* indices;
	unsigned int* compare_indices;
	double target_pos_edge_0[3];
	double target_pos_edge_1[3];
	double target_pos_compare_edge_0[3];
	double target_pos_compare_edge_1[3];
	double stiffness;
	double stiffness_initial = collision_stiffness[SELF_EDGE_EDGE];
	std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	std::vector<double>* record_stiffness = target_pos->stiffness.data();
	bool** vetex_need_update = target_pos->need_update;
	unsigned int* pair;
	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		stiffness = stiffness_initial;
		pair = pair_ + i;
		indices = edge_vertices[*(pair + 1)] + ((*(pair)) << 1);
		compare_indices = edge_vertices[*(pair + 3)] + ((*(pair + 2)) << 1);
		//if (collision_constraint.edgeEdgeResponse(target_pos_edge_0, target_pos_edge_1, target_pos_compare_edge_0, target_pos_compare_edge_1,
		//	vertex_position[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
		//	vertex_for_render[*(pair + 1)][*indices].data(), vertex_for_render[*(pair + 1)][*(indices + 1)].data(),
		//	vertex_position[*(pair + 3)][*compare_indices].data(), vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
		//	vertex_for_render[*(pair + 3)][*compare_indices].data(), vertex_for_render[*(pair + 3)][*(compare_indices + 1)].data(),
		//	d_hat_, stiffness, epsilon_,
		//	mass[*(pair + 1)][*indices], mass[*(pair + 1)][*(indices + 1)],
		//	mass[*(pair + 3)][*compare_indices], mass[*(pair + 3)][*(compare_indices + 1)]))
		if (dcd.edgeEdge(target_pos_edge_0, target_pos_edge_1, target_pos_compare_edge_0, target_pos_compare_edge_1,
			vertex_position[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
			vertex_for_render[*(pair + 1)][*indices].data(), vertex_for_render[*(pair + 1)][*(indices + 1)].data(),
			vertex_position[*(pair + 3)][*compare_indices].data(), vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
			vertex_for_render[*(pair + 3)][*compare_indices].data(), vertex_for_render[*(pair + 3)][*(compare_indices + 1)].data(),
			tolerance, 
			//1.0, 1.0, 1.0, 1.0
			mass[*(pair + 1)][*indices], mass[*(pair + 1)][*(indices + 1)],	mass[*(pair + 3)][*compare_indices], mass[*(pair + 3)][*(compare_indices + 1)]
		))
		{
			addTargetPosToSystemTotal(vertex_b_sum[*(pair + 1)][*indices].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*indices].data(),
				target_pos_edge_0, stiffness, record_stiffness[*(pair + 1)][*indices], vetex_need_update[*(pair + 1)][*indices]);
			addTargetPosToSystemTotal(vertex_b_sum[*(pair + 1)][*(indices + 1)].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*(indices + 1)].data(),
				target_pos_edge_1, stiffness, record_stiffness[*(pair + 1)][*(indices + 1)], vetex_need_update[*(pair + 1)][*(indices + 1)]);
			addTargetPosToSystemTotal(vertex_b_sum[*(pair + 3)][*compare_indices].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][*compare_indices].data(),
				target_pos_compare_edge_0, stiffness, record_stiffness[*(pair + 3)][*compare_indices], vetex_need_update[*(pair + 3)][*compare_indices]);
			addTargetPosToSystemTotal(vertex_b_sum[*(pair + 3)][*(compare_indices + 1)].data(), target_pos->collision_energy, vertex_position[*(pair + 3)][*(compare_indices + 1)].data(),
				target_pos_compare_edge_1, stiffness, record_stiffness[*(pair + 3)][*(compare_indices + 1)], vetex_need_update[*(pair + 3)][*(compare_indices + 1)]);
			memcpy(target_pos_index_, pair, 16);
			target_pos_index_ += 4;

			//*(target_pos_index_++) = *indices;
			//*(target_pos_index_++) = *(pair + 1);
			//memcpy(target_pos_, target_pos_edge_0, 24);
			//target_pos_ += 3;
			//*(target_pos_++) = stiffness;
			//*(target_pos_index_++) = *(indices + 1);
			//*(target_pos_index_++) = *(pair + 1);
			//memcpy(target_pos_, target_pos_edge_1, 24);
			//target_pos_ += 3;
			//*(target_pos_++) = stiffness;
			//*(target_pos_index_++) = *compare_indices;
			//*(target_pos_index_++) = *(pair + 3);
			//memcpy(target_pos_, target_pos_compare_edge_0, 24);
			//target_pos_ += 3;
			//*(target_pos_++) = stiffness;
			//*(target_pos_index_++) = *(compare_indices + 1);
			//*(target_pos_index_++) = *(pair + 3);
			//memcpy(target_pos_, target_pos_compare_edge_1, 24);
			//target_pos_ += 3;
			//*(target_pos_++) = stiffness;
			//test_edge_edge_record_true_number[thread_No]++;
		}
	}
	edge_edge_target_pos_index[thread_No][0] = target_pos_index_ - edge_edge_target_pos_index[thread_No].data() - 1;
	//target_position_index[thread_No][0] = (target_pos_index_ - target_position_index[thread_No].data() - 1);
}
void Collision::edgeEdgeColliderResponse(unsigned int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, TargetPosition* target_pos)
{
	////unsigned int* target_pos_index_ = target_position_index[thread_No].data() + 1 + (target_position_index[thread_No][0] << 1);
	//unsigned int* target_pos_index_ = edge_edge_collider_target_pos_index[thread_No].data() + 1 + edge_edge_collider_target_pos_index[thread_No][0];
	////double* target_pos_ = target_position_and_stiffness[thread_No].data() + (target_position_index[thread_No][0] << 2);
	//unsigned int* pair_ = spatial_hashing.edge_edge_pair_collider[pair_thread_No] + 1;
	////double d_hat_ = d_hat;
	//double tolerance = tolerance_radius[BODY_POINT_TRIANGLE];
	//unsigned int* indices;
	//unsigned int* compare_indices;
	//double target_pos_edge_0[3];
	//double target_pos_edge_1[3];
	//double stiffness;
	//double stiffness_initial = collision_stiffness[BODY_POINT_TRIANGLE];
	//std::vector<std::array<double, 3>>* vertex_b_sum = target_pos->b_sum.data();
	//std::vector<double>* record_stiffness = target_pos->stiffness.data();
	//bool** vetex_need_update = target_pos->need_update;
	//unsigned int* pair;
	//stiffness = stiffness_initial;
	//for (int i = start_pair_index; i < end_pair_index; i += 4) {		
	//	pair = pair_ + i;
	//	indices = edge_vertices[*(pair + 1)]+ ((*pair) << 1);
	//	compare_indices = collider_edge_vertices[*(pair + 3)] + (2*(*(pair + 2)));
	//	if (dcd.edgeEdgeCollider(target_pos_edge_0, target_pos_edge_1,
	//		vertex_position[*(pair + 1)][*indices].data(), vertex_position[*(pair + 1)][*(indices + 1)].data(),
	//		vertex_for_render[*(pair + 1)][*indices].data(), vertex_for_render[*(pair + 1)][*(indices + 1)].data(),
	//		vertex_position_collider[*(pair + 3)][*compare_indices].data(), vertex_position_collider[*(pair + 3)][*(compare_indices + 1)].data(),
	//		vertex_for_render_collider[*(pair + 3)][*compare_indices].data(), vertex_for_render_collider[*(pair + 3)][*(compare_indices + 1)].data(),
	//		tolerance, mass[*(pair + 1)][*indices],mass[*(pair + 1)][*(indices + 1)]))
	//	{
	//		addTargetPosToSystemTotal(vertex_b_sum[*(pair + 1)][*indices].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*indices].data(),
	//			target_pos_edge_0, stiffness, record_stiffness[*(pair + 1)][*indices], vetex_need_update[*(pair + 1)][*indices]);
	//		addTargetPosToSystemTotal(vertex_b_sum[*(pair + 1)][*(indices + 1)].data(), target_pos->collision_energy, vertex_position[*(pair + 1)][*(indices + 1)].data(),
	//			target_pos_edge_1, stiffness, record_stiffness[*(pair + 1)][*(indices + 1)], vetex_need_update[*(pair + 1)][*(indices + 1)]);
	//		memcpy(target_pos_index_, pair, 16);
	//		target_pos_index_ += 4;
	//	}
	//}
	//edge_edge_collider_target_pos_index[thread_No][0] = target_pos_index_ - edge_edge_collider_target_pos_index[thread_No].data() - 1;
}


//COLLISION_CONSTRAINT_IPC
//void Collision::collisionConstraint(int thread_No)
//{
//	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
//	target_pos->initial();
//	TriangleMeshStruct* mesh_struct;
//	std::vector<std::vector<int>>* neighbor_primitve;
//	int index_begin;
//	int index_end;
//	double* mass;
//	std::array<double, 3>* current_pos;
//	std::array<double, 3>* initial_pos;
//	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
//		mesh_struct = &(*cloth)[cloth_No].mesh_struct;
//		index_begin = mesh_struct->vertex_index_begin_per_thread[thread_No];
//		index_end = mesh_struct->vertex_index_begin_per_thread[thread_No + 1];
//		neighbor_primitve = (*cloth)[cloth_No].vertex_neighbor_obj_triangle.data();
//		mass = mesh_struct->mass.data();
//		current_pos = mesh_struct->vertex_position.data();
//		initial_pos = mesh_struct->vertex_for_render.data();
//		for (int i = index_begin; i < index_end; ++i) {
//			pointSelfTriangleClose(neighbor_primitve[i].data(), initial_pos[i].data(), current_pos[i].data(), i, cloth_No, mass[i], target_pos);
//		}
//	}
//	std::array<int, 3>* triangle_vertex_index;
//	std::vector<double*> triangle_pos(3);
//	std::array<double, 3>* triangle_normal;
//	for (int collider_No = 0; collider_No < collider->size(); ++collider_No) {
//		mesh_struct = &(*collider)[collider_No].mesh_struct;
//		index_begin = mesh_struct->face_index_begin_per_thread[thread_No];
//		index_end = mesh_struct->face_index_begin_per_thread[thread_No + 1];
//		current_pos = mesh_struct->vertex_position.data();
//		neighbor_primitve = (*collider)[collider_No].triangle_neighbor_obj_vertex.data();
//		triangle_vertex_index = mesh_struct->triangle_indices.data();
//		triangle_normal = mesh_struct->face_normal.data();
//		for (int i = index_begin; i < index_end; ++i) {
//			for (int j = 0; j < 3; ++j) {
//				triangle_pos[j] = current_pos[triangle_vertex_index[i][j]].data();
//			}
//			pointColliderTriangleClose(triangle_vertex_index[i].data(), neighbor_primitve[i].data(), triangle_pos, triangle_normal[i].data(), target_pos);
//		}
//	}
//	unsigned int* edge_vertex_index;
//	double mass_[4];
//	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
//		mesh_struct = &(*cloth)[cloth_No].mesh_struct;
//		index_begin = mesh_struct->edge_index_begin_per_thread[thread_No];
//		index_end = mesh_struct->edge_index_begin_per_thread[thread_No + 1];
//		neighbor_primitve = (*cloth)[cloth_No].edge_neighbor_obj_edge.data();
//		initial_pos = mesh_struct->vertex_for_render.data();
//		current_pos = mesh_struct->vertex_position.data();
//		for (int i = index_begin; i < index_end; ++i) {
//			edge_vertex_index = mesh_struct->edge_vertices.data() + (i << 1);// edges[i].vertex;
//			mass_[0] = mesh_struct->mass[edge_vertex_index[0]];
//			mass_[1] = mesh_struct->mass[edge_vertex_index[1]];
//			edgeEdgeClose(neighbor_primitve[i].data(), initial_pos[edge_vertex_index[0]].data(), initial_pos[edge_vertex_index[1]].data(),
//				current_pos[edge_vertex_index[0]].data(), current_pos[edge_vertex_index[1]].data(), cloth_No, edge_vertex_index[0],
//				edge_vertex_index[1], mass_, target_pos);
//		}
//	}
//}


void Collision::vertexColliderTriangleCollisionTime(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, double& collision_time)
{
	//double conservative_rescaling_ = conservative_rescaling;
	//unsigned int* pair = spatial_hashing.vertex_collider_triangle_obj_pair[pair_thread_No] + 1;
	//double time;
	//int* indices;
	//for (int i = start_pair_index; i < end_pair_index; i += 4) {
	//	indices = triangle_indices[pair[i + 3]][pair[i + 2]].data();
	//	if (approx_CCD.vertexTriangleCCD(time, vertex_for_render_collider[pair[i + 1]][pair[i]].data(), vertex_position_collider[pair[i + 1]][pair[i]].data(),
	//		vertex_for_render[pair[i + 3]][indices[0]].data(), vertex_position[pair[i + 3]][indices[0]].data(), vertex_for_render[pair[i + 3]][indices[1]].data(),
	//		vertex_position[pair[i + 3]][indices[1]].data(), vertex_for_render[pair[i + 3]][indices[2]].data(), vertex_position[pair[i + 3]][indices[2]].data(),
	//		conservative_rescaling_)) {
	//		if (time < collision_time) {
	//			collision_time = time;
	//		}
	//	}
	//	//time = CCD::pointTriangleCcd(vertex_for_render_collider[pair[i + 1]][pair[i]].data(),
	//	//	vertex_for_render[pair[i + 3]][indices[0]].data(),
	//	//	vertex_for_render[pair[i + 3]][indices[1]].data(),
	//	//	vertex_for_render[pair[i + 3]][indices[2]].data(),
	//	//	vertex_position_collider[pair[i + 1]][pair[i]].data(),
	//	//	vertex_position[pair[i + 3]][indices[0]].data(),
	//	//	vertex_position[pair[i + 3]][indices[1]].data(),
	//	//	vertex_position[pair[i + 3]][indices[2]].data(), eta, tolerance);
	//	//if (time < collision_time) {
	//	//	collision_time = time;
	//	//}
	//}
}


void Collision::vertexTriangleColliderCollisionTime(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, double& collision_time)
{
	double conservative_rescaling_ = conservative_rescaling;
	unsigned int* pair = spatial_hashing.vertex_obj_triangle_collider_pair[pair_thread_No] + 1;
	double time;
	int* indices;

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		indices = triangle_indices_collider[pair[i + 3]][pair[i + 2]].data();

		//if (approx_CCD.pointTriangleCollisionTime(time, vertex_for_render[pair[i + 1]][pair[i]].data(), vertex_position[pair[i + 1]][pair[i]].data(),
		//	vertex_for_render_collider[pair[i + 3]][indices[0]].data(), vertex_position_collider[pair[i + 3]][indices[0]].data(), vertex_for_render_collider[pair[i + 3]][indices[1]].data(),
		//	vertex_position_collider[pair[i + 3]][indices[1]].data(), vertex_for_render_collider[pair[i + 3]][indices[2]].data(), vertex_position_collider[pair[i + 3]][indices[2]].data(),
		//	triangle_normal_render_not_normalized_collider[pair[i + 3]][pair[i + 2]].data(), triangle_normal_not_normalized_collider[pair[i + 3]][pair[i + 2]].data(),
		//	cross_for_approx_CCD_collider[pair[i + 3]][pair[i + 2]].data(), tolerance_)) {
		//	if (time < collision_time) {
		//		collision_time = time;
		//	}
		//}
		//if (approx_CCD.vertexTriangleCCD(time, vertex_for_render[pair[i + 1]][pair[i]].data(), vertex_position[pair[i + 1]][pair[i]].data(),
		//	vertex_for_render_collider[pair[i + 3]][indices[0]].data(), vertex_position_collider[pair[i + 3]][indices[0]].data(), vertex_for_render_collider[pair[i + 3]][indices[1]].data(),
		//	vertex_position_collider[pair[i + 3]][indices[1]].data(), vertex_for_render_collider[pair[i + 3]][indices[2]].data(), vertex_position_collider[pair[i + 3]][indices[2]].data(), conservative_rescaling_)) {
		//	if (time < collision_time) {
		//		collision_time = time;
		//	}
		//}
		time = CCD::pointTriangleCcd(vertex_for_render[pair[i + 1]][pair[i]].data(),
			vertex_for_render_collider[pair[i + 3]][indices[0]].data(),
			vertex_for_render_collider[pair[i + 3]][indices[1]].data(),
			vertex_for_render_collider[pair[i + 3]][indices[2]].data(),
			vertex_position[pair[i + 1]][pair[i]].data(),
			vertex_position_collider[pair[i + 3]][indices[0]].data(),
			vertex_position_collider[pair[i + 3]][indices[1]].data(),
			vertex_position_collider[pair[i + 3]][indices[2]].data(), eta, tolerance);
		//if (pair[i] == chosen_show_vertex) {
		//	//std::cout << "collision time " << time << std::endl;
		//}
		if (time < collision_time) {
			collision_time = time;
		}
	}
}


void Collision::findVertexEdgePair(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, CollisionIndicateType type)
{
	unsigned int* pair;
	unsigned int* indices;
	int edge_num;
	unsigned int* vertex_edge_pair_ ;

	unsigned int** face_edges_;
	unsigned int** representative_edge_num_;
	std::array<double, 6>** vertex_aabb_;
	std::array<double, 6>** edge_aabb_;	
	

	if (type == VERTEX_EDGE) {
		pair = spatial_hashing.vertex_triangle_pair[pair_thread_No] + 1;
		vertex_edge_pair_ = vertex_edge_pair[thread_No] + 1 + vertex_edge_pair[thread_No][0];
		representative_edge_num_= representative_edge_num.data();
		face_edges_ = face_edges.data();
		vertex_aabb_ = vertex_aabb.data();
		edge_aabb_ = edge_aabb.data();

	}
	else{// if(type == VERTEX_OBJ_EDGE_COLLIDER){
		pair = spatial_hashing.vertex_obj_triangle_collider_pair[pair_thread_No] + 1;
		vertex_edge_pair_ = vertex_obj_edge_collider_pair[thread_No] + 1 + vertex_obj_edge_collider_pair[thread_No][0];
		representative_edge_num_ = representative_edge_num_collider.data();
		face_edges_ = collider_face_edges.data();
		vertex_aabb_ = vertex_aabb.data();
		edge_aabb_ = edge_aabb_collider.data();
	}
	//else {
		//pair = spatial_hashing.vertex_collider_triangle_obj_pair[pair_thread_No] + 1;
		//vertex_edge_pair_ = vertex_collider_edge_obj_pair[thread_No] + 1 + vertex_collider_edge_obj_pair[thread_No][0];
		//representative_edge_num_ = representative_edge_num.data();
		//face_edges_ = face_edges.data();
		//vertex_aabb_ = vertex_aabb_collider.data();
		//edge_aabb_ = edge_aabb.data();
	//}

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		edge_num = representative_edge_num_[pair[i + 3]][pair[i + 2]];
		indices = face_edges_[pair[i + 3]] + 3 * pair[i + 2];
		for (int j = 0; j < edge_num; ++j) {
			if (AABB::AABB_intersection(vertex_aabb_[pair[i + 1]][pair[i]].data(), edge_aabb_[pair[i + 3]][indices[j]].data())) {
				memcpy(vertex_edge_pair_, pair + i, 8);
				vertex_edge_pair_ += 2;
				*vertex_edge_pair_ = indices[j];
				vertex_edge_pair_++;
				*vertex_edge_pair_ = pair[i + 3];
				vertex_edge_pair_++;
			}
		}
	}

	if (type == VERTEX_EDGE) {
		vertex_edge_pair[thread_No][0] = vertex_edge_pair_ - vertex_edge_pair[thread_No] - 1;
	}
	else if (type == VERTEX_OBJ_EDGE_COLLIDER) {
		vertex_obj_edge_collider_pair[thread_No][0] = vertex_edge_pair_ - vertex_obj_edge_collider_pair[thread_No] - 1;
	}
	else {
		vertex_collider_edge_obj_pair[thread_No][0] = vertex_edge_pair_ - vertex_collider_edge_obj_pair[thread_No] - 1;
	}
	
	
}


void Collision::findVertexVertexPair(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, CollisionIndicateType type)
{
	unsigned int* pair; 
	int* indices;
	int vertex_num;
	unsigned int* vertex_vertex_pair_;
	unsigned int** representative_vertex_num_;
	std::array<double, 6>** vertex_aabb_;
	std::array<int, 3>** triangle_index_in_order_;
	if (type == VERTEX_VERTEX) {
		vertex_vertex_pair_ = vertex_vertex_pair[thread_No] + 1 + vertex_vertex_pair[thread_No][0];
		representative_vertex_num_ = representative_vertex_num.data();
		vertex_aabb_ = vertex_aabb.data();
		triangle_index_in_order_ = triangle_index_in_order.data();
		pair = spatial_hashing.vertex_triangle_pair[pair_thread_No] + 1;		
	}
	else {
		vertex_vertex_pair_ = vertex_vertex_pair_collider[thread_No] + 1 + vertex_vertex_pair_collider[thread_No][0];
		representative_vertex_num_ = representative_vertex_num_collider.data();
		vertex_aabb_ = vertex_aabb_collider.data();
		triangle_index_in_order_ = triangle_index_in_order_collider.data();
		pair = spatial_hashing.vertex_obj_triangle_collider_pair[pair_thread_No] + 1;
	}

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		vertex_num = representative_vertex_num_[pair[i + 3]][pair[i + 2]];
		indices = triangle_index_in_order_[pair[i + 3]][pair[i + 2]].data();
		for (int j = 0; j < vertex_num; ++j) {
			if (AABB::AABB_intersection(vertex_aabb[pair[i + 1]][pair[i]].data(), vertex_aabb_[pair[i + 3]][indices[j]].data())) {
				memcpy(vertex_vertex_pair_, pair + i, 8);
				vertex_vertex_pair_ += 2;
				*vertex_vertex_pair_ = indices[j];
				vertex_vertex_pair_++;
				*vertex_vertex_pair_ = pair[i + 3];
				vertex_vertex_pair_++;
			}
		}
	}

	if (type == VERTEX_VERTEX) {
		vertex_vertex_pair[thread_No][0] = vertex_vertex_pair_ - vertex_vertex_pair[thread_No] - 1;
	}
	else {
		vertex_vertex_pair_collider[thread_No][0] = vertex_vertex_pair_ - vertex_vertex_pair_collider[thread_No] - 1;
	}

}

void Collision::vertexEdgeCollisionTime(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, double& collision_time, CollisionIndicateType type)
{
	double conservative_rescaling_ = conservative_rescaling;
	unsigned int* pair;
	double time;
	unsigned int** edge_vertices_;
	std::array<double, 3>** vertex_for_render_vertex_;
	std::array<double, 3>** vertex_position_vertex_;
	std::array<double, 3>** vertex_for_render_edge_;
	std::array<double, 3>** vertex_position_edge_;
	if (type == VERTEX_EDGE) {
		pair = vertex_edge_pair[pair_thread_No] + 1;
		vertex_for_render_vertex_ = vertex_for_render.data();
		vertex_position_vertex_ = vertex_position.data();
		vertex_for_render_edge_ = vertex_for_render.data();
		vertex_position_edge_ = vertex_position.data();
		edge_vertices_ = edge_vertices.data();
	}
	else if(type==VERTEX_COLLIDER_EDGE_OBJ) {
		pair = vertex_collider_edge_obj_pair[pair_thread_No] + 1;
		vertex_for_render_vertex_ = vertex_for_render_collider.data();
		vertex_position_vertex_ = vertex_position_collider.data();
		vertex_for_render_edge_ = vertex_for_render.data();
		vertex_position_edge_ = vertex_position.data();
		edge_vertices_ = edge_vertices.data();
	}
	else {
		pair = vertex_obj_edge_collider_pair[pair_thread_No] + 1;
		vertex_for_render_edge_ = vertex_for_render_collider.data();
		vertex_position_edge_ = vertex_position_collider.data();
		vertex_for_render_vertex_ = vertex_for_render.data();
		vertex_position_vertex_ = vertex_position.data();
		edge_vertices_ = collider_edge_vertices.data();
	}
	unsigned int* indices;

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		indices = edge_vertices_[pair[i + 3]] + (pair[i + 2] << 1);
		if (approx_CCD.vertexEdgeCCD(time, vertex_for_render_vertex_[pair[i + 1]][pair[i]].data(), vertex_for_render_edge_[pair[i + 3]][indices[0]].data(), 
			vertex_for_render_edge_[pair[i + 3]][indices[1]].data(),vertex_position_vertex_[pair[i + 1]][pair[i]].data(),
			vertex_position_edge_[pair[i + 3]][indices[0]].data(), vertex_position_edge_[pair[i + 3]][indices[1]].data(), conservative_rescaling_)) {
			if (time < collision_time) {
				collision_time = time;
			}
		}
	}

}


void Collision::vertexEdgeCollisionTimeCompare(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, double& collision_time, CollisionIndicateType type)
{
	double conservative_rescaling_ = conservative_rescaling;
	unsigned int* pair;
	double time;
	unsigned int** edge_vertices_;
	std::vector<Eigen::Vector3d>* vertex_for_render_vertex_;
	std::vector<Eigen::Vector3d>* vertex_position_vertex_;
	std::vector<Eigen::Vector3d>* vertex_for_render_edge_;
	std::vector<Eigen::Vector3d>* vertex_position_edge_;

	pair = vertex_edge_pair[pair_thread_No] + 1;
	vertex_for_render_vertex_ = vertex_for_render_eigen.data();
	vertex_position_vertex_ = vertex_position_eigen.data();
	vertex_for_render_edge_ = vertex_for_render_eigen.data();
	vertex_position_edge_ = vertex_position_eigen.data();
	edge_vertices_ = edge_vertices.data();

	unsigned int* indices;

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		indices = edge_vertices_[pair[i + 3]] + (pair[i + 2] << 1);
		if (approx_CCD.vertexEdgeCCD(time, vertex_for_render_vertex_[pair[i + 1]][pair[i]], vertex_for_render_edge_[pair[i + 3]][indices[0]],
			vertex_for_render_edge_[pair[i + 3]][indices[1]], vertex_position_vertex_[pair[i + 1]][pair[i]],
			vertex_position_edge_[pair[i + 3]][indices[0]], vertex_position_edge_[pair[i + 3]][indices[1]], conservative_rescaling_)) {
			if (time < collision_time) {
				collision_time = time;
			}
		}
	}

}



void Collision::vertexVertexCollisionTime(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, double& collision_time, CollisionIndicateType type)
{
	double conservative_rescaling_ = conservative_rescaling;
	unsigned int* pair;
	double time;
	std::array<double, 3>** vertex_for_render_;
	std::array<double, 3>** vertex_position_;
	if (type == VERTEX_VERTEX) {
		pair = vertex_vertex_pair[pair_thread_No] + 1;
		vertex_for_render_ = vertex_for_render.data();
		vertex_position_ = vertex_position.data();
	}
	else {
		pair = vertex_vertex_pair_collider[pair_thread_No] + 1;
		vertex_for_render_ = vertex_for_render_collider.data();
		vertex_position_ = vertex_position_collider.data();
	}
	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		if (approx_CCD.vertexVertexCCD(time, vertex_for_render[pair[i + 1]][pair[i]].data(), vertex_position[pair[i + 1]][pair[i]].data(),
			vertex_for_render_[pair[i + 3]][pair[i+2]].data(), vertex_position_[pair[i + 3]][pair[i + 2]].data(),conservative_rescaling_)) {
			if (time < collision_time) {
				collision_time = time;
			}
		}
	}
}



void Collision::vertexVertexCollisionTimeCompare(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, double& collision_time, CollisionIndicateType type)
{
	double conservative_rescaling_ = conservative_rescaling;
	unsigned int* pair;
	double time;
	std::vector<Eigen::Vector3d>* vertex_for_render_;
	std::vector<Eigen::Vector3d>* vertex_position_;
	pair = vertex_vertex_pair[pair_thread_No] + 1;
	vertex_for_render_ = vertex_for_render_eigen.data();
	vertex_position_ = vertex_position_eigen.data();

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		if (approx_CCD.vertexVertexCCD(time, vertex_for_render_eigen[pair[i + 1]][pair[i]], vertex_position_eigen[pair[i + 1]][pair[i]],
			vertex_for_render_[pair[i + 3]][pair[i + 2]], vertex_position_[pair[i + 3]][pair[i + 2]], conservative_rescaling_)) {
			if (time < collision_time) {
				collision_time = time;
			}
		}
	}
}


void Collision::vertexTriangleCollisionTimeCompare(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, double& collision_time)
{
	double conservative_rescaling_ = conservative_rescaling;
	unsigned int* pair = spatial_hashing.vertex_triangle_pair[pair_thread_No] + 1;
	double time;
	int* indices;
	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		indices = triangle_indices[pair[i + 3]][pair[i + 2]].data();
		if (approx_CCD.vertexTriangleCCD(time, vertex_for_render_eigen[pair[i + 1]][pair[i]], vertex_position_eigen[pair[i + 1]][pair[i]],
			vertex_for_render_eigen[pair[i + 3]][indices[0]], vertex_position_eigen[pair[i + 3]][indices[0]], vertex_for_render_eigen[pair[i + 3]][indices[1]],
			vertex_position_eigen[pair[i + 3]][indices[1]], vertex_for_render_eigen[pair[i + 3]][indices[2]], vertex_position_eigen[pair[i + 3]][indices[2]], conservative_rescaling_)) {
			if (time < collision_time) {
				collision_time = time;
			}
		}

	}
}


void Collision::EECollisionTimeOneEdgeAll(double* initial_pos_a0, double* initial_pos_a1, double* current_pos_a0, double* current_pos_a1, double& collision_time,
	unsigned int num, unsigned int* edge_indices,
	std::array<double, 3>** vertex_for_render, std::array<double, 3>** vertex_pos, unsigned int** compare_edge_vertices, int move_size)
{
	double time;
	unsigned int* indices;

	for (int i = 0; i < num; i += move_size) {
		indices = compare_edge_vertices[edge_indices[i]] + (edge_indices[i + 1] << 1);
		time = CCD::edgeEdgeCcd(initial_pos_a0, initial_pos_a1, vertex_for_render[edge_indices[i]][*indices].data(),
			vertex_for_render[edge_indices[i]][*(indices + 1)].data(),
			current_pos_a0, current_pos_a1, vertex_pos[edge_indices[i]][*indices].data(),
			vertex_pos[edge_indices[i]][*(indices + 1)].data(), eta, tolerance);
		if (time < collision_time) {
			collision_time = time;
		}		
	}
}



void Collision::EECollisionTimeOneEdge(double* initial_pos_a0, double* initial_pos_a1, double* current_pos_a0, double* current_pos_a1, double& collision_time,
	unsigned int edge_index, unsigned int edge_obj_No, unsigned int num, unsigned int* edge_indices,
	std::array<double, 3>** vertex_for_render, std::array<double, 3>** vertex_pos, bool is_self, unsigned int** edge_vertices)
{
	double time;
	unsigned int* indices;
	if (is_self) {
		for (int i = 0; i < num; i += 3) {
			if (edge_obj_No < edge_indices[i]) {
				indices = edge_vertices[edge_indices[i]] + (edge_indices[i + 1] << 1);
				time = CCD::edgeEdgeCcd(initial_pos_a0, initial_pos_a1, vertex_for_render[edge_indices[i]][*indices].data(),
					vertex_for_render[edge_indices[i]][*(indices + 1)].data(),
					current_pos_a0, current_pos_a1, vertex_pos[edge_indices[i]][*indices].data(),
					vertex_pos[edge_indices[i]][*(indices + 1)].data(), eta, tolerance);
				if (time < collision_time) {
					collision_time = time;
				}
			}
			else if (edge_obj_No == edge_indices[i] && edge_index < edge_indices[i + 1]) {
				indices = edge_vertices[edge_indices[i]] + (edge_indices[i + 1] << 1);
				time = CCD::edgeEdgeCcd(initial_pos_a0, initial_pos_a1, vertex_for_render[edge_indices[i]][*indices].data(),
					vertex_for_render[edge_indices[i]][*(indices + 1)].data(),
					current_pos_a0, current_pos_a1, vertex_pos[edge_indices[i]][*indices].data(),
					vertex_pos[edge_indices[i]][*(indices + 1)].data(), eta, tolerance);
				if (time < collision_time) {
					collision_time = time;
				}
			}
		}
	}
	else {
		for (int i = 0; i < num; i += 2) {
			indices = edge_vertices[edge_indices[i]] + (edge_indices[i + 1] << 1);
			time = CCD::edgeEdgeCcd(initial_pos_a0, initial_pos_a1, vertex_for_render[edge_indices[i]][*indices].data(),
				vertex_for_render[edge_indices[i]][*(indices + 1)].data(),
				current_pos_a0, current_pos_a1, vertex_pos[edge_indices[i]][*indices].data(),
				vertex_pos[edge_indices[i]][*(indices + 1)].data(), eta, tolerance);
			if (time < collision_time) {
				collision_time = time;
			}
		}
	}

}


void Collision::EECollisionTimeOneEdge(double* initial_pos_a0, double* initial_pos_a1, double* current_pos_a0, double* current_pos_a1, double& collision_time, 
	unsigned int edge_index, unsigned int edge_obj_No,  unsigned int num, unsigned int* edge_indices, 
	std::array<double, 3>** vertex_for_render, std::array<double, 3>** vertex_pos, bool is_self, unsigned int** edge_vertices, bool* is_used)
{
	double time;
	unsigned int* indices;
	if (is_self) {
		for (int i = 0; i < num; i += 2) {
			if (!(*is_used)) {
				//if (edge_obj_No < edge_indices[i]) {
					indices = edge_vertices[edge_indices[i]] + (edge_indices[i + 1] << 1);
					time = CCD::edgeEdgeCcd(initial_pos_a0, initial_pos_a1, vertex_for_render[edge_indices[i]][*indices].data(),
						vertex_for_render[edge_indices[i]][*(indices + 1)].data(),
						current_pos_a0, current_pos_a1, vertex_pos[edge_indices[i]][*indices].data(),
						vertex_pos[edge_indices[i]][*(indices + 1)].data(), eta, tolerance);
					if (time < collision_time) {
						collision_time = time;
					}
				//}
				//else if (edge_obj_No == edge_indices[i] && edge_index < edge_indices[i + 1]) {
				//	indices = edge_vertices[edge_indices[i]] + (edge_indices[i + 1] << 1);
				//	time = CCD::edgeEdgeCcd(initial_pos_a0, initial_pos_a1, vertex_for_render[edge_indices[i]][*indices].data(),
				//		vertex_for_render[edge_indices[i]][*(indices + 1)].data(),
				//		current_pos_a0, current_pos_a1, vertex_pos[edge_indices[i]][*indices].data(),
				//		vertex_pos[edge_indices[i]][*(indices + 1)].data(), eta, tolerance);
				//	if (time < collision_time) {
				//		collision_time = time;
				//	}
				//}
			}
		}
		is_used++;
	}
	else {
		for (int i = 0; i < num; i += 2) {
			if (!(*is_used)) {
				indices = edge_vertices[edge_indices[i]] + (edge_indices[i + 1] << 1);
				time = CCD::edgeEdgeCcd(initial_pos_a0, initial_pos_a1, vertex_for_render[edge_indices[i]][*indices].data(),
					vertex_for_render[edge_indices[i]][*(indices + 1)].data(),
					current_pos_a0, current_pos_a1, vertex_pos[edge_indices[i]][*indices].data(),
					vertex_pos[edge_indices[i]][*(indices + 1)].data(), eta, tolerance);
				if (time < collision_time) {
					collision_time = time;
				}
			}
			is_used++;
		}
	}

	memset(is_used - (num >> 1), 0, num >> 1);
}





bool Collision::floorCollisionTime(double* initial_position, double* current_pos, unsigned int dimension, bool direction,
	double floor_value, double& collision_time, double tolerance)
{
	if (direction) {
		floor_value += tolerance;
		if (current_pos[dimension] > floor_value) {
			return false;
		}
		floor_value += tolerance;
	}
	else {
		floor_value -= tolerance;
		if (current_pos[dimension] < floor_value) {
			return false;
		}
		
	}
	double time = (initial_position[dimension] - floor_value) / (initial_position[dimension] - current_pos[dimension]);
	if (collision_time > time) {
		collision_time = time;
	}
	if (time < 0) {
		std::cout << "error collision time of floor is negative" << std::endl;
	}

	if (direction) {
		if (initial_position[dimension] < floor_value) {
			std::cout << "error exceed floor" << std::endl;
		}
	}

	return true;
}


void Collision::TVCollisionTimeOneTriangle(double* initial_pos_0, double* initial_pos_1, double* initial_pos_2,
	double* current_pos_0, double* current_pos_1, double* current_pos_2,
	double& collision_time, unsigned int num,
	unsigned int* vertex_index, std::array<double, 3>** initial_vertex, std::array<double, 3>** current_vertex, bool* is_used, int size_of_a_pair)
{
	double time;
	int* indices;
	for (int i = 0; i < num; i += size_of_a_pair) {
		if (!(* is_used)) {
			time = CCD::pointTriangleCcd(initial_vertex[vertex_index[i]][vertex_index[i + 1]].data(),
				initial_pos_0, initial_pos_1, initial_pos_2,
				current_vertex[vertex_index[i]][vertex_index[i + 1]].data(),
				current_pos_0, current_pos_1, current_pos_2, eta, tolerance);
			if (time < collision_time) {
				collision_time = time;
			}
		}
		is_used++;
	}
	memset(is_used - (num / size_of_a_pair), 0, num/ size_of_a_pair);
}

void Collision::TVCollisionTimeOneTriangleSelfColor(double* initial_pos_0, double* initial_pos_1, double* initial_pos_2,
	double* current_pos_0, double* current_pos_1, double* current_pos_2,
	double& collision_time, unsigned int num,
	unsigned int* vertex_index, std::array<double, 3>** initial_vertex, std::array<double, 3>** current_vertex, int size_of_a_pair)
{
	double time;
	for (int i = 0; i < num; i += size_of_a_pair) {
		//if (!vertex_belong_to_color_group[vertex_index[i]][vertex_index[i + 1]]) {
			time = CCD::pointTriangleCcd(initial_vertex[vertex_index[i]][vertex_index[i + 1]].data(),
				initial_pos_0, initial_pos_1, initial_pos_2,
				current_vertex[vertex_index[i]][vertex_index[i + 1]].data(),
				current_pos_0, current_pos_1, current_pos_2, eta, tolerance);
			if (time < collision_time) {
				collision_time = time;
			}
		//}
	}	
}


void Collision::TVCollisionTimeOneTriangle(double* initial_pos_0, double* initial_pos_1, double* initial_pos_2,
	double* current_pos_0, double* current_pos_1, double* current_pos_2,
	double& collision_time, unsigned int num,
	unsigned int* vertex_index, std::array<double, 3>** initial_vertex, std::array<double, 3>** current_vertex, int size_of_a_pair)
{
	double time;
	for (int i = 0; i < num; i += size_of_a_pair) {
		time = CCD::pointTriangleCcd(initial_vertex[vertex_index[i]][vertex_index[i + 1]].data(),
			initial_pos_0, initial_pos_1, initial_pos_2,
			current_vertex[vertex_index[i]][vertex_index[i + 1]].data(),
			current_pos_0, current_pos_1, current_pos_2, eta, tolerance);
		if (time < collision_time) {
			collision_time = time;
		}
	}
}



void Collision::VTCollisionTimeOneVertex(double* initial_pos, double* current_pos, double& collision_time, unsigned int num,
	unsigned int* triangle_index, std::array<double, 3>** initial_vertex, std::array<double, 3>** current_vertex, std::array<int, 3>** triangle_indices, bool* is_used)
{
	double time;
	int* indices;
	for (int i = 0; i < num; i += 2) {
		if (!(*is_used)) {
			indices = triangle_indices[triangle_index[i]][triangle_index[i + 1]].data();
			time = CCD::pointTriangleCcd(initial_pos,
				initial_vertex[triangle_index[i]][indices[0]].data(),
				initial_vertex[triangle_index[i]][indices[1]].data(),
				initial_vertex[triangle_index[i]][indices[2]].data(),
				current_pos,
				current_vertex[triangle_index[i]][indices[0]].data(),
				current_vertex[triangle_index[i]][indices[1]].data(),
				current_vertex[triangle_index[i]][indices[2]].data(), eta, tolerance);
			if (time < collision_time) {
				collision_time = time;
			}
		}
		is_used++;
	}
	memset(is_used - (num >> 1), 0, num >> 1);
}



void Collision::VTCollisionTimeOneVertex(double* initial_pos, double* current_pos, double& collision_time, unsigned int num,
	unsigned int* triangle_index, std::array<double, 3>** initial_vertex, std::array<double, 3>** current_vertex, std::array<int, 3>** triangle_indices, bool is_collider, int vertex_index)
{
	double time;
	int* indices;
	for (int i = 0; i < num; i += 2) {
		indices = triangle_indices[triangle_index[i]][triangle_index[i + 1]].data();
		time = CCD::pointTriangleCcd(initial_pos,
			initial_vertex[triangle_index[i]][indices[0]].data(),
			initial_vertex[triangle_index[i]][indices[1]].data(),
			initial_vertex[triangle_index[i]][indices[2]].data(),
			current_pos,
			current_vertex[triangle_index[i]][indices[0]].data(),
			current_vertex[triangle_index[i]][indices[1]].data(),
			current_vertex[triangle_index[i]][indices[2]].data(), eta, tolerance);
		if (time < collision_time) {
			collision_time = time;

			if (time == 0.0) {
				std::cout << vertex_index<<" "<< triangle_index[i] << " " << triangle_index[i + 1] << std::endl;
			}

		}		
	}
}


void Collision::edgeEdgeCollisionTimePair(int start_pair_index,
	int end_pair_index, double& collision_time, unsigned int* pair, std::array<double, 3>** vertex_for_render_0, std::array<double, 3>** vertex_position_0,
	std::array<double, 3>** vertex_for_render_1, std::array<double, 3>** vertex_position_1, 
	std::array<double, 3>** ori_pos_0, std::array<double, 3>** ori_pos_1,
	unsigned int** edge_vertices_0, unsigned int** edge_vertices_1, std::vector<unsigned int>* record_index,
	std::vector<double>* record_d_hat, unsigned int* hash_size_record, unsigned int* hash_record,
	unsigned int* edge_0_index_prefix_sum, unsigned int* edge_1_index_prefix_sum)
{
	//unsigned int* pair = spatial_hashing.edge_edge_pair[pair_thread_No] + 1;
	double time;
	unsigned int* indices_1;
	unsigned int* indices_0;

	double distance = d_hat_2;
	unsigned int hash_value;
	unsigned int actual_ele_0, actual_ele_1;
	unsigned int* address;

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		indices_0 = edge_vertices_0[pair[i + 1]] + (pair[i] << 1);
		indices_1 = edge_vertices_1[pair[i + 3]] + (pair[i + 2] << 1);

		time = CCD::edgeEdgeCcd(vertex_for_render_0[pair[i + 1]][indices_0[0]].data(),
			vertex_for_render_0[pair[i + 1]][indices_0[1]].data(),
			vertex_for_render_1[pair[i + 3]][indices_1[0]].data(),
			vertex_for_render_1[pair[i + 3]][indices_1[1]].data(),
			vertex_position_0[pair[i + 1]][indices_0[0]].data(),
			vertex_position_0[pair[i + 1]][indices_0[1]].data(),
			vertex_position_1[pair[i + 3]][indices_1[0]].data(),
			vertex_position_1[pair[i + 3]][indices_1[1]].data(), eta, tolerance);

		if (time < 1.0) {
			actual_ele_0 = edge_0_index_prefix_sum[pair[i + 1]] + pair[i];
			actual_ele_1 = edge_1_index_prefix_sum[pair[i + 3]] + pair[i + 2];
			hash_value = ((actual_ele_0 * P1) ^ (P2 * actual_ele_1)) % pair_hash_table_size;

			if (!hash_size_record[hash_value]) {
				distance = CCD::internal::edgeEdgeDistanceUnclassified(ori_pos_0[pair[i + 1]][indices_0[0]].data(),
					ori_pos_0[pair[i + 1]][indices_0[1]].data(),
					ori_pos_1[pair[i + 3]][indices_1[0]].data(),					
					ori_pos_1[pair[i + 3]][indices_1[1]].data());

				hash_size_record[hash_value] = 2;
				hash_record[hash_value * pair_hash_table_cell_size] = actual_ele_0;
				hash_record[hash_value * pair_hash_table_cell_size + 1] = actual_ele_1;

				record_index->emplace_back(pair[i + 1]);
				record_index->emplace_back(pair[i]);
				record_index->emplace_back(pair[i + 3]);
				record_index->emplace_back(pair[i + 2]);
				record_d_hat->emplace_back((std::max)(distance, d_hat_2));
			}
			else {
				address = hash_record + hash_value * pair_hash_table_cell_size;
				for (unsigned int i = 0; i < hash_size_record[hash_value]; i += 2) {
					if (address[i] == actual_ele_0 && address[i + 1] == actual_ele_1) {
						goto not_add_value;
					}
				}

				distance = CCD::internal::edgeEdgeDistanceUnclassified(ori_pos_0[pair[i + 1]][indices_0[0]].data(),
					ori_pos_0[pair[i + 1]][indices_0[1]].data(),
					ori_pos_1[pair[i + 3]][indices_1[0]].data(),
					ori_pos_1[pair[i + 3]][indices_1[1]].data());

				address[hash_size_record[hash_value]] = actual_ele_0;
				address[hash_size_record[hash_value] + 1] = actual_ele_1;
				hash_size_record[hash_value] += 2;

				record_index->emplace_back(pair[i + 1]);
				record_index->emplace_back(pair[i]);
				record_index->emplace_back(pair[i + 3]);
				record_index->emplace_back(pair[i + 2]);
				record_d_hat->emplace_back((std::max)(distance, d_hat_2));

			}
		not_add_value:;

			if (time < collision_time) {
				collision_time = time;

			}
		}
	}
}


void Collision::vertexTriangleCollisionTimePair(int start_pair_index,
	int end_pair_index, double& collision_time, std::array<int,3>** triangle_indices,
	std::array<double, 3>** vertex_for_render_0, std::array<double, 3>** vertex_for_render_1,
	std::array<double, 3>** vertex_position_0,
	std::array<double, 3>** vertex_position_1, unsigned int* pair, std::vector<unsigned int>* record_index,
	std::vector<double>* record_d_hat, unsigned int* hash_size_record, unsigned int* hash_record,
	unsigned int* vertex_index_prefix_sum, unsigned int* triangle_index_prefix_sum,
	std::array<double, 3>** ori_pos_0, std::array<double, 3>**ori_pos_1)
{
	//unsigned int* pair = spatial_hashing.vertex_triangle_pair[pair_thread_No] + 1;
	double time;
	int* indices;
	double distance = d_hat_2;
	unsigned int hash_value;
	unsigned int actual_ele_0, actual_ele_1;

	unsigned int* address;

		//need to check if the pair already exist in the list
	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		indices = triangle_indices[pair[i + 3]][pair[i + 2]].data();
		time = CCD::pointTriangleCcd(vertex_for_render_0[pair[i + 1]][pair[i]].data(),
			vertex_for_render_1[pair[i + 3]][indices[0]].data(),
			vertex_for_render_1[pair[i + 3]][indices[1]].data(),
			vertex_for_render_1[pair[i + 3]][indices[2]].data(),
			vertex_position_0[pair[i + 1]][pair[i]].data(),
			vertex_position_1[pair[i + 3]][indices[0]].data(),
			vertex_position_1[pair[i + 3]][indices[1]].data(),
			vertex_position_1[pair[i + 3]][indices[2]].data(), eta, tolerance);
		if (time < 1.0) {

			actual_ele_0 = vertex_index_prefix_sum[pair[i + 1]] + pair[i];
			actual_ele_1 = triangle_index_prefix_sum[pair[i + 3]] + pair[i + 2];

			hash_value = ((actual_ele_0 * P1) ^ (P2 * actual_ele_1)) % pair_hash_table_size;

			if (!hash_size_record[hash_value]) {
				distance = CCD::internal::pointTriangleDistanceUnclassified(ori_pos_0[pair[i + 1]][pair[i]].data(), 
					ori_pos_1[pair[i + 3]][indices[0]].data(), ori_pos_1[pair[i + 3]][indices[1]].data(), ori_pos_1[pair[i + 3]][indices[2]].data());

				hash_size_record[hash_value] = 2;
				hash_record[hash_value * pair_hash_table_cell_size] = actual_ele_0;
				hash_record[hash_value * pair_hash_table_cell_size +1 ] = actual_ele_1;

				record_index->emplace_back(pair[i + 1]);
				record_index->emplace_back(pair[i]);
				record_index->emplace_back(pair[i + 3]);
				record_index->emplace_back(pair[i + 2]);
				record_d_hat->emplace_back((std::max)(distance, d_hat_2));
			}
			else {
				address = hash_record + hash_value * pair_hash_table_cell_size;
				for (unsigned int i = 0; i < hash_size_record[hash_value]; i += 2) {
					if (address[i] == actual_ele_0 && address[i + 1] == actual_ele_1) {							
						goto not_add_value;
					}
				}

				distance = CCD::internal::pointTriangleDistanceUnclassified(ori_pos_0[pair[i + 1]][pair[i]].data(),
					ori_pos_1[pair[i + 3]][indices[0]].data(), ori_pos_1[pair[i + 3]][indices[1]].data(), ori_pos_1[pair[i + 3]][indices[2]].data());

				address[hash_size_record[hash_value]] = actual_ele_0;
				address[hash_size_record[hash_value] + 1] = actual_ele_1;
				hash_size_record[hash_value] += 2;

				record_index->emplace_back(pair[i + 1]);
				record_index->emplace_back(pair[i]);
				record_index->emplace_back(pair[i + 3]);
				record_index->emplace_back(pair[i + 2]);
				record_d_hat->emplace_back((std::max)(distance, d_hat_2));
			}
		not_add_value:;
			if (time < collision_time) {
				collision_time = time;
			}
		}
	}
	
}


void Collision::vertexTriangleCollisionTime(int thread_No, unsigned int pair_thread_No, int start_pair_index, 
	int end_pair_index, double& collision_time)
{
	unsigned int* pair = spatial_hashing.vertex_triangle_pair[pair_thread_No] + 1;
	double time;
	int* indices;

	//double* record_collision_time = record_VT_collision_time[pair_thread_No].data();

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		indices = triangle_indices[pair[i + 3]][pair[i + 2]].data();
		//if (approx_CCD.vertexTriangleCCD(time, vertex_for_render[pair[i + 1]][pair[i]].data(), vertex_position[pair[i + 1]][pair[i]].data(),
		//	vertex_for_render[pair[i + 3]][indices[0]].data(), vertex_position[pair[i + 3]][indices[0]].data(), vertex_for_render[pair[i + 3]][indices[1]].data(),
		//	vertex_position[pair[i + 3]][indices[1]].data(), vertex_for_render[pair[i + 3]][indices[2]].data(), vertex_position[pair[i + 3]][indices[2]].data(), conservative_rescaling_)) {
		//	if (time < collision_time) {
		//		collision_time = time;
		//	}
		//}
		//if (approx_CCD.pointTriangleCollisionTime(time, vertex_for_render[pair[i + 1]][pair[i]].data(), vertex_position[pair[i + 1]][pair[i]].data(),
		//	vertex_for_render[pair[i + 3]][indices[0]].data(), vertex_position[pair[i + 3]][indices[0]].data(), vertex_for_render[pair[i + 3]][indices[1]].data(),
		//	vertex_position[pair[i + 3]][indices[1]].data(), vertex_for_render[pair[i + 3]][indices[2]].data(), vertex_position[pair[i + 3]][indices[2]].data(),
		//	triangle_normal_render_not_normalized[pair[i + 3]][pair[i + 2]].data(), triangle_normal_not_normalized[pair[i + 3]][pair[i + 2]].data(),
		//	cross_for_approx_CCD[pair[i + 3]][pair[i + 2]].data(), tolerance_)) {
		//	if (time < collision_time) {
		//		collision_time = time;
		//	}
		//}
		time = CCD::pointTriangleCcd(vertex_for_render[pair[i + 1]][pair[i]].data(),
			vertex_for_render[pair[i + 3]][indices[0]].data(),
			vertex_for_render[pair[i + 3]][indices[1]].data(),
			vertex_for_render[pair[i + 3]][indices[2]].data(),
			vertex_position[pair[i + 1]][pair[i]].data(),
			vertex_position[pair[i + 3]][indices[0]].data(),
			vertex_position[pair[i + 3]][indices[1]].data(),
			vertex_position[pair[i + 3]][indices[2]].data(), eta, tolerance);

		//if (*(pair + i + 3) == 0 && *(pair + i + 2) == 37 && *(pair + i) == 43 && *(pair + 1 + i) == 1) {
		//	double test_time = 1.0;
		//	if (approx_CCD.vertexTriangleCCD(test_time, vertex_for_render[pair[i + 1]][pair[i]].data(), vertex_position[pair[i + 1]][pair[i]].data(),
		//		vertex_for_render[pair[i + 3]][indices[0]].data(), vertex_position[pair[i + 3]][indices[0]].data(), vertex_for_render[pair[i + 3]][indices[1]].data(),
		//		vertex_position[pair[i + 3]][indices[1]].data(), vertex_for_render[pair[i + 3]][indices[2]].data(), vertex_position[pair[i + 3]][indices[2]].data(), conservative_rescaling_)) {
		//		std::cout << "check collision " << test_time << std::endl;
		//	}
		//	else{
		//		std::cout<< "check collision not collide " << test_time << std::endl;
		//	}
		//	std::cout << "collision_time-- " << time << std::endl;
		//}
		//	double test_time=1.0;
		//	if (approx_CCD.vertexTriangleCCD(test_time, vertex_for_render[pair[i + 1]][pair[i]].data(), vertex_position[pair[i + 1]][pair[i]].data(),
		//		vertex_for_render[pair[i + 3]][indices[0]].data(), vertex_position[pair[i + 3]][indices[0]].data(), vertex_for_render[pair[i + 3]][indices[1]].data(),
		//		vertex_position[pair[i + 3]][indices[1]].data(), vertex_for_render[pair[i + 3]][indices[2]].data(), vertex_position[pair[i + 3]][indices[2]].data(), conservative_rescaling_)) {
		//		//std::cout << "check collision " << test_time << std::endl;
		//	}
		//	else{
		//		std::cout<< "check collision not collide " << test_time << std::endl;
		//	}
		//}
		//if (pair[i + 1] == 0 && pair[i] == 1 && pair[i+2]==0 && pair[i+3]==1) {
		//	//std::cout << "show collision time "<< time << std::endl;
		//	double triangle_normal[3];
		//	getTriangleNormal(vertex_for_render[pair[i + 3]][indices[0]].data(),
		//		vertex_for_render[pair[i + 3]][indices[1]].data(),
		//		vertex_for_render[pair[i + 3]][indices[2]].data(), triangle_normal);
		//	double barycentric[3];
		//	double d_2 = CCD::internal::pointTriangleNearestDistance(vertex_for_render[pair[i + 1]][pair[i]].data(),
		//		vertex_for_render[pair[i + 3]][indices[0]].data(),
		//		vertex_for_render[pair[i + 3]][indices[1]].data(),
		//		vertex_for_render[pair[i + 3]][indices[2]].data(), triangle_normal, barycentric, d_hat_2);
		//	//std::cout << "distance render " << d_2<<" "<<d_hat_2 << std::endl;
		//	getTriangleNormal(vertex_position[pair[i + 3]][indices[0]].data(),
		//		vertex_position[pair[i + 3]][indices[1]].data(),
		//		vertex_position[pair[i + 3]][indices[2]].data(), triangle_normal);
		//	d_2= CCD::internal::pointTriangleNearestDistance(vertex_position[pair[i + 1]][pair[i]].data(),
		//		vertex_position[pair[i + 3]][indices[0]].data(),
		//		vertex_position[pair[i + 3]][indices[1]].data(),
		//		vertex_position[pair[i + 3]][indices[2]].data(), triangle_normal, barycentric, d_hat_2);
		//	//std::cout << "distance current " << d_2 <<" "<<d_hat_2<< std::endl;
		//	//std::cout << "position " << vertex_for_render[pair[i + 1]][pair[i]][0] << " " << vertex_for_render[pair[i + 1]][pair[i]][1] << " " << vertex_for_render[pair[i + 1]][pair[i]][2] << std::endl;
		//	//std::cout << vertex_for_render[pair[i + 3]][indices[0]][0] << " " << vertex_for_render[pair[i + 3]][indices[0]][1] << " " << vertex_for_render[pair[i + 3]][indices[0]][2] << std::endl;
		//	//std::cout << vertex_for_render[pair[i + 3]][indices[1]][0] << " " << vertex_for_render[pair[i + 3]][indices[1]][1] << " " << vertex_for_render[pair[i + 3]][indices[1]][2] << std::endl;
		//	//std::cout << vertex_for_render[pair[i + 3]][indices[2]][0] << " " << vertex_for_render[pair[i + 3]][indices[2]][1] << " " << vertex_for_render[pair[i + 3]][indices[2]][2] << std::endl;
		//	//std::cout << "== " << std::endl;
		//	//std::cout << vertex_position[pair[i + 1]][pair[i]][0] << " " << vertex_position[pair[i + 1]][pair[i]][1] << " " << vertex_position[pair[i + 1]][pair[i]][2] << std::endl;
		//	//std::cout << vertex_position[pair[i + 3]][indices[0]][0] << " " << vertex_position[pair[i + 3]][indices[0]][1] << " " << vertex_position[pair[i + 3]][indices[0]][2] << std::endl;
		//	//std::cout << vertex_position[pair[i + 3]][indices[1]][0] << " " << vertex_position[pair[i + 3]][indices[1]][1] << " " << vertex_position[pair[i + 3]][indices[1]][2] << std::endl;
		//	//std::cout << vertex_position[pair[i + 3]][indices[2]][0] << " " << vertex_position[pair[i + 3]][indices[2]][1] << " " << vertex_position[pair[i + 3]][indices[2]][2] << std::endl;
		//}
		//record_collision_time[i >> 2] = time;
		//////std::cout << "time "<< pair[i + 1]<<" "<< pair[i] << " " << pair[i + 3]<<" "<< pair[i + 2]<<" " << time << std::endl;
		//////std::cout << "index info " << pair[i + 3] << " " << indices[0] << " " << indices[1] << " " << indices[2] << std::endl;
		//////std::cout << vertex_for_render[pair[i + 1]][pair[i]][1] << " " << vertex_position[pair[i + 1]][pair[i]][1] <<
		//	" " << vertex_for_render[pair[i + 3]][indices[0]][1] << " " << vertex_position[pair[i + 3]][indices[0]][1] << std::endl;

		if (time < collision_time) {
			collision_time = time;
		}		
	}
}

void Collision::edgeEdgeCollisionTimeCompare(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, double& collision_time)
{
	unsigned int* pair = spatial_hashing.edge_edge_pair[pair_thread_No] + 1;
	double conservative_rescaling_ = conservative_rescaling;
	double time;
	unsigned int* indices_1;
	unsigned int* indices_0;

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		indices_0 = edge_vertices[pair[i + 1]] + (pair[i] << 1);
		indices_1 = edge_vertices[pair[i + 3]] + (pair[i + 2] << 1);
		if (approx_CCD.edgeEdgeCCD(time, vertex_position_eigen[pair[i + 1]][indices_0[0]],
			vertex_position_eigen[pair[i + 1]][indices_0[1]], vertex_for_render_eigen[pair[i + 1]][indices_0[0]],
			vertex_for_render_eigen[pair[i + 1]][indices_0[1]], vertex_position_eigen[pair[i + 3]][indices_1[0]],
			vertex_position_eigen[pair[i + 3]][indices_1[1]], vertex_for_render_eigen[pair[i + 3]][indices_1[0]],
			vertex_for_render_eigen[pair[i + 3]][indices_1[1]], conservative_rescaling_)) {
			if (time < collision_time) {
				collision_time = time;
			}
		}
	}
}





void Collision::edgeEdgeCollisionTime(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, double& collision_time)
{
	unsigned int* pair = spatial_hashing.edge_edge_pair[pair_thread_No] + 1;
	double conservative_rescaling_ = conservative_rescaling;
	double time;
	unsigned int* indices_1;
	unsigned int* indices_0;

	for (int i = start_pair_index; i < end_pair_index; i += 4) {
		indices_0 = edge_vertices[pair[i + 1]] + (pair[i]<<1);
		indices_1 = edge_vertices[pair[i + 3]] + (pair[i + 2]<<1);
		//if (approx_CCD.edgeEdgeCCD(time, vertex_position[pair[i + 1]][indices_0[0]].data(),
		//	vertex_position[pair[i + 1]][indices_0[1]].data(), vertex_for_render[pair[i + 1]][indices_0[0]].data(),
		//	vertex_for_render[pair[i + 1]][indices_0[1]].data(), vertex_position[pair[i + 3]][indices_1[0]].data(),
		//	vertex_position[pair[i + 3]][indices_1[1]].data(), vertex_for_render[pair[i + 3]][indices_1[0]].data(),
		//	vertex_for_render[pair[i + 3]][indices_1[1]].data(), conservative_rescaling_)) {
		//	if (time < collision_time) {
		//		collision_time = time;
		//	}
		//}


		time = CCD::edgeEdgeCcd(vertex_for_render[pair[i + 1]][indices_0[0]].data(),
			vertex_for_render[pair[i + 1]][indices_0[1]].data(),
			vertex_for_render[pair[i + 3]][indices_1[0]].data(),
			vertex_for_render[pair[i + 3]][indices_1[1]].data(),
			vertex_position[pair[i + 1]][indices_0[0]].data(),
			vertex_position[pair[i + 1]][indices_0[1]].data(),
			vertex_position[pair[i + 3]][indices_1[0]].data(),
			vertex_position[pair[i + 3]][indices_1[1]].data(), eta, tolerance);


		//if (pair[i + 1] == 0 && pair[i] == 0 && pair[i + 2] == 3 && pair[i + 3] == 1) {

		//	double test_time=2.0;
		//	if (approx_CCD.edgeEdgeCCD(test_time, vertex_position[pair[i + 1]][indices_0[0]].data(),
		//		vertex_position[pair[i + 1]][indices_0[1]].data(), vertex_for_render[pair[i + 1]][indices_0[0]].data(),
		//		vertex_for_render[pair[i + 1]][indices_0[1]].data(), vertex_position[pair[i + 3]][indices_1[0]].data(),
		//		vertex_position[pair[i + 3]][indices_1[1]].data(), vertex_for_render[pair[i + 3]][indices_1[0]].data(),
		//		vertex_for_render[pair[i + 3]][indices_1[1]].data(), conservative_rescaling_)) {
		//	}

		//	//std::cout << "print this pair " << time<<" "<< test_time << std::endl;
		//}

		////std::cout << "ee collision " << time << std::endl;

		if (time < collision_time) {
			collision_time = time;

		}
	}
}

void Collision::edgeEdgeColliderCollisionTime(int thread_No, unsigned int pair_thread_No, int start_pair_index,
	int end_pair_index, double& collision_time)
{
	//unsigned int* pair = spatial_hashing.edge_edge_pair_collider[pair_thread_No] + 1;
	//double conservative_rescaling_ = conservative_rescaling;
	//double time;
	//unsigned int* indices_1;
	//unsigned int* indices_0;
	//for (int i = start_pair_index; i < end_pair_index; i += 4) {
	//	indices_0 = edge_vertices[pair[i + 1]] + (pair[i] << 1);
	//	indices_1 = collider_edge_vertices[pair[i + 3]] + (pair[i + 2] << 1);
	//	if (approx_CCD.edgeEdgeCCD(time, vertex_position[pair[i + 1]][indices_0[0]].data(),
	//		vertex_position[pair[i + 1]][indices_0[1]].data(), vertex_for_render[pair[i + 1]][indices_0[0]].data(),
	//		vertex_for_render[pair[i + 1]][indices_0[1]].data(), vertex_position_collider[pair[i + 3]][indices_1[0]].data(),
	//		vertex_position_collider[pair[i + 3]][indices_1[1]].data(), vertex_for_render_collider[pair[i + 3]][indices_1[0]].data(),
	//		vertex_for_render_collider[pair[i + 3]][indices_1[1]].data(), conservative_rescaling_)) {
	//		if (time < collision_time) {
	//			collision_time = time;
	//		}
	//	}
	//	//time = CCD::edgeEdgeCcd(vertex_for_render[pair[i + 1]][indices_0[0]].data(),
	//	//	vertex_for_render[pair[i + 1]][indices_0[1]].data(),
	//	//	vertex_for_render_collider[pair[i + 3]][indices_1[0]].data(),
	//	//	vertex_for_render_collider[pair[i + 3]][indices_1[1]].data(),
	//	//	vertex_position[pair[i + 1]][indices_0[0]].data(),
	//	//	vertex_position[pair[i + 1]][indices_0[1]].data(),
	//	//	vertex_position_collider[pair[i + 3]][indices_1[0]].data(),
	//	//	vertex_position_collider[pair[i + 3]][indices_1[1]].data(), eta, tolerance);
	//	//if (time < collision_time) {
	//	//	collision_time = time;
	//	//}
	//}
}

//FIND_VERTEX_VERTEX_VERTEX_EDGE_PAIRS
void Collision::findAllVertexVertexEdgePairs(int thread_No)
{
	findAllVertexVertexPairs(thread_No);
	findAllVertexEdgePairs(thread_No);
}

void Collision::findAllVertexVertexPairs(int thread_No)
{
	vertex_vertex_pair[thread_No][0] = 0;
	vertex_vertex_pair_collider[thread_No][0] = 0;

	if (vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1] >
		vertex_triangle_pair_index_start_per_thread[thread_No << 1]) {
		findVertexVertexPair(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]][0], VERTEX_VERTEX);
		for (int i = vertex_triangle_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			findVertexVertexPair(thread_No, i, 0,
				spatial_hashing.vertex_triangle_pair[i][0], VERTEX_VERTEX);
		}
		findVertexVertexPair(thread_No, vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], VERTEX_VERTEX);
	}
	else {
		findVertexVertexPair(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], VERTEX_VERTEX);
	}

	if (has_collider) {
		if (vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]) {
			findVertexVertexPair(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
				spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]][0], VERTEX_VERTEX_COLLIDER);
			for (int i = vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1] + 1;
				i < vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
				findVertexVertexPair(thread_No, i, 0,
					spatial_hashing.vertex_obj_triangle_collider_pair[i][0], VERTEX_VERTEX_COLLIDER);
			}
			findVertexVertexPair(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
				vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3], VERTEX_VERTEX_COLLIDER);
		}
		else {
			findVertexVertexPair(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
				vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3], VERTEX_VERTEX_COLLIDER);
		}
	}

}


void Collision::findAllVertexEdgePairs(int thread_No)
{
	vertex_edge_pair[thread_No][0] = 0;
	vertex_obj_edge_collider_pair[thread_No][0] = 0;
	vertex_collider_edge_obj_pair[thread_No][0] = 0;

	if (vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1] >
		vertex_triangle_pair_index_start_per_thread[thread_No << 1]) {
		findVertexEdgePair(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]][0], VERTEX_EDGE);
		for (int i = vertex_triangle_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			findVertexEdgePair(thread_No, i, 0,
				spatial_hashing.vertex_triangle_pair[i][0], VERTEX_EDGE);
		}
		findVertexEdgePair(thread_No, vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], VERTEX_EDGE);
	}
	else {
		findVertexEdgePair(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], VERTEX_EDGE);
	}


	if (has_collider) {
		//vertex_collider_triangle
		//if (vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1]) {
		//	findVertexEdgePair(thread_No, vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1], vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 1],
		//		spatial_hashing.vertex_collider_triangle_obj_pair[vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1]][0], VERTEX_COLLIDER_EDGE_OBJ);
		//	for (int i = vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1] + 1;
		//		i < vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
		//		findVertexEdgePair(thread_No, i, 0,
		//			spatial_hashing.vertex_collider_triangle_obj_pair[i][0], VERTEX_COLLIDER_EDGE_OBJ);
		//	}
		//	findVertexEdgePair(thread_No, vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
		//		vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 3], VERTEX_COLLIDER_EDGE_OBJ);
		//}
		//else {
		//	findVertexEdgePair(thread_No, vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1], vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 1],
		//		vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 3], VERTEX_COLLIDER_EDGE_OBJ);
		//}

		if (vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]) {
			findVertexEdgePair(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
				spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]][0], VERTEX_OBJ_EDGE_COLLIDER);
			for (int i = vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1] + 1;
				i < vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
				findVertexEdgePair(thread_No, i, 0,
					spatial_hashing.vertex_obj_triangle_collider_pair[i][0], VERTEX_OBJ_EDGE_COLLIDER);
			}
			findVertexEdgePair(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
				vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3], VERTEX_OBJ_EDGE_COLLIDER);
		}
		else {
			findVertexEdgePair(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
				vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3], VERTEX_OBJ_EDGE_COLLIDER);
		}	
	}
}


//
void Collision::collisionTimeCompare(int thread_No)
{
	double collision_time = 2.0;
	//vertex_triangle
	if (vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_triangle_pair_index_start_per_thread[thread_No << 1]) {
		vertexTriangleCollisionTimeCompare(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]][0], collision_time);
		for (int i = vertex_triangle_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			vertexTriangleCollisionTimeCompare(thread_No, i, 0,
				spatial_hashing.vertex_triangle_pair[i][0], collision_time);
		}
		vertexTriangleCollisionTimeCompare(thread_No, vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
	}
	else {
		vertexTriangleCollisionTimeCompare(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
	}

	//edge edge
	if (edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1] > edge_edge_pair_index_start_per_thread[thread_No << 1]) {
		edgeEdgeCollisionTimeCompare(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.edge_edge_pair[edge_edge_pair_index_start_per_thread[thread_No << 1]][0], collision_time);
		for (int i = edge_edge_pair_index_start_per_thread[thread_No << 1] + 1;
			i < edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			edgeEdgeCollisionTimeCompare(thread_No, i, 0,
				spatial_hashing.edge_edge_pair[i][0], collision_time);
		}
		edgeEdgeCollisionTimeCompare(thread_No, edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
	}
	else {
		edgeEdgeCollisionTimeCompare(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
	}

	// vertex edge
	if (vertex_edge_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_edge_pair_index_start_per_thread[thread_No << 1]) {
		vertexEdgeCollisionTimeCompare(thread_No, vertex_edge_pair_index_start_per_thread[thread_No << 1], vertex_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_edge_pair[vertex_edge_pair_index_start_per_thread[thread_No << 1]][0], collision_time, VERTEX_EDGE);
		for (int i = vertex_edge_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_edge_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			vertexEdgeCollisionTimeCompare(thread_No, i, 0,
				vertex_edge_pair[i][0], collision_time, VERTEX_EDGE);
		}
		vertexEdgeCollisionTimeCompare(thread_No, vertex_edge_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_edge_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_EDGE);
	}
	else {
		vertexEdgeCollisionTimeCompare(thread_No, vertex_edge_pair_index_start_per_thread[thread_No << 1], vertex_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_edge_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_EDGE);
	}
	////vertex vertex
	if (vertex_vertex_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_vertex_pair_index_start_per_thread[thread_No << 1]) {
		vertexVertexCollisionTimeCompare(thread_No, vertex_vertex_pair_index_start_per_thread[thread_No << 1], vertex_vertex_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_vertex_pair[vertex_vertex_pair_index_start_per_thread[thread_No << 1]][0], collision_time, VERTEX_VERTEX);
		for (int i = vertex_vertex_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_vertex_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			vertexVertexCollisionTimeCompare(thread_No, i, 0,
				vertex_vertex_pair[i][0], collision_time, VERTEX_VERTEX);
		}
		vertexVertexCollisionTimeCompare(thread_No, vertex_vertex_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_vertex_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_VERTEX);
	}
	else {
		vertexVertexCollisionTimeCompare(thread_No, vertex_vertex_pair_index_start_per_thread[thread_No << 1], vertex_vertex_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_vertex_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_VERTEX);
	}
}


//void Collision::collisionTimeSingleVertex(unsigned int obj_index, unsigned int vertex_index, unsigned int vertex_index_on_surface,
//	std::array<double, 3>* initial_pos, std::array<double, 3>* current_pos)
//{
//	double collision_time = 2.0;
//	VTCollisionTimeOneVertex(initial_pos[vertex_index].data(), current_pos[vertex_index].data(), collision_time, 
//		spatial_hashing.vertex_triangle_pair_num_record[obj_index][vertex_index_on_surface],
//		spatial_hashing.vertex_triangle_pair_by_vertex[obj_index] + estimate_coeff_for_vt_pair_num * vertex_index_on_surface,
//		vertex_collision_free.data(), vertex_position.data(),triangle_indices.data());
//	
//	std::vector<unsigned int>* element = &mesh_struct[obj_index]->vertices[vertex_index].edge;
//	unsigned int* element_;
//	for (auto i = element->begin(); i != element->end(); ++i) {
//		element_ = edge_vertices[obj_index] + ((*i) << 1);
//		EECollisionTimeOneEdge(initial_pos[*element_].data(), initial_pos[*(element_ + 1)].data(),
//			current_pos[*element_].data(), current_pos[*(element_ + 1)].data(),
//			collision_time, *i, obj_index, spatial_hashing.edge_edge_pair_num_record[obj_index][*i],
//			spatial_hashing.edge_edge_pair_by_edge[obj_index] + estimate_coeff_for_ee_pair_num * (*i), vertex_collision_free.data(),
//			vertex_position.data(), true, edge_vertices.data());
//	}
//	//element = &mesh_struct[obj_index]->vertices[vertex_index].face;
//	//int* triangle_;
//	//for (auto i = element->begin(); i != element->end(); ++i) {
//	//	triangle_ = triangle_indices[obj_index][*i].data();
//	//	TVCollisionTimeOneVertex(initial_pos[*triangle_].data(), initial_pos[*(triangle_ + 1)].data(),
//	//		initial_pos[*(triangle_ + 2)].data(), current_pos[*triangle_].data(), current_pos[*(triangle_ + 1)].data(),
//	//		current_pos[*(triangle_ + 2)].data(), collision_time,
//	//		spatial_hashing.triangle_vertex_pair_num_record[obj_index][*i],
//	//		spatial_hashing.triangle_vertex_pair_by_triangle[obj_index] + estimate_coeff_for_tv_pair_num * (*i),
//	//		vertex_collision_free.data(), vertex_position.data());
//	//}
//
//
//}


// lack of tv collider
void Collision::collisionFreeOneVertex(unsigned int obj_No, unsigned int vertex_No, unsigned int vertex_index_on_surface, double* initial_vertex_pos, double* current_vertex_pos,
	std::array<double, 3>* initial_pos_this_obj, std::array<double, 3>* current_pos_this_obj,
	std::array<double, 3>**current_pos)
{
	double collision_time = 1.0;
	std::vector<unsigned int>* element;

	//VT
	VTCollisionTimeOneVertex(initial_vertex_pos, current_vertex_pos, collision_time,
		vertex_triangle_pair_num_record[obj_No][vertex_index_on_surface],
		vertex_triangle_pair_by_vertex[obj_No] + close_vt_pair_num * vertex_index_on_surface, current_pos,
		current_pos,triangle_indices.data(),false,vertex_No);
	
	//TV
	element = &mesh_struct[obj_No]->vertices[vertex_No].face;
	int* triangle_;
	for (auto i = element->begin(); i != element->end(); ++i) {
		triangle_ = triangle_indices[obj_No][*i].data();
		TVCollisionTimeOneTriangle(initial_pos_this_obj[*triangle_].data(), initial_pos_this_obj[*(triangle_ + 1)].data(),
			initial_pos_this_obj[*(triangle_ + 2)].data(), current_pos_this_obj[*triangle_].data(), 
			current_pos_this_obj[*(triangle_ + 1)].data(),
			current_pos_this_obj[*(triangle_ + 2)].data(), collision_time,
			triangle_vertex_pair_num_record[obj_No][*i],
			triangle_vertex_pair_by_triangle[obj_No] + close_tv_pair_num * (*i),
			current_pos, current_pos,3);
	}
	element = &mesh_struct[obj_No]->vertices[vertex_No].edge;
	unsigned int* element_;
	for (auto i = element->begin(); i != element->end(); ++i) {
		element_ = edge_vertices[obj_No] + ((*i) << 1);
		EECollisionTimeOneEdge(initial_pos_this_obj[*element_].data(), initial_pos_this_obj[*(element_ + 1)].data(),
			current_pos_this_obj[*element_].data(), current_pos_this_obj[*(element_ + 1)].data(),
			collision_time, *i, obj_No, edge_edge_pair_number_record[obj_No][*i],
			edge_edge_pair_by_edge[obj_No] + close_ee_pair_num * (*i), current_pos, current_pos,true,edge_vertices.data());
	}
	//VT collider
	if (has_collider) {
		VTCollisionTimeOneVertex(initial_vertex_pos, current_vertex_pos, collision_time, 
			vertex_obj_triangle_collider_num_record[obj_No][vertex_index_on_surface],
			vertex_obj_triangle_collider_pair_by_vertex[obj_No] + close_vt_collider_pair_num * vertex_index_on_surface,
			vertex_position_collider.data(), vertex_position_collider.data(),triangle_indices_collider.data(),true, vertex_No);
	}
	//floor
	if (floor->exist) {	
		floorCollisionTime(initial_vertex_pos, current_vertex_pos, floor->dimension,
			floor->normal_direction, floor->value, collision_time, tolerance);		
	}

	if (collision_time < 1.0) {
		collision_time *= 0.9;
		current_vertex_pos[0] = initial_vertex_pos[0] + collision_time * (current_vertex_pos[0] - initial_vertex_pos[0]);
		current_vertex_pos[1] = initial_vertex_pos[1] + collision_time * (current_vertex_pos[1] - initial_vertex_pos[1]);
		current_vertex_pos[2] = initial_vertex_pos[2] + collision_time * (current_vertex_pos[2] - initial_vertex_pos[2]);
	}
}


//COLOR_COLLISION_TIME
void Collision::colorCollisionTime(int thread_No, int color_No)
{
	int i;
	unsigned int end_per_thread;
	unsigned int* index_of_a_tet_color;
	int* vertex_index_on_surface;

	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* ini_vertex_pos;


	unsigned int* edge_vertices;



	double* edge_length;
	int surface_index;


	unsigned int* edge_vertex;

	double collision_time = 1.0;

	bool has_collider = this->has_collider;

	bool has_floor = floor->exist;

	int color_group_index;
	for (int tet_obj_no = 0; tet_obj_no < tetrahedron->size(); ++tet_obj_no) {
		i = tet_obj_no + cloth->size();
		color_group_index = *inner_iteration_number % tet_color_groups[i]->size();

		//if (color_No >=tet_color_groups[i]->data()[color_group_index].size()-1) {
		//	continue;
		//}
	
		vertex_pos = vertex_position[i];
		ini_vertex_pos = vertex_record_for_this_color[i];

		index_of_a_tet_color = surface_vertex_index_of_a_tet_color_group[i][color_group_index][color_No].data();
		vertex_index_on_surface = this->general_index_to_surface_index[i];
		end_per_thread = surface_vertex_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color_No][thread_No + 1];
		//vt vt_collider floor
		for (int j = surface_vertex_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color_No][thread_No]; j < end_per_thread; ++j) {
			surface_index = vertex_index_on_surface[index_of_a_tet_color[j]];

			VTCollisionTimeOneVertex(ini_vertex_pos[index_of_a_tet_color[j]].data(),
				vertex_pos[index_of_a_tet_color[j]].data(), collision_time,
				vertex_triangle_pair_num_record[i][surface_index],
				vertex_triangle_pair_by_vertex[i] + close_vt_pair_num * surface_index,
				vertex_record_for_this_color.data(), vertex_position.data(), triangle_indices.data(),false, index_of_a_tet_color[j]);

			if (has_collider) {
				VTCollisionTimeOneVertex(ini_vertex_pos[index_of_a_tet_color[j]].data(),
					vertex_pos[index_of_a_tet_color[j]].data(), collision_time,
					vertex_obj_triangle_collider_num_record[i][surface_index],
					vertex_obj_triangle_collider_pair_by_vertex[i] + close_vt_collider_pair_num * surface_index,
					vertex_position_collider.data(), vertex_position_collider.data(), triangle_indices_collider.data(),true, index_of_a_tet_color[j]);
			}
			if (has_floor) {
				floorCollisionTime(ini_vertex_pos[index_of_a_tet_color[j]].data(),
					vertex_pos[index_of_a_tet_color[j]].data(), floor->dimension,
					floor->normal_direction, floor->value, collision_time, tolerance);
			}
			
		}
		//ee ee_collider
		end_per_thread = edge_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color_No][thread_No + 1];
		index_of_a_tet_color = edge_index_of_a_tet_color_group[i][color_group_index][color_No].data();
		edge_vertices = this->edge_vertices[i];
		edge_length = rest_edge_length[i];
		for (int j = edge_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color_No][thread_No]; j < end_per_thread; ++j) {
			edge_vertex = edge_vertices + (index_of_a_tet_color[j] << 1);
			EECollisionTimeOneEdgeAll(ini_vertex_pos[*edge_vertex].data(), ini_vertex_pos[*(edge_vertex + 1)].data(),
				vertex_pos[*edge_vertex].data(), vertex_pos[*(edge_vertex + 1)].data(), collision_time,
				edge_edge_pair_number_record[i][index_of_a_tet_color[j]],
				edge_edge_pair_by_edge[i] + close_ee_pair_num * index_of_a_tet_color[j],
				vertex_record_for_this_color.data(), vertex_position.data(), this->edge_vertices.data(),3);

			if (has_collider) {
				EECollisionTimeOneEdgeAll(ini_vertex_pos[*edge_vertex].data(), ini_vertex_pos[*(edge_vertex + 1)].data(),
					vertex_pos[*edge_vertex].data(), vertex_pos[*(edge_vertex + 1)].data(), collision_time,
					edge_edge_collider_pair_num_record[i][index_of_a_tet_color[j]],
					edge_edge_collider_pair_by_edge[i] + close_ee_collider_pair_num * index_of_a_tet_color[j],
					vertex_position_collider.data(), vertex_position_collider.data(), this->collider_edge_vertices.data(),2);

			}
		}

		//TV TV collider
		end_per_thread = triangle_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color_No][thread_No + 1];
		index_of_a_tet_color = triangle_index_of_a_tet_color_group[i][color_group_index][color_No].data();
		int* triangle_vertex;

		for (int j = triangle_index_of_a_tet_color_per_thread_start_group[i][color_group_index][color_No][thread_No]; j < end_per_thread; ++j) {
			triangle_vertex = triangle_indices[i][index_of_a_tet_color[j]].data();
			TVCollisionTimeOneTriangleSelfColor(ini_vertex_pos[triangle_vertex[0]].data(), ini_vertex_pos[triangle_vertex[1]].data(),
				ini_vertex_pos[triangle_vertex[2]].data(),
				vertex_pos[triangle_vertex[0]].data(), vertex_pos[triangle_vertex[1]].data(),
				vertex_pos[triangle_vertex[2]].data(), collision_time,
				triangle_vertex_pair_num_record[i][index_of_a_tet_color[j]],
				triangle_vertex_pair_by_triangle[i] + close_tv_pair_num * index_of_a_tet_color[j], vertex_record_for_this_color.data(),
				vertex_position.data(), 3);
			if (has_collider) {
				TVCollisionTimeOneTriangle(ini_vertex_pos[triangle_vertex[0]].data(), ini_vertex_pos[triangle_vertex[1]].data(),
					ini_vertex_pos[triangle_vertex[2]].data(),
					vertex_pos[triangle_vertex[0]].data(), vertex_pos[triangle_vertex[1]].data(),
					vertex_pos[triangle_vertex[2]].data(), collision_time,
					triangle_vertex_collider_pair_num_record[i][index_of_a_tet_color[j]],
					triangle_vertex_collider_pair_by_triangle[i] + close_tv_collider_pair_num * index_of_a_tet_color[j],
					vertex_position_collider.data(), vertex_position_collider.data(), 2);
			}
		}
	}
	collision_time_thread[thread_No] = collision_time;
}


//CLOSE_PAIR_COLLISION_TIME
void Collision::collisionTimeAllClosePair(int thread_No)
{
	double collision_time = 1.0;
	//vt:
	auto start = vt_pair_compressed_record.begin() + vt_per_thread_start_index[thread_No];
	auto end = vt_pair_compressed_record.begin()+ vt_per_thread_start_index[thread_No + 1];
	int* indices;
	double time;


	double vt_collision_time = 1.0;

	for (auto i = start; i < end; i+=5) {
		indices = triangle_indices[*(i+2)][*(i+3)].data();
		time = CCD::pointTriangleCcd(vertex_record_for_this_color[*i][*(i+1)].data(),
			vertex_record_for_this_color[*(i + 2)][indices[0]].data(),
			vertex_record_for_this_color[*(i + 2)][indices[1]].data(),
			vertex_record_for_this_color[*(i + 2)][indices[2]].data(),
			vertex_position[*i][*(i + 1)].data(),
			vertex_position[*(i + 2)][indices[0]].data(),
			vertex_position[*(i + 2)][indices[1]].data(),
			vertex_position[*(i + 2)][indices[2]].data(),  eta, tolerance);
		if (time < vt_collision_time) {
			vt_collision_time = time;
		}
	}

	if (vt_collision_time < collision_time) {
		collision_time = vt_collision_time;
	}


	if (vt_collision_time == 0) {
		std::cout << "self vt collision color time zero " << std::endl;
	}

	//ee
	start = ee_pair_compressed_record.begin() + ee_per_thread_start_index[thread_No];
	end = ee_pair_compressed_record.begin() + ee_per_thread_start_index[thread_No + 1];
	unsigned int* edge_0_vertex;
	unsigned int* edge_1_vertex;
	for (auto i = start; i < end; i += 5) {
		edge_0_vertex = edge_vertices[*i] + ((*(i + 1)) << 1);
		edge_1_vertex = edge_vertices[*(i + 2)] + ((*(i + 3)) << 1);
		time = CCD::edgeEdgeCcd(vertex_record_for_this_color[*i][*edge_0_vertex].data(),
			vertex_record_for_this_color[*i][*(edge_0_vertex+1)].data(), vertex_record_for_this_color[*(i + 2)][*edge_1_vertex].data(),
			vertex_record_for_this_color[*(i+2)][*(edge_1_vertex+1)].data(),
			vertex_position[*i][*edge_0_vertex].data(),
			vertex_position[*i][*(edge_0_vertex + 1)].data(), vertex_position[*(i + 2)][*edge_1_vertex].data(),
			vertex_position[*(i + 2)][*(edge_1_vertex + 1)].data(), eta, tolerance);
		if (time < collision_time) {
			collision_time = time;
		}
	}


	if (has_collider) {

		std::array <double, 3>* initial_pos;
		std::array <double, 3>* current_pos;

		unsigned int* element;
		unsigned int* element_num;

		int* normal_to_surface;


		double vt_collider_collision_time = 1.0;
		double tv_collider_collision_time = 1.0;

		for (int i = 0; i < total_obj_num; ++i) {
			initial_pos = vertex_record_for_this_color[i];
			current_pos = vertex_position[i];

			// vt collider
			element = vertex_obj_triangle_collider_pair_by_vertex[i];
			element_num = vertex_obj_triangle_collider_num_record[i];

			start = vertex_index_collide_with_collider[i].begin() + vertex_index_collide_with_collider_start_per_thread[i * (thread_num + 1) + thread_No];
			end = vertex_index_collide_with_collider[i].begin() + vertex_index_collide_with_collider_start_per_thread[i * (thread_num + 1) + thread_No + 1];
			if (i < cloth->size()) {
				for (auto j = start; j < end; ++j) {
					VTCollisionTimeOneVertex(initial_pos[*j].data(), current_pos[*j].data(), vt_collider_collision_time, element_num[*j],
						element + close_vt_collider_pair_num * (*j), vertex_for_render_collider.data(),
						vertex_position_collider.data(), triangle_indices_collider.data(),true,*j);
				}
			}
			else {
				unsigned int j;
				normal_to_surface = general_index_to_surface_index[i];
				for (auto j = start; j < end; ++j) {
					VTCollisionTimeOneVertex(initial_pos[*j].data(), current_pos[*j].data(), vt_collider_collision_time, element_num[normal_to_surface[*j]],
						element + close_vt_collider_pair_num * normal_to_surface[*j], vertex_for_render_collider.data(),
						vertex_position_collider.data(), triangle_indices_collider.data(),false,*j);
				}
			}

			// ee collider
			element = edge_edge_collider_pair_by_edge[i];
			element_num = edge_edge_collider_pair_num_record[i];

			start = edge_index_collide_with_collider[i].begin() + edge_index_collide_with_collider_start_per_thread[i * (thread_num + 1) + thread_No];
			end = edge_index_collide_with_collider[i].begin() + edge_index_collide_with_collider_start_per_thread[i * (thread_num + 1) + thread_No + 1];

			for (auto j = start; j < end; ++j) {
				edge_0_vertex = edge_vertices[i] + ((*j) << 1);
				EECollisionTimeOneEdge(initial_pos[*edge_0_vertex].data(), initial_pos[*(edge_0_vertex + 1)].data(),
					current_pos[*edge_0_vertex].data(), current_pos[*(edge_0_vertex + 1)].data(),
					collision_time, *j, i, element_num[*j],
					element + close_ee_collider_pair_num * (*j), vertex_for_render_collider.data(), vertex_position_collider.data(),
					false, collider_edge_vertices.data());				

			}

			// tv collider
			element =triangle_vertex_collider_pair_by_triangle[i];
			element_num = triangle_vertex_collider_pair_num_record[i];
			start = triangle_index_collide_with_collider[i].begin() + triangle_index_collide_with_collider_start_per_thread[i * (thread_num + 1) + thread_No];
			end = triangle_index_collide_with_collider[i].begin() + triangle_index_collide_with_collider_start_per_thread[i * (thread_num + 1) + thread_No + 1];

			for (auto j = start; j < end; ++j) {
				indices = triangle_indices[i][*j].data();
				TVCollisionTimeOneTriangle(initial_pos[indices[0]].data(), initial_pos[indices[1]].data(), initial_pos[indices[2]].data(),
					current_pos[indices[0]].data(), current_pos[indices[1]].data(), current_pos[indices[2]].data(), tv_collider_collision_time,
					element_num[*j], element + close_tv_collider_pair_num * (*j), vertex_for_render_collider.data(), vertex_position_collider.data(), 2);
			}
		}
		if (tv_collider_collision_time < collision_time) {
			collision_time = tv_collider_collision_time;
		}
		if (vt_collider_collision_time < collision_time) {
			collision_time = vt_collider_collision_time;
		}

		//if (tv_collider_collision_time == 0) {
		//	std::cout << "tv collider collision color time zero " << std::endl;
		//}
		//if (vt_collider_collision_time == 0) {
		//	std::cout << "vt collider collision color time zero " << std::endl;
		//}
	}




	//floor:
	if (floor->exist) {
		std::array<double, 3>* q_pre;
		std::array<double, 3>* q_end;
		int* involved;
		int index_end;
		unsigned int* surface_to_global;
		unsigned int dimension; bool direction;	double floor_value;
		for (unsigned int i = 0; i < total_obj_num; ++i) {
			q_end = vertex_position[i];
			q_pre = vertex_record_for_this_color[i];
			involved = indicate_if_involved_in_last_color[i][0].data();
			index_end = vertex_index_start_per_thread[i][thread_No + 1];
			surface_to_global = vertex_index_on_surface[i + cloth->size()];
			dimension = floor->dimension;
			direction = floor->normal_direction;
			floor_value = floor->value;

			for (unsigned int j = vertex_index_start_per_thread[i][thread_No]; j < index_end; ++j) {
				if (involved[surface_to_global[j]]) {
					floorCollisionTime(q_pre[surface_to_global[j]].data(), q_end[surface_to_global[j]].data(), dimension, direction, floor_value,
						collision_time, d_hat_2);
				}
			}
		}
	}
	
	collision_time_thread[thread_No] = collision_time;
}





void Collision::collisionTimeByElement(int thread_No)
{
	double collision_time = 2.0;
	unsigned int element_end;
	std::array<double, 3>* initial_pos;
	std::array<double, 3>* current_pos;

	unsigned int* element;
	unsigned int* element_num;

	unsigned int* surface_to_normal;

	double record_previous_collision_time = 1.0;
	//VT
	bool* need_to_check;
	for (int i = 0; i < total_obj_num; ++i) {
		element_end = vertex_index_start_per_thread[i][thread_No + 1];
		initial_pos = vertex_collision_free[i];
		current_pos = vertex_position[i];
		element = spatial_hashing.vertex_triangle_pair_by_vertex[i];
		element_num = spatial_hashing.vertex_triangle_pair_num_record[i];

		need_to_check = spatial_hashing.is_used_vertex_triangle_pair_by_vertex[i];
		if(i<cloth->size()){
			for (int j = vertex_index_start_per_thread[i][thread_No]; j < element_end; ++j) {
				VTCollisionTimeOneVertex(initial_pos[j].data(), current_pos[j].data(), collision_time, element_num[j],
					element + estimate_coeff_for_vt_pair_num * j, vertex_collision_free.data(), vertex_position.data(),
					triangle_indices.data(), need_to_check + estimate_coeff_for_vt_pair_num_exist * j);
			}
		}
		else {
			unsigned int j;
			surface_to_normal = vertex_index_on_surface[i];
			for (int k = vertex_index_start_per_thread[i][thread_No]; k < element_end; ++k) {
				j = surface_to_normal[k];
				VTCollisionTimeOneVertex(initial_pos[j].data(), current_pos[j].data(), collision_time, element_num[k],
					element + estimate_coeff_for_vt_pair_num * k, vertex_collision_free.data(), vertex_position.data(), triangle_indices.data(),
					need_to_check + estimate_coeff_for_vt_pair_num_exist * k);
				if (collision_time == 0) {
					std::cout << "sefl vt collision time zero" << std::endl;
				}

			}

		}
	}
	record_previous_collision_time = collision_time;

	////EE
	unsigned int* edge_vertex;
	for (int i = 0; i < total_obj_num; ++i) {
		element_end = edge_index_start_per_thread[i][thread_No + 1];
		initial_pos = vertex_collision_free[i];
		current_pos = vertex_position[i];
		element = spatial_hashing.edge_edge_pair_by_edge[i];
		element_num = spatial_hashing.edge_edge_pair_num_record[i];		
		need_to_check = spatial_hashing.is_used_edge_edge_pair_by_edge[i];
		for (int j = edge_index_start_per_thread[i][thread_No]; j < element_end; ++j) {
			edge_vertex=edge_vertices[i]+(j<<1);
			EECollisionTimeOneEdge(initial_pos[*edge_vertex].data(),initial_pos[*(edge_vertex+1)].data(),
				current_pos[*edge_vertex].data(), current_pos[*(edge_vertex+1)].data(), 
				collision_time,j,i, element_num[j],
				element + estimate_coeff_for_ee_pair_num * j, vertex_collision_free.data(),vertex_position.data(),true,edge_vertices.data(),
				need_to_check + estimate_coeff_for_ee_pair_num_exist * j);
		}
	}
	
	if (has_collider) {
		//VT collider
		for (int i = 0; i < total_obj_num; ++i) {
			element_end = vertex_index_start_per_thread[i][thread_No + 1];
			initial_pos = vertex_collision_free[i];
			current_pos = vertex_position[i];
			element = spatial_hashing.vertex_obj_triangle_collider_pair_by_vertex[i];
			element_num = spatial_hashing.vertex_obj_triangle_collider_num_record[i];
			need_to_check = spatial_hashing.is_used_vertex_obj_triangle_collider_pair_by_vertex[i];
			if (i < cloth->size()) {
				for (int j = vertex_index_start_per_thread[i][thread_No]; j < element_end; ++j) {
					VTCollisionTimeOneVertex(initial_pos[j].data(), current_pos[j].data(), collision_time, element_num[j],
						element + estimate_coeff_for_vt_collider_pair_num * j, vertex_for_render_collider.data(), 
						vertex_position_collider.data(), triangle_indices_collider.data(), need_to_check + estimate_coeff_for_vt_collider_pair_num_exist * j);
				}
			}
			else {
				unsigned int j;
				surface_to_normal = vertex_index_on_surface[i];
				for (int k = vertex_index_start_per_thread[i][thread_No]; k < element_end; ++k) {
					j = surface_to_normal[k];
					VTCollisionTimeOneVertex(initial_pos[j].data(), current_pos[j].data(), collision_time, element_num[k],
						element + estimate_coeff_for_vt_collider_pair_num * k, vertex_for_render_collider.data(), 
						vertex_position_collider.data(),triangle_indices_collider.data(),
						need_to_check + estimate_coeff_for_vt_collider_pair_num_exist * k);
					if (collision_time == 0 && record_previous_collision_time!=0.0) {
						std::cout << "vt collider collision time zero" << std::endl;
					}
				}
			}
		}

		record_previous_collision_time = collision_time;

		//EE collider
		for (int i = 0; i < total_obj_num; ++i) {
			element_end = edge_index_start_per_thread[i][thread_No + 1];
			initial_pos = vertex_collision_free[i];
			current_pos = vertex_position[i];
			element = spatial_hashing.edge_obj_edge_collider_pair_by_edge[i];
			element_num = spatial_hashing.edge_obj_edge_collider_num_record[i];
			need_to_check = spatial_hashing.is_used_edge_obj_edge_collider_pair_by_edge[i];
			for (int j = edge_index_start_per_thread[i][thread_No]; j < element_end; ++j) {
				edge_vertex = edge_vertices[i] + (j << 1);
				EECollisionTimeOneEdge(initial_pos[*edge_vertex].data(), initial_pos[*(edge_vertex + 1)].data(),
					current_pos[*edge_vertex].data(), current_pos[*(edge_vertex + 1)].data(),
					collision_time, j, i, element_num[j],
					element + estimate_coeff_for_ee_collider_pair_num * j, vertex_for_render_collider.data(), vertex_position_collider.data(),
					false, collider_edge_vertices.data(), need_to_check+ estimate_coeff_for_ee_collider_pair_num_exist * j);
			}
		}

		//TV collider
		int* triangle_vertex;
		for (int i = 0; i < total_obj_num; ++i) {
			element_end = triangle_index_start_per_thread[i][thread_No + 1];
			initial_pos = vertex_collision_free[i];
			current_pos = vertex_position[i];
			element = spatial_hashing.triangle_obj_vertex_collider_pair_by_triangle[i];
			element_num = spatial_hashing.triangle_obj_vertex_collider_num_record[i];
			need_to_check = spatial_hashing.is_used_triangle_obj_vertex_collider_pair_by_triangle[i];
			for (int j = triangle_index_start_per_thread[i][thread_No]; j < element_end; ++j) {
				triangle_vertex = triangle_indices[i][j].data();
				TVCollisionTimeOneTriangle(initial_pos[*triangle_vertex].data(), initial_pos[*(triangle_vertex + 1)].data(),
					initial_pos[*(triangle_vertex + 2)].data(),
					current_pos[*triangle_vertex].data(), current_pos[*(triangle_vertex + 1)].data(),
					current_pos[*(triangle_vertex + 2)].data(),
					collision_time, element_num[j],
					element + estimate_coeff_for_tv_collider_pair_num * j, vertex_for_render_collider.data(), vertex_position_collider.data(),
					need_to_check + estimate_coeff_for_tv_collider_pair_num_exist * j,2);
				if (collision_time == 0 && record_previous_collision_time!=0.0) {
					std::cout << "tv collider collision time zero" << std::endl;
				}

			}
		}

	}



	//floor
	if (floor->exist) {
		for (int i = 0; i < total_obj_num; ++i) {
			element_end = vertex_index_start_per_thread[i][thread_No + 1];
			initial_pos = vertex_collision_free[i];
			current_pos = vertex_position[i];
			element = spatial_hashing.vertex_obj_triangle_collider_pair_by_vertex[i];
			element_num = spatial_hashing.vertex_obj_triangle_collider_num_record[i];
			if (i < cloth->size()) {
				for (int j = vertex_index_start_per_thread[i][thread_No]; j < element_end; ++j) {
					floorCollisionTime(initial_pos[j].data(), current_pos[j].data(), floor->dimension,
						floor->normal_direction, floor->value, collision_time,tolerance);
				}
			}
			else {
				unsigned int j;
				surface_to_normal = vertex_index_on_surface[i];
				for (int k = vertex_index_start_per_thread[i][thread_No]; k < element_end; ++k) {
					j = surface_to_normal[k];
					floorCollisionTime(initial_pos[j].data(), current_pos[j].data(), floor->dimension,
						floor->normal_direction, floor->value, collision_time,tolerance);
				}
			}

		}
	}

	collision_time_thread[thread_No] = collision_time;

	

}


//GLOBAL_COLLISION_TIME_ADD_PAIR
void Collision::collisionTimeWithPair(int thread_No)
{
	//VT
	if (vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_triangle_pair_index_start_per_thread[thread_No << 1]) {
		vertexTriangleCollisionTimePair(	vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]][0], collision_time,
			triangle_indices.data(),vertex_collision_free.data(), vertex_collision_free.data(),
			vertex_position.data(), vertex_position.data(),
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]]+1,
			&record_vt_pair[thread_No], &record_vt_pair_d_hat[thread_No], vt_hash_size_record.data(),
			vt_hash_record.data(), vertex_index_prefix_sum_obj.data(), triangle_index_prefix_sum_obj.data(),vertex_for_render.data(),
			vertex_for_render.data());
		for (int i = vertex_triangle_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			vertexTriangleCollisionTimePair(0,
				spatial_hashing.vertex_triangle_pair[i][0], collision_time, 
				triangle_indices.data(), vertex_collision_free.data(), vertex_collision_free.data(),
				vertex_position.data(), vertex_position.data(),
				spatial_hashing.vertex_triangle_pair[i] + 1,
				&record_vt_pair[thread_No], &record_vt_pair_d_hat[thread_No], vt_hash_size_record.data(),
				vt_hash_record.data(), vertex_index_prefix_sum_obj.data(), triangle_index_prefix_sum_obj.data(), vertex_for_render.data(),
				vertex_for_render.data());
		}
		vertexTriangleCollisionTimePair(0,
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time,
			triangle_indices.data(), vertex_collision_free.data(), vertex_collision_free.data(),
			vertex_position.data(), vertex_position.data(),
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1]] + 1,
			&record_vt_pair[thread_No], &record_vt_pair_d_hat[thread_No], vt_hash_size_record.data(),
			vt_hash_record.data(), vertex_index_prefix_sum_obj.data(), triangle_index_prefix_sum_obj.data(), vertex_for_render.data(),
			vertex_for_render.data());
	}
	else {
		vertexTriangleCollisionTimePair(vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time,
			triangle_indices.data(), vertex_collision_free.data(), vertex_collision_free.data(),
			vertex_position.data(), vertex_position.data(),
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]] + 1,
			&record_vt_pair[thread_No], &record_vt_pair_d_hat[thread_No], vt_hash_size_record.data(),
			vt_hash_record.data(), vertex_index_prefix_sum_obj.data(), triangle_index_prefix_sum_obj.data(), vertex_for_render.data(),
			vertex_for_render.data());
	}

	//EE
	if (edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1] > edge_edge_pair_index_start_per_thread[thread_No << 1]) {
		edgeEdgeCollisionTimePair(edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.edge_edge_pair[edge_edge_pair_index_start_per_thread[thread_No << 1]][0], collision_time,
			spatial_hashing.edge_edge_pair[edge_edge_pair_index_start_per_thread[thread_No << 1]] + 1,
			vertex_collision_free.data(), vertex_position.data(), vertex_collision_free.data(), vertex_position.data(),
			vertex_for_render.data(), vertex_for_render.data(), edge_vertices.data(), edge_vertices.data(),
			&record_ee_pair[thread_No], &record_ee_pair_d_hat[thread_No], ee_hash_size_record.data(),
			ee_hash_record.data(), edge_index_prefix_sum_obj.data(), edge_index_prefix_sum_obj.data());

		for (int i = edge_edge_pair_index_start_per_thread[thread_No << 1] + 1;
			i < edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			edgeEdgeCollisionTimePair(0, spatial_hashing.edge_edge_pair[i][0], collision_time,
				spatial_hashing.edge_edge_pair[i] + 1,
				vertex_collision_free.data(), vertex_position.data(), vertex_collision_free.data(), vertex_position.data(),
				vertex_for_render.data(), vertex_for_render.data(), edge_vertices.data(), edge_vertices.data(),
				&record_ee_pair[thread_No], &record_ee_pair_d_hat[thread_No], ee_hash_size_record.data(),
				ee_hash_record.data(), edge_index_prefix_sum_obj.data(), edge_index_prefix_sum_obj.data());
		}
		edgeEdgeCollisionTimePair(0,edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time,
			spatial_hashing.edge_edge_pair[edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1]] + 1,
			vertex_collision_free.data(), vertex_position.data(), vertex_collision_free.data(), vertex_position.data(),
			vertex_for_render.data(), vertex_for_render.data(), edge_vertices.data(), edge_vertices.data(),
			&record_ee_pair[thread_No], &record_ee_pair_d_hat[thread_No], ee_hash_size_record.data(),
			ee_hash_record.data(), edge_index_prefix_sum_obj.data(), edge_index_prefix_sum_obj.data());
	}
	else {
		edgeEdgeCollisionTimePair( edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time,
			spatial_hashing.edge_edge_pair[edge_edge_pair_index_start_per_thread[thread_No << 1]] + 1,
			vertex_collision_free.data(), vertex_position.data(), vertex_collision_free.data(), vertex_position.data(),
			vertex_for_render.data(), vertex_for_render.data(), edge_vertices.data(), edge_vertices.data(),
			&record_ee_pair[thread_No], &record_ee_pair_d_hat[thread_No], ee_hash_size_record.data(),
			ee_hash_record.data(), edge_index_prefix_sum_obj.data(), edge_index_prefix_sum_obj.data());
	}

	if (has_collider) {
		//tv collider
		if (vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1]) {
			vertexTriangleCollisionTimePair(vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 1],
				spatial_hashing.vertex_collider_triangle_obj_pair[vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1]][0], collision_time,
				triangle_indices.data(), vertex_for_render_collider.data(), vertex_collision_free.data(),
				vertex_position_collider.data(), vertex_position.data(),
				spatial_hashing.vertex_collider_triangle_obj_pair[vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1]] + 1,
				&record_tv_collider_pair[thread_No], &record_tv_collider_pair_d_hat[thread_No], tv_collider_hash_size_record.data(),
				tv_collider_hash_record.data(), vertex_index_prefix_sum_obj_collider.data(), triangle_index_prefix_sum_obj.data(), vertex_for_render_collider.data(),
				vertex_for_render.data());
			for (int i = vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1] + 1;
				i < vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
				vertexTriangleCollisionTimePair(0,	spatial_hashing.vertex_collider_triangle_obj_pair[i][0], collision_time,
					triangle_indices.data(), vertex_for_render_collider.data(), vertex_collision_free.data(),
					vertex_position_collider.data(), vertex_position.data(),
					spatial_hashing.vertex_collider_triangle_obj_pair[i] + 1,
					&record_tv_collider_pair[thread_No], &record_tv_collider_pair_d_hat[thread_No], tv_collider_hash_size_record.data(),
					tv_collider_hash_record.data(), vertex_index_prefix_sum_obj_collider.data(), triangle_index_prefix_sum_obj.data(), vertex_for_render_collider.data(),
					vertex_for_render.data());

			}
			vertexTriangleCollisionTimePair(0, vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time,
				triangle_indices.data(), vertex_for_render_collider.data(), vertex_collision_free.data(),
				vertex_position_collider.data(), vertex_position.data(),
				spatial_hashing.vertex_collider_triangle_obj_pair[vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1]] + 1,
				&record_tv_collider_pair[thread_No], &record_tv_collider_pair_d_hat[thread_No], tv_collider_hash_size_record.data(),
				tv_collider_hash_record.data(), vertex_index_prefix_sum_obj_collider.data(), triangle_index_prefix_sum_obj.data(), vertex_for_render_collider.data(),
				vertex_for_render.data());
		}
		else {
			vertexTriangleCollisionTimePair(vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 1],
				vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time,
				triangle_indices.data(), vertex_for_render_collider.data(), vertex_collision_free.data(),
				vertex_position_collider.data(), vertex_position.data(),
				spatial_hashing.vertex_collider_triangle_obj_pair[vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1]] + 1,
				&record_tv_collider_pair[thread_No], &record_tv_collider_pair_d_hat[thread_No], tv_collider_hash_size_record.data(),
				tv_collider_hash_record.data(), vertex_index_prefix_sum_obj_collider.data(), triangle_index_prefix_sum_obj.data(), vertex_for_render_collider.data(),
				vertex_for_render.data());
		}

		//vt collider
		if (vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]) {
			vertexTriangleCollisionTimePair(vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
				spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]][0], 
				collision_time,
				triangle_indices_collider.data(), vertex_collision_free.data(), vertex_for_render_collider.data(),
				vertex_position.data(), vertex_position_collider.data(),
				spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]] + 1,
				&record_vt_collider_pair[thread_No], &record_vt_collider_pair_d_hat[thread_No], vt_collider_hash_size_record.data(),
				vt_collider_hash_record.data(), vertex_index_prefix_sum_obj.data(), triangle_index_prefix_sum_obj_collider.data(), vertex_for_render.data(),
				vertex_for_render_collider.data());


			for (int i = vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1] + 1;
				i < vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {

				vertexTriangleCollisionTimePair(0,spatial_hashing.vertex_obj_triangle_collider_pair[i][0],collision_time,
					triangle_indices_collider.data(), vertex_collision_free.data(), vertex_for_render_collider.data(),
					vertex_position.data(), vertex_position_collider.data(),
					spatial_hashing.vertex_obj_triangle_collider_pair[i] + 1,
					&record_vt_collider_pair[thread_No], &record_vt_collider_pair_d_hat[thread_No], vt_collider_hash_size_record.data(),
					vt_collider_hash_record.data(), vertex_index_prefix_sum_obj.data(), triangle_index_prefix_sum_obj_collider.data(), vertex_for_render.data(),
					vertex_for_render_collider.data());

			}
			vertexTriangleCollisionTimePair(0, vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time,
				triangle_indices_collider.data(), vertex_collision_free.data(), vertex_for_render_collider.data(),
				vertex_position.data(), vertex_position_collider.data(),
				spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1]] + 1,
				&record_vt_collider_pair[thread_No], &record_vt_collider_pair_d_hat[thread_No], vt_collider_hash_size_record.data(),
				vt_collider_hash_record.data(), vertex_index_prefix_sum_obj.data(), triangle_index_prefix_sum_obj_collider.data(), vertex_for_render.data(),
				vertex_for_render_collider.data());

		}
		else {
			vertexTriangleCollisionTimePair(vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
				vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3],
				collision_time,
				triangle_indices_collider.data(), vertex_collision_free.data(), vertex_for_render_collider.data(),
				vertex_position.data(), vertex_position_collider.data(),
				spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]] + 1,
				&record_vt_collider_pair[thread_No], &record_vt_collider_pair_d_hat[thread_No], vt_collider_hash_size_record.data(),
				vt_collider_hash_record.data(), vertex_index_prefix_sum_obj.data(), triangle_index_prefix_sum_obj_collider.data(), vertex_for_render.data(),
				vertex_for_render_collider.data());
		}

		//ee collider
		if (edge_edge_pair_collider_index_start_per_thread[(thread_No + 1) << 1] > edge_edge_pair_collider_index_start_per_thread[thread_No << 1]) {
			
			edgeEdgeCollisionTimePair(edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 1],
				spatial_hashing.edge_edge_pair_collider[edge_edge_pair_collider_index_start_per_thread[thread_No << 1]][0], collision_time,
				spatial_hashing.edge_edge_pair_collider[edge_edge_pair_collider_index_start_per_thread[thread_No << 1]] + 1,
				vertex_collision_free.data(), vertex_position.data(), vertex_for_render_collider.data(), 
				vertex_position_collider.data(),
				vertex_for_render.data(), vertex_for_render_collider.data(), edge_vertices.data(), collider_edge_vertices.data(),
				&record_ee_collider_pair[thread_No], &record_ee_collider_pair_d_hat[thread_No], ee_collider_hash_size_record.data(),
				ee_collider_hash_record.data(), edge_index_prefix_sum_obj.data(), edge_index_prefix_sum_obj_collider.data());


			for (int i = edge_edge_pair_collider_index_start_per_thread[thread_No << 1] + 1;
				i < edge_edge_pair_collider_index_start_per_thread[(thread_No + 1) << 1]; ++i) {

				edgeEdgeCollisionTimePair(0,spatial_hashing.edge_edge_pair_collider[i][0], collision_time,
					spatial_hashing.edge_edge_pair_collider[i] + 1,
					vertex_collision_free.data(), vertex_position.data(), vertex_for_render_collider.data(),
					vertex_position_collider.data(),
					vertex_for_render.data(), vertex_for_render_collider.data(), edge_vertices.data(), collider_edge_vertices.data(),
					&record_ee_collider_pair[thread_No], &record_ee_collider_pair_d_hat[thread_No], ee_collider_hash_size_record.data(),
					ee_collider_hash_record.data(), edge_index_prefix_sum_obj.data(), edge_index_prefix_sum_obj_collider.data());

			}
			edgeEdgeCollisionTimePair(0, edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 3], collision_time,
				spatial_hashing.edge_edge_pair_collider[edge_edge_pair_collider_index_start_per_thread[(thread_No + 1) << 1]] + 1,
				vertex_collision_free.data(), vertex_position.data(), vertex_for_render_collider.data(),
				vertex_position_collider.data(),
				vertex_for_render.data(), vertex_for_render_collider.data(), edge_vertices.data(), collider_edge_vertices.data(),
				&record_ee_collider_pair[thread_No], &record_ee_collider_pair_d_hat[thread_No], ee_collider_hash_size_record.data(),
				ee_collider_hash_record.data(), edge_index_prefix_sum_obj.data(), edge_index_prefix_sum_obj_collider.data());

		}
		else {
			edgeEdgeCollisionTimePair(edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 1],
				edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 3], collision_time,
				spatial_hashing.edge_edge_pair_collider[edge_edge_pair_collider_index_start_per_thread[thread_No << 1]] + 1,
				vertex_collision_free.data(), vertex_position.data(), vertex_for_render_collider.data(),
				vertex_position_collider.data(),
				vertex_for_render.data(), vertex_for_render_collider.data(), edge_vertices.data(), collider_edge_vertices.data(),
				&record_ee_collider_pair[thread_No], &record_ee_collider_pair_d_hat[thread_No], ee_collider_hash_size_record.data(),
				ee_collider_hash_record.data(), edge_index_prefix_sum_obj.data(), edge_index_prefix_sum_obj_collider.data());
		}
	}
	
	if (floor->exist) {
		unsigned int element_end;
		std::array<double, 3>* initial_pos;
		std::array<double, 3>* current_pos;

		unsigned int* element;
		unsigned int* element_num;

		unsigned int* surface_to_normal;
		char* with_floor;
		std::vector<double>* d_hat_with_floor = &record_vertex_collide_with_floor_d_hat[thread_No];

		std::vector<unsigned int>* record_with_floor= &record_vertex_collide_with_floor[thread_No];

		for (int i = 0; i < total_obj_num; ++i) {
			element_end = vertex_index_start_per_thread[i][thread_No + 1];
			initial_pos = vertex_collision_free[i];
			current_pos = vertex_position[i];
			with_floor = indicate_vertex_collide_with_floor[i].data();
			if (i < cloth->size()) {
				for (int j = vertex_index_start_per_thread[i][thread_No]; j < element_end; ++j) {
					if(floorCollisionTime(initial_pos[j].data(), current_pos[j].data(), floor->dimension,
						floor->normal_direction, floor->value, collision_time, tolerance)) {
						if (!with_floor[j]) {
							with_floor[j] = '\1';
							record_with_floor->emplace_back(i);
							record_with_floor->emplace_back(j);
							d_hat_with_floor->emplace_back((vertex_for_render[i][j][floor->dimension] - floor->value)*
								vertex_for_render[i][j][floor->dimension] - floor->value);
						}
					}
				}
			}
			else {
				unsigned int j;
				surface_to_normal = vertex_index_on_surface[i];
				for (int k = vertex_index_start_per_thread[i][thread_No]; k < element_end; ++k) {
					j = surface_to_normal[k];
					floorCollisionTime(initial_pos[j].data(), current_pos[j].data(), floor->dimension,
						floor->normal_direction, floor->value, collision_time, tolerance);
					if (!with_floor[k]) {
						with_floor[k] = '\1';
						record_with_floor->emplace_back(i);
						record_with_floor->emplace_back(j);
						d_hat_with_floor->emplace_back((std::max)((vertex_for_render[i][j][floor->dimension] - floor->value) *
							(vertex_for_render[i][j][floor->dimension] - floor->value),d_hat_2));
					}
				}
			}
		}
	}

	collision_time_thread[thread_No] = collision_time;
}


void Collision::collisionTimeByPair(int thread_No)
{
	double collision_time = 2.0;
	//vertex_triangle
	if (vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_triangle_pair_index_start_per_thread[thread_No << 1]) {
		vertexTriangleCollisionTime(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.vertex_triangle_pair[vertex_triangle_pair_index_start_per_thread[thread_No << 1]][0], collision_time);
		for (int i = vertex_triangle_pair_index_start_per_thread[thread_No << 1] + 1;
			i < vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			vertexTriangleCollisionTime(thread_No, i, 0,
				spatial_hashing.vertex_triangle_pair[i][0], collision_time);
		}
		vertexTriangleCollisionTime(thread_No, vertex_triangle_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
	}
	else {
		vertexTriangleCollisionTime(thread_No, vertex_triangle_pair_index_start_per_thread[thread_No << 1], vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 1],
			vertex_triangle_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
	}
	//edge edge
	if (edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1] > edge_edge_pair_index_start_per_thread[thread_No << 1]) {
		edgeEdgeCollisionTime(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			spatial_hashing.edge_edge_pair[edge_edge_pair_index_start_per_thread[thread_No << 1]][0], collision_time);
		for (int i = edge_edge_pair_index_start_per_thread[thread_No << 1] + 1;
			i < edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
			edgeEdgeCollisionTime(thread_No, i, 0,
				spatial_hashing.edge_edge_pair[i][0], collision_time);
		}
		edgeEdgeCollisionTime(thread_No, edge_edge_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
	}
	else {
		edgeEdgeCollisionTime(thread_No, edge_edge_pair_index_start_per_thread[thread_No << 1], edge_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
			edge_edge_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
	}

	//// vertex edge
	//if (vertex_edge_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_edge_pair_index_start_per_thread[thread_No << 1]) {
	//	vertexEdgeCollisionTime(thread_No, vertex_edge_pair_index_start_per_thread[thread_No << 1], vertex_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
	//		vertex_edge_pair[vertex_edge_pair_index_start_per_thread[thread_No << 1]][0], collision_time, VERTEX_EDGE);
	//	for (int i = vertex_edge_pair_index_start_per_thread[thread_No << 1] + 1;
	//		i < vertex_edge_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
	//		vertexEdgeCollisionTime(thread_No, i, 0,
	//			vertex_edge_pair[i][0], collision_time, VERTEX_EDGE);
	//	}
	//	vertexEdgeCollisionTime(thread_No, vertex_edge_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
	//		vertex_edge_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_EDGE);
	//}
	//else {
	//	vertexEdgeCollisionTime(thread_No, vertex_edge_pair_index_start_per_thread[thread_No << 1], vertex_edge_pair_index_start_per_thread[(thread_No << 1) + 1],
	//		vertex_edge_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_EDGE);
	//}
	//////vertex vertex
	//if (vertex_vertex_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_vertex_pair_index_start_per_thread[thread_No << 1]) {
	//	vertexVertexCollisionTime(thread_No, vertex_vertex_pair_index_start_per_thread[thread_No << 1], vertex_vertex_pair_index_start_per_thread[(thread_No << 1) + 1],
	//		vertex_vertex_pair[vertex_vertex_pair_index_start_per_thread[thread_No << 1]][0], collision_time, VERTEX_VERTEX);
	//	for (int i = vertex_vertex_pair_index_start_per_thread[thread_No << 1] + 1;
	//		i < vertex_vertex_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
	//		vertexVertexCollisionTime(thread_No, i, 0,
	//			vertex_vertex_pair[i][0], collision_time, VERTEX_VERTEX);
	//	}
	//	vertexVertexCollisionTime(thread_No, vertex_vertex_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
	//		vertex_vertex_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_VERTEX);
	//}
	//else {
	//	vertexVertexCollisionTime(thread_No, vertex_vertex_pair_index_start_per_thread[thread_No << 1], vertex_vertex_pair_index_start_per_thread[(thread_No << 1) + 1],
	//		vertex_vertex_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_VERTEX);
	//}

	if (has_collider) {
		//vertex_collider_triangle
		if (vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1]) {
			vertexColliderTriangleCollisionTime(thread_No, vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1], vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 1],
				spatial_hashing.vertex_collider_triangle_obj_pair[vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1]][0], collision_time);
			for (int i = vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1] + 1;
				i < vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
				vertexColliderTriangleCollisionTime(thread_No, i, 0,
					spatial_hashing.vertex_collider_triangle_obj_pair[i][0], collision_time);
			}
			vertexColliderTriangleCollisionTime(thread_No, vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
				vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
		}
		else {
			vertexColliderTriangleCollisionTime(thread_No, vertex_collider_triangle_obj_pair_index_start_per_thread[thread_No << 1], vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 1],
				vertex_collider_triangle_obj_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
		}

		//vertex_triangle_collider

		if (vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]) {
			vertexTriangleColliderCollisionTime(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
				spatial_hashing.vertex_obj_triangle_collider_pair[vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1]][0], collision_time);
			for (int i = vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1] + 1;
				i < vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
				vertexTriangleColliderCollisionTime(thread_No, i, 0,
					spatial_hashing.vertex_obj_triangle_collider_pair[i][0], collision_time);
			}
			vertexTriangleColliderCollisionTime(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
				vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
		}
		else {
			vertexTriangleColliderCollisionTime(thread_No, vertex_obj_triangle_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
				vertex_obj_triangle_collider_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time);
		}

		//edge_edge_collider
		if (edge_edge_pair_collider_index_start_per_thread[(thread_No + 1) << 1] > edge_edge_pair_collider_index_start_per_thread[thread_No << 1]) {
			edgeEdgeColliderCollisionTime(thread_No, edge_edge_pair_collider_index_start_per_thread[thread_No << 1], edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 1],
				spatial_hashing.edge_edge_pair_collider[edge_edge_pair_collider_index_start_per_thread[thread_No << 1]][0], collision_time);
			for (int i = edge_edge_pair_collider_index_start_per_thread[thread_No << 1] + 1;
				i < edge_edge_pair_collider_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
				edgeEdgeColliderCollisionTime(thread_No, i, 0,
					spatial_hashing.edge_edge_pair_collider[i][0], collision_time);
			}
			edgeEdgeColliderCollisionTime(thread_No, edge_edge_pair_collider_index_start_per_thread[(thread_No + 1) << 1], 0,
				edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 3], collision_time);
		}
		else {
			edgeEdgeColliderCollisionTime(thread_No, edge_edge_pair_collider_index_start_per_thread[thread_No << 1], edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 1],
				edge_edge_pair_collider_index_start_per_thread[(thread_No << 1) + 3], collision_time);
		}

		////vertex_collider_edge_obj
		//if (vertex_collider_edge_obj_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_collider_edge_obj_pair_index_start_per_thread[thread_No << 1]) {
		//	vertexEdgeCollisionTime(thread_No, vertex_collider_edge_obj_pair_index_start_per_thread[thread_No << 1], vertex_collider_edge_obj_pair_index_start_per_thread[(thread_No << 1) + 1],
		//		vertex_collider_edge_obj_pair[vertex_collider_edge_obj_pair_index_start_per_thread[thread_No << 1]][0], collision_time, VERTEX_COLLIDER_EDGE_OBJ);
		//	for (int i = vertex_collider_edge_obj_pair_index_start_per_thread[thread_No << 1] + 1;
		//		i < vertex_collider_edge_obj_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
		//		vertexEdgeCollisionTime(thread_No, i, 0,
		//			vertex_collider_edge_obj_pair[i][0], collision_time, VERTEX_COLLIDER_EDGE_OBJ);
		//	}
		//	vertexEdgeCollisionTime(thread_No, vertex_collider_edge_obj_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
		//		vertex_collider_edge_obj_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_COLLIDER_EDGE_OBJ);
		//}
		//else {
		//	vertexEdgeCollisionTime(thread_No, vertex_collider_edge_obj_pair_index_start_per_thread[thread_No << 1], vertex_collider_edge_obj_pair_index_start_per_thread[(thread_No << 1) + 1],
		//		vertex_collider_edge_obj_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_COLLIDER_EDGE_OBJ);
		//}
		////vertex_obj_edge_collider
		//if (vertex_obj_edge_collider_pair_index_start_per_thread[(thread_No + 1) << 1] > vertex_obj_edge_collider_pair_index_start_per_thread[thread_No << 1]) {
		//	vertexEdgeCollisionTime(thread_No, vertex_obj_edge_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_edge_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
		//		vertex_obj_edge_collider_pair[vertex_obj_edge_collider_pair_index_start_per_thread[thread_No << 1]][0], collision_time, VERTEX_OBJ_EDGE_COLLIDER);
		//	for (int i = vertex_obj_edge_collider_pair_index_start_per_thread[thread_No << 1] + 1;
		//		i < vertex_obj_edge_collider_pair_index_start_per_thread[(thread_No + 1) << 1]; ++i) {
		//		vertexEdgeCollisionTime(thread_No, i, 0,
		//			vertex_obj_edge_collider_pair[i][0], collision_time, VERTEX_OBJ_EDGE_COLLIDER);
		//	}
		//	vertexEdgeCollisionTime(thread_No, vertex_obj_edge_collider_pair_index_start_per_thread[(thread_No + 1) << 1], 0,
		//		vertex_obj_edge_collider_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_OBJ_EDGE_COLLIDER);
		//}
		//else {
		//	vertexEdgeCollisionTime(thread_No, vertex_obj_edge_collider_pair_index_start_per_thread[thread_No << 1], vertex_obj_edge_collider_pair_index_start_per_thread[(thread_No << 1) + 1],
		//		vertex_obj_edge_collider_pair_index_start_per_thread[(thread_No << 1) + 3], collision_time, VERTEX_OBJ_EDGE_COLLIDER);
		//}

	}

	if (floor->exist) {
		unsigned int element_end;
		std::array<double, 3>* initial_pos;
		std::array<double, 3>* current_pos;

		unsigned int* element;
		unsigned int* element_num;

		unsigned int* surface_to_normal;

		for (int i = 0; i < total_obj_num; ++i) {
			element_end = vertex_index_start_per_thread[i][thread_No + 1];
			initial_pos = vertex_for_render[i];
			current_pos = vertex_position[i];
			if (i < cloth->size()) {
				for (int j = vertex_index_start_per_thread[i][thread_No]; j < element_end; ++j) {
					floorCollisionTime(initial_pos[j].data(), current_pos[j].data(), floor->dimension,
						floor->normal_direction, floor->value, collision_time,tolerance);
				}
			}
			else {
				unsigned int j;
				surface_to_normal = vertex_index_on_surface[i];
				for (int k = vertex_index_start_per_thread[i][thread_No]; k < element_end; ++k) {
					j = surface_to_normal[k];
					floorCollisionTime(initial_pos[j].data(), current_pos[j].data(), floor->dimension,
						floor->normal_direction, floor->value, collision_time,tolerance);
				}
			}

		}
		////std::cout << "F collision time " << collision_time << std::endl;
	}

	collision_time_thread[thread_No] = collision_time;
}




//GLOBAL_COLLISION_TIME
void Collision::collisionTime(int thread_No)
{
	if (record_pair_by_element) {
		collisionTimeByElement(thread_No);
	}
	else {
		collisionTimeByPair(thread_No);
	}
	
}


//void Collision::initialCollisionTimeRecord()
//{
//	for (unsigned int i = 0; i < thread_num; ++i) {
//		initialCollisionTimeRecord(i);	
//	}
//}

//void Collision::initialCollisionTimeRecord(int thread_No)
//{
//	record_VT_collision_time[thread_No].resize(spatial_hashing.vertex_triangle_pair[thread_No][0]>>2);
//	record_EE_collision_time[thread_No].resize(spatial_hashing.edge_edge_pair[thread_No][0]>>2);
//	record_VTCollider_collision_time[thread_No].resize(spatial_hashing.vertex_obj_triangle_collider_pair[thread_No][0]>>2);
//}

void Collision::collisionTime(int thread_No)
{
	int thread_test = 2;
	double* collision_time = &collision_time_thread[thread_No];
	(*collision_time) = 2.0;
	TriangleMeshStruct* mesh_struct;
	std::vector<std::vector<int>>* neighbor_primitve;
	int index_begin;
	int index_end;
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		mesh_struct = &(*cloth)[cloth_No].mesh_struct;
		index_begin = mesh_struct->vertex_index_begin_per_thread[thread_No];
		index_end = mesh_struct->vertex_index_begin_per_thread[thread_No + 1];
		neighbor_primitve = (*cloth)[cloth_No].vertex_neighbor_obj_triangle.data();
		for (int i = index_begin; i < index_end; ++i) {
			pointSelfTriangleCollisionTime(collision_time, neighbor_primitve[i].data(), mesh_struct->vertex_for_render[i].data(), mesh_struct->vertex_position[i].data(), i);
		}
	}
	std::array<double, 3>* initial_pos;
	std::array<double, 3>* current_pos;
	for (int collider_No = 0; collider_No < collider->size(); ++collider_No) {
		mesh_struct = &(*collider)[collider_No].mesh_struct;
		index_begin = mesh_struct->face_index_begin_per_thread[thread_No];
		index_end = mesh_struct->face_index_begin_per_thread[thread_No + 1];
		initial_pos = mesh_struct->vertex_for_render.data();
		current_pos = mesh_struct->vertex_position.data();
		neighbor_primitve = (*collider)[collider_No].triangle_neighbor_obj_vertex.data();
		for (int i = index_begin; i < index_end; ++i) {
			pointColliderTriangleCollisionTime(collision_time, mesh_struct->triangle_indices[i].data(), neighbor_primitve[i].data(), initial_pos, current_pos,
				mesh_struct->ori_face_normal_for_render[i].data(), mesh_struct->ori_face_normal[i].data(), mesh_struct->cross_for_approx_CCD[i].data(),
				mesh_struct->f_face_normal_for_render[i].data(), mesh_struct->f_face_normal[i].data(), mesh_struct->f_cross_for_approx_CCD[i].data(), i);

		}

	}

	unsigned int* edge_vertex_index;
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		mesh_struct = &(*cloth)[cloth_No].mesh_struct;
		index_begin = mesh_struct->edge_index_begin_per_thread[thread_No];
		index_end = mesh_struct->edge_index_begin_per_thread[thread_No + 1];
		neighbor_primitve = (*cloth)[cloth_No].edge_neighbor_obj_edge.data();
		initial_pos = mesh_struct->vertex_for_render.data();
		current_pos = mesh_struct->vertex_position.data();
		for (int i = index_begin; i < index_end; ++i) {
			edge_vertex_index = mesh_struct->edge_vertices.data() + (i << 1);//edges[i].vertex;
			edgeEdgeCollisionTime(collision_time, neighbor_primitve[i].data(), initial_pos[edge_vertex_index[0]].data(), initial_pos[edge_vertex_index[1]].data(),
				current_pos[edge_vertex_index[0]].data(), current_pos[edge_vertex_index[1]].data());
		}
	}
}


void Collision::testCollision()
{
	//int k = 0;
	//for (int i = 0; i < (*cloth)[0].triangle_neighbor_obj_triangle.size(); ++i) {
	//	//(*cloth)[0].triangle_neighbor_cloth_triangle[i][0].push_back(i);
	//	k += (*cloth)[0].triangle_neighbor_obj_triangle[i][0].size();
	//}
	//////std::cout << k + (*cloth)[0].triangle_neighbor_cloth_triangle.size()<<std::endl;

	//for (int i = 0; i < (*cloth)[0].mesh_struct.vertices.size(); ++i) {
	//for (int j = 0; j < (*cloth)[0].vertex_neighbor_collider_triangle[i][0].size(); ++j) {
	//	//////std::cout << (*cloth)[0].vertex_AABB[i].min[0] << " " << (*cloth)[0].vertex_AABB[i].min[1] << " " << (*cloth)[0].vertex_AABB[i].min[2] << " "
	//		<< (*cloth)[0].vertex_AABB[i].max[0] << " " << (*cloth)[0].vertex_AABB[i].max[1] << " " << (*cloth)[0].vertex_AABB[i].max[2] << std::endl;
	//	//////std::cout << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].min[0] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].min[1] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].min[2] << " "
	//		<< (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].max[0] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].max[1] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].max[2] << std::endl;
	//	//////std::cout << i << " " << (*cloth)[0].vertex_neighbor_collider_triangle[i][0][j] << std::endl;
	//	/*	if (!(*cloth)[0].vertex_AABB[i].AABB_intersection((*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]])) {
	//			//////std::cout << i << " " << (*cloth)[0].vertex_neighbor_collider_triangle[i][0][j] << std::endl;
	//		}*/
	//}

	//if (!(*cloth)[0].vertex_neighbor_collider_triangle[i][0].empty()) {
	//	system("pause");
	//}
	//}
	//for (int i = 0; i < (*collider)[0].mesh_struct.triangle_indices.size(); ++i) {
	//	for (int j = 0; j < (*collider)[0].triangle_neighbor_obj_vertex[i][0].size(); ++j) {
	//		if ((*collider)[0].triangle_neighbor_obj_vertex[i][0][j] == 0) {
	//			if (!(*collider)[0].triangle_AABB[i].AABB_intersection((*cloth)[0].vertex_AABB[(*collider)[0].triangle_neighbor_obj_vertex[i][0][j]])) {
	//				//////std::cout << (*collider)[0].triangle_neighbor_obj_vertex[i][0][j] << " " << i << std::endl;
	//			}
	//			system("pause");
	//		}
	//	}
	//}
	//if (!(*cloth)[0].vertex_neighbor_collider_triangle[i][0].empty()) {
	//	
	//}
	//for (int i = 0; i < (*cloth)[0].mesh_struct.triangle_indices.size(); ++i) {
	//	for (int j = 0; j < (*cloth)[0].triangle_neighbor_collider_triangle[i][0].size(); ++j) {
	//		if (!(*cloth)[0].triangle_AABB[i].AABB_intersection((*collider)[0].triangle_AABB[(*cloth)[0].triangle_neighbor_collider_triangle[i][0][j]])) {
	//			//////std::cout << i << " " << (*cloth)[0].triangle_neighbor_collider_triangle[i][0][j] << std::endl;
	//		}
	//		
	//	}
	//}

}

void Collision::updateCollisionPosition()
{
	for (int i = 0; i < cloth->size(); ++i) {
		(*cloth)[i].mesh_struct.getNormal();
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		(*tetrahedron)[i].mesh_struct.getNormal();
	}
	thread->assignTask(this, RE_DETECTION);
	resumTargetPosition();
}



void Collision::resumTargetPosition()
{
	obj_target_pos.partialInitial();
	thread->assignTask(this, RESUM_TARGET_POSITION);
	for (int i = 0; i < thread_num; ++i) {
		obj_target_pos.collision_energy += obj_target_pos_per_thread[i].collision_energy;
	}
}

void Collision::sumTargetPosition()
{
	obj_target_pos.initial();
	for (int i = 0; i < thread_num; ++i) {
		obj_target_pos.collision_energy += obj_target_pos_per_thread[i].collision_energy;
	}

	thread->assignTask(this, SUM_TARGET_POSITION);

}


//SUM_TARGET_POSITION
void Collision::sumTargetPositionPerThread(int thread_id)
{
	unsigned int* index_begin;
	bool* need_update;
	std::vector<std::array<double, 3>>* b_sum;
	std::vector<std::array<double, 3>>* b_sum_per_thread;
	bool* global_need_update;
	double* stiffness;
	double* global_stiffness;

	for (int j = 0; j < tetrahedron_begin_obj_index; ++j) {
		global_need_update = obj_target_pos.need_update[j];
		b_sum = &obj_target_pos.b_sum[j];
		index_begin = (*cloth)[j].mesh_struct.vertex_index_begin_per_thread.data();
		global_stiffness = obj_target_pos.stiffness[j].data();
		for (int i = 0; i < thread_num; ++i) {
			need_update = obj_target_pos_per_thread[i].need_update[j];
			b_sum_per_thread = &obj_target_pos_per_thread[i].b_sum[j];
			stiffness = obj_target_pos_per_thread[i].stiffness[j].data();
			for (int k = index_begin[thread_id]; k < index_begin[thread_id + 1]; ++k) {
				if (need_update[k]) {
					global_need_update[k] = true;
					SUM_((*b_sum)[k], (*b_sum_per_thread)[k]);
					global_stiffness[k] += stiffness[k];
				}
			}
		}
	}

	unsigned int* surface_vertex_index;
	unsigned int vertex_index;

	for (int j = tetrahedron_begin_obj_index; j < total_obj_num; ++j) {
		global_need_update = obj_target_pos.need_update[j];
		b_sum = &obj_target_pos.b_sum[j];
		index_begin = (*tetrahedron)[j - tetrahedron_begin_obj_index].mesh_struct.vertex_index_on_surface_begin_per_thread.data();
		surface_vertex_index = vertex_index_on_surface[j];
		global_stiffness = obj_target_pos.stiffness[j].data();
		for (int i = 0; i < thread_num; ++i) {
			need_update = obj_target_pos_per_thread[i].need_update[j];
			b_sum_per_thread = &obj_target_pos_per_thread[i].b_sum[j];
			stiffness = obj_target_pos_per_thread[i].stiffness[j].data();
			for (int k = index_begin[thread_id]; k < index_begin[thread_id + 1]; ++k) {
				vertex_index = surface_vertex_index[k];
				if (need_update[vertex_index]) {
					global_need_update[vertex_index] = true;
					SUM_((*b_sum)[vertex_index], (*b_sum_per_thread)[vertex_index]);
					global_stiffness[vertex_index] += stiffness[vertex_index];
				}
			}
		}
	}
	////////std::cout << "collision vertex ";
	//for (int i = 0; i < (*cloth)[0].ori_vertices.size(); ++i) {
	//	if (obj_target_pos.need_update[0][i]) {
	//		//////std::cout << i << " ";
	//	}
	//}
	////////std::cout << std::endl;
}
//RESUM_TARGET_POSITION
void Collision::resumTargetPositionPerThread(int thread_id)
{
	unsigned int* index_begin;
	bool* need_update;
	std::vector<std::array<double, 3>>* b_sum;
	std::vector<std::array<double, 3>>* b_sum_per_thread;

	for (int j = 0; j < tetrahedron_begin_obj_index; ++j) {
		b_sum = &obj_target_pos.b_sum[j];
		index_begin = (*cloth)[j].mesh_struct.vertex_index_begin_per_thread.data();
		need_update = obj_target_pos.need_update[j];
		for (int i = 0; i < thread_num; ++i) {
			b_sum_per_thread = &obj_target_pos_per_thread[i].b_sum[j];
			for (int k = index_begin[thread_id]; k < index_begin[thread_id + 1]; ++k) {
				if (need_update[k]) {
					SUM_((*b_sum)[k], (*b_sum_per_thread)[k]);
				}
			}
		}
	}

	unsigned int* surface_vertex_index;
	unsigned int vertex_index;
	for (int j = tetrahedron_begin_obj_index; j < total_obj_num; ++j) {
		b_sum = &obj_target_pos.b_sum[j];
		index_begin = (*tetrahedron)[j - tetrahedron_begin_obj_index].mesh_struct.vertex_index_on_surface_begin_per_thread.data();
		surface_vertex_index = vertex_index_on_surface[j];
		need_update = obj_target_pos.need_update[j];
		for (int i = 0; i < thread_num; ++i) {
			b_sum_per_thread = &obj_target_pos_per_thread[i].b_sum[j];
			for (int k = index_begin[thread_id]; k < index_begin[thread_id + 1]; ++k) {
				vertex_index = surface_vertex_index[k];
				if (need_update[vertex_index]) {
					SUM_((*b_sum)[vertex_index], (*b_sum_per_thread)[vertex_index]);
				}
			}
		}
	}

}
//RE_DETECTION
void Collision::collisionReDetection(int thread_No)
{
	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
	target_pos->partialInitial();
	unsigned int* index_begin;

	std::vector<int>* vertex_neighbor_obj_triangle;
	std::vector<int>* collide_vertex_obj_triangle;
	std::vector<int>* vertex_neighbor_collider_triangle;
	std::vector<int>* collide_vertex_collider_triangle;
	int end;
	int obj_No;
	MeshStruct* mesh_struct;

	double PC_radius0;
	double PC_radius1;
	double vertex_collision_stiffness0;
	double vertex_collision_stiffness1;


	//if (use_BVH) {
	//	for (int j = 0; j < total_obj_num; ++j) {
	//		if (j<tetrahedron_begin_obj_index) {
	//			index_begin = (*cloth)[j].mesh_struct.vertex_index_begin_per_thread.data();
	//			end = index_begin[thread_No + 1];
	//			mesh_struct = &(*cloth)[j].mesh_struct;
	//			vertex_collision_stiffness0 = (*cloth)[j].collision_stiffness[SELF_POINT_TRIANGLE];
	//			vertex_collision_stiffness1 = (*cloth)[j].collision_stiffness[BODY_POINT_TRIANGLE];
	//		}
	//		else {
	//			obj_No = j - tetrahedron_begin_obj_index;
	//			index_begin = (*tetrahedron)[obj_No].mesh_struct.vertex_index_on_surface_begin_per_thread.data();
	//			end = index_begin[thread_No + 1];
	//			mesh_struct = &(*tetrahedron)[obj_No].mesh_struct;
	//			vertex_collision_stiffness0 = (*tetrahedron)[obj_No].collision_stiffness[SELF_POINT_TRIANGLE];
	//			vertex_collision_stiffness1 = (*tetrahedron)[obj_No].collision_stiffness[BODY_POINT_TRIANGLE];
	//		}	
	//		for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
	//			if (j < tetrahedron_begin_obj_index) {
	//				vertex_neighbor_obj_triangle = (*cloth)[j].vertex_neighbor_obj_triangle[i].data();
	//				collide_vertex_obj_triangle = (*cloth)[j].collide_vertex_obj_triangle[i].data();
	//				PC_radius0 = (*cloth)[j].PC_radius[SELF_POINT_TRIANGLE];
	//				PC_radius1 = (*cloth)[j].PC_radius[BODY_POINT_TRIANGLE];
	//				vertex_neighbor_collider_triangle = (*cloth)[j].vertex_neighbor_collider_triangle[i].data();
	//				collide_vertex_collider_triangle = (*cloth)[j].collide_vertex_collider_triangle[i].data();
	//			}
	//			else {
	//				vertex_neighbor_obj_triangle = (*tetrahedron)[obj_No].surface_vertex_neighbor_obj_triangle[i].data();
	//				collide_vertex_obj_triangle = (*tetrahedron)[obj_No].collide_vertex_obj_triangle[i].data();
	//				PC_radius0 = (*tetrahedron)[obj_No].PC_radius[SELF_POINT_TRIANGLE];
	//				PC_radius1 = (*tetrahedron)[obj_No].PC_radius[BODY_POINT_TRIANGLE];
	//				vertex_neighbor_collider_triangle = (*tetrahedron)[obj_No].surface_vertex_neighbor_collider_triangle[i].data();
	//				collide_vertex_collider_triangle = (*tetrahedron)[obj_No].collide_vertex_collider_triangle[i].data();
	//			}
	//			pointSelfTriangleCollisionReDetection(thread_No, i, j, collide_vertex_obj_triangle, &(*cloth)[j].mesh_struct,
	//				(*cloth)[j].PC_radius[SELF_POINT_TRIANGLE], (*cloth)[j].collision_stiffness[SELF_POINT_TRIANGLE], target_pos);
	//			//pointColliderTriangleCollisionReDetection(thread_No, i, j, &(*cloth)[j].collide_vertex_collider_triangle[i], &(*cloth)[j].mesh_struct,
	//			//	(*cloth)[j].PC_radius[BODY_POINT_TRIANGLE], (*cloth)[j].collision_stiffness[BODY_POINT_TRIANGLE], target_pos);
	//		}
	//	}
	//}
	//else {
	//	for (int j = 0; j < cloth->size(); ++j) {
	//		index_begin = (*cloth)[j].mesh_struct.vertex_index_begin_per_thread.data();
	//		for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
	//			pointSelfTriangleCollisionReDetection(thread_No, i, j, &(*cloth)[j].collide_vertex_cloth_triangle[i], &(*cloth)[j].mesh_struct,
	//				(*cloth)[j].PC_radius[SELF_POINT_TRIANGLE], &(*cloth)[j].collision_stiffness[SELF_POINT_TRIANGLE], target_pos);
	//		}
	//	}
	//	for (int collider_No = 0; collider_No < collider->size(); ++collider_No) {
	//		index_begin = (*collider)[collider_No].mesh_struct.face_index_begin_per_thread.data();
	//		for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
	//			colliderTriangleVertexCollisionReDetection(thread_No, i, collider_No, (*collider)[collider_No].collider_triangle_cloth_vertex[i].data(),
	//				&(*collider)[collider_No].mesh_struct, (*collider)[collider_No].tolerance, target_pos);
	//		}
	//	}
	//}
	//for (int j = 0; j < cloth->size(); ++j) {
	//	index_begin = (*cloth)[j].mesh_struct.edge_index_begin_per_thread.data();
	//	for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
	//		edgeSelfEdgeCollisionReDetection(thread_No, i, j, &(*cloth)[j].collide_edge_cloth_edge[i], &(*cloth)[j].mesh_struct,
	//			(*cloth)[j].PC_radius[SELF_EDGE_EDGE], (*cloth)[j].collision_stiffness[SELF_EDGE_EDGE], target_pos);
	//	}
	//}
}
//GLOBAL_COLLISION_DETECTION
void Collision::collisionDetection(int thread_No)
{
	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
	target_pos->initial();
	unsigned int* index_begin;

	std::vector<int>* vertex_neighbor_obj_triangle;
	std::vector<int>* collide_vertex_obj_triangle;
	std::vector<int>* vertex_neighbor_collider_triangle;
	std::vector<int>* collide_vertex_collider_triangle;
	unsigned int end;
	unsigned int obj_No;
	MeshStruct* mesh_struct;

	double PC_radius0;
	double PC_radius1;
	double vertex_collision_stiffness0;
	double vertex_collision_stiffness1;
	//if (use_BVH) {
	//	for (int j = 0; j < total_obj_num; ++j) {
	//		if (j < tetrahedron_begin_obj_index) {
	//			index_begin = (*cloth)[j].mesh_struct.vertex_index_begin_per_thread.data();
	//			end = index_begin[thread_No + 1];
	//			mesh_struct = &(*cloth)[j].mesh_struct;
	//			vertex_collision_stiffness0 = (*cloth)[j].collision_stiffness[SELF_POINT_TRIANGLE];
	//			vertex_collision_stiffness1 = (*cloth)[j].collision_stiffness[BODY_POINT_TRIANGLE];
	//		}
	//		else {
	//			obj_No = j - tetrahedron_begin_obj_index;
	//			index_begin = (*tetrahedron)[obj_No].mesh_struct.vertex_index_on_surface_begin_per_thread.data();
	//			end = index_begin[thread_No + 1];
	//			mesh_struct = &(*tetrahedron)[obj_No].mesh_struct;
	//			vertex_collision_stiffness0 = (*tetrahedron)[obj_No].collision_stiffness[SELF_POINT_TRIANGLE];
	//			vertex_collision_stiffness1 = (*tetrahedron)[obj_No].collision_stiffness[BODY_POINT_TRIANGLE];
	//		}			
	//		
	//		for (int i = index_begin[thread_No]; i < end; ++i) {
	//			if (j < tetrahedron_begin_obj_index) {
	//				vertex_neighbor_obj_triangle = (*cloth)[j].vertex_neighbor_obj_triangle[i].data();
	//				collide_vertex_obj_triangle = (*cloth)[j].collide_vertex_obj_triangle[i].data();					
	//				PC_radius0 = (*cloth)[j].PC_radius[SELF_POINT_TRIANGLE];
	//				PC_radius1 = (*cloth)[j].PC_radius[BODY_POINT_TRIANGLE];
	//				vertex_neighbor_collider_triangle = (*cloth)[j].vertex_neighbor_collider_triangle[i].data();
	//				collide_vertex_collider_triangle = (*cloth)[j].collide_vertex_collider_triangle[i].data();
	//			}
	//			else {
	//				vertex_neighbor_obj_triangle = (*tetrahedron)[obj_No].surface_vertex_neighbor_obj_triangle[i].data();
	//				collide_vertex_obj_triangle = (*tetrahedron)[obj_No].collide_vertex_obj_triangle[i].data();
	//				PC_radius0 = (*tetrahedron)[obj_No].PC_radius[SELF_POINT_TRIANGLE];
	//				PC_radius1 = (*tetrahedron)[obj_No].PC_radius[BODY_POINT_TRIANGLE];
	//				vertex_neighbor_collider_triangle = (*tetrahedron)[obj_No].surface_vertex_neighbor_collider_triangle[i].data();
	//				collide_vertex_collider_triangle = (*tetrahedron)[obj_No].collide_vertex_collider_triangle[i].data();
	//			}
	//			pointSelfTriangleCollisionDetection(thread_No, i, j, vertex_neighbor_obj_triangle,
	//				collide_vertex_obj_triangle, mesh_struct, PC_radius0, target_pos, vertex_collision_stiffness0);
	//			pointColliderTriangleCollisionDetection(thread_No, i, j, vertex_neighbor_collider_triangle,
	//				collide_vertex_collider_triangle, mesh_struct, PC_radius1, target_pos, vertex_collision_stiffness1);
	//		}
	//	}
	//}
	//else {
	for (int j = 0; j < total_obj_num; ++j) {
		if (j < tetrahedron_begin_obj_index) {
			index_begin = (*cloth)[j].mesh_struct.vertex_index_begin_per_thread.data();
			end = index_begin[thread_No + 1];
			mesh_struct = &(*cloth)[j].mesh_struct;
			vertex_collision_stiffness0 = (*cloth)[j].collision_stiffness[SELF_POINT_TRIANGLE];
		}
		else {
			obj_No = j - tetrahedron_begin_obj_index;
			index_begin = (*tetrahedron)[obj_No].mesh_struct.vertex_index_on_surface_begin_per_thread.data();
			end = index_begin[thread_No + 1];
			mesh_struct = &(*tetrahedron)[obj_No].mesh_struct;
			vertex_collision_stiffness0 = (*tetrahedron)[obj_No].collision_stiffness[SELF_POINT_TRIANGLE];
		}
		for (int i = index_begin[thread_No]; i < end; ++i) {
			if (j < tetrahedron_begin_obj_index) {
				vertex_neighbor_obj_triangle = (*cloth)[j].vertex_neighbor_obj_triangle[i].data();
				collide_vertex_obj_triangle = (*cloth)[j].collide_vertex_obj_triangle[i].data();
				PC_radius0 = (*cloth)[j].PC_radius[SELF_POINT_TRIANGLE];
				vertex_neighbor_collider_triangle = (*cloth)[j].vertex_neighbor_collider_triangle[i].data();
			}
			else {
				vertex_neighbor_obj_triangle = (*tetrahedron)[obj_No].surface_vertex_neighbor_obj_triangle[i].data();
				collide_vertex_obj_triangle = (*tetrahedron)[obj_No].collide_vertex_obj_triangle[i].data();
				PC_radius0 = (*tetrahedron)[obj_No].PC_radius[SELF_POINT_TRIANGLE];
				vertex_neighbor_collider_triangle = (*tetrahedron)[obj_No].surface_vertex_neighbor_collider_triangle[i].data();
			}
			pointSelfTriangleCollisionDetection(thread_No, i, j, vertex_neighbor_obj_triangle,
				collide_vertex_obj_triangle, mesh_struct, PC_radius0, target_pos, vertex_collision_stiffness0);
		}
	}
	for (int collider_No = 0; collider_No < collider->size(); ++collider_No) {
		index_begin = (*collider)[collider_No].mesh_struct.face_index_begin_per_thread.data();
		for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
			colliderTriangleVertexCollisionDetection(thread_No, i, collider_No, &(*collider)[collider_No].triangle_neighbor_obj_vertex[i],
				&(*collider)[collider_No].collider_triangle_obj_vertex[i], &(*collider)[collider_No].mesh_struct, (*collider)[collider_No].tolerance, target_pos);
		}
	}
	//}

	std::vector<int>* edge_neighbor_obj_edge;
	std::vector<int>* collide_edge_obj_edge;

	for (int j = 0; j < total_obj_num; ++j) {
		if (j < tetrahedron_begin_obj_index) {
			mesh_struct = &(*cloth)[j].mesh_struct;
			PC_radius0 = (*cloth)[j].PC_radius[SELF_EDGE_EDGE];
			vertex_collision_stiffness0 = (*cloth)[j].collision_stiffness[SELF_EDGE_EDGE];
			index_begin = (*cloth)[j].mesh_struct.edge_index_begin_per_thread.data();
		}
		else {
			obj_No = j - tetrahedron_begin_obj_index;
			mesh_struct = &(*tetrahedron)[obj_No].mesh_struct;
			PC_radius0 = (*tetrahedron)[obj_No].PC_radius[SELF_EDGE_EDGE];
			vertex_collision_stiffness0 = (*tetrahedron)[obj_No].collision_stiffness[SELF_EDGE_EDGE];
			index_begin = (*tetrahedron)[obj_No].mesh_struct.vertex_index_begin_per_thread.data();
		}

		for (int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
			if (j < tetrahedron_begin_obj_index) {
				vertex_neighbor_obj_triangle = (*cloth)[j].edge_neighbor_obj_edge[i].data();
				collide_edge_obj_edge = (*cloth)[j].collide_edge_obj_edge[i].data();
			}
			else {
				vertex_neighbor_obj_triangle = (*tetrahedron)[obj_No].edge_neighbor_obj_edge[i].data();
				collide_edge_obj_edge = (*tetrahedron)[obj_No].collide_edge_obj_edge[i].data();
			}
			edgeSelfEdgeCollisionDetection(thread_No, i, j, vertex_neighbor_obj_triangle, collide_edge_obj_edge,
				mesh_struct, PC_radius0, target_pos, vertex_collision_stiffness0);
		}
	}
}

void Collision::colliderTriangleVertexCollisionDetection(int thread_No, int triangle_index, int collider_No,
	std::vector<std::vector<int>>* triangle_neighbor_vertex, std::vector<std::vector<int>>* collide_triangle_vertex, MeshStruct* triangle_mesh, double radius0,
	TargetPosition* target_pos)
{
	MeshStruct* vertex_mesh;
	std::vector<int>* neighbor_vertex;
	std::vector<int>* collide_vertex;
	double radius1;
	double stiffness;
	int obj_No;
	for (int i = 0; i < total_obj_num; ++i) {
		if (i < tetrahedron_begin_obj_index) {
			vertex_mesh = &(*cloth)[i].mesh_struct;
			radius1 = (*cloth)[i].PC_radius[BODY_POINT_TRIANGLE];
			stiffness = (*cloth)[i].collision_stiffness[BODY_POINT_TRIANGLE];
		}
		else {
			obj_No = i - tetrahedron_begin_obj_index;
			vertex_mesh = &(*tetrahedron)[obj_No].mesh_struct;
			radius1 = (*tetrahedron)[obj_No].PC_radius[BODY_POINT_TRIANGLE];
			stiffness = (*tetrahedron)[obj_No].collision_stiffness[BODY_POINT_TRIANGLE];
		}
		neighbor_vertex = &(*triangle_neighbor_vertex)[i];
		collide_vertex = &(*collide_triangle_vertex)[i];
		collide_vertex->clear();
		collide_vertex->reserve(neighbor_vertex->size());

		for (int k = 0; k < neighbor_vertex->size(); ++k) {
			if (checkPointColliderTriangleCollision(vertex_mesh, triangle_mesh, radius0 + radius1, (*neighbor_vertex)[k], triangle_index, i,
				collider_No, target_pos, true, stiffness)) {
				collide_vertex->push_back((*neighbor_vertex)[k]);
			}
		}
	}
}







void Collision::pointColliderTriangleCollisionDetection(int thread_No, int vertex_index, int cloth_No,
	std::vector<int>* vertex_neighbor_triangle, std::vector<int>* collide_vertex_triangle, MeshStruct* vertex_mesh, double radius0,
	TargetPosition* target_pos, double vertex_collision_stiffness)
{
	MeshStruct* triangle_mesh;
	std::vector<int>* neighbor_triangle;
	std::vector<int>* collide_triangle;
	double radius1;
	for (int i = 0; i < collider->size(); ++i) {
		triangle_mesh = &(*collider)[i].mesh_struct;
		neighbor_triangle = &vertex_neighbor_triangle[i];
		collide_triangle = &collide_vertex_triangle[i];
		collide_triangle->clear();
		collide_triangle->reserve(neighbor_triangle->size());
		radius1 = (*collider)[i].tolerance;
		for (int k = 0; k < neighbor_triangle->size(); ++k) {
			if (checkPointColliderTriangleCollision(vertex_mesh, triangle_mesh, radius0 + radius1, vertex_index, (*neighbor_triangle)[k], cloth_No,
				i, target_pos, true, vertex_collision_stiffness)) {
				collide_triangle->push_back((*neighbor_triangle)[k]);
			}
		}
	}
}



void Collision::pointSelfTriangleClose(std::vector<int>* vertex_neighbor_triangle, double* initial_vertex_pos, double* current_vertex_pos,
	int vertex_index, int cloth_No, double mass, TargetPosition* target_position)
{
	TriangleMeshStruct* triangle_mesh;
	std::vector<int>* neighbor_triangle;
	int* triangle_vertex_index;
	std::array<double, 3>* current_position;
	std::array<double, 3>* initial_position;
	std::array<int, 3>* triangle_indices;
	std::array<double, 3>* initial_face_normal;
	int triangle_index;

	std::vector<double*>triangle_initial_pos(3);
	std::vector<double*>triangle_current_pos(3);

	double target_pos[3];
	std::vector<std::array<double, 3>> triangle_target_pos(3);
	double stiffness;
	double triangle_mass[3];
	double record_stiffness;
	double* vertex_b_sum = target_position->b_sum[cloth_No][vertex_index].data();
	double* vertex_stiffness = &target_position->stiffness[cloth_No][vertex_index];
	bool* vetex_need_update = &target_position->need_update[cloth_No][vertex_index];

	for (int i = 0; i < cloth->size(); ++i) {
		triangle_mesh = &(*cloth)[i].mesh_struct;
		neighbor_triangle = &(vertex_neighbor_triangle[i]);
		current_position = triangle_mesh->vertex_position.data();
		initial_position = triangle_mesh->vertex_for_render.data();
		triangle_indices = triangle_mesh->triangle_indices.data();
		initial_face_normal = triangle_mesh->face_normal_for_render.data();
		record_stiffness = (*cloth)[i].collision_stiffness[SELF_POINT_TRIANGLE];
		for (int k = 0; k < neighbor_triangle->size(); ++k) {
			triangle_index = (*neighbor_triangle)[k];
			triangle_vertex_index = triangle_indices[triangle_index].data();

			for (int j = 0; j < 3; ++j) {
				triangle_initial_pos[j] = initial_position[triangle_vertex_index[j]].data();
				triangle_current_pos[j] = current_position[triangle_vertex_index[j]].data();
				triangle_mass[j] = triangle_mesh->mass[triangle_vertex_index[j]];
			}
			stiffness = record_stiffness;
			if (collision_constraint.pointSelfTriangle(initial_vertex_pos, current_vertex_pos, triangle_initial_pos, triangle_current_pos,
				initial_face_normal[triangle_index].data(), target_pos, triangle_target_pos.data(), d_hat, stiffness, mass, triangle_mass)) {
				addTargetPosToSystemTotal(vertex_b_sum, target_position->collision_energy, initial_vertex_pos,
					target_pos, stiffness, *vertex_stiffness, *vetex_need_update);
				int triangle_vertex;
				for (int j = 0; j < 3; ++j) {
					triangle_vertex = triangle_vertex_index[j];
					addTargetPosToSystemTotal(target_position->b_sum[i][triangle_vertex].data(), target_position->collision_energy, initial_position[triangle_vertex].data(),
						triangle_target_pos[j].data(), stiffness,
						target_position->stiffness[i][triangle_vertex], target_position->need_update[i][triangle_vertex]);
					//if (vertex_index == 1) {
					//	//////std::cout << triangle_initial_pos[j][0] << " " << triangle_initial_pos[j][1] << " " << triangle_initial_pos[j][2] << std::endl;
					//	//////std::cout << triangle_current_pos[j][0] << " " << triangle_current_pos[j][1] << " " << triangle_current_pos[j][2] << std::endl;
					//	//////std::cout << triangle_target_pos[j][0] << " " << triangle_target_pos[j][1] << " " << triangle_target_pos[j][2] << std::endl;
					//	//////std::cout << initial_face_normal[triangle_index][0] << " " << initial_face_normal[triangle_index][1] << " " << initial_face_normal[triangle_index][2] << std::endl;
					//}
					//
					////////std::cout << target_position->collision_energy[i] << std::endl;
				}

			}
		}
	}
}

void Collision::edgeEdgeClose(std::vector<int>* edge_neighbor_edge, double* initial_edge_vertex_0, double* initial_edge_vertex_1, double* current_edge_vertex_0, double* current_edge_vertex_1,
	int cloth_No, int edge_vertex_index_0, int edge_vertex_index_1, double* mass, TargetPosition* target_position)
{
	TriangleMeshStruct* compare_mesh;
	std::vector<int>* neighbor_edge;
	double current_collision_time;
	std::array<double, 3>* current_position;
	std::array<double, 3>* initial_position;
	std::vector<std::array<double, 3>> target_pos(2);
	std::vector<std::array<double, 3>> compare_target_pos(2);
	int compare_edge_index_0;
	int compare_edge_index_1;
	double stiffness;
	double record_stiffness;
	for (int i = 0; i < cloth->size(); ++i) {
		compare_mesh = &(*cloth)[i].mesh_struct;
		neighbor_edge = &(edge_neighbor_edge[i]);
		current_position = compare_mesh->vertex_position.data();
		initial_position = compare_mesh->vertex_for_render.data();
		record_stiffness = (*cloth)[i].collision_stiffness[SELF_EDGE_EDGE];
		for (int k = 0; k < neighbor_edge->size(); ++k) {
			compare_edge_index_0 = compare_mesh->edge_vertices[(*neighbor_edge)[k] << 1];//   edges[(*neighbor_edge)[k]].vertex[0];
			compare_edge_index_1 = compare_mesh->edge_vertices[((*neighbor_edge)[k] << 1) + 1];// edges[(*neighbor_edge)[k]].vertex[1];
			mass[2] = compare_mesh->mass[compare_edge_index_0];
			mass[3] = compare_mesh->mass[compare_edge_index_1];
			stiffness = record_stiffness;
			if (collision_constraint.edgeEdgeCollision(target_pos, compare_target_pos, current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0,
				initial_edge_vertex_1, current_position[compare_edge_index_0].data(), current_position[compare_edge_index_1].data(),
				initial_position[compare_edge_index_0].data(), initial_position[compare_edge_index_1].data(), mass, d_hat, stiffness)) {
				addTargetPosToSystemTotal(target_position->b_sum[cloth_No][edge_vertex_index_0].data(),
					target_position->collision_energy, initial_edge_vertex_0,
					target_pos[0].data(), stiffness, target_position->stiffness[cloth_No][edge_vertex_index_0], target_position->need_update[cloth_No][edge_vertex_index_0]);
				addTargetPosToSystemTotal(target_position->b_sum[cloth_No][edge_vertex_index_1].data(),
					target_position->collision_energy, initial_edge_vertex_1,
					target_pos[1].data(), stiffness, target_position->stiffness[cloth_No][edge_vertex_index_1], target_position->need_update[cloth_No][edge_vertex_index_1]);
				addTargetPosToSystemTotal(target_position->b_sum[i][compare_edge_index_0].data(),
					target_position->collision_energy, initial_position[compare_edge_index_0].data(),
					compare_target_pos[0].data(), stiffness, target_position->stiffness[i][compare_edge_index_0], target_position->need_update[i][compare_edge_index_0]);
				addTargetPosToSystemTotal(target_position->b_sum[i][compare_edge_index_1].data(),
					target_position->collision_energy, initial_position[compare_edge_index_1].data(),
					compare_target_pos[1].data(), stiffness, target_position->stiffness[i][compare_edge_index_1], target_position->need_update[i][compare_edge_index_1]);
			}
		}
	}
}



void Collision::pointColliderTriangleClose(int* triangle_vertex_index, std::vector<int>* triangle_neighbor_vertex,
	std::vector<double*>& current_position, double* current_face_normal, TargetPosition* target_position)
{
	MeshStruct* vertex_mesh;
	std::vector<int>* neighbor_vertex;
	double current_collision_time;
	double target_pos[3];
	double stiffness;
	double record_stiffness;
	std::array<double, 3>* vertex_b_sum;
	double* vertex_stiffness;
	bool* vertex_need_update;
	int vertex_index;
	double* energy;
	for (int i = 0; i < cloth->size(); ++i) {
		vertex_mesh = &(*cloth)[i].mesh_struct;
		neighbor_vertex = &(triangle_neighbor_vertex[i]);
		vertex_b_sum = target_position->b_sum[i].data();
		vertex_stiffness = target_position->stiffness[i].data();
		vertex_need_update = target_position->need_update[i];
		record_stiffness = (*cloth)[i].collision_stiffness[BODY_POINT_TRIANGLE];
		energy = &target_position->collision_energy;
		for (int k = 0; k < neighbor_vertex->size(); ++k) {
			vertex_index = (*neighbor_vertex)[k];
			stiffness = record_stiffness;
			if (collision_constraint.pointColliderTriangle(vertex_mesh->vertex_for_render[vertex_index].data(), vertex_mesh->vertex_position[vertex_index].data(),
				current_position, current_face_normal, target_pos, d_hat, stiffness)) {
				addTargetPosToSystemTotal(vertex_b_sum[vertex_index].data(),
					*energy, vertex_mesh->vertex_for_render[vertex_index].data(),
					target_pos, stiffness,
					vertex_stiffness[vertex_index], vertex_need_update[vertex_index]);
				////////std::cout << vertex_index << " " << stiffness << std::endl;
			}
		}
	}
}



void Collision::pointSelfTriangleCollisionTime(double* collision_time, std::vector<int>* vertex_neighbor_triangle, double* initial_vertex_pos, double* current_vertex_pos, int vertex_index)
{
	MeshStruct* triangle_mesh;
	std::vector<int>* neighbor_triangle;
	double current_collision_time;
	int* triangle_vertex_index;
	std::array<double, 3>* current_position;
	std::array<double, 3>* initial_position;
	std::array<int, 3>* triangle_indices;
	std::array<double, 3>* current_ori_face_normal;
	std::array<double, 3>* initial_ori_face_normal;
	std::array<double, 3>* cross_for_CCD;
	std::array<floating, 3>* f_current_face_normal;
	std::array<floating, 3>* f_initial_face_normal;
	std::array<floating, 3>* f_cross_for_CCD;

	int triangle_index;
	for (int i = 0; i < cloth->size(); ++i) {
		triangle_mesh = &(*cloth)[i].mesh_struct;
		neighbor_triangle = &(vertex_neighbor_triangle[i]);
		current_position = triangle_mesh->vertex_position.data();
		initial_position = triangle_mesh->vertex_for_render.data();
		triangle_indices = triangle_mesh->triangle_indices.data();
		current_ori_face_normal = triangle_mesh->ori_face_normal.data();
		initial_ori_face_normal = triangle_mesh->ori_face_normal_for_render.data();
		cross_for_CCD = triangle_mesh->cross_for_approx_CCD.data();

		f_current_face_normal = triangle_mesh->f_face_normal.data();
		f_initial_face_normal = triangle_mesh->f_face_normal_for_render.data();
		f_cross_for_CCD = triangle_mesh->f_cross_for_approx_CCD.data();

		for (int k = 0; k < neighbor_triangle->size(); ++k) {
			triangle_index = (*neighbor_triangle)[k];
			triangle_vertex_index = triangle_indices[triangle_index].data();

			current_collision_time = CCD::pointTriangleCcd(initial_vertex_pos, initial_position[triangle_vertex_index[0]].data(), initial_position[triangle_vertex_index[1]].data(), initial_position[triangle_vertex_index[2]].data(),
				current_vertex_pos, current_position[triangle_vertex_index[0]].data(), current_position[triangle_vertex_index[1]].data(), current_position[triangle_vertex_index[2]].data(), eta, tolerance);
			if ((*collision_time) > current_collision_time) {
				(*collision_time) = current_collision_time;
			}
			//if (approx_CCD.pointTriangleCollisionTime(current_collision_time, initial_vertex_pos, current_vertex_pos, initial_position[triangle_vertex_index[0]].data(), current_position[triangle_vertex_index[0]].data(),
			//	initial_position[triangle_vertex_index[1]].data(), current_position[triangle_vertex_index[1]].data(), initial_position[triangle_vertex_index[2]].data(), current_position[triangle_vertex_index[2]].data(),
			//	initial_ori_face_normal[triangle_index].data(),
			//	current_ori_face_normal[triangle_index].data(), cross_for_CCD[triangle_index].data(), tolerance_2,
			//	f_initial_face_normal[triangle_index].data(), f_current_face_normal[triangle_index].data(),f_cross_for_CCD[triangle_index].data(),vertex_index)) {//
			//	//if (current_collision_time < 1e-4) {
			//	//	////std::cout << current_collision_time << " " << vertex_index << " " << triangle_vertex_index[0] << " " << triangle_vertex_index[1] << " " << triangle_vertex_index[2] << std::endl;
			//	//}
			//	if ((*collision_time) > current_collision_time) {
			//		(*collision_time) = current_collision_time;
			//	}
			//}
		}
	}
}

void Collision::pointColliderTriangleCollisionTime(double* collision_time, int* triangle_vertex_index, std::vector<int>* triangle_neighbor_vertex,
	std::array<double, 3>* initial_position, std::array<double, 3>* current_position,
	double* initial_ori_face_normal, double* current_ori_face_normal, double* cross_for_CCD,
	floating* f_initial_normal, floating* f_current_normal, floating* f_cross_for_CCD, int triangle_index)
{
	MeshStruct* vertex_mesh;
	std::vector<int>* neighbor_vertex;
	double current_collision_time;
	for (int i = 0; i < cloth->size(); ++i) {
		vertex_mesh = &(*cloth)[i].mesh_struct;
		neighbor_vertex = &(triangle_neighbor_vertex[i]);
		//if (*time_stamp == 15 && triangle_index == 9104)
		//{
		//	////std::cout << "k11" << std::endl;
		//}
		for (int k = 0; k < neighbor_vertex->size(); ++k) {
			//if (*time_stamp == 15 && triangle_index == 9104 && k==3)
			//{
			//	////std::cout << vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]][0]<<" "<< vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]][1]<<" "
			//		<< vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]][2] << std::endl;
			//	////std::cout << vertex_mesh->vertex_position[(*neighbor_vertex)[k]][0] << " " << vertex_mesh->vertex_position[(*neighbor_vertex)[k]][1] << " "
			//		<< vertex_mesh->vertex_position[(*neighbor_vertex)[k]][2] << std::endl;
			//	////std::cout << initial_position[triangle_vertex_index[0]][0] << " " << initial_position[triangle_vertex_index[0]][1] << " " <<
			//		initial_position[triangle_vertex_index[0]][2] << std::endl;
			//	////std::cout << initial_position[triangle_vertex_index[1]][0] << " " << initial_position[triangle_vertex_index[1]][1] << " " <<
			//		initial_position[triangle_vertex_index[1]][2] << std::endl;
			//	////std::cout << initial_position[triangle_vertex_index[2]][0] << " " << initial_position[triangle_vertex_index[2]][1] << " " <<
			//		initial_position[triangle_vertex_index[2]][2] << std::endl;
			//	////std::cout << current_position[triangle_vertex_index[0]][0] << " " << current_position[triangle_vertex_index[0]][1] << " " <<
			//		current_position[triangle_vertex_index[0]][2] << std::endl;
			//	////std::cout << current_position[triangle_vertex_index[1]][0] << " " << current_position[triangle_vertex_index[1]][1] << " " <<
			//		current_position[triangle_vertex_index[1]][2] << std::endl;
			//	////std::cout << current_position[triangle_vertex_index[2]][0] << " " << current_position[triangle_vertex_index[2]][1] << " " <<
			//		current_position[triangle_vertex_index[2]][2] << std::endl;
			//}
			current_collision_time = CCD::pointTriangleCcd(vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]].data(), initial_position[triangle_vertex_index[0]].data(), initial_position[triangle_vertex_index[1]].data(), initial_position[triangle_vertex_index[2]].data(),
				vertex_mesh->vertex_position[(*neighbor_vertex)[k]].data(), current_position[triangle_vertex_index[0]].data(), current_position[triangle_vertex_index[1]].data(), current_position[triangle_vertex_index[2]].data(), eta, tolerance);
			if ((*collision_time) > current_collision_time) {
				(*collision_time) = current_collision_time;
			}
			//if (*time_stamp == 15 && triangle_index == 9104)
			//{
			//	////std::cout << k << " " << (*neighbor_vertex)[k] << std::endl;
			//}
			/*if (approx_CCD.pointTriangleCollisionTime(current_collision_time, vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]].data(), vertex_mesh->vertex_position[(*neighbor_vertex)[k]].data(), initial_position[triangle_vertex_index[0]].data(), current_position[triangle_vertex_index[0]].data(),
				initial_position[triangle_vertex_index[1]].data(), current_position[triangle_vertex_index[1]].data(),
				initial_position[triangle_vertex_index[2]].data(), current_position[triangle_vertex_index[2]].data(),
				initial_ori_face_normal, current_ori_face_normal, cross_for_CCD, tolerance_2,
				f_initial_normal, f_current_normal, f_cross_for_CCD, (*neighbor_vertex)[k])) {
				if ((*collision_time) > current_collision_time) {
					(*collision_time) = current_collision_time;
				}
			}*/
		}
		//if (*time_stamp == 15 && triangle_index == 9104)
		//{
		//	////std::cout << "k22" << std::endl;
		//}
	}
}

void Collision::edgeEdgeCollisionTime(double* collision_time, std::vector<int>* edge_neighbor_edge, double* initial_edge_vertex_0, double* initial_edge_vertex_1, double* current_edge_vertex_0, double* current_edge_vertex_1)
{
	TriangleMeshStruct* compare_mesh;
	std::vector<int>* neighbor_edge;
	double current_collision_time;
	std::array<double, 3>* current_position;
	std::array<double, 3>* initial_position;
	unsigned int* edge_vertex_index;
	for (int i = 0; i < cloth->size(); ++i) {
		compare_mesh = &(*cloth)[i].mesh_struct;
		neighbor_edge = &(edge_neighbor_edge[i]);
		current_position = compare_mesh->vertex_position.data();
		initial_position = compare_mesh->vertex_for_render.data();
		for (int k = 0; k < neighbor_edge->size(); ++k) {
			edge_vertex_index = compare_mesh->edge_vertices.data() + ((*neighbor_edge)[k] << 1);//edges[(*neighbor_edge)[k]].vertex;

			current_collision_time = CCD::edgeEdgeCcd(initial_edge_vertex_0, initial_edge_vertex_1, initial_position[edge_vertex_index[0]].data(), initial_position[edge_vertex_index[1]].data(),
				current_edge_vertex_0, current_edge_vertex_1, current_position[edge_vertex_index[0]].data(), current_position[edge_vertex_index[1]].data(), eta, tolerance);
			
			
			if ((*collision_time) > current_collision_time) {
				(*collision_time) = current_collision_time;
			}

			/*if (approx_CCD.edgeEdgeCollisionTime(current_collision_time, current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0, initial_edge_vertex_1, current_position[edge_vertex_index[0]].data(), current_position[edge_vertex_index[1]].data(),
				initial_position[edge_vertex_index[0]].data(), initial_position[edge_vertex_index[1]].data(), tolerance_2)) {
				if ((*collision_time) > current_collision_time) {
					(*collision_time) = current_collision_time;
				}
			}*/
		}
	}
}


void Collision::pointSelfTriangleCollisionDetection(int thread_No, int vertex_index, int cloth_No,
	std::vector<int>* vertex_neighbor_triangle, std::vector<int>* collide_vertex_triangle, MeshStruct* vertex_mesh, double radius0,
	TargetPosition* target_pos, double vertex_collision_stiffness)
{
	MeshStruct* triangle_mesh;
	std::vector<int>* neighbor_triangle;
	std::vector<int>* collide_triangle;
	double radius1;
	double triangle_collision_stiffness;

	int tet_No;

	for (int i = 0; i < total_obj_num; ++i) {
		if (i < tetrahedron_begin_obj_index) {
			triangle_mesh = &(*cloth)[i].mesh_struct;
			radius1 = (*cloth)[i].PC_radius[SELF_POINT_TRIANGLE];
			triangle_collision_stiffness = (*cloth)[i].collision_stiffness[SELF_POINT_TRIANGLE];
		}
		else {
			tet_No = i - tetrahedron_begin_obj_index;
			triangle_mesh = &(*tetrahedron)[tet_No].mesh_struct;
			radius1 = (*tetrahedron)[tet_No].PC_radius[SELF_POINT_TRIANGLE];
			triangle_collision_stiffness = (*tetrahedron)[tet_No].collision_stiffness[SELF_POINT_TRIANGLE];
		}
		neighbor_triangle = &vertex_neighbor_triangle[i];
		collide_triangle = &collide_vertex_triangle[i];
		collide_triangle->clear();
		collide_triangle->reserve(neighbor_triangle->size());
		for (int k = 0; k < neighbor_triangle->size(); ++k) {
			if (checkPointTriangleCollision(vertex_mesh, triangle_mesh, radius0 + radius1, vertex_index, (*neighbor_triangle)[k], cloth_No,
				i, target_pos, true, vertex_collision_stiffness, triangle_collision_stiffness)) {
				collide_triangle->push_back((*neighbor_triangle)[k]);
			}
		}
	}
}

void Collision::pointColliderTriangleCollisionReDetection(int thread_No, int vertex_index, int cloth_No, std::vector<std::vector<int>>* collide_vertex_triangle, MeshStruct* vertex_mesh,
	double radius0, std::vector<double>* collision_stiffness, TargetPosition* target_postion_)
{
	MeshStruct* triangle_mesh;
	std::vector<int>* collide_triangle;
	double radius1;
	for (int i = 0; i < collider->size(); ++i) {
		//triangle_mesh = &(*collider)[i].mesh_struct;
		//collide_triangle = &(*collide_vertex_triangle)[i];
		//radius1 = (*collider)[i].tolerance;
		//for (int k = 0; k < collide_triangle->size(); ++k) {
		//	if (!checkPointColliderTriangleCollision(vertex_mesh, triangle_mesh, radius0 + radius1, vertex_index, (*collide_triangle)[k], cloth_No,
		//		i, target_postion_, false, (*collision_stiffness))) {
		//		addTargetPosToSystem(target_postion_->b_sum[cloth_No][vertex_index].data(),
		//			target_postion_->collision_energy[cloth_No], vertex_mesh->vertex_position[vertex_index].data(), vertex_mesh->vertex_position[vertex_index].data(), (*collision_stiffness)[vertex_index]);
		//	}
		//}
	}
}

void Collision::pointSelfTriangleCollisionReDetection(int thread_No, int vertex_index, int cloth_No, std::vector<int>* collide_vertex_triangle, MeshStruct* vertex_mesh,
	double radius0, double collision_stiffness, TargetPosition* target_postion_)
{
	MeshStruct* triangle_mesh;
	std::vector<int>* collide_triangle;
	double radius1;
	int* triangle_vertex_index;
	int triangle_vertex;
	std::vector<double>* triangle_collision_stiffness;
	for (int i = 0; i < cloth->size(); ++i) {
		//triangle_mesh = &(*cloth)[i].mesh_struct;
		//collide_triangle = &collide_vertex_triangle[i];
		//radius1 = (*cloth)[i].PC_radius[SELF_POINT_TRIANGLE]+ radius0;
		//triangle_collision_stiffness = &(*cloth)[i].collision_stiffness[SELF_POINT_TRIANGLE];
		//for (int k = 0; k < collide_triangle->size(); ++k) {
		//	if (!checkPointTriangleCollision(vertex_mesh, triangle_mesh, radius1, vertex_index, (*collide_triangle)[k], cloth_No,
		//		i, target_postion_, false, (*collision_stiffness), (*triangle_collision_stiffness))) {
		//		addTargetPosToSystem(target_postion_->b_sum[cloth_No][vertex_index].data(),
		//			target_postion_->collision_energy[cloth_No], vertex_mesh->vertex_position[vertex_index].data(), vertex_mesh->vertex_position[vertex_index].data(), (*collision_stiffness)[vertex_index]);
		//		
		//		triangle_vertex_index = triangle_mesh->triangle_indices[(*collide_triangle)[k]].data();
		//		for (int j = 0; j < 3; ++j) {
		//			triangle_vertex = triangle_vertex_index[j];
		//			addTargetPosToSystem(target_postion_->b_sum[i][triangle_vertex].data(),
		//				target_postion_->collision_energy[i], triangle_mesh->vertex_position[triangle_vertex].data(), triangle_mesh->vertex_position[triangle_vertex].data(), (*triangle_collision_stiffness)[triangle_vertex]);
		//		}
		//	}
		//}
	}
}

void Collision::colliderTriangleVertexCollisionReDetection(int thread_No, int triangle_index, int collider_No,
	std::vector<int>* collide_triangle_vertex, MeshStruct* triangle_mesh, double radius0, TargetPosition* target_pos)
{
	MeshStruct* vertex_mesh;
	std::vector<int>* collide_vertex;
	double radius1;
	std::vector<double>* collision_stiffness;
	for (int i = 0; i < cloth->size(); ++i) {
		//vertex_mesh = &(*cloth)[i].mesh_struct;
		//collide_vertex = &collide_triangle_vertex[i];
		//radius1 = (*cloth)[i].PC_radius[BODY_POINT_TRIANGLE]+radius0;
		//collision_stiffness = &(*cloth)[i].collision_stiffness[BODY_POINT_TRIANGLE];
		//for (int k = 0; k < collide_vertex->size(); ++k) {
		//	if (!checkPointColliderTriangleCollision(vertex_mesh, triangle_mesh, radius1, (*collide_vertex)[k], triangle_index, i,
		//		collider_No, target_pos, false, (*collision_stiffness))) {
		//		addTargetPosToSystem(target_pos->b_sum[i][(*collide_vertex)[k]].data(),
		//			target_pos->collision_energy[i], vertex_mesh->vertex_position[(*collide_vertex)[k]].data(), vertex_mesh->vertex_position[(*collide_vertex)[k]].data(), (*collision_stiffness)[(*collide_vertex)[k]]);
		//	}
		//}
	}
}


void Collision::edgeSelfEdgeCollisionReDetection(int thread_No, int edge_index, int cloth_No, std::vector<std::vector<int>>* collide_edge_edge, TriangleMeshStruct* edge_mesh,
	double radius0, std::vector<double>& collision_stiffness, TargetPosition* target_postion_)
{
	MeshStruct* compare_edge_mesh;
	std::vector<int>* collide_edge;
	double radius1;
	std::vector<double>* compare_collision_stiffness;
	unsigned int* edge_vertex = edge_mesh->edge_vertices.data() + (edge_index << 1);// edges[edge_index].vertex;
	int* compare_edge_index;
	for (int i = 0; i < cloth->size(); ++i) {
		//compare_edge_mesh = &(*cloth)[i].mesh_struct;
		//collide_edge = &(*collide_edge_edge)[i];
		//radius1 = (*cloth)[i].PC_radius[SELF_EDGE_EDGE];
		//compare_collision_stiffness= &(*cloth)[i].collision_stiffness[SELF_EDGE_EDGE];
		//for (int k = 0; k < collide_edge->size(); ++k) {
		//	compare_edge_index = compare_edge_mesh->edges[(*collide_edge)[k]].vertex;
		//	if (!checkEdgeEdgeCollision(edge_mesh, compare_edge_mesh, radius0 + radius1, edge_index, (*collide_edge)[k], cloth_No,
		//		i, target_postion_, false, collision_stiffness, (*compare_collision_stiffness))) {
		//		for (int j = 0; j < 2; ++j) {
		//			addTargetPosToSystem(target_postion_->b_sum[cloth_No][edge_vertex[j]].data(),target_postion_->collision_energy[cloth_No], 
		//				edge_mesh->vertex_position[edge_vertex[j]].data(), edge_mesh->vertex_position[edge_vertex[j]].data(), collision_stiffness[edge_vertex[j]]);
		//			addTargetPosToSystem(target_postion_->b_sum[i][compare_edge_index[j]].data(),target_postion_->collision_energy[i],
		//				compare_edge_mesh->vertex_position[compare_edge_index[j]].data(), compare_edge_mesh->vertex_position[compare_edge_index[j]].data(), (*compare_collision_stiffness)[compare_edge_index[j]]);
		//		}
		//	}
		//}
	}
}

void Collision::edgeSelfEdgeCollisionDetection(int thread_No, int edge_index, int cloth_No,
	std::vector<int>* edge_neighbor_edge, std::vector<int>* collide_edge_edge, MeshStruct* edge_mesh, double radius0,
	TargetPosition* target_pos, double stiffness_0)
{
	MeshStruct* compare_edge_mesh;
	std::vector<int>* neighbor_edge;
	std::vector<int>* collide_edge;
	double radius1;
	double stiffness;
	for (int i = 0; i < total_obj_num; ++i) {
		if (i < tetrahedron_begin_obj_index) {
			compare_edge_mesh = &(*cloth)[i].mesh_struct;
			stiffness = (*cloth)[i].collision_stiffness[SELF_EDGE_EDGE];
			radius1 = (*cloth)[i].PC_radius[SELF_EDGE_EDGE];
		}
		else {
			compare_edge_mesh = &(*tetrahedron)[i - tetrahedron_begin_obj_index].mesh_struct;
			stiffness = (*tetrahedron)[i - tetrahedron_begin_obj_index].collision_stiffness[SELF_EDGE_EDGE];
			radius1 = (*tetrahedron)[i - tetrahedron_begin_obj_index].PC_radius[SELF_EDGE_EDGE];
		}
		neighbor_edge = &edge_neighbor_edge[i];
		collide_edge = &collide_edge_edge[i];
		collide_edge->clear();
		collide_edge->reserve(neighbor_edge->size());
		for (int k = 0; k < neighbor_edge->size(); ++k) {
			if (checkEdgeEdgeCollision(edge_mesh, compare_edge_mesh, radius0 + radius1, edge_index, (*neighbor_edge)[k], cloth_No,
				i, target_pos, true, stiffness_0, stiffness)) {
				collide_edge->push_back((*neighbor_edge)[k]);
			}
		}
	}
}


bool Collision::checkPointTriangleCollision(MeshStruct* vertex_mesh, MeshStruct* triangle_mesh,
	double radius, int vertex_index, int triangle_index, int vertex_cloth_No, int triangle_cloth_No, TargetPosition* target_position, bool new_collision_registration,
	double vertex_collision_stiffness, double triangle_collision_stiffness)
{
	std::vector<double*>initial_triangle_pos(3); std::vector<double*>current_triangle_pos(3);
	double triangle_mass[3];
	int* triangle_vertex_index = triangle_mesh->triangle_indices[triangle_index].data();
	for (int i = 0; i < 3; ++i) {
		initial_triangle_pos[i] = triangle_mesh->vertex_for_render[triangle_vertex_index[i]].data();
		current_triangle_pos[i] = triangle_mesh->vertex_position[triangle_vertex_index[i]].data();
		triangle_mass[i] = triangle_mesh->mass[triangle_vertex_index[i]];
	}
	double vertex_target_pos[3];
	std::vector<std::array<double, 3>> triangle_target_pos;
	if (predictive_contact.pointTriangleCollision(vertex_mesh->vertex_for_render[vertex_index].data(), vertex_mesh->vertex_position[vertex_index].data(), initial_triangle_pos, current_triangle_pos,
		triangle_mesh->face_normal_for_render[triangle_index].data(), triangle_mesh->face_normal[triangle_index].data(), vertex_target_pos, triangle_target_pos, radius, triangle_mesh->triangle_normal_magnitude_reciprocal[triangle_index],
		vertex_mesh->mass[vertex_index], triangle_mass)) {
		//if (vertex_index == 0) {
		//	//////std::cout << "++" << std::endl;
		//	//////std::cout << vertex_index << " " << triangle_index << " " << vertex_mesh->vertex_for_render[vertex_index][0] << " " << vertex_mesh->vertex_for_render[vertex_index][1] << " " << vertex_mesh->vertex_for_render[vertex_index][2] << " " << std::endl;
		//	for (int i = 0; i < triangle_target_pos.size(); ++i) {
		//		//////std::cout << vertex_index << " " << triangle_index << " " << initial_triangle_pos[i][0] << " " << initial_triangle_pos[i][1] << " " << initial_triangle_pos[i][2] << " " << std::endl;
		//	}
		//	//////std::cout << vertex_mesh->vertex_position[vertex_index][0] << " " << vertex_mesh->vertex_position[vertex_index][1] << " " << vertex_mesh->vertex_position[vertex_index][2] << " " << std::endl;
		//	for (int i = 0; i < triangle_target_pos.size(); ++i) {
		//		//////std::cout << current_triangle_pos[i][0] << " " << current_triangle_pos[i][1] << " " << current_triangle_pos[i][2] << " " << std::endl;
		//	}
		//	//////std::cout << vertex_index << " " << triangle_index <<" "<< vertex_target_pos[0] << " " << vertex_target_pos[1] << " " << vertex_target_pos[2] << " " << std::endl;
		//	for (int i = 0; i < triangle_target_pos.size(); ++i) {
		//		//////std::cout << vertex_index<<" "<< triangle_index << " " << triangle_target_pos[i][0] << " " << triangle_target_pos[i][1] << " " << triangle_target_pos[i][2] << " " << std::endl;
		//	}
		//	//////std::cout << "++" << std::endl;
		//}

		addTargetPosToSystem(target_position->b_sum[vertex_cloth_No][vertex_index].data(),
			target_position->collision_energy, vertex_mesh->vertex_position[vertex_index].data(), vertex_target_pos, vertex_collision_stiffness);
		int triangle_vertex;
		for (int i = 0; i < 3; ++i) {
			triangle_vertex = triangle_vertex_index[i];
			addTargetPosToSystem(target_position->b_sum[triangle_cloth_No][triangle_vertex].data(),
				target_position->collision_energy, triangle_mesh->vertex_position[triangle_vertex].data(), triangle_target_pos[i].data(), triangle_collision_stiffness);
		}
		if (new_collision_registration) {
			target_position->stiffness[vertex_cloth_No][vertex_index] += vertex_collision_stiffness;
			target_position->need_update[vertex_cloth_No][vertex_index] = true;
			for (int i = 0; i < 3; ++i) {
				target_position->stiffness[triangle_cloth_No][triangle_vertex_index[i]] += triangle_collision_stiffness;
				target_position->need_update[triangle_cloth_No][triangle_vertex_index[i]] = true;
			}
		}

		return true;
	}
	return false;
}

bool Collision::checkEdgeEdgeCollision(MeshStruct* edge_mesh, MeshStruct* compare_mesh,
	double radius, int edge_index, int compare_edge_index, int edge_cloth_No, int compare_cloth_No, TargetPosition* target_position, bool new_collision_registration,
	double collision_stiffness, double compare_collision_stiffness)
{
	double mass[4];
	unsigned int* edge_vertex_index = edge_mesh->edge_vertices.data() + (edge_index << 1);// edges[edge_index].vertex;
	unsigned int* compare_edge_vertex_index = compare_mesh->edge_vertices.data() + (compare_edge_index << 1);// edges[compare_edge_index].vertex;
	for (int i = 0; i < 2; ++i) {
		mass[i] = edge_mesh->mass[edge_vertex_index[i]];
		mass[2 + i] = edge_mesh->mass[compare_edge_vertex_index[i]];
	}
	std::vector<std::array<double, 3>> target_pos;
	std::vector<std::array<double, 3>> compare_target_pos;
	if (predictive_contact.edgeEdgeCollision(target_pos, compare_target_pos, radius, edge_mesh->vertex_position[edge_vertex_index[0]].data(), edge_mesh->vertex_position[edge_vertex_index[1]].data(),
		edge_mesh->vertex_for_render[edge_vertex_index[0]].data(), edge_mesh->vertex_for_render[edge_vertex_index[1]].data(), compare_mesh->vertex_position[compare_edge_vertex_index[0]].data(),
		compare_mesh->vertex_position[compare_edge_vertex_index[1]].data(), compare_mesh->vertex_for_render[compare_edge_vertex_index[0]].data(), compare_mesh->vertex_for_render[compare_edge_vertex_index[1]].data(),
		mass)) {
		for (int i = 0; i < 2; ++i) {
			addTargetPosToSystem(target_position->b_sum[edge_cloth_No][edge_vertex_index[i]].data(),
				target_position->collision_energy, edge_mesh->vertex_position[edge_vertex_index[i]].data(), target_pos[i].data(), collision_stiffness);
			addTargetPosToSystem(target_position->b_sum[compare_cloth_No][compare_edge_vertex_index[i]].data(),
				target_position->collision_energy, compare_mesh->vertex_position[compare_edge_vertex_index[i]].data(), compare_target_pos[i].data(), compare_collision_stiffness);
		}
		if (new_collision_registration) {
			for (int i = 0; i < 2; ++i) {
				target_position->stiffness[edge_cloth_No][edge_vertex_index[i]] += collision_stiffness;
				target_position->need_update[edge_cloth_No][edge_vertex_index[i]] = true;
				target_position->stiffness[compare_cloth_No][compare_edge_vertex_index[i]] += compare_collision_stiffness;
				target_position->need_update[edge_cloth_No][compare_edge_vertex_index[i]] = true;
			}
		}
		return true;
	}
	return false;
}

bool Collision::checkPointColliderTriangleCollision(MeshStruct* vertex_mesh, MeshStruct* triangle_mesh,
	double radius, int vertex_index, int triangle_index, int vertex_cloth_No, int triangle_obj_No, TargetPosition* target_position, bool new_collision_registration,
	double vertex_collision_stiffness)
{
	std::vector<double*>initial_triangle_pos(3); std::vector<double*>current_triangle_pos(3);
	int* triangle_vertex_index = triangle_mesh->triangle_indices[triangle_index].data();
	for (int i = 0; i < 3; ++i) {
		initial_triangle_pos[i] = triangle_mesh->vertex_for_render[triangle_vertex_index[i]].data();
		current_triangle_pos[i] = triangle_mesh->vertex_position[triangle_vertex_index[i]].data();
	}
	double vertex_target_pos[3];
	if (predictive_contact.pointColliderTriangleCollision(vertex_mesh->vertex_for_render[vertex_index].data(), vertex_mesh->vertex_position[vertex_index].data(), initial_triangle_pos, current_triangle_pos,
		triangle_mesh->face_normal_for_render[triangle_index].data(), triangle_mesh->face_normal[triangle_index].data(), vertex_target_pos, radius)) {
		if (new_collision_registration) {
			target_position->stiffness[vertex_cloth_No][vertex_index] += vertex_collision_stiffness;
			target_position->need_update[vertex_cloth_No][vertex_index] = true;
			//if (vertex_index == 0) {
			//	//////std::cout << vertex_target_pos[0] << " " << vertex_target_pos[1] << " " << vertex_target_pos[2] << " "
			//		<< vertex_mesh->vertex_position[vertex_index][0] << " " << vertex_mesh->vertex_position[vertex_index][1] << " " << vertex_mesh->vertex_position[vertex_index][2] << std::endl;
			//	//////std::cout << current_triangle_pos[0][1] << std::endl;
			//}
		}
		else {
			/*if (vertex_index == 0) {
				//////std::cout << vertex_target_pos[0] << " " << vertex_target_pos[1] << " " << vertex_target_pos[2] << " "
					<< vertex_mesh->vertex_position[vertex_index][0] << " " << vertex_mesh->vertex_position[vertex_index][1] << " " << vertex_mesh->vertex_position[vertex_index][2] << std::endl;
			}*/
		}

		addTargetPosToSystem(target_position->b_sum[vertex_cloth_No][vertex_index].data(),
			target_position->collision_energy, vertex_mesh->vertex_position[vertex_index].data(), vertex_target_pos, vertex_collision_stiffness);
		return true;
	}
	return false;
}


void Collision::construct_b_sum(double* b_sum, double* target_pos, double stiffness)
{
	b_sum[0] += target_pos[0] * stiffness;
	b_sum[1] += target_pos[1] * stiffness;
	b_sum[2] += target_pos[2] * stiffness;
}

void Collision::addTargetPosToSystem(double* b_sum, double& energy, double* current_pos, double* target_pos, double stiffness)
{
	energy += 0.5 * stiffness * EDGE_LENGTH(current_pos, target_pos);
	b_sum[0] += target_pos[0] * stiffness;
	b_sum[1] += target_pos[1] * stiffness;
	b_sum[2] += target_pos[2] * stiffness;
}

void Collision::addTargetPosToSystemTotal(double* b_sum, double& energy, double* current_pos, double* target_pos, double stiffness, double& sum_stiffness, bool& update)
{
	energy += 0.5 * stiffness * EDGE_LENGTH(current_pos, target_pos);
	b_sum[0] += target_pos[0] * stiffness;
	b_sum[1] += target_pos[1] * stiffness;
	b_sum[2] += target_pos[2] * stiffness;
	update = true;
	sum_stiffness += stiffness;
}

void Collision::test()
{
	//findAllNeighborPairs();

	//////std::cout << "vertex " << std::endl;
	//for (int i = 0; i < (*cloth)[0].vertex_neighbor_obj_triangle[0][0].size(); ++i) {
	//	//////std::cout << (*cloth)[0].vertex_neighbor_obj_triangle[0][0][i] << " ";
	//}
	////////std::cout << std::endl;
	////////std::cout << "triangle " << std::endl;
	//for (int i = 0; i < (*cloth)[0].triangle_neighbor_cloth_triangle[0][0].size(); ++i) {
	//	//////std::cout << (*cloth)[0].triangle_neighbor_cloth_triangle[0][0][i] << " ";
	//}
	////////std::cout << std::endl;
	////////std::cout << "edge " << std::endl;
	//for (int i = 0; i < (*cloth)[0].edge_neighbor_obj_edge[0][0].size(); ++i) {
	//	//////std::cout << (*cloth)[0].edge_neighbor_obj_edge[0][0][i] << " ";
	//}
	//////std::cout << std::endl;
}


void Collision::searchPatch(double* aabb, unsigned int compare_index, unsigned int obj_No, bool& intersect)
{
	//for (int i = 0; i < total_obj_num; ++i) {
	//	if (obj_BVH[i].searchIfPatchIntersect(aabb, compare_index, i == obj_No, 1, 0, mesh_patch.patch_AABB[i].size())) {
	//		intersect = true;
	//	}	
	//}
}

void Collision::searchTriangle(double* aabb, unsigned int compare_index, unsigned int obj_No, std::vector<std::vector<unsigned int>>* obj_neighbor_index,
	std::vector<std::vector<unsigned int>>* collider_neighbor_index, bool is_collider)
{

	for (int i = 0; i < total_obj_num; ++i) {
		(*obj_neighbor_index)[i].clear();
		if (i < tetrahedron_begin_obj_index) {
			obj_BVH[i].search(aabb, compare_index, i == obj_No, &((*obj_neighbor_index)[i]), 1, 0, (*cloth)[i].triangle_AABB.size());
		}
		else {
			obj_BVH[i].search(aabb, compare_index, i == obj_No, &((*obj_neighbor_index)[i]), 1, 0, (*tetrahedron)[i - tetrahedron_begin_obj_index].triangle_AABB.size());
		}
	}
	for (int i = 0; i < collider->size(); ++i) {
		(*collider_neighbor_index)[i].clear();
		collider_BVH[i].search(aabb, compare_index, false, &((*collider_neighbor_index)[i]), 1, 0, (*collider)[i].triangle_AABB.size());
	}
}

////FIND_PATCH_PAIRS
//void Collision::findAllPatchPairs(int thread_No)
//{
//	unsigned int* thread_begin;
//	unsigned int end;
//	unsigned int obj_No;
//
//	for (int i = 0; i < total_obj_with_collider; ++i) {
//		thread_begin = mesh_patch.patch_index_start_per_thread[i].data();
//		end = thread_begin[thread_No + 1];
//		for (int j = thread_begin[thread_No]; j < end; ++j) {
//			searchPatch(mesh_patch.patch_AABB[i][j].data(), j, i, mesh_patch.patch_is_intersect[i][j]);
//		}
//	}
//}




//FIND_TRIANGLE_PAIRS
void Collision::findAllTrianglePairs(int thread_No)
{
	unsigned int* thread_begin;
	unsigned int end;
	unsigned int obj_No;
	for (int i = 0; i < cloth->size(); ++i) {
		thread_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
		end = thread_begin[thread_No + 1];
		for (int j = thread_begin[thread_No]; j < end; ++j) {
			searchTriangle((*cloth)[i].triangle_AABB[j].data(), j, i, &(*cloth)[i].triangle_neighbor_obj_triangle[j],
				&(*cloth)[i].triangle_neighbor_collider_triangle[j], false);
		}
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		thread_begin = (*tetrahedron)[i].mesh_struct.face_index_begin_per_thread.data();
		end = thread_begin[thread_No + 1];
		for (int j = thread_begin[thread_No]; j < end; ++j) {
			searchTriangle((*tetrahedron)[i].triangle_AABB[j].data(), j, i, &(*tetrahedron)[i].triangle_neighbor_obj_triangle[j],
				&(*tetrahedron)[i].triangle_neighbor_collider_triangle[j], false);
		}
	}
	//else {
	//	std::vector<std::vector<int>> test_neighbor_triangle;
	//	test_neighbor_triangle.resize(cloth->size() + tetrahedron->size());
	//	for (int i = 0; i < cloth->size(); ++i) {
	//		thread_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
	//		end = thread_begin[thread_No + 1];
	//		for (int j = thread_begin[thread_No]; j < end; ++j) {
	//			spatial_hashing.searchTriangle((*cloth)[i].triangle_AABB[j].data(),i, j,(*cloth)[i].triangle_neighbor_obj_triangle[j].data(),
	//				false,thread_No);
	//			//searchTriangle((*cloth)[i].triangle_AABB[j], j, i, &test_neighbor_triangle,
	//			//	&test_neighbor_triangle, false);
	//			//testTwoVectorsAreSame((*cloth)[i].triangle_neighbor_obj_triangle[j], test_neighbor_triangle,i,j);
	//			
	//		}
	//	}
	//	for (int i = 0; i < collider->size(); ++i) {
	//		thread_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread.data();
	//		end = thread_begin[thread_No + 1];
	//		for (int j = thread_begin[thread_No]; j < end; ++j) {
	//			spatial_hashing.searchTriangle((*collider)[i].triangle_AABB[j].data(), i, j, (*collider)[i].triangle_neighbor_obj_triangle[j].data(),
	//				true, thread_No);			
	//			//searchTriangle((*collider)[i].triangle_AABB[j], j, i, &test_neighbor_triangle,
	//			//	&test_neighbor_triangle, true);
	//			//if (thread_No == 0) {
	//				//////std::cout <<"after func"<< (*collider)[i].triangle_neighbor_obj_triangle[j][0].size() << std::endl;
	//				//testTwoVectorsAreSame((*collider)[i].triangle_neighbor_obj_triangle[j], test_neighbor_triangle, i, j);
	//			//}
	//		}
	//	}
	//	for (int i = 0; i < tetrahedron->size(); ++i) {
	//		thread_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread.data();
	//		end = thread_begin[thread_No + 1];
	//		obj_No = tetrahedron_begin_obj_index + i;
	//		for (int j = thread_begin[thread_No]; j < end; ++j) {
	//			spatial_hashing.searchTriangle((*tetrahedron)[i].triangle_AABB[j].data(), obj_No, j, (*tetrahedron)[i].triangle_neighbor_obj_triangle[j].data(),
	//				true, thread_No);
	//		}
	//	}
	//}
}


void Collision::testTwoVectorsAreSame(std::vector<std::vector<int>>& vec1, std::vector<std::vector<int>>& vec2, unsigned int obj_index,
	unsigned int triangle_index)
{
	/*for (int i = 0; i < vec1.size(); ++i) {
		if (vec1[i].size() != vec2[i].size()) {
			////std::cout << "error " << obj_index << " " << triangle_index << " size: " << vec1[i].size() << " " << vec2[i].size() << std::endl;
			if (!vec1[i].empty()) {
				////std::cout << vec1[i][0] << std::endl;
				////std::cout << "is intersect " << (*collider)[obj_index].triangle_AABB[triangle_index].AABB_intersection((*cloth)[i].triangle_AABB[vec1[i][0]]);
			}

			break;
		}
		std::sort(vec1[i].begin(), vec1[i].end());
		std::sort(vec2[i].begin(), vec2[i].end());
		for (int j = 0; j < vec1[i].size(); ++j) {
			if (vec1[i][j] != vec2[i][j]) {
				//std::cout<<"error element " << obj_index << " " << triangle_index << " element: " << vec1[i][j] << " " << vec2[i][j] << std::endl;
				break;
			}
		}
	}*/
}

//FIND_PRIMITIVE_AROUND
void Collision::findPrimitivesAround(int thread_No)
{
	findObjTriangleAroundVertex(thread_No);
	//if (use_BVH) {
	//	findColliderTriangleAroundVertex(thread_No);
	//}
	//else {
	//	findVertexAroundColliderTriangle(thread_No);
	//}
	findEdgeAroundEdge(thread_No);
}

void Collision::findVertexAroundColliderTriangle(int thread_No)
{
	//unsigned int* thread_begin;
	//std::vector<std::vector<std::vector<unsigned int>>>* triangle_neighbor_triangle;
	//std::vector<std::vector<std::vector<int>>>* triangle_neighbor_vertex; //triangle index near vertex
	//int* rep_vertex_num; // record of vertex's representative triangle index

	//std::vector<int>* triangle_neighbor_this_obj_vertex; //vector to record the triangle index near vertex which triangle and vertex in same cloth
	//std::vector<unsigned int>* triangle_neighbor_this_obj_triangle;//vector to record the triangle index near triangle in same cloth

	//std::array<double,6>* triangle_aabb;
	//std::array<double, 6>* compared_vertex_aabb;
	//int cloth_triangle_index;

	//std::array<int,3>* faces;
	//unsigned int obj_No;
	//unsigned int start, end;
	//for (int i = 0; i < collider->size(); ++i) {
	//	thread_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread.data();
	//	triangle_aabb = (*collider)[i].triangle_AABB.data();
	//	triangle_neighbor_vertex = &(*collider)[i].triangle_neighbor_obj_vertex;
	//	triangle_neighbor_triangle = &(*collider)[i].triangle_neighbor_obj_triangle;
	//	start = thread_begin[thread_No];
	//	end = thread_begin[thread_No+1];
	//	for (int k = 0; k < total_obj_num; ++k) {
	//		if (k < tetrahedron_begin_obj_index) {
	//			compared_vertex_aabb = (*cloth)[k].vertex_AABB.data();
	//			rep_vertex_num = (*cloth)[k].representative_vertex_num.data();
	//			faces = (*cloth)[k].mesh_struct.surface_triangle_index_in_order.data();
	//		}
	//		else {
	//			obj_No = k - tetrahedron_begin_obj_index;
	//			compared_vertex_aabb = (*tetrahedron)[obj_No].vertex_AABB.data();
	//			rep_vertex_num = (*tetrahedron)[obj_No].representative_vertex_num.data();
	//			faces = (*tetrahedron)[obj_No].mesh_struct.surface_triangle_index_in_order.data();
	//		}
	//		
	//		for (int j = start; j < end; ++j) {
	//			triangle_neighbor_this_obj_vertex = &(*triangle_neighbor_vertex)[j][k];
	//			triangle_neighbor_this_obj_triangle = &(*triangle_neighbor_triangle)[j][k];
	//			triangle_neighbor_this_obj_vertex->clear();				
	//			triangle_neighbor_this_obj_vertex->reserve(triangle_neighbor_this_obj_triangle->size());
	//			for (int m = 0; m < triangle_neighbor_this_obj_triangle->size(); ++m) {
	//				cloth_triangle_index = (*triangle_neighbor_this_obj_triangle)[m];
	//				for (int n = 0; n < rep_vertex_num[cloth_triangle_index]; ++n) {
	//					if (AABB::AABB_intersection(triangle_aabb[j].data(),
	//						compared_vertex_aabb[faces[cloth_triangle_index][n]].data())) {
	//						triangle_neighbor_this_obj_vertex->push_back(faces[cloth_triangle_index][n]);
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
}

void Collision::findObjTriangleAroundVertex(int thread_No)
{
	unsigned int* thread_begin;
	std::vector<std::vector<std::vector<unsigned int>>>* triangle_neighbor_triangle;
	std::vector<std::vector<std::vector<int>>>* vertex_neighbor_triangle; //triangle index near vertex
	std::vector<int>* vertex_from_rep_triangle_index; // record of vertex's representative triangle index
	std::vector<std::array<int, 3>>* face_indices;//the record of every triangle's index
	std::vector<int>* vertex_neighbor_this_obj_triangle; //vector to record the triangle index near vertex which triangle and vertex in same cloth
	std::vector<unsigned int>* triangle_neighbor_this_obj_triangle;//vector to record the triangle index near triangle in same cloth
	std::array<double, 6>* vertex_aabb;
	std::array<double, 6>* compared_triangle_aabb;
	unsigned int obj_No;
	unsigned int start, end;

	unsigned int vertex_global_index;

	for (int i = 0; i < total_obj_num; ++i) {
		if (i < tetrahedron_begin_obj_index) {
			thread_begin = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread.data();
			triangle_neighbor_triangle = &(*cloth)[i].triangle_neighbor_obj_triangle;
			vertex_neighbor_triangle = &(*cloth)[i].vertex_neighbor_obj_triangle;
			vertex_from_rep_triangle_index = &(*cloth)[i].vertex_from_rep_triangle_index;
			vertex_aabb = (*cloth)[i].vertex_AABB.data();
		}
		else {
			obj_No = i - tetrahedron_begin_obj_index;
			thread_begin = (*tetrahedron)[obj_No].mesh_struct.vertex_index_on_surface_begin_per_thread.data();
			triangle_neighbor_triangle = &(*tetrahedron)[obj_No].triangle_neighbor_obj_triangle;
			vertex_neighbor_triangle = &(*tetrahedron)[obj_No].surface_vertex_neighbor_obj_triangle;
			vertex_from_rep_triangle_index = &(*tetrahedron)[obj_No].surface_vertex_from_rep_triangle_index;
			vertex_aabb = (*tetrahedron)[obj_No].vertex_AABB.data();
		}
		start = thread_begin[thread_No];
		end = thread_begin[thread_No+1];
		for (int k = 0; k < total_obj_num; ++k) {
			if (k < tetrahedron_begin_obj_index) {
				compared_triangle_aabb = (*cloth)[k].triangle_AABB.data();
			}
			else {
				compared_triangle_aabb = (*tetrahedron)[k-tetrahedron_begin_obj_index].triangle_AABB.data();
			}
			if (i != k) {
				for (int j = start; j < end; ++j) {
					vertex_neighbor_this_obj_triangle = &(*vertex_neighbor_triangle)[j][k];
					triangle_neighbor_this_obj_triangle = &(*triangle_neighbor_triangle)[(*vertex_from_rep_triangle_index)[j]][k];//vertex's represent triangle index = (*vertex_from_rep_triangle_index)[j];
					vertex_neighbor_this_obj_triangle->clear();
					vertex_neighbor_this_obj_triangle->reserve(triangle_neighbor_this_obj_triangle->size());
					for (int m = 0; m < triangle_neighbor_this_obj_triangle->size(); ++m) {						
						if (AABB::AABB_intersection(vertex_aabb[j].data(), 
							compared_triangle_aabb[(*triangle_neighbor_this_obj_triangle)[m]].data())) {
							vertex_neighbor_this_obj_triangle->push_back((*triangle_neighbor_this_obj_triangle)[m]);
						}
					}
				}
			}
			else {
				if (k < tetrahedron_begin_obj_index) {
					face_indices = &(*cloth)[k].mesh_struct.triangle_indices;
				}
				else {
					face_indices = &(*tetrahedron)[k - tetrahedron_begin_obj_index].mesh_struct.surface_triangle_index_in_order;
				}
				for (int j = start; j < end; ++j) {
					vertex_neighbor_this_obj_triangle = &(*vertex_neighbor_triangle)[j][k];
					triangle_neighbor_this_obj_triangle = &(*triangle_neighbor_triangle)[(*vertex_from_rep_triangle_index)[j]][k];//vertex's represent triangle index = (*vertex_from_rep_triangle_index)[j];
					vertex_neighbor_this_obj_triangle->clear();
					vertex_neighbor_this_obj_triangle->reserve(triangle_neighbor_this_obj_triangle->size());
					
					vertex_global_index = vertex_index_on_surface[k][j];

					for (int m = 0; m < triangle_neighbor_this_obj_triangle->size(); ++m) {
						if (!vertexInTriangle((*face_indices)[(*triangle_neighbor_this_obj_triangle)[m]].data(), vertex_global_index)) {
							if (AABB::AABB_intersection(vertex_aabb[j].data(),
								compared_triangle_aabb[(*triangle_neighbor_this_obj_triangle)[m]].data())) {
								vertex_neighbor_this_obj_triangle->push_back((*triangle_neighbor_this_obj_triangle)[m]);
							}
						}
					}
				}
			}
		}
	}
}


void Collision::findColliderTriangleAroundVertex(int thread_No)
{
	unsigned int* thread_begin;
	std::vector<std::vector<std::vector<unsigned int>>>* triangle_neighbor_triangle;
	std::vector<std::vector<std::vector<int>>>* vertex_neighbor_triangle; //triangle index near vertex
	std::vector<int>* vertex_from_rep_triangle_index; // record of vertex's representative triangle index
	std::vector<int>* vertex_neighbor_this_obj_triangle; //vector to record the triangle index near vertex which triangle and vertex in same cloth
	std::vector<unsigned int>* triangle_neighbor_this_obj_triangle;//vector to record the triangle index near triangle in same cloth
	std::array<double, 6>* vertex_aabb;
	std::array<double, 6>* compared_triangle_aabb;
	unsigned int obj_No;
	unsigned int start, end;
	for (int i = 0; i < total_obj_num; ++i) {
		if (i < tetrahedron_begin_obj_index) {
			thread_begin = (*cloth)[i].mesh_struct.vertex_index_begin_per_thread.data();
			vertex_from_rep_triangle_index = &(*cloth)[i].vertex_from_rep_triangle_index;
			vertex_aabb = (*cloth)[i].vertex_AABB.data();
			vertex_neighbor_triangle = &(*cloth)[i].vertex_neighbor_collider_triangle;
			triangle_neighbor_triangle = &(*cloth)[i].triangle_neighbor_collider_triangle;
		}
		else {
			obj_No = i - tetrahedron_begin_obj_index;
			thread_begin = (*tetrahedron)[obj_No].mesh_struct.vertex_index_on_surface_begin_per_thread.data();
			vertex_from_rep_triangle_index = &(*tetrahedron)[obj_No].surface_vertex_from_rep_triangle_index;
			vertex_aabb = (*tetrahedron)[obj_No].vertex_AABB.data();
			vertex_neighbor_triangle = &(*tetrahedron)[obj_No].surface_vertex_neighbor_collider_triangle;
			triangle_neighbor_triangle = &(*tetrahedron)[obj_No].triangle_neighbor_collider_triangle;
		}		
		start = thread_begin[thread_No];
		end = thread_begin[thread_No + 1];
		for (int k = 0; k < collider->size(); ++k) {
			compared_triangle_aabb = (*collider)[k].triangle_AABB.data();
			for (int j = start; j < end; ++j) {
				vertex_neighbor_this_obj_triangle = &(*vertex_neighbor_triangle)[j][k];
				triangle_neighbor_this_obj_triangle = &(*triangle_neighbor_triangle)[(*vertex_from_rep_triangle_index)[j]][k];//vertex's represent triangle index = (*vertex_from_rep_triangle_index)[j];
				vertex_neighbor_this_obj_triangle->clear();
				vertex_neighbor_this_obj_triangle->reserve(triangle_neighbor_this_obj_triangle->size());
				for (int m = 0; m < triangle_neighbor_this_obj_triangle->size(); ++m) {
					if (AABB::AABB_intersection(vertex_aabb[j].data(),
						compared_triangle_aabb[(*triangle_neighbor_this_obj_triangle)[m]].data())) {
						(*vertex_neighbor_this_obj_triangle).push_back((*triangle_neighbor_this_obj_triangle)[m]);
					}
				}
			}
		}
	}
}


void Collision::findEdgeAroundEdge(int thread_No)
{
	std::vector<std::array<int, 3>>* face_indices;//the record of every triangle's index
	unsigned int* thread_begin;
	std::vector<std::vector<std::vector<unsigned int>>>* triangle_neighbor_triangle;
	std::vector<std::vector<std::vector<int>>>* edge_neighbor_edge; //edge index near edge
	std::vector<int>* edge_from_rep_triangle_index; // record of edge's representative triangle index
	std::vector<int>* edge_neighbor_one_obj_edge; //vector to record the edge index near edge which triangle and vertex in same cloth
	std::vector<unsigned int>* triangle_neighbor_one_obj_triangle;//vector to record the triangle index near triangle in same cloth
	std::vector<unsigned int>* representative_edge_num;
	std::vector<MeshStruct::Face>* face;
	int face_index;
	std::vector<MeshStruct::Edge>* edge;
	std::array<double, 6>* edge_aabb;
	std::array<double, 6>* compared_edge_aabb;
	int compare_edge_index;
	unsigned int obj_No, obj_No2;

	unsigned int* face_edges;
	unsigned int* edge_vertex_index;

	for (int i = 0; i < total_obj_num; ++i) {
		if (i < tetrahedron_begin_obj_index) {
			thread_begin = (*cloth)[i].mesh_struct.edge_index_begin_per_thread.data();
			triangle_neighbor_triangle = &(*cloth)[i].triangle_neighbor_obj_triangle;
			edge_neighbor_edge = &(*cloth)[i].edge_neighbor_obj_edge;
			edge_from_rep_triangle_index = &(*cloth)[i].edge_from_rep_triangle_index;
			edge = &(*cloth)[i].mesh_struct.edges;
			edge_aabb = (*cloth)[i].edge_AABB.data();

			edge_vertex_index = (*cloth)[i].mesh_struct.edge_vertices.data();

		}
		else {
			obj_No = i - tetrahedron_begin_obj_index;
			thread_begin = (*tetrahedron)[obj_No].mesh_struct.edge_index_begin_per_thread.data();
			triangle_neighbor_triangle = &(*tetrahedron)[obj_No].triangle_neighbor_obj_triangle;
			edge_neighbor_edge = &(*tetrahedron)[obj_No].edge_neighbor_obj_edge;
			edge_from_rep_triangle_index = &(*tetrahedron)[obj_No].edge_from_rep_triangle_index;
			edge = &(*tetrahedron)[obj_No].mesh_struct.edges;
			edge_aabb = (*tetrahedron)[obj_No].edge_AABB.data();
			edge_vertex_index = (*tetrahedron)[obj_No].mesh_struct.edge_vertices.data();
		}
		for (int k = i; k < total_obj_num; ++k) {
			if (k < tetrahedron_begin_obj_index) {
				compared_edge_aabb = (*cloth)[k].edge_AABB.data();
				representative_edge_num = &(*cloth)[k].representative_edge_num;
				face = &(*cloth)[k].mesh_struct.faces;
				face_edges = (*cloth)[k].mesh_struct.face_edges.data();
			}
			else {
				obj_No2 = k - tetrahedron_begin_obj_index;
				compared_edge_aabb = (*tetrahedron)[obj_No2].edge_AABB.data();
				representative_edge_num = &(*tetrahedron)[obj_No2].representative_edge_num;
				face = &(*tetrahedron)[obj_No2].mesh_struct.faces;
				face_edges = (*tetrahedron)[obj_No2].mesh_struct.face_edges.data();
			}
			if (i != k) {
				for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
					edge_neighbor_one_obj_edge = &(*edge_neighbor_edge)[j][k];
					triangle_neighbor_one_obj_triangle = &(*triangle_neighbor_triangle)[(*edge_from_rep_triangle_index)[j]][k];
					edge_neighbor_one_obj_edge->clear();
					edge_neighbor_one_obj_edge->reserve(triangle_neighbor_one_obj_triangle->size());					
					for (int m = 0; m < triangle_neighbor_one_obj_triangle->size(); ++m) {
						face_index = (*triangle_neighbor_one_obj_triangle)[m];
						for (int n = 0; n < (*representative_edge_num)[face_index]; ++n) {
							if (AABB::AABB_intersection(edge_aabb[j].data(), 								
								compared_edge_aabb[face_edges[3*face_index+n]].data())) {
								(*edge_neighbor_one_obj_edge).push_back(face_edges[3 * face_index + n]);
							}
						}						
					}
				}
			}
			else {
				for (int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
					edge_neighbor_one_obj_edge = &(*edge_neighbor_edge)[j][k];
					triangle_neighbor_one_obj_triangle = &(*triangle_neighbor_triangle)[(*edge_from_rep_triangle_index)[j]][k];
					edge_neighbor_one_obj_edge->clear();
					edge_neighbor_one_obj_edge->reserve(triangle_neighbor_one_obj_triangle->size());
					for (int m = 0; m < triangle_neighbor_one_obj_triangle->size(); ++m) {
						face_index = (*triangle_neighbor_one_obj_triangle)[m];
						for (int n = 0; n < (*representative_edge_num)[face_index]; ++n) {
							compare_edge_index = face_edges[3 * face_index + n];// (*face)[face_index].edge[n];
							if (j < compare_edge_index) {
								if (!edgeEdgeconnected(edge_vertex_index+2*j, edge_vertex_index + 2*compare_edge_index)) {
									if (AABB::AABB_intersection(edge_aabb[j].data(), 
										compared_edge_aabb[compare_edge_index].data())) {
										(*edge_neighbor_one_obj_edge).push_back(compare_edge_index);
									}
								}
							}							
						}
					}
				}
			}
		}
	}
}


void Collision::initialVolume()
{
	VT_volume.resize(total_obj_num);
	TV_volume.resize(total_obj_num);
	EE_volume.resize(total_obj_num);
	VT_collider_volume.resize(total_obj_num);
	VT_start_index = new unsigned int* [total_obj_num];
	TV_start_index = new unsigned int* [total_obj_num];
	EE_start_index = new unsigned int* [total_obj_num];
	VT_collider_start_index = new unsigned int* [total_obj_num];
	for (int i = 0; i < total_obj_num; ++i) {
		VT_volume[i].reserve(mesh_struct[i]->vertex_position.size());
		TV_volume[i].reserve(mesh_struct[i]->vertex_position.size());
		EE_volume[i].reserve(mesh_struct[i]->vertex_position.size());
		VT_start_index[i] = new unsigned int[vertex_index_start_per_thread[i][thread_num] + 1];
		TV_start_index[i] = new unsigned int[mesh_struct[i]->triangle_indices.size()+1];
		EE_start_index[i] = new unsigned int[mesh_struct[i]->edge_vertices.size()/2+1];
		if (!collider->empty()) {
			VT_collider_volume[i].reserve(mesh_struct[i]->vertex_position.size());
			VT_collider_start_index[i] = new unsigned int[vertex_index_start_per_thread[i][thread_num] + 1];
		}
	}
}



void Collision::initialHessianRecord()
{
	int total_vertex_num=0;

	vertex_num_on_surface_prefix_sum.resize(total_obj_num+1);
	vertex_num_on_surface_prefix_sum[0] = 0;
	for (int i = 0; i < total_obj_num; ++i) {
		total_vertex_num += vertex_index_start_per_thread[i][thread_num];
	}
	for (int i = 1; i <= total_obj_num; ++i) {
		vertex_num_on_surface_prefix_sum[i] = vertex_num_on_surface_prefix_sum[i - 1] + vertex_index_start_per_thread[i - 1][thread_num];
	}

	vt_hessian_record_index.reserve(total_vertex_num * 50);
	vt_hessian_record.reserve(total_vertex_num * 1440);
	vt_grad_record.reserve(total_vertex_num * 120);

	ee_hessian_record_index.reserve(total_vertex_num * 50);
	ee_hessian_record.reserve(total_vertex_num * 1440);
	ee_grad_record.reserve(total_vertex_num * 120);

	if (!collider->empty()) {
		vt_colldier_hessian_record.resize(total_vertex_num*9);
		vt_colldier_grad_record.resize(total_vertex_num*3);
		//vt_colldier_hessian_record_is_not_empty=new bool*[total_obj_num];

		ee_collider_hessian_record_index.reserve(total_vertex_num*30);
		ee_collider_hessian_record.reserve(total_vertex_num*360);
		ee_collider_grad_record.reserve(total_vertex_num*60);

		tv_colldier_hessian_record_index.reserve(total_vertex_num*40);
		tv_colldier_hessian_record.reserve(total_vertex_num*810);
		tv_colldier_grad_record.reserve(total_vertex_num*90);
	}




	vertex_belong_to_color_group = new bool* [total_obj_num];
	for (int i = 0; i < total_obj_num; ++i) {
		vertex_belong_to_color_group[i] = new bool[mesh_struct[i]->vertex_position.size()];
	}
}




void Collision::initialPairRecord()
{
	vertex_triangle_pair_by_vertex = new unsigned int* [total_obj_num];
	triangle_vertex_pair_by_triangle = new unsigned int* [total_obj_num];
	edge_edge_pair_by_edge = new unsigned int* [total_obj_num];
	vertex_obj_triangle_collider_pair_by_vertex = new unsigned int* [total_obj_num];

	vertex_triangle_pair_num_record = new unsigned int* [total_obj_num];
	triangle_vertex_pair_num_record = new unsigned int* [total_obj_num];
	edge_edge_pair_number_record = new unsigned int* [total_obj_num];
	vertex_obj_triangle_collider_num_record = new unsigned int* [total_obj_num];


	triangle_vertex_collider_pair_by_triangle= new unsigned int* [total_obj_num]; 
	triangle_vertex_collider_pair_num_record = new unsigned int* [total_obj_num];
	edge_edge_collider_pair_by_edge = new unsigned int* [total_obj_num];
	edge_edge_collider_pair_num_record = new unsigned int* [total_obj_num];


	vertex_triangle_pair_num_record_prefix_sum = new unsigned int* [total_obj_num];
	edge_edge_pair_num_record_prefix_sum = new unsigned int* [total_obj_num];
	triangle_vertex_collider_num_record_prefix_sum = new unsigned int* [total_obj_num];
	edge_edge_collider_num_record_prefix_sum = new unsigned int* [total_obj_num];


	triangle_index_collide_with_collider.resize(total_obj_num);
	edge_index_collide_with_collider.resize(total_obj_num);
	vertex_index_collide_with_collider.resize(total_obj_num);


	//triangle_index_collide_with_collider_prefix_sum.resize(total_obj_num + 1,0);
	//edge_index_collide_with_collider_prefix_sum.resize(total_obj_num + 1,0);
	//vertex_index_collide_with_collider_prefix_sum.resize(total_obj_num + 1,0);

	unsigned int vertex_num, triangle_num, edge_num;

	for (int i = 0; i <total_obj_num; ++i) {
		edge_edge_pair_by_edge[i] = new unsigned int[close_ee_pair_num * edge_index_start_per_thread[i][thread_num]];// 
		memset(edge_edge_pair_by_edge[i], 0, 4 * close_ee_pair_num * edge_index_start_per_thread[i][thread_num]);

		vertex_triangle_pair_by_vertex[i] = new unsigned int[close_vt_pair_num * vertex_index_start_per_thread[i][thread_num]];// 
		memset(vertex_triangle_pair_by_vertex[i], 0, 4 * (close_vt_pair_num * vertex_index_start_per_thread[i][thread_num]));
		vertex_triangle_pair_num_record[i] = new unsigned int[vertex_index_start_per_thread[i][thread_num]];// 
		memset(vertex_triangle_pair_num_record[i], 0, 4 * vertex_index_start_per_thread[i][thread_num]);
	
		triangle_vertex_pair_by_triangle[i] = new unsigned int[close_tv_pair_num * mesh_struct[i]->triangle_indices.size()];// 
		memset(triangle_vertex_pair_by_triangle[i], 0, 4 * (close_tv_pair_num * mesh_struct[i]->triangle_indices.size()));

		edge_edge_pair_number_record[i] = new unsigned int[mesh_struct[i]->edge_length.size()];// 
		memset(edge_edge_pair_number_record[i], 0, 4 * mesh_struct[i]->edge_length.size());
		triangle_vertex_pair_num_record[i] = new unsigned int[mesh_struct[i]->triangle_indices.size()];// 
		memset(triangle_vertex_pair_num_record[i], 0, 4 * mesh_struct[i]->triangle_indices.size());

	

		edge_edge_pair_num_record_prefix_sum[i] = new unsigned int[mesh_struct[i]->edge_length.size()+1];// 
		vertex_triangle_pair_num_record_prefix_sum[i] = new unsigned int[vertex_index_start_per_thread[i][thread_num]+1];// 

		if (has_collider) {
			triangle_index_collide_with_collider[i].reserve(mesh_struct[i]->triangle_indices.size() / 4);
			edge_index_collide_with_collider[i].reserve(mesh_struct[i]->edge_vertices.size() / 8);
			vertex_index_collide_with_collider[i].reserve(mesh_struct[i]->vertex_position.size() / 4);


			vertex_obj_triangle_collider_pair_by_vertex[i] = new unsigned int[close_vt_collider_pair_num * vertex_index_start_per_thread[i][thread_num]];
			memset(vertex_obj_triangle_collider_pair_by_vertex[i], 0, 4 * (close_vt_collider_pair_num * vertex_index_start_per_thread[i][thread_num]));

			vertex_obj_triangle_collider_num_record[i] = new unsigned int[vertex_index_start_per_thread[i][thread_num]];// 
			memset(vertex_obj_triangle_collider_num_record[i], 0, 4 * vertex_index_start_per_thread[i][thread_num]);


			triangle_vertex_collider_pair_by_triangle[i] = new unsigned int[close_tv_collider_pair_num * triangle_index_start_per_thread[i][thread_num]];
			memset(triangle_vertex_collider_pair_by_triangle[i], 0, 4 * (close_tv_collider_pair_num * triangle_index_start_per_thread[i][thread_num]));

			triangle_vertex_collider_pair_num_record[i] = new unsigned int[triangle_index_start_per_thread[i][thread_num]];// 
			memset(triangle_vertex_collider_pair_num_record[i], 0, 4 * triangle_index_start_per_thread[i][thread_num]);

			edge_edge_collider_pair_by_edge[i] = new unsigned int[close_ee_collider_pair_num * edge_index_start_per_thread[i][thread_num]];
			memset(edge_edge_collider_pair_by_edge[i], 0, 4 * (close_ee_collider_pair_num * edge_index_start_per_thread[i][thread_num]));

			edge_edge_collider_pair_num_record[i] = new unsigned int[edge_index_start_per_thread[i][thread_num]];// 
			memset(edge_edge_collider_pair_num_record[i], 0, 4 * edge_index_start_per_thread[i][thread_num]);


		
			triangle_vertex_collider_num_record_prefix_sum[i] = new unsigned int[mesh_struct[i]->triangle_indices.size()+1];

			edge_edge_collider_num_record_prefix_sum[i] = new unsigned int[edge_index_start_per_thread[i][thread_num]+1];


		}
		else {
			//edge_edge_pair_collider[i] = new unsigned int[1];
			vertex_obj_triangle_collider_pair_by_vertex[i] = new unsigned int[1];
			vertex_obj_triangle_collider_num_record[i] = new unsigned int[1];
			//vertex_collider_triangle_obj_pair[i] = new unsigned int[1];

			triangle_vertex_collider_pair_by_triangle[i] = new unsigned int[1];
			triangle_vertex_collider_pair_num_record[i] = new unsigned int[1];
			edge_edge_collider_pair_by_edge[i] = new unsigned int[1]; 
			edge_edge_collider_pair_num_record[i] = new unsigned int[1];

			triangle_vertex_collider_num_record_prefix_sum[i]= new unsigned int[1];
			edge_edge_collider_num_record_prefix_sum[i]= new unsigned int[1];

		}
	}
	initialVolume();
}


void Collision::initialPair()
{
	unsigned int total_triangle_num = 0;
	for (int i = 0; i < cloth->size(); ++i) {
		total_triangle_num += cloth->data()[i].mesh_struct.triangle_indices.size();
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		total_triangle_num += tetrahedron->data()[i].mesh_struct.triangle_indices.size();
	}

	for (int i = 0; i < collider->size(); ++i) {
		total_triangle_num += collider->data()[i].mesh_struct.triangle_indices.size();
	}
	
	if (CCD_compare) {
		vertex_edge_pair = new unsigned int* [thread_num];
		vertex_obj_edge_collider_pair = new unsigned int* [thread_num];
		vertex_collider_edge_obj_pair = new unsigned int* [thread_num];

		vertex_vertex_pair = new unsigned int* [thread_num];
		vertex_vertex_pair_collider = new unsigned int* [thread_num];

		for (int i = 0; i < thread_num; ++i) {
			vertex_edge_pair[i] = new unsigned int[estimate_coeff_for_vt_pair_num * total_triangle_num];// 
			memset(vertex_edge_pair[i], 0, 4 * (estimate_coeff_for_vt_pair_num * total_triangle_num));
			vertex_vertex_pair[i] = new unsigned int[estimate_coeff_for_vt_pair_num * total_triangle_num / 4];// 
			memset(vertex_vertex_pair[i], 0, 4 * (estimate_coeff_for_vt_pair_num * total_triangle_num / 4));

			if (has_collider) {
				vertex_obj_edge_collider_pair[i] = new unsigned int[estimate_coeff_for_vt_pair_num * total_triangle_num];
				memset(vertex_obj_edge_collider_pair[i], 0, 4 * (estimate_coeff_for_vt_pair_num * total_triangle_num));

				vertex_collider_edge_obj_pair[i] = new unsigned int[estimate_coeff_for_vt_pair_num * total_triangle_num];
				memset(vertex_collider_edge_obj_pair[i], 0, 4 * (estimate_coeff_for_vt_pair_num * total_triangle_num));

				vertex_vertex_pair_collider[i] = new unsigned int[estimate_coeff_for_vt_pair_num * total_triangle_num / 4];
				memset(vertex_vertex_pair_collider[i], 0, 4 * (estimate_coeff_for_vt_pair_num * total_triangle_num / 4));

				
			}
			else {
				vertex_obj_edge_collider_pair[i] = new unsigned int[1];
				vertex_collider_edge_obj_pair[i] = new unsigned int[1];
				vertex_vertex_pair_collider[i] = new unsigned int[1];
			}
		}
	}
	unsigned int pair_num = 0;
	for (int i = 0; i < thread_num; ++i) {
		pair_num += estimate_coeff_for_vt_pair_num * total_triangle_num / 2;
		if (has_collider) {
			pair_num += estimate_coeff_for_vt_pair_num * total_triangle_num / 4;
		}
	}
	//target_position_and_stiffness.resize(thread_num);
//target_position_index.resize(thread_num);

	//should not appear when using PT

	point_triangle_target_pos_index.resize(thread_num);
	//point_collider_triangle_target_pos_index.resize(thread_num);
	point_triangle_collider_target_pos_index.resize(thread_num);
	edge_edge_target_pos_index.resize(thread_num);
	//edge_edge_collider_target_pos_index.resize(thread_num);

	point_triangle_target_pos_record.resize(thread_num);
	//point_collider_triangle_target_pos_index.resize(thread_num);
	point_triangle_collider_target_pos_record.resize(thread_num);
	edge_edge_target_pos_record.resize(thread_num);


	for (int i = 0; i < thread_num; ++i) {
		//target_position_and_stiffness[i].resize(pair_num * 4);
		//target_position_index[i].resize(pair_num * 4 + 1);
		point_triangle_target_pos_index[i].resize(pair_num * 4 + 1);				
		edge_edge_target_pos_index[i].resize(pair_num * 4 + 1);

		point_triangle_target_pos_record[i].resize(pair_num * 13);
		edge_edge_target_pos_record[i].resize(pair_num * 13);

		if (has_collider) {
			//point_collider_triangle_target_pos_index[i].resize(pair_num * 4 + 1);
			//edge_edge_collider_target_pos_index[i].resize(pair_num * 4 + 1);
			point_triangle_collider_target_pos_index[i].resize(pair_num * 4 + 1);
			point_triangle_collider_target_pos_record[i].resize(pair_num * 4);
		}

	}
	
}


void Collision::testPointPair()
{
	test_triangle_index.clear();
	 chosen_show_vertex = 350;
	unsigned int pair_num;
	for (unsigned int i = 0; i < thread_num; ++i) {
		pair_num = spatial_hashing.vertex_obj_triangle_collider_pair[i][0];
		for (unsigned int j = 0; j < pair_num; j += 4) {
			if (*(spatial_hashing.vertex_obj_triangle_collider_pair[i] + 1 + j) == chosen_show_vertex) {
				test_triangle_index.push_back(*(spatial_hashing.vertex_obj_triangle_collider_pair[i] + 3 + j));
			}
		}		
	}
}






//PREFIX_SUM_ALL_PAIRS
void Collision::prefixSumAllPair(int thread_No)
{
	unsigned int vt_start = 0, ee_start = 0, tv_collider_start = 0, ee_collider_start = 0;
	if (thread_num == 1) {		
		for (unsigned int i = 0; i < total_obj_num; ++i) {
			prefixSumRecordPairNum(vertex_triangle_pair_num_record[i], vertex_triangle_pair_num_record_prefix_sum[i], vertex_index_start_per_thread[i][thread_num],vt_start,2);
			prefixSumRecordPairNum(edge_edge_pair_number_record[i], edge_edge_pair_num_record_prefix_sum[i], edge_index_start_per_thread[i][thread_num],ee_start,3);
			if (has_collider) {
				prefixSumRecordPairNum(triangle_vertex_collider_pair_num_record[i], triangle_vertex_collider_num_record_prefix_sum[i], triangle_index_start_per_thread[i][thread_num], tv_collider_start,2);
				prefixSumRecordPairNum(edge_edge_collider_pair_num_record[i], edge_edge_collider_num_record_prefix_sum[i], edge_index_start_per_thread[i][thread_num], ee_collider_start,2);
			}
		}		
	}
	else if (thread_num > 3) {
		switch (thread_No)
		{
		case 0:
			for (unsigned int i = 0; i < total_obj_num; ++i) {
				prefixSumRecordPairNum(vertex_triangle_pair_num_record[i], vertex_triangle_pair_num_record_prefix_sum[i], vertex_index_start_per_thread[i][thread_num], vt_start,2);
			}
			break;
		case 1:
			for (unsigned int i = 0; i < total_obj_num; ++i) {
				prefixSumRecordPairNum(edge_edge_pair_number_record[i], edge_edge_pair_num_record_prefix_sum[i], edge_index_start_per_thread[i][thread_num], ee_start,3);
			}
			break;
		case 2:
			if (has_collider) {
				for (unsigned int i = 0; i < total_obj_num; ++i) {
					prefixSumRecordPairNum(triangle_vertex_collider_pair_num_record[i], triangle_vertex_collider_num_record_prefix_sum[i], triangle_index_start_per_thread[i][thread_num], tv_collider_start,2);
				}
			}
			break;
		case 3:
			if (has_collider) {
				for (unsigned int i = 0; i < total_obj_num; ++i) {
					prefixSumRecordPairNum(edge_edge_collider_pair_num_record[i], edge_edge_collider_num_record_prefix_sum[i], edge_index_start_per_thread[i][thread_num], ee_collider_start,2);
				}
			}
			break;
		}
	}
	else {
		switch (thread_No)
		{
		case 0:
			for (unsigned int i = 0; i < total_obj_num; ++i) {
				prefixSumRecordPairNum(vertex_triangle_pair_num_record[i], vertex_triangle_pair_num_record_prefix_sum[i], vertex_index_start_per_thread[i][thread_num], vt_start,2);
				if (has_collider) {
					prefixSumRecordPairNum(edge_edge_collider_pair_num_record[i], edge_edge_collider_num_record_prefix_sum[i], edge_index_start_per_thread[i][thread_num], ee_collider_start,2);
				}
			}
			break;
		case 1:
			for (unsigned int i = 0; i < total_obj_num; ++i) {
				prefixSumRecordPairNum(edge_edge_pair_number_record[i], edge_edge_pair_num_record_prefix_sum[i], edge_index_start_per_thread[i][thread_num], ee_start,3);
				if (has_collider) {
					prefixSumRecordPairNum(triangle_vertex_collider_pair_num_record[i], triangle_vertex_collider_num_record_prefix_sum[i], triangle_index_start_per_thread[i][thread_num], tv_collider_start,2);
				}
			}
			break;
		}
	}
}


void Collision::prefixSumRecordPairNum(unsigned int* num_record, unsigned int* prefix_sum, int num, unsigned int& start_index, int move_size)
{
	prefix_sum[0] = start_index;
	for (int i = 1; i <= num; ++i) {
		prefix_sum[i] = prefix_sum[i - 1] + (num_record[i - 1]/move_size);
	}
	start_index = prefix_sum[num];
}



void Collision::testNearestPoint()
{
	double point[3] = {0.0,1.0,0.0 };
	double triangle_point[3][3] = { {1.0,1.0,0.0},{1.0,1.0,1.0},{0.0,1.0,1.0} };
	double normal[3] = { 0.0,1.0,0.0 };
	double bary[3];
	CCD::internal::pointTriangleNearestPoint(point, triangle_point[0], triangle_point[1], triangle_point[2], normal, bary);
	////std::cout <<"barycentric "<< bary[0] << " " << bary[1] << " " << bary[2] << std::endl;
}




void Collision::setSelfCollisionPairPerThread()
{
	//collision_pair_around_pair.resize(prefix_sum_of_different_type_pair[5] * collision_pair_around_pair_size_per_pair);
	//collision_pair_around_pair_size.resize(prefix_sum_of_different_type_pair[5]);
	//memset(collision_pair_around_pair_size.data(), 0, collision_pair_around_pair_size.size() << 2);
	vt_per_thread_start_index.resize(thread_num + 1,0);
	ee_per_thread_start_index.resize(thread_num + 1,0);
	//vt_collider_per_thread_start_index.resize(thread_num + 1,0);
	//tv_collider_per_thread_start_index.resize(thread_num + 1,0);
	//ee_collider_per_thread_start_index.resize(thread_num + 1,0);

	memset(vt_per_thread_start_index.data(), 0, vt_per_thread_start_index.size() << 2);
	memset(ee_per_thread_start_index.data(), 0, ee_per_thread_start_index.size() << 2);
	//memset(vt_collider_per_thread_start_index.data(), 0, vt_collider_per_thread_start_index.size() << 2);
	//memset(tv_collider_per_thread_start_index.data(), 0, tv_collider_per_thread_start_index.size() << 2);
	//memset(ee_collider_per_thread_start_index.data(), 0, ee_collider_per_thread_start_index.size() << 2);
	
	//setPairStartPerThread(vertex_triangle_pair_num_record_prefix_sum, vertex_index_start_per_thread.data(),
	//	vt_per_thread_start_index.data());

	//setPairStartPerThread(edge_edge_pair_num_record_prefix_sum, edge_index_start_per_thread.data(),
	//	ee_per_thread_start_index.data());
	arrangeIndex(thread_num, vt_pair_compressed_record.size() /5, vt_per_thread_start_index.data());
	arrangeIndex(thread_num, ee_pair_compressed_record.size() /5, ee_per_thread_start_index.data());
	for (int i = 0; i <= thread_num; ++i) {
		vt_per_thread_start_index[i] *= 5;
		ee_per_thread_start_index[i] *= 5;
	}
	//arrangeIndex(thread_num, triangle_index_collide_with_collider_prefix_sum[total_obj_num], tv_collider_per_thread_start_index.data());
	//arrangeIndex(thread_num, edge_index_collide_with_collider_prefix_sum[total_obj_num], ee_collider_per_thread_start_index.data());
	//arrangeIndex(thread_num, vertex_index_collide_with_collider_prefix_sum[total_obj_num], vt_collider_per_thread_start_index.data());



}


void Collision::addFlagToColorGroup(std::vector<unsigned int>& pair_compress_record, std::vector<unsigned int>** tet_around_element0, std::vector<unsigned int>** tet_around_element1)
{
	int color_size;
	int* address_of_tet_order_in_group;
	std::vector<unsigned int>::iterator end;
	int group_size;
	std::vector<std::vector<char>>* tet_color_groups_label_;
	std::vector<unsigned int>* tet_involved_in_collision_;

	for (auto i = pair_compress_record.begin(); i < pair_compress_record.end(); i += 3) { // here add another 2 in the loop
		if (*i >= cloth->size()) {
			end = tet_around_element0[*i][*(i + 1)].end();
			group_size = tet_color_groups[*i]->size();
			tet_color_groups_label_ = tet_color_groups_label[*i];
			tet_involved_in_collision_ = tet_involved_in_collision[*i].data();
			for (auto j = tet_around_element0[*i][*(i + 1)].begin(); j < end; ++j) {
				address_of_tet_order_in_group = tet_order_in_color_group[*i] + (group_size <<1) * (*j);	
				for (int k = 0; k < group_size; ++k) {
					if (*(address_of_tet_order_in_group + (k << 1)) < tet_color_groups[*i]->data()[k].size() - 1) {
						if (!tet_color_groups_label_[k][*(address_of_tet_order_in_group + (k << 1))][*(address_of_tet_order_in_group + (k << 1) + 1)]) {
							tet_color_groups_label_[k][*(address_of_tet_order_in_group + (k << 1))][*(address_of_tet_order_in_group + (k << 1) + 1)] = '\1';
							tet_involved_in_collision_[k].emplace_back(*j);
						}
					}				
				}
			}
		}

		i += 2;

		if (*i >= cloth->size()) {
			end = tet_around_element1[*i][*(i + 1)].end();
			group_size = tet_color_groups[*i]->size();
			tet_color_groups_label_ = tet_color_groups_label[*i];
			tet_involved_in_collision_ = tet_involved_in_collision[*i].data();
			for (auto j = tet_around_element1[*i][*(i + 1)].begin(); j < end; ++j) {
				address_of_tet_order_in_group = tet_order_in_color_group[*i] + (group_size << 1) * (*j);
				for (int k = 0; k < group_size; ++k) {
					if (*(address_of_tet_order_in_group + (k << 1))< tet_color_groups[*i]->data()[k].size()-1) {
						if (!tet_color_groups_label_[k][*(address_of_tet_order_in_group + (k << 1))][*(address_of_tet_order_in_group + (k << 1) + 1)]) {
							tet_color_groups_label_[k][*(address_of_tet_order_in_group + (k << 1))][*(address_of_tet_order_in_group + (k << 1) + 1)] = '\1';
							tet_involved_in_collision_[k].emplace_back(*j);
						}
					}					
				}
			}
		}
	}
}


void Collision::addFlagToColorGroup(std::vector<unsigned int>*element_collide_with_collider, std::vector<unsigned int>** tet_around_element)
{
	int group_size;
	int i;
	int* address_of_tet_order_in_group;
	std::vector<std::vector<char>>* tet_color_groups_label_;
	std::vector<unsigned int>* tet_around_element_;
	std::vector<unsigned int>* tet_involved_in_collision_;

	for (int tet_No = 0; tet_No < tetrahedron->size(); ++tet_No) {
		i = tet_No + cloth->size();
		group_size = tet_color_groups[i]->size();
		tet_color_groups_label_ = tet_color_groups_label[i];
		tet_around_element_ = tet_around_element[i];
		tet_involved_in_collision_ = tet_involved_in_collision[i].data();
		for (auto j = element_collide_with_collider[i].begin(); j < element_collide_with_collider[i].end(); ++j) {
			for (auto k = tet_around_element_[*j].begin(); k < tet_around_element_[*j].end(); ++k) {
				address_of_tet_order_in_group = tet_order_in_color_group[i] + (group_size << 1) * (*k);
				for (int m = 0; m < group_size; ++m) {
					if (*(address_of_tet_order_in_group + (m << 1)) < tet_color_groups[i]->data()[m].size() - 1) {
						if (!tet_color_groups_label_[m][*(address_of_tet_order_in_group + (m << 1))][*(address_of_tet_order_in_group + (m << 1) + 1)]) {
							tet_color_groups_label_[m][*(address_of_tet_order_in_group + (m << 1))][*(address_of_tet_order_in_group + (m << 1) + 1)] = '\1';
							tet_involved_in_collision_[m].emplace_back(*k);
						}
					}
				}
			}
		}
	}
}

void Collision::addFlagToColorGroup()
{
	//initial label
	int i;
	std::vector < std::vector<char>>* label;
	for (int tet_No = 0; tet_No < tetrahedron->size(); ++tet_No) {
		i = cloth->size() + tet_No;
		for (int j = 0; j < tet_color_groups[i]->size(); ++j) {
			label = &tet_color_groups_label[i][j];
			for (auto k = label->begin(); k < label->end(); ++k) {
				memset(k->data(), 0, k->size());
			}
		}
	}
	for (auto i = tet_involved_in_collision.begin(); i < tet_involved_in_collision.end(); ++i) {
		for (auto j = i->begin(); j < i->end(); ++j) {
			j->clear();
		}
	}

	//add label
	// vt pair
	addFlagToColorGroup(vt_pair_compressed_record, tet_around_vertex.data(), tet_around_triangle.data());
	//EE pair
	addFlagToColorGroup(ee_pair_compressed_record, tet_around_edge.data(), tet_around_edge.data());
	if (has_collider) {
		// vt collider
		addFlagToColorGroup(vertex_index_collide_with_collider.data(), tet_around_vertex.data());
		// tv collider
		addFlagToColorGroup(triangle_index_collide_with_collider.data(), tet_around_triangle.data());
		//ee collider
		addFlagToColorGroup(edge_index_collide_with_collider.data(), tet_around_edge.data());
	}

	for (int tet_No = 0; tet_No < tetrahedron->size(); ++tet_No) {
		i = cloth->size() + tet_No;
		for (int j = 0; j < tet_involved_in_collision_start_per_thread[i].size(); ++j) {
			arrangeIndex(thread_num, tet_involved_in_collision[i][j].size(), tet_involved_in_collision_start_per_thread[i][j].data());
			//std::cout << "tet tet_involved_in_collision_start_per_thread:  ";
			//for (int tt = 0; tt < thread_num; ++tt) {
			//	std::cout << tet_involved_in_collision_start_per_thread[i][j][tt + 1] << " ";
			//}
			//std::cout << std::endl;
		}
	}
}



void Collision::setPairStartPerThread(unsigned int** prefix_sum_record, unsigned int**element_index_start_per_thread,
	int* pair_per_thread_start_index)
{
	std::vector<int> start_per_thread(thread_num + 1, 0);
	arrangeIndex(thread_num, 
		prefix_sum_record[total_obj_num-1][element_index_start_per_thread[total_obj_num-1][thread_num]], start_per_thread);
	
	int start_index;

	for (unsigned int j =0; j < total_obj_num; ++j) {
		for (int k = 0; k < element_index_start_per_thread[j][thread_num]; ++k) {
			if (prefix_sum_record[j][k] != 0) {
				if (k != 0) {
					pair_per_thread_start_index[0] = j;
					pair_per_thread_start_index[1] = k - 1;
					pair_per_thread_start_index[2] = 0;
				}
				else {
					pair_per_thread_start_index[0] = j-1;
					pair_per_thread_start_index[1] = element_index_start_per_thread[j-1][thread_num];
					pair_per_thread_start_index[2] = 0;
				}
				goto finish_loop_in_setPairStartPerThread;
			}
		}
	}
	finish_loop_in_setPairStartPerThread:


	for (unsigned int i = 1; i <= thread_num; ++i) {
		for (unsigned int j = pair_per_thread_start_index[3*(i-1)]; j < total_obj_num; ++j) {
			if (j == pair_per_thread_start_index[3 * (i - 1)]) {
				start_index = pair_per_thread_start_index[3 * (i - 1) + 1];
				if (start_index == 0) {
					start_index = 1;
				}
			}
			else {
				start_index =1;
			}
			for (int k = start_index; k <= element_index_start_per_thread[j][thread_num]; ++k) {
				if (prefix_sum_record[j][k] == start_per_thread[i]) {
					if (j < total_obj_num - 1) {
						if (k != element_index_start_per_thread[j][thread_num]) {
							pair_per_thread_start_index[3 * i] = j;
							pair_per_thread_start_index[3 * i + 1] = k;
							pair_per_thread_start_index[3 * i + 2] = 0;							
						}
					}
					else{
						if (k != element_index_start_per_thread[j][thread_num]) {
							pair_per_thread_start_index[3 * i] = j;
							pair_per_thread_start_index[3 * i + 1] = k;
							pair_per_thread_start_index[3 * i + 2] = 0;
						}
						else {
							pair_per_thread_start_index[3 * i] = j;
							pair_per_thread_start_index[3 * i + 1] = k-1;
							pair_per_thread_start_index[3 * i + 2] = start_per_thread[i] - prefix_sum_record[j][pair_per_thread_start_index[3 * i + 1]];
						}
					}				
					goto next_loop;
				}
				else if (prefix_sum_record[j][k] > start_per_thread[i]) {
					if (k == 0) {
						pair_per_thread_start_index[3 * i] = j-1;
						pair_per_thread_start_index[3 * i + 1] = element_index_start_per_thread[j-1][thread_num]-1;
						pair_per_thread_start_index[3 * i + 2] = start_per_thread[i] - prefix_sum_record[j-1][pair_per_thread_start_index[3 * i + 1]];
					}
					else {
						pair_per_thread_start_index[3 * i] = j;
						pair_per_thread_start_index[3 * i + 1] = k - 1;
						pair_per_thread_start_index[3 * i + 2] = start_per_thread[i] - prefix_sum_record[j][k - 1];
					}					
					goto next_loop;
				}
			}
		}
	next_loop:;
	}
}




//SET_ELEMENT_COLLIDE_WITH_COLLIDER
void Collision::setElementCollideWithCollider(int thread_No)
{
	if (thread_num > 2) {
		switch (thread_No)
		{
		case 0:
			setTriangleCollideWithCollider();
			break;
		case 1:
			setEdgeCollideWithCollider();
			break;
		case 2:
			setVertexCollideWithCollider();
			break;
		}
	}
	else if (thread_num > 1) {
		switch (thread_No)
		{
		case 0:
			setTriangleCollideWithCollider();
			setVertexCollideWithCollider();
			break;
		case 1:
			setEdgeCollideWithCollider();
			break;
		}
	}
	else {
		setTriangleCollideWithCollider();
		setVertexCollideWithCollider();
		setEdgeCollideWithCollider();
	}
}

void Collision::findMinMaxDegreeOfCollisionPair(unsigned int& max_degree, unsigned int& min_degree)
{
	max_degree = 0;
	min_degree = UINT_MAX;
	//check VT
	unsigned int* vt_num;
	unsigned int* ee_num;
	unsigned int* vertex_pair;
	int cloth_size = cloth->size();
	
	unsigned int* vertex_index_on_surface_;

	int total_pair_num;
	int* triangle_index;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vt_num = vertex_triangle_pair_num_record[i];
		vertex_pair = vertex_triangle_pair_by_vertex[i];
		for (int j = 0; j < vertex_index_start_per_thread[i][thread_num]; ++i) {
			if (vt_num[j] != 0) {
				//add vt node
				total_pair_num = vt_num[j] - 1;
				for (int k = 0; k < vt_num[j]; k += 2) {
					triangle_index = triangle_indices[vertex_pair[k]][vertex_pair[k + 1]].data();
					if (vertex_pair[k] >= cloth_size) {
						total_pair_num += vertex_triangle_pair_num_record[vertex_pair[k]][general_index_to_surface_index[vertex_pair[k]][triangle_index[0]]];
						total_pair_num += vertex_triangle_pair_num_record[vertex_pair[k]][general_index_to_surface_index[vertex_pair[k]][triangle_index[1]]];
						total_pair_num += vertex_triangle_pair_num_record[vertex_pair[k]][general_index_to_surface_index[vertex_pair[k]][triangle_index[2]]];
					}
					else {
						total_pair_num += vertex_triangle_pair_num_record[vertex_pair[k]][triangle_index[0]];
						total_pair_num += vertex_triangle_pair_num_record[vertex_pair[k]][triangle_index[1]];
						total_pair_num += vertex_triangle_pair_num_record[vertex_pair[k]][triangle_index[2]];
					}					
				}				
				//add ee node

			}
		}
	}
}


void Collision::setTriangleCollideWithCollider()
{
	unsigned int* tri_pair;
	std::vector<unsigned int>* tri_index_reocrd;
	unsigned int size;
	for (int i = 0; i < total_obj_num; ++i) {
		size = mesh_struct[i]->triangle_indices.size();
		tri_index_reocrd = &triangle_index_collide_with_collider[i];
		tri_index_reocrd->clear();
		tri_pair = triangle_vertex_collider_pair_num_record[i];
		for (unsigned int j = 0; j < size;++j) {
			if (tri_pair[j] != 0) {
				tri_index_reocrd->emplace_back(j);
			}
		}
		arrangeIndex(thread_num, tri_index_reocrd->size(), triangle_index_collide_with_collider_start_per_thread.data() + i * (thread_num + 1));
	}
}

void Collision::setEdgeCollideWithCollider()
{
	unsigned int* edge_pair;
	std::vector<unsigned int>* edge_index_reocrd;
	unsigned int size;
	for (int i = 0; i < total_obj_num; ++i) {
		edge_index_reocrd = &edge_index_collide_with_collider[i];
		edge_index_reocrd->clear();
		edge_pair = edge_edge_collider_pair_num_record[i];
		size = mesh_struct[i]->edge_vertices.size() >> 1;
		for (unsigned int j = 0; j < size; ++j) {
			if (edge_pair[j] != 0) {
				edge_index_reocrd->emplace_back(j);
			}
		}
		arrangeIndex(thread_num, edge_index_reocrd->size(), edge_index_collide_with_collider_start_per_thread.data() + i * (thread_num + 1));
	}
}


void Collision::setVertexCollideWithCollider()
{
	unsigned int* vertex_pair;
	std::vector<unsigned int>* vertex_index_reocrd;
	unsigned int size;
	unsigned int* vertex_index_on_surface_;
	for (int i = 0; i < total_obj_num; ++i) {
		vertex_index_reocrd = &vertex_index_collide_with_collider[i];
		vertex_index_reocrd->clear();
		vertex_pair = vertex_obj_triangle_collider_num_record[i];
		size = vertex_index_start_per_thread[i][thread_num];
		if (i < cloth->size()) {
			for (unsigned int j = 0; j < size; ++j) {
				if (vertex_pair[j] != 0) {
					vertex_index_reocrd->emplace_back(j);
				}
			}
		}
		else {
			vertex_index_on_surface_ = vertex_index_on_surface[i];
			for (unsigned int j = 0; j < size; ++j) {
				if (vertex_pair[j] != 0) {
					vertex_index_reocrd->emplace_back(vertex_index_on_surface_[j]);
				}
			}
		}
		arrangeIndex(thread_num, vertex_index_reocrd->size(), vertex_index_collide_with_collider_start_per_thread.data() + i * (thread_num + 1));
	}
}

////FIND_PRIMITIVE_AROUND
//void Collision::findPointTriangleEdgeEdgePair(int thread_No)
//{
//	unsigned int total_triangle_pair_num = spatial_hashing.triangle_pair[thread_No][0];	
//	unsigned int* triangle_pair_ = spatial_hashing.triangle_pair[thread_No]+1;	
//	unsigned int* point_triangle_pair_= point_triangle_pair[thread_No] + 1;	
//
//	unsigned int* edge_edge_pair_ = edge_edge_pair[thread_No] + 1;
//	unsigned int obj0_index, triangle0_index, obj1_index, triangle1_index;
//	bool check_aabb;
//	unsigned int vertex_index;
//
//	unsigned int edge_edge_count_ = 0;
//	unsigned int vertex_triangle_count_ = 0;
//
//	for (int i = 0; i < total_triangle_pair_num; i += 4) {
//		triangle0_index = triangle_pair_[i];
//		obj0_index = triangle_pair_[i + 1];
//		triangle1_index = triangle_pair_[i + 2];
//		obj1_index = triangle_pair_[i + 3];
//		if (obj0_index != obj1_index) {
//			for (int j = 0; j < representative_vertex_num[obj0_index][triangle0_index]; ++j) {
//				//vertex_triangle_count_++;
//				if (AABB::AABB_intersection(vertex_aabb[obj0_index][triangle_index_in_order[obj0_index]->data()[3 * triangle0_index + j]].data(),
//					obj_tri_aabb[obj1_index][triangle1_index].data())) {
//					*(point_triangle_pair_++) = triangle_index_in_order[obj0_index]->data()[3 * triangle0_index + j];
//					memcpy(point_triangle_pair_, triangle_pair_ + i + 1, 12);
//					point_triangle_pair_ += 3;
//				}
//			}
//		}
//		else {
//			for (int j = 0; j < representative_vertex_num[obj0_index][triangle0_index]; ++j) {
//				//vertex_triangle_count_++;
//				if (AABB::AABB_intersection(vertex_aabb[obj0_index][triangle_index_in_order[obj0_index]->data()[3 * triangle0_index + j]].data(),
//					obj_tri_aabb[obj1_index][triangle1_index].data())) {
//					if (!vertexInTriangle(triangle_index_in_order[obj1_index][triangle1_index].data(),
//						triangle_index_in_order[obj0_index]->data()[3 * triangle0_index + j])) {
//						*(point_triangle_pair_++) = triangle_index_in_order[obj0_index]->data()[3 * triangle0_index + j];
//						memcpy(point_triangle_pair_, triangle_pair_ + i + 1, 12);
//						point_triangle_pair_ += 3;
//						//*(point_triangle_pair_++) = obj0_index;
//						//*(point_triangle_pair_++) = triangle1_index;
//						//*(point_triangle_pair_++) = obj1_index;
//					}
//				}
//			}
//		}
//	}
//	for (int i = 0; i < total_triangle_pair_num; i += 4) {
//		triangle0_index = triangle_pair_[i];
//		obj0_index = triangle_pair_[i + 1];
//		triangle1_index = triangle_pair_[i + 2];
//		obj1_index = triangle_pair_[i + 3];
//		if (obj0_index != obj1_index) {
//			for (int j = 0; j < representative_vertex_num[obj1_index][triangle1_index]; ++j) {
//				//vertex_triangle_count_++;
//				if (AABB::AABB_intersection(vertex_aabb[obj1_index][triangle_index_in_order[obj1_index]->data()[3 * triangle1_index + j]].data(),
//					obj_tri_aabb[obj0_index][triangle0_index].data())) {
//					*(point_triangle_pair_++) = triangle_index_in_order[obj1_index]->data()[3 * triangle1_index + j];
//					*(point_triangle_pair_++) = obj1_index;
//					*(point_triangle_pair_++) = triangle0_index;
//					*(point_triangle_pair_++) = obj0_index;
//				}
//			}
//		}
//		else {
//			for (int j = 0; j < representative_vertex_num[obj1_index][triangle1_index]; ++j) {
//				//vertex_triangle_count_++;
//				if (AABB::AABB_intersection(vertex_aabb[obj1_index][triangle_index_in_order[obj1_index]->data()[3 * triangle1_index + j]].data(),
//					obj_tri_aabb[obj0_index][triangle0_index].data())) {
//					if (!vertexInTriangle(triangle_index_in_order[obj0_index][triangle0_index].data(),
//						triangle_index_in_order[obj1_index]->data()[3 * triangle1_index + j])) {
//						*(point_triangle_pair_++) = triangle_index_in_order[obj1_index]->data()[3 * triangle1_index + j];
//						*(point_triangle_pair_++) = obj1_index;
//						*(point_triangle_pair_++) = triangle0_index;
//						*(point_triangle_pair_++) = obj0_index;
//					}
//				}
//			}
//		}
//	}
//	for (int i = 0; i < total_triangle_pair_num; i += 4) {
//		triangle0_index = triangle_pair_[i];
//		obj0_index = triangle_pair_[i + 1];
//		triangle1_index = triangle_pair_[i + 2];
//		obj1_index = triangle_pair_[i + 3];
//		if (obj0_index != obj1_index) {
//			for (int j = 0; j < representative_edge_num[obj0_index][triangle0_index]; ++j) {
//				for (int k = 0; k < representative_edge_num[obj1_index][triangle1_index]; ++k) {
//					//edge_edge_count_++;
//					if (AABB::AABB_intersection(edge_aabb[obj0_index][face_edges[obj0_index][3*triangle0_index+j]].data(),
//						edge_aabb[obj1_index][face_edges[obj1_index][3*triangle1_index+k]].data())) {
//						*(edge_edge_pair_++) = face_edges[obj0_index][3 * triangle0_index + j];
//						*(edge_edge_pair_++) = obj0_index;
//						*(edge_edge_pair_++) = face_edges[obj1_index][3 * triangle1_index + k];
//						*(edge_edge_pair_++) = obj1_index;
//					}
//				}
//			}
//		}
//		else {
//			for (int j = 0; j < representative_edge_num[obj0_index][triangle0_index]; ++j) {
//				for (int k = 0; k < representative_edge_num[obj1_index][triangle1_index]; ++k) {
//					//edge_edge_count_++;
//					if (AABB::AABB_intersection(edge_aabb[obj0_index][face_edges[obj0_index][3 * triangle0_index + j]].data(),
//						edge_aabb[obj1_index][face_edges[obj1_index][3 * triangle1_index + k]].data())) {
//						if (!edgeEdgeconnected(edge_vertices[obj0_index] + (face_edges[obj0_index][3 * triangle0_index + j] << 1),
//							edge_vertices[obj1_index] + (face_edges[obj1_index][3 * triangle1_index + k] << 1))) {
//							*(edge_edge_pair_++) = face_edges[obj0_index][3 * triangle0_index + j];
//							*(edge_edge_pair_++) = obj0_index;
//							*(edge_edge_pair_++) = face_edges[obj1_index][3 * triangle1_index + k];
//							*(edge_edge_pair_++) = obj1_index;
//						}
//					}
//				}
//			}
//		}	
//	}
//	point_triangle_pair[thread_No][0] = point_triangle_pair_ - point_triangle_pair[thread_No] - 1;
//	edge_edge_pair[thread_No][0] = edge_edge_pair_ - edge_edge_pair[thread_No] - 1;
//	//edge_edge_count[thread_No] = edge_edge_count_;
//	//vertex_triangle_count[thread_No] = vertex_triangle_count_;
//}


//FIND_PRIMITIVE_AROUND
void Collision::findPointTriangleEdgeEdgePair(int thread_No)
{
	//unsigned int total_triangle_pair_num = spatial_hashing.triangle_pair[thread_No][0];
	//unsigned int* triangle_pair_ = spatial_hashing.triangle_pair[thread_No] + 1;
	//unsigned int* point_triangle_pair_ = point_triangle_pair[thread_No] + 1;
	//unsigned int* edge_edge_pair_ = edge_edge_pair[thread_No] + 1;
	//unsigned int obj0_index, triangle0_index, obj1_index, triangle1_index;
	//bool check_aabb;
	//unsigned int vertex_index;
	//unsigned int edge_edge_count_ = 0;
	//unsigned int vertex_triangle_count_ = 0;
	//for (int i = 0; i < total_triangle_pair_num; i += 4) {
	//	triangle0_index = triangle_pair_[i];
	//	obj0_index = triangle_pair_[i + 1];
	//	triangle1_index = triangle_pair_[i + 2];
	//	obj1_index = triangle_pair_[i + 3];
	//	if (obj0_index != obj1_index) {
	//		for (int j = 0; j < representative_vertex_num[obj0_index][triangle0_index]; ++j) {
	//			//vertex_triangle_count_++;
	//			if (AABB::AABB_intersection(vertex_aabb[obj0_index][triangle_index_in_order[obj0_index]->data()[3 * triangle0_index + j]].data(),
	//				obj_tri_aabb[obj1_index][triangle1_index].data())) {
	//				*(point_triangle_pair_++) = triangle_index_in_order[obj0_index]->data()[3 * triangle0_index + j];
	//				memcpy(point_triangle_pair_, triangle_pair_ + i + 1, 12);
	//				point_triangle_pair_ += 3;
	//				//*(point_triangle_pair_++) = obj0_index;
	//				//*(point_triangle_pair_++) = triangle1_index;
	//				//*(point_triangle_pair_++) = obj1_index;
	//			}
	//		}
	//		for (int j = 0; j < representative_vertex_num[obj1_index][triangle1_index]; ++j) {
	//			//vertex_triangle_count_++;
	//			if (AABB::AABB_intersection(vertex_aabb[obj1_index][triangle_index_in_order[obj1_index]->data()[3 * triangle1_index + j]].data(),
	//				obj_tri_aabb[obj0_index][triangle0_index].data())) {
	//				*(point_triangle_pair_++) = triangle_index_in_order[obj1_index]->data()[3 * triangle1_index + j];
	//				*(point_triangle_pair_++) = obj1_index;
	//				*(point_triangle_pair_++) = triangle0_index;
	//				*(point_triangle_pair_++) = obj0_index;
	//				//_mm_stream_si32((int*)point_triangle_pair_++, *(int*)&triangle_index_in_order[obj1_index][triangle1_index][j]);
	//				//_mm_stream_si32((int*)point_triangle_pair_++, *(int*)&obj1_index);
	//				//_mm_stream_si32((int*)point_triangle_pair_++, *(int*)&triangle0_index);
	//				//_mm_stream_si32((int*)point_triangle_pair_++, *(int*)&obj0_index);
	//			}
	//		}
	//		for (int j = 0; j < representative_edge_num[obj0_index][triangle0_index]; ++j) {
	//			for (int k = 0; k < representative_edge_num[obj1_index][triangle1_index]; ++k) {
	//				//edge_edge_count_++;
	//				if (AABB::AABB_intersection(edge_aabb[obj0_index][face_edges[obj0_index][3 * triangle0_index + j]].data(),
	//					edge_aabb[obj1_index][face_edges[obj1_index][3 * triangle1_index + k]].data())) {
	//					*(edge_edge_pair_++) = face_edges[obj0_index][3 * triangle0_index + j];
	//					*(edge_edge_pair_++) = obj0_index;
	//					*(edge_edge_pair_++) = face_edges[obj1_index][3 * triangle1_index + k];
	//					*(edge_edge_pair_++) = obj1_index;
	//				}
	//			}
	//		}
	//	}
	//	else {
	//		for (int j = 0; j < representative_vertex_num[obj0_index][triangle0_index]; ++j) {
	//			//vertex_triangle_count_++;
	//			if (AABB::AABB_intersection(vertex_aabb[obj0_index][triangle_index_in_order[obj0_index]->data()[3 * triangle0_index + j]].data(),
	//				obj_tri_aabb[obj1_index][triangle1_index].data())) {
	//				if (!vertexInTriangle(triangle_index_in_order[obj1_index][triangle1_index].data(),
	//					triangle_index_in_order[obj0_index]->data()[3 * triangle0_index + j])) {
	//					*(point_triangle_pair_++) = triangle_index_in_order[obj0_index]->data()[3 * triangle0_index + j];
	//					memcpy(point_triangle_pair_, triangle_pair_ + i + 1, 12);
	//					point_triangle_pair_ += 3;
	//					//*(point_triangle_pair_++) = obj0_index;
	//					//*(point_triangle_pair_++) = triangle1_index;
	//					//*(point_triangle_pair_++) = obj1_index;
	//				}
	//			}
	//		}
	//		for (int j = 0; j < representative_vertex_num[obj1_index][triangle1_index]; ++j) {
	//			//vertex_triangle_count_++;
	//			if (AABB::AABB_intersection(vertex_aabb[obj1_index][triangle_index_in_order[obj1_index]->data()[3 * triangle1_index + j]].data(),
	//				obj_tri_aabb[obj0_index][triangle0_index].data())) {
	//				if (!vertexInTriangle(triangle_index_in_order[obj0_index][triangle0_index].data(),
	//					triangle_index_in_order[obj1_index]->data()[3 * triangle1_index + j])) {
	//					*(point_triangle_pair_++) = triangle_index_in_order[obj1_index]->data()[3 * triangle1_index + j];
	//					*(point_triangle_pair_++) = obj1_index;
	//					*(point_triangle_pair_++) = triangle0_index;
	//					*(point_triangle_pair_++) = obj0_index;
	//				}
	//			}
	//		}
	//		for (int j = 0; j < representative_edge_num[obj0_index][triangle0_index]; ++j) {
	//			for (int k = 0; k < representative_edge_num[obj1_index][triangle1_index]; ++k) {
	//				//edge_edge_count_++;
	//				if (AABB::AABB_intersection(edge_aabb[obj0_index][face_edges[obj0_index][3 * triangle0_index + j]].data(),
	//					edge_aabb[obj1_index][face_edges[obj1_index][3 * triangle1_index + k]].data())) {
	//					if (!edgeEdgeconnected(edge_vertices[obj0_index] + (face_edges[obj0_index][3 * triangle0_index + j] << 1),
	//						edge_vertices[obj1_index] + (face_edges[obj1_index][3 * triangle1_index + k] << 1))) {
	//						*(edge_edge_pair_++) = face_edges[obj0_index][3 * triangle0_index + j];
	//						*(edge_edge_pair_++) = obj0_index;
	//						*(edge_edge_pair_++) = face_edges[obj1_index][3 * triangle1_index + k];
	//						*(edge_edge_pair_++) = obj1_index;
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
	//if (has_collider) {
	//	unsigned int total_triangle_pair_num_with_collider = spatial_hashing.triangle_pair_with_collider[thread_No][0];
	//	unsigned int* point_obj_triangle_collider_pair_ = point_obj_triangle_collider_pair[thread_No] + 1;
	//	unsigned int* triangle_pair_collider_ = spatial_hashing.triangle_pair_with_collider[thread_No] + 1;
	//	unsigned int* edge_edge_collider_pair_ = edge_edge_collider_pair[thread_No] + 1;
	//	unsigned int* point_collider_triangle_obj_pair_ = point_collider_triangle_obj_pair[thread_No] + 1;
	//	for (int i = 0; i < total_triangle_pair_num_with_collider; i += 4) {
	//		triangle0_index = triangle_pair_collider_[i];
	//		obj0_index = triangle_pair_collider_[i + 1];
	//		triangle1_index = triangle_pair_collider_[i + 2];
	//		obj1_index = triangle_pair_collider_[i + 3];
	//		for (int j = 0; j < representative_vertex_num[obj0_index][triangle0_index]; ++j) {
	//			if (AABB::AABB_intersection(vertex_aabb[obj0_index][triangle_index_in_order[obj0_index][triangle0_index][j]].data(),
	//				obj_tri_aabb_collider[obj1_index][triangle1_index].data())) {
	//				*(point_obj_triangle_collider_pair_++) = triangle_index_in_order[obj0_index][triangle0_index][j];
	//				memcpy(point_obj_triangle_collider_pair_, triangle_pair_collider_ + i + 1, 12);
	//				//*(point_obj_triangle_collider_pair_++) = obj0_index;
	//				//*(point_obj_triangle_collider_pair_++) = triangle1_index;
	//				//*(point_obj_triangle_collider_pair_++) = obj1_index;
	//				point_obj_triangle_collider_pair_ += 3;
	//			}
	//		}
	//		for (int j = 0; j < representative_vertex_num_collider[obj1_index][triangle1_index]; ++j) {
	//			if (AABB::AABB_intersection(vertex_aabb_collider[obj1_index][triangle_index_in_order_collider[obj1_index][triangle1_index][j]].data(),
	//				obj_tri_aabb[obj0_index][triangle0_index].data())) {
	//				*(point_collider_triangle_obj_pair_++) = triangle_index_in_order_collider[obj1_index][triangle1_index][j];
	//				*(point_collider_triangle_obj_pair_++) = obj1_index;
	//				*(point_collider_triangle_obj_pair_++) = triangle0_index;
	//				*(point_collider_triangle_obj_pair_++) = obj0_index;
	//			}
	//		}
	//		for (int j = 0; j < representative_edge_num[obj0_index][triangle0_index]; ++j) {
	//			for (int k = 0; k < representative_edge_num_collider[obj1_index][triangle1_index]; ++k) {
	//				if (AABB::AABB_intersection(edge_aabb[obj0_index][face_edges[obj0_index][3 * triangle0_index + j]].data(),
	//					edge_aabb_collider[obj1_index][collider_face_edges[obj1_index][3 * triangle1_index + k]].data())) {
	//					*(edge_edge_collider_pair_++) = face_edges[obj0_index][3 * triangle0_index + j];
	//					*(edge_edge_collider_pair_++) = obj0_index;
	//					*(edge_edge_collider_pair_++) = collider_face_edges[obj1_index][3 * triangle1_index + k];
	//					*(edge_edge_collider_pair_++) = obj1_index;
	//				}
	//			}
	//		}
	//	}
	//	edge_edge_collider_pair[thread_No][0] = edge_edge_collider_pair_ - edge_edge_collider_pair[thread_No] - 1;
	//	point_obj_triangle_collider_pair[thread_No][0] = point_obj_triangle_collider_pair_ - point_obj_triangle_collider_pair[thread_No] - 1;
	//	point_collider_triangle_obj_pair[thread_No][0] = point_collider_triangle_obj_pair_ - point_collider_triangle_obj_pair[thread_No] - 1;
	//}
	//point_triangle_pair[thread_No][0] = point_triangle_pair_ - point_triangle_pair[thread_No] - 1;
	//edge_edge_pair[thread_No][0] = edge_edge_pair_ - edge_edge_pair[thread_No] - 1;
	////edge_edge_count[thread_No] = edge_edge_count_;
	////vertex_triangle_count[thread_No] = vertex_triangle_count_;
}
