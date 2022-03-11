#include"collision.h"

void Collision::initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider,
	std::vector<Tetrahedron>* tetrahedron, Thread* thread, double* tolerance_ratio)
{
	collision_time = 1.0;
	has_collider = !collider->empty();
	tetrahedron_begin_obj_index = cloth->size();
	total_obj_num = cloth->size() + tetrahedron->size();
	total_obj_with_collider = total_obj_num + collider->size();
	this->cloth = cloth;
	this->collider = collider;
	this->tetrahedron = tetrahedron;
	this->thread = thread;

	max_index_number_in_one_cell = 600;
	max_index_number_in_one_cell_collider = 300;
	estimate_coeff_for_pair_num = 20;

	draw_culling.initial(cloth, collider, tetrahedron, thread);
	//findPatchOfObjects();
	initialBVH(cloth, collider, tetrahedron, thread);
	initialTargetPos(cloth, tetrahedron, thread);
	initialSpatialHashing(cloth, collider, tetrahedron, thread, tolerance_ratio);
	//use_BVH = true;
	initialNeighborPrimitive();
	collision_time_thread.resize(thread->thread_num);

	initialCollidePairInfo();
	reorganzieDataOfObjects();
	//collision_constraint.testPT();
	//approx_CCD.test();
	//CCD::test();
	//std::cout <<"floor coordinate "<< (*collider)[0].ori_vertices[0][1] << std::endl;

	//draw_culling.setInSpatialHashingValue(spatial_hashing.spatial_hashing_value,
	//	spatial_hashing.spatial_hashing_triangle_index, spatial_hashing.spatial_hashing_value_collider,
	//	spatial_hashing.spatial_hashing_triangle_index_collider,
	//	&spatial_hashing.prefix_sum, &spatial_hashing.prefix_sum_collider, &spatial_hashing.cell_begin_per_thread);
	//the above last input variable should be actual_exist_cell_begin_per_thread(sorting) /cell_begin_per_thread(unsorting)
	//draw_culling.setInSpatialHashingValue(spatial_hashing.spatial_hashing_cell, spatial_hashing.spatial_hashing_cell_collider,
	//	spatial_hashing.hash_cell_count);
	draw_culling.vertex_tet_pair = spatial_hashing.vertex_tet_pair.data();


	edge_edge_count.resize(thread_num, 0);
	vertex_triangle_count.resize(thread_num, 0);

}


void Collision::initialCollidePairInfo()
{
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



	//for (unsigned int i = 0; i < thread_num; ++i) {
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

	for (unsigned int i = 0; i < cloth->size(); ++i) {
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
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
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

		for (unsigned int i = 0; i < collider->size(); ++i) {
			obj_tri_aabb_collider[i] = collider->data()[i].triangle_AABB.data();
			vertex_aabb_collider[i] = collider->data()[i].vertex_AABB.data();
			edge_aabb_collider[i] = collider->data()[i].edge_AABB.data();
			collider_face_edges[i] = collider->data()[i].mesh_struct.face_edges.data();
			collider_edge_vertices[i] = collider->data()[i].mesh_struct.edge_vertices.data();
			representative_vertex_num_collider[i] = collider->data()[i].representative_vertex_num.data();
			representative_edge_num_collider[i] = collider->data()[i].representative_edge_num.data();
			triangle_index_in_order_collider[i] = collider->data()[i].mesh_struct.surface_triangle_index_in_order.data();
		}
	}

}

void Collision::initialDHatTolerance(double ave_edge_length)
{
	d_hat = 1e-2 * ave_edge_length;
	eta = 0.01;
	tolerance = 1e-3 * d_hat;
	tolerance_2 = tolerance * tolerance;
	//std::cout << "d_hat_2 " << d_hat_2 << std::endl;

}


void Collision::findPatchOfObjects()
{
	//getAABBWithoutTolerance();
	//mesh_patch.initialPatch(cloth, collider, tetrahedron, thread);
	//mesh_patch.setBuffer(0,tetrahedron_begin_obj_index);

}

//void Collision::drawMeshPatch(Camera* camera)
//{
//	mesh_patch.draw(camera);
//}

void Collision::initialNeighborPrimitive()
{
	for (int i = 0; i < cloth->size(); ++i) {
		(*cloth)[i].initialNeighborPrimitiveRecording(cloth->size(), tetrahedron->size(), collider->size(), true);
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		(*tetrahedron)[i].initialNeighborPrimitiveRecording(cloth->size(), tetrahedron->size(), collider->size(), true);
	}
}


//void Collision::testBVHUpdate()
//{
	//BVH test_bvh;
	//test_bvh.init((*tetrahedron)[0].mesh_struct.triangle_indices.size(), (*tetrahedron)[0].mesh_struct.face_index_begin_per_thread, thread);	
	//test_bvh.updateBVH(&(*tetrahedron)[0].triangle_AABB);
	//int j = 0;
	////for (int i = 0; i < test_bvh.triangle_node_index.size(); ++i) {
	////	std::cout << test_bvh.triangle_node_index[i] << std::endl;
	////}
	////for (int i = 0; i < (*tetrahedron)[0].mesh_struct.triangle_indices.size(); ++i) {
	////	if (!(test_bvh.aabb_list[test_bvh.triangle_node_index[i]] == obj_BVH[0].aabb_list[test_bvh.triangle_node_index[i]])) {
	////	std::cout << test_bvh.aabb_list[i].max[0] << " " << test_bvh.aabb_list[i].max[1] << " " << test_bvh.aabb_list[i].max[2] << " "
	////		<< test_bvh.aabb_list[i].min[0] << " " << test_bvh.aabb_list[i].min[1] << " " << test_bvh.aabb_list[i].min[2] << std::endl;
	////	std::cout << obj_BVH[0].aabb_list[i].max[0] << " " << obj_BVH[0].aabb_list[i].max[1] << " " << obj_BVH[0].aabb_list[i].max[2] << " "
	////		<< obj_BVH[0].aabb_list[i].min[0] << " " << obj_BVH[0].aabb_list[i].min[1] << " " << obj_BVH[0].aabb_list[i].min[2] << std::endl;
	////		j++;
	////	}
	////}
	//for (int i = 1; i < test_bvh.aabb_list.size(); ++i) {
	//	if (!(test_bvh.aabb_list[i] == obj_BVH[0].aabb_list[i])) {
	//		std::cout << i << std::endl;
	//		std::cout << test_bvh.aabb_list[i].max[0] << " " << test_bvh.aabb_list[i].max[1] << " " << test_bvh.aabb_list[i].max[2] << " "
	//			<< test_bvh.aabb_list[i].min[0] << " " << test_bvh.aabb_list[i].min[1] << " " << test_bvh.aabb_list[i].min[2] << std::endl;
	//		std::cout << obj_BVH[0].aabb_list[i].max[0] << " " << obj_BVH[0].aabb_list[i].max[1] << " " << obj_BVH[0].aabb_list[i].max[2] << " "
	//			<< obj_BVH[0].aabb_list[i].min[0] << " " << obj_BVH[0].aabb_list[i].min[1] << " " << obj_BVH[0].aabb_list[i].min[2] << std::endl;
	//		j++;
	//	}
	//}
	//std::cout <<j<<" "<< test_bvh.aabb_list.size() << " run test bvh update" << std::endl;
	//for (unsigned int i = 0; i < cloth->size(); ++i) {
	//	obj_BVH[i].test(cloth->data()[0].triangle_AABB.data());
	//}
	//for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
	//	obj_BVH[i - tetrahedron_begin_obj_index].test(tetrahedron->data()[0].triangle_AABB.data());
	//}
//}

void Collision::initialTargetPos(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, Thread* thread)
{
	thread_num = thread->thread_num;
	obj_target_pos_per_thread.resize(thread_num);
	for (int i = 0; i < obj_target_pos_per_thread.size(); ++i) {
		obj_target_pos_per_thread[i].initialSet(total_obj_num);
		for (int j = 0; j < total_obj_num; ++j) {
			if (j < tetrahedron_begin_obj_index) {
				obj_target_pos_per_thread[i].initialSet2(j, (*cloth)[j].ori_vertices.size());
			}
			else {
				obj_target_pos_per_thread[i].initialSet2(j, (*tetrahedron)[j - tetrahedron_begin_obj_index].mesh_struct.vertex_index_on_sureface.size());
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
			obj_target_pos.initialSet2(j, (*tetrahedron)[j - tetrahedron_begin_obj_index].mesh_struct.vertex_index_on_sureface.size());
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
	spatial_hashing.setInObject(cloth, collider, tetrahedron, thread, tolerance_ratio, 8, false, 
		max_index_number_in_one_cell, max_index_number_in_one_cell_collider, estimate_coeff_for_pair_num);
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
//std::cout << "build " << clock() - t1 << std::endl;
//t1 = clock();
//for (int i = 0; i < 100; ++i) {
//thread->assignTask(this, FIND_TRIANGLE_PAIRS);
//}
//std::cout << "find triangle pair " << clock() - t1 << std::endl;

////testCollision();

//t1 = clock();
//for (int i = 0; i < 100; ++i) {
	//thread->assignTask(this, FIND_PRIMITIVE_AROUND);

	//}
	//std::cout << "find around primitive " << clock() - t1 << std::endl;
	////std::cout << "search " << clock() - t1 << std::endl;




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

	for (unsigned int i = 0; i < thread_num; ++i) {
		vertex_triangle_count_total_final += spatial_hashing.vertex_triangle_pair[i][0] >> 2;
		edge_edge_count_count_total_final += spatial_hashing.edge_edge_pair[i][0] >> 2;
		//triangle_triangle_count_total += spatial_hashing.triangle_pair[i][0] >> 2;
	}

	//std::cout << "triangle triangle pair " << triangle_triangle_count_total << std::endl;
	std::cout << " vertex triangle after cut " << vertex_triangle_count_total_final << std::endl;
	std::cout << " edge edge after cut " << edge_edge_count_count_total_final << std::endl;

}

void Collision::testRepeatability()
{
	unsigned int pair_num = 0;
	unsigned int edge_pair_num = 0;
	for (unsigned int i = 0; i < thread_num; ++i) {
		pair_num += spatial_hashing.vertex_triangle_pair[i][0] / 4;
		edge_pair_num += spatial_hashing.edge_edge_pair[i][0] / 4;
	}
	std::vector<TriangleElementPair> triangle_pair(pair_num);
	std::vector<TriangleElementPair> edge_pair(edge_pair_num);
	unsigned int index = 0;
	unsigned int index_edge = 0;

	unsigned int* tri_pair_;
	unsigned int* edge_pair_;

	for (unsigned int i = 0; i < thread_num; ++i) {
		tri_pair_ = spatial_hashing.vertex_triangle_pair[i] + 1;
		edge_pair_ = spatial_hashing.edge_edge_pair[i] + 1;
		for (unsigned int j = 0; j < spatial_hashing.vertex_triangle_pair[i][0]; j += 4) {
			memcpy(triangle_pair[index].index, tri_pair_ + j, 16);
			index++;
		}
		for (unsigned int j = 0; j < spatial_hashing.edge_edge_pair[i][0]; j += 4) {
			//if (edge_pair_[j] == 12106) {
			//	std::cout <<".."<< edge_pair_[j + 2] << std::endl;
			//}
			//else if (edge_pair_[j + 2] == 12106) {
			//	std::cout << edge_pair_[j] << std::endl;
			//}

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
	for (unsigned int i = 1; i < triangle_pair.size(); ++i) {
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
	for (unsigned int i = 1; i < edge_pair.size(); ++i) {
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
	for (unsigned int i = 0; i < count.size(); ++i) {
		if (count[i] == 2) {
			count_2++;
		}
		else if (count[i] > 2) {
			count_2_larger++;
		}
	}
	for (unsigned int i = 0; i < count_edge.size(); ++i) {
		if (count_edge[i] == 2) {
			count_2_collider++;
		}
		else if (count_edge[i] > 2) {
			count_2_larger_collider++;
		}
	}



	std::cout << "vertex triangle Repeatability 2: " << (double)count_2 / (double)count.size() << " larger than 2 " << (double)count_2_larger / (double)count.size() << std::endl;
	std::cout << count_2 << " " << count_2_larger << " " << count.size() << std::endl;
	std::cout << "total VT pair " << triangle_pair.size() << std::endl;
	std::cout << "edge edge Repeatability 2: " << (double)count_2_collider / (double)count_edge.size() << " larger than 2 " << (double)count_2_larger_collider / (double)count_edge.size() << std::endl;
	std::cout << count_2_collider << " " << count_2_larger_collider << " " << count_edge.size() << std::endl;

	std::cout << "total edge pair " << edge_pair.size() << std::endl;
	//std::vector<std::vector<unsigned int>> prfix_sum_vertex(total_obj_num);
	//std::vector<std::vector<unsigned int>> prfix_sum_edge(total_obj_num);	
	//for (unsigned int i = 0; i < total_obj_num; ++i) {
	//	if (i < cloth->size()) {
	//		prfix_sum_vertex[i].resize(cloth->data()[i].mesh_struct.vertex_position.size() + 1, 0);
	//		prfix_sum_edge[i].resize(cloth->data()[i].mesh_struct.edges.size() + 1, 0);
	//	}
	//	else {
	//		prfix_sum_vertex[i].resize(tetrahedron->data()[i- cloth->size()].mesh_struct.vertex_position.size() + 1, 0);
	//		prfix_sum_edge[i].resize(tetrahedron->data()[i- cloth->size()].mesh_struct.edges.size() + 1, 0);
	//	}
	//}

	//for (unsigned int i = 0; i < triangle_pair.size(); ++i) {
	//	prfix_sum_vertex[triangle_pair[i].index[1]][triangle_pair[i].index[0] + 1]++;
	//}

	//for (unsigned int i = 0; i < edge_pair.size(); ++i) {
	//	prfix_sum_edge[edge_pair[i].index[1]][edge_pair[i].index[0] + 1]++;
	//}
	////for (unsigned int i = 0; i < total_obj_num; ++i) {
	//	for (unsigned int j = 1; j < prfix_sum_vertex[0].size() + 1; ++j) {
	//		prfix_sum_vertex[0][j] += prfix_sum_vertex[0][j - 1];
	//	}
	//	for (unsigned int j = 1; j < prfix_sum_edge[0].size() + 1; ++j) {
	//		prfix_sum_edge[0][j] += prfix_sum_edge[0][j - 1];
	//	}
	////}
	//	for (unsigned int i = 1; i < total_obj_num; ++i) {
	//		prfix_sum_vertex[i][0] = prfix_sum_vertex[i - 1].back();
	//		prfix_sum_edge[i][0] = prfix_sum_edge[i - 1].back();

	//		for (unsigned int j = 1; j < prfix_sum_vertex[i].size() + 1; ++j) {
	//			prfix_sum_vertex[i][j] += prfix_sum_vertex[i][j - 1];
	//		}
	//		for (unsigned int j = 1; j < prfix_sum_edge[i].size() + 1; ++j) {
	//			prfix_sum_edge[i][j] += prfix_sum_edge[i][j - 1];
	//		}
	//	}



	//int surface_vertex_index_on_global;
	//std::vector<std::vector<std::vector<int>>>* tri_obj_tri;
	//std::vector<std::vector<std::vector<int>>>* edge_obj_edge;
	//bool need_break;
	//unsigned int tri_pair_num_ = 0;
	//unsigned int edge_pair_num_ = 0;
	//for (unsigned int m = cloth->size(); m < total_obj_num; ++m) { //obj 0
	//	if (m < cloth->size()) {
	//		tri_obj_tri = &cloth->data()[m].vertex_neighbor_obj_triangle;
	//	}
	//	else {
	//		tri_obj_tri = &tetrahedron->data()[m-cloth->size()].surface_vertex_neighbor_obj_triangle;
	//	}		
	//	for (int i = 0; i < tri_obj_tri->size(); ++i) {		//vertex 0			
	//		surface_vertex_index_on_global = tetrahedron->data()[m - cloth->size()].mesh_struct.vertex_index_on_sureface[i];
	//		for (int j = 0; j < tri_obj_tri->data()[i].size(); ++j) {		//obj 1		
	//			for (int k = 0; k < tri_obj_tri->data()[i][j].size(); ++k) {  // triangle 1
	//				tri_pair_num_++;
	//				//obj_index, i, j, tri_obj->data()[i][j][k]
	//				//need_break = false;
	//				//for (unsigned int l = prfix_sum_vertex[m][surface_vertex_index_on_global]; l < prfix_sum_vertex[m][surface_vertex_index_on_global + 1]; ++l) {
	//				//	if (triangle_pair[l].index[3] == j && triangle_pair[l].index[2] == tri_obj_tri->data()[i][j][k]) {
	//				//		need_break = true;
	//				//		break;
	//				//	}
	//				//}
	//				//if (!need_break) {
	//				//	std::cout << "does not find VT pair: " << m << " " << surface_vertex_index_on_global << " " << j << " " << tri_obj_tri->data()[i][j][k] << std::endl;
	//				//}
	//			}
	//		}
	//	}
	//	edge_obj_edge = &tetrahedron->data()[m - cloth->size()].edge_neighbor_obj_edge;
	//	for (int i = 0; i < edge_obj_edge->size(); ++i) {		//edge 0			
	//		for (int j = 0; j < edge_obj_edge->data()[i].size(); ++j) {		//obj 1		
	//			for (int k = 0; k < edge_obj_edge->data()[i][j].size(); ++k) {  // edge 1
	//				//obj_index, i, j, tri_obj->data()[i][j][k]
	//				edge_pair_num_++;
	//				//need_break = false;
	//				//for (unsigned int l = prfix_sum_edge[m][i]; l < prfix_sum_edge[m][i + 1]; ++l) {
	//				//	if (edge_pair[l].index[3] == j && edge_pair[l].index[2] == edge_obj_edge->data()[i][j][k]) {
	//				//		need_break = true;
	//				//		break;
	//				//	}
	//				//}
	//				//if (!need_break) {
	//				//	std::cout << "does not find EE pair: " << m << " " << i << " " << j << " " << edge_obj_edge->data()[i][j][k] << std::endl;
	//				//}
	//			}
	//		}
	//	}
	//}
	//std::cout << "edge pair num " << edge_pair_num_ << std::endl;
	//std::cout << "triangle pair num " << tri_pair_num_ << std::endl;
	//std::cout << "SH is correct " << std::endl;

	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	for (unsigned int j = 0; j < spatial_hashing.triangle_pair[i].size(); j+=4) {
	//		if (spatial_hashing.triangle_pair[i][j] == 1 || spatial_hashing.triangle_pair[i][j+2]==1) {
	//			std::cout <<i<<" "<< spatial_hashing.triangle_pair[i][j + 1] << " " << spatial_hashing.triangle_pair[i][j] << " "
	//				<< spatial_hashing.triangle_pair[i][j + 3] << " " << spatial_hashing.triangle_pair[i][j + 2] << std::endl;
	//		}
	//	}
	//}
	//std::cout << cloth->data()[0].triangle_neighbor_obj_triangle[0][0].size() << std::endl;
	//for (unsigned int i = 0; i < cloth->data()[0].triangle_neighbor_obj_triangle[0][0].size(); ++i) {
	//	std::cout << cloth->data()[0].triangle_neighbor_obj_triangle[0][0][i]<<" " << std::endl;
	//}
	//std::cout <<(int) AABB::AABB_intersection(obj_tri_aabb[0][0].data(), obj_tri_aabb[0][6].data()) << std::endl;

	

	//for (unsigned int i = 0; i < triangle_pair.size(); ++i) {
	//	std::cout << triangle_pair[i].index[1] << " " << triangle_pair[i].index[0] << " " << triangle_pair[i].index[3] << " " << triangle_pair[i].index[2] << std::endl;
	//}
	//TriangleElementPair a;
	//std::vector<unsigned int> count;
	//std::vector<unsigned int> count_collider;
	//a = triangle_pair[0];
	//count.push_back(1);
	//for (unsigned int i = 1; i < triangle_pair.size(); ++i) {
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
	//	for (unsigned int i = 1; i < triangle_pair_collider.size(); ++i) {
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
	//for (unsigned int i = 0; i < count.size(); ++i) {
	//	if (count[i] == 2) {
	//		count_2++;
	//	}
	//	else if (count[i] > 2) {
	//		count_2_larger++;
	//	}
	//}
	//if (!triangle_pair_collider.empty()) {
	//	for (unsigned int i = 0; i < count_collider.size(); ++i) {
	//		if (count_collider[i] == 2) {
	//			count_2_collider++;
	//		}
	//		else if (count_collider[i] > 2) {
	//			count_2_larger_collider++;
	//		}
	//	}
	//}
	//std::cout << "Repeatability 2: " << (double)count_2 / (double)count.size() << " larger than 2 " << (double)count_2_larger / (double)count.size() << std::endl;
	//std::cout << count_2 << " " << count_2_larger << " " << count.size() << std::endl;
	//std::cout << "Collider Repeatability 2: " << (double)count_2_collider / (double)count_collider.size() << " " << (double)count_2_larger_collider / (double)count_collider.size() << std::endl;


}


void Collision::testIfSPRight()
{
	unsigned int* vertex_tri_pair;
	for (unsigned int i = 0; i < thread_num; ++i) {
		vertex_tri_pair = spatial_hashing.vertex_triangle_pair[i] + 1;
		for (unsigned int j = 0; j < vertex_tri_pair[0]; j += 4) {
			findVertexTriangleInBVH(vertex_tri_pair[j + 1], vertex_tri_pair[j], vertex_tri_pair[j + 3], vertex_tri_pair[j + 2]);
		}
	}
	unsigned int* edge_edge_pair;
	for (unsigned int i = 0; i < thread_num; ++i) {
		edge_edge_pair = spatial_hashing.edge_edge_pair[i] + 1;
		for (unsigned int j = 0; j < edge_edge_pair[0]; j += 4) {
			findEdgeEdgeInBVH(edge_edge_pair[j + 1], edge_edge_pair[j], edge_edge_pair[j + 3], edge_edge_pair[j + 2]);
		}
	}
	

	//std::vector<unsigned int>* tri_tri_pair = spatial_hashing.triangle_pair;
	//std::vector<unsigned int>* tri_tri_pair_collider = spatial_hashing.triangle_pair_with_collider;
	//for (unsigned int i = 0; i < thread_num; ++i) {
	//	for (unsigned int j = 0; j < tri_tri_pair[i].size(); j+=4) {
	//		findInBVH(tri_tri_pair[i][j + 1], tri_tri_pair[i][j], tri_tri_pair[i][j + 3], tri_tri_pair[i][j + 2],false);
	//	}
	//	for (unsigned int j = 0; j < tri_tri_pair_collider[i].size(); j += 4) {
	//		findInBVH(tri_tri_pair_collider[i][j + 1], tri_tri_pair_collider[i][j], tri_tri_pair_collider[i][j + 3], tri_tri_pair_collider[i][j + 2], true);
	//	}
	//}

	//std::cout << "find pair in BVH " << std::endl;

	//std::vector<std::vector<std::vector<unsigned int>>>* tri_obj_tri;
	//for (unsigned int i = 0; i < total_obj_num; ++i) {
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
	////for (unsigned int i = 0; i < thread_num; ++i) {
	////	for (unsigned int j = 0; j < spatial_hashing.triangle_pair[i].size(); j+=4) {
	////		if (spatial_hashing.triangle_pair[i][j] == 1 || spatial_hashing.triangle_pair[i][j+2]==1) {
	////			std::cout <<i<<" "<< spatial_hashing.triangle_pair[i][j + 1] << " " << spatial_hashing.triangle_pair[i][j] << " "
	////				<< spatial_hashing.triangle_pair[i][j + 3] << " " << spatial_hashing.triangle_pair[i][j + 2] << std::endl;
	////		}
	////	}
	////}
	////std::cout << cloth->data()[0].triangle_neighbor_obj_triangle[0][0].size() << std::endl;
	////for (unsigned int i = 0; i < cloth->data()[0].triangle_neighbor_obj_triangle[0][0].size(); ++i) {
	////	std::cout << cloth->data()[0].triangle_neighbor_obj_triangle[0][0][i]<<" " << std::endl;
	////}
	////std::cout <<(int) AABB::AABB_intersection(obj_tri_aabb[0][0].data(), obj_tri_aabb[0][6].data()) << std::endl;
	std::cout << "test is right" << std::endl;


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
	//			for (unsigned int l = 0; l < thread_num; ++l) {
	//				for (unsigned int m = 0; m < triangle_pair[l].size(); m += 4) {
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
	//				std::cout << "does not find the collider pair: " << obj_index << " " << i << " " << j << " " << tri_obj->data()[i][j][k] << std::endl;
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
	//			for (unsigned int l = 0; l < thread_num; ++l) {
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
	//				std::cout << "does not find the pair: " << obj_index << " " << i << " " << j << " " << tri_obj->data()[i][j][k] << std::endl;
	//			}
	//		}
	//	}
	//}
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
	for (unsigned int i = 0; i < neighbor->size(); ++i) {
		if (neighbor->data()[i] == edge_index_1) {
			find = true;
		}
	}
	if (!find) {
		std::cout << "does not find edge edge pair: " << obj_0 << " " << edge_index_0 << " " << obj_1 << " " << edge_index_1 << std::endl;
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
		//std::cout << vertex_index_0 << " " << tetrahedron->data()[0].mesh_struct.vertex_surface_index.size() << std::endl;
		int vertex_index_in_surface = tetrahedron->data()[obj_0 - cloth->size()].mesh_struct.vertex_surface_index[vertex_index_0];
		//std::cout << vertex_index_in_surface << std::endl;
		neighbor = &tetrahedron->data()[obj_0 - cloth->size()].surface_vertex_neighbor_obj_triangle[vertex_index_in_surface][obj_1];
	}
	for (unsigned int i = 0; i < neighbor->size(); ++i) {
		if (neighbor->data()[i] == tri_index1) {
			find = true;
		}
	}
	if (!find) {
		std::cout << "does not find vertex triangle pair: " << obj_0 << " " << vertex_index_0 << " " << obj_1 << " " << tri_index1 << std::endl;
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
		for (unsigned int i = 0; i < neighbor->size(); ++i) {
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
		for (unsigned int i = 0; i < neighbor->size(); ++i) {
			if (neighbor->data()[i] == tri_index_0) {
				find = true;
			}
		}
		if (!find) {
			std::cout << "does not find the pair: " << obj_0 << " " << tri_index_0 << " " << obj_1 << " " << tri_index1 << std::endl;
		}
	}
	else {
		if (obj_0 < cloth->size()) {
			neighbor = &cloth->data()[obj_0].triangle_neighbor_collider_triangle[tri_index_0][obj_1];
		}
		else {
			neighbor = &tetrahedron->data()[obj_0 - cloth->size()].triangle_neighbor_collider_triangle[tri_index_0][obj_1];
		}
		for (unsigned int i = 0; i < neighbor->size(); ++i) {
			if (neighbor->data()[i] == tri_index1) {
				find = true;
			}
		}
		if (!find) {
			std::cout << "does not find the collider pair: " << obj_0 << " " << tri_index_0 << " " << obj_1 << " " << tri_index1 << std::endl;
		}
	}

}

void Collision::testCulling()
{
	for (int i = 0; i < (*collider)[0].triangle_AABB.size(); ++i) {
		for (int j = 0; j < (*collider)[0].triangle_neighbor_obj_vertex[i][0].size(); ++j) {
			if ((*collider)[0].triangle_neighbor_obj_vertex[i][0][j] == 13) {
				std::cout << i << std::endl;
			}
		}
	}
}

void Collision::collisionCulling()
{
	time_t t = clock();
	time_t t1 = clock();
	t = clock();
	for (unsigned int i = 0; i < 100; ++i) {
	getAABB();
	//std::cout << "end here " << std::endl;
	getSceneAABB();
	//std::cout << "end here 2 " << std::endl;
	}
	t1 = clock();
	std::cout << "AABB " << t1 - t << std::endl;


	buildBVH();

	spatial_hashing.buildSpatialHashing(scene_aabb);

	//thread->assignTask(this, FIND_TRIANGLE_PAIRS);
	//t = clock();
	//for (unsigned int i = 0; i < 10; ++i) {
	//thread->assignTask(this, FIND_PRIMITIVE_AROUND);
	//}
	//t1 = clock();
	//std::cout << "find primitive around multi thread" << t1 - t << std::endl;
	//t = clock();
	//for (unsigned int i = 0; i < 10; ++i) {
	//	for (unsigned int j = 0; j < thread_num; ++j) {
	//		findPointTriangleEdgeEdgePair(j);
	//	}	
	//}
	//t1 = clock();
	//std::cout << "find primitive around single thread" << t1 - t << std::endl;

	//record_time1.push_back(t1 - t);
	//if (record_time.size() == 5) {
	//	time_t t_0 = 0, t_1=0;
	//	for (unsigned int i = 0; i < 5; ++i) {
	//		t_0 += record_time[i];
	//		t_1 += record_time1[i];
	//	}
	//	std::cout << "ave 5 AABB " << (double)t_0 / 5000.0 << std::endl;
	//	std::cout << "ave 5 find primitive around " << (double)t_1 / 500.0 << std::endl;
	//}
	//testIfSPRight();

	testRepeatability();


	//testCulling();
	//draw_culling.drawAABBIntersectBetweenObjects();
	//draw_culling.setCellData(&spatial_hashing.hash_value_for_test, spatial_hashing.cell_length, spatial_hashing.cell_number,
	//	spatial_hashing.scene_aabb);
	//draw_culling.setSingleCellData(spatial_hashing.cell_length, spatial_hashing.cell_number, spatial_hashing.scene_aabb);
	//draw_culling.setTetrahedronVertex();
	//totalCount();
}





void Collision::getSceneAABB()
{
	memset(scene_aabb + 3, 0xFE, 24); //set double to -5.31401e+303
	memset(scene_aabb, 0x7F, 24); //set double to 1.38242e+306
	double* aabb_;
	for (unsigned int i = 0; i < total_obj_with_collider; ++i) {
		if (i < tetrahedron_begin_obj_index) {
			aabb_ = cloth->data()[i].obj_aabb;
		}
		else if (i < total_obj_num) {
			aabb_ = tetrahedron->data()[i - tetrahedron_begin_obj_index].obj_aabb;
		}
		else {
			aabb_ = collider->data()[i - total_obj_num].obj_aabb;
		}
		for (unsigned int j = 0; j < 3; ++j) {
			if (scene_aabb[j] > aabb_[j]) {
				scene_aabb[j] = aabb_[j];
			}
		}
		for (unsigned int j = 3; j < 6; ++j) {
			if (scene_aabb[j] < aabb_[j]) {
				scene_aabb[j] = aabb_[j];
			}
		}
	}
	for (int i = 0; i < 3; ++i) {
		scene_aabb[i + 3] += 0.1;
		scene_aabb[i] -= 0.1;
	}
}

void Collision::globalCollisionTime()
{
	for (int i = 0; i < cloth->size(); ++i) {
		thread->assignTask(&(*cloth)[i].mesh_struct, FACE_NORMAL);
	}
	for (int i = 0; i < tetrahedron->size(); ++i) {
		thread->assignTask(&(*tetrahedron)[i].mesh_struct, FACE_NORMAL);
	}

	thread->assignTask(this, GLOBAL_COLLISION_TIME);
	//std::cout << "finish" << std::endl;
	collision_time = collision_time_thread[0];
	for (int i = 1; i < thread_num; ++i) {
		if (collision_time > collision_time_thread[i]) {
			collision_time = collision_time_thread[i];
		}
	}
	collision_time *= 0.9;

	if (collision_time > 1.0) {
		collision_time = 1.0;
	}

	//std::cout <<"collision time "<< collision_time << std::endl;
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

	//testIfBuildCollisionConstraint();
}


void Collision::testIfBuildCollisionConstraint()
{
	for (int i = 0; i < obj_target_pos.b_sum[0].size(); ++i) {
		if (obj_target_pos.need_update[0][i]) {
			//std::cout << "build collision constraint "<<i << std::endl;
		}
	}

}

//COLLISION_CONSTRAINT
void Collision::collisionConstraint(int thread_No)
{
	TargetPosition* target_pos = &obj_target_pos_per_thread[thread_No];
	target_pos->initial();
	TriangleMeshStruct* mesh_struct;
	std::vector<std::vector<int>>* neighbor_primitve;
	int index_begin;
	int index_end;
	double* mass;
	std::array<double, 3>* current_pos;
	std::array<double, 3>* initial_pos;
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		mesh_struct = &(*cloth)[cloth_No].mesh_struct;
		index_begin = mesh_struct->vertex_index_begin_per_thread[thread_No];
		index_end = mesh_struct->vertex_index_begin_per_thread[thread_No + 1];
		neighbor_primitve = (*cloth)[cloth_No].vertex_neighbor_obj_triangle.data();
		mass = mesh_struct->mass.data();
		current_pos = mesh_struct->vertex_position.data();
		initial_pos = mesh_struct->vertex_for_render.data();
		for (int i = index_begin; i < index_end; ++i) {
			pointSelfTriangleClose(neighbor_primitve[i].data(), initial_pos[i].data(), current_pos[i].data(), i, cloth_No, mass[i], target_pos);
		}
	}
	std::array<int, 3>* triangle_vertex_index;
	std::vector<double*> triangle_pos(3);
	std::array<double, 3>* triangle_normal;
	for (int collider_No = 0; collider_No < collider->size(); ++collider_No) {
		mesh_struct = &(*collider)[collider_No].mesh_struct;
		index_begin = mesh_struct->face_index_begin_per_thread[thread_No];
		index_end = mesh_struct->face_index_begin_per_thread[thread_No + 1];
		current_pos = mesh_struct->vertex_position.data();
		neighbor_primitve = (*collider)[collider_No].triangle_neighbor_obj_vertex.data();
		triangle_vertex_index = mesh_struct->triangle_indices.data();
		triangle_normal = mesh_struct->face_normal.data();
		for (int i = index_begin; i < index_end; ++i) {
			for (int j = 0; j < 3; ++j) {
				triangle_pos[j] = current_pos[triangle_vertex_index[i][j]].data();
			}
			pointColliderTriangleClose(triangle_vertex_index[i].data(), neighbor_primitve[i].data(), triangle_pos, triangle_normal[i].data(), target_pos);
		}
	}
	unsigned int* edge_vertex_index;
	double mass_[4];
	for (int cloth_No = 0; cloth_No < cloth->size(); ++cloth_No) {
		mesh_struct = &(*cloth)[cloth_No].mesh_struct;
		index_begin = mesh_struct->edge_index_begin_per_thread[thread_No];
		index_end = mesh_struct->edge_index_begin_per_thread[thread_No + 1];
		neighbor_primitve = (*cloth)[cloth_No].edge_neighbor_obj_edge.data();
		initial_pos = mesh_struct->vertex_for_render.data();
		current_pos = mesh_struct->vertex_position.data();
		for (int i = index_begin; i < index_end; ++i) {
			edge_vertex_index = mesh_struct->edge_vertices.data() + (i << 1);// edges[i].vertex;
			mass_[0] = mesh_struct->mass[edge_vertex_index[0]];
			mass_[1] = mesh_struct->mass[edge_vertex_index[1]];
			edgeEdgeClose(neighbor_primitve[i].data(), initial_pos[edge_vertex_index[0]].data(), initial_pos[edge_vertex_index[1]].data(),
				current_pos[edge_vertex_index[0]].data(), current_pos[edge_vertex_index[1]].data(), cloth_No, edge_vertex_index[0],
				edge_vertex_index[1], mass_, target_pos);
		}
	}
}




//GLOBAL_COLLISION_TIME
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
	int k = 0;
	for (int i = 0; i < (*cloth)[0].triangle_neighbor_obj_triangle.size(); ++i) {
		//(*cloth)[0].triangle_neighbor_cloth_triangle[i][0].push_back(i);
		k += (*cloth)[0].triangle_neighbor_obj_triangle[i][0].size();
	}
	//std::cout << k + (*cloth)[0].triangle_neighbor_cloth_triangle.size()<<std::endl;

	//for (int i = 0; i < (*cloth)[0].mesh_struct.vertices.size(); ++i) {
	//int i = 0;

	//for (int j = 0; j < (*cloth)[0].vertex_neighbor_collider_triangle[i][0].size(); ++j) {
	//	//std::cout << (*cloth)[0].vertex_AABB[i].min[0] << " " << (*cloth)[0].vertex_AABB[i].min[1] << " " << (*cloth)[0].vertex_AABB[i].min[2] << " "
	//		<< (*cloth)[0].vertex_AABB[i].max[0] << " " << (*cloth)[0].vertex_AABB[i].max[1] << " " << (*cloth)[0].vertex_AABB[i].max[2] << std::endl;
	//	//std::cout << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].min[0] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].min[1] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].min[2] << " "
	//		<< (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].max[0] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].max[1] << " " << (*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]].max[2] << std::endl;
	//	//std::cout << i << " " << (*cloth)[0].vertex_neighbor_collider_triangle[i][0][j] << std::endl;
	//	/*	if (!(*cloth)[0].vertex_AABB[i].AABB_intersection((*collider)[0].triangle_AABB[(*cloth)[0].vertex_neighbor_collider_triangle[i][0][j]])) {
	//			//std::cout << i << " " << (*cloth)[0].vertex_neighbor_collider_triangle[i][0][j] << std::endl;
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
	//				//std::cout << (*collider)[0].triangle_neighbor_obj_vertex[i][0][j] << " " << i << std::endl;
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
	//			//std::cout << i << " " << (*cloth)[0].triangle_neighbor_collider_triangle[i][0][j] << std::endl;
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

	int tet_index;
	for (int i = 0; i < thread_num; ++i) {
		for (int j = 0; j < total_obj_num; ++j) {
			obj_target_pos.collision_energy[j] += obj_target_pos_per_thread[i].collision_energy[j];
		}
	}
}

void Collision::sumTargetPosition()
{
	obj_target_pos.initial();
	int tet_index;
	for (int i = 0; i < thread_num; ++i) {
		for (int j = 0; j < total_obj_num; ++j) {
			obj_target_pos.collision_energy[j] += obj_target_pos_per_thread[i].collision_energy[j];
		}
	}

	thread->assignTask(this, SUM_TARGET_POSITION);

}

//RESUM_TARGET_POSITION
void Collision::resumTargetPositionPerThread(int thread_id)
{
	unsigned int* index_begin;
	bool* need_update;
	std::vector<std::array<double, 3>>* b_sum;
	std::vector<std::array<double, 3>>* b_sum_per_thread;

	for (unsigned int j = 0; j < total_obj_num; ++j) {
		b_sum = &obj_target_pos.b_sum[j];
		if (j < tetrahedron_begin_obj_index) {
			index_begin = (*cloth)[j].mesh_struct.vertex_index_begin_per_thread.data();
		}
		else {
			index_begin = (*tetrahedron)[j - tetrahedron_begin_obj_index].mesh_struct.vertex_index_on_surface_begin_per_thread.data();
		}
		for (unsigned int i = 0; i < thread_num; ++i) {
			need_update = obj_target_pos_per_thread[i].need_update[j];
			b_sum_per_thread = &obj_target_pos_per_thread[i].b_sum[j];
			for (unsigned int k = index_begin[thread_id]; k < index_begin[thread_id + 1]; ++k) {
				if (need_update[k]) {
					SUM((*b_sum)[k], (*b_sum)[k], (*b_sum_per_thread)[k]);
				}
			}
		}
	}



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

	for (unsigned int j = 0; j < total_obj_num; ++j) {
		global_need_update = obj_target_pos.need_update[j];
		b_sum = &obj_target_pos.b_sum[j];
		if (j < tetrahedron_begin_obj_index) {
			index_begin = (*cloth)[j].mesh_struct.vertex_index_begin_per_thread.data();
		}
		else {
			index_begin = (*tetrahedron)[j - tetrahedron_begin_obj_index].mesh_struct.vertex_index_on_surface_begin_per_thread.data();
		}
		global_stiffness = obj_target_pos.stiffness[j].data();
		for (unsigned int i = 0; i < thread_num; ++i) {
			need_update = obj_target_pos_per_thread[i].need_update[j];
			b_sum_per_thread = &obj_target_pos_per_thread[i].b_sum[j];
			stiffness = obj_target_pos_per_thread[i].stiffness[j].data();
			for (unsigned int k = index_begin[thread_id]; k < index_begin[thread_id + 1]; ++k) {
				if (need_update[k]) {
					global_need_update[k] = true;
					SUM_((*b_sum)[k], (*b_sum_per_thread)[k]);
					global_stiffness[k] += stiffness[k];
				}
			}
		}
	}
	////std::cout << "collision vertex ";
	//for (int i = 0; i < (*cloth)[0].ori_vertices.size(); ++i) {
	//	if (obj_target_pos.need_update[0][i]) {
	//		//std::cout << i << " ";
	//	}
	//}
	////std::cout << std::endl;
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
	//	for (unsigned int j = 0; j < total_obj_num; ++j) {
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
	//		for (unsigned int i = index_begin[thread_No]; i < index_begin[thread_No + 1]; ++i) {
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
	//	for (unsigned int j = 0; j < total_obj_num; ++j) {
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
	//		for (unsigned int i = index_begin[thread_No]; i < end; ++i) {
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
		record_stiffness = (*cloth)[i].collision_stiffness_initial[SELF_POINT_TRIANGLE];
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
				addTargetPosToSystemTotal(vertex_b_sum, target_position->collision_energy[cloth_No], initial_vertex_pos,
					target_pos, stiffness, *vertex_stiffness, *vetex_need_update);
				int triangle_vertex;
				for (int j = 0; j < 3; ++j) {
					triangle_vertex = triangle_vertex_index[j];
					addTargetPosToSystemTotal(target_position->b_sum[i][triangle_vertex].data(), target_position->collision_energy[i], initial_position[triangle_vertex].data(),
						triangle_target_pos[j].data(), stiffness,
						target_position->stiffness[i][triangle_vertex], target_position->need_update[i][triangle_vertex]);
					//if (vertex_index == 1) {
					//	//std::cout << triangle_initial_pos[j][0] << " " << triangle_initial_pos[j][1] << " " << triangle_initial_pos[j][2] << std::endl;
					//	//std::cout << triangle_current_pos[j][0] << " " << triangle_current_pos[j][1] << " " << triangle_current_pos[j][2] << std::endl;
					//	//std::cout << triangle_target_pos[j][0] << " " << triangle_target_pos[j][1] << " " << triangle_target_pos[j][2] << std::endl;
					//	//std::cout << initial_face_normal[triangle_index][0] << " " << initial_face_normal[triangle_index][1] << " " << initial_face_normal[triangle_index][2] << std::endl;
					//}
					//
					////std::cout << target_position->collision_energy[i] << std::endl;
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
		record_stiffness = (*cloth)[i].collision_stiffness_initial[SELF_EDGE_EDGE];
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
					target_position->collision_energy[cloth_No], initial_edge_vertex_0,
					target_pos[0].data(), stiffness, target_position->stiffness[cloth_No][edge_vertex_index_0], target_position->need_update[cloth_No][edge_vertex_index_0]);
				addTargetPosToSystemTotal(target_position->b_sum[cloth_No][edge_vertex_index_1].data(),
					target_position->collision_energy[cloth_No], initial_edge_vertex_1,
					target_pos[1].data(), stiffness, target_position->stiffness[cloth_No][edge_vertex_index_1], target_position->need_update[cloth_No][edge_vertex_index_1]);
				addTargetPosToSystemTotal(target_position->b_sum[i][compare_edge_index_0].data(),
					target_position->collision_energy[i], initial_position[compare_edge_index_0].data(),
					compare_target_pos[0].data(), stiffness, target_position->stiffness[i][compare_edge_index_0], target_position->need_update[i][compare_edge_index_0]);
				addTargetPosToSystemTotal(target_position->b_sum[i][compare_edge_index_1].data(),
					target_position->collision_energy[i], initial_position[compare_edge_index_1].data(),
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
		record_stiffness = (*cloth)[i].collision_stiffness_initial[BODY_POINT_TRIANGLE];
		energy = &target_position->collision_energy[i];
		for (int k = 0; k < neighbor_vertex->size(); ++k) {
			vertex_index = (*neighbor_vertex)[k];
			stiffness = record_stiffness;
			if (collision_constraint.pointColliderTriangle(vertex_mesh->vertex_for_render[vertex_index].data(), vertex_mesh->vertex_position[vertex_index].data(),
				current_position, current_face_normal, target_pos, d_hat, stiffness)) {
				addTargetPosToSystemTotal(vertex_b_sum[vertex_index].data(),
					*energy, vertex_mesh->vertex_for_render[vertex_index].data(),
					target_pos, stiffness,
					vertex_stiffness[vertex_index], vertex_need_update[vertex_index]);
				////std::cout << vertex_index << " " << stiffness << std::endl;
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
			//	//	std::cout << current_collision_time << " " << vertex_index << " " << triangle_vertex_index[0] << " " << triangle_vertex_index[1] << " " << triangle_vertex_index[2] << std::endl;
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
		//	std::cout << "k11" << std::endl;
		//}
		for (int k = 0; k < neighbor_vertex->size(); ++k) {
			//if (*time_stamp == 15 && triangle_index == 9104 && k==3)
			//{
			//	std::cout << vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]][0]<<" "<< vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]][1]<<" "
			//		<< vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]][2] << std::endl;
			//	std::cout << vertex_mesh->vertex_position[(*neighbor_vertex)[k]][0] << " " << vertex_mesh->vertex_position[(*neighbor_vertex)[k]][1] << " "
			//		<< vertex_mesh->vertex_position[(*neighbor_vertex)[k]][2] << std::endl;
			//	std::cout << initial_position[triangle_vertex_index[0]][0] << " " << initial_position[triangle_vertex_index[0]][1] << " " <<
			//		initial_position[triangle_vertex_index[0]][2] << std::endl;
			//	std::cout << initial_position[triangle_vertex_index[1]][0] << " " << initial_position[triangle_vertex_index[1]][1] << " " <<
			//		initial_position[triangle_vertex_index[1]][2] << std::endl;
			//	std::cout << initial_position[triangle_vertex_index[2]][0] << " " << initial_position[triangle_vertex_index[2]][1] << " " <<
			//		initial_position[triangle_vertex_index[2]][2] << std::endl;
			//	std::cout << current_position[triangle_vertex_index[0]][0] << " " << current_position[triangle_vertex_index[0]][1] << " " <<
			//		current_position[triangle_vertex_index[0]][2] << std::endl;
			//	std::cout << current_position[triangle_vertex_index[1]][0] << " " << current_position[triangle_vertex_index[1]][1] << " " <<
			//		current_position[triangle_vertex_index[1]][2] << std::endl;
			//	std::cout << current_position[triangle_vertex_index[2]][0] << " " << current_position[triangle_vertex_index[2]][1] << " " <<
			//		current_position[triangle_vertex_index[2]][2] << std::endl;
			//}
			current_collision_time = CCD::pointTriangleCcd(vertex_mesh->vertex_for_render[(*neighbor_vertex)[k]].data(), initial_position[triangle_vertex_index[0]].data(), initial_position[triangle_vertex_index[1]].data(), initial_position[triangle_vertex_index[2]].data(),
				vertex_mesh->vertex_position[(*neighbor_vertex)[k]].data(), current_position[triangle_vertex_index[0]].data(), current_position[triangle_vertex_index[1]].data(), current_position[triangle_vertex_index[2]].data(), eta, tolerance);
			if ((*collision_time) > current_collision_time) {
				(*collision_time) = current_collision_time;
			}
			//if (*time_stamp == 15 && triangle_index == 9104)
			//{
			//	std::cout << k << " " << (*neighbor_vertex)[k] << std::endl;
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
		//	std::cout << "k22" << std::endl;
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
		//	//std::cout << "++" << std::endl;
		//	//std::cout << vertex_index << " " << triangle_index << " " << vertex_mesh->vertex_for_render[vertex_index][0] << " " << vertex_mesh->vertex_for_render[vertex_index][1] << " " << vertex_mesh->vertex_for_render[vertex_index][2] << " " << std::endl;
		//	for (int i = 0; i < triangle_target_pos.size(); ++i) {
		//		//std::cout << vertex_index << " " << triangle_index << " " << initial_triangle_pos[i][0] << " " << initial_triangle_pos[i][1] << " " << initial_triangle_pos[i][2] << " " << std::endl;
		//	}
		//	//std::cout << vertex_mesh->vertex_position[vertex_index][0] << " " << vertex_mesh->vertex_position[vertex_index][1] << " " << vertex_mesh->vertex_position[vertex_index][2] << " " << std::endl;
		//	for (int i = 0; i < triangle_target_pos.size(); ++i) {
		//		//std::cout << current_triangle_pos[i][0] << " " << current_triangle_pos[i][1] << " " << current_triangle_pos[i][2] << " " << std::endl;
		//	}
		//	//std::cout << vertex_index << " " << triangle_index <<" "<< vertex_target_pos[0] << " " << vertex_target_pos[1] << " " << vertex_target_pos[2] << " " << std::endl;
		//	for (int i = 0; i < triangle_target_pos.size(); ++i) {
		//		//std::cout << vertex_index<<" "<< triangle_index << " " << triangle_target_pos[i][0] << " " << triangle_target_pos[i][1] << " " << triangle_target_pos[i][2] << " " << std::endl;
		//	}
		//	//std::cout << "++" << std::endl;
		//}

		addTargetPosToSystem(target_position->b_sum[vertex_cloth_No][vertex_index].data(),
			target_position->collision_energy[vertex_cloth_No], vertex_mesh->vertex_position[vertex_index].data(), vertex_target_pos, vertex_collision_stiffness);
		int triangle_vertex;
		for (int i = 0; i < 3; ++i) {
			triangle_vertex = triangle_vertex_index[i];
			addTargetPosToSystem(target_position->b_sum[triangle_cloth_No][triangle_vertex].data(),
				target_position->collision_energy[triangle_cloth_No], triangle_mesh->vertex_position[triangle_vertex].data(), triangle_target_pos[i].data(), triangle_collision_stiffness);
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
				target_position->collision_energy[edge_cloth_No], edge_mesh->vertex_position[edge_vertex_index[i]].data(), target_pos[i].data(), collision_stiffness);
			addTargetPosToSystem(target_position->b_sum[compare_cloth_No][compare_edge_vertex_index[i]].data(),
				target_position->collision_energy[compare_cloth_No], compare_mesh->vertex_position[compare_edge_vertex_index[i]].data(), compare_target_pos[i].data(), compare_collision_stiffness);
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
			//	//std::cout << vertex_target_pos[0] << " " << vertex_target_pos[1] << " " << vertex_target_pos[2] << " "
			//		<< vertex_mesh->vertex_position[vertex_index][0] << " " << vertex_mesh->vertex_position[vertex_index][1] << " " << vertex_mesh->vertex_position[vertex_index][2] << std::endl;
			//	//std::cout << current_triangle_pos[0][1] << std::endl;
			//}
		}
		else {
			/*if (vertex_index == 0) {
				//std::cout << vertex_target_pos[0] << " " << vertex_target_pos[1] << " " << vertex_target_pos[2] << " "
					<< vertex_mesh->vertex_position[vertex_index][0] << " " << vertex_mesh->vertex_position[vertex_index][1] << " " << vertex_mesh->vertex_position[vertex_index][2] << std::endl;
			}*/
		}

		addTargetPosToSystem(target_position->b_sum[vertex_cloth_No][vertex_index].data(),
			target_position->collision_energy[vertex_cloth_No], vertex_mesh->vertex_position[vertex_index].data(), vertex_target_pos, vertex_collision_stiffness);
		return true;
	}
	return false;
}

void Collision::addTargetPosToSystem(double* b_sum, double& energy, double* current_pos, double* target_pos, double stiffness)
{
	double temp[3];
	SUB(temp, current_pos, target_pos);
	energy += 0.5 * stiffness * DOT(temp, temp);
	b_sum[0] += target_pos[0] * stiffness;
	b_sum[1] += target_pos[1] * stiffness;
	b_sum[2] += target_pos[2] * stiffness;
}

void Collision::addTargetPosToSystemTotal(double* b_sum, double& energy, double* current_pos, double* target_pos, double stiffness, double& sum_stiffness, bool& update)
{
	double temp[3];
	SUB(temp, current_pos, target_pos);
	energy += 0.5 * stiffness * DOT(temp, temp);
	b_sum[0] += target_pos[0] * stiffness;
	b_sum[1] += target_pos[1] * stiffness;
	b_sum[2] += target_pos[2] * stiffness;
	update = true;
	sum_stiffness += stiffness;
}

void Collision::test()
{
	//findAllNeighborPairs();

	//std::cout << "vertex " << std::endl;
	//for (int i = 0; i < (*cloth)[0].vertex_neighbor_obj_triangle[0][0].size(); ++i) {
	//	//std::cout << (*cloth)[0].vertex_neighbor_obj_triangle[0][0][i] << " ";
	//}
	////std::cout << std::endl;
	////std::cout << "triangle " << std::endl;
	//for (int i = 0; i < (*cloth)[0].triangle_neighbor_cloth_triangle[0][0].size(); ++i) {
	//	//std::cout << (*cloth)[0].triangle_neighbor_cloth_triangle[0][0][i] << " ";
	//}
	////std::cout << std::endl;
	////std::cout << "edge " << std::endl;
	//for (int i = 0; i < (*cloth)[0].edge_neighbor_obj_edge[0][0].size(); ++i) {
	//	//std::cout << (*cloth)[0].edge_neighbor_obj_edge[0][0][i] << " ";
	//}
	//std::cout << std::endl;
}


void Collision::searchPatch(double* aabb, unsigned int compare_index, unsigned int obj_No, bool& intersect)
{
	//for (unsigned int i = 0; i < total_obj_num; ++i) {
	//	if (obj_BVH[i].searchIfPatchIntersect(aabb, compare_index, i == obj_No, 1, 0, mesh_patch.patch_AABB[i].size())) {
	//		intersect = true;
	//	}	
	//}
}

void Collision::searchTriangle(double* aabb, unsigned int compare_index, unsigned int obj_No, std::vector<std::vector<unsigned int>>* obj_neighbor_index,
	std::vector<std::vector<unsigned int>>* collider_neighbor_index, bool is_collider)
{

	for (unsigned int i = 0; i < total_obj_num; ++i) {
		(*obj_neighbor_index)[i].clear();
		if (i < tetrahedron_begin_obj_index) {
			obj_BVH[i].search(aabb, compare_index, i == obj_No, &((*obj_neighbor_index)[i]), 1, 0, (*cloth)[i].triangle_AABB.size());
		}
		else {
			obj_BVH[i].search(aabb, compare_index, i == obj_No, &((*obj_neighbor_index)[i]), 1, 0, (*tetrahedron)[i - tetrahedron_begin_obj_index].triangle_AABB.size());
		}
	}
	for (unsigned int i = 0; i < collider->size(); ++i) {
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
//	for (unsigned int i = 0; i < total_obj_with_collider; ++i) {
//		thread_begin = mesh_patch.patch_index_start_per_thread[i].data();
//		end = thread_begin[thread_No + 1];
//		for (unsigned int j = thread_begin[thread_No]; j < end; ++j) {
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
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		thread_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
		end = thread_begin[thread_No + 1];
		for (unsigned int j = thread_begin[thread_No]; j < end; ++j) {
			searchTriangle((*cloth)[i].triangle_AABB[j].data(), j, i, &(*cloth)[i].triangle_neighbor_obj_triangle[j],
				&(*cloth)[i].triangle_neighbor_collider_triangle[j], false);
		}
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		thread_begin = (*tetrahedron)[i].mesh_struct.face_index_begin_per_thread.data();
		end = thread_begin[thread_No + 1];
		for (unsigned int j = thread_begin[thread_No]; j < end; ++j) {
			searchTriangle((*tetrahedron)[i].triangle_AABB[j].data(), j, i, &(*tetrahedron)[i].triangle_neighbor_obj_triangle[j],
				&(*tetrahedron)[i].triangle_neighbor_collider_triangle[j], false);
		}
	}
	//else {
	//	std::vector<std::vector<int>> test_neighbor_triangle;
	//	test_neighbor_triangle.resize(cloth->size() + tetrahedron->size());
	//	for (unsigned int i = 0; i < cloth->size(); ++i) {
	//		thread_begin = (*cloth)[i].mesh_struct.face_index_begin_per_thread.data();
	//		end = thread_begin[thread_No + 1];
	//		for (unsigned int j = thread_begin[thread_No]; j < end; ++j) {
	//			spatial_hashing.searchTriangle((*cloth)[i].triangle_AABB[j].data(),i, j,(*cloth)[i].triangle_neighbor_obj_triangle[j].data(),
	//				false,thread_No);
	//			//searchTriangle((*cloth)[i].triangle_AABB[j], j, i, &test_neighbor_triangle,
	//			//	&test_neighbor_triangle, false);
	//			//testTwoVectorsAreSame((*cloth)[i].triangle_neighbor_obj_triangle[j], test_neighbor_triangle,i,j);
	//			
	//		}
	//	}
	//	for (unsigned int i = 0; i < collider->size(); ++i) {
	//		thread_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread.data();
	//		end = thread_begin[thread_No + 1];
	//		for (int j = thread_begin[thread_No]; j < end; ++j) {
	//			spatial_hashing.searchTriangle((*collider)[i].triangle_AABB[j].data(), i, j, (*collider)[i].triangle_neighbor_obj_triangle[j].data(),
	//				true, thread_No);			
	//			//searchTriangle((*collider)[i].triangle_AABB[j], j, i, &test_neighbor_triangle,
	//			//	&test_neighbor_triangle, true);
	//			//if (thread_No == 0) {
	//				//std::cout <<"after func"<< (*collider)[i].triangle_neighbor_obj_triangle[j][0].size() << std::endl;
	//				//testTwoVectorsAreSame((*collider)[i].triangle_neighbor_obj_triangle[j], test_neighbor_triangle, i, j);
	//			//}
	//		}
	//	}
	//	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
	//		thread_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread.data();
	//		end = thread_begin[thread_No + 1];
	//		obj_No = tetrahedron_begin_obj_index + i;
	//		for (unsigned int j = thread_begin[thread_No]; j < end; ++j) {
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
			std::cout << "error " << obj_index << " " << triangle_index << " size: " << vec1[i].size() << " " << vec2[i].size() << std::endl;
			if (!vec1[i].empty()) {
				std::cout << vec1[i][0] << std::endl;
				std::cout << "is intersect " << (*collider)[obj_index].triangle_AABB[triangle_index].AABB_intersection((*cloth)[i].triangle_AABB[vec1[i][0]]);
			}

			break;
		}
		std::sort(vec1[i].begin(), vec1[i].end());
		std::sort(vec2[i].begin(), vec2[i].end());
		for (int j = 0; j < vec1[i].size(); ++j) {
			if (vec1[i][j] != vec2[i][j]) {
				std::cout<<"error element " << obj_index << " " << triangle_index << " element: " << vec1[i][j] << " " << vec2[i][j] << std::endl;
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
	//for (unsigned int i = 0; i < collider->size(); ++i) {
	//	thread_begin = (*collider)[i].mesh_struct.face_index_begin_per_thread.data();
	//	triangle_aabb = (*collider)[i].triangle_AABB.data();
	//	triangle_neighbor_vertex = &(*collider)[i].triangle_neighbor_obj_vertex;
	//	triangle_neighbor_triangle = &(*collider)[i].triangle_neighbor_obj_triangle;
	//	start = thread_begin[thread_No];
	//	end = thread_begin[thread_No+1];
	//	for (unsigned int k = 0; k < total_obj_num; ++k) {
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
	//		for (unsigned int j = start; j < end; ++j) {
	//			triangle_neighbor_this_obj_vertex = &(*triangle_neighbor_vertex)[j][k];
	//			triangle_neighbor_this_obj_triangle = &(*triangle_neighbor_triangle)[j][k];
	//			triangle_neighbor_this_obj_vertex->clear();				
	//			triangle_neighbor_this_obj_vertex->reserve(triangle_neighbor_this_obj_triangle->size());
	//			for (unsigned int m = 0; m < triangle_neighbor_this_obj_triangle->size(); ++m) {
	//				cloth_triangle_index = (*triangle_neighbor_this_obj_triangle)[m];
	//				for (unsigned int n = 0; n < rep_vertex_num[cloth_triangle_index]; ++n) {
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

	for (unsigned int i = 0; i < total_obj_num; ++i) {
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
		for (unsigned int k = 0; k < total_obj_num; ++k) {
			if (k < tetrahedron_begin_obj_index) {
				compared_triangle_aabb = (*cloth)[k].triangle_AABB.data();
			}
			else {
				compared_triangle_aabb = (*tetrahedron)[k-tetrahedron_begin_obj_index].triangle_AABB.data();
			}
			if (i != k) {
				for (unsigned int j = start; j < end; ++j) {
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
				for (unsigned int j = start; j < end; ++j) {
					vertex_neighbor_this_obj_triangle = &(*vertex_neighbor_triangle)[j][k];
					triangle_neighbor_this_obj_triangle = &(*triangle_neighbor_triangle)[(*vertex_from_rep_triangle_index)[j]][k];//vertex's represent triangle index = (*vertex_from_rep_triangle_index)[j];
					vertex_neighbor_this_obj_triangle->clear();
					vertex_neighbor_this_obj_triangle->reserve(triangle_neighbor_this_obj_triangle->size());
					
					vertex_global_index = (*tetrahedron)[k - tetrahedron_begin_obj_index].mesh_struct.vertex_index_on_sureface[j];

					for (unsigned int m = 0; m < triangle_neighbor_this_obj_triangle->size(); ++m) {
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
	for (unsigned int i = 0; i < total_obj_num; ++i) {
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
		for (unsigned int k = 0; k < collider->size(); ++k) {
			compared_triangle_aabb = (*collider)[k].triangle_AABB.data();
			for (unsigned int j = start; j < end; ++j) {
				vertex_neighbor_this_obj_triangle = &(*vertex_neighbor_triangle)[j][k];
				triangle_neighbor_this_obj_triangle = &(*triangle_neighbor_triangle)[(*vertex_from_rep_triangle_index)[j]][k];//vertex's represent triangle index = (*vertex_from_rep_triangle_index)[j];
				vertex_neighbor_this_obj_triangle->clear();
				vertex_neighbor_this_obj_triangle->reserve(triangle_neighbor_this_obj_triangle->size());
				for (unsigned int m = 0; m < triangle_neighbor_this_obj_triangle->size(); ++m) {
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

	int* face_edges;
	unsigned int* edge_vertex_index;

	for (unsigned int i = 0; i < total_obj_num; ++i) {
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
		for (unsigned int k = i; k < total_obj_num; ++k) {
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
				for (unsigned int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
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
				for (unsigned int j = thread_begin[thread_No]; j < thread_begin[thread_No + 1]; ++j) {
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
//	for (unsigned int i = 0; i < total_triangle_pair_num; i += 4) {
//		triangle0_index = triangle_pair_[i];
//		obj0_index = triangle_pair_[i + 1];
//		triangle1_index = triangle_pair_[i + 2];
//		obj1_index = triangle_pair_[i + 3];
//		if (obj0_index != obj1_index) {
//			for (unsigned int j = 0; j < representative_vertex_num[obj0_index][triangle0_index]; ++j) {
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
//			for (unsigned int j = 0; j < representative_vertex_num[obj0_index][triangle0_index]; ++j) {
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
//	for (unsigned int i = 0; i < total_triangle_pair_num; i += 4) {
//		triangle0_index = triangle_pair_[i];
//		obj0_index = triangle_pair_[i + 1];
//		triangle1_index = triangle_pair_[i + 2];
//		obj1_index = triangle_pair_[i + 3];
//		if (obj0_index != obj1_index) {
//			for (unsigned int j = 0; j < representative_vertex_num[obj1_index][triangle1_index]; ++j) {
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
//			for (unsigned int j = 0; j < representative_vertex_num[obj1_index][triangle1_index]; ++j) {
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
//	for (unsigned int i = 0; i < total_triangle_pair_num; i += 4) {
//		triangle0_index = triangle_pair_[i];
//		obj0_index = triangle_pair_[i + 1];
//		triangle1_index = triangle_pair_[i + 2];
//		obj1_index = triangle_pair_[i + 3];
//		if (obj0_index != obj1_index) {
//			for (unsigned int j = 0; j < representative_edge_num[obj0_index][triangle0_index]; ++j) {
//				for (unsigned int k = 0; k < representative_edge_num[obj1_index][triangle1_index]; ++k) {
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
//			for (unsigned int j = 0; j < representative_edge_num[obj0_index][triangle0_index]; ++j) {
//				for (unsigned int k = 0; k < representative_edge_num[obj1_index][triangle1_index]; ++k) {
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
	//for (unsigned int i = 0; i < total_triangle_pair_num; i += 4) {
	//	triangle0_index = triangle_pair_[i];
	//	obj0_index = triangle_pair_[i + 1];
	//	triangle1_index = triangle_pair_[i + 2];
	//	obj1_index = triangle_pair_[i + 3];
	//	if (obj0_index != obj1_index) {
	//		for (unsigned int j = 0; j < representative_vertex_num[obj0_index][triangle0_index]; ++j) {
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
	//		for (unsigned int j = 0; j < representative_vertex_num[obj1_index][triangle1_index]; ++j) {
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
	//		for (unsigned int j = 0; j < representative_edge_num[obj0_index][triangle0_index]; ++j) {
	//			for (unsigned int k = 0; k < representative_edge_num[obj1_index][triangle1_index]; ++k) {
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
	//		for (unsigned int j = 0; j < representative_vertex_num[obj0_index][triangle0_index]; ++j) {
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
	//		for (unsigned int j = 0; j < representative_vertex_num[obj1_index][triangle1_index]; ++j) {
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
	//		for (unsigned int j = 0; j < representative_edge_num[obj0_index][triangle0_index]; ++j) {
	//			for (unsigned int k = 0; k < representative_edge_num[obj1_index][triangle1_index]; ++k) {
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
	//	for (unsigned int i = 0; i < total_triangle_pair_num_with_collider; i += 4) {
	//		triangle0_index = triangle_pair_collider_[i];
	//		obj0_index = triangle_pair_collider_[i + 1];
	//		triangle1_index = triangle_pair_collider_[i + 2];
	//		obj1_index = triangle_pair_collider_[i + 3];
	//		for (unsigned int j = 0; j < representative_vertex_num[obj0_index][triangle0_index]; ++j) {
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
	//		for (unsigned int j = 0; j < representative_vertex_num_collider[obj1_index][triangle1_index]; ++j) {
	//			if (AABB::AABB_intersection(vertex_aabb_collider[obj1_index][triangle_index_in_order_collider[obj1_index][triangle1_index][j]].data(),
	//				obj_tri_aabb[obj0_index][triangle0_index].data())) {
	//				*(point_collider_triangle_obj_pair_++) = triangle_index_in_order_collider[obj1_index][triangle1_index][j];
	//				*(point_collider_triangle_obj_pair_++) = obj1_index;
	//				*(point_collider_triangle_obj_pair_++) = triangle0_index;
	//				*(point_collider_triangle_obj_pair_++) = obj0_index;
	//			}
	//		}
	//		for (unsigned int j = 0; j < representative_edge_num[obj0_index][triangle0_index]; ++j) {
	//			for (unsigned int k = 0; k < representative_edge_num_collider[obj1_index][triangle1_index]; ++k) {
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