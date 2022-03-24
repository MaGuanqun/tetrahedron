#include"ApproxCCD.h"
#include"../basic/global.h"
#include <iostream>
#include <iomanip>

#undef OBTAIN_CURRENT_POS
#define OBTAIN_CURRENT_POS(dest, v0,v1,e0,e1,t)\
dest[0] = v1[0] - v0[0] + t * (e1[0] - e0[0]);\
dest[1] = v1[1] - v0[1] + t * (e1[1] - e0[1]);\
dest[2] = v1[2] - v0[2] + t * (e1[2] - e0[2]);

bool ApproxCCD::edgeEdgeCCD(double& t, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
	double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
	double* initial_compare_edge_vertex_1, double conservative_rescaling)
{
	//t = 1.0;
	double initial_distance_2 =
		CCD::internal::edgeEdgeDistanceUnclassified(initial_edge_vertex_0,initial_edge_vertex_1,initial_compare_edge_vertex_0,
			initial_compare_edge_vertex_1);
	if (initial_distance_2 == 0.0) {
		t = 0.0;
		return true;
	}	
	double min_distance_2 = (1.0 - conservative_rescaling) * (1.0 - conservative_rescaling) * initial_distance_2;
	double time_ = 1.0;
	//std::cout << "min_distance_2 " << min_distance_2 << std::endl;

	bool is_impacting = edgeEdgeCollisionTime(time_, current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0, initial_edge_vertex_1,
		current_compare_edge_vertex_0, current_compare_edge_vertex_1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1,
		min_distance_2);
	if (is_impacting && time_ < 1e-6) {
		is_impacting = edgeEdgeCollisionTime(time_, current_edge_vertex_0, current_edge_vertex_1, initial_edge_vertex_0, initial_edge_vertex_1,
			current_compare_edge_vertex_0, current_compare_edge_vertex_1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, 0.0);
		if (is_impacting) {
			time_ *= conservative_rescaling;
		}
	}

	if (is_impacting) {
		t = time_;
	}

	//bool is_collide=false;
	//time_ = 1.0;
	//is_collide = vertexEdgeCCD(time_, initial_edge_vertex_0, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, current_edge_vertex_0,
	//	current_compare_edge_vertex_0, current_compare_edge_vertex_1, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = false;
	//is_collide = vertexEdgeCCD(time_, initial_edge_vertex_1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, current_edge_vertex_1,
	//	current_compare_edge_vertex_0, current_compare_edge_vertex_1, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = false;
	//is_collide = vertexEdgeCCD(time_, initial_compare_edge_vertex_0, initial_edge_vertex_0, initial_edge_vertex_1, current_compare_edge_vertex_0,
	//	current_edge_vertex_0, current_edge_vertex_1, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = false;
	//is_collide = vertexEdgeCCD(time_, initial_compare_edge_vertex_1, initial_edge_vertex_0, initial_edge_vertex_1, current_compare_edge_vertex_1,
	//	current_edge_vertex_0, current_edge_vertex_1, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = vertexVertexCCD(time_, initial_compare_edge_vertex_1, current_compare_edge_vertex_1,
	//	initial_edge_vertex_0, current_edge_vertex_0, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = vertexVertexCCD(time_, initial_compare_edge_vertex_1, current_compare_edge_vertex_1,
	//	initial_edge_vertex_1, current_edge_vertex_1, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = vertexVertexCCD(time_, initial_compare_edge_vertex_0, current_compare_edge_vertex_0,
	//	initial_edge_vertex_0, current_edge_vertex_0, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = vertexVertexCCD(time_, initial_compare_edge_vertex_0, current_compare_edge_vertex_0,
	//	initial_edge_vertex_1, current_edge_vertex_1, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//Vec3d initial_edge_0 = Vec3d(initial_edge_vertex_0[0], initial_edge_vertex_0[1], initial_edge_vertex_0[2]);
	//Vec3d current_edge_0 = Vec3d(current_edge_vertex_0[0], current_edge_vertex_0[1], current_edge_vertex_0[2]);
	//Vec3d initial_edge_1 = Vec3d(initial_edge_vertex_1[0], initial_edge_vertex_1[1], initial_edge_vertex_1[2]);
	//Vec3d current_edge_1 = Vec3d(current_edge_vertex_1[0], current_edge_vertex_1[1], current_edge_vertex_1[2]);
	//Vec3d initial_comapre_edge_0 = Vec3d(initial_compare_edge_vertex_0[0], initial_compare_edge_vertex_0[1], initial_compare_edge_vertex_0[2]);
	//Vec3d current_comapre_edge_0 = Vec3d(current_compare_edge_vertex_0[0], current_compare_edge_vertex_0[1], current_compare_edge_vertex_0[2]);
	//Vec3d initial_comapre_edge_1 = Vec3d(initial_compare_edge_vertex_1[0], initial_compare_edge_vertex_1[1], initial_compare_edge_vertex_1[2]);
	//Vec3d current_comapre_edge_1 = Vec3d(current_compare_edge_vertex_1[0], current_compare_edge_vertex_1[1], current_compare_edge_vertex_1[2]);
	//rootparity::RootParityCollisionTest root_parity(initial_edge_0, initial_edge_1, initial_comapre_edge_0, initial_comapre_edge_1,
	//	current_edge_0, current_edge_1, current_comapre_edge_0, current_comapre_edge_1, true);
	//if (root_parity.edge_edge_collision())
	//{
	//	if (!is_impacting) {
	//		std::cout << "EE actually collide but does not test out " << std::endl;
	//		//std::cout << "ee not equation coeff " << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
	//	}
	//	else {
	//		std::cout << "successfully test EE" << std::endl;
	//	}
	//}
	//if (is_impacting) {
		//testEdgeInsideOutside(t, initial_edge_vertex_0, v0, initial_edge_vertex_1, v2, initial_compare_edge_vertex_0, v2, initial_compare_edge_vertex_1, v3);
		//std::cout << "CCD successfully test EE" << std::endl;
	//}
	//if (time_ < 1e-15) {
	//	std::cout << initial_edge_vertex_0[0] << " " << initial_edge_vertex_0[1] << " " << initial_edge_vertex_0[2] << std::endl;
	//	std::cout << initial_edge_vertex_1[0] << " " << initial_edge_vertex_1[1] << " " << initial_edge_vertex_1[2] << std::endl;
	//	std::cout << initial_compare_edge_vertex_0[0] << " " << initial_compare_edge_vertex_0[1] << " " << initial_compare_edge_vertex_0[2] << std::endl;
	//	std::cout << initial_compare_edge_vertex_1[0] << " " << initial_compare_edge_vertex_1[1] << " " << initial_compare_edge_vertex_1[2] << std::endl;		
	//}
	return is_impacting;
}



bool ApproxCCD::edgeEdgeCCD(double& t, Eigen::Vector3d& current_edge_vertex_0, Eigen::Vector3d& current_edge_vertex_1, Eigen::Vector3d& initial_edge_vertex_0, Eigen::Vector3d& initial_edge_vertex_1,
	Eigen::Vector3d& current_compare_edge_vertex_0, Eigen::Vector3d& current_compare_edge_vertex_1, Eigen::Vector3d& initial_compare_edge_vertex_0,
	Eigen::Vector3d& initial_compare_edge_vertex_1, double conservative_rescaling)
{
	//t = 1.0;
	double initial_distance_2 =
		CCD::internal::edgeEdgeDistanceUnclassified(initial_edge_vertex_0.data(), initial_edge_vertex_1.data(), initial_compare_edge_vertex_0.data(),
			initial_compare_edge_vertex_1.data());
	if (initial_distance_2 == 0.0) {
		t = 0.0;
		return true;
	}
	double min_distance_2 = (1.0 - conservative_rescaling) * (1.0 - conservative_rescaling) * initial_distance_2;
	double time_ = 1.0;
	//std::cout << "min_distance_2 " << min_distance_2 << std::endl;

	bool is_impacting = ctcd.edgeEdgeCTCD(initial_edge_vertex_1, initial_edge_vertex_0, initial_compare_edge_vertex_1,
		initial_compare_edge_vertex_0,	current_edge_vertex_1, current_edge_vertex_0, current_compare_edge_vertex_1, current_compare_edge_vertex_0,
		min_distance_2, time_);
	if (is_impacting && time_ < 1e-6) {
		is_impacting = ctcd.edgeEdgeCTCD(initial_edge_vertex_1, initial_edge_vertex_0, initial_compare_edge_vertex_1,
			initial_compare_edge_vertex_0, current_edge_vertex_1, current_edge_vertex_0, current_compare_edge_vertex_1, current_compare_edge_vertex_0,
			0.0, time_);
		if (is_impacting) {
			time_ *= conservative_rescaling;
		}
	}

	if (is_impacting) {
		t = time_;
	}
	return is_impacting;
}



bool ApproxCCD::vertexTriangleCCD(double& t, double* initial_position, double* current_position,
	double* initial_triangle_1, double* current_triangle_1, double* initial_triangle_2, double* current_triangle_2, double* initial_triangle_3, double* current_triangle_3,
	double conservative_rescaling)
{
	double initial_distance_2 = CCD::internal::pointTriangleDistanceUnclassified(initial_position, initial_triangle_1, initial_triangle_2, initial_triangle_3);
	if (initial_distance_2 == 0.0) {
		t = 0.0;
		return true;
	}
	double min_distance_2 = (1.0 - conservative_rescaling)* (1.0 - conservative_rescaling) * initial_distance_2;	
	//t = 1.0;

	double time_ = 1.0;

	bool is_impacting = vertexTriangleCollisionTime(time_, initial_position, current_position, initial_triangle_1, current_triangle_1,
		initial_triangle_2, current_triangle_2, initial_triangle_3, current_triangle_3, min_distance_2);

	if (is_impacting && time_ < 1e-6) {
		is_impacting = vertexTriangleCollisionTime(time_, initial_position, current_position, initial_triangle_1, current_triangle_1,
			initial_triangle_2, current_triangle_2, initial_triangle_3, current_triangle_3, 0.0);
		if (is_impacting) {
			time_ *= conservative_rescaling;
		}
	}

	if (is_impacting) {
		t = time_;
	}	
	//bool is_collide=false;
	//
	//is_collide = vertexEdgeCCD(time_, initial_position, initial_triangle_1, initial_triangle_2, current_position,
	//	current_triangle_1, current_triangle_2, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;

	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//	if (time_ < 1e-15) {
	//		std::cout << "time_t " << std::endl;
	//		std::cout << initial_position[0] << " " << initial_position[1] << " " << initial_position[2] << std::endl;

	//		std::cout << initial_triangle_1[0] << " " << initial_triangle_1[1] << " " << initial_triangle_1[2] << std::endl;
	//		std::cout << initial_triangle_2[0] << " " << initial_triangle_2[1] << " " << initial_triangle_2[2] << std::endl;
	//		std::cout << initial_triangle_3[0] << " " << initial_triangle_3[1] << " " << initial_triangle_3[2] << std::endl;
	//		std::cout << current_position[0] << " " << current_position[1] << " " << current_position[2] << std::endl;
	//		std::cout << current_triangle_1[0] << " " << current_triangle_1[1] << " " << current_triangle_1[2] << std::endl;
	//		std::cout << current_triangle_2[0] << " " << current_triangle_2[1] << " " << current_triangle_2[2] << std::endl;
	//		std::cout << current_triangle_3[0] << " " << current_triangle_3[1] << " " << current_triangle_3[2] << std::endl;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = vertexEdgeCCD(time_, initial_position, initial_triangle_1, initial_triangle_3, current_position,
	//	current_triangle_1, current_triangle_3, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = vertexEdgeCCD(time_, initial_position, initial_triangle_3, initial_triangle_2, current_position,
	//	current_triangle_3, current_triangle_2, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = vertexVertexCCD(time_, initial_position, current_position,
	//	initial_triangle_1, current_triangle_1, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = vertexVertexCCD(time_, initial_position, current_position,
	//	initial_triangle_2, current_triangle_2, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//time_ = 1.0;
	//is_collide = vertexVertexCCD(time_, initial_position, current_position,
	//	initial_triangle_3, current_triangle_3, conservative_rescaling);
	//is_impacting = is_collide || is_impacting;
	//if (is_collide) {
	//	if (time_ < t) {
	//		t = time_;
	//	}
	//}
	//if (time_ < 1e-15) {
	//	std::cout << "t" << std::endl;
	//	std::cout << initial_position[0] << " " << initial_position[1] << " " << initial_position[2] << std::endl;
	//	
	//	std::cout << initial_triangle_1[0] << " " << initial_triangle_1[1] << " " << initial_triangle_1[2] << std::endl;
	//	std::cout << initial_triangle_2[0] << " " << initial_triangle_2[1] << " " << initial_triangle_2[2] << std::endl;
	//	std::cout << initial_triangle_3[0] << " " << initial_triangle_3[1] << " " << initial_triangle_3[2] << std::endl;
	//	std::cout << current_position[0] << " " << current_position[1] << " " << current_position[2] << std::endl;
	//	std::cout << current_triangle_1[0] << " " << current_triangle_1[1] << " " << current_triangle_1[2] << std::endl;
	//	std::cout << current_triangle_2[0] << " " << current_triangle_2[1] << " " << current_triangle_2[2] << std::endl;
	//	std::cout << current_triangle_3[0] << " " << current_triangle_3[1] << " " << current_triangle_3[2] << std::endl;
	//}
	//	//We add an extral collision test for CCD 
	//Vec3d initial_pos_a = Vec3d(initial_position[0], initial_position[1], initial_position[2]);
	//Vec3d current_pos_a = Vec3d(current_position[0], current_position[1], current_position[2]);
	//Vec3d initial_triangle_0_a = Vec3d(initial_triangle_1[0], initial_triangle_1[1], initial_triangle_1[2]);
	//Vec3d initial_triangle_1_a = Vec3d(initial_triangle_2[0], initial_triangle_2[1], initial_triangle_2[2]);
	//Vec3d initial_triangle_2_a = Vec3d(initial_triangle_3[0], initial_triangle_3[1], initial_triangle_3[2]);
	//Vec3d current_triangle_0_a = Vec3d(current_triangle_1[0], current_triangle_1[1], current_triangle_1[2]);
	//Vec3d current_triangle_1_a = Vec3d(current_triangle_2[0], current_triangle_2[1], current_triangle_2[2]);
	//Vec3d current_triangle_2_a = Vec3d(current_triangle_3[0], current_triangle_3[1], current_triangle_3[2]);
	//rootparity::RootParityCollisionTest root_parity(initial_pos_a, initial_triangle_0_a, initial_triangle_1_a, initial_triangle_2_a,
	//	current_pos_a, current_triangle_0_a, current_triangle_1_a, current_triangle_2_a, false);
	//if (root_parity.point_triangle_collision())
	//{
	//	if (!is_impacting) {
	//		std::cout << "PT actually collide but does not test out " << std::endl;
	//		Eigen::Vector3d initial_pos_a_ = Eigen::Vector3d(initial_position[0], initial_position[1], initial_position[2]);
	//		Eigen::Vector3d current_pos_a_ = Eigen::Vector3d(current_position[0], current_position[1], current_position[2]);
	//		Eigen::Vector3d initial_triangle_0_a_ = Eigen::Vector3d(initial_triangle_1[0], initial_triangle_1[1], initial_triangle_1[2]);
	//		Eigen::Vector3d initial_triangle_1_a_ = Eigen::Vector3d(initial_triangle_2[0], initial_triangle_2[1], initial_triangle_2[2]);
	//		Eigen::Vector3d initial_triangle_2_a_ = Eigen::Vector3d(initial_triangle_3[0], initial_triangle_3[1], initial_triangle_3[2]);
	//		Eigen::Vector3d current_triangle_0_a_ = Eigen::Vector3d(current_triangle_1[0], current_triangle_1[1], current_triangle_1[2]);
	//		Eigen::Vector3d current_triangle_1_a_ = Eigen::Vector3d(current_triangle_2[0], current_triangle_2[1], current_triangle_2[2]);
	//		Eigen::Vector3d current_triangle_2_a_ = Eigen::Vector3d(current_triangle_3[0], current_triangle_3[1], current_triangle_3[2]);
	//		//std::cout << "====" << std::endl;
	//		std::cout << initial_pos_a_ << std::endl;
	//		std::cout << initial_triangle_0_a_ << std::endl;
	//		std::cout << initial_triangle_1_a_ << std::endl;
	//		std::cout << initial_triangle_2_a_ << std::endl;
	//		//double initial_distance_2 =
	//		//	CCD::internal::pointTriangleDistanceUnclassified(initial_position, initial_triangle_1, initial_triangle_2,
	//		//		initial_triangle_3);
	//		//CTCD ctcd;
	//		//double tt;
	//		//if (ctcd.vertexFaceCTCD(initial_pos_a_, initial_triangle_0_a_, initial_triangle_1_a_, initial_triangle_2_a_,
	//		//	current_pos_a_, current_triangle_0_a_, current_triangle_1_a_, current_triangle_2_a_, tolerance_2, tt)) {
	//		//	std::cout << tolerance_2 << " ground truth ctcd collide " << tt << std::endl;
	//		//}
	//		//else {
	//		//	std::cout << tolerance_2 << " ground truth ctcd not collide " << std::endl;
	//		//}
	//	}
	//	else {
	//		std::cout << "successfully test PT" << std::endl;
	//	}
	//}
	//if (is_collide) {
	//	std::cout << "CCD successfully test PT" << std::endl;
	//	std::cout << "true " << std::endl;
	//	if (!testInsideOutside(t, initial_position, v0, initial_triangle_1, v1, initial_triangle_2, v2, initial_triangle_3, v3)) {
	//		std::cout << mint << " " << tolerance_2 << std::endl;
	//		
	//	
	//	}
	//}

	
	return is_impacting;
}

bool ApproxCCD::vertexTriangleCCD(double& t, Eigen::Vector3d& initial_position, Eigen::Vector3d& current_position,
	Eigen::Vector3d& initial_triangle_1, Eigen::Vector3d& current_triangle_1, Eigen::Vector3d& initial_triangle_2, Eigen::Vector3d& current_triangle_2, Eigen::Vector3d& initial_triangle_3, Eigen::Vector3d& current_triangle_3,
	double conservative_rescaling)
{
	double initial_distance_2 = CCD::internal::pointTriangleDistanceUnclassified(initial_position.data(), initial_triangle_1.data(), 
		initial_triangle_2.data(), initial_triangle_3.data());
	if (initial_distance_2 == 0.0) {
		t = 0.0;
		return true;
	}
	double min_distance_2 = (1.0 - conservative_rescaling) * (1.0 - conservative_rescaling) * initial_distance_2;	//t = 1.0;

	double time_ = 1.0;

	bool is_impacting = ctcd.vertexFaceCTCD(initial_position, initial_triangle_1, initial_triangle_2, initial_triangle_3,
		current_position,  current_triangle_1, current_triangle_2, current_triangle_3, min_distance_2, time_);

	if (is_impacting && time_ < 1e-6) {
		is_impacting = ctcd.vertexFaceCTCD(initial_position, initial_triangle_1, initial_triangle_2, initial_triangle_3,
			current_position, current_triangle_1, current_triangle_2, current_triangle_3, 0.0, time_);
		if (is_impacting) {
			time_ *= conservative_rescaling;
		}
	}

	if (is_impacting) {
		t = time_;
	}
	return is_impacting;
}



bool ApproxCCD::vertexEdgeCCD(double& t, double* q0_initial, double* q1_initial, double* q2_initial,
	double* q0_current, double* q1_current, double* q2_current, double conservative_rescaling)
{
	double initial_distance_2 = CCD::internal::pointEdgeDistanceUnclassified(q0_initial, q1_initial, q2_initial);
	if (initial_distance_2 == 0.0) {
		t = 0.0;
		return true;
	}

	double time_ = 1.0;
	double min_distance_2 = (1.0 - conservative_rescaling) * (1.0 - conservative_rescaling) * initial_distance_2;
	bool is_impacting = vertexEdgeCollisionTime(time_, q0_initial, q1_initial, q2_initial, q0_current,
		q1_current, q2_current, min_distance_2);

	if (is_impacting && time_ < 1e-6) {
		is_impacting = vertexEdgeCollisionTime(time_, q0_initial, q1_initial, q2_initial, q0_current,
			q1_current, q2_current, 0.0);
		if (is_impacting) {
			time_ *= conservative_rescaling;
		}
	}

	if (is_impacting) {
		t = time_;
	}

	return is_impacting;
}


bool ApproxCCD::vertexEdgeCCD(double& t, Eigen::Vector3d& q0_initial, Eigen::Vector3d& q1_initial, Eigen::Vector3d& q2_initial,
	Eigen::Vector3d& q0_current, Eigen::Vector3d& q1_current, Eigen::Vector3d& q2_current, double conservative_rescaling)
{
	double initial_distance_2 = CCD::internal::pointEdgeDistanceUnclassified(q0_initial.data(), q1_initial.data(), q2_initial.data());
	if (initial_distance_2 == 0.0) {
		t = 0.0;
		return true;
	}

	double time_ = 1.0;
	double min_distance_2 = (1.0 - conservative_rescaling) * (1.0 - conservative_rescaling) * initial_distance_2;
	bool is_impacting = ctcd.vertexEdgeCTCD(q0_initial, q1_initial, q2_initial, q0_current,
		q1_current, q2_current, min_distance_2, time_);

	if (is_impacting && time_ < 1e-6) {
		is_impacting = ctcd.vertexEdgeCTCD(q0_initial, q1_initial, q2_initial, q0_current,
			q1_current, q2_current, 0.0, time_);
		if (is_impacting) {
			time_ *= conservative_rescaling;
		}
	}

	if (is_impacting) {
		t = time_;
	}

	return is_impacting;
}


bool ApproxCCD::vertexVertexCCD(double& t, double* initial_position, double* current_position,
	double* initial_position_1, double* current_position_1, double conservative_rescaling)
{
	double initial_distance_2 = CCD::internal::pointPointDistance(initial_position, initial_position_1);
	if (initial_distance_2 == 0.0) {
		t = 0.0;
		return true;
	}
	double time_ = 1.0;
	double min_distance_2 = (1.0 - conservative_rescaling) * (1.0 - conservative_rescaling) * initial_distance_2;
	bool is_impacting = vertexVertexCollisionTime(time_, initial_position, current_position, initial_position_1, current_position_1, min_distance_2);

	if (is_impacting && time_ < 1e-6) {
		is_impacting = vertexVertexCollisionTime(time_, initial_position, current_position, initial_position_1, current_position_1, 0.0);
		if (is_impacting) {
			time_ *= conservative_rescaling;
		}
	}

	if (is_impacting) {
		t = time_;
	}
	return is_impacting;
}


bool ApproxCCD::vertexVertexCCD(double& t, Eigen::Vector3d& initial_position, Eigen::Vector3d& current_position,
	Eigen::Vector3d& initial_position_1, Eigen::Vector3d& current_position_1, double conservative_rescaling)
{
	double initial_distance_2 = CCD::internal::pointPointDistance(initial_position.data(), initial_position_1.data());
	if (initial_distance_2 == 0.0) {
		t = 0.0;
		return true;
	}
	double time_ = 1.0;
	double min_distance_2 = (1.0 - conservative_rescaling) * (1.0 - conservative_rescaling) * initial_distance_2;
	bool is_impacting = ctcd.vertexVertexCTCD(initial_position, initial_position_1,  current_position, current_position_1, min_distance_2,time_);

	if (is_impacting && time_ < 1e-6) {
		is_impacting = ctcd.vertexVertexCTCD(initial_position, initial_position_1, current_position, current_position_1, 0.0, time_);
		if (is_impacting) {
			time_ *= conservative_rescaling;
		}
	}

	if (is_impacting) {
		t = time_;
	}
	return is_impacting;
}


bool ApproxCCD::edgeEdgeCollisionTime(double& t, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
	double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
	double* initial_compare_edge_vertex_1, double tolerance_2)
{
	//std::cout << initial_edge_vertex_0[0] << " " << initial_edge_vertex_0[1] << " " << initial_edge_vertex_0[2] << std::endl;
	//std::cout << initial_edge_vertex_1[0] << " " << initial_edge_vertex_1[1] << " " << initial_edge_vertex_1[2] << std::endl;
	//std::cout << initial_compare_edge_vertex_0[0] << " " << initial_compare_edge_vertex_0[1] << " " << initial_compare_edge_vertex_0[2] << std::endl;
	//std::cout << initial_compare_edge_vertex_1[0] << " " << initial_compare_edge_vertex_1[1] << " " << initial_compare_edge_vertex_1[2] << std::endl;


	std::vector<TimeInterval> rawcoplane, a0, a1, b0, b1;
	double v0[3], v1[3], v2[3], v3[3];
	SUB(v0, current_edge_vertex_0, initial_edge_vertex_0);
	SUB(v1, current_edge_vertex_1, initial_edge_vertex_1);
	SUB(v2, current_compare_edge_vertex_0, initial_compare_edge_vertex_0);
	SUB(v3, current_compare_edge_vertex_1, initial_compare_edge_vertex_1);

	double x10[3], x20[3], x30[3];
	double v10[3], v20[3], v30[3];

	SUB(x10, initial_edge_vertex_0, initial_compare_edge_vertex_0);
	SUB(v10, v0, v2);
	SUB(x20, initial_edge_vertex_0, initial_edge_vertex_1);
	SUB(v20, v0, v1);
	SUB(x30, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
	SUB(v30, v2, v3);

	distancePoly3D(x10, x20, x30, v10, v20, v30, tolerance_2, rawcoplane);

	// check for parallel edges
	std::vector<TimeInterval> coplane;
	//coplane.reserve(rawcoplane.size());
	std::vector<TimeInterval> parallel;
	//parallel.reserve(rawcoplane.size());
	double temp[3];
	for (unsigned int i = 0; i < rawcoplane.size(); i++)
	{
		double midt = (rawcoplane[i].u + rawcoplane[i].l) / 2;
		OBTAIN_CURRENT_POS(x10, initial_edge_vertex_0, initial_edge_vertex_1, v0, v1, midt);
		OBTAIN_CURRENT_POS(x20, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1, v2, v3, midt);
		CROSS(temp, x10, x20);
		if (DOT(temp,temp) < 1e-16)
		{
			parallel.emplace_back(rawcoplane[i]);
		}
		else
		{
			coplane.emplace_back(rawcoplane[i]);
		}
	}
	if (coplane.empty())
		return false;

	SUB(x10, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
	SUB(v10, v2, v3);
	SUB(x20, initial_edge_vertex_0, initial_edge_vertex_1);
	SUB(v20, v0, v1);
	SUB(x30, initial_edge_vertex_1, initial_compare_edge_vertex_1);
	SUB(v30, v1, v3);

	barycentricPoly3D(x10, x20, x30, v10, v20, v30, a0);
	if (a0.empty())
		return false;

	SUB(x20, initial_edge_vertex_1, initial_edge_vertex_0);
	SUB(v20, v1, v0);
	SUB(x30, initial_edge_vertex_0, initial_compare_edge_vertex_1);
	SUB(v30, v0, v3);
	barycentricPoly3D(x10, x20, x30, v10, v20, v30, a1);
	if (a1.empty())
		return false;

	SUB(x10, initial_edge_vertex_0, initial_edge_vertex_1);
	SUB(v10, v0, v1);
	SUB(x20, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1);
	SUB(v20, v2, v3);
	SUB(x30, initial_compare_edge_vertex_1, initial_edge_vertex_1);
	SUB(v30, v3, v1);
	barycentricPoly3D(x10, x20, x30, v10, v20, v30, b0);
	if (b0.empty())
		return false;

	SUB(x20, initial_compare_edge_vertex_1, initial_compare_edge_vertex_0);
	SUB(v20, v3, v2);
	SUB(x30, initial_compare_edge_vertex_0, initial_edge_vertex_1);
	SUB(v30, v2, v1);
	barycentricPoly3D(x10, x20, x30, v10, v20, v30, b1);
	if (b1.empty())
		return false;

	// check intervals for overlap
	bool col = false;
	double mint = 1.0;
	std::vector<TimeInterval> intervals(5);
	std::vector<TimeInterval> pcheck(2);
	TimeInterval isect;
	for (unsigned int i = 0; i < coplane.size(); i++){
		for (unsigned int j = 0; j < a0.size(); j++){
			for (unsigned int k = 0; k < a1.size(); k++){
				for (unsigned int l = 0; l < b0.size(); l++){
					for (unsigned int m = 0; m < b1.size(); m++){						
						intervals[0]=coplane[i];
						intervals[1]=a0[j];
						intervals[2]=a1[k];
						intervals[3]=b0[l];
						intervals[4]=b1[m];
						if (TimeInterval::overlap(intervals)){
							isect = TimeInterval::intersect(intervals);
							bool skip = false;
							for (unsigned int p = 0; p < parallel.size(); p++){								
								pcheck[0]=isect;
								pcheck[1]=parallel[p];
								if (TimeInterval::overlap(pcheck)){
									skip = true;
									break;
								}
							}
							if (!skip){
								mint = (std::min)(mint, isect.l);
								col = true;
							}
						}
					}
				}
			}
		}
	}


	//std::cout << "co--------" << std::endl;
	//for (unsigned int i = 0; i < coplane.size(); i++) {
	//	std::cout << coplane[i].l << " " << coplane[i].u << std::endl;
	//}
	//std::cout << "a0--------" << std::endl;
	//for (unsigned int i = 0; i < a0.size(); i++) {
	//	std::cout << a0[i].l << " " << a0[i].u << std::endl;
	//}
	//std::cout << "a1--------" << std::endl;
	//for (unsigned int i = 0; i < a1.size(); i++) {
	//	std::cout << a1[i].l << " " << a1[i].u << std::endl;
	//}
	//std::cout << "b0-------" << std::endl;
	//for (unsigned int i = 0; i < b0.size(); i++) {
	//	std::cout << b0[i].l << " " << b0[i].u << std::endl;
	//}
	//std::cout << "b1--------" << std::endl;
	//for (unsigned int i = 0; i < b1.size(); i++) {
	//	std::cout << b1[i].l << " " << b1[i].u << std::endl;
	//}

	bool is_collide = false;

	if (col){
		t = mint;
		is_collide = true;
		//return true;
	}
	
	return is_collide;

}


bool ApproxCCD::vertexEdgeCollisionTime(double& t, double* q0_initial, double* q1_initial, double* q2_initial,
	double* q0_current, double* q1_current, double* q2_current, double tolerance_2)
{
	double op[5];
	double v0[3], v1[3], v2[3];
	SUB(v0, q0_current, q0_initial);
	SUB(v1, q1_current, q1_initial);
	SUB(v2, q2_current, q2_initial);

	std::vector<TimeInterval> colin, e1, e2;
	double ab[3], ac[3], cb[3];
	double vab[3], vac[3], vcb[3];

	SUB(ab, q2_initial, q1_initial);
	SUB(ac, q0_initial, q1_initial);
	SUB(cb, q2_initial, q0_initial);
	SUB(vab, v2, v1);
	SUB(vac, v0, v1);
	SUB(vcb, v2, v0);

	double c = DOT(ab, ac);
	double b = DOT(ac, vab) + DOT(ab, vac);
	double a = DOT(vab, vac);
	op[0] = a;
	op[1] = b;
	op[2] = c;
	findIntervals(op, 2, e1, true);
	if (e1.empty())
		return false;
	c = DOT(ab, cb);
	b = DOT(cb,vab) + DOT(ab,vcb);
	a = DOT(vab,vcb);
	op[0] = a;
	op[1] = b;
	op[2] = c;
	findIntervals(op, 2, e2, true);
	if (e2.empty())
		return false;
	double A = DOT(ab,ab);
	double B = 2.0 * DOT(ab,vab);
	double C = DOT(vab, vab);
	double D = DOT(ac, ac);
	double E = 2.0 * DOT(ac,vac);
	double F = DOT(vac,vac);
	double G = DOT(ac, ab);
	double H = DOT(vab, ac) + DOT(vac, ab);
	double I = DOT(vab, vac);
	op[4] = A * D - G * G - tolerance_2 * A;
	op[3] = B * D + A * E - 2 * G * H - tolerance_2 * B;
	op[2] = B * E + A * F + C * D - H * H - 2 * G * I - tolerance_2 * C;
	op[1] = B * F + C * E - 2 * H * I;
	op[0] = C * F - I * I;
	findIntervals(op, 4, colin, false);
	if (colin.empty())
		return false;

	double mint = 1.0;
	bool col = false;
	std::vector<TimeInterval> intervals(3);
	for (unsigned int i = 0; i < colin.size(); i++){
		for (unsigned int j = 0; j < e1.size(); j++){
			for (unsigned int k = 0; k < e2.size(); k++){				
				intervals[0]=colin[i];
				intervals[1]=e1[j];
				intervals[2]=e2[k];
				if (TimeInterval::overlap(intervals)){
					mint = (std::min)(TimeInterval::intersect(intervals).l, mint);
					col = true;
				}
			}
		}
	}
	if (col)
	{
		t = mint;
		return true;
	}
	return false;
}

bool ApproxCCD::vertexVertexCollisionTime(double& t, double* initial_position, double* current_position,
	double* initial_position_1, double* current_position_1, double tolerance_2)
{
	int roots = 0;
	double t1 = 0, t2 = 0;
	double v1[3], v2[3];
	SUB(v1, current_position, initial_position);
	SUB(v2, current_position_1, initial_position_1);
	
	// t^2 term
	double a = DOT(v1,v1) + DOT(v2,v2) - 2.0 * DOT(v1,v2);
	// t term
	double b = 2.0 * (DOT(v1, initial_position) - DOT(v2, initial_position) - DOT(v1, initial_position_1) + DOT(v2, initial_position_1));
	// current distance - min_d
	double c = DOT(initial_position, initial_position) + DOT(initial_position_1, initial_position_1) - 2.0 * DOT(initial_position, initial_position_1) - tolerance_2;
	if (a != 0){
		roots = getQuadRoots(a, b, c, t1, t2);
	}
	else if (b != 0){
		t1 = -c / b;
		roots = 1;
	}
	else{
		if (c <= 0)	{
			t = 0;
			return true;
		}
		return false;
	}

	double op[3];
	op[0] = a;
	op[1] = b;
	op[2] = c;
	std::vector<TimeInterval> interval;
	if (roots == 2){
		checkInterval(0, t1, op, 2, interval, false);
		if (!interval.empty())
		{
			t = 0;
			return true;
		}
		checkInterval(t1, t2, op, 2, interval, false);
		if (!interval.empty())
		{
			t = t1;
			return true;
		}
		checkInterval(t2, 1.0, op, 2, interval, false);
		if (!interval.empty())
		{
			t = t2;
			return true;
		}
		return false;
	}
	else if (roots == 1){
		checkInterval(0, t1, op, 2, interval, false);
		if (!interval.empty())
		{
			t = 0;
			return true;
		}
		checkInterval(t1, 1.0, op, 2, interval, false);
		if (!interval.empty())
		{
			t = t1;
			return true;
		}
		return false;
	}
	checkInterval(0, 1.0, op, 2, interval, false);
	if (!interval.empty()){
		t = 0;
		return true;
	}
	return false;
}

bool ApproxCCD::vertexTriangleCollisionTime(double& t, double* initial_position, double* current_position,
	double* initial_triangle_1, double* current_triangle_1, double* initial_triangle_2, double* current_triangle_2, double* initial_triangle_3, double* current_triangle_3,
	double tolerance_2)
{
	double v0[3], v1[3], v2[3], v3[3];
	SUB(v0, current_position, initial_position);
	SUB(v1, current_triangle_1, initial_triangle_1);
	SUB(v2, current_triangle_2, initial_triangle_2);
	SUB(v3, current_triangle_3, initial_triangle_3);

	double x10[3], x20[3], x30[3];
	double v10[3], v20[3], v30[3];

	std::vector<TimeInterval> coplane, e1, e2, e3;

	double temp[3];
	//here, the plane is constructed by v1v3 and normal,
	//if there is a collision, the collision point must on the positve (front) side of the plane,
	SUB(x10, initial_position, initial_triangle_1);
	SUB(v10, v0, v1);
	SUB(x30, initial_triangle_3, initial_triangle_1);
	SUB(v30, v3, v1);
	SUB(temp, initial_triangle_2, initial_triangle_1);
	CROSS(x20, x30, temp);
	SUB(temp, v2, v1);
	CROSS(v20, v30, temp);
	
	planePoly3D(x10, x20, x30, v10, v20, v30, e1);

	if (e1.empty()) {
		return false;
	}

	SUB(x10, initial_position, initial_triangle_2);
	SUB(v10, v0, v2);
	SUB(x30, initial_triangle_1, initial_triangle_2);
	SUB(v30, v1, v2);
	SUB(temp, initial_triangle_3, initial_triangle_2);
	CROSS(x20, x30, temp);
	SUB(temp, v3, v2);
	CROSS(v20, v30, temp);
	planePoly3D(x10, x20, x30, v10, v20, v30, e2);

	if (e2.empty())
		return false;

	SUB(x10, initial_position, initial_triangle_3);
	SUB(v10, v0, v3);
	SUB(x30, initial_triangle_2, initial_triangle_3);
	SUB(v30, v2, v3);
	SUB(temp, initial_triangle_1, initial_triangle_3);
	CROSS(x20, x30, temp);
	SUB(temp, v1, v3);
	CROSS(v20, v30, temp);

	planePoly3D(x10, x20, x30, v10, v20, v30, e3);

	if (e3.empty())
		return false;

	SUB(x10, initial_position, initial_triangle_1);
	SUB(v10, v0, v1);
	SUB(x30, initial_triangle_3, initial_triangle_1);
	SUB(v30, v3, v1);
	SUB(x20, initial_triangle_2, initial_triangle_1);
	SUB(v20, v2, v1);

	distancePoly3D(x10, x20, x30, v10, v20, v30, tolerance_2, coplane);

	if (coplane.empty())
		return false;

	bool col = false;
	double mint = 1.0;


	std::vector<TimeInterval> intervals(4);
	for (unsigned int i = 0; i < coplane.size(); i++){
		for (unsigned int j = 0; j < e1.size(); j++){
			for (unsigned int k = 0; k < e2.size(); k++){
				for (unsigned int l = 0; l < e3.size(); l++){
					intervals[0]=coplane[i];
					intervals[1]=e1[j];
					intervals[2]=e2[k];
					intervals[3]=e3[l];
					if (TimeInterval::overlap(intervals))
					{
						mint = (std::min)(TimeInterval::intersect(intervals).l, mint);
						col = true;
					}
				}
			}
		}
	}

	bool is_collide = false;
	if (col){
		t = mint;
		is_collide = true;
	}


	return is_collide;

}






bool ApproxCCD::TimeInterval::overlap(const TimeInterval& t1, const TimeInterval& t2)
{
	return !(t1.l > t2.u || t2.l > t1.u);
}


bool ApproxCCD::TimeInterval::overlap(const std::vector<TimeInterval>& intervals)
{
	for (std::vector<TimeInterval>::const_iterator it1 = intervals.begin(); it1 != intervals.end(); ++it1)
	{		std::vector<TimeInterval>::const_iterator it2 = it1;
		for (++it2; it2 != intervals.end(); ++it2)
			if (it1->l > it2->u || it2->l > it1->u)
				return false;
	}
	return true;
}

ApproxCCD::TimeInterval ApproxCCD::TimeInterval::intersect(const std::vector<TimeInterval>& intervals)
{
	TimeInterval isect(0.0, 1.0);
	for (std::vector<TimeInterval>::const_iterator it = intervals.begin(); it != intervals.end(); ++it)
	{
		if (isect.l < it->l) {
			isect.l = it->l;
		}
		if (isect.u > it->u) {
			isect.u = it->u;
		}
	}
	return isect;
}

void ApproxCCD::distancePoly3D(double* x10, double* x20, double* x30, double* v10, double* v20, double* v30, double minDSquared,
	std::vector<TimeInterval>& result)
{
	double temp_v20_v30[3], temp_v20_x30_x20_v30[3], temp_x20_v30[3];
	double temp_x20_x30[3];
	CROSS(temp_v20_v30, v20, v30);
	CROSS(temp_v20_x30_x20_v30, v20, x30);
	CROSS(temp_x20_v30, x20, v30);
	SUM_(temp_v20_x30_x20_v30, temp_x20_v30);
	double A = DOT(v10, temp_v20_v30);
	double B = DOT(x10, temp_v20_v30) + DOT(v10, temp_v20_x30_x20_v30);
	CROSS(temp_x20_x30, x20, x30);
	double C = DOT(x10, temp_v20_x30_x20_v30) + DOT(v10, temp_x20_x30);
	double D = DOT(x10, temp_x20_x30);
	double op[7];
	op[0] = A * A;
	op[1] = 2.0 * A * B;
	op[2] = B * B + 2.0 * A * C - DOT(temp_v20_v30, temp_v20_v30) * minDSquared;
	op[3] = 2.0 * A * D + 2.0 * B * C - 2.0 * DOT(temp_v20_v30, temp_v20_x30_x20_v30) * minDSquared;
	op[4] = 2.0 * B * D + C * C - (2.0 * DOT(temp_v20_v30, temp_x20_x30) + DOT(temp_v20_x30_x20_v30, temp_v20_x30_x20_v30)) * minDSquared;
	op[5] = 2.0 * C * D - 2.0 * DOT(temp_v20_x30_x20_v30, temp_x20_x30) * minDSquared;
	op[6] = D * D - DOT(temp_x20_x30, temp_x20_x30) * minDSquared;
	//result.reserve(6);
	findIntervals(op, 6, result, false);
}




void ApproxCCD::planePoly3D(double* x10, double* x20, double* x30, double* v10, double* v20, double* v30,
	std::vector<TimeInterval>& result)
{
	double op[4];
	double temp_2030[3], temp_v20_x30[3], temp_x20_v30[3];
	CROSS(temp_2030, v20, v30);
	CROSS(temp_v20_x30, v20, x30);
	CROSS(temp_x20_v30, x20, v30);
	op[0] = DOT(v10, temp_2030);
	op[1] = DOT(x10, temp_2030) + DOT(v10, temp_x20_v30) + DOT(v10,temp_v20_x30);
	CROSS(temp_2030, x20, x30);
	op[2] = DOT(x10, temp_x20_v30) + DOT(x10, temp_v20_x30) + DOT(v10, temp_2030);
	op[3] = DOT(x10, temp_2030);

	//result.reserve(3);
	findIntervals(op, 3, result, true);
}




void ApproxCCD::barycentricPoly3D(double* x10, double* x20, double* x30, double* v10, double* v20, double* v30,
	std::vector<TimeInterval>& result)
{
	double A = DOT(x10, x10);
	double B = 2.0 * DOT(x10, v10);
	double C = DOT(v10, v10);

	double D = DOT(x20, x10);
	double E = DOT(x20, v10) + DOT(v20, x10);
	double F = DOT(v20, v10);

	double G = DOT(x30, x20);
	double H = DOT(x30, v20) + DOT(v30, x20);
	double I = DOT(v30, v20);

	double J = DOT(x30, x10);
	double K = DOT(x30, v10) + DOT(v30, x10);
	double L = DOT(v30, v10);

	double op[5];
	op[0] = F * L - C * I;
	op[1] = F * K + E * L - C * H - B * I;
	op[2] = F * J + D * L + E * K - C * G - A * I - B * H;
	op[3] = D * K + E * J - A * H - B * G;
	op[4] = D * J - A * G;
	//result.reserve(4);
	findIntervals(op, 4, result, true);
}


//bool ApproxCCD::pointTriangleCollisionTime(double& t, double* initial_position, double* current_position,
//	double* initial_triangle_0, double* current_triangle_0, double* initial_triangle_1, double* current_triangle_1, double* initial_triangle_2, double* current_triangle_2,
//	double* initial_normal_not_normalized, double* current_normal_not_normalized, double* cross_for_CCD, double tolerance_2)//floating* f_initial_normal, floating* f_current_normal, floating* f_cross_for_CCD, 
//{
//	double e_0[3], e_1[3], e_2[3];
//	SUB(e_0, initial_position, initial_triangle_0);
//	SUB(e_1, initial_position, initial_triangle_1);
//	SUB(e_2, initial_position, initial_triangle_2);
//
//	double d = DOT(e_0, initial_normal_not_normalized);
//
//	double u[3], u_0[3], u_1[3], u_2[3];
//	SUB(u, current_position, initial_position);
//	SUB(u_0, current_triangle_0, initial_triangle_0);
//	SUB(u_1, current_triangle_1, initial_triangle_1);
//	SUB(u_2, current_triangle_2, initial_triangle_2);
//
//	double qp_1[3];
//	SUB(qp_1, current_position, current_triangle_0);
//	double b1 = DOT(qp_1, current_normal_not_normalized);
//	double b2 = DOT(qp_1, initial_normal_not_normalized);
//	double b3 = DOT(e_0, cross_for_CCD);
//	double b4 = DOT(qp_1, cross_for_CCD);
//	double b5 = DOT(e_0, current_normal_not_normalized);
//
//	double a3 = -d + b1 + b2 + b3 - b4 - b5;
//	double a2 = 3 * d + b5 - 2 * b3 + b4 - 2 * b2;
//	double a1 = -3 * d + b2 + b3;
//
//	std::vector<double>time;
//	time.reserve(7);
//	double t0 = 2.0; double t1 = 2.0; double t2 = 2.0;
//	bool is_collide = false;
//
//	if(a3 > -NEAR_ZERO && a3<NEAR_ZERO
//		&& a2 > -NEAR_ZERO && a2<NEAR_ZERO
//		&& a1 > -NEAR_ZERO && a1 < NEAR_ZERO){
//		if (d<NEAR_ZERO && d>-NEAR_ZERO) {
//			double e_[3];
//			SUB(e_, initial_triangle_1, initial_triangle_0);
//			if (pointEdgeCollisionTime(t, u, u_0, u_1, e_, e_0, e_1, tolerance_2)) {
//				time.push_back(t);
//			}
//			SUB(e_, initial_triangle_2, initial_triangle_1);
//			if (pointEdgeCollisionTime(t, u, u_1, u_2, e_, e_1, e_2, tolerance_2)) {
//				time.push_back(t);
//			}
//			SUB(e_, initial_triangle_2, initial_triangle_0);
//			if (pointEdgeCollisionTime(t, u, u_0, u_2, e_, e_0, e_2, tolerance_2)) {
//				time.push_back(t);
//			}
//			if (pointPointCollisionTime(t, e_0, u_0, u, tolerance_2)) {
//				time.push_back(t);
//			}
//			if (pointPointCollisionTime(t, e_1, u_1, u, tolerance_2)) {
//				time.push_back(t);
//			}
//			if (pointPointCollisionTime(t, e_2, u_2, u, tolerance_2)) {
//				time.push_back(t);
//			}
//			//double e_0_1[3], e_1_2[3], e_2_0[3];
//			//CROSS(e_0_1, e_0, e_1);
//			//CROSS(e_1_2, e_1, e_2);
//			//CROSS(e_2_0, e_2, e_0);
//			//if (DOT(e_0_1, e_1_2) > 0 && DOT(e_1_2, e_2_0) > 0 && DOT(e_0_1, e_2_0) > 0) {
//			//	t = 0;
//			//	return true;
//			//}		
//		}
//		else {
//			is_collide = false;
//		}
//	}
//	else {
//		if (solveEquation(t0, t1, t2, a3, a2, a1, d)) {
//			if (checkInside(t0, initial_position, initial_triangle_0, initial_triangle_1, initial_triangle_2,
//				u, u_0, u_1, u_2)) {
//				time.push_back(t0);
//				/*if (t < 1e-11) {
//					std::cout << "point-triangle-inside "<<t << std::endl;
//				}*/
//			}
//			else {
//				if (t1 < 1.5) {
//					if (checkInside(t1, initial_position, initial_triangle_0, initial_triangle_1, initial_triangle_2,
//						u, u_0, u_1, u_2)) {
//						time.push_back(t1);
//					}
//					else {
//						if (t2 < 1.5) {
//							if (checkInside(t2, initial_position, initial_triangle_0, initial_triangle_1, initial_triangle_2,
//								u, u_0, u_1, u_2)) {
//								time.push_back(t2);
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//
//
//
//	if (!time.empty()) {
//		t = time[0];
//		for (int i = 1; i < time.size(); ++i) {
//			if (t > time[i]) {
//				t = time[i];
//			}
//		}
//		t *= 0.8;
//
//		is_collide = true;
//	}
//	
//
//	//We add an extral collision test for CCD 
//	Vec3d initial_pos_a = Vec3d(initial_position[0], initial_position[1], initial_position[2]);
//	Vec3d current_pos_a = Vec3d(current_position[0], current_position[1], current_position[2]);
//	Vec3d initial_triangle_0_a = Vec3d(initial_triangle_0[0], initial_triangle_0[1], initial_triangle_0[2]);
//	Vec3d initial_triangle_1_a = Vec3d(initial_triangle_1[0], initial_triangle_1[1], initial_triangle_1[2]);
//	Vec3d initial_triangle_2_a = Vec3d(initial_triangle_2[0], initial_triangle_2[1], initial_triangle_2[2]);
//	Vec3d current_triangle_0_a = Vec3d(current_triangle_0[0], current_triangle_0[1], current_triangle_0[2]);
//	Vec3d current_triangle_1_a = Vec3d(current_triangle_1[0], current_triangle_1[1], current_triangle_1[2]);
//	Vec3d current_triangle_2_a = Vec3d(current_triangle_2[0], current_triangle_2[1], current_triangle_2[2]);
//	rootparity::RootParityCollisionTest root_parity(initial_pos_a, initial_triangle_0_a, initial_triangle_1_a, initial_triangle_2_a,
//				current_pos_a, current_triangle_0_a, current_triangle_1_a, current_triangle_2_a,false);
//	if (root_parity.point_triangle_collision())
//	{
//		if (!is_collide) {
//			std::cout << "PT actually collide but does not test out " << std::endl;
//			std::cout << "vt not equation coeff " << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
//		}
//	}	
//	testInsideOutside(t, initial_position, u, initial_triangle_0, u_0, initial_triangle_1, u_1, initial_triangle_2, u_2);
//
//
//	if (is_collide) {
//		std::cout << "vt equation coeff " << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
//	}
//
//	return is_collide;
//}



//bool ApproxCCD::edgeEdgeCollisionTime(double& t, double* current_edge_vertex_0, double* current_edge_vertex_1, double* initial_edge_vertex_0, double* initial_edge_vertex_1,
//	double* current_compare_edge_vertex_0, double* current_compare_edge_vertex_1, double* initial_compare_edge_vertex_0,
//	double* initial_compare_edge_vertex_1, double tolerance_2)
//{
//	double e0[3], e0_[3], e1[3], e1_[3], e2[3], e2_[3];
//	//double n0[3];
//	//floating f_e0[3], f_e0_[3], f_e1[3], f_e1_[3], f_e2[3], f_e2_[3];
//	//floating f_current_edge_vertex_0[3], f_current_edge_vertex_1[3], f_initial_edge_vertex_0[3], f_initial_edge_vertex_1[3],
//	//	f_current_compare_edge_vertex_0[3], f_current_compare_edge_vertex_1[3], f_initial_compare_edge_vertex_0[3],
//	//	f_initial_compare_edge_vertex_1[3];
//	//make_vector(current_edge_vertex_0, f_current_edge_vertex_0);
//	//make_vector(current_edge_vertex_1, f_current_edge_vertex_1);
//	//make_vector(initial_edge_vertex_0, f_initial_edge_vertex_0);
//	//make_vector(initial_edge_vertex_1, f_initial_edge_vertex_1);
//	//make_vector(initial_compare_edge_vertex_0, f_initial_compare_edge_vertex_0);
//	//make_vector(initial_compare_edge_vertex_1, f_initial_compare_edge_vertex_1);
//	//make_vector(current_compare_edge_vertex_0, f_current_compare_edge_vertex_0);
//	//make_vector(current_compare_edge_vertex_1, f_current_compare_edge_vertex_1);
//
//
//	SUB(e0, initial_compare_edge_vertex_1, initial_edge_vertex_0);//d-a
//	SUB(e0_, current_compare_edge_vertex_1, current_edge_vertex_0);
//	SUB(e1, initial_edge_vertex_1, initial_edge_vertex_0);//b-a
//	SUB(e1_, current_edge_vertex_1, current_edge_vertex_0);
//	SUB(e2, initial_compare_edge_vertex_0, initial_edge_vertex_0);//c-a
//	SUB(e2_, current_compare_edge_vertex_0, current_edge_vertex_0);
//
//	double a[3];
//	CROSS(a, e1, e2);
//	double d = DOT(e0, a);
//
//	double b[3], c[3];
//	CROSS(c, e1, e2_);
//	CROSS(b, e1_, e2);
//	SUM_(b, c);
//	CROSS(c, e1_, e2_);
//	double b1 = DOT(e0_, c);
//	double b2 = DOT(e0_, a);
//	double b3 = DOT(e0, b);
//	double b4 = DOT(e0_, b);
//	double b5 = DOT(e0, c);
//
//
//	double a3 = -d + b1 + b2 + b3 - b4 - b5;
//	double a2 = 3 * d + b5 - 2 * b3 + b4 - 2 * b2;
//	double a1 = -3 * d + b2 + b3;
//
//	std::vector<double>time;
//	time.reserve(7);
//
//	double d_c[3];
//	SUB(d_c, initial_compare_edge_vertex_1, initial_compare_edge_vertex_0);//d-c
//
//	double u_0_0[3], u_0_1[3], u_1_0[3], u_1_1[3];
//	SUB(u_0_0, current_edge_vertex_0, initial_edge_vertex_0);
//	SUB(u_0_1, current_edge_vertex_1, initial_edge_vertex_1);
//	SUB(u_1_0, current_compare_edge_vertex_0, initial_compare_edge_vertex_0);
//	SUB(u_1_1, current_compare_edge_vertex_1, initial_compare_edge_vertex_1);
//
//	bool is_collide = false;
//
//
//	if (a3 > -NEAR_ZERO && a3<NEAR_ZERO
//		&& a2 > -NEAR_ZERO && a2<NEAR_ZERO
//		&& a1 > -NEAR_ZERO && a1 < NEAR_ZERO) {
//		if (d<NEAR_ZERO && d>-NEAR_ZERO) { //the twos edge always coplanar
//			double v_0[3], v_1[3];
//			//initial_edge_vertex_0
//			SUB(v_0, initial_edge_vertex_0, initial_compare_edge_vertex_0);
//			SUB(v_1, initial_edge_vertex_0, initial_compare_edge_vertex_1);
//			if (pointEdgeCollisionTime(t, u_0_0, u_1_0, u_1_1, d_c, v_0, v_1, tolerance_2)) {
//				time.push_back(t);
//			}
//			if (pointPointCollisionTime(t, v_0, u_1_0, u_0_0, tolerance_2)) {
//				time.push_back(t);
//			}
//			if (pointPointCollisionTime(t, v_1, u_1_1, u_0_0, tolerance_2)) {
//				time.push_back(t);
//			}
//
//			//initial_edge_vertex_1
//			SUB(v_0, initial_edge_vertex_1, initial_compare_edge_vertex_0);
//			SUB(v_1, initial_edge_vertex_1, initial_compare_edge_vertex_1);
//			if (pointEdgeCollisionTime(t, u_0_1, u_1_0, u_1_1, d_c, v_0, v_1, tolerance_2)) {
//				time.push_back(t);
//			}
//			if (pointPointCollisionTime(t, v_0, u_1_0, u_0_1, tolerance_2)) {
//				time.push_back(t);
//			}
//			if (pointPointCollisionTime(t, v_1, u_1_1, u_0_1, tolerance_2)) {
//				time.push_back(t);
//			}
//
//			//initial_compare_edge_vertex_0
//			SUB(v_0, initial_compare_edge_vertex_0, initial_edge_vertex_0);
//			SUB(v_1, initial_compare_edge_vertex_0, initial_edge_vertex_1);
//			if (pointEdgeCollisionTime(t, u_1_0, u_0_0, u_0_1, e1, v_0, v_1, tolerance_2)) {
//				time.push_back(t);
//			}
//			//initial_compare_edge_vertex_1
//			SUB(v_0, initial_compare_edge_vertex_1, initial_edge_vertex_0);
//			SUB(v_1, initial_compare_edge_vertex_1, initial_edge_vertex_1);
//			if (pointEdgeCollisionTime(t, u_1_1, u_0_0, u_0_1, e1, v_0, v_1, tolerance_2)) {
//				time.push_back(t);
//			}
//			////check when t=0,the two edges are coplanar if the two edges are collide
//			//double e_temp0[3], compare_direction0[3], compare_direction1[3];
//			//SUB(e_temp0, initial_compare_edge_vertex_0, initial_edge_vertex_1);//c-b		
//			//CROSS(compare_direction0, d_c, e_temp0);
//			//CROSS(compare_direction1, d_c, e2);
//			//double d2 = DOT(compare_direction0, compare_direction1);
//			////SUB(e_temp0, initial_compare_edge_vertex_1, initial_edge_vertex_0);//d-a
//			//CROSS(compare_direction0, e1, e0);
//			//CROSS(compare_direction1, e1, e2);
//			//double d1 = DOT(compare_direction0, compare_direction1);		
//			//if ( d1 < 0 && d2<0) { //this means the vertices of one edge is on different sides of the other edge
//			//	t = 0;
//			//	return true;
//			//}
//			//// initial_compare_edge_vertex_0  initial_edge
//			//if (pointEdgeIsClose(e2, e1, tolerance_2)) {
//			//	t = 0;
//			//	return true;
//			//}
//			//// initial_compare_edge_vertex_1  initial_edge
//			//if (pointEdgeIsClose(e0, e1, tolerance_2)) {
//			//	t = 0;
//			//	return true;
//			//}
//			////initial_edge_vertex_0 initial_compare_edge
//			//SUB(e_temp0, initial_edge_vertex_0, initial_compare_edge_vertex_0);//a-c
//			//if (pointEdgeIsClose(e_temp0, d_c, tolerance_2)) {
//			//	t = 0;
//			//	return true;
//			//}
//			//SUB(e_temp0, initial_edge_vertex_1, initial_compare_edge_vertex_0);//b-c
//			//if (pointEdgeIsClose(e_temp0, d_c, tolerance_2)) {
//			//	t = 0;
//			//	return true;
//			//}
//		}
//		else {
//			is_collide = false;
//		}
//
//	}
//	else {
//		double t0 = 2.0; double t1 = 2.0; double t2 = 2.0;
//		if (solveEquation(t0,t1,t2, a3, a2, a1, d)) {
//			//if(estimate){
//			//	if (tight_CCD.insideTest(f_initial_edge_vertex_0, f_initial_edge_vertex_1, f_initial_compare_edge_vertex_0,
//			//		f_initial_compare_edge_vertex_1, f_current_edge_vertex_0, f_current_edge_vertex_1, f_current_compare_edge_vertex_0,
//			//		f_current_compare_edge_vertex_1, f_n0, f_n1, f_cross_for_CCD, true)) {
//			//		time.push_back(t);
//			//		/*if (t < 1e-11) {
//			//			std::cout << "edge-edge approx " << t << std::endl;
//			//			std::cout << std::setprecision(20) << initial_edge_vertex_0[0] << ", " << initial_edge_vertex_0[1] << ", " << initial_edge_vertex_0[2] << std::endl;
//			//			std::cout << std::setprecision(20) << initial_edge_vertex_1[0] << ", " << initial_edge_vertex_1[1] << ", " << initial_edge_vertex_1[2] << std::endl;
//			//			std::cout << std::setprecision(20) << initial_compare_edge_vertex_0[0] << ", " << initial_compare_edge_vertex_0[1] << ", " << initial_compare_edge_vertex_0[2] << std::endl;
//			//			std::cout << std::setprecision(20) << initial_compare_edge_vertex_1[0] << ", " << initial_compare_edge_vertex_1[1] << ", " << initial_compare_edge_vertex_1[2] << std::endl;
//			//			std::cout << std::setprecision(20) << current_edge_vertex_0[0] << ", " << current_edge_vertex_0[1] << ", " << current_edge_vertex_0[2] << std::endl;
//			//			std::cout << std::setprecision(20) << current_edge_vertex_1[0] << ", " << current_edge_vertex_1[1] << ", " << current_edge_vertex_1[2] << std::endl;
//			//			std::cout << std::setprecision(20) << current_compare_edge_vertex_0[0] << ", " << current_compare_edge_vertex_0[1] << ", " << current_compare_edge_vertex_0[2] << std::endl;
//			//			std::cout << std::setprecision(20) << current_compare_edge_vertex_1[0] << ", " << current_compare_edge_vertex_1[1] << ", " << current_compare_edge_vertex_1[2] << std::endl;
//			//		}*/
//			//	}
//			//}
//			//else {
//			if (edgeEdgeCheckInside(t0, initial_edge_vertex_0, initial_edge_vertex_1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1,
//				u_0_0, u_0_1, u_1_0, u_1_1)) {
//				time.push_back(t0);
//				/*if (t < 1e-11) {
//					std::cout << "edge-edge inside " << t << std::endl;
//				}*/
//			}
//			else {
//				if (t1 < 1.5) {
//					if (edgeEdgeCheckInside(t1, initial_edge_vertex_0, initial_edge_vertex_1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1,
//						u_0_0, u_0_1, u_1_0, u_1_1)) {
//						time.push_back(t1);
//					}
//					else {
//						if (t2 < 1.5) {
//							if (edgeEdgeCheckInside(t2, initial_edge_vertex_0, initial_edge_vertex_1, initial_compare_edge_vertex_0, initial_compare_edge_vertex_1,
//								u_0_0, u_0_1, u_1_0, u_1_1)) {
//								time.push_back(t2);
//							}
//						}
//					}
//				}
//			}
//			//}
//		}
//	}
//
//	
//
//	if (!time.empty()) {
//		t = time[0];
//		for (int i = 1; i < time.size(); ++i) {
//			if (t > time[i]) {
//				t = time[i];
//			}
//		}
//		t *= 0.8;
//		is_collide = true;
//	}
//
//	Vec3d initial_edge_0 = Vec3d(initial_edge_vertex_0[0], initial_edge_vertex_0[1], initial_edge_vertex_0[2]);
//	Vec3d current_edge_0 = Vec3d(current_edge_vertex_0[0], current_edge_vertex_0[1], current_edge_vertex_0[2]);
//	Vec3d initial_edge_1 = Vec3d(initial_edge_vertex_1[0], initial_edge_vertex_1[1], initial_edge_vertex_1[2]);
//	Vec3d current_edge_1 = Vec3d(current_edge_vertex_1[0], current_edge_vertex_1[1], current_edge_vertex_1[2]);
//	Vec3d initial_comapre_edge_0 = Vec3d(initial_compare_edge_vertex_0[0], initial_compare_edge_vertex_0[1], initial_compare_edge_vertex_0[2]);
//	Vec3d current_comapre_edge_0 = Vec3d(current_compare_edge_vertex_0[0], current_compare_edge_vertex_0[1], current_compare_edge_vertex_0[2]);
//	Vec3d initial_comapre_edge_1 = Vec3d(initial_compare_edge_vertex_1[0], initial_compare_edge_vertex_1[1], initial_compare_edge_vertex_1[2]);
//	Vec3d current_comapre_edge_1 = Vec3d(current_compare_edge_vertex_1[0], current_compare_edge_vertex_1[1], current_compare_edge_vertex_1[2]);
//	rootparity::RootParityCollisionTest root_parity(initial_edge_0, initial_edge_1, initial_comapre_edge_0, initial_comapre_edge_1,
//		current_edge_0, current_edge_1, current_comapre_edge_0, current_comapre_edge_1, true);
//	if (root_parity.edge_edge_collision())
//	{
//		if (!is_collide) {
//			std::cout << "EE actually collide but does not test out " << std::endl;
//			std::cout << "ee not equation coeff " << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
//		}
//	}
//
//	testEdgeInsideOutside(t, initial_edge_vertex_0, u_0_0, initial_edge_vertex_1, u_0_1, initial_compare_edge_vertex_0, u_1_0, initial_compare_edge_vertex_1, u_1_1);
//
//	if (is_collide) {
//		std::cout << "ee equation coeff " << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
//	}
//
//	return is_collide;
//}


bool ApproxCCD::checkInside(double t, double* v0, double* v1, double* v2, double* v3,
	double* e0, double* e1, double* e2, double* e3)
{
	double current_v0[3], current_v1[3], current_v2[3];

	double e0_1[3], e0_2[3], e0_3[3];
	//for (unsigned int i = 0; i < 3; ++i) {
	//	e0_1[i] = v1[i] - v0[i] + t * (e1[i] - e0[i]);
	//	e0_2[i] = v2[i] - v0[i] + t * (e2[i] - e0[i]);
	//	e0_3[i] = v3[i] - v0[i] + t * (e3[i] - e0[i]);
	//}

	OBTAIN_CURRENT_POS(e0_1, v0, v1, e0, e1, t);
	OBTAIN_CURRENT_POS(e0_2, v0, v2, e0, e2, t);
	OBTAIN_CURRENT_POS(e0_3, v0, v3, e0, e3, t);

	CROSS(current_v0, e0_1, e0_2);//here we use current_vi to store the norm to save storage
	CROSS(current_v1, e0_2, e0_3);

	if (DOT(current_v0, current_v1) < 0) {
		return false;
	}

	CROSS(current_v2, e0_3, e0_1);
	
	if (DOT(current_v1, current_v2) < 0) {
		return false;
	}
	if (DOT(current_v2, current_v0) < 0) {
		return false;
	}
	return true;
}



bool ApproxCCD::edgeEdgeCheckInside(double t, double* v0, double* v1, double* v2, double* v3,
	double* e0, double* e1, double* e2, double* e3)
{
	//double current_v0[3], current_v1[3], current_v2[3], current_v3[3];
	//for (int i = 0; i < 3; ++i) {
	//	current_v0[i] = v0[i] + t * e0[i];
	//	current_v1[i] = v1[i] + t * e1[i];
	//	current_v2[i] = v2[i] + t * e2[i];
	//	current_v3[i] = v3[i] + t * e3[i];
	//}
	double e_0_2[3], e_1_2[3], e_0_3[3], e_1_3[3], direct0[3], direct1[3];

	OBTAIN_CURRENT_POS(e_0_2, v2, v0, e2, e0, t);
	OBTAIN_CURRENT_POS(e_1_2, v2, v1, e2, e1, t);
	OBTAIN_CURRENT_POS(e_0_3, v3, v0, e3, e0, t);
	OBTAIN_CURRENT_POS(e_1_3, v3, v1, e3, e1, t);
	//SUB(e_0, current_v0, current_v2);
	//SUB(e_1, current_v1, current_v2);
	CROSS(direct0, e_0_2, e_1_2);
	//SUB(e_0, current_v0, current_v3);
	//SUB(e_1, current_v1, current_v3);
	CROSS(direct1, e_0_3, e_1_3);

	if (DOT(direct0, direct1) > 0) {
		return false;
	}
	//SUB(e_0, current_v2, current_v0);
	//SUB(e_1, current_v3, current_v0);
	CROSS(direct0, e_0_2, e_0_3);
	//SUB(e_0, current_v2, current_v1);
	//SUB(e_1, current_v3, current_v1);
	CROSS(direct1, e_1_2, e_1_3);
	if (DOT(direct0, direct1) > 0) {
		return false;
	}
	return true;
}


bool ApproxCCD::solveCubicEquation(double a, double b, double c, double d, double& t0, double& t1, double& t2)
{
	b = b / a;
	c = c / a;
	d = d / a;

	double alpha = (-2.0 * b * b * b + 9.0 * b * c - 27.0 * d) / 54.0;
	double beta = (b * b - 3.0 * c) / 9.0;
	double dum1 = beta * beta * beta;
	double delta = alpha * alpha - dum1;
	b /= 3.0;

	std::cout << "delta " << delta<< std::endl;

	if (delta > 0) {
		delta = sqrt(delta);
		//double R1 = ;
		//double R2 = ;
		t0 = -b + cbrt(alpha + delta) + cbrt(alpha - delta);
		std::cout << "t " << t0 << std::endl;
		if (t0 >= 0.0 && t0 <= 1.0) {
			return true;
		}
		else {
			t0 = 2.0;
		}
	}
	else if (delta == 0.0) {	//all roots are real, at least two are equal.
		alpha = cbrt(alpha);
		if ((-b + alpha + alpha) >= 0.0 && (-b + alpha + alpha) <= 1.0)
		{
			t0 = -b + alpha + alpha;
		}
		if ((alpha + b) <= 0.0 && (alpha + b) >= -1.0)
		{
			t1 = -alpha - b;
		}
		if (t1 < t0) {
			d = t1;
			t1 = t0;
			t0 = d;
		}
		if (t0 < 1.5) {
			return true;
		}
	}
	else { //all roots are real and different
		dum1 = acos(alpha / sqrt(dum1));
		alpha = 2.0 * std::sqrt(beta);
		t0 = -b + alpha * cos(dum1 / 3.0);
		t1 = -b + alpha * cos((dum1 + 2.0 * M_PI) / 3.0);
		t2 = -b + alpha * cos((dum1 + 4.0 * M_PI) / 3.0);

		if (t0 < 0.0 || t0>1.0) {
			t0 = 2.0;
		}
		if (t1 < 0.0 || t1>1.0) {
			t1 = 2.0;
		}
		if (t2 < 0.0 || t2>1.0) {
			t2 = 2.0;
		}
		sortABC(t0, t1, t2);
		if (t0 < 1.5) {
			return true;
		}
	}
	return false;

}


inline void ApproxCCD::sortABC(double& a, double& b, double& c)
{
	double d;
	if (a > b) { d = a; a = b; b = d; }
	if (a > c) { d = a; a = c; c = d; }
	if (b > c) { d = b; b = c; c = d; }
}

//t
bool ApproxCCD::solveEquation(double& t0, double& t1, double& t2, double a3, double a2, double a1, double d)
{
	//std::cout << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
	if (d<NEAR_ZERO && d>-NEAR_ZERO) {
		return solveQuadraticEquation(t0,t1, a3, a2, a1);
	}
	if (a3<NEAR_ZERO && a3>-NEAR_ZERO) {
		return solveQuadraticEquation(t0,t1, a2, a1, d);
	}
	return solveCubicEquation(a3, a2, a1, d, t0, t1, t2);
}

//bool ApproxCCD::solveEquation(double&t, double a3, double a2, double a1, double d)
//{
//	//std::cout << a3 << " " << a2 << " " << a1 << " " << d << std::endl;
//	if (d<NEAR_ZERO2 && d>-NEAR_ZERO2) {
//		return solveQuadraticEquation(t, a3, a2, a1);
//	}	
//	if (a3<NEAR_ZERO2 && a3>-NEAR_ZERO2) {
//		return solveQuadraticEquation(t, a2, a1, d);
//	}
//	double a_2 = 3 * a3;
//	double a_1 = 2 * a2;
//	double check_root = a_1 * a_1 - 4 * a_2 * a1;
//	if (check_root > 0) {
//		check_root = sqrt(check_root);
//		if (d < 0) {
//			a3 = -a3;
//			a2 = -a2;
//			a1 = -a1;
//			d = -d;
//			a_2 = -a_2;
//			a_1 = -a_1;
//		}
//		double i_0, i_1;
//		if (a3 > 0) {
//			i_0 = (-a_1 - check_root) / (2 * a_2);
//			i_1 = (-a_1 + check_root) / (2 * a_2);
//			//std::cout << "i_0 " << i_0 << " " << i_1 << std::endl;
//			if (i_1 < 0) {
//				return false;
//			}
//			if (i_0 > 1) {
//				return false;
//			}
//			if (i_0 < 0) {
//				if (i_1 > 1) {//interval [0, 1]
//					if (d * (a3 + a2 + a1 + d) < 0) {
//						double t_s = -a2 / a_2;
//						if (t_s < 0) {
//							t = -d / a1;
//						}
//						else {
//							if (t_s < 1) {
//								t = -d / (a_2 * t_s * t_s + a_1 * t_s + a1);
//							}
//							else {
//								t = -d / (a_2 + a_1 + a1);
//							}
//						}						
//						return true;
//					}						
//				}
//				else {//interval [0, i_1]
//					if (d * (a3 * i_1 * i_1 * i_1 + a2 * i_1 * i_1 + a1 * i_1 + d) < 0) {
//						double t_s = -a2 / a_2;
//						if (t_s < 0) {
//							t = -d / a1;
//						}
//						else {
//							t = -d / (a_2 * t_s * t_s + a_1 * t_s + a1);
//						}						
//						return true;
//					}
//				}
//				return false;
//			}		
//			double r_0 = a3 * i_0 * i_0 * i_0 + a2 * i_0 * i_0 + a1 * i_0 + d;
//			if (i_1 > 1) {//interval [i_0, 1]
//				if (r_0 * (a3 + a2 + a1 + d) < 0) {
//					double t_s = -a2 / a_2;
//					if (t_s > 1) {
//						t = -d / (a_2 + a_1 + a1);
//					}
//					else {
//						t = i_0 - r_0 / (a_2 * t_s * t_s + a_1 * t_s + a1);
//					}					
//					return true;
//				}						
//			}
//			else {//interval [i_0, i_1]
//				if (r_0 * (a3 * i_1 * i_1 * i_1 + a2 * i_1 * i_1 + a1 * i_1 + d) < 0) {
//					double t_s = -a2 / a_2;
//					t = i_0 - r_0 / (a_2 * t_s * t_s + a_1 * t_s + a1);
//					return true;
//				}
//			}
//			return false;			
//		}
//		else {
//			i_0 = (-a_1 + check_root) / (2 * a_2);
//			i_1 = (-a_1 - check_root) / (2 * a_2);
//			if (i_0 > 0) {
//				if (i_0 > 1) {//interval [0, 1]
//					if (d * (a3 + a2 + a1 + d) < 0) {
//						t = -d / a1;
//						return true;
//					}
//				}
//				else {//interval [0, i_0]
//					if (d * (a3 * i_0 * i_0 * i_0 + a2 * i_0 * i_0 + a1 * i_0 + d) < 0) {
//						t = -d / a1;
//						return true;
//					}
//				}
//				return false;
//			}
//			if (i_1 > 1) {
//				return false;
//			}
//			if (i_1 < 0) { //interval [0,1]
//				if (d * (a3 + a2 + a1 + d) < 0) {
//					t = -d / (a_2 + a_1 + a1);
//					return true;
//				}
//			}
//			else { //interval [i_1,1]
//				double r_i_1 = a3 * i_1 * i_1 * i_1 + a2 * i_1 * i_1 + a1 * i_1 + d;
//				if (r_i_1 * (a3 + a2 + a1 + d) < 0) {
//					t = i_1 - r_i_1 / (a_2 + a_1 + a1);
//					return true;
//				}
//			}
//			return false;		
//		}
//	}
//	else {
//		if (d * (a3 + a2 + a1 + d) < 0) {
//			double r_1 = a_2 + a_1 + a1;
//			if (fabs(a1) > fabs(r_1)) {
//				t = -d / a1;
//			}
//			else {
//				t = -d / r_1;
//			}
//			return true;
//		}
//		return false;
//	}	
//}


bool ApproxCCD::pointPointCollisionTime(double& t, double* e_1_0, double* e_0, double* e_1, double tolerance_2)
{
	std::vector<double> t_list;
	t_list.reserve(3);
	double t_temp;
	//x
	t_temp = e_1_0[0] / (e_0[0] - e_1[0]);
	if (t_temp > 0 && t_temp <= 1.0) {
		t_list.emplace_back(t_temp);
	}
	//y
	t_temp = e_1_0[1] / (e_0[1] - e_1[1]);
	if (t_temp > 0 && t_temp <= 1.0) {
		t_list.emplace_back(t_temp);
	}
	//z
	t_temp = e_1_0[2] / (e_0[2] - e_1[2]);
	if (t_temp > 0 && t_temp <= 1.0) {
		t_list.emplace_back(t_temp);
	}

	if (t_list.empty()) {
		return false;
	}
	if (t_list.size() == 1) {
		if (pointPointIsClose(t_list[0], e_1_0, e_0, e_1, tolerance_2)) {
			t = t_list[0];
			return true;
		}
		return false;
	}
	std::sort(t_list.begin(), t_list.end());
	for (int i = 0; i < t_list.size(); ++i) {
		if (pointPointIsClose(t_list[i], e_1_0, e_0, e_1, tolerance_2)) {
			t = t_list[i];
			return true;
		}
	}
	return false;
}


void ApproxCCD::pointEdge2D(std::vector<double>&t, double* u, double* u0, double* u1, double* e_1_0, double* e0, double* e1, 
	double* u_0, double* u_1, double* u_10)
	// u is v(t1)-v(t0), e_1_0 is e_1(t0)-e_0(t0), e0 is v(t0)-e_0(t0)	u_0=v(t1)-v(t0)-(e_0(t1)-e_0(t0))
	//u_1_0 is e_1(t1)-e_1(t0)-(e_0(t1)-e_0(t0))
{

	double a, b, c;
	double t_temp[2];
	int root_number;
	double v_0[2], v_1[2], v_0_1[2];
	double x0t_x01t;
	double x01t_x01t;
	 

	//XY
	a = u_0[0] * u_10[1] - u_0[1] * u_10[0];
	b = u_0[0] * e_1_0[1] - u_0[1] * e_1_0[0]
		+ e0[0] * u_10[1] - e0[1] * u_10[0];
	c = e0[0] * e_1_0[1] - e0[1] * e_1_0[0];



	solveQuadratic(t_temp, a, b, c, root_number);
	for (int i = 0; i < root_number; ++i) {
		v_0[0] = e0[0] + t_temp[i] * u_0[0];
		v_0[1] = e0[1] + t_temp[i] * u_0[1];
		v_0_1[0] = e_1_0[0] + t_temp[i] * u_10[0];
		v_0_1[1] = e_1_0[1] + t_temp[i] * u_10[1];
		x0t_x01t = v_0[0] * v_0_1[0] + v_0[1] * v_0_1[1];
		x01t_x01t= v_0_1[0] * v_0_1[0] + v_0_1[1] * v_0_1[1];
		if (x0t_x01t >= 0 && x0t_x01t <= x01t_x01t) {
			t.push_back(t_temp[i]);
		}
	}	

	//YZ
	a = u_0[2] * u_10[1] - u_0[1] * u_10[2];
	b = u_0[2] * e_1_0[1] - u_0[1] * e_1_0[2]
		+ e0[2] * u_10[1] - e0[1] * u_10[2];
	c = e0[2] * e_1_0[1] - e0[1] * e_1_0[2];
	solveQuadratic(t_temp, a, b, c, root_number);
	for (int i = 0; i < root_number; ++i) {
		v_0[0] = e0[1] + t_temp[i] * u_0[1];
		v_0[1] = e0[2] + t_temp[i] * u_0[2];
		v_0_1[0] = e_1_0[1] + t_temp[i] * u_10[1];
		v_0_1[1] = e_1_0[2] + t_temp[i] * u_10[2];
		x0t_x01t = v_0[0] * v_0_1[0] + v_0[1] * v_0_1[1];
		x01t_x01t = v_0_1[0] * v_0_1[0] + v_0_1[1] * v_0_1[1];
		if (x0t_x01t >= 0 && x0t_x01t <= x01t_x01t) {
			t.push_back(t_temp[i]);
		}
	}

	//XZ
	a = u_0[2] * u_10[0] - u_0[0] * u_10[2];
	b = u_0[2] * e_1_0[0] - u_0[0] * e_1_0[2]
		+ e0[2] * u_10[0] - e0[0] * u_10[2];
	c = e0[2] * e_1_0[0] - e0[0] * e_1_0[2];
	solveQuadratic(t_temp, a, b, c, root_number);
	for (int i = 0; i < root_number; ++i) {
		v_0[0] = e0[0] + t_temp[i] * u_0[0];
		v_0[1] = e0[2] + t_temp[i] * u_0[2];
		v_0_1[0] = e_1_0[0] + t_temp[i] * u_10[0];
		v_0_1[1] = e_1_0[2] + t_temp[i] * u_10[2];
		x0t_x01t = v_0[0] * v_0_1[0] + v_0[1] * v_0_1[1];
		x01t_x01t = v_0_1[0] * v_0_1[0] + v_0_1[1] * v_0_1[1];
		if (x0t_x01t >= 0 && x0t_x01t <= x01t_x01t) {
			t.push_back(t_temp[i]);
		}
	}
}

void ApproxCCD::solveQuadratic(double* t, double a, double b, double c, int& root_number)
{
	root_number = 0;
	if (a<NEAR_ZERO2 && a>-NEAR_ZERO2) {
		if (b<NEAR_ZERO2 && b>-NEAR_ZERO2) {
			return;			
		}
		double time = -c / b;
		if (0 < time && time <= 1) {
			t[root_number++] = time;
		}		
		return;
	}
	if (a < 0) { a = -a; b = -b; c = -c; }	
	double delta = b * b - 4 * a * c;
	if (delta <= 0){
		if (-b > 0 && -b < 2 * a)  t[root_number++] = -b / (2 * a);
		return;
	}
	if (b <= 0){
		double temp = -b + sqrt(delta);
		double twice_c = 2 * c;
		double twice_a = 2 * a;
		if (twice_c > 0 && twice_c < temp)	t[root_number++] = twice_c / temp;
		if (temp < twice_a)				t[root_number++] = temp / twice_a;
	}
	else{
		double temp = -b - sqrt(delta);
		double twice_c = 2 * c;
		double twice_a = 2 * a;
		if (twice_a < temp)				t[root_number++] = temp / twice_a;
		if (twice_c < 0 && temp < twice_c)	t[root_number++] = twice_c / temp;
	}
}


unsigned int ApproxCCD::getQuadRoots(double a, double b, double c, double& t0, double& t1)
{
	unsigned int roots = 0;
	int sign = 1;

	if (b < 0) {
		sign = -1;
	}
	double D= b * b - 4.0 * a * c;
	if (D >= 0) {
		roots = 2;
		double q = -0.5 * (b + sign * sqrt(D)); //here we try to use sum for two positive values, avoid use subtract
		t0 = q / a;
		t1 = c / q;
		if (t0 > t1) {
			a = t1;
			t1 = t0;
			t0 = a;
		}			
	}
	return roots;
}

bool ApproxCCD::solveQuadraticEquation(double& t0, double& t1, double a2, double a1, double a0)
{
	if (a0<NEAR_ZERO2 && a0>-NEAR_ZERO2) {
		if (a2<NEAR_ZERO2 && a2>-NEAR_ZERO2) {
			return false;
			//std::cout << "all zero "<<a2<<" "<<a1<<" "<<a0 << std::endl;
		}
		t0 = -a1 / a2;
		if (t0 > 0.0 && t0 <= 1.0) {
			return true;
		}
		else {
			return false;
		}
	}
	if (a2<NEAR_ZERO2 && a2>-NEAR_ZERO2) {
		if (a1<NEAR_ZERO2 && a1>-NEAR_ZERO2) {
			return false;
		}
		t0 = -a0 / a1;
		if (t0 > 0.0 && t0 <= 1.0) {
			return true;
		}
		else {
			return false;
		}
	}
	else {
		double temp = a1 * a1 - 4 * a2 * a0;
		if (temp < 0) {
			return false;
		}
		else {
			temp = sqrt(temp);
			if (a2 < 0) {
				a2 = -a2;
				a1 = -a1;
				a0 = -a0;
			}
			t0 = (-a1 - temp) / (a2 + a2);
			if (t0 > 1.0) {
				return false;
			}
			else if (t0 < 0.0) {
				t0 = 2.0;
			}
			t1 = (-a1 + temp) / (a2 + a2);
			if (t1 < 0.0) {
				return false;
			}
			if (t1 > 1.0) {
				t1 = 2.0;
			}
			if (t1 < t0) {
				temp = t1;
				t1 = t0;
				t0 = temp;
			}
			if (t0 <= 1.0 && t0 >= 0.0) {
				return true;
			}
			return false;
		}
	}
}


//log:
//Here, we first test vertex-edge and vertex-vertex. 
//



bool ApproxCCD::testInsideOutside(double t, double* initial_pos, double* u, double* initial_triangle_0,
	double* u_0, double* initial_triangle_1, double* u_1, double* initial_triangle_2, double* u_2)
{
	double result_v[3];
	double result_p0[3];
	double result_p1[3];
	double result_p2[3];
	POINT_ON_EDGE(result_v, 1.0, t, initial_pos, u);
	POINT_ON_EDGE(result_p0, 1.0, t, initial_triangle_0, u_0);
	POINT_ON_EDGE(result_p1, 1.0, t, initial_triangle_1, u_1);
	POINT_ON_EDGE(result_p2, 1.0, t, initial_triangle_2, u_2);
	double e_result1[3], e_result2[3];
	SUB(e_result1, result_p0, result_p2);
	SUB(e_result2, result_p1, result_p2);
	double norm_result[3];
	CROSS(norm_result, e_result1, e_result2);
	//normalize(norm_result);
	SUB(e_result1, initial_triangle_0, initial_triangle_2);
	SUB(e_result2, initial_triangle_1, initial_triangle_2);



	double norm_initial[3];
	CROSS(norm_initial, e_result1, e_result2);
	//normalize(norm_initial);
	double e_result[3];
	SUB(e_result, result_v, result_p0);
	double e_initial[3];
	SUB(e_initial, initial_pos, initial_triangle_0);	
	if (DOT(e_initial, norm_initial) * DOT(e_result, norm_result) < 0) {
		t = 0.553698;;
		POINT_ON_EDGE(result_v, 1.0, t, initial_pos, u);
		POINT_ON_EDGE(result_p0, 1.0, t, initial_triangle_0, u_0);
		POINT_ON_EDGE(result_p1, 1.0, t, initial_triangle_1, u_1);
		POINT_ON_EDGE(result_p2, 1.0, t, initial_triangle_2, u_2);		
		SUB(e_result1, result_p0, result_p2);
		SUB(e_result2, result_p1, result_p2);
		CROSS(norm_result, e_result1, e_result2);
		SUB(e_result, result_v, result_p1);

		std::cout << "vt side indicate " << DOT(e_initial, norm_initial) << " " << DOT(e_result, norm_result) << " " << DOT(norm_result, norm_result) << std::endl;
		std::cout << "vt side has changed " << std::endl;

		

		//POINT_ON_EDGE(result_v, 1.0, 1.0, initial_pos, u);
		//POINT_ON_EDGE(result_p0, 1.0, 1.0, initial_triangle_0, u_0);
		//POINT_ON_EDGE(result_p1, 1.0, 1.0, initial_triangle_1, u_1);
		//POINT_ON_EDGE(result_p2, 1.0, 1.0, initial_triangle_2, u_2);


		double e[3];
		SUB(e, result_p0, result_p1);
		double norm_test[3];
		CROSS(norm_test, norm_result, e);
		normalize(norm_test);
		SUB(e_result, result_v, result_p1);
		std::cout << "distace0 " << DOT(norm_test, e_result) << std::endl;

		SUB(e, initial_triangle_0, initial_triangle_1);
		CROSS(norm_test, norm_initial, e);
		normalize(norm_test);
		std::cout << "distace1 " << DOT(norm_test, e_initial) << std::endl;

		return false;
	}
	return true;
}

void ApproxCCD::testEdgeInsideOutside(double t, double* initial_pos_0_0, double* u_0_0, double* initial_pos_0_1,
	double* u_0_1, double* initial_triangle_1_0, double* u_1_0, double* initial_triangle_1_1, double* u_1_1)
{
	double result_v[3];
	double result_p0[3];
	double result_p1[3];
	double result_p2[3];
	POINT_ON_EDGE(result_v, 1.0, t, initial_pos_0_0, u_0_0);
	POINT_ON_EDGE(result_p0, 1.0, t, initial_pos_0_1, u_0_1);
	POINT_ON_EDGE(result_p1, 1.0, t, initial_triangle_1_0, u_1_0);
	POINT_ON_EDGE(result_p2, 1.0, t, initial_triangle_1_1, u_1_1);
	double e_result1[3], e_result2[3];
	SUB(e_result1, result_v, result_p0);
	SUB(e_result2, result_p1, result_p2);
	double norm_result[3];
	CROSS(norm_result, e_result1, e_result2);
	normalize(norm_result);
	SUB(e_result1, initial_pos_0_0, initial_pos_0_1);
	SUB(e_result2, initial_triangle_1_0, initial_triangle_1_1);
	double norm_initial[3];
	CROSS(norm_initial, e_result1, e_result2);
	normalize(norm_initial);
	double e_result[3];
	SUB(e_result, result_v, result_p1);
	double e_initial[3];
	SUB(e_initial, initial_pos_0_0, initial_triangle_1_0);	
	if (DOT(e_initial, norm_initial) * DOT(e_result, norm_result) < 0) {
		std::cout << "ee side indicate " << DOT(e_initial, norm_initial) << " " << DOT(e_result, norm_result) << " " << DOT(norm_result, norm_result) << std::endl;
		std::cout << "ee side has changed " << std::endl;
	}
}


//edge-edge
void ApproxCCD::test()
{
	//double a = (double)rand() / RAND_MAX * 10.0;
	//double b = (double)rand() / RAND_MAX;
	//double c = (double)rand() / RAND_MAX;
	//double d = (double)rand() / RAND_MAX * 0.02;
	////double a = 10.0;
	////double b = 35.0;
	////double c = 25.0;
	////double d = 4.0;
	//double x1, x2, x3;
	//double x11, x12, x13;
	//for (unsigned int i = 0; i < 10; ++i) {
	//	a = (double)rand() / RAND_MAX * 10.0;
	//	b = (double)rand() / RAND_MAX;
	//	c = (double)rand() / RAND_MAX;
	//	d = (double)rand() / RAND_MAX * 0.02;
	//	std::cout << "a: " << a << ", b: " << b << ", c: " << c << ", d: " << d << std::endl;
	//	solveCubicEquation(a, b, c, d, x1, x2, x3);
	//	std::cout << "x1: " << x1 << ", x2: " << x2 << ", x3: " << x3 << std::endl;
	//	cubicsolve(a, b, c, d, x11, x12, x13);
	//	std::cout << "x11: " << x11 << ", x12: " << x12 << ", x13: " << x13 << std::endl;
	//	std::cout << "===" << std::endl;
	//}
	//double xx1 = a * x1 * x1 * x1 + b * x1 * x1 + c * x1 + d;
	//double xx2 = a * x2 * x2 * x2 + b * x2 * x2 + c * x2 + d;
	//double xx3 = a * x3 * x3 * x3 + b * x3 * x3 + c * x3 + d;
	//std::cout << "error1: " << std::abs(xx1) << ", error2: " << std::abs(xx2) << ", error3: " << std::abs(xx3) << std::endl;
	//tight_CCD.testInside();
	////edge 1
	//double initial_pos[3] = { 1.0,1.0,1.0 };
	//double initial_triangle_0[3] = { -1.0,0.0,-1.0 };
	//
	////edge 2
	//double initial_triangle_1[3] = { -1.0,0.0,1.0 };
	//double initial_triangle_2[3] = { 1.0,0.0,-1.0 };
	//double current_pos[3] = { 1.0,0.3,1.0 };
	//double current_triangle_0[3] = { -1.0,-0.1,-1.0 };
	//double current_triangle_1[3] = { -0.5,0.5,1.0 };
	//double current_triangle_2[3] = { 0.5,0.7,-1.0 };
	//edge 1
	double initial_pos[3] = { 1.0,1.0,1.0 };
	double initial_triangle_0[3] = { -1.0,0.0,-1.0 };
	//edge 2
	double initial_triangle_1[3] = { -1,-0.5,0.3 };
	double initial_triangle_2[3] = { 0.0,-0.5,-0.8 };
	double current_pos[3] = { -0.5,0.3,1.0 };
	double current_triangle_0[3] = { -1.0,-0.1,-0.5 };
	double current_triangle_1[3] = { -0.5,0.5,0.0 };
	double current_triangle_2[3] = { 0.5,0.5,-1.0 };
	double t;
	double e_0[3], e_1[3], e_2[3];
	SUB(e_0, initial_pos, initial_triangle_0);
	SUB(e_1, initial_pos, initial_triangle_1);
	SUB(e_2, initial_pos, initial_triangle_2);
	double u[3], u_0[3], u_1[3], u_2[3];
	SUB(u, current_pos, initial_pos);
	SUB(u_0, current_triangle_0, initial_triangle_0);
	SUB(u_1, current_triangle_1, initial_triangle_1);
	SUB(u_2, current_triangle_2, initial_triangle_2);
	double e_[3];
	SUB(e_, initial_triangle_1, initial_triangle_0);
	double tolerance_2 = 0.9;
	
	double time1 = 1.0;
	time_t t1 = clock();
	//for (unsigned int i = 0; i < 100; ++i) {
		for (unsigned int j = 0; j < 112580; ++j) {
			if (edgeEdgeCCD(t, current_pos, current_triangle_0, initial_pos, initial_triangle_0,
				current_triangle_1, current_triangle_2, initial_triangle_1, initial_triangle_2, tolerance_2)) {
				if (time1 > t) {
					time1 = t;
				}
				//std::cout << t << std::endl;
				//t *= 0.9;
			/*	double result_v[3];
				double result_p0[3];
				double result_p1[3];
				double result_p2[3];
				POINT_ON_EDGE(result_v, 1.0, t, initial_pos, u);
				POINT_ON_EDGE(result_p0, 1.0, t, initial_triangle_0, u_0);
				POINT_ON_EDGE(result_p1, 1.0, t, initial_triangle_1, u_1);
				POINT_ON_EDGE(result_p2, 1.0, t, initial_triangle_2, u_2);
				double e_result1[3], e_result2[3];
				SUB(e_result1, result_v, result_p0);
				SUB(e_result2, result_p1, result_p2);
				double norm_result[3];
				CROSS(norm_result, e_result1, e_result2);
				SUB(e_result1, initial_pos, initial_triangle_0);
				SUB(e_result2, initial_triangle_1, initial_triangle_2);
				double norm_initial[3];
				CROSS(norm_initial, e_result1, e_result2);
				double e_result[3];
				SUB(e_result, result_v, result_p1);
				double e_initial[3];
				SUB(e_initial, initial_pos, initial_triangle_1);
				std::cout << "side indicate " << DOT(e_initial, norm_initial) << " " << DOT(e_result, norm_result)<<" "<< DOT(norm_result, norm_result) << std::endl;
				if (DOT(e_initial, norm_initial) * DOT(e_result, norm_result) < 0) {
					std::cout << "side has changed " << std::endl;
				}*/
			}
		}
	//}
	t1 = clock() - t1;
	std::cout <<"test time " << t1 << std::endl;
	//else {
	//	//std::cout << "does not collide " << std::endl;
	//}
	Vec3d initial_pos_a = Vec3d(initial_pos[0], initial_pos[1], initial_pos_a[2]);
	Vec3d current_pos_a = Vec3d(current_pos[0], current_pos[1], current_pos[2]);
	Vec3d initial_triangle_0_a = Vec3d(initial_triangle_0[0], initial_triangle_0[1], initial_triangle_0[2]);
	Vec3d initial_triangle_1_a = Vec3d(initial_triangle_1[0], initial_triangle_1[1], initial_triangle_1[2]);
	Vec3d initial_triangle_2_a = Vec3d(initial_triangle_2[0], initial_triangle_2[1], initial_triangle_2[2]);
	Vec3d current_triangle_0_a = Vec3d(current_triangle_0[0], current_triangle_0[1], current_triangle_0[2]);
	Vec3d current_triangle_1_a = Vec3d(current_triangle_1[0], current_triangle_1[1], current_triangle_1[2]);
	Vec3d current_triangle_2_a = Vec3d(current_triangle_2[0], current_triangle_2[1], current_triangle_2[2]);
	rootparity::RootParityCollisionTest root_parity(initial_pos_a, initial_triangle_0_a, initial_triangle_1_a, initial_triangle_2_a,
		current_pos_a, current_triangle_0_a, current_triangle_1_a, current_triangle_2_a, true);
	if (root_parity.edge_edge_collision())
	{
		std::cout << "actually collide " << std::endl;
	}

	Eigen::Vector3d initial_pos_a_ = Eigen::Vector3d(initial_pos[0], initial_pos[1], initial_pos[2]);
	Eigen::Vector3d current_pos_a_ = Eigen::Vector3d(current_pos[0], current_pos[1], current_pos[2]);
	Eigen::Vector3d initial_triangle_0_a_ = Eigen::Vector3d(initial_triangle_0[0], initial_triangle_0[1], initial_triangle_0[2]);
	Eigen::Vector3d initial_triangle_1_a_ = Eigen::Vector3d(initial_triangle_1[0], initial_triangle_1[1], initial_triangle_1[2]);
	Eigen::Vector3d initial_triangle_2_a_ = Eigen::Vector3d(initial_triangle_2[0], initial_triangle_2[1], initial_triangle_2[2]);
	Eigen::Vector3d current_triangle_0_a_ = Eigen::Vector3d(current_triangle_0[0], current_triangle_0[1], current_triangle_0[2]);
	Eigen::Vector3d current_triangle_1_a_ = Eigen::Vector3d(current_triangle_1[0], current_triangle_1[1], current_triangle_1[2]);
	Eigen::Vector3d current_triangle_2_a_ = Eigen::Vector3d(current_triangle_2[0], current_triangle_2[1], current_triangle_2[2]);

	double initial_distance_2 =
		CCD::internal::edgeEdgeDistanceUnclassified(initial_pos, initial_triangle_0, initial_triangle_1,
			initial_triangle_2);
	double min_distance_2 = (1.0 - tolerance_2) * (1.0 - tolerance_2) * initial_distance_2;

	CTCD ctcd;
	double tt;
	if (ctcd.edgeEdgeCTCD(initial_triangle_0_a_, initial_pos_a_, initial_triangle_2_a_, initial_triangle_1_a_,
		current_triangle_0_a_, current_pos_a_, current_triangle_2_a_, current_triangle_1_a_, min_distance_2, tt)) {
		std::cout << "ground truth ctcd collide " << tt << std::endl;
	}
	else {
		std::cout << "ground truth ctcd not collide " << std::endl;
	}

}

//vertex-triangle
//void ApproxCCD::test()
//{
//	//RootFinder rf;
//	//RootFinder_ rf_;
//	//double time[6];
//	//// We don't care one bit about these imaginary roots
//	//double zeroi[6];
//	//double time_[6];
//	//double zeroi_[6];
//	////memset()
//	//int degree = 6;
//	//double op[6];
//	//for (unsigned int k = 0; k < 100; ++k) {
//	//	for (unsigned int i = 0; i < 6; ++i) {
//	//		op[i] = (double)rand() / (double)RAND_MAX;
//	//	}
//	////memset(op, 0, 48);
//	//int root, root_;
//	////op[0] = 1;
//	////op[5] = -1;
//	//root=rf.rpoly(op, degree, time, zeroi);
//	//root_=rf_.rpoly(op, degree, time_, zeroi_);
//	//	//for (unsigned int i = 0; i < 6; ++i) {
//	//	//	if (abs(time[i] - time_[i]) > 1e-10) {
//	//std::cout << "===== " << std::endl;
//	//std::cout << root << " " << root_ << std::endl;
//	//for (unsigned int j = 0; j < 6; ++j) {
//	//	std::cout << time[j]- time_[j]<<" "<< time[j] << " " << time_[j] << " " << std::endl;
//	//}
//	//			//break;
//	//		//}
//	//	//}
//	//}	
//
//
//	double initial_pos[3] = { 0.5,0.0,0.3 };
//	double current_pos[3] = { 0.5,0.0,-0.3 };
//	double initial_triangle_0[3] = { 1.0,0.0,0.0 };
//	double initial_triangle_1[3] = { 0.0,1.0,0.0 };
//	double initial_triangle_2[3] = { 0.0,-1.0,0.0 };
//	double current_triangle_0[3] = { 1.0,0.0,0.0 };
//	double current_triangle_1[3] = { 0.0,1.0,0.0 };
//	double current_triangle_2[3] = { 0.0,-1.0,0.0 };
//	//double current_triangle_1[3] = { -1.0,0.1,-0.8 };
//	//double current_triangle_2[3] = { -1.0,-0.1,0.8 };
//	double initial_tri_normal[3]; double tri_normal[3];
//	double cross_for_CCD[3]; double temp[3];
//	double e1[3], e2[3];
//	double e3[3], e4[3];
//	SUB(e1, initial_triangle_1, initial_triangle_0);
//	SUB(e2, initial_triangle_2, initial_triangle_0);
//	CROSS(initial_tri_normal, e1, e2);
//	SUB(e3, current_triangle_1, current_triangle_0);
//	SUB(e4, current_triangle_2, current_triangle_0);
//	CROSS(tri_normal, e3, e4);
//	CROSS(cross_for_CCD, e1, e4);
//	CROSS(temp, e3, e2);
//	SUM_(cross_for_CCD, temp);
//	double t;
//	double e_0[3], e_1[3], e_2[3];
//	SUB(e_0, initial_pos, initial_triangle_0);
//	SUB(e_1, initial_pos, initial_triangle_1);
//	SUB(e_2, initial_pos, initial_triangle_2);
//	double u[3], u_0[3], u_1[3], u_2[3];
//	SUB(u, current_pos, initial_pos);
//	SUB(u_0, current_triangle_0, initial_triangle_0);
//	SUB(u_1, current_triangle_1, initial_triangle_1);
//	SUB(u_2, current_triangle_2, initial_triangle_2);
//	double e_[3];
//	SUB(e_, initial_triangle_1, initial_triangle_0);
//	double tolerance_2 = 0.9;
//	vertexTriangleCCD(t, initial_pos, current_pos, initial_triangle_0, current_triangle_0,
//		initial_triangle_1, current_triangle_1, initial_triangle_2, current_triangle_2, tolerance_2);
//	if (vertexTriangleCCD(t, initial_pos, current_pos, initial_triangle_0, current_triangle_0,
//		initial_triangle_1, current_triangle_1, initial_triangle_2, current_triangle_2, tolerance_2)) {
//		std::cout << t << std::endl;
//		//t *= 0.9;
//		double result_v[3];
//		double result_p0[3];
//		double result_p1[3];
//		double result_p2[3];
//		POINT_ON_EDGE(result_v, 1.0, t, initial_pos, u);
//		POINT_ON_EDGE(result_p0, 1.0, t, initial_triangle_0, u_0);
//		POINT_ON_EDGE(result_p1, 1.0, t, initial_triangle_1, u_1);
//		POINT_ON_EDGE(result_p2, 1.0, t, initial_triangle_2, u_2);
//		double e_result1[3], e_result2[3];
//		SUB(e_result1, result_p1, result_p0);
//		SUB(e_result2, result_p2, result_p0);
//		double norm_result[3];
//		CROSS(norm_result, e_result1, e_result2);
//		double e_result[3];
//		SUB(e_result, result_v, result_p0);
//		std::cout <<"side indicate "<< DOT(e_0, initial_tri_normal) << " " << DOT(e_result, norm_result) << std::endl;
//		if (DOT(e_0, initial_tri_normal) * DOT(e_result, norm_result) < 0) {			
//			std::cout << "side has changed " << std::endl;
//		}
//	}
//	else {
//		std::cout << "does not collide " << std::endl;
//	}
//	Vec3d initial_pos_a = Vec3d(initial_pos[0], initial_pos[1], initial_pos_a[2]);
//	Vec3d current_pos_a = Vec3d(current_pos[0], current_pos[1], current_pos[2]);
//	Vec3d initial_triangle_0_a = Vec3d(initial_triangle_0[0], initial_triangle_0[1], initial_triangle_0[2]);
//	Vec3d initial_triangle_1_a = Vec3d(initial_triangle_1[0], initial_triangle_1[1], initial_triangle_1[2]);
//	Vec3d initial_triangle_2_a = Vec3d(initial_triangle_2[0], initial_triangle_2[1], initial_triangle_2[2]);
//	Vec3d current_triangle_0_a = Vec3d(current_triangle_0[0], current_triangle_0[1], current_triangle_0[2]);
//	Vec3d current_triangle_1_a = Vec3d(current_triangle_1[0], current_triangle_1[1], current_triangle_1[2]);
//	Vec3d current_triangle_2_a = Vec3d(current_triangle_2[0], current_triangle_2[1], current_triangle_2[2]);
//	rootparity::RootParityCollisionTest root_parity(initial_pos_a, initial_triangle_0_a, initial_triangle_1_a, initial_triangle_2_a,
//		current_pos_a, current_triangle_0_a, current_triangle_1_a, current_triangle_2_a,false);
//	if (root_parity.point_triangle_collision())
//	{
//		std::cout << "actually collide " << std::endl;
//	}
//}

//void ApproxCCD::test()
//{
//	double initial_pos[3] = { -1.0001,0.1,0.0 };
//	double current_pos[3] = { -1.0001,-0.1,0.0 };
//	double initial_triangle_0[3]= { 1.0,0.0,0.0 };
//	double initial_triangle_1[3] = { -1.0,-0.1,-1.0 };
//	double initial_triangle_2[3] = { -1.0,0.1,1.0 };
//	double current_triangle_0[3]= { 1.0,0.0,0.0 };
//	double current_triangle_1[3]= { -1.0,0.1,-0.8 };
//	double current_triangle_2[3]= { -1.0,-0.1,0.8 };
//	double initial_tri_normal[3]; double tri_normal[3];
//	double cross_for_CCD[3]; double temp[3];
//	double e1[3], e2[3];
//	double e3[3], e4[3];
//	SUB(e1, initial_triangle_1, initial_triangle_0);
//	SUB(e2, initial_triangle_2, initial_triangle_0);
//	CROSS(initial_tri_normal, e1, e2);
//	SUB(e3, current_triangle_1, current_triangle_0);
//	SUB(e4, current_triangle_2, current_triangle_0);
//	CROSS(tri_normal, e3, e4);
//	CROSS(cross_for_CCD, e1, e4);
//	CROSS(temp, e3, e2);
//	SUM_(cross_for_CCD, temp);
//	double t;
//	double e_0[3], e_1[3], e_2[3];
//	SUB(e_0, initial_pos, initial_triangle_0);
//	SUB(e_1, initial_pos, initial_triangle_1);
//	SUB(e_2, initial_pos, initial_triangle_2);
//	double u[3], u_0[3], u_1[3], u_2[3];
//	SUB(u, current_pos, initial_pos);
//	SUB(u_0, current_triangle_0, initial_triangle_0);
//	SUB(u_1, current_triangle_1, initial_triangle_1);
//	SUB(u_2, current_triangle_2, initial_triangle_2);
//	double e_[3];
//	SUB(e_, initial_triangle_1, initial_triangle_0);
//	tolerance_2 = 1e-6;
//	if (pointEdgeCollisionTime(t, u, u_0, u_1, e_, e_0, e_1)) {
//		std::cout << "v e01 " << t<<std::endl;
//	}
//	SUB(e_, initial_triangle_2, initial_triangle_1);
//	if (pointEdgeCollisionTime(t, u, u_1, u_2, e_, e_1, e_2)) {
//		std::cout << "v e12 " << t << std::endl;
//	}
//	SUB(e_, initial_triangle_2, initial_triangle_0);
//	if (pointEdgeCollisionTime(t, u, u_0, u_2, e_, e_0, e_2)) {
//		std::cout << "v e02 " << t << std::endl;
//	}
//	if (pointPointCollisionTime(t, e_0, u_0, u)) {
//		std::cout << "v v0 " << t << std::endl;
//	}
//	if (pointPointCollisionTime(t, e_1, u_1, u)) {
//		std::cout << "v v1 " << t << std::endl;
//	}
//	if (pointPointCollisionTime(t, e_2, u_2, u)) {
//		std::cout << "v v2 " << t << std::endl;
//	}
//}

bool ApproxCCD::pointEdgeCollisionTime(double& t, double* u, double* u0, double* u1, double* e_1_0,
	double* e0, double* e1, double tolerance_2)
	// u is v(t1)-v(t0), e_1_0 is e_1(t0)-e_0(t0), e0 is v(t0)-e_0(t0)	//u_0=v(t1)-v(t0)-(e_0(t1)-e_0(t0))

{
	double u_0[3], u_1[3], u_10[3];
	SUB(u_0, u, u0);
	SUB(u_1, u, u1);
	SUB(u_10, u1, u0);


	std::vector<double>time;
	pointEdge2D(time, u, u0, u1, e_1_0, e0, e1, u_0, u_1, u_10);
	if (time.empty()) {
		return false;
	}
	std::sort(time.begin(), time.end());
	double v_0[3], v_1_0[3];
	for (int i = 0; i < time.size(); ++i) {
		t = time[i];
		v_0[0] = e0[0] + t * u_0[0];
		v_0[1] = e0[1] + t * u_0[1];
		v_0[2] = e0[2] + t * u_0[2];	
		v_1_0[0] = e_1_0[0] + t * u_10[0];
		v_1_0[1] = e_1_0[1] + t * u_10[1];
		v_1_0[2] = e_1_0[2] + t * u_10[2];		
		if (pointEdgeIsClose(v_0,v_1_0,tolerance_2)) {
			return true;
		}
	}
	return false;
}


bool ApproxCCD::pointEdgeIsClose(double* v_0, double* v_1_0, double tolerance_2)
//v_0 is v-e_0, v_1_0 is e_1-e_0
{
	double l_v_0 = DOT(v_0, v_0);
	double l_v_0_ = DOT(v_0, v_1_0);
	double l_v_1_0 = DOT(v_1_0, v_1_0);
	if (l_v_0 - l_v_0_ * l_v_0_ / l_v_1_0 > tolerance_2) {//use Pythagorean Theorem d^2=v_0^2-(v_0*v_1_0/||v_1_0||)^2
		return false;
	}	
	if (l_v_0_ > l_v_1_0 || l_v_0_ < 0)
		return false;
	return true;
}

bool ApproxCCD::pointPointIsClose(double t, double* e10_0, double* e_0, double* e_1, double tolerance_2)
{
	double distance[3];
	for (int i = 0; i < 3; ++i) {
		distance[i] = e10_0[i] + t * (e_1[i] - e_0[i]);
	}
	if (DOT(distance, distance) > tolerance_2) {
		return false;
	}
	return true;
}

inline void ApproxCCD::make_vector(double* v, floating* out)
{
	for (int i = 0; i < 3; i++) {
		out[i] = floating(v[i], 0);
	}
}

void ApproxCCD::cubicsolve(const double& a, const double& b, const double& c, const double& d, double& x1, double& x2, double& x3)
{
	if (a < 0.000001) { // end if a == 0 
		assert("The coefficient of the cube of x is 0. Please use the utility for a SECOND degree quadratic. No further action taken.");
		return;
	}
	if (d < 0.000001) { // End if d == 0, d is definitely a root of the function
		assert("One root is 0. Now divide through by x and use the utility for a SECOND degree quadratic to solve the resulting equation for the other two roots. No further action taken.");
		return;
	}
	double bb = b / a;
	double cc = c / a;
	double dd = d / a;
	double disc, q, r, dum1, s, t, term1, r13;
	q = (3.0 * cc - (bb * bb)) / 9.0;
	r = -(27.0 * dd) + bb * (9.0 * cc - 2.0 * (bb * bb));
	r /= 54.0;
	std::cout << "q: " << q << ", r: " << r << std::endl;
	disc = q * q * q + r * r;
	x1 = 0.0;			// The first root is always real.
	term1 = (bb / 3.0);
	if (disc > 0)       // one root real, two are complex
	{
		std::cout << "one root real, two are complex" << std::endl;
		s = r + std::sqrt(disc);
		s = ((s < 0) ? -std::pow(-s, (1.0 / 3.0)) : std::pow(s, (1.0 / 3.0)));
		t = r - std::sqrt(disc);
		t = ((t < 0) ? -std::pow(-t, (1.0 / 3.0)) : std::pow(t, (1.0 / 3.0)));
		x1 = -term1 + s + t;
		term1 += (s + t) / 2.0;
		x3 = x2 = -term1;
		term1 = std::sqrt(3.0) * (-t + s) / 2;
		x2 = term1;
		x3 = -term1;
		return;
	}
	// End if (disc > 0)
	// The remaining options are all real
	x3 = x2 = 0.0;
	if (disc == 0) // All roots real, at least two are equal.
	{
		std::cout << "All roots real, at least two are equal" << std::endl;
		r13 = ((r < 0) ? -std::pow(-r, (1.0 / 3.0)) : std::pow(r, (1.0 / 3.0)));
		x1 = -term1 + 2.0 * r13;
		x3 = x2 = -(r13 + term1);
		return;
	} // End if (disc == 0)
	// Only option left is that all roots are real and unequal (to get here, q < 0)
	q = -q;
	dum1 = q * q * q;
	dum1 = std::acos(r / std::sqrt(dum1));
	r13 = 2.0 * std::sqrt(q);
	x1 = -term1 + r13 * std::cos(dum1 / 3.0);
	x2 = -term1 + r13 * std::cos((dum1 + 2.0 * M_PI) / 3.0);
	x3 = -term1 + r13 * std::cos((dum1 + 4.0 * M_PI) / 3.0);
	std::cout << "All things are ok" << std::endl;
	return;
}


void ApproxCCD::checkInterval(double t1, double t2, double* op, int degree, std::vector<TimeInterval>& intervals, bool pos)
{
	// clamp values
	if (t1 < 0.0) {
		t1 = 0.0;
	}
	if (t2 < 0.0) {
		t2 = 0.0;
	}
	if (t1 > 1.0) {
		t1 = 1.0;
	}
	if (t2 > 1.0) {
		t2 = 1.0;
	}
	//this is to compute ((ax+b)x+c)x+....
	double tmid = 0.5 * (t2 + t1);
	double f = op[0];
	for (unsigned int i = 1; i <= degree; i++)
	{
		f *= tmid;
		f += op[i];
	}
	//for pos =true, that is to check the intersection with plane constructed by edge and triangle normal. 
	//the sign of f will be consistent in [t1,t2], thus we use midpoint.
	//if the vertex actually intersect with triangle, the sign must be positive.
	//here is to exclude when the sign is negative,
	if (pos && f >= 0)
		intervals.emplace_back(TimeInterval(t1, t2));
	else if (!pos && f <= 0)
		intervals.emplace_back(TimeInterval(t1, t2));
}

void ApproxCCD::findIntervals(double* op, unsigned int n, std::vector<TimeInterval>& intervals, bool pos)
{
	unsigned int roots = 0;
	unsigned int reducedDegree = n;
	double time[6];
	double zeroi[6];
	//std::uninitialized_fill(time, time + 6, 2.0);
	
	//normalized coeffs
	double maxval = 0;
	for (unsigned int i = 0; i <= n; i++) {
		if (maxval < fabs(op[i])) {
			maxval = fabs(op[i]);
		}
	}
	if (maxval != 0) {
		for (int i = 0; i <= n; i++) {
			op[i] /= maxval;
		}			
	}
	for (unsigned int i = 0; i < n; i++){
		if (op[i] == 0)
			reducedDegree--;
		else
			break;
	}	
	if (reducedDegree < n) {
		//move lower term coeff to higher term, this is for convenience
		for (int i = 0; i <= reducedDegree; i++) {
			op[i] = op[i + n - reducedDegree];
		}			
	}
	if (reducedDegree > 2) {
		if (!couldHaveRoots(op, reducedDegree, pos))
			return;

		roots = root_finder.rpoly(op,reducedDegree,time, zeroi);
	}
	else if (reducedDegree == 2) {
		roots = getQuadRoots(op[0], op[1], op[2], time[0], time[1]);
	}
	else if (reducedDegree == 1) {
		time[0] = -op[1] / op[0];
		roots = 1;
	}
	else { //reduced degree =0, for example: the vertex moving paralle to the plane
		// if  pos=true, the only coeff, coeff[0] must be positive (on the front side of the plane)
		if (!pos && op[0] <= 0 || (pos && op[0] >= 0)) {
			intervals.emplace_back(TimeInterval(0, 1.0));
		}
		return;
	}
	// need to check intervals
	//because time is the root, when pos =true, the vertex position in a time range won't change sideness of the plane
	if (roots > 0) {
		std::sort(time, time + roots);
		if (time[0] >= 0) {
			checkInterval(0, time[0], op, reducedDegree, intervals, pos);
		}
		for (unsigned int i = 0; i < roots - 1; ++i) {
			if (!((time[i] < 0 && time[i + 1] < 0) || (time[i] > 1.0 && time[i + 1] > 1.0))) {
				checkInterval(time[i], time[i + 1], op, reducedDegree, intervals, pos);
			}
		}
		if (time[roots - 1] <= 1.0)
			checkInterval(time[roots - 1], 1.0, op, reducedDegree, intervals, pos);
	}
	else {
		checkInterval(0.0, 1.0, op, reducedDegree, intervals, pos);
	}
}

bool ApproxCCD::couldHaveRoots(double* op, int degree, bool pos) {
	double result = 0;

	if (pos) {
		if (op[0] > 0) {
			result = op[0];
		}
		for (unsigned int i = 1; i < degree; i++)
		{
			if (op[i] > 0) {
				result += op[i];
			}
		}
		result += op[degree];
		return result >= 0;
	}
	else {
		if (op[0] < 0) {
			result = op[0];
		}
		for (unsigned int i = 1; i < degree; i++)
		{
			if (op[i] < 0) {
				result += op[i];
			}
		}
		result += op[degree];
		return result <= 0;
	}
}