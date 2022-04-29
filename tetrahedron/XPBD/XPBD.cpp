#include"XPBD.h"



XPBD::XPBD()
{
	gravity_ = 9.8;
	total_thread_num = std::thread::hardware_concurrency();
	sub_step_num = 1;
	iteration_number = 50;
	time_step = 1.0 / 100.0;
	damping_coe = 0.0;
}


void XPBD::initial()
{
	resetExternalForce();
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
	}
}

void XPBD::reset()
{
	resetExternalForce();
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
	}
}



void XPBD::setForXPBD(std::vector<Cloth>* cloth, std::vector<Tetrahedron>* tetrahedron, std::vector<Collider>* collider,
	Thread* thread, double* tolerance_ratio, DrawCulling* draw_culling_)
{
	this->cloth = cloth;
	this->tetrahedron = tetrahedron;
	this->collider = collider;
	sub_time_step = time_step / (double)sub_step_num;
	this->thread = thread;

	//collision.initial(cloth, collider, tetrahedron, thread, tolerance_ratio);

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
	gravity[0] = 0;
	gravity[1] = -gravity_;
	gravity[2] = 0;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		f_ext[i].resize(mesh_struct[i]->vertex_position.size());
		velocity[i].resize(mesh_struct[i]->vertex_position.size());
	}
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(velocity[i][0].data(), 0, 24 * velocity[i].size());
		memset(f_ext[i][0].data(), 0, 24 * f_ext[i].size());
	}
}

void XPBD::reorganzieDataOfObjects()
{
	vertex_position.resize(total_obj_num);
	initial_vertex_position.resize(total_obj_num);
	mesh_struct.resize(total_obj_num);
	vertex_index_begin_per_thread.resize(total_obj_num);
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		vertex_position[i]= cloth->data()[i].mesh_struct.vertex_position.data();
		initial_vertex_position[i]= cloth->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i]= &cloth->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i]= cloth->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
	}
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_position.data();
		initial_vertex_position[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_for_render.data();
		mesh_struct[i + cloth->size()] = &tetrahedron->data()[i].mesh_struct;
		vertex_index_begin_per_thread[i + cloth->size()] = tetrahedron->data()[i].mesh_struct.vertex_index_begin_per_thread.data();
	}
}

void XPBD::PBDsolve()
{
	memset(lambda.data(), 0, 8*lambda.size());
	memset(lambda_collision.data(), 0, 8*lambda_collision.size());
	thread->assignTask(this, SET_POS_PREDICT);
	for (unsigned int i = 0; i < iteration_number; ++i) {
		solveConstraint();
	}
	thread->assignTask(this, XPBD_VELOCITY);	
	updatePosition();
	updateNormal();
}



void XPBD::updatePosition()
{
	
	unsigned int vertex_num;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_num = mesh_struct[i]->vertex_position.size();
		memcpy(initial_vertex_position[i][0].data(), vertex_position[i][0].data(), 24 * vertex_num);
	}
}

void XPBD::solveConstraint()
{
	solveBendingConstraint();
	solveEdgeLengthConstraint();
	solveTetStrainConstraint();
}



void XPBD::updateNormal()
{
	TriangleMeshStruct* mesh_struct;
	for (unsigned int j = 0; j < cloth->size(); ++j) {
		mesh_struct = &(*cloth)[j].mesh_struct;
		thread->assignTask(mesh_struct, FACE_NORMAL_RENDER);
		thread->assignTask(mesh_struct, VERTEX_NORMAL_RENDER);
	}
	for (unsigned int j = 0; j < collider->size(); ++j) {
		mesh_struct = &(*collider)[j].mesh_struct;
		thread->assignTask(mesh_struct, FACE_NORMAL_RENDER);
		thread->assignTask(mesh_struct, VERTEX_NORMAL_RENDER);
	}
	TetrahedronMeshStruct* mesh_struct_;
	for (unsigned int j = 0; j < tetrahedron->size(); ++j) {
		mesh_struct_ = &(*tetrahedron)[j].mesh_struct;
		thread->assignTask(mesh_struct_, FACE_NORMAL_RENDER);
		thread->assignTask(mesh_struct_, VERTEX_NORMAL_RENDER);
	}

}


void XPBD::solveEdgeLengthConstraint()
{
	unsigned int size;
	MeshStruct* mesh_struct_;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	double stiffness;
	unsigned int* edge_vertex_index;
	double* mass_inv;
	//double delta_t = sub_time_step;
	double* lambda_ = lambda.data()+ constraint_index_start[1];
	//double damping_stiffness = damping_coe;
	//for (unsigned int i = 0; i < mesh_struct[0]->vertex_position.size(); ++i) {
	//	std::cout << vertex_position[0][i].data()[0]<<" "<< vertex_position[0][i][0] << std::endl;
	//}

	for (unsigned int i = 0; i < cloth->size(); ++i) {
		mesh_struct_ = mesh_struct[i];
		size = mesh_struct_->edge_length.size();
		vertex_pos = vertex_position[i];
		initial_vertex_pos = initial_vertex_position[i];

		//std::cout << vertex_position[i] << " " << vertex_pos << std::endl;

		stiffness = cloth->data()[i].length_stiffness[0];
		edge_vertex_index = mesh_struct_->edge_vertices.data();
		mass_inv = mesh_struct_->mass_inv.data();
		for (unsigned int j = 0; j < size; ++j) {			
			//std::cout << *lambda_ << " " << vertex_position[i][edge_vertex_index[j << 1]].data()[0] << " " << vertex_position[i][edge_vertex_index[(j << 1) + 1]].data()[0] << std::endl;
			XPBD_constraint.solveEdgeLengthConstraint(vertex_pos[edge_vertex_index[j << 1]].data(),
				vertex_pos[edge_vertex_index[(j << 1) + 1]].data(), mesh_struct_->edge_length[j], stiffness, sub_time_step, mass_inv[edge_vertex_index[j << 1]],
				mass_inv[edge_vertex_index[(j << 1) + 1]], *lambda_, damping_coe, initial_vertex_pos[edge_vertex_index[j << 1]].data(),
				initial_vertex_pos[edge_vertex_index[(j << 1)+1]].data());
			//std::cout << *lambda_ << " " << vertex_position[i][edge_vertex_index[j << 1]].data()[0] << " " << vertex_position[i][edge_vertex_index[(j << 1) + 1]].data()[0] << std::endl;
			//std::cout << vertex_pos[edge_vertex_index[j << 1]].data()[0]<<" "<< vertex_pos[edge_vertex_index[(j << 1) + 1]].data()[0]<<" "<< edge_vertex_index[j << 1] << " " << edge_vertex_index[(j << 1) + 1] << std::endl;		
			lambda_++;
			
		}
	}
}



void XPBD::solveTetStrainConstraint()
{
	unsigned int size;
	std::array<int, 4>* indices;
	MeshStruct* mesh_struct_;
	double* volume;
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	double* mass_inv;
	double stiffness;
	Matrix<double, 3, 4>* A;
	double* lambda_ = lambda.data() + constraint_index_start[2];
	double* sigma_limit;
	double youngs_modulus, poisson_ratio;
	for (unsigned int i = 0; i < tetrahedron->size(); ++i) {
		mesh_struct_ = mesh_struct[i+cloth->size()];
		size =tetrahedron->data()[i].mesh_struct.indices.size();
		indices =tetrahedron->data()[i].mesh_struct.indices.data();
		volume = tetrahedron->data()[i].mesh_struct.volume.data();
		vertex_pos = vertex_position[i+cloth->size()];
		initial_vertex_pos = initial_vertex_position[i + cloth->size()];
		stiffness = tetrahedron->data()[i].ARAP_stiffness;
		A = tetrahedron->data()[i].mesh_struct.A.data();
		sigma_limit = (*tetrahedron)[i].sigma_limit;
		mass_inv = mesh_struct_->mass_inv.data();
		youngs_modulus = tetrahedron->data()[i].youngs_modulus;
		poisson_ratio = tetrahedron->data()[i].poisson_ratio;


		//std::cout << "ARAP " << stiffness << " yong " << youngs_modulus << " " << "poiss " << poisson_ratio << std::endl;

		//for (unsigned int j = 0; j < mesh_struct_->vertex_position.size(); ++j) {
		//	std::cout << mass_inv[j] << std::endl;
		//}


		for (unsigned int j = 0; j < size; ++j) {
			XPBD_constraint.solveTetStrainConstraint(vertex_pos, initial_vertex_pos, stiffness, sub_time_step, A[j], indices[j].data(), mass_inv,
				*lambda_, damping_coe, volume[j], youngs_modulus, poisson_ratio);
			lambda_ ++;
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
	std::array<double, 3>* initial_vertex_pos;
	double stiffness;
	double* lambda_ = lambda.data();
	//double delta_t = sub_time_step;
	for (unsigned int i = 0; i < cloth->size(); ++i) {
		stiffness = cloth->data()[i].bend_stiffness;
		mesh_struct_ = mesh_struct[i];
		mass_inv = mesh_struct_->mass_inv.data();
		size = mesh_struct_->vertex_position.size();
		rest_mean_curvature_norm_ = rest_mean_curvature_norm[i].data();
		lbo_weight_ = lbo_weight[i].data();
		vertex_lbo_ = vertex_lbo[i].data();
		vertex_pos = vertex_position[i];
		initial_vertex_pos = initial_vertex_position[i];
		for (unsigned int j = 0; j < size; ++j) {
			XPBD_constraint.solveBendingConstraint(vertex_pos[j].data(), mass_inv[j], vertex_pos, mesh_struct_->vertices[j].neighbor_vertex,
				rest_mean_curvature_norm_[j], lbo_weight_[j], vertex_lbo_[j], stiffness, sub_time_step, mass_inv, *lambda_,damping_coe,
				initial_vertex_pos[j].data(),initial_vertex_pos);
			lambda_++;
		}		
	}
}

//XPBD_VELOCITY
void XPBD::computeVelocity(int thread_No)
{
	std::array<double, 3>* vertex_pos;
	std::array<double, 3>* initial_vertex_pos;
	unsigned int vertex_end = 0;
	std::array<double, 3>* velocity_;
	double delta_t = sub_time_step;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		initial_vertex_pos = initial_vertex_position[i];
		vertex_end = vertex_index_begin_per_thread[i][thread_No + 1];
		velocity_ = velocity[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < vertex_end; ++j) {
			for (unsigned int k = 0; k < 3; ++k) {
				velocity_[j][k] = (vertex_pos[j][k] - initial_vertex_pos[j][k]) / delta_t;
			}
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
	double gravity__[3];
	memcpy(gravity__, gravity, 24);
	double delta_t = sub_time_step;
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		vertex_pos = vertex_position[i];
		vertex_end = vertex_index_begin_per_thread[i][thread_No + 1];
		mass_inv = mesh_struct[i]->mass_inv.data();
		f_ext_ = f_ext[i].data();
		velocity_ = velocity[i].data();
		for (unsigned int j = vertex_index_begin_per_thread[i][thread_No]; j < vertex_end; ++j) {
			if (mass_inv[j] == 0) {
				for (unsigned int k = 0; k < 3; ++k) {
					vertex_pos[j][k] += delta_t * velocity_[j][k] + delta_t_2 * (mass_inv[j] * f_ext_[j][k]);
				}
			}
			else {
				for (unsigned int k = 0; k < 3; ++k) {
					vertex_pos[j][k] += delta_t * velocity_[j][k] + delta_t_2 * (mass_inv[j] * f_ext_[j][k] + gravity__[k]);
				}
			}		
		}
	}
}


void XPBD::resetExternalForce()
{
	for (unsigned int i = 0; i < total_obj_num; ++i) {
		memset(f_ext[i][0].data(), 0, 24 * f_ext[i].size());
	}
}

void XPBD::initialDHatTolerance(double ave_edge_length)
{
	collision.initialDHatTolerance(ave_edge_length);
}

void XPBD::updateTetrahedronAnchorVertices()
{
	double* mass_inv;
	int* anchor_vertex;
	unsigned int anchor_vertex_size;
	for (unsigned int i = cloth->size(); i < cloth->size()+tetrahedron->size(); ++i) {
		mass_inv = mesh_struct[i]->mass_inv.data();
		anchor_vertex_size = mesh_struct[i]->anchor_vertex.size();
		anchor_vertex = mesh_struct[i]->anchor_vertex.data();
		mesh_struct[i]->mass_inv = mesh_struct[i]->initial_mass_inv;
		for (unsigned int j = 0; j < anchor_vertex_size; ++j) {
			mass_inv[anchor_vertex[j]] = 0.0;
		}
	}
}

void XPBD::addExternalForce(double* neighbor_vertex_force_direction, std::vector<double>& coe, std::vector<int>& neighbor_vertex, int obj_No)
{
	if (!coe.empty()) {
		for (unsigned int i = 0; i < coe.size(); ++i) {
			for (unsigned int j = 0; j < 3; ++j) {
				f_ext[obj_No][neighbor_vertex[i]][j] += coe[i] * neighbor_vertex_force_direction[j];
			}
		}
	}
}