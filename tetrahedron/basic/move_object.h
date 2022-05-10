#pragma once
#include"../object/cloth.h"
#include"../object/collider.h"
#include"../object/tetrahedron.h"
#include"../thread.h"

class MoveObject
{
public:
	void initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
		Thread* thread);

	void move(unsigned int obj_No, double* displacement, bool only_move_vertex_pos);
	void move(int thread_No, unsigned int obj_No);
	void moveDiffInitialCurrent(int thread_No, unsigned int obj_No);


	//type=1 for dinosaur, type=2 for dragon, type=3 cloth
	void moveScript(unsigned int type);


	int select_object_index;

	void rotation(double angle_move, unsigned int dimension, bool only_move_vertex_pos);
	void rotateAroundAxis(int thread_No, unsigned int obj_No);

	unsigned int select_dimension;

private:
	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;

	unsigned int total_obj_num; //include collider
	unsigned int tetrahedron_end_index;
	unsigned int thread_num;
	double* displacement;

	std::vector<std::array<double, 3>> total_displacement;

	std::vector<std::array<double, 3>> displacement_test;
	std::vector<std::array<double, 3>> displacement_test_render;



	double angle_move;

	double center[3];

};
