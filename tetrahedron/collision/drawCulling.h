#pragma once
#include"spatial_hashing.h"

class DrawCulling
{
public:
	void initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
		Thread* thread);//, std::vector<std::vector<std::vector<unsigned int>>>* triangle_patch,std::vector<std::vector<std::vector<unsigned int>>>* patch_vertex

	void draw(Camera* camera, float& far_plane);
	void move(unsigned int obj_No, double* displacement);
	void move(int thread_No, unsigned int obj_No);

	void moveDiffInitialCurrent(int thread_No, unsigned int obj_No);

	void setThreadDataTogether(int thread_No);
	void setAllTriangle(int thread_No);
	void drawAABBIntersectBetweenObjects();

	void setInSpatialHashingValue(unsigned int** spatial_hashing_value,
		unsigned int** spatial_hashing_triangle_index, unsigned int** spatial_hashing_value_collider, 
		unsigned int** spatial_hashing_triangle_index_collider, std::vector<std::vector<unsigned int>>* prefix_sum,
		std::vector<std::vector<unsigned int>>* prefix_sum_collider, std::vector<unsigned int>* actual_exist_cell_begin_per_thread);

	//type=1 for dinosaur, type=2 for dragon, type=3 cloth
	void moveScript(unsigned int type);
	size_t time_step = 0;


	void drawCell(Camera* camera, Shader* shader);
	void setCellData(std::vector<unsigned int>* hash, double cell_length, unsigned int* cell_number,
		double* initial_aabb);

	void setInSpatialHashingValue(std::vector<unsigned int>** spatial_hashing_cell, std::vector<unsigned int>** spatial_hashing_cell_collider,
		unsigned int hash_cell_count);

	void setSingleCellData(double cell_length, unsigned int* cell_number,
		double* initial_aabb);

	void setTetrahedronVertex();
	void drawTetTriangle(Camera* camera, Shader* object_shader_front, Light& light, float& far_plane);

	std::vector<unsigned int>* vertex_tet_pair;


private:
	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;

	unsigned int total_obj_num; //include collider
	unsigned int tetrahedron_end_index;
	unsigned int thread_num;

	double* displacement;

	unsigned int** spatial_hashing_value;
	//unsigned int** spatial_hashing_obj_index;
	unsigned int** spatial_hashing_triangle_index;

	unsigned int** spatial_hashing_value_collider;
	//unsigned int** spatial_hashing_obj_index_collider;
	unsigned int** spatial_hashing_triangle_index_collider;

	std::vector<std::vector<unsigned int>>* prefix_sum;
	std::vector<std::vector<unsigned int>>* prefix_sum_collider;
	std::vector<unsigned int>* actual_exist_cell_begin_per_thread;

	bool*** obj_is_used0;
	bool*** obj_is_used1;
	bool*** collider_is_used0;
	void initialUsedIndicator();

	std::vector<std::array<double, 3>> pos;
	std::vector<std::array<double, 3>> color;
	std::vector<std::array<double, 3>> normal;

	std::vector<std::vector<std::array<double, 3>>> pos_thread;
	std::vector<std::vector<std::array<double, 3>>> color_thread;
	std::vector<std::vector<std::array<double, 3>>> normal_thread;

	bool checkIfExistTwoObject(unsigned int cell_index);
	bool checkIfExistTwoObjectSP(unsigned int cell_index);

	std::vector<std::array<double, 3>*> vertex_position;
	std::vector<std::array<double, 3>*> vertex_normal_render;
	std::vector<std::array<int, 3>*> triangle_indices;

	std::vector<std::array<double, 3>*> collider_vertex_position;
	std::vector<std::array<double, 3>*> collider_vertex_normal_render;
	std::vector<std::array<int, 3>*> collider_triangle_indices;

	std::vector<std::array<double, 3>> palette;
	void setPalette();

	unsigned int VAO;
	unsigned int VBO[3];


	unsigned int VAO1;
	unsigned int VBO1;
	unsigned int EBO1;

	unsigned int VAO2;
	unsigned int VBO2;
	unsigned int EBO2;

	void genBuffer();
	Shader* shader;

	void setThreadDataTogether();
	std::vector<unsigned int> pos_start_thread;//the pos & color start index for every thread
	Light light;

	void recordDisplacement();
	std::vector<std::array<double, 3>> total_displacement;


	std::vector<std::array<double, 3>> displacement_test;
	std::vector<std::array<double, 3>> displacement_test_render;

	std::vector<unsigned int>* hash;

	
	void reverseHash(unsigned int hash_index, double cell_length, unsigned int* cell_number,
		double* initial_aabb, double* point);

	std::vector<double> point_value;
	std::vector<unsigned int> edge_index;

	unsigned int cell_num0_cell_num1;

	unsigned int cell_index_basic[24];
	void setCellIndexBasic();

	std::vector<unsigned int>** spatial_hashing_cell;//in every vector, we store triangle index,obj_index respectively
	std::vector<unsigned int>** spatial_hashing_cell_collider;//in every vector, we store triangle index, obj_index respectively
	unsigned int hash_cell_count;

	void findTheMaxNumberOfTriangleInCell();
	unsigned int hash_index;


	

	std::vector<unsigned int> tet_triangle_indices;

};
