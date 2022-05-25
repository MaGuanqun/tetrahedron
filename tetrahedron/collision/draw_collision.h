#pragma once

#include"../object/cloth.h"
#include"../object/collider.h"
#include"../object/tetrahedron.h"
#include"../thread.h"
#include"../external/shader.h"
#include"../basic/draw_vertex.h"


class DrawCollision
{
public:
	void initial(std::vector<Cloth>* cloth, std::vector<Collider>* collider, std::vector<Tetrahedron>* tetrahedron,
		Thread* thread);

	void setInPairInfo(std::vector<unsigned int>* point_triangle, std::vector<unsigned int>* point_triangle_collider, std::vector<unsigned int>* edge_edge);

	void setElementIndices();

	void drawCollision(bool draw_VT, Light& light, Camera* camera, Shader* object_shader_front, std::vector<std::vector<bool>>& show_collision_element);

	void setElementInOneCell(std::vector<std::vector<unsigned int>>& vertex_index, std::vector<std::vector<unsigned int>>& triangle_index,
		std::vector<std::vector<unsigned int>>& edge_index);

	void drawCollisionCell(bool draw_VT, Light& light, Camera* camera, Shader* object_shader_front, std::vector<std::vector<bool>>& show_collision_element,
		bool show_all_element);

private:

	std::vector<unsigned int>* point_triangle_target_pos_index;
	std::vector<unsigned int>* point_triangle_collider_target_pos_index;
	std::vector<unsigned int>* edge_edge_target_pos_index;

	DrawVertex draw_vertex;


	std::vector<Cloth>* cloth;
	std::vector<Tetrahedron>* tetrahedron;
	std::vector<Collider>* collider;
	Thread* thread;

	unsigned int total_obj_num; //include collider
	unsigned int tetrahedron_end_index;
	unsigned int thread_num;

	std::vector<unsigned int>VT_VAO;
	std::vector<unsigned int>VT_VBO;
	std::vector<unsigned int>VT_EBO;

	std::vector<unsigned int>EE_VAO;
	std::vector<unsigned int>EE_VBO;
	std::vector<unsigned int>EE_EBO;


	std::vector<unsigned int>VT_VAO_collide_in_a_cell;
	std::vector<unsigned int>VT_VBO_collide_in_a_cell;
	std::vector<unsigned int>VT_EBO_collide_in_a_cell;


	std::vector<unsigned int>EE_VAO_collide_in_a_cell;
	std::vector<unsigned int>EE_VBO_collide_in_a_cell;
	std::vector<unsigned int>EE_EBO_collide_in_a_cell;


	std::vector<unsigned int>VT_VAO_all_cell;
	std::vector<unsigned int>VT_VBO_all_cell;
	std::vector<unsigned int>VT_EBO_all_cell;

	std::vector<unsigned int>EE_VAO_all_cell;
	std::vector<unsigned int>EE_VBO_all_cell;
	std::vector<unsigned int>EE_EBO_all_cell;



	std::vector<std::vector<int>> triangle_vertex_index;
	std::vector<std::vector<int>> collider_triangle_vertex_index;
	std::vector<std::vector<unsigned int>> edge_vertex_index;

	std::vector<std::vector<unsigned int>> vertex_index;


	std::vector<std::vector<int>> triangle_vertex_index_in_a_cell;
	std::vector<std::vector<int>> collider_triangle_vertex_index_in_a_cell;
	std::vector<std::vector<unsigned int>> edge_vertex_index_in_a_cell;
	std::vector<std::vector<unsigned int>> vertex_index_in_a_cell;



	std::vector<std::vector<int>> triangle_vertex_index_all_cell;
	std::vector<std::vector<int>> collider_triangle_vertex_index_all_cell;
	std::vector<std::vector<unsigned int>> edge_vertex_index_all_cell;
	std::vector<std::vector<unsigned int>> vertex_index_all_cell;


	Shader* shader;
	Light light;

	void genBuffer();

	void setIndicesSize();
	void setTriangleIndices();
	void reorganzieDataOfObjects();
	std::vector<std::array<int, 3>*> triangle_indices;
	std::vector<std::array<int, 3>*> collider_triangle_indices;
	std::vector<unsigned int*>edge_indices;
	void setColliderTriangleIndices();
	void setEdgeIndices();
	void setVertexIndices();
	void setVertexTriangleBuffer(unsigned int obj_index);
	void setEdgeEdgeBuffer(unsigned int obj_index);
	std::vector<unsigned int>vertex_number;
	std::vector<std::array<double, 3>*> vertex_for_render;
	std::vector<std::array<double, 3>*> vertex_normal_for_render;
	void setBuffer();
	void setColliderBuffer(unsigned int obj_index);

	std::vector<std::vector<bool>> obj_is_used;

	void initialBoolean();
	void resetBooleanVector();

	void drawVertex(Camera* camera, std::vector<std::vector<bool>>& show_collision_element);

	void drawVertexCell(Camera* camera, std::vector<std::vector<bool>>& show_collision_element, bool show_all_element);


	void drawVT_triangle(Light& light, Camera* camera, Shader* object_shader_front, std::vector<std::vector<bool>>& show_collision_element,
		unsigned int* VT_VAO, std::vector<int>* triangle_vertex_index, std::vector<int>* collider_triangle_vertex_index);

	void drawEdge(Light& light, Camera* camera, Shader* object_shader_front, std::vector<std::vector<bool>>& show_collision_element,
		unsigned int* EE_VAO, std::vector<unsigned int>* edge_vertex_index);


	void setTriangleIndicesInOneCell(std::vector<std::vector<unsigned int>>& triangle_index);
	void setEdgeIndicesInOneCell(std::vector<std::vector<unsigned int>>& edge_index);
	void setVertexIndicesInOneCell(std::vector<std::vector<unsigned int>>& vertex_index_);
};
