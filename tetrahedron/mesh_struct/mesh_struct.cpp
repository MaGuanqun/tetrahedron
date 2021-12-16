#include"mesh_struct.h"
void MeshStruct::initialNormalSize()
{
	face_normal_for_render.resize(triangle_indices.size(), {0.0,0.0,0.0});
	vertex_normal_for_render.resize(vertex_position.size());
	triangle_normal_magnitude_reciprocal.resize(triangle_indices.size());
	ori_face_normal_for_render.resize(triangle_indices.size(), { 0.0,0.0,0.0 });
	cross_for_approx_CCD.resize(triangle_indices.size(), { 0.0,0.0,0.0 });
	face_normal = face_normal_for_render;
	ori_face_normal = ori_face_normal_for_render;
	vertex_normal = vertex_normal_for_render;
	
	f_face_normal_for_render.resize(triangle_indices.size());
	f_face_normal.resize(triangle_indices.size());
	f_cross_for_approx_CCD.resize(triangle_indices.size());
}


void MeshStruct::setAnchorPosition()
{
	anchor_position.resize(anchor_vertex.size());
	for (int i = 0; i < anchor_vertex.size(); ++i) {
		anchor_position[i] = vertex_position[anchor_vertex[i]];
	}
}




