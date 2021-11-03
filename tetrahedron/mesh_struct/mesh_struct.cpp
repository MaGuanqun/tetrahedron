#include"mesh_struct.h"
void MeshStruct::initialNormalSize()
{
	face_normal_for_render.resize(triangle_indices.size(), {0.0,0.0,0.0});
	vertex_normal_for_render.resize(vertex_position.size());
	triangle_normal_magnitude_reciprocal.resize(triangle_indices.size());

}


void MeshStruct::setAnchorPosition()
{
	anchor_position.resize(anchor_vertex.size());
	for (int i = 0; i < anchor_vertex.size(); ++i) {
		anchor_position[i] = vertex_position[anchor_vertex[i]];
	}
}




