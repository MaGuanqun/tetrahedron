#include"mesh_struct.h"
void MeshStruct::initialNormalSize()
{
	face_norm_for_render.resize(triangle_indices.size() / 3, {0.0,0.0,0.0});
	temp_face_norm.resize(triangle_indices.size() / 3, { 0.0,0.0,0.0 });
	vertex_norm_for_render.resize(vertex_position.size());
}

void MeshStruct::arrangeIndex(int total_thread_num, int total_num, std::vector<int>& begin) {
	int interval1 = total_num / total_thread_num;
	int resi1 = total_num % total_thread_num;
	if (resi1 == 0) {
		for (int i = 0; i < total_thread_num; ++i) {
			begin[i] = interval1 * i;
		}
	}
	else {
		for (int i = 0; i < total_thread_num; ++i) {
			if (i < resi1) {
				begin[i] = (interval1 + 1) * i;
			}
			else {
				begin[i] = (interval1 + 1) * resi1 + interval1 * (i - resi1);
			}
		}		
	}
	begin[total_thread_num] = total_num;
}

void MeshStruct::setAnchorPosition()
{
	for (int i = 0; i < anchor_vertex.size(); ++i) {
		anchor_position.push_back(vertex_position[anchor_vertex[i]]);

	}

}




