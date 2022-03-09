#include"object.h"

void Object::genBuffer()
{
	glGenVertexArrays(1, &VAO);
	glGenBuffers(3, VBO);
	glGenBuffers(1, &EBO);
}


void Object::combineObjAABB()
{
	memcpy(obj_aabb, obj_aabb_per_thread[0].data(), 48);
	for (unsigned int i = 1; i < total_thread_num; ++i) {
		for (unsigned int j = 0; j < 3; ++j) {
			if (obj_aabb[j] > obj_aabb_per_thread[i][j]) {
				obj_aabb[j] = obj_aabb_per_thread[i][j];
			}
		}
		for (unsigned int j = 3; j < 6; ++j) {
			if (obj_aabb[j] < obj_aabb_per_thread[i][j]) {
				obj_aabb[j] = obj_aabb_per_thread[i][j];
			}
		}
	}
}

void Object::setWireframwColor(double* color)
{
	wireframe_color = glm::vec3(color[0], color[1], color[2]);
}

void Object::getAABB(double* target, double* aabb0, double* aabb1, double* aabb2)
{
	for (int i = 0; i < 3; ++i) {
		target[i] = myMin(aabb0[i], aabb1[i]);
		if (target[i] > aabb2[i]) {
			target[i] = aabb2[i];
		}
	}
	for (int i = 3; i < 6; ++i) {
		target[i] = myMax(aabb0[i], aabb1[i]);
		if (target[i] < aabb2[i]) {
			target[i] = aabb2[i];
		}
	}
}


void Object::getAABB(double* target, double* aabb0, double* aabb1, double* aabb2, double* aabb3)
{
	for (int i = 0; i < 3; ++i) {
		target[i] = myMin(aabb0[i], aabb1[i]);
		if (target[i] > aabb2[i]) {
			target[i] = aabb2[i];
		}
		if (target[i] > aabb3[i]) {
			target[i] = aabb3[i];
		}
	}
	for (int i = 3; i < 6; ++i) {
		target[i] = myMax(aabb0[i], aabb1[i]);
		if (target[i] < aabb2[i]) {
			target[i] = aabb2[i];
		}
		if (target[i] < aabb3[i]) {
			target[i] = aabb3[i];
		}
	}
}

void Object::getAABB(double* target, double* aabb0, double* aabb1)
{
	for (int i = 0; i < 3; ++i) {
		target[i] = myMin(aabb0[i], aabb1[i]);
	}
	for (int i = 3; i < 6; ++i) {
		target[i] = myMax(aabb0[i], aabb1[i]);
	}
}

void Object::getAABB(double* target, double* aabb0, double* aabb1, double radius)
{
	for (int i = 0; i < 3; ++i) {
		target[i] = myMin(aabb0[i], aabb1[i]) - radius;
	}
	for (int i = 3; i < 6; ++i) {
		target[i] = myMax(aabb0[i], aabb1[i]) + radius;
	}
}

void Object::setOrder(bool* in_this_triangle, unsigned int count, int* index)
{
	if (count == 1) {
		if (in_this_triangle[1]) {
			int temp = index[0];
			index[0] = index[1];
			index[1] = index[2];
			index[2] = temp;
		}
		else if (in_this_triangle[2]) {
			int temp = index[0];
			index[0] = index[2];
			index[2] = index[1];
			index[1] = temp;
		}
	}
	else if (count == 2) {
		if (!in_this_triangle[0]) {
			int temp = index[0];
			index[0] = index[1];
			index[1] = index[2];
			index[2] = temp;
		}
		else if (!in_this_triangle[1]) {
			int temp = index[0];
			index[0] = index[2];
			index[2] = index[1];
			index[1] = temp;
		}
	}
}

void Object::setOrderEdge(bool* in_this_triangle, unsigned int count, int* index)
{
	if (count == 1) {
		if (in_this_triangle[1]) {
			int temp = index[0];
			index[0] = index[1];
			index[1] = temp;
		}
		else if (in_this_triangle[2]) {
			int temp = index[0];
			index[0] = index[2];
			index[2] = temp;
		}
	}
	else if (count == 2) {
		if (!in_this_triangle[0]) {
			int temp = index[0];
			index[0] = index[2];
			index[2] = temp;
		}
		else if (!in_this_triangle[1]) {
			int temp = index[2];
			index[2] = index[1];
			index[1] = temp;
		}
	}
}




void Object::setRepresentativeEdge(int* face_edges, int face_size, int edge_size)
{
	unsigned int count;
	bool in_this_triangle[3];
	std::vector<bool> is_used(edge_size, false);
	for (int i = 0; i < face_size; ++i) {
		count = 0;
		memset(in_this_triangle, 0, 3);
		for (int j = 0; j < 3; ++j) {
			if (!is_used[face_edges[3 * i + j]]) {
				count++;
				is_used[face_edges[3 * i + j]] = true;
				in_this_triangle[j] = true;
			}
		}
		representative_edge_num[i] = count;
		setOrderEdge(in_this_triangle, count, face_edges + 3 * i);
	}
}


void Object::setRepresentativeVertex(std::vector<std::array<int, 3>>& face, std::vector<MeshStruct::Vertex>& vertex)
{
	unsigned int count;
	bool in_this_triangle[3];
	std::vector<bool> is_used(vertex.size(), false);
	for (int i = 0; i < face.size(); ++i) {
		count = 0;
		memset(in_this_triangle, 0, 3);
		for (int j = 0; j < 3; ++j) {
			if (!is_used[face[i][j]]) {
				count++;
				is_used[face[i][j]] = true;
				in_this_triangle[j] = true;
			}
		}
		representative_vertex_num[i] = count;
		setOrder(in_this_triangle, count, face[i].data());
	}
}