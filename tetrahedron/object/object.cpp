#include"object.h"

void Object::genBuffer()
{
	glGenVertexArrays(1, &VAO);
	glGenBuffers(3, VBO);
	glGenBuffers(1, &EBO);
}


void Object::setWireframwColor(double* color)
{
	wireframe_color = glm::vec3(color[0], color[1], color[2]);
}

void Object::getAABB(AABB& target, AABB& aabb0, AABB& aabb1, AABB& aabb2)
{
	for (int i = 0; i < 3; ++i) {
		target.min[i] = myMin(aabb0.min[i], aabb1.min[i]);
		target.max[i] = myMax(aabb0.max[i], aabb1.max[i]);
		if (target.min[i] > aabb2.min[i]) {
			target.min[i] = aabb2.min[i];
		}
		if (target.max[i] < aabb2.max[i]) {
			target.max[i] = aabb2.max[i];
		}
	}	
}

void Object::getAABB(AABB& target, AABB& aabb0, AABB& aabb1)
{
	for (int i = 0; i < 3; ++i) {
		target.min[i] = myMin(aabb0.min[i], aabb1.min[i]);
		target.max[i] = myMax(aabb0.max[i], aabb1.max[i]);
	}
}

void Object::getAABB(AABB& target, AABB& aabb0, AABB& aabb1, double radius)
{
	for (int i = 0; i < 3; ++i) {
		target.min[i] = myMin(aabb0.min[i], aabb1.min[i])-radius;
		target.max[i] = myMax(aabb0.max[i], aabb1.max[i])+radius;
	}
}

void Object::setOrder(bool* in_this_triangle, int count, int* index)
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

void Object::setOrderEdge(bool* in_this_triangle, int count, int* index)
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




void Object::setRepresentativeEdge(std::vector<MeshStruct::Face>& face, std::vector<MeshStruct::Edge>& edge)
{
	int count;
	bool in_this_triangle[3];
	std::vector<bool> is_used(edge.size(), false);
	for (int i = 0; i < face.size(); ++i) {
		count = 0;
		memset(in_this_triangle, 0, 3);
		for (int j = 0; j < 3; ++j) {
			if (!is_used[face[i].edge[j]]) {
				count++;
				is_used[face[i].edge[j]] = true;
				in_this_triangle[j] = true;
			}
		}
		representative_edge_num[i] = count;
		setOrderEdge(in_this_triangle, count, face[i].edge.data());
	}
}


void Object::setRepresentativeVertex(std::vector<std::array<int,3>>& face, std::vector<MeshStruct::Vertex>& vertex)
{
	int count;
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