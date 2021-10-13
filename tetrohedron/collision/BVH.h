#pragma once
#include"../global_struct.h"


class BVH
{
public:
	void init_BVH_recursive(std::vector<AABB>& triangle_AABB, int node_index, int b, int e);

private:

};
