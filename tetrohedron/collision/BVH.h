#pragma once
#include"../basic/aabb.h"
#include"morton.h"
#include<array>
#include"../thread.h"

class BVH
{
public:
	struct SortStruct {
		int order;
		MortonCode64 morton;
	};

	void init(int triangle_num, std::vector<int>& triangle_index_begin_per_thread, Thread* thread);
	void buildBVH(std::vector<AABB>* triangle_AABB);
	void calCenterPerThread(int thread_No);
	void calMortonCode(int thread_No);
	void search(AABB& aabb, int compare_index, bool search_same_object, std::vector<int>* neighbor_list, int n, int b, int e);private:
	
	std::vector<std::array<double, 3>> aabb_center;
	std::vector<SortStruct> list;
	std::vector<int>triangle_index_begin_per_thread;
	std::vector<std::array<double, 3>> aabb_min_per_thread;
	std::vector<std::array<double, 3>> aabb_max_per_thread;
	void setMortonCode();
	
	int total_thread_num;
	std::vector<AABB>* triangle_AABB;
	double center[3];
	const int multi = 10000;
	Thread* thread;
	std::vector<int> new2old;
	int maxNodeIndex(int node_index, int b, int e);
	std::vector<AABB> aabb_list;
	void initBVHRecursive(std::vector<AABB>* triangle_AABB, int node_index, int b, int e);
};
