#pragma once
#include"../basic/aabb.h"
#include<array>
#include"../thread.h"
#include"parallel_radix_sort.h"

class BVH
{
public:
	//struct SortStruct {
	//	int order;
	//	MortonCode64 morton;
	//};
	~BVH();

	void init(int triangle_num, std::vector<unsigned int>& triangle_index_begin_per_thread, Thread* thread);
	void buildBVH(std::vector<AABB>* triangle_AABB);
	void calCenterPerThread(int thread_No);
	void calMortonCode(int thread_No);
	void search(AABB& aabb, unsigned int compare_index, bool search_same_object, std::vector<unsigned int>* neighbor_list, unsigned int n, unsigned int b, unsigned int e);
	void updateBVH(std::vector<AABB>* aabb);
	void updateNodeValue(int thread_No);
	void updateNodeValueLastLayer(int thread_No);

	std::vector<AABB> aabb_list;


	std::vector<unsigned int> triangle_node_index; //triangle -> node index

private:

	std::vector<std::array<double, 3>> aabb_center;
	//std::vector<SortStruct> list;

	std::vector<uint64_t> morton_list;
	std::vector<unsigned int>new2old;

	std::vector<unsigned int>triangle_index_begin_per_thread;
	std::vector<std::array<double, 3>> aabb_min_per_thread;
	std::vector<std::array<double, 3>> aabb_max_per_thread;
	void setMortonCode();
	
	int total_thread_num;
	std::vector<AABB>* triangle_AABB;
	//double center[3];
	double min[3];
	const int multi = 10000;
	Thread* thread;	
	std::vector<unsigned int> old2new;
	int maxNodeIndex(int node_index, int b, int e);
	
	void initBVHRecursive(std::vector<AABB>* triangle_AABB, int node_index, int b, int e);
	RadixSort radix_sort;

	uint64_t splitBy3Bits21(uint32_t x);
	uint64_t mortonCode64(uint32_t x, uint32_t y, uint32_t z);
	uint64_t calMaxMortonCode(double* max, double* min);

	void test();

	unsigned int assumed_thread_num;
	std::vector<unsigned int>start_leaf_node_per_thread;
	unsigned int triangle_num;
	unsigned int last_layer_start_index;
	void setLeafNode();
	unsigned int total_layers;
	unsigned int final_multi_thread_layers;

	void obtainAABB(AABB& result, AABB& left, AABB& right);
	void recursiveUpdate(unsigned int start_index, unsigned int end_index);

	
	bool* is_node_leaf;

	void recordNodeInfo(int node_index, int b, int e);
};
