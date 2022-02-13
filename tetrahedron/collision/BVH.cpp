#include"BVH.h"

//const auto morton_compare(const BVH::SortStruct& a, const BVH::SortStruct& b) {
//	return (a.morton < b.morton);
//};
BVH::~BVH()
{
	delete[] is_node_leaf;
}


void BVH::init(int triangle_num, std::vector<unsigned int>&triangle_index_begin_per_thread, Thread* thread)
{
	this->triangle_num = triangle_num;
	aabb_center.resize(triangle_num);
	//list.resize(triangle_num);
	morton_list.resize(triangle_num);
	this->triangle_index_begin_per_thread = triangle_index_begin_per_thread;
	total_thread_num = std::thread::hardware_concurrency();
	aabb_min_per_thread.resize(total_thread_num);
	aabb_max_per_thread.resize(total_thread_num);
	this->thread = thread;
	new2old.resize(triangle_num);
	old2new.resize(triangle_num);	

	aabb_list.resize(maxNodeIndex(1, 0, triangle_num) + 1);
	is_node_leaf = new bool[aabb_list.size()];
	memset(is_node_leaf, 0,sizeof(bool)*aabb_list.size());
	triangle_node_index.resize(triangle_num);
	recordNodeInfo(1, 0, triangle_num);

	radix_sort.initial(thread);
	radix_sort.initialMortonArray(triangle_num);
	

	//here we assume the thread_num is 2^x.
	
	//std::cout << "BVH construct " << aabb_list.size()<<" "<< new2old.size() << std::endl;
	setLeafNode();
	//unsigned int a = 7 >> 1;
	//unsigned int b = 6 >> 1;
	//a <<= 1;
	//b <<= 1;
	//std::cout << "test bit operator" << a << " " << b << std::endl;
}


void BVH::setLeafNode()
{
	final_multi_thread_layers = floor(log2(total_thread_num));
	assumed_thread_num =pow(2, final_multi_thread_layers);
	total_layers = ceil(log2(triangle_num));
	
	start_leaf_node_per_thread.resize(assumed_thread_num + 1);
	last_layer_start_index = pow(2, (int)ceil(log2(triangle_num)));
	arrangeIndex(total_thread_num, last_layer_start_index, start_leaf_node_per_thread.data());
	if (assumed_thread_num < total_thread_num) {
		start_leaf_node_per_thread.resize(total_thread_num + 1);
		for (int i = assumed_thread_num + 1; i < start_leaf_node_per_thread.size(); ++i) {
			start_leaf_node_per_thread[i] = start_leaf_node_per_thread[assumed_thread_num];
		}
	}
	
	for (int i = 0; i < start_leaf_node_per_thread.size(); ++i) {
		start_leaf_node_per_thread[i] += last_layer_start_index;
	}
	
}

void BVH::buildBVH(std::vector<AABB>* triangle_AABB)
{
	this->triangle_AABB = triangle_AABB;
	setMortonCode();
	initBVHRecursive(triangle_AABB, 1, 0, triangle_AABB->size());
}


void BVH::updateBVH(std::vector<AABB>* aabb)
{
	this->triangle_AABB = aabb;
	setMortonCode();

	if (triangle_AABB->size() <= assumed_thread_num >> 1) {
		initBVHRecursive(triangle_AABB, 1, 0, triangle_AABB->size());
	}
	else {
		thread->assignTask(this, UPDATE_LAST_LAYER_NODE_VALUE);
		if (total_layers >= final_multi_thread_layers) {
			thread->assignTask(this, UPDATE_NODE_VALUE);
		}
		unsigned int start, end;
		start = assumed_thread_num >> 1;
		end = assumed_thread_num;
		for (int j = final_multi_thread_layers - 1; j >= 0; --j) {
			for (int i = start; i < end; ++i) {
				obtainAABB(aabb_list[i], aabb_list[2 * i], aabb_list[2 * i + 1]);
			}
			start >>= 1;
			end >>= 1;
		}
	}
	
}


//UPDATE_LAST_LAYER_NODE_VALUE
void BVH::updateNodeValueLastLayer(int thread_No)
{
	for (unsigned int i = triangle_index_begin_per_thread[thread_No]; i < triangle_index_begin_per_thread[thread_No + 1]; ++i) {
		aabb_list[triangle_node_index[i]] = triangle_AABB->data()[new2old[i]];
	}
}

//UPDATE_NODE_VALUE
void BVH::updateNodeValue(int thread_No)
{	
	
	unsigned int end_index = start_leaf_node_per_thread[thread_No + 1] >>1;
	for (unsigned int i = start_leaf_node_per_thread[thread_No] >>1; i < end_index; ++i) {
		if (!is_node_leaf[i]) {				
			obtainAABB(aabb_list[i], aabb_list[2 * i], aabb_list[2 * i + 1]);
		}
	}
	//recursiveUpdate(start_leaf_node_per_thread[thread_No] / 4, start_leaf_node_per_thread[thread_No + 1] / 4);
	unsigned int start_index = start_leaf_node_per_thread[thread_No] >> 2;
	end_index = start_leaf_node_per_thread[thread_No + 1] >> 2;
	while (true)
	{
		if (start_index == end_index) {
			break;
		}
		for (int i = start_index; i < end_index; ++i) {
			obtainAABB(aabb_list[i], aabb_list[2 * i], aabb_list[2 * i + 1]);
		}
		start_index >>= 1;
		end_index >>= 1;
	}
	

}



void BVH::obtainAABB(AABB& result, AABB& left, AABB& right)
{
	//for (unsigned int i = 0; i < 3; ++i) {
	//	result.min[i] = std::min(left.min[i], right.min[i]);
	//	result.max[i] = std::min(left.max[i], right.max[i]);
	//}
	result = left;
	for (unsigned int i = 0; i < 3; ++i) {
		if (result.min[i] > right.min[i]) {
			result.min[i] = right.min[i];
		}
		if (result.max[i] < right.max[i]) {
			result.max[i] = right.max[i];
		}
	}
}

void BVH::recursiveUpdate(unsigned int start_index, unsigned int end_index)
{
	if (start_index == end_index) {
		return;
	}
	for (int i = start_index; i < end_index; ++i) {
		obtainAABB(aabb_list[i], aabb_list[2 * i], aabb_list[2 * i + 1]);
	}

	recursiveUpdate(start_index / 2, end_index / 2);
}


void BVH::search(AABB& aabb, unsigned int compare_index, bool search_same_object, std::vector<unsigned int>* neighbor_list, unsigned int n, unsigned int b, unsigned int e)
{
	if (!aabb.AABB_intersection(aabb_list[n])) {
		return;
	}
	if (e == b + 1) {
		if(search_same_object){
			if (compare_index != new2old[b]) {
				neighbor_list->push_back(new2old[b]);
			}
		}
		else {
			neighbor_list->push_back(new2old[b]);
		}
		return;
	}
	unsigned int m = b + (e - b) / 2;
	unsigned int child_left = 2 * n;
	unsigned int child_right = 2 * n + 1;

	search(aabb, compare_index, search_same_object, neighbor_list, child_left, b, m);
	search(aabb, compare_index, search_same_object, neighbor_list, child_right, m, e);
}


//CAL_CENTER
void BVH::calCenterPerThread(int thread_No)
{
	aabb_min_per_thread[thread_No] = std::array{ DBL_MAX,DBL_MAX ,DBL_MAX };
	aabb_max_per_thread[thread_No] = std::array{ DBL_MIN,DBL_MIN ,DBL_MIN };
	double* center; double* max; double* min; double* aabb_min; double* aabb_max;
	aabb_min = aabb_min_per_thread[thread_No].data();
	aabb_max = aabb_max_per_thread[thread_No].data();
	for (int i = triangle_index_begin_per_thread[thread_No]; i < triangle_index_begin_per_thread[thread_No + 1]; ++i) {
		center = aabb_center[i].data();
		min = triangle_AABB->data()[i].min;
		max = triangle_AABB->data()[i].max;		
		TWO_POINTS_CENTER(center, min, max);
		for (int j = 0; j < 3; ++j) {
			if (aabb_min[j] > min[j]) {
				aabb_min[j] = min[j];
			}
			if (aabb_max[j] < center[j]) {
				aabb_max[j] = center[j];
			}
		}		
	}
}

//CAL_MORTON
void BVH::calMortonCode(int thread_No)
{
	double* center_;
	for (unsigned int i = triangle_index_begin_per_thread[thread_No]; i < triangle_index_begin_per_thread[thread_No + 1]; ++i) {
		center_ = aabb_center[i].data();
		SUB(center_, center_, min);
		MULTI(center_, center_, multi);
		//list[i].morton = MortonCode64(unsigned int(center_[0]), unsigned int(center_[1]), unsigned int(center_[2]));
		//list[i].order = i;
		morton_list[i] = mortonCode64(unsigned int(center_[0]), unsigned int(center_[1]), unsigned int(center_[2]));
		new2old[i] = i;
	}
}

void BVH::setMortonCode()
{
	thread->assignTask(this, CAL_CENTER);

	memcpy(min, aabb_min_per_thread[0].data(), 24);
	memcpy(v_max, aabb_max_per_thread[0].data(), 24);
	for (int i = 1; i < total_thread_num; ++i) {
		for (int j = 0; j < 3; ++j) {
			min[j] = myMin(min[j], aabb_min_per_thread[i][j]);
			v_max[j] = myMax(v_max[j], aabb_max_per_thread[i][j]);
		}
	}
	for (int i = 0; i < 3; ++i) {
		min[i] -= 1e-7;
	}
	//TWO_POINTS_CENTER(center, min, v_max);
	thread->assignTask(this, CAL_MORTON);
	//std::sort(list.begin(), list.end(), morton_compare);
	radix_sort.radixSort(calMaxMortonCode(v_max, min), &morton_list, new2old.data());
	//for (int i = 1; i < morton_list.size(); ++i) {
	//	//std::cout << morton_list[i] << std::endl;
	//	if (morton_list[i] < morton_list[i - 1]) {
	//		std::cout << "error morton code" << std::endl;
	//	}
	//}
	for (unsigned int i = 0; i < new2old.size(); ++i) {
		old2new[new2old[i]] = i;
	}
}


void BVH::initBVHRecursive(std::vector<AABB>* triangle_AABB, int node_index, int b, int e)
{
	if (b + 1 == e) {
		aabb_list[node_index] = triangle_AABB->data()[new2old[b]];
		return;
	}
	int m = b + (e - b) / 2;
	int child_left = 2 * node_index;
	int child_right = 2 * node_index + 1;

	initBVHRecursive(triangle_AABB, child_left, b, m);
	initBVHRecursive(triangle_AABB, child_right, m, e);
	for (int i = 0; i < 3; ++i) {
		aabb_list[node_index].min[i] = std::min(aabb_list[child_left].min[i], aabb_list[child_right].min[i]);
		aabb_list[node_index].max[i] = std::max(aabb_list[child_left].max[i], aabb_list[child_right].max[i]);
	}
}


int BVH::maxNodeIndex(int node_index, int b, int e)
{
	if (b + 1 == e) {
		return node_index;
	}


	int m = b + (e - b) / 2;
	int child_left = 2 * node_index;
	int child_right = 2 * node_index + 1;
	return std::max(maxNodeIndex(child_left, b, m), maxNodeIndex(child_right, m, e));
}

void BVH::recordNodeInfo(int node_index, int b, int e)
{
	if (b + 1 == e) {
		is_node_leaf[node_index] = true;
		triangle_node_index[b] = node_index;
		return;
	}


	int m = b + (e - b) / 2;
	int child_left = 2 * node_index;
	int child_right = 2 * node_index + 1;
	recordNodeInfo(child_left, b, m);
	recordNodeInfo(child_right, m, e);
}

uint64_t BVH::splitBy3Bits21(uint32_t x)
{
	uint64_t r = x & 0x1fffff;//only use first 21
	//uint64_t r = x;

	r = (r | r << 32) & 0x1f00000000ffff;//0000000000011111000000000000000000000000000000001111111111111111
	r = (r | r << 16) & 0x1f0000ff0000ff;//0000000000011111000000000000000011111111000000000000000011111111
	r = (r | r << 8) & 0x100f00f00f00f00f;//0001000000001111000000001111000000001111000000001111000000001111
	r = (r | r << 4) & 0x10c30c30c30c30c3;//0001000011000011000011000011000011000011000011000011000011000011
	r = (r | r << 2) & 0x1249249249249249;//0001001001001001001001001001001001001001001001001001001001001001
	return r;
}

uint64_t BVH::mortonCode64(uint32_t x, uint32_t y, uint32_t z)
{
	return (splitBy3Bits21(x) | splitBy3Bits21(y) << 1 | splitBy3Bits21(z) << 2);
	//data |= 0x7000000000000000;
}


uint64_t BVH::calMaxMortonCode(double* max, double* min)
{
	double max_[3];
	SUB(max_, max, min);
	MULTI(max_, max_, multi);
	return mortonCode64(unsigned int(max_[0]), unsigned int(max_[1]), unsigned int(max_[2]));
}


void BVH::test(std::vector<AABB>* aabb)
{
	//std::vector<SortStruct> temp_list;
//std::vector<uint64_t> temp_morton_list;
//std::vector<int>temp_triangle_list;
//time_t t = clock();
//for (int i = 0; i < 1000; ++i) {
//	//temp_morton_list = morton_list;
//	//temp_triangle_list = triangle_list;
//	temp_list = list;
//	//radix_sort.radixSort(calMaxMortonCode(v_max, min), &temp_morton_list, &temp_triangle_list);
//	std::sort(temp_list.begin(), temp_list.end(), morton_compare);
//}
//time_t t2 = clock() - t;
//int k = 0;
//t = clock();
//for (int i = 0; i < 1000; ++i) {
///*	temp_morton_list = morton_list;
//	temp_triangle_list = triangle_list;*/
//	temp_list = list;
//	k++;
//}
//time_t t3 = clock() - t;
//std::cout << temp_morton_list.size() << " " << t2 << " " << t3 << " " << t2 - t3 << std::endl;

	//for (int i = 1; i < morton_list.size(); ++i) {
	//	if (morton_list[i] - morton_list[i - 1] < 0) {
	//		std::cout << "error " << std::endl;
	//	}
	//}
	std::vector<AABB>test_aabb = *triangle_AABB;
	this->triangle_AABB = triangle_AABB;
	int itr_num = 1000;
	int k = 0;
	time_t t = clock();
	for (int i = 0; i < itr_num; ++i) {
		setMortonCode();
		k++;
	}
	std::cout << "morton code1 1 " << clock() - t << std::endl;
	std::cout << k << std::endl;
	k = 0;
	t = clock();
	for (int i = 0; i < itr_num; ++i) {
		setMortonCode();
		test_buildBVH(aabb);
		k++;
	}
	std::cout << "total update 1 " << clock() - t << std::endl;
	std::cout << k << std::endl;
	k = 0;
	t = clock();
	for (int i = 0; i < itr_num; ++i) {
		setMortonCode();
		test_updateBVH(aabb);
		k++;
	}
	std::cout << "total update 1 " << clock() - t << std::endl;
	std::cout << k << std::endl;
	std::cout << "size(triangle num) " << triangle_num << std::endl;
	k = 0;
	t = clock();
	for (int i = 0; i < itr_num; ++i) {
		//setMortonCode();
		//setMortonCode();
		//setMortonCode2();
		//setMortonCode3();
		test_buildBVH(aabb);
		k++;
	}
	std::cout << "single thread " << clock() - t << std::endl;
	std::cout << k << std::endl;
	k = 0;
	t = clock();
	for (int i = 0; i < itr_num; ++i) {
		setMortonCode();
		k++;
	}
	std::cout << "morton code 2 " << clock() - t << std::endl;
	std::cout << k << std::endl;
	k = 0;
	t = clock();
	for (int i = 0; i < itr_num; ++i) {
		//setMortonCode();
		//setMortonCode2();
		//setMortonCode3();
		test_updateBVH(aabb);
		k++;
	}
	std::cout << "multi thread " << clock() - t << std::endl;
	std::cout << k << std::endl;
	std::cout <<"assume "<< (assumed_thread_num >> 1) << " " << assumed_thread_num << std::endl;

}

void BVH::test_buildBVH(std::vector<AABB>* triangle_AABB)
{
	initBVHRecursive(triangle_AABB, 1, 0, triangle_AABB->size());
}
void BVH::test_updateBVH(std::vector<AABB>* triangle_AABB)
{
	if (triangle_AABB->size() <= (assumed_thread_num >> 1)) {
		initBVHRecursive(triangle_AABB, 1, 0, triangle_AABB->size());
	}
	else {
		thread->assignTask(this, UPDATE_LAST_LAYER_NODE_VALUE);
		if (total_layers >= final_multi_thread_layers) {
			thread->assignTask(this, UPDATE_NODE_VALUE);
		}
		unsigned int start, end;
		start = assumed_thread_num >> 1;
		end = assumed_thread_num;
		for (int j = final_multi_thread_layers - 1; j >= 0; --j) {
			for (int i = start; i < end; ++i) {
				obtainAABB(aabb_list[i], aabb_list[2 * i], aabb_list[2 * i + 1]);
			}
			start >>= 1;
			end >>= 1;
		}
	}
}