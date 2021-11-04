#include"BVH.h"

const auto morton_compare(const BVH::SortStruct& a, const BVH::SortStruct& b) {
	return (a.morton < b.morton);
};

void BVH::init(int triangle_num, std::vector<int>&triangle_index_begin_per_thread, Thread* thread)
{
	aabb_center.resize(triangle_num);
	list.resize(triangle_num);
	this->triangle_index_begin_per_thread = triangle_index_begin_per_thread;
	total_thread_num = std::thread::hardware_concurrency();
	aabb_min_per_thread.resize(total_thread_num);
	aabb_max_per_thread.resize(total_thread_num);
	this->thread = thread;
	new2old.resize(triangle_num);
	aabb_list.resize(maxNodeIndex(1, 0, triangle_num) + 1);
}

void BVH::buildBVH(std::vector<AABB>* triangle_AABB)
{
	this->triangle_AABB = triangle_AABB;
	setMortonCode();
	initBVHRecursive(triangle_AABB, 1, 0, triangle_AABB->size());
}

void BVH::search(AABB& aabb, int compare_index, bool search_same_object, std::vector<int>* neighbor_list, int n, int b, int e)
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
	int m = b + (e - b) / 2;
	int child_left = 2 * n;
	int child_right = 2 * n + 1;

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
		min = (*triangle_AABB)[i].min;
		max = (*triangle_AABB)[i].max;		
		TWO_POINTS_CENTER(center, min, max);
		for (int j = 0; j < 3; ++j) {
			if (aabb_min[j] > center[j]) {
				aabb_min[j] = center[j];
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
	for (int i = triangle_index_begin_per_thread[thread_No]; i < triangle_index_begin_per_thread[thread_No + 1]; ++i) {
		center_ = aabb_center[i].data();
		SUB(center_, center_, center);
		MULTI(center_, center_, multi);
		list[i].morton = MortonCode64(int(center_[0]), int(center_[1]), int(center_[2]));
		list[i].order = i;
	}
}

void BVH::setMortonCode()
{
	thread->assignTask(this, CAL_CENTER);
	double v_max[3]; double v_min[3];
	memcpy(v_min, aabb_min_per_thread[0].data(), 24);
	memcpy(v_max, aabb_max_per_thread[0].data(), 24);
	for (int i = 1; i < total_thread_num; ++i) {
		for (int j = 0; j < 3; ++j) {
			v_min[j] = myMin(v_min[j], aabb_min_per_thread[i][j]);
			v_max[j] = myMax(v_max[j], aabb_max_per_thread[i][j]);
		}
	}
	TWO_POINTS_CENTER(center, v_min, v_max);
	thread->assignTask(this, CAL_MORTON);
	std::sort(list.begin(), list.end(), morton_compare);
	for (int i = 0; i < new2old.size(); ++i) {
		new2old[i] = list[i].order;
	}
}


void BVH::initBVHRecursive(std::vector<AABB>* triangle_AABB, int node_index, int b, int e)
{
	if (b + 1 == e) {
		aabb_list[node_index] = (*triangle_AABB)[new2old[b]];
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
