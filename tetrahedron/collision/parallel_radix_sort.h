//we use key with 4 bit for each bucket

#pragma once
#include"../thread.h"
#include<array>
#include"../basic/global.h"

class RadixSort
{
public:
	void initial(Thread* thread);
	
	void radixSort(int spatial_hashing_index_size, std::vector<int>* hash_value, std::vector<int>* hash_triangle_index, std::vector<int>* hash_cloth_No);
	void reorder(int thread_No, int key_id);
	void setCountBucket(int thread_No, int key_id);
private:
	Thread* thread;
	void findHighestBit(unsigned int spatial_hashing_index_size, int& key_num);
	int leading_zeros(unsigned int value);
	int key_num;
	void addCount(int size, int* count_bucket, int move_byte, std::array<int, 3>* array);
	void addCount(int array_begin, int array_end, int* count_bucket, int move_byte, int* array);
	void reorder(std::array<int, 3>* array, std::array<int, 3>* stack, int move_byte, int* index_bucket, int size);
	void lsdSort(std::vector<int>* hash_value, std::vector<int>* hash_triangle_index, std::vector<int>* hash_cloth_No);
	std::vector<int> array_index_begin;
	std::vector<int>* hash_value;
	std::vector<int>* hash_triangle_index;
	std::vector<int>* hash_cloth_No;
	int thread_num;
	std::vector<std::array<int, 0x100>> histogram;

	std::vector<int> stack_hash_value;
	std::vector<int> stack_hash_triangle_index;
	std::vector<int> stack_hash_cloth_No;
	void reorder(int* hash_value, int* stack_hash_value, int* hash_triangle_index, int* stack_hash_triangle_index, int* hash_cloth_No, int* stack_hash_cloth_No, int move_byte, int* index_bucket, int array_index_start, int array_index_end);
};

