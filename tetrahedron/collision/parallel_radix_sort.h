//we use key with 4 bit for each bucket

#pragma once
#include"../thread.h"
#include<array>
#include"../basic/global.h"

class RadixSort
{
public:
	void initial(Thread* thread);
	
	void radixSort(int size, std::vector<std::array<int, 3>>* array);
	void reorder(int thread_No, int key_id);
	void setCountBucket(int thread_No, int key_id);
private:
	Thread* thread;
	void findHighestBit(unsigned int spatial_hashing_index_size, int& key_num);
	int leading_zeros(unsigned int value);
	int key_num;
	void addCount(int size, int* count_bucket, int move_byte, std::array<int, 3>* array);
	void addCount(int array_begin, int array_end, int* count_bucket, int move_byte, std::array<int, 3>* array);
	void reorder(std::array<int, 3>* array, std::array<int, 3>* stack, int move_byte, int* index_bucket, int size);
	void lsdSort(std::vector<std::array<int, 3>>* array);
	std::vector<int> array_index_begin;
	std::vector<std::array<int, 3>>* array;
	int thread_num;
	std::vector<std::array<int, 0x100>> histogram;

	std::vector<std::array<int, 3>> stack;
	void reorder(std::array<int, 3>* array, std::array<int, 3>* stack, int move_byte, int* index_bucket, int array_index_start, int array_index_end);
};

