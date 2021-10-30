//we use key with 4 bit for each bucket

#pragma once
#include"../thread.h"
#include<array>

class RadixSort
{
public:
	void initial(Thread* thread);
	
	void radixSort(int size, std::vector<std::array<int, 2>>& array);

	
private:
	Thread* thread;
	void findHighestBit(unsigned int spatial_hashing_index_size, int& highest_bit, int& key_num);
	int leading_zeros(unsigned int value);
	int highest_bit;
	int key_num;
	void addCount(int size, int* count_bucket, int move_byte, std::vector<std::array<int, 2>>& array);
	void reorder(std::vector<std::array<int, 2>>& array, std::vector<std::array<int, 2>>& stack, int move_byte, int* index_bucket, int size);
	void lsdSort(std::vector<std::array<int, 2>>& array);
};

