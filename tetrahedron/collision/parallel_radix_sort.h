//we use key with 4 bit for each bucket

#pragma once
#include"../thread.h"
#include<array>
#include"../basic/global.h"

class RadixSort
{
public:
	void initial(Thread* thread);
	
	void radixSort(unsigned int spatial_hashing_index_size, unsigned int* value, unsigned int* triangle_index,
		unsigned int* hash_cloth_No, unsigned int list_size, unsigned int& largest_count);
	void reorder(int thread_No, int key_id);
	void setCountBucket(int thread_No, int key_id);
	void setCountBucketMorton(int thread_No, int key_id);

	void radixSort(uint64_t max_morton_code, std::vector<uint64_t>* morton_value, unsigned int* triangle_index);
	void reorderMorton(int thread_No, int key_id);
	void initialArray(unsigned int max_length);
	void deleteArray();
	void initialMortonArray(unsigned int max_length);
	void deleteMortonArray();
	void copyArray(int thread_No);
	void copyArrayMorton(int thread_No);


private:
	Thread* thread;
	void findHighestBit(unsigned int size, int& key_num);
	void findHighestBit(uint64_t size, int& key_num);
	int leading_zeros(unsigned int value);
	int leading_zeros(uint64_t value);
	int key_num;
	void addCount(int size, unsigned int* count_bucket, unsigned int move_byte, std::array<int, 3>* array);
	void addCount(unsigned int array_begin, unsigned int array_end, unsigned int* count_bucket, unsigned int move_byte, unsigned int* array);
	void addCount(unsigned int array_begin, int array_end, unsigned int* count_bucket, int move_byte, uint64_t* array);
	void reorder(std::array<int, 3>* array, std::array<int, 3>* stack, int move_byte, int* index_bucket, int size);
	void lsdSort(unsigned int* value, unsigned int* triangle_index, unsigned int* hash_cloth_No, unsigned int list_size,
		unsigned int& largest_count);
	void lsdSort(std::vector<uint64_t>* value, unsigned int* triangle_index);
	std::vector<unsigned int> array_index_begin;
	unsigned int* value;
	std::vector<uint64_t>* morton_value;
	unsigned int* triangle_index;
	unsigned int* hash_cloth_No;
	int thread_num;
	unsigned int** histogram;

	unsigned int* stack_value;
	unsigned int* stack_triangle_index;
	unsigned int* stack_hash_cloth_No;

	uint64_t* stack_morton_value;

	void reorder(unsigned int* value, unsigned int* stack_value, unsigned int* triangle_index, unsigned int* stack_triangle_index, unsigned int* hash_cloth_No, unsigned int* stack_hash_cloth_No, unsigned int move_byte, unsigned int* index_bucket, unsigned int array_index_start, unsigned int array_index_end);
	void reorder(uint64_t* value, uint64_t* stack_value, unsigned int* triangle_index, unsigned int* stack_triangle_index, unsigned int move_byte, unsigned int* index_bucket, unsigned int array_index_start, unsigned int array_index_end);

	
};

