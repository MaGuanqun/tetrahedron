#include"parallel_radix_sort.h"

#define BIT_PER_KEY 8

void RadixSort::initial(Thread* thread)
{
	this->thread = thread;
    thread_num = thread->thread_num;
    array_index_begin.resize(thread_num+1);
    histogram.resize(thread_num);
}

void RadixSort::findHighestBit(unsigned int size, int& key_num)
{
	int highest_bit = 32 - leading_zeros(size);
    key_num = highest_bit / BIT_PER_KEY;
    if (highest_bit % BIT_PER_KEY > 0) {
        key_num += 1;
    }
}

void RadixSort::radixSort(int spatial_hashing_index_size, std::vector<std::array<int, 3>>* array)
{
    arrangeIndex(thread->thread_num, array->size(), array_index_begin);
	findHighestBit(spatial_hashing_index_size, key_num);
    this->array = array;
    lsdSort(array);   
}


void RadixSort::lsdSort(std::vector<std::array<int, 3>>* array)
{
    stack.resize(array->size());
    int s, t;
    for (int j = 0; j < key_num; ++j) {
        thread->assignTask(this, SET_COUNT_BUCKET,j);
        s = 0;
        for (int i = 0; i < 0x100; ++i) {
            for (int k = 0; k < thread_num; ++k) {
                t = s + histogram[k][i];
                histogram[k][i] = s;
                s = t;
            }
        }
        thread->assignTask(this, REORDER, j);
    }
    if (key_num % 2 == 1) {
        (*array) = stack;
    }
}

//void RadixSort::lsdSort(std::vector<std::array<int, 3>>* array)
//{
//    std::vector<std::array<int, 0x100>> index(key_num);
//    std::vector<std::array<int, 0x100>> count(key_num, { 0 });
//    std::vector<std::array<int, 3>> stack_(array->size());
//    int* index_bucket;
//    int* count_bucket;
//    for (int j = 0; j < key_num; ++j) {
//        count_bucket = count[j].data();
//        addCount(array->size(), count_bucket, 8 * j, array->data());
//        index_bucket = index[j].data();
//        count_bucket = count[j].data();
//        index_bucket[0] = 0;
//        for (int i = 1; i < 0xff; ++i) {
//            index_bucket[i] = index_bucket[i - 1] + count_bucket[i - 1];
//        }
//        if (j % 2 == 0) {
//            reorder(array->data(), stack_.data(), 8 * j, index[j].data(), array->size());
//        }
//        else {
//            reorder(stack_.data(), array->data(), 8 * j, index[j].data(), array->size());
//        }
//    }
//    if (key_num % 2 == 1) {
//        (*array) = stack;
//    }
//}


//SET_COUNT_BUCKET
void RadixSort::setCountBucket(int thread_No, int key_id)
{
    int* count_bucket;
    count_bucket = histogram[thread_No].data();
    memset(count_bucket, 0, 4 * 0x100);
    if (key_id % 2 == 0) {
        addCount(array_index_begin[thread_No], array_index_begin[thread_No + 1], count_bucket, 8 * key_id, array->data());
    }
    else {
        addCount(array_index_begin[thread_No], array_index_begin[thread_No + 1], count_bucket, 8 * key_id, stack.data());
    }
}



//REORDER
void RadixSort::reorder(int thread_No, int key_id)
{
    if (key_id % 2 == 0) {
        reorder(array->data(), stack.data(), 8 * key_id, histogram[thread_No].data(), array_index_begin[thread_No], array_index_begin[thread_No+1]);
    }
    else {
        reorder(stack.data(), array->data(), 8 * key_id, histogram[thread_No].data(), array_index_begin[thread_No], array_index_begin[thread_No + 1]);
    }
}



void RadixSort::reorder(std::array<int, 3>* array, std::array<int, 3>* stack, int move_byte, int* index_bucket, int array_index_start, int array_index_end)
{
    if (move_byte == 0) {
        for (int i = array_index_start; i < array_index_end; ++i) {
            stack[index_bucket[array[i][2] & 0xff]++] = array[i];
        }
    }
    else {
        for (int i = array_index_start; i < array_index_end; ++i) {
            stack[index_bucket[(array[i][2] >> move_byte) & 0xff]++] = array[i];
        }
    }
}

void RadixSort::reorder(std::array<int, 3>* array, std::array<int, 3>* stack, int move_byte, int* index_bucket, int size)
{
    if (move_byte == 0) {
        for (int i = 0; i < size; ++i) {
            stack[index_bucket[array[i][2] & 0xff]++] = array[i];
        }
    }
    else {
        for (int i = 0; i < size; ++i) {
            stack[index_bucket[(array[i][2]>>move_byte) & 0xff]++] = array[i];
        }
    }
}

void RadixSort::addCount(int array_begin, int array_end, int* count_bucket, int move_byte, std::array<int, 3>* array)
{
    if (move_byte == 0) {
        for (int i = array_begin; i < array_end; ++i) {
            count_bucket[array[i][2] & 0xff]++;
        }
    }
    else {
        for (int i = array_begin; i < array_end; ++i) {
            count_bucket[(array[i][2] >> move_byte) & 0xff]++;
        }
    }
}

void RadixSort::addCount(int size, int* count_bucket, int move_byte, std::array<int, 3>* array)
{
    if (move_byte == 0) {
        for (int i = 0; i < size; ++i) {
            count_bucket[array[i][2] & 0xff]++;
        }
    }
    else {
        for (int i = 0; i < size; ++i) {
            count_bucket[(array[i][2]>>move_byte) & 0xff]++;
        }
    }
}



int RadixSort::leading_zeros(unsigned int value) {
    int count = 0;
    if ((value & 0xffff0000u) == 0) {
        count += 16;
        value <<= 16;
    }
    if ((value & 0xff000000u) == 0) {
        count += 8;
        value <<= 8;
    }
    if ((value & 0xf0000000u) == 0) {
        count += 4;
        value <<= 4;
    }
    if ((value & 0xc0000000u) == 0) {
        count += 2;
        value <<= 2;
    }
    if ((value & 0x80000000u) == 0) {
        count += 1;
    }
    return count;
}