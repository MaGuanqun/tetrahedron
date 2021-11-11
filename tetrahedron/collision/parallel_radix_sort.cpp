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

void RadixSort::findHighestBit(uint64_t size, int& key_num)
{
    int highest_bit = 64 - leading_zeros(size);
    key_num = highest_bit / BIT_PER_KEY;
    if (highest_bit % BIT_PER_KEY > 0) {
        key_num += 1;
    }
}

void RadixSort::radixSort(uint64_t max_morton_code, std::vector<uint64_t>* morton_value, std::vector<int>* triangle_index)
{
    arrangeIndex(thread->thread_num, morton_value->size(), array_index_begin);
    findHighestBit(max_morton_code, key_num);
    this->morton_value = morton_value;
    this->triangle_index = triangle_index;
    lsdSort(morton_value, triangle_index);
}

void RadixSort::radixSort(unsigned int spatial_hashing_index_size, std::vector<int>* value, std::vector<int>* triangle_index, std::vector<int>* hash_cloth_No)
{
    arrangeIndex(thread->thread_num, value->size(), array_index_begin);
	findHighestBit(spatial_hashing_index_size, key_num);
    this->value = value;
    this->triangle_index = triangle_index;
    this->hash_cloth_No = hash_cloth_No;
    lsdSort(value, triangle_index, hash_cloth_No);
}


void RadixSort::lsdSort(std::vector<int>* value, std::vector<int>* triangle_index, std::vector<int>* hash_cloth_No)
{
    stack_value.resize(value->size());
    stack_triangle_index.resize(value->size());
    stack_hash_cloth_No.resize(value->size());
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
        (*value) = stack_value;
        (*triangle_index) = stack_triangle_index;
        (*hash_cloth_No) = stack_hash_cloth_No;
    }
}


void RadixSort::lsdSort(std::vector<uint64_t>* value, std::vector<int>* triangle_index)
{
    stack_morton_value.resize(value->size());
    stack_triangle_index.resize(value->size());
    int s, t;
    for (int j = 0; j < key_num; ++j) {
        thread->assignTask(this, SET_COUNT_BUCKET_MORTON, j);
        s = 0;
        for (int i = 0; i < 0x100; ++i) {
            for (int k = 0; k < thread_num; ++k) {
                t = s + histogram[k][i];
                histogram[k][i] = s;
                s = t;
            }
        }
        thread->assignTask(this, MORTON_REORDER, j);
    }
    if (key_num % 2 == 1) {
        (*value) = stack_morton_value;
        (*triangle_index) = stack_triangle_index;      
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
        addCount(array_index_begin[thread_No], array_index_begin[thread_No + 1], count_bucket, 8 * key_id, value->data());
    }
    else {
        addCount(array_index_begin[thread_No], array_index_begin[thread_No + 1], count_bucket, 8 * key_id, stack_value.data());
    }
}

//SET_COUNT_BUCKET_MORTON
void RadixSort::setCountBucketMorton(int thread_No, int key_id)
{
    int* count_bucket;
    count_bucket = histogram[thread_No].data();
    memset(count_bucket, 0, 4 * 0x100);
    if (key_id % 2 == 0) {
        addCount(array_index_begin[thread_No], array_index_begin[thread_No + 1], count_bucket, 8 * key_id, morton_value->data());
    }
    else {
        addCount(array_index_begin[thread_No], array_index_begin[thread_No + 1], count_bucket, 8 * key_id, stack_morton_value.data());
    }
}

//REORDER
void RadixSort::reorder(int thread_No, int key_id)
{
    if (key_id % 2 == 0) {
        reorder(value->data(), stack_value.data(), triangle_index->data(), stack_triangle_index.data(), hash_cloth_No->data(), stack_hash_cloth_No.data(),  8 * key_id, histogram[thread_No].data(), array_index_begin[thread_No], array_index_begin[thread_No+1]);
    }
    else {
        reorder(stack_value.data(), value->data(), stack_triangle_index.data(), triangle_index->data(), stack_hash_cloth_No.data(), hash_cloth_No->data(), 8 * key_id, histogram[thread_No].data(), array_index_begin[thread_No], array_index_begin[thread_No + 1]);
    }
}


//MORTON_REORDER
void RadixSort::reorderMorton(int thread_No, int key_id)
{
    if (key_id % 2 == 0) {
        reorder(morton_value->data(), stack_morton_value.data(), triangle_index->data(), stack_triangle_index.data(), 8 * key_id, histogram[thread_No].data(), array_index_begin[thread_No], array_index_begin[thread_No + 1]);
    }
    else {
        reorder(stack_morton_value.data(), morton_value->data(), stack_triangle_index.data(), triangle_index->data(), 8 * key_id, histogram[thread_No].data(), array_index_begin[thread_No], array_index_begin[thread_No + 1]);
    }
}


void RadixSort::reorder(uint64_t* value, uint64_t* stack_value, int* triangle_index, int* stack_triangle_index, int move_byte, int* index_bucket, int array_index_start, int array_index_end)
{
    int* index;
    if (move_byte == 0) {
        for (int i = array_index_start; i < array_index_end; ++i) {
            index = &index_bucket[value[i] & 0xff];
            stack_value[*index] = value[i];
            stack_triangle_index[*index] = triangle_index[i];
            (*index)++;
        }
    }
    else {
        for (int i = array_index_start; i < array_index_end; ++i) {
            index = &index_bucket[(value[i] >> move_byte) & 0xff];
            stack_value[*index] = value[i];
            stack_triangle_index[*index] = triangle_index[i];
            (*index)++;
        }
    }
}

void RadixSort::reorder(int* value, int* stack_value, int* triangle_index, int* stack_triangle_index,int* hash_cloth_No, int* stack_hash_cloth_No, int move_byte, int* index_bucket, int array_index_start, int array_index_end)
{
    int* index;
    if (move_byte == 0) {
        for (int i = array_index_start; i < array_index_end; ++i) {
            index = &index_bucket[value[i] & 0xff];
            stack_value[*index] = value[i];
            stack_triangle_index[*index] = triangle_index[i];
            stack_hash_cloth_No[*index] = hash_cloth_No[i];
            (*index)++;
        }
    }
    else {
        for (int i = array_index_start; i < array_index_end; ++i) {
            index = &index_bucket[(value[i] >> move_byte) & 0xff];
            stack_value[*index] = value[i];
            stack_triangle_index[*index] = triangle_index[i];
            stack_hash_cloth_No[*index] = hash_cloth_No[i];
            (*index)++;            
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

void RadixSort::addCount(int array_begin, int array_end, int* count_bucket, int move_byte, int* array)
{
    if (move_byte == 0) {
        for (int i = array_begin; i < array_end; ++i) {
            count_bucket[array[i] & 0xff]++;
        }
    }
    else {
        for (int i = array_begin; i < array_end; ++i) {
            count_bucket[(array[i] >> move_byte) & 0xff]++;
        }
    }
}

void RadixSort::addCount(int array_begin, int array_end, int* count_bucket, int move_byte, uint64_t* array)
{
    if (move_byte == 0) {
        for (int i = array_begin; i < array_end; ++i) {
            count_bucket[array[i] & 0xff]++;
        }
    }
    else {
        for (int i = array_begin; i < array_end; ++i) {
            count_bucket[(array[i] >> move_byte) & 0xff]++;
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

int RadixSort::leading_zeros(uint64_t value) {
    int count = 0;
    if ((value & 0xffffffff00000000u) == 0) {
        count += 32;
        value <<= 32;
    }
    if ((value & 0xffff000000000000u) == 0) {
        count += 16;
        value <<= 16;
    }
    if ((value & 0xff00000000000000u) == 0) {
        count += 8;
        value <<= 8;
    }
    if ((value & 0xf000000000000000u) == 0) {
        count += 4;
        value <<= 4;
    }
    if ((value & 0xc000000000000000u) == 0) {
        count += 2;
        value <<= 2;
    }
    if ((value & 0x8000000000000000u) == 0) {
        count += 1;
    }
    return count;
}