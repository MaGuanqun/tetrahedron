#include"parallel_radix_sort.h"

#define BIT_PER_KEY 8

void RadixSort::initial(Thread* thread, bool is_single_thread)
{
    this->is_single_thread = is_single_thread;
    if (is_single_thread) {
        histogram_single_thread = new unsigned int[256];
    }
    else {
        this->thread = thread;
        thread_num = thread->thread_num;
        array_index_begin.resize(thread_num + 1);
        histogram = new unsigned int* [thread_num];
        for (int i = 0; i < thread_num; ++i) {
            histogram[i] = new unsigned int[256];
        }
    }
}

void RadixSort::findHighestBit(unsigned int size, unsigned  int& key_num)
{
	unsigned int highest_bit = 32 - leading_zeros(size);
    key_num = highest_bit / BIT_PER_KEY;
    if (highest_bit % BIT_PER_KEY > 0) {
        key_num += 1;
    }
}

void RadixSort::findHighestBit(uint64_t size, unsigned int& key_num)
{
    unsigned int highest_bit = 64 - leading_zeros(size);
    key_num = highest_bit / BIT_PER_KEY;
    if (highest_bit % BIT_PER_KEY > 0) {
        key_num += 1;
    }
}



void RadixSort::initialArray(unsigned int max_length)
{
    stack_value = new unsigned int[max_length];
    stack_triangle_index= new unsigned int[max_length];
    stack_hash_cloth_No = new unsigned int[max_length];
    this->array_size = max_length;
}

void RadixSort::deleteArray()
{
    delete[] stack_value;
    delete[] stack_triangle_index;
    delete[] stack_hash_cloth_No;
    if (is_single_thread) {
        delete[] histogram_single_thread;
    }
    else {
        delete[] histogram;
    }  
}

void RadixSort::initialMortonArray(unsigned int max_length)
{
    stack_morton_value = new uint64_t[max_length];
    stack_triangle_index = new unsigned int[max_length];
    this->array_size = max_length;
}

void RadixSort::deleteMortonArray()
{
    delete[] stack_morton_value;
    delete[] stack_triangle_index;
    if (is_single_thread) {
        delete[] histogram_single_thread;
    }
    else {
        delete[] histogram;
    }
}

void RadixSort::radixSort(uint64_t max_morton_code, std::vector<uint64_t>* morton_value, unsigned int* triangle_index)
{
    arrangeIndex(thread->thread_num, morton_value->size(), array_index_begin.data());
    findHighestBit(max_morton_code, key_num);
    this->morton_value = morton_value;
    this->triangle_index = triangle_index;
    lsdSort(morton_value, triangle_index);
}

void RadixSort::radixSort(unsigned int spatial_hashing_index_size, unsigned int* value, unsigned int* triangle_index, 
    unsigned int* hash_cloth_No, unsigned int list_size)
    //largest count is the count of elements which has is 11111111 in the largest part
{
    findHighestBit(spatial_hashing_index_size, key_num);
   
    this->value = value;
    this->triangle_index = triangle_index;
    this->hash_cloth_No = hash_cloth_No;

    if (is_single_thread) {
        lsdSortSingleThread(value, triangle_index, hash_cloth_No, list_size, array_size);
    }
    else {
        arrangeIndex(thread->thread_num, list_size, array_index_begin.data());
        lsdSort(value, triangle_index, hash_cloth_No, list_size, array_size);
    }
}


void RadixSort::lsdSortSingleThread(unsigned int* value, unsigned int* triangle_index, unsigned int* hash_cloth_No, unsigned int list_size,
    unsigned int& largest_count)
{
    int s, t;
    for (unsigned int j = 0; j < key_num; ++j) {
        setCountBucket(j);        
        if (j == key_num - 1) {
            largest_count = histogram_single_thread[255];
        }
        s = 0;
        for (int i = 0; i < 0x100; ++i) {
            t = s + histogram_single_thread[i];
            histogram_single_thread[i] = s;
            s = t;
        }
        reorder(j);
    }
    if (key_num % 2 == 1) {
        memcpy(value, stack_value, 4 * list_size);
        memcpy(triangle_index, stack_triangle_index, 4 * list_size);
        memcpy(hash_cloth_No, stack_hash_cloth_No, 4 * list_size);
    }
}


void RadixSort::lsdSort(unsigned int* value, unsigned int* triangle_index, unsigned int* hash_cloth_No, unsigned int list_size,
    unsigned int& largest_count)
{
    int s, t;
    for (unsigned int j = 0; j < key_num; ++j) {
        thread->assignTask(this, SET_COUNT_BUCKET,j);
        s = 0;
        if (j == key_num - 1) {
            largest_count = histogram[0][255];
            for (int i = 1; i < thread_num; ++i) {
                largest_count+= histogram[i][255];
            }
        }
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
        thread->assignTask(this, COPY_ARRAY, 0);
        //memcpy(value, stack_value, 4 * list_size);
        //memcpy(triangle_index, stack_triangle_index, 4 * list_size);
        //memcpy(hash_cloth_No, stack_hash_cloth_No, 4 * list_size);

    }
}

//COPY_ARRAY
void RadixSort::copyArray(int thread_No)
{
    memcpy(value + array_index_begin[thread_No], stack_value + array_index_begin[thread_No], 4 * (array_index_begin[thread_No + 1] - array_index_begin[thread_No]));
    memcpy(triangle_index + array_index_begin[thread_No], stack_triangle_index + array_index_begin[thread_No], 4 * (array_index_begin[thread_No + 1] - array_index_begin[thread_No]));
    memcpy(hash_cloth_No + array_index_begin[thread_No], stack_hash_cloth_No + array_index_begin[thread_No], 4 * (array_index_begin[thread_No + 1] - array_index_begin[thread_No]));
}




void RadixSort::lsdSort(std::vector<uint64_t>* value, unsigned int* triangle_index)
{
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
        memcpy(value->data(), stack_morton_value, 8 * value->size());
        memcpy(triangle_index, stack_triangle_index, 4 * value->size());
    }
}
//void RadixSort::lsdSort(std::vector<int>* value, std::vector<int>* triangle_index, std::vector<int>* hash_cloth_No)
//{
//    std::vector<int> stack_(value->size());
//    std::vector<int> stack_triangle_index(value->size());
//    std::vector<int> stack_cloth_No(value->size());
//    int s, t;
//    for (int j = 0; j < key_num; ++j) {
//        memset(histogram[0].data(), 0, 4 * 0x100);
//        addCount(0,value->size(), histogram[0].data(), 8 * j, value->data());      
//        s = 0;
//        for (int i = 0; i < 0x100; ++i) {
//            t = s + histogram[0][i];
//            histogram[0][i] = s;
//            s = t;
//        }
//        if (j % 2 == 0) {
//            reorder(value->data(), stack_.data(), triangle_index->data(), stack_triangle_index.data(), hash_cloth_No->data(), stack_cloth_No.data(), 8 * j, histogram[0].data(),0, value->size());
//        }
//        else {
//            reorder(stack_.data(), value->data(), stack_triangle_index.data(), triangle_index->data(),  stack_cloth_No.data(), hash_cloth_No->data(), 8 * j, histogram[0].data(), 0, value->size());
//        }
//    }
//    if (key_num % 2 == 1) {
//        (*value) = stack_;
//        (*triangle_index) = stack_triangle_index;
//        (*hash_cloth_No) = stack_cloth_No;
//    }
//}

void RadixSort::setCountBucket(unsigned int key_id)
{
    memset(histogram_single_thread, 0, 4 * 0x100);
    if (key_id % 2 == 0) {
        addCount(0, array_size, histogram_single_thread, 8 * key_id, value);
    }
    else {
        addCount(0, array_size, histogram_single_thread, 8 * key_id, stack_value);
    }
}


//SET_COUNT_BUCKET
void RadixSort::setCountBucket(int thread_No, unsigned int key_id)
{
    unsigned int* count_bucket;
    count_bucket = histogram[thread_No];
    memset(count_bucket, 0, 4 * 0x100);
    if (key_id % 2 == 0) {
        addCount(array_index_begin[thread_No], array_index_begin[thread_No + 1], count_bucket, 8 * key_id, value);
    }
    else {
        addCount(array_index_begin[thread_No], array_index_begin[thread_No + 1], count_bucket, 8 * key_id, stack_value);
    }
}

//SET_COUNT_BUCKET_MORTON
void RadixSort::setCountBucketMorton(int thread_No, unsigned int key_id)
{
    unsigned int* count_bucket;
    count_bucket = histogram[thread_No];
    memset(count_bucket, 0, 4 * 0x100);
    if (key_id % 2 == 0) {
        addCount(array_index_begin[thread_No], array_index_begin[thread_No + 1], count_bucket, 8 * key_id, morton_value->data());
    }
    else {
        addCount(array_index_begin[thread_No], array_index_begin[thread_No + 1], count_bucket, 8 * key_id, stack_morton_value);
    }
}

void RadixSort::reorder(unsigned int key_id)
{
    if (key_id % 2 == 0) {
        reorder(value, stack_value, triangle_index, stack_triangle_index, hash_cloth_No, stack_hash_cloth_No, 8 * key_id, histogram_single_thread, 0, array_size);
    }
    else {
        reorder(stack_value, value, stack_triangle_index, triangle_index, stack_hash_cloth_No, hash_cloth_No, 8 * key_id, histogram_single_thread, 0, array_size);
    }
}


//REORDER
void RadixSort::reorder(int thread_No, unsigned int key_id)
{
    if (key_id % 2 == 0) {
        reorder(value, stack_value, triangle_index, stack_triangle_index, hash_cloth_No, stack_hash_cloth_No,  8 * key_id, histogram[thread_No], array_index_begin[thread_No], array_index_begin[thread_No+1]);
    }
    else {
        reorder(stack_value, value, stack_triangle_index, triangle_index, stack_hash_cloth_No, hash_cloth_No, 8 * key_id, histogram[thread_No], array_index_begin[thread_No], array_index_begin[thread_No + 1]);
    }
}


//MORTON_REORDER
void RadixSort::reorderMorton(int thread_No, unsigned int key_id)
{
    if (key_id % 2 == 0) {
        reorder(morton_value->data(), stack_morton_value, triangle_index, stack_triangle_index, 8 * key_id, histogram[thread_No], array_index_begin[thread_No], array_index_begin[thread_No + 1]);
    }
    else {
        reorder(stack_morton_value, morton_value->data(), stack_triangle_index, triangle_index, 8 * key_id, histogram[thread_No], array_index_begin[thread_No], array_index_begin[thread_No + 1]);
    }
}


void RadixSort::reorder(uint64_t* value, uint64_t* stack_value, unsigned int* triangle_index, unsigned int* stack_triangle_index, unsigned int move_byte, unsigned int* index_bucket, unsigned int array_index_start, unsigned int array_index_end)
{
    unsigned int* index;
    if (move_byte == 0) {
        for (unsigned int i = array_index_start; i < array_index_end; ++i) {
            index = &index_bucket[value[i] & 0xff];
            stack_value[*index] = value[i];
            stack_triangle_index[*index] = triangle_index[i];
            (*index)++;
        }
    }
    else {
        for (unsigned int i = array_index_start; i < array_index_end; ++i) {
            index = &index_bucket[(value[i] >> move_byte) & 0xff];
            stack_value[*index] = value[i];
            stack_triangle_index[*index] = triangle_index[i];
            (*index)++;
        }
    }
}

void RadixSort::reorder(unsigned int* value, unsigned int* stack_value, unsigned int* triangle_index, unsigned int* stack_triangle_index,unsigned int* hash_cloth_No, unsigned int* stack_hash_cloth_No, unsigned int move_byte, unsigned int* index_bucket, unsigned int array_index_start, unsigned int array_index_end)
{
    unsigned int* index;
    if (move_byte == 0) {
        for (unsigned int i = array_index_start; i < array_index_end; ++i) {
            index = &index_bucket[value[i] & 0xff];
            stack_value[*index] = value[i];
            stack_triangle_index[*index] = triangle_index[i];
            stack_hash_cloth_No[*index] = hash_cloth_No[i];
            (*index)++;
        }
    }
    else {
        for (unsigned int i = array_index_start; i < array_index_end; ++i) {
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

void RadixSort::addCount(unsigned int array_begin, unsigned int array_end, unsigned int* count_bucket, unsigned int move_byte, unsigned int* array)
{
    if (move_byte == 0) {
        for (unsigned int i = array_begin; i < array_end; ++i) {
            count_bucket[array[i] & 0xff]++;
        }
    }
    else {
        for (unsigned int i = array_begin; i < array_end; ++i) {
            count_bucket[(array[i] >> move_byte) & 0xff]++;
        }
    }
}

void RadixSort::addCount(unsigned int array_begin, int array_end, unsigned int* count_bucket, int move_byte, uint64_t* array)
{
    if (move_byte == 0) {
        for (unsigned int i = array_begin; i < array_end; ++i) {
            count_bucket[array[i] & 0xff]++;
        }
    }
    else {
        for (unsigned int i = array_begin; i < array_end; ++i) {
            count_bucket[(array[i] >> move_byte) & 0xff]++;
        }
    }
}

void RadixSort::addCount(int size, unsigned int* count_bucket, unsigned int move_byte, std::array<int, 3>* array)
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