#ifndef THREAD_H
#define THREAD_H
#include <thread>
#include <iostream>
#include <sstream>
#include <future>
#include <queue>
#include <condition_variable>
#include <mutex>
#include"basic/enum_setting.h"


class ProjectDynamic;
class TriangleMeshStruct;
class Cloth;
class Collider;
class Tetrohedron;
class BVH;

using job = std::packaged_task<void()>;

struct ThreadData
{
    int id; // Could use thread::id, but this is filled before the thread is started
    std::thread t; // The thread object
    std::queue<job> jobs; // The job queue
    std::condition_variable cv; // The condition variable to wait for threads
    std::mutex m; // Mutex used for avoiding data races
    bool stop = false; // When set, this flag tells the thread that it should exit
};

class Thread
{
public:
    Thread();
    void initial();
	~Thread();
    void assignTask(ProjectDynamic* func, PDFuncSendToThread taskType);
    //void assignTask(SpatialHashing* func, SpatialHashingFuncSendToThread taskType, int cloth_No, int compare_cloth_No);
    void assignTask(TriangleMeshStruct* func, MeshStructFuncSendToThread taskType);
    //void assignTask(SpatialHashing* func, SpatialHashingFuncSendToThread taskType);
    void assignTask(Cloth* func, ObjectFunc taskType);
    void assignTask(Collider* func, ObjectFunc taskType);
    void assignTask(BVH* func, BVHFunc taskType);
   
    int thread_num;
private:
    ThreadData* threads;
    std::vector<std::future<void>> futures;
    
    void thread_func(ThreadData* pData);
    job create_task(ProjectDynamic* func, int thread_id, PDFuncSendToThread function_type);
    //job create_task(SpatialHashing* func, int thread_id, SpatialHashingFuncSendToThread function_type, int cloth_No, int compare_cloth_No);
    job create_task(TriangleMeshStruct* func, int thread_id, MeshStructFuncSendToThread function_type);// int jobNumber
    //job create_task(SpatialHashing* func, int thread_id, SpatialHashingFuncSendToThread function_type);
    job create_task(Cloth* func, int thread_id, ObjectFunc function_type);
    job create_task(Collider* func, int thread_id, ObjectFunc function_type);
    job create_task(BVH* func, int thread_id, BVHFunc function_type);
};


#endif // !THREAD_H

#pragma once
