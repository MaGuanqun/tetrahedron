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
class TetrahedronMeshStruct;
class Cloth;
class Collider;
class Tetrahedron;
class BVH;
class Collision;
class SpatialHashing;
class RadixSort;
class IterationMethod;

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
private:
    ThreadData* threads;
    std::vector<std::future<void>> futures;

    void thread_func(ThreadData* pData);
    job create_task(ProjectDynamic* func, int thread_id, PDFuncSendToThread function_type);
    job create_task(TriangleMeshStruct* func, int thread_id, MeshStructFuncSendToThread function_type);// int jobNumber
    job create_task(SpatialHashing* func, int thread_id, SpatialHashingFuncSendToThread function_type);
    job create_task(Cloth* func, int thread_id, ObjectFunc function_type);
    job create_task(Collider* func, int thread_id, ObjectFunc function_type);
    job create_task(BVH* func, int thread_id, BVHFunc function_type);
    job create_task(TetrahedronMeshStruct* func, int thread_id, MeshStructFuncSendToThread function_type);
    job create_task(Collision* func, int thread_id, CollisionFuncSendToThread function_type);
    job create_task(RadixSort* func, int thread_id, RadixSortFunc function_type, int key_id);
    job create_task(IterationMethod* func, int thread_id, IterationMethodFunc function_type);

public:
    Thread();
    void initial();
	~Thread();
    template <class T, typename U>
    void assignTask(T* func, U taskType)
    {
        for (int i = 0; i < thread_num; ++i)
        {
            // std::cout << threads[i].id << std::endl;
            job j = create_task(func, threads[i].id, taskType);
            futures.push_back(j.get_future());
            std::unique_lock<std::mutex> l(threads[i].m);
            threads[i].jobs.push(std::move(j));
            // Notify the thread that there is work do to...
            threads[i].cv.notify_one();
        }
        for (auto& f : futures) { f.wait(); }
        futures.clear();
    }
    template <class T, typename U>
    void assignTask(T* func, U taskType, int key_id)
    {
        for (int i = 0; i < thread_num; ++i)
        {
            // std::cout << threads[i].id << std::endl;
            job j = create_task(func, threads[i].id, taskType, key_id);
           
            futures.push_back(j.get_future());
            std::unique_lock<std::mutex> l(threads[i].m);
            threads[i].jobs.push(std::move(j));
            // Notify the thread that there is work do to...
            threads[i].cv.notify_one();
        }
        for (auto& f : futures) { f.wait(); }
        futures.clear();
    }  
    int thread_num;
};


#endif // !THREAD_H

#pragma once
