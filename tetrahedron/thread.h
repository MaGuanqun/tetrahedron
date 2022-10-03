#pragma once
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
#include"external/Eigen/Dense"
#include"external/Eigen/Sparse"

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
//class DrawCulling;
class Scene;
class XPBD;
class XPBD_IPC;
class MoveObject;
class NewtonMethod;
class SecondOrderLargeSystem;

//class MeshPatch;

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
    void assignTask(T* func, U taskType, unsigned int key_id)
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

    void assignTask(IterationMethod* func, IterationMethodFunc function_type, Eigen::VectorXd* u, Eigen::VectorXd* b,
        double* residual_norm, double omega_chebyshev, Eigen::VectorXd* u_last, Eigen::VectorXd* u_previous);
    void assignTask(IterationMethod* func, IterationMethodFunc function_type, int* vertex_index, double* coefficient, int* vertex_index_start,
        Eigen::VectorXd* x, Eigen::VectorXd* b, Eigen::VectorXd* result, double* residual_norm, int* vertex_index_thread_begin,
        Eigen::VectorXd* u_last, Eigen::VectorXd* u_previous);

    void assignTask(IterationMethod* func, std::vector<int>* vertex_index, std::vector<double>* coefficient,
        double* x, double* b, double* result, int* vertex_index_thread_begin, int sys_size);

    int thread_num;

private:
    ThreadData* threads;
    std::vector<std::future<void>> futures;

    void thread_func(ThreadData* pData);
    job create_task(TriangleMeshStruct* func, int thread_id, MeshStructFuncSendToThread function_type);// int jobNumber
    job create_task(SpatialHashing* func, int thread_id, SpatialHashingFuncSendToThread function_type);
    job create_task(SpatialHashing* func, int thread_id, SpatialHashingFuncSendToThread function_type, unsigned int key_id);
    job create_task(Cloth* func, int thread_id, ObjectFunc function_type);
    job create_task(Collider* func, int thread_id, ObjectFunc function_type);
    job create_task(BVH* func, int thread_id, BVHFunc function_type);
    job create_task(TetrahedronMeshStruct* func, int thread_id, MeshStructFuncSendToThread function_type);
    job create_task(Collision* func, int thread_id, CollisionFuncSendToThread function_type);
    job create_task(RadixSort* func, int thread_id, RadixSortFunc function_type, unsigned int key_id);

    job create_task(IterationMethod* func, int thread_id, IterationMethodFunc function_type);
    job create_task(IterationMethod* func, int thread_id, std::vector<int>* vertex_index, std::vector<double>* coefficient,
        double* x, double* b, double* result, int* vertex_index_thread_begin, int sys_size);
    job create_task(IterationMethod* func, IterationMethodFunc function_type, int thread_id, int* vertex_index, double* coefficient, int* vertex_index_start,
        Eigen::VectorXd* x, Eigen::VectorXd* b, Eigen::VectorXd* result, double* residual_norm, int* vertex_index_thread_begin,
        Eigen::VectorXd* u_last, Eigen::VectorXd* u_previous);
    job create_task(IterationMethod* func, int thread_id, IterationMethodFunc function_type, Eigen::VectorXd* u, Eigen::VectorXd* b, double* residual_norm,
        double omega_chebyshev, Eigen::VectorXd* u_last, Eigen::VectorXd* u_previous);
    job create_task(ProjectDynamic* func, int thread_id, PDFuncSendToThread function_type);
    job create_task(Tetrahedron* func, int thread_id, ObjectFunc function_type);// int jobNumber
    job create_task(Scene* func, int thread_id, SceneFuc function_type);// int jobNumber
    //job create_task(MeshPatch* func, int thread_id, MeshPatchFunc function_type);
    //job create_task(DrawCulling* func, int thread_id, DrawCullingFunc function_type, unsigned int key_id);
    job create_task(XPBD* func, int thread_id, XPBDFunc function_type);
    job create_task(XPBD_IPC* func, int thread_id, XPBD_IPC_Func function_type);
    job create_task(NewtonMethod* func, int thread_id, NewtonMethodFunc function_type);
    job create_task(SecondOrderLargeSystem* func, int thread_id, NewtonMethodFunc function_type);
    job create_task(MoveObject* func, int thread_id, MoveObjectFunc function_type, unsigned int key_id);
};


#endif // !THREAD_H


