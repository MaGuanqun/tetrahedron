#include"thread.h"
#include"project_dynamic.h"
#include"collision/spatial_hashing.h"
#include"./mesh_struct/triangle_mesh_struct.h"
#include"./object/cloth.h"
#include"./object/collider.h"
#include"./object/tetrahedron.h"
#include"collision/BVH.h"
#include"collision/collision.h"
#include"collision/parallel_radix_sort.h"


Thread::Thread()
{
    thread_num = std::thread::hardware_concurrency();
    initial();
}

Thread::~Thread()
{
    for (int i = 0; i < thread_num; ++i)
    {
        std::unique_lock<std::mutex> l(threads[i].m);
        threads[i].stop = true;
        threads[i].cv.notify_one();
    }
    // Join all the threads
    for (int i = 0; i < thread_num; ++i) { threads[i].t.join(); }
}

void Thread::initial()
{
   // std::cout << thread_num;
    threads = new ThreadData[thread_num];
   // threads.push_back(a);
	//threads.resize(thread_num);
    int tdi = 0;
    for(int i=0;i< thread_num;i++)
    {
        threads[i].id = tdi++;
        threads[i].t = std::thread(&Thread::thread_func,this, &threads[i]);
    }
}

job Thread::create_task(ProjectDynamic* func, int thread_id, PDFuncSendToThread function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case LOCAL_PROJECTION:
        k = job([func, thread_id]() {func->localProjectionPerThread(thread_id); });
        break;
    case SOLVE_SYSYTEM:
        k = job([func, thread_id]() {func->solveSystemPerThead(thread_id); });
        break;
    //case UPDATE_COLLISION_STIFFNESS:
    //    k = job([func, thread_id]() {func->updateCollisionStiffnessClothPerThread(thread_id); });
    //    break;
    case UPDATE_UV:
        k = job([func, thread_id]() {func->updateUVPerThread(thread_id); });
        break;
    case MATRIX_DECOMPOSITION:
        k = job([func, thread_id]() {func->matrixDecomposition(thread_id); });
        break;
    case UPDATE_MATRIX: {
        k = job([func, thread_id]() {func->updateMatrixPerThread(thread_id); });
        break;
    }
    //case VIRTUAL_LOCAL_PROJECTION: {
    //    k = job([func, thread_id]() {func->virtualLocalProjectionPerThread(thread_id); });
    //}
    //    break;
    //case PREPARE_WOODBURY: {
    //    k = job([func, thread_id]() {func->prepareWoodbury(thread_id); });
    //    break;
    //}
    //case SET_K_COLUMN: {
    //    k = job([func, thread_id]() {func->setKColumnPerThread(thread_id); });
    //    break;
    //}
    //case FACTORIZE_WOODBURY_K: {
    //    k = job([func, thread_id]() {func->factorizeWoodburyK(thread_id); });
    //    break;
    //}
    //case GET_FRICTION:
    //    k= job([func, thread_id]() {func->getFrictionPerThread(thread_id); });
    //    break;
    }
    return k;
}

job Thread::create_task(TriangleMeshStruct* func, int thread_id, MeshStructFuncSendToThread function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case FACE_NORMAL:
        k = job([func, thread_id]() {func->getFaceNormalPerThread(thread_id); });
        break;
    case FACE_NORMAL_RENDER:
        k = job([func, thread_id]() {func->getRenderFaceNormalPerThread(thread_id); });
        break;
    case VERTEX_NORMAL_RENDER:
        k = job([func, thread_id]() {func->getRenderVertexNormalPerThread(thread_id); });
        break;
    case VERTEX_NORMAL:
        k = job([func, thread_id]() {func->getVertexNormalPerThread(thread_id); });
        break;     
    }
    return k;
}

job Thread::create_task(TetrahedronMeshStruct* func, int thread_id, MeshStructFuncSendToThread function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case FACE_NORMAL_RENDER:
        k = job([func, thread_id]() {func->getRenderFaceNormalPerThread(thread_id); });
        break;
    case VERTEX_NORMAL_RENDER:
        k = job([func, thread_id]() {func->getRenderVertexNormalPerThread(thread_id); });
        break;
    case SET_VOLUME:
        k = job([func, thread_id]() {func->setVolume(thread_id); });
        break;
    }
    return k;
}

job Thread::create_task(Cloth* func, int thread_id, ObjectFunc function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case TRIANGLE_AABB:
        k = job([func, thread_id]() {func->getTriangleAABBPerThread(thread_id); });
        break;
    case EDGE_AABB:
        k = job([func, thread_id]() {func->getEdgeAABBPerThread(thread_id); });
        break;
    case VERTEX_AABB:
        k = job([func, thread_id]() {func->getVertexAABBPerThread(thread_id); });
        break;
    }
    return k;
}


job Thread::create_task(Collider* func, int thread_id, ObjectFunc function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case TRIANGLE_AABB:
        k = job([func, thread_id]() {func->getTriangleAABBPerThread(thread_id); });
        break;
    }
    return k;
}

job Thread::create_task(BVH* func, int thread_id, BVHFunc function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case CAL_CENTER:
        k = job([func, thread_id]() {func->calCenterPerThread(thread_id); });
        break;
    case CAL_MORTON:
        k = job([func, thread_id]() {func->calMortonCode(thread_id); });
        break;
    }
    return k;
}

job Thread::create_task(Collision* func, int thread_id, CollisionFuncSendToThread function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case RE_DETECTION:
        k = job([func, thread_id]() {func->collisionReDetection(thread_id); });
        break;
    case RESUM_TARGET_POSITION:
        k = job([func, thread_id]() {func->collisionReDetection(thread_id); });
        break;
    case FIND_TRIANGLE_PAIRS:
        k = job([func, thread_id]() {func->findAllTrianglePairs(thread_id); });
        break;
    case  FIND_PRIMITIVE_AROUND:
        k = job([func, thread_id]() {func->findPrimitivesAround(thread_id); });
        break;
    case GLOBAL_COLLISION_DETECTION:
        k = job([func, thread_id]() {func->collisionDetection(thread_id); });
        break;
    case SUM_TARGET_POSITION:
        k = job([func, thread_id]() {func->sumTargetPositionPerThread(thread_id); });
        break;
    }
    return k;
}



job Thread::create_task(SpatialHashing* func, int thread_id, SpatialHashingFuncSendToThread function_type)//
{
    job k;
    switch (function_type)
    {
    case TRIANGLE_HASHING:
        k = job([func, thread_id]() {func->triangleHashing(thread_id); });
        break;
    case SCENE_AABB:
        k = job([func, thread_id]() {func->getSceneAABB(thread_id); });
        break;
    }
    return k;
}


job Thread::create_task(RadixSort* func, int thread_id, RadixSortFunc function_type, int key_id)
{
    job k;
    switch (function_type)
    {
    case SET_COUNT_BUCKET:
        k = job([func, thread_id, key_id]() {func->setCountBucket(thread_id,key_id); });
        break;
    case REORDER:
        k = job([func, thread_id, key_id]() {func->reorder(thread_id, key_id); });
    }
    return k;
}


void Thread::assignTask(SpatialHashing* func, SpatialHashingFuncSendToThread taskType)
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


void Thread::assignTask(TriangleMeshStruct* func, MeshStructFuncSendToThread taskType)
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


void Thread::assignTask(TetrahedronMeshStruct* func, MeshStructFuncSendToThread taskType)
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


void Thread::assignTask(ProjectDynamic* func, PDFuncSendToThread taskType)
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

void Thread::assignTask(Cloth* func, ObjectFunc taskType)
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

void Thread::assignTask(Collider* func, ObjectFunc taskType)
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



void Thread::assignTask(BVH* func, BVHFunc taskType)
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


void Thread::assignTask(Collision* func, CollisionFuncSendToThread taskType)
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

void Thread::assignTask(RadixSort* func, RadixSortFunc taskType, int key_id)
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


void Thread::thread_func(ThreadData* pData)
{
    std::unique_lock<std::mutex> l(pData->m, std::defer_lock);
    while (true)
    {
        l.lock();
        // Wait until the queue won't be empty or stop is signaled
        pData->cv.wait(l, [pData]() {
            return (pData->stop || !pData->jobs.empty());
            });
        // Stop was signaled, let's exit the thread
        if (pData->stop) { return; }
        // Pop one task from the queue...
        job j = std::move(pData->jobs.front());
        pData->jobs.pop();
        l.unlock();
        // Execute the task!
        j();
    }
}