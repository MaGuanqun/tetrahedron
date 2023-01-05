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
#include"iteration_method.h"
//#include"collision/drawCulling.h"
#include"scene.h"
#include"./XPBD/XPBD.h"
#include"./basic/move_object.h"
#include"newton_method.h"
#include"XPBD_large_system.h"
#include"XPBD_IPC.h"
//#include"collision/mesh_patch.h"

Thread::Thread()
{
    thread_num = 1;// std::thread::hardware_concurrency();
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
    for (int i = 0; i < thread_num; i++)
    {
        threads[i].id = tdi++;
        threads[i].t = std::thread(&Thread::thread_func, this, &threads[i]);
    }
}

job Thread::create_task(ProjectDynamic* func, int thread_id, PDFuncSendToThread function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case LOCAL_PROJECTION:
        k = job([func, thread_id]() {func->localProjectionPerThread(thread_id, true); });
        break;
    case LOCAL_PROJECTION_WITHOUT_ENERGY:
        k = job([func, thread_id]() {func->localProjectionPerThread(thread_id, false); });
        break;
    case COLLISION_FREE_POSITION:
        k = job([func, thread_id]() {func->computeCollisionFreePosition(thread_id); });
        break;
    case CONSTRUCT_B:
        k = job([func, thread_id]() {func->constructbPerThead(thread_id, true); });
        break;
    case CONSTRUCT_B_WITHOUT_COLLISION:
        k = job([func, thread_id]() {func->constructbPerThead(thread_id, false); });
        break;
    case SOLVE_WITH_COLLISION:
        k = job([func, thread_id]() {func->solveSystemPerThread(thread_id, true); });
        break;
    case SOLVE_WITHOUT_COLLISION:
        k = job([func, thread_id]() {func->solveSystemPerThread(thread_id, false); });
        break;
    case COMPUTE_DISPLACEMENT:
        k = job([func, thread_id]() {func->computeDisplacement(thread_id); });
        break;
    case COMPUTE_ENERGY:
        k = job([func, thread_id]() {func->computeEnergyPerThread(thread_id); });
        break;
        //case UPDATE_COLLISION_STIFFNESS:
        //    k = job([func, thread_id]() {func->updateCollisionStiffnessClothPerThread(thread_id); });
        //    break;
    case UPDATE_UV:
        k = job([func, thread_id]() {func->updateUVPerThread(thread_id); });
        break;
    case UPDATE_MATRIX: {
        k = job([func, thread_id]() {func->updateMatrixPerThread(thread_id); });
        break;
    }
    case UPDATE_DIAGONAL: {
        k = job([func, thread_id]() {func->updateDiagonalPerThread(thread_id); });
        break;
    }
                        //case LOCAL_EDGE_LENGTH_PROJECTION:
                        //    k = job([func, thread_id]() {func->localEdgeLengthProjectionPerThread(thread_id); });
                        //    break;
    case TEST_LOCAL_PROJECTION: {
        k = job([func, thread_id]() {func->testLocalProjectionPerThread(thread_id); });
        break;
    }
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
    case VERTEX_NORMAL_FROM_RENDER:
        k = job([func, thread_id]() {func->getVertexNormalFromRenderPerThread(thread_id); });
        break;
    case SORT_TRIANGLE_EDGE_AROUND_TRIANGLE_EDGE:
        k = job([func, thread_id]() {func->setFaceEdgeAroundFace(thread_id); });
        break;
    case SORT_TRIANGLE_AROUND_VERTEX_EDGE:
        k = job([func, thread_id]() {func->sortTriangleAroundVertexEdge(thread_id); });
        break;
    }
    return k;
}


job Thread::create_task(Scene* func, int thread_id, SceneFuc function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case TEST_ARRAY:
        k = job([func, thread_id]() {func->testForWritetToArray(thread_id); });
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
    case FACE_NORMAL:
        k = job([func, thread_id]() {func->getFaceNormalPerThread(thread_id); });
        break;
    case VERTEX_NORMAL:
        k = job([func, thread_id]() {func->getVertexNormalPerThread(thread_id); });
        break;
    case SET_VOLUME:
        k = job([func, thread_id]() {func->setVolume(thread_id); });
        break;
    case VERTEX_NORMAL_FROM_RENDER:
        k = job([func, thread_id]() {func->getVertexNormalFromRenderPerThread(thread_id); });
        break;
    case TET_NEIGHBOR_TET_VERTEX_INDEX:
        k = job([func, thread_id]() {func->updateTetNeighborTetVertexIndex(thread_id); });
        break;
    case SORT_TRIANGLE_EDGE_AROUND_TRIANGLE_EDGE:
        k = job([func, thread_id]() {func->setFaceEdgeAroundFace(thread_id);
        func->setTetAroundFace(thread_id); });
        break;
    case SORT_TRIANGLE_AROUND_VERTEX_EDGE:
        k = job([func, thread_id]() {func->sortTriangleAroundVertexEdge(thread_id);
        func->sortTetAroundVertexEdge(thread_id); });
        break;
    case TET_AROUND_TET_COLOR_GROUP:
        k = job([func, thread_id]() {func->setTetAroundTetColor(thread_id); });
        break;
    }
    return k;
}

job Thread::create_task(Cloth* func, int thread_id, ObjectFunc function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case EDGE_TRIANGLE_AABB:
        k = job([func, thread_id]() {func->getEdgeTriangleAABBPerThread(thread_id); });
        break;
        //case TRIANGLE_AABB:
        //    
        //    break;
        //case EDGE_AABB:
        //    k = job([func, thread_id]() {func->getEdgeAABBPerThread(thread_id); });
        //    break;
    case VERTEX_AABB:
        k = job([func, thread_id]() {func->getVertexAABBPerThread(thread_id, true); });
        break;
    case VERTEX_AABB_WITHOUT_TOLERANCE:
        k = job([func, thread_id]() {func->getVertexAABBPerThread(thread_id, false); });
        break;
    case CURRENT_AABB:
        k = job([func, thread_id]() {func->getCurrentPosAABB(thread_id); });
        break;
    }
    return k;
}


job Thread::create_task(Tetrahedron* func, int thread_id, ObjectFunc function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case EDGE_TRIANGLE_AABB:
        k = job([func, thread_id]() {func->getEdgeTriangleAABBPerThread(thread_id); });
        break;
        //case TRIANGLE_AABB:
        //    
        //    break;
        //case EDGE_AABB:
        //    k = job([func, thread_id]() {func->getEdgeAABBPerThread(thread_id); });
        //    break;
    case VERTEX_AABB:
        k = job([func, thread_id]() {func->getVertexAABBPerThread(thread_id, true); });
        break;
    case VERTEX_AABB_WITHOUT_TOLERANCE:
        k = job([func, thread_id]() {func->getVertexAABBPerThread(thread_id, false); });
        break;
    case TETRAHEDRON_AABB:
        k = job([func, thread_id]() {func->getTetAABBPerThread(thread_id); });
        break;
    case CURRENT_AABB:
        k = job([func, thread_id]() {func->getCurrentPosAABB(thread_id); });
        break;
    case FIND_NEIGHBOR_VERTEX:
        k = job([func, thread_id]() {func->findAllNeighborVertex(thread_id); });
        break;
    }
    return k;
}


job Thread::create_task(Collider* func, int thread_id, ObjectFunc function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
    case EDGE_TRIANGLE_AABB:
        k = job([func, thread_id]() {func->getEdgeTriangleAABBPerThread(thread_id); });
        break;
        //case TRIANGLE_AABB:
        //    k = job([func, thread_id]() {func->getTriangleAABBPerThread(thread_id); });
        //    break;
    case VERTEX_AABB:
        k = job([func, thread_id]() {func->getVertexAABBPerThread(thread_id, true); });
        break;
    case VERTEX_AABB_WITHOUT_TOLERANCE:
        k = job([func, thread_id]() {func->getVertexAABBPerThread(thread_id, false); });
        break;
    case CURRENT_AABB:
        k = job([func, thread_id]() {func->getCurrentPosAABB(thread_id); });
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
    case UPDATE_NODE_VALUE:
        k = job([func, thread_id]() {func->updateNodeValue(thread_id); });
        break;
    case UPDATE_LAST_LAYER_NODE_VALUE:
        k = job([func, thread_id]() {func->updateNodeValueLastLayer(thread_id); });
        break;
    }

    return k;
}


job Thread::create_task(Collision* func, int thread_id, CollisionFuncSendToThread function_type, int para)
{
    job k;
    switch (function_type)
    {
    case UPDATE_COLLISION_HESSIAN_COLOR:
        k = job([func, thread_id, para]() {func->computeHessianPerThread(para,thread_id); });
        break;      
    case COLOR_COLLISION_TIME:
        k = job([func, thread_id, para]() {func->colorCollisionTime(thread_id, para); });
        break;
    case UPDATE_COLOR_POSITION:
        k = job([func, thread_id, para]() {func->updatePositionColor(thread_id, para); });
        break;
    case SUM_COLLISION_HESSIAN:
        k = job([func, thread_id, para]() {func->sumAllCollisionHessian(para, thread_id); });
        break;
    //case FLOOR_COLLISION_TIME:
    //    k = job([func, thread_id, para]() {func->floorCollisionTime(thread_id, para); });
    //    break;
    //case UPDATE_POSITION_FOR_FLOOR_COLLISION:
    //    k = job([func, thread_id, para]() {func->updatePositionForFloor(thread_id, para); });
    //    break;
    }
    return k;
}



job Thread::create_task(Collision* func, int thread_id, CollisionFuncSendToThread function_type)// int jobNumber
{
    job k;
    switch (function_type)
    {
        //case FIND_PATCH_PAIRS:
        //    k = job([func, thread_id]() {func->findAllPatchPairs(thread_id); });
        //    break;  
    case GLOBAL_COLLISION_TIME_ADD_PAIR:
        k = job([func, thread_id]() {func->collisionTimeWithPair(thread_id); });
        break;
    case FIND_COLLISION_PAIR:
        k = job([func, thread_id]() {func->getCollisionPair(thread_id); });
        break;
    case  FIND_PRIMITIVE_AROUND:
        k = job([func, thread_id]() {func->findPrimitivesAround(thread_id); });
        break;
    case GLOBAL_COLLISION_TIME:
        k = job([func, thread_id]() {func->collisionTime(thread_id); });
        break;
    case COMPUTE_VOLUME:
        k = job([func, thread_id]() {func->computeVolume(thread_id); });
        break;
    case FIND_VERTEX_VERTEX_VERTEX_EDGE_PAIRS:
        k = job([func, thread_id]() {func->findAllVertexVertexEdgePairs(thread_id); });
        break;
    case COMPUTE_COLLISION_ENERGY:
        k = job([func, thread_id]() {func->collisionEnergy(thread_id); });
        break;
    case COLLISION_CONSTRAINT_IPC:
        k = job([func, thread_id]() {func->collisionConstraintIPC(thread_id); });
        break;
    case RE_COLLISION_CONSTRAINT_IPC:
        k = job([func, thread_id]() {func->re_collisionConstraintIPC(thread_id); });
        break;
    case COLLISION_CONSTRAINT:
        k = job([func, thread_id]() {func->collisionConstraint(thread_id); });
        break;
    case RE_DETECTION:
        k = job([func, thread_id]() {func->collisionReDetection(thread_id); });
        break;
    case RESUM_TARGET_POSITION:
        k = job([func, thread_id]() {func->resumTargetPositionPerThread(thread_id); });
        break;
    case GLOBAL_COLLISION_DETECTION:
        k = job([func, thread_id]() {func->collisionDetection(thread_id); });
        break;
    case SUM_TARGET_POSITION:
        k = job([func, thread_id]() {func->sumTargetPositionPerThread(thread_id); });
        break;
    case FIND_TRIANGLE_PAIRS:
        k = job([func, thread_id]() {func->findAllTrianglePairs(thread_id); });
        break;
    case RE_COLLISION_CONSTRAINT:
        k = job([func, thread_id]() {func->re_collisionConstraint(thread_id); });
        break;
    case PREFIX_SUM_ALL_PAIRS:
        k = job([func, thread_id]() {func->prefixSumAllPair(thread_id); });
        break;
    case SET_ELEMENT_COLLIDE_WITH_COLLIDER:
        k = job([func, thread_id]() {func->setElementCollideWithCollider(thread_id); });
        break;
    case UPDATE_RECORD_VERTEX_POSITION:
        k = job([func, thread_id]() {func->updateVertexRecordForColor(thread_id); });
        break;
    //case RECORD_VT_PAIR_COMPRESS:
    //    k = job([func, thread_id]() {func->recordVTPairCompress(thread_id); });
    //    break;
    //case RECORD_EE_PAIR_COMPRESS:
    //    k = job([func, thread_id]() {func->recordEEPairCompress(thread_id); });
    //    break;
    //case RECORD_TRIANGLE_HAS_COLLISION_PAIR:
    //    k = job([func, thread_id]() {func->recordTriangleHasTVPair(thread_id); });
    //    break;
    //case CLOSE_PAIR_COLLISION_TIME:
    //    k = job([func, thread_id]() {func->collisionTimeAllClosePair(thread_id); });
    //    break;
    case COLLISION_FREE_POSITION_LAST_COLOR:
        k = job([func, thread_id]() {func->computeCollisionFreePositionForColor(thread_id); });
        break;
    case FIND_CLOSE_PAIR:
        k = job([func, thread_id]() {func->findClosePair(thread_id); });
        break;
    case COMPUTE_D_HAT_WITH_INDEX_IN_GLOBAL:
        k = job([func, thread_id]() {func->computeDhatWithIndexInGlobal(thread_id); });
        break;
    case EXTRACT_ELEMENTS_WITH_COLLIDER:
        k = job([func, thread_id]() {func->extractElementCollideWithCollider(thread_id); });
        break;
    case ADD_TET_IN_COLLISION:
        k = job([func, thread_id]() {func->addTetInvolvedInCollision(thread_id); });
        break;
    case SET_PAIR_BY_ELEMENT:
        k = job([func, thread_id]() {func->setPairByElement(thread_id); });
        break;
    }
    return k;
}



job Thread::create_task(SpatialHashing* func, int thread_id, SpatialHashingFuncSendToThread function_type)//
{
    job k;
    switch (function_type)
    {
    case SCENE_AABB:
        k = job([func, thread_id]() {func->getSceneAABB(thread_id); });
        break;
    case TRIANGLE_HASHING_SMALLER_HASH_TABLE:
        k = job([func, thread_id]() {func->triangleHashingSmallerHashTable(thread_id); });
        break;
    case RECORD_NONEMPTY_CELL:
        k = job([func, thread_id]() {func->recordNonEmptyCell(thread_id); });
        break;
    case OBTAIN_PAIR_COUNT:
        k = job([func, thread_id]() {func->obtainPairCount(thread_id); });
        break;
    case SET_HASH_CELL_PAIR_NUM_PREFIX_SUM_TOGETHER:
        k = job([func, thread_id]() {func->setHashCellPairNumPrefixSumTogether(thread_id); });
        break;
    case FIND_ALL_PAIRS_HASH_TABLE:
        k = job([func, thread_id]() {func->findAllPairsHashTable(thread_id); });
        break;
    case SET_PAIR_AVE:
        k = job([func, thread_id]() {func->setPairAveInThread(thread_id); });
        break;
    case FIND_ALL_TRIANGLE_PAIRS_HASH_TABLE_ELEMENTWISE:
        k = job([func, thread_id]() {func->findAllPairsHashTableElementwise(thread_id); });
        break;
    case ORI_TRIANGLE_HASHING:
        k = job([func, thread_id]() {func->oriTriangleHashing(thread_id); });
        break;
    case COMBINE_HASH_TABLE:
        k = job([func, thread_id]() {func->combineHashTable(thread_id); });
        break;
    case FIND_ALL_PAIRS_HASH_TABLE_BY_ELEMENT:
        k = job([func, thread_id]() {func->findAllPairsByPrimitive(thread_id); });
        break;
        //case SET_HASH_TOGETHER:
        //    k = job([func, thread_id]() {func->setHashTogether(thread_id); });
        //    break;
        //case PREPARE_FOR_ACTUAL_HASH_VALUE_COUNT_THREAD: {
        //    k = job([func, thread_id]() {func->prepareForActualHashValueCountThread(thread_id); });
        //    break;
        //}
        //case ADD_COUNT_FOR_PRIFIX_SUM: {
        //    k = job([func, thread_id]() {func->prifixSum1(thread_id); });
        //    break;
        //}
        //case PREFIX_SUM_THREAD_1:
        //    k = job([func, thread_id]() {func->prifixSum2(thread_id); });
        //    break;
        //case PREFIX_SUM_THREAD_2:
        //    k = job([func, thread_id]() {func->prifixSum3(thread_id); });
        //    break;
        //case MEMSET_PREFIX:
        //    k = job([func, thread_id]() {func->memsetThread(thread_id); });
        //    break;
    }
    return k;
}


job Thread::create_task(SpatialHashing* func, int thread_id, SpatialHashingFuncSendToThread function_type, unsigned int key_id)
{
    job k;
    switch (function_type)
    {
    case 0:

        break;
        //case PREFIX_SUM_UP:
        //    k = job([func, thread_id, key_id]() {func->prefixSumParallelUp(thread_id, key_id); });
            //break;
        //case PREFIX_SUM_DOWN:
        //    k = job([func, thread_id, key_id]() {func->prefixSumParallelDown(thread_id, key_id); });
            //break;
    }
    return k;
}

job Thread::create_task(RadixSort* func, int thread_id, RadixSortFunc function_type, unsigned int key_id)
{
    job k;
    switch (function_type)
    {
    case SET_COUNT_BUCKET:
        k = job([func, thread_id, key_id]() {func->setCountBucket(thread_id, key_id); });
        break;
    case REORDER:
        k = job([func, thread_id, key_id]() {func->reorder(thread_id, key_id); });
        break;
    case MORTON_REORDER:
        k = job([func, thread_id, key_id]() {func->reorderMorton(thread_id, key_id); });
        break;
    case SET_COUNT_BUCKET_MORTON:
        k = job([func, thread_id, key_id]() {func->setCountBucketMorton(thread_id, key_id); });
        break;
    case COPY_ARRAY:
        k = job([func, thread_id]() {func->copyArray(thread_id); });
        break;
    }

    return k;
}

job Thread::create_task(SecondOrderLargeSystem* func, int thread_id, NewtonMethodFunc function_type)
{
    job k;
    switch (function_type)
    {
    case UPDATE_HESSIAN_FIXED_STRUCTURE:
        k = job([func, thread_id]() {func->updateHessianFixedStructure(thread_id); });
        break;
    case UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE:
        k = job([func, thread_id]() {func->setHessianDiagonalFixedStructure(thread_id); });
        break;
    case UPDATE_INTERNAL_FORCE:
        k = job([func, thread_id]() {func->updateInternalForce(thread_id); });
        break;
    case SUM_B:
        k = job([func, thread_id]() {func->sumB(thread_id); });
        break;
    case SET_S_N:
        k = job([func, thread_id]() {func->setSn(thread_id); });
        break;
    case UPDATE_POSITION_NEWTON:
        k = job([func, thread_id]() {func->updatePosition(thread_id); });
        break;
    case UPDATE_POSITION_NEWTON_FROM_ORI:
        k = job([func, thread_id]() {func->updatePositionFromOri(thread_id); });
        break;
    case VELOCITY_NEWTON:
        k = job([func, thread_id]() {func->updateVelocity(thread_id); });
        break;
    case NEWTON_METHOD_ENERGY:
        k = job([func, thread_id]() {func->computeEnergy(thread_id); });
        break;
    case SET_MASS_SPRING:
        k = job([func, thread_id]() {func->massSpring(thread_id); });
        break;
    case SET_HESSIAN_DIAGONAL:
        k = job([func, thread_id]() {func->setHessianDiagonal(thread_id); });
        break;
    case GET_COEFF_ADDRESS:
        k = job([func, thread_id]() {func->hessianCoeffAddress(thread_id); });
        break;
    case UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE_INITIAL_STIFFNESS:
        k = job([func, thread_id]() {func->setHessianDiagonalFixedStructureInitialStiffness(thread_id); });
        break;
    }
    return k;
}

job Thread::create_task(NewtonMethod* func, int thread_id, NewtonMethodFunc function_type)
{
    job k;
    switch (function_type)
    {
    case UPDATE_HESSIAN_FIXED_STRUCTURE:
        k = job([func, thread_id]() {func->updateHessianFixedStructure(thread_id); });
        break;
    case UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE:
        k = job([func, thread_id]() {func->setHessianDiagonalFixedStructure(thread_id); });
        break;
    case UPDATE_INTERNAL_FORCE:
        k = job([func, thread_id]() {func->updateInternalForce(thread_id); });
        break;
    case SUM_B:
        k = job([func, thread_id]() {func->sumB(thread_id); });
        break;
    case SET_S_N:
        k = job([func, thread_id]() {func->setSn(thread_id); });
        break;
    case SET_B_N:
        k = job([func, thread_id]() {func->setBn(thread_id); });
        break;
    case UPDATE_ANCHOR_POINT_HESSIAN:
        k = job([func, thread_id]() {func->updateHessianForFixPoint(thread_id); });
        break;
    case UPDATE_POSITION_NEWTON:
        k = job([func, thread_id]() {func->updatePosition(thread_id); });
        break;
    case UPDATE_POSITION_NEWTON_FROM_ORI:
        k = job([func, thread_id]() {func->updatePositionFromOri(thread_id); });
        break;
    case VELOCITY_NEWTON:
        k = job([func, thread_id]() {func->updateVelocity(thread_id); });
        break;
    case VELOCITY_NEWTON_2:
        k = job([func, thread_id]() {func->updateVelocity2(thread_id); });
        break;
    case UPDATEVELOCITY_ACCELERATION_NEWMARK:
        k = job([func, thread_id]() {func->updateVelocityAccelerationNewMark(thread_id); });
        break;
    case NEWTON_METHOD_ENERGY:
        k = job([func, thread_id]() {func->computeEnergy(thread_id); });
        break;
    case UPDATE_DAMP:
        k = job([func, thread_id]() {func->updateHessianForDamp(thread_id); });
        break;
    case SET_MASS_SPRING:
        k = job([func, thread_id]() {func->massSpring(thread_id); });
        break;
    case SET_HESSIAN_DIAGONAL:
        k = job([func, thread_id]() {func->setHessianDiagonal(thread_id); });
        break;
    case GET_COEFF_ADDRESS:
        k = job([func, thread_id]() {func->hessianCoeffAddress(thread_id); });
        break;
    case UPDATE_DIAGONAL_HESSIAN_FIXED_STRUCTURE_INITIAL_STIFFNESS:
        k = job([func, thread_id]() {func->setHessianDiagonalFixedStructureInitialStiffness(thread_id); });
        break;
    }
    return k;
}


job Thread::create_task(IterationMethod* func, int thread_id, IterationMethodFunc function_type)
{
    job k;
    switch (function_type)
    {
    case UPDATE_JACOBI_OPERATOR:
        k = job([func, thread_id]() {func->updateJacobiOperator(thread_id); });
        break;
    case UPDATE_2_A_JACOBI_ITR_MATRIX:
        k = job([func, thread_id]() {func->update2AJaocbiIterationMatrix(thread_id); });
        break;
    case UPDATE_3_A_JACOBI_ITR_MATRIX:
        k = job([func, thread_id]() {func->update3AJaocbiIterationMatrix(thread_id); });
        break;
    }
    return k;
}


job Thread::create_task(IterationMethod* func, int thread_id, IterationMethodFunc function_type, Eigen::VectorXd* u, Eigen::VectorXd* b, double* residual_norm,
    double omega_chebyshev, Eigen::VectorXd* u_last, Eigen::VectorXd* u_previous)
{
    job k;
    switch (function_type)
    {
    case JACOBI_ITR:
        k = job([func, thread_id, u, b, residual_norm]() {func->JacobiIterationPerThread(thread_id, u, b, residual_norm); });
        break;
    case A_JACOBI_2_ITR:
        k = job([func, thread_id, u, b, residual_norm]() {func->SuperJacobi2IterationPerThread(thread_id, u, b, residual_norm); });
        break;
    case A_JACOBI_3_ITR:
        k = job([func, thread_id, u, b, residual_norm]() {func->SuperJacobi3IterationPerThread(thread_id, u, b, residual_norm); });
        break;
    case CHEBYSHEV_JACOBI_ITR:
        k = job([func, thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous]()
            {func->ChebyshevSemiIterativeJacobiIterationPerThread(thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous); });
        break;
    case CHEBYSHEV_A_JACOBI_2_ITR:
        k = job([func, thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous]()
            {func->ChebyshevSemiIterativeAJacobi2IterationPerThread(thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous); });
        break;
    case CHEBYSHEV_A_JACOBI_3_ITR:
        k = job([func, thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous]()
            {func->ChebyshevSemiIterativeAJacobi3IterationPerThread(thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous); });
        break;
    case PCG_ITR1:
        k = job([func, thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous]()
            {func->PCGIterationPerThread1(thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous); });
        break;
    case PCG_ITR2:
        k = job([func, thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous]()
            {func->PCGIterationPerThread2(thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous); });
        break;
    case GAUSS_SEIDEL_ITR:
        k = job([func, thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous]()
            {func->GaussSeidelIterationPerThread(thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous); });
        break;
    case CHEBYSHEV_GAUSS_SEIDEL_ITR:
        k = job([func, thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous]()
            {func->ChebyshevSemiIterativeGaussSeidelIterationPerThread(thread_id, u, b, residual_norm, omega_chebyshev, u_last, u_previous); });
        break;
    }

    return k;
}


job Thread::create_task(IterationMethod* func, int thread_id, std::vector<int>* vertex_index, std::vector<double>* coefficient,
    double* x, double* b, double* result, int* vertex_index_thread_begin, int sys_size)
{
    job k;
    k = job([func, thread_id, vertex_index, coefficient, x, b, result, vertex_index_thread_begin, sys_size]()
        {func->RMultiXPlusb(vertex_index, coefficient, x, b, result, vertex_index_thread_begin[thread_id],
            vertex_index_thread_begin[thread_id + 1], sys_size); });
    return k;
}


job Thread::create_task(IterationMethod* func, IterationMethodFunc function_type, int thread_id, int* vertex_index, double* coefficient, int* vertex_index_start,
    Eigen::VectorXd* x, Eigen::VectorXd* b, Eigen::VectorXd* result, double* residual_norm, int* vertex_index_thread_begin,
    Eigen::VectorXd* u_last, Eigen::VectorXd* u_previous)
{
    job k;
    switch (function_type)
    {
    case R_MULTIPLY_X_PLUS_B:
        k = job([func, thread_id, vertex_index, coefficient, vertex_index_start, x, b, result, vertex_index_thread_begin]()
            {func->RMultiXPlusb(vertex_index, coefficient, vertex_index_start, x, b, result, vertex_index_thread_begin[thread_id],
                vertex_index_thread_begin[thread_id + 1]); });
        break;
    case W_R_MULTIPLY_X_PLUS_B_1_W_X:
        k = job([func, thread_id, vertex_index, coefficient, vertex_index_start, x, b, result, vertex_index_thread_begin]()
            {func->A_jacobi_RMultiXPlusb(vertex_index, coefficient, vertex_index_start, x, b, result, vertex_index_thread_begin[thread_id],
                vertex_index_thread_begin[thread_id + 1]); });
        break;
    case COMPUTE_RESIDUAL:
        k = job([func, thread_id, vertex_index, coefficient, vertex_index_start, x, b, result, residual_norm, vertex_index_thread_begin]()
            {func->computeResidual(vertex_index, coefficient, vertex_index_start, x, b, result, &residual_norm[thread_id], vertex_index_thread_begin[thread_id],
                vertex_index_thread_begin[thread_id + 1]); });
        break;
    case CHEBYSHEV_A_JACOBI_ITERATION: {
        k = job([func, thread_id, vertex_index, coefficient, vertex_index_start, x, b, result, residual_norm, vertex_index_thread_begin,
            u_last, u_previous]()
            {func->ChebyshevAJacobiIterationPerThread(vertex_index, coefficient, vertex_index_start, x, b, result, residual_norm, vertex_index_thread_begin[thread_id],
                vertex_index_thread_begin[thread_id + 1], u_last, u_previous); });
        break;
    }
    case ESTIMATE_A_JACOBI_2_EIGEN_VALUE:
        k = job([func, thread_id, coefficient, residual_norm, x, vertex_index_thread_begin]()
            {func->estimateAJacobi2EigenValue(x, &coefficient[thread_id], &residual_norm[thread_id], vertex_index_thread_begin[thread_id],
                vertex_index_thread_begin[thread_id + 1]); });
        break;
    case ESTIMATE_A_JACOBI_3_EIGEN_VALUE:
        k = job([func, thread_id, coefficient, residual_norm, x, vertex_index_thread_begin]()
            {func->estimateAJacobi3EigenValue(x, &coefficient[thread_id], &residual_norm[thread_id], vertex_index_thread_begin[thread_id],
                vertex_index_thread_begin[thread_id + 1]); });
        break;
    }
    return k;
}


job Thread::create_task(MoveObject* func, int thread_id, MoveObjectFunc function_type, unsigned int key_id)
{
    job k;
    switch (function_type)
    {
    case MOVE_OBJECT:
        k = job([func, thread_id, key_id]() {func->move(thread_id, key_id); });
        break;
    case MOVE_OBJECT2:
        k = job([func, thread_id, key_id]() {func->moveDiffInitialCurrent(thread_id, key_id); });
        break;
    case ROTATE_AROUND_AXIS:
        k = job([func, thread_id, key_id]() {func->rotateAroundAxis(thread_id, key_id); });
        break;
    }
    return k;
}

//job Thread::create_task(DrawCulling* func, int thread_id, DrawCullingFunc function_type, unsigned int key_id)
//{
//    job k;
//    switch (function_type)
//    {
//    case SET_POSITION_COLOR:
//        k = job([func, thread_id, key_id]() {func->setAllTriangle(thread_id); });
//        break;
//    case SET_DATA_TOGETHER:
//        k = job([func, thread_id, key_id]() {func->setThreadDataTogether(thread_id); });
//        break;
//    }
//    return k;
//}

job Thread::create_task(XPBD* func, int thread_id, XPBDFunc function_type)
{
    job k;
    switch (function_type)
    {
    case SET_POS_PREDICT:
        k = job([func, thread_id]() {func->setPosPredict(thread_id); });
        break;
    case SET_POS_PREDICT_SUB_TIME_STEP:
        k = job([func, thread_id]() {func->setPosPredictSubTimeStep(thread_id,false); });
        break;
    case SET_POS_PREDICT_SUB_TIME_STEP_FOR_CULLING:
        k = job([func, thread_id]() {func->setPosPredictSubTimeStep(thread_id, true); });
        break;
    case XPBD_VELOCITY:
        k = job([func, thread_id]() {func->computeVelocity(thread_id); });
        break;
    }
    return k;
}



job Thread::create_task(XPBD_IPC* func, int thread_id, XPBD_IPC_Func function_type, unsigned int para)
{
    job k;
    switch (function_type)
    {
    case SOLVE_TET_BLOCK:
        //k = job([func, thread_id, para]() {func->newtonCDTetBlockAGroupTest(thread_id, para);});
        k = job([func, thread_id, para]() {func->newtonCDTetBlockAGroup(thread_id, para); });
         break;
    case FIRST_COLOR_ARAP_ENERGY:
        k = job([func, thread_id, para]() {func->computePreviousColorARAPEnergy(thread_id, para); });
        break;
    case PREVIOUS_COLOR_INERTIAL_ENERGY:
        k = job([func, thread_id, para]() {func->computePreviousColorInertialEnergy(thread_id, para); });
        break;
    case UPDATE_TET_GRAD_SHARED:
        k = job([func, thread_id, para]() {func->tetGradForColor(thread_id, para); });
        break;
    //case UPDATE_TET_GRAD_SHARED_COLLISION:
    //    k = job([func, thread_id, para]() {func->tetGradForColorCollision(thread_id, para); });
    //    break;
    //case SOLVE_TET_BLOCK_COLLISION:
    //    k = job([func, thread_id, para]() {func->newtonCDTetBlockAGroupCollision(thread_id, para); });
    //    break;
    case UPDATE_TET_GRAD_SHARED_COLLISION_NEIGHBOR:
        k = job([func, thread_id, para]() {func->tetGradForColorCollisionNeighbor(thread_id, para); });
        break;
    }
    return k;
}

job Thread::create_task(XPBD_IPC* func, int thread_id, XPBD_IPC_Func function_type)
{
    job k;
    switch (function_type)
    {
    case SET_POS_PREDICT_:
        k = job([func, thread_id]() {func->setPosPredict(thread_id); });
        break;
    case COLLISION_FREE_POSITION_:
        k = job([func, thread_id]() {func->computeCollisionFreePosition(thread_id); });
        break;
    case XPBD_IPC_VELOCITY:
        k = job([func, thread_id]() {func->computeVelocity(thread_id); });
        break;
    case SUM_ALL_GRAD:
        k = job([func, thread_id]() {func->sumAllGrad(thread_id); });
        break;
    //case UPDATE_TET_HESSIAN:
    //    k = job([func, thread_id]() {func->tetHessian(thread_id); });
    //    break;
    case UPDATE_POSITION_AVERAGE:
        k = job([func, thread_id]() {func->updatePositionAverage(thread_id); });
       break;
    case UPDATE_LAST_COLOR_VERTEX_BELONG:
        k = job([func, thread_id]() {func->lastColorVertexBelongToGroup(thread_id); });
        break;
    case COLLISION_FREE_POSITION_FROM_RECORD:
        k = job([func, thread_id]() {func->computeCollisionFreePositionFromRecord(thread_id); });
        break;
    case INERTIAL_ENERGY:
        k = job([func, thread_id]() {func->inertialEnergyPerThread(thread_id); });
        break;
    case ARAP_ENERGY:
        k = job([func, thread_id]() {func->computeARAPEnergyPerThread(thread_id); });
        break;
    case COMPUTE_BARRIER_ENERGY:
        k = job([func, thread_id]() {func->computeBarrierEnergy(thread_id); });
        break;
    case COMPUTE_RPREVIOUS_COLOR_BARRIER_ENERGY:
        k = job([func, thread_id]() {func->computePreviousColorCollisionEnergy(thread_id); });
        break;
    case LAST_COLOR_INERTIAL_ENERGY:
        k = job([func, thread_id]() {func->computeColorInertialEnergy(thread_id); });
        break;
    case LAST_COLOR_ARAP_ENERGY:
        k = job([func, thread_id]() {func->computeLastColorARAPEnergy(thread_id); });
        break;
    }
    return k;
}

//job Thread::create_task(MeshPatch* func, int thread_id, MeshPatchFunc function_type)
//{
//    job k;
//    switch (function_type)
//    {  
//    case PATCH_AABB:
//        k = job([func, thread_id]() {func->obtainAABB(thread_id); });
//        break;
//    case DECIDE_TRIANGLE_INDEX_SIZE:
//        k = job([func, thread_id]() {func->decideTriangleIndexSize(thread_id); });
//        break;
//    case FIND_TRIANGLE_INDEX:
//        k = job([func, thread_id]() {func->findNotedTrueTriangleIndex(thread_id); });
//        break;
//    case FIND_VERTEX:
//        k = job([func, thread_id]() {func->findVertex(thread_id); });
//        break;   
//    }
//    return k;
//}

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

void Thread::assignTask(IterationMethod* func, std::vector<int>* vertex_index, std::vector<double>* coefficient,
    double* x, double* b, double* result, int* vertex_index_thread_begin, int sys_size)
{
    for (int i = 0; i < thread_num; ++i)
    {
        // std::cout << threads[i].id << std::endl;
        job j = create_task(func, threads[i].id, vertex_index, coefficient, x, b, result, vertex_index_thread_begin, sys_size);
        futures.push_back(j.get_future());
        std::unique_lock<std::mutex> l(threads[i].m);
        threads[i].jobs.push(std::move(j));
        // Notify the thread that there is work do to...
        threads[i].cv.notify_one();
    }
    for (auto& f : futures) { f.wait(); }
    futures.clear();
}

void Thread::assignTask(IterationMethod* func, IterationMethodFunc function_type, Eigen::VectorXd* u, Eigen::VectorXd* b,
    double* residual_norm, double omega_chebyshev, Eigen::VectorXd* u_last, Eigen::VectorXd* u_previous)
{
    for (int i = 0; i < thread_num; ++i)
    {
        // std::cout << threads[i].id << std::endl;
        job j = create_task(func, threads[i].id, function_type, u, b, residual_norm, omega_chebyshev, u_last, u_previous);
        futures.push_back(j.get_future());
        std::unique_lock<std::mutex> l(threads[i].m);
        threads[i].jobs.push(std::move(j));
        // Notify the thread that there is work do to...
        threads[i].cv.notify_one();
    }
    for (auto& f : futures) { f.wait(); }
    futures.clear();
}

void Thread::assignTask(IterationMethod* func, IterationMethodFunc function_type, int* vertex_index, double* coefficient, int* vertex_index_start,
    Eigen::VectorXd* x, Eigen::VectorXd* b, Eigen::VectorXd* result, double* residual_norm, int* vertex_index_thread_begin,
    Eigen::VectorXd* u_last, Eigen::VectorXd* u_previous)
{
    for (int i = 0; i < thread_num; ++i)
    {
        // std::cout << threads[i].id << std::endl;
        job j = create_task(func, function_type, threads[i].id, vertex_index, coefficient, vertex_index_start, x, b, result,
            residual_norm, vertex_index_thread_begin, u_last, u_previous);
        futures.push_back(j.get_future());
        std::unique_lock<std::mutex> l(threads[i].m);
        threads[i].jobs.push(std::move(j));
        // Notify the thread that there is work do to...
        threads[i].cv.notify_one();
    }
    for (auto& f : futures) { f.wait(); }
    futures.clear();
}