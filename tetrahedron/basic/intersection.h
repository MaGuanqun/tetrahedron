#ifndef INTERSECTION_H
#define INTERSECTION_H

class Intersection
{
public:
    Intersection() {
        happened = false;
        happened_include_collider = false;
        face_index = -1;
        obj_No = -1;
    }
    bool happened_include_collider;
    bool happened;
    int face_index;
    int obj_No;
    bool first_intersection_frame;

    void initialIntersection() {
        happened = false;
        happened_include_collider = false;
        face_index = -1;
        obj_No = -1;
        first_intersection_frame = false;
    }

    void setIntersection(int* triangle_index) {
        happened_include_collider = true;
        face_index = triangle_index[0];
        obj_No = triangle_index[1];
    }

};


#endif // !INTERSECTION_H

#pragma once
