#ifndef INTERSECTION_H
#define INTERSECTION_H

class Intersection
{
public:
	Intersection() {
        happened = false;
        face_index = -1;
       obj_No = -1;
	}

    bool happened;
    int face_index;
    int obj_No;
    bool is_cloth;
    void initialIntersection() {
        happened = false;
        face_index = -1;
        obj_No = -1;
        is_cloth = true;
    }

    void setIntersection(int* triangle_index, bool is_cloth) {
        happened = true;
        face_index = triangle_index[0];
        obj_No = triangle_index[1];
        this->is_cloth = is_cloth;
    }

};


#endif // !INTERSECTION_H

#pragma once
