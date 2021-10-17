#ifndef INTERSECTION_H
#define INTERSECTION_H

class Intersection
{
public:
	Intersection() {
        happened = false;
        face_index = -1;
        cloth_No = -1;
	}

    bool happened;
    int face_index;
    int cloth_No;

    void initialIntersection() {
        happened = false;
        face_index = -1;
        cloth_No = -1;
    }

    void setIntersection(int* triangle_index) {
        happened = true;
        face_index = triangle_index[0];
        cloth_No = triangle_index[1];
    }

};


#endif // !INTERSECTION_H

#pragma once
