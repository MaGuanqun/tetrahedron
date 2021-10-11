#ifndef INTERSECTION_H
#define INTERSECTION_H

class Intersection
{
public:
	Intersection() {
        happened = false;
        distance = DBL_MAX;
        baryCenCoord[0] = 0.0;
        baryCenCoord[1] = 0.0;
        face_index = -1;
        cloth_No = -1;

	}

    bool happened;
    double distance;
    int face_index;
    int cloth_No;
    std::vector<int>vertex_index;
    //Vector3d emission;
    double baryCenCoord[2];
    double position[3];

    void initialIntersection() {
        happened = false;
        distance = DBL_MAX;
        baryCenCoord[0] = 0.0;
        baryCenCoord[1] = 0.0;
        face_index = -1;
        cloth_No = -1;
    }

};


#endif // !INTERSECTION_H

#pragma once
