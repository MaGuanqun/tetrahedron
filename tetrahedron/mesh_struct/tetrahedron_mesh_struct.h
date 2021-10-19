#pragma once
#include"mesh_struct.h"
#include<map>
class TetrahedronMeshStruct:public MeshStruct
{
public:
	struct TetrahedronVertex
	{
		std::vector<int>face;
	};
	std::vector<double>mass;
	std::vector<TetrahedronVertex> vertices;
	void setVertex();
	void findSurface();
	void getRenderVertexNormalPerThread(int thread_id);
	void setThreadIndex(int total_thread_num_);
	void getRenderNormal();
	void getRenderFaceNormalPerThread(int thread_id);
	std::vector<std::array<int, 4>> indices;
	std::vector<double>volume;

	std::vector<int> tetrahedron_index_begin_per_thread;
	void setVolume(int thread_No);
	double setVolumeMass(double density);
	std::vector<bool>vertex_on_surface;
private:
	struct TetrahedronFace {
		std::array<int,3> index;
		int sorted_index[3];
		TetrahedronFace(int v0, int v1, int v2) {
			index.data()[0] = v0;
			index.data()[1] = v1;
			index.data()[2] = v2;
			memcpy(sorted_index, index.data(), 12);
			std::sort(sorted_index, sorted_index + 3);
		}
		bool operator<(const TetrahedronFace& t1) const
		{
			if (sorted_index[0] < t1.sorted_index[0])
				return true;
			else if (sorted_index[0] == t1.sorted_index[0]) {
				if (sorted_index[1] < t1.sorted_index[1])
					return true;
				else if (sorted_index[1] == t1.sorted_index[1]) {
					if (sorted_index[2] < t1.sorted_index[2])
						return true;
				}
			}
			return false;
		}
	};
	void buildMap(std::map<TetrahedronFace, int>& face_in_tet, int v0, int v1, int v2);
	double getTetrahedronVolume(double* v1, double* v2, double* v3, double* v4);
};

