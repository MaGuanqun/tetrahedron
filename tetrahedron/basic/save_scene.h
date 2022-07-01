#pragma once
#include <io.h>
#include <direct.h>
#include <fstream>
#include <iostream>
#include<vector>
#include<array>
#include"../mesh_struct/mesh_struct.h"



class SaveScene
{
public:
	//save time stamp, simulate_scene_indicator, position velocity of all objects, position of all colliders
	void save_scene_XPBD(size_t time_stamp, unsigned int simulate_scene_indicator, std::vector<MeshStruct*>& obj_mesh_struct,
		std::vector<std::vector<std::array<double, 3>>>* velocity, std::vector<MeshStruct*>& collider_mesh_struct);

	bool read_scene_XPBD(const char* file_name, size_t* time_stamp, unsigned int* simulate_scene_indicator, std::vector<MeshStruct*>& obj_mesh_struct,
		std::vector<std::vector<std::array<double, 3>>>* velocity, std::vector<MeshStruct*>& collider_mesh_struct);

private:
	size_t record_time_stamp= ULLONG_MAX;
	template<typename T>
	void write_binary_trunc(const char* filename, T* data, unsigned int size)
	{
		std::ofstream input_file(filename, std::ios::out | std::ios::binary | std::ios::trunc);
		//input_file.precision(64);
		input_file.write((char*)&size, sizeof(unsigned int));
		input_file.write((char*)data, size * sizeof(T));
		input_file.close();
	}
	template<typename T>
	void write_binary_add(std::ofstream& input_file, const char* filename, T* data, unsigned int size)
	{
		//std::ofstream input_file(filename, std::ios::out | std::ios::binary | std::ios::app);
		//input_file.precision(64);
		input_file.write((char*)&size, sizeof(unsigned int));
		input_file.write((char*)data, size * sizeof(T));
		//input_file.close();
	}

	template<typename T>
	bool read_binary(std::ifstream& in, T* data, unsigned int ori_size)
	{
		unsigned int size=0;
		in.read((char*)(&size), sizeof(unsigned int));
		if (ori_size != size) {
			std::cout << ori_size << " " << size << std::endl;
			std::cout << "the vertex number does not match" << std::endl;
			return false;
		}
		in.read((char*)data, size * sizeof(T));
		return true;
	}
};
