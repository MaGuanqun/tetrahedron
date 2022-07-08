#ifndef WRITE_OBJ_H
#define WRITE_OBJ_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include <iomanip>
#include <direct.h>
#include <io.h>
#include<array>

class WriteObj
{
public:
	template <class T>
	void write(std::vector<std::array<T, 3>>& position, std::vector<std::array<int, 3>>& vertex_indices, int precision, int time_stamp, int cloth_index, std::string name_, MeshMaterial& material)
	{
		std::string prefix = "./obj_record/";
		if (_access(prefix.c_str(), 0) == -1)
			_mkdir(prefix.c_str());
		std::string name;
		std::string basic_name = prefix + name_;
		name = basic_name + std::to_string(cloth_index);
		name += "_";
		name += std::to_string(time_stamp);
		std::string mtl_name = prefix + name_ + std::to_string(cloth_index) +"_"+ std::to_string(time_stamp);
		writeToFile(mtl_name,position, vertex_indices, name, precision, material);
		
	}

private:
	template <class T>
	void writeToFile(std::string& mtl_name,  std::vector<std::array<T,3>>& position, std::vector<std::array<int, 3>>& vertex_indices,  std::string& name, int precision, MeshMaterial& material)
	{
		std::ofstream input_file;
		input_file.precision(precision);
		std::string obj_name = name + ".obj";
		std::string material_name = mtl_name + ".mtl";
		input_file.open(obj_name.c_str(), std::ios::trunc);
		input_file << "mtllib " << material_name << "\n";
		for (int i = 0; i < position.size(); ++i) {
			input_file << "v " << position[i][0] << " " << position[i][1] << " " << position[i][2] << "\n";
		}
		//for (int i = 0; i < model.texture.size(); ++i) {
		//	input_file << "vt " << model.texture[i][0] << " " << model.texture[i][1] <<"\n";
		//}
		//for (int i = 0; i < model.face_groups.size(); ++i) {
			//input_file << "g " << model.face_groups[i].group_name << "\n";
			input_file << "usemtl front_material" << "\n";
			for (int j = 0; j < vertex_indices.size(); ++j) {
				input_file << "f " << vertex_indices[j][0]+1 << " " << vertex_indices[j][1]+1 << " " << vertex_indices[j][2]+1 << "\n";
			}
		//}
		input_file.close();
		input_file.open(material_name.c_str(), std::ios::trunc);
		//for (int i = 0; i < model.face_groups.size(); ++i) {
		input_file << "newmtl front_material" << "\n";
		input_file << "illum " << material.illum << "\n";
		input_file << "Kd " << material.Kd[0] << " " << material.Kd[1] << " " << material.Kd[2] << "\n";
		input_file << "Ka " << material.Ka[0] << " " << material.Ka[1] << " " << material.Ka[2] << "\n";
		input_file << "Ks " << material.Ks[0] << " " << material.Ks[1] << " " << material.Ks[2] << "\n";
		input_file << "Tf " << material.Tf[0] << " " << material.Tf[1] << " " << material.Tf[2] << "\n";
		input_file << "Ni " << material.Ni << "\n";
		input_file << "Ns " << material.Ns << "\n";
		input_file << "\n";
		//}
		input_file.close();
	}

};


#endif // !WRITE_OBJ_H

#pragma once
