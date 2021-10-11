#ifndef WRITE_OBJ_H
#define WRITE_OBJ_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include <iomanip>
#include"../object/cloth.h"
#include <direct.h>
#include <io.h>

class WriteObj
{
public:
	void write(std::vector<OriMesh>& clothes, int precision, int time_stamp)
	{
		std::string prefix = "./obj_record/";
		if (_access(prefix.c_str(), 0) == -1)
			_mkdir(prefix.c_str());
		std::string name;
		std::string basic_name = prefix+"cloth";
		for (int i = 0; i < clothes.size(); ++i) {
			name = basic_name + std::to_string(i);
			name += "_";
			name+= std::to_string(time_stamp);
			writeToFile(clothes[i], name, precision);
		}
	}

	void write(std::vector<Cloth>& clothes, int precision, int time_stamp)
	{
		std::string prefix = "./obj_record/";
		if (_access(prefix.c_str(), 0) == -1)
			_mkdir(prefix.c_str());
		std::string name;
		std::string basic_name = prefix + "cloth";
		for (int i = 0; i < clothes.size(); ++i) {
			name = basic_name + std::to_string(i);
			name += "_";
			name += std::to_string(time_stamp);
			writeToFile(clothes[i], name, precision);
		}
	}

private:


	void writeToFile(Cloth& model, std::string& name, int precision)
	{
		std::ofstream input_file;
		input_file.precision(precision);
		std::string obj_name = name + ".obj";
		std::string material_name = name + ".mtl";
		input_file.open(obj_name.c_str(), std::ios::trunc);
		std::vector<std::array<double, 3>>* position = &model.mesh_struct.vertex_position;
		std::vector<int>* indices = &model.mesh_struct.triangle_indices;
		for (int i = 0; i < position->size(); ++i) {
			input_file << "v " << (*position)[i][0] << " " << (*position)[i][1] << " " << (*position)[i][2] << "\n";
		}	
		for (int j = 0; j < indices->size() / 3; ++j) {
			input_file << "f " << (*indices)[3 * j] + 1 << " " << (*indices)[3 * j + 1] + 1 << " " << (*indices)[3 * j + 2] + 1 << "\n";
		}
		input_file.close();
	}

	void writeToFile(OriMesh& model, std::string& name, int precision)
	{
		std::ofstream input_file;
		input_file.precision(precision);
		std::string obj_name = name + ".obj";
		std::string material_name = name + ".mtl";
		input_file.open(obj_name.c_str(), std::ios::trunc);
		input_file << "mtllib " << material_name << "\n";
		for (int i = 0; i < model.vertices.size(); ++i) {
			input_file << "v " << model.vertices[i][0] << " " << model.vertices[i][1] << " " << model.vertices[i][2] << "\n";
		}
		//for (int i = 0; i < model.texture.size(); ++i) {
		//	input_file << "vt " << model.texture[i][0] << " " << model.texture[i][1] <<"\n";
		//}
	

		input_file << "usemtl " << model.front_material.material_name << "\n";
		for (int j = 0; j < model.indices.size() / 3; ++j) {
			input_file << "f " << model.indices[3 * j]+1 << " " << model.indices[3 * j + 1]+1 << " " << model.indices[3 * j + 2]+1 << "\n";
		}
		
		input_file.close();
		input_file.open(material_name.c_str(), std::ios::trunc);
	
		input_file << "newmtl " << model.front_material.material_name << "\n";
		input_file << "illum " << model.front_material.illum << "\n";
		input_file << "Kd " << model.front_material.Kd[0] << " " << model.front_material.Kd[1] << " " << model.front_material.Kd[2] << "\n";
		input_file << "Ka " << model.front_material.Ka[0] << " " << model.front_material.Ka[1] << " " << model.front_material.Ka[2] << "\n";
		input_file << "Ks " << model.front_material.Ks[0] << " " << model.front_material.Ks[1] << " " << model.front_material.Ks[2] << "\n";
		input_file << "Tf " << model.front_material.Tf[0] << " " << model.front_material.Tf[1] << " " << model.front_material.Tf[2] << "\n";
		input_file << "Ni " << model.front_material.Ni << "\n";
		input_file << "Ns " << model.front_material.Ns << "\n";
		input_file << "\n";
		
		input_file.close();
	}

};


#endif // !WRITE_OBJ_H

#pragma once
