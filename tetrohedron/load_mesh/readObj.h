#ifndef READFILE_H
#define READFILE_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include <iomanip>
#include<time.h>
#include"../basic/global.h"
#include<array>
#include"../global_struct.h"

class ReadObj {
public:
	void load(const char* name, OriMesh& mesh)
	{
		std::vector<std::vector<std::string>> readString;
		std::vector<std::string> temString;
		std::string mtl_file_name;
		loadFile(name, temString);
		splitting(temString, readString);
		transferToMesh(readString, mesh, mtl_file_name);
		if (!mtl_file_name.empty()) {
			std::string obj_path(name);
			std::string mtlpath;
			splitLastCharacter(obj_path, mtlpath);
			mtlpath += mtl_file_name;

			std::vector<std::string> tem_mtl_string;
			std::vector<std::vector<std::string>> read_mtl_string;
			loadFile(mtlpath.c_str(), tem_mtl_string);
			splitting(tem_mtl_string, read_mtl_string);
			transferMaterial(read_mtl_string, mesh);
		}
		//refineMesh(mesh);
		//time_t t2 = clock();
		//std::cout << "read file cost " << (double)(t2 - t1) / (double)CLOCKS_PER_SEC << std::endl;
		//test(mesh);
	}
private:

	void refineMesh(OriMesh& mesh)//let the center of mesh be origin point
	{
		double xMin, xMax, yMin, yMax, zMin, zMax;
		xMin = xMax = mesh.vertices[0][0];
		yMin = yMax = mesh.vertices[0][1];
		zMin = zMax = mesh.vertices[0][2];
		for (int i = 1; i < mesh.vertices.size(); ++i) {
			xMin = myMin(xMin, mesh.vertices[i][0]);
			xMax = myMax(xMax, mesh.vertices[i][0]);
			yMin = myMin(yMin, mesh.vertices[i][1]);
			yMax = myMax(yMax, mesh.vertices[i][1]);
			zMin = myMin(zMin, mesh.vertices[i][2]);
			zMax = myMax(zMax, mesh.vertices[i][2]);
		}
		double xMid = 0.5 * (xMin + xMax);
		double yMid = 0.5 * (yMin + yMax);
		double zMid = 0.5 * (zMin + zMax);
		for (int i = 0; i < mesh.vertices.size(); ++i) {
			mesh.vertices[i][0] -= xMid;
			mesh.vertices[i][1] -= yMid;
			mesh.vertices[i][2] -= zMid;
		}
	}

	void transferMaterial(std::vector<std::vector<std::string>>& result, OriMesh& mesh)
	{
		bool start_record = false;
		for (int i = 0; i < result.size(); i++) {
			if (result[i][0] == "newmtl") {
				if (result[i][1] == mesh.front_material.material_name) {
					start_record = true;
				}
				else {
					start_record = false;
				}
			}
			if (start_record) {
				if (result[i][0] == "illum") {
					mesh.front_material.illum = std::stoi(result[i][1]);
				}
				else if (result[i][0] == "Kd") {
					mesh.front_material.Kd[0] = std::stof(result[i][1]);
					mesh.front_material.Kd[1] = std::stof(result[i][2]);
					mesh.front_material.Kd[2] = std::stof(result[i][3]);
				}
				else if (result[i][0] == "Ka") {
					mesh.front_material.Ka[0] = std::stof(result[i][1]);
					mesh.front_material.Ka[1] = std::stof(result[i][2]);
					mesh.front_material.Ka[2] = std::stof(result[i][3]);
				}
				else if (result[i][0] == "Ks") {
					mesh.front_material.Ks[0] = std::stof(result[i][1]);
					mesh.front_material.Ks[1] = std::stof(result[i][2]);
					mesh.front_material.Ks[2] = std::stof(result[i][3]);
				}
				else if (result[i][0] == "Tf") {
					mesh.front_material.Tf[0] = std::stof(result[i][1]);
					mesh.front_material.Tf[1] = std::stof(result[i][2]);
					mesh.front_material.Tf[2] = std::stof(result[i][3]);
				}
				else if (result[i][0] == "Ni") {
					mesh.front_material.Ni = std::stof(result[i][1]);
				}
				else if (result[i][0] == "Ns") {
					mesh.front_material.Ns = std::stof(result[i][1]);
				}
			}
		}
	
	}

	void transferToMesh(std::vector<std::vector<std::string>>& result, OriMesh& mesh, std::string& mtl_file_name)
	{
		int g_;
		std::vector<int>temp_face_indices;
		temp_face_indices.resize(4);
		mesh.vertices.reserve(1000);
		mesh.texture.reserve(1000);
		for (int i = 0; i < result.size(); i++) {
			if (result[i][0] == "v") {
				mesh.vertices.push_back({ stod(result[i][1]),stod(result[i][2]),stod(result[i][3]) });
			}
			else if (result[i][0] == "vt") {
				mesh.texture.push_back({ stod(result[i][1]),stod(result[i][2]) });
			}
			else if (result[i][0] == "f") {
				std::vector<std::vector<std::string>> split_f(4);//to record the splitting of f

				switch (result[i].size())
				{
				case 4:
					for (int j = 0; j < 3; j++) {
						split(result[i][j + 1], "/", split_f[j]);
						mesh.indices.push_back(stoi(split_f[j][0]) - 1);
					}
					break;
				case 5:
					for (int j = 0; j < 4; j++) {
						split(result[i][j + 1], "/", split_f[j]);
						temp_face_indices[j] = stoi(split_f[j][0]) - 1;
					}
					mesh.indices.push_back(temp_face_indices[0]);
					mesh.indices.push_back(temp_face_indices[1]);
					mesh.indices.push_back(temp_face_indices[2]);
					mesh.indices.push_back(temp_face_indices[0]);
					mesh.indices.push_back(temp_face_indices[2]);
					mesh.indices.push_back(temp_face_indices[3]);

					break;
				}
			}
			else if (result[i][0] == "mtllib") {
				mtl_file_name = result[i][1];
			}
			else if (result[i][0] == "usemtl") {
				mesh.front_material.material_name = result[i][1];
			}
		}
	}

	void loadFile(const char* name, std::vector<std::string>& temString)
	{
		std::string line;
		std::ifstream fin;
		fin.open(name);
		if (!fin.is_open()) {
			std::cout << "error opening " << *name << std::endl;
		}
		else
		{
			while (std::getline(fin, line)) {

				//std::stringstream ss(line);
				//std::string token;
				//std::vector<std::string> lineString;
				//while (ss>>token)
				////while (std::getline(ss,token,'\t'))
				//{
				//	lineString.push_back(token);
				//}
				//readString.push_back(lineString);
				temString.push_back(line);
			}
			fin.close();
		}
	}

	void splitting(std::vector<std::string>& ori, std::vector<std::vector<std::string>>& result)
	{
		for (int i = 0; i < ori.size(); i++)
		{
			std::vector<std::string>tem;
			split(ori[i], " ", tem);
			result.push_back(tem);
		}
	}

	void test(OriMesh& mesh)
	{
		std::cout << std::endl;
		//std::cout << mesh.face_groups[0].group_name << " " << mesh.face_groups[0].front_material.Kd[0] << " " << mesh.face_groups[0].front_material.Kd[1] << " " << mesh.face_groups[0].front_material.Kd[2] << " " << mesh.face_groups[0].front_material.Ns;
		std::cout << std::endl;
		/*		for (int i = 0; i < mesh.face_groups.size(); i++) {
					for(int j=0;j<20;j++)
						std::cout << mesh.face_groups[i].indices[j] << " ";
					std::cout << std::endl;
				}	*/
	}

	void split(std::string& str, const std::string& delim, std::vector<std::string>& dest)
	{
		std::string::size_type start = 0, index;
		if (str.find_last_of(delim) == str.length() - 1) {
			while (start != std::string::npos)
			{
				index = str.find_first_of(delim, start);
				if (index == std::string::npos) {
					dest.push_back(str.substr(start, str.length() - 1 - start));
					return;
				}
				else {
					dest.push_back(str.substr(start, index - start));
				}
				start = str.find_first_not_of(delim, index);
			}
		}
		else {
			while (start != std::string::npos)
			{
				index = str.find_first_of(delim, start);
				if (index == std::string::npos) {
					dest.push_back(str.substr(start, str.length() - start));
					return;
				}
				else {
					dest.push_back(str.substr(start, index - start));
				}
				start = str.find_first_not_of(delim, index);
			}
		}
	}
	void splitLastCharacter(std::string& str, std::string& dest)//get the substring before last /
	{
		std::string::size_type index = str.find_last_of("\\");
		if (index == std::string::npos) {
			index = str.find_last_of("/");
			if (index == std::string::npos) {
				dest = str;
			}
			else {
				dest = str.substr(0, index + 1);
			}
		}
		else {
			dest = str.substr(0, index + 1);
		}
	}


};

#endif // !READFILE_H

#pragma once