#pragma once

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include <iomanip>
#include <direct.h>
#include <io.h>
#include"enum_setting.h"
#include<array>
//we need to record object number, path, anchor vertex of every object, stiffness, method(PD, XPBD, newton), time step, 
//for XPBD, we need to record substep number, iteration number in every substep.
//for PD, we need to record the convergence rate for both inner and outer iteration

namespace SaveParameter{

	inline void split(std::string& str, const std::string& delim, std::vector<std::string>& dest)
	{
		dest.clear();
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


	inline void readFile(std::string& scene_path, std::vector<std::string>& obj_path, std::vector<std::string>& collider_path,
		std::vector<std::vector<double>>& obj_stiffness, std::vector<std::vector<double>>& collide_stiffness,
		std::vector<std::vector<int>>& anchor_vertex, double& time_step, unsigned int& use_method, unsigned int& sub_step_num, unsigned int& iteration_num,
		double& local_conv_rate, double& outer_conv_rate, double& cloth_density, double& tet_density, double& velocity_damp, double* friction_coe,
		unsigned int& sub_step_per_detection, bool& floor_exist, unsigned int& floor_dimension, bool& floor_normal_direction, double& floor_value)
	{

		std::string line;
		std::ifstream in(scene_path.c_str(), std::ios::in);
		if (!in.good())
		{
			std::cout << "file not open" << std::endl;
			return;
		}
		std::getline(in, line);
		if (line != "object") {
			std::cout << "error read object " << std::endl;
			return;
		}
		std::getline(in, line);
		obj_path.resize(std::stoi(line));
		for (unsigned int i = 0; i < obj_path.size(); ++i) {
			std::getline(in, line);
			obj_path[i] = line;
		}
		obj_stiffness.resize(obj_path.size());
		collide_stiffness.resize(obj_path.size());
		std::getline(in, line);
		if (line != "collider") {
			std::cout << "error read collider " << std::endl;
			return;
		}
		std::getline(in, line);
		if (std::stoi(line) > 0) {
			collider_path.resize(std::stoi(line));
			for (unsigned int i = 0; i < collider_path.size(); ++i) {
				std::getline(in, line);
				collider_path[i] = line;
			}
		}
		std::getline(in, line);
		if (line != "stiffness") {
			std::cout << "error read stiffness " << std::endl;
			return;
		}
		std::vector<std::string> vec;
		for (unsigned int i = 0; i < obj_path.size(); ++i) {
			obj_stiffness[i].reserve(6);
			collide_stiffness[i].reserve(8);
			std::getline(in, line);
			split(line, " ", vec);
			for (unsigned int j = 0; j < vec.size(); ++j) {
				obj_stiffness[i].emplace_back(stod(vec[j]));
			}
			std::getline(in, line);
			split(line, " ", vec);
			for (unsigned int j = 0; j < vec.size(); ++j) {
				collide_stiffness[i].emplace_back(stod(vec[j]));
			}
		}
		std::getline(in, line);
		if (line != "anchor_vertex") {
			std::cout << "error read anchor vertex index " << std::endl;
			return;
		}
		anchor_vertex.resize(obj_path.size());
		for (unsigned int i = 0; i < obj_path.size(); ++i) {
			std::getline(in, line);
			std::getline(in, line);
			split(line, " ", vec);
			if (std::stoi(vec[0]) > 0) {
				if (std::stoi(vec[0]) != vec.size() - 1) {
					std::cout << "the anchor vertex number is not equal to the record number" << std::endl;
					return;
				}
				anchor_vertex[i].reserve(std::stoi(vec[0]));
				for (unsigned int j = 1; j < vec.size(); ++j) {
					anchor_vertex[i].emplace_back(stod(vec[j]));
				}
			}
		}
		std::getline(in, line);
		if (line != "time_step") {
			std::cout << "error read time step" << std::endl;
			return;
		}
		std::getline(in, line);
		time_step = std::stod(line);
		std::getline(in, line);
		if (line != "cloth_density") {
			std::cout << "error read cloth density" << std::endl;
			return;
		}
		std::getline(in, line);
		cloth_density = std::stod(line);
		std::getline(in, line);
		if (line != "tet_density") {
			std::cout << "error read tet density" << std::endl;
			return;
		}
		std::getline(in, line);
		tet_density = std::stod(line);
		std::getline(in, line);
		if (line != "velocity_damp") {
			std::cout << "error read velocity damp" << std::endl;
			return;
		}
		std::getline(in, line);
		velocity_damp = std::stod(line);
		std::getline(in, line);
		if (line != "friction_coe") {
			std::cout << "error read friction coe" << std::endl;
			return;
		}
		std::getline(in, line);
		split(line, " ", vec);
		for (unsigned int i = 0; i < vec.size(); ++i) {
			friction_coe[i] = std::stod(vec[i]);
		}		
		std::getline(in, line);
		if (line != "floor_dimension") {
			std::cout << "error read floor dimension" << std::endl;
			return;
		}
		std::getline(in, line);
		int dimension = std::stoi(line);
		if (dimension==4) {
			floor_exist = false;
		}
		else {
			floor_exist = true;
			if (dimension > 0) {
				floor_normal_direction = true;
			}
			else {
				floor_normal_direction = false;
			}
			floor_dimension = std::abs(dimension) - 1;
		}
		std::getline(in, line);
		if (line != "floor_value") {
			std::cout << "error read floor value" << std::endl;
			return;
		}
		std::getline(in, line);
		floor_value = std::stod(line);
		std::getline(in, line);
		if (line == "XPBD") {
			use_method = XPBD_;
		}
		else if (line == "PD") {
			use_method = PD_;
		}
		else if (line == "newton") {
			use_method = NEWTON_;
		}
		else  if(line=="second_order_large_XPBD") {
			use_method = XPBD_SECOND_ORDER_LARGE_;
		}
		else if (line == "XPBD_IPC") {
			use_method = XPBD_IPC_;
		}
		else {
			std::cout << "error reading the simulation method" << std::endl;
			return;
		}
		switch (use_method)
		{
		case XPBD_IPC_:
		{
			std::getline(in, line);
			if (line != "substep_num") {
				std::cout << "error read substep num" << std::endl;
				return;
			}
			std::getline(in, line);
			sub_step_num = std::stoi(line);
			std::getline(in, line);
			if (line != "iteration_num") {
				std::cout << "error read iteration num" << std::endl;
				return;
			}
			std::getline(in, line);
			iteration_num = std::stoi(line);
			std::getline(in, line);
			if (line != "sub_step_per_detection") {
				std::cout << "error read sub_step per detection" << std::endl;
			}
			std::getline(in, line);
			sub_step_per_detection = std::stoi(line);
		}
		break;
		case XPBD_:
		{
			std::getline(in, line);
			if (line != "substep_num") {
				std::cout << "error read substep num" << std::endl;
				return;
			}
			std::getline(in, line);
			sub_step_num = std::stoi(line);
			std::getline(in, line);
			if (line != "iteration_num") {
				std::cout << "error read iteration num" << std::endl;
				return;
			}
			std::getline(in, line);
			iteration_num = std::stoi(line);
			std::getline(in, line);
			if (line != "sub_step_per_detection") {
				std::cout << "error read sub_step per detection" << std::endl;
			}
			std::getline(in, line);
			sub_step_per_detection=std::stoi(line);
		}
			break;
		case PD_:
		{
			std::getline(in, line);
			if (line != "convergence_rate") {
				std::cout << "error read convergence rate" << std::endl;
				return;
			}
			std::getline(in, line);
			local_conv_rate = std::stod(line);
			std::getline(in, line);
			outer_conv_rate = std::stod(line);
		}
		break;
		}
	}


	

	inline void writeBasicPara(std::ofstream& input_file, std::vector<std::string>& path, std::vector<std::string>& collider_path,
		unsigned int use_method, std::vector<std::vector<int>*>& anchor_veretx, double time_step,
		std::vector<std::array<double, 6>>& cloth_stiffness, std::vector<std::array<double, 6>>& tet_stiffness,
		std::vector<std::array<double, 8>>& cloth_collision_stiffness, std::vector<std::array<double, 8>>& tet_collision_stiffness, double cloth_density, double tet_density,
		double velocity_damp, double* friction_coe, bool floor_exist, int floor_dimension, bool floor_normal_direction, double floor_value)
	{
		//input_file.precision(64);
		input_file << "object" << "\n";
		input_file << path.size() << "\n";
		for (unsigned int i = 0; i < path.size(); ++i) {
			input_file << path[i] << "\n";
		}
		input_file << "collider" << "\n";
		input_file << collider_path.size() << "\n";
		for (unsigned int i = 0; i < collider_path.size(); ++i) {
			input_file << collider_path[i] << "\n";
		}
		input_file << "stiffness" << "\n";
		for (unsigned int i = 0; i < cloth_stiffness.size(); ++i) {
			for (unsigned int j = 0; j < cloth_stiffness.data()[i].size(); ++j) {
				input_file << cloth_stiffness.data()[i][j] << " ";
			}
			input_file << "\n";
			for (unsigned int j = 0; j < cloth_collision_stiffness.data()[i].size(); ++j) {
				input_file << cloth_collision_stiffness.data()[i][j] << " ";
			}
			input_file << "\n";
		}
		for (unsigned int i = 0; i < tet_stiffness.size(); ++i) {
			for (unsigned int j = 0; j < tet_stiffness.data()[i].size(); ++j) {
				input_file << tet_stiffness.data()[i][j] << " ";
			}
			input_file << "\n";
		}
		for (unsigned int i = 0; i < tet_collision_stiffness.size(); ++i) {
			for (unsigned int j = 0; j < tet_collision_stiffness.data()[i].size(); ++j) {
				input_file << tet_collision_stiffness.data()[i][j] << " ";
			}
			input_file << "\n";
		}
		input_file << "anchor_vertex" << "\n";
		for (unsigned int i = 0; i < anchor_veretx.size(); ++i) {		
			input_file << i << "\n";
			input_file << anchor_veretx[i]->size() << " ";
			if (!anchor_veretx[i]->empty()) {
				for (unsigned int j = 0; j < anchor_veretx[i]->size(); ++j) {
					input_file << anchor_veretx[i]->data()[j] << " ";
				}		
			}		
			input_file << "\n";
		}		
		input_file << "time_step" << "\n";
		input_file << time_step << "\n";
		input_file<<"cloth_density" << "\n";
		input_file << cloth_density << "\n";
		input_file << "tet_density" << "\n";
		input_file << tet_density << "\n";
		input_file << "velocity_damp" << "\n";
		input_file << velocity_damp << "\n";
		input_file << "friction_coe" << "\n";
		input_file << friction_coe[0]<<" "<<friction_coe[1]<<" "<< friction_coe[2] << "\n";
		input_file << "floor_dimension" << "\n";
		if (!floor_exist) {
			input_file << 4 << "\n";
		}
		else {
			if (floor_normal_direction) {
				input_file << floor_dimension + 1 << std::endl;
			}
			else {
				input_file << -1*(floor_dimension + 1) << std::endl;
			}
		}
		input_file << "floor_value" << "\n";
		input_file << floor_value<< "\n";
	}


	inline void writeParameter(std::vector<std::string>& path, std::vector<std::string>& collider_path, std::vector<std::array<double,6>>& cloth_stiffness, 
		std::vector<std::array<double, 6>>& tet_stiffness, std::vector<std::array<double, 8>>& cloth_collision_stiffness, std::vector<std::array<double, 8>>& tet_collision_stiffness,
		unsigned int use_method, std::vector<std::vector<int>*>&anchor_veretx, double time_step, double outer_convergence_rate,
		double local_convergence_rate,
		unsigned int sub_step_num, unsigned int iteration_num, double cloth_density, double tet_density, double velocity_damp,
		double* friction_coe, unsigned int sub_step_per_detection, bool floor_exist, int floor_dimension, bool floor_normal_direction, double floor_value)
	{
		std::string prefix = "./record_scene_data/";
		if (_access(prefix.c_str(), 0) == -1)
			_mkdir(prefix.c_str());

		std::ofstream input_file;
		std::string obj_name = prefix + "scene.scene";
		input_file.open(obj_name.c_str(), std::ios::trunc);
		writeBasicPara(input_file, path, collider_path, use_method, anchor_veretx, time_step,cloth_stiffness,tet_stiffness, cloth_collision_stiffness, 
			tet_collision_stiffness,cloth_density,tet_density,velocity_damp,friction_coe,floor_exist,floor_dimension,floor_normal_direction,floor_value);
		switch (use_method)
		{
		case XPBD_:
			input_file << "XPBD"<< "\n";
			input_file << "substep_num"<< "\n";
			input_file << sub_step_num << "\n";
			input_file << "iteration_num" << "\n";
			input_file << iteration_num << "\n";
			input_file << "sub_step_per_detection" << "\n";
			input_file << sub_step_per_detection << "\n";
			break;
		case PD_:
			input_file << "PD" << "\n";
			input_file << "convergence_rate" << "\n";
			input_file << local_convergence_rate << "\n";
			input_file << outer_convergence_rate << "\n";
			break;
		case NEWTON_:
			input_file << "newton" << "\n";
			break;
		case XPBD_SECOND_ORDER_LARGE_:
			input_file << "second_order_large_XPBD" << "\n";
			break;
		case XPBD_IPC_:
			input_file << "XPBD_IPC" << "\n";
			input_file << "substep_num" << "\n";
			input_file << sub_step_num << "\n";
			input_file << "iteration_num" << "\n";
			input_file << iteration_num << "\n";
			input_file << "sub_step_per_detection" << "\n";
			input_file << sub_step_per_detection << "\n";
			break;
		}		
	}





};


