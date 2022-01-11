#pragma once

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include <iomanip>
#include <direct.h>
#include <io.h>

namespace WriteTxt {
	void writeTxt(std::vector<double>& result, std::vector<int>& time, int precision, std::string name, std::string first_line)
	{
		std::ofstream input_file;
		input_file.precision(precision);
		std::string obj_name = name + ".txt";
		//std::string material_name = name + ".mtl";
		input_file.open(obj_name.c_str(), std::ios::trunc);
		//input_file << first_line << "\t";
		for (int i = 0; i < result.size(); ++i) {
			input_file << result[i]<<"\t";
		}
		input_file << "\n";
		//input_file << "time ";
		for (int i = 0; i < time.size(); ++i) {
			input_file << time[i] << " ";
		}
		input_file << "\n";
		input_file.close();
		std::cout << "write "<<obj_name<< result.size() << std::endl;
	}


	void addToTxt(std::vector<double>& result, std::vector<int>& time, int precision, std::string name, std::string first_line)
	{
		std::ofstream input_file;
		input_file.precision(precision);
		std::string obj_name = name + ".txt";
		input_file.open(obj_name.c_str(), std::ios::app);
		//input_file << first_line << " ";
		for (int i = 0; i < result.size(); ++i) {
			input_file << result[i] << "\t";
		}
		input_file << "\n";
		//input_file << "time ";
		for (int i = 0; i < time.size(); ++i) {
			input_file << time[i] << "\t";
		}
		input_file << "\n";
		input_file.close();
	}
}


