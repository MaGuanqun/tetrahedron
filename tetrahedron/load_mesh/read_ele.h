#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include <iomanip>
#include<time.h>
#include"../basic/global.h"
#include<array>
#include"../global_struct.h"

class ReadEle
{
public:
	void load(const char* name, OriMesh& mesh);
 

private:
	void loadFile(const char* name, std::vector<std::string>& temString);
	void splitting(std::vector<std::string>& ori, std::vector<std::vector<std::string>>& result);
	void split(std::string& str, const std::string& delim1, const std::string& delim2, std::vector<std::string>& dest);
	void saveTetrahedronIndex(std::vector<std::vector<std::string>>& result, OriMesh& mesh);
	void setNodePath(std::string& node_path, std::string& ele_path);
	void saveTetrahedronVertex(std::vector<std::vector<std::string>>& result, OriMesh& mesh);
};

