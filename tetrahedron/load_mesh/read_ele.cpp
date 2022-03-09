#include"read_ele.h"

void ReadEle::load(const char* name, OriMesh& mesh)
{
	std::vector<std::vector<std::string>> read_string;
	std::vector<std::string> tem_string;
	std::string node_path;
	std::string tetro_path(name);
	setNodePath(node_path, tetro_path);
	loadFile(name, tem_string);
	splitting(tem_string, read_string);
	tem_string.clear();
	saveTetrahedronIndex(read_string, mesh);
	read_string.clear();
	loadFile(node_path.c_str(), tem_string);
	splitting(tem_string, read_string);
	tem_string.clear();
	saveTetrahedronVertex(read_string, mesh);

	read_string.clear();
}

void ReadEle::loadFile(const char* name, std::vector<std::string>& temString)
{
	std::string line;
	std::ifstream fin;
	fin.open(name);
	if (!fin.is_open()) {
		std::cout << "error opening " << *name << std::endl;
	}
	else {
		temString.reserve(2000);
		while (std::getline(fin, line)) {
			temString.push_back(line);
		}
		fin.close();
	}
}

void ReadEle::splitting(std::vector<std::string>& ori, std::vector<std::vector<std::string>>& result)
{
	result.reserve(ori.size());

	for (int i = 0; i < ori.size(); i++)
	{
		std::vector<std::string>tem;
		split(ori[i], " ", "\t", tem);
		result.push_back(tem);
	}
}


void ReadEle::split(std::string& str, const std::string& delim1, const std::string& delim2, std::vector<std::string>& dest)
{
	std::string::size_type start = str.find_first_not_of(delim1, 0), index;
	std::string::size_type start2 = str.find_first_not_of(delim2, 0), index2;
	if (start < start2) {
		start = start2;
	}
	if (str.find_last_of(delim1) == str.length() - 1 || str.find_last_of(delim2) == str.length() - 1) {
		while (start != std::string::npos)
		{
			index = str.find_first_of(delim1, start);
			index2 = str.find_first_of(delim2, start);
			if (index > index2) {
				index = index2;
			}
			if (index == std::string::npos) {
				return;
			}
			else {
				dest.push_back(str.substr(start, index - start));
			}
			start = str.find_first_not_of(delim1, index);
			start2 = str.find_first_not_of(delim2, index);
			if (start < start2 && str[start] == delim2.c_str()[0]) {
				start = start2;
			}
			else if (start > start2 && str[start2] != delim1.c_str()[0]) {
				start = start2;
			}
		}
	}
	else {
		while (start != std::string::npos)
		{
			index = str.find_first_of(delim1, start);
			index2 = str.find_first_of(delim2, start);
			if (index > index2) {
				index = index2;
			}
			if (index == std::string::npos) {
				dest.push_back(str.substr(start, str.length() - start));
				return;
			}
			else {
				dest.push_back(str.substr(start, index - start));
			}
			start = str.find_first_not_of(delim1, index);
			start2 = str.find_first_not_of(delim2, index);
			if (start < start2 && str[start] == delim2.c_str()[0]) {
				start = start2;
			}
			else if (start > start2 && str[start2] != delim1.c_str()[0]) {
				start = start2;
			}
		}
	}
}

void ReadEle::saveTetrahedronIndex(std::vector<std::vector<std::string>>& result, OriMesh& mesh)
{
	bool start_from_zero = false;
	int tetrahedron_num = stoi(result[0][0]);
	mesh.indices.reserve(4 * tetrahedron_num);
	for (int i = 1; i < result.size() - 1; i++) {
		for (int j = 1; j < 5; j++) {
			if (result[i][j] == "0") {
				start_from_zero = true;
			}
			mesh.indices.push_back(stoi(result[i][j]));
		}
	}
	if (result[result.size() - 1][0] != "#") {
		if (result[result.size() - 1].size() > 1) {
			for (int j = 1; j < 5; j++) {
				if (result[result.size() - 1][j] == "0") {
					start_from_zero = true;
				}
				mesh.indices.push_back(stoi(result[result.size() - 1][j]));
			}
		}
	}
	if (!start_from_zero) {
		for (int i = 0; i < mesh.indices.size(); ++i) {
			mesh.indices[i]--;
		}
	}
}

void ReadEle::saveTetrahedronVertex(std::vector<std::vector<std::string>>& result, OriMesh& mesh)
{
	int tetrahedron_vertex_num = stoi(result[0][0]);
	mesh.vertices.reserve(3 * tetrahedron_vertex_num);
	for (int i = 1; i < result.size() - 1; ++i) {
		mesh.vertices.push_back({ stod(result[i][1]),stod(result[i][2]),stod(result[i][3]) });
	}
	if (result[result.size() - 1][0] != "#") {
		mesh.vertices.push_back({ stod(result[result.size() - 1][1]),stod(result[result.size() - 1][2]),stod(result[result.size() - 1][3]) });
	}
	if (mesh.vertices.size() != tetrahedron_vertex_num) {
		std::cout << "error: the vertex num is not compatible." << std::endl;
	}
}

void ReadEle::setNodePath(std::string& node_path, std::string& ele_path)
{
	std::string::size_type index;
	index = ele_path.find_last_of(".");
	node_path = ele_path.substr(0, index) + ".node";
}
