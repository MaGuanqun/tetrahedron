#include"save_scene.h"
void SaveScene::save_scene_XPBD(size_t time_stamp, unsigned int simulate_scene_indicator, std::vector<MeshStruct*>& obj_mesh_struct,
	std::vector<std::vector<std::array<double, 3>>>* velocity, std::vector<MeshStruct*>& collider_mesh_struct)
{
	if (record_time_stamp == time_stamp) {
		return;
	}
	record_time_stamp = time_stamp;

	std::string prefix = "./record_simulation_data/";
	if (_access(prefix.c_str(), 0) == -1)
		_mkdir(prefix.c_str());
	std::string file_name = prefix + "obj_";
	std::string final_file_name;
	final_file_name = file_name+std::to_string(time_stamp) + ".dat";
	std::ofstream input_file(final_file_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
	//input_file.precision(64);
	input_file.write((char*)&time_stamp, sizeof(size_t));
	input_file.write((char*)&simulate_scene_indicator, sizeof(unsigned int));
	unsigned int obj_num = obj_mesh_struct.size();
	input_file.write((char*)&obj_num, sizeof(unsigned int));
	obj_num = collider_mesh_struct.size();
	input_file.write((char*)&obj_num, sizeof(unsigned int));
	for (unsigned int i = 0; i < obj_mesh_struct.size(); ++i) {		
		unsigned int vertex_num = obj_mesh_struct[i]->vertex_position.size();	
		write_binary_add(input_file, final_file_name.c_str(), obj_mesh_struct[i]->vertex_position[0].data(), 3 * vertex_num);
		write_binary_add(input_file, final_file_name.c_str(), velocity->data()[i][0].data(), 3 * vertex_num);
	}

	for (unsigned int i = 0; i < collider_mesh_struct.size(); ++i) {
		unsigned int vertex_num = collider_mesh_struct[i]->vertex_position.size();
		write_binary_add(input_file, final_file_name.c_str(), collider_mesh_struct[i]->vertex_position[0].data(), 3 * vertex_num);
	}
	input_file.close();

}



bool SaveScene::read_scene_XPBD(const char* file_name, size_t* time_stamp, unsigned int* simulate_scene_indicator, std::vector<MeshStruct*>& obj_mesh_struct,
	std::vector<std::vector<std::array<double, 3>>>* velocity, std::vector<MeshStruct*>& collider_mesh_struct)
{
	std::ifstream in(file_name, std::ios::in | std::ios::binary);
	if (!in.good())
	{
		std::cout << "file not open" << std::endl;
		return false;
	}

	in.read((char*)time_stamp, sizeof(size_t));
	in.read((char*)simulate_scene_indicator, sizeof(unsigned int));

	unsigned int obj_num = 0;
	unsigned int collider_num=0;

	in.read((char*)&obj_num, sizeof(unsigned int));
	in.read((char*)&collider_num, sizeof(unsigned int));




	if (obj_num != obj_mesh_struct.size()) {
		std::cout << "the object count does not match" << std::endl;
		return false;
	}
	if (collider_num != collider_mesh_struct.size()) {
		std::cout << "the collider count does not match" << std::endl;
		return false;
	}

	for (unsigned int i = 0; i < obj_num; ++i) {
		if (!read_binary(in, obj_mesh_struct[i]->vertex_position[0].data(), 3*obj_mesh_struct[i]->vertex_position.size())) {
			return false;
		}
		if (!read_binary(in, velocity->data()[i][0].data(), 3*obj_mesh_struct[i]->vertex_position.size())) {
			return false;
		}
	}

	for (unsigned int i = 0; i < collider_num; ++i) {
		std::cout << "coll " << std::endl;
		if (!read_binary(in, collider_mesh_struct[i]->vertex_position[0].data(), 3*collider_mesh_struct[i]->vertex_position.size())) {
			return false;
		}
	}




	return true;

}