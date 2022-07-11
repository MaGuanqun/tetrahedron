#ifndef CREATE_MESH_H
#define CREATE_MESH_H
#include"../load_mesh/readObj.h"
#include"../external/glm/glm.hpp"
class CreateMesh
{
public:


	CreateMesh(OriMesh& mesh)
	{


	}

	void setMaterial1(OriMesh& mesh)
	{



		//only one cloth
		//int col_num = 15;
		//int row_num = 15;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{ 0.35 * (double)j / (double)(row_num - 1) -0.175, 0.15 * (double)i / (double)(col_num - 1)+0.3,
		//			0.175 - 0.35 * (double)i / (double)(col_num - 1) });// 
		//	}
		//}
		//int col_num = 6;
		//int row_num = 6;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{ 0.14 * (double)j / (double)(row_num - 1) - 0.07, -0.38,
		//			0.07 - 0.14 * (double)i / (double)(col_num - 1) });// 
		//	}
		//}

		//int col_num = 2;
		//int row_num = 2;

		//mesh.vertices.push_back({ -0.025,0.2,0.0 });
		//mesh.vertices.push_back({ 0.0,0.2,-0.025 });
		//mesh.vertices.push_back({ 0.0,0.2,0.0125 });
		//mesh.vertices.push_back({ 0.025,0.2,0.0 });

		//mesh.vertices.push_back({ -0.025,0.2,0.0 });
		//mesh.vertices.push_back({ -0.04,0.2,-0.04 });
		//mesh.vertices.push_back({ 0.04,0.2,0.04 });
		//mesh.vertices.push_back({ 0.025,0.2,0.0 });

		//mesh.vertices.push_back({ -0.125,0.2,-0.125 });
		////mesh.vertices.push_back({ -0.125,0.2,0.125 });
		//mesh.vertices.push_back({ 0.125,0.2,-0.125 });
		//mesh.vertices.push_back({ 0.125,0.2,0.125 });

		//mesh.vertices.push_back({ 0.0,0.2,0.0 });
		//mesh.vertices.push_back({ 0.02,0.2,0.02 });
		//mesh.vertices.push_back({ -0.02,0.2,0.02 });
		//mesh.vertices.push_back({ -0.02,0.2,-0.02 });
		//mesh.vertices.push_back({ 0.02,0.2,-0.02 });


		//mesh.vertices.push_back({ 0.0,0.125,0.0 });
		//mesh.vertices.push_back({ 0.125,0.0,0.0 });
		//mesh.vertices.push_back({ -0.125,0.0,0.0 });
		//mesh.vertices.push_back({ 0.0,-0.2,0.0 });
		//mesh.vertices.push_back({ 0.0,0.0,0.0 });
		//mesh.vertices.push_back({ 0.02,-0.02,0.0 });
		//mesh.vertices.push_back({ -0.02,-0.02,0.0 });
		//mesh.vertices.push_back({ -0.02,0.2,0.02 });
		//mesh.vertices.push_back({ -0.02,0.2,-0.02 });

		//indicesForTwoTriangle(mesh);

		//int col_num = 70;
		//int row_num = 70;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{1.4 * (double)i / (double)(col_num - 1) - 0.7,  //+0.9
		//			1.4 * (double)j / (double)(row_num - 1) - 0.7, 0.2 * (double)j / (double)(row_num - 1) - 0.1});//
		//	}
		//}

		//int col_num = 3;
		//int row_num = 3;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		//mesh.vertices.push_back(std::array<double, 3>{ 0.7 * (double)j / (double)(row_num - 1) + 0.575, 0.35 * (double)i / (double)(col_num - 1),//+0.9
		//		//mesh.vertices.push_back(std::array<double, 3>{0.7* (double)i / (double)(col_num - 1) - 0.35,  //+0.9
		//		//	0.4 * (double)j / (double)(row_num - 1) - 0.285, 0.7 * (double)j / (double)(row_num - 1) - 0.35});//
		//		mesh.vertices.push_back(std::array<double, 3>{0.07 * (double)i / (double)(col_num - 1) - 0.035,  //+0.9
		//			0.07 * (double)j / (double)(row_num - 1) - 0.6, 0.05 * (double)j / (double)(row_num - 1) - 0.01});//
		//	}
		//}

		//int col_num = 50;
		//int row_num = 50;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		//mesh.vertices.push_back(std::array<double, 3>{ 0.7 * (double)j / (double)(row_num - 1) + 0.575, 0.35 * (double)i / (double)(col_num - 1),//+0.9
		//		//mesh.vertices.push_back(std::array<double, 3>{0.7* (double)i / (double)(col_num - 1) - 0.35,  //+0.9
		//		//	0.4 * (double)j / (double)(row_num - 1) - 0.285, 0.7 * (double)j / (double)(row_num - 1) - 0.35});//
		//		mesh.vertices.push_back(std::array<double, 3>{1.0 * (double)i / (double)(col_num - 1),  // - 0.35+0.9
		//			0.0,1.0 * (double)j / (double)(row_num - 1)}); //-0.585//0.05 * (double)j / (double)(row_num - 1) - 0.1
		//	}
		//}

		//int col_num = 75;
		//int row_num = 75;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		//mesh.vertices.push_back(std::array<double, 3>{ 0.7 * (double)j / (double)(row_num - 1) + 0.575, 0.35 * (double)i / (double)(col_num - 1),//+0.9
		//		//mesh.vertices.push_back(std::array<double, 3>{0.7* (double)i / (double)(col_num - 1) - 0.35,  //+0.9
		//		//	0.4 * (double)j / (double)(row_num - 1) - 0.285, 0.7 * (double)j / (double)(row_num - 1) - 0.35});//
		//		mesh.vertices.push_back(std::array<double, 3>{2.0 * (double)i / (double)(col_num - 1)-1.0,  // - 0.35+0.9
		//			0.0 * (double)j / (double)(row_num - 1)+0.6,2.0 * (double)j / (double)(row_num - 1)-1.0}); //-0.585//0.05 * (double)j / (double)(row_num - 1) - 0.1
		//	}
		//}

		//int col_num = 100;
		//int row_num = 100;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{2.0 * (double)i / (double)(col_num - 1) - 1.0,  //+0.9
		//			1.8 * (double)j / (double)(row_num - 1) - 0.6, 0.4 * (double)j / (double)(row_num - 1) - 0.2});//
		//	}
		//}

		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{1.0 * (double)i / (double)(col_num - 1) - 0.45,  //+0.9
		//			0.954 * (double)j / (double)(row_num - 1) - 0.6, 0.3 * (double)j / (double)(row_num - 1) - 0.6});// 
		//	}
		//}

		//this is for adding force
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{1.0 * (double)i / (double)(col_num - 1) - 0.45,  //+0.9
		//			1.0 * (double)j / (double)(row_num - 1) - 0.6, 0.2 * (double)j / (double)(row_num - 1) - 0.6});// 
		//	}
		//}

		//band falling
		int col_num = 71;
		int row_num = 3;
		for (int j = 0; j < row_num; j++) {
			for (int i = 0; i < col_num; i++) {		
				mesh.vertices.push_back(std::array<double, 3>{0.05* (double)i / (double)(col_num - 1)-0.025,  //+0.9
					1.4 * (double)i / (double)(col_num - 1) - 1.0, 0.04* (double)j / (double)(row_num - 1) - 0.07});//
			}
		}



		//this is for sphere rotating
		//int col_num = 120;
		//int row_num = 120;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {		
		//		mesh.vertices.push_back(std::array<double, 3>{1.77 * (double)i / (double)(col_num - 1) - 1.0,  //+0.9
		//			0.25, 1.77* (double)j / (double)(row_num - 1) - 0.5});//
		//	}
		//}

		//int col_num = 10;
		//int row_num = 10;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		//mesh.vertices.push_back(std::array<double, 3>{ 0.7 * (double)j / (double)(row_num - 1) + 0.575, 0.35 * (double)i / (double)(col_num - 1),//+0.9
		//		//mesh.vertices.push_back(std::array<double, 3>{0.7* (double)i / (double)(col_num - 1) - 0.35,  //+0.9
		//		//	0.4 * (double)j / (double)(row_num - 1) - 0.285, 0.7 * (double)j / (double)(row_num - 1) - 0.35});//
		//		mesh.vertices.push_back(std::array<double, 3>{0.7 * (double)i / (double)(col_num - 1) - 0.35,  //+0.9
		//			0.7 * (double)j / (double)(row_num - 1) - 0.585, 0.1 * (double)j / (double)(row_num - 1) - 0.1});//
		//	}
		//}		
		//int col_num = 70;
		//int row_num = 70;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		//mesh.vertices.push_back(std::array<double, 3>{ 0.7 * (double)j / (double)(row_num - 1) + 0.575, 0.35 * (double)i / (double)(col_num - 1),//+0.9
		//		mesh.vertices.push_back(std::array<double, 3>{1.64* (double)i / (double)(col_num - 1) - 0.82,  //+0.9
		//			1.64 * (double)j / (double)(row_num - 1) - 0.42, 0.4 * (double)j / (double)(row_num - 1) - 0.2});//
		//	}
		//}
		//int col_num = 2;
		//int row_num = 2;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.14 * (double)i / (double)(col_num - 1) - 0.07,  //+0.9
		//			0.2, 0.14 * (double)j / (double)(row_num - 1)-0.07});//0.1 * (double)j / (double)(row_num - 1) - 0.05
		//	}
		//}
		//setBandRotate();
		//double cloth_center[3] = { 0.75,1.075,0.8 };
		//rotateCloth(cloth_center, mesh);


		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{ 0.0, 0.5 * (double)j / (double)(row_num - 1) + 1.25,
		//			0.5 * (double)i / (double)(col_num - 1) + 1.5 });// 
		//	}
		//}

		//for (int j = 0; j < col_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.0 * (double)i / (double)(col_num - 1) + 0.3,
		//			0.8 * (double)i / (double)(col_num - 1) + 0.2, 0.7 * (double)j / (double)(col_num - 1) - 0.35});
		//	}
		//}
		//col_num=44;
		//for (int j = 0; j < col_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{ 0.1 * (double)i / (double)(col_num - 1) -0.15,
		//			0.7 * (double)i / (double)(col_num - 1)+0.55, 0.7 * (double)j / (double)(col_num - 1) -0.35});	
		//	}
		//}

		//move capsule
		//int col_num = 41;
		//int row_num = 41;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.35 * (double)j / (double)(row_num - 1) - 0.175,
		//			0.7 * (double)i / (double)(col_num - 1)-0.7, -0.3 * (double)i / (double)(col_num - 1) + 0.15});
		//	}
		//}
		//////band
		//int col_num = 115;
		//int row_num = 9;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.25 * (double)j / (double)(row_num - 1) -0.125,
		//			-0.2,2.9 * (double)i / (double)(col_num - 1) -1.45});
		//	}
		//}
		// 
		// 
		//int col_num = 91;
		//int row_num = 3;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.05 * (double)j / (double)(row_num - 1), 1.5 * (double)i / (double)(col_num - 1) + 0.05,
		//			0.05 * (double)i / (double)(col_num - 1)});//
		//	}
		//}


		//int col_num = 10;
		//int row_num = 10;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.2 * (double)j / (double)(row_num - 1)+1.5, 0.2 * (double)i / (double)(col_num - 1) + 0.3,
		//			1.5 + 0.005 * (double)i / (double)(col_num - 1)});//
		//	}
		//}

		//int col_num = 2;
		//int row_num = 2;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.04 * (double)j / (double)(row_num - 1) + 1.5, 0.04 * (double)i / (double)(col_num - 1) + 0.3,
		//			1.5 + 0.005 * (double)i / (double)(col_num - 1)});//
		//	}
		//}

		//int col_num = 91;
		//int row_num = 3;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.05 * (double)j / (double)(row_num - 1) + 1.5, 1.9 * (double)i / (double)(col_num - 1) + 0.3,
		//			1.5+ 0.02 * (double)i / (double)(col_num - 1) });
		//	}
		//}
		//sphere
		//int col_num = 60;
		//int row_num = 60;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{1.5 * (double)j / (double)(row_num - 1) - 0.75,
		//			0.3, 1.5 * (double)i / (double)(col_num - 1) - 0.75});
		//	}
		//}
		//cloth for vertically moving cloth
		//int col_num = 28;
		//int row_num = 14;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.0,
		//			0.68 * (double)i / (double)(col_num - 1) - 0.4, 0.28 * (double)j / (double)(row_num - 1) - 0.14});
		//	}
		//}

		//indices(mesh, col_num, row_num);
		meshIndices(mesh, col_num, row_num);

		//mesh.vertices.push_back({ 0.0125,0.2,0.0125 });		
		//mesh.vertices.push_back({ 0.0125,0.2,-0.0125 });
		//mesh.vertices.push_back({ -0.0125,0.2,0.0125 });
		//indicesForTwoTriangle(mesh);
		//setVirtualIndices(mesh,2, col_num, row_num);

		initialMaterial(mesh, 1);


	}



	void setMaterial3(OriMesh& mesh)
	{

		mesh.vertices.push_back({ 0.15,-0.15,0.03 });
		initialMaterial(mesh, 2);

	}


	void setMaterial2(OriMesh& mesh)
	{
		//int col_num = 40;
		//int row_num = 40;
		//for (int j = 0; j < col_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.15 * (double)i / (double)(col_num - 1),
		//			0.7 * (double)i / (double)(col_num - 1)+0.3, 0.8 * (double)j / (double)(col_num - 1)-0.3});
		//	}
		//}
		// 
		//mesh.vertices.push_back({ 0.0,-0.2,0.3 });
		//mesh.vertices.push_back({ 0.0,0.2,0.0 });
		//mesh.vertices.push_back({ -0.02,0.2,-0.02 });
		//mesh.vertices.push_back({ 0.02,0.2,-0.02 });
	// 

		// 
		// 
		//mesh.vertices.push_back({ 0.125,0.2,0.125 });		
		//mesh.vertices.push_back({ 0.125,0.2,-0.125 });
		//mesh.vertices.push_back({ -0.125,0.2,0.125 });
		//indicesForTwoTriangle(mesh);
		//
		//cloth falls on other cloth
		int col_num = 6;
		int row_num = 11;
		for (int j = 0; j < row_num; j++) {
			for (int i = 0; i < col_num; i++) {
				//mesh.vertices.push_back(std::array<double, 3>{ 0.7 * (double)j / (double)(col_num - 1) -0.35, 0.2,// 0.5 * (double)i / (double)(col_num - 1) + 0.6,
				//	0.7 * (double)i / (double)(col_num - 1) -0.35});// +  //, 1.35 - 0.2 * (double)i / (double)(col_num - 1)
				mesh.vertices.push_back(std::array<double, 3>{0.4 * (double)i / (double)(col_num - 1) - 0.2,  //+0.9
					0.7 * (double)j / (double)(row_num - 1) - 0.35,0.0});// 0. * (double)j / (double)(col_num - 1) - 0.8
			}
		}



		//for (int j = 0; j < col_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{ 0.7, 0.7 * (double)i / (double)(col_num - 1) + 0.1,// 0.5 * (double)i / (double)(col_num - 1) + 0.6,
		//			0.7 * (double)j / (double)(col_num - 1) + 0.4});// +  //, 1.35 - 0.2 * (double)i / (double)(col_num - 1)
		//	}
		//}
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{ 0.35 * (double)j / (double)(row_num - 1) + 0.575, 0.35 * (double)i / (double)(col_num - 1) + 0.9,
		//			0.8 });// 
		//	}
		//}
		//meshIndices(mesh, col_num, col_num);

		//col_num *= 0.5;
		//for (int j = 0; j < col_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		//mesh.vertices.push_back(std::array<double, 3>{0.05 * (double)i / (double)(col_num - 1) + 0.3,
		//		mesh.vertices.push_back(std::array<double, 3>{0.07,
		//			1.0 * (double)i / (double)(col_num - 1)+0.3, 1.0 * (double)j / (double)(col_num - 1) - 0.35});
		//	}
		//}
		//cloth for vertically moving cloth
		//int col_num = 21;
		//int row_num = 7;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{-0.025,
		//			0.51 * (double)i / (double)(col_num - 1) - 0.23, 0.14 * (double)j / (double)(row_num - 1) - 0.07});
		//	}
		//}

		indices(mesh, col_num, row_num);
		//setVirtualIndices(mesh, 2, col_num, row_num);
		//indicesForTwoTriangle(mesh);
		initialMaterial(mesh, 2);

	}

	void setFloor(OriMesh& mesh)
	{
		//int col_num = 30;
		//for (int j = 0; j <col_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{2.8*(double)i / (double)(col_num - 1) - 1.4,
		//			-0.7, 2.8*(double)j / (double)(col_num - 1) - 1.4});;
		//	}
		//}

		////for two prisms.
		int col_num = 40;
		int row_num = 40;
		for (int j = 0; j < row_num; j++) {
			for (int i = 0; i < col_num; i++) {
				mesh.vertices.push_back(std::array<double, 3>{1.6 * (double)i / (double)(col_num - 1) - 0.8,  //+0.9
					-0.35, 1.6 * (double)j / (double)(row_num - 1) - 0.8});//-0.525, - 0.425
			}
		}
		////for two dragons.
		//int col_num = 50;
		//int row_num = 100;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{1.0 * (double)i / (double)(col_num - 1) - 0.5,  //+0.9
		//			-0.525, 2.0 * (double)j / (double)(row_num - 1) - 1.0});//-0.525,
		//	}
		//}


		//int col_num = 100;
		//int row_num =50;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{2.0 * (double)i / (double)(col_num - 1) - 1.0,  //+0.9
		//			-0.525, 1.0 * (double)j / (double)(row_num - 1) - 0.5});//-0.525,
		//	}
		//}


		//int col_num = 100;
		//int row_num = 100;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{
		//			2.0 * (double)i / (double)(col_num - 1) - 1.0,
		//			-0.525,
		//			2.0 * (double)j / (double)(row_num - 1) - 1.0
		//		});//-0.525,
		//	}
		//}


		//int col_num = 100;
		//int row_num = 100;
		//for (int j = 0; j < row_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.525,				
		//			2.0 * (double)i / (double)(col_num - 1) - 0.75,
		//			2.0 * (double)j / (double)(row_num - 1) - 1.3
		//			});//-0.525,
		//	}
		//}


		//int col_num = 2;
		//for (int j = 0; j <col_num; j++) {
		//	for (int i = 0; i < col_num; i++) {
		//		mesh.vertices.push_back(std::array<double, 3>{0.5*(double)i / (double)(col_num - 1)-0.3,
		//			-0.5, 0.5*(double)j / (double)(col_num - 1)-0.1});
		//	}
		//}
		//mesh.vertices.push_back({-0.3,-0.5,-0.1});
		//mesh.vertices.push_back({0.2,-0.5,0.4});
		//mesh.vertices.push_back({0.2,-0.5,-0.1});
		//mesh.face_groups.resize(1);
		//for (int i = 0; i < 3; ++i) {
		//	mesh.indices.push_back(i);
		//}

		indices(mesh, col_num, row_num);
		floorMaterial(mesh);
		//for (int i = 0; i < mesh.vertices.size(); ++i) {
		//	mesh.vertices[i][1] -= 1.0;
		//}
	}

	void setCylinder(OriMesh& mesh, double r, double length, int col_num, int row_num, int type)
	{
		for (int j = 0; j < row_num; j++) {
			for (int i = 0; i < col_num; i++) {
				mesh.vertices.push_back(std::array<double, 3>{
					r* (1.0 + 1.5 * length * (1.0 - (double)j / (double)(row_num - 1)))* sin(2.0 * M_PI * (double)i / (double)col_num),
						0.5 - length * (double)(row_num - 1 - j) / (double)(row_num - 1),// - 0.5* length,				
						r* (1.0 + 1.5 * length * (1.0 - (double)j / (double)(row_num - 1)))* cos(2.0 * M_PI * (double)i / (double)col_num)});
			}
		}
		cylinderIndices(mesh, 0, col_num * row_num, col_num, row_num);
		initialMaterial(mesh, type);
	}



	void setCapsule(OriMesh& mesh)
	{
		int col_num = 34;
		int row_num = 34;
		double r = 0.12;
		int globe_row = 9;
		int globe_num_2 = 17;

		double half_capsule_length = 0.4;
		double capsule_length = 2.0 * half_capsule_length;

		mesh.vertices.push_back(std::array<double, 3>{half_capsule_length + r, 0.0, 0.0});
		for (int j = 1; j < globe_row; j++) {
			for (int i = 0; i < col_num; i++) {
				mesh.vertices.push_back(std::array<double, 3>{
					r* cos(M_PI* (double)j / (double)globe_num_2) + half_capsule_length,
						r* sin(M_PI* (double)j / (double)globe_num_2)* sin(2.0 * M_PI * (double)i / (double)col_num),
						r* sin(M_PI* (double)j / (double)globe_num_2)* cos(2.0 * M_PI * (double)i / (double)col_num)});
			}
		}

		int first_half_sphere = mesh.vertices.size() - col_num;

		for (int j = row_num - 2; j > 0; j--) {
			for (int i = 0; i < col_num; i++) {
				mesh.vertices.push_back(std::array<double, 3>{capsule_length* (double)j / (double)(row_num - 1) - half_capsule_length,
					r* sin(2.0 * M_PI * (double)i / (double)col_num),
					r* cos(2.0 * M_PI * (double)i / (double)col_num)});
			}
		}

		int before_last_half_sphere = mesh.vertices.size();
		int total_num = row_num * col_num;

		for (int j = globe_row - 1; j > 0; j--) {
			for (int i = 0; i < col_num; i++) {
				mesh.vertices.push_back(std::array<double, 3>{
					-r * cos(M_PI * (double)j / (double)globe_num_2) - half_capsule_length,
						r* sin(M_PI* (double)j / (double)globe_num_2)* sin(2.0 * M_PI * (double)i / (double)col_num),
						r* sin(M_PI* (double)j / (double)globe_num_2)* cos(2.0 * M_PI * (double)i / (double)col_num)});
			}
		}
		mesh.vertices.push_back(std::array<double, 3>{-half_capsule_length - r, 0.0, 0.0});



		globeIndices(globe_row, col_num, mesh, before_last_half_sphere);
		cylinderIndices(mesh, first_half_sphere, total_num, col_num, row_num);

	}

	void setSphere(OriMesh& mesh)
	{
		int col_num = 100;
		double r = 0.2;
		int globe_num_2 = 50;
		//double radius = 0.3;		
		for (int j = 1; j < globe_num_2; j++) {
			for (int i = 0; i < col_num; i++) {
				mesh.vertices.push_back(std::array<double, 3>{
					r* cos(M_PI* (double)j / (double)globe_num_2),
						r* sin(M_PI* (double)j / (double)globe_num_2)* sin(2.0 * M_PI * (double)i / (double)col_num),
						r* sin(M_PI* (double)j / (double)globe_num_2)* cos(2.0 * M_PI * (double)i / (double)col_num)});
			}
		}
		mesh.vertices.push_back(std::array<double, 3>{r, 0.0, 0.0});
		mesh.vertices.push_back(std::array<double, 3>{-r, 0.0, 0.0});
		sphereIndices(globe_num_2, col_num, mesh);

	}


private:



	glm::mat3x3 band_rotate;
	glm::mat3x3 band_rotate_reverse;
	glm::mat3x3 capsule_rotate;
	glm::mat3x3 capsule_rotate_reverse;

	void setBandRotate()
	{
		glm::vec3 rotate_axe = glm::vec3(0.0, 1.0, 0.0);
		glm::mat3x3 ux = glm::mat3x3(0, -rotate_axe.z, rotate_axe.y, rotate_axe.z, 0, -rotate_axe.x, -rotate_axe.y, rotate_axe.x, 0);
		glm::mat3x3 uxu = glm::mat3x3(rotate_axe.x * rotate_axe.x, rotate_axe.x * rotate_axe.y, rotate_axe.x * rotate_axe.z,
			rotate_axe.x * rotate_axe.y, rotate_axe.y * rotate_axe.y, rotate_axe.y * rotate_axe.z,
			rotate_axe.x * rotate_axe.z, rotate_axe.y * rotate_axe.z, rotate_axe.z * rotate_axe.z);
		float angle = 0.02 * M_PI;
		band_rotate = cos(angle) * glm::mat3(1.0f) + sin(angle) * ux + (1 - cos(angle)) * uxu;
		angle = -angle;
		band_rotate_reverse = cos(angle) * glm::mat3(1.0f) + sin(angle) * ux + (1 - cos(angle)) * uxu;


		rotate_axe = glm::vec3(0.0, 0.0, 1.0);
		ux = glm::mat3x3(0, -rotate_axe.z, rotate_axe.y, rotate_axe.z, 0, -rotate_axe.x, -rotate_axe.y, rotate_axe.x, 0);
		uxu = glm::mat3x3(rotate_axe.x * rotate_axe.x, rotate_axe.x * rotate_axe.y, rotate_axe.x * rotate_axe.z,
			rotate_axe.x * rotate_axe.y, rotate_axe.y * rotate_axe.y, rotate_axe.y * rotate_axe.z,
			rotate_axe.x * rotate_axe.z, rotate_axe.y * rotate_axe.z, rotate_axe.z * rotate_axe.z);
		angle = 0.2 * M_PI;
		capsule_rotate = cos(angle) * glm::mat3(1.0f) + sin(angle) * ux + (1 - cos(angle)) * uxu;
		angle = -angle;
		capsule_rotate_reverse = cos(angle) * glm::mat3(1.0f) + sin(angle) * ux + (1 - cos(angle)) * uxu;

	}


	void rotateCloth(double* body_center, OriMesh& mesh) {
		double move[3];
		memcpy(move, body_center, 24);
		glm::vec3 body_pos;
		for (int i = 0; i < mesh.vertices.size(); ++i) {
			body_pos = glm::vec3(mesh.vertices[i][0] - move[0], mesh.vertices[i][1] - move[1], mesh.vertices[i][2] - move[2]);// -body_center_y
			body_pos = capsule_rotate * body_pos;
			mesh.vertices[i][0] = body_pos.x + move[0];
			mesh.vertices[i][1] = body_pos.y + move[1];// +body_center_y;
			mesh.vertices[i][2] = body_pos.z + move[2];
		}

	}

	void initialMaterial(OriMesh& mesh, int type)
	{
		if (type == 1) {
			mesh.front_material.Kd[0] = 0.95;
			mesh.front_material.Kd[1] = 0.95;
			mesh.front_material.Kd[2] = 0.0;
			mesh.front_material.Ka[0] = 0.2;
			mesh.front_material.Ka[1] = 0.2;
			mesh.front_material.Ka[2] = 0.0;
			mesh.front_material.Ks[0] = 0.8;
			mesh.front_material.Ks[1] = 0.8;
			mesh.front_material.Ks[2] = 0.0;
			mesh.front_material.Tf[0] = 1.0;
			mesh.front_material.Tf[1] = 1.0;
			mesh.front_material.Tf[2] = 1.0;
			mesh.front_material.Ni = 1.0;
			mesh.front_material.Ns = 28.76;
		}
		else if (type == 2) {
			mesh.front_material.Kd[0] = 0.12;
			mesh.front_material.Kd[1] = 0.9;
			mesh.front_material.Kd[2] = 0.12;

			mesh.front_material.Ka[0] = 0.04;
			mesh.front_material.Ka[1] = 0.3;
			mesh.front_material.Ka[2] = 0.04;

			mesh.front_material.Ks[0] = 0.08;
			mesh.front_material.Ks[1] = 0.6;
			mesh.front_material.Ks[2] = 0.08;

			mesh.front_material.Tf[0] = 1.0;
			mesh.front_material.Tf[1] = 1.0;
			mesh.front_material.Tf[2] = 1.0;
			mesh.front_material.Ni = 1.0;
			mesh.front_material.Ns = 28.76;
		}
	}



	//int col_num = 40;
	//int row_num = 20;
	void indicesForTwoTriangle(OriMesh& mesh)
	{
		mesh.indices.push_back(0);
		mesh.indices.push_back(1);
		mesh.indices.push_back(2);
		//mesh.indices.push_back(0);
		//mesh.indices.push_back(3);
		//mesh.indices.push_back(4);
	}

	void indices(OriMesh& mesh, int col_num, int row_num)
	{
		int totalNum = mesh.vertices.size();
		for (int i = 0; i < totalNum; i++) {
			if ((i + 1) % col_num != 0 && i < totalNum - col_num) {
				mesh.indices.push_back(i);
				mesh.indices.push_back(i + col_num + 1);
				mesh.indices.push_back(i + 1);
				mesh.indices.push_back(i);
				mesh.indices.push_back(i + col_num);
				mesh.indices.push_back(i + col_num + 1);
			}
		}
		
		//mesh.virtual_face_indices = mesh.indices;
	}


	//void setVirtualIndices(OriMesh& mesh, int total_coarse_grid_num, int col_num, int row_num)
	//{
	//	mesh.virtual_face_indices.resize(total_coarse_grid_num);
	//	mesh.not_select_as_virtual_vertex.resize(total_coarse_grid_num);
	//	mesh.neighbor_virtual_vertex.resize(total_coarse_grid_num);
	//	mesh.barycentric.resize(total_coarse_grid_num);
	//	int interval[2] = { 4,3 };
	//	int coarse_grid_row_num;
	//	int coarse_grid_col_num;
	//	int upper_grid_col_num = col_num;
	//	int upper_grid_row_num = row_num;
	//	for (int i = 0; i < total_coarse_grid_num; ++i) {
	//		virtualIndices(mesh, upper_grid_col_num, upper_grid_row_num, i, coarse_grid_col_num, coarse_grid_row_num, interval[i]);
	//		upper_grid_col_num = coarse_grid_col_num;
	//		upper_grid_row_num = coarse_grid_row_num;
	//		std::cout << upper_grid_col_num << " " << upper_grid_row_num << std::endl;
	//	}
	//}

	//void virtualIndices(OriMesh& mesh, int col_num, int row_num, int coarse_grid_no, int &use_col_num,int &use_row_num, int interval)
	//{		
	//	if ((col_num-1) % interval == 0) {
	//		use_col_num = (col_num-1) / interval + 1;
	//	}
	//	else {
	//		use_col_num = (col_num - 1) / interval + 2;
	//	}
	//	if ((row_num-1) % interval == 0) {
	//		use_row_num = (row_num-1) / interval + 1;
	//	}
	//	else {
	//		use_row_num = (row_num-1) / interval + 2;
	//	}
	//	int virtual_total_num = use_col_num* use_row_num;
	//	mesh.virtual_face_indices[coarse_grid_no].reserve(2 * virtual_total_num);
	//	for (int i = 0; i < virtual_total_num; i++) {
	//		if ((i + 1) % use_col_num != 0 && i < virtual_total_num - use_col_num) {
	//			mesh.virtual_face_indices[coarse_grid_no].push_back(i);
	//			mesh.virtual_face_indices[coarse_grid_no].push_back(i + use_col_num + 1);
	//			mesh.virtual_face_indices[coarse_grid_no].push_back(i + 1);
	//			mesh.virtual_face_indices[coarse_grid_no].push_back(i);
	//			mesh.virtual_face_indices[coarse_grid_no].push_back(i + use_col_num);
	//			mesh.virtual_face_indices[coarse_grid_no].push_back(i + use_col_num + 1);
	//		}
	//	}
	//	std::vector<int>actual_index;
	//	actual_index.reserve(col_num*row_num);
	//	for (int i = 0; i < use_row_num-1; ++i) {
	//		for (int j = 0; j < use_col_num-1; ++j) {
	//			actual_index.push_back(col_num * (i * interval) + j * interval);
	//		}
	//		actual_index.push_back(col_num * (i * interval + 1) - 1);
	//	}
	//	for (int j = 0; j < use_col_num - 1; ++j) {
	//		actual_index.push_back(col_num * (row_num-1) + j * interval);
	//	}
	//	actual_index.push_back(col_num*row_num-1);
	//	
	//	for (int i = 0; i < mesh.virtual_face_indices[coarse_grid_no].size(); ++i) {
	//		mesh.virtual_face_indices[coarse_grid_no][i] = actual_index[mesh.virtual_face_indices[coarse_grid_no][i]];
	//	}
	//	mesh.not_select_as_virtual_vertex[coarse_grid_no].reserve(col_num * row_num);
	//	mesh.barycentric[coarse_grid_no].reserve(col_num * row_num);
	//	mesh.neighbor_virtual_vertex[coarse_grid_no].reserve(col_num * row_num);
	//	int virtual_col_no;
	//	int virtual_row_no;

	//	int chosen_row_num;
	//	int chosen_col_num;

	//	chosen_col_num = (use_col_num - 2) * interval;
	//	chosen_row_num = (use_row_num - 2) * interval;

	//	double barycentric0,barycentric1;
	//	for (int i = 0; i < chosen_row_num; ++i) {
	//		virtual_row_no = i / interval;
	//		if (i % interval == 0) {
	//			for (int j = 0; j < chosen_col_num; ++j) {
	//				virtual_col_no = j / interval;
	//				if (j % interval != 0) {
	//					mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i * col_num + j);
	//					mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no* use_col_num + virtual_col_no,(virtual_row_no + 1) * use_col_num + (virtual_col_no + 1), virtual_row_no * use_col_num + (virtual_col_no + 1)});
	//					mesh.barycentric[coarse_grid_no].push_back({ 1.0-(double)(j % interval) / (double)interval, 0.0,(double)(j % interval) / (double)interval});
	//				}
	//			}
	//			for (int j = chosen_col_num + 1; j < col_num-1; ++j) {
	//				virtual_col_no = j / interval;
	//				mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i * col_num + j);
	//				mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no ,(virtual_row_no + 1) * use_col_num + (virtual_col_no + 1), virtual_row_no * use_col_num+ (virtual_col_no + 1) });
	//				mesh.barycentric[coarse_grid_no].push_back({ 1.0-(double)(j % interval) / (double)(col_num-1- virtual_col_no * interval),0.0,(double)(j % interval) / (double)(col_num - 1 - virtual_col_no * interval) });
	//			}
	//		}
	//		else {
	//			for (int j = 0; j < chosen_col_num; ++j) {
	//				virtual_col_no = j / interval;
	//				mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i * col_num + j);
	//				if (j%interval >= i%interval) {						
	//					mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no,(virtual_row_no + 1) * use_col_num + (virtual_col_no + 1) , virtual_row_no * use_col_num + (virtual_col_no + 1) });
	//					barycentric0 = 1.0 - (double)(j % interval) / (double)interval;
	//					barycentric1 = 1.0 - barycentric0 - (double)(i % interval) / (double)interval;
	//					mesh.barycentric[coarse_grid_no].push_back({ barycentric0,1.0 - barycentric0 - barycentric1,barycentric1});
	//				}
	//				else {
	//					mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no, (virtual_row_no+1) * use_col_num + virtual_col_no , (virtual_row_no + 1) * use_col_num + (virtual_col_no + 1) });
	//					barycentric0 = 1.0 - (double)(i % interval) / (double)interval;
	//					barycentric1 = 1.0 - (double)(j % interval) / (double)interval-barycentric0;						
	//					mesh.barycentric[coarse_grid_no].push_back({ barycentric0,barycentric1,1.0 - barycentric0 - barycentric1 });
	//				}
	//			}
	//			if (col_num % interval == 1) {
	//				for (int j = chosen_col_num; j < col_num-1; ++j) {
	//					virtual_col_no = j / interval;
	//					double coe = (double)(col_num - 1 - virtual_col_no * interval) / (double)interval;
	//					mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i * col_num + j);
	//					if ((double)(j % interval) / (double)(i % interval) > coe) {
	//						mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no, (virtual_row_no + 1) * use_col_num + (virtual_col_no + 1) , virtual_row_no * use_col_num + (virtual_col_no + 1) });
	//						barycentric0 = 1.0 - (double)(j % interval) / (double)(col_num - 1 - virtual_col_no * interval);
	//						barycentric1 = 1.0 - barycentric0 - (double)(i % interval) / (double)interval;
	//						mesh.barycentric[coarse_grid_no].push_back({ barycentric0,1.0 - barycentric0 - barycentric1, barycentric1 });
	//					}
	//					else {
	//						mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no, (virtual_row_no + 1) * use_col_num + virtual_col_no,(virtual_row_no + 1) * use_col_num + (virtual_col_no + 1) });
	//						barycentric0 = 1.0 - (double)(i % interval) / (double)interval;
	//						barycentric1 = 1.0 - (double)(j % interval) / (double)(col_num - 1 - virtual_col_no * interval) - barycentric0;
	//						mesh.barycentric[coarse_grid_no].push_back({ barycentric0,barycentric1,1.0 - barycentric0 - barycentric1 });
	//					}
	//				}
	//				virtual_col_no = col_num / interval;
	//				mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i* col_num + col_num-1);
	//				mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no, virtual_row_no * use_col_num + (virtual_col_no - 1) , (virtual_row_no+1) * use_col_num + virtual_col_no});
	//				barycentric0 = 1.0 - (double)(i % interval) / (double)(interval);
	//				barycentric1 = 1.0 - barycentric0;
	//				mesh.barycentric[coarse_grid_no].push_back({ barycentric0,0.0, barycentric1 });
	//			}
	//			else {
	//				for (int j = chosen_col_num; j < col_num; ++j) {
	//					virtual_col_no = j / interval;
	//					double coe = (double)(col_num - 1 - virtual_col_no * interval) / (double)interval;
	//					mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i * col_num + j);
	//					if ((double)(j % interval) / (double)(i % interval) > coe) {
	//						mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no, (virtual_row_no + 1) * use_col_num + (virtual_col_no + 1) , virtual_row_no * use_col_num + (virtual_col_no + 1) });
	//						barycentric0 = 1.0 - (double)(j % interval) / (double)(col_num - 1 - virtual_col_no * interval);
	//						barycentric1 = 1.0 - barycentric0 - (double)(i % interval) / (double)interval;
	//						mesh.barycentric[coarse_grid_no].push_back({ barycentric0,1.0 - barycentric0 - barycentric1, barycentric1 });
	//					}
	//					else {
	//						mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no, (virtual_row_no + 1) * use_col_num + virtual_col_no,(virtual_row_no + 1) * use_col_num + (virtual_col_no + 1) });
	//						barycentric0 = 1.0 - (double)(i % interval) / (double)interval;
	//						barycentric1 = 1.0 - (double)(j % interval) / (double)(col_num - 1 - virtual_col_no * interval) - barycentric0;
	//						mesh.barycentric[coarse_grid_no].push_back({ barycentric0,barycentric1,1.0 - barycentric0 - barycentric1 });
	//					}
	//				}
	//			}
	//			
	//		}
	//	}
	//	for (int i = chosen_row_num; i < row_num-1; ++i) {
	//		virtual_row_no = i / interval;
	//		if (i % interval == 0) {
	//			for (int j = 0; j < chosen_col_num; ++j) {
	//				virtual_col_no = j / interval;
	//				if (j % interval != 0) {
	//					mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i * col_num + j);
	//					mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no,(virtual_row_no + 1) * use_col_num + (virtual_col_no + 1) , virtual_row_no * use_col_num + (virtual_col_no + 1) });
	//					barycentric0 = 1.0 - (double)(j % interval) / (double)interval;
	//					barycentric1 = 1.0 - barycentric0 - (double)(i % interval) / (double)(row_num-1- virtual_row_no*interval);
	//					mesh.barycentric[coarse_grid_no].push_back({ barycentric0,1.0 - barycentric0 - barycentric1, barycentric1});
	//				}
	//			}
	//			for (int j = chosen_col_num + 1; j < col_num-1; ++j) {
	//				virtual_col_no = j / interval;
	//				mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i * col_num + j);
	//				mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no,(virtual_row_no + 1) * use_col_num + (virtual_col_no + 1) , virtual_row_no * use_col_num + (virtual_col_no + 1) });
	//				barycentric0 = 1.0 - (double)(j % interval) / (double)(col_num - 1 - virtual_col_no * interval);
	//				barycentric1 = 1.0 - barycentric0 - (double)(i % interval) / (double)(row_num - 1 - virtual_row_no * interval);
	//				mesh.barycentric[coarse_grid_no].push_back({ barycentric0,1.0 - barycentric0 - barycentric1,barycentric1});
	//			}
	//		}
	//		else {
	//			double coe;
	//			for (int j = 0; j < chosen_col_num; ++j) {
	//				virtual_col_no = j / interval;
	//				coe= (double)interval / (double)(row_num - 1 - virtual_row_no * interval);
	//				mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i * col_num + j);
	//				if ((double)(j % interval) /(double) (i % interval)>coe) {
	//					mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no, (virtual_row_no + 1) * use_col_num + (virtual_col_no + 1), virtual_row_no * use_col_num + (virtual_col_no + 1) });
	//					barycentric0 = 1.0 - (double)(j % interval) / (double)interval;
	//					barycentric1 = 1.0 - barycentric0 - (double)(i % interval) / (double)(row_num - 1 - virtual_row_no * interval);
	//					mesh.barycentric[coarse_grid_no].push_back({ barycentric0,1.0 - barycentric0 - barycentric1,barycentric1 });
	//				}
	//				else {
	//					mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no,(virtual_row_no + 1) * use_col_num + virtual_col_no,  (virtual_row_no + 1) * use_col_num + (virtual_col_no + 1)});
	//					barycentric0 = 1.0 - (double)(i % interval) / (double)(row_num - 1 - virtual_row_no * interval);
	//					barycentric1 = 1.0 - (double)(j % interval) / (double)interval - barycentric0;						
	//					mesh.barycentric[coarse_grid_no].push_back({ barycentric0,barycentric1,1.0 - barycentric0 - barycentric1 });
	//				}
	//			}			
	//			if (col_num % interval == 1) {
	//				for (int j = chosen_col_num; j < col_num - 1; ++j) {
	//					virtual_col_no = j / interval;
	//					coe = (double)(col_num - 1 - virtual_col_no * interval) / (double)(row_num - 1 - virtual_row_no * interval);
	//					mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i* col_num + j);
	//					if ((double)(j % interval) / (double)(i % interval) >= coe) {
	//						mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no, (virtual_row_no + 1) * use_col_num + (use_col_num - 1) , virtual_row_no * use_col_num + (use_col_num - 1) });
	//						barycentric0 = 1.0 - (double)(j % interval) / (double)(col_num - 1 - virtual_col_no * interval);
	//						barycentric1 = 1.0 - barycentric0 - (double)(i % interval) / (double)(row_num - 1 - virtual_row_no * interval);
	//						mesh.barycentric[coarse_grid_no].push_back({ barycentric0,1.0 - barycentric0 - barycentric1,barycentric1 });
	//					}
	//					else {
	//						mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no,(virtual_row_no + 1) * use_col_num + virtual_col_no,(virtual_row_no + 1) * use_col_num + (use_col_num - 1) });
	//						barycentric0 = 1.0 - (double)(i % interval) / (double)(row_num - 1 - virtual_row_no * interval);
	//						barycentric1 = 1.0 - (double)(j % interval) / (double)(col_num - 1 - virtual_col_no * interval) - barycentric0;
	//						mesh.barycentric[coarse_grid_no].push_back({ barycentric0,barycentric1,1.0 - barycentric0 - barycentric1 });
	//					}
	//				}
	//				virtual_col_no = col_num / interval;
	//				mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i * col_num + col_num - 1);
	//				mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no, virtual_row_no * use_col_num + (virtual_col_no - 1) , (virtual_row_no + 1) * use_col_num + virtual_col_no });
	//				barycentric0 = 1.0 - (double)(i % interval) / (double)(row_num - 1 - virtual_row_no * interval);
	//				barycentric1 = 1.0 - barycentric0;
	//				mesh.barycentric[coarse_grid_no].push_back({ barycentric0,0.0, barycentric1 });
	//			}
	//			else {
	//				for (int j = chosen_col_num; j < col_num; ++j) {
	//					virtual_col_no = j / interval;
	//					coe = (double)(col_num - 1 - virtual_col_no * interval) / (double)(row_num - 1 - virtual_row_no * interval);
	//					mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back(i * col_num + j);

	//					if ((double)(j % interval) / (double)(i % interval) >= coe) {

	//						mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no, (virtual_row_no + 1) * use_col_num + (use_col_num - 1) , virtual_row_no * use_col_num + (use_col_num - 1) });
	//						barycentric0 = 1.0 - (double)(j % interval) / (double)(col_num - 1 - virtual_col_no * interval);
	//						barycentric1 = 1.0 - barycentric0 - (double)(i % interval) / (double)(row_num - 1 - virtual_row_no * interval);
	//						mesh.barycentric[coarse_grid_no].push_back({ barycentric0,1.0 - barycentric0 - barycentric1,barycentric1 });
	//					}
	//					else {

	//						mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ virtual_row_no * use_col_num + virtual_col_no,(virtual_row_no + 1) * use_col_num + virtual_col_no,(virtual_row_no + 1) * use_col_num + (use_col_num - 1) });
	//						barycentric0 = 1.0 - (double)(i % interval) / (double)(row_num - 1 - virtual_row_no * interval);
	//						barycentric1 = 1.0 - (double)(j % interval) / (double)(col_num - 1 - virtual_col_no * interval) - barycentric0;
	//						mesh.barycentric[coarse_grid_no].push_back({ barycentric0,barycentric1,1.0 - barycentric0 - barycentric1 });
	//					}
	//				}
	//			}
	//		}
	//	}
	//	for (int j = 0; j < chosen_col_num; ++j) {
	//		virtual_col_no = j / interval;
	//		if (j % interval != 0) {
	//			mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back((row_num - 1) * col_num + j);
	//			mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ (use_row_num-1) * use_col_num + virtual_col_no, (use_row_num - 1) * use_col_num + (virtual_col_no + 1) , (use_row_num - 1) * use_col_num + (virtual_col_no + 1)});
	//			barycentric0 = 1.0 - (double)(j % interval) / (double)interval;
	//			barycentric1 = 1.0 - barycentric0;
	//			mesh.barycentric[coarse_grid_no].push_back({ barycentric0,0.0,barycentric1});
	//		}
	//	}
	//	for (int j = chosen_col_num + 1; j < col_num-1; ++j) {
	//		virtual_col_no = j / interval;
	//		mesh.not_select_as_virtual_vertex[coarse_grid_no].push_back((row_num - 1) * col_num + j);
	//		mesh.neighbor_virtual_vertex[coarse_grid_no].push_back({ (use_row_num - 1) * use_col_num + virtual_col_no, (use_row_num - 1) * use_col_num + (use_col_num - 1)  , (use_row_num - 1) * use_col_num + (use_col_num - 1)});
	//		barycentric0 = 1.0 - (double)(j % interval) / (double)(col_num - 1 - virtual_col_no * interval);
	//		barycentric1 = 1.0 - barycentric0;
	//		mesh.barycentric[coarse_grid_no].push_back({ barycentric0,0.0 ,barycentric1});
	//	}
	//	//for (int i = 0; i < mesh.not_select_as_virtual_vertex.size(); ++i) {
	//	//	//if (mesh.not_select_as_virtual_vertex[i] == 24) {
	//	//		std::cout << "+++"<<mesh.not_select_as_virtual_vertex[i] << std::endl;
	//	//	//}
	//	//}
	//	//for (int i = 0; i < mesh.not_select_as_virtual_vertex.size(); ++i) {
	//	//	std::cout << i << " " << mesh.neighbor_virtual_vertex[i][0] << " " << mesh.neighbor_virtual_vertex[i][1] << " " << mesh.neighbor_virtual_vertex[i][2]<<" "<< mesh.barycentric[i][0]<<" "<< mesh.barycentric[i][1]<<" "<< mesh.barycentric[i][2]<< std::endl;
	//	//}
	//}

	void meshIndices(OriMesh& mesh, int col_num, int row_num)
	{

		for (int i = 0; i < row_num - 1; i++)
		{
			for (int j = 0; j < col_num - 1; j++)
			{
				int helper = 0;
				if (i % 2 == j % 2)
					helper = 1;
				mesh.indices.emplace_back(i * col_num + j);
				mesh.indices.emplace_back(i * col_num + j + 1);
				mesh.indices.emplace_back((i + 1) * col_num + j + helper);
				mesh.indices.emplace_back((i + 1) * col_num + j + 1);
				mesh.indices.emplace_back((i + 1) * col_num + j);
				mesh.indices.emplace_back(i * col_num + j + 1 - helper);
			}
		}
		//assert(col_num % 2 == 0);
		//assert(row_num % 2 == 0);
		//int totalNum = mesh.vertices.size();
		//int index;
		//for (int j = 0; j < row_num - 1; j += 2) {
		//	for (int i = 0; i < col_num - 1; i += 2) {
		//		index = col_num * j + i;
		//		mesh.indices.push_back(index);
		//		mesh.indices.push_back(index + col_num);
		//		mesh.indices.push_back(index + 1);
		//		mesh.indices.push_back(index + 1);
		//		mesh.indices.push_back(index + col_num);
		//		mesh.indices.push_back(index + col_num + 1);
		//		mesh.indices.push_back(index + 1);
		//		mesh.indices.push_back(index + col_num + 1);
		//		mesh.indices.push_back(index + col_num + 2);
		//		mesh.indices.push_back(index + col_num);
		//		mesh.indices.push_back(index + 2 * col_num + 1);
		//		mesh.indices.push_back(index + col_num + 1);
		//		mesh.indices.push_back(index + col_num + 1);
		//		mesh.indices.push_back(index + 2 * col_num + 1);
		//		mesh.indices.push_back(index + col_num + 2);
		//		mesh.indices.push_back(index + col_num);
		//		mesh.indices.push_back(index + 2 * col_num);
		//		mesh.indices.push_back(index + 2 * col_num + 1);
		//		mesh.indices.push_back(index + col_num + 2);
		//		mesh.indices.push_back(index + 2 * col_num + 1);
		//		mesh.indices.push_back(index + 2 * col_num + 2);
		//		mesh.indices.push_back(index + 1);
		//		mesh.indices.push_back(index + col_num + 2);
		//		mesh.indices.push_back(index + 2);
		//	}
		//}

	}

	void cylinderIndices(OriMesh& mesh, int first_half_sphere, int total_num, int col_num, int row_num)
	{

		for (int i = 0; i < total_num; i++) {
			if ((i + 1) % col_num != 0 && i < total_num - col_num) {
				mesh.indices.push_back(i + first_half_sphere);
				mesh.indices.push_back(i + 1 + first_half_sphere);
				mesh.indices.push_back(i + col_num + 1 + first_half_sphere);
				mesh.indices.push_back(i + first_half_sphere);
				mesh.indices.push_back(i + col_num + 1 + first_half_sphere);
				mesh.indices.push_back(i + col_num + first_half_sphere);
			}
		}
		for (int i = 0; i < row_num - 1; ++i) {
			mesh.indices.push_back(col_num * (i + 1) + first_half_sphere);
			mesh.indices.push_back(col_num * (i + 1) - 1 + first_half_sphere);
			mesh.indices.push_back(col_num * i + first_half_sphere);
			mesh.indices.push_back(col_num * (i + 1) - 1 + first_half_sphere);
			mesh.indices.push_back(col_num * (i + 1) + first_half_sphere);
			mesh.indices.push_back(col_num * (i + 2) - 1 + first_half_sphere);
		}
	}

	void sphereIndices(int row, int col, OriMesh& mesh)
	{
		std::vector<int>globe_indices;
		for (int j = 1; j < row - 1; ++j) {
			for (int i = 0; i < col - 1; ++i) {
				globe_indices.push_back(col * (j - 1) + i);
				globe_indices.push_back(col * (j - 1) + i + 1);
				globe_indices.push_back(col * j + i);
				globe_indices.push_back(col * j + i);
				globe_indices.push_back(col * (j - 1) + i + 1);
				globe_indices.push_back(col * j + i + 1);
			}
			globe_indices.push_back(col * (j - 1) + col - 1);
			globe_indices.push_back(col * (j - 1));
			globe_indices.push_back(col * j + col - 1);
			globe_indices.push_back(col * (j - 1));
			globe_indices.push_back(col * j);
			globe_indices.push_back(col * j + col - 1);
		}
		int vertex_num = (row - 1) * col;
		int pre_vertex_num = (row - 2) * col;
		for (int i = 0; i < col - 1; ++i) {
			globe_indices.push_back(vertex_num);
			globe_indices.push_back(i + 1);
			globe_indices.push_back(i);

			globe_indices.push_back(pre_vertex_num + i + 1);
			globe_indices.push_back(vertex_num + 1);
			globe_indices.push_back(pre_vertex_num + i);

		}
		globe_indices.push_back(vertex_num);
		globe_indices.push_back(0);
		globe_indices.push_back(col - 1);
		globe_indices.push_back(pre_vertex_num);
		globe_indices.push_back(vertex_num + 1);
		globe_indices.push_back(pre_vertex_num + col - 1);
		mesh.indices.insert(mesh.indices.end(), globe_indices.begin(), globe_indices.end());
	}


	void globeIndices(int row, int col, OriMesh& mesh, int cylinder_vertex_num)
	{
		std::vector<int>globe_indices;
		for (int j = 1; j < row - 1; ++j) {
			for (int i = 0; i < col - 1; ++i) {
				globe_indices.push_back(col * (j - 1) + i);
				globe_indices.push_back(col * (j - 1) + i + 1);
				globe_indices.push_back(col * j + i);
				globe_indices.push_back(col * j + i);
				globe_indices.push_back(col * (j - 1) + i + 1);
				globe_indices.push_back(col * j + i + 1);
			}
			globe_indices.push_back(col * (j - 1) + col - 1);
			globe_indices.push_back(col * (j - 1));
			globe_indices.push_back(col * j + col - 1);
			globe_indices.push_back(col * (j - 1));
			globe_indices.push_back(col * j);
			globe_indices.push_back(col * j + col - 1);
		}
		//	mesh.indices.insert(mesh.indices.end(), globe_indices.begin(), globe_indices.end());
		for (int i = 0; i < globe_indices.size() / 3; ++i) {
			mesh.indices.push_back(globe_indices[3 * i] + 1);
			mesh.indices.push_back(globe_indices[3 * i + 1] + 1);
			mesh.indices.push_back(globe_indices[3 * i + 2] + 1);
			mesh.indices.push_back(globe_indices[3 * i] + cylinder_vertex_num);
			mesh.indices.push_back(globe_indices[3 * i + 1] + cylinder_vertex_num);
			mesh.indices.push_back(globe_indices[3 * i + 2] + cylinder_vertex_num);
		}

		globe_indices.clear();
		for (int i = 0; i < col - 1; ++i) {
			globe_indices.push_back(0);
			globe_indices.push_back(i + 2);
			globe_indices.push_back(i + 1);
		}
		globe_indices.push_back(0);
		globe_indices.push_back(1);
		globe_indices.push_back(col);

		mesh.indices.insert(mesh.indices.end(), globe_indices.begin(), globe_indices.end());
		int vertex_num = (row - 1) * col;
		globe_indices.clear();
		for (int i = 0; i < col - 1; ++i) {
			globe_indices.push_back((row - 2) * col + i);
			globe_indices.push_back((row - 2) * col + i + 1);
			globe_indices.push_back(vertex_num);
		}
		globe_indices.push_back((row - 1) * col - 1);
		globe_indices.push_back((row - 2) * col);
		globe_indices.push_back(vertex_num);
		for (int i = 0; i < globe_indices.size() / 3; ++i) {
			mesh.indices.push_back(globe_indices[3 * i] + cylinder_vertex_num);
			mesh.indices.push_back(globe_indices[3 * i + 1] + cylinder_vertex_num);
			mesh.indices.push_back(globe_indices[3 * i + 2] + cylinder_vertex_num);
		}
	}


	void floorMaterial(OriMesh& mesh)
	{
		mesh.front_material.Kd[0] = 0.8;
		mesh.front_material.Kd[1] = 0.8;
		mesh.front_material.Kd[2] = 0.8;
		mesh.front_material.Ka[0] = 0.2136;
		mesh.front_material.Ka[1] = 0.3136;
		mesh.front_material.Ka[2] = 0.1136;
		mesh.front_material.Ks[0] = 0.0136;
		mesh.front_material.Ks[1] = 0.0136;
		mesh.front_material.Ks[2] = 0.0136;
		mesh.front_material.Tf[0] = 1.0;
		mesh.front_material.Tf[1] = 1.0;
		mesh.front_material.Tf[2] = 1.0;
		mesh.front_material.Ni = 1.0;
		mesh.front_material.Ns = 28.76;
	}
};


#endif // !CREATE_MESH_H



#pragma once

