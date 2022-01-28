#ifndef SAVE_IMAGE_H
#define SAVE_IMAGE_H
#include <io.h>
#include <direct.h>
#include"../basic/global.h"
//#include"../external/FreeImage.h"
#define STB_IMAGE_IMPLEMENTATION
#include "../external/stb-master/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../external/stb-master/stb_image_write.h"
class SaveImage 
{
public:
	SaveImage()
	{
		std::string prefix = "./screen_record/";
		if (_access(prefix.c_str(), 0) == -1)
			_mkdir(prefix.c_str());

		total_bit_size= 3 * SCR_WIDTH * SCR_HEIGHT;
		pixels = new BYTE[total_bit_size];
		initialHeader();	
	}

	//void outputImage()
	//{
	//	//glPixelStorei(GL_PACK_ALIGNMENT, 1);

	//	glReadPixels(0,0, SCR_WIDTH, SCR_HEIGHT, GL_BGR, GL_UNSIGNED_BYTE, pixels);
	//	//glReadPixels(0,0, SCR_WIDTH, SCR_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	//	FIBITMAP* image = FreeImage_ConvertFromRawBits(pixels, SCR_WIDTH, SCR_HEIGHT, 3 * SCR_WIDTH, 24, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, false);
	//	std::string name = "./screen_record/SCR_" + std::to_string(num);
	//	name+=".jpeg";
	//	FreeImage_Save(FIF_JPEG, image, name.c_str(),0);
	//	//cout << "Image Saved."<<endl;
	//	//    // Free resources
	//	FreeImage_Unload(image);
	//	num++;
	//	//delete[] pixels;
	//}

	void outputImage(int num)
	{
		glPixelStorei(GL_PACK_ALIGNMENT, 1);

		//glReadPixels(0, 0, SCR_WIDTH, SCR_HEIGHT, GL_BGR, GL_UNSIGNED_BYTE, pixels);
		glReadPixels(0,0, SCR_WIDTH, SCR_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, pixels);		
		std::string name = "./screen_record/SCR_" + std::to_string(num);
		name += ".png";
		stbi_flip_vertically_on_write(true);
		stbi_write_png(name.c_str(), SCR_WIDTH, SCR_HEIGHT, 3, pixels, SCR_WIDTH * 3);
		//name += ".bmp";
		//FILE* fp = fopen(name.c_str(), "wb+");
		//if (fp == NULL) {
		//	std::cout << "Could not open file: %s" << std::endl;
		//}
		//fwrite(header, sizeof(unsigned char), 54, fp);
		//fwrite(pixels, 1, total_bit_size, fp);
		//fclose(fp);
		//delete[] pixels;
	}
private:
	
	BYTE* pixels;

	unsigned char header[54];

	int total_bit_size;

	void writBmpheader(unsigned char* inputImg, int offset, int bytes, int value) {
		int i;
		for (i = 0; i < bytes; i++)
			inputImg[offset + i] = (value >> (i << 3)) & 0xFF;
	}
	void initialHeader()
	{
		header[0] = 'B';
		header[1] = 'M';
		writBmpheader(header, 2, 4, total_bit_size + 54); //whole file size
		writBmpheader(header, 0xA, 4, 54); //offset before bitmap raw data
		writBmpheader(header, 0xE, 4, 40); //length of bitmap info header
		writBmpheader(header, 0x12, 4, SCR_WIDTH); //width
		writBmpheader(header, 0x16, 4, SCR_HEIGHT); //height
		writBmpheader(header, 0x1A, 2, 1);
		writBmpheader(header, 0x1C, 2, 24); //bit per pixel
		writBmpheader(header, 0x1E, 4, 0); //compression
		writBmpheader(header, 0x22, 4, total_bit_size); //size of bitmap raw data
		for (int i = 0x26; i < 0x36; i++)
			header[i] = 0;
	}

};
#endif // !SCENE_H
#pragma once
