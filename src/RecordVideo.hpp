#pragma once
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cassert>
#include <vector>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

using namespace std;

class ImageRecorder {
public:

	string path = "../out/";
	int frameNumber = 0;
	string filename = "";
	vector<unsigned char> pixels;
	vector<unsigned char> rowData;

	void writeCurrentFrameToFile(GLFWwindow* window) {
		// get the image
		int comp = 3; // RGB
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		if (pixels.size() != width * height * comp) {
			pixels.resize(width * height * comp);
			rowData.resize(width * comp);
		}
		glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, &pixels[0]);
		// gross, need to flip the buffer
		for (int i = 0; i < height / 2; i++) {
			int r1 = i * width * comp;;
			int r2 = (height - 1 - i) * width * comp;
			memcpy(&rowData[0], &pixels[r1], rowData.size());
			memcpy(&pixels[r1], &pixels[r2], rowData.size());
			memcpy(&pixels[r2], &rowData[0], rowData.size());
		}
		// make the filename
		std::ostringstream oss;
		oss << path << "image" << std::setfill('0') << std::setw(5) << frameNumber++ << ".png";
		filename = oss.str();
		// write to file
		int stride_in_bytes = width * comp * sizeof(unsigned char);
		int rc = stbi_write_png(filename.c_str(), width, height, comp, &pixels[0], stride_in_bytes);
		if (rc == 0) {
			cout << "Couldn't write to " << filename << endl;
		}
	}
};