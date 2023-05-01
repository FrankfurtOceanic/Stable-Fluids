#include <iostream>
#include <vector>
#include <chrono>
#include <limits>
#include <cmath>

#define GLEW_STATIC

#define DPI_SCALE 1.0f
#define MSAA_SAMPLES 16

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include "ShaderProgram.h"
#include "Buffer.h"
#include "GridSolver.h"
#include "Vertex.h"
#include "RecordVideo.hpp"

int winWidth = 1280, winHeight = 960;
bool simulating = false;
int fps = 30; // dipslay FPS (not simulation)

ImageRecorder imageRecorder;
int speed = 4;
bool recordFrames = false; // toggled by keyboard
bool recordFrame = false; // set to record when stepped in display`
bool showFilament = true;
bool showGrid = true;
bool showVel = false;
float velFactor = 1000.f;
float srcTemp = 50.f;
float pointedTemp = std::numeric_limits<double>::quiet_NaN();

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::shared_ptr;
using glm::mat4;
using glm::vec3;
using std::chrono::steady_clock;

void keyboardCB(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (action == GLFW_PRESS) {
		switch (key) {
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose(window, GLFW_TRUE);
			break;
		case GLFW_KEY_S:
			step();
			//updateGrid();
			if (recordFrames)
				recordFrame = true;
			break;
		case GLFW_KEY_R:
			initGrid();
			break;
		case GLFW_KEY_V:
			recordFrames = !recordFrames;
			break;
		case GLFW_KEY_SPACE:
			simulating = !simulating;
			break;
		}
	}
}

void getIJ(double x, double y, int &i, int &j) {
	y = winHeight - y;
	double x_offset, y_offset;
	x_offset = (winWidth - frameWidth) / 2;
	y_offset = (winHeight - frameHeight) / 2;
	j = (int)((x - x_offset) / cellSize) + 1;
	i = (int)((y - y_offset) / cellSize) + 1;
}

void cursorPosCB(GLFWwindow* window, double xpos, double ypos) {
	static double xprv = 0, yprv = 0;
	if (!ImGui::GetIO().WantCaptureMouse) {
		int i, j;
		getIJ(xpos, ypos, i, j);
		if (0 < i && i <= rows && 0 < j && j <= cols)
			pointedTemp = tpSrc[i][j];
		else
			pointedTemp = std::numeric_limits<double>::quiet_NaN();
		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
			double dx = xpos - xprv;
                        double dy = ypos - yprv;
			if (0 < i && i <= rows && 0 < j && j <= cols) {
				double vx = velFactor * +dx;
				double vy = velFactor * -dy;
				vxSrc[i][j] = vx;
				vySrc[i][j] = vy;
			}
		}
	}
	xprv = xpos;
	yprv = ypos;
}

void mouseButtonCB(GLFWwindow* window, int button, int action, int mods)
{
	if (!ImGui::GetIO().WantCaptureMouse) {
		if (button == GLFW_MOUSE_BUTTON_RIGHT) {
			if (action == GLFW_PRESS) {
				double x, y;
				int i, j;
				glfwGetCursorPos(window, &x, &y);
				getIJ(x, y, i, j);
				//cout << "row : " << i << "  col : " << j << endl; debugging cells
				if (0 < i && i <= rows && 0 < j && j <= cols) {
					if (tpSrc[i][j] == 0)
						tpSrc[i][j] = srcTemp;
					else
						tpSrc[i][j] = 0;
				}
			}
		} else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
			if (action == GLFW_PRESS) {
				double x, y;
				glfwGetCursorPos(window, &x, &y);
				x = x - (winWidth-frameWidth)/2;
				y = (frameHeight+winHeight)/2 - y;
				if (0 < x && x < frameWidth && 0 < y && y < frameHeight)
					addFilament(x, y);
			}
		}
	}
}

void scrollCB(GLFWwindow* window, double xoffset, double yoffset) {
	if (!ImGui::GetIO().WantCaptureMouse) {
		double x, y;
		int i, j;
		glfwGetCursorPos(window, &x, &y);
		getIJ(x, y, i, j);
		if (0 < i && i <= rows && 0 < j && j <= cols) {
			tpSrc[i][j] += 0.1 * yoffset;
			pointedTemp = tpSrc[i][j];
		}
	}
}

int main(void) {
	GLFWwindow* window;

	if (!glfwInit()) {
		cerr << "glfwInit() failed" << endl;
		exit(EXIT_FAILURE);
	}

	// GL 4.0
	// GLSL 400
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

	glfwWindowHint(GLFW_SAMPLES, MSAA_SAMPLES);
	window = glfwCreateWindow(winWidth, winHeight, "Fluid", NULL, NULL);
	if (!window) {
		cerr << "glfwCreateWindow() failed" << endl;
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwSetKeyCallback(window, keyboardCB);
	glfwSetCursorPosCallback(window, cursorPosCB);
	glfwSetMouseButtonCallback(window, mouseButtonCB);
	glfwSetScrollCallback(window, scrollCB);
	glfwMakeContextCurrent(window);

	GLenum errn = glewInit();
	if (GLEW_OK != errn) {
		cerr << "glewInit() failed: " << glewGetErrorString(errn) << endl;
		exit(EXIT_FAILURE);
	}

	glEnable(GL_MULTISAMPLE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	ShaderProgram shader("../resources/shaders/simple.vs", "../resources/shaders/simple.fs");

	vector<size_t> attrSizes{ 2, 4 };

	initGrid();

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	// Setup Dear ImGui style
	ImGui::StyleColorsDark();
	//ImGui::StyleColorsLight();

	// Setup Platform/Renderer backends
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 400");
#ifdef DPI_SCALE
	ImFontConfig font_cfg;
	font_cfg.SizePixels = (int)(13 * DPI_SCALE);
	ImGui::GetIO().Fonts->AddFontDefault(&font_cfg);
	ImGui::GetStyle().ScaleAllSizes(DPI_SCALE);
#endif

	Vertex frameVertices[2 * (rows + cols + 2)];
	size_t k = 0;
	const float c = 0.2f;
	for (size_t i = 0; i <= cols; ++i) {
		frameVertices[k++] = Vertex((float)i * cellSize, (float)0, 0.2f, 1.0f, 0.2f, c);
		frameVertices[k++] = Vertex((float)i * cellSize, (float)frameHeight, 1.0f, 1.0f, 1.0f, c);
	}
	for (size_t i = 0; i <= rows; ++i) {
		frameVertices[k++] = Vertex((float)0, (float)i * cellSize, 0.2f, 1.0f, 0.2f, c);
		frameVertices[k++] = Vertex((float)frameWidth, (float)i * cellSize, 1.0f, 1.0f, 1.0f, c);
	}
	Buffer frameBuffer(attrSizes, (GLfloat*)frameVertices, sizeof(frameVertices));
	const float dotRad = 2.5f;
	GLfloat rd[] = {
		-dotRad, -dotRad, 1, 0, 0, 1,
		 dotRad, -dotRad, 1, 0, 0, 1,
		 dotRad,  dotRad, 1, 0, 0, 1,
		-dotRad, -dotRad, 1, 0, 0, 1,
		 dotRad,  dotRad, 1, 0, 0, 1,
		-dotRad,  dotRad, 1, 0, 0, 1
	};
	Buffer redDotBuffer(attrSizes, rd, sizeof(rd));
	GLfloat bd[] = {
		-dotRad, -dotRad, 0, 0, 1, 1,
		 dotRad, -dotRad, 0, 0, 1, 1,
		 dotRad,  dotRad, 0, 0, 1, 1,
		-dotRad, -dotRad, 0, 0, 1, 1,
		 dotRad,  dotRad, 0, 0, 1, 1,
		-dotRad,  dotRad, 0, 0, 1, 1
	};
	Buffer blueDotBuffer(attrSizes, bd, sizeof(bd));

	auto t_now = steady_clock::now();
	while (!glfwWindowShouldClose(window)) {

		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		ImGui::Begin("Sim Setting");
		ImGui::Text("Shortcut:\n\tspace - run/stop\n\ts     - step\n\tr     - restart\n\tv     - record\n\tESC   - exit");
		ImGui::Text("\tright button     - add/delete heat source\n\tmiddle button    - add filament\n\tleft button drag - add force\n\tscrolling        - temperature control");
		ImGui::Text("simTime: %6.2f", simTime);
		ImGui::Text("srcTemp: %6.2f", pointedTemp);
		ImGui::Text("avgTemp: %6.2f", (float)getReferenceTemperature(cur_tp));
		ImGui::SliderInt("speed", &speed, 1, 10);
		ImGui::Checkbox("recordFrames", &recordFrames);
		ImGui::SliderFloat("dt", &dt, 0.0001f, 0.1f);
		ImGui::SliderFloat("dx", &dx, 0.01f, 5.0f);
		ImGui::SliderFloat("nv", &nv, 0.f, 20.f);
		ImGui::SliderFloat("kappa", &kappa, 0.f, 5.f);
		ImGui::SliderFloat("buoyancy", &buoyancy, 0.f, 5.f);
		ImGui::SliderFloat("srcTemp", &srcTemp, -500.f, 500.f);
		ImGui::SliderFloat("velFactor", &velFactor, 10.f, 5000.f);
		ImGui::SliderFloat("maxAge", &maxAge, 1.0f, 50.0f);
		ImGui::SliderFloat("refinementThreshold", &refinementThreshold, 1.0f, 50.0f);
		ImGui::SliderInt("iterations", &iterations, 1, 100);
		ImGui::End();

		ImGui::Begin("Viz Setting");
		ImGui::Checkbox("showFilament", &showFilament);
		ImGui::Checkbox("showGrid", &showGrid);
		ImGui::Checkbox("showVel", &showVel);
		ImGui::SliderFloat("velScale", &velScale, 1.f, 50.f);
		ImGui::End();

		ImGui::Render();

		glfwGetFramebufferSize(window, &winWidth, &winHeight);

		glViewport(0, 0, winWidth, winHeight);
		glClear(GL_COLOR_BUFFER_BIT);

		mat4 pv = glm::ortho(-(float)winWidth / 2, (float)winWidth / 2, -(float)winHeight / 2, (float)winHeight / 2); // Note that glm::ortho is a template function. Must cast the arguments to float to use the function with the correct type!
		mat4 m = glm::translate(vec3(-frameWidth/2, -frameHeight/2, 0.f));
		shader.use();
		shader.set_mat4("MVP", pv * m);

		if (simulating && std::chrono::duration_cast<std::chrono::milliseconds>(steady_clock::now() - t_now).count() > 1000 / fps) {
			for (int i = 0; i < speed; ++i)
				step();
			updateGrid();
			t_now = steady_clock::now();
			if (recordFrames)
				recordFrame = true;
		}
		//updateGrid();

		Buffer gridBuffer(attrSizes, (GLfloat*)gridVertices, sizeof(gridVertices), gridIndices, sizeof(gridIndices));
		Buffer velBuffer(attrSizes, (GLfloat*)velVertices, sizeof(velVertices));
		gridBuffer.draw(GL_TRIANGLES);

		if (showVel)
			velBuffer.draw(GL_LINES);

		if (showFilament) {
			for (shared_ptr<vector<Vertex>> &filament: filaments) {
				Buffer buffer(attrSizes, (GLfloat*)filament->data(), filament->size() * sizeof(Vertex));
				buffer.draw(GL_LINE_STRIP);
			}
		}

		if (showGrid)
			frameBuffer.draw(GL_LINES);

		for (size_t i = 1; i <= rows; ++i)
			for (size_t j = 1; j <= cols; ++j)
				if (tpSrc[i][j] > 0) {
					mat4 m = glm::translate(vec3(-frameWidth/2 + (j-0.5) * cellSize, -frameHeight/2 + (i-0.5) * cellSize, 0.f));
					shader.set_mat4("MVP", pv * m);
					redDotBuffer.draw(GL_TRIANGLES);
				} else if (tpSrc[i][j] < 0) {
					mat4 m = glm::translate(vec3(-frameWidth/2 + (j-0.5) * cellSize, -frameHeight/2 + (i-0.5) * cellSize, 0.f));
					shader.set_mat4("MVP", pv * m);
					blueDotBuffer.draw(GL_TRIANGLES);
				}

		if (recordFrame)
		{
			recordFrame = false;
			imageRecorder.writeCurrentFrameToFile(window);
		}

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);

	glfwTerminate();
	exit(EXIT_SUCCESS);
}
