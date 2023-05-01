#include <iostream>
#include <fstream>
#include <sstream>
#include <GL/glew.h>
#include <glm/gtc/type_ptr.hpp>

#include "ShaderProgram.h"

using namespace std;

string readFile(string fname);

ShaderProgram::ShaderProgram(string vsFilename, string fsFilename) {
	int status;
	const size_t infoLen = 1024;
	char info[infoLen];

	// Vertex shader
	vShader = glCreateShader(GL_VERTEX_SHADER);
	string v_prg = readFile(vsFilename);
	const char* vs = v_prg.c_str();
	glShaderSource(vShader, 1, &vs, NULL);
	glCompileShader(vShader);
	glGetShaderiv(vShader, GL_COMPILE_STATUS, &status);
	if (!status) {
		glGetShaderInfoLog(vShader, infoLen, NULL, info);
		cerr << "Compile shader " << vsFilename << " failed: " << info << endl;
		exit(EXIT_FAILURE);
	}

	// Fragment shader
	fShader = glCreateShader(GL_FRAGMENT_SHADER);
	string f_prg = readFile(fsFilename);
	const char* fs = f_prg.c_str();
	glShaderSource(fShader, 1, &fs, NULL);
	glCompileShader(fShader);
	glGetShaderiv(fShader, GL_COMPILE_STATUS, &status);
	if (!status) {
		glGetShaderInfoLog(fShader, infoLen, NULL, info);
		cerr << "Compile shader " << fsFilename << " failed: " << info << endl;
		exit(EXIT_FAILURE);
	}

	program = glCreateProgram();
	glAttachShader(program, vShader);
	glAttachShader(program, fShader);
	glLinkProgram(program);
	glGetProgramiv(program, GL_LINK_STATUS, &status);
	if (!status) {
		glGetProgramInfoLog(program, infoLen, NULL, info);
		cerr << "Link program failed: " << info << endl;
	}
	glDeleteShader(vShader);
	glDeleteShader(fShader);
}

void ShaderProgram::use() {
	glUseProgram(program);
}

void ShaderProgram::set_int(string name, int i) {
	glUniform1i(glGetUniformLocation(program, name.c_str()), i);
}

void ShaderProgram::set_float(string name, float f) {
	glUniform1i(glGetUniformLocation(program, name.c_str()), f);
}

void ShaderProgram::set_vec2(string name, glm::vec2 v) {
	glUniform2fv(glGetUniformLocation(program, name.c_str()), 1, glm::value_ptr(v));
}

void ShaderProgram::set_vec3(string name, glm::vec3 v) {
	glUniform3fv(glGetUniformLocation(program, name.c_str()), 1, glm::value_ptr(v));
}

void ShaderProgram::set_vec4(string name, glm::vec4 v) {
	glUniform4fv(glGetUniformLocation(program, name.c_str()), 1, glm::value_ptr(v));
}

void ShaderProgram::set_mat2(string name, glm::mat2 v) {
	glUniformMatrix2fv(glGetUniformLocation(program, name.c_str()), 1, GL_FALSE, glm::value_ptr(v));
}

void ShaderProgram::set_mat3(string name, glm::mat3 v) {
	glUniformMatrix3fv(glGetUniformLocation(program, name.c_str()), 1, GL_FALSE, glm::value_ptr(v));
}

void ShaderProgram::set_mat4(string name, glm::mat4 v) {
	glUniformMatrix4fv(glGetUniformLocation(program, name.c_str()), 1, GL_FALSE, glm::value_ptr(v));
}

string readFile(string fname) {
	ifstream fin(fname);
	if (fin.good()) {
		stringstream ss;
		ss << fin.rdbuf();
		return ss.str();
	} else {
		cerr << "Cannot open " << fname << endl;
		exit(EXIT_FAILURE);
	}
}
