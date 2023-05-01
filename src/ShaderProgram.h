#ifndef SHADER_PROGRAM_H
#define SHADER_PROGRAM_H

#include <string>
#include <glm/glm.hpp>

class ShaderProgram
{
public:
    ShaderProgram(std::string vsFileame, std::string fsFilename);
    void use();
    void set_int(std::string name, int i);
    void set_float(std::string name, float f);
    void set_vec2(std::string name, glm::vec2 v);
    void set_vec3(std::string name, glm::vec3 v);
    void set_vec4(std::string name, glm::vec4 v);
    void set_mat2(std::string name, glm::mat2 m);
    void set_mat3(std::string name, glm::mat3 m);
    void set_mat4(std::string name, glm::mat4 m);

private:
    int program, vShader, fShader;
};

#endif
