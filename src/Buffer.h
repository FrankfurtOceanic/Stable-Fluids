#ifndef BUFFER_H
#define BUFFER_H

#include <vector>

#define GLEW_STATIC

#include <GL/glew.h>

class Buffer
{
public:
    Buffer(const std::vector<size_t> &attrSizes, const GLfloat *vertices, size_t vtxSize, const GLuint *indices = nullptr, size_t idxSize = 0);
    void bind();
    void unbind();
    void draw(GLenum mode = GL_TRIANGLES);
    ~Buffer();

private:
    GLuint vao, vbo, ebo;
    size_t numVert, numIdx;
    bool useElement;
};

#endif
