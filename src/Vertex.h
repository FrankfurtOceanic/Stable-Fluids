#ifndef VERTEX_H
#define VERTEX_H

#define GLEW_STATIC

#include <GL/glew.h>

/*
typedef struct SPosition
{
    GLfloat x, y;
} Position;

typedef struct SColor
{
    GLfloat r, g, b;
} Color;
*/

typedef struct SVertex
{
    /*
    Position position;
    Color cololr;
    */
    GLfloat x, y, r, g, b, a;
    SVertex(): x(0), y(0), r(0), g(0), b(0), a(0) {}
    SVertex(float x, float y, float r, float g, float b, float a): x(x), y(y), r(r), g(g), b(b), a(a) {}
} Vertex;

#endif
