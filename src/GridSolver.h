#ifndef GRID_SOLVER_H
#define GRID_SOLVER_H

#define BIG_GRID

#define GLEW_STATIC

#include <vector>
#include <deque>
#include <memory>
#include <GL/glew.h>

#include "Vertex.h"

#ifndef BIG_GRID
const int cols = 5, rows = 6;
const int cellSize = 40; // cell size in display
#else
const int cols = 50, rows = 40;
const int cellSize = 20; // cell size in display
#endif
const int frameWidth = cols * cellSize, frameHeight = rows * cellSize;

extern float nv, kappa, buoyancy, dt, velScale;
extern float simTime;
extern float dx;
extern float refinementThreshold;
extern int iterations;

extern double tpSrc[rows + 2][cols + 2];
extern double vxSrc[rows + 2][cols + 2];
extern double vySrc[rows + 2][cols + 2];

extern double (*cur_tp)[cols + 2];

extern std::deque<std::shared_ptr<std::vector<Vertex>>> filaments;
extern std::deque<double> ages;
extern float maxAge;

extern Vertex gridVertices[rows + 1][cols + 1];
extern GLuint gridIndices[6 * rows * cols];
extern Vertex velVertices[rows][cols][2];

void initGrid();
double getReferenceTemperature(double tp[rows + 2][cols + 2]);
void step();
void updateGrid();
void addFilament(double x, double y);
void updateFilament();

#endif
