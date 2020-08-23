#ifndef SOURCE_H
#define SOURCE_H

#include "BOV.h"
#include "./neighborhood_search.h"
#include "./point.h"
#include "kernel.h"
#include <time.h>
#include <math.h>

extern int NPTS;
extern double MASSE;
extern double SOURCE_TEMP;
extern double INIT_TEMP;
extern double dt;
extern double vmax;

// see stringification process
#define xstr(s) str(s)
#define str(s) #s

// Structure to represent a point

void updateDensity(point* points, node_t** voisinage, double kh);

void setDensity(point* points, node_t** voisinage, double kh);

void updatePressure(point* points, node_t** voisinage, double kh);

void updateVelocity(point* points, node_t** voisinage, double kh, float* vx, float* vy);

void updatePosition(point* points, int direction);

#endif