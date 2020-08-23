#ifndef POINT_H
#define POINT_H

#include "BOV.h"
#include <time.h>
#include <math.h>

extern int NPTS;
extern int NPTS_boundaries;
extern double vmax;
extern double MASSE;
extern double SOURCE_TEMP;
extern double INIT_TEMP;

// see stringification process
#define xstr(s) str(s)
#define str(s) #s
#define M_PI 3.14159265358979323846

// Structure to represent a point
typedef struct point {
	float x;
	float y;
	float vx;
	float vy;
	float old_density;
	float density;
	float old_dpdx;
	float old_dpdy;
	float dpdx;
	float dpdy;
	float pressure;
	int type; //1: normal point, 2: boundary point
	/*
	double val;
	double div;
	double grad_x;
	double grad_y;
	double lapl;
	*/
}point;

void fillPointsRand(point* points);

void updateData(point* points, GLfloat(*data)[8]);

static void colormap(float v, float color[3]);

typedef struct voisin {
	int id;
	double distance;
	double dx;
	double dy;
}voisin;

typedef struct node
{
	voisin actu;
	struct node * next;
}node_t;

void update_neighboor(point* points, node_t** voisinage, double kh);
double compute_kh(double NPTS);
#endif