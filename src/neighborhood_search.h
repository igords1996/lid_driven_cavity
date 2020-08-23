#ifndef NEIGHBORHOOD_SEARCH_H
#define NEIGHBORHOOD_SEARCH_H

#include "BOV.h"
#include "point.h"
#include <time.h>
#include <math.h>

extern int NPTS;

// malloc verification
// To use after each call to malloc, calloc and realloc
#define CHECK_MALLOC(ptr) if((ptr)==NULL) { \
		BOV_ERROR_LOG(BOV_OUT_OF_MEM_ERROR, "Memory allocation failed"); \
		exit(EXIT_FAILURE); }


// Structure to represent a neighbours, which is a node of the neighborhoods linked lists
// index : the index in the data table, supposed to be available everywhere it is needed
// distance : distance between the particle at the index in data and the particle that owns the neighborhood
// next : pointer to the next neighbours of the linked list
typedef struct neighbours {
	int index;
	double distance;
	double distX;
	double distY;
	double W;
	double grad_W;
	struct neighbours* next;
}neighbours;

// Structure to represent the neighborhood of a particle
// index : the index in the data table, supposed to be available everywhere it is needed
// nNeighbours : number of neighbours in the linked list list
// list : linked list of the actual neighbours of the particle at the index index in the data table
typedef struct neighborhood {
	int index;
	int nNeighbours;
	neighbours* list;
}neighborhood;

double compute_KH();

neighborhood* createNeighborhood(double kh, point* points);

void deleteNeighborhood(neighborhood* nh);

#endif
