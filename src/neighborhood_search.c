#include "neighborhood_search.h"
#include "point.h"
#include <math.h>


// Structure to represent a node of the linked list ResidentList of a cell, used to represent the particles contained in a cell
// index : the index in the data table, supposed to be available everywhere it is needed
// next : pointer to the next particle of the linked list
typedef struct node {
	int index;
	struct node* next;
}node;

// Structure to represent a cell and the particles it contains
// nResident : number of particles contained in the linked list ResidentList
// ResidentList : linked list used to contain the index of the particles that are in the current cell
typedef struct cell {
	int nResident;
	node* ResidentList;
}cell;

// funtion to compute the radius kh of the circle of influence of a particle
double compute_KH() {
	double target = 21.0 / NPTS;
	if (NPTS < 21)
		return sqrt(2.0);
	else if (target <= M_PI / 4)
		return sqrt(4 * target / M_PI);
	double tolerance = 0.0000000001;
	double kh_min = 1.0;
	double kh_max = sqrt(2.0);
	while (fabs(kh_max - kh_min) >= tolerance) {
		kh_max = kh_min;
		kh_min = sqrt((target - sin(acos(1.0 / kh_min)) * kh_min) / ((M_PI / 4 - acos(1.0 / kh_min))));
	}
	return kh_min;
}

void checkCellVois(neighborhood* nh, cell* cells, double kh, point* points, int A, int B) {
	node *pointA = cells[A].ResidentList;
	for (int i = 0; i < cells[A].nResident; i++) {
		node* pointB = cells[B].ResidentList;
		neighbours* pointeur;
		if (nh[pointA->index].nNeighbours > 0) {
			pointeur = nh[pointA->index].list;
			for (int k = 0; k < nh[pointA->index].nNeighbours-1; k++) {
				pointeur = pointeur->next;
			}
		}
		for (int j = 0; j < cells[B].nResident; j++) {
			double distX = fabs(points[pointA->index].x - points[pointB->index].x);
			double distY = fabs(points[pointA->index].y - points[pointB->index].y);
			double dist = pow(pow(distX, 2) + pow(distY, 2), 0.5);

			if (dist < kh && !(A == B && i == j)) {
				neighbours *neigh=calloc(1,sizeof(neighbours));
				neigh->index = pointB->index;
				neigh->distance = dist;
				neigh->distX = distX;
				neigh->distY = distY;
				neigh->W;//To determine
				neigh->grad_W;//To determine
				nh[pointA->index].nNeighbours++;
				if (nh[pointA->index].nNeighbours == 1) {
					nh[pointA->index].list = neigh;
				}
				else {
					pointeur->next = neigh;
				}
				pointeur = neigh;
			}
			pointB = pointB->next;
		}

		pointA = pointA->next;
	}
}

neighborhood *createNeighborhood(double kh, point *points){
	neighborhood* nh = calloc(NPTS, sizeof(neighborhood));
	CHECK_MALLOC(nh);

	const int nCells = 200 / kh;
	cell* cells = calloc(nCells*nCells,sizeof(cell));
	node** pointeur = calloc(nCells * nCells,sizeof(point));
	node* noeuds = calloc(NPTS,sizeof(node));
	for (int i = 0; i < NPTS; i++) {
		double x = points[i].x;
		double y = points[i].y;
		nh[i].index = i;
		noeuds[i].index = i;
		int cellNumber = floor((x + 100.0) * nCells / 200.0) + floor((y + 100.0) * nCells / 200.0 ) * nCells;
		cells[cellNumber].nResident++;
		if (cells[cellNumber].nResident == 1) {
			cells[cellNumber].ResidentList = &noeuds[i];
		}
		else {
			pointeur[cellNumber]->next = &noeuds[i];
		}
		pointeur[cellNumber] = &noeuds[i];
	}

	free(pointeur);

	for (int iCell = 0; iCell < nCells * nCells; iCell++) {

		if (iCell==0) {//Cell en b à g
			checkCellVois(nh, cells, kh, points, iCell, iCell + 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells + 1);
		}
		else if (iCell == nCells-1) {//Cell en b à d
			checkCellVois(nh, cells, kh, points, iCell, iCell - 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells - 1);
		}
		else if (iCell == nCells*(nCells-1)) {//Cell en h à g
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells + 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + 1);
		}
		else if (iCell == nCells * nCells - 1) {//Cell en h à d
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells - 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell - 1);
		}
		else if (iCell > 0 && iCell < nCells-1) {//Cell en b
			checkCellVois(nh, cells, kh, points, iCell, iCell - 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells - 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells + 1);
		}
		else if (iCell > nCells* (nCells - 1) && iCell < nCells * nCells - 1) {//Cell en h
			checkCellVois(nh, cells, kh, points, iCell, iCell - 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells - 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells + 1);
		}
		else if (iCell%nCells ==0) {//Cell à g
			checkCellVois(nh, cells, kh, points, iCell, iCell + 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells + 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells + 1);
		}
		else if (iCell%nCells == nCells-1) {//Cell à d
			checkCellVois(nh, cells, kh, points, iCell, iCell - 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells - 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells - 1);
		}
		else {
			checkCellVois(nh, cells, kh, points, iCell, iCell - 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells - 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell - nCells + 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells - 1);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells);
			checkCellVois(nh, cells, kh, points, iCell, iCell + nCells + 1);
		}
		checkCellVois(nh, cells, kh, points, iCell, iCell);
	}
	free(cells);
	free(noeuds);
	return nh;
}

void deleteNeighborhood(neighborhood *nh){
	for (int i = 0; i < NPTS; i++) {
		neighbours* point = nh[i].list;
		for (int j = 0; j < nh[i].nNeighbours/*-1*/; j++) {
			neighbours* temp = point;
			point = point->next;
			free(temp);
		}
	}
	free(nh);
}