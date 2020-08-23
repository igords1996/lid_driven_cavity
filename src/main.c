#include "point.h"
#include "kernel.h"
#include "source.h"
#include <math.h>
#include <stdio.h>
#include <time.h>
//
#include <omp.h>
int NPTS = 1000;// +25 * 5 * 4 + 20;
int NPTS_boundaries = 0; // 5 * 25 + 5;// 24;
double vmax = 1;
double MASSE = 1;
int NUMBER_ITERATIONS = 50000;
double dt = 0.001;//0.0005
int main(){
	// Creation des points
	clock_t begin = clock();
	printf("max = %d \n", omp_get_max_threads());
	omp_set_num_threads(omp_get_max_threads());
	bov_window_t* window = bov_window_new(1024, 780, "ANM Project: SPH");
	bov_window_set_color(window, (GLfloat[]) { 0.9f, 0.85f, 0.8f, 0.0f });
	point* points = malloc(sizeof(point) * NPTS);
	fillPointsRand(points);
	GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
	CHECK_MALLOC(data);

	updateData(points, data);
	/* send data to GPU, and receive reference to those data in a points object */
	bov_points_t* particles = bov_particles_new(data, NPTS, GL_STATIC_DRAW);

	/* setting particles appearance */
	bov_points_set_width(particles, 0.01);
	bov_points_set_outline_width(particles, 0.00125);

	/* fit particles in a [-0.8 0.8]x[-0.8 0.8] square */
	bov_points_scale(particles, (GLfloat[2]) { 0.008, 0.008 });
	bov_points_set_pos(particles, (GLfloat[2]) { 0.0, -0.1 });
	bov_text_t* msg = bov_text_new(
		(GLubyte[]) {
		"Rendering " xstr(NPTS) " particles"
	},
		GL_STATIC_DRAW);
	bov_text_set_pos(msg, (GLfloat[2]) { -0.95, 0.82 });
	bov_text_set_fontsize(msg, 0.1);
	// Seed the random
	time_t seed = time(NULL);
	//printf(" %u \n", seed);
	srand(seed);
	//neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
	//neighborhood* nh = options->nh;
	double kh = compute_kh(NPTS)/2; 
	//	point* points = malloc(sizeof(point) * NPTS);
	float* vx = malloc(sizeof(vx)*NPTS); 
	float* vy = malloc(sizeof(vy) * NPTS);
	int direction = 1;
	int iteration = 1;
	while (!bov_window_should_close(window)&& iteration<NUMBER_ITERATIONS) 
	{
		if (GetAsyncKeyState(VK_RETURN) & 0x8000) {
			fillPointsRand(points);
			updateData(points, data);
			iteration = 1;
		}
		node_t** voisinage = malloc(sizeof(node_t) * NPTS); // Initialisation neighbour
		update_neighboor(points, voisinage, kh);
		setDensity(points, voisinage, kh);
		//updateDensity(points, voisinage, kh); 
		updatePressure(points,voisinage, kh);
		updateVelocity(points, voisinage, kh, vx, vy);
		if(GetAsyncKeyState(VK_DOWN) & 0x8000) 
		{
			direction = 1;
		}
		if(GetAsyncKeyState(VK_RIGHT) & 0x8000)
		{
			direction = 2;
		}
		if(GetAsyncKeyState(VK_LEFT) & 0x8000)
		{
			direction = 3;
		}
		if(GetAsyncKeyState(VK_UP) & 0x8000)
		{
			direction = 4;
		}
		if (GetAsyncKeyState(VK_DELETE) & 0x8000)
		{
			direction = 5;
		}
		updatePosition(points, direction);
		updateData(points, data);
		free(voisinage);
		bov_text_draw(window, msg);
		bov_particles_draw(window, particles, 0, BOV_TILL_END);
		bov_particles_update(particles, data, NPTS);
		char s[60];
		sprintf(s, "Number of iterations: %d", iteration);
		bov_text_update(msg, s);
		iteration += 1;
		bov_window_update(window); //don't wait for events => bov_window_update(window)
	}
	//deleteNeighborhood(nh);

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("timer start = %f \n", time_spent);
	free(data);
	free(points);
	free(vx);
	free(vy);
	return EXIT_SUCCESS;
}
