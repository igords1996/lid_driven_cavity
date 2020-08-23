#include "./source.h"
#include <math.h>
void setDensity(point* points, node_t** voisinage, double kh)
{
	int i;
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for ( i = 0; i < NPTS; i++)
		{
			double density = 0;
			// double density = MASS * w_lucy(0, kh);
			//neighbours* point = nh[i].list;
			//for (int j = 0; j < (int)nh[i].nNeighbours; j++)
			node_t* current = voisinage[i];
			while (current != NULL)
			{
				density += MASSE * w_lucy(current->actu.distance, kh);
				current = current->next;
			}
			points[i].old_density = density;
		}
	}
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for (i = 0; i < NPTS; i++)
		{
			double num = 1;
			double denom = 1 / points[i].old_density;
			node_t* current = voisinage[i];
			while (current != NULL)
			{
				num += w_lucy(current->actu.distance, kh);
				denom += w_lucy(current->actu.distance, kh) / points[current->actu.id].old_density;
				current = current->next;
			}
			/*for (int j = 0; j < NPTS; j++)
			{
				double dx = ((double)points[i].x - (double)points[j].x) / 100.0;
				double dy = ((double)points[i].y - (double)points[j].y) / 100.0;
				double distance = sqrt(dx * dx + dy * dy);
				if (distance < kh)
				{
					num += w_lucy(distance, kh);
					denom += w_lucy(distance, kh) / points[j].old_density;
				}

			}*/
			points[i].density = num / denom;
		}
	}
}

// DENSITY SET OK
void updateDensity(point* points, node_t** voisinage, double kh)
{
	int i;
	#pragma omp parallel
	{	
		#pragma omp for schedule(dynamic)

		for (i = 0; i < NPTS - 4 * NPTS_boundaries; i++)
		{
			double der_density = 0;


			double vxi = (double)points[i].vx;
			double vyi = (double)points[i].vy;
			double speed_i = sqrt(vxi * vxi + vyi * vyi);
			node_t* current = voisinage[i]->next;
			while (current != NULL)
			{
				int j = current->actu.id;
				double vxj = (double)points[j].vx;
				double vyj = (double)points[j].vy;
				double speed_j = sqrt(vxj * vxj + vyj * vyj);
				double der_densityx = MASSE * (vxi - vxj) * grad_w_lucy1D(current->actu.dx, current->actu.distance, kh);
				double der_densityy = MASSE * (vyi - vyj) * grad_w_lucy1D(current->actu.dy, current->actu.distance, kh);
				der_density += der_densityx + der_densityy;
				//der_density += MASSE * (speed_i - speed_j) * grad_w_lucy(current->actu.distance, kh);
				current = current->next;
			}
			points[i].density = points[i].density + dt * der_density;
		}
	}
}

void updatePressure(point* points, node_t** voisinage, double kh)
{
	double rho_0= 1000; // kg/m^3 => minimum value with 400 particles
	//rho_0 = 100;
	double B = 10130;//10130; // Pa
	double gamma = 7;
	int i;
	#pragma omp parallel
	{
	#pragma omp for schedule(dynamic)
		for ( i = 0; i < NPTS; i++)
		{
			double dens = points[i].density;
			int mindensity = 1000;
			if (points[i].density < mindensity)
			{
				dens = mindensity;
			}


			points[i].pressure = B * (pow(dens / rho_0, gamma) - 1);
		}
	}
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for ( i = 0; i < NPTS - 4 * NPTS_boundaries; i++)
		{
			points[i].old_dpdx = 0;
			points[i].old_dpdy = 0; // Need to initialise the variable first
			node_t* current = voisinage[i]->next; // on veut pas le i=j
			while (current != NULL)
			{
				int j = current->actu.id;
				points[i].old_dpdx -= MASSE * (points[i].pressure / (points[i].density * points[i].density) + points[j].pressure / (points[j].density * points[j].density)) * grad_w_lucy1D(current->actu.dx, current->actu.distance, kh);
				points[i].old_dpdy -= MASSE * (points[i].pressure / (points[i].density * points[i].density) + points[j].pressure / (points[j].density * points[j].density)) * grad_w_lucy1D(current->actu.dy, current->actu.distance, kh);
				current = current->next;
			}

			/*for (int j = 0; j < NPTS ; j++)
			{
				double dx = ((double)points[i].x - (double)points[j].x) / 100.0;
				double dy = ((double)points[i].y - (double)points[j].y) / 100.0;
				double distance = sqrt(dx * dx + dy * dy);
				if (distance < kh && i != j)
				{
					points[i].old_dpdx -= MASSE * (points[i].pressure / (points[i].density * points[i].density) + points[j].pressure / (points[j].density * points[j].density)) * grad_w_lucy1D(dx, distance, kh);
					points[i].old_dpdy -= MASSE * (points[i].pressure / (points[i].density * points[i].density) + points[j].pressure / (points[j].density * points[j].density)) * grad_w_lucy1D(dy, distance, kh);
				}
			}*/
		}
	}
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for (i = 0; i < NPTS - 4 * NPTS_boundaries; i++)
		{
			double numx = points[i].old_dpdx;
			double numy = points[i].old_dpdy;
			double denomx = 0;
			double denomy = 0;
			node_t* current = voisinage[i]->next;
			while (current != NULL)
			{
				if (points[current->actu.id].type == 1)
				{
					denomx -= grad_w_lucy1D(current->actu.dx, current->actu.distance, kh) * current->actu.dx * MASSE / points[current->actu.id].density;
					denomy -= grad_w_lucy1D(current->actu.dy, current->actu.distance, kh) * current->actu.dy * MASSE / points[current->actu.id].density;
				}
				current = current->next;
			}
			/*for (int j = 0; j < NPTS - 4 * NPTS_boundaries; j++)
			{
				double dx = ((double)points[i].x - (double)points[j].x) / 100.0;
				double dy = ((double)points[i].y - (double)points[j].y) / 100.0;
				double distance = sqrt(dx * dx + dy * dy);
				if (distance < kh && i!=j )
				{
					denomx -= grad_w_lucy1D(dx , distance, kh) * dx * MASSE / points[j].density;
					denomy -= grad_w_lucy1D(dy, distance, kh) * dy * MASSE / points[j].density;
				}
			}*/
			points[i].dpdx = numx / denomx;
			points[i].dpdy = numy / denomy;
		}
	}
}

void updateVelocity(point* points, node_t** voisinage, double kh, float* vx, float* vy)
{
	double kin_v = 0.1;// kinematic viscosity m^2/s
	double g = -9.81;// 0.005; // gravity m/s^2
	int i;
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for ( i = 0; i < NPTS - 4 * NPTS_boundaries; i++)
		{
			double cut_off = kh;//0.05;kh =0.25
			// PRESSURE PART
			float dvx = points[i].dpdx;//
			float dvy = points[i].dpdy;///g
			node_t* current = voisinage[i]->next; // on évite le i =j 
			while (current != NULL)
			{
				int j = current->actu.id;
				dvx += 2 * kin_v * MASSE / (points[j].density) * ((double)points[i].vx - (double)points[j].vx) * grad_w_lucy1D(current->actu.dy, current->actu.distance, kh) / (current->actu.dy + 0.000001);// *dx / (dx * dx);
				dvy += 2 * kin_v * MASSE / (points[j].density) * ((double)points[i].vy - (double)points[j].vy) * grad_w_lucy1D(current->actu.dx, current->actu.distance, kh) / (current->actu.dx + 0.000001); //*dy / (dy * dy);
				current = current->next;
			}
			/*
			for (int j = 0; j < NPTS-4*NPTS_boundaries; j++)
			{
				double dx = ((double)points[i].x - (double)points[j].x) / 100.0;
				double dy = ((double)points[i].y - (double)points[j].y) / 100.0;
				double distance = sqrt(dx * dx + dy * dy);
				if (distance < kh && i!=j)
				{
					// PRESSURE PART
					//dvx -= 1 / points[i].density * MASSE * ((double)points[i].pressure / ((double)points[i].density * (double)points[i].density) + (double)points[j].pressure / ((double)points[j].density * (double)points[j].density)) * grad_w_lucy(distance, kh) * dx / distance;
					//dvy -= 1 / points[i].density * MASSE * ((double)points[i].pressure / ((double)points[i].density * (double)points[i].density) + (double)points[j].pressure / ((double)points[j].density * (double)points[j].density)) * grad_w_lucy(distance, kh) * dy / distance;
					// VELOCITY PART
					dvx += 2 * kin_v * MASSE / (points[j].density) * ((double)points[i].vx - (double)points[j].vx) * grad_w_lucy1D(dy,distance, kh) /(dy+0.0001);// *dx / (dx * dx);
					dvy += 2 * kin_v * MASSE / (points[j].density) * ((double)points[i].vy - (double)points[j].vy) * grad_w_lucy1D(dx,distance, kh) /(dx+0.0001); //*dy / (dy * dy);
					// type Lennard-Jones
					if (distance < cut_off)
					{
						//double damper = 0.00000005;
						//dvx += damper*0.01 * dx / distance * (pow(cut_off / distance, 12) - pow(cut_off / distance, 4)) / MASSE;
						//dvy += damper*0.01 * dy / distance * (pow(cut_off / distance, 12) - pow(cut_off / distance, 4)) / MASSE;
					}
					// type coulomb
					//double k = 100;
					//dvx += k / points[i].density * (1/(double)points[i].density + 1/(double)points[i].density) / distance * dx / distance;
					//dvy += k / points[i].density * (1/(double)points[i].pressure + 1/(double)points[i].pressure) / distance * dy / distance;
					// loi parabolique
					//dvx += k*points[j].density*1/distance*dx / distance;
					//dvy +=k* points[j].density * 1 / distance * dy / distance;
					// d rho => f =>
					//double k = 0.000000001;
					//dvx += k * points[j].density * dx / (distance+0.001) * ((double)points[i].density - (double)points[j].density);
					//dvy += k * points[j].density * dy / (distance+0.001) * ((double)points[i].density - (double)points[j].density);
				}
			}*/
			vx[i] = points[i].vx + dvx * dt; // euler method
			vy[i] = points[i].vy + dvy * dt;
		}
	}
	#pragma omp parallel
	{
	#pragma omp for schedule(dynamic)
		for ( i = 0; i < NPTS - 4 * NPTS_boundaries; i++)
		{
			points[i].vx = vx[i];
			points[i].vy = vy[i];
		}
	}

}

void updatePosition(point* points, int direction)
{
	int i;
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for (i = 0; i < NPTS - 4 * NPTS_boundaries; i++)
		{
			//points[i].x += points[i].vx * dt; //euler
			//points[i].y += points[i].vy * dt;
			double cut_off = 0.05;
			if (points[i].x > 1 && direction == 2)
			{
				points[i].vy = 5;
			}
			if (points[i].x < 0 && direction == 3)
			{
				points[i].vy = 5;
			}
			if (points[i].y > 1 && direction == 4)
			{
				points[i].vx = 5;
			}
			if (points[i].y < 0 && direction == 1)
			{
				points[i].vx = 5;
			}
			if (points[i].x > 1)
			{
				points[i].vx -= dt*(pow(cut_off / (cut_off - ((double)points[i].x - 1)), 12) - pow(cut_off / (cut_off - ((double)points[i].x - 1)), 4)) / MASSE;
				/*points[i].x = 100 - (points[i].x - 100);
				points[i].vx = - 0.8*points[i].vx;
				//points[i].vy = - 1;
				if (points[i].x < 0)
				{
					points[i].x = i/2;
					points[i].vx = 0;
					points[i].vy = 0;
				}*/
			}
			else if (points[i].x < 0)
			{
				points[i].vx += dt*(pow(cut_off / (cut_off + points[i].x), 12) - pow(cut_off / (cut_off + points[i].x), 4)) / MASSE;
				/*points[i].x = -5-(points[i].x+5);
				points[i].vx = -0.8*points[i].vx;
				points[i].vy = 1;
				if (points[i].x > 100)
				{
					points[i].x = i/2;
					points[i].vx = 0;
					points[i].vy = 0;
				}*/
			}
			if (points[i].y > 1)
			{
				points[i].vy -= dt*(pow(cut_off / (cut_off - ((double)points[i].y - 1)), 12) - pow(cut_off / (cut_off - ((double)points[i].y - 1)), 4)) / MASSE;
				/*
				points[i].y = 100 - (points[i].y - 100);
				points[i].vy = -0.8*points[i].vy;
				points[i].vx =  1;
				if (points[i].y < -0)
				{
					points[i].y = i/2;
					points[i].vx = 0;
					points[i].vy = 0;
				}*/

			}
			else if (points[i].y < 0)
			{
				points[i].vy += dt*(pow(cut_off / (cut_off + points[i].y), 12) - pow(cut_off / (cut_off + points[i].y), 4)) / MASSE;
				/*
				points[i].y = -5 - (points[i].y + 5);
				points[i].vy = -0.8*points[i].vy;
				points[i].vx = - 1;
				if (points[i].y > 100)
				{
					points[i].y = i/2;
					points[i].vx = 0;
					points[i].vy = 0;
				}*/
			}
			
			double speed = sqrt(points[i].vx * points[i].vx + points[i].vy * points[i].vy);
			if (speed > 8) // Help limit error propagations
			{
				points[i].vx = points[i].vx / speed * 8;
				points[i].vy = points[i].vy / speed * 8;
			}
			points[i].x += points[i].vx * dt; //euler
			points[i].y += points[i].vy * dt;
		}
	}
}
