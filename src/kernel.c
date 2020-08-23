#include "kernel.h"
//#include "neighborhood_search_for_mac.h"
#include <math.h>


// Implementation of the kernel cubic function and return the weight regarding the distance and the radius of the circle
/*
void kernel(GLfloat(*data)[14], GLfloat(*coord)[2], neighborhood* nh, double kh) {
    for (int i = 0; i < NPTS; i++) {
        double val_node_x = data[i][8];
        double val_node_y = data[i][9];
        int nNeigh = nh[i].nNeighbours;
        double val_div = 0;
        double val_grad_x = 0;
        double val_grad_y = 0;
        double val_lapl = 0;
        double dens2 = pow(DENSITY, 2);
        neighbours* List = nh[i].list;
        if (nNeigh > 0) {
            for (int j = 0; j < nNeigh; j++) {
                int index_node2 = List->index;
                double distance = List->distance;
                double d_x = data[index_node2][0] - data[i][0];
                double d_y = data[index_node2][1] - data[i][1];
                
                
                // You can choose here the desired kernel function for your code.
                 
                
                //double weight_x = grad_w_cubic(distance, kh, d_x);
                //double weight_y = grad_w_cubic(distance, kh, d_y);
                
                double weight_x = grad_w_lucy(distance, kh, d_x);
                double weight_y = grad_w_lucy(distance, kh, d_y);
                
                //double weight_x = grad_w_newquartic(distance, kh, d_x);
                //double weight_y = grad_w_newquartic(distance, kh, d_y);
                
                //double weight_x = grad_w_quinticspline(distance, kh, d_x);
                //double weight_y = grad_w_quinticspline(distance, kh, d_y);
                
                val_div += -MASS / DENSITY * ((data[index_node2][8] - val_node_x) * weight_x + (data[index_node2][9] - val_node_y) * weight_y);
                val_grad_x += -DENSITY * MASS * ((val_node_x / dens2) + (data[index_node2][8] / dens2)) * weight_x;
                val_grad_y += -DENSITY * MASS * ((val_node_x / dens2) + (data[index_node2][8] / dens2)) * weight_y;
                val_lapl += 2.0 * MASS / DENSITY * (val_node_x - data[index_node2][8]) * (d_x * weight_x + d_y * weight_y) / pow(distance,2);
                
                List = List->next;
            }
        }
        // All the values of the divergent gradient and laplacien are stored in the data table
        data[i][10] = val_div;
        data[i][11] = val_grad_x;
        data[i][12] = val_grad_y;
        data[i][13] = val_lapl;
    }
    
    //Computation of the error based on the already know function.
    for (int j = 0; j < NPTS; j++) {
        double exact = 3 * pow(data[j][0], 2);
        double error = exact - data[j][10];
    }
}
*/
double grad_w_lucy(double distance, double kh)
{
    double q = distance / kh;
    double alpha = (5 / (M_PI * kh * kh));
    double grad_w = 0;
    if (q >= 0 && q <= 1)
    {
        grad_w = alpha * (-12 * q + 24 * q * q - 12 * q * q * q);
    }
    return  grad_w;
}

double w_lucy(double distance,double kh)
{
    double q = distance / kh;
    double alpha = (5 / (M_PI * kh * kh));
    double W = 0;
    if (distance < kh)
    {
        W = alpha*(1 + 3 * q) * (1 - q) * (1 - q) * (1 - q);
    }
    return W;
}


double grad_w_lucy1D(double dx, double distance, double kh)
{
    double q = distance / kh;
    double alpha = (5 / (M_PI * kh * kh));
    double grad_w = 0;
    if (q >= 0 && q <= 1)
    {
        grad_w = alpha * (-12 * q + 24 * q * q - 12 * q * q * q)*(dx)/(distance*kh);
    }
    return  grad_w;
}














