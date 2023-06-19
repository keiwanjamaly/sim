#ifndef COMPUTE_PHYSICS_H
#define COMPUTE_PHYSICS_H
#include "computation_data.h"

double cal_k(double t, struct computation_data *data);

void add_diffusion(double t, double k, double *grid, double *u, double *u_dot,
                   int N, struct computation_data *input_data);

void compute_source(double t, double k, double *grid, double *u_dot, int N,
                    struct computation_data *input_data);

#endif // !COMPUTE_PHYSICS_H
