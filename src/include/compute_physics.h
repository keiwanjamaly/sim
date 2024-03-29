#ifndef COMPUTE_PHYSICS_H
#define COMPUTE_PHYSICS_H
#include "computation_data.h"

double cal_k(double t, struct computation_data *data);

extern "C" void add_diffusion(double t, double k, double *u, double *u_dot,
                              struct computation_data *input_data);

void compute_source(double t, double k, double *u_dot,
                    struct computation_data *input_data);

#endif // !COMPUTE_PHYSICS_H
