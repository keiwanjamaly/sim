#ifndef PHYSICS_H
#define PHYSICS_H
#include <stdbool.h>

#include "grid.h"

struct physics_data;

typedef struct physics_data PhysicsData;

// Diffusion Flux (t, k, ux)
double Q(double t, double k, double ux, struct physics_data *data);

// Source (t, k, x)
double S(double t, double k, double x, struct physics_data *data);

// initial condition
double initial_condition(double x, struct physics_data *data);

bool activate_diffusion(struct physics_data *data);

#endif
