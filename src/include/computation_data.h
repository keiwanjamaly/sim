#ifndef COMPUTATION_DATA_H
#define COMPUTATION_DATA_H

#include "physics.h"

typedef struct computation_data
{
    double Lambda;
    double kir;
    double tir;
    double tolerances;
    struct grid *computation_grid;
    struct physics_data *data;
} ComputationData;

extern "C" ComputationData *initialize_computation_data(double Lambda, double kir, Grid *computation_grid,
                                             PhysicsData *data, double tolerances);
extern "C" void destroy_computation_data(ComputationData *computation_data_to_be_freed);

#endif
