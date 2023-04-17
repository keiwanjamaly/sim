#ifndef COMPUTATION_DATA_H
#define COMPUTATION_DATA_H

#include "physics.h"

struct computation_data
{
    double Lambda;
    double kir;
    double tir;
    double tolerances;
    struct grid *computation_grid;
    struct physics_data *data;
};

#endif