#include <stdlib.h>
#include <math.h>

#include "computation_data.h"

ComputationData *initialize_computation_data(double Lambda, double kir, Grid *computation_grid,
                                             PhysicsData *data, double tolerances)
{
    ComputationData *new_computation_data = (ComputationData *)malloc(sizeof(ComputationData));
    new_computation_data->Lambda = Lambda;
    new_computation_data->kir = kir;
    new_computation_data->tir = -log(kir / Lambda);
    new_computation_data->tolerances = tolerances;
    new_computation_data->computation_grid = computation_grid;
    new_computation_data->data = data;

    return new_computation_data;
}

void destroy_computation_data(ComputationData *computation_data_to_be_freed)
{
    free(computation_data_to_be_freed);
}
