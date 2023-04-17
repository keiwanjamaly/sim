#include <stdlib.h>
#include <math.h>

#include "computation_data.h"

struct computation_data *initialize_computation_data(double Lambda, double kir, struct grid *computation_grid,
                                                     struct physics_data *data, double tolerances)
{
    struct computation_data *new_computation_data = malloc(sizeof(struct computation_data));
    new_computation_data->Lambda = Lambda;
    new_computation_data->kir = kir;
    new_computation_data->tir = -log(kir / Lambda);
    new_computation_data->tolerances = tolerances;
    new_computation_data->computation_grid = computation_grid;
    new_computation_data->data = data;

    return new_computation_data;
}

void free_computation_data(struct computation_data *computation_data_to_be_freed)
{
    free(computation_data_to_be_freed);
}