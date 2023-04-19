#ifndef RETURN_DATA_H
#define RETURN_DATA_H

#include "grid.h"

typedef struct return_data
{
    int grid_size;
    int samples;
    double *grid;
    double **solution_y;
    double *solution_time;
} ReturnData;

ReturnData *create_return_data(int, Grid *);
void destroy_return_data(ReturnData *);

void save_step(ReturnData *return_data_to_be_saved_to, int index, double *y, double time, double left_point, double right_point);
#endif