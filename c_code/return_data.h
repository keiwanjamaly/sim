#ifndef RETURN_DATA_H
#define RETURN_DATA_H

struct return_data
{
    int grid_size;
    int samples;
    double *grid;
    double **solution_y;
    double *solution_time;
};

void save_step(struct return_data *return_data_to_be_saved_to, int index, double *y, double time, double left_point, double right_point);

#endif