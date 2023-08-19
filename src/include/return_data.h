#ifndef RETURN_DATA_H
#define RETURN_DATA_H

#include "grid.h"

typedef struct return_data {
  int grid_size;
  int samples;
  double *grid;
  double **solution_y;
  double *solution_time;
  double **diffusion;
  double **source;
} ReturnData;

extern "C" ReturnData *create_return_data(int, Grid *);
extern "C" void destroy_return_data(ReturnData *);

void save_step(ReturnData *return_data_to_be_saved_to, int index, double *y,
               double time, double left_point, double right_point,
               double *source, double *diffusion);
#endif
