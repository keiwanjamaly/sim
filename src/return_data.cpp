#include "return_data.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" ReturnData *create_return_data(int samples, Grid *computation_grid) {
  ReturnData *new_return_data = (ReturnData *)malloc(sizeof(ReturnData));

  // set samples and grid size
  new_return_data->samples = samples;
  new_return_data->grid_size = computation_grid->N + 2;

  // copy over grid points
  new_return_data->grid =
      (double *)malloc(new_return_data->grid_size * sizeof(double));
  new_return_data->grid[0] = computation_grid->grid_left;
  for (int i = 1; i < new_return_data->grid_size - 1; i++) {
    new_return_data->grid[i] = computation_grid->grid_points[i - 1];
  }
  new_return_data->grid[new_return_data->grid_size - 1] =
      computation_grid->grid_right;

  // allocate memory for solution_y
  new_return_data->solution_y = (double **)malloc(samples * sizeof(double *));
  for (int i = 0; i < samples; i++) {
    new_return_data->solution_y[i] =
        (double *)malloc(new_return_data->grid_size * sizeof(double));
  }

  // allocate memory for source
  new_return_data->source = (double **)malloc(samples * sizeof(double *));
  for (int i = 0; i < samples; i++) {
    new_return_data->source[i] =
        (double *)malloc(new_return_data->grid_size * sizeof(double));
  }

  // allocate memory for diffusion
  new_return_data->diffusion = (double **)malloc(samples * sizeof(double *));
  for (int i = 0; i < samples; i++) {
    new_return_data->diffusion[i] =
        (double *)malloc(new_return_data->grid_size * sizeof(double));
  }

  new_return_data->solution_time = (double *)malloc(samples * sizeof(double));

  // return
  return new_return_data;
}

extern "C" void destroy_return_data(ReturnData *return_data) {
  // free solution_y, source, diffusion
  for (int i = 0; i < return_data->samples; i++) {
    free(return_data->solution_y[i]);
    free(return_data->diffusion[i]);
    free(return_data->source[i]);
  }
  free(return_data->solution_y);
  free(return_data->diffusion);
  free(return_data->source);

  // free time array
  free(return_data->solution_time);
  // free grid
  free(return_data->grid);

  // free return_data structure itself
  free(return_data);
}

void save_step(ReturnData *return_data_to_be_saved_to, int index, double *y,
               double time, double left_point, double right_point,
               double *source, double *diffusion) {
  // Save the left point, y values, and right point to the solution_y array at
  // the specified index
  return_data_to_be_saved_to->solution_y[index][0] = left_point;
  memcpy(return_data_to_be_saved_to->solution_y[index] + 1, y,
         (return_data_to_be_saved_to->grid_size - 2) * sizeof(double));
  return_data_to_be_saved_to
      ->solution_y[index][return_data_to_be_saved_to->grid_size - 1] =
      right_point;

  memcpy(return_data_to_be_saved_to->source[index], source,
         (return_data_to_be_saved_to->grid_size) * sizeof(double));

  memcpy(return_data_to_be_saved_to->diffusion[index], diffusion,
         return_data_to_be_saved_to->grid_size * sizeof(double));

  // Save the time value to the corresponding time array
  return_data_to_be_saved_to->solution_time[index] = time;
}
