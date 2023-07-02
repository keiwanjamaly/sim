#ifndef GRID_H
#define GRID_H

#include <stdlib.h>

typedef struct grid {
  int N;
  double c_1, c_2;
  double grid_left;
  double grid_right;
  double *grid_points;
  double *dx;
  double *dx_midpoints;
} Grid;

extern "C" Grid *create_grid(int N, double *points);
double right_boundary(Grid *computation_grid, double *u);
double left_boundary(Grid *_, double *u);
extern "C" void destroy_grid(Grid *grid_to_be_freed);

#endif // GRID_H
