#ifndef GRID_H
#define GRID_H

#include <stdlib.h>

typedef struct grid {
  int N;
  double grid_left;
  double grid_right;
  double *grid_points;
  double *dx;
  double *dx_midpoints;
} Grid;

extern "C" Grid *create_grid(int N, double *points);
void copy_internal_grid_points(Grid *grid, int N, double *points);
void compute_dx(Grid *grid, int N, double *points);
void compute_grid_midpoints(Grid *grid, int N, double *points,
                            double *grid_midpoints);
void compute_dx_midpoints(Grid *grid, int N, double *grid_midpoints);
extern "C" void destroy_grid(Grid *grid_to_be_freed);

#endif // GRID_H
