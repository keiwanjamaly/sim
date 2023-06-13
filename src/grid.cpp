#include "grid.h"
#include <stdio.h>

extern "C" Grid *create_grid(int N, double *points) {
  Grid *new_grid = (Grid *)malloc(sizeof(Grid));
  new_grid->N = N - 2;
  new_grid->grid_points = (double *)malloc(new_grid->N * sizeof(double));

  double *grid_midpoints = (double *)malloc((new_grid->N + 1) * sizeof(double));
  new_grid->dx = (double *)malloc((new_grid->N + 1) * sizeof(double));
  new_grid->dx_midpoints = (double *)malloc(new_grid->N * sizeof(double));

  copy_internal_grid_points(new_grid, N, points);
  compute_dx(new_grid, N, points);
  compute_grid_midpoints(new_grid, N, points, grid_midpoints);
  compute_dx_midpoints(new_grid, N, grid_midpoints);

  free(grid_midpoints);
  return new_grid;
}

void copy_internal_grid_points(Grid *grid, int N, double *points) {
  grid->grid_left = points[0];
  for (int i = 1; i < (N - 1); i++) {
    grid->grid_points[i - 1] = points[i];
  }
  grid->grid_right = points[N - 1];
}

void compute_dx(Grid *grid, int N, double *points) {
  grid->dx[0] = grid->grid_points[0] - grid->grid_left;
  for (int i = 1; i < grid->N; i++) {
    grid->dx[i] = grid->grid_points[i] - grid->grid_points[i - 1];
  }
  grid->dx[grid->N] = grid->grid_right - grid->grid_points[grid->N - 1];
}

void compute_grid_midpoints(Grid *grid, int N, double *points,
                            double *grid_midpoints) {
  grid_midpoints[0] = (grid->grid_points[0] + grid->grid_left) / 2;
  for (int i = 1; i < grid->N; i++) {
    grid_midpoints[i] = (grid->grid_points[i] + grid->grid_points[i - 1]) / 2;
  }
  grid_midpoints[grid->N] =
      (grid->grid_right + grid->grid_points[grid->N - 1]) / 2;
}

void compute_dx_midpoints(Grid *grid, int N, double *grid_midpoints) {
  for (int i = 0; i < grid->N; i++) {
    grid->dx_midpoints[i] = grid_midpoints[i + 1] - grid_midpoints[i];
  }
}

extern "C" void destroy_grid(Grid *grid_to_be_freed) {
  free(grid_to_be_freed->grid_points);
  free(grid_to_be_freed->dx);
  free(grid_to_be_freed->dx_midpoints);
  free(grid_to_be_freed);
}
