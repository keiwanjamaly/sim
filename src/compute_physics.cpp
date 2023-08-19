#include "compute_physics.h"
#include <alloca.h>
#include <math.h>

double cal_k(double t, struct computation_data *data) {
  return data->Lambda * exp(-t);
}

void add_diffusion(double t, double k, double *u, double *u_dot,
                   struct computation_data *input_data) {
  double u_boundary;
  int N = input_data->computation_grid->N;
  double *dx = input_data->computation_grid->dx;
  double *dx_midpoints = input_data->computation_grid->dx_midpoints;
  double *u_x = (double *)alloca((N + 1) * sizeof(double));
  double *Q_cal = (double *)alloca((N + 1) * sizeof(double));
  // double *u_x = (double *)malloc((N + 1) * sizeof(double));
  // double *Q_cal = (double *)malloc((N + 1) * sizeof(double));

  // compute u_x
  // handle left boundary
  u_boundary = left_boundary(input_data->computation_grid, u);
  u_x[0] = (u[0] - u_boundary) / dx[0];
  for (int i = 1; i < N; i++) {
    u_x[i] = (u[i] - u[i - 1]) / dx[i];
  }
  // handle right boundary
  u_boundary = right_boundary(input_data->computation_grid, u);
  u_x[N] = (u_boundary - u[N - 1]) / dx[N];

  for (int i = 0; i < N + 1; i++) {
    Q_cal[i] = Q(t, k, u_x[i], input_data->data);
  }

  for (int i = 0; i < N; i++) {
    u_dot[i] += (Q_cal[i + 1] - Q_cal[i]) / dx_midpoints[i];
  }

  // free(u_x);
  // free(Q_cal);
}

void compute_source(double t, double k, double *u_dot,
                    struct computation_data *input_data) {
  int N = input_data->computation_grid->N;
  double *grid = input_data->computation_grid->grid_points;
  for (int i = 0; i < N; i++) {
    u_dot[i] = S(t, k, grid[i], input_data->data);
  }
}
