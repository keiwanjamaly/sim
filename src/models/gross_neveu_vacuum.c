#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "helpers.h"
#include "physics.h"

struct physics_data {
  double h;
  double one_over_g2;
  int dimension;
  int dimension_gamma;
  double A_d;
  double N;
  double one_over_N;
};

struct physics_data *initialize_physics_data(double h, double one_over_g2,
                                             int dimension, int dimension_gamma,
                                             double N) {
  struct physics_data *new_physics_data = malloc(sizeof(struct physics_data));
  new_physics_data->h = h;
  new_physics_data->dimension = dimension;
  new_physics_data->dimension_gamma = dimension_gamma;
  new_physics_data->N = N;
  new_physics_data->one_over_N = 1 / N;
  new_physics_data->one_over_g2 = one_over_g2;
  new_physics_data->A_d = calculate_prefactor(dimension);
  return new_physics_data;
}

void free_physics_data(struct physics_data *physics_data_to_be_freed) {
  free(physics_data_to_be_freed);
}

bool activate_diffusion(struct physics_data *data) { return !isinf(data->N); }

double initial_condition(double x, struct physics_data *data) {
  double one_over_g2 = data->one_over_g2;
  double h = data->h;

  return one_over_g2 * x * pow(h, 2);
}

double Q(double t, double k, double ux, struct physics_data *data) {
  int d = data->dimension;
  int d_gamma = data->dimension_gamma;
  double A_d = data->A_d;
  double h = data->h;
  double e = e_b(k, ux);
  double one_over_N = data->one_over_N;

  double return_value = -A_d * one_over_N * pow(k, d + 2) / (2.0 * e);

  return return_value;
}

double S(double t, double k, double x, struct physics_data *data) {
  int d = data->dimension;
  int d_gamma = data->dimension_gamma;
  double A_d = data->A_d;
  double h = data->h;
  double e = e_f(k, x);

  double return_value =
      -x * pow(h, 2) * pow(k, d + 2) * A_d * d_gamma / (2.0 * pow(e, 3));

  return return_value;
}

double left_boundary(Grid *_, double *u, PhysicsData *data) { return 0; }

double right_boundary(Grid *computation_grid, double *u, PhysicsData *data) {
  int N = computation_grid->N;
  return 2 * u[N - 1] - u[N - 2];
}
