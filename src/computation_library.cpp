#include <alloca.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunmatrix/sunmatrix_band.h>
#include <time.h>

#include "computation_data.h"
#include "compute_physics.h"
#include "grid.h"
#ifdef ACTIVATE_LIVE_PLOTTING
#include "live_plotting.h"
#endif
#include "physics.h"
#include "return_data.h"

int f_without_diffusion(double t, N_Vector u, N_Vector udot, void *input_data) {
  struct computation_data *user_data = (struct computation_data *)input_data;
  int N = user_data->computation_grid->N;
  double k = cal_k(t, user_data);
  double *grid = user_data->computation_grid->grid_points;
  double *udot_data = N_VGetArrayPointer(udot);

  compute_source(t, k, udot_data, user_data);
  return 0;
}

int f_with_diffusion(double t, N_Vector u, N_Vector udot, void *input_data) {
  struct computation_data *user_data = (struct computation_data *)input_data;
  int N = user_data->computation_grid->N;
  double k = cal_k(t, user_data);
  double *grid = user_data->computation_grid->grid_points;
  double *udot_data = N_VGetArrayPointer(udot);
  double *u_data = N_VGetArrayPointer(u);

  compute_source(t, k, udot_data, user_data);
  add_diffusion(t, k, u_data, udot_data, user_data);
  return 0;
}

double compute_time_from_start(clock_t start_clock) {
  return (double)(clock() - start_clock) / CLOCKS_PER_SEC;
}

void save_computation_step(struct computation_data *data,
                           double *u_output_pointer, double t, int i,
                           struct return_data *return_struct) {
  double left_point = left_boundary(data->computation_grid, u_output_pointer);
  double right_point = right_boundary(data->computation_grid, u_output_pointer);
  int N = return_struct->grid_size;
  double *source = (double *)malloc(N * sizeof(double));
  double *diffusion = (double *)calloc(N, sizeof(double));
  double k = cal_k(t, data);
  source[0] = S(t, k, return_struct->grid[0], data->data);
  compute_source(t, k, source + 1, data);
  add_diffusion(t, k, u_output_pointer, diffusion + 1, data);
  diffusion[N - 2] = diffusion[N - 3];
  diffusion[N - 1] = diffusion[N - 2];
  source[N - 1] = S(t, k, return_struct->grid[N - 1], data->data);

  save_step(return_struct, i, u_output_pointer, t, left_point, right_point,
            source, diffusion);

  free(source);
  free(diffusion);
}

extern "C" int compute(struct computation_data *data,
                       struct return_data *return_struct) {
  int steps_to_save = return_struct->samples;
  double t_final = data->tir;
  double t_out;
  double t_now = 0.0;
  double dt = t_final / (double)(steps_to_save - 1);
  double reltol = data->tolerances;
  double abstol = data->tolerances;
  SUNContext sunctx;
  void *package_mem;
  int status;
  N_Vector u0, u_out;
  SUNMatrix jacobi_matrix;
  SUNLinearSolver lin_sol;
  SUNContext_Create(NULL, &sunctx);

  // setup plotting library
#ifdef ACTIVATE_LIVE_PLOTTING
  LivePlottingData *plotting_data = setup_live_plotting(data);
#endif // DEBUG

  u_out = N_VNew_Serial(data->computation_grid->N, sunctx);
  u0 = N_VNew_Serial(data->computation_grid->N, sunctx);
  double *u0_c_pointer = N_VGetArrayPointer(u0);
  double *u_output_pointer = N_VGetArrayPointer(u_out);
  for (int i = 0; i < data->computation_grid->N; i++) {
    u0_c_pointer[i] =
        initial_condition(data->computation_grid->grid_points[i], data->data);
    u_output_pointer[i] = u0_c_pointer[i];
  }

  package_mem = CVodeCreate(CV_BDF, sunctx);

  if (activate_diffusion(data->data)) {
    CVodeInit(package_mem, f_with_diffusion, t_now, u0);
  } else {
    CVodeInit(package_mem, f_without_diffusion, t_now, u0);
  }

  // set user data
  CVodeSetUserData(package_mem, data);

  CVodeSStolerances(package_mem, reltol, abstol);

#ifdef ACTIVATE_LIVE_PLOTTING
  CVodeSetMaxNumSteps(package_mem, 10);
#endif // ACTIVATE_LIVE_PLOTTING
#ifndef ACTIVATE_LIVE_PLOTTING
  CVodeSetMaxNumSteps(package_mem, 1000);
#endif // !ACTIVATE_LIVE_PLOTTING

  jacobi_matrix = SUNBandMatrix(data->computation_grid->N, 1, 1, sunctx);

  lin_sol = SUNLinSol_Band(u0, jacobi_matrix, sunctx);

  CVodeSetLinearSolver(package_mem, lin_sol, jacobi_matrix);

  CVodeSetMaxErrTestFails(package_mem, 1000);

  CVodeSetErrFile(package_mem, NULL);

  clock_t start_clock = clock();

  // save zeroth step
  int i = 0;
  t_now = 0.0;
  // save_step(return_struct, 0, u_output_pointer, t_now, left_point,
  // right_point);
  save_computation_step(data, u_output_pointer, t_now, i, return_struct);
  for (i = 1; i < steps_to_save; i++) {
    t_out = dt * i;
    status = CV_TOO_MUCH_WORK;
    while (t_now < t_out && status == CV_TOO_MUCH_WORK) {
      status = CVode(package_mem, t_out, u_out, &t_now, CV_NORMAL);
#ifdef ACTIVATE_LIVE_PLOTTING
      draw_frame(plotting_data, t_now, u_output_pointer,
                 compute_time_from_start(start_clock));
#endif // DEBUG
    }

    if (status != CV_SUCCESS) {
      printf("Error: something went wrong! CVODE error code %d\n", status);
      return status;
    }
    // save after each step
    save_computation_step(data, u_output_pointer, t_now, i, return_struct);
  }

  // destory plotting_library
#ifdef ACTIVATE_LIVE_PLOTTING
  tear_down_live_plotting(plotting_data);
#endif

  // Free stuff
  N_VDestroy(u0);
  N_VDestroy(u_out);
  SUNMatDestroy(jacobi_matrix);
  SUNLinSolFree(lin_sol);
  CVodeFree(&package_mem);
  SUNContext_Free(&sunctx);
  return 0;
}
