#include "u_plotting.h"
#include "compute_physics.h"
#include "implot.h"
#include "physics.h"
#include <algorithm>
#include <stdlib.h>

void plot_u(LivePlottingData *plotting_data, double *u) {
  Grid *computation_grid = plotting_data->data->computation_grid;
  if (ImPlot::BeginPlot("Function Plot")) {
    int N = computation_grid->N;
    double *x_data = computation_grid->grid_points;
    double *y_data = u;

    ImPlot::PlotLine("u", x_data, y_data, N);
    ImPlot::EndPlot();
  }
}

void plot_u_within_zero_to_one(LivePlottingData *plotting_data, double *u) {
  Grid *computation_grid = plotting_data->data->computation_grid;
  int N = 0;
  double max_u = 0.0;
  double min_u = 0.0;
  while (computation_grid->grid_points[N] <= 1.1)
  {
      N++;
      max_u = std::max(max_u, u[N]);
      min_u = std::min(min_u, u[N]);
  }
  max_u *= 1.1;
  if (min_u >= 0) {
      min_u *= 0.9;
  } else {
      min_u *= 1.3;
  }
  if (ImPlot::BeginPlot("Function Plot")) {
    ImPlot::SetupAxis(ImAxis_X1, "sigma");
    ImPlot::SetupAxis(ImAxis_Y1, "");
    ImPlot::SetupAxisLimits(ImAxis_Y1, min_u, max_u,
                            ImPlotCond_Always);
    double *x_data = computation_grid->grid_points;
    double *y_data = u;

    ImPlot::PlotLine("u", x_data, y_data, N);
    ImPlot::EndPlot();
  }
}

void plot_Q(double *u, double t, LivePlottingData *plotting_data) {
  Grid *computation_grid = plotting_data->data->computation_grid;

  double k = cal_k(t, plotting_data->data);
  int N = computation_grid->N;
  double *Q = (double *)malloc(N * sizeof(double));
  for (int i = 0; i < N; i++) {
    Q[i] = 0.0;
  }
  add_diffusion(t, k, computation_grid->grid_points, u, Q, N,
                plotting_data->data);

  // calculate the maximum
  double maximum_of_Q = 0.0;
  for (int i = 0; i < N; i++) {
    plotting_data->max_Q[i] = std::max(Q[i], plotting_data->max_Q[i]);
    maximum_of_Q = std::max(maximum_of_Q, plotting_data->max_Q[i]);
  }

  // ImPlot::SetNextAxisLimits(1, 0, maximum_of_Q);
  if (ImPlot::BeginPlot("Function Plot")) {
    ImPlot::SetupAxis(ImAxis_X1, "sigma");
    ImPlot::SetupAxis(ImAxis_Y1, "");
    ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, maximum_of_Q * 1.1,
                            ImPlotCond_Always);
    double *x_data = computation_grid->grid_points;
    double *y_data = Q;

    ImPlot::PlotLine("max Diffusion", x_data, plotting_data->max_Q, N);
    ImPlot::PlotLine("Diffusion", x_data, y_data, N);
    ImPlot::EndPlot();
  }
}
