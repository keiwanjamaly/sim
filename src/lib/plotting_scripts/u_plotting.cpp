#include "u_plotting.h"
#include "compute_physics.h"
#include "implot.h"
#include "physics.h"
#include <algorithm>
#include <stdlib.h>

void plot_u(LivePlottingData *plotting_data, double *u) {
  Grid *computation_grid = plotting_data->data->computation_grid;
  ImGui::SetNextWindowPos(ImVec2(200, 0));
  ImGui::Begin("Plotting");
  if (ImPlot::BeginPlot("Function Plot")) {
    int N = computation_grid->N;
    double *x_data = computation_grid->grid_points;
    double *y_data = u;

    ImPlot::PlotLine("u", x_data, y_data, N);
    ImPlot::EndPlot();
  }
  ImGui::End();
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

  ImGui::SetNextWindowPos(ImVec2(200, 0));
  ImGui::SetNextWindowSize(ImVec2(1000, 1000));
  ImGui::Begin("Diffusion Plot");
  // ImPlot::SetNextAxisLimits(1, 0, maximum_of_Q);
  if (ImPlot::BeginPlot("Function Plot")) {
    ImPlot::SetupAxis(ImAxis_X1, "sigma");
    ImPlot::SetupAxis(ImAxis_Y1, "");
    ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, maximum_of_Q * 1.1, ImPlotCond_Always);
    double *x_data = computation_grid->grid_points;
    double *y_data = Q;

    ImPlot::PlotLine("max Diffusion", x_data, plotting_data->max_Q, N);
    ImPlot::PlotLine("Diffusion", x_data, y_data, N);
    ImPlot::EndPlot();
  }
  ImGui::End();
}
