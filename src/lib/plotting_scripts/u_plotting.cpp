#include "u_plotting.h"
#include "physics.h"

void plot_u(Grid *computation_grid, double *u) {
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

void plot_Q(Grid *computation_grid, double *u) {
    
}
