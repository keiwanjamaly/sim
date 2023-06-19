#ifndef U_PLOTTING_H
#define U_PLOTTING_H
#include "imgui.h"
#include "implot.h"
#include "grid.h"
#include "live_plotting.h"

void plot_u(LivePlottingData *plotting_data, double *u);

void plot_u_within_zero_to_one(LivePlottingData *plotting_data, double *u);

void plot_Q(double *u, double t, LivePlottingData *plotting_data);

#endif // !U_PLOTTING_H
