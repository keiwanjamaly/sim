#ifndef LIVE_PLOTTING_H
#define LIVE_PLOTTING_H

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <stdio.h>
#define GL_SILENCE_DEPRECATION
#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <GLES2/gl2.h>
#endif
#include "computation_data.h"
#include "grid.h"
#include <GLFW/glfw3.h> // Will drag system OpenGL headers

typedef struct live_plotting_data {
  GLFWwindow *window;
  // ImGuiIO &io;
  ComputationData *data;
  double *max_Q;
} LivePlottingData;

LivePlottingData *setup_live_plotting(ComputationData *data);
void tear_down_live_plotting(LivePlottingData *data);
void draw_frame(LivePlottingData *plotting_data, double, double *, double);

#endif // !LIVE_PLOTTING_H
