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
#include <GLFW/glfw3.h> // Will drag system OpenGL headers

GLFWwindow *setup_live_plotting();
void tear_down_live_plotting(GLFWwindow *);
void draw_frame(GLFWwindow *, double);

#endif // !LIVE_PLOTTING_H
