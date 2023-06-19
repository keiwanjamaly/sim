#include "live_plotting.h"
#include "computation_data.h"
#include "compute_physics.h"
#include "imgui.h"
#include "implot.h"
#include "physics.h"
#include "u_plotting.h"
#include <filesystem>
#include <stdlib.h>
#include <unistd.h>

static void glfw_error_callback(int error, const char *description) {
  fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

LivePlottingData *setup_live_plotting(ComputationData *data) {
  glfwSetErrorCallback(glfw_error_callback);

  LivePlottingData *new_plotting_data =
      (LivePlottingData *)malloc(sizeof(LivePlottingData));

  new_plotting_data->data = data;
  int N = data->computation_grid->N;
  new_plotting_data->max_Q = (double *)malloc(N * sizeof(double));
  for (int i = 0; i < N; i++) {
    new_plotting_data->max_Q[i] = 0.0;
  }

  std::filesystem::path cwd = std::filesystem::current_path();
  if (!glfwInit())
    exit(1);
  std::filesystem::current_path(cwd);

  // Decide GL+GLSL versions
#if defined(IMGUI_IMPL_OPENGL_ES2)
  // GL ES 2.0 + GLSL 100
  const char *glsl_version = "#version 100";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#elif defined(__APPLE__)
  // GL 3.2 + GLSL 150
  const char *glsl_version = "#version 150";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // 3.2+ only
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);           // Required on Mac
#else
  // GL 3.0 + GLSL 130
  const char *glsl_version = "#version 130";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+
  // only glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // 3.0+ only
#endif

  // Create window with graphics context
  GLFWwindow *window =
      glfwCreateWindow(1280, 720, "FLOW Debug Window", nullptr, nullptr);
  new_plotting_data->window = window;
  if (window == nullptr)
    exit(1);
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImPlot::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  // new_plotting_data->io = io;
  (void)io;
  io.ConfigFlags |=
      ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
  io.ConfigFlags |=
      ImGuiConfigFlags_NavEnableGamepad; // Enable Gamepad Controls

  // Setup Dear ImGui style
  ImGui::StyleColorsDark();
  // ImGui::StyleColorsLight();

  // Setup Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);

  // Our state
  bool show_demo_window = true;
  bool show_another_window = false;
  ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
  return new_plotting_data;
}

void tear_down_live_plotting(LivePlottingData *data) {
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImPlot::DestroyContext();
  ImGui::DestroyContext();

  glfwDestroyWindow(data->window);
  glfwTerminate();
  free(data->max_Q);
  free(data);
}

void draw_frame(LivePlottingData *plotting_data, double time, double *u, double simulation_time) {
  if (glfwWindowShouldClose(plotting_data->window)) {
    printf("Program exited due to user interrupt!\n");
    exit(1);
  }
  // Poll and handle events (inputs, window resize, etc.)
  // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to
  // tell if dear imgui wants to use your inputs.
  // - When io.WantCaptureMouse is true, do not dispatch mouse input data to
  // your main application, or clear/overwrite your copy of the mouse data.
  // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input
  // data to your main application, or clear/overwrite your copy of the
  // keyboard data. Generally you may always pass all inputs to dear imgui,
  // and hide them from your application based on those two flags.
  glfwPollEvents();

  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  {
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImVec2(200, 100));
    if (ImGui::Begin("Simulation stats")) {
      ImGui::Text("RG time is: %.3f", time);
      ImGui::Text("k time is: %.4e", cal_k(time, plotting_data->data));
      ImGui::Text("Simulation time is: %.1f", simulation_time);
      ImGui::End();
    }

    plot_u(plotting_data, u);

    if (activate_diffusion(plotting_data->data->data))
      plot_Q(u, time, plotting_data);
  }

  // Rendering
  ImGui::Render();
  ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
  int display_w, display_h;
  glfwGetFramebufferSize(plotting_data->window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w,
               clear_color.z * clear_color.w, clear_color.w);
  glClear(GL_COLOR_BUFFER_BIT);
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

  glfwSwapBuffers(plotting_data->window);
}
