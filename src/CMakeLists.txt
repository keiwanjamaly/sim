set(LIBRARY_SOURCES 
  grid.cpp
  return_data.cpp
  computation_data.cpp
  computation_library.cpp
  helpers.cpp
)

set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/imgui)
set(IMGUI_INCLUDE ${IMGUI_DIR} ${IMGUI_DIR}/backends)
set(IMGUI_SOURCE 
    ${IMGUI_DIR}/imgui.cpp 
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    )

set(IMPLOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/implot)
set(IMPLOT_INCLUDE ${IMPLOT_DIR})
set(IMPLOT_SOURCE 
    ${IMPLOT_DIR}/implot.cpp
    ${IMPLOT_DIR}/implot_items.cpp
    )

add_subdirectory(lib/plotting_scripts)

message(LIVE_PLOTTING_INCLUDE=${LIVE_PLOTTING_INCLUDE})

add_library(sim STATIC ${LIBRARY_SOURCES})
target_link_libraries(sim SUNDIALS::cvode m)
target_include_directories(sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if (ACTIVATE_LIVE_PLOTTING)
    find_package(glfw3 3.3 REQUIRED)
    find_package(OpenGL REQUIRED)
    target_sources(sim PRIVATE ${LIVE_PLOTTING_SOURCE} ${IMGUI_SOURCE} ${IMPLOT_SOURCE})
    target_include_directories(sim PUBLIC ${IMGUI_INCLUDE} ${IMPLOT_INCLUDE} ${LIVE_PLOTTING_INCLUDE})
    target_link_libraries(sim glfw OpenGL::GL)
endif()

file(GLOB MODELS ${CMAKE_CURRENT_SOURCE_DIR}/models/*)
foreach(FILE ${MODELS})
  get_filename_component(FILENAME ${FILE} NAME)
  get_filename_component(LIBNAME ${FILE} NAME_WE)
  message("File: ${FILENAME}")
  message("Libname: lib${LIBNAME}")
  add_library(${LIBNAME} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/models/${FILENAME} ${LIBRARY_SOURCES})
  target_link_libraries(${LIBNAME} PUBLIC m SUNDIALS::cvode sim)
  target_include_directories(${LIBNAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
  if (ACTIVATE_LIVE_PLOTTING)
      target_sources(${LIBNAME} PRIVATE ${LIVE_PLOTTING_SOURCE} ${IMGUI_SOURCE} ${IMPLOT_SOURCE})
      target_include_directories(${LIBNAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${IMGUI_INCLUDE} ${IMPLOT_INCLUDE} ${LIVE_PLOTTING_INCLUDE})
      target_compile_definitions(${LIBNAME} PUBLIC ACTIVATE_LIVE_PLOTTING)
      target_compile_definitions(${LIBNAME} PRIVATE ACTIVATE_LIVE_PLOTTING)
  endif()
endforeach()
