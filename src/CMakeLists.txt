add_library(grid_code STATIC grid.c)
target_include_directories(grid_code PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})