# cmake_minimum_required(VERSION 3.22)

# project(dif_lib VERSION 1.0.0)

# include(CheckIPOSupported)
# check_ipo_supported(RESULT ipo_supported)
# if(ipo_supported)
#         set(CMAKE_INTERPROCEDURAL_OPTIMIZATION True)
# endif()

# set(CMAKE_C_STANDARD 11)

# find_package(SUNDIALS REQUIRED)

# set(MODEL_FILE "gross_neveu.c" CACHE STRING "" FORCE)

# add_library(sim SHARED
#         c_code/computation_library.c
#         c_code/computation_data.c
#         c_code/grid.c
#         c_code/return_data.c
#         c_code/${MODEL_FILE})

# target_link_libraries(sim SUNDIALS::cvode m)

cmake_minimum_required(VERSION 3.14)
project(my_project)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Enable double assertions
find_package(unity REQUIRED)

# Include directories
include_directories(src)

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)


