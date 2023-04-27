cmake_minimum_required(VERSION 3.22)

# project(dif_lib VERSION 1.0.0)

set(SUNDIALS_INSTALL_PREFIX "~/install/sundials" CACHE PATH "" FORCE)
set(UNITY_INSTALL_PREFIX "~/install/unity" CACHE PATH "" FORCE)

find_package(SUNDIALS REQUIRED)

cmake_minimum_required(VERSION 3.14)
project(my_project)

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_C_STANDARD 11)
# set(CMAKE_C_STANDARD_REQUIRED ON)

# Enable double assertions
find_package(unity REQUIRED)

# Include directories
include_directories(src)

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)


