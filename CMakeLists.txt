cmake_minimum_required(VERSION 3.22)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O1")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -fPIC -c -soname")
set(CMAKE_VERBOSE_MAKEFILE ON)
# set(CMAKE_C_COMPILER "gcc")
# set(CMAKE_CXX_COMPILER "/usr/bin/g++")
# project(dif_lib VERSION 1.0.0)

set(SUNDIALS_INSTALL_PREFIX "~/install/sundials" CACHE PATH "")
set(ACTIVATE_LIVE_PLOTTING ON CACHE BOOL "")
set(ACTIVATE_TESTS ON CACHE BOOL "")

find_package(SUNDIALS REQUIRED)

cmake_minimum_required(VERSION 3.14)
project(my_project)


set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_C_STANDARD 11)
# set(CMAKE_C_STANDARD_REQUIRED ON)

# Enable double assertions

# Include directories
# include_directories(src)

# Add subdirectories
add_subdirectory(src)
if(ACTIVATE_TESTS)
    set(UNITY_INSTALL_PREFIX "~/install/unity" CACHE PATH "")
    find_package(unity REQUIRED)
    add_subdirectory(tests)
endif()


