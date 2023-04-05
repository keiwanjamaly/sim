cmake_minimum_required(VERSION 3.22)
project(sim C)

set(CMAKE_C_STANDARD 11)

message(CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}")

find_package(SUNDIALS REQUIRED)

add_executable(sim
        computation_library.c)

target_link_libraries(sim SUNDIALS::cvode m)
