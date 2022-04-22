

cmake_minimum_required(VERSION 3.21.1)
project(black_hole)

add_subdirectory("ray_tracing")

add_executable(black_hole main.cpp)


target_link_libraries(black_hole PUBLIC ray_tracing)
