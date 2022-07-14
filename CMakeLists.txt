cmake_minimum_required(VERSION 3.21.1)
project(black_hole)

set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP REQUIRED)


add_subdirectory(src)


include_directories(src)

add_executable(black_hole src/main.cpp)

add_subdirectory("src/ray_tracing")

target_link_libraries(black_hole PUBLIC ray_tracing)

target_link_libraries(black_hole PRIVATE OpenMP::OpenMP_CXX)
