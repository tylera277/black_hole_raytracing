

add_library(ray_tracing STATIC ray_tracing.cpp ray_tracing.hpp)

add_subdirectory("math_formula")

target_link_libraries(ray_tracing PUBLIC math_formula)

include_directories("/usr/local/include/eigen3")

