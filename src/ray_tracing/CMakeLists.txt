

add_library(ray_tracing STATIC ray_tracing.cpp ray_tracing.hpp)

add_subdirectory("math_formula")
add_subdirectory("ray_equations")


target_link_libraries(ray_tracing PUBLIC math_formula)
target_link_libraries(ray_tracing PUBLIC ray_equations)


