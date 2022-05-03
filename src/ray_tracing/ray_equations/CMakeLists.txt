
add_library(ray_equations STATIC ray_equations.cpp ray_equations.hpp)

add_subdirectory("equations_of_motion")

target_link_libraries(ray_equations PUBLIC equations_of_motion)
