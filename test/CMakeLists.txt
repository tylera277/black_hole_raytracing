

set(CATCH_INCLUDE_DIR </catch.hpp>)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

# Add test executable
add_executable(tests_raytracing test_raytracing.cpp ../src/ray_tracing/ray_tracing.cpp ../src/ray_tracing/ray_tracing.hpp ../src/ray_tracing/math_formula/math_formula.cpp ../src/ray_tracing/math_formula/math_formula.hpp)

add_executable(tests_mathformula  test_mathformula.cpp ../src/ray_tracing/ray_tracing.cpp ../src/ray_tracing/ray_tracing.hpp ../src/ray_tracing/math_formula/math_formula.cpp ../src/ray_tracing/math_formula/math_formula.hpp)



target_link_libraries(tests_raytracing Catch )
target_link_libraries(tests_mathformula Catch)
