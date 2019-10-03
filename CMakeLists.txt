cmake_minimum_required(VERSION 2.8)

project(learn-numerical-method Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "-ffree-line-length-none -g -fbacktrace")
endif ()

set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

set (sources
  src/euler_method_mod.F90
  src/improved_euler_method_mod.F90
)

add_library (numerical_method SHARED ${sources})
target_link_libraries (numerical_method)

add_executable (solve_ode src/solve_ode.F90)
target_link_libraries (solve_ode numerical_method)
