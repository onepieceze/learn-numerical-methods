cmake_minimum_required(2.8)

project(learn-numerical-method Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "-ffree-line-length-none -g -fbacktrace")
endif

set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set (source
  src/euler_method_mod.F90
)

add_library (numerical_method SHARED ${source})

add_executable (test_euler_method src/test_euler_method.F90)
target_link_libraries (test_euler_method, numerical_method)
