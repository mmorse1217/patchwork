cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
#include_directories($ENV{BLENDSURF_DIR})
# export MPI_DIR=$PETSC_DIR/$PETSC_ARCH
#include_directories($ENV{MPI_DIR}/include)
#include_directories($ENV{P4EST_DIR}/include)
# export OPENGL_DIR=/usr/X11R6/include on mac
#include_directories($ENV{OPENGL_DIR})
#include_directories(${Patchwork_SOURCE_DIR}/src)
#include_directories(${Patchwork_SOURCE_DIR}/tests)
#include_directories(/usr/local/include /usr/include)
#find_package(OpenMP REQUIRED)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()
#find_package(lapack REQUIRED)
#find_package(png REQUIRED)
set(PATCHWORK_TEST_SRC
    test_benchmark.cpp
    test_bezier.cpp
    test_function.cpp
    test_patchwork.hpp
    utils.hpp
    test_geometry_error.cpp
    test_refinement.cpp
    test_regression.cpp
    test_second_derivatives.cpp
    test_spline_error_convergence.cpp
    utils.cpp
)
set(PATCHWORK_TEST_HEADER
    test_patchwork.hpp
    utils.hpp
)
#move executable to build/bin
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
add_executable(test_patchwork test_patchwork.cpp ${PATCHWORK_TEST_SRC})
#target_compile_options(test_patchwork PUBLIC ${PATCHWORK_TEST_SRC})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lstdc++")
target_link_libraries(test_patchwork patchwork)
