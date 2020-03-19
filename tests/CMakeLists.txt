cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)


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

target_link_libraries(test_patchwork patchwork
    ${MPI_LIB}
    ${BLENDSURF_LIB} ${P4EST_LIB} ${SC_LIB} ${LAPACK_LIBRARIES} patchwork) 
