cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
include_directories(${Patchwork_SOURCE_DIR}/src)
#move libraries to build/lib
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )


set(Patchwork_src
    analytic_evaluator.cpp 
    bezier.cpp 
    blended_evaluator.cpp 
    bspline.cpp 
    bspline_codegen.cpp 
    bspline_precomp.cpp 
    diff_geom.cpp 
    face_map.cpp 
    face_map_obj.cpp 
    function.cpp
    p4est_interface.cpp
    polynomial_evaluator.cpp
    sampling.cpp
)
set(Patchwork_header
    analytic_evaluator.hpp 
    bezier.hpp 
    blended_evaluator.hpp 
    bspline.hpp 
    bspline_codegen.hpp 
    bspline_precomp.hpp 
    diff_geom.hpp 
    face_map.hpp 
    face_map_obj.hpp 
    function.hpp
    p4est_interface.hpp
    polynomial_evaluator.hpp
    sampling.hpp
)
add_library(patchwork STATIC ${Patchwork_src} ${Patchwork_header})
#target_link_libraries(patchwork PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(patchwork PUBLIC ${BLENDSURF_LIB} ${P4EST_LIB} ${SC_LIB} 
      ${MPI_LIB} png lapack blas )
 #target_compile_options(patchwork PUBLIC ${OpenMP_CXX_FLAGS})
