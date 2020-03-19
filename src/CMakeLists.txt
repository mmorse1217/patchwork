cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)

set(Patchwork_SRC
    analytic_evaluator.cpp 
    bezier.cpp 
    blended_evaluator.cpp 
    bspline.cpp 
    bspline_codegen.cpp 
    bspline_precomp.cpp 
    diff_geom.cpp 
    face_map.cpp 
    #face_map_obj.cpp 
    function.cpp
    p4est_interface.cpp
    polynomial_evaluator.cpp
    sampling.cpp
)



PREPEND(Patchwork_SRC)
set(Patchwork_SRC ${Patchwork_SRC}  PARENT_SCOPE)
