set(Patchwork_INC
    analytic_evaluator.hpp 
    bezier.hpp 
    blended_evaluator.hpp 
    bspline.hpp 
    bspline_codegen.hpp 
    bspline_precomp.hpp 
    diff_geom.hpp 
    face_map.hpp 
    #face_map_obj.hpp 
    function.hpp
    p4est_interface.hpp
    p4est_interface.tpp
    polynomial_evaluator.hpp
    surface_evaluator.hpp
    sampling.hpp
)

PREPEND(Patchwork_INC)
set(Patchwork_INC ${Patchwork_INC}  PARENT_SCOPE)
