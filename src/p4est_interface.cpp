#include "p4est_interface.hpp"

// ----------------------------------------------------------------------------
// Simple functions that just pass along a function handle to p4est_iterate.
// ----------------------------------------------------------------------------
void update_refinement_criteria_on_quads(p4est* p4est, 
        SurfaceEvaluator* evaluator,
        p4est_iter_volume_t is_quad_valid){

    p4est_iterate(p4est,
            NULL,
            evaluator,
            is_quad_valid,
            NULL,
            NULL);
}

void update_quad_user_data(p4est* p4est, 
        p4est_iter_volume_t quad_updator){

    p4est_iterate(p4est,
            NULL,
            NULL,
            quad_updator,
            NULL,
            NULL);
}



double p4est_quad_to_face_map_coords(int32_t quad_coord, int8_t quad_level){
    int face_map_coord_as_int = quad_coord >> (P4EST_MAXLEVEL - quad_level);
    double face_map_coord = double(face_map_coord_as_int)/double(pow2(quad_level));
    return face_map_coord;
}

int32_t face_map_to_p4est_quad_coords(double x, int8_t level){
    int32_t inv = (int(x*pow2(level) ) << (P4EST_MAXLEVEL - level) );
    return inv;
}

void quad_to_face_map_uv_bounds(p4est_quadrant_t* quad, Rectangle& domain, int level_override){

    int8_t level;
    if(level_override >= 0){
        level = level_override;
    } else {
        level= quad->level;
    }
    double quad_width = 1./double(pow2(level));

    double x_min = p4est_quad_to_face_map_coords(quad->x, level);
    double y_min = p4est_quad_to_face_map_coords(quad->y, level);
    
    double x_max = x_min + quad_width;
    double y_max = y_min + quad_width;

    domain.first.first  = x_min;
    domain.first.second = x_max;

    domain.second.first  = y_min;
    domain.second.second = y_max;

}

