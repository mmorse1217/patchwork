#ifndef __P4EST_INTERFACE_HPP__
#define __P4EST_INTERFACE_HPP__
#include <gpmesh.hpp>
#include "surface_evaluator.hpp"
// IMPORTANT mpi must be included outside of the extern "C" block, before p4est
// is included. otherwise, the include searches for non-name-mangled functions
// since mpi is included inside the extern "C" and expects C functions, but 
// only finds the c++ MPI, and gets upset.
#include <mpi.h>
extern "C" {
#include <p4est.h>
#include <p4est_iterate.h>
#include <p4est_extended.h>
}


/**
 * struct stored in each p4est_quadrant_t->p.user_data. Keeps a reference to the
 * associated FaceMapPatch object, marks the quad for refinement, and counts the
 * number of refinement steps to arrive at the current quad.
 */
// TODO refactor into a standalone templated library that takes an arbitrary
// patch type. Use this:
// https://stackoverflow.com/questions/26174510/how-to-make-a-function-with-c-linkage-from-template
template<class PatchType>
struct RefinementData {
    bool refine;
    PatchType* patch;
    size_t num_refinement_steps;
};

typedef void (*refinement_callback)(p4est_t* p4est,
        p4est_topidx_t which_tree, 
        int num_outgoing,
        p4est_quadrant_t* outgoing[],
        int num_incoming,
        p4est_quadrant_t* incoming[]);

/**
 * Construct a forest of quadtrees based on a user-provided base quad mesh. 
 * The trees are unrefined i.e. one root-level quad per mesh face. Refinement
 * takes place elsewhere.
 * @param   GpMesh          mesh            user-defined mesh                          
 * @return  p4est_t*                        p4est based on the connectivity
 *                                          structure of GpMesh
 */
template<class PatchType>
p4est_t* build_p4est_from_mesh(GpMesh* mesh);

/**
 * Mark all leaf patches in p4est for refinement. Can mark all as valid as well,
 * but mostly needed for debugging purposes.
 * @param   p4est_t*    p4est        p4est to mark 
 * @param   bool        to_refine    optional bool to mark all quads as valid
 * @return                           p4est with quads marked for refinement
 *                                   according to to_refine
 */
template<class PatchType>
void mark_all_patches_for_refinement(p4est_t* p4est, bool to_refine=true);
/**
 * Check if any of the leaves in p4est need more refinement.
 * @param   p4est_t*    p4est        p4est to inspect 
 * @return  bool                     indicating whether or not any of the p4est
 *                                   quads needs additional refinement
 */
template<class PatchType>
bool check_refinement_criteria(p4est* p4est);

/**
 * Updates the (quad->p.user_data)->refine variable on each leaf quad, according
 * to the external rule is_quad_valid. It is assumed that the refine variable is
 * updated within is_quad_valid; this merely passes the callback along.
 * @param   p4est_t*    p4est                       p4est to inspect 
 * @param   p4est_iter_volume_t* is_quad_valid      rule with which each quad is
 *                                                  updated
 * @return  void                     the final tree will have newly updated 
 *                                   refinement values based on the rule passed
 *                                   by the user
 */
void update_refinement_criteria_on_quads(p4est* p4est, 
        SurfaceEvaluator* evaluator,
        p4est_iter_volume_t is_quad_valid);
/**
 * the same as above, without the evaluator argument.
 */
void update_quad_user_data(p4est* p4est, 
        p4est_iter_volume_t quad_updator);

/**
 * Aggregate void* from all leaf quads in p4est. User needs to recast to desired
 * patch type as needed.
 * @param   p4est_t*    p4est                forest of quadtrees 
 * @return  vector<void*>                    list of "patches" stored on the
 *                                           forest leaves. patches can be
 *                                           anything so long as refinement is
 *                                           well defined
 */
template<class PatchType>
vector<PatchType*> collect_patches(p4est* p4est);

/**
 * Refine p4est uniformly "level" number of times
 * @param   p4est_t*    p4est                forest of quadtrees 
 * @param   int*        level                number of refinement levels to
 *                                           compute
 * @param   refinement_callback callback     callback function passed to
 *                                           p4est_refine_ext; this should
 *                                           populate p.user_data appropriately
 *                                           based on the parent's p.user_data
 * @return  void                     the final tree will have "level" levels of
 *                                   uniform quadtree subdivision
 * */
template<class PatchType>
void uniform_refinement(p4est_t* p4est, int level, refinement_callback callback);

/**
 * Callback passed to p4est_refine_ext to determine whether or not to split a
 * quad.
 * @param   p4est_t*                        unused
 * @param   p4est_topidx_t*                 unused
 * @param   p4est_quadrant_t*   quadrant    current quad of interest
 * @return  int                             whether the quad should be split
 */
template<class PatchType>
int is_quad_valid(p4est_t*, 
        p4est_topidx_t, 
        p4est_quadrant_t* quadrant);

void quad_to_face_map_uv_bounds(p4est_quadrant_t* quad, Rectangle& domain,
        int level_override =-1);
double p4est_quad_to_face_map_coords(int32_t quad_coord, int8_t quad_level);
int32_t face_map_to_p4est_quad_coords(double x, int8_t level);

#include "p4est_interface.tpp"
#endif
