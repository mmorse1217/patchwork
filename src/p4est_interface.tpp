#include "p4est_interface.hpp"

template<class PatchType>
p4est_t* build_p4est_from_mesh(GpMesh* mesh){

    p4est_topidx_t num_vertices = mesh->numVertices();
    p4est_topidx_t num_trees = mesh->numFaces();

    // MJM TODO BUG DOESN'T HANDLE DOMAINS WITH CORNERS
    p4est_topidx_t num_ctt = 0;

    const int dim = 3;
    const int quad_size = 4;
    double vertices[num_vertices * dim];
    p4est_topidx_t tree_to_vertex[num_trees * quad_size];
    p4est_topidx_t tree_to_tree[num_trees * quad_size];
    int8_t        tree_to_face[num_trees * quad_size];

    // Copy vertices from mesh array into p4est array
    for(int v = 0; v < num_vertices; v++){
        Point3 vertex = mesh->vpoint(v);
        for(int d =0; d < dim; d++){
            vertices[v*dim + d] = vertex(d);
        }
    }

    // GpMesh vertices are ordered clockwise around a face:
    //  1 ---> 2
    //  ^      |
    //  |      |
    //  |      \/
    //  0 <--- 3
    // (0, 1, 2, 3)  

    // P4est faces are ordered in z-ordering:
    //  1 ---> 2 
    //    ^      
    //     |   
    //      | 
    //  0 ---> 3
    // (0,3,1,2)
    // The following map will map local vertex id's in clockwise order to local
    // vertex id's in z-order (Figure 2-left in p4est paper).
    vector<int> clockwise_to_z_order_vertex(quad_size, -1);
    clockwise_to_z_order_vertex[0] = 0;
    clockwise_to_z_order_vertex[1] = 2;
    clockwise_to_z_order_vertex[2] = 3;
    clockwise_to_z_order_vertex[3] = 1;

    // For each quad, indicate the global vertex id's that bound the quad (and
    // hence the tree that lives on that quad)
    for(int F = 0; F < num_trees; F++){
        for(int v = 0; v < quad_size; v++){ 
            int z_order_v = clockwise_to_z_order_vertex[v];
            pair<int, int> ids = mesh->Fv2Vf(F, v);
            
            int V = ids.first; // global vertex id
            int index = F*4 + z_order_v;
            tree_to_vertex[index] = V;
        }
    }

    p4est_connectivity_t* conn = p4est_connectivity_new_copy (
            num_vertices, num_trees, 0,
            vertices, tree_to_vertex,
            tree_to_tree, tree_to_face,
            NULL, &num_ctt, NULL, NULL);

    // fill in remaining arrays because who knows what tree_to_face is suppposed
    // to be
    p4est_connectivity_complete(conn);
    p4est_t* p4est = 
        p4est_new(MPI_COMM_WORLD, conn, sizeof(RefinementData<PatchType>*), NULL, NULL);
    return p4est;


}

// ----------------------------------------------------------------------------
template<class PatchType>
void mark_all_quads(p4est_iter_volume_info_t * info,
        void *user_data){
    p4est_quadrant_t* quad = info->quad;
    RefinementData<PatchType>* r = static_cast<RefinementData<PatchType>*>(quad->p.user_data);
    bool* to_refine = (bool*) user_data;
    r->refine = *to_refine;
}

template<class PatchType>
void mark_all_patches_for_refinement(p4est_t* p4est, bool to_refine){
    p4est_iterate(p4est,
            NULL,
            &to_refine,
            mark_all_quads<PatchType>,
            NULL,
            NULL);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template<class PatchType>
struct CheckRefinementCriteria{
    public:
        // store the result
        bool _further_refinement_needed;

        CheckRefinementCriteria(): 
            _further_refinement_needed(false) {}// initially, no refinement

        void operator()(p4est_iter_volume_info_t * info,
                void *user_data){
            p4est_quadrant_t* quad = info->quad;
            RefinementData<PatchType>* r = static_cast<RefinementData<PatchType>*>(quad->p.user_data);

            // if any of the quads we touch need refinement, we need another
            // pass at the loop
            _further_refinement_needed = _further_refinement_needed || r->refine;
        }
};
/*
// used to allow for templated functions to be passed as an extern C callback...
// maybe
extern "C" typedef void p4est_iterate_quad(
        p4est_iter_volume_info_t * info,
        void *user_data);
template <class PatchType> 
p4est_iterate_quad check_refinement_wrapper;

template <class PatchType> 
*/

template<class PatchType>
void check_refinement_wrapper(p4est_iter_volume_info_t * info,
                    void *user_data){
    CheckRefinementCriteria<PatchType>* check_ref = 
        static_cast<CheckRefinementCriteria<PatchType>*>(user_data);
    return(*check_ref)(info, user_data);
}

template<class PatchType>
bool check_refinement_criteria(p4est* p4est){
   CheckRefinementCriteria<PatchType> check_refinement;

    p4est_iterate(p4est,
            NULL,
            &check_refinement,
            check_refinement_wrapper<PatchType>,
            NULL,
            NULL);
    return check_refinement._further_refinement_needed;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// TODO Not parallel-ready; assumes forest is contained on one machine. need to
// communicate patches across nodes to synchronize after refinement terminates
template<class PatchType>
struct PatchAggregator{
    public:
        // store the result
        vector<PatchType*> _patches;

        PatchAggregator(){}

        void operator()(p4est_iter_volume_info_t * info,
                void* /* user_data*/){
            p4est_quadrant_t* quad = info->quad;
            RefinementData<PatchType>* r = static_cast<RefinementData<PatchType>*>(quad->p.user_data);

            _patches.push_back(r->patch);
        }
};

template<class PatchType>
void patch_aggregator_wrapper(p4est_iter_volume_info_t * info,
                    void *user_data){
    PatchAggregator<PatchType>* patch_aggregator= 
        static_cast<PatchAggregator<PatchType>*>(user_data);
    return(*patch_aggregator)(info, user_data);
}

template<class PatchType>
vector<PatchType*>  collect_patches(p4est* p4est){
    PatchAggregator<PatchType> patch_aggregator;

    p4est_iterate(p4est,
            NULL,
            &patch_aggregator,
            patch_aggregator_wrapper<PatchType>,
            NULL,
            NULL);
    return patch_aggregator._patches;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template<class PatchType>
void uniform_refinement(p4est_t* p4est, int level, refinement_callback callback){
    //assert(level > 0);
    mark_all_patches_for_refinement<PatchType>(p4est);
    for(int l =0; l < level; l++){
        p4est_refine_ext(p4est, 0, -1, is_quad_valid<PatchType>, NULL, callback);
    }
}

template<class PatchType>
int is_quad_valid(p4est_t*, 
        p4est_topidx_t, 
        p4est_quadrant_t* quadrant){
  return ((RefinementData<PatchType>*)quadrant->p.user_data)->refine ;
}




