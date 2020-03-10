#include "face_map.hpp"
#include <set>
#include <queue>
#include "blended_evaluator.hpp"
#include <limits.h>
#include <omp.h>
#include <algorithm>
#include "sampling.hpp"
#include "polynomial_evaluator.hpp"
#include "diff_geom.hpp"
typedef pair<int, int> intpair;


FaceMapSurf::FaceMapSurf(SurfaceEvaluator* evaluator):
    _fitting_error(1e-4),
    _evaluator(evaluator), 
    _use_functional_rep(false){
    _mesh = *(evaluator->mesh());
    _legacy = true;
    
    _p4est = build_p4est_from_mesh<FaceMapPatch>(&_mesh);
}

void FaceMapSurf::setup(){

    int num_faces= _mesh.numFaces();

    // TODO might be able to remove this now
    if(_adaptive){
        _refinement_factor = 0;
    }

    // Precompute (n choose i) for i=0, ..., n
    _binomial_coefs_array = compute_binomial_coefs_array(_patch_order);
    _basis_type = Common::BEZIER;
    //_basis_type = Common::BSPLINE;
    // list of invalid patches that need to be subdivided further
    queue<FaceMapPatch*> invalid_patches;

    vector<FaceMapPatch*> valid_patches;
    
    // group ids of parent surface (from wrl file)
    vector<int> blended_group_ids = gpmesh().get_face_group_ids();
    vector<int> boundary_orientations= gpmesh().get_intpt_orientation();
   
    // 3d samples used to fit each patches to surface
    // outer vector is size num_faces; inner vectors are size (num_samples^2)
    _samples = new vector<vector<Point3> >(0, vector<Point3>() );
    cout << "gpmesh faces: " << num_faces << endl;

    // Initial pass: all top level patches derived directly from a mesh will be
    // split once
    int num_samples_per_dim = 4*_patch_order;
    double ctrlstep = 1.0/float(num_samples_per_dim-1);
    if(!_legacy){
        NumMatrix uv_samples(2, num_samples_per_dim*num_samples_per_dim);
        for(int j = 0; j < num_samples_per_dim; j++){
            for(int i = 0; i < num_samples_per_dim; i++){
                // Equal sized steps in each direction over [0,1] x [0,1]
                double cd[2];
                cd[0] = i*ctrlstep;
                cd[1] = j*ctrlstep;

                for(int d =0; d < uv_samples.m(); d++)
                    uv_samples(d,j*(num_samples_per_dim)+i) = cd[d];
                //Function::compute_basis_function_value_cache(2, _patch_order, *_binomial_coefs, NumVector(2, false, cd), this->cache);
            }
        }

        if(_basis_type == Common::BSPLINE){
            int N = 2;
            int M = 3; 
            int spline_piece_degree = 5;
            Function f(N, M, spline_piece_degree, _patch_order+1, _basis_type);
            _pseudo_inverse = f.compute_pseudo_inverse(uv_samples);

        } else if (_basis_type == Common::BEZIER){
        // For bezier...
        _pseudo_inverse = 
            Function::compute_pseudo_inverse(_patch_order, 
                    _patch_order+1, 
                    *_binomial_coefs_array, 
                    _basis_type, 
                    uv_samples);
        }
    }

    for(int F = 0; F < num_faces; F++){

        // make sure the quadtree forest is initially unrefined for all
        // quadtrees

        p4est_tree_t* tree = p4est_tree_array_index(_p4est->trees, F);
        //cout << "num quads in tree " << F << ": " << tree->quadrants.elem_count << endl;
        //assert(tree->quadrants.elem_count == 1);
        if(tree->quadrants.elem_count > 0){
        p4est_quadrant_t* quad = p4est_quadrant_array_index(&tree->quadrants, 0);
        
        quad->p.user_data = new RefinementData<FaceMapPatch>(); 
        // Sub-intverals of parametric domain of blendsurf patch for this patch
        // to sample and fit to
        Interval x_int(0., 1.);
        Interval y_int(0., 1.);

        // we're on the root (0) level w.r.t. to the forest of quadtrees
        int base_level = 0; 
        
        // Build a patch
        FaceMapPatch* patch;
        if(!_legacy){
            patch = new FaceMapPatch(F, base_level, x_int, y_int, _basis_type,this);
        } else {
            patch = new FaceMapPatch(F, base_level, x_int, y_int, this);
            patch->set_cache(&_pseudo_inverse, NULL);
        }

        patch->_group_id = blended_group_ids[F];
        patch->_orientation= boundary_orientations[patch->_group_id];

        // Let the quad know which patch its associated with and vice versa
        RefinementData<FaceMapPatch>* r = 
            static_cast<RefinementData<FaceMapPatch>*>(quad->p.user_data);
        r->patch = patch;
        patch->_quad = quad;

        // Store it for processing later
        invalid_patches.push(patch);
        }
    }
    // adaptive geometric refinement: fit face-map patches to underlying
    // mesh-based surface
    if(_adaptive){
        // split patches until surface is resolve to _fitting_error  rel. err.

        clock_t eval_time =omp_get_wtime(); 
        resolve_function_tree(_evaluator, 
                            _fitting_error, 
                            _fitting_error);
        valid_patches = _patches;
        /*
        //valid_patches = extract_patches_from_quadtrees();
        valid_patches = resolve_function(_evaluator, 
                            invalid_patches, 
                            _fitting_error, 
                            _fitting_error);
        */
        eval_time = omp_get_wtime() - eval_time;
        cout << endl << "Adaptive geometry refinement time: " << eval_time << endl << endl; 
    } else {
        // uniform geometric refinement: split all top-level patches
        // _refinement_factor times into (# num top-level
        // patches)*4^_refinement_factor) new patches
        uniform_refinement<FaceMapPatch>(_p4est, _refinement_factor, build_child_patches_from_parent);
        _patches = collect_patches<FaceMapPatch>(_p4est);;
        valid_patches = _patches;
    }

    // Resize the interface to blendsurf and the legacy plotting code to handle
    // the refined number of patches
    _samples->resize(valid_patches.size());

    cout << "num patches: " << _samples->size() << endl;
    // Update internal data structures
    for(int i =0; i < valid_patches.size(); i++){
        FaceMapPatch* patch = valid_patches[i];
    //cout << "mean_curvature: " << Differential::mean_curvature(Point2(.5, .5), this, i) << endl;
        (*_samples)[i] = (*patch)._samples;
    }
    _patches = valid_patches;
}
FaceMapPatch::FaceMapPatch(int F, int level, 
        Interval x_int, Interval y_int, 
        Common::BasisType basis_type, 
        FaceMapSurf* surface):
    // assumes function from R^2 \to \R^3 if not specified, see function.cpp
    //Function(surface->_patch_order, basis_type){
    Function(2, 3, 
            basis_type == Common::BSPLINE? 5: surface->_patch_order, 
            surface->_patch_order+1, 
            basis_type),
    _level_of_refinement(level),
    _blended_face_to_sample(F), 
    _surface(surface),
    _x_int(x_int),
    _y_int(y_int),
    _patch_order(surface->_patch_order){

         // TODO THIS IS BUG PRONE FIX BY FACTORING OUT INTO OPTION   
        int num_samples_per_dim =  2*surface->_patch_order;

        _xy_coords_mat = NumMatrix(_domain_dim, num_samples_per_dim*num_samples_per_dim);
        _samples_mat = NumMatrix(_range_dim, num_samples_per_dim*num_samples_per_dim);
        // sample the surface
        surface->compute_xy_coordinates_and_surface_positions(
                this,
                _xy_coords_mat, 
                _samples_mat);
        
        set_cache(&(surface->_pseudo_inverse), NULL);
        //_basis_function_values = &(surface->cache);
        // approximate in least squares sense
        this->compute_polynomial_approximation(_xy_coords_mat, _samples_mat);
        
        _derivatives.resize(5);

        // Compute derivatives
        Function dp_du = Function::d_dx(Common::U,*this);
        Function dp_dv = Function::d_dx(Common::V,*this);
        Function d2p_du2 = Function::d_dx(Common::U, dp_du);
        Function d2p_dudv = Function::d_dx(Common::V, dp_du);
        Function d2p_dv2  = Function::d_dx(Common::V,dp_dv);

        _derivatives[0] = dp_du;
        _derivatives[1] = dp_dv;
        _derivatives[2] = d2p_du2;
        _derivatives[3] = d2p_dudv;
        _derivatives[4] = d2p_dv2;
        // support legacy interface to viewer... 
        _samples.resize(_samples_mat.n());
        for(int i = 0; i < _samples_mat.n(); i++){
            Point3 sample(_samples_mat.clmdata(i));
            _samples[i] = sample;
        }

        //enforce_periodicity();
        Vector<Point3>* coefficients = new Vector<Point3>(_polynomial_coefficents.n());
        for(int i = 0; i < _polynomial_coefficents.n(); i++){
            Point3 coeff(_polynomial_coefficents.clmdata(i));
            (*coefficients)(i) = coeff;
        }
        _polynomial_patch = coefficients;
}

FaceMapPatch::FaceMapPatch(int F,
        int level,
        Interval x_int,
        Interval y_int,
        FaceMapSurf* surface):
    _level_of_refinement(level),
    _blended_face_to_sample(F), 
    _surface(surface),
    _x_int(x_int),
    _y_int(y_int),
    _patch_order(surface->_patch_order){

    assert(surface != NULL);
    
    // Uniformly sample the surface face 
    surface->compute_xy_coordinates_and_blended_positions(
            this,
            _xy_coords, 
            _samples);
    PolynomialEvaluator* e = dynamic_cast<PolynomialEvaluator*>(surface->_evaluator);
    if(e && !_surface->_adaptive){
        // copy coefficients directly
        Vector<Point3>* explicit_coeffs = e->polynomial_patch(F);
        Vector<Point3>* copied_coeffs = new Vector<Point3>(explicit_coeffs->length());
        for(int i = 0; i < explicit_coeffs->length(); i++){
            Point3 coeff = (*explicit_coeffs)(i);
            for(int d = 0; d < 3; d++)
            (*copied_coeffs)(i)(d) = coeff(d);
        }
        _polynomial_patch = copied_coeffs;
    } else {

    // Fit a polynomial to the samples on the surface face
    _polynomial_patch = surface->construct_bernstein_polynomial_patch(
            _xy_coords,
            _samples);
    }
}


vector<FaceMapPatch*> FaceMapPatch::split(){
    Interval x_int = this->_x_int;
    Interval y_int = this->_y_int;

    // Compute the x- and y-midpoints of the current patch. These along with the
    // original parent domain bounds define the upper/lower bounds of the
    // subpatches
    double x_mid = (x_int.second + x_int.first)/2.;
    double y_mid = (y_int.second + y_int.first)/2.;

    vector<FaceMapPatch*> subpatches;

    for(int i =0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            Interval x_split;
            Interval y_split;
            
            // Enumerate the bounds of the four subpatches
            //  ------------------ <-y_int.second
            //  |        |       | 
            //  |        |       | 
            //  |        |       | 
            //  ------------------ <- y_mid 
            //  |        |       | 
            //  |        |       | 
            //  |        |       | 
            //  ------------------ <-y_int.first
            //  ^        ^-x_mid ^-x_int.second
            //  |-x_int.first
            
            if (i == 0){
                // lower-left patch
                x_split.first = x_int.first;
                x_split.second= x_mid;
            } else {
                // lower right patch 
                x_split.first = x_mid; 
                x_split.second= x_int.second;
            }
            
            if (j == 0){
                // upper-left patch
                y_split.first = y_int.first;
                y_split.second= y_mid;
            } else {
                // upper-right patch
                y_split.first = y_mid;
                y_split.second= y_int.second;
            }
            
            // Make actually make new patch on the next level down defined on
            // x_split x y_split
            FaceMapPatch* patch;
            if(!_surface->_legacy){
                patch = new FaceMapPatch(
                        this->_blended_face_to_sample, // sample blended face as parent
                        this->_level_of_refinement+1, // next level of refinement
                        x_split,
                        y_split, 
                        _surface->_basis_type,
                        this->_surface);
            } else {
                patch = new FaceMapPatch(
                        this->_blended_face_to_sample, // sample blended face as parent
                        this->_level_of_refinement+1, // next level of refinement
                        x_split,
                        y_split, 
                        //Common::BEZIER,
                        this->_surface);
            }
            patch->_group_id = _group_id;
            patch->_orientation = _orientation;
            subpatches.push_back(patch);

        }
    }


    return subpatches;

}

vector<FaceMapPatch*> FaceMapPatch::split(vector<p4est_quadrant_t*> child_quads){
    assert(child_quads.size() == 4);

    vector<FaceMapPatch*> subpatches;

    for(int qi = 0; qi < child_quads.size(); qi++){        
        p4est_quadrant_t* child_quad = child_quads[qi];

        // Generate the new child's subdomain of the parent directly from the
        // quad data
        Rectangle child_domain;
        quad_to_face_map_uv_bounds(child_quad, child_domain);
        // Make actually make new patch on the next level down defined on that
        // domain
        FaceMapPatch* patch;
        if(!_surface->_legacy){
            patch = new FaceMapPatch(
                    this->_blended_face_to_sample, // sample blended face as parent
                    this->_level_of_refinement+1, // next level of refinement
                    child_domain.first,
                    child_domain.second, 
                    _surface->_basis_type,
                    this->_surface);
        } else {
            patch = new FaceMapPatch(
                    this->_blended_face_to_sample, // sample blended face as parent
                    this->_level_of_refinement+1, // next level of refinement
                    child_domain.first,
                    child_domain.second, 
                    //Common::BEZIER,
                    this->_surface);
        }
        patch->_group_id = _group_id;
        patch->_orientation = _orientation;
        subpatches.push_back(patch);

    }


    return subpatches;

}





bool FaceMapPatch::is_patch_valid(SurfaceEvaluator* evaluator,
        double eps_pos, double eps_normal){

    // we need to check the patch value at at least _surface->patch_order+1
    // points to resolve the polynomial, but we don't want to check at the same
    // points we fit to. So we add two extra points 
    int num_samples_per_dim = 2*_surface->_patch_order;


    // We need the surface bounding box in order to properly normalize the
    // distance from a point on a face-map patch to the original surface 
    Point3 bbox_max;
    Point3 bbox_min;
    //_surface->_evaluator->bounding_box(bbox_min, bbox_max);
    evaluator->bounding_box(bbox_min, bbox_max);

    // buffer the size of the bounding box to prevent numerical issues for
    // control points defined on the boundary of [-1,1]^3
    Point3 bbox_diameter = bbox_max - bbox_min+ .001*Point3(1,1,1);

    // relative distance is determined by dividing absolute distance by the 
    // longest edge of the bounding box.
    double normalization_factor = bbox_diameter.x();
    normalization_factor = max(normalization_factor, bbox_diameter.y());
    normalization_factor = max(normalization_factor, bbox_diameter.z());

    // for each tensor product sample, compute the error of the polynomial
    // approximation to the underlying surface representation at each point.
    double position_error = 0.;
    double normal_error = 0.;
    
    BlendedEvaluator* blended_evaluator = dynamic_cast<BlendedEvaluator*>(evaluator);
    bool is_blended_evaluator = blended_evaluator != NULL;
    
    Rectangle subpatch_domain(_x_int, _y_int);
    
    // Sample the patch you're going to split the current patch into 4 patches,
    // defined on [0,1].
    NumMatrix face_map_equispaced_samples = 
        Sampling::sample_2d<Sampling::equispaced>(num_samples_per_dim, Sampling::base_domain);
    
    assert(face_map_equispaced_samples.m() == Sampling::parameter_dim);
    assert(face_map_equispaced_samples.n() == num_samples_per_dim*num_samples_per_dim);
    
    // sample the four subdomains of the quad of interest
    for(int i =0 ; i < num_samples_per_dim; i++){
        for(int j =0 ; j < num_samples_per_dim; j++){
            int index = i*num_samples_per_dim + j;
            
            Point2 xy(face_map_equispaced_samples.clmdata(index));

            // translate uv coordinate in patch domain [0,1]^2 to the mesh 
            // domain containing it
            Point2 cd = Sampling::affine_coordinate_transform(
                    xy,
                    Sampling::base_domain, 
                    subpatch_domain);
            
            // Evaluate patch
            Point3 face_map_ret[3];
            if(!_surface->_legacy){
                eval_legacy(xy.array(), EVAL_VALUE|EVAL_1ST_DERIV, face_map_ret);
            }else {
                FaceMapSurf::evaluate_bernstein_polynomial_patch(
                        _polynomial_patch,
                        _surface->_binomial_coefs_array,
                        _surface->_patch_order,
                        EVAL_VALUE|EVAL_1ST_DERIV, 
                        xy.array(),
                        face_map_ret);
            } 
            
            // Scale positions
            Point3 face_map_position = face_map_ret[0];
            face_map_position /= normalization_factor;
            Point3 face_map_n(cross(face_map_ret[1], face_map_ret[2]));
            assert(face_map_n.l2() > 0.);
            face_map_n /= face_map_n.l2();


            int eval_type = EVAL_VALUE;
           
            // blended surfaces can't evaluate derivatives at (x,y) =(0,0) or
            // (1,1) for some amazing reason; so we need to check if we're with
            // \eps_mach any of the corners.
            bool blendsurf_derivative_eval = is_blended_evaluator && 
                    (fabs(cd(0)) >= DBL_EPSILON && fabs(cd(1)) >= DBL_EPSILON) &&
                    (fabs(cd(0)) <= 1. -DBL_EPSILON && fabs(cd(1)) <= 1. -DBL_EPSILON);
            if(blendsurf_derivative_eval){
                eval_type = eval_type | EVAL_1ST_DERIV;
            }

            vector<Point3> position_and_derivs = 
                evaluator->evaluate(this->_blended_face_to_sample, 
                        eval_type, cd); 

            Point3* blended_ret = position_and_derivs.data();

            Point3 blended_position = blended_ret[0];
            blended_position /= normalization_factor;

            position_error = max(position_error, (face_map_position - blended_position).l2());
            
            //Obviously we can only compare values we've computed.
            if(blendsurf_derivative_eval){
                assert(position_and_derivs.size() == 3); 

                Point3 blended_n(cross(blended_ret[1], blended_ret[2]));
                blended_n /= blended_n.l2();

                normal_error = max(normal_error, (face_map_n - blended_n).l2());
            }
        }
    }
    return position_error < eps_pos && normal_error < eps_normal;

}
double FaceMapPatch::compute_jacobian(double* xy){
    Point3 ret[3];
    FaceMapSurf::evaluate_bernstein_polynomial_patch(
            _polynomial_patch,
            _surface->_binomial_coefs_array,
            _surface->_patch_order,
            EVAL_VALUE|EVAL_1ST_DERIV,
            xy,
            ret);
    Point3 n(cross(ret[1], ret[2]));
    return n.l2();
}

double FaceMapSurf::compute_jacobian(int F, double* xy){
    return _patches[F]->compute_jacobian(xy);
}

void FaceMapPatch::eval_legacy(double* xy, int flags, Point3* ret){
  NumVector uv_coords(2, false, xy);
  FaceMapPatch* patch = this;
  if(flags & EVAL_VALUE){
      NumVector value = (*patch)(uv_coords);

      for(int d =0; d < 3; d++) 
          ret[0](d) = value(d);
  }
  if(flags & EVAL_1ST_DERIV){
      NumVector dp_du= patch->_derivatives[0](uv_coords);
      NumVector dp_dv= patch->_derivatives[1](uv_coords);

      for(int d =0; d < 3; d++){
          ret[1](d) = dp_du(d);
          ret[2](d) = dp_dv(d);
      }
  }
  if(flags & EVAL_2ND_DERIV){
      NumVector d2p_du2  = patch->_derivatives[2](uv_coords);
      NumVector d2p_dudv = patch->_derivatives[3](uv_coords);
      NumVector d2p_dv2  = patch->_derivatives[4](uv_coords);

      for(int d =0; d < 3; d++){
          ret[3](d) = d2p_du2(d);
          ret[4](d) = d2p_dudv(d);
          ret[5](d) = d2p_dv2(d);
      }
  }

}

extern "C" void update_quad_validity(p4est_iter_volume_info_t * info, void *user_data){
    RefinementData<FaceMapPatch>* r = ((RefinementData<FaceMapPatch>*)info->quad->p.user_data);
    FaceMapPatch* patch = static_cast<FaceMapPatch*>(r->patch);
    assert(patch != NULL);
    double fit_error = patch->_surface->_fitting_error;
    r->refine = !(patch->is_patch_valid( static_cast<SurfaceEvaluator*>(user_data), fit_error, fit_error));
}


void FaceMapSurf::build_child_patches_from_parent(p4est_t* /* p4est */ ,
        p4est_topidx_t /*which_tree*/, 
        int num_outgoing,
        p4est_quadrant_t* outgoing[],
        int num_incoming,
        p4est_quadrant_t* incoming[]){

    assert(num_incoming == 4);
    assert(num_outgoing == 1);
    
    p4est_quadrant_t* parent = outgoing[0];
    
    RefinementData<FaceMapPatch>* parent_ref_data = 
        static_cast<RefinementData<FaceMapPatch>*>(parent->p.user_data);
    FaceMapPatch* patch = parent_ref_data->patch;
    
    vector<p4est_quadrant_t*> incoming_quads(incoming, incoming + num_incoming);
    vector<FaceMapPatch*> child_patches = patch->split(incoming_quads);

    for(int qi = 0; qi < num_incoming; qi++){
        p4est_quadrant_t* quad = incoming_quads[qi];
        quad->p.user_data = new RefinementData<FaceMapPatch>();
        RefinementData<FaceMapPatch>* child_ref_data = 
            static_cast<RefinementData<FaceMapPatch>*>(quad->p.user_data);
        
        child_ref_data->refine = parent_ref_data->refine;
        child_ref_data->patch= child_patches[qi];


    }
}


void FaceMapSurf::resolve_function_tree(SurfaceEvaluator* evaluator, 
        double eps_a, double eps_r){

    // mark all patches in the quadtree forest for refinement
    mark_all_patches_for_refinement<FaceMapPatch>(_p4est);
    bool patches_still_need_refinement = true;
    // while there are still quads that need refinement
    while(patches_still_need_refinement){
        // iterate over the leaves, check if each one is valid via callbacks
        update_refinement_criteria_on_quads(_p4est, evaluator, update_quad_validity);

        // if so, continue
        p4est_refine_ext(_p4est, 0, -1, is_quad_valid<FaceMapPatch>, NULL, build_child_patches_from_parent);

        // else, split it using split() inside of the refine_ext callback
        patches_still_need_refinement = check_refinement_criteria<FaceMapPatch>(_p4est);
    }
    // extract final set of patches of the leaves of the refined forest.
    vector<FaceMapPatch*> patches= collect_patches<FaceMapPatch>(_p4est);
    _patches = collect_patches<FaceMapPatch>(_p4est);

}

vector<FaceMapPatch*> FaceMapSurf::resolve_function(SurfaceEvaluator* evaluator, 
        queue<FaceMapPatch*>& invalid_patches, double eps_a, double eps_r){

    vector<FaceMapPatch*> valid_patches;
    
    while(!invalid_patches.empty() ){
        FaceMapPatch* invalid_patch = invalid_patches.front();
        invalid_patches.pop();

        // Split the patch into four new ones
        vector<FaceMapPatch*> subpatches = invalid_patch->split();

        // Check to see if the new patches are valid or not; if so, keep
        // them for later; if not, add to the back of the queue
        for(int i = 0; i < subpatches.size(); i++){
            FaceMapPatch* subpatch = subpatches[i];

            if(subpatch->is_patch_valid(evaluator, eps_a, eps_r)){
                valid_patches.push_back(subpatch);
            } else {
                invalid_patches.push(subpatch);
            }
        }

        cout << invalid_patches.size() << endl;
    }
    return valid_patches;

}

// Sample and evaluate the blended C^\inf surface
void FaceMapSurf::compute_xy_coordinates_and_blended_positions(
        //int F, // glocal face number F
        FaceMapPatch* patch,
        vector<Point2>& xy_coords, // xy coordinates
        vector<Point3>& limit_positions){ // surface limit positions coordinates

    // Choose the 1st vertex cooresponding to global face F as defined by GpMesh
    // we don't care which point this is so long it is uniquely determined;
    // the other corners of the face will be mapped to corners of unit square
    // accordingly
    
    // MJM TODO BUG CHECK no reference to vertex id; this could be a problem,
    // but maybe not, since we assume we're in the face-centered coordinates
    // below
    int ctrlsize = pow2(patch->_level_of_refinement);
    int num_samples_per_dim = 4*_patch_order;
    // compute number of samples
    //int ss = (int)ceil( _evalUpperBnd * (ctrlsize) ); assert(ss==ctrlsize); 

    /*
     * We subdivide patches in the following fashion:
     *
     *  ------------        --------------------
     *  |          |        | 12 | 13 | 14 | 15 | 
     *  |          |        --------------------- 
     *  |          |  --->  | 8  |  9 | 10 | 11 | 
     *  |          |        --------------------- 
     *  |          |        | 4  |  5 | 6  |  7 | 
     *  ------------        --------------------- 
     *                      | 0  |  1 |  2 |  3 | 
     *                      --------------------- 
     *                       
     *  As a result, we can determine the sampling interval upper and lower
     *  bounds for subpatch i (= F % num_blendsurf_faces) with refinment factor
     *  r by:
     *  x_lower = (1/2)^r * i mod 2^r
     *  x_upper = (1/2)^r * i mod 2^r + (1/2)^r
     *  y_lower = (1/2)^r * i / 2^r
     *  y_upper = (1/2)^r * i / 2^r + (1/2)^r
     *
     * where "/" corresponds to integer division.
     */
   
    double x_lower_bnd = patch->_x_int.first;
    double x_upper_bnd = patch->_x_int.second;
    double y_lower_bnd = patch->_y_int.first;
    double y_upper_bnd = patch->_y_int.second;

    int blended_id = patch->_blended_face_to_sample;
    // stepsize in parameter domain
    double ctrlstep = 1.0/float(num_samples_per_dim-1);

    for(int j = 0; j < num_samples_per_dim; j++){
        for(int i = 0; i < num_samples_per_dim; i++){
            // Equal sized steps in each direction over [0,1] x [0,1]
            double cd[2];
            cd[0] = i*ctrlstep/float(ctrlsize) + x_lower_bnd;
            cd[1] = j*ctrlstep/float(ctrlsize) + y_lower_bnd;
            
            
            vector<Point3> position_and_derivs = 
                _evaluator->evaluate(blended_id, EVAL_VALUE, Point2(cd)); 

            assert(position_and_derivs.size() == 1); 
            Point3 ret = position_and_derivs[0];
            
            // rescale points to [0,1]x[0,1] for bernstein eval
           cd[0] = 1./(x_upper_bnd - x_lower_bnd)* (cd[0] - x_lower_bnd);
           cd[1] = 1./(y_upper_bnd - y_lower_bnd)* (cd[1] - y_lower_bnd);

            // Save em
            xy_coords.push_back(Point2(cd));
            limit_positions.push_back(ret);
        }
    }
}
// Sample and evaluate the underlying surface we'd like to approximate
void FaceMapSurf::compute_xy_coordinates_and_surface_positions(
        //int F, // glocal face number F
        FaceMapPatch* patch,
        NumMatrix& xy_coords, // xy coordinates
        NumMatrix& surface_positions){ // surface positions coordinates

    // Choose the 1st vertex cooresponding to global face F as defined by GpMesh
    // we don't care which point this is so long it is uniquely determined;
    // the other corners of the face will be mapped to corners of unit square
    // accordingly
    
    // MJM TODO BUG CHECK no reference to vertex id; this could be a problem,
    // but maybe not, since we assume we're in the face-centered coordinates
    // below
    int ctrlsize = pow2(patch->_level_of_refinement);
    // TODO BUG PRONE NEEDS TO MATCH 2*_patch_order LISTED ABOVE FACTOR OUT INTO
    // AN OPTION!!!
    int num_samples_per_dim = 2*_patch_order;
    // compute number of samples
    //int ss = (int)ceil( _evalUpperBnd * (ctrlsize) ); assert(ss==ctrlsize); 

   
    double x_lower_bnd = patch->_x_int.first;
    double x_upper_bnd = patch->_x_int.second;
    double y_lower_bnd = patch->_y_int.first;
    double y_upper_bnd = patch->_y_int.second;

    int blended_id = patch->_blended_face_to_sample;
    // stepsize in parameter domain
    double ctrlstep = 1.0/float(num_samples_per_dim-1);

    for(int j = 0; j < num_samples_per_dim; j++){
        for(int i = 0; i < num_samples_per_dim; i++){
            // Equal sized steps in each direction over [0,1] x [0,1]
            double cd[2];
            cd[0] = i*ctrlstep/float(ctrlsize) + x_lower_bnd;
            cd[1] = j*ctrlstep/float(ctrlsize) + y_lower_bnd;
            
            vector<Point3> position_and_derivs = 
                _evaluator->evaluate(blended_id, EVAL_VALUE, Point2(cd)); 

            assert(position_and_derivs.size() == 1); 
            Point3 ret = position_and_derivs[0];
            
            // rescale points to [0,1]x[0,1] for bernstein eval
           cd[0] = 1./(x_upper_bnd - x_lower_bnd)* (cd[0] - x_lower_bnd);
           cd[1] = 1./(y_upper_bnd - y_lower_bnd)* (cd[1] - y_lower_bnd);

            // Save em
            for(int d =0; d < xy_coords.m(); d++)
                xy_coords(d,j*(num_samples_per_dim)+i) = cd[d];

            for(int d =0; d < surface_positions.m(); d++)
                surface_positions(d, j*(num_samples_per_dim)+i) = position_and_derivs[0](d);
        }
    }
}


Vector<Point3>* FaceMapSurf::construct_bernstein_polynomial_patch(
        //int F, // glocal face number F
        //FaceMapPatch* patch,
        vector<Point2>& xy_coords, // xy coordinates
        vector<Point3>& limit_positions){ // surface limit positions coordinates
    //NumMatrix normal_eqs(num_rows, num_cols);
    //clear(normal_eqs);
    //cout << "blended face fit: " << patch->_blended_face_to_sample << endl;
    int poly_degree = _patch_order;
    //int ctrlsize = pow2(_subdivCtrlLevel);  //assert(_subdivCtrlLevel==2);
    //int ctrlsize = pow2(_refinement_factor);
    //int ctrlsize = 4;
    int ctrlsize = int(sqrt(xy_coords.size()))-1;
    int stride = ctrlsize +1;
    
    int num_rows = xy_coords.size();
    int num_cols = (poly_degree+1)*(poly_degree+1);

    // Form the least-squares overdetermined system matrix A
    // Legacy: adapted from blendsurf vertex-centered patch fitting
    // system matrix A of size (num 2D samples) x (tensor poly basis size)
    // A_ij := coefficient of the jth tensor product polynomial evaluated
    // at the ith 2D sample point 
    NumMatrix A(num_rows, num_cols);
    clear(A);
    for(int i = 0; i < num_rows; i++){
        double xi = xy_coords[i](0);
        double yi = xy_coords[i](1);

        // TODO switch to std::vector
        //double x_poly[poly_degree+1];
        //double y_poly[poly_degree+1];
        NumVector x_poly(poly_degree+1);
        NumVector y_poly(poly_degree+1);
        //setvalue(x_poly,0.);
        //setvalue(y_poly,0.);

        init_bernstein_polynomial(poly_degree, xi, _binomial_coefs_array, x_poly);
        init_bernstein_polynomial(poly_degree, yi, _binomial_coefs_array, y_poly);

        // TODO refactor
        int j = 0;
        for(int di = 0; di <= poly_degree; di++){
            for(int dj = 0; dj <= poly_degree; dj++){
                A(i,j++) = x_poly(di) * y_poly(dj);
            }
        }
        assert(j == num_cols);
    }


    //NumMatrix pseudo_inverse(num_cols, num_rows); //???
    //pinv(normal_eqs, 1e-10, pseudo_inverse);
    //updatePolyControlPoint(F, pseudo_inverse, limit_positions);
    
    NumMatrix pseudo_inverse((poly_degree+1) * (poly_degree+1),xy_coords.size() ); 
    setvalue(pseudo_inverse, 0.);
    // -----------------------------------------------------------------------
    // Compute 1d least squares poly on each edge first and store coefficients
    // -----------------------------------------------------------------------
    
    // TODO remove unused index vectors 
    vector<double> edge_points_x_is_0;
    vector<double> edge_points_x_is_1;
    vector<double> edge_points_y_is_0;
    vector<double> edge_points_y_is_1;
    
    vector<Point3> limit_points_x_is_0;
    vector<Point3> limit_points_x_is_1;
    vector<Point3> limit_points_y_is_0;
    vector<Point3> limit_points_y_is_1;
    
    vector<Point2> interior_points;
    vector<Point3> interior_limit_points;

    // Construct index sets for boundary points
    vector<int> x_is_0_idx;
    vector<int> x_is_1_idx;
    vector<int> y_is_0_idx;
    vector<int> y_is_1_idx;

    // Construct index sets for the non-zero 1d poly's and points on each boundary
    vector<int> x_is_0_poly_idx;
    vector<int> x_is_1_poly_idx;
    vector<int> y_is_0_poly_idx;
    vector<int> y_is_1_poly_idx;
    set<int> boundary_set;
    set<int> boundary_poly_set;
    set<int> interior_set;
    set<int> interior_poly_set;

    for(int i = 0; i <= ctrlsize; i++){
        x_is_0_idx.push_back(i*stride);
        y_is_0_idx.push_back(i);
        x_is_1_idx.push_back(i*stride + ctrlsize);
        y_is_1_idx.push_back(stride*ctrlsize + i);
        boundary_set.insert(i*stride);
        boundary_set.insert(i);
        boundary_set.insert(i*stride + ctrlsize);
        boundary_set.insert(stride*ctrlsize + i);
    }
   
    int deg = poly_degree+ 1;
    for(int di = 0; di < deg; di++){
        x_is_0_poly_idx.push_back(poly_degree*deg + di);
        y_is_0_poly_idx.push_back(di*deg + poly_degree);
        x_is_1_poly_idx.push_back(di);
        y_is_1_poly_idx.push_back(di*deg);
        boundary_poly_set.insert(poly_degree*deg + di);
        boundary_poly_set.insert(di*deg + poly_degree);
        boundary_poly_set.insert(di);
        boundary_poly_set.insert(di*deg);
    }
    
    
    for(int di = 1; di < poly_degree; di++){
        for(int dj = 1; dj < poly_degree; dj++){
            interior_poly_set.insert(di*(poly_degree+1) +dj);
        }
    }
    for(int i = 1; i < ctrlsize; i++){
        for(int j = 1; j < ctrlsize; j++){
            interior_set.insert(i*stride + j);
        }
    }

    vector<int> boundary_idx(boundary_set.begin(), boundary_set.end());
    vector<int> boundary_poly_idx(boundary_poly_set.begin(), boundary_poly_set.end());
    vector<int> interior_idx(interior_set.begin(), interior_set.end());
    vector<int> interior_poly_idx(interior_poly_set.begin(), interior_poly_set.end());
    std::sort(interior_poly_idx.begin(), interior_poly_idx.end());
    std::sort(boundary_poly_idx.begin(), boundary_poly_idx.end());
   
    // split  A into A = ( A' | A_b ) where A_b is the columns  of A 
    // corresponding to the basis polynomials whose coefficients are completely
    // determined by the by the points on the boundary, and A' is the remaining
    // basis polynomials
    
    // we're solving A'*x' + A_b*x_b = b in the least squares sense.
    // MJM BUG TODO 3*deg is a consequence of use cubic poly's as the basis
    // degree. For high-order patches, revisit and refactor code
    NumMatrix A_b(xy_coords.size(), 4*poly_degree); // TODO BUG hard coded for cubic
    NumMatrix A_prime(xy_coords.size(), deg*deg - 4*poly_degree); // TODO BUG hard coded for cubic
    NumVector temp(xy_coords.size());
    
    for(int i = 0; i < boundary_poly_idx.size(); i++){
        A.getColumn(boundary_poly_idx[i], temp);
        A_b.setColumn(i, temp);
        
    }
    for(int i = 0; i < interior_poly_idx.size(); i++){
        A.getColumn(interior_poly_idx[i], temp);
        A_prime.setColumn(i, temp);
        
    }

    // Construct the least-squares system for the points on the boundary by
    // extracting out rows of A_b corresponding to the boundary points only
    // This produces a system of that characterizes the exact constraints of the
    // boundary points on each other
    NumMatrix boundary_system(4*ctrlsize, 4*poly_degree); // TODO check that number of rows is correct
    
    NumVector temp2(4*poly_degree);
    for(int i = 0; i < boundary_idx.size(); i++){
        A_b.getRow(boundary_idx[i], temp2);
        boundary_system.setRow(i, temp2);
    }

    // Now solve the boundary system in the least-squares sense
    NumMatrix boundary_system_pseudo_inv(4*poly_degree, 4*ctrlsize);
    pinv(boundary_system, 1e-14, boundary_system_pseudo_inv); // computing pseudo-inverse via truncated SVD


    // Compute hard boundary constraints given the limit positions there
    vector<Point3> x_b(4*poly_degree);
    for(int di = 0; di < 4*poly_degree; di++){ // TODO check correctness
        x_b[di] = Point3(0.0);
        for(int i = 0; i < 4*ctrlsize; i++){ // 4 because we're using quads
            // compute x_b = C^+*b
            x_b[di] += 
                boundary_system_pseudo_inv(di,i)*limit_positions[boundary_idx[i]];
        }
    }
    
    // Now we're computing b - A_b*x_b, i.e. the contribution of the boundary
    // fitting on the limit positions of all points on the patch
    // x_b was computed just above, A_b was formed earlier as a sample of
    // columns of A. This will remove from b the component contributed by the
    // boundary, which we solved for above
    vector<Point3> corrected_rhs(xy_coords.size());
    // Copy the original rhs into corrected_rhs
    for(int i = 0; i < xy_coords.size(); i++){
        corrected_rhs[i] = limit_positions[i];
    }
    
    // Subtract off the boundary contribution
    for(int i = 0; i < xy_coords.size(); i++){
        for(int di = 0; di < 4*poly_degree; di++){ // TODO check correctness 
            corrected_rhs[i] -=  A_b(i,di)*x_b[di];
        }
    }

    // Find the remain coefficients of the polynomial patch
    NumMatrix A_prime_pinv(A_prime.n(), A_prime.m());
    setvalue(A_prime_pinv, 0.);
    pinv(A_prime, 1e-10, A_prime_pinv);

    // Apply it to the remaining limit points
    //vector<Point3> x_i(deg); // hard coded for cubics
    vector<Point3> x_i(deg*deg - 4*poly_degree); // hard coded for cubics
    for(int di = 0; di < x_i.size(); di++){ // TODO check correctness
        x_i[di] = Point3(0.0);
        for(int i = 0; i < A_prime_pinv.n(); i++){ 
            x_i[di] += 
                A_prime_pinv(di,i)*corrected_rhs[i];
        }
    }

    // Store the final results
  Vector<Point3>* coefficients = new Vector<Point3>((poly_degree+1)*(poly_degree+1));
  for(int i = 0; i < (poly_degree+1)*(poly_degree+1); i++)
      (*coefficients)(i) = Point3(0.0);

  for(int i = 0; i < boundary_poly_idx.size(); i++)
      (*coefficients)(boundary_poly_idx[i]) = x_b[i];
  
  for(int i = 0; i < interior_poly_idx.size(); i++)
      (*coefficients)(interior_poly_idx[i]) = x_i[i];
   
  //_controls[F] = coefficients;
    return coefficients;

}

void FaceMapSurf::evaluate_bernstein_polynomial_patch(
        Vector<Point3>* _polynomial_patch,
        IntMatrix* _binomial_coefs,
        int poly_degree,
        int flags, 
        double* xy, Point3* ret){     
    double init=  omp_get_wtime();

    //Vector<Point3>* ivp = _controls[F];
    assert( _polynomial_patch!=NULL );
    double start = omp_get_wtime();
    Vector<Point3> vp = *_polynomial_patch;
    IntMatrix binomial_coefs = *_binomial_coefs;
    //cout << "copy time : " << start - omp_get_wtime() << endl;
    double x = xy[0];	 double y = xy[1];


    //int poly_degree = _surface->_patch_order;

    int cn = (poly_degree+1)*(poly_degree+1);
    assert(cn==vp.length());

    NumVector x_poly(poly_degree+1);
    NumVector y_poly(poly_degree+1);
    //setvalue(x_poly,0.);
    //setvalue(y_poly,0.);
    //double x_poly[poly_degree+1];
    //double y_poly[poly_degree+1];

    // Fill x_poly and y_poly with x_poly[i] = x^i, y_poly[i] = y^i for
    // i = 0, ..., poly_degree
    
    start = omp_get_wtime();
    init_bernstein_polynomial(poly_degree, x, _binomial_coefs, x_poly );
    init_bernstein_polynomial(poly_degree, y, _binomial_coefs, y_poly );
    //cout << "basis func eval time : " << start - omp_get_wtime() << endl;


    if(flags & EVAL_VALUE){
    start = omp_get_wtime();
        NumVector tensor_berstein_poly(cn);
        //setvalue(tensor_berstein_poly, 0.);
        for(int di=0; di < poly_degree+1; di++){
            for(int dj=0; dj < poly_degree+1; dj++){
                tensor_berstein_poly(di*(poly_degree+1) +dj) = x_poly(di)*y_poly(dj);
            }
        }
        Point3 position(0.,0.,0.);
        for(int i = 0; i < cn; i++){
            position += tensor_berstein_poly(i)*vp(i);
        }

        ret[0] = position;
    //cout << "func eval time : " << start - omp_get_wtime() << endl;
    }
    if(flags & EVAL_1ST_DERIV){
    start = omp_get_wtime();
        int n = poly_degree;

        NumVector x_prime(poly_degree+1);
        NumVector y_prime(poly_degree+1);
        //setvalue(x_prime, 0.);
        //setvalue(y_prime, 0.);

        x_prime(0) = double(n)*pow(x, n-1);
        y_prime(0) = double(n)*pow(y, n-1);

        x_prime(n) = -double(n)*pow(1.-x, n-1);
        y_prime(n) = -double(n)*pow(1.-y, n-1);

        for(int i = 1; i < n; i++){
            //x_prime(i) = double(n)*( _surface->B(i, n-1, x) - _surface->B(i-1, n-1, x) );
            //y_prime(i) = double(n)*( _surface->B(i, n-1, y) - _surface->B(i-1, n-1, y) );
            x_prime(i) = double(n)*( B(i, n-1, x, binomial_coefs(n-1,i)) - B(i-1, n-1, x, binomial_coefs(n-1,i-1)) );
            y_prime(i) = double(n)*( B(i, n-1, y, binomial_coefs(n-1,i)) - B(i-1, n-1, y, binomial_coefs(n-1,i-1)) );
        }


        NumVector tensor_berstein_derivs((poly_degree+1)*(poly_degree+1));
        for(int di=0; di < poly_degree+1; di++){
            for(int dj=0; dj < poly_degree+1; dj++){
                tensor_berstein_derivs(di*(poly_degree+1) +dj) = x_prime(di)*y_poly(dj);
            }
        }

        Point3 dfdx(0.);
        for(int i = 0; i < cn; i++){
            dfdx += tensor_berstein_derivs(i)*vp(i);
        }

        ret[1] = dfdx;

        //NumVector tensor_berstein_derivs(poly_degree*poly_degree);
        //setvalue(tensor_berstein_derivs, 0.);
        for(int di=0; di < poly_degree+1; di++){
            for(int dj=0; dj < poly_degree+1; dj++){
                tensor_berstein_derivs(di*(poly_degree+1) +dj) = x_poly(di)*y_prime(dj);
            }
        }
        Point3 dfdy(0.);
        for(int i = 0; i < cn; i++){
            dfdy += tensor_berstein_derivs(i)*vp(i);
        }
        ret[2] = dfdy;

    //cout << "1st deriveval time : " << start - omp_get_wtime() << endl;
    }
    if(flags & EVAL_2ND_DERIV){
    start = omp_get_wtime();
        int n = poly_degree;

        NumVector dx(poly_degree+1);
        NumVector dy(poly_degree+1);

        NumVector dx2(poly_degree+1);
        NumVector dy2(poly_degree+1);

        //setvalue(dx, 0.);
        //setvalue(dy, 0.);
        //setvalue(dx2, 0.);
        //setvalue(dy2, 0.);

        // --------------------------------------------------------------------
        dx(0) = double(n)*pow(x, n-1);
        dy(0) = double(n)*pow(y, n-1);

        dx(n) = -double(n)*pow(1.-x, n-1);
        dy(n) = -double(n)*pow(1.-y, n-1);

        for(int i = 1; i < n; i++){
            dx(i) = double(n)*( B(i, n-1, x, binomial_coefs(n-1,i)) - B(i-1, n-1, x, binomial_coefs(n-1,i-1)) );
            dy(i) = double(n)*( B(i, n-1, y, binomial_coefs(n-1,i)) - B(i-1, n-1, y, binomial_coefs(n-1,i-1)) );
        }
        // --------------------------------------------------------------------
        dx2(0) = double(n*(n-1))*pow(x, n-2);
        dy2(0) = double(n*(n-1))*pow(y, n-2);
        
        dx2(n) = double(n*(n-1))*pow(1.-x, n-2);
        dy2(n) = double(n*(n-1))*pow(1.-y, n-2);
       dx2(1)   = double((n-1)*(n-2))*pow(x, n-3) - double(n*(n-1))*pow(x, n-2);
       dy2(1)   = double((n-1)*(n-2))*pow(y, n-3) - double(n*(n-1))*pow(y, n-2);
       dx2(n-1) = double((n-1)*(n-2))*x*pow(1.-x, n-3)-2.*double(n-1)*pow(1.-x,n-2);
       dy2(n-1) = double((n-1)*(n-2))*y*pow(1.-y, n-3)-2.*double(n-1)*pow(1.-y,n-2);
      /* 
        dx2(1)   = double((n-1)*(n-2))*pow(x, n-3) - double(n*(n-1))*pow(x, n-2);
        dy2(1)   = double((n-1)*(n-2))*pow(y, n-3) - double(n*(n-1))*pow(y, n-2);
        
        dx2(n-1) = double((n-1)*(n-2))*pow(1.-x, n-3);
        dy2(n-1) = double((n-1)*(n-2))*pow(1.-y, n-3);
        */
        //dx2(1)   = (1.-x)*double((n-1)*(n-2))*pow(x,n-3) - 2.*double(n-1)*pow(x,n-2);
        //dy2(1)   = (1.-y)*double((n-1)*(n-2))*pow(y,n-3) - 2.*double(n-1)*pow(y,n-2);
        //dx2(1)   = -double(n-1)*(n*(x-1.) +2)*pow(x,n-3);
        //dy2(1)   = -double(n-1)*(n*(y-1.) +2)*pow(y,n-3);
        //dx2(1)   = -double(n*(n-1))*(n*x-n+2.)*pow(x, n-3);
        //dy2(1)   = -double(n*(n-1))*(n*y-n+2.)*pow(y, n-3);
        
        //dx2(n-1)   = -double(n-1)*(n*x -2.)*pow(1.-x,n-3);
        //dy2(n-1)   = -double(n-1)*(n*y -2.)*pow(1.-y,n-3);
        //dx2(n-1)   = -double(n*(n-1))*(n*x -2.)*pow(1.-x,n-3);
        //dy2(n-1)   = -double(n*(n-1))*(n*y -2.)*pow(1.-y,n-3);

        //dx2(n-1) = (x)*double((n-1)*(n-2))*pow(1.-x,n-3) - 2.*double(n-1)*pow(1.-x,n-2);
        //dy2(n-1) = (y)*double((n-1)*(n-2))*pow(1.-y,n-3) - 2.*double(n-1)*pow(1.-y,n-2);

        double n_choose_one = binomial_coefs(n,1);
        double n_choose_n_minus_one = binomial_coefs(n,n-1);
        dx2(1) *= n_choose_one;
        dy2(1) *= n_choose_one;
        
        dx2(n-1) *= n_choose_n_minus_one;
        dy2(n-1) *= n_choose_n_minus_one;

        for(int i = 2; i < n-1; i++){
            dx2(i) = double(n*(n-1))*( 
                        B(i,   n-2, x, binomial_coefs(n-2,i)  ) 
                    - 2*B(i-1, n-2, x, binomial_coefs(n-2,i-1)) 
                      + B(i-2, n-2, x, binomial_coefs(n-2,i-2)) 
                    );
            dy2(i) = double(n*(n-1))*( 
                        B(i,   n-2, y, binomial_coefs(n-2,i)  ) 
                    - 2*B(i-1, n-2, y, binomial_coefs(n-2,i-1)) 
                      + B(i-2, n-2, y, binomial_coefs(n-2,i-2)) 
                    );
        }

        NumVector tensor_berstein_derivs((poly_degree+1)*(poly_degree+1));
        for(int di=0; di < poly_degree+1; di++){
            for(int dj=0; dj < poly_degree+1; dj++){
                tensor_berstein_derivs(di*(poly_degree+1) +dj) = dx2(di)*y_poly(dj);
            }
        }

        Point3 d2f_dx2(0.);
        for(int i = 0; i < cn; i++){
            d2f_dx2 += tensor_berstein_derivs(i)*vp(i);
        }

        ret[3] = d2f_dx2;

        //NumVector tensor_berstein_derivs(poly_degree*poly_degree);
        //setvalue(tensor_berstein_derivs, 0.);
        for(int di=0; di < poly_degree+1; di++){
            for(int dj=0; dj < poly_degree+1; dj++){
                tensor_berstein_derivs(di*(poly_degree+1) +dj) = dx(di)*dy(dj);
            }
        }
        Point3 d2f_dydx(0.);
        for(int i = 0; i < cn; i++){
            d2f_dydx+= tensor_berstein_derivs(i)*vp(i);
        }
        ret[4] = d2f_dydx;
        //ret[5] = d2f_dydx;

        //NumVector tensor_berstein_derivs(poly_degree*poly_degree);
        //setvalue(tensor_berstein_derivs, 0.);
        for(int di=0; di < poly_degree+1; di++){
            for(int dj=0; dj < poly_degree+1; dj++){
                tensor_berstein_derivs(di*(poly_degree+1) +dj) = x_poly(di)*dy2(dj);
            }
        }
        Point3 d2f_dy2(0.);
        for(int i = 0; i < cn; i++){
            d2f_dy2+= tensor_berstein_derivs(i)*vp(i);
        }
        ret[5] = d2f_dy2;
    //cout << "2nd deriveval time : " << start - omp_get_wtime() << endl;
    }
//cout << "total eval time : " << init - omp_get_wtime() << endl;
}

IntMatrix* FaceMapSurf::compute_binomial_coefs_array(int _patch_order){
    // Compute binomial coefficients by dynamic programming with the relation
    // (n choose k) = (n-1 choose k) + (n-1 choose k-1)

    // Base cases: (n choose 0) and (n choose n) is 1 for any n;
    //             (n choose 1) = n for any n
    IntMatrix* binomial_coefs_array = new IntMatrix(_patch_order+1, _patch_order+1);

    (*binomial_coefs_array)(0,0) = 1;
   

    for(int n = 1; n <= _patch_order; n++){
        (*binomial_coefs_array)(n,0) = 1;
        (*binomial_coefs_array)(n,1) = n;
        (*binomial_coefs_array)(n,n) = 1;
    }

    // Can be further accelerated because (n choose k) = (n choose n-k)
    // I don't think we're doing high enough order for this to matter
    // practically, since this is a precomputation anyways
    for(int n = 2; n <= _patch_order; n++){
        for(int k=2; k < n; k++){
            int n_minus_one_choose_k = (*binomial_coefs_array)(n-1,k);
            int n_minus_one_choose_k_minus_one = (*binomial_coefs_array)(n-1,k-1);
            
            // (n choose k) = (n-1 choose k) + (n-1 choose k-1)
            (*binomial_coefs_array)(n,k) =
                n_minus_one_choose_k + n_minus_one_choose_k_minus_one;

        }
    }
    return binomial_coefs_array;
}

double pow_sq(double x, int n){
    assert(0); // this is broken. fix it
    if( n < 0){
        x = 1./x;
        n *= -1.;
    }
    if( n == 0 ){
        return 1;
    }
    double y = 1;
    while(n > 1){
        if(n % 2){
            x *= x;
            n /= 2;
        } else {
            y = x*y;
            x = x*x;
            n = (n-1)/2;
        }
    }
    return x*y; 
    
    /*
    if( n < 0){
        return pow_sq(1./x, -n);
    } else if( n == 0){
        return 1;
    } else if (n == 1){
        return x;
    } else if (n % 2){
        return pow_sq(x*x, n/2);
    } else if( n % 2 == 1){
        return pow_sq(x*x, (n-1)/2);
    }*/
}


double FaceMapSurf::B(int i, int n, double x, int binomial_coef){

    double p_x;
    if(i == 0){
        p_x= pow(x, n);

    } else if(i == n){
        p_x= pow(1. - x, n);

    } else {
        p_x= pow(x, n-i)*pow(1. - x, i);
    }

    // Multiply by the proper binomial coefficient pg. 79 of Curves and
    // Surface for Computer Aided Graphical Design by Gerald Farin
    return binomial_coef*p_x;
}


// ith n-degree Berstein polynomial evaluated at x
// B_i^n = x^{n-i}*(1-x)^i
double FaceMapSurf::B(int i, int n, double x, IntMatrix* binomial_coefs){

    double p_x;
    if(i == 0){
        p_x= pow(x, n);

    } else if(i == n){
        p_x= pow(1. - x, n);

    } else {
        p_x= pow(x, n-i)*pow(1. - x, i);
    }

    // Multiply by the proper binomial coefficient pg. 79 of Curves and
    // Surface for Computer Aided Graphical Design by Gerald Farin
    return (*binomial_coefs)(n,i)*p_x;
}

void FaceMapSurf::init_bernstein_polynomial(int degree, double x,
        IntMatrix* binomial_coefs,  NumVector& x_poly){
    //assert(degree <= 3);

    // evaluate x^(3-i) * (1-x)^i
    for(int i = 0; i <= degree; i++){
        // Evaluate ith Bernstein polynomial
        x_poly(i) = B(i, degree, x, (*binomial_coefs)(degree,i));
    }
}

int FaceMapSurf::eval(int flags, int F, double* xy, Point3* ret){
    // BLENDSURF
  assert( 
          (flags== EVAL_VALUE) || 
          (flags== (EVAL_VALUE|EVAL_1ST_DERIV)) ||
          (flags== (EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV))
        );

  double reprLB = 0;
  double reprUB = 1;
  double* xyi = xy;
  NumVector uv_coords(2, false, xy);
  FaceMapPatch* patch = _patches[F];
  if(!_legacy){
      patch->eval_legacy(xy, flags, ret);
  } else {
      // MJM TODO check  possibly redundant
      if(abs(xyi[0])<=reprUB && abs(xyi[1])<=reprUB &&
              abs(xyi[0])>=reprLB && abs(xyi[1])>=reprLB
        )  
          evaluate_bernstein_polynomial_patch(
                  _patches[F]->_polynomial_patch,
                  _binomial_coefs_array,
                  _patch_order,
                  flags, 
                  xy,
                  ret);
  }
  return 0;
}


// Needed for interface with the renderer...
int FaceMapSurf::evalall(int _lvl, int _gen){
  
  int TTL = pow2(_lvl);
  int NS = TTL-1;
  for(int V=0; V<1; V++) {
      //int K = _bdvalence(V);
    // BLENDSURF
      //int K = _bdsurf->gpmesh().valence(V);
      int K = _mesh.valence(V);
      if(K==0) continue;
      if (K==4 &&_gen == 2) continue;
      for(int f=0; f<K; f++) {
          double step = 1.0/double(TTL);      
          for(int j=0; j<=NS; j++)
              for(int i=0; i<=NS; i++) {
                  double cd[2];  
                  cd[0] = i*step; cd[1] = j*step;
                  //double xy[2];
                  Point3 ret[6];
                  // BLENDSURF
                  eval( EVAL_VALUE|EVAL_1ST_DERIV, V, cd, ret);

              }
      }
  }

  return 0;
}

vector<Point3> FaceMapSurf::control_points_on_subdomain(
        int F, 
        pair<double, double> x_interval,
        pair<double, double> y_interval){
    vector<Point3> control_point_copy;
    Vector<Point3> control_point_vec = control_points(F);
    
    for(int i=0; i < control_point_vec.length(); i++)
        control_point_copy.push_back(control_point_vec(i));

    vector<Point3> subdomain_control_points = 
        Bezier::compute_control_points_on_subdomain(
            control_point_copy,
            _patch_order,
                x_interval,
                y_interval);
    return control_point_copy;
}
void FaceMapSurf::bounding_box( int F, Point3& bounding_box_min, Point3& bounding_box_max){
    // Compute bounding boxes on [0,.5]^2, [.5,1]^2, [0,.5] x [.5, 1], and
    // [.5, 1] x [0,.5] via subdivision. 
    vector<Point3> subdomain_mins(4, Point3());
    vector<Point3> subdomain_maxes(4, Point3());
    
    vector<pair<double, double> > subdomains(2, pair<double, double>());
    subdomains[0].first = 0.;
    subdomains[0].second = .5;
    subdomains[1].first = .5;
    subdomains[1].second = 1.;
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            int idx = 2*i+j;
            auto x_int = subdomains[i];
            auto y_int = subdomains[j];
            bounding_box_on_subdomain(F, x_int, y_int, 
                    subdomain_mins[idx], subdomain_maxes[idx]);
            
        }
    }
    
    bounding_box_min = Point3(DBL_MAX, DBL_MAX, DBL_MAX);
    bounding_box_max = -Point3(DBL_MAX, DBL_MAX, DBL_MAX);

    // compute the bounding box of the bounding boxes
    for(const auto subdomain_min : subdomain_mins){ 
        bounding_box_min.x() = min(subdomain_min.x(), bounding_box_min.x());
        bounding_box_min.y() = min(subdomain_min.y(), bounding_box_min.y());
        bounding_box_min.z() = min(subdomain_min.z(), bounding_box_min.z());
    }
    for(const auto subdomain_max : subdomain_maxes){ 
        bounding_box_max.x() = max(subdomain_max.x(), bounding_box_max.x());
        bounding_box_max.y() = max(subdomain_max.y(), bounding_box_max.y());
        bounding_box_max.z() = max(subdomain_max.z(), bounding_box_max.z());
    }
        
}


void FaceMapSurf::bounding_box_on_subdomain(
        int F, 
        pair<double, double> x_interval,
        pair<double, double> y_interval,
        Point3& bounding_box_min,
        Point3& bounding_box_max){
    // TODO address this stupid vector copying
    vector<Point3> control_point_copy;
    Vector<Point3> control_point_vec = control_points(F);
    
    for(int i=0; i < control_point_vec.length(); i++)
        control_point_copy.push_back(control_point_vec(i));

    vector<Point3> subdomain_control_points = 
        Bezier::compute_control_points_on_subdomain(
            control_point_copy,
            _patch_order,
                x_interval,
                y_interval);

    bounding_box_min = Point3(DBL_MAX, DBL_MAX, DBL_MAX);
    bounding_box_max = -Point3(DBL_MAX, DBL_MAX, DBL_MAX);

    for(unsigned i = 0; i < subdomain_control_points.size(); i++){
        Point3 control_point = subdomain_control_points[i];
        bounding_box_min.x() = min(control_point.x(), bounding_box_min.x());
        bounding_box_min.y() = min(control_point.y(), bounding_box_min.y());
        bounding_box_min.z() = min(control_point.z(), bounding_box_min.z());
        
        bounding_box_max.x() = max(control_point.x(), bounding_box_max.x());
        bounding_box_max.y() = max(control_point.y(), bounding_box_max.y());
        bounding_box_max.z() = max(control_point.z(), bounding_box_max.z());
    }
    // Buffer by 10*machine epsilon to ensure all points are actually contained
    // inside (addressing roundoff in polynomial evaluation)
    bounding_box_max += Point3(2e-15);
    bounding_box_min -= Point3(2e-15);

}

vector<FaceMapPatch*> FaceMapSurf::cast_to_face_map_patches(p4est_t* p4est){
    vector<FaceMapPatch*> patches= collect_patches<FaceMapPatch>(p4est);
    return patches;
}


