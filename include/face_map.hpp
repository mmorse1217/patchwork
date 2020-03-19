#ifndef __FACE_MAP_HPP__
#define __FACE_MAP_HPP__

#include "bdsurf.hpp"
#include "surface_evaluator.hpp"
#include "vecmatop.hpp"
#include "function.hpp"
#include <map>
#include "common.hpp"
#include "p4est_interface.hpp"

enum { 
    EVAL_VALUE=1,  
    EVAL_1ST_DERIV=2,    
    EVAL_2ND_DERIV=4  
};

class FaceMapSurf;

class FaceMapPatch: public Function {
    public:
        FaceMapPatch(int F, int level, Interval x_int, Interval y_int, FaceMapSurf* surface);
        FaceMapPatch(int F, int level, Interval x_int, Interval y_int, 
                Common::BasisType basis_type, FaceMapSurf* surface);
        ~FaceMapPatch(){
            delete _polynomial_patch;
        }
        int _group_id;
        int _orientation;
        int _level_of_refinement;
        int _blended_face_to_sample;
        FaceMapSurf* _surface;
        int _global_refined_patch_id;
        Interval _x_int;
        Interval _y_int;
        vector<Function> _derivatives;
        p4est_quadrant_t* _quad;
        int _patch_order;

        // TODO REFACTOR use Function
        // these should be matrices of size 2 x (_patch_order+1)^2 and 3 x
        // (_patch_order+1)^2
        vector<Point3> _samples;
        vector<Point2> _xy_coords;
        Vector<Point3>* _polynomial_patch;

        // New FaceMapPatch Interface
        NumMatrix _samples_mat;
        NumMatrix _xy_coords_mat;
        NumMatrix _polynomial_patch_mat;
        // TODO REFACTOR use Function
        bool static is_patch_valid(
                SurfaceEvaluator* evaluator,
                int patch_order,
                double eps_pos, double eps_normal);
        bool is_patch_valid(SurfaceEvaluator* evaluator,
                double eps_pos, double eps_normal);
        // TODO REFACTOR use Function
        vector<FaceMapPatch*> split();
        vector<FaceMapPatch*> split(vector<p4est_quadrant_t*> child_quads);

        void eval_legacy(double* xy, int flags, Point3* ret);
        // TODO REFACTOR use Function d/dx
        double compute_jacobian(double* xy);

};



class FaceMapSurf{
    private:

        static void build_child_patches_from_parent(p4est_t* p4est,
                p4est_topidx_t which_tree, 
                int num_outgoing,
                p4est_quadrant_t* outgoing[],
                int num_incoming,
                p4est_quadrant_t* incoming[]);
    public:
        bool _legacy;
        double _fitting_error;
        int _refinement_factor;
        
        SurfaceEvaluator* _evaluator;
        p4est_t* _p4est;
        vector<int> _group_ids;
        
        // TODO REFACTOR maintain legacy interface for renderer
        vector<vector<Point3> >* _samples;
        vector<FaceMapPatch*> _patches;


        int _patch_order;
        bool _adaptive;
        bool _use_functional_rep;
        Common::BasisType _basis_type;

        NumMatrix _pseudo_inverse;
        // TODO REFACTOR remove pointer... this is insane
        IntMatrix* _binomial_coefs_array;

        // New interface
        std::map<Index, int> _binomial_coefs_stack;
        GpMesh _mesh;
            

        FaceMapSurf(SurfaceEvaluator* evaluator);
        ~FaceMapSurf(){
            if(_evaluator != NULL){
                //delete _evaluator;
            }
            if(_samples != NULL){
                delete _samples;
            }
            for(unsigned i = 0; i < _patches.size(); i++){
                // MJM TODO MEMORY LEAK make sure old patches are deleted when
                // doing refinement!!!
                delete _patches[i];
            }
            if(_p4est != NULL){
                //p4est_destroy(_p4est);
            }
        }

        vector<FaceMapPatch*> resolve_function(SurfaceEvaluator* evaluator, 
                queue<FaceMapPatch*>& invalid_patches, 
                double eps_a=1e-4, 
                double eps_r=1e-4);
        
        void resolve_function_tree(SurfaceEvaluator* evaluator, 
                double eps_a=1e-4, 
                double eps_r=1e-4);


        void setup();
        
        // TODO REFACTOR should just pass along to a call to Function()
        int eval(int flags, int V, double* xy, Point3* ret);
        // TODO maintain for render interface
        int evalall(int _lvl, int _gen);

        // TODO REFACTOR change to matrices for Function interface
        void compute_xy_coordinates_and_blended_positions(
                FaceMapPatch* patch,
                vector<Point2>& xy_coords, // xy coordinates
                vector<Point3>& limit_positions); // surface limit positions coordinates

        void compute_xy_coordinates_and_surface_positions(
                FaceMapPatch* patch,
                NumMatrix& xy_coords, // 2 x(_patch_order+1)^2 xy coordinates
                NumMatrix& surface_positions); // 3 x(_patch_order+1)^2 surface positions coordinates
        // Using Bernstein polynomials as basis for least-squares
        // TODO REFACTOR factor out C^0 enforcing and make it a switch to flip in options
        Vector<Point3>* construct_bernstein_polynomial_patch(
                vector<Point2>& xy_coords, 
                vector<Point3>& limit_positions); 

        NumMatrix construct_bernstein_polynomial_patch(
                NumMatrix& xy_coords, 
                NumMatrix& limit_positions); 

        
        // TODO REFACTOR to use basis functions directly instead of evaluating each time
        static void init_bernstein_polynomial(int degree, double x, 
                map<Index, int>* binomial_coefs, NumVector& x_poly);
        static void init_bernstein_polynomial(int degree, double x, 
                IntMatrix* binomial_coefs, NumVector& x_poly);


        // TODO remove depreacted
        static void evaluate_bernstein_polynomial_patch(
                Vector<Point3>* _polynomial_patch,
                IntMatrix* _binomial_coefs,
                int poly_degree,
                int flags, 
                double* xy, Point3* ret);

        int get_group_id(int F){
            return _patches[F]->_group_id;
        }
        vector<int> get_group_ids(){
            return _group_ids;
        }
        int get_orientation(int F){
            return _patches[F]->_orientation;
        }

        // Passing along calls to bdsurf variables for the sake of rendering
      void set_params( int refinement_factor, int patch_order,
            bool adaptive, double fitting_accuracy=1e-4){
          _patch_order = patch_order;
          _refinement_factor = refinement_factor;
          _adaptive = adaptive;
          _fitting_error = fitting_accuracy;
      }
      
      static IntMatrix* compute_binomial_coefs_array(int patch_order);
      
      // ith n-degree Berstein polynomial evaluated at x
      // B_i^n = x^{n-i}*(1-x)^i
      static double B(int i, int n, double x, int binomial_coefs);
      static double B(int i, int n, double x, IntMatrix* binomial_coefs);
        
        Point3 ctr() {
            return gpmesh().ctr();
        }

        void bounding_box(int F, Point3& min, Point3& max);

        GpMesh& gpmesh(){
            return _mesh;
        }
        SurfaceEvaluator* evaluator(){
            return _evaluator;
        }
        
        void bbox(Point3& bbmin, Point3& bbmax){
            gpmesh().bbox(bbmin, bbmax);
        }
        int numVertices(){
            return gpmesh().numVertices();
        }
        int numFaces(){
            return _patches.size();
        }
        int valence(int V){
            return gpmesh().valence(V);

        }
        
        vector<Point3> samples(int F){
            return (*_samples)[F];
        }
        Vector<Point3> control_points(int F){
            return (*_patches[F]->_polynomial_patch);
        }
        Vector<Point3>* control_points_write(int F){
            return &(*_patches[F]->_polynomial_patch);
        }

        pair<int, int> Fv2Vf(int F, int v){
            return gpmesh().Fv2Vf(F,v);
        }
        
        pair<int, int> Vf2Fv(int V, int f){
            return gpmesh().Vf2Fv(V,f);
        }
        double compute_jacobian(int F, double* xy);
        
        void bounding_box_on_subdomain(
                int F, 
                pair<double, double> x_interval,
                pair<double, double> y_interval,
                Point3& min,
                Point3& max);
        vector<Point3> control_points_on_subdomain(
                int F, 
                pair<double, double> x_interval,
                pair<double, double> y_interval);
        static vector<FaceMapPatch*> cast_to_face_map_patches(p4est_t* p4est);
};

#endif
