#include "catch.hpp"
#include "test_patchwork.hpp"
#include <face_map.hpp>
#include <polynomial_evaluator.hpp>
#include <analytic_evaluator.hpp>
#include <blended_evaluator.hpp>
#include <sampling.hpp>

using Sampling::base_domain;
using Sampling::equispaced;
using Sampling::sample_2d;
using Sampling::affine_coordinate_transform;
struct FaceMapParameters{
    int refinement_factor;
    int patch_order;
    bool  adaptive;
    double fit_accuracy;
};

unique_ptr<FaceMapSurf> setup_face_map(const unique_ptr<SurfaceEvaluator>& evaluator, 
        FaceMapParameters parameters){

    auto face_map = unique_ptr<FaceMapSurf>(new FaceMapSurf(evaluator.get()));

    //face_map->_legacy = true;
    face_map->set_params(parameters.refinement_factor, 
            parameters.patch_order, 
            parameters.adaptive, 
            parameters.fit_accuracy);
    face_map->setup();
    return face_map;
}

void check_face_map_fit_accuracy(const unique_ptr<SurfaceEvaluator>& evaluator,
        const unique_ptr<FaceMapSurf>& face_map, double abs_eps){
    int num_samples = min(face_map->_patch_order*3,30);
    NumMatrix test_samples = sample_2d<equispaced>(num_samples, base_domain);
    for(const auto& patch : face_map->_patches){
        // subdomain of evaluator patch that patch is approximating
        Rectangle patch_subdomain(patch->_x_int, patch->_y_int); 
        for(int i =0; i < test_samples.n(); i++){
            // xy coordinates of patch
            Point2 xy_patch(test_samples.clmdata(i));

            // xy coordinates of evaluator 
            Point2 xy_evaluator = affine_coordinate_transform(
                    xy_patch,
                    base_domain, 
                    patch_subdomain);

            // evaluate face-map patch
            Point3 face_map_position[3];
            if(!face_map->_legacy){
                patch->eval_legacy(xy_patch.array(), EVAL_VALUE|EVAL_1ST_DERIV, face_map_position);
            } else {
            FaceMapSurf::evaluate_bernstein_polynomial_patch(
                    patch->_polynomial_patch,
                    patch->_surface->_binomial_coefs_array,
                    patch->_surface->_patch_order,
                    EVAL_VALUE|EVAL_1ST_DERIV, 
                    xy_patch.array(),
                    face_map_position);
            }
            // evaluate evalutor
            vector<Point3> evaluator_position= 
                evaluator->evaluate(patch->_blended_face_to_sample, 
                        EVAL_VALUE, xy_evaluator); 
            // make sure they match
            CHECK((evaluator_position[0] - face_map_position[0]).length() <=abs_eps);
        }
    }
}

TEST_CASE("", "[explicit][face-map][integration]"){
    SECTION("exact flat patch"){

        string poly_file = "wrl_meshes/poly/flat_patch.poly";
        string patch_file = "wrl_meshes/wrl/flat_patch.wrl";
        unique_ptr<SurfaceEvaluator> polynomial_evaluator = 
            unique_ptr<PolynomialEvaluator>(new PolynomialEvaluator(poly_file, patch_file));

        FaceMapParameters p;
        p.refinement_factor = 0;
        p.adaptive = false; // patch is explicit
        p.patch_order = 3;  // cubic patches
        p.fit_accuracy = 1.;
        auto face_map = setup_face_map(polynomial_evaluator, p);
        check_face_map_fit_accuracy(polynomial_evaluator, face_map, 1e-10);
    }
    SECTION("exact torus"){

        string poly_file = "wrl_meshes/poly/explicit_torus_patches.poly";
        string patch_file = "wrl_meshes/wrl/newtorus.wrl";
        unique_ptr<SurfaceEvaluator> polynomial_evaluator = 
            unique_ptr<PolynomialEvaluator>(new PolynomialEvaluator(poly_file, patch_file));

        FaceMapParameters p;
        p.refinement_factor = 0;
        p.adaptive = false; // patch is explicit
        p.patch_order = 3;  // cubic patches
        p.fit_accuracy = 1.;
        auto face_map = setup_face_map(polynomial_evaluator, p);
        check_face_map_fit_accuracy(polynomial_evaluator, face_map, 1e-14);
    }
}
TEST_CASE("", "[analytic][face-map][integration]"){
    SECTION("single gaussian patch"){
        string patch_file = "wrl_meshes/wrl/flat_patch.wrl";
        unique_ptr<SurfaceEvaluator> polynomial_evaluator = 
            unique_ptr<AnalyticEvaluator>(new AnalyticEvaluator(patch_file, &gaussian));

        FaceMapParameters p;
        p.refinement_factor = 0;
        p.adaptive = true; 
        p.patch_order = 10;  
        p.fit_accuracy = 1e-7;
        auto face_map = setup_face_map(polynomial_evaluator, p);
        check_face_map_fit_accuracy(polynomial_evaluator, face_map, p.fit_accuracy*10);
    }
}
TEST_CASE("", "[blended][face-map][integration][cube]"){
    string patch_file = "wrl_meshes/wrl/cube.wrl";

    unique_ptr<SurfaceEvaluator> polynomial_evaluator = 
        unique_ptr<BlendedEvaluator>(new BlendedEvaluator(patch_file));

    FaceMapParameters p;
    p.refinement_factor = 0;
    SECTION("blended cube, no adaptivity"){
        p.adaptive = false; 
        p.patch_order = 30;  
        p.fit_accuracy = 1e-7;
        auto face_map = setup_face_map(polynomial_evaluator, p);
        check_face_map_fit_accuracy(polynomial_evaluator, face_map, p.fit_accuracy*10);
    }
    SECTION("blended cube, with adaptivity"){
        p.adaptive = true; // patch is explicit
        p.patch_order = 10;  // cubic patches
        p.fit_accuracy = 1e-3;
        auto face_map = setup_face_map(polynomial_evaluator, p);
        check_face_map_fit_accuracy(polynomial_evaluator, face_map, 1e-4);
    }
}
TEST_CASE("", "[blended][face-map][integration][prop]"){
    string patch_file = "wrl_meshes/wrl/ppp.wrl";

    unique_ptr<SurfaceEvaluator> polynomial_evaluator = 
        unique_ptr<BlendedEvaluator>(new BlendedEvaluator(patch_file));

    FaceMapParameters p;
    p.refinement_factor = 0;
    SECTION("blended propeller, no adaptivity"){
        p.adaptive = false; 
        p.patch_order = 30;  
        p.fit_accuracy = 1e-5;
        auto face_map = setup_face_map(polynomial_evaluator, p);
        check_face_map_fit_accuracy(polynomial_evaluator, face_map, p.fit_accuracy*10);
    }
    SECTION("blended propeller, with adaptivity"){
        p.adaptive = true; // patch is explicit
        p.patch_order = 6;  // cubic patches
        p.fit_accuracy = 1e-3;
        auto face_map = setup_face_map(polynomial_evaluator, p);
        check_face_map_fit_accuracy(polynomial_evaluator, face_map, 1e-4);
    }
}
