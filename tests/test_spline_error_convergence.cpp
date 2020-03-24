#include "catch.hpp"
#include "test_patchwork.hpp"
#include "analytic_evaluator.hpp"
#include "face_map.hpp"
TEST_CASE("Test spline face-map patch convergence under uniform refinement", "[converge][refine]"){
    string filename = "wrl_meshes/wrl/flat_patch.wrl";
    AnalyticEvaluator* evaluator = new AnalyticEvaluator(filename, &exp_sin_cos);

    int upper = 24;
    //int upper = 22;
    //int upper = 10;
    int lower = 6;
    int step = 4;
    int iteration_count = (upper - lower)/step;
    int count = 0;
    NumMatrix relative_error(iteration_count+1, 6+1);
    for(int patch_order = lower; patch_order <= upper; patch_order += step){
        FaceMapSurf* face_map = new FaceMapSurf(evaluator);
        //int patch_order = 20;
        int refinement_factor = 0;
        int adaptive = 0;
        face_map->set_params(refinement_factor, patch_order, adaptive, 1e-4);
        face_map->setup();
        assert(face_map->_patches.size() == 1);
        //FaceMapSurf::FaceMapPatch* patch = face_map->_patches[0];

        int num_samples = 2*(patch_order+1);
        double step = 1./(num_samples-1);

        vector<double> abs_error_inf(6, -DBL_MAX);
        vector<double> inf_norm(6, -DBL_MAX);
        for(int i = 0; i < num_samples; i++){
            for(int j = 0; j < num_samples; j++){
                Point2 uv(i*step, j*step);
                uv *= .5;
                uv += .25;
                vector<Point3> patch_values(6, Point3());

                vector<Point3> surface_values = 
                    evaluator->evaluate(0, EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV, uv);

                face_map->eval( EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV, 0, uv.array(), patch_values.data());
                for(int d =0 ; d < 6; d++){
                    //cout << "d: " << d << ", " << surface_values[d] << " =?= " << patch_values[d] << endl;
                    double ith_abs_error = (surface_values[d] - patch_values[d]).linfty();
                    double ith_abs = (surface_values[d]).linfty();
                    abs_error_inf[d] = max(ith_abs_error ,abs_error_inf[d]);
                    inf_norm[d] = max(ith_abs , inf_norm[d]);
                }

            }
        }
        vector<double> rel_error(6,0.);
            relative_error(count, 0) = patch_order;
        for(int d =0 ; d < 6; d++){

            rel_error[d] = inf_norm[d] > 1e-13 ? abs_error_inf[d]/inf_norm[d] : 1e-13;
            relative_error(count, d+1) = rel_error[d];
            //relative_error(count, d+1) = abs_error_inf[d] -inf_norm[d]*1e-15;
            //relative_error(count, d+1) = abs_error_inf[d];
        }
        count++;

        delete face_map;
    }
//    ofstream f;
//    f.open("tests/bspline_relative_error_exp_sin_cos.csv");
//        cout << relative_error << endl;
//        f << relative_error;
//
//        f.close();


}
