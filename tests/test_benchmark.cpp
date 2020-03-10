#include "catch.hpp"
#include "test_patchwork.hpp"
#include <face_map.hpp>
#include <polynomial_evaluator.hpp>
#include <omp.h>

void benchmark_face_map_evaluation(FaceMapSurf* face_map, vector<Point2> samples);


TEST_CASE("Benchmark FaceMap evaluations", "[bench]"){

    char filec[300];
    int refinement_factor;
    int patch_order;
    bool  adaptive;
    char filepoly[300];

    load_options(filec, refinement_factor, patch_order, adaptive, filepoly);
    refinement_factor = 0; 
    adaptive = 0;
    patch_order = 11;

    PolynomialEvaluator* polynomial_evaluator = 
        new PolynomialEvaluator(string("../wrl_files/poly/flat_patch.poly"), string("../wrl_files/flat_patch.wrl"));

    FaceMapSurf* face_map = new FaceMapSurf(polynomial_evaluator);
    face_map->set_params(refinement_factor, patch_order, adaptive, 1e-4);
    face_map->setup();
    
    FaceMapSurf* new_face_map = new FaceMapSurf(polynomial_evaluator);
    new_face_map->set_params(refinement_factor, patch_order, adaptive, 1e-4);
    new_face_map->_legacy = false;
    new_face_map->setup();

    // Generate equispaced samples in [0,1]^2 to evaluate the bezier patch at 
    int num_samples = 10;
    double step = 1./double(num_samples-1);
    vector<Point2> samples(num_samples*num_samples);
    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            int index = i*num_samples + j;

            Point2 uv(i*step, j*step);
            samples[index] = uv;
        }
    }
    SECTION("Benchmark legacy face-map evaluation vs. new Function evaluation"){
        cout << "Legacy FaceMap: " << endl;
        benchmark_face_map_evaluation(face_map, samples);
        //cout << "Function-based FaceMap: " << endl;
        //benchmark_face_map_evaluation(new_face_map, samples);
    }
    
}

void benchmark_face_map_evaluation(FaceMapSurf* face_map, vector<Point2> samples){
    int N = 1e3;
    int num_samples = samples.size();
    double eval_time = omp_get_wtime(); 
    for(int i = 0; i < N; i++){
        for(int si = 0; si < samples.size(); si++){
            Point2 uv = samples[si]; 
            vector<Point3> results(6,Point3(0.));
            face_map->eval(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV,
                    0, uv.array(), results.data());
        }
    }
    eval_time = omp_get_wtime() - eval_time;
    cout << "number of iterations: " << N*num_samples<< endl;
    cout << "total evaluation time: " << eval_time << endl;
    cout << "mean evaluation time: " << eval_time/double(N*num_samples) << endl;

}
