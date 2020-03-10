
#include "catch.hpp"
#include "test_patchwork.hpp"
#include <face_map.hpp>
#include <diff_geom.hpp>
#include <analytic_evaluator.hpp>
#include <polynomial_evaluator.hpp>
#include <sampling.hpp>


TEST_CASE("Test face-map approximation of geometric quantities", "[diff-geom]"){
    char filec[300];
    int refinement_factor;
    int patch_order;
    bool  adaptive;
    char filepoly[300];

    load_options(filec, refinement_factor, patch_order, adaptive, filepoly);
    patch_order = 15;
    refinement_factor = 0;
    //AnalyticEvaluator* polynomial_evaluator = 
        //new AnalyticEvaluator( string("../wrl_files/flat_patch.wrl"), &gaussian);
    PolynomialEvaluator* polynomial_evaluator = 
       new PolynomialEvaluator( string("../wrl_files/poly/parabaloid.poly"), string("../wrl_files/parabaloid.wrl") );

    //FaceMapSurf* face_map = new FaceMapSurf(polynomial_evaluator);
    //face_map->set_params(refinement_factor, patch_order, 1); // adaptive
    //face_map->setup();

    FaceMapSurf* reference_face_map = new FaceMapSurf(polynomial_evaluator);
    reference_face_map->set_params(refinement_factor, patch_order, 0, 1e-4); // base polynomial for comparison
    reference_face_map->_legacy = true;
    reference_face_map->setup();
    assert(reference_face_map->_patches.size() == 1);
    SECTION("test determinant of first fundamental form"){
        int pi = 0;

        int num_samples = 50;
        NumMatrix samples = Sampling::sample_2d<Sampling::equispaced>(num_samples, Sampling::base_domain);
        for(int i = 0; i < samples.n(); i++){
            Point2 uv(samples.clmdata(i));
            vector<Point3> position_and_derivs(6,Point3());
            reference_face_map->eval(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV, pi, uv, position_and_derivs.data());

            Point3 P_u = position_and_derivs[1];
            Point3 P_v = position_and_derivs[2];
            Point3 P_uu = position_and_derivs[3];
            Point3 P_uv = position_and_derivs[4];
            Point3 P_vv = position_and_derivs[5];
            Point3 normal = cross(P_u, P_v);
            double normal_magnitude = normal.length();
            normal /= normal_magnitude;
            vector<double> coefs = Differential::first_fundamental_form(P_u, P_v);
            CHECK(fabs(normal_magnitude - Differential::determinant_of_first_fundamental_form(coefs)) <= 1e-14);
            vector<double> second_form = Differential::second_fundamental_form(P_u,P_v,P_uu, P_uv, P_vv);
            vector<double> first_form = Differential::first_fundamental_form(P_u,P_v);
            double E = first_form[0];
            double F = first_form[1];
            double G = first_form[2];

            double L = second_form[0];
            double M = second_form[1];
            double N = second_form[2];
            CHECK(fabs(L - dot(P_uu, normal)) <=1e-12); // L
            CHECK(fabs(M - dot(P_uv, normal)) <=1e-12); // L
            CHECK(fabs(N - dot(P_vv, normal)) <=1e-12); // L
            double gaussian_curvature = Differential::gaussian_curvature(uv, reference_face_map,pi);
            CHECK(fabs(gaussian_curvature - (L*N-M*M)/(E*G-F*F)) <=1e-10);
            double mean_curvature = Differential::mean_curvature(uv, reference_face_map,pi);

            Point3 du = cross(P_uu, P_v) + cross(P_u, P_uv); // d/du(P_u \cross P_v)
            Point3 dv = cross(P_uv, P_v) + cross(P_u, P_vv);// d/dv(P_u \cross P_v)
            Point3 n_u = (du - dot(normal,du)/normal_magnitude*normal)/normal_magnitude;
            Point3 n_v = (dv - dot(normal,dv)/normal_magnitude*normal)/normal_magnitude;
            //n_u/=normal_magnitude;
            //n_v/=normal_magnitude;
            double surface_div = Differential::surface_divergence(uv, n_u, n_v,reference_face_map, pi);
            CHECK(fabs(2*mean_curvature + surface_div) <=1e-10);
            //CHECK(fabs(surface_div) <=1e-14);


        }

    }
    /*
    SECTION("test determinant of first fundamental form"){
        //for(int pi =0; pi < reference_face_map->_patches.size(); pi++){
        int pi = 0;

        int num_samples = 50;
        double step_size = 1./double(num_samples);
        //FaceMapSurf::FaceMapPatch* patch = reference_face_map->_patches[pi];
        Function patch_du = reference_face_map->_patches[pi]->_derivatives[0];
        Function patch_dv = reference_face_map->_patches[pi]->_derivatives[1];
        for(int i = 0; i <= num_samples; i++){
            for(int j = 0; j <= num_samples; j++){

                Point2 uv(i*step_size, j*step_size);
                Point3 position_and_derivs[3];
                reference_face_map->eval(EVAL_VALUE|EVAL_1ST_DERIV, pi, uv, position_and_derivs);

                double normal_magnitude = cross(position_and_derivs[1], position_and_derivs[2]).length();
                //vector<double> coefs = Differential::first_fundamental_form(position_and_derivs[1], position_and_derivs[2]);
                vector<double> coefs = Differential::first_fundamental_form(patch_du, patch_dv, uv);
                CHECK(fabs(normal_magnitude - Differential::determinant_of_first_fundamental_form(coefs)) <= 1e-14);

            }
        }
    } SECTION("Test grad, div, laplace-beltrami on sphere with vector fields that produce zero"){

        // This test is sketchy; getting some excellent O(1) error at the
        // "patch" boundaries, but the patch is fit to a sphere, so it's kind of
        // a hack as it is. moreover, increasing polynomial order increases
        // error, which is quite excellent. Ask Denis/Abtin why.
        //
        // Anyway, low order approximations and skipping the boundaries is a
        // reasonable test for surface grad for now.

        int patch_order = 4;
        int num_samples = 2*patch_order;
        const int DIM = 3;
        double step_size = 1./double(num_samples-1);

        NumMatrix uv_coords(2, num_samples*num_samples);
        NumMatrix scalar_function_values(1, num_samples*num_samples);
        NumMatrix surface_values(3, num_samples*num_samples);
        NumMatrix div_free_values(3, num_samples*num_samples);
        for(int i = 0; i < num_samples; i++){
            for(int j = 0; j < num_samples; j++){
                int index = i*num_samples +j;
                Point2 uv(i*step_size, j*step_size);
                double u = uv.x();
                double v = uv.y();
                Point3 sphere_sample(
                        cos(u/(2.*M_PI))*sin(v/M_PI),
                        sin(u/(2.*M_PI))*sin(v/M_PI),
                        cos(v/M_PI));

                // partial derivatives will be divergence-free on-surface, so
                // that's a quick and dirty test
                Point3 sphere_du_exact(
                        -1./(2.*M_PI)*sin(u/(2.*M_PI))*sin(v/M_PI),
                        1./(2.*M_PI)*cos(u/(2.*M_PI))*sin(v/M_PI),
                        0.);

                double x = sphere_sample.x();
                double y = sphere_sample.y();
                double z = sphere_sample.z();

                CHECK(fabs (sphere_sample.length() -1.) <=1e-15);

                scalar_function_values(0,index) = .5*(x*x + y*y + z*z);
                for(int d = 0; d < DIM; d++){
                    surface_values(d,index) = sphere_sample(d);

                    // V(x,y,z) = (z, x, y)
                    div_free_values(d, index) = sphere_du_exact(d);
                }
                for(int d = 0; d < 2; d++)
                    uv_coords(d,index) = uv(d);
            }
        }

        Function scalar_func = Function::construct_function(
                uv_coords,
                scalar_function_values,
                Common::BEZIER,
                patch_order);
        Function sphere = Function::construct_function(
                uv_coords,
                surface_values,
                Common::BEZIER,
                patch_order);
        Function div_free = Function::construct_function(
                uv_coords,
                div_free_values,
                Common::BEZIER,
                patch_order);

        assert(scalar_func.domain_dim() == 2);
        assert(sphere.domain_dim() == 2);
        assert(div_free.domain_dim() == 2);
        assert(scalar_func.range_dim() == 1);
        assert(sphere.range_dim() == 3);
        assert(div_free.range_dim() == 3);

        Function sphere_du = Function::d_dx(0, sphere);
        Function sphere_dv = Function::d_dx(1, sphere);
        Function sphere_dudu = Function::d_dx(0, sphere_du);
        Function sphere_dudv = Function::d_dx(1, sphere_du);
        Function sphere_dvdv = Function::d_dx(1, sphere_dv);

        Function scalar_func_du = Function::d_dx(0, scalar_func);
        Function scalar_func_dv = Function::d_dx(1, scalar_func);
        Function scalar_func_dudu = Function::d_dx(0, scalar_func_du);
        Function scalar_func_dudv = Function::d_dx(1, scalar_func_du);
        Function scalar_func_dvdv = Function::d_dx(1, scalar_func_dv);

        Function div_free_du = Function::d_dx(0, div_free);
        Function div_free_dv = Function::d_dx(1, div_free);


        vector<Function> scalar_derivs(5);
        vector<Function> sphere_derivs(5);
        scalar_derivs[0] = scalar_func_du;
        scalar_derivs[1] = scalar_func_dv;
        scalar_derivs[2] = scalar_func_dudu;
        scalar_derivs[3] = scalar_func_dudv;
        scalar_derivs[4] = scalar_func_dvdv;

        sphere_derivs[0] = sphere_du;
        sphere_derivs[1] = sphere_dv;
        sphere_derivs[2] = sphere_dudu;
        sphere_derivs[3] = sphere_dudv;
        sphere_derivs[4] = sphere_dvdv;

        num_samples = 4*patch_order;
        step_size = 1./double(num_samples-1);

        for(int i = 1; i < num_samples-1; i++){
            for(int j = 1; j < num_samples-1; j++){
                int index = i*num_samples +j;
                Point2 uv(i*step_size, j*step_size);
                Point3 gradient = Differential::surface_gradient(
                        sphere_du, sphere_dv, 
                        scalar_func_du, scalar_func_dv,
                        uv);
                CHECK((gradient).length() <= 1e-11);

                double divergence = Differential::surface_divergence(
                        sphere_du, sphere_dv, 
                        div_free_du, div_free_dv, 
                        uv);

                // the vector field is exactly along surface normal and is
                CHECK(fabs(divergence) <= 1e-13);
                double surface_laplacian= Differential::laplace_beltrami(
                        sphere_derivs,
                        scalar_derivs,
                        uv);
                CHECK(fabs(surface_laplacian) <= 5e-8);

            }
        }


    }
    SECTION("test analytic solutions"){

    }*/
                    /*
    SECTION("test laplace-beltrami surface operator"){
        for(int pi =0; pi < reference_face_map->_patches.size(); pi++){
            Function patch_du = reference_face_map->_patches[pi]->_derivatives[0];
            Function patch_dv = reference_face_map->_patches[pi]->_derivatives[1];

            int num_samples = 8;
            double step_size = 1./double(num_samples);
            for(int i = 0; i <= num_samples; i++){
                for(int j = 0; j <= num_samples; j++){

                    Point2 uv(i*step_size, j*step_size);
                    Point3 position_and_derivs[6];
                    Point3 laplace_beltrami = 
                        Differential::laplace_beltrami(
                                    reference_face_map->_patches[pi]->_derivatives,
                                uv);
                    reference_face_map->eval(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV, 
                            pi, uv, position_and_derivs);

                    Point3 n = cross(position_and_derivs[1], position_and_derivs[2]).dir();
                    double H = Differential::mean_curvature(uv, reference_face_map, pi);
                    Point3 laplace_beltrami = 
                        Differential::laplace_beltrami(uv, 
                                position_and_derivs[1],
                                position_and_derivs[2],
                                position_and_derivs[3],
                                position_and_derivs[4],
                                position_and_derivs[5],
                                reference_face_map, pi);
                    Point3 true_laplace_beltrami = -2*H*n;
                    cout << laplace_beltrami << ", " << true_laplace_beltrami << endl;
                    CHECK((laplace_beltrami - true_laplace_beltrami).length()/true_laplace_beltrami.length() <= 1e-8);

                }
            }
        }
    }

    SECTION(""){

        int num_samples = 3;
        double step_size = 1./double(num_samples);

        map<pair<int, double>, double> gaussian_curvature_error;
        map<pair<int, double>, double> mean_curvature_error;
        for(int patch_order = 4; patch_order <=20; patch_order+=4){
            for(int digits = 3; digits<=7; digits++){
                double rel_error = pow(10, -digits);
                double mean_curvature_inf_norm = DBL_MIN;
                double gaussian_curvature_inf_norm = DBL_MIN;

                double mean_curvature_abs_error= DBL_MIN;
                double gaussian_curvature_abs_error= DBL_MIN;


                FaceMapSurf* face_map = new FaceMapSurf(polynomial_evaluator);
                face_map->set_params(refinement_factor, patch_order, 1); // adaptive
                face_map->_fitting_error = rel_error;
                face_map->setup();
                for(int pi =0; pi < face_map->_patches.size(); pi++){
                    Interval x_int = face_map->_patches[pi]->_x_int;
                    Interval y_int = face_map->_patches[pi]->_y_int;
                    for(int i = 0; i <= num_samples; i++){
                        for(int j = 0; j <= num_samples; j++){
                            //cout << mean_curvature_inf_norm << ", " << gaussian_curvature_inf_norm << ", " << mean_curvature_abs_error << ", " << gaussian_curvature_abs_error << endl;
                            // sample the face-map patch
                            Point2 subpatch_xy(i*step_size, j*step_size);

                            // sample the original surface we're approximating
                            Point2 parent_patch_xy(
                                    (x_int.second - x_int.first)*(subpatch_xy.x()) + x_int.first, 
                                    (y_int.second - y_int.first)*(subpatch_xy.y()) + y_int.first);

                            //cout << subpatch_xy << ", " << parent_patch_xy << endl;
                            // evaluate geometry quantities on face-map and "true"
                            // surface

                            // face-map
                            double face_map_gaussian_curvature =
                                Differential::gaussian_curvature(subpatch_xy, face_map, pi);
                            double face_map_mean_curvature =
                                Differential::mean_curvature(subpatch_xy, face_map, pi);

                            // true surface
                            double true_mean_curvature =
                                Differential::mean_curvature(parent_patch_xy, reference_face_map, 0);
                            double true_gaussian_curvature =
                                Differential::gaussian_curvature(parent_patch_xy, reference_face_map, 0);
                            //cout << "gaussian: " << face_map_gaussian_curvature << ", " << true_gaussian_curvature << endl;
                            //cout << "mean: " << face_map_mean_curvature << ", " << true_mean_curvature << endl << endl;


                            // compare error
                            // gaussian
                            gaussian_curvature_abs_error = 
                                max(gaussian_curvature_abs_error, 
                                        fabs(face_map_gaussian_curvature - true_gaussian_curvature ));
                            // mean
                            mean_curvature_abs_error = 
                                max(mean_curvature_abs_error, 
                                        fabs(face_map_mean_curvature - true_mean_curvature ));

                            // also keep track of the infinity norm of the true
                            // curvature
                            // gaussian
                            gaussian_curvature_inf_norm= 
                                max(gaussian_curvature_inf_norm, 
                                        fabs(true_gaussian_curvature ));
                            // mean
                            mean_curvature_inf_norm = 
                                max(mean_curvature_inf_norm, 
                                        fabs(true_mean_curvature ));
                        }
                    }
                }
                delete face_map;
            cout << "gaussian abs error: " << gaussian_curvature_abs_error << endl;
            cout << "true gaussian inf norm: " << gaussian_curvature_inf_norm << endl;
            cout << "gaussian curvature relative error: " << gaussian_curvature_abs_error/gaussian_curvature_inf_norm << endl << endl;
            cout << "mean abs error: " << mean_curvature_abs_error << endl;;
            cout << "true mean inf norm: " << mean_curvature_inf_norm << endl;
            cout << "mean curvature relative error: " << mean_curvature_abs_error/mean_curvature_inf_norm << endl;
                gaussian_curvature_error[pair<int,double>(patch_order, rel_error)] = gaussian_curvature_abs_error/gaussian_curvature_inf_norm;
                mean_curvature_error[pair<int,double>(patch_order, rel_error)] = mean_curvature_abs_error/mean_curvature_inf_norm;
            }
        }
            cout << "gaussian curvature error" << endl;
            for(map<pair<int,double>, double>::iterator it = gaussian_curvature_error.begin();
                    it != gaussian_curvature_error.end();
                    it++){
                pair<int, double> key = it->first;
                cout << "patch_order: "  << key.first << ", rel_error: " << key.second << "; gaussian_curvature_error: " << it->second <<endl ;
            }
            cout << "mean curvature error" << endl;
            for(map<pair<int,double>, double>::iterator it = mean_curvature_error.begin();
                    it != mean_curvature_error.end();
                    it++){
                pair<int, double> key = it->first;
                cout << "patch_order: "  << key.first << ", rel_error: " << key.second << "; mean_curvature_error: " << it->second<<endl;
            }
    }
    SECTION(""){
        double mean_curvature_inf_norm = DBL_MIN;
        double gaussian_curvature_inf_norm = DBL_MIN;

        double mean_curvature_abs_error= DBL_MIN;
        double gaussian_curvature_abs_error= DBL_MIN;

        int num_samples = 3;
        double step_size = 1./double(num_samples);
        for(int pi =0; pi < face_map->_patches.size(); pi++){
            Interval x_int = face_map->_patches[pi]->_x_int;
            Interval y_int = face_map->_patches[pi]->_y_int;
            for(int i = 0; i <= num_samples; i++){
                for(int j = 0; j <= num_samples; j++){
                    //cout << mean_curvature_inf_norm << ", " << gaussian_curvature_inf_norm << ", " << mean_curvature_abs_error << ", " << gaussian_curvature_abs_error << endl;
                    // sample the face-map patch
                    Point2 subpatch_xy(i*step_size, j*step_size);

                    // sample the original surface we're approximating
                    Point2 parent_patch_xy(
                            (x_int.second - x_int.first)*(subpatch_xy.x()) + x_int.first, 
                            (y_int.second - y_int.first)*(subpatch_xy.y()) + y_int.first);

                    //cout << subpatch_xy << ", " << parent_patch_xy << endl;
                    // evaluate geometry quantities on face-map and "true"
                    // surface

                    // face-map
                    double face_map_gaussian_curvature =
                        Differential::gaussian_curvature(subpatch_xy, face_map, pi);
                    double face_map_mean_curvature =
                        Differential::mean_curvature(subpatch_xy, face_map, pi);

                    // true surface
                    double true_mean_curvature =
                        Differential::mean_curvature(parent_patch_xy, reference_face_map, 0);
                    double true_gaussian_curvature =
                        Differential::gaussian_curvature(parent_patch_xy, reference_face_map, 0);
                    //cout << "gaussian: " << face_map_gaussian_curvature << ", " << true_gaussian_curvature << endl;
                    //cout << "mean: " << face_map_mean_curvature << ", " << true_mean_curvature << endl << endl;


                    // compare error
                    // gaussian
                    gaussian_curvature_abs_error = 
                        max(gaussian_curvature_abs_error, 
                                fabs(face_map_gaussian_curvature - true_gaussian_curvature ));
                    // mean
                    mean_curvature_abs_error = 
                        max(mean_curvature_abs_error, 
                                fabs(face_map_mean_curvature - true_mean_curvature ));
                    
                    // also keep track of the infinity norm of the true
                    // curvature
                    // gaussian
                    gaussian_curvature_inf_norm= 
                        max(gaussian_curvature_inf_norm, 
                                fabs(true_gaussian_curvature ));
                    // mean
                    mean_curvature_inf_norm = 
                        max(mean_curvature_inf_norm, 
                                fabs(true_mean_curvature ));

                }
            }
        }
        cout << "gaussian abs error: " << gaussian_curvature_abs_error << endl;
        cout << "true gaussian inf norm: " << gaussian_curvature_inf_norm << endl;
        cout << "gaussian curvature relative error: " << gaussian_curvature_abs_error/gaussian_curvature_inf_norm << endl << endl;
        cout << "mean abs error: " << mean_curvature_abs_error << endl;;
        cout << "true mean inf norm: " << mean_curvature_inf_norm << endl;
        cout << "mean curvature relative error: " << mean_curvature_abs_error/mean_curvature_inf_norm << endl;

    }*/
}
