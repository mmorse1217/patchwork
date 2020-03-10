#include "catch.hpp"
#include "test_patchwork.hpp"
#include <face_map.hpp>
#include <polynomial_evaluator.hpp>
#include <blended_evaluator.hpp>
TEST_CASE("Test bezier related functions", "[bezier]"){
    char filec[300];
    int refinement_factor;
    int patch_order;
    bool  adaptive;
    char filepoly[300];

    load_options(filec, refinement_factor, patch_order, adaptive, filepoly);
    PolynomialEvaluator* polynomial_evaluator = 
        new PolynomialEvaluator(string("../wrl_files/poly/flat_patch.poly"), string("../wrl_files/flat_patch.wrl"));

    FaceMapSurf* face_map = new FaceMapSurf(polynomial_evaluator);
    face_map->set_params(refinement_factor, patch_order, adaptive, 1e-4);
    face_map->setup();

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
    
    vector<Point3> control_points( (patch_order+1)*(patch_order+1));
    Vector<Point3> control_points_vec = face_map->control_points(0);

    for(int i= 0 ; i <  (patch_order+1)*(patch_order+1); i++)
            control_points[i] = control_points_vec(i);


    SECTION("Test 1d de Casteljau"){
        int poly_degree = 3;
        Vector<Point3> control_points_vec(poly_degree+1);

        control_points_vec(0) = Point3(-1.,-1.,-1.);
        control_points_vec(1) = Point3(-2./3.,-2./3.,-2./3.);
        control_points_vec(2) = Point3(0.,0., 0.);
        control_points_vec(3) = Point3(0.,0., 0.);

        vector<Point3> control_points(poly_degree+1);

        for(int i= 0 ; i <  poly_degree+1; i++)
            control_points[i] = control_points_vec(i);


        for(int si=0; si < num_samples; si++){
            double t = samples[si].y();

            cout << "start1" << endl;
            Point3 blossom_position = Bezier::de_casteljau_recursive(control_points, t);
            cout << "start2" << endl;
            Point3 curve_position(0.,0.,0.);                    

            for(int i = 0; i < poly_degree+1; i++){
                curve_position += 
                    control_points[i]*FaceMapSurf::B(i, poly_degree, t, face_map->_binomial_coefs_array);
            }
            CHECK( (blossom_position - curve_position).length() <1e-12);

        }
    }
    SECTION("Test 1d de Casteljau backed by iterative blossom"){
        int poly_degree = 3;
        Vector<Point3> control_points_vec(poly_degree+1);

        control_points_vec(0) = Point3(-1.,-1.,-1.);
        control_points_vec(1) = Point3(-2./3.,-2./3.,-2./3.);
        control_points_vec(2) = Point3(0.,0., 0.);
        control_points_vec(3) = Point3(0.,0., 0.);

        vector<Point3> control_points(poly_degree+1);

        for(int i= 0 ; i <  poly_degree+1; i++)
            control_points[i] = control_points_vec(i);


        for(int si=0; si < num_samples; si++){
            double t = samples[si].y();

            cout << "start" << endl;
            Point3 blossom_position = Bezier::de_casteljau(control_points, t);
            cout << "done" << endl;
            Point3 curve_position(0.,0.,0.);                    

            for(int i = 0; i < poly_degree+1; i++){
                curve_position += 
                    control_points[i]*FaceMapSurf::B(i, poly_degree, t, face_map->_binomial_coefs_array);
            }
            CHECK( (blossom_position - curve_position).length() <1e-12);

        }
    }

    SECTION("Test tensor product blossom evaluation"){
        // Blossoming with all equal t-values should produce the function
        // value on the patch
        for(unsigned si=0; si < samples.size(); si++){
            int patch_order = face_map->_patch_order;
            Point2 sample = samples[si];

            // Evaluate by blossiming (all t-values are equal)
            vector<Point2> t_values(patch_order, sample);

            // change to Vector
            //vector<Point3> control_points( (patch_order+1)*(patch_order+1));
            //Vector<Point3> control_points_vec = face_map->control_points(0);

            //for(int i= 0 ; i <  (patch_order+1)*(patch_order+1); i++)
                //control_points[i] = control_points_vec(i);

            Point3 blossoming_position = 
                Bezier::tensor_product_blossom(
                        control_points,
                        t_values);


            // Evaluate the polynomial directly
            Point3 position;

            FaceMapSurf::evaluate_bernstein_polynomial_patch(
                    &control_points_vec,
                    face_map->_binomial_coefs_array,
                    patch_order,
                    EVAL_VALUE,
                    sample.array(), 
                    &position);

            // Should be the same
            CHECK( ( blossoming_position - position).length() <= 1e-12);
        }


    }
    
    SECTION("Test tensor product blossom to recover control points"){
        // blossoming can be used to recover the ith Bezier control point if t_i =
        // 1, and t_j = 0 for j \not= i. Here we check that we can actually do
        // this reliably.

        for(int di = 0; di <= patch_order; di++){
            for(int dj = 0; dj <= patch_order; dj++){
                // In 1d, a blossom will compute the ith control point in of a
                // bezier curve when evaluated using n-i 0's and i 1's as
                // t-values (order of values doesn't matter due to symmetry). 
                // Here, we enumerate the tuples (0, ..., 0), (1, 0, ..., 0),
                // etc. for both x and y values 
                vector<Point2> t_values(patch_order, Point2(0.));
                for(int i = 0; i < patch_order-di; i++){
                    t_values[i].x() = 1.;
                }
                for(int j = 0; j < patch_order-dj; j++){
                    t_values[j].y() = 1.;
                }
                
                // This should recover the (i,j)th control point
                Point3 ijth_control_point = control_points[di*(patch_order+1)+dj];
                
                // Actually do the blossom
                Point3 blossomed_value = 
                    Bezier::tensor_product_blossom_recursive(
                            control_points,
                            t_values);
                // Make sure it's correct
                CHECK( (ijth_control_point - blossomed_value).length() <=1e-15);
            }
        }
    }
    SECTION("Test iterative tensor product blossom to recover control points"){
        // blossoming can be used to recover the ith Bezier control point if t_i =
        // 1, and t_j = 0 for j \not= i. Here we check that we can actually do
        // this reliably.

        for(int di = 0; di <= patch_order; di++){
            for(int dj = 0; dj <= patch_order; dj++){
                // In 1d, a blossom will compute the ith control point in of a
                // bezier curve when evaluated using n-i 0's and i 1's as
                // t-values (order of values doesn't matter due to symmetry). 
                // Here, we enumerate the tuples (0, ..., 0), (1, 0, ..., 0),
                // etc. for both x and y values 
                vector<Point2> t_values(patch_order, Point2(0.));
                for(int i = 0; i < patch_order-di; i++){
                    t_values[i].x() = 1.;
                }
                for(int j = 0; j < patch_order-dj; j++){
                    t_values[j].y() = 1.;
                }
                
                // This should recover the (i,j)th control point
                Point3 ijth_control_point = control_points[di*(patch_order+1)+dj];
                
                // Actually do the blossom
                Point3 blossomed_value = 
                    Bezier::tensor_product_blossom(
                            control_points,
                            t_values);
                // Make sure it's correct
                CHECK( (ijth_control_point - blossomed_value).length() <=1e-15);
            }
        }
    }
    SECTION("Test tensor product blossom is a multilinear function"){
        // If one blossoms using the values (a*r+b*s, t_2, ..., t_n), producing
        // the value B[a*r+b*s, t_2, ..., t_n], this is equal to 
        // a*B[r, t_2, ..., t_n] + b*B[s, t_2, ..., t_n]. This is called
        // multilinearity. We'll test this with random a,b,r,s, and t_i values.
        

        for(int iter= 0; iter < 10; iter++){
            double a = drand48();
            double b = 1.-a; // enforce that a+b=1
            
            // Generate random values on [0,1]^2 to evaluate the blossom at

            vector<Point2> t_values(patch_order, Point2(0.));
            for(int i =1; i < patch_order; i++){
                t_values[i] = Point2(drand48(), drand48());
            }
            // Save the chosen (but random) values of r and s to test
            Point2 r(drand48(), drand48()); 
            Point2 s(drand48(), drand48()); 

            // First compute B[ar+bs, t_1, ...]
            t_values[0] = a*r+b*s;

            Point3 ar_plus_bs_blossom = 
                Bezier::tensor_product_blossom(
                        control_points, 
                        t_values);

            // Then compute B[r, ...] and B[s, ...]
            t_values[0] = r;
            Point3 r_blossom = 
                Bezier::tensor_product_blossom(
                        control_points, 
                        t_values);

            t_values[0] = s;
            Point3 s_blossom = 
                Bezier::tensor_product_blossom(
                        control_points, 
                        t_values);
            
            // Check that they're the same 
            CHECK((ar_plus_bs_blossom - (a*r_blossom + b*s_blossom)).length() <=1e-14);
        }
    }

    SECTION("Test tensor product blossom is symmetric"){
        // Test that B[t_1, ..., t_n] = B[t_2, t_1, ..., t_n]
        // We'll just test that this holds for regenerating control points i.e.
        // random permutations of a shuffle of [0, ..., 0, 1, ..., 1] for both x
        // and y coordinates of the "t"-values. 

        for(int di = 0; di <= patch_order; di++){
            for(int dj = 0; dj <= patch_order; dj++){
                for(int iter = 0; iter < 10; iter++){
                    vector<Point2> t_values(patch_order, Point2(0.));
                    for(int i = 0; i < patch_order-di; i++){
                        t_values[i].x() = 1.;
                    }
                    random_shuffle(t_values.begin(), t_values.end());
                    for(int j = 0; j < patch_order-dj; j++){
                        t_values[j].y() = 1.;
                    }
                    random_shuffle(t_values.begin(), t_values.end());

                    // This should recover the (i,j)th control point with each
                    // permuatation
                    Point3 ijth_control_point = control_points[di*(patch_order+1)+dj];

                    // Actually do the blossom on the permuted values
                    Point3 blossomed_value = 
                        Bezier::tensor_product_blossom(
                                control_points,
                                t_values);
                    // Make sure it's correct
                    CHECK( (ijth_control_point - blossomed_value).length() <=1e-15);
                }
            }
        }
    }
    SECTION("Test blossoming to compute control points on [a,b]^2 "){
        // First, compute the blossom B[a,..., a, b, ..., b] (n-i a's and i b's)
        // to find the control points for the bezier patch on [a,b]^2 \subset
        // [0,1]^2. Then test that patch evaluation using the original control 
        // points and the newly computed control points are equal.

        // Arbitrarily look at the subdomain [1/3, .45]^2
        double a = 1./3;
        double b = .45;
        int num_samples = 10;
        vector<Point3> subdomain_control_points((patch_order+1)*(patch_order+1));
        for(int di = 0; di <= patch_order; di++){
            for(int dj = 0; dj <= patch_order; dj++){
                // set everything to (a, a)
                vector<Point2> t_values(patch_order, Point2(a));

                // fill in the first i values with b (implies the remaining n-i
                // values are a's) Recall order doesn't matter by symmetry
                for(int i = 0; i < patch_order-di; i++){
                    t_values[i].x() = b;
                }
                for(int i = 0; i < patch_order-dj; i++){
                    t_values[i].y() = b;
                }

                //Evaluate the blossom and save the control point
                Point3 control_point_on_subdomain = 
                    Bezier::tensor_product_blossom(
                            control_points, 
                            t_values);
                
                //int index = (patch_order-di)*(patch_order+1)+(patch_order-dj);
                int index = di*(patch_order+1)+dj;
                subdomain_control_points[index] = control_point_on_subdomain;

            }
        }
        // generate samples on [a,b]^2 to evaluate the original bezier patch
        // and samples on [0,1]^2 to evaluate the new bezier patch
        double step = (b-a)/double(num_samples-1);
        vector<Point2> old_samples(num_samples*num_samples);
        vector<Point2> new_samples(num_samples*num_samples);
        for(int i = 0; i < num_samples; i++){
            for(int j = 0; j < num_samples; j++){
                int index = i*num_samples + j;

                Point2 uv_new(i/double(num_samples-1), j/double(num_samples-1));
                Point2 uv_old = (b-a)*uv_new + Point2(a);
                
                old_samples[index] = uv_old;
                new_samples[index] = uv_new;
            }
        }
        
        // for each sample..
        for(unsigned i =0; i < num_samples*num_samples; i++){

            vector<Point2> t_values_old(patch_order, old_samples[i]);
            // First evaluate the bezier patch using the original control
            // points
            Point3 true_blossomed_position = 
                Bezier::tensor_product_blossom(
                        control_points,
                        t_values_old);

            vector<Point2> t_values_new(patch_order, new_samples[i]);
            // then with the newly computed control points
            Point3 computed_blossomed_position = 
                Bezier::tensor_product_blossom(
                        subdomain_control_points,
                        t_values_new);

            // and make sure the values are the same
            CHECK((computed_blossomed_position - true_blossomed_position).length() <=1e-14);

        }

    }
    
    SECTION("Test bounding box on subdomain"){
        Point3 bounding_box_min;
        Point3 bounding_box_max;
        int F = 0;

        double a = .1;
        double b = .62;
        int num_samples = 10;
        vector<Point3> subdomain_control_points((patch_order+1)*(patch_order+1));
        for(int di = 0; di <= patch_order; di++){
            for(int dj = 0; dj <= patch_order; dj++){
                // set everything to (a, a)
                vector<Point2> t_values(patch_order, Point2(a));

                // fill in the first i values with b (implies the remaining n-i
                // values are a's) Recall order doesn't matter by symmetry
                for(int i = 0; i < patch_order-di; i++){
                    t_values[i].x() = b;
                }
                for(int i = 0; i < patch_order-dj; i++){
                    t_values[i].y() = b;
                }

                //Evaluate the blossom and save the control point
                Point3 control_point_on_subdomain = 
                    Bezier::tensor_product_blossom(
                            control_points, 
                            t_values);
                
                int index = di*(patch_order+1)+dj;
                subdomain_control_points[index] = control_point_on_subdomain;

            }
        }

        // Compute the bounding box
        Point3 min, max;
        face_map->bounding_box_on_subdomain(F, Interval(a,b), Interval(a,b), 
                min, max);
        // Check that each sample point on the subdomain [a,b]^2 is contained in
        // the computed bounding box
        double step = (b-a)/double(num_samples-1);
        vector<Point2> samples(num_samples*num_samples);
        for(int i = 0; i < num_samples; i++){
            for(int j = 0; j < num_samples; j++){
                int index = i*num_samples + j;

                Point2 uv(step*i + a, step*j+ a);
                Point3 position; 
                //samples[index] = uv;
                FaceMapSurf::evaluate_bernstein_polynomial_patch(
                        &control_points_vec,
                        face_map->_binomial_coefs_array,
                        patch_order,
                        EVAL_VALUE,
                        uv.array(), 
                        &position);
                CHECK(position.x() > min.x()); 
                CHECK(position.x() < max.x());
                CHECK(position.y() > min.y()); 
                CHECK(position.y() < max.y());
                CHECK(position.z() > min.z()); 
                CHECK(position.z() < max.z());
            }
        }

    }

}


TEST_CASE("Test face-map basis funcs.", "[basis]"){

    int poly_degree = 6;
    double x = .5;
    double y = .5;
    IntMatrix* _binomial_coefs = FaceMapSurf::compute_binomial_coefs_array(poly_degree);
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
        dx(i) = double(n)*( FaceMapSurf::B(i, n-1, x, (*_binomial_coefs)(n-1,i)) - FaceMapSurf::B(i-1, n-1, x, (*_binomial_coefs)(n-1,i-1)) );
        dy(i) = double(n)*( FaceMapSurf::B(i, n-1, y, (*_binomial_coefs)(n-1,i)) - FaceMapSurf::B(i-1, n-1, y, (*_binomial_coefs)(n-1,i-1)) );
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

       //dx2(n-1) = double((n-1)*(n-2))*pow(1.-x, n-3);
       //dy2(n-1) = double((n-1)*(n-2))*pow(1.-y, n-3);
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

    double n_choose_one = (*_binomial_coefs)(n,1);
    double n_choose_n_minus_one = (*_binomial_coefs)(n,n-1);
    dx2(1) *= n_choose_one;
    dy2(1) *= n_choose_one;

    dx2(n-1) *= n_choose_n_minus_one;
    dy2(n-1) *= n_choose_n_minus_one;

    for(int i = 2; i < n-1; i++){
        dx2(i) = double(n*(n-1))*( 
                FaceMapSurf::B(i,   n-2, x, (*_binomial_coefs)(n-2,i)  ) 
                - 2*FaceMapSurf::B(i-1, n-2, x, (*_binomial_coefs)(n-2,i-1)) 
                + FaceMapSurf::B(i-2, n-2, x, (*_binomial_coefs)(n-2,i-2)) 
                );
        dy2(i) = double(n*(n-1))*( 
                FaceMapSurf::B(i,   n-2, y, (*_binomial_coefs)(n-2,i)  ) 
                - 2*FaceMapSurf::B(i-1, n-2, y, (*_binomial_coefs)(n-2,i-1)) 
                + FaceMapSurf::B(i-2, n-2, y, (*_binomial_coefs)(n-2,i-2)) 
                );
    }
    cout << "d/dx" << endl;
    cout << dx << endl;
    cout << "d/dy" << endl;
    cout << dy << endl;
    cout << "d2/dx2" << endl;
    cout << dx2 << endl;
    cout << "d2/dy2" << endl;
    cout << dy2 << endl;


}

/*TEST_CASE("Test bezier related functions on blended surface", "[bezier]"){
    char filec[300];
    int refinement_factor;
    int patch_order;
    bool  adaptive;
    char filepoly[300];








    load_options(filec, refinement_factor, patch_order, adaptive, filepoly);
    BlendedEvaluator* blended_evaluator = 
        new BlendedEvaluator("../wrl_files/cube.wrl");

    FaceMapSurf* face_map = new FaceMapSurf(blended_evaluator);
    adaptive = 1;
    patch_order = 6;
    face_map->set_params(refinement_factor, patch_order, adaptive);
    face_map->setup();

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
    
    SECTION("Test bounding box on subdomain"){
        Point3 bounding_box_min;
        Point3 bounding_box_max;
        int F = 0;
        for(int pi = 0; pi < face_map->_patches.size(); pi++){
            vector<Point3> control_points( (patch_order+1)*(patch_order+1));
            Vector<Point3> control_points_vec = face_map->control_points(pi);

            for(int i= 0 ; i <  (patch_order+1)*(patch_order+1); i++)
                control_points[i] = control_points_vec(i);

            double a = .62;
            double b = .62000001;
            int num_samples = 10;
            vector<Point3> subdomain_control_points((patch_order+1)*(patch_order+1));
              for(int di = 0; di <= patch_order; di++){
              for(int dj = 0; dj <= patch_order; dj++){
            // set everything to (a, a)
            vector<Point2> t_values(patch_order, Point2(a));

            // fill in the first i values with b (implies the remaining n-i
            // values are a's) Recall order doesn't matter by symmetry
            for(int i = 0; i < patch_order-di; i++){
            t_values[i].x() = b;
            }
            for(int i = 0; i < patch_order-dj; i++){
            t_values[i].y() = b;
            }

            //Evaluate the blossom and save the control point
            Point3 control_point_on_subdomain = 
            Bezier::tensor_product_blossom(
            control_points, 
            t_values);

            int index = di*(patch_order+1)+dj;
            subdomain_control_points[index] = control_point_on_subdomain;

            }
            }

            // Compute the bounding box
            Point3 min, max;
            face_map->bounding_box_on_subdomain(pi, Interval(a,b), Interval(a,b), 
                    min, max);
            // Check that each sample point on the subdomain [a,b]^2 is contained in
            // the computed bounding box
            double step = (b-a)/double(num_samples-1);
            vector<Point2> samples(num_samples*num_samples);
            for(int i = 0; i < num_samples; i++){
                for(int j = 0; j < num_samples; j++){
                    int index = i*num_samples + j;

                    Point2 uv(step*i + a, step*j+ a);
                    Point3 position; 
                    //samples[index] = uv;
                    FaceMapSurf::evaluate_bernstein_polynomial_patch(
                            &control_points_vec,
                            face_map->_binomial_coefs_array,
                            patch_order,
                            EVAL_VALUE,
                            uv.array(), 
                            &position);
                    CHECK(position.x() > min.x()); 
                    CHECK(position.x() < max.x());
                    CHECK(position.y() > min.y()); 
                    CHECK(position.y() < max.y());
                    CHECK(position.z() > min.z()); 
                    CHECK(position.z() < max.z());
                }
            }
        }
    }

}*/
