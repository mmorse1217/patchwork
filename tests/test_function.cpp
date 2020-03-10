#include "catch.hpp"
#include "test_patchwork.hpp"
#include "blended_evaluator.hpp"
#include <polynomial_evaluator.hpp>
#include "function.hpp"
#include "utils.hpp"
#include <fstream>


void test_function_values_1d(NumMatrix samples, NumMatrix values, 
        NumMatrix test_samples, NumMatrix test_values,
        Common::BasisType basis_type, 
        int num_control_points,
        int basis_bidegree, double eps_a, double eps_r){
    assert(samples.m() == 1);
    assert(values.m() == 1);
    Function func= Function::construct_spline_function(
            samples, values, num_control_points, basis_bidegree);

    // should interpolate 1 exactly.
    int num_samples = samples.n();
    for(int i = 0; i < num_samples; i++){

        // get the ith sample
        NumVector coords(func.domain_dim(), false, test_samples.clmdata(i));
        // evaluate function 
        NumVector evaluated_values = func(coords);

        //  check the size to be sure
        //CHECK(evaluated_values.length() == func.range_dim());
        for(int d =0 ; d < func.range_dim(); d++){
            CHECK(fabs(test_values(d, i) - evaluated_values(d)) <= eps_a + fabs(test_values(d, i))*eps_r);
        }
    }
}


void test_function_derivatives_1d(NumMatrix samples, NumMatrix values, 
        NumMatrix derivative_samples, NumMatrix derivative_values,
        Common::BasisType basis_type, 
        int num_control_points,
        int basis_bidegree, double eps_a, double eps_r){
    assert(samples.m() == 1);
    assert(values.m() == 1);

    Function func= Function::construct_spline_function(
            samples, values, num_control_points, basis_bidegree);


    // should interpolate 1 exactly.
    int num_samples = int(sqrt(samples.n()));
    
    // evaluate function 
    Function df_dx;
    df_dx = Function::d_dx(0, func);
    for(int i = 0; i < num_samples; i++){

        // get the ith sample
        NumVector coords(func.domain_dim(), false, derivative_samples.clmdata(i));
        NumVector evaluated_values = df_dx(coords);

        //  check the size to be sure
        CHECK(evaluated_values.length() == func.range_dim());

        for(int d =0 ; d < func.range_dim(); d++){
            CHECK(fabs(derivative_values(d, i) - evaluated_values(d)) <= eps_a + fabs(derivative_values(d, i))*eps_r);
        }
    }
}





void test_function_values(NumMatrix samples, NumMatrix values, 
        NumMatrix test_samples, NumMatrix test_values,
        Common::BasisType basis_type, 
        int basis_bidegree, double eps_a, double eps_r){
    Function func;
    if(basis_type == Common::BEZIER){
        func= Function::construct_function(
                samples, values, basis_type, basis_bidegree);
    } else if (basis_type == Common::BSPLINE){
        func= Function::construct_spline_function(
                samples, values, basis_bidegree,5);
        
    }

    int num_samples = int(sqrt(samples.n()));
    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            int index = i*num_samples + j;

            // get the ith sample
            NumVector coords(func.domain_dim(), false, test_samples.clmdata(index));
            // evaluate function 
            NumVector evaluated_values = func(coords);

            //  check the size to be sure
            //CHECK(evaluated_values.length() == func.range_dim());
            for(int d =0 ; d < func.range_dim(); d++){
                CHECK(fabs(test_values(d, index) - evaluated_values(d)) <= eps_a + fabs(test_values(d, index))*eps_r);
            }
        }
    }
}

void test_function_derivative(NumMatrix samples, NumMatrix values, 
        NumMatrix derivative_samples, NumMatrix derivative_values,
        Common::BasisType basis_type, 
        int basis_bidegree, int xi, double eps_a, double eps_r){

    //Function func= Function::construct_function(
    //samples, values, Common::BEZIER, basis_bidegree);
    Function func;
    if(basis_type == Common::BEZIER){
        func= Function::construct_function(
                samples, values, basis_type, basis_bidegree);
    } else if (basis_type == Common::BSPLINE){
        func= Function::construct_spline_function(
                samples, values, basis_bidegree,5);

    }


    // should interpolate 1 exactly.
    int num_samples = int(sqrt(samples.n()));
    int k = xi;
    // evaluate function 
    Function df_dx;
    df_dx = Function::d_dx(k, func);
    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            int index = i*num_samples + j;

            // get the ith sample
            NumVector coords(func.domain_dim(), false, derivative_samples.clmdata(index));
            NumVector evaluated_values = df_dx(coords);

            //  check the size to be sure
            CHECK(evaluated_values.length() == func.range_dim());

            for(int d =0 ; d < func.range_dim(); d++){
                //cout << derivative_values(d, index) << " =?= " <<  evaluated_values(d) << endl;
                CHECK(fabs(derivative_values(d, index) - evaluated_values(d)) <= eps_a + fabs(derivative_values(d, index))*eps_r);
            }
        }
    }
}

void test_function_second_derivative(NumMatrix samples, NumMatrix values, 
        NumMatrix derivative_samples, NumMatrix derivative_values,
        Common::BasisType basis_type, 
        int basis_bidegree, int xi_1, int xi_2,double eps_a, double eps_r){

    //Function func= Function::construct_function(
    //samples, values, Common::BEZIER, basis_bidegree);
    Function func;
    if(basis_type == Common::BEZIER){
        func= Function::construct_function(
                samples, values, basis_type, basis_bidegree);
    } else if (basis_type == Common::BSPLINE){
        func= Function::construct_spline_function(
                samples, values, basis_bidegree,5);

    }


    // should interpolate 1 exactly.
    int num_samples = int(sqrt(samples.n()));
    
    // evaluate function 
    //Function df_dx1;
    //df_dx1 = Function::d_dx(xi_1, func);
    //Function df_dx12;
    //df_dx12 = Function::d_dx(xi_2, df_dx1);
    Function df_dx12;
    df_dx12 = Function::d_dx(xi_2, Function::d_dx(xi_1, func));
    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            int index = i*num_samples + j;

            // get the ith sample
            NumVector coords(func.domain_dim(), false, derivative_samples.clmdata(index));
            NumVector evaluated_values = df_dx12(coords);

            //  check the size to be sure
            CHECK(evaluated_values.length() == func.range_dim());

            for(int d =0 ; d < func.range_dim(); d++){
                //cout << "d: " << d << "; " << derivative_values(d, index) << " =?= " <<  evaluated_values(d) << endl;
                CHECK(fabs(derivative_values(d, index) - evaluated_values(d)) <= eps_a + fabs(derivative_values(d, index))*eps_r);
            }
        }
    }
}


void test_periodic_function_values(NumMatrix spline_coefficients,
        NumMatrix test_samples, NumMatrix test_values,
        Common::BasisType basis_type, 
        int basis_bidegree, double eps_a, double eps_r){
    assert(basis_type == Common::BSPLINE);

    int N = test_samples.m();
    int M = spline_coefficients.m();
    int n = 3;
    int num_coefficients = int(sqrt(spline_coefficients.n()));
    int num_samples = int(sqrt(test_samples.n()));
    int periodic_padding_size = (n-1)/2;
    //int periodic_padding_size = ;
    //int stride = num_coefficients + 2*periodic_padding_size;
    int stride = num_coefficients;// + 2*periodic_padding_size;
cout << "func " << endl;
    Function func(N, M, n, num_coefficients, Common::BSPLINE);

    NumMatrix samples(N, spline_coefficients.n());
    for(int i =0; i < spline_coefficients.n(); i++)
        for(int d =0; d < N; d++)
            samples(d,i) = spline_coefficients(d,i);
    cout << samples << endl;
    func.compute_polynomial_approximation(samples, spline_coefficients);
    NumMatrix padded_coeffs(M,pow(stride,2));
    setvalue(padded_coeffs, 0.);
    //for(int i = n; i < num_coefficients-n; i++){
        //for(int j = n; j < num_coefficients-n; j++){
        /*
    for(int i = 0; i < num_coefficients; i++){
        for(int j = 0; j < num_coefficients; j++){
            int index = (i+periodic_padding_size)*(stride) +(j+periodic_padding_size);
            for(int d = 0; d< M; d++){
                padded_coeffs(d, index) = func._polynomial_coefficents(d,i*num_coefficients+j);
            }
        }
    }*/


    Function periodic_func(N, M, n, stride, Common::BSPLINE);
    //cout << padded_coeffs << endl;


    //func._polynomial_coefficents = spline_coefficients;
    //periodic_func._polynomial_coefficents = padded_coeffs;
   periodic_func._polynomial_coefficents = func._polynomial_coefficents;
    periodic_func.set_periodicity(Common::PERIODIC_XY, 1., true);


    //cout << func._polynomial_coefficents << endl;
    for(int d =0; d < M; d++){
        for(int i= 0; i < stride; i++){
            for(int j =0; j < stride; j++){
                cout << periodic_func._polynomial_coefficents(d, i*stride + j) << ", ";
            }
            cout << endl;
        }
        cout << endl;
    }
    {
    ofstream out("tests/spline_coeffs.csv");
    //out << spline_values;
    out << periodic_func._polynomial_coefficents;
    //out << spline_coefficients;
    //out << padded_coeffs;
    out.close();
    }


    //func._polynomial_coefficents = padded_coeffs;
    
    //Function func= Function::construct_spline_function(
    //            samples, values, basis_bidegree,5);
    cout << "finished setup" << endl;
    //func.set_periodicity(Common::PERIODIC_XY, 1.);
    cout << "finished periodicity" << endl;
    NumMatrix spline_values(M, num_samples*num_samples);

    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            int index = i*num_samples + j;

            // get the ith sample
            NumVector coords(periodic_func.domain_dim(), false, test_samples.clmdata(index));
            
            // evaluate spline approximation 
            NumVector evaluated_values = periodic_func(coords);

            //  check the size to be sure
            //CHECK(evaluated_values.length() == func.range_dim());
                //cout << evaluated_values << endl;
            for(int d =0 ; d < M; d++){
                spline_values(d,index) = evaluated_values(d);
                //cout << test_values(d, index) << ", " << evaluated_values(d) << endl;
                //cout << fabs(test_values(d, index) - evaluated_values(d)) << endl;
                CHECK(fabs(test_values(d, index) - evaluated_values(d)) <= eps_a + fabs(test_values(d, index))*eps_r);
            }
        }
    }
    ofstream out("tests/spline_points.csv");
    out << spline_values;
    //out << func._polynomial_coefficents;
    //out << spline_coefficients;
    //out << padded_coeffs;
    out.close();
}

void test_periodic_function_derivative(NumMatrix spline_coefficients, 
        NumMatrix derivative_samples, NumMatrix derivative_values,
        Common::BasisType basis_type, 
        int basis_bidegree, int xi, double eps_a, double eps_r){
    assert(basis_type == Common::BSPLINE);
    assert(derivative_samples.n() == derivative_values.n());

    int N = derivative_samples.m();
    int M = spline_coefficients.m();
    int n = 3;
    int num_coefficients = int(sqrt(spline_coefficients.n()));
    int num_samples = int(sqrt(derivative_samples.n()));
    int periodic_padding_size = (n-1)/2;
    //int periodic_padding_size = (n);
    //int periodic_padding_size = n;
    int stride = num_coefficients + 2*periodic_padding_size;
cout << "func " << endl;
    Function func(N, M, n, stride, Common::BSPLINE);
    NumMatrix padded_coeffs(M,pow(stride,2));
    setvalue(padded_coeffs, 0.);
    //for(int i = n; i < num_coefficients-n; i++){
        //for(int j = n; j < num_coefficients-n; j++){
    for(int i = 0; i < num_coefficients; i++){
        for(int j = 0; j < num_coefficients; j++){
            int index = (i+periodic_padding_size)*(stride) +(j+periodic_padding_size);
            for(int d = 0; d< M; d++){
                padded_coeffs(d, index) = spline_coefficients(d,i*num_coefficients+j);
            }
        }
    }


    //cout << padded_coeffs << endl;


    //func._polynomial_coefficents = spline_coefficients;
    func._polynomial_coefficents = padded_coeffs;
    func._periodic = Common::PERIODIC_XY;

func.set_periodicity(Common::PERIODIC_XY, 1.);

    //cout << func._polynomial_coefficents << endl;
    for(int d =0; d < M; d++){
        for(int i= 0; i < stride; i++){
            for(int j =0; j < stride; j++){
                cout << func._polynomial_coefficents(d, i*stride + j) << ", ";
            }
            cout << endl;
        }
        cout << endl;
    }
    {
    ofstream out("tests/spline_coeffs.csv");
    //out << spline_values;
    out << func._polynomial_coefficents;
    //out << spline_coefficients;
    //out << padded_coeffs;
    out.close();
    }
    
    //Function func= Function::construct_spline_function(
    //            samples, values, basis_bidegree,5);
    //func.enforce_periodicity();
    //func.enforce_periodicity(Common::PERIODIC_XY, 1.);
    //func.set_periodicity(Common::PERIODIC_XY, 1.);



    // should interpolate 1 exactly.
    //int num_samples = int(sqrt(samples.n()));
    int k = xi;
    // evaluate function 
    Function df_dx;
    df_dx = Function::d_dx(k, func);
    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            int index = i*num_samples + j;

            // get the ith sample
            NumVector coords(func.domain_dim(), false, derivative_samples.clmdata(index));
            NumVector evaluated_values = df_dx(coords);

            //  check the size to be sure
            CHECK(evaluated_values.length() == func.range_dim());

            for(int d =2 ; d < 3; d++){
                //cout << derivative_values(d, index) << " =?= " <<  evaluated_values(d) << endl;
                CHECK(fabs(derivative_values(d-2, index) - evaluated_values(d)) <= eps_a + fabs(derivative_values(d-2, index))*eps_r);
            }
        }
    }
}




void sample_function_R1_to_R1( double (*scalar_func)(double), 
       NumMatrix samples, NumMatrix& values) {
    
    assert(samples.n() == values.n());
    assert(samples.m() == 1);
    assert(values.m() == 1);

    int num_samples = samples.n();
    for(int i = 0; i < num_samples; i++){
        double xi =  samples(0,i);
        //values(0,index) = pow(xi, 6)*pow(yj, 6)+7.;
        values(0,i) = scalar_func(xi);
    }
}



void sample_function_R2_to_R1( double (*scalar_func)(double, double), 
       NumMatrix samples, NumMatrix& values) {
    
    assert(samples.n() == values.n());
    assert(samples.m() == 2);
    assert(values.m() == 1);

    int num_samples = int(sqrt(samples.n()));
    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            int index = i*num_samples + j;
            double xi =  samples(0,index);
            double yj =  samples(1,index);
            //values(0,index) = pow(xi, 6)*pow(yj, 6)+7.;
            values(0,index) = scalar_func(xi, yj);
        }
    }


}


void sample_function_R2_to_R3( Point3 (*vector_func)(double, double), 
       NumMatrix samples, NumMatrix& values) {
    
    assert(samples.n() == values.n());
    assert(samples.m() == 2);
    assert(values.m() == 3);

    int num_samples = int(sqrt(samples.n()));
    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            int index = i*num_samples + j;
            double xi =  samples(0,index);
            double yj =  samples(1,index);

            Point3 value = vector_func(xi,yj);

            for(int d = 0; d < values.m(); d++){
                values(d,index) = value(d);
            }
        }
    }


}


TEST_CASE("","[function]"){
    int basis_bidegree = 12;
    int num_samples = 2*(basis_bidegree+1);
    //int num_samples = 20;
    double step_size = 1./(num_samples-1);

    const int N = 2;

    // uniform samples on [0,1]^2
    NumMatrix equispaced_samples(N, num_samples*num_samples);
    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            int index = i*num_samples + j;

            equispaced_samples(0,index) = i*step_size;
            equispaced_samples(1,index) = j*step_size;
        }
    }
    // uniform samples in [.25, .75]^2. By design half of the points will be 
    // unused in the original fitting
    NumMatrix off_node_samples(N, num_samples*num_samples);
    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            int index = i*num_samples + j;

            off_node_samples(0,index) = i*step_size/2. +.25;
            off_node_samples(1,index) = j*step_size/2.+ .25;
        }
    }
    SECTION("Test constant function"){

        int basis_bidegree = 6;
        int num_samples = 2*basis_bidegree;
        double step_size = 1./(num_samples-1);

        const int N = 2;

        // uniform samples on [0,1]^2
        NumMatrix equispaced_samples(N, num_samples*num_samples);
        for(int i = 0; i < num_samples; i++){
            for(int j = 0; j < num_samples; j++){
                int index = i*num_samples + j;

                equispaced_samples(0,index) = i*step_size;
                equispaced_samples(1,index) = j*step_size;
            }
        }
        // uniform samples in [.25, .75]^2. By design half of the points will be 
        // unused in the original fitting
        NumMatrix off_node_samples(N, num_samples*num_samples);
        for(int i = 0; i < num_samples; i++){
            for(int j = 0; j < num_samples; j++){
                int index = i*num_samples + j;

                off_node_samples(0,index) = i*step_size/2. +.25;
                off_node_samples(1,index) = j*step_size/2.+ .25;
            }
        }
        const int M = 1;

        NumMatrix constant(M, num_samples*num_samples);
        setvalue(constant, 1.); 
        test_function_values(equispaced_samples, constant, equispaced_samples, constant, Common::BEZIER, basis_bidegree, 1e-14, 1e-14);
        test_function_values(equispaced_samples, constant, off_node_samples, constant, Common::BEZIER, basis_bidegree, 1e-14, 1e-14);
        
        NumMatrix zero(M, num_samples*num_samples);
        setvalue(zero, 0.); 
        test_function_derivative(equispaced_samples, zero, equispaced_samples, zero,Common::BEZIER,  basis_bidegree,0, 1e-14, 1e-14);
        test_function_derivative(equispaced_samples, zero, off_node_samples, zero, Common::BEZIER, basis_bidegree, 0, 1e-14, 1e-14);
        test_function_derivative(equispaced_samples, zero, equispaced_samples, zero, Common::BEZIER, basis_bidegree,1, 1e-14, 1e-14);
        test_function_derivative(equispaced_samples, zero, off_node_samples, zero, Common::BEZIER, basis_bidegree, 1, 1e-14, 1e-14);


    }
    SECTION("Test constant function, R^2 -> R^10"){
        int basis_bidegree = 6;
        int num_samples = 2*basis_bidegree;
        double step_size = 1./(num_samples-1);

        const int N = 2;

        // uniform samples on [0,1]^2
        NumMatrix equispaced_samples(N, num_samples*num_samples);
        for(int i = 0; i < num_samples; i++){
            for(int j = 0; j < num_samples; j++){
                int index = i*num_samples + j;

                equispaced_samples(0,index) = i*step_size;
                equispaced_samples(1,index) = j*step_size;
            }
        }
        // uniform samples in [.25, .75]^2. By design half of the points will be 
        // unused in the original fitting
        NumMatrix off_node_samples(N, num_samples*num_samples);
        for(int i = 0; i < num_samples; i++){
            for(int j = 0; j < num_samples; j++){
                int index = i*num_samples + j;

                off_node_samples(0,index) = i*step_size/2. +.25;
                off_node_samples(1,index) = j*step_size/2.+ .25;
            }
        }
        const int M = 10;

        NumMatrix constant(M, num_samples*num_samples);
        setvalue(constant, 1.);

        test_function_values(equispaced_samples, constant, equispaced_samples, constant, Common::BEZIER, basis_bidegree, 1e-14, 1e-14);
        test_function_values(equispaced_samples, constant, off_node_samples, constant, Common::BEZIER, basis_bidegree, 1e-14, 1e-14);

        NumMatrix zero(M, num_samples*num_samples);
        setvalue(zero, 0.); 
        test_function_derivative(equispaced_samples, zero, equispaced_samples, zero, Common::BEZIER, basis_bidegree,0, 1e-14, 1e-14);
        test_function_derivative(equispaced_samples, zero, off_node_samples, zero, Common::BEZIER, basis_bidegree, 0, 1e-14, 1e-14);
        test_function_derivative(equispaced_samples, zero, equispaced_samples, zero, Common::BEZIER, basis_bidegree,1, 1e-14, 1e-14);
        test_function_derivative(equispaced_samples, zero, off_node_samples, zero, Common::BEZIER, basis_bidegree, 1, 1e-14, 1e-14);

    }
    SECTION("Test degree six polynomial"){
        const int M = 1;

        NumMatrix polynomial(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sixth_degree_poly, equispaced_samples, polynomial);
        // check that approximated values at the nodes should be exact
        test_function_values(equispaced_samples, polynomial, equispaced_samples, polynomial, Common::BEZIER, basis_bidegree, 1e-11, 1e-11);

        NumMatrix off_node_polynomial(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sixth_degree_poly, off_node_samples, off_node_polynomial);
        
        // should be able to resolve the polynomial accurately
        test_function_values(equispaced_samples, polynomial, off_node_samples, off_node_polynomial, Common::BEZIER, basis_bidegree, 1e-11, 1e-11);

        // same test for derivative
        NumMatrix poly_dx(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sixth_degree_poly_dx, off_node_samples, poly_dx);

        // test computing derivatives properly
        test_function_derivative(equispaced_samples, polynomial, off_node_samples, poly_dx, Common::BEZIER, basis_bidegree, 0, 1e-10, 1e-10);

        NumMatrix poly_dy(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sixth_degree_poly_dy, off_node_samples, poly_dy);

        test_function_derivative(equispaced_samples, polynomial, off_node_samples, poly_dy, Common::BEZIER, basis_bidegree, 1, 1e-10, 1e-10);
    }
    SECTION("Test trig function approximation"){
        const int M = 1;

        NumMatrix trig_func(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sin2cos2, equispaced_samples, trig_func);

        // values at the nodes should be close  but not exact
        test_function_values(equispaced_samples, trig_func, equispaced_samples, trig_func, Common::BEZIER, basis_bidegree, 0, 1e-8);

        NumMatrix off_node_trig_func(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sin2cos2, off_node_samples, off_node_trig_func);

        test_function_values(equispaced_samples, trig_func, off_node_samples, off_node_trig_func, Common::BEZIER, basis_bidegree, 0, 1e-6);;

        NumMatrix trig_values_dx(M, num_samples*num_samples);
        sample_function_R2_to_R1(&trig_dx, off_node_samples, trig_values_dx);

        // test computing derivatives properly
        test_function_derivative(equispaced_samples, trig_func, off_node_samples, trig_values_dx, Common::BEZIER, basis_bidegree, 0, 1e-5, 1e-5);

        NumMatrix trig_values_dy(M, num_samples*num_samples);
        sample_function_R2_to_R1(&trig_dy, off_node_samples, trig_values_dy);

        test_function_derivative(equispaced_samples, trig_func, off_node_samples, trig_values_dy, Common::BEZIER, basis_bidegree, 1, 1e-5, 1e-5);

    }
    SECTION("Test R^2 \to R^3 approximation"){
        const int M = 3;
        //basis_bidegree = 20;

        NumMatrix func_values(M, num_samples*num_samples);
        sample_function_R2_to_R3(&paraboloid, equispaced_samples, func_values);

        // values at the nodes should be exact
        //test_function_values(equispaced_samples, func_values, equispaced_samples, func_values, Common::BEZIER, basis_bidegree, 1e-8, 1e-8);

        NumMatrix off_node_func_values(M, num_samples*num_samples);
        sample_function_R2_to_R3(&paraboloid, off_node_samples, off_node_func_values);

        //test_function_values(equispaced_samples, func_values, off_node_samples, off_node_func_values, Common::BEZIER, basis_bidegree, 1e-6, 1e-6);;
        NumMatrix func_values_dx(M, num_samples*num_samples);
        sample_function_R2_to_R3(&paraboloid_dx, off_node_samples, func_values_dx);

        // test computing derivatives properly
        test_function_derivative(equispaced_samples, func_values, off_node_samples, func_values_dx, Common::BEZIER, basis_bidegree, 0, 1e-7, 1e-7);

        NumMatrix func_values_dy(M, num_samples*num_samples);
        sample_function_R2_to_R3(&paraboloid_dy, off_node_samples, func_values_dy);

        test_function_derivative(equispaced_samples, func_values, off_node_samples, func_values_dy, Common::BEZIER, basis_bidegree, 1, 1e-7, 1e-7);
    } /*
         SECTION("Test original legacy face-map implementation vs. function based version"){
        char filea[300]; char fileb[300]; char filec[300]; char filepoly[300];
        int cl; int pc;	int ct;   
        int bt;  int pp; double lb; double ub;int gl;double fs; int cot;  
        int sd; int poudeg; int lvl; int gen; int refinement_factor; 
        int patch_order; bool adaptive; 

        load_blended_options(
                filea,  fileb, filec,  cl,  pc,	 ct,   
                bt,  pp, lb,  ub, gl, fs,  cot,  
                sd, poudeg, lvl,  gen,  refinement_factor, 
                patch_order,  adaptive,  filepoly);
        ifstream tmpa(filea);
        CCSubMatLib* submatlib = new CCSubMatLib();
        cerr<<"submatlib setup"<<endl;
        submatlib->setup(tmpa);
        
        ifstream tmpb(fileb);
        BdULib* bdulib = new BdULib();
        cerr<<"bdulib setup"<<endl;
        bdulib->setup(tmpb);
        ifstream tmpc(filec);
        BdSurf* bdsurf = new BdSurf();
        bdsurf->setParams(cl,pc,ct,bt,pp, lb,ub, gl, fs, cot, sd ,poudeg, refinement_factor, 0);
        bdsurf->setup(tmpc,0, Point3(1.0), (void*)submatlib, (void*)bdulib);

        BlendedEvaluator* blended_evaluator = new BlendedEvaluator(bdsurf);
        FaceMapSurf* legacy_face_map = new FaceMapSurf(blended_evaluator);
        FaceMapSurf* function_face_map = new FaceMapSurf(blended_evaluator);
        
        legacy_face_map->_legacy = true;
        function_face_map->_legacy = false;
        
        patch_order = 6;
        refinement_factor = 0;
        adaptive = 1;
        
        legacy_face_map->set_params(refinement_factor, patch_order, adaptive);
        function_face_map->set_params(refinement_factor, patch_order, adaptive);
        
        legacy_face_map->setup();
        function_face_map->setup();

        // should have the same number of patches
        CHECK(legacy_face_map->_patches.size() == function_face_map->_patches.size());
        // fail if not :(
        assert(legacy_face_map->_patches.size() == function_face_map->_patches.size());
        int num_patches = legacy_face_map->_patches.size();
        for(int pi = 0; pi < num_patches; pi++){
            int num_samples = 2*patch_order;
            double step = 1./double(num_samples-1);
            for(int i = 0; i < num_samples; i++){
                for(int j = 0; j < num_samples; j++){
                    Point2 uv(i*step, j*step);
                    
                    vector<Point3> legacy_values(6);
                    vector<Point3> function_values(6);
                    
                    legacy_face_map->eval(
                            EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV,
                            pi,
                            uv.array(),
                            legacy_values.data());
                    
                    function_face_map->eval(
                            EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV,
                            pi,
                            uv.array(),
                            function_values.data());
                    for(int k =0; k < 6; k++){
                        // not enforcing C^0 on functional face-map, so this
                        // should tighten the abs error down
                        CHECK((legacy_values[k] - function_values[k]).l2() <=1e-3);
                    }
                }
            }

        }
    }*/
}
/*TEST_CASE("","[function]"){
    SECTION("Test degree 0-15 polynomial b-spline 1d"){
        // should be able to exactly recover sampled degree n polynomial function exactly
        // with n+1 control points

        for(int deg = 0; deg < 10; deg++){
            // at least 10 to avoid possible sampling issues
            int num_samples = max(3*(deg+1), 10);
            double step_size = 1./(num_samples-1);

            const int N = 1;

            // uniform samples on [0,1] to be used for function fitting
            NumMatrix equispaced_samples(N, num_samples);
            for(int i = 0; i < num_samples; i++){
                equispaced_samples(0,i) = i*step_size;
            }
            // uniform samples in [.25, .75]^2. 
            NumMatrix off_node_samples(N, num_samples);
            for(int i = 0; i < num_samples; i++){
                off_node_samples(0,i) = i*step_size/2. +.25;
            }
            const int M = 1;


            // evaluate polynomial to approximate at equispaced nodes
            NumMatrix polynomial(M, num_samples);
            setvalue(polynomial,0.);

            for(int i = 0; i < num_samples; i++){
                double xi =  equispaced_samples(0,i);
                for(int k = 0; k < deg; k++)
                    polynomial(0,i) += k*pow(xi, deg);
            }

            // check that approximated values at the nodes is exact
            test_function_values_1d(equispaced_samples, polynomial, equispaced_samples, polynomial, Common::BSPLINE, deg+1, deg, 1e-13, 1e-13);
            // evaluate derivative of polynomial at equispaced nodes
            NumMatrix polynomial_dx(M, num_samples);
            setvalue(polynomial_dx,0.);

            for(int i = 0; i < num_samples; i++){
                double xi =  equispaced_samples(0,i);
                for(int k = 0; k < deg; k++)
                    polynomial_dx(0,i) += double(deg)*k*pow(xi, deg-1);
            }
            //test_function_derivatives_1d(equispaced_samples, polynomial, equispaced_samples, polynomial_dx, Common::BSPLINE, deg+1, deg, 1e-12, 1e-12);

            // evaluate polynomial at target points to check fit accuracy
            NumMatrix off_node_polynomial(M, num_samples);
            setvalue(off_node_polynomial,0.);
            for(int i = 0; i < num_samples; i++){
                double xi =  off_node_samples(0,i);
                for(int k = 0; k < deg; k++)
                    off_node_polynomial(0,i) += k*pow(xi, deg);
            }

            // should be able to exactly resolve the polynomial accurately away from the
            // nodes used to fit
            test_function_values_1d(equispaced_samples, polynomial, off_node_samples,  off_node_polynomial, Common::BSPLINE, deg+1, deg, 1e-13, 1e-13);
            // same test for derivative

            // TODO 
        } 
    }
}*/
TEST_CASE("","[function]"){

        //int basis_bidegree = 25;
        int basis_bidegree = 15;
        int num_samples = 2*basis_bidegree; // is this under sampled? might need 4x degree to resolve cubics per patch
        //int num_samples = 10;
        double step_size = 1./(num_samples-1);

        const int N = 2;

        // uniform samples on [0,1]^2
        NumMatrix equispaced_samples(N, num_samples*num_samples);
        for(int i = 0; i < num_samples; i++){
            for(int j = 0; j < num_samples; j++){
                int index = i*num_samples + j;

                equispaced_samples(0,index) = i*step_size;
                equispaced_samples(1,index) = j*step_size;
            }
        }
        // uniform samples in [.25, .75]^2. By design half of the points will be 
        // unused in the original fitting
        NumMatrix off_node_samples(N, num_samples*num_samples);
        for(int i = 0; i < num_samples; i++){
            for(int j = 0; j < num_samples; j++){
                int index = i*num_samples + j;

                off_node_samples(0,index) = i*step_size/2. +.25;
                off_node_samples(1,index) = j*step_size/2.+ .25;
            }
        }
        const int M = 1;

SECTION("Test constant function b spline 2d"){
        NumMatrix constant(M, num_samples*num_samples);
        setvalue(constant, 1.); 
        test_function_values(equispaced_samples, constant, equispaced_samples, constant, Common::BSPLINE, basis_bidegree, 1e-12, 1e-12);
        test_function_values(equispaced_samples, constant, off_node_samples, constant, Common::BSPLINE, basis_bidegree, 1e-12, 1e-12);
        
        NumMatrix zero(M, num_samples*num_samples);
        setvalue(zero, 0.); 
        test_function_derivative(equispaced_samples, zero, equispaced_samples, zero,Common::BSPLINE,  basis_bidegree,0, 1e-14, 1e-14);
        test_function_derivative(equispaced_samples, zero, off_node_samples, zero, Common::BSPLINE, basis_bidegree, 0, 1e-14, 1e-14);
        test_function_derivative(equispaced_samples, zero, equispaced_samples, zero, Common::BSPLINE, basis_bidegree,1, 1e-14, 1e-14);
        test_function_derivative(equispaced_samples, zero, off_node_samples, zero, Common::BSPLINE, basis_bidegree, 1, 1e-14, 1e-14);


    }
    SECTION("Test degree six polynomial bspline 2d"){
        //const int M = 1;

        NumMatrix polynomial(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sixth_degree_poly, equispaced_samples, polynomial);
        // check that approximated values at the nodes should be exact
        test_function_values(equispaced_samples, polynomial, equispaced_samples, polynomial, Common::BSPLINE, basis_bidegree, 1e-7, 1e-7);

        NumMatrix off_node_polynomial(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sixth_degree_poly, off_node_samples, off_node_polynomial);
        
        // should be able to resolve the polynomial accurately
        test_function_values(equispaced_samples, polynomial, off_node_samples, off_node_polynomial, Common::BSPLINE, basis_bidegree, 1e-6, 1e-6);
        // same test for derivative
        NumMatrix poly_dx(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sixth_degree_poly_dx, off_node_samples, poly_dx);

        // test computing derivatives properly
        test_function_derivative(equispaced_samples, polynomial, off_node_samples, poly_dx, Common::BSPLINE, basis_bidegree, 0, 1e-6, 1e-6);

        NumMatrix poly_dy(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sixth_degree_poly_dy, off_node_samples, poly_dy);

        test_function_derivative(equispaced_samples, polynomial, off_node_samples, poly_dy, Common::BSPLINE, basis_bidegree, 1, 1e-6, 1e-6);
    }
    SECTION("Test trig function approximation b spline 2d"){
        const int M = 1;

        NumMatrix trig_func(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sin2cos2, equispaced_samples, trig_func);

        // values at the nodes should be close  but not exact
        test_function_values(equispaced_samples, trig_func, equispaced_samples, trig_func, Common::BSPLINE, basis_bidegree, 1e-5, 1e-5);

        NumMatrix off_node_trig_func(M, num_samples*num_samples);
        sample_function_R2_to_R1(&sin2cos2, off_node_samples, off_node_trig_func);

        test_function_values(equispaced_samples, trig_func, off_node_samples, off_node_trig_func, Common::BSPLINE, basis_bidegree, 1e-5, 1e-5);
        NumMatrix trig_values_dx(M, num_samples*num_samples);
        sample_function_R2_to_R1(&trig_dx, off_node_samples, trig_values_dx);

        // test computing derivatives properly
        test_function_derivative(equispaced_samples, trig_func, off_node_samples, trig_values_dx, Common::BSPLINE, basis_bidegree, 0, 1e-5, 1e-5);

        NumMatrix trig_values_dy(M, num_samples*num_samples);
        sample_function_R2_to_R1(&trig_dy, off_node_samples, trig_values_dy);

        test_function_derivative(equispaced_samples, trig_func, off_node_samples, trig_values_dy, Common::BSPLINE, basis_bidegree, 1, 1e-5, 1e-5);

    }
    SECTION("Test R^2 \to R^3 approximation"){
        const int M = 3;
        //basis_bidegree = 20;

        NumMatrix func_values(M, num_samples*num_samples);
        sample_function_R2_to_R3(&paraboloid, equispaced_samples, func_values);

        // values at the nodes should be exact
        test_function_values(equispaced_samples, func_values, equispaced_samples, func_values, Common::BSPLINE, basis_bidegree, 1e-8, 1e-8);

        NumMatrix off_node_func_values(M, num_samples*num_samples);
        sample_function_R2_to_R3(&paraboloid, off_node_samples, off_node_func_values);

        test_function_values(equispaced_samples, func_values, off_node_samples, off_node_func_values, Common::BSPLINE, basis_bidegree, 1e-7, 1e-7);;
        NumMatrix func_values_dx(M, num_samples*num_samples);
        sample_function_R2_to_R3(&paraboloid_dx, off_node_samples, func_values_dx);

        // test computing derivatives properly
        test_function_derivative(equispaced_samples, func_values, off_node_samples, func_values_dx, Common::BSPLINE, basis_bidegree, 0, 1e-7, 1e-7);

        NumMatrix func_values_dy(M, num_samples*num_samples);
        sample_function_R2_to_R3(&paraboloid_dy, off_node_samples, func_values_dy);

        test_function_derivative(equispaced_samples, func_values, off_node_samples, func_values_dy, Common::BSPLINE, basis_bidegree, 1, 1e-7, 1e-7);

        NumMatrix paraboloid_dxdy(M, num_samples*num_samples);
        setvalue(paraboloid_dxdy, 0.);
        //sample_function_R2_to_R3(&zero, off_node_samples, paraboloid_dxdy);
        //sample_function_R2_to_R3(&paraboloid_dxdx, off_node_samples, paraboloid_dxdy);
        
        test_function_second_derivative(equispaced_samples, func_values, 
                off_node_samples, paraboloid_dxdy,
                Common::BSPLINE, 
                basis_bidegree, 0, 1, 1e-6, 1e-6);


    }
}
TEST_CASE("","[periodic]"){

    //int basis_bidegree = 25;
    int basis_bidegree = 3; 
    //int num_samples = 2*basis_bidegree+1; // is this under sampled? might need 4x degree to resolve cubics per patch
    int num_samples = 10;
    //int num_samples = 10;
    double step_size = 1./(num_samples-1);

    const int N = 2;

    // uniform samples on [0,1]^2
    NumMatrix equispaced_samples(N, num_samples*num_samples);
    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            int index = i*num_samples + j;

            equispaced_samples(0,index) = i*step_size;
            equispaced_samples(1,index) = j*step_size;
        }
    }
    int upsample_rate = 2;
    double upsample_step = 1./(upsample_rate*num_samples-1);
    // uniform samples in [.25, .75]^2. By design half of the points will be 
    // unused in the original fitting
    NumMatrix off_node_samples(N, upsample_rate*upsample_rate*num_samples*num_samples);
    for(int i = 0; i < upsample_rate*num_samples; i++){
        for(int j = 0; j < upsample_rate*num_samples; j++){
            int index = i*(upsample_rate*num_samples) + j;

            off_node_samples(0,index) = i*upsample_step;
            off_node_samples(1,index) = j*upsample_step;
        }
    }

    SECTION("Test constant function b spline 2d"){
        const int M = 3;
        NumMatrix constant(M, num_samples*num_samples);
        NumMatrix off_node_func_values(M, upsample_rate*num_samples*upsample_rate*num_samples);
        
        sample_function_R2_to_R3(&constant_plane, equispaced_samples, constant);
        sample_function_R2_to_R3(&constant_plane, off_node_samples, off_node_func_values);
        
        //test_periodic_function_values(constant, equispaced_samples, constant, Common::BSPLINE, basis_bidegree, 1e-13, 1e-13);
        //test_periodic_function_values(constant, off_node_samples, off_node_func_values, Common::BSPLINE, basis_bidegree, 1e-13, 1e-13);
        /*
        sample_function_R2_to_R3(&constant_plane_dx, equispaced_samples, constant);
        sample_function_R2_to_R3(&constant_plane_dx, off_node_samples, off_node_func_values);
        //test_periodic_function_derivative(constant, equispaced_samples, constant,Common::BSPLINE,  basis_bidegree,0, 1e-13, 1e-13);
        test_periodic_function_derivative(equispaced_samples, zero, off_node_samples, zero, Common::BSPLINE, basis_bidegree, 0, 1e-13, 1e-13);
        test_periodic_function_derivative(equispaced_samples, zero, equispaced_samples, zero, Common::BSPLINE, basis_bidegree,1, 1e-13, 1e-13);
        test_periodic_function_derivative(equispaced_samples, zero, off_node_samples, zero, Common::BSPLINE, basis_bidegree, 1, 1e-13, 1e-13);
        */ 

    }
    SECTION("Test sin(2pi*x)*sin(2pi*y)"){
        const int M = 3;
        //basis_bidegree = 20;

        NumMatrix func_values(M, num_samples*num_samples);
        sample_function_R2_to_R3(&sinsin, equispaced_samples, func_values);
        cout << "function samples: " << endl;
        // values at the nodes should be exact
        NumMatrix off_node_func_values(M, upsample_rate*num_samples*upsample_rate*num_samples);
        sample_function_R2_to_R3(&sinsin, off_node_samples, off_node_func_values);
    ofstream out("tests/func_values.csv");
    //out << off_node_func_values;
    out << func_values;
    //out << constant;
    out.close();
        test_periodic_function_values( func_values, equispaced_samples, func_values, Common::BSPLINE, basis_bidegree, 1e-3, 1e-3);

        //test_periodic_function_values(func_values, off_node_samples, off_node_func_values, Common::BSPLINE, basis_bidegree, 1e-3, 1e-3);;
    }/*
    SECTION("Test bump function b spline 2d"){
        NumMatrix constant(M, num_samples*num_samples);
        setvalue(constant, 1.); 
        test_periodic_function_values(equispaced_samples, constant, equispaced_samples, constant, Common::BSPLINE, basis_bidegree, 1e-13, 1e-13);
        test_periodic_function_values(equispaced_samples, constant, off_node_samples, constant, Common::BSPLINE, basis_bidegree, 1e-13, 1e-13);

        NumMatrix zero(M, num_samples*num_samples);
        setvalue(zero, 0.); 
        test_periodic_function_derivative(equispaced_samples, zero, equispaced_samples, zero,Common::BSPLINE,  basis_bidegree,0, 1e-13, 1e-13);
        test_periodic_function_derivative(equispaced_samples, zero, off_node_samples, zero, Common::BSPLINE, basis_bidegree, 0, 1e-13, 1e-13);
        test_periodic_function_derivative(equispaced_samples, zero, equispaced_samples, zero, Common::BSPLINE, basis_bidegree,1, 1e-13, 1e-13);
        test_periodic_function_derivative(equispaced_samples, zero, off_node_samples, zero, Common::BSPLINE, basis_bidegree, 1, 1e-13, 1e-13);


    }*/
}

