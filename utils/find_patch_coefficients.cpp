#include <iostream>
#include <face_map.hpp>
#include <nummatrix.hpp>
#include <numvector.hpp>
using namespace std;


void dump_to_file(Vector<Point3> coeffs, int num_patches, 
        int poly_bidegree, string output_filepath){
    ofstream file_writer(output_filepath.c_str());
    int num_coeffs = num_patches*(poly_bidegree+1)*(poly_bidegree+1);
    assert(coeffs.length() == num_coeffs);

    file_writer << num_patches;
    file_writer << " ";
    file_writer << poly_bidegree;
    file_writer << endl;
    for(int i = 0; i < num_coeffs; i++){
        for(int d = 0; d < 3; d++){
            file_writer << coeffs(i)(d);
            file_writer << " ";
        }
        file_writer << endl;
    }
    file_writer.close();

}

void to_bernstein_coeffs_2d(){
    int polynomial_degree = 3;
    int num_samples = 3*polynomial_degree+1; 
    //map<Index, int>* binomial_coefs = 
        //FaceMapSurf::compute_binomial_coefs(polynomial_degree);
    IntMatrix* binomial_coefs = 
        FaceMapSurf::compute_binomial_coefs_array(polynomial_degree);
    
    Vector<Point3> control_points(num_samples*num_samples);
    Vector<Point2> sample_points(num_samples*num_samples);
   
    NumMatrix least_squares( num_samples*num_samples, (polynomial_degree+1)*(polynomial_degree+1));
    string output_file = "../wrl_files/poly/two_parallel_curved_patches.poly";
    cout.precision(16);
    for(int i = 0; i < num_samples; i++){
        for(int j = 0; j < num_samples; j++){
            double step = 1./double(num_samples -1);
            
            int index = i*num_samples + j;
            
            double xi =  double(i)*step;
            double yi =  double(j)*step;
            //cout << xi<< "," <<  yi << endl;

            sample_points(index) = Point2(xi,yi);
            // flat patch
            //control_points(index) = Point3(-1 + 2*xi,-1 +2*yi,0);
            
            // smaller flat patch
            //control_points(index) = Point3(-.1 + .2*xi,-.1 +.2*yi,-.00005);
            //control_points(index) = Point3(.1 + .2*xi,.3 +.2*yi,0.);
            
            // slanted flat patches sharing an edge
            //control_points(index) = Point3(-.1 + .1*xi, -.1 + .2*xi,-.1 +.2*yi);
            //control_points(index) = Point3(.1 - .1*xi, -.1 + .2*xi,-.1 +.2*yi);
            
            // slanted flat patches sharing a vertex
            //control_points(index) = Point3(-.1 + .1*xi, -.1 + .2*xi,-.2 +.2*yi);
            //control_points(index) = Point3(.1 - .1*xi, -.1+  .2*xi,0 +.2*yi);
            
            // slanted flat patches intersecting
            //control_points(index) = Point3(-.1 + .2*xi, -.1 + .2*xi,-.1 +.2*yi);
            //control_points(index) = Point3(.1 - .2*xi, -.1 + .2*xi,-.1 +.2*yi);
            
            // curved, parallel close patches 
            double alpha = .2;
            //num_patches = 2;
            //control_points(index) = Point3(-.1 + .2*yi, -.1 + .2*xi, alpha*(pow(xi, 2)+pow(yi,2)));
           //control_points(index) = Point3(-.1 + .2*xi, -.1 + .2*yi, alpha*(pow(xi, 2)+pow(yi,2))+.1);
           
           // parabaloid, max at origin
           //double u = -.1 + .2*xi;
           //double v = -.1 + .2*yi;
           double u = -1 + 2*xi;
           double v = -1 + 2*yi;
           control_points(index) = Point3(u, v, -(pow(u, 2)+pow(v,2)));

            
            // Gaussian, mean mu and standard deviation sigma
            /*
            Point3 mu(.5,.5,0);
            Point3 sigma(1,1,1);
            sigma *= .01;
            Point3 sigma_squared = sigma*sigma;

            Point3 gaussian(
                    2*xi - 1., 2*yi-1.,
                    exp(-(
                              (xi - mu.x())*(xi - mu.x())/(2.*sigma_squared.x())
                            + (yi - mu.y())*(yi - mu.y())/(2.*sigma_squared.y())
                         )
                       ));
            control_points(index) = gaussian;*/
            //control_points(index) = Point3(-1 + 2*xi, -1 + 2*xi,-1 +2*yi);
            //control_points(index) = Point3(1 - 2*xi, -1 + 2*xi,-1 +2*yi);
            //control_points(index) = Point3(-.5 + yi, -.5 + xi, -1e-1);
            //control_points(index) = Point3(-.5 + xi, -.5 + yi, 1e-1);
            
            NumVector x_poly(polynomial_degree+1);
            NumVector y_poly(polynomial_degree+1);

            FaceMapSurf::init_bernstein_polynomial(polynomial_degree, xi, binomial_coefs, x_poly);
            FaceMapSurf::init_bernstein_polynomial(polynomial_degree, yi, binomial_coefs, y_poly);

            for(int di = 0; di <= polynomial_degree; di++){
                for(int dj = 0; dj <= polynomial_degree; dj++){
                    int poly_index = di*(polynomial_degree+1) + dj;
                    least_squares(index, poly_index) = x_poly(di)*y_poly(dj);
                }
            }
        }
    }
    NumMatrix pseudo_inverse( (polynomial_degree+1)*(polynomial_degree+1), num_samples*num_samples);
    setvalue(pseudo_inverse, 0.);
    pinv(least_squares, 1e-10, pseudo_inverse);

    Vector<Point3> bernstein_coeffs(pow(polynomial_degree+1,2));
    for(int i = 0; i < pseudo_inverse.m(); i++){
        for(int j = 0; j < pseudo_inverse.n(); j++){
           bernstein_coeffs(i)  +=  pseudo_inverse(i, j)*control_points(j);
        }
        cout << bernstein_coeffs(i) << endl;
    }
    //dump_to_file(bernstein_coeffs, 1, polynomial_degree, output_file);


}

void to_bernstein_coeffs_1d(){
    int polynomial_degree = 3;
    int num_samples = 5; 
    //map<Index, int>* binomial_coefs = 
     //   FaceMapSurf::compute_binomial_coefs(polynomial_degree);
    IntMatrix* binomial_coefs = 
        FaceMapSurf::compute_binomial_coefs_array(polynomial_degree);

    Vector<Point3> control_points(num_samples);
    Vector<Point2> sample_points(num_samples);

    NumMatrix least_squares( num_samples, (polynomial_degree+1));
    cout.precision(16);
    for(int i = 0; i < num_samples; i++){
        double step = 1./double(num_samples -1);


        double xi =  double(i)*step;

        sample_points(i) = xi;
        control_points(i) = xi*xi*xi - 2*xi*xi;

        NumVector x_poly(polynomial_degree+1);

        FaceMapSurf::init_bernstein_polynomial(polynomial_degree, xi, binomial_coefs, x_poly);

        for(int di = 0; di <= polynomial_degree; di++){
            least_squares(i, di) = x_poly(di);
        }
    }
    NumMatrix pseudo_inverse( (polynomial_degree+1), num_samples);
    setvalue(pseudo_inverse, 0.);
    pinv(least_squares, 1e-10, pseudo_inverse);

    Vector<Point3> bernstein_coeffs(polynomial_degree+1);
    for(int i = 0; i < pseudo_inverse.m(); i++){
        for(int j = 0; j < pseudo_inverse.n(); j++){
            bernstein_coeffs(i)  +=  pseudo_inverse(i, j)*control_points(j);
        }
        cout << bernstein_coeffs(i) << endl;
    }

}

/*

2 3
0.1 0.1 0.4 
0.1 0.0333333 0.2 
0.1 -0.0333333 0.2 
0.1 -0.1 0.2 
0.0333333 0.1 0.2 
0.0333333 0.0333333 2.77556e-16 
0.0333333 -0.0333333 -1.11022e-16 
0.0333333 -0.1 1.52656e-16 
-0.0333333 0.1 0.2 
-0.0333333 0.0333333 -3.33067e-16 
-0.0333333 -0.0333333 5.55112e-17 
-0.0333333 -0.1 1.73472e-17 
-0.1 0.1 0.2 
-0.1 0.0333333 -9.02056e-17 
-0.1 -0.0333333 3.1572e-16 
-0.1 -0.1 3.43475e-16 
0.09999999999999992 0.1 0.4400000000000002
0.09999999999999996 0.0333333333333334 0.2399999999999997
0.09999999999999978 -0.03333333333333346 0.2399999999999997
0.1000000000000002 -0.1 0.2399999999999999
0.03333333333333355 0.09999999999999956 0.2399999999999998
0.03333333333333349 0.03333333333333365 0.04000000000000004
0.03333333333333331 -0.03333333333333358 0.04000000000000009
0.03333333333333345 -0.09999999999999989 0.04000000000000015
-0.03333333333333369 0.1000000000000003 0.2399999999999999
-0.03333333333333335 0.03333333333333288 0.03999999999999981
-0.03333333333333313 -0.03333333333333303 0.03999999999999983
-0.03333333333333337 -0.1000000000000001 0.04000000000000004
-0.09999999999999996 0.09999999999999994 0.2399999999999998
-0.09999999999999998 0.03333333333333346 0.03999999999999981
-0.09999999999999989 -0.03333333333333355 0.04000000000000042
-0.09999999999999998 -0.09999999999999985 0.04000000000000035
 */



int main(){
    to_bernstein_coeffs_2d();
}
