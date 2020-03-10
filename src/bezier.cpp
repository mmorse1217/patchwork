#include "bezier.hpp"

void Bezier::set_value_edge_case(size_t n, 
        size_t i, 
        size_t j,  
        size_t coeff_dim,   
        Common::Direction dx_i,
        NumMatrix original_coeffs, NumMatrix& derivative){
    //int index = first_term.first*n + first_term.second;
    //int index_plus_one = first_term.first*n + first_term.second;
    int index, index_plus_one;
    if(dx_i == Common::U){
        if(i == 0){
            index = i*(n+1) + j;
            index_plus_one = (i+1)*(n+1) + j;

        } else if (i == n){
            index = (i-1)*(n+1) + j;
            index_plus_one = i*(n+1) + j;
        } else {
            // should only be called for i == 0 or i == n
            assert(0);
        }
    } else if(dx_i == Common::V){
        if(j == 0){
            index = i*(n+1) + j;
            index_plus_one = i*(n+1) + j+1;

        } else if (j == n){
            index = i*(n+1) + j-1;
            index_plus_one = i*(n+1) + j;
        } else {
            // should only be called for j == 0 or j == n
            assert(0);
        }
    }
    for(int k = 0; k < coeff_dim; k++){
        double a = original_coeffs(k,index);
        double a_plus_one = original_coeffs(k,index_plus_one); 
        //derivative(k, index) = double(n)*(a_plus_one - a);
        derivative(k, i*(n+1)+j) = double(n)*(a_plus_one - a);
    }


}

void Bezier::set_value(size_t n, 
        size_t i, 
        size_t j,  
        size_t coeff_dim,   
        Common::Direction dx_i,
        NumMatrix original_coeffs, NumMatrix& derivative){
    assert((Common::V == dx_i && (j != 0 && j != n)) | 
            (Common::U == dx_i && ( i != 0 && i != n)));
    int index, index_plus_one, index_minus_one;
    double coeff, coeff_plus_one, coeff_minus_one;
    
    double i_double = double(i);
    double j_double = double(j);

    if(dx_i == Common::U){
        index_minus_one = (i-1)*(n+1) + j;
        index           = i*(n+1) + j;
        index_plus_one  = (i+1)*(n+1) + j;
        
        coeff_minus_one = -i_double;
        coeff           = 2*i_double - n;
        coeff_plus_one  = n-i_double;

    } else if(dx_i == Common::V){
        index_minus_one = i*(n+1) + (j-1);
        index           = i*(n+1) + j;
        index_plus_one  = i*(n+1) + (j+1);

        coeff_minus_one = -j_double;
        coeff           = 2*j_double - n;
        coeff_plus_one  = n-j_double;
    }
    for(int k = 0; k < coeff_dim; k++){
        double a_minus_one  = original_coeffs(k,index_minus_one); 
        double a            = original_coeffs(k,index);
        double a_plus_one   = original_coeffs(k,index_plus_one); 
        /*
        cout << "a_{i-1}: " << a_minus_one << endl;
        cout << "a_{i}: " << a<< endl;
        cout << "a_{i+1}: " << a_plus_one << endl;
        cout << "c_{i-1}: " << coeff_minus_one<< endl;
        cout << "c_{i}: " << coeff<< endl;
        cout << "c_{i+1}: " << coeff_plus_one<< endl;
       */ 
        derivative(k, index) = coeff_plus_one*a_plus_one
                             + coeff*a
                             + coeff_minus_one*a_minus_one;

    }
}


void Bezier::derivative_of_tensor_product_bezier(
       NumMatrix tensor_product_coeffs,
       Common::Direction dx_i,
       size_t coeff_dim,
       size_t basis_bidegree, NumMatrix& derivative_coeffs){
    assert(tensor_product_coeffs.m() == coeff_dim);
    assert(tensor_product_coeffs.n() == pow(basis_bidegree+1, 2));
    setvalue(derivative_coeffs, 0.); 
    //NumMatrix derivative_coeffs(tensor_product_coeffs.m(), tensor_product_coeffs.n());
    assert(derivative_coeffs.m() == tensor_product_coeffs.m());
    assert(derivative_coeffs.n() == tensor_product_coeffs.n());
    int n = basis_bidegree; // I'm lazy
    // base cases i = 0; j=0; i = basis_bidegree; j= basis_bidegree;
    // just differentiate the original bezier curve with no elevation i.e. 
    // C^n_{i,j} = n*(a_i+1 - a_i) or something to this effect
    for(size_t k = 0; k <= n; k++){
        // j= 0
        if(dx_i == Common::V){

            set_value_edge_case(n, k, 0, coeff_dim, dx_i, 
                    tensor_product_coeffs, derivative_coeffs);
            // j = n
            set_value_edge_case(n, k, n, coeff_dim, dx_i, 
                    tensor_product_coeffs, derivative_coeffs);
        }
        if(dx_i == Common::U){
            // i = 0
            set_value_edge_case(n, 0, k, coeff_dim, dx_i, 
                    tensor_product_coeffs, derivative_coeffs);
            // i = n
            set_value_edge_case(n, n, k, coeff_dim, dx_i, 
                    tensor_product_coeffs, derivative_coeffs);
        }

    }
    if(dx_i == Common::U){
        for(size_t i = 1; i < n; i++){
            for(size_t j = 0; j <= n; j++){
                set_value(n, i, j, coeff_dim, dx_i, 
                        tensor_product_coeffs, derivative_coeffs);
            }
        }
    }
    if(dx_i == Common::V){
        for(size_t i = 0; i <= n; i++){
            for(size_t j = 1; j < n; j++){
                set_value(n, i, j, coeff_dim, dx_i, 
                        tensor_product_coeffs, derivative_coeffs);
            }
        }
    }
}

// ---------------------------------------------------------------------------
Point3 Bezier::de_casteljau_recursive(vector<Point3> control_points, 
        double t){
    // Note: assuming 1d polynomials only
    int poly_degree = control_points.size()-1;

    vector<double> blossom_eval_points(poly_degree, t);
    return blossom_recursive(control_points, blossom_eval_points, poly_degree, 0);
}

Point3 Bezier::blossom_recursive(vector<Point3> control_points, 
        vector<double> t_values, 
        int r,
        int i){
        if(r == 0){
            return control_points[i];
        } else {
            double ti = t_values[t_values.size()-1];
            vector<double> remaining_t_values(t_values.size()-1);
            for(int k =0; k < t_values.size()-1; k++){
                remaining_t_values[k] = t_values[k];
            }

            return (1. - ti)* blossom_recursive(
                                    control_points,
                                    remaining_t_values, 
                                    r-1,
                                    i+1) 
                         + ti*blossom_recursive(
                                  control_points,
                                  remaining_t_values,
                                  r-1,
                                  i);
        }
}

Point3 Bezier::de_casteljau(vector<Point3> control_points, 
        double t){
    // Note: assuming 1d polynomials only
    int poly_degree = control_points.size()-1;

    vector<double> blossom_eval_points(poly_degree, t);
    return blossom(control_points, blossom_eval_points, poly_degree, 0);
}

Point3 Bezier::blossom(vector<Point3> control_points, 
        vector<double> t_values, 
        int r,
        int i){

    int n = r;
    int stack_size = n+1;
    assert(control_points.size() == stack_size);
    assert(t_values.size() == n);
    vector<Point3> stack = control_points;

    int k = 0;
    //reverse(t_values.begin(), t_values.end()); 
    while( k<= n){
        stack_size -= 1;
        for(int j = 0; j < stack_size; j++){

            //int index = stack_size-j-1;
            int index = stack_size-1;
            assert(index >= 0);
            double ti = t_values[index];
            stack[j] = ti*stack[j] + (1. - ti)*stack[j+1];
        }
        k += 1;

        // contract the size of the stack after the reduction

    }
    // all done
    return stack[0];

}


Point3 Bezier::tensor_product_blossom_recursive(vector<Point3> control_points, vector<Point2> t){

    int tensor_product_degree = control_points.size();
    int poly_degree = int(sqrt(tensor_product_degree)) - 1;
    // samples are (x_i, y_j)
    // coefficients are x-major, then yy
    // first, do a blossom evaluation on each of the n lines of "y" control points, and
    // store the result
    //assert(t.size() == 2*poly_degree);

    vector<Point3> second_blossom_ctrl_points_in_x_coords(poly_degree+1, Point3(0.));;
    vector<double> t_values_along_line(poly_degree, 0.);

    for(int i = 0; i < poly_degree; i++){
        Point2 t_value  = t[i];
        t_values_along_line[i] = t_value.y();
    }

    for(int di = 0; di <= poly_degree; di++){
        vector<Point3> first_blossom_ctrl_points_in_y_coords(poly_degree+1, Point3(0.));

        // yank out control points along a line x=c in parameteric space
        for(int dj = 0; dj <= poly_degree; dj++){
            Point3 control_point_y = control_points[di*(poly_degree+1)+dj];
            first_blossom_ctrl_points_in_y_coords[dj] = control_point_y;
        }
        
        // do a 1d blossom along these points, and save them
        Point3 x_control_point = blossom_recursive(first_blossom_ctrl_points_in_y_coords,
                                        t_values_along_line, poly_degree, 0);
        second_blossom_ctrl_points_in_x_coords[di] = x_control_point;
    }
    for(int i = 0; i < poly_degree; i++){
        Point2 t_value  = t[i];
        t_values_along_line[i] = t_value.x();
    }
    return blossom_recursive(second_blossom_ctrl_points_in_x_coords, t_values_along_line, poly_degree, 0);

}



Point3 Bezier::tensor_product_blossom(vector<Point3> control_points, vector<Point2> t){

    int tensor_product_degree = control_points.size();
    int poly_degree = int(sqrt(tensor_product_degree)) - 1;
    // samples are (x_i, y_j)
    // coefficients are x-major, then yy
    // first, do a blossom evaluation on each of the n lines of "y" control points, and
    // store the result
    //assert(t.size() == 2*poly_degree);

    vector<Point3> second_blossom_ctrl_points_in_x_coords(poly_degree+1, Point3(0.));;
    vector<double> t_values_along_line(poly_degree, 0.);

    for(int i = 0; i < poly_degree; i++){
        Point2 t_value  = t[i];
        t_values_along_line[i] = t_value.y();
    }

    for(int di = 0; di <= poly_degree; di++){
        vector<Point3> first_blossom_ctrl_points_in_y_coords(poly_degree+1, Point3(0.));

        // yank out control points along a line x=c in parameteric space
        for(int dj = 0; dj <= poly_degree; dj++){
            Point3 control_point_y = control_points[di*(poly_degree+1)+dj];
            first_blossom_ctrl_points_in_y_coords[dj] = control_point_y;
        }
        
        // do a 1d blossom along these points, and save them
        Point3 x_control_point = blossom(first_blossom_ctrl_points_in_y_coords,
                                        t_values_along_line, poly_degree, 0);
        second_blossom_ctrl_points_in_x_coords[di] = x_control_point;
    }
    for(int i = 0; i < poly_degree; i++){
        Point2 t_value  = t[i];
        t_values_along_line[i] = t_value.x();
    }
    return blossom(second_blossom_ctrl_points_in_x_coords, t_values_along_line, poly_degree, 0);

}

vector<Point3> Bezier::compute_control_points_on_subdomain(
        vector<Point3> control_points,
        int patch_order,
        pair<double, double> x_interval,
        pair<double, double> y_interval){

    vector<Point3> subdomain_control_points((patch_order+1)*(patch_order+1));
    for(int di = 0; di <= patch_order; di++){
        for(int dj = 0; dj <= patch_order; dj++){
            // set everything to (x_0, y_0)
            vector<Point2> t_values(patch_order, 
                    Point2(x_interval.first, y_interval.first));

            // fill in the first i values with (x_1, y_1) 
            // (implies the remaining n-i values are (x_0, y_0)'s)
            // Recall order doesn't matter by symmetry
            for(int i = 0; i < patch_order-di; i++){
                t_values[i].x() = x_interval.second;
            }
            for(int i = 0; i < patch_order-dj; i++){
                t_values[i].y() = y_interval.second;
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
    return subdomain_control_points;
}

// ---------------------------------------------------------------------------




/*
Point3 Bezier::de_casteljau(vector<Point3> control_points, 
        double t){
    // Note: assuming 1d polynomials only
    int poly_degree = control_points.size()-1;

    vector<double> blossom_eval_points(poly_degree, t);
    return blossom(control_points, blossom_eval_points, poly_degree, 0);
}

Point3 Bezier::blossom(vector<Point3> control_points, 
        vector<double> t_values, 
        int r,
        int i){
        if(r == 0){
            return control_points[i];
        } else {
            double ti = t_values[t_values.size()-1];
            vector<double> remaining_t_values(t_values.size()-1);
            for(int k =0; k < t_values.size()-1; k++){
                remaining_t_values[k] = t_values[k];
            }

            return (1. - ti)* blossom(
                                    control_points,
                                    remaining_t_values, 
                                    r-1,
                                    i+1) 
                         + ti*blossom(
                                  control_points,
                                  remaining_t_values,
                                  r-1,
                                  i);
        }
}

Point3 Bezier::tensor_product_blossom(vector<Point3> control_points, vector<Point2> t){

    int tensor_product_degree = control_points.size();
    int poly_degree = int(sqrt(tensor_product_degree)) - 1;
    // samples are (x_i, y_j)
    // coefficients are x-major, then yy
    // first, do a blossom evaluation on each of the n lines of "y" control points, and
    // store the result
    //assert(t.size() == 2*poly_degree);

    vector<Point3> second_blossom_ctrl_points_in_x_coords(poly_degree+1, Point3(0.));;
    vector<double> t_values_along_line(poly_degree, 0.);

    for(int i = 0; i < poly_degree; i++){
        Point2 t_value  = t[i];
        t_values_along_line[i] = t_value.y();
    }

    for(int di = 0; di <= poly_degree; di++){
        vector<Point3> first_blossom_ctrl_points_in_y_coords(poly_degree+1, Point3(0.));

        // yank out control points along a line x=c in parameteric space
        for(int dj = 0; dj <= poly_degree; dj++){
            Point3 control_point_y = control_points[di*(poly_degree+1)+dj];
            first_blossom_ctrl_points_in_y_coords[dj] = control_point_y;
        }
        
        // do a 1d blossom along these points, and save them
        Point3 x_control_point = blossom(first_blossom_ctrl_points_in_y_coords,
                                        t_values_along_line, poly_degree, 0);
        second_blossom_ctrl_points_in_x_coords[di] = x_control_point;
    }
    for(int i = 0; i < poly_degree; i++){
        Point2 t_value  = t[i];
        t_values_along_line[i] = t_value.x();
    }
    return blossom(second_blossom_ctrl_points_in_x_coords, t_values_along_line, poly_degree, 0);

}

vector<Point3> Bezier::compute_control_points_on_subdomain(
        vector<Point3> control_points,
        pair<double, double> x_interval,
        pair<double, double> y_interval){

    vector<Point3> subdomain_control_points((_patch_order+1)*(_patch_order+1));
    for(int di = 0; di <= _patch_order; di++){
        for(int dj = 0; dj <= _patch_order; dj++){
            // set everything to (x_0, y_0)
            vector<Point2> t_values(_patch_order, 
                    Point2(x_interval.first, y_interval.first));

            // fill in the first i values with (x_1, y_1) 
            // (implies the remaining n-i values are (x_0, y_0)'s)
            // Recall order doesn't matter by symmetry
            for(int i = 0; i < _patch_order-di; i++){
                t_values[i].x() = x_interval.second;
            }
            for(int i = 0; i < _patch_order-dj; i++){
                t_values[i].y() = y_interval.second;
            }

            //Evaluate the blossom and save the control point
            Point3 control_point_on_subdomain = 
                Bezier::tensor_product_blossom(
                        control_points, 
                        t_values);

            int index = di*(_patch_order+1)+dj;
            subdomain_control_points[index] = control_point_on_subdomain;

        }
    }
    return subdomain_control_points;
}

*/
