#include "bspline.hpp"

using namespace BSpline;


BSpline::KnotVector::KnotVector(size_t spline_degree, size_t num_control_points, SplineType spline_type){
        int n = spline_degree;
    if(spline_type == PERIODIC){
        for(int i = -n; i < num_control_points + 2*n - 1; i++){
            this->push_back(pair<double, int>(double(i), 1));
        }
    } else if (spline_type == OPEN){
        this->resize(num_control_points+n+1);
        for(int i = -n; i < num_control_points  + 1; i++){
            this->at(i) = pair<double, int>(double(i), 1);
        }
        
        /*
        //for(int i = 0; i < n; i++)
        this->push_back(pair<double, int>(double(0),n));

        for(int i = 0; i < num_control_points - n + 1; i++)
            this->push_back(pair<double, int>(double(i), 1));

        //for(int i = 0; i < n; i++)
        this->push_back(pair<double, int>(double(num_control_points- n), n));*/
    } else {
        assert(0);
    }





    _flattened_knots  = generate_bspline_knots(spline_degree, num_control_points, spline_type);
}

BSpline::KnotVector::KnotVector(vector<double> knots, vector<int> multiplicities){
    assert(knots.size() == multiplicities.size());
    int n = knots.size();
    for(int i =0; i < n; i++){
        this->push_back(pair<double, int>(knots[i], multiplicities[i]));
    }
    _flattened_knots = this->flatten();
}

int BSpline::KnotVector::knot_multiplicity(int index){
    return this->at(index).second;
}
double BSpline::KnotVector::knot_value(int index){
    return this->at(index).first;
}

void BSpline::KnotVector::set_knot_multiplicity(int index, int multiplicty){
    this->at(index).second = multiplicty;
}
vector<double> BSpline::KnotVector::flatten(){
    vector<double> flattened_knots;
    
    for(int i = 0; i < this->size(); i++){
        
        double  knot_value = this->knot_value(i);
        int multiplicty = this->knot_multiplicity(i);
        
        for(int mi = 0; mi < multiplicty; mi++){
            flattened_knots.push_back(knot_value);
        }
    }
    
    return flattened_knots;
}



vector<double> BSpline::generate_bspline_knots(int n, int num_control_points, 
        SplineType spline_type){
    vector<double> knots;

    if(spline_type == PERIODIC){
        assert(0); // WARNING this probably doesn't work
        for(int i = -n; i < num_control_points + 2*n - 1; i++){
            knots.push_back(i);
        }
    } else if (spline_type == OPEN){
        for(int i = -n; i < num_control_points  + 1; i++){
            knots.push_back(i);
        }
        /*
        for(int i = 0; i < n; i++)
            knots.push_back(0);

        for(int i = 0; i < num_control_points - n + 1; i++)
            knots.push_back(i);

        for(int i = 0; i < n; i++)
            knots.push_back(num_control_points- n);
        */
    } else {
        assert(0);
    }

    return knots;
}

double BSpline::rescale_to_bspline_interval(double x, int n, int num_control_points,
        SplineType spline_type){

    // TODO double check that this work or periodic splines
        return (num_control_points-n)*x;
        /*
    if(spline_type == OPEN){
        return (num_control_points-n)*x;
        //return x;
    } else if (spline_type == PERIODIC){

        // yanked from
        // https://stackoverflow.com/questions/34803197/fast-b-spline-algorithm-with-numpy-scipy/35007804
        // this is some scary stuff. not really sure why it works
        double result = double(num_control_points)*x - .5*double(n-1); 
        return result < 0 ? result + num_control_points : result;
    } else {
        assert(0);
    }*/
}

double BSpline::de_boor(int i, int k, double x, vector<double> knots){
    /**
     * base case
     * k = 0; if t \in [u_i, u_{i+1}] (the ith knot span), the value of the basis
     * function is 1, otherwise it's zero. The ith basis function is also only
     * non-zero on [u_i, u_{i+n+1}], so we'll only look at these knot spans
     */
    if( k == 0){
        if(x <= knots[i+1] && knots[i] <= x){
            return 1.;
        } else {
            return 0.;
        }
    } else {
        // compute interpolation ratio denominators
        double first_diff = knots[i+k+1] - knots[i+1];
        double  second_diff = knots[i+k] - knots[i];
        
        // check for division by zero; if there is one, zero out that term
        double A = fabs(first_diff) > 0 ? (knots[i+k+1] - x)/first_diff : 0;
        double B = fabs(second_diff) > 0 ? (x - knots[i])/second_diff : 0;

        // recursively evaluate k-1 B-Splines and linearly inteporlate between
        // the two.
        return A*de_boor(i+1, k-1, x, knots) + B*de_boor(i, k-1, x, knots);
    }
}


double BSpline::de_boor_iterative(int i, int n, double x, vector<double> knots){
    /**
     * base case
     * k = 0; if t \in [u_i, u_{i+1}] (the ith knot span), the value of the basis
     * function is 1, otherwise it's zero. The ith basis function is also only
     * non-zero on [u_i, u_{i+n+1}], so we'll only look at these knot spans
     */
    int stack_size = n+1;
    vector<double> stack(stack_size, 0.);

    for(int k = 0; k < stack_size; k++){
        // if the evaluation target is inside the knot span of the current basis
        // function...
        if(x <= knots[i+k+1] &&  knots[i+k] <= x){
            // give it a weight of 1.
            stack[k] = 1.;
        }
        // otherwise, it remains 0. as per default initialization in constructor
    }

    /* the first step of the while loop is the first step up from the base case
     * of the recursion, so we will linearly interpolate between the previously
     * computed basis function values
     */
    int k = 1;

    /* stack size is proportional to the number of contributing k-degree basis
     * functions to the value of the k+1 degree basis function. Hence, if the
     * stack size is one, there's no more values to interpolate between, so we're
     * finished.
     *
     * We'll maintain the "actual" size of the stack we need to look at
     * ourselves to avoid reallocating + copying or vector resizing at each
     * iteration.
     */
    while( k<= n){
        for(int j = 0; j < stack_size-1; j++){
            /*  j is an index relative to the compact support of the kth degree basis
             *  function (j=0,..., k), which is stored in the stack. To index into
             *  the appropriate piece of the knot vector, we shift the index by
             *  the basis function index.
             */
            int jj = j+i;

            // compute interpolation ratio denominators
            double first_diff = knots[jj+k+1] - knots[jj+1];
            double  second_diff = knots[jj+k] - knots[jj];

            // check for division by zero; if there is one, zero out that term
            double A = fabs(first_diff) > 0 ? (knots[jj+k+1] - x)/first_diff : 0;
            double B = fabs(second_diff) > 0 ? (x - knots[jj])/second_diff : 0;

            // linear interpolation, de Boor-style
            stack[j] = A*stack[j+1] + B*stack[j];
        }
        // repeat, interpolating the degree k+1 basis functions
        k += 1;

        // contract the size of the stack after the reduction
        stack_size -= 1;

    }
    // all done
    return stack[0];
}



double BSpline::evaluate_bspline(int i, int n, int num_control_points, double x, SplineType spline_type, EvalType eval_type, vector<double> knots){
    //if(i == num_control_points)
        //return 0.;
    double x_rescaled = BSpline::rescale_to_bspline_interval(x, n, num_control_points, spline_type);
    //cout << "unscaled: " << x << "; scaled x: " << x_rescaled << endl;
    //double b_spline_value = BSpline::de_boor(i, n, x_rescaled, knots);
    //double b_spline_value = BSpline::de_boor_iterative(i, n, x_rescaled, knots);
    double b_spline_value = BSpline::bspline_precomp(i, n, x_rescaled, knots, eval_type);

    if(std::isnan(b_spline_value)){
        return 0.;
    } else {
        return b_spline_value;
    }
}


double BSpline::bspline_precomp(int i, int n, double x, vector<double> knots, EvalType eval_type){
    switch(eval_type){
        case VALUE:
            return bspline_value( x ,  i, n, knots);
            break;
        case DERIV_1ST:
            return bspline_first_deriv( x , i, n, knots);
            break;
        case DERIV_2ND:
            return bspline_second_deriv( x , i, n, knots);
            break;
        case DERIV_3RD:
            return bspline_third_deriv( x ,  n, knots);
            break;
        case DERIV_4TH:
            return bspline_fourth_deriv( x ,  n, knots);
            break;
        default:
            assert(0);
    }

}
/*
double BSpline::evaluate_bspline_degree_elevated(int i, int n, int num_control_points, double x, SplineType spline_type, vector<double> knots){
    if(i == num_control_points)
        return 0.;
    
    double b_spline_value = 0.;
    double x_rescaled = BSpline::rescale_to_bspline_interval(x, n, num_control_points, spline_type);
    for(int j = 0; j < n; j++){
        vector<double> elevated_knots(knots.size(), 0.);
        
    }
    double b_spline_value = BSpline::de_boor(i, n, x_rescaled, knots);

    if(std::isnan(b_spline_value) | b_spline_value >=1e2){
        return 0.;
    } else {
        return b_spline_value;
    }

}*/

double BSpline::evaluate_bspline(int i, int n, int num_control_points, double x, SplineType spline_type){
    if(i == num_control_points)
        return 0.;
    
    vector<double> knots = BSpline::generate_bspline_knots(n, num_control_points, spline_type);
    
    double x_rescaled = BSpline::rescale_to_bspline_interval(x, n, num_control_points, spline_type);
    double b_spline_value = BSpline::de_boor(i, n, x_rescaled, knots);

    if(std::isnan(b_spline_value) || b_spline_value >=1e2){
        return 0.;
    } else {
        return b_spline_value;
    }
}

