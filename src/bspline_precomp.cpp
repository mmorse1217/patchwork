#include "bspline_precomp.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "bspline_codegen.hpp"
using namespace std;

int find_knot_interval(double x, int n, vector<double> knots){
    vector<double> knots_unique(knots);
    int interval = -1;
    knots.erase( std::unique( knots.begin(), knots.end() ), knots.end() );
    for(int i = 0; i < knots.size()-1; i++){
        if(x <= knots[i+1] && knots[i] <= x && fabs(knots[i] -knots[i+1]) > 1e-14){
            interval = i;
            //return interval;
        }
    }
    //return interval != -1 ? interval - n : interval;
    return interval;
}


// NOTE to future poor soul: values here were computed with 
// sympy.functions.special.bsplines.bspline_basis(5, range(7), 0, x)
//
double BSpline::bspline_value(double x, int i, int n, vector<double> knots){
    assert(n == 5 || n == 3);
    x -= i;
    switch(n){

        case 3:
            return bspline_deg3_deriv0(x);
        case 5:
            return bspline_deg5_deriv0(x);
        default:
            assert(0);

    }
}

double BSpline::bspline_first_deriv(double x , int i, int n, vector<double>  knots){
    //assert(n == 5);
    assert(n == 5 || n == 3);
    int interval = find_knot_interval(x, n,  knots);
    assert(interval != -1);
    x -= i;
    switch(n){

        case 3:
            return bspline_deg3_deriv1(x);
        case 5:
            return bspline_deg5_deriv1(x);
        default:
            assert(0);

    }
}

double BSpline::bspline_second_deriv(double x, int i, int n, vector<double>  knots){
    //assert(n == 5);
    assert(n == 5 || n == 3);
    int interval = find_knot_interval(x, n,  knots);
    assert(interval != -1);
    x -= i;
    switch(n){

        case 3:
            return bspline_deg3_deriv2(x);
        case 5:
            return bspline_deg5_deriv2(x);
        default:
            assert(0);

    }
}

double BSpline::bspline_third_deriv(double x, int n, vector<double>  knots){
    assert(n == 5);
    int interval = find_knot_interval(x, n,  knots);
    assert(interval != -1);
    double value = 0.;
    
    double x2 = pow(x, 2);

    switch(interval){
        case 0: // x \in [0,1]
            value = 1./2.*x2;
            break; 
        case 1:// x \in [1,2]
            value = -5./2.*x2 + 6.*x - 3.;
            break; 
        case 2:// x \in [2,3]
            value = 5.*x2 - 24.*x + 27.;
            break; 
        case 3:// x \in [3,4]
            value = -5.*x2 + 36.*x - 63.;
            break; 
        case 4:// x \in [4,5]
            value = 5./2.*x2 - 24.*x + 57.;
            break; 
        case 5:// x \in [5,6]
            value = -1./2.*x2 + 6.*x - 18.;
            break; 
        default: // interval > 6
            value = 0.;
            break; 
    }   
    
    return value;
}

double BSpline::bspline_fourth_deriv(double x, int n, vector<double>  knots){
    assert(n == 5);
    int interval = find_knot_interval(x, n,  knots);
    assert(interval != -1);
    double value = 0.;
    
    switch(interval){
        case 0: // x \in [0,1]
            value = x;
            break; 
        case 1:// x \in [1,2]
            value = -5.*x + 6.;
            break; 
        case 2:// x \in [2,3]
            value = 10.*x - 24.;
            break; 
        case 3:// x \in [3,4]
            value = -10.*x + 36.;
            break; 
        case 4:// x \in [4,5]
            value = 5.*x - 24.;
            break; 
        case 5:// x \in [5,6]
            value = -x + 6.;
            break; 
        default: // interval > 6
            value = 0.;
            break; 
    }   
    

    return value;
}
