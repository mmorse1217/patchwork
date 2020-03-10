#ifndef __BSPLINE_PRECOMP_HPP__
#define __BSPLINE_PRECOMP_HPP__

#include <vector>
#include "assert.h"
using namespace std;
namespace BSpline{
    enum EvalType{
        VALUE = 0,
        DERIV_1ST = 1,
        DERIV_2ND = 2,
        DERIV_3RD = 3,
        DERIV_4TH = 4
    };

    double bspline_value(double x, int i, int n, vector<double> knots);
    double bspline_first_deriv(double x ,int i,  int n, vector<double>  knots);
    double bspline_second_deriv(double x, int i, int n, vector<double>  knots);
    double bspline_third_deriv(double x, int n, vector<double>  knots);
    double bspline_fourth_deriv(double x, int n, vector<double>  knots);
};
#endif
