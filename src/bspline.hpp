#ifndef __BSPLINE_HPP__
#define __BSPLINE_HPP__
#include "nummatrix.hpp"
#include "common.hpp"
#include "bspline_precomp.hpp"


namespace BSpline {

    enum SplineType {
        OPEN = 0,
        PERIODIC = 1
    };

    class KnotVector: public vector<pair<double, int> >{
        vector<double> _flattened_knots;
        public:
        KnotVector(): vector<pair<double, int> >(){;}
        KnotVector(vector<double> knots, vector<int> multiplicities);
        KnotVector(size_t spline_degree, size_t num_control_points, SplineType spline_type);

        void set_knot_multiplicity(int index, int multiplicty);
        int knot_multiplicity(int index);
        double knot_value(int index);

        double operator[](int i){
            assert(_flattened_knots.size() > 0);
            return _flattened_knots[i];
        }
        size_t knot_count(){
            return _flattened_knots.size();
        }
        vector<double> flatten();
        vector<double> flattened_knots(){return _flattened_knots;}

    };
    
    double evaluate_bspline(int i, int n,  
            int num_control_points, double x, SplineType spline_type=OPEN);
    double evaluate_bspline(int i, int n, int num_control_points, double x, SplineType spline_type, EvalType eval_type, vector<double> knots);

    vector<double> generate_bspline_knots(int n, int num_control_points, 
            SplineType spline_type);

    double rescale_to_bspline_interval(double x, int n, int num_control_points,
            SplineType spline_type);
    /**
     * Recursive implementation of de Boor's algorithm. 
     * TODO elaborate
     *
     * @param int               i           which B-spline B_k(x - i) to evalaute
     * @param int               k           degree of B-spline
     * @param double            x           evaluation point
     * @param vector<double>    knots       list of knot values for the B-spline
     *
     */
    double de_boor(int i, int k, double x, vector<double> knots);
    /**
     * Iterative implementation of de Boor's algorithm. 
     * TODO elaborate
     *
     * @param int               i           which B-spline B_k(x - i) to evalaute
     * @param int               k           degree of B-spline
     * @param double            x           evaluation point
     * @param vector<double>    knots       list of knot values for the B-spline
     *
     */
    double de_boor_iterative(int i, int k, double x, vector<double> knots);
    
    double bspline_precomp(int i, int n, double x, vector<double> knots, EvalType eval_type);

};

#endif
