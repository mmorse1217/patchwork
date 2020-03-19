#ifndef __FUNCTION_HPP__
#define __FUNCTION_HPP__
#include <numvector.hpp>
#include <nummatrix.hpp>
//#include "face_map.hpp"
#include "bezier.hpp"
#include "bspline.hpp"
#include "common.hpp"
#include "vec2t.hpp"

class Function {
    
    class BasisFunction {
        protected:
            size_t _i;
            size_t _basis_degree;
        public:
        BSpline::EvalType _eval_type;
            BasisFunction(size_t i, size_t basis_degree):
                _i(i), _basis_degree(basis_degree),
            _eval_type(BSpline::VALUE){;}
            virtual double evaluate(double x) = 0;
    };

    class BezierBasisFunction: public BasisFunction{
        public:
        IntMatrix* _binomial_coeffs;
        BezierBasisFunction(size_t i, 
                size_t basis_degree, 
                IntMatrix* binomial_coeffs):
            BasisFunction(i, basis_degree), 
            _binomial_coeffs(binomial_coeffs) {;}
        double evaluate(double x);
    };

    class BSplineBasisFunction: public BasisFunction{
        size_t _num_basis_functions;
        public:
        BSpline::KnotVector _knots;
        BSplineBasisFunction(size_t i,
                size_t basis_degree, 
                size_t num_basis_functions, 
                BSpline::KnotVector knots):

            BasisFunction(i, basis_degree), 
            _num_basis_functions(num_basis_functions), 
            _knots(knots)
            {;}
        double evaluate(double x);
    };

    // Generic Function class for F: R^n \to R^m 
    protected:
    size_t _domain_dim; // n
    size_t _range_dim;  // m
    size_t _basis_bidegree; // polynomial basis order
    size_t _num_basis_functions; // actual number of basis functions per dimension to evaluate
    // for Bezier, this is _basis_bidegree +1,
    // for B-Splines, this is # of knots
    // points == # num ctrl points + spline degree + 1
    //IntVector _basis_degrees;
    Common::BasisType _basis_type;
    bool _cached_matrices;
    NumMatrix* _pseudo_inverse;
    NumMatrix* _values_at_basis_functions;
    BSpline::KnotVector _knots;
    BSpline::SplineType _spline_type;


        NumVector (*_eval)(NumVector);
        
        
    public:
        // _range_dim x (basis_bidegree)^_domain_dim matrix of polynomial
        // coefficients. 
        // Each column is a coefficient for the corresponding
        // tensor product basis function.
        NumMatrix _polynomial_coefficents;
        IntMatrix _binomial_coefficients_array;
        vector<BasisFunction*> _basis_functions;
        Common::Periodicity _periodic;
        double _period_width;
        // temp fix until we can match poly degree for splines
        vector<BasisFunction*> _basis_functions_x;
        vector<BasisFunction*> _basis_functions_y;
        // to do load options from options file... or set_from_opts function
        Function();
        Function(size_t domain_dim, size_t range_dim, 
                size_t basis_bidegree, Common::BasisType basis_type);
        Function(size_t basis_bidegree, Common::BasisType basis_type);
        Function(size_t domain_dim, size_t range_dim, 
                size_t basis_bidegree, size_t num_control_points,
                Common::BasisType basis_type);
        ~Function();
        //Function(const Function& func);
        vector<BasisFunction*> generate_basis_set(size_t basis_bidegree, 
                size_t num_basis_functions, 
                Common::BasisType basis_type);

        Function& operator=(const Function& func);
        NumVector operator()(NumVector x);

        int domain_dim(){
            return _domain_dim;
        }

        int range_dim(){
            return _range_dim;
        }
        /*void set_polynomial_basis(BasisType basis_type, size_t basis_bidegree){
            _basis_bidegree = basis_bidegree;
            _basis_type = basis_type;
            _binomial_coefficients = 
                FaceMapSurf::compute_binomial_coefs(_basis_bidegree);
        }*/

        /**
         * Construct a functional representation of data given a matrix of
         * coordinates in R^n and values in R^m. The Function is represented as
         * an n-fold polynomial tensor product polynomial approximation to the
         * input data. Only tested for n=2.
         * @param NumMatrix     coordinates         n x (num_points) matrix of
         *                                          coordinates. the ith column of
         *                                          "coordinates" is a member of
         *                                          R^n and is mapped to the ith
         *                                          column of "values" under the
         *                                          action of F
         * @param NumMatrix     values              m x (num_points) matrix of
         *                                          function values. the ith column of
         *                                          "values" is a member of
         *                                          R^m and is the result of
         *                                          apply the target function F
         *                                          to the ith column of 
         *                                          "coordinates".
         * @param BasisType     basis_type          basis function type
         * @param size_t        basis_bidegree      bidegree of basis functions
         * @return Function                         a continuous functional 
         *                                          representation of the input
         *                                          data
         */
        Function static construct_function(NumMatrix coordinates, 
                NumMatrix values, Common::BasisType basis_type, size_t basis_order);
        
        Function static construct_spline_function(NumMatrix coordinates, NumMatrix values, 
                size_t num_control_points, size_t basis_order);
        /**
         * Differentiate func with respect to the variable x_i. This ultimately
         * just computes the derivative of the underlying polynomial approximation.
         * It's worth noting that high derivatives will likely give the user
         * trouble if they aren't careful.
         * @param size_t        i                   Which variable we should be
         *                                          differentiating with respect
         *                                          to (i.e. d/dx_i( * ) 
         * @param Fucntion      func                function we'd like the
         *                                          derivative of
         * @return Function                         the derivative of func 
         *                                          w.r.t x_i
         */
        Function static d_dx(size_t i, Function func);

        
        // TODO REFACTOR and optionally cache basis function values at sample
        // points and the least squares matrix to accelerate
        NumMatrix compute_polynomial_approximation(NumMatrix coordinates, 
                NumMatrix values);

        NumMatrix static compute_polynomial_approximation(NumMatrix coordinates, 
                NumMatrix values, Common::BasisType basis_type, size_t basis_order);
        

        double basis_function(IntVector basis_index, NumVector x);
        
        NumMatrix static compute_pseudo_inverse(
                int basis_bidegree,
                IntMatrix binomial_coefficients,
                Common::BasisType basis_type,
                NumMatrix coordinates){
            
            return compute_pseudo_inverse(basis_bidegree, 
                    basis_bidegree+1, 
                    binomial_coefficients, 
                    basis_type, 
                    coordinates);
        }

        NumMatrix compute_pseudo_inverse(
                vector<BasisFunction*> basis_functions,
                NumMatrix coordinates);

        NumMatrix static compute_pseudo_inverse(
                vector<BasisFunction*> basis_functions_x,
                vector<BasisFunction*> basis_functions_y,
                NumMatrix coordinates);


        NumMatrix static compute_pseudo_inverse(
                int basis_bidegree,
                int num_basis_functions,
                IntMatrix binomial_coefficients,
                Common::BasisType basis_type,
                NumMatrix coordinates);

        NumMatrix static compute_basis_function_values(
                size_t domain_dim,
                vector<BasisFunction*>& basis_functions_x,
                vector<BasisFunction*>& basis_functions_y,
                NumVector x);
        
        NumMatrix static compute_basis_function_values(
                size_t domain_dim,
                vector<BasisFunction*> basis_functions,
                NumVector x);


        NumMatrix static compute_basis_function_values(
                size_t domain_dim,
                size_t basis_bidegree,
                size_t num_basis_functions,
                IntMatrix binomial_coefficients,
                Common::BasisType basis_type,
                NumVector x);


        // cached version of above
        NumMatrix compute_basis_function_values(NumVector x);
        NumMatrix compute_pseudo_inverse(
                NumMatrix coordinates);
                
        
        void set_cache( NumMatrix* pseudo_inverse, NumMatrix* basis_function_values_1d){
            _pseudo_inverse = pseudo_inverse;
            _values_at_basis_functions = basis_function_values_1d;
            _cached_matrices = true;
        }


        void set_periodicity(Common::Periodicity periodic_type, double periodic_width, bool unshift=false);

};
        long long double_to_index(double x);
#endif
