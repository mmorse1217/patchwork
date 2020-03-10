#include "function.hpp"
#include "face_map.hpp"
#include "bezier.hpp"
#include "bspline.hpp"

using Common::BasisType;
using Common::BEZIER;
using Common::BSPLINE;
using BSpline::KnotVector;
using BSpline::VALUE;
using BSpline::DERIV_1ST;
using BSpline::DERIV_2ND;
using BSpline::DERIV_3RD;
using BSpline::DERIV_4TH;

Function::Function():
    _domain_dim(0), 
    _range_dim(0), 
    _basis_bidegree(0), 
    _num_basis_functions(0),
    _basis_type(Common::BEZIER),
    _cached_matrices(false),
    _pseudo_inverse(NULL),
    _knots(KnotVector(_basis_bidegree, _num_basis_functions, BSpline::OPEN))
{

}
Function::Function(size_t basis_bidegree, Common::BasisType basis_type): 
    _domain_dim(2), 
    _range_dim(3), 
    _basis_bidegree(basis_bidegree), 
    _basis_type(basis_type),
    _cached_matrices(false),
    _pseudo_inverse(NULL)
    //_knots(KnotVector(_basis_bidegree, _num_basis_functions, BSpline::OPEN))
{
    if(basis_type == Common::BEZIER){
        _num_basis_functions = _basis_bidegree + 1;
    } else if(basis_type == Common::BSPLINE){
        // assumes _basis_bidegree == num control points for open,
        //                         == num_control_points + spline deg for
        //                              periodic
        _num_basis_functions = _basis_bidegree;
        _knots = KnotVector(_basis_bidegree, _num_basis_functions, BSpline::OPEN);
    }
    _binomial_coefficients_array = 
        *FaceMapSurf::compute_binomial_coefs_array(_basis_bidegree);
    _basis_functions = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
        _basis_functions_x = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
        _basis_functions_y = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);

}

Function::Function(size_t domain_dim, size_t range_dim, 
        size_t basis_bidegree, Common::BasisType basis_type): 
    _domain_dim(domain_dim), 
    _range_dim(range_dim), 
    _basis_bidegree(basis_bidegree), 
    _basis_type(basis_type),
    _cached_matrices(false),
    _pseudo_inverse(NULL)
{
    if(basis_type == Common::BEZIER){
        _num_basis_functions = _basis_bidegree + 1;
    } else if(basis_type == Common::BSPLINE){
        // assumes _basis_bidegree == num control points for open,
        //                         == num_control_points + spline deg for
        //                              periodic
        _num_basis_functions = _basis_bidegree;
        _knots = KnotVector(_basis_bidegree, _num_basis_functions, BSpline::OPEN);
    }

    _binomial_coefficients_array = 
        *FaceMapSurf::compute_binomial_coefs_array(_basis_bidegree);
    //if(domain_dim == 1){
    _basis_functions = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
    /*} else if(domain_dim == 2){

        _basis_functions_x = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
        _basis_functions_y = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
    }*/
        _basis_functions_x = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
        _basis_functions_y = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
}

Function::Function(size_t domain_dim, size_t range_dim, 
        size_t basis_bidegree, size_t num_control_points, Common::BasisType basis_type): 
    _domain_dim(domain_dim), 
    _range_dim(range_dim), 
    _basis_bidegree(basis_bidegree), 
    _num_basis_functions(num_control_points),
    _basis_type(basis_type),
    _cached_matrices(false),
    _pseudo_inverse(NULL)
    //_knots(KnotVector(_basis_bidegree, _num_basis_functions, BSpline::OPEN))
{

    _binomial_coefficients_array = 
        *FaceMapSurf::compute_binomial_coefs_array(_basis_bidegree);
    _knots = KnotVector(_basis_bidegree, _num_basis_functions, BSpline::OPEN);
    //if(domain_dim == 1){
    _basis_functions = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
    /*} else if(domain_dim == 2){

        _basis_functions_x = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
        _basis_functions_y = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
    }*/
    //generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);

        _basis_functions_x = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
        _basis_functions_y = generate_basis_set(_basis_bidegree, _num_basis_functions, basis_type);
}

vector<Function::BasisFunction*> Function::generate_basis_set(size_t basis_bidegree, 
        size_t num_basis_functions, 
        Common::BasisType basis_type){
    // sanity checks
    if(basis_type == Common::BEZIER){
        // need exactly bidegree + 1 polynomials
        assert(num_basis_functions == basis_bidegree+1);
    } else if(basis_type == Common::BSPLINE){
        // need at least bidegree + 1 control points (i.e. Bsplines) to
        // represent a (biegree) curve 
        assert(basis_bidegree < num_basis_functions);
        assert(_knots.knot_count() == basis_bidegree + num_basis_functions + 1);
    } else {
        assert(basis_type == Common::BEZIER || basis_type == Common::BSPLINE);
    }
    vector<BasisFunction*> basis_functions(num_basis_functions); 
    for(size_t i = 0; i < num_basis_functions; i++){
        if(basis_type == Common::BEZIER){
            basis_functions[i] = new BezierBasisFunction(i, basis_bidegree, &_binomial_coefficients_array);
        } else if(basis_type == Common::BSPLINE){
            basis_functions[i] = new BSplineBasisFunction(i, basis_bidegree, num_basis_functions, _knots);
        }

    }
    return basis_functions;
}

double Function::BSplineBasisFunction::evaluate(double x){
    return BSpline::evaluate_bspline(_i, _basis_degree, _num_basis_functions, x, BSpline::OPEN, _eval_type, _knots.flattened_knots());
}

double Function::BezierBasisFunction::evaluate(double x){
    // I'm a horrifying person.
    // TODO REMOVE DUPLICATE B() FUNCTION FOR POINTER AND NON-POINTER BINOMIAL
    // COEFFS. 
    int binomial_coeff = (*_binomial_coeffs)(_basis_degree, _i);
    return FaceMapSurf::B(_i, _basis_degree, x, binomial_coeff);
}


Function::~Function(){
    // add deletion of basis functions...
    //delete _binomial_coefficients;

}
/*
   Function::Function(const Function& func):
   _domain_dim(func._domain_dim), 
   _range_dim(func._range_dim), 
   _basis_bidegree(func._basis_bidegree), 
   _basis_type(func._basis_type),
   _polynomial_coefficents(func._polynomial_coefficents){

   }
   */
Function& Function::operator=(const Function& func){
    _domain_dim = func._domain_dim;
    _range_dim = func._range_dim;
    _basis_bidegree = func._basis_bidegree;
    _basis_type = func._basis_type;
    _num_basis_functions = func._num_basis_functions;

    // need it's own copy of binomial coeffs; if we're initializing the funciton
    // now, we should create new ones
    // using func._binomial_coefficients is made because func is a copy that
    // goes out-of-scope after this function call, so don't do that.
    //if(_binomial_coefficients != NULL){
    //delete _binomial_coefficients;
    //}
    //_binomial_coefficients = FaceMapSurf::compute_binomial_coefs_stack(_basis_bidegree);
    _binomial_coefficients_array = func._binomial_coefficients_array;
    _polynomial_coefficents = func._polynomial_coefficents;
    _pseudo_inverse = func._pseudo_inverse;
    _basis_functions = func._basis_functions;
    _basis_functions_x = func._basis_functions_x;
    _basis_functions_y = func._basis_functions_y;
    _knots = func._knots;
    if(_basis_type == Common::BEZIER){
        for(int bi = 0; bi < _num_basis_functions; bi++){
            ((BezierBasisFunction*)_basis_functions[bi])->_binomial_coeffs = &_binomial_coefficients_array;
        }
    }
    //generate_basis_set(_basis_bidegree, _num_basis_functions, _basis_type);
    //_basis_function_values= func._basis_function_values;
    return *this;
}

NumVector Function::operator()(NumVector x){
    assert(x.length() == _domain_dim);
    assert(_domain_dim == 1 || _domain_dim == 2);
    NumMatrix function_value(_range_dim, 1);
    
    NumMatrix values_of_basis_functions(pow(_num_basis_functions, _domain_dim),1);
    NumMatrix basis_function_values_1d = compute_basis_function_values(x);
    if (_domain_dim == 2){
        for(int j = 0; j < _num_basis_functions; j++){
            for(int k = 0; k < _num_basis_functions; k++){
                int index = j*(_num_basis_functions) + k;
                values_of_basis_functions(index,0) = basis_function_values_1d(0,j)*basis_function_values_1d(1,k);
            }
        } 
    } else if (_domain_dim == 1){
        for(int j = 0; j <_num_basis_functions; j++){
            values_of_basis_functions(j,0) = basis_function_values_1d(0,j);
        } 
    }
    if(_periodic == Common::PERIODIC_XY){
        //NumMatrix unshifted_polynomial_coefficents = _polynomial_coefficents;
        //this->set_periodicity(_periodic, 1.);
        NumMatrixMult(_polynomial_coefficents, values_of_basis_functions, function_value);
        //this->set_periodicity(_periodic, 1., true);
        //_polynomial_coefficents = unshifted_polynomial_coefficents;
        
    } else {
        NumMatrixMult(_polynomial_coefficents, values_of_basis_functions, function_value);
    }
    //NumMatrixMult(_polynomial_coefficents, values_of_basis_functions, function_value);
    NumVector function_value_vector(function_value.m());

    for(int i = 0; i < function_value.m(); i++){
        function_value_vector(i) = function_value(i,0);
    }

    return function_value_vector;

}


// implement the Bezier forward differencing operator given on page 49 of CAGD
// by Gerald Farin.
Function Function::d_dx(size_t i, Function f){
    //assert(f._basis_type == Common::BEZIER);

    // copy the contents of f into the derivative
    Function df_dxi = f;
    df_dxi._pseudo_inverse = f._pseudo_inverse;
    //df_dxi._basis_function_values= f._basis_function_values;
    if(f._basis_type == Common::BEZIER){
        Bezier::derivative_of_tensor_product_bezier(
                f._polynomial_coefficents,
                i == 0 ? Common::U  : Common::V,
                f._range_dim,
                f._basis_bidegree,
                df_dxi._polynomial_coefficents
                );
    } else if(f._basis_type == Common::BSPLINE){
        // compute new coefficients of degree (n-1) x (n-1) b-spline basis
        /*
           BSpline::derivative_of_tensor_product_bspline(
           f._polynomial_coefficents,
           i == 0 ? Common::U  : Common::V,
           f._range_dim,
           f._basis_bidegree,
           df_dxi._polynomial_coefficents
           );*/

        // create smaller, lower order basis set 
        //df_dxi._basis_functions = df_dxi.generate_basis_set(f._basis_bidegree-1, f._num_basis_functions-1, Common::BSPLINE);
        /*
           df_dxi.generate_basis_set(f._basis_bidegree-1, f._num_basis_functions-1, Common::BSPLINE);
           int remaining_basis_functions = f._num_basis_functions - df_dxi._basis_functions.size();
           for(int i = 0; i < remaining_basis_functions; i++){
           df_dxi._basis_functions.push_back(new Function::ZeroBasisFunction());
           }
           df_dxi._basis_bidegree = f._basis_bidegree-1;
           */
        //df_dxi.generate_basis_set(f._basis_bidegree, f._num_basis_functions, 
        df_dxi._basis_functions_x = df_dxi.generate_basis_set(f._basis_bidegree, f._num_basis_functions, f._basis_type);
        df_dxi._basis_functions_y = df_dxi.generate_basis_set(f._basis_bidegree, f._num_basis_functions, f._basis_type);
        int num_basis_functions = df_dxi._basis_functions_x.size();
        df_dxi._polynomial_coefficents = f._polynomial_coefficents;
        multvalue( df_dxi._polynomial_coefficents, double(df_dxi._num_basis_functions - df_dxi._basis_bidegree));
        //multvalue( df_dxi._polynomial_coefficents, .5);

        for(int j =0; j < num_basis_functions; j++){
            df_dxi._basis_functions_x[j]->_eval_type = f._basis_functions_x[j]->_eval_type;
            df_dxi._basis_functions_y[j]->_eval_type = f._basis_functions_y[j]->_eval_type;
        }
        for(int j =0; j < num_basis_functions; j++){
            //df_dxi._basis_functions_x[j]->_eval_type = f._basis_functions_x[j]->_eval_type;
            //df_dxi._basis_functions_y[j]->_eval_type = f._basis_functions_y[j]->_eval_type;

            BasisFunction* dfdxi_basis_function;
            BasisFunction* f_basis_function;
            if(i == 0){
                dfdxi_basis_function= df_dxi._basis_functions_x[j];
                f_basis_function= f._basis_functions_x[j];
            } else if (i == 1){
                dfdxi_basis_function = df_dxi._basis_functions_y[j];
                f_basis_function = f._basis_functions_y[j];

            }

            //cout << "before: " << i <<  endl;
            //cout << dfdxi_basis_function->_eval_type << endl;
            switch(f_basis_function->_eval_type){
                case VALUE:
                    dfdxi_basis_function->_eval_type = DERIV_1ST;
                    break;
                case DERIV_1ST:
                    //basis_function->_eval_type = DERIV_1ST;
                    dfdxi_basis_function->_eval_type = DERIV_2ND;
                    break;
                case DERIV_2ND:
                    dfdxi_basis_function->_eval_type = DERIV_3RD;
                    break;
                case DERIV_3RD:
                    dfdxi_basis_function->_eval_type = DERIV_4TH;
                    break;
                case DERIV_4TH:
                    assert(0);
                    break;
                default:
                    assert(0);
            }
            //cout << "after" << endl;
            //cout << dfdxi_basis_function->_eval_type << endl;

        }
    }
    return df_dxi;
}


NumMatrix Function::compute_basis_function_values(
        size_t domain_dim,
        vector<BasisFunction*>& basis_functions_x,
        vector<BasisFunction*>& basis_functions_y,
        NumVector x){

    assert(domain_dim == 2);
    assert(x.length() == domain_dim);
    assert(basis_functions_y.size() == basis_functions_x.size());
    int num_basis_functions = basis_functions_x.size();
    NumMatrix basis_function_values_1d(domain_dim, num_basis_functions);
    setvalue(basis_function_values_1d, 0.);
    for(int j = 0; j < num_basis_functions; j++){
        basis_function_values_1d(0, j) = basis_functions_x[j]->evaluate(x(0));
        basis_function_values_1d(1, j) = basis_functions_y[j]->evaluate(x(1));
    }
    return basis_function_values_1d;
}


NumMatrix Function::compute_basis_function_values(
        size_t domain_dim,
        vector<BasisFunction*> basis_functions,
        NumVector x){

    assert(x.length() == domain_dim);
    int num_basis_functions = basis_functions.size();
    NumMatrix basis_function_values_1d(domain_dim, num_basis_functions);
    setvalue(basis_function_values_1d, 0.);

    for(int j = 0; j < num_basis_functions; j++){
        BasisFunction* basis_function = basis_functions[j];
        for(int d = 0; d < domain_dim; d++){
            basis_function_values_1d(d, j) = basis_function->evaluate(x(d));
        }
    }
    return basis_function_values_1d;
}

NumMatrix Function::compute_basis_function_values(
        size_t domain_dim,
        size_t basis_bidegree,
        size_t num_basis_functions,
        IntMatrix binomial_coefficients,
        Common::BasisType basis_type,
        NumVector x){

    assert(x.length() == domain_dim);
    NumMatrix basis_function_values_1d(domain_dim, num_basis_functions);
    //setvalue(basis_function_values_1d, 0.);

    for(int j = 0; j < num_basis_functions; j++){
        for(int d = 0; d < domain_dim; d++){
            if(basis_type == BEZIER){
                int binomial_coefficient = binomial_coefficients(basis_bidegree, j);
                basis_function_values_1d(d, j) = 
                    FaceMapSurf::B(j, basis_bidegree, x(d), binomial_coefficient);

            } else if(basis_type == BSPLINE){
                basis_function_values_1d(d, j) = BSpline::evaluate_bspline(j, basis_bidegree, num_basis_functions, x(d));
            }
        }
    }
    return basis_function_values_1d;
}

NumMatrix Function::compute_basis_function_values(NumVector x){
    assert(x.length() == _domain_dim);

    //if(_cached_matrices && _values_at_basis_functions != NULL){
    if(false){ // TODO cache for evaluation...?
        return *_values_at_basis_functions;
    } else {
        if(_domain_dim == 2 && _basis_type == Common::BSPLINE){
            return compute_basis_function_values(
                    _domain_dim,
                    _basis_functions_x, 
                    _basis_functions_y, 
                    x);
        } else{
            return compute_basis_function_values(
                    _domain_dim,
                    _basis_functions, 
                    x);
        }
    }
}

NumMatrix Function::compute_pseudo_inverse(
        vector<BasisFunction*> basis_functions_x,
        vector<BasisFunction*> basis_functions_y,
        NumMatrix coordinates){
    const int N = coordinates.m();
    int num_basis_functions_x = basis_functions_x.size();
    int num_basis_functions_y = basis_functions_y.size();
    //const int M = values.m();

    assert(N == 2);
    assert(dynamic_cast<BSplineBasisFunction*>(basis_functions_x[0])
        && dynamic_cast<BSplineBasisFunction*>(basis_functions_y[0])
          );
    int num_samples = coordinates.n();
    // possible b spline bug in pseudo inverse size
    NumMatrix A(num_samples, num_basis_functions_x*num_basis_functions_y);
    setvalue(A, 0.);


    for(int i = 0; i < num_samples; i++){
        NumVector x(N, false, coordinates.clmdata(i));
        // TODO REIMPLEMENT FOR ARBITRARY DOMAINDIMENSION
        NumMatrix basis_function_values_1d = 
            compute_basis_function_values(N, basis_functions_x, basis_functions_y, x);
        for(int j = 0; j < num_basis_functions_x; j++){
            for(int k = 0; k < num_basis_functions_y; k++){

                A(i,j*(num_basis_functions_x) +k) = basis_function_values_1d(0,j)*basis_function_values_1d(1,k);
            }
        }
    }
    
    // maps function values to polynomial coefficients
    NumMatrix pseudo_inverse(num_basis_functions_x*num_basis_functions_y, num_samples);
    setvalue(pseudo_inverse, 0.);
    pinv(A, 1e-10, pseudo_inverse);
    return pseudo_inverse;

}



NumMatrix Function::compute_pseudo_inverse(
        vector<BasisFunction*> basis_functions,
        NumMatrix coordinates){
    const int N = coordinates.m();
    int num_basis_functions = basis_functions.size();
    //const int M = values.m();

    assert(N == 2 || N == 1);
    int num_samples = coordinates.n();
    // possible b spline bug in pseudo inverse size
    NumMatrix A(num_samples, pow(num_basis_functions, N));
    setvalue(A, 0.);


    for(int i = 0; i < num_samples; i++){
        NumVector x(N, false, coordinates.clmdata(i));
        // TODO REIMPLEMENT FOR ARBITRARY DOMAINDIMENSION
        //NumMatrix basis_function_values_1d = 
            //compute_basis_function_values(N, basis_bidegree, num_basis_functions, binomial_coefficients, basis_type, x);
        NumMatrix basis_function_values_1d = 
            compute_basis_function_values(N, basis_functions, x);
            //compute_basis_function_values(N, basis_bidegree, num_basis_functions, binomial_coefficients, basis_type, x);
        //Function::compute_basis_function_values_cache(N, basis_bidegree, binomial_coefficients, x);
        if(N == 1){
            for(int j = 0; j < num_basis_functions; j++){
                A(i,j) = basis_function_values_1d(0,j);
            }
        } else if (N == 2){
            for(int j = 0; j < num_basis_functions; j++){
                for(int k = 0; k < num_basis_functions; k++){

                    A(i,j*(num_basis_functions) +k) = basis_function_values_1d(0,j)*basis_function_values_1d(1,k);
                }
            }
        } else { 
            assert(0); 
        }
    }
    
    // maps function values to polynomial coefficients
    NumMatrix pseudo_inverse(pow(num_basis_functions, N), num_samples);
    setvalue(pseudo_inverse, 0.);
    pinv(A, 1e-10, pseudo_inverse);
    return pseudo_inverse;

}



NumMatrix Function::compute_pseudo_inverse(
        int basis_bidegree,
        int num_basis_functions,
        IntMatrix binomial_coefficients,
        BasisType basis_type,
        NumMatrix coordinates
        ){
        //NumMatrix values){
    const int N = coordinates.m();
    //const int M = values.m();

    assert(N == 2 || N == 1);
    int num_samples = coordinates.n();
    // possible b spline bug in pseudo inverse size
    NumMatrix A(num_samples, pow(num_basis_functions, N));
    setvalue(A, 0.);


    for(int i = 0; i < num_samples; i++){
        NumVector x(N, false, coordinates.clmdata(i));
        // TODO REIMPLEMENT FOR ARBITRARY DOMAINDIMENSION
        //NumMatrix basis_function_values_1d = 
            //compute_basis_function_values(N, basis_bidegree, num_basis_functions, binomial_coefficients, basis_type, x);
        NumMatrix basis_function_values_1d = 
            compute_basis_function_values(N, basis_bidegree, num_basis_functions, binomial_coefficients, basis_type, x);
        //Function::compute_basis_function_values_cache(N, basis_bidegree, binomial_coefficients, x);
        if(N == 1){
            for(int j = 0; j < num_basis_functions; j++){
                A(i,j) = basis_function_values_1d(0,j);
            }
        } else if (N == 2){
            for(int j = 0; j < num_basis_functions; j++){
                for(int k = 0; k < num_basis_functions; k++){

                    A(i,j*(num_basis_functions) +k) = basis_function_values_1d(0,j)*basis_function_values_1d(1,k);
                }
            }
        } else { 
            assert(0); 
        }
    }
    
    // maps function values to polynomial coefficients
    NumMatrix pseudo_inverse(pow(num_basis_functions, N), num_samples);
    setvalue(pseudo_inverse, 0.);
    pinv(A, 1e-10, pseudo_inverse);
    return pseudo_inverse;

}


NumMatrix Function::compute_pseudo_inverse(
        NumMatrix coordinates){
    if(_cached_matrices && _pseudo_inverse != NULL){
        //cout << "cached pinv" << endl;
        return *_pseudo_inverse;
    } else {
        if(_domain_dim == 2 && _basis_type == Common::BSPLINE){
            return compute_pseudo_inverse(
                    _basis_functions_x, 
                    _basis_functions_y, 
                    coordinates);
        } else {
            return compute_pseudo_inverse(
                    _basis_functions,
                    coordinates);
        }
    }
}

NumMatrix Function::compute_polynomial_approximation(NumMatrix coordinates, 
        NumMatrix values){
    assert(coordinates.n() == values.n()); // we need num coordindates == num values
    //assert(this->_basis_type == Common::BEZIER);
    const int N = coordinates.m();
    const int M = values.m();
    assert(N == 2 || N == 1); // 
    
    
    NumMatrix pseudo_inverse = compute_pseudo_inverse(coordinates);
    NumMatrix values_t(values.n(), values.m());
    NumMatrixTranspose(values, values_t);
    int num_basis_functions;
    //if(_basis_type == Common::BSPLINE && N == 2)
        //num_basis_functions = _basis_functions_x.size()*_basis_functions_y.size();
    //else 
        num_basis_functions = pow(_num_basis_functions, N);
    //NumMatrix basis_coefficients_t(pow(_num_basis_functions, N), M);
    NumMatrix basis_coefficients_t(num_basis_functions, M);
    setvalue(basis_coefficients_t, 0.);
    NumMatrixMult(pseudo_inverse, values_t, basis_coefficients_t);
    
    NumMatrix basis_coefficients(basis_coefficients_t.n(), basis_coefficients_t.m());
    NumMatrixTranspose(basis_coefficients_t, basis_coefficients);
    
    _polynomial_coefficents = basis_coefficients;
    return basis_coefficients;
    
}

Function Function::construct_function(NumMatrix coordinates, NumMatrix values, 
        Common::BasisType basis_type, size_t basis_order){
    // check domain/range sizes
    assert(coordinates.n() == values.n()); // we need num coordindates == num values

    int N = coordinates.m();
    int M = values.m();
    Function func(N, M, basis_order, basis_type);
    //func.set_polynomial_basis(basis_type, basis_order);
    
    // least square fit for basis function coefficients
    //NumMatrix polynomial_coeffs = 
    func.compute_polynomial_approximation(coordinates, values);
    // store coefficients in func

    return func;
}

Function Function::construct_spline_function(NumMatrix coordinates, NumMatrix values, 
        size_t num_control_points, size_t basis_order){
    // check domain/range sizes
    assert(coordinates.n() == values.n()); // we need num coordindates == num values

    int N = coordinates.m();
    int M = values.m();
    Function func(N, M, basis_order, num_control_points, Common::BSPLINE);
    //func.set_polynomial_basis(basis_type, basis_order);
    
    // least square fit for basis function coefficients
    //NumMatrix polynomial_coeffs = 
    func.compute_polynomial_approximation(coordinates, values);
    // store coefficients in func

    return func;
}

void reassign(int i, int j, int i_periodic_index, int j_periodic_index, 
        int _num_basis_functions, Point3 coeff_shift, NumMatrix _polynomial_coefficents, NumMatrix&periodic_coefficients){
    /*
    cout << "num_control_points: " << _num_basis_functions << endl;
        cout << "i: " << i << endl;
        cout << "j: " << j << endl;
        cout << "i_periodic_index: " << i_periodic_index << endl;
        cout << "j_periodic_index: " << j_periodic_index << endl;*/
    assert(j_periodic_index >=0 && j_periodic_index < _num_basis_functions);

    assert(i_periodic_index >=0 && i_periodic_index < _num_basis_functions);
    int  index = i*_num_basis_functions + j;
    int periodic_index = i_periodic_index*_num_basis_functions+j_periodic_index;
    for(int d = 0; d < _polynomial_coefficents.m(); d++){
        periodic_coefficients(d,periodic_index) = _polynomial_coefficents(d, index) + coeff_shift(d);
    }
}

void Function::set_periodicity(Common::Periodicity periodic_type, double period_width, bool unshift){
    assert(_basis_type == Common::BSPLINE);
    if(!unshift){
        int n = _basis_bidegree;
        assert(n % 2 == 1); // need odd degree b-spline
        int padding_size = (n-1)/2;
        //int padding_size = 0;
        int interval_size = _num_basis_functions - 2*padding_size-1;
        assert(periodic_type == Common::PERIODIC_XY);

        _periodic = periodic_type;
        _period_width = period_width;
        NumMatrix periodic_coefficients(_polynomial_coefficents.m(), _polynomial_coefficents.n());

        setvalue(periodic_coefficients, 0.);
        //NumMatrix periodic_coefficients = _polynomial_coefficents;
        for(int i = 0; i < _polynomial_coefficents.m(); i++){
            for(int j = 0; j < _polynomial_coefficents.n(); j++){
                periodic_coefficients(i,j) = _polynomial_coefficents(i,j);
            }
        }
        for(int i =0; i <_num_basis_functions; i++){
            for(int j =0; j <_num_basis_functions; j++){
                Point3 coeff_shift(0.);

                int i_periodic_index;
                int j_periodic_index;
                if(i < padding_size){
                    i_periodic_index = i + interval_size;
                    coeff_shift(0) = -_period_width;
                } else if ( i >= _num_basis_functions - padding_size){
                    i_periodic_index = i - interval_size;
                    coeff_shift(0) = _period_width;
                } else {
                    i_periodic_index = i;
                    coeff_shift(0) = 0.;
                }
                if(j < padding_size){
                    j_periodic_index = j + interval_size;
                    coeff_shift(1) = -_period_width;
                } else if ( j >= _num_basis_functions - padding_size){
                    j_periodic_index = j - interval_size;
                    coeff_shift(1) = _period_width;
                } else {
                    j_periodic_index = j;
                    coeff_shift(1) = 0.;
                }
                //cout << "padding_size: " << padding_size << endl;
                reassign( i_periodic_index, j_periodic_index, i,j,
                        _num_basis_functions, coeff_shift, _polynomial_coefficents, periodic_coefficients);

            }
        }
    
    for(int i = 0; i < _polynomial_coefficents.m(); i++){
        for(int j = 0; j < _polynomial_coefficents.n(); j++){
            _polynomial_coefficents(i,j) = periodic_coefficients(i,j);
        }
    }
    } else{
        int n = _basis_bidegree;
        assert(n % 2 == 1); // need odd degree b-spline
        int padding_size = (n-1)/2;
        //int padding_size = 0;
        int interval_size = _num_basis_functions - 2*padding_size-1;
        assert(periodic_type == Common::PERIODIC_XY);

        _periodic = periodic_type;
        _period_width = period_width;
        NumMatrix periodic_coefficients(_polynomial_coefficents.m(), _polynomial_coefficents.n());

        setvalue(periodic_coefficients, 0.);
        //NumMatrix periodic_coefficients = _polynomial_coefficents;
        for(int i = 0; i < _polynomial_coefficents.m(); i++){
            for(int j = 0; j < _polynomial_coefficents.n(); j++){
                periodic_coefficients(i,j) = _polynomial_coefficents(i,j);
            }
        }
        for(int i =0; i <_num_basis_functions; i++){
            for(int j =0; j <_num_basis_functions; j++){
                Point3 coeff_shift(0.);

                int i_periodic_index;
                int j_periodic_index;
                if(i < padding_size){
                    i_periodic_index = i + interval_size;
                    coeff_shift(0) = -_period_width;
                } else if ( i >= _num_basis_functions - padding_size){
                    i_periodic_index = i - interval_size;
                    coeff_shift(0) = _period_width;
                } else {
                    i_periodic_index = i;
                    coeff_shift(0) = 0.;
                }
                if(j < padding_size){
                    j_periodic_index = j + interval_size;
                    coeff_shift(1) = -_period_width;
                } else if ( j >= _num_basis_functions - padding_size){
                    j_periodic_index = j - interval_size;
                    coeff_shift(1) = _period_width;
                } else {
                    j_periodic_index = j;
                    coeff_shift(1) = 0.;
                }
                //cout << "padding_size: " << padding_size << endl;
                reassign( i_periodic_index, j_periodic_index, i,j,
                        _num_basis_functions, coeff_shift, _polynomial_coefficents, periodic_coefficients);

            }
        }
    
    for(int i = 0; i < _polynomial_coefficents.m(); i++){
        for(int j = 0; j < _polynomial_coefficents.n(); j++){
            _polynomial_coefficents(i,j) = periodic_coefficients(i,j);
        }
    }
    }
    
}

