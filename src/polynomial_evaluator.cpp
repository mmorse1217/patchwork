#include "polynomial_evaluator.hpp"
#include <string>

PolynomialEvaluator::PolynomialEvaluator(string polynomial_filename, 
        string mesh_filename){
    ifstream polynomials(polynomial_filename.c_str());
    ifstream mesh(mesh_filename.c_str());

    _mesh =  new GpMesh(); 
    _mesh->setup(mesh, 0, Point3(1.));
    
    polynomials >> _num_patches;
    polynomials >> _poly_degree;
    cout << "concrete poly degree: " << _poly_degree << endl;

    int poly_order = (_poly_degree+1)*(_poly_degree+1);
    
    _polynomial_patches = new vector<Vector<Point3> >(
                        _num_patches, 
                        Vector<Point3>(poly_order));

    for(int pi = 0; pi < _num_patches; pi++){
        Vector<Point3> polynomial_coeffs(poly_order);

        for(int ci = 0; ci < poly_order; ci++){
            Point3 coeff;
            polynomials >> coeff;

            polynomial_coeffs(ci) = coeff;
        }
        _polynomial_patches->at(pi) = polynomial_coeffs;
    /*
    double dx = 0.673469;
    double dy = -0.632653;
    Point3 deltax(0.,0.,dx);
    Point3 deltay(0.,0.,dy);
    for (int i = 0; i < 4; i+=3) {
        for (int j = 1; j < 3; j++) {
            _polynomial_patches->at(pi)(i*4 + j) += deltax/2.; 
        }
    }
    for (int i = 1; i < 3; i++) {
        for (int j = 0; j < 4; j+=3) {
            _polynomial_patches->at(pi)(i*4 + j) += deltay/2.; 
        }
    }
    for (int i = 0; i < 4; i+=3) {
        for (int j = 0; j < 4; j+=3) {
            _polynomial_patches->at(pi)(i*4 + j) += (deltay+deltax)*.5; 
            
        }
    }
    */
    
    //for (int i = 0; i < 4; i=i+3) {
    //    for (int j = 0; j < 4; j=j+3) {
    //        _polynomial_patches->at(pi)(i*4 + j) += delta/8.; 
    //    }
    //}
    // move center coefficients by t/2
    /*for (int i = 1; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            _polynomial_patches->at(pi)(i*4+j) -= delta/3.; 
        }
    }*/
    /*for (int i = 0; i < 4; i=i+3) {
        for (int j = 1; j < 3; j++) {
            _polynomial_patches->at(pi)(i*4+j) -= delta/2.; 
            _polynomial_patches->at(pi)(j*4+i) -= delta/2.; 
        }
    }*/
    }
    polynomials.close();
    _binomial_coefs = FaceMapSurf::compute_binomial_coefs_array(_poly_degree);
}

PolynomialEvaluator::~PolynomialEvaluator(){
    delete _polynomial_patches;
    delete _binomial_coefs;
}

vector<Point3> PolynomialEvaluator::evaluate(int patch_id, int evaluation_flags, Point2 uv_coords){
    
    int return_value_size;
    if(evaluation_flags & EVAL_1ST_DERIV){
        return_value_size = 6;
    } else if (evaluation_flags & EVAL_VALUE){
        return_value_size = 2;
    } else {
        assert(0); // second derivatives not yet implementated
    }
    
    vector<Point3> positions_and_derivs(return_value_size/2, Point3(0.));
    Vector<Point3>* polynomial_coefs = &((*_polynomial_patches)[patch_id]);

    FaceMapSurf::evaluate_bernstein_polynomial_patch(
            polynomial_coefs,
            _binomial_coefs,
            _poly_degree,
            evaluation_flags,
            uv_coords.array(),
            positions_and_derivs.data());
    return positions_and_derivs;
}
void PolynomialEvaluator::setup(){
    assert(0);
}

void PolynomialEvaluator::bounding_box(Point3& min, Point3& max){
    _mesh->bbox(min, max);
}
