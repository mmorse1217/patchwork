#ifndef __BEZIER_HPP__
#define __BEZIER_HPP__
#include "nummatrix.hpp"
#include "common.hpp"
#include "vec3t.hpp"
#include "vec2t.hpp"
namespace Bezier {
    
    /**
     * Derived from Gerald Farin's Curves and Surface for CAGD (4th ed. Chapter
     * 15.5 and 15.6)
     * Suppose a function 
     
     * f(u,v) = \sum_i=0^n\sum_j=0^n a_ij*B_i^n(u)*B_j^n(v). 
     *
     * then, by the first eq of 15.6, we have
     *
     * df/du = \sum_i=0^{n-1}\sum_j=0^n (a_{i+1,j} - a_ij)*B_i^{n-1}(u)*B_j^n(v).
     *
     * This is a bummer, because the basis elements changed from B_i^n to
     * B_i^{n-1}, and so the number of basis elements changed. To fix this, we
     * elevate the degree of the bezier approximation of f w.r.t u, which is
     * determined by a linear interpolation of the n-1 terms in the parameter
     * i/n (Chapter 5.1/eq 15.8). 
     *
     * The ith coefficient of the elevated representation is given by:
     * 
     * b_{i,j}^{(1,0)} = i/n*b_{i-1,j} + (1-i/n)*b_{i,j}, i,j=0,...,n
     * where b_ij are control points of the original un-elevated function 
     *
     * This gives us:
     *
     * df/du = \sum_i=0^{n-1}\sum_j=0^n C^n_{ij}*B_i^{n-1}(u)*B_j^n(v).
     * where 
     * C^n_{ij} = n*(i/n(a_{i,j} - a_{i-1,j}) + (1-i/n)*(a_{i+1,j} - a_{ij}))
     *          = i(a_{i,j} - a_{i-1,j}) + (n-i)*(a_{i+1,j} - a_{ij}))
     *          = (n-i)*(a_{i+1,j} + (2i-n)*a_{ij}) - i*a_{i-1,j}
     * 
     * This gives us a tri-diagonal matrix that maps n^2 coefficients of a
     * tensor product bezier patch to n^2 coefficients of it's derivative.
     *
     * @param   NumMatrix tensor_product_coeffs   list of 
     *                                            basis_bidegree^2 x coeff_dim
     *                                            control points for a bezier
     *                                            patch. (i*basis_bidegree +j)th 
     *                                            row is the i,jth coefficient
     *                                            and the kth column is the kth
     *                                            component of the coefficient, 
     *                                            k=1,... coeff_dim
     * @param Direction    dx_i                   which variable to
     *                                            differentiate w.r.t 
     * @param size_t       coeff_dim              dimension of the space the
     *                                            polynomial maps to
     * @param size_t       basis_bidegree         number of basis elements per
     *                                            dimension
     *
     *
     */
   void derivative_of_tensor_product_bezier(
           NumMatrix tensor_product_coeffs,
           Common::Direction dx_i,
           size_t coeff_dim,
           size_t basis_bidegree,
            NumMatrix& derivative_coeffs
           );
    /**
     * Helper function to initialize the derivative coefficient special cases of
     * i = 0, j = 0, i = n, j = n. For example for d/dv, this sets either
     *  C_i,0 = n*(a_i,1 - a_i,0) or,
     *  C_i,n = n*(a_i,n - a_i,n-1); similarly for i
     *
     * @param   size_t      n                   polynomial bidegree
     * @param   size_t      i                   polynomial index w.r.t u
     * @param   size_t      j                   polynomial index w.r.t v
     * @param   size_t      coeff_dim           dimension of range poly is
     *                                          mapping to
     * @param   Direction   dx_i                which variable we're
     *                                          differentiating w.r.t
     * @param   NumMatrix   original_coeffs     original function coefficients
     *                                          in bernstein basis
     * @param   NumMatrix&  derivative_coeffs   where to store the derivative
     *                                          coefficients
     */
    void set_value_edge_case(size_t n, // polynomial bidegree
        size_t i,  
        size_t j,  // polynomial index w.r.t v
        size_t coeff_dim,   //degree of space poly's are mapping to (usually R^3)
        Common::Direction dx_i,
        NumMatrix original_coeffs, NumMatrix& derivative);

    /**
     * Helper function to initialize the derivative coefficient for the general
     * i,j case, i,j \not = 0, n
     * For d/du, this sets 
     *  C_i,j = (n-i)*(a_{i+1,j} + (2i-n)*a_{ij}) - i*a_{i-1,j}
     *
     * @param   size_t      n                   polynomial bidegree
     * @param   size_t      i                   polynomial index w.r.t u
     * @param   size_t      j                   polynomial index w.r.t v
     * @param   size_t      coeff_dim           dimension of range poly is
     *                                          mapping to
     * @param   Direction   dx_i                which variable we're
     *                                          differentiating w.r.t
     * @param   NumMatrix   original_coeffs     original function coefficients
     *                                          in bernstein basis
     * @param   NumMatrix&  derivative_coeffs   where to store the derivative
     *                                          coefficients
     */
    void set_value(size_t n, 
        size_t i,  
        size_t j,  
        size_t coeff_dim,   
        Common::Direction dx_i,
        NumMatrix original_coeffs, NumMatrix& derivative);

    // TODO 
    // also change interface to use Function() and accept matrix of
    // coefficients instead of vector of Point3's :(
    Point3 blossom(vector<Point3> control_points, 
            vector<double> parameter_values, 
            int n, 
            int i);
    // blossom with parameter_values = [t, t, t, ..., t]
    Point3 de_casteljau(vector<Point3> control_points,
            double t);
    Point3 blossom_recursive(vector<Point3> control_points,
            vector<double> parameter_values,
            int n,
            int i);
    // blossom with parameter_values = [t, t, t, ..., t]
    Point3 de_casteljau_recursive(vector<Point3> control_points,
            double t);
    Point3 tensor_product_blossom(vector<Point3> control_points,
            vector<Point2> parameter_values);
    Point3 tensor_product_blossom_recursive(vector<Point3> control_points,
            vector<Point2> parameter_values);

    vector<Point3> compute_control_points_on_subdomain(
            vector<Point3> control_points,
            int patch_order,
            pair<double, double> x_interval,
            pair<double, double> y_interval);

};
#endif
