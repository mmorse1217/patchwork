#ifndef __DIFF_GEOM_HPP__
#define __DIFF_GEOM_HPP__
#include "face_map.hpp"
#include "vec3t.hpp"
#include "nummatrix.hpp"
#include "function.hpp"

namespace Differential{
    // Note: for a function f: R^n  \to R^m, with x = (x_1, ..., x_2) \in R^n, 
    // we write df/dx_i = f_{x_i}. if f(x,y,z) = y \in R^m, we'll just write
    // f_x, f_y, f_z. 
    /**
     * For a point on a parametrized surface X(u,v), compute the elements of the
     * matrix
     * [E  F]
     * [F  G],
     * where E = X_u \dot X_u, F = X_u \dot X_v, G = X_v \dot X_v
     * @param Point3    X_u         derivative at X(u,v) w.r.t. u
     * @param Point3    X_v         derivative at X(u,v) w.r.t. v
     * @ret   vector<double>        list of three doubles corresponding to E,F,G
     *                              in that order
     */
    vector<double> first_fundamental_form(Point3 X_u, Point3 X_v); 
    vector<double> first_fundamental_form(Function X_u, Function X_v, Point2 uv_coords); 

    /**
     * Let n = X_u \cross X_v / |X_u \cross X_v|. For a point on a parametrized 
     * surface X(u,v), compute the elements of the matrix
     * [L  M]
     * [M  N],
     * where L = X_uu \dot n, F = X_uv \dot n, G = X_vv \dot n 
     * @param Point3    X_u         first dersvative at X(u,v) w.r.t. u
     * @param Point3    X_v         first derivative at X(u,v) w.r.t. v
     * @param Point3    X_uu        second derivative at X(u,v) w.r.t. u twice
     * @param Point3    X_uv        second derivative at X(u,v) w.r.t. u, then v
     *                              (equal to deriv w.r.t v then u)
     * @param Point3    X_vv        derivative at X(u,v) w.r.t. v twice 
     * @ret   vector<double>        list of three doubles corresponding to L,M,N
     *                              in that order
     */
    vector<double> second_fundamental_form(Point3 X_u, Point3 X_v, 
            Point3 X_uu, Point3 X_uv, Point3 X_vv);
    vector<double> second_fundamental_form(Function X_u, Function X_v, 
            Function X_uu, Function X_uv, Function X_vv, Point2 uv_coords);

    vector<double> first_fundamental_form_partial_derivatives(
            Function X_u, Function X_v, 
            Function X_uu, Function X_uv, Function X_vv,
            Point2 uv_coords);
    vector<double> first_fundamental_form_partial_derivatives(Point3 X_u, Point3 X_v, 
            Point3 X_uu, Point3 X_uv, Point3 X_vv);


    /**
     * Return sqrt(EG*F^2)
     */

    double determinant_of_first_fundamental_form(vector<double> first_fundamental_form_coefs);
    
    vector<double> derivatives_of_det_of_first_form(
        vector<double> first_form_coefs,
        vector<double> first_form_derivatives);

    double gaussian_curvature(Point3 X_u, Point3 X_v, 
            Point3 X_uu, Point3 X_uv, Point3 X_vv);

double mean_curvature(Function X_u, Function X_v, 
        Function X_uu, Function X_uv, Function X_vv, Point2 uv);


    double mean_curvature(Point3 X_u, Point3 X_v, 
            Point3 X_uu, Point3 X_uv, Point3 X_vv);

    Point3 surface_gradient(Function X_u, Function X_v, 
        Function f_u, Function f_v, Point2 uv_coords);
    
    Point3 surface_gradient(Point3 X_u, Point3 X_v, 
            double f_u, double f_v);

    NumMatrix surface_gradient(Point3 X_u, Point3 X_v, 
            Point3 V_u, Point3 V_v);

    double surface_divergence(Function X_u, Function X_v, 
            Function V_u, Function V_v, Point2 uv_coords);

    double surface_divergence(Point3 X_u, Point3 X_v, 
            Point3 V_u, Point3 V_v);


    double laplace_beltrami(
            vector<Function> X_derivatives,
            vector<Function> f_derivatives,
            Point2 uv_coords);

    double laplace_beltrami(Point3 X_u, Point3 X_v, 
        Point3 X_uu, Point3 X_uv, Point3 X_vv,
        double f_u, double f_v,
        double f_uu, double f_uv, double f_vv);

    Point3 laplace_beltrami(
        Point3 X_u, Point3 X_v, Point3 X_uu, Point3 X_uv, Point3 X_vv,
        Point3 V_u, Point3 V_v, Point3 V_uu, Point3 V_uv, Point3 V_vv);


    // above functions re-defined per patch as a functino of uv-coordindates
    double gaussian_curvature(Point2 uv, FaceMapSurf* surface, int patch_id);
    double mean_curvature(Point2 uv, FaceMapSurf* surface, int patch_id);
    
    Point3 surface_gradient(Point2 uv, double f_u, double f_v,
            FaceMapSurf* surface, int patch_id);

    double surface_divergence(Point2 uv, Point3 V_u, Point3 V_v,
            FaceMapSurf* surface, int patch_id);
    
    double laplace_beltrami(Point2 uv, double f_u, double f_v,
        double f_uu, double f_uv, double f_vv,
            FaceMapSurf* surface, int patch_id);
    
    Point3 laplace_beltrami(Point2 uv, Point3 V_u, Point3 V_v,
        Point3 V_uu, Point3 V_uv, Point3 V_vv,
            FaceMapSurf* surface, int patch_id);


};
#endif
