#include "utils.hpp"

double sixth_degree_poly(double x, double y){
    return pow(x,6)*pow(y,6) + pow(x,5)*pow(y,5);
}

double sixth_degree_poly_dx(double x, double y){
    return 6*pow(x,5)*pow(y,6)+ 5*pow(x,4)*pow(y,5);
}
double sixth_degree_poly_dy(double x, double y){
    return 6*pow(x,6)*pow(y,5)+ 5.*pow(x,5)*pow(y,4);
}
double constant(double, double){
    return 1.;
}

double sin7(double x){
    return pow(sin(x), 7);
}
double sixth_degree_poly_1d(double x){
    return pow(x, 6);
}

double sin2cos2(double x, double y){
    return sin(x)*sin(x)*cos(y)*cos(y) +1.;
}
double trig_dx(double x, double y){
    return sin(2*x)*cos(y)*cos(y);
}
double trig_dy(double x, double y){
    //return -2*sin(x)*sin(x)*cos(y)*sin(y);
    return -sin(x)*sin(x)*sin(2*y);
}

Point3 paraboloid(double x, double y){
    return Point3(x, y, x*x - y*y);
}

Point3 paraboloid_dx(double x, double ){
    return Point3(1., 0., 2.*x);
}

Point3 paraboloid_dy(double , double y){
    return Point3(0., 1., -2.*y);
}
Point3 zero(double, double){
    return Point3(0.);
}
Point3 paraboloid_dxdx(double, double){
    return Point3(0., 0., 2.);
}

Point3 constant_plane(double x, double y){
    return Point3(x, y, 1.);
}
Point3 constant_plane_dx(double, double ){
    return Point3(1, 0, 0.);
}
Point3 constant_plane_dy(double, double){
    return Point3(0., 1, 0.);
}

Point3 sinsin(double x, double y){
    //return Point3(x, y, sin(2*M_PI*y));
    return Point3(x, y, .25*sin(2*M_PI*x)*sin(2*M_PI*y));
}
Point3 sinsin_dx(double x, double y){
    //return Point3(1, 0, 2*M_PI*cos(2*M_PI*x));
    return Point3(1, 0, .25*2*M_PI*cos(2*M_PI*x)*sin(2*M_PI*y));
}
Point3 sinsin_dy(double x, double y){
    //return Point3(1, 0, 0.);
    return Point3(0, 1, .25*2*M_PI*sin(2*M_PI*x)*cos(2*M_PI*y));
}
