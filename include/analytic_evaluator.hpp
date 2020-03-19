#ifndef __ANALYTIC_EVALUATOR_HPP__
#define __ANALYTIC_EVALUATOR_HPP__

#include "surface_evaluator.hpp"
#include "face_map.hpp"
#include "numvector.hpp"

struct Constants {
    // gaussian related
    Point2 mean;
    Point2 std_dev;
    double amp;
};
typedef void (*AnalyticFunction)(Constants, Point2, int, vector<Point3>&) ;
class AnalyticEvaluator: public SurfaceEvaluator {
    private:
        vector<Vector<Point3> >* _polynomial_patches;
        void (*_func)(Constants, Point2, int, vector<Point3>&);
        map<Index, int>* _binomial_coefs;
        int _poly_degree;
        int _num_patches;
    public: 
        //AnalyticEvaluator(string polynomial_filename, string mesh_filename);
        AnalyticEvaluator(string mesh_filename, void (*_func)(Constants, Point2, int, vector<Point3>&));
        ~AnalyticEvaluator();

        vector<Point3> evaluate(int patch_id, int evaluation_flags, Point2 uv_coords);
        void setup();
        void bounding_box(Point3& min, Point3& max);
        static AnalyticFunction select_function(int analytic_func);
};

void gaussian(Constants constants, Point2 uv, 
        int flags, vector<Point3>& functions_values);

void sin3(Constants constants, Point2 uv, 
        int flags, vector<Point3>& functions_values);

void sin_plus_sin(Constants constants, Point2 uv, int flags, vector<Point3>& values);
void exp_sin_cos(Constants constants, Point2 uv, int flags, vector<Point3>& values);
void torus(Constants constants, Point2 uv, int flags, vector<Point3>& values);
#endif
