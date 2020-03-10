#ifndef __POLYNOMIAL_EVALUATOR_HPP__
#define __POLYNOMIAL_EVALUATOR_HPP__

#include "surface_evaluator.hpp"
#include "face_map.hpp"
#include "numvector.hpp"
class PolynomialEvaluator: public SurfaceEvaluator {
    private:
        vector<Vector<Point3> >* _polynomial_patches;
        IntMatrix* _binomial_coefs;
        int _poly_degree;
        int _num_patches;
    public: 
        PolynomialEvaluator(string polynomial_filename, string mesh_filename);
        ~PolynomialEvaluator();

        vector<Point3> evaluate(int patch_id, int evaluation_flags, Point2 uv_coords);
        void setup();
        void bounding_box(Point3& min, Point3& max);
        Vector<Point3>* polynomial_patch(int F){
            assert(0 <= F && F < _polynomial_patches->size());
            return &_polynomial_patches->at(F);
        }
};

#endif
