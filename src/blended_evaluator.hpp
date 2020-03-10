#ifndef __BLENDED_EVALUATOR_HPP__
#define __BLENDED_EVALUATOR_HPP__
#include "surface_evaluator.hpp"
#include <bdsurf.hpp>

class BlendedEvaluator: public SurfaceEvaluator {
    public: 
        BdSurf* _blended_surface;
        BlendedEvaluator(string wrl_filename);
        BlendedEvaluator(BdSurf* bdsurf){
            _blended_surface = bdsurf;
            _mesh = &bdsurf->gpmesh();
            }
        ~BlendedEvaluator(){
            delete _blended_surface;
        }

        vector<Point3> evaluate(int patch_id, int evaluation_flags, Point2 uv_coords);
        void setup();
        void bounding_box(Point3& min, Point3& max);
};

#endif
