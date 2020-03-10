#ifndef __EVALUATOR_HPP__
#define __EVALUATOR_HPP__

#include <vec2t.hpp>
#include <vec3t.hpp>
#include <gpmesh.hpp>
#include "function.hpp"

// virtual class describing the interface needed to construct a face-map surface
// from another mesh- or patch-based surface
class SurfaceEvaluator: public Function {
    protected:
        GpMesh* _mesh;
    public:
        virtual ~SurfaceEvaluator(){;} // lol wut

        virtual vector<Point3> evaluate(int patch_id, int evaluation_flags, Point2 uv_coords)=0;
        virtual void setup()=0;
        virtual void bounding_box(Point3& min, Point3& max)=0;
        const GpMesh* mesh(){
            return _mesh;
        }

};
#endif
