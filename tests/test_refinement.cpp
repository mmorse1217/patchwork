
#include "catch.hpp"
#include "test_patchwork.hpp"
#include <face_map.hpp>
#include <polynomial_evaluator.hpp>
#include <blended_evaluator.hpp>
TEST_CASE("test parallel face map setup", "[face-map][setup]"){
    BlendedEvaluator* blended_evaluator = 
        new BlendedEvaluator("wrl_meshes/wrl/cube.wrl");

    FaceMapSurf* face_map = new FaceMapSurf(blended_evaluator);
    int adaptive = 1;
    int patch_order = 6;
    int refinement_factor = 2;
    face_map->set_params(refinement_factor, patch_order, adaptive);
    face_map->setup();

}
