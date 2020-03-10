#include "catch.hpp"
#include <polynomial_evaluator.hpp>
#include <blended_evaluator.hpp>
#include <sampling.hpp>
#include <face_map.hpp>
#include "test_patchwork.hpp"

using Sampling::sample_2d;
using Sampling::base_domain;
using Sampling::equispaced;

Point3 first_derivative_finite_diff(Point3 x_minus_h, Point3 x_plus_h, double h){
    return (x_plus_h - x_minus_h)/(2.*h);
}
Point3 second_derivative_finite_diff(Point3 x_minus_h, Point3 x, Point3 x_plus_h, double h){
    return (x_plus_h- 2.*x + x_minus_h )/(h*h);
}

TEST_CASE("Test second derivative patch evaluation known solution", "[derivative]"){


    SECTION("test flat patch has small second derivatives."){
        char filec[300];
        int refinement_factor;
        int patch_order;
        bool  adaptive;
        char filepoly[300];

        load_options(filec, refinement_factor, patch_order, adaptive, filepoly);
        PolynomialEvaluator* polynomial_evaluator = new PolynomialEvaluator(string(filepoly), string(filec));
        FaceMapSurf* face_map = new FaceMapSurf(polynomial_evaluator);

        face_map->set_params(refinement_factor, patch_order, adaptive, 1e-4);
        face_map->setup();
        int num_samples = 6;
        double step = 1./double(num_samples-1);
        for(int i = 0; i < num_samples; i++){
            for(int j = 0; j < num_samples; j++){
                double xy[2];
                xy[0] = i*step;
                xy[1] = j*step;
                Point3 results[6];
                face_map->eval(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV, 0, xy, results);
                //for(int di =3; di < 6; di++)
                //cout << results[di] << endl;
                Point3 Xuu = results[3];
                Point3 Xuv = results[4];
                Point3 Xvv = results[5];
                CHECK(Xuu.length() <=1e-12);
                CHECK(Xuv.length() <=1e-12);
                CHECK(Xvv.length() <=1e-12);
            }
        }
    }
}
TEST_CASE("Test second derivative patch evaluation approx", "[derivative]"){


    char filea[300]; char fileb[300]; char filec[300]; int cl; int pc;	int ct;   
    int bt;  int pp; double lb; double ub;int gl;double fs; int cot;  
    int sd; int poudeg; int lvl; int gen; int refinement_factor; 
    int patch_order; bool adaptive; char filepoly[300];

    load_blended_options(
            filea, fileb, filec, cl, pc,	ct,   
            bt, pp, lb, ub,gl,fs, cot,  
            sd,  poudeg, lvl, gen, refinement_factor, 
            patch_order, adaptive, filepoly);
    
    ifstream tmpa(filea);
    CCSubMatLib* submatlib = new CCSubMatLib();
    submatlib->setup(tmpa);

    ifstream tmpb(fileb);
    BdULib* bdulib = new BdULib();
    bdulib->setup(tmpb);
    ifstream tmpc(string("../wrl_files/cube.wrl").c_str());
    
    BdSurf* bdsurf = new BdSurf();
    bdsurf->setParams(cl,pc,ct,bt,pp, lb,ub, gl, fs, cot, sd ,poudeg, refinement_factor, 0);
    bdsurf->setup(tmpc,0, Point3(1.0), (void*)submatlib, (void*)bdulib);
    
    tmpa.close();
    tmpb.close();
    tmpc.close();
    BlendedEvaluator* blended_evaluator = new BlendedEvaluator(bdsurf);
    FaceMapSurf* face_map = new FaceMapSurf(blended_evaluator);
    face_map->set_params(0, 6, 0, 1e-4);
    face_map->_legacy = true;
    face_map->setup();
    SECTION("Test second derivatives on patches approximating a blended surface"){
        double eps = 1e-3;
        int num_samples = 10;
        NumMatrix samples = sample_2d<equispaced>(num_samples,base_domain);
        for(int pi =0; pi < face_map->_patches.size(); pi++){
            for(int i =1; i < samples.n()-1; i++){
                Point2 s(samples.clmdata(i));

                vector<Point3> position_and_derivs(6,Point3());
                
                cout << i << endl;
                // evaluate second derivative at sample
                face_map->eval(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV,
                        pi, s.array(), position_and_derivs.data());
                
                // eval with finite differences
                Point2 s_plus_eps_x(s.x()+eps, s.y());
                Point2 s_minus_eps_x (s.x()-eps, s.y());
                Point2 s_plus_eps_y(s.x(), s.y()+eps);
                Point2 s_minus_eps_y(s.x(), s.y()-eps);

                Point3 p = position_and_derivs[0];
                Point3 p_plus_eps_x  ;
                Point3 p_minus_eps_x ;
                Point3 p_plus_eps_y  ;
                Point3 p_minus_eps_y ;

                face_map->eval(EVAL_VALUE,
                        pi, s_plus_eps_x.array(), &p_plus_eps_x);
                face_map->eval(EVAL_VALUE,
                        pi, s_minus_eps_x.array(), &p_minus_eps_x);
                face_map->eval(EVAL_VALUE,
                        pi, s_plus_eps_y.array(), &p_plus_eps_y);
                face_map->eval(EVAL_VALUE,
                        pi, s_minus_eps_y.array(), &p_minus_eps_y);

                Point3 dx = first_derivative_finite_diff(p_minus_eps_x, p_plus_eps_x, eps);
                Point3 dy = first_derivative_finite_diff(p_minus_eps_y, p_plus_eps_y, eps);
                Point3 ddx = second_derivative_finite_diff(p_minus_eps_x, p, p_plus_eps_x, eps);
                Point3 ddy = second_derivative_finite_diff(p_minus_eps_y, p, p_plus_eps_y, eps);

                CHECK((dx - position_and_derivs[1]).length() <= 1e-4);
                CHECK((dy - position_and_derivs[2]).length() <= 1e-4);
                CHECK((ddx - position_and_derivs[3]).length() <= 1e-4);
                CHECK((ddy - position_and_derivs[5]).length() <= 1e-4);

            }

        }

    delete face_map;
    }

}
