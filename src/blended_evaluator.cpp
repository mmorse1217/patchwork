#include "blended_evaluator.hpp"
#include "face_map_obj.hpp"
#include <iomanip>

BlendedEvaluator::BlendedEvaluator(string wrl_filename){
    map<string,string> opts;
    opts.clear();
    ifstream fin("visoption3d"); assert(fin.good());
    string name;  fin>>name;
    while(fin.good()) {
        char cont[300];   fin.getline(cont, 299);
        opts[name] = string(cont);
        fin>>name;
    }
    fin.close();

    map<string,string>::iterator mi;

    int objtype;
    mi = opts.find("-objtype"); assert(mi!=opts.end());  { istringstream ss((*mi).second); ss>>objtype; }


    //------------------------
    char filea[300];	 mi = opts.find("-bdsurf_submatlibfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filea; }
    ifstream tmpa(filea);
    CCSubMatLib* submatlib = new CCSubMatLib();
    cerr<<"submatlib setup"<<endl;
    submatlib->setup(tmpa);
    char fileb[300];	 mi = opts.find("-bdsurf_bdulibfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>fileb; }
    ifstream tmpb(fileb);
    BdULib* bdulib = new BdULib();
    cerr<<"bdulib setup"<<endl;
    bdulib->setup(tmpb);
    char filec[300];	 mi = opts.find("-bdsurf_meshfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filec; }
    ifstream tmpc(wrl_filename.c_str());

    int cl;     mi = opts.find("-bdsurf_ctrllvl");   assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>cl; }// assert(cl>=0 && cl<=2);
    int pc;	    mi = opts.find("-bdsurf_pouctrl");   assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>pc; } assert(pc==0 || pc==1 || pc==2 || pc == 3);

    int ct;     mi = opts.find("-bdsurf_chttyp");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>ct; } 
    int bt;     mi = opts.find("-bdsurf_bsstyp");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>bt; }
    int pp;	    mi = opts.find("-bdsurf_stpppg");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>pp; }
    double lb;	mi = opts.find("-bdsurf_lb");        assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>lb; }
    double ub;	mi = opts.find("-bdsurf_ub");        assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>ub; }
    int gl;	    mi = opts.find("-bdsurf_indepboun"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>gl; }
    double fs;	mi = opts.find("-bdsurf_flats");     assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>fs; }
    int cot;    mi = opts.find("-bdsurf_concat");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>cot; }
    int sd;	    mi = opts.find("-bdsurf_spdeg");   assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>sd; }
    int poudeg; mi = opts.find("-bdsurf_poubsdeg");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>poudeg; }
    int lvl; mi = opts.find("-bdsurf_renderlvl");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>lvl; }
    int gen; mi = opts.find("-bdsurf_rendergen");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>gen; }
    int refinement_factor; mi = opts.find("-bdsurf_refinement_factor");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>refinement_factor; }
    int patch_order; mi = opts.find("-bdsurf_patch_order");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>patch_order; }
    bool  adaptive; mi = opts.find("-bdsurf_adaptive");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>adaptive; }
    double spacing; mi = opts.find("-bdsurf_interpolant_spacing");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>spacing; }
    bool interpolate; mi = opts.find("-bdsurf_interpolate");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>interpolate; }
    _blended_surface = new BdSurf();
    _blended_surface->setParams(cl,pc,ct,bt,pp, lb,ub, gl, fs, cot, sd ,poudeg,spacing, interpolate, refinement_factor,  0);
    //cerr<<"bdsurf setup "<< filec << endl;
    _blended_surface->setup(tmpc,0, Point3(1.0), (void*)submatlib, (void*)bdulib);

            _mesh = &_blended_surface->gpmesh();

}
vector<Point3> BlendedEvaluator::evaluate(int patch_id, int evaluation_flags, Point2 uv_coords){
    int local_vertex = 0;
    int global_vertex;
    int local_face;
    

    int return_value_size;

    if(evaluation_flags & BdSurf::EVAL_1ST_DERIV){
        return_value_size = 6;
    } else if (evaluation_flags & BdSurf::EVAL_VALUE){
        return_value_size = 2;
    } else {
        assert(0); // second derivatives not yet implementated
    }
    
    double vfcd_coords[return_value_size];
    double vxy_coords[return_value_size];
    vector<Point3> positions_and_derivs(return_value_size/2, Point3(0.));

    _blended_surface->Fcd2Vfcd(
            evaluation_flags,
            patch_id,
            uv_coords.array(), // coordinates in Fcd parameter space
            local_vertex,
            global_vertex,
            local_face,
            vfcd_coords);
    
    _blended_surface->Vfcd2Vxy(
            evaluation_flags,
            global_vertex,
            local_face,
            vfcd_coords,
            vxy_coords);
    
    _blended_surface->eval(
            evaluation_flags, 
            global_vertex,
            vxy_coords,
            positions_and_derivs.data());

    return positions_and_derivs;
}
void BlendedEvaluator::setup(){
    assert(0);
}

void BlendedEvaluator::bounding_box(Point3& min, Point3& max){
    _blended_surface->bbox(min, max);
}
