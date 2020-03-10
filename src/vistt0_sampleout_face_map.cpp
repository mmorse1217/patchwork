#include "ccsubmatlib.hpp"
#include "ccsurf.hpp"
#include "face_map.hpp"
#include "blended_evaluator.hpp"
#include "polynomial_evaluator.hpp"

#include "ballviewer.hpp"
#include "ccsurfobj.hpp"
#include "bdsurfobj.hpp"
#include "face_map_obj.hpp"

int WriteAndQuit = 0;

int optionsCreate(const char* optfile, map<string,string>& options)
{
    options.clear();
    ifstream fin(optfile); assert(fin.good());
    string name;  fin>>name;
    while(fin.good()) {
        char cont[300];   fin.getline(cont, 299);
        options[name] = string(cont);
        fin>>name;
    }
    fin.close();
    return 0;
}





// ---------------------------------------------------------------------- 
int main(int argc, char** argv)
{


    assert( argc == 2 );

    map<string,string> opts;
    optionsCreate(argv[1], opts);

    map<string,string>::iterator mi;

    int objtype;
    mi = opts.find("-objtype"); assert(mi!=opts.end());  { istringstream ss((*mi).second); ss>>objtype; }
    GeoObject* obj = NULL;


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
    ifstream tmpc(filec);

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
    double fit_accuracy; mi = opts.find("-bdsurf_fit_accuracy");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>adaptive; }




    char filepoly[300];	 mi = opts.find("-poly_coeffs_file"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filepoly; }
    ifstream tmppoly(filepoly);

    BdSurf* bdsurf = new BdSurf();
    cout << "refinement factor: " << refinement_factor << endl;
    //bdsurf->set_params(cl,pc,ct,bt,pp, lb,ub, gl, fs, cot, sd ,poudeg,refinement_factor,patch_order, adaptive, 0);
    //bdsurf->setup(tmpc,0,Point3(1.0),(void*)submatlib,(void*)bdulib);
    //ofstream fout("bddump.dat"); bdsurf->dump(fout); fout.close();
    
    //PolynomialEvaluator* polynomial_evaluator = new PolynomialEvaluator(string(filepoly), string(filec));
    //FaceMapSurf* face_map = new FaceMapSurf(polynomial_evaluator);
    
    bdsurf->setParams(cl,pc,ct,bt,pp, lb,ub, gl, fs, cot, sd ,poudeg,.01, 0);
    cerr<<"bdsurf setup "<< filec << endl;
    bdsurf->setup(tmpc,0, Point3(1.0), (void*)submatlib, (void*)bdulib);
    BlendedEvaluator* blended_evaluator = new BlendedEvaluator(bdsurf);
    FaceMapSurf* face_map = new FaceMapSurf(blended_evaluator);

    face_map->set_params(refinement_factor, patch_order, adaptive, fit_accuracy);
    face_map->setup();

    int ctrlsize = pow2(2);
    int ss = (int)ceil( 1.0 * (ctrlsize) ) + 1;
    double ctrlstep = 1.0/ctrlsize;




    // ---------------------------------------------------
    ofstream fout_int("sample_pts_interior.dat");
    fout_int.precision(16);
    for (int K=3; K<=13; K++){
        cerr << "int" << K << endl;
        fout_int<<K<<endl;            
        for(int f=0; f<K; f++) {
            for(int j=0; j<ss; j++){
                for(int i=0; i<ss; i++) {
                    if(i==0 && j==0 && f!=0) continue; 
                    if(i==0 && j>0) continue;
                    double cd[2];             
                    cd[0] = i*ctrlstep;  cd[1] = j*ctrlstep; //V position
                    double xy[2];    
                    xy[0] = cd[0];
                    xy[1] = cd[1];
                    //bdsurf->Vfcd2Vxy_val(1, K,0, GpMesh::INTERIOR_VERTEX, f, cd, xy); 
                    fout_int<<xy[0]<<" "<<xy[1]<<endl;
                }
            }
        }
    }
    fout_int.close();
    // ---------------------------------------------------
    /*
       ofstream fout_boun("sample_pts_boun.dat");
       fout_boun.precision(16);
       for (int K=2; K<=13; K++){
       cerr << "boun" << K << endl;
       fout_boun<<K<<endl;            
       for(int f=0; f<K; f++) {
       for(int j=0; j<ss; j++){
       for(int i=0; i<ss; i++) {
       if(i==0 && j==0 && f!=0) continue; 
       if(i==0 && j>0) continue;
       double cd[2];             
       cd[0] = i*ctrlstep;  cd[1] = j*ctrlstep; //V position
       double xy[2];    
       bdsurf->Vfcd2Vxy_val(1, K,1, GpMesh::CREASE_VERTEX, f, cd, xy); 
       fout_boun<<xy[0]<<" "<<xy[1]<<endl;                    
       }
       }
       }
       int f = K-1;
       int i=0;
       for (int j=1; j<ss; j++){    
       double cd[2];             
       cd[0] = i*ctrlstep;  cd[1] = j*ctrlstep; //V position
       double xy[2];          
       bdsurf->Vfcd2Vxy_val(1, K, 1, GpMesh::CREASE_VERTEX, f, cd, xy); 
       fout_boun<<xy[0]<<" "<<xy[1]<<endl;
       }
       }
       fout_boun.close();
    // ---------------------------------------------------
    ofstream fout_vex("sample_pts_convex.dat");
    fout_vex.precision(16);
    for (int K=1; K<=8; K++){
    cerr << "convex" << K << endl;
    fout_vex<<K<<endl;            	 
    for(int f=0; f<K; f++) {
    for(int j=0; j<ss; j++){
    for(int i=0; i<ss; i++) {
    if(i==0 && j==0 && f!=0) continue; 
    if(i==0 && j>0) continue;
    double cd[2];             
    cd[0] = i*ctrlstep;  cd[1] = j*ctrlstep; //V position
    double xy[2];    
    bdsurf->Vfcd2Vxy_val(1, K,1, GpMesh::CONVEX_VERTEX, f, cd, xy); 
    fout_vex<<xy[0]<<" "<<xy[1]<<endl;
    }
    }
    }
    int f = K-1;
    int i=0;
    for (int j=1; j<ss; j++){    
    double cd[2];             
    cd[0] = i*ctrlstep;  cd[1] = j*ctrlstep; //V position
    double xy[2];    
    bdsurf->Vfcd2Vxy_val(1, K, 1, GpMesh::CONVEX_VERTEX, f, cd, xy); 
    fout_vex<<xy[0]<<" "<<xy[1]<<endl;
    }
    }
    fout_vex.close();
    // ---------------------------------------------------
    ofstream fout_cave("sample_pts_concave.dat");
    fout_cave.precision(16);

    for (int K=2; K<=8; K++){
    cerr << "concave" << K << endl;
    fout_cave<<K<<endl;            
    for(int f=0; f<K; f++) {
    for(int j=0; j<ss; j++){
    for(int i=0; i<ss; i++) {
    if(i==0 && j==0 && f!=0) continue; 
    if(i==0 && j>0) continue;
    double cd[2];             
    cd[0] = i*ctrlstep;  cd[1] = j*ctrlstep; //V position
    double xy[2];    
    bdsurf->Vfcd2Vxy_val(1, K,1, GpMesh::CONCAVE_VERTEX, f, cd, xy); 
    fout_cave<<xy[0]<<" "<<xy[1]<<endl;                    
}
}
}
int f = K-1;
int i=0;
for (int j=1; j<ss; j++){    
    double cd[2];             
    cd[0] = i*ctrlstep;  cd[1] = j*ctrlstep; //V position
    double xy[2];    
    bdsurf->Vfcd2Vxy_val(1, K, 1, GpMesh::CONCAVE_VERTEX, f, cd, xy); 
    fout_cave<<xy[0]<<" "<<xy[1]<<endl;                    
}
}
fout_cave.close();
*/

cerr<<"bdsurfobj setup"<<endl;
char filetex[100];
mi = opts.find("-bdsurfobj_texfile"); 
assert(mi!=opts.end());  
{ istringstream ss((*mi).second); ss>>filetex; }
//int lvl = 2; 
//int gen = 0;
int alt = 0;
obj = new FaceMapSurfObj(face_map, lvl, gen, alt);
face_map->evalall(lvl, gen);
//2. viewer

char vrname[100];  mi = opts.find("-vrname"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>vrname; }
cout << "setup viewer" << endl;
Viewer::initGL(&argc, argv);
BallViewer *ballviewer = new BallViewer(vrname,600,600);
ballviewer->setObject(obj);
cout << "finished viewer setup " << endl;

glutMainLoop();
glCheck();


}
