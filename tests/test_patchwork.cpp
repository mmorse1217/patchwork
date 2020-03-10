//#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "test_patchwork.hpp"
#include <mpi.h>
#include <p4est.h>
#include <sc.h>
using namespace std;


void optionsCreate(const char* optfile, map<string,string>& options)
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
}

void load_options(char* filec, int& refinement_factor, int& patch_order, 
        bool& adaptive, char* filepoly){
    map<string,string> opts;
    string opts_file = "visoption3d";
    optionsCreate(opts_file.c_str(), opts);
    map<string,string>::iterator mi;

     mi = opts.find("-bdsurf_meshfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filec; }
    ifstream tmpc(filec);
    mi = opts.find("-bdsurf_refinement_factor");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>refinement_factor; }
     mi = opts.find("-bdsurf_patch_order");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>patch_order; }
     mi = opts.find("-bdsurf_adaptive");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>adaptive; }

     mi = opts.find("-poly_coeffs_file"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filepoly; }

}

void load_blended_options(
        char* filea, char* fileb, char* filec, int& cl, int& pc,	int& ct,   
        int& bt,  int& pp, double& lb, double& ub,int& gl,double& fs, int& cot,  
        int& sd, int& poudeg, int& lvl, int& gen, int& refinement_factor, 
        int& patch_order, bool& adaptive, char* filepoly){
    map<string,string> opts;
    string opts_file = "visoption3d";
    optionsCreate(opts_file.c_str(), opts);
    map<string,string>::iterator mi;

    mi = opts.find("-bdsurf_submatlibfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filea; }
    mi = opts.find("-bdsurf_bdulibfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>fileb; }
    mi = opts.find("-bdsurf_meshfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filec; }
    mi = opts.find("-bdsurf_ctrllvl");   assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>cl; }// assert(cl>=0 && cl<=2);
    mi = opts.find("-bdsurf_pouctrl");   assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>pc; } assert(pc==0 || pc==1 || pc==2 || pc == 3);
    mi = opts.find("-bdsurf_chttyp");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>ct; } 
    mi = opts.find("-bdsurf_bsstyp");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>bt; }
    mi = opts.find("-bdsurf_stpppg");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>pp; }
    mi = opts.find("-bdsurf_lb");        assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>lb; }
    mi = opts.find("-bdsurf_ub");        assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>ub; }
    mi = opts.find("-bdsurf_indepboun"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>gl; }
    mi = opts.find("-bdsurf_flats");     assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>fs; }
    mi = opts.find("-bdsurf_concat");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>cot; }
    mi = opts.find("-bdsurf_spdeg");   assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>sd; }
    mi = opts.find("-bdsurf_poubsdeg");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>poudeg; }
    mi = opts.find("-bdsurf_renderlvl");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>lvl; }
    mi = opts.find("-bdsurf_rendergen");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>gen; }
    mi = opts.find("-bdsurf_refinement_factor");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>refinement_factor; }
    mi = opts.find("-bdsurf_patch_order");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>patch_order; }
    mi = opts.find("-bdsurf_adaptive");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>adaptive; }
    mi = opts.find("-poly_coeffs_file"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filepoly; }

}

int main( int argc, char**  argv ) {

    // Load options from options file

    cerr << "starting tests... " << endl;
    MPI_Init(NULL,NULL);

  //int mpiret = sc_MPI_Init ();
  //SC_CHECK_MPI (mpiret);
  //sc_MPI_Comm mpicomm = MPI_COMM_WORLD;

  /* These functions are optional.  If called they store the MPI rank as a
   * static variable so subsequent global p4est log messages are only issued
   * from processor zero.  Here we turn off most of the logging; see sc.h. */
  //sc_init (MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  //p4est_init (NULL, SC_LP_PRODUCTION);
    int result = Catch::Session().run( argc, argv );
    MPI_Finalize();
    cerr << "finished tests! " << endl;

    return result;
}
