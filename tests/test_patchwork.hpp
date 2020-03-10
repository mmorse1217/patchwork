#ifndef __TEST_FACE_MAP__
#define __TEST_FACE_MAP__
#include <iostream>
#include <map>
#include <fstream>
#include <assert.h>
using namespace std;
void optionsCreate(const char* optfile, map<string,string>& options);

void load_options(char* filec, int& refinement_factor, int& patch_order, 
        bool& adaptive, char* filepoly);
void load_blended_options(
        char* filea, char* fileb, char* filec, int& cl, int& pc,	int& ct,   
        int& bt,  int& pp, double& lb, double& ub,int& gl,double& fs, int& cot,  
        int& sd, int& poudeg, int& lvl, int& gen, int& refinement_factor, 
        int& patch_order, bool& adaptive, char* filepoly);

#endif
