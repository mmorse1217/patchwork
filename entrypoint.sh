#!/bin/bash

# code to execute during CI build
mkdir -p /patchwork/build/ 
cd /patchwork/build/
export BLENDSURF_DIR=/libs/blendsurf 
export P4EST_DIR=/libs/p4est-1.1
cmake -DCMAKE_MODULE_PATH=/usr/share/cmake-3.10/Modules/ -DCOMPILE_RENDERER=True ..  
make 
cd ../
bin/test_patchwork [bezier],[diff-geom],[derivative]
