cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)

#sets compiler to match $HEDGEHOG_CC/CXX
execute_process(COMMAND which $ENV{HEDGEHOG_CC} OUTPUT_VARIABLE cc)
execute_process(COMMAND which $ENV{HEDGEHOG_CXX} OUTPUT_VARIABLE cpp)
message(STATUS "$ENV{CC}")
set(ENV{CC}  $ENV{HEDGEHOG_CC} )
set(ENV{CXX} $ENV{HEDGEHOG_CXX})
message(STATUS "$ENV{CC}")

project(Patchwork)

find_library(MPI_LIB
    NAMES  libmpi.dylib libmpich.dylib
    PATHS $ENV{PETSC_DIR}/$ENV{PETSC_ARCH} $ENV{MPI_HOME} 
    PATH_SUFFIXES lib
    SHARED
    )
find_library(P4EST_LIB
    NAMES  libp4est.dylib 
    PATHS $ENV{P4EST_DIR}
    PATH_SUFFIXES lib
    SHARED
    )
find_library(SC_LIB
    NAMES  libsc.dylib 
    PATHS $ENV{P4EST_DIR}
    PATH_SUFFIXES lib
    SHARED
    )
find_library(BLENDSURF_LIB
    NAMES  libblend.a 
    PATHS $ENV{BLENDSURF_DIR}
    STATIC
    )
find_library(OPENGL_LIB
    NAMES  libGL.dylib
    PATHS $ENV{OPENGL_DIR}../lib
    STATIC
    )
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
find_package(lapack REQUIRED)
find_package(blas REQUIRED)
find_package(png REQUIRED)

include_directories($ENV{BLENDSURF_DIR})
# export MPI_DIR=$PETSC_DIR/$PETSC_ARCH
include_directories($ENV{MPI_DIR}/include)
include_directories($ENV{P4EST_DIR}/include)
# export OPENGL_DIR=/usr/X11R6/include on mac
include_directories($ENV{OPENGL_DIR})
include_directories(${Patchwork_SOURCE_DIR}/src)
include_directories(${Patchwork_SOURCE_DIR}/tests)


add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(utils)

    

