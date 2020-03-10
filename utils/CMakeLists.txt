find_library(GLUT_LIB
    NAMES  libglut.dylib
    PATHS $ENV{OPENGL_DIR}../lib
    STATIC
    )
find_library(OPENGLU_LIB
    NAMES  libGLU.dylib
    PATHS $ENV{OPENGL_DIR}../lib
    STATIC
    )
#move executable to build/bin
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
add_executable(find_patch_coefficients find_patch_coefficients.cpp)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lstdc++")

target_link_libraries(find_patch_coefficients 
    ${MPI_LIB}
    ${BLENDSURF_LIB} ${P4EST_LIB} ${SC_LIB} lapack png patchwork) 
add_executable(face_map_viewer face_map_viewer.cpp)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lstdc++")

target_link_libraries(face_map_viewer patchwork 
    ${OPENGL_LIB} ${OPENGLU_LIB} ${GLUT_LIB} png  ) 
