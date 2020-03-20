find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

add_executable(find_patch_coefficients find_patch_coefficients.cpp)

target_link_libraries(find_patch_coefficients patchwork)
#${MPI_LIB}
#    ${BLENDSURF_LIB} ${P4EST_LIB} ${SC_LIB} ${LAPACK_LIBRARIES} patchwork) 
if (COMPILE_RENDERER)
add_executable(face_map_viewer face_map_viewer.cpp face_map_obj.cpp)
target_include_directories(face_map_viewer PUBLIC 
    ${Patchwork_SOURCE_DIR}/include
    #${Patchwork_SOURCE_DIR}/utils 
    ${Blendsurf_render_INCLUDE_DIR})

target_link_libraries(face_map_viewer ${Blendsurf_render_LIBRARIES}
    OpenGL::GL
    OpenGL::GLU
    GLUT::GLUT 
    patchwork
    )
endif()
