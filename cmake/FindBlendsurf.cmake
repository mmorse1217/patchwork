# - Try to find the Blendsurf library
# Once done this will define
#
#  Blendsurf_FOUND - system has Blendsurf
#  Blendsurf_INCLUDE_DIR - Blendsurf include directory
#  Blendsurf_LIB - Blendsurf library directory
#  Blendsurf_LIBRARIES - Blendsurf libraries to link

if(Blendsurf_FOUND)
    return()
endif()

find_path(Blendsurf_INCLUDE_DIR bdsurf.hpp
    HINTS
        ${Blendsurf_DIR}
        ${BLENDSURF_DIR}
        $ENV{Blendsurf_DIR}
        $ENV{BLENDSURF_DIR}
        ENV BLENDSURF_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/blendsurf
        ${CMAKE_SOURCE_DIR}/../blendsurf
        ${CMAKE_SOURCE_DIR}/../../blendsurf
        /usr
        /usr/local
    PATH_SUFFIXES include
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Blendsurf
    "\nBlendsurf not found --- You can download it using:\n\tgit clone https://github.com/mmorse1217/blendsurf\n and setting the BLENDSURF_DIR environment variable accordingly"
    Blendsurf_INCLUDE_DIR)
mark_as_advanced(Blendsurf_INCLUDE_DIR)

set(Blendsurf_LIB "${Blendsurf_INCLUDE_DIR}/../lib/"  )
set(Blendsurf_LIBRARIES "${Blendsurf_INCLUDE_DIR}/../lib/libblendsurf.a"  )
set(Blendsurf_render_LIBRARIES "${Blendsurf_INCLUDE_DIR}/../lib/libblendsurf_render.a"  )
set(Blendsurf_render_INCLUDE_DIR "${Blendsurf_INCLUDE_DIR}/../vis/"  )
list(APPEND CMAKE_MODULE_PATH "${Blendsurf_INCLUDE_DIR}/../cmake")
