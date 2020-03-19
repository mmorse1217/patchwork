# - Try to find the p4est library
# Once done this will define
#
#  p4est_FOUND - system has p4est
#  p4est_INCLUDE_DIR - p4est include directory
#  p4est_LIB - p4est library directory
#  p4est_LIBRARIES - p4est libraries to link

if(p4est_FOUND)
    return()
endif()

find_path(p4est_INCLUDE_DIR p4est.h
    HINTS
        ${p4est_DIR}
        ${P4EST_DIR}
        $ENV{p4est_DIR}
        $ENV{P4EST_DIR}
        ENV P4EST_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/p4est
        ${CMAKE_SOURCE_DIR}/../p4est
        ${CMAKE_SOURCE_DIR}/../../p4est
        /usr
        /usr/local
    PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(p4est
    "\np4est not found --- You can download it using:\n\tgit clone https://github.com/cburstedde/p4est\n compiling the library and setting the P4EST_DIR environment variable accordingly"
    p4est_INCLUDE_DIR)
mark_as_advanced(p4est_INCLUDE_DIR)

set(p4est_LIB "${p4est_INCLUDE_DIR}/../lib/")
set(p4est_LIBRARIES "${p4est_INCLUDE_DIR}/../lib/libp4est.a" "${p4est_INCLUDE_DIR}/../lib/libsc.a" )
list(APPEND CMAKE_MODULE_PATH "${p4est_INCLUDE_DIR}/../cmake")



