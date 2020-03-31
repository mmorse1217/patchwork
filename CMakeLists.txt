cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)


project(Patchwork)
set(Patchwork_VERSION 0.1)
list(APPEND CMAKE_MODULE_PATH "${Patchwork_SOURCE_DIR}/cmake")

FUNCTION(PREPEND var )
   SET(listVar "")
   FOREACH(f ${${var}})
       LIST(APPEND listVar "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${Patchwork_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${Patchwork_SOURCE_DIR}/lib")

find_package(OpenMP REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Blendsurf REQUIRED)
find_package(MPI REQUIRED)
find_package(p4est REQUIRED)
find_package(ZLIB REQUIRED)

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(tests)

add_library(patchwork STATIC ${Patchwork_SRC})
target_include_directories(patchwork PUBLIC 
    ${Blendsurf_INCLUDE_DIR}
    ${p4est_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_DIRS} 
    )

# 

if(OpenMP_CXX_FOUND)
    target_link_libraries(patchwork OpenMP::OpenMP_CXX)
endif()

target_link_libraries(patchwork ${Blendsurf_LIBRARIES})
target_link_libraries(patchwork -L${p4est_LIB} ${p4est_LIBRARIES})
target_link_libraries(patchwork ZLIB::ZLIB)
target_link_libraries(patchwork ${LAPACK_LIBRARIES} )
target_link_libraries(patchwork ${MPI_LIBRARIES})
target_link_libraries(patchwork stdc++)

add_subdirectory(utils)



set(Patchwork_INCLUDE_DEST "include/patchwork-${Patchwork_VERSION}")
set(Patchwork_LIB_DEST "lib/patchwork-${Patchwork_VERSION}")


# generator expressions are needed for the include directories, since installing headers changes the include path
target_include_directories(patchwork PUBLIC
   # including headers when building from source
   $<BUILD_INTERFACE:${Patchwork_SOURCE_DIR}/include> 
   $<BUILD_INTERFACE:${Patchwork_BINARY_DIR}/include> 
   # including headers when installing 
   $<INSTALL_INTERFACE:include> 
   )
target_link_libraries(patchwork ${LAPACK_LIBRARIES})
target_link_libraries(patchwork stdc++)

# compile the renderer
if (COMPILE_RENDERER)
    add_subdirectory(vis)
    # TODO add install support for renderer...?
endif()

# Install patchwork in CMAKE_INSTALL_PREFIX (defaults to /usr/local on linux)
install(
    # targets to install
    TARGETS patchwork 
    # name of the CMake "export group" containing the targets we want to install
    EXPORT patchworkTargets
    # Dynamic, static library and include destination locations after running
    # "make install"
    LIBRARY DESTINATION ${Patchwork_LIB_DEST}
    ARCHIVE DESTINATION ${Patchwork_LIB_DEST} 
    INCLUDES DESTINATION ${Patchwork_INCLUDE_DEST}
    )
install(
    # export group we want to save (matches name defined above containing the
    # install targets)
    EXPORT patchworkTargets
    # CMake file in which to store the export group info
    FILE  patchworkTargets.cmake
    # namespace prepends all targets in the export group
    NAMESPACE Patchwork::
    # where to place the resulting file (here we're putting it with the library)
    DESTINATION ${Patchwork_LIB_DEST}
    )

# Copies the files in ${Patchwork_INC} to ${Patchwork_INCLUDE_DEST}
install(FILES ${Patchwork_INC} DESTINATION ${Patchwork_INCLUDE_DEST})

# Write a "version file" in case some wants to only load a particular version of
# patchwork @BUG this fails at the moment.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    patchworkConfigVersion.cmake
    VERSION ${Patchwork_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

# Copies the resulting CMake config files to the installed library directory
install(FILES 
    "cmake/patchwork-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/patchworkConfigVersion.cmake"
    DESTINATION ${Patchwork_LIB_DEST}
    )

    

