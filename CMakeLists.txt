cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

project(Patchwork)

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

message(STATUS "env var " ${P4EST_DIR})
message(STATUS "p4est include dir" ${p4est_INCLUDE_DIR})
add_subdirectory(src)
add_subdirectory(tests)

add_library(patchwork STATIC ${Patchwork_SRC})
target_include_directories(patchwork PUBLIC 
    ${Patchwork_SOURCE_DIR}/include
    ${Blendsurf_INCLUDE_DIR}
    ${p4est_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_DIRS} )

# 

if(OpenMP_CXX_FOUND)
    target_link_libraries(patchwork OpenMP::OpenMP_CXX)
endif()
target_link_libraries(patchwork -L${Blendsurf_LIB} -lblendsurf)
target_link_libraries(patchwork -L${p4est_LIB} ${p4est_LIBRARIES})
target_link_libraries(patchwork ${ZLIB_LIBRARIES})
target_link_libraries(patchwork ${LAPACK_LIBRARIES} )
target_link_libraries(patchwork ${MPI_LIBRARIES})
target_link_libraries(patchwork stdc++)

add_subdirectory(utils)

    

