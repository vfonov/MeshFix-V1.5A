cmake_minimum_required (VERSION 2.8)
project(MeshFix C CXX)

ENABLE_TESTING()
INCLUDE(CTest)

add_subdirectory(contrib/JMeshLib)
add_subdirectory(3rd_party)

include_directories(
    include
    contrib/JMeshLib/include
    contrib/OpenNL3.2.1/src
    contrib/jrs_predicates
    3rd_party/superlu/SRC
)

link_directories(
    ${LINK_DIRECTORIES}
)

set(meshfix_h 
    include/component.h
    include/detectIntersections.h
    include/epsilonSampling.h
    include/exttrimesh.h
    include/holeFilling.h
    include/simplification.h
    include/sparseLSystem.h)


set(meshfix_src
    src/meshfix.cpp 
    contrib/jrs_predicates/jrs_predicates.c
    src/detectIntersections.cpp
    src/holeFilling.cpp
    src/sparseLSystem.cpp
    src/mwExtensions.cpp
    src/smoothing.cpp
    contrib/OpenNL3.2.1/src/nl_single_file.c
    src/epsilonSampling.cpp
    src/simplification.cpp
    src/uniform.cpp
    src/cleaning.cpp
)

if(CMAKE_CXX_SIZEOF_DATA_PTR EQUAL 8)
    ADD_DEFINITIONS(-DIS64BITPLATFORM)  ## needed for 64 bit 
endif()

ADD_DEFINITIONS(-DNL_USE_SUPERLU)
if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_EXE_LINKER_FLAGS "-Wl,--allow-multiple-definition")
endif()

message("BLAS_LIBRARY=${BLAS_LIBRARY}")
message("SUPERLU_LIBRARY=${SUPERLU_LIBRARY}")

add_executable(meshfix ${meshfix_src} ${meshfix_h} )
target_link_libraries(meshfix jmesh ${BLAS_LIBRARY} ${SUPERLU_LIBRARY} ) 

IF(BUILD_TESTING)
    ADD_SUBDIRECTORY(test)
ENDIF(BUILD_TESTING)


