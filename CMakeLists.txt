cmake_minimum_required (VERSION 2.6)
project (MeshFix)
add_subdirectory(contrib/JMeshLib)

# Commented out in OSX
#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET(BUILD_SHARED_LIBRARIES OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

# OSX: Downloaded and compiled superlu from the official website
add_library(superlu STATIC IMPORTED SHARED)
find_library(SUPERLU_LIBRARY superlu /usr/lib)
MESSAGE(STATUS "SUPERLU: " ${SUPERLU_LIBRARY})

# NOTE: installed the BLAS library from homebrew
add_library(blas STATIC IMPORTED)
find_library(BLAS_LIBRARY blas) # usr/lib) Commented out in OSX
MESSAGE(STATUS "BLAS: " ${BLAS_LIBRARY})


include_directories(
    include
    /usr/include/superlu
    contrib/JMeshLib/include
    contrib/OpenNL3.2.1/src
    contrib/jrs_predicates
)

link_directories(
    ${LINK_DIRECTORIES}
    ${CMAKE_CURRENT_SOURCE_DIR}/contrib/JMeshLib/lib
)
file(GLOB meshfix_h include/*.h)
set(meshfix_src
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
add_custom_command(
    TARGET meshfix
    PRE_LINK
    COMMAND make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/contrib/JMeshLib/
)
if(CMAKE_CXX_SIZEOF_DATA_PTR EQUAL 8)
    ADD_DEFINITIONS(-DIS64BITPLATFORM)  ## needed for 64 bit 
endif()
ADD_DEFINITIONS(-DNL_USE_SUPERLU)
if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_EXE_LINKER_FLAGS "-Wl,--allow-multiple-definition")
endif()
add_executable(meshfix ${meshfix_src} ${meshfix_h} meshfix.cpp)
target_link_libraries(meshfix jmesh ${SUPERLU_LIBRARY} ${BLAS_LIBRARY}) 
