# ===================================================================
# System configuration file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 05 June 2015
# normally no need to alter any thing in this file
# please DO NOT modify unless if it is necessary
# beginners are advised NOT to alter anything the below contents
# ========================================================================

cmake_minimum_required (VERSION 2.8)

project (ParMooN)

include(${CMAKE_SOURCE_DIR}/UserConfig.cmake)

set (ParMooN_VERSION_MAJOR 1)
set (ParMooN_VERSION_MINOR 0)

# set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/lib)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} $ENV{HOME}/lib)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} $ENV{HOME}/software/lib)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} $ENV{HOME}/Software/lib)

set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /opt/local/include)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /opt/include)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} $ENV{HOME}/include)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} $ENV{HOME}/software/include)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} $ENV{HOME}/Software/include)

# set(ENV{PKG_CONFIG_PATH} $ENV{PKG_CONFIG_PATH}:/opt/local/lib/pkgconfig )
set(PARMOON_EXTLIB_PATH ${PROJECT_SOURCE_DIR}/EXT_LIB)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


# search the blas in the system, if available it will be default
message("==================Searching BLAS in the System - begin===================================================")
find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_DIRS})
if(CMAKE_VERBOSE_MAKEFILE)
  if(BLAS_FOUND)
    message("Available system BLAS_LIB =${BLAS_LIBRARIES}")
  elseif(NOT BLAS_FOUND)
      message(FATAL_ERROR "BLAS_NOT_FOUND in your system")
  endif()
endif()
message("==================Searching BLAS in the System - end=====================================================")

# search the LAPACK in the system, if available it will be default
message("==================Searching LAPACK in the System - begin===================================================")
find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIRS})
if(CMAKE_VERBOSE_MAKEFILE)
  if(LAPACK_FOUND)
    message("Available system LAPACK_LIB =${LAPACK_LIBRARIES}")
  elseif(NOT LAPACK_FOUND)
      message(FATAL_ERROR "LAPACK_NOT_FOUND in your system")
  endif()
endif()
message("==================Searching LAPACK in the System - end=====================================================")


# search the UMFPACK in the system, if available it wil be default
message("==================Searching UMFPACK in the System - begin===================================================")
find_package(UMFPACK COMPONENTS SCALAPACK LAPACK)
include_directories(${UMFPACK_INCLUDE_DIRS})
if(CMAKE_VERBOSE_MAKEFILE)
  if(UMFPACK_FOUND)
    message("Available system UMFPACK_LIB =${UMFPACK_LIBRARIES}")
  elseif(NOT UMFPACK_FOUND)
    message(FATAL_ERROR "UMFPACK_NOT_FOUND in your system or in the ParMooN")
  endif()
endif()
message("==================Searching UMFPACK in the System - end=====================================================")


# search the MUMPS in the system, if available it wil be default
message("==================Searching MUMPS in the System - begin=====================================================")
find_package(MUMPS)
include_directories(${MUMPS_INCLUDE_DIRS})
if(CMAKE_VERBOSE_MAKEFILE)
  if(MUMPS_FOUND)
    message("Available system MUMPS_LIB =${MUMPS_LIBRARIES}")
  elseif(NOT MUMPS_FOUND)
    message(FATAL_ERROR "MUMPS_NOT_FOUND in your system or in the ParMooN")
  endif()
endif()
message("==================Searching MUMPS in the System - end=======================================================")


# search the TRIANGLE in the system, if available it wil be default
message("==================Searching TRIANGLE in the System - begin==================================================")
find_package(TRIANGLE)
include_directories(${TRIANGLE_INCLUDE_DIRS})
if(CMAKE_VERBOSE_MAKEFILE)
  if(TRIANGLE_FOUND)
    message("Available system TRIANGLE_LIB =${TRIANGLE_LIBRARIES}")
  elseif(NOT TRIANGLE_FOUND)
    message(FATAL_ERROR "TRIANGLE_NOT_FOUND in your system or in the ParMooN")
  endif()
endif()
message("==================Searching TRIANGLE in the System - end====================================================")

# search the TETGEN in the system, if available it wil be default
if("${GEO}" STREQUAL "3D")
  message("==================Searching TETGEN in the System - begin==================================================")
  find_package(TETGEN)
  include_directories(${TETGEN_INCLUDE_DIRS})
  if(CMAKE_VERBOSE_MAKEFILE)
    if(TETGEN_FOUND)
      message("Available system TETGEN_LIB =${TETGEN_LIBRARIES}")
    elseif(NOT TETGEN_FOUND)
      message(FATAL_ERROR "TETGEN_NOT_FOUND in your system or in the ParMooN")
    endif()
  endif()
message("==================Searching TETGEN in the System - end======================================================")
endif("${GEO}" STREQUAL "3D")

if (NOT CONFIGURED_ONCE)
    set(CMAKE_CXX_FLAGS "${warnings}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_C_FLAGS   "${warnings}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif()


# select the compiler type,
# CMAKE_BUILD_TYPE [DEBUG|RELEASE|RELWITHDEBINFO|MINSIZEREL]

set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "options")
set(PARMOON_INCLUDE )

include (CMakeForceCompiler)
  set(CMAKE_SYSTEM_NAME ${ARCH})
  
if("${ARCH}" STREQUAL "MAC64")
  CMAKE_FORCE_C_COMPILER   (clang Clang)
  CMAKE_FORCE_CXX_COMPILER (clang++ Clang)
  set(CMAKE_C_FLAGS "-O3 -fapple-pragma-pack")
  set(CMAKE_C_FLAGS_DEBUG "-g -fapple-pragma-pack")
  set(CMAKE_C_FLAGS_RELEASE "O3 -fapple-pragma-pack")
  set(CMAKE_CXX_FLAGS "-O3 ${PARMOON_CXX_DEF} -DREDUCED -DNO_TIMER -DMKL_ILP64 -m64 -fapple-pragma-pack")
  set(CMAKE_CXX_FLAGS_DEBUG "${PARMOON_CXX_DEF} -DREDUCED -DNO_TIMER -DMKL_ILP64 -m64 -fapple-pragma-pack")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${PARMOON_CXX_DEF} -DREDUCED -DNO_TIMER -DMKL_ILP64 -m64 -fapple-pragma-pack") 
elseif("${ARCH}" STREQUAL "LINUX64")
  if("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
    CMAKE_FORCE_C_COMPILER   (mpicc)
    CMAKE_FORCE_CXX_COMPILER (mpicxx)
  elseif("${PARALLEL_TYPE}" STREQUAL "SEQUENTIAL"  OR "${PARALLEL_TYPE}" STREQUAL "OMPONLY")
    CMAKE_FORCE_C_COMPILER   (gcc)
    CMAKE_FORCE_CXX_COMPILER (g++)
  endif()  
  
  set(CMAKE_C_FLAGS "-D__$(GEO)__  -std=c++11 -fopenmp")
  set(CMAKE_C_FLAGS_DEBUG "-g -D__$(GEO)__ $(INC) -std=c++11 -fopenmp")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -D__$(GEO)__ -s $(INC) -std=c++11 -fopenmp")
  set(CMAKE_CXX_FLAGS "-D__$(GEO)__ -DTRILIBRARY -DREDUCED -DNO_TIMER $(INC) -DMKL_ILP64 -m64 -fopenmp")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -D__$(GEO)__ -DTRILIBRARY -DREDUCED -DNO_TIMER $(INC) -DMKL_ILP64 -m64  -fopenmp")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -D__$(GEO)__ -s -DTRILIBRARY -DREDUCED -DNO_TIMER $(INC) -DMKL_ILP64 -m64 -fopenmp") 
  
endif()




include_directories("${CMAKE_SOURCE_DIR}/include/AMG")
include_directories("${CMAKE_SOURCE_DIR}/include/Refinement")
include_directories("${CMAKE_SOURCE_DIR}/include/QuadFormulas")
include_directories("${CMAKE_SOURCE_DIR}/include/Geometry")
include_directories("${CMAKE_SOURCE_DIR}/include/General")
include_directories("${CMAKE_SOURCE_DIR}/include/FE")
include_directories("${CMAKE_SOURCE_DIR}/include/System")
include_directories("${CMAKE_SOURCE_DIR}/include/PBE")
include_directories("${CMAKE_SOURCE_DIR}/Examples")


# message( "INCLUDE_DIRECTORIES = ${INCLUDE_DIRECTORIES}")


add_subdirectory(src/AMG) 

add_subdirectory(src/General)  
add_subdirectory(src/Geometry)          

# link_directories(${CMAKE_BINARY_DIR}/src/AMG/obj3D)
# link_directories(${CMAKE_BINARY_DIR}/src/General/obj3D)



add_executable(parmoon  ${MODEL})

target_link_libraries(parmoon PRIVATE "amg" ${BLAS_LIBRARIES})

# target_link_libraries(parmoon PRIVATE "gen" "geo" )


message("============================================================================================================"  )        
message("Configuration for the ParMooN Version 1.1")  
message("Written by Sashikumaar Ganesan, SERC, IISc Bangalore, India")  
message("============================================================================================================"  ) 


# FOREACH(var ${CMAKE_SYSTEM_PREFIX_PATH})
#     MESSAGE(STATUS "CMAKE_SYSTEM_PREFIX_PATH = " ${var})
# ENDFOREACH(var)
# MESSAGE(STATUS)

#         

# if (USE_MUMPS)
#   include_directories ("${PROJECT_SOURCE_DIR}/EXT_LIB/MUMPS/")
#   add_subdirectory (EXT_LIB/MUMPS/)
#   set (EXTRA_LIBS ${EXTRA_LIBS} EXT_LIB/MUMPS/)
# endif (USE_MUMPS)        
        


#  add_library(mylib src1.c src2.c)
#  add_executable(myexe main1.c)
#   set_property(
#    DIRECTORY
#    PROPERTY COMPILE_DEFINITIONS A AV=1
#    )
#  set_property(
#    TARGET mylib
#    PROPERTY COMPILE_DEFINITIONS B BV=2
#    )
#  set_property(
#    SOURCE src1.c
#    PROPERTY COMPILE_DEFINITIONS C CV=3
#    )
#  set_property(
#    SOURCE src1.c
#    APPEND PROPERTY COMPILE_DEFINITIONS D DV=4
#    )
# will build the source files with these definitions:
#  src1.c:   -DA -DAV=1 -DB -DBV=2 -DC -DCV=3
#  src2.c:   -DA -DAV=1 -DB -DBV=2
#  main2.c:  -DA -DAV=1


# find_package (BZip2)
# if (BZIP2_FOUND)
#   include_directories(${BZIP_INCLUDE_DIRS})
#   target_link_libraries (helloworld ${BZIP2_LIBRARIES})
# endif (BZIP2_FOUND)

