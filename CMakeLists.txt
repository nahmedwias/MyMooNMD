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

if(USE_SYSTEM_BLAS)
  find_package(BLAS)
  include_directories(${BLAS_INCLUDE_DIRS})
  if(CMAKE_VERBOSE_MAKEFILE)
    if(BLAS_FOUND)
      message("Available system BLAS_LIB =${BLAS_LIBRARIES}")
    elseif(NOT BLAS_FOUND)
      message("BLAS_NOT_FOUND in your system, Searching in ParMooN")
    endif()
  endif()
endif(USE_SYSTEM_BLAS)

if(NOT BLAS_FOUND)
    message("BLAS_NOT_FOUND in your system, Searching in ParMooN")
    find_package(ACMLBLAS REQUIRED)
    include_directories(${ACMLBLAS_INCLUDE_DIRS})
    if(CMAKE_VERBOSE_MAKEFILE)
      if(ACMLBLAS_FOUND)
        message("Available system ACMLBLAS_LIB =${ACMLBLAS_LIBRARIES}")
      elseif(NOT ACMLBLAS_FOUND)
        message("ACMLBLAS_NOT_FOUND in ParMooN")
      endif(ACMLBLAS_FOUND)
    endif(CMAKE_VERBOSE_MAKEFILE)        
endif(NOT BLAS_FOUND)

message("==================Searching BLAS in the System - end=====================================================")


# search the UMFPACK in the system, if available it wil be default
message("==================Searching UMFPACK in the System - begin===================================================")
find_package(UMFPACK)
include_directories(${UMFPACK_INCLUDE_DIRS})
if(CMAKE_VERBOSE_MAKEFILE)
  if(UMFPACK_FOUND)
    message("Available system UMFPACK_LIB =${UMFPACK_LIBRARIES}")
  elseif(NOT UMFPACK_FOUND)
    message(FATAL_ERROR "UMFPACK_NOT_FOUND in your system or in the ParMooN")
  endif()
endif()
message("==================Searching UMFPACK in the System - end=====================================================")

if("${PARALLEL_TYPE}" STREQUAL "MPI" OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
# search the LAPACK in the system, if available it will be default
  message("==================Searching LAPACK in the System - begin===================================================")
 if(USE_SYSTEM_LAPACK)  
  find_package(LAPACK)
  include_directories(${LAPACK_INCLUDE_DIRS})
  if(CMAKE_VERBOSE_MAKEFILE)
    if(LAPACK_FOUND)
      message("Available system LAPACK_LIB =${LAPACK_LIBRARIES}")
    elseif(NOT LAPACK_FOUND)
        message(FATAL_ERROR "LAPACK_NOT_FOUND in your system")
    endif(LAPACK_FOUND)
  endif(CMAKE_VERBOSE_MAKEFILE)
 endif(USE_SYSTEM_LAPACK)
  
  if(NOT LAPACK_FOUND)
    message("LAPACK_NOT_FOUND in your system, Searching in ParMooN")
    find_package(LAPACKLIB REQUIRED)
#     include_directories(${ACMLBLAS_INCLUDE_DIRS})
    if(CMAKE_VERBOSE_MAKEFILE)
    if(LAPACK_FOUND)
      message("Available system ACMLBLAS_LIB =${ACMLBLAS_LIBRARIES}")
    elseif(NOT LAPACK_FOUND)
      message("LAPACK_NOT_FOUND in ParMooN")
    endif(LAPACK_FOUND)
   endif(CMAKE_VERBOSE_MAKEFILE)   
  endif(NOT LAPACK_FOUND)  
message("==================Searching LAPACK in the System - end=====================================================")  

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
endif("${PARALLEL_TYPE}" STREQUAL "MPI" OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")

# search the GRIDGEN in the system, if available it wil be default
message("==================Searching GRIDGEN in the System - begin==================================================")
find_package(GRIDGEN)
include_directories(${GRIDGEN_INCLUDE_DIRS})
if(CMAKE_VERBOSE_MAKEFILE)
  if(GRIDGEN_FOUND)
    message("Available system GRIDGEN_LIB =${GRIDGEN_LIBRARIES}")
  elseif(NOT GRIDGEN_FOUND)
    message(FATAL_ERROR "GRIDGEN_NOT_FOUND in your system or in the ParMooN")
  endif()
endif()
message("==================Searching GRIDGEN in the System - end====================================================")

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

# search the TECPLOT in the system, if available it wil be default
message("==================Searching TECPLOT in the System - begin==================================================")
find_package(TECPLOT)
include_directories(${TECPLOT_INCLUDE_DIRS})
if(CMAKE_VERBOSE_MAKEFILE)
  if(TECPLOT_FOUND)
    message("Available system TECPLOT_LIB =${TECPLOT_LIBRARIES}")
  elseif(NOT TECPLOT_FOUND)
    message(FATAL_ERROR "TECPLOT_NOT_FOUND in your system or in the ParMooN")
  endif()
endif()
message("==================Searching TECPLOT in the System - end====================================================")


# select the compiler type,
# CMAKE_BUILD_TYPE [DEBUG|RELEASE|RELWITHDEBINFO|MINSIZEREL]
set(CMAKE_BUILD_TYPE DEBUG CACHE STRING "options")
set(PARMOON_INCLUDE )

include (CMakeForceCompiler)
  set(CMAKE_SYSTEM_NAME ${ARCH})
  
if("${ARCH}" STREQUAL "MAC64")

  if("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
    find_path(MPI_CPATH   mpicc-mpich-mp)
    set(CMAKE_C_COMPILER ${MPI_CPATH}/mpicc-mpich-mp)
    find_path(MPI_CXXPATH   mpicxx-mpich-mp)  
    set(CMAKE_CXX_COMPILER ${MPI_CXXPATH}/mpicxx-mpich-mp)      
  elseif("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
    find_path(CPATH  clang)
    set(CMAKE_C_COMPILER ${CPATH}/clang)
    find_path(CXXPATH Clang)  
    set(CMAKE_CXX_COMPILER ${CXXPATH}/Clang)  
  endif("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
 
elseif("${ARCH}" STREQUAL "LINUX64")

  if("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
    find_path(MPI_CPATH   mpicc)
    set(CMAKE_C_COMPILER ${MPI_CPATH}/mpicc)
    find_path(MPI_CXXPATH   mpicxx)  
    set(CMAKE_CXX_COMPILER ${MPI_CXXPATH}/mpicxx)      
  elseif("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
    find_path(CPATH  gcc)
    set(CMAKE_C_COMPILER ${CPATH}/gcc)
    find_path(CXXPATH g++)  
    set(CMAKE_CXX_COMPILER ${CXXPATH}/g++)  
  endif("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
  
elseif("${ARCH}" STREQUAL "CRAY64")

    find_path(CPATH cc)
    set(CMAKE_C_COMPILER ${CPATH}/cc)
    find_path(CXXPATH CC)  
    set(CMAKE_CXX_COMPILER ${CXXPATH}/CC)  

endif()




  set(CMAKE_C_FLAGS "${PARMOON_CXX_DEF}")
  set(CMAKE_C_FLAGS_DEBUG "-g ${PARMOON_CXX_DEF}")
  set(CMAKE_C_FLAGS_RELEASE "O3 -s  ${PARMOON_CXX_DEF}")
  set(CMAKE_CXX_FLAGS "${PARMOON_CXX_DEF} ${PARMOON_CXX_DEF}" )
  set(CMAKE_CXX_FLAGS_DEBUG "${PARMOON_CXX_DEF} ${PARMOON_CXX_DEF}")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -s ${PARMOON_CXX_DEF} ${PARMOON_CXX_DEF}") 
 
include_directories("${CMAKE_SOURCE_DIR}/include/General")
include_directories("${CMAKE_SOURCE_DIR}/include/FE")
include_directories("${CMAKE_SOURCE_DIR}/include/Geometry")
include_directories("${CMAKE_SOURCE_DIR}/include/QuadFormulas")
include_directories("${CMAKE_SOURCE_DIR}/include/FE1D")
include_directories("${CMAKE_SOURCE_DIR}/include/FE2D")
include_directories("${CMAKE_SOURCE_DIR}/include/FE3D")
include_directories("${CMAKE_SOURCE_DIR}/include/AMG")
include_directories("${CMAKE_SOURCE_DIR}/include/PBE")
include_directories("${CMAKE_SOURCE_DIR}/include/Refinement")
include_directories("${CMAKE_SOURCE_DIR}/include/System")
include_directories("${CMAKE_SOURCE_DIR}/include/Parallel")
include_directories("${CMAKE_SOURCE_DIR}/Examples")


add_subdirectory(src/AMG) 
add_subdirectory(src/Refinement) 
add_subdirectory(src/QuadFormulas)  
add_subdirectory(src/Geometry)          
add_subdirectory(src/General) 
add_subdirectory(src/FE) 
add_subdirectory(src/System) 
# add_subdirectory(src/PBE) 

if("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
   add_subdirectory(src/Parallel) 
endif("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")

add_executable(parmoon_${PARALLEL_TYPE}.exe  ${MODEL})

if( "${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID" )
   target_link_libraries(parmoon_${PARALLEL_TYPE}.exe PRIVATE "par" "geo"  "fe" "gen"  "sys"
                      "amg"  "quad" "ref"  ${UMFPACK_LIBRARIES} ${TECPLOT_LIBRARIES}
                      ${BLAS_LIBRARIES} ${ACMLBLAS_LIBRARIES} ${TETGEN_LIBRARIES} ${GRIDGEN_LIBRARIES} ${MUMPS_LIBRARIES} "mpi" "mpifort" "mpicxx")  
elseif( "${PARALLEL_TYPE}" STREQUAL "SEQUENTIAL"  OR "${PARALLEL_TYPE}" STREQUAL "OMPONLY")
   target_link_libraries(parmoon_${PARALLEL_TYPE}.exe PRIVATE  "geo"  "fe" "gen"  "sys"
                      "amg"  "quad" "ref"  ${UMFPACK_LIBRARIES} ${TECPLOT_LIBRARIES}
                      ${BLAS_LIBRARIES} ${ACMLBLAS_LIBRARIES} ${TETGEN_LIBRARIES} ${GRIDGEN_LIBRARIES} )
endif("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")

message("============================================================================================================"  )  
message("GEO=${GEO}")
message("MODEL=${MODEL}")
message("ARCH=${ARCH}")
message("PARALLEL_TYPE=${PARALLEL_TYPE}")
message("C_COMPILER=${CMAKE_C_COMPILER}")
message("CXX_COMPILER=${CMAKE_CXX_COMPILER}")
message("PARMOON_CXX_DEF=${PARMOON_CXX_DEF}")
message("EXE_FILE=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/parmoon_${PARALLEL_TYPE}.exe")
message("============================================================================================================"  )        
message("Configuration for the ParMooN Version 1.1")  
message("Author : Sashikumaar Ganesan")
message("         Numerical Mathematics and Scientific Computing")
message("         Supercomputer Edu. and Research Centre")
message("         Indian Institute of Science")
message("         Bangalore, India")  
message("Question: www.serc.iisc.ernet.in/~sashi/") 
message("============================================================================================================"  ) 




