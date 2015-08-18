# ========================================================================
#     This file is part of ParMooN.
# 
#     ParMooN is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     ParMooN is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with ParMooN.  If not, see <http://www.gnu.org/licenses/>.
# =======================================================================
# PURPOSE :  configuration file for the ParMooN Version 1.1
# AUTHOR  : Sashikumaar Ganesan, SERC, IISc Bangalore, India
# DATE    : 05 June 2015
# ========================================================================

# normally no need to alter any thing in this file
# please DO NOT modify unless if it is necessary
# beginners are advised NOT to alter anything the below contents
#
# Change History 
# 2015/08/18 Clemens Bartsch 

cmake_minimum_required(VERSION 2.8)

project(ParMooN)
set(ParMooN_VERSION_MAJOR 1)
set(ParMooN_VERSION_MINOR 1)

include(${CMAKE_SOURCE_DIR}/UserConfig.cmake)

set(PARMOON_EXTLIB_PATH ${PROJECT_SOURCE_DIR}/EXT_LIB)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

# Search and include BLAS library.
message("==================Searching BLAS - begin===================================================")

if(USE_SYSTEM_BLAS)
message("Using default CMake module to find BLAS.")
  find_package(BLAS MODULE)
  if(BLAS_FOUND)
    message("BLAS found with default CMake module.")
    if(CMAKE_VERBOSE_MAKEFILE)
  	  message("Available BLAS_LIB =${BLAS_LIBRARIES}")
  	endif(CMAKE_VERBOSE_MAKEFILE)
  	include_directories(${BLAS_INCLUDE_DIRS})
  endif()
endif(USE_SYSTEM_BLAS)

if(NOT BLAS_FOUND)
    message("Using ParMooN CMake module to find BLAS.")
    find_package(ACMLBLAS REQUIRED MODULE)
    if(BLAS_FOUND)
      message("BLAS found with ParMooN CMake module.")
      if(CMAKE_VERBOSE_MAKEFILE)
  	    message("Available BLAS_LIB =${BLAS_LIBRARIES}")
  	  endif(CMAKE_VERBOSE_MAKEFILE)
      include_directories(${BLAS_INCLUDE_DIRS})
    else()
      # This case should never be reached, because of the REQUIRED option.
      message("BLAS not found. Terminating.")
    endif() 
endif(NOT BLAS_FOUND)

message("==================Searching BLAS - end=====================================================")


# Search and include UMFPACK library.
message("==================Searching UMFPACK - begin================================================")
message("Using ParMooN CMake module to find UMFPACK.")
find_package(UMFPACK MODULE REQUIRED)
if(UMFPACK_FOUND)
  message("UMFPACK found with ParMooN CMake module.")
  include_directories(${UMFPACK_INCLUDE_DIRS})
  if(CMAKE_VERBOSE_MAKEFILE)
    message("Available UMFPACK_LIB =${UMFPACK_LIBRARIES}")
  endif(CMAKE_VERBOSE_MAKEFILE)
else()
  # This case should never be reached, because of the REQUIRED option.
  message("UMFPACK not found. Terminating.")
endif() 
message("==================Searching UMFPACK - end==================================================")

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
#     include_directories(${BLAS_INCLUDE_DIRS})
    if(CMAKE_VERBOSE_MAKEFILE)
    if(LAPACK_FOUND)
      message("Available system LAPACK_LIB =${LAPACK_LIBRARIES}")
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
message("cmake system name " ${CMAKE_SYSTEM_NAME})
set(CMAKE_SYSTEM_NAME ${ARCH})
message("cmake system name " ${CMAKE_SYSTEM_NAME})

if("${ARCH}" STREQUAL "MAC64")

  if("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
    find_path(MPI_CPATH   mpicc-mpich-mp)
    set(CMAKE_C_COMPILER ${MPI_CPATH}/mpicc-mpich-mp)
    find_path(MPI_CXXPATH   mpicxx-mpich-mp)  
    set(CMAKE_CXX_COMPILER ${MPI_CXXPATH}/mpicxx-mpich-mp)      
  elseif("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
    find_path(CPATH  clang)
    set(CMAKE_C_COMPILER ${CPATH}/clang  )
    find_path(CXXPATH Clang)  
    set(CMAKE_CXX_COMPILER ${CXXPATH}/Clang)  
  endif("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
 
elseif("${ARCH}" STREQUAL "LINUX64")

  if("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
    find_path(MPI_CPATH   mpicc)
    set(CMAKE_C_COMPILER ${MPI_CPATH}/mpicc)
    find_path(MPI_CXXPATH   mpicxx)  
    set(CMAKE_CXX_COMPILER ${MPI_CXXPATH}/mpicxx)      
  elseif("${PARALLEL_TYPE}" STREQUAL "SEQUENTIAL"  OR "${PARALLEL_TYPE}" STREQUAL "OMPONLY")
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


# set the define for compiler  options
  set(CMAKE_C_FLAGS "${PARMOON_C_DEF}")
  set(CMAKE_C_FLAGS_DEBUG "-g ${PARMOON_C_DEF}")
  set(CMAKE_C_FLAGS_RELEASE "O3 -s  ${PARMOON_C_DEF}")
  set(CMAKE_CXX_FLAGS "${PARMOON_CXX_DEF}" )
  set(CMAKE_CXX_FLAGS_DEBUG "${PARMOON_CXX_DEF} ")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -s ${PARMOON_CXX_DEF}") 
 
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
add_subdirectory(Examples)
# add_subdirectory(src/PBE) 

if("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
   add_subdirectory(src/Parallel) 
endif("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")

add_executable(parmoon_${GEO}_${PARALLEL_TYPE}.exe  ${MODEL})

if( "${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID" )
   target_link_libraries(parmoon_${GEO}_${PARALLEL_TYPE}.exe "par" "geo"  "fe" "gen"  "sys" "examples"
                      "amg"  "quad" "ref"  ${UMFPACK_LIBRARIES} ${TECPLOT_LIBRARIES}
                      ${BLAS_LIBRARIES} ${TETGEN_LIBRARIES} ${GRIDGEN_LIBRARIES} ${MUMPS_LIBRARIES} "mpi" "mpifort" "mpicxx")  
elseif( "${PARALLEL_TYPE}" STREQUAL "SEQUENTIAL"  OR "${PARALLEL_TYPE}" STREQUAL "OMPONLY")
   target_link_libraries(parmoon_${GEO}_${PARALLEL_TYPE}.exe "geo"  "fe" "gen" "geo"  "sys" "examples"
                      "amg"  "quad" "ref" 
                      "geo"  "fe" "gen" "geo"  "sys" "examples"
                      "amg"  "quad" "ref"  ${UMFPACK_LIBRARIES} ${TECPLOT_LIBRARIES}
                      ${BLAS_LIBRARIES} ${TETGEN_LIBRARIES} ${GRIDGEN_LIBRARIES} )
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
message("         Supercomputer Edu. and Research Centre, Indian Institute of Science")
message("         Bangalore, India")  
message("Contact: mailt0:sashi@serc.iisc.in") 
message("============================================================================================================"  ) 




