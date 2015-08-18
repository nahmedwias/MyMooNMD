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
# 2015/08/18 Clemens Bartsch: Introduced distinction between 
# FIND_USER_{PACKAGENAME} and USE_SYSTEM_{PACKAGENAME}. Tidied up output
# a bit. First rework of search modules in directory cmake.
# Moved several variables between CMakeLists.txt and UserConfig.cmake.
# ========================================================================

cmake_minimum_required(VERSION 2.8)

project(ParMooN)
set(ParMooN_VERSION_MAJOR 1)
set(ParMooN_VERSION_MINOR 1)

# Include the user specific settings.
include(${CMAKE_SOURCE_DIR}/UserConfig.cmake)

# Set path to the external libraries which ship with ParMooN
set(PARMOON_EXTLIB_PATH ${PROJECT_SOURCE_DIR}/EXT_LIB)

# Set path to ParMooN package search modules.
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

# Search and include BLAS library.
message("==================Searching BLAS - begin===================================================")

	if(USE_SYSTEM_BLAS)
	  message("Using default CMake module to find BLAS.")
	    find_package(BLAS MODULE)
	    if(BLAS_FOUND)
	      message("BLAS found with default CMake module.")
		  message(STATUS "Available BLAS_LIB =${BLAS_LIBRARIES}")
	  	  include_directories(${BLAS_INCLUDE_DIRS})
	    endif(BLAS_FOUND)
	endif(USE_SYSTEM_BLAS)
	
	if(NOT BLAS_FOUND)
	  message("Using ParMooN CMake module to find BLAS.")
	  find_package(ACMLBLAS REQUIRED MODULE)
	  # Proceeding form here only when found - due to REQUIRED keyword.
	  message("BLAS found with ParMooN CMake module.")
	  message(STATUS "Available BLAS_LIB =${BLAS_LIBRARIES}")
	  include_directories(${BLAS_INCLUDE_DIRS})
	endif(NOT BLAS_FOUND)

message("==================Searching BLAS - end=====================================================")


# Search and include UMFPACK library.
message("==================Searching UMFPACK - begin================================================")
	message("Using ParMooN CMake module to find UMFPACK.")
	find_package(UMFPACK MODULE REQUIRED)
	message("UMFPACK found with ParMooN CMake module.")
	include_directories(${UMFPACK_INCLUDE_DIRS})
	message(STATUS "Available UMFPACK_LIB =${UMFPACK_LIBRARIES}")
message("==================Searching UMFPACK - end==================================================")

# The parallel solvers MUMPS and LAPACK will only be searched for, when PARALLEL_TYPE is MPI or HYBRID.
if("${PARALLEL_TYPE}" STREQUAL "MPI" OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
  # Search and include LAPACK library.
  message("==================Searching LAPACK - begin=================================================")
  
	  if(USE_SYSTEM_LAPACK)
	    message("Using default CMake module to find LAPACK.")
	    find_package(LAPACK MODULE)
	    if(LAPACK_FOUND)
	      message("LAPACK found with default CMake module.")
		  message(STATUS "Available LAPACK_LIB =${LAPACK_LIBRARIES}")
	  	  include_directories(${LAPACK_INCLUDE_DIRS})
	    endif(LAPACK_FOUND)  
	  endif(USE_SYSTEM_LAPACK)
	  
	  if(NOT LAPACK_FOUND)
	    message("Using ParMooN CMake module to find LAPACK.")
	    find_package(LAPACKLIB REQUIRED MODULE)
	    message("LAPACK found with ParMooN CMake module.")
	    message(STATUS "Available LAPACK_LIB =${LAPACK_LIBRARIES}")
	    include_directories(${LAPACK_INCLUDE_DIRS})
	  endif(NOT LAPACK_FOUND)
  
  message("==================Searching LAPACK - end===================================================")  

  # Search and include MUMPS library.
  message("==================Searching MUMPS - begin==================================================")
	  message("Using ParMooN CMake module to find MUMPS.")
	  find_package(MUMPS REQUIRED MODULE)
	  message("MUMPS found with ParMooN CMake module.")
	  include_directories(${MUMPS_INCLUDE_DIRS})
	  message(STATUS "Available MUMPS_LIB =${MUMPS_LIBRARIES}")
  message("==================Searching MUMPS - end====================================================")
  
endif("${PARALLEL_TYPE}" STREQUAL "MPI" OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")

# Search and include GRIDGEN library.
message("==================Searching GRIDGEN - begin================================================")
	message("Using ParMooN CMake module to find GRIDGEN.")
	find_package(GRIDGEN REQUIRED MODULE)
	message("GRIDGEN found with ParMooN CMake module.")
	include_directories(${GRIDGEN_INCLUDE_DIRS})
	message(STATUS "Available system GRIDGEN_LIB =${GRIDGEN_LIBRARIES}")
message("==================Searching GRIDGEN - end==================================================")

# Search and include TETGEN library, only when GEO=3D.
if("${GEO}" STREQUAL "3D")
  message("==================Searching TETGEN - begin================================================")
	  message("Using ParMooN CMake module to find TETGEN.")
	  find_package(TETGEN REQUIRED MODULE)
	  message("TETGEN found with ParMooN CMake module.")
	  include_directories(${TETGEN_INCLUDE_DIRS}) 
	  message(STATUS "Available system TETGEN_LIB =${TETGEN_LIBRARIES}")
  message("==================Searching TETGEN - end==================================================")
endif("${GEO}" STREQUAL "3D")

# TODO CB 2015/08/18 Without setting CONFIGURED_ONCE an internal cache flag
# or fill the variable warnings with compiler flags, this code is useless.
#if (NOT CONFIGURED_ONCE)
#    set(CMAKE_CXX_FLAGS "${warnings}"
#        CACHE STRING "Flags used by the compiler during all build types." FORCE)
#    set(CMAKE_C_FLAGS   "${warnings}"
#        CACHE STRING "Flags used by the compiler during all build types." FORCE)
#endif()

# Search and include TECPLOT library.
message("==================Searching TECPLOT - begin================================================")
	message("Using ParMooN CMake module to find TECPLOT.")
	find_package(TECPLOT REQUIRED MODULE)
	message("TECPLOT found with ParMooN CMake module.")
	include_directories(${TECPLOT_INCLUDE_DIRS})
	message(STATUS "Available system TECPLOT_LIB =${TECPLOT_LIBRARIES}")
message("==================Searching TECPLOT - end==================================================")

# select the compiler type,
# CMAKE_BUILD_TYPE [DEBUG|RELEASE|RELWITHDEBINFO|MINSIZEREL]
set(CMAKE_BUILD_TYPE DEBUG CACHE STRING "options")
set(PARMOON_INCLUDE )

include (CMakeForceCompiler)
message("What is CMAKE_SYSTEM? " ${CMAKE_SYSTEM})
message("What is CMAKE_SYSTEM_NAME? " ${CMAKE_SYSTEM_NAME})
#TODO CB I don't think it's a good idea to change the system name (which we never use anyway) to fit our custom ARCH.
#set(CMAKE_SYSTEM_NAME ${ARCH})
message("What is CMAKE_SYSTEM? " ${CMAKE_SYSTEM})
message("What is CMAKE_SYSTEM_NAME? " ${CMAKE_SYSTEM_NAME})


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
    message("What is MPI_CPATH?" ${MPI_CPATCH})
    set(CMAKE_C_COMPILER ${MPI_CPATH}/mpicc)
    find_path(MPI_CXXPATH   mpicxx)
    message("What is MPI_CXXPATH?" ${MPI_CXXPATCH}) 
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




