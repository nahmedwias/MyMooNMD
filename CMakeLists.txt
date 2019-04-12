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
#			Clemens Bartsch, WIAS Berlin
# DATE    : 05 June 2015
# ========================================================================

# normally no need to alter any thing in this file
# please DO NOT modify unless if it is necessary
# beginners are advised NOT to alter anything the below contents
# ========================================================================

cmake_minimum_required(VERSION 3.5)

project(ParMooN)
enable_language(Fortran)
include(FortranCInterface) 
#get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} gfortran)

# determine the parallel type (see https://stackoverflow.com/a/39671055)
list(APPEND PARMOON_PARALLEL_TYPE_LIST "MPI" "OMPONLY" "HYBRID" "SEQUENTIAL")
set(PARMOON_PARALLEL_TYPE "SEQUENTIAL" 
    CACHE STRING "The used parallelism. Choose between ${PARMOON_PARALLEL_TYPE_LIST}.")
# set the allowed strings as property. This way (graphical) editors for the 
# cache (CMakeCache.txt) can select only the allowed values
set_property(CACHE PARMOON_PARALLEL_TYPE PROPERTY STRINGS ${PARMOON_PARALLEL_TYPE_LIST})
list(FIND PARMOON_PARALLEL_TYPE_LIST ${PARMOON_PARALLEL_TYPE} _VALID_PAR)
if(_VALID_PAR EQUAL -1)
  message(FATAL_ERROR "Unknown PARMOON_PARALLEL_TYPE ${PARMOON_PARALLEL_TYPE}. "
                      "Choose between ${PARMOON_PARALLEL_TYPE_LIST}")
endif()
message(STATUS "Chosen parallel type is " ${PARMOON_PARALLEL_TYPE}.)

# Set two auxiliary variables, which make the file better readable imho
if("${PARMOON_PARALLEL_TYPE}" STREQUAL "MPI" OR "${PARMOON_PARALLEL_TYPE}" STREQUAL "HYBRID")
  set(PARMOON_USING_MPI TRUE)
endif()
if("${PARMOON_PARALLEL_TYPE}" STREQUAL "OMPONLY" OR "${PARMOON_PARALLEL_TYPE}" STREQUAL "HYBRID")
  set(PARMOON_USING_OMP TRUE)
endif()

# include all external libraries
include(cmake/external_libraries.cmake)

# include BOOST - needed for ROM/POD
find_package(Boost 1.47.0 REQUIRED COMPONENTS regex filesystem system program_options)
IF(Boost_FOUND)
    MESSAGE("Boost includes located in: ${Boost_INCLUDE_DIRS}")
    MESSAGE("Boost libraries located in: ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIRS})
    list(APPEND PARMOON_EXTERNAL_LIBRARIES ${Boost_LIBRARIES})
ENDIF(Boost_FOUND)

# turn on all warnings
add_compile_options(
  # clang/GCC warnings
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-pedantic>
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wno-attributes>
  # MSVC warnings
  $<$<CXX_COMPILER_ID:MSVC>:"/W4">)

add_compile_options($<$<BOOL:${PARMOON_USING_MPI}>:-D_MPI>
                    $<$<BOOL:${PARMOON_USING_OMP}>:-D_OMP>)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_compile_options(-D__MAC64__) # how can we get rid of this?
endif()


# If we want global compile flags to reach the subdirectories, we have to 
# run the following includes AFTER setting global flags.
include(cmake/sources_assembleroutines.cmake)
include(cmake/sources_assembleroutinessaddle.cmake)
include(cmake/sources_examples.cmake)
include(cmake/sources_fe.cmake)
include(cmake/sources_general.cmake)
include(cmake/sources_geometry.cmake)
include(cmake/sources_matrix.cmake)
include(cmake/sources_multigrid.cmake)
include(cmake/sources_parallel.cmake)
include(cmake/sources_quadformulas.cmake)
include(cmake/sources_refinement.cmake)
include(cmake/sources_rom.cmake)
include(cmake/sources_solver.cmake)
include(cmake/sources_system.cmake)
include(cmake/sources_timediscroutines.cmake)
# include(cmake/sources_pbe.cmake)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Define two static libraries. Each has its own precompiler flag (2D/3D).
add_library(parmoon_2d_${PARMOON_PARALLEL_TYPE} STATIC ${PARMOON_SOURCES_2D})
target_compile_definitions(parmoon_2d_${PARMOON_PARALLEL_TYPE} PUBLIC -D__2D__)
#target_compile_features(parmoon_2d_${PARMOON_PARALLEL_TYPE} PUBLIC cxx_std_11) # cmake 3.8

add_library(parmoon_3d_${PARMOON_PARALLEL_TYPE} STATIC ${PARMOON_SOURCES_3D} dummy.f)
target_compile_definitions(parmoon_3d_${PARMOON_PARALLEL_TYPE} PUBLIC -D__3D__)
#target_compile_features(parmoon_3d_${PARMOON_PARALLEL_TYPE} PUBLIC cxx_std_11) # cmake 3.8


##############################################################
# Start adding targets.

include(cmake/main_programs.cmake)
include(cmake/documentation_target.cmake)
include(cmake/get_repository_information.cmake)
###############################################################################
# parmoon defines 
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/all_defines_external_libraries.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/config/all_defines_external_libraries.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/config/)

###############################################################################
enable_testing()

###############################################################################
# the CMakeLists.txt for the user_projects may define more targets, tests, ...
add_subdirectory(user_projects)

###############################################################################
# tests 
add_subdirectory(test EXCLUDE_FROM_ALL)
####################################################
#                THE CHECK TARGET                  #
#   making this target results in building and     #
#   running all tests corresponding to the current #
#   parallel type. Use it excessively!             #
####################################################
add_custom_target(build_all_tests DEPENDS ${parmoon_tests})
add_custom_target(build_all_main DEPENDS ${parmoon_main_programs})
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --exclude-regex all
                  DEPENDS build_all_tests build_all_main)
add_custom_target(check_all COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure 
                  DEPENDS build_all_tests build_all_main)

include(cmake/additional_helper_targets.cmake)

############
# End adding executable targets.
############

message("============================================================================================================"  )        
message("Configuration for ParMooN, Version " ${PARMOON_VERSION})  
message("Author : Sashikumaar Ganesan")
message("         Numerical Mathematics and Scientific Computing")
message("         Supercomputer Edu. and Research Centre, Indian Institute of Science")
message("         Bangalore, India")  
message("Contact: mailt0:sashi@serc.iisc.in") 
message("============================================================================================================"  ) 




