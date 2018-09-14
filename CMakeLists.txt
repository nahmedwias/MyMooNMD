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

cmake_minimum_required(VERSION 2.8.12)


project(ParMooN)
set(ParMooN_VERSION_MAJOR 1)
set(ParMooN_VERSION_MINOR 1)

# Enable Fortran compiling, linking and mangling.
#enable_language(Fortran)
#include(FortranCInterface) 
#get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} /Users/blank/ParMooN/ParMooN/gfortran64)

# Determine ParMoon-internal architecture variable "_ARCH" from the used system.
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    set(_ARCH "LINUX64")
elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    set(_ARCH "MAC64")
elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
    set(_ARCH "INTEL64") #TODO Is this what should be used on windows?
else()
    #TODO Add the right handling for Tyrone and Clay cluster here.
    message(FATAL_ERROR "Unknown CMAKE_SYSTEM_NAME: " ${CMAKE_SYSTEM_NAME} 
        ". Rework CMakeLists.txt to fit your system.")
endif()
message(STATUS "Used system name is " ${CMAKE_SYSTEM_NAME} ", " 
    "ParMooN internal architecture name is " ${_ARCH} "."  )

#Set PARMOON_... cache variables which are designated for change by the user,
#and process them immediately.
set(PARMOON_OUTPUT_DIRECTORY CACHE STRING "The directory where ParMooN executables will be put.")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PARMOON_OUTPUT_DIRECTORY})

set(PARMOON_PARALLEL_TYPE "SEQUENTIAL" CACHE STRING "The used parallelism. Choose between MPI, OMPONLY, HYBRID and SEQUENTIAL.")
list(APPEND _PAR_LIST "MPI" "OMPONLY" "HYBRID" "SEQUENTIAL")
list(FIND _PAR_LIST ${PARMOON_PARALLEL_TYPE} _VALID_PAR)
if(_VALID_PAR EQUAL -1)
    message(FATAL_ERROR "Unknown PARMOON_PARALLEL_TYPE '" ${PARMOON_PARALLEL_TYPE} "'. Choose between MPI, OMPONLY, HYBRID and SEQUENTIAL.")
endif()
message(STATUS "Chosen parallel type is " ${PARMOON_PARALLEL_TYPE}.)

# Set two auxiliary variables, which make the file better readable imho
if("${PARMOON_PARALLEL_TYPE}" STREQUAL "MPI" OR "${PARMOON_PARALLEL_TYPE}" STREQUAL "HYBRID")
  set(_USING_MPI TRUE)
endif()
if("${PARMOON_PARALLEL_TYPE}" STREQUAL "OMPONLY" OR "${PARMOON_PARALLEL_TYPE}" STREQUAL "HYBRID")
  set(_USING_OMP TRUE)
endif()

# include all external libraries
include(cmake/external_libraries.cmake)

# turn on all warnings
IF(CMAKE_BUILD_TOOL MATCHES "make")
  IF(NOT CMAKE_CXX_FLAGS MATCHES "-Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-attributes")
  ENDIF(NOT CMAKE_CXX_FLAGS MATCHES "-Wall")
  IF(NOT CMAKE_C_FLAGS MATCHES "-Wall")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-attributes")
  ENDIF(NOT CMAKE_C_FLAGS MATCHES "-Wall")
ENDIF(CMAKE_BUILD_TOOL MATCHES "make")



# Determine precompiler definitions ("PCD") depending on the chosen _PARALLEL_TYPE
if("${PARMOON_PARALLEL_TYPE}" STREQUAL "MPI")
    set(_PCD_PARALLEL_TYPE -D_PAR -D_MPIONLY -D_MPI)
elseif("${PARMOON_PARALLEL_TYPE}" STREQUAL "OMPONLY")
    set(_PCD_PARALLEL_TYPE -D_PAR -D_OMPONLY -D_OMP)
elseif("${PARMOON_PARALLEL_TYPE}" STREQUAL "HYBRID")
    set(_PCD_PARALLEL_TYPE -D_PAR -D_MPI -D_HYBRID)
elseif("${PARMOON_PARALLEL_TYPE}" STREQUAL "SEQUENTIAL")
    set(_PCD_PARALLEL_TYPE -D_SEQ)
endif()
add_definitions(${_PCD_PARALLEL_TYPE})

# Enable the right C++ 11 compiler flag.
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" _COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" _COMPILER_SUPPORTS_CXX0X)
if(_COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(_COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(WARNING "The compiler ${CMAKE_CXX_COMPILER} does neither support -std=c++11 nor "
    "-std=c++0x compiler flag. Are you sure C++11 is supported?")
endif()

# Compiler flags depending on system architecture.
if("${_ARCH}" STREQUAL "LINUX64")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMKL_ILP64 ")
elseif("${_ARCH}" STREQUAL "MAC64")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMKL_ILP64 -m64 -fapple-pragma-pack ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -fapple-pragma-pack -Wdeprecated-register ")
elseif("${_ARCH}" STREQUAL "INTEL64")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMKL_ILP64 -m64  -fapple-pragma-pack") 
elseif("${_ARCH}" STREQUAL "TYRONE64")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMPICH_IGNORE_CXX_SEEK")
elseif("${_ARCH}" STREQUAL "CRAY64")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_IGNORE_CXX_SEEK")
endif()

# TODO The PRIVATE flag is used at very few places only - remove it from global
# and e.g. attach to those targets, which make use of it (TNSE2D/SystemMatTNSE2D/SystemTNSE2D)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__${_ARCH}__   -D__PRIVATE__") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__${_ARCH}__   -D__PRIVATE__")


# If we want global compile flags to reach the subdirectories, we have to 
# run the following includes AFTER setting global flags.
include(cmake/sources_assembleroutines.cmake)
include(cmake/sources_assembleroutinessaddle.cmake)
include(cmake/sources_examples.cmake)
include(cmake/sources_fe.cmake)
include(cmake/sources_general.cmake)
include(cmake/sources_geometry.cmake)
include(cmake/sources_multigrid.cmake)
include(cmake/sources_parallel.cmake)
include(cmake/sources_quadformulas.cmake)
include(cmake/sources_refinement.cmake)
include(cmake/sources_solver.cmake)
include(cmake/sources_system.cmake)
include(cmake/sources_matrix.cmake)
include(cmake/sources_timediscroutines.cmake)
# include(cmake/sources_pbe.cmake)

# Define two static libraries. Each has its own precompiler flag (2D/3D).
# Adding dummy.f to the libraries is a hack - no automatic fortran linking (lapack!) without fortran compiling.
add_library(parmoon_2d_${PARMOON_PARALLEL_TYPE} STATIC ${PARMOON_SOURCES_2D} dummy.f)
target_compile_definitions(parmoon_2d_${PARMOON_PARALLEL_TYPE} PUBLIC -D__2D__)

add_library(parmoon_3d_${PARMOON_PARALLEL_TYPE} STATIC ${PARMOON_SOURCES_3D} dummy.f)
target_compile_definitions(parmoon_3d_${PARMOON_PARALLEL_TYPE} PUBLIC -D__3D__)


##############################################################
# Start adding executable targets.
##############################################################

### Standard 2D Programs ###

### CD2D ###
# Standard executable target.
add_executable(cd2d ${PROJECT_SOURCE_DIR}/2DPrograms/CD2D_ParMooN.C)
set_target_properties(cd2d PROPERTIES OUTPUT_NAME parMooN_cd2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(cd2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(cd2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### Darcy2D ###
# Standard executable target.
add_executable(darcy2d ${PROJECT_SOURCE_DIR}/2DPrograms/Darcy2D_ParMooN.C)
set_target_properties(darcy2d PROPERTIES OUTPUT_NAME parMooN_darcy2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(darcy2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(darcy2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### NSE2D ###
# Standard executable target.
add_executable(nse2d ${PROJECT_SOURCE_DIR}/2DPrograms/NSE2D_ParMooN.C)
set_target_properties(nse2d PROPERTIES OUTPUT_NAME parMooN_nse2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(nse2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(nse2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})


### Brinkman2D ###
# Standard executable target.
add_executable(brinkman2d ${PROJECT_SOURCE_DIR}/2DPrograms/Brinkman2D_ParMooN.C)
set_target_properties(brinkman2d PROPERTIES OUTPUT_NAME parMooN_brinkman2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(brinkman2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(brinkman2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### Brinkman3D ###
# Standard executable target.
add_executable(brinkman3d ${PROJECT_SOURCE_DIR}/3DPrograms/Brinkman3D_ParMooN.C)
set_target_properties(brinkman3d PROPERTIES OUTPUT_NAME parMooN_brinkman3d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(brinkman3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(brinkman3d parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES}) 

### TCD2D ###
# Standard executable target.
add_executable(tcd2d ${PROJECT_SOURCE_DIR}/2DPrograms/TCD2D_ParMooN.C)
set_target_properties(tcd2d PROPERTIES OUTPUT_NAME parMooN_tcd2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tcd2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(tcd2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### TNSE2D ###
# Standard executable target.
add_executable(tnse2d ${PROJECT_SOURCE_DIR}/2DPrograms/TNSE2D_ParMooN.C)
set_target_properties(tnse2d PROPERTIES OUTPUT_NAME parMooN_tnse2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tnse2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(tnse2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### Thermohydraulic Brinkman2D ###
# Standard executable target.
add_executable(thermohydraulic_brinkman2d ${PROJECT_SOURCE_DIR}/2DPrograms/Thermohydraulic_Brinkman2D_ParMooN.C)
set_target_properties(thermohydraulic_brinkman2d PROPERTIES OUTPUT_NAME parMooN_thermohydraulic_brinkman2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(thermohydraulic_brinkman2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(thermohydraulic_brinkman2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### TThermohydraulic Brinkman2D ###
# Standard executable target.
add_executable(tthermohydraulic_brinkman2d ${PROJECT_SOURCE_DIR}/2DPrograms/TThermohydraulic_Brinkman2D_ParMooN.C)
set_target_properties(tthermohydraulic_brinkman2d PROPERTIES OUTPUT_NAME parMooN_tthermohydraulic_brinkman2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tthermohydraulic_brinkman2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(tthermohydraulic_brinkman2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### MESH ###
# Standard executable target.
add_executable(mesh ${PROJECT_SOURCE_DIR}/2DPrograms/mesh_ParMooN.C)
set_target_properties(mesh PROPERTIES OUTPUT_NAME parMooN_mesh_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(mesh PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(mesh parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### MESH2GEO(2D) ###
# Standard executable target.
add_executable(geo2mesh2d ${PROJECT_SOURCE_DIR}/2DPrograms/geo2mesh2d_ParMooN.C)
set_target_properties(geo2mesh2d PROPERTIES OUTPUT_NAME parMooN_geo2mesh2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(geo2mesh2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(geo2mesh2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### Standard 3D Programs ###

### CD3D ###
# Standard executable target.
add_executable(cd3d ${PROJECT_SOURCE_DIR}/3DPrograms/CD3D_ParMooN.C)
set_target_properties(cd3d PROPERTIES OUTPUT_NAME parMooN_cd3d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(cd3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(cd3d parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### Darcy3D ###
add_executable(darcy3d ${PROJECT_SOURCE_DIR}/3DPrograms/Darcy3D_ParMooN.C)
set_target_properties(darcy3d PROPERTIES OUTPUT_NAME parMooN_darcy3d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(darcy3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(darcy3d parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### TCD3D ###
add_executable(tcd3d ${PROJECT_SOURCE_DIR}/3DPrograms/TCD3D_ParMooN.C)
set_target_properties(tcd3d PROPERTIES OUTPUT_NAME parMooN_tcd3d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tcd3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(tcd3d parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### NSE3D ###
# Standard executable target.
add_executable(nse3d ${PROJECT_SOURCE_DIR}/3DPrograms/NSE3D_ParMooN.C)
set_target_properties(nse3d PROPERTIES OUTPUT_NAME parMooN_nse3d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(nse3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(nse3d parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### TNSE3D ###
# Standard executable target.
add_executable(tnse3d ${PROJECT_SOURCE_DIR}/3DPrograms/TNSE3D_ParMooN.C)
set_target_properties(tnse3d PROPERTIES OUTPUT_NAME parMooN_tnse3d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tnse3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(tnse3d parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

### MESH(3D) ###
# Standard executable target.
add_executable(mesh3d ${PROJECT_SOURCE_DIR}/3DPrograms/mesh3d_ParMooN.C)
set_target_properties(mesh3d PROPERTIES OUTPUT_NAME parMooN_mesh3d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(mesh3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(mesh3d parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

###############################################################################
# building the documentation using doxygen (define target 'doc')
find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(PARMOON_DOCUMENTATION_DIRECTORY "" CACHE PATH 
      "The directory for the Doxygen documentation.")
  if(NOT PARMOON_DOCUMENTATION_DIRECTORY)
    add_custom_target(doc ${CMAKE_COMMAND} -E cmake_echo_color --red 
                      "Please set PARMOON_DOCUMENTATION_DIRECTORY in your "
                      "CMakeCache.txt file. No documentation is created!")
  else()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/documentation/Doxyfile.in
                   ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
                      ${DOXYGEN_EXECUTABLE} 
                      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating documentation with Doxygen" VERBATIM)
  endif(NOT PARMOON_DOCUMENTATION_DIRECTORY)
else()
  add_custom_target(doc ${CMAKE_COMMAND} -E cmake_echo_color --red 
                    "It seems cmake could not find a doxygen executable. "
                    "Therefore, no documentation is created!")
endif(DOXYGEN_FOUND)

include(cmake/get_repository_information.cmake)

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
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure 
                  DEPENDS ${parmoon_tests})


############
# End adding executable targets.
############

message("============================================================================================================"  )        
message("CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message("Configuration for the ParMooN Version 1.1")  
message("Author : Sashikumaar Ganesan")
message("         Numerical Mathematics and Scientific Computing")
message("         Supercomputer Edu. and Research Centre, Indian Institute of Science")
message("         Bangalore, India")  
message("Contact: mailt0:sashi@serc.iisc.in") 
message("============================================================================================================"  ) 




