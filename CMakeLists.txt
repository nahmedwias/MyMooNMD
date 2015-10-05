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
#
# 2015/08/18 Clemens Bartsch: Introduced distinction between 
# FIND_USER_{PACKAGENAME} and DEFAULT_MODULE_{PACKAGENAME}. Tidied up output
# a bit. First rework of search modules in directory cmake.
# Moved several variables between CMakeLists.txt and UserConfig.cmake.
#
# 2015/08/25 Clemens Bartsch: Major rework, of find modules as well as the
# (portable) adding of compiler and linker options. Also changed the executable
# target system. Now all default targets exist at the same time, but only
# for the currently chosen PARALLEL_TYPE. To get other PARALLEL_TYPE
# targets, change UserConfig.cmake accordingly and rerun cmake.
#
# 2015/08/25 Clemens Bartsch:
#    - renamed gridgen to triangle in the whole project
#    - temporarily disabled default find of lapack and blas library, which continuously made trouble
#    - moved -DTRILIBRARY and -DREDUCED flag to triangle.c
# 2015/08/28 Clemens Bartsch:
#    - stopped writing a history into the files, for that's what hg is for
# ========================================================================

cmake_minimum_required(VERSION 2.8.12)


project(ParMooN)
set(ParMooN_VERSION_MAJOR 1)
set(ParMooN_VERSION_MINOR 1)

# Enable Fortran compiling, linking and mangling.
enable_language(Fortran)
include(FortranCInterface) 

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

# Set path to the external libraries which ship with ParMooN
set(PARMOON_EXTLIB_PATH ${CMAKE_SOURCE_DIR}/EXT_LIB)
# Set path to ParMooN package search modules.
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

# Search and include LAPACK library, which also includes searching for BLAS. 
message("==================Searching LAPACK - begin=================================================")
    find_package(LAPACK REQUIRED) #Complains and exits when LAPACK or BLAS are not found.
    # Process what has been found - add linker flags.
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")  
message("==================Searching LAPACK - end===================================================")  

# When using OpenMP, find it and use the module's output.
if(_USING_OMP)
  message("==================Searching OPENMP - begin=================================================")
      # CB Kick off the default FindOpenMP module.
      find_package(OpenMP REQUIRED)
      if(OPENMP_FOUND)
        # If OpenMP is used, add the found compiler flags for c and c++ globally.
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      endif(OPENMP_FOUND)
  message("==================Searching OPENMP - end===================================================")
endif(_USING_OMP)

# When using MPI, find mpi and some other libraries which use mpi.
if(_USING_MPI) 
  #Kickoff the default MPI search module and use its output.
  message("==================Searching MPI - begin====================================================")
  find_package(MPI REQUIRED)
  #set path to the mpi header files.
  include_directories(${MPI_INCLUDE_PATH})
  #set compile and link flags
  if(MPI_COMPILE_FLAGS)
      message(STATUS "There are mpi compile flags!")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
      message(STATUS "There are mpi linker flags!")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
  endif() 
  #Combine MPI libraries for all 3 languages to one list used for linking.
  list(APPEND _MPI_LIBRARIES ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
  message("==================Searching MPI - end======================================================")
  
  #Find a thread library.
  message("==================Searching Thread Library - begin=========================================")
    set(CMAKE_THREAD_PREFER_PTHREAD ON)
    find_package(Threads REQUIRED)
  message("==================Searching Thread Library - end===========================================")
  
  #When using MPI we're also using Parmetis. Include it.
  message("==================Including Parmetis - begin===============================================")
      find_path(PARMETIS_INCLUDE_DIR  parmetis.h PATHS ${PARMOON_EXTLIB_PATH}/Metis/Lib NO_DEFAULT_PATH)
      find_library(PARMETIS_LIBRARY NAMES parmetis_${_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Metis NO_DEFAULT_PATH)
      include_directories(BEFORE ${PARMETIS_INCLUDE_DIRS})
      
      set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR})
      set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY})
      
      if(NOT (PARMETIS_LIBRARIES AND  PARMETIS_INCLUDE_DIRS)) #Mimic "REQUIRED" behaviour
        message(FATAL_ERROR "Parmetis is missing from ParMooN.")
    else()
        message(STATUS "Done.")
    endif()
  message("==================Including Parmetis - end=================================================")
   
  #When using MPI we're also using MUMPS solver. Include it.
  message("==================Including MUMPS - begin==================================================")
    
    find_path(MUMPS_INCLUDE_DIR  mumps_compat.h PATHS ${PARMOON_EXTLIB_PATH}/MUMPS/include NO_DEFAULT_PATH)
    set(_MUMPS_LIBS_PATH ${PARMOON_EXTLIB_PATH}/MUMPS/lib)
    # Find mumps and its dependent libraries.
    find_library(MUMPS_LIBRARY NAMES dmumps_${_ARCH} PATHS ${_MUMPS_LIBS_PATH} NO_DEFAULT_PATH)
    # mumps_common library
    find_library(MUMPS_LIBRARY_COMMON NAMES mumps_common_${_ARCH} PATHS ${_MUMPS_LIBS_PATH} NO_DEFAULT_PATH)
    # esmumps/dmumps/tesmumps
    find_library(MUMPS_ESMUMPS_LIBRARY NAMES esmumps_TYRONE64 PATHS ${_MUMPS_LIBS_PATH} NO_DEFAULT_PATH)
    find_library(MUMPS_DMUMPS_LIBRARY NAMES dmumps_${_ARCH} PATHS ${_MUMPS_LIBS_PATH} NO_DEFAULT_PATH)
    find_library(MUMPS_TESMUMPS_LIBRARY NAMES ptesmumps_${_ARCH} PATHS ${_MUMPS_LIBS_PATH} NO_DEFAULT_PATH)
    # scotch and pord libraries
    find_library(MUMPS_PORD_LIBRARY NAMES pord_${_ARCH} PATHS ${_MUMPS_LIBS_PATH} NO_DEFAULT_PATH)
    find_library(MUMPS_SCOTCH_LIBRARY NAMES scotch_${_ARCH} PATHS ${_MUMPS_LIBS_PATH} NO_DEFAULT_PATH)
    find_library(MUMPS_SCOTCHERR_LIBRARY NAMES scotcherr_TYRONE64 PATHS ${_MUMPS_LIBS_PATH} NO_DEFAULT_PATH)
    #scalapack library 
    find_library(MUMPS_SCALAPACK_LIBRARY NAMES scalapack_${_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS NO_DEFAULT_PATH) 
     
    set(MUMPS_LIBRARIES ${MUMPS_LIBRARY} ${MUMPS_DMUMPS_LIBRARY} ${MUMPS_ESMUMPS_LIBRARY} ${MUMPS_TESMUMPS_LIBRARY}
                        ${MUMPS_LIBRARY_COMMON} ${MUMPS_SCOTCH_LIBRARY} ${MUMPS_SCOTCHERR_LIBRARY} ${MUMPS_PORD_LIBRARY}
                        ${MUMPS_SCALAPACK_LIBRARY})
    set(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})
    
    if(NOT (MUMPS_LIBRARIES AND  MUMPS_INCLUDE_DIRS)) #Mimic "REQUIRED" behaviour
        message(FATAL_ERROR "Mumps is missing from ParMooN.")
    else()
        message(STATUS "Done.")
    endif()
    
    include_directories(${MUMPS_INCLUDE_DIRS})  
    
  message("==================Including MUMPS - end====================================================")
endif(_USING_MPI)

#Include Metis
message("==================Including METIS - begin==================================================")
    find_path(METIS_INCLUDE_DIR  metis.h PATHS ${PARMOON_EXTLIB_PATH}/Metis/Lib NO_DEFAULT_PATH)
    find_library(METIS_LIBRARY NAMES metis_${_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Metis NO_DEFAULT_PATH)
    
    set(METIS_LIBRARIES ${METIS_LIBRARY} )
    set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR} )
    if(NOT (METIS_LIBRARIES AND  METIS_INCLUDE_DIRS)) #Mimic "REQUIRED" behaviour
        message(FATAL_ERROR "Metis is missing from ParMooN.")
    else()
        message(STATUS "Done.")
    endif()
    include_directories(${METIS_INCLUDE_DIRS})
message("==================Indcluding METIS - end====================================================")

# Search and include TRIANGLE library.
message("==================Including TRIANGLE - begin================================================")
    find_path(TRIANGLE_INCLUDE_DIR  triangle.h PATHS ${PARMOON_EXTLIB_PATH}/Triangle NO_DEFAULT_PATH)
    find_library(TRIANGLE_LIBRARY NAMES triangle_${_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Triangle NO_DEFAULT_PATH)
    
    set(TRIANGLE_LIBRARIES ${TRIANGLE_LIBRARY})
    set(TRIANGLE_INCLUDE_DIRS ${TRIANGLE_INCLUDE_DIR})
    if(NOT (TRIANGLE_LIBRARIES AND  TRIANGLE_INCLUDE_DIRS)) #Mimic "REQUIRED" behaviour
        message(FATAL_ERROR "Triangle is mising from ParMooN.")
    else()
        message(STATUS "Done.")
    endif()
	include_directories(${TRIANGLE_INCLUDE_DIRS})
message("==================Including TRIANGLE - end==================================================")

# Search and include TETGEN library.
message("==================Including TETGEN - begin================================================")
    find_path(TETGEN_INCLUDE_DIR  tetgen.h PATHS ${PARMOON_EXTLIB_PATH}/tetgen)
    find_library(TETGEN_LIBRARY NAMES tet_${_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/tetgen)
    
    set(TETGEN_LIBRARIES ${TETGEN_LIBRARY})
    set(TETGEN_INCLUDE_DIRS ${TETGEN_INCLUDE_DIR})
    if(NOT (TETGEN_LIBRARIES AND  TETGEN_INCLUDE_DIRS)) #Mimic "REQUIRED" behaviour
        message(FATAL_ERROR "Tetgen is missing from ParMooN.")
    else()
        message(STATUS "Done.")
    endif()
    include_directories(${TETGEN_INCLUDE_DIRS}) 
message("==================Including TETGEN - end==================================================")

# Include TECPLOT library.
message("==================Including TECPLOT - begin================================================")
    find_path(TECPLOT_INCLUDE_DIR  TECIO.h PATHS ${PARMOON_EXTLIB_PATH}/tecplot/include NO_DEFAULT_PATH)
    find_library(TECPLOT_LIBRARY NAMES tecio_${_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/tecplot/lib NO_DEFAULT_PATH)
    
    set(TECPLOT_LIBRARIES ${TECPLOT_LIBRARY})
    set(TECPLOT_INCLUDE_DIRS ${TECPLOT_INCLUDE_DIR})
    if(NOT (TECPLOT_LIBRARIES AND  TECPLOT_INCLUDE_DIRS)) #Mimic "REQUIRED" behaviour
        message(FATAL_ERROR "Tecplot is missing from ParMooN.")
    else()
        message(STATUS "Done.")
    endif()
	include_directories(${TECPLOT_INCLUDE_DIRS})
message("==================Including TECPLOT - end==================================================")

# Include UMFPACK library.
message("==================Including UMFPACK - begin================================================")
    find_path(UMFPACK_INCLUDE_DIR  umfpack.h PATHS ${PARMOON_EXTLIB_PATH}/UMFPACK/Include NO_DEFAULT_PATH)
    find_library(UMFPACK_LIBRARY NAMES umfpack_${_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/UMFPACK/Lib NO_DEFAULT_PATH) 
    get_filename_component(_UMFPACK_LIBDIR ${UMFPACK_LIBRARY} PATH)
    find_library(UMFPACK_SUITESE_LIBRARY NAMES suitesparseconfig_${_ARCH} PATHS ${_UMFPACK_LIBDIR} NO_DEFAULT_PATH)
    find_library(UMFPACK_AMD_LIBRARY NAMES amd_${_ARCH} PATHS ${_UMFPACK_LIBDIR} NO_DEFAULT_PATH)
    
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARY} ${UMFPACK_SUITESE_LIBRARY} ${UMFPACK_AMD_LIBRARY})
    set(UMFPACK_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIR})
    if(NOT (UMFPACK_LIBRARIES AND  UMFPACK_INCLUDE_DIRS)) #Mimic "REQUIRED" behaviour
        message(FATAL_ERROR "Umfpack is missing from ParMooN.")
    else()
        message(STATUS "Done.")
    endif()
    include_directories(${UMFPACK_INCLUDE_DIRS})
message("==================Including UMFPACK - end==================================================")

# turn on all warnings
IF(CMAKE_BUILD_TOOL MATCHES "make")
  IF(NOT CMAKE_CXX_FLAGS MATCHES "-Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  ENDIF(NOT CMAKE_CXX_FLAGS MATCHES "-Wall")
  IF(NOT CMAKE_C_FLAGS MATCHES "-Wall")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
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
include(cmake/sources_amg.cmake)
include(cmake/sources_examples.cmake)
include(cmake/sources_fe.cmake)
include(cmake/sources_general.cmake)
include(cmake/sources_geometry.cmake)
include(cmake/sources_parallel.cmake)
include(cmake/sources_quadformulas.cmake)
include(cmake/sources_refinement.cmake)
include(cmake/sources_system.cmake)
include(cmake/sources_matrix.cmake)
# include(cmake/sources_pbe.cmake)

# Define two static libraries. Each has its own precompiler flag (2D/3D). 
add_library(parmoon_2d STATIC ${PARMOON_SOURCES_2D})
target_compile_definitions(parmoon_2d PUBLIC -D__2D__)

add_library(parmoon_3d STATIC ${PARMOON_SOURCES_3D})
target_compile_definitions(parmoon_3d PUBLIC -D__3D__)

# Set up list of extern libraries, alphabetically. Mutual dependencies are listed locally.

LIST(APPEND _EXTERN_LIBRARIES ${BLAS_LIBRARIES})
LIST(APPEND _EXTERN_LIBRARIES ${TRIANGLE_LIBRARIES})
LIST(APPEND _EXTERN_LIBRARIES ${LAPACK_LIBRARIES})
LIST(APPEND _EXTERN_LIBRARIES ${METIS_LIBRARIES})
LIST(APPEND _EXTERN_LIBRARIES ${MUMPS_LIBRARIES})
LIST(APPEND _EXTERN_LIBRARIES ${PARMETIS_LIBRARIES})
LIST(APPEND _EXTERN_LIBRARIES ${TECPLOT_LIBRARIES})
LIST(APPEND _EXTERN_LIBRARIES ${TETGEN_LIBRARIES})
LIST(APPEND _EXTERN_LIBRARIES ${UMFPACK_LIBRARIES})
if(_USING_MPI)
    LIST(APPEND _EXTERN_LIBRARIES ${_MPI_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endif(_USING_MPI)


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
target_link_libraries(cd2d parmoon_2d ${_EXTERN_LIBRARIES})

### Darcy2D ###
# Standard executable target.
add_executable(darcy2d ${PROJECT_SOURCE_DIR}/2DPrograms/Darcy2D_ParMooN.C)
set_target_properties(darcy2d PROPERTIES OUTPUT_NAME parMooN_darcy2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(darcy2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(darcy2d parmoon_2d ${_EXTERN_LIBRARIES})

### NSE2D ###
# Standard executable target.
add_executable(nse2d ${PROJECT_SOURCE_DIR}/2DPrograms/NSE2D_ParMooN.C)
set_target_properties(nse2d PROPERTIES OUTPUT_NAME parMooN_nse2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(nse2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(nse2d parmoon_2d ${_EXTERN_LIBRARIES})

### TCD2D ###
# Standard executable target.
add_executable(tcd2d ${PROJECT_SOURCE_DIR}/2DPrograms/TCD2D_ParMooN.C)
set_target_properties(tcd2d PROPERTIES OUTPUT_NAME parMooN_tcd2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tcd2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(tcd2d parmoon_2d ${_EXTERN_LIBRARIES})

### TNSE2D ###
# Standard executable target.
add_executable(tnse2d ${PROJECT_SOURCE_DIR}/2DPrograms/TNSE2D_ParMooN.C)
set_target_properties(tnse2d PROPERTIES OUTPUT_NAME parMooN_tnse2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tnse2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(tnse2d parmoon_2d ${_EXTERN_LIBRARIES})

### Standard 3D Programs ###

### CD3D ###
# Standard executable target.
add_executable(cd3d ${PROJECT_SOURCE_DIR}/3DPrograms/CD3D_ParMooN.C)
set_target_properties(cd3d PROPERTIES OUTPUT_NAME parMooN_cd3d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(cd3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(cd3d parmoon_3d ${_EXTERN_LIBRARIES})
# Set a fortran linker flag TODO This is architecture specific, find a portable workaround!
set_target_properties(cd3d PROPERTIES LINK_FLAGS "-lgfortran")


### NSE3D ###
# Standard executable target.
add_executable(nse3d ${PROJECT_SOURCE_DIR}/3DPrograms/NSE3D_ParMooN.C)
set_target_properties(nse3d PROPERTIES OUTPUT_NAME parMooN_nse3d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(nse3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(nse3d parmoon_3d ${_EXTERN_LIBRARIES})
# Set a fortran linker flag TODO This is architecture specific, find a portable workaround!
set_target_properties(nse3d PROPERTIES LINK_FLAGS "-lgfortran")

### TCD3D ###
# Standard executable target.
add_executable(tcd3d ${PROJECT_SOURCE_DIR}/3DPrograms/TCD3D_ParMooN.C)
set_target_properties(tcd3d PROPERTIES OUTPUT_NAME parMooN_tcd3d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tcd3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(tcd3d parmoon_3d ${_EXTERN_LIBRARIES})
# Set a fortran linker flag TODO This is architecture specific, find a portable workaround!
set_target_properties(tcd3d PROPERTIES LINK_FLAGS "-lgfortran")


###############################################################################
# building the documentation using doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/documentation/Doxyfile.in
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
                    ${DOXYGEN_EXECUTABLE} 
                    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generating documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)


add_subdirectory(user_projects)

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




