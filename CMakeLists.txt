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

include(${CMAKE_SOURCE_DIR}/UserConfig.cmake)

project(ParMooN)
set(ParMooN_VERSION_MAJOR 1)
set(ParMooN_VERSION_MINOR 1)

# Enable Fortran compiling, linking and mangling.
enable_language(Fortran)
include(FortranCInterface) 

# Set two auxiliary variables, which make the file better readable imho
if("${_PARALLEL_TYPE}" STREQUAL "MPI" OR "${_PARALLEL_TYPE}" STREQUAL "HYBRID")
  set(_USING_MPI TRUE)
endif()
if("${_PARALLEL_TYPE}" STREQUAL "OMPONLY" OR "${_PARALLEL_TYPE}" STREQUAL "HYBRID")
  set(_USING_OMP TRUE)
endif()

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

# Set output path for binaries. Make sure the directory exists!
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${_OUTPUT_DIRECTORY})

# Set path to the external libraries which ship with ParMooN
set(PARMOON_EXTLIB_PATH ${PROJECT_SOURCE_DIR}/EXT_LIB CACHE PATH "Directory path to ParMooN extern libraries.")

# Set path to ParMooN package search modules.
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/ CACHE PATH "Directory path to ParMooN CMake modules.")

# Search and include LAPACK library. 
# Note: The default lapack search module searches for a BLAS library, too, and eventually sets BLAS_FOUND
# to false, even if our module (below) already found a BLAS library, which leads to nasty linker errors.
# Therefore LAPACK is searched before BLAS, which is basically a hack and should be changed someday TODO.
message("==================Searching LAPACK - begin=================================================")
  
	  #if(DEFAULT_MODULE_LAPACK)
	  #  message("Using default CMake module to find LAPACK.")
	  #  find_package(LAPACK MODULE)
	  #  if(LAPACK_FOUND)
	  #    message("LAPACK found with default CMake module.")
	  #	  include_directories(${LAPACK_INCLUDE_DIRS})
	  #  endif(LAPACK_FOUND)  
	  #endif(DEFAULT_MODULE_LAPACK)
	  
	  if(NOT LAPACK_FOUND)
	    message("Using ParMooN CMake module to find LAPACK.")
	    find_package(LAPACKLIB REQUIRED MODULE)
	    message("LAPACK found with ParMooN CMake module.")
	    #include_directories(${LAPACK_INCLUDE_DIRS}) #CB Our module does not find lapack headers!
	  endif(NOT LAPACK_FOUND)
  
message("==================Searching LAPACK - end===================================================")  

# Search and include LAPACK library.
message("==================Searching BLAS - begin===================================================")

	#if(DEFAULT_MODULE_BLAS)
	#  message("Using default CMake module to find BLAS.")
	#    find_package(BLAS MODULE)
	#    if(BLAS_FOUND)
	#      message("BLAS found with default CMake module.")
	#  	  include_directories(${BLAS_INCLUDE_DIRS})
	#    endif(BLAS_FOUND)
	#endif(DEFAULT_MODULE_BLAS)
	
	if(NOT BLAS_FOUND)
	  message("Using ParMooN CMake module to find BLAS.")
	  find_package(ACMLBLAS REQUIRED MODULE)
	  message("BLAS found with ParMooN CMake module.")
	  include_directories(${BLAS_INCLUDE_DIRS})
	endif(NOT BLAS_FOUND)

message("==================Searching BLAS - end=====================================================")

# When using MPI, find mpi
if(_USING_MPI)  
  #Kickoff the default MPI search module and use its output.
  message("==================Searching MPI - begin====================================================")
  find_package(MPI MODULE REQUIRED)
  if(MPI_FOUND)
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
    #Combine MPI libraries for all 3 languages to one list.
    list(APPEND _MPI_LIBRARIES ${MPI_Fortran_LIBRARIES})
    list(APPEND _MPI_LIBRARIES ${MPI_C_LIBRARIES})
    list(APPEND _MPI_LIBRARIES ${MPI_CXX_LIBRARIES})
  endif(MPI_FOUND)
  message("==================Searching MPI - end======================================================")
  
  # Find a thread library, too.
  message("==================Searching Thread Library - begin=========================================")
  set(CMAKE_THREAD_PREFER_PTHREAD ON)
  find_package(Threads REQUIRED)
  message("==================Searching Thread Library - end===========================================")
  
  # When using MPI we're also using Parmetis.
  message("==================Searching Parmetis - begin===============================================")
      find_package(PARMETIS MODULE REQUIRED)
      include_directories(${PARMETIS_INCLUDE_DIRS})
   message("==================Searching Parmetis - end=================================================")
   
  # When using MPI we're also using MUMPS solver.
  message("==================Searching MUMPS - begin==================================================")
    message("Using ParMooN CMake module to find MUMPS.")
    find_package(MUMPS REQUIRED MODULE)
	message("MUMPS found with ParMooN CMake module.")
	include_directories(${MUMPS_INCLUDE_DIRS})
  message("==================Searching MUMPS - end====================================================")
endif(_USING_MPI)

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

# Search and include Metis
message("==================Searching METIS - begin==================================================")
    find_package(METIS MODULE REQUIRED)
    include_directories(${METIS_INCLUDE_DIRS})
 message("==================Searching METIS - end====================================================")
 
# Search and include TRIANGLE library (which is actually Shewchuk's "Triangle").
message("==================Searching TRIANGLE - begin================================================")
	message("Using ParMooN CMake module to find TRIANGLE.")
	find_package(TRIANGLE REQUIRED MODULE)
	message("TRIANGLE found with ParMooN CMake module.")
	include_directories(${TRIANGLE_INCLUDE_DIRS})
message("==================Searching TRIANGLE - end==================================================")

# Search and include TETGEN library.
message("==================Searching TETGEN - begin================================================")
    message("Using ParMooN CMake module to find TETGEN.")
	find_package(TETGEN REQUIRED MODULE)
	message("TETGEN found with ParMooN CMake module.")
    include_directories(${TETGEN_INCLUDE_DIRS}) 
message("==================Searching TETGEN - end==================================================")

# Search and include TECPLOT library.
message("==================Searching TECPLOT - begin================================================")
	message("Using ParMooN CMake module to find TECPLOT.")
	find_package(TECPLOT REQUIRED MODULE)
	message("TECPLOT found with ParMooN CMake module.")
	include_directories(${TECPLOT_INCLUDE_DIRS})
message("==================Searching TECPLOT - end==================================================")

# Search and include UMFPACK library.
message("==================Searching UMFPACK - begin================================================")
	message("Using ParMooN CMake module to find UMFPACK.")
	find_package(UMFPACK MODULE REQUIRED)
	message("UMFPACK found with ParMooN CMake module.")
	include_directories(${UMFPACK_INCLUDE_DIRS})
message("==================Searching UMFPACK - end==================================================")


# Include directories which contain header files.
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

# Determine precompiler definitions ("PCD") depending on the chosen _PARALLEL_TYPE
if("${_PARALLEL_TYPE}" STREQUAL "MPI")
    set(_PCD_PARALLEL_TYPE -D_PAR -D_MPIONLY -D_MPI)
elseif("${_PARALLEL_TYPE}" STREQUAL "OMPONLY")
    set(_PCD_PARALLEL_TYPE -D_PAR -D_OMPONLY -D_OMP)
elseif("${_PARALLEL_TYPE}" STREQUAL "HYBRID")
    set(_PCD_PARALLEL_TYPE -D_PAR -D_MPI -D_HYBRID)
elseif("${_PARALLEL_TYPE}" STREQUAL "SEQUENTIAL")
    set(_PCD_PARALLEL_TYPE -D_SEQ)
endif()
# Use the directory-based "add_definitions" command to enable those
# ParMooN-specific precompiler definitions.
# NOTE: Should precompiler flags of the same name be used by extern
# libraries which are not in subdirectories, the won't be passed there this way!
add_definitions(${_PCD_PARALLEL_TYPE})

# Determine compiler options depending on _ARCH.
# Is there a way to better localize them?

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
# run the add_subdirectory AFTER setting global flags.
add_subdirectory(Examples)
add_subdirectory(src/AMG)
add_subdirectory(src/General)
add_subdirectory(src/Geometry)          
add_subdirectory(src/FE)
add_subdirectory(src/Parallel)
add_subdirectory(src/QuadFormulas)
add_subdirectory(src/Refinement) 
add_subdirectory(src/System)
# add_subdirectory(src/PBE)

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

# Set up lists of the 2D libraries ready for linking.
LIST(APPEND _PARMOON_LIBRARIES_2D "examples_2d")
LIST(APPEND _PARMOON_LIBRARIES_2D "fe_2d")
LIST(APPEND _PARMOON_LIBRARIES_2D "gen_2d")
LIST(APPEND _PARMOON_LIBRARIES_2D "geo_2d")
LIST(APPEND _PARMOON_LIBRARIES_2D "par_2d")
LIST(APPEND _PARMOON_LIBRARIES_2D "quad_2d")
LIST(APPEND _PARMOON_LIBRARIES_2D "ref_2d")
LIST(APPEND _PARMOON_LIBRARIES_2D "sys_2d")
LIST(APPEND _PARMOON_LIBRARIES_2D "amg")
# Double everything, for the libraries are heavily interdependent.
LIST(APPEND _PARMOON_LIBRARIES_2D ${_PARMOON_LIBRARIES_2D})

# Set up lists of the 3D libraries ready for linking.
LIST(APPEND _PARMOON_LIBRARIES_3D "examples_3d")
LIST(APPEND _PARMOON_LIBRARIES_3D "fe_3d")
LIST(APPEND _PARMOON_LIBRARIES_3D "gen_3d")
LIST(APPEND _PARMOON_LIBRARIES_3D "geo_3d")
LIST(APPEND _PARMOON_LIBRARIES_3D "par_3d")
LIST(APPEND _PARMOON_LIBRARIES_3D "quad_3d")
LIST(APPEND _PARMOON_LIBRARIES_3D "ref_3d")
LIST(APPEND _PARMOON_LIBRARIES_3D "sys_3d")
LIST(APPEND _PARMOON_LIBRARIES_3D "amg")
# Double everything, for the libraries are heavily interdependent.
LIST(APPEND _PARMOON_LIBRARIES_3D ${_PARMOON_LIBRARIES_3D})

##############################################################
# Start adding executable targets.
##############################################################

### Standard 2D Programs ###

### CD2D ###
# Standard executable target.
add_executable(cd2d ${PROJECT_SOURCE_DIR}/2DPrograms/CD2D_ParMooN.C)
set_target_properties(cd2d PROPERTIES OUTPUT_NAME parMooN_cd2d_${_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(cd2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(cd2d ${_PARMOON_LIBRARIES_2D} ${_EXTERN_LIBRARIES})
if(_USING_MPI)
    target_link_libraries(cd2d ${_MPI_LIBRARIES})
endif(_USING_MPI)

### Darcy2D ###
# Standard executable target.
add_executable(darcy2d ${PROJECT_SOURCE_DIR}/2DPrograms/Darcy2D_ParMooN.C)
set_target_properties(darcy2d PROPERTIES OUTPUT_NAME parMooN_darcy2d_${_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(darcy2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(darcy2d ${_PARMOON_LIBRARIES_2D} ${_EXTERN_LIBRARIES})
if(_USING_MPI)
    target_link_libraries(darcy2d ${_MPI_LIBRARIES})
endif(_USING_MPI)

### NSE2D ###
# Standard executable target.
add_executable(nse2d ${PROJECT_SOURCE_DIR}/2DPrograms/NSE2D_ParMooN.C)
set_target_properties(nse2d PROPERTIES OUTPUT_NAME parMooN_nse2d_${_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(nse2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(nse2d ${_PARMOON_LIBRARIES_2D} ${_EXTERN_LIBRARIES})
if(_USING_MPI)
    target_link_libraries(nse2d ${_MPI_LIBRARIES})
endif(_USING_MPI)

### TCD2D ###
# Standard executable target.
add_executable(tcd2d ${PROJECT_SOURCE_DIR}/2DPrograms/TCD2D_ParMooN.C)
set_target_properties(tcd2d PROPERTIES OUTPUT_NAME parMooN_tcd2d_${_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tcd2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(tcd2d ${_PARMOON_LIBRARIES_2D} ${_EXTERN_LIBRARIES})
if(_USING_MPI)
    target_link_libraries(tcd2d ${_MPI_LIBRARIES})
endif(_USING_MPI)

### TNSE2D ###
# Standard executable target.
add_executable(tnse2d ${PROJECT_SOURCE_DIR}/2DPrograms/TNSE2D_ParMooN.C)
set_target_properties(tnse2d PROPERTIES OUTPUT_NAME parMooN_tnse2d_${_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tnse2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(tnse2d ${_PARMOON_LIBRARIES_2D} ${_EXTERN_LIBRARIES})
if(_USING_MPI)
    target_link_libraries(tnse2d ${_MPI_LIBRARIES})
endif(_USING_MPI)

### Standard 3D Programs ###

### CD3D ###
# Standard executable target.
add_executable(cd3d ${PROJECT_SOURCE_DIR}/3DPrograms/CD3D_ParMooN.C)
set_target_properties(cd3d PROPERTIES OUTPUT_NAME parMooN_cd3d_${_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(cd3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(cd3d ${_PARMOON_LIBRARIES_3D} ${_EXTERN_LIBRARIES})
if(_USING_MPI)
    # Link in MPI libraries.
    target_link_libraries(cd3d ${_MPI_LIBRARIES})
    # Link in thread libraries
    target_link_libraries(cd3d ${CMAKE_THREAD_LIBS_INIT})
endif(_USING_MPI)
# Set a fortran linker flag TODO This is architecture specific, find a portable workaround!
set_target_properties(cd3d PROPERTIES LINK_FLAGS "-lgfortran")

### NSE3D ###
# Standard executable target.
add_executable(nse3d ${PROJECT_SOURCE_DIR}/3DPrograms/NSE3D_ParMooN.C)
set_target_properties(nse3d PROPERTIES OUTPUT_NAME parMooN_nse3d_${_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(nse3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(nse3d ${_PARMOON_LIBRARIES_3D} ${_EXTERN_LIBRARIES})
if(_USING_MPI)
    # Link in MPI libraries.
    target_link_libraries(nse3d ${_MPI_LIBRARIES})
    # Link in thread libraries
    target_link_libraries(nse3d ${CMAKE_THREAD_LIBS_INIT})
endif(_USING_MPI)
# Set a fortran linker flag TODO This is architecture specific, find a portable workaround!
set_target_properties(nse3d PROPERTIES LINK_FLAGS "-lgfortran")

### TCD3D ###
# Standard executable target.
add_executable(tcd3d ${PROJECT_SOURCE_DIR}/3DPrograms/TCD3D_ParMooN.C)
set_target_properties(tcd3d PROPERTIES OUTPUT_NAME parMooN_tcd3d_${_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tcd3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(tcd3d ${_PARMOON_LIBRARIES_3D} ${_EXTERN_LIBRARIES})
if(_USING_MPI)
    # Link in MPI libraries.
    target_link_libraries(tcd3d ${_MPI_LIBRARIES})
    # Link in thread libraries
    target_link_libraries(tcd3d ${CMAKE_THREAD_LIBS_INIT})
endif(_USING_MPI)
    # Set a fortran linker flag TODO This is architecture specific, find a portable workaround!
    set_target_properties(tcd3d PROPERTIES LINK_FLAGS "-lgfortran")


############
# End adding executable targets.
############

message("============================================================================================================"  )        
message("Configuration for the ParMooN Version 1.1")  
message("Author : Sashikumaar Ganesan")
message("         Numerical Mathematics and Scientific Computing")
message("         Supercomputer Edu. and Research Centre, Indian Institute of Science")
message("         Bangalore, India")  
message("Contact: mailt0:sashi@serc.iisc.in") 
message("============================================================================================================"  ) 




