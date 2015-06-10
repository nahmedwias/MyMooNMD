# ===================================================================
# This is a user configuration file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 05 June 2015
# ===================================================================

# controlling the output messages
set(CMAKE_VERBOSE_MAKEFILE FALSE)

# selection of dimension (2D 3D)
set(GEO "3D" CACHE STRING "Change GEO, to select the Dimensio of the problem")

# select this line accordingly to include your main program
# set(MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/CD2D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model") 
set(MODEL "${PROJECT_SOURCE_DIR}/3DPrograms/CD3D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model") 

# selection of architect type (LINUX64 MAC64 INTEL64 TYRONE64 CRAY64)
set(ARCH "MAC64" CACHE STRING "select the machine type")  

#  selection of program type (SEQUENTIAL MPI OMPONLY HYBRID)
set(PARALLEL_TYPE "MPI" CACHE STRING "select the parallel type")

#  set MORTAR, if needed
set(MORTAR " ")

# set the path to save the exe file
set(DESTDIR "${CMAKE_SOURCE_DIR}/../ParMooN_Output/cd3d" CACHE STRING "select the model")

# ========================================================================================================================
# no need to change anyting after this line
# used only when ccmake or cmake-gui is used
# ========================================================================================================================
set_property(CACHE GEO PROPERTY STRINGS 1D 2D 3D ) 

# selection of all main programs
if("${GEO}" STREQUAL "2D")
  file(GLOB_RECURSE MAIN_SOURCES "${PROJECT_SOURCE_DIR}/2DPrograms/*.C")
  set_property(CACHE MODEL PROPERTY STRINGS  ${MAIN_SOURCES})   
elseif("${GEO}" STREQUAL "3D")
  file(GLOB_RECURSE MAIN_SOURCES "${PROJECT_SOURCE_DIR}/3DPrograms/*.C")
  set_property(CACHE MODEL PROPERTY STRINGS  ${MAIN_SOURCES})  
endif()  

# selection of all architects
set_property(CACHE ARCH PROPERTY STRINGS LINUX64 MAC64 INTEL64 TYRONE64 CRAY64)

# selection of all program types
set_property(CACHE PARALLEL_TYPE PROPERTY STRINGS SEQUENTIAL MPI OMPONLY HYBRID)

# ======================================================================
# general settings
# ======================================================================
 if("${PARALLEL_TYPE}" STREQUAL "MPI")
   set(PARMOON_PRG_DEFINE "-D_PAR -D_MPIONLY -D_MPI")
 elseif("${PARALLEL_TYPE}" STREQUAL "OMPONLY")
   set(PARMOON_PRG_DEFINE "-D_PAR -D_OMPONLY -D_OMP")
 elseif("${PARALLEL_TYPE}" STREQUAL "HYBRID")
   set(PARMOON_PRG_DEFINE "-D_PAR -D_MPI -D_HYBRID")
 elseif("${PARALLEL_TYPE}" STREQUAL "SEQUENTIAL")
   set(PARMOON_PRG_DEFINE "-D_SEQ")
 endif()

 if("${ARCH}" STREQUAL "LINUX64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF}  -DTRILIBRARY -DREDUCED -DNO_TIMER -DMKL_ILP64 -m64 -fopenmp")
 elseif("${ARCH}" STREQUAL "MAC64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF}  -DREDUCED -DNO_TIMER -DMKL_ILP64 -m64 -fapple-pragma-pack -Wmacro-redefined -Wdangling-else
                         -Wcomment -Wparentheses-equality -Wdelete-non-virtual-dtor -Wnull-conversion")
 elseif("${ARCH}" STREQUAL "INTEL64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF}  ")
 elseif("${ARCH}" STREQUAL "TYRONE64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} -DTRILIBRARY -DREDUCED -DNO_TIMER -DMPICH_IGNORE_CXX_SEEK")
 elseif("${ARCH}" STREQUAL "CRAY64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} ")  
 endif()
 
 
set(PARMOON_CXX_DEF "-D__${GEO}__ -D__${ARCH}__ -DTETLIBRARY -D__PRIVATE__ ${MORTAR} ${PARMOON_PRG_DEFINE}")
