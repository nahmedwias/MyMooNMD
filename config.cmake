# ===================================================================
# This is basic configuration file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 05 June 2015
# please do not modify unless if it necessary
# run ccmake or cmake-gui to select the options
# any variable value can be overwritten if it is not available in drop-down menu
# to add value other than the one in drop-down menu use: cmake -DMACHINE:STRING=BLUEGENE 
# ===================================================================

# selection of dimension
set(GEO "2D" CACHE STRING 
   "If you change GEO, configure immediately [c] before proceeding further as all other variables depends on this choice !!!")
set_property(CACHE GEO PROPERTY STRINGS 1D 2D 3D) 

# selection of architect type
set(ARCH "MAC64" CACHE STRING "select the machine type")
set_property(CACHE ARCH PROPERTY STRINGS LINUX64 MAC64 INTEL64 TYRONE64 CRAY64)       
 
#  selection of program type
set(PARALLEL_TYPE "SEQUENTIAL" CACHE STRING "select the parallel type")
set_property(CACHE PARALLEL_TYPE PROPERTY STRINGS SEQUENTIAL MPI OMPONLY HYBRID )

# selection of main program
if("${GEO}" STREQUAL "2D")
  set(MODEL " " CACHE STRING "Enter to select the Main file of the model")
  file(GLOB_RECURSE MAIN_SOURCES "${PROJECT_SOURCE_DIR}/2DPrograms/*.C")
  set_property(CACHE MODEL PROPERTY STRINGS  ${MAIN_SOURCES})   
elseif("${GEO}" STREQUAL "3D")
  set(MODEL "${PROJECT_SOURCE_DIR}/3DPrograms/CD3D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")
  file(GLOB_RECURSE MAIN_SOURCES "${PROJECT_SOURCE_DIR}/3DPrograms/*.C")
  set_property(CACHE MODEL PROPERTY STRINGS  ${MAIN_SOURCES})  
endif()  

# selection of BLAS
set(EXTLIB "${PROJECT_SOURCE_DIR}/EXT_LIB")
set(BLAS_OWN " " CACHE STRING "Enter your BLAS lib, it will appear in BLAS_USED option after Pressing [c]: Eg. PathToBlasLibFolder/ -lbblas")
set(ACML_PATH "${EXTLIB}/ACML/gfortran64/lib")

if("${ARCH}" STREQUAL "MAC64")
  set(AVAILABLE_BLAS "-framework Accelerate -lpthread" "BLAS for  ${ARCH} Not Available in ParMooN, use system or enter your BLAS in BLAS_OWN" ${BLAS_OWN})
  set(BLAS_USED "-framework Accelerate -lpthread" CACHE STRING "Select the BLAS")
  set_property(CACHE BLAS_USED PROPERTY STRINGS  ${AVAILABLE_BLAS})   
elseif("${ARCH}" STREQUAL "LINUX64")
  set(AVAILABLE_BLAS "-L${MKLROOT}/lib/intel64 -lmkl_rt -lgfortran   -lpthread" "-L${ACML_PATH} -lacml -lacml_mp -lgfortran  -lpthread" ${BLAS_OWN})
  set(BLAS_USED "-L$(ACML_PATH) -lacml -lacml_mp -lgfortran  -lpthread" CACHE STRING "Select the BLAS")
  set_property(CACHE BLAS_USED PROPERTY STRINGS  ${AVAILABLE_BLAS})  
endif()  


# selection of Triangle lib mesh generator
set(TRIANGLE_OWN " " CACHE STRING "Enter your Triangle lib, it will appear in TRIANGLE_USED option after Pressing [c]: Eg. PathToLibFolder/gridgen_${ARCH}.o")
set(AVAILABLE_TRILIBS "${EXTLIB}/GridGen/gridgen_${ARCH}.o" ${TRIANGLE_OWN})
set(TRIANGLE_USED "${EXTLIB}/GridGen/gridgen_${ARCH}.o" CACHE STRING "Select the Triangle lib")
set_property(CACHE TRIANGLE_USED PROPERTY STRINGS  ${AVAILABLE_TRILIBS})   

# selection of Triangle lib mesh generator
set(TETGEN_OWN " " CACHE STRING "Enter your Triangle lib, it will appear in TETGEN_USED option after Pressing [c]: Eg. -L/PathToLibFolder/ -ltet_${ARCH}")
set(AVAILABLE_TETGENLIBS "-L${EXTLIB}/tetgen -ltet_${ARCH}" ${TETGEN_OWN})
set(TETGEN_USED "-L${EXTLIB}/tetgen -ltet_${ARCH}" CACHE STRING "Select the Tetgen lib")
set_property(CACHE TETGEN_USED PROPERTY STRINGS  ${AVAILABLE_TETGENLIBS})   


# selection of UMFPACK lib mesh generator
set(UMFPACK_OWN " " CACHE STRING "Enter your UMFPACK lib, it will appear in UMFPACK_USED option after Pressing [c]: Eg. -L/PathToLibFolder/ -lumfpack_${ARCH} ...")
set(UMPACK_PARMOON "-L${EXTLIB}/UMFPACK/Lib -lumfpack_${ARCH} -lamd_${ARCH} -lsuitesparseconfig_${ARCH} -lcholmod_${ARCH} -lcolamd_${ARCH} -lrt")
set(AVAILABLE_UMFPACKLIBS  ${UMPACK_PARMOON} ${UMFPACK_OWN})
set(UMFPACK_USED ${UMPACK_PARMOON} CACHE STRING "Select the UMFPACK lib")
set_property(CACHE UMFPACK_USED PROPERTY STRINGS  ${AVAILABLE_UMFPACKLIBS})  

# selection of Teclot lib
set(TECPLOT_OWN " " CACHE STRING "Enter your Tecplot lib, it will appear in TECPLOT_USED option after Pressing [c]: Eg. -L/PathToLibFolder/ -lumfpack_${ARCH}")
set(TECPLOT_PARMOON "-L${EXTLIB}/tecplot/lib -ltecio_${ARCH}")
set(AVAILABLE_TECPLOTLIBS  ${TECPLOT_PARMOON} ${TECPLOT_OWN})
set(TECPLOT_USED ${TECPLOT_PARMOON} CACHE STRING "Select the TECPLOT lib")
set_property(CACHE TECPLOT_USED PROPERTY STRINGS  ${AVAILABLE_TECPLOTLIBS}) 


# selecting the MUMPS lib

if("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
# first set the blacs for MUMPS_PARMOON
if("${ARCH}" STREQUAL "MAC64")
 set(MPIBLACS_${ARCH} "-L${EXTLIB}/MPIBLACS -lscalapack_${ARCH} -llapack_${ARCH}")
elseif("${ARCH}" STREQUAL "LINUX64")
 set(MPIBLACS_${ARCH} "-L${EXTLIB}/MPIBLACS   -lscalapack_${ARCH}  -lblacs_MPI_${ARCH}-0 -lblacsF77init_MPI_${ARCH}-0  -lblacs_MPI_${ARCH}-0    -L/opt/mpich/lib64  -lmpi -lmpifort -lmpicxx")
endif()

set(MUMPS_OWN " " CACHE STRING "Enter your MUMPS lib including all blacs/Metis/scotch etc, it will appear in MUMPS_USED option after Pressing [c]")
set(MUMPS_PARMOON_${ARCH} "-L${EXTLIB}/MUMPS/lib -ldmumps_${ARCH} -lmumps_common_${ARCH} -L${EXTLIB}/Metis  -lparmetis_${ARCH} -lmetis_${ARCH} -L${EXTLIB}/MUMPS/lib -lpord_${ARCH} -lptesmumps_${ARCH} -lptscotch_${ARCH} -lscotch_${ARCH} -lptscotcherrexit_${ARCH} -lptscotchparmetis_${ARCH} ${MPIBLACS_${ARCH}}")

set(AVAILABLE_MUMPSLIBS  ${MUMPS_PARMOON_${ARCH}} ${MUMPS_OWN})
set(MUMPS_USED ${MUMPS_PARMOON_${ARCH}} CACHE STRING "Select the MUMPS lib")
set_property(CACHE MUMPS_USED PROPERTY STRINGS  ${AVAILABLE_MUMPSLIBS}) 

elseif("${PARALLEL_TYPE}" STREQUAL "SEQUENTIAL"  OR "${PARALLEL_TYPE}" STREQUAL "OMPONLY")

set(MUMPS_USED " " CACHE STRING "No need to select the MUMPS lib")
set_property(CACHE MUMPS_USED PROPERTY STRINGS  " ") 

endif()



set(DESTDIR "${CMAKE_SOURCE_DIR}/ParMooN/OutPut" CACHE STRING "select the model")


# select the compiler type,
# CMAKE_BUILD_TYPE [DEBUG|RELEASE|RELWITHDEBINFO|MINSIZEREL]
set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "options")

include (CMakeForceCompiler)
  set(CMAKE_SYSTEM_NAME ${ARCH})
  
if("${ARCH}" STREQUAL "MAC64")

  CMAKE_FORCE_C_COMPILER   (clang Clang)
  CMAKE_FORCE_CXX_COMPILER (clang++ Clang)
  set(CMAKE_C_FLAGS "-O3 -D${GEO} -DTRIMACHINE -DREDUCED -DNO_TIMER $(INC) -DMKL_ILP64 -m64  -framework -fapple-pragma-pack   -Wno-unused-parameter")
  set(CMAKE_C_FLAGS_DEBUG "-D${GEO} -Weverything  -Wno-unused-parameter")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -D${GEO} -DTRIMACHINE -DREDUCED -DNO_TIMER $(INC) -DMKL_ILP64 -m64 -fapple-pragma-pack  -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS "-D${GEO} -fapple-pragma-pack  -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_DEBUG " -D${GEO} -Weverything -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -D${GEO} -fapple-pragma-pack -Wno-unused-parameter") 
  
elseif("${ARCH}" STREQUAL "LINUX64")

if("${PARALLEL_TYPE}" STREQUAL "MPI"  OR "${PARALLEL_TYPE}" STREQUAL "HYBRID")
  CMAKE_FORCE_C_COMPILER   (mpicc)
  CMAKE_FORCE_CXX_COMPILER (mpicxx)
elseif("${PARALLEL_TYPE}" STREQUAL "SEQUENTIAL"  OR "${PARALLEL_TYPE}" STREQUAL "OMPONLY")
  CMAKE_FORCE_C_COMPILER   (gcc)
  CMAKE_FORCE_CXX_COMPILER (g++)
endif()  
  
  set(CMAKE_C_FLAGS "-D${GEO}  -std=c++11 -fopenmp")
  set(CMAKE_C_FLAGS_DEBUG "-g -D${GEO} $(INC) -std=c++11 -fopenmp")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -D${GEO} -s $(INC) -std=c++11 -fopenmp")
  set(CMAKE_CXX_FLAGS "-D${GEO} -DTRILIBRARY -DREDUCED -DNO_TIMER $(INC) -DMKL_ILP64 -m64 -fopenmp")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -D${GEO} -DTRILIBRARY -DREDUCED -DNO_TIMER $(INC) -DMKL_ILP64 -m64  -fopenmp")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -D${GEO} -s -DTRILIBRARY -DREDUCED -DNO_TIMER $(INC) -DMKL_ILP64 -m64 -fopenmp") 
  
endif()

 
#  SET(PART_NAME SAMPLE_PART)
# SET(PART_TARGET_NAME samplepart)
# 
# SET(BUILD_${PART_NAME} off CACHE BOOL "Build part ${PART_NAME}")
# SET(${PART_NAME}_FILES some_file.c another_file.c)
# 
# IF(NOT BUILD_${PART_NAME})
# SET(EXCLUDE_${PART_NAME} EXCLUDE_FROM_ALL)
# ENDIF(NOT BUILD_${PART_NAME})
# 
# ADD_LIBRARY(${PART_TARGET_NAME} SHARED ${EXCLUDE_SAMPLE_PART} ${PART_NAME}_FILES)
# TARGET_LINK_LIBRARY(${PART_TARGET_NAME} some_other_target)

 