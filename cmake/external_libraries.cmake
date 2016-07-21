# Set path to the external libraries which ship with ParMooN
set(PARMOON_EXTLIB_PATH ${CMAKE_SOURCE_DIR}/EXT_LIB)
# Set path to ParMooN package search modules.
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/PETSc)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules/)


# find PETSc
# what we want is 'find_package(PETSc REQUIRED)' but with a custom error 
# message. To do so we omit the 'REQUIRED' and then check for success by hand.
set(PETSC_DIR CACHE STRING "The path to the directory to find PETSc in.")
set(PETSC_ARCH CACHE STRING
    "The directory (under PETSC_DIR) to find the compiled PETSc in.")
find_package(PETSc)
if(NOT PETSC_FOUND)
  message(WARNING "could not find PETSc. This is essential for ParMooN.")
  message(STATUS "If you don't have PETSc installed, then you can download "
          "PETSc via ")
  message(STATUS "    git clone -b maint https://bitbucket.org/petsc/petsc "
          "petsc")
  message(STATUS "into the subdirectory petsc (preferably outside of the ")
  message(STATUS "ParMooN directory). Then in that directory you have to run ")
  message(STATUS "the configure script. As an example try ")
  message(STATUS "    ./configure --with-cc=gcc --with-cxx=g++ \
--with-fc=gfortran \
--with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native' \
CXXOPTFLAGS='-O3 -march=native -mtune=native' \
FOPTFLAGS='-O3 -march=native -mtune=native' \
--download-fblaslapack --download-openmpi --download-scalapack \
--download-mumps --download-metis --download-suitesparse --download-parmetis")
  message(STATUS "In the end, the output of the configure script will tell you")
  message(STATUS "the command you should call to compile all of PETSc together")
  message(STATUS "with all its dependecies. Then you have to set PETSC_DIR and")
  message(STATUS "PETSC_ARCH correctly within the CMakeCache of your ParMooN")
  message(STATUS "build directory. They are now: ")
  message(STATUS ${PETSC_DIR} " and " ${PETSC_ARCH})
  message(STATUS "Please also read the documentation of PETSc for further")
  message(STATUS "details")
  message(FATAL_ERROR "ParMooN can not be compiled without PETSc.")
endif(NOT PETSC_FOUND)
set(PETSC_LIBS_PATH ${PETSC_DIR}/${PETSC_ARCH}/lib/)
# the following makes all calls to 'include_directories' unnecessary for
# libraries provided by PETSc. Also PETSc needs to find the header file 'mpi.h'
# even in sequential mode (because it is compiled for MPI). With this line it 
# can find it.
include_directories(${PETSC_INCLUDES})


function(check_if_lib_exists lib_name petsc_configure_flag)
  string(FIND "${PETSC_LIBRARIES}" "${lib_name}" found)
  if(found EQUAL -1)
    message(FATAL_ERROR 
            "${lib_name} not found in the PETSc libraries. You have to "
            "configure PETSc with ${lib_name}, e.g.: "
            "--download-${petsc_configure_flag}")
  endif()
endfunction(check_if_lib_exists)


check_if_lib_exists("blas" "fblaslapack")
check_if_lib_exists("lapack" "fblaslapack")
if(_USING_MPI)
  check_if_lib_exists("mpi" "openmpi")
  check_if_lib_exists("metis" "metis")
  check_if_lib_exists("parmetis" "parmetis")
  check_if_lib_exists("dmumps" "mumps")
  check_if_lib_exists("mumps_common" "mumps")
  # Find a thread library. (Mumps needs this)
  set(CMAKE_THREAD_PREFER_PTHREAD ON)
  find_package(Threads REQUIRED)
  # pord library (the sequential ordering tool which ships with mumps)
  check_if_lib_exists("pord" "mumps")
  check_if_lib_exists("scalapack" "scalapack")
  check_if_lib_exists("umfpack" "suitesparse")
  check_if_lib_exists("suitesparseconfig" "suitesparse")
  check_if_lib_exists("amd" "suitesparse")
  check_if_lib_exists("cholmod" "suitesparse")
  check_if_lib_exists("colamd" "suitesparse")
endif(_USING_MPI)



# When using OpenMP, find it and use the module's output.
if(_USING_OMP)
  # CB Kick off the default FindOpenMP module.
  find_package(OpenMP REQUIRED)
  if(OPENMP_FOUND)
    # If OpenMP is used, add the found compiler flags for c and c++ globally.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif(OPENMP_FOUND)
endif(_USING_OMP)

# Search and include TRIANGLE library.
find_path(TRIANGLE_INCLUDE_DIR triangle.h PATHS ${PARMOON_EXTLIB_PATH}/Triangle
          NO_DEFAULT_PATH)
find_library(TRIANGLE_LIBRARY NAMES triangle_${_ARCH} PATHS
             ${PARMOON_EXTLIB_PATH}/Triangle NO_DEFAULT_PATH)
# Mimic "REQUIRED" behaviour
if(NOT (TRIANGLE_LIBRARY AND TRIANGLE_INCLUDE_DIR))
  message(FATAL_ERROR "Could not find Triangle. It should be in the repository")
endif()
include_directories(${TRIANGLE_INCLUDE_DIR})



# Search and include TETGEN library.
find_path(TETGEN_INCLUDE_DIR tetgen.h PATHS ${PARMOON_EXTLIB_PATH}/tetgen)
find_library(TETGEN_LIBRARY NAMES tet_${_ARCH} PATHS
             ${PARMOON_EXTLIB_PATH}/tetgen)
# Mimic "REQUIRED" behaviour
if(NOT (TETGEN_LIBRARY AND TETGEN_INCLUDE_DIR))
  message(FATAL_ERROR "Could not find Tetgen. It should be in the repository")
endif()
include_directories(${TETGEN_INCLUDE_DIR})



# Include pardiso library
# TODO this should be made optional, so it should be possible to compile ParMooN
# even without pardiso but with openmp.
if(_USING_OMP)
  # for some reason the pardiso library is in the mumps directory
  find_library(PARDISO_LIBRARY
               NAMES libpardisoj83.2.9.x86_omp_gfortran_gcc43.a
               PATHS ${PARMOON_EXTLIB_PATH}/pardiso NO_DEFAULT_PATH)
  # Mimic "REQUIRED" behaviour
  if(NOT (PARDISO_LIBRARY))
    message(FATAL_ERROR "pardiso is missing from ParMooN.")
  endif()
endif(_USING_OMP)




# Set up list of external libraries. If library A depends on library B then A
# should be listed before B. For example the Blas library is the last one 
# because it depends on no other library but others depend on it.
list(APPEND _EXTERN_LIBRARIES ${TRIANGLE_LIBRARY})
list(APPEND _EXTERN_LIBRARIES ${TETGEN_LIBRARY})
if(_USING_OMP)
  list(APPEND _EXTERN_LIBRARIES ${PARDISO_LIBRARY})
endif(_USING_OMP)
list(APPEND _EXTERN_LIBRARIES ${PETSC_LIBRARIES})
