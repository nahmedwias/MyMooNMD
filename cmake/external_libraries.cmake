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
include_directories(${PETSC_DIR}/include)
find_library(PETSC_LIBRARIES NAMES petsc PATHS ${PETSC_LIBS_PATH} 
             NO_DEFAULT_PATH)

# Search and include LAPACK library
# In fact we would like to use 'find_package(LAPACK REQUIRED)' to find LAPACK 
# and BLAS at the same time. However this will not use the PETSc LAPACK and 
# BLAS, so we have to include these by hand.
find_library(LAPACK_LIBRARIES NAMES libflapack.a PATHS ${PETSC_LIBS_PATH}
             NO_DEFAULT_PATH)
# Mimic "REQUIRED" behaviour
if(NOT LAPACK_LIBRARIES)
  message(FATAL_ERROR "Could not find LAPACK. You have to configure PETSc to "
          "include LAPACK")
endif()
find_library(BLAS_LIBRARIES NAMES libfblas.a PATHS ${PETSC_LIBS_PATH}
             NO_DEFAULT_PATH)
# Mimic "REQUIRED" behaviour
if(NOT BLAS_LIBRARIES)
  message(FATAL_ERROR "Could not find BLAS. You have to configure PETSc to "
          "include BLAS")
endif()


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


# When using MPI, find mpi and some other libraries which use mpi.
if(_USING_MPI)
  # ======================== find mpi library =================================
  # we want to use "find_package(MPI)" to find the mpi libraries, compilers, 
  # and includes. However one has to use the same ones which PETSc uses, so it
  # is done by hand here. The find_package command is then used to get correct
  # compiler and linker flags.
  set(MPI_C_COMPILER ${PETSC_DIR}/${PETSC_ARCH}/bin/mpicc)
  set(MPI_CXX_COMPILER ${PETSC_DIR}/${PETSC_ARCH}/bin/mpic++)
  set(MPI_Fortran_COMPILER ${PETSC_DIR}/${PETSC_ARCH}/bin/mpifort)
  set(MPI_C_INCLUDE_PATH ${PETSC_DIR}/${PETSC_ARCH}/include/)
  set(MPI_CXX_INCLUDE_PATH ${PETSC_DIR}/${PETSC_ARCH}/include/)
  set(MPI_Fortran_INCLUDE_PATH ${PETSC_DIR}/${PETSC_ARCH}/include/)
  find_library(MPI_C_LIBRARIES NAMES mpi PATHS PETSC_LIBS_PATH NO_DEFAULT_PATH)
  find_library(MPI_CXX_LIBRARIES NAMES mpi_cxx PATHS PETSC_LIBS_PATH 
               NO_DEFAULT_PATH)
  find_library(MPI_Fortran_LIBRARIES NAMES mpifort PATHS PETSC_LIBS_PATH 
               NO_DEFAULT_PATH)
  
  # call to find_package, this now looks for appropriate compiler/linker flags
  find_package(MPI REQUIRED)
  # set path to the mpi header files.
  include_directories(${MPI_INCLUDE_PATH})
  # set compile and link flags
  if(MPI_COMPILE_FLAGS)
    message(STATUS "There are mpi compile flags!")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} 
        ${MPI_Fortran_COMPILE_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    message(STATUS "There are mpi linker flags!")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
  endif() 
  # Combine MPI libraries for all 3 languages to one list used for linking.
  list(APPEND _MPI_LIBRARIES ${MPI_Fortran_LIBRARIES})
  list(APPEND _MPI_LIBRARIES ${MPI_C_LIBRARIES})
  list(APPEND _MPI_LIBRARIES ${MPI_CXX_LIBRARIES})
  
  
  
  # Find a thread library. (Mumps needs this)
  set(CMAKE_THREAD_PREFER_PTHREAD ON)
  find_package(Threads REQUIRED)
  
  
  
  # Include Metis
  find_path(METIS_INCLUDE_DIR metis.h PATHS ${PETSC_INCLUDES} NO_DEFAULT_PATH)
  find_library(METIS_LIBRARY NAMES metis PATHS ${PETSC_LIBS_PATH} 
               NO_DEFAULT_PATH)
  # Mimic "REQUIRED" behaviour
  if(NOT (METIS_LIBRARY AND METIS_INCLUDE_DIR))
    message(FATAL_ERROR "Could not find Metis. You have to configure PETSc to "
            "include metis")
  endif()
  include_directories(${METIS_INCLUDE_DIR})
  
  
  
  #When using MPI we're also using Parmetis. Include it.
  find_path(PARMETIS_INCLUDE_DIR parmetis.h PATHS ${PETSC_INCLUDES} 
            NO_DEFAULT_PATH)
  find_library(PARMETIS_LIBRARY NAMES parmetis PATHS ${PETSC_LIBS_PATH} 
               NO_DEFAULT_PATH)
  # Mimic "REQUIRED" behaviour
  if(NOT (PARMETIS_LIBRARY AND PARMETIS_INCLUDE_DIR))
    message(FATAL_ERROR "Could not find Parmetis. You have to configure PETSc "
            "to include Parmetis")
  endif()
  include_directories(BEFORE ${PARMETIS_INCLUDE_DIR})
  
  
  
  # When using MPI we're also using the MUMPS solver. Take the supplied dmumps,
  # mumps_common and scalapack dependencies.
  find_path(MUMPS_INCLUDE_DIR mumps_compat.h PATHS ${PETSC_INCLUDES}
            NO_DEFAULT_PATH)
  # Find mumps and its dependent libraries.
  find_library(MUMPS_LIBRARY NAMES dmumps PATHS ${PETSC_LIBS_PATH}
               NO_DEFAULT_PATH)
  # mumps_common library
  find_library(MUMPS_LIBRARY_COMMON NAMES mumps_common PATHS ${PETSC_LIBS_PATH}
               NO_DEFAULT_PATH)
  # pord library (the sequential ordering tool which ships with mumps)
  find_library(PORD_LIBRARY NAMES pord PATHS ${PETSC_LIBS_PATH} NO_DEFAULT_PATH)
  # scalapack library
  find_library(SCALAPACK_LIBRARY NAMES scalapack PATHS ${PETSC_LIBS_PATH}
               NO_DEFAULT_PATH)
  set(MUMPS_LIBRARIES ${MUMPS_LIBRARY} ${MUMPS_LIBRARY_COMMON} 
      ${SCALAPACK_LIBRARY} ${PORD_LIBRARY} gfortran)
  # Mimic "REQUIRED" behaviour
  if(NOT (MUMPS_LIBRARIES AND MUMPS_INCLUDE_DIR)) 
    message(FATAL_ERROR "Mumps or one of its dependencies could not be found."
            " - check CMakeCache.txt.")
  endif()
  include_directories(${MUMPS_INCLUDE_DIR})
endif(_USING_MPI)



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



# Include TECPLOT library.
find_path(TECPLOT_INCLUDE_DIR TECIO.h PATHS
          ${PARMOON_EXTLIB_PATH}/tecplot/include NO_DEFAULT_PATH)
find_library(TECPLOT_LIBRARY NAMES tecio_${_ARCH} PATHS
             ${PARMOON_EXTLIB_PATH}/tecplot/lib NO_DEFAULT_PATH)
# Mimic "REQUIRED" behaviour
if(NOT (TECPLOT_LIBRARY AND TECPLOT_INCLUDE_DIR)) 
  message(FATAL_ERROR "Could not find Tecplot. It should be in the repository")
endif()
include_directories(${TECPLOT_INCLUDE_DIR})



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



# Include UMFPACK library.
find_path(UMFPACK_INCLUDE_DIR umfpack.h PATHS ${PETSC_INCLUDES} NO_DEFAULT_PATH)
find_library(UMFPACK_LIBRARY NAMES umfpack 
             PATHS ${PETSC_DIR}/${PETSC_ARCH}/lib/ NO_DEFAULT_PATH)
find_library(UMFPACK_SUITESE_LIBRARY NAMES suitesparseconfig
             PATHS ${PETSC_DIR}/${PETSC_ARCH}/lib/ NO_DEFAULT_PATH)
find_library(UMFPACK_AMD_LIBRARY NAMES amd
             PATHS ${PETSC_DIR}/${PETSC_ARCH}/lib/ NO_DEFAULT_PATH)
find_library(UMFPACK_CHOLMOD_LIBRARY NAMES cholmod
             PATHS ${PETSC_DIR}/${PETSC_ARCH}/lib/ NO_DEFAULT_PATH)
find_library(UMFPACK_COLAMD_LIBRARY NAMES colamd
             PATHS ${PETSC_DIR}/${PETSC_ARCH}/lib/ NO_DEFAULT_PATH)
set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARY} ${UMFPACK_SUITESE_LIBRARY}
    ${UMFPACK_AMD_LIBRARY} ${UMFPACK_CHOLMOD_LIBRARY} ${UMFPACK_COLAMD_LIBRARY})
include_directories(${UMFPACK_INCLUDE_DIR})



# Set up list of external libraries. If library A depends on library B then A
# should be listed before B. For example the Blas library is the last one 
# because it depends on no other library but others depend on it.
list(APPEND _EXTERN_LIBRARIES ${PETSC_LIBRARIES})
list(APPEND _EXTERN_LIBRARIES ${TRIANGLE_LIBRARY})
list(APPEND _EXTERN_LIBRARIES ${TECPLOT_LIBRARY})
list(APPEND _EXTERN_LIBRARIES ${TETGEN_LIBRARY})
if(_USING_OMP)
  list(APPEND _EXTERN_LIBRARIES ${PARDISO_LIBRARY})
endif(_USING_OMP)
list(APPEND _EXTERN_LIBRARIES ${UMFPACK_LIBRARIES})
if(_USING_MPI)
    list(APPEND _EXTERN_LIBRARIES ${MUMPS_LIBRARIES})
    list(APPEND _EXTERN_LIBRARIES ${PARMETIS_LIBRARY})
    list(APPEND _EXTERN_LIBRARIES ${METIS_LIBRARY})
    list(APPEND _EXTERN_LIBRARIES ${_MPI_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endif(_USING_MPI)
list(APPEND _EXTERN_LIBRARIES ${LAPACK_LIBRARIES})
list(APPEND _EXTERN_LIBRARIES ${BLAS_LIBRARIES})
