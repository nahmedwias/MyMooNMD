# Set path to ParMooN package search modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

# mandatory packages
set(PARMOON_WITH_BLAS false)
set(PARMOON_WITH_LAPACK false)
# optional packages
set(PARMOON_WITH_UMFPACK false)
set(PARMOON_WITH_PARDISO false)
set(PARMOON_WITH_MPI false)
set(PARMOON_WITH_METIS false)
set(PARMOON_WITH_PARMETIS false)
set(PARMOON_WITH_MUMPS false)
set(PARMOON_WITH_SCALAPACK false)
set(PARMOON_WITH_PETSC false)



find_package(PETSc)
if(${PETSC_FOUND})
  list(APPEND PARMOON_EXTERNAL_LIBRARIES ${PETSC_LIBRARIES})
  include_directories(${PETSC_INCLUDES})
  set(PARMOON_WITH_PETSC true)
  
  # check for more libraries in petsc:
  if("${PETSC_LIBRARIES}" MATCHES blas)
    set(PARMOON_WITH_BLAS true)
  endif("${PETSC_LIBRARIES}" MATCHES blas)
  if("${PETSC_LIBRARIES}" MATCHES lapack)
    set(PARMOON_WITH_LAPACK true)
  endif("${PETSC_LIBRARIES}" MATCHES lapack)
  if("${PETSC_LIBRARIES}" MATCHES umfpack
     AND "${PETSC_LIBRARIES}" MATCHES suitesparseconfig
     AND "${PETSC_LIBRARIES}" MATCHES amd
     AND "${PETSC_LIBRARIES}" MATCHES cholmod
     AND "${PETSC_LIBRARIES}" MATCHES colamd)
    set(PARMOON_WITH_UMFPACK true)
  endif()
  if("${PETSC_LIBRARIES}" MATCHES mpi)
    set(PARMOON_WITH_MPI true)
    set(MPIEXEC ${PETSC_MPIEXEC})
  endif("${PETSC_LIBRARIES}" MATCHES mpi)
  if("${PETSC_LIBRARIES}" MATCHES metis)
    set(PARMOON_WITH_METIS true)
  endif("${PETSC_LIBRARIES}" MATCHES metis)
  if("${PETSC_LIBRARIES}" MATCHES dmumps 
     AND "${PETSC_LIBRARIES}" MATCHES mumps_common
     AND "${PETSC_LIBRARIES}" MATCHES pord)
    set(PARMOON_WITH_MUMPS true)
  endif()
  if("${PETSC_LIBRARIES}" MATCHES scalapack)
    set(PARMOON_WITH_SCALAPACK true)
  endif("${PETSC_LIBRARIES}" MATCHES scalapack)
endif(${PETSC_FOUND})

message(STATUS "done with FindPETSc, "
               "blas: ${PARMOON_WITH_BLAS}, lapack: ${PARMOON_WITH_LAPACK}, "
               "umfpack: ${PARMOON_WITH_UMFPACK}, mpi: ${PARMOON_WITH_MPI}, "
               "metis: ${PARMOON_WITH_METIS}, mumps: ${PARMOON_WITH_MUMPS}, "
               "scalapack: ${PARMOON_WITH_SCALAPACK}")

if(NOT PARMOON_WITH_BLAS)
  find_package(BLAS REQUIRED)
  list(APPEND PARMOON_EXTERNAL_LIBRARIES ${BLAS_LIBRARIES})
  set(PARMOON_WITH_BLAS true)
  message(STATUS "blas found  ${BLAS_LIBRARIES}" )
endif(NOT PARMOON_WITH_BLAS)

if(NOT PARMOON_WITH_LAPACK)
  find_package(LAPACK REQUIRED)
  list(APPEND PARMOON_EXTERNAL_LIBRARIES ${LAPACK_LIBRARIES})
  set(PARMOON_WITH_LAPACK true)
  message(STATUS "lapack found  ${LAPACK_LIBRARIES}" )
endif(NOT PARMOON_WITH_LAPACK)

if(NOT PARMOON_WITH_MPI)
  find_package(MPI)
  list(APPEND PARMOON_EXTERNAL_LIBRARIES ${MPI_CXX_LIBRARIES})
  include_directories(${MPI_CXX_INCLUDE_PATH})
  set(PARMOON_WITH_MPI true)
  message(STATUS "blas found  ${MPI_CXX_INCLUDE_PATH}  ${MPI_CXX_LIBRARIES}" )
endif(NOT PARMOON_WITH_MPI)

if(NOT PARMOON_WITH_UMFPACK)
  find_package(UMFPACK)
  if(${UMFPACK_FOUND})
    list(APPEND PARMOON_EXTERNAL_LIBRARIES ${UMFPACK_LIBRARIES})
    include_directories(${UMFPACK_INCLUDES})
    set(PARMOON_WITH_UMFPACK true)
    message(STATUS "Umfpack found  ${UMFPACK_INCLUDES}  ${UMFPACK_LIBRARIES}" )
  endif(${UMFPACK_FOUND})
endif(NOT PARMOON_WITH_UMFPACK)

if(NOT PARMOON_WITH_METIS)
  find_package(METIS)
  if(${METIS_FOUND})
    list(APPEND PARMOON_EXTERNAL_LIBRARIES ${METIS_LIBRARIES})
    include_directories(${METIS_INCLUDE_DIRS})
    set(PARMOON_WITH_METIS true)
    message(STATUS "Metis found  ${METIS_INCLUDE_DIRS}  ${METIS_LIBRARIES}" )
  endif(${METIS_FOUND})
endif(NOT PARMOON_WITH_METIS)

if(NOT PARMOON_WITH_MUMPS)
  find_package(MUMPS)
  if(${MUMPS_FOUND})
    list(APPEND PARMOON_EXTERNAL_LIBRARIES ${MUMPS_LIBRARIES})
    include_directories(${MUMPS_INCLUDE_DIRS})
    set(PARMOON_WITH_MUMPS true)
    message(STATUS "Mumps found  ${MUMPS_INCLUDE_DIRS}  ${MUMPS_LIBRARIES}" )
  endif(${MUMPS_FOUND})
endif(NOT PARMOON_WITH_MUMPS)

if(NOT PARMOON_WITH_SCALAPACK)
  find_package(SCALAPACK)
  if(${SCALAPACK_FOUND})
    list(APPEND PARMOON_EXTERNAL_LIBRARIES ${SCALAPACK_LIBRARIES})
    include_directories(${SCALAPACK_INCLUDES})
    set(PARMOON_WITH_SCALAPACK true)
  endif(${SCALAPACK_FOUND})
endif(NOT PARMOON_WITH_SCALAPACK)

# sanity checks:
if(PARMOON_USING_MPI AND NOT PARMOON_WITH_MPI)
  message(FATAL_ERROR
          "ParMooN needs an mpi library in mpi mode but could not find one.")
endif(PARMOON_USING_MPI AND NOT PARMOON_WITH_MPI)
