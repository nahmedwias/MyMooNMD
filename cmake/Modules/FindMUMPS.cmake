# ===================================================================
# This is FindMUMPS.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a MUMPS lib in the system 
# if found, this will define
#  MUMPS_FOUND - System has MUMPS
#  MUMPS_INCLUDE_DIRS - The MUMPS include directories
#  MUMPS_LIBRARIES - The libraries needed to use MUMPS
# History:
# 2015/08/18 Clemens Bartsch: Introduced variable FIND_USER_MUMPS 
#	for activating non-ParMooN search. Added NO_DEFAULT_PATH
#	keyword. Slightly changed the output. Added comments.
# 2015/08/20 Clemens Bartsch: Placed standard handling of arguments
#         outside the if(..._library) block - that was a bug.
# 2015/08/24 Clemens Bartsch: Gave the package a complete rework,
#            because not resolving MUMPS dependencies is one of the main
#            reasons for MPI compiling not working.
#
#            TODO The esmumps and scotcherr libraries are not supplied for
#            all ${_ARCH}-types.
#
#            TODO Handling of the variable MUMPS_LIBRARY and MUMPS_LIBRARIES
#            - which should appear in the cache, and in which form?
#
# ===================================================================
if (MUMPS_LIBRARIES AND MUMPS_INCLUDE_DIRS)
    set(MUMPS_FIND_QUIETLY TRUE)
endif(MUMPS_LIBRARIES AND MUMPS_INCLUDE_DIRS)

if (NOT MUMPS_FOUND)
  #Make sure the dependent libraries are there.
  #LAPACK Note: The default FindLAPACK also muddles up the BLAS!
  if(NOT LAPACK_LIBRARIES) #Lapack used without header files...
      find_package(LAPACK REQUIRED)
  endif()
  #BLAS
  if(NOT (BLAS_LIBRARIES AND BLAS_INCLUDE_DIRS))
      find_package(BLAS REQUIRED)
  endif()
  #Parmetis
  if(NOT (PARMETIS_LIBRARIES AND PARMETIS_INCLUDE_DIRS))
      find_package(PARMETIS REQUIRED)
  endif()
      
  # Search for the Library in standard non-ParMooN paths and in those specified after PATHS.
  if(FIND_USER_MUMPS) 
    message("Find user MUMPS not implemented yet.")
  endif(FIND_USER_MUMPS)
  
  # Search for the library exclusively in the ParMooN EXT_LIB path.
    message("Searching in ParMooN EXT_LIB path. Selected architecture _ARCH=${_ARCH}")
    find_path(MUMPS_INCLUDE_DIR  mumps_compat.h PATHS ${PARMOON_EXTLIB_PATH}/MUMPS/include NO_DEFAULT_PATH)
    find_library(MUMPS_LIBRARY NAMES dmumps_${_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/MUMPS/lib NO_DEFAULT_PATH)
    get_filename_component(_MUMPS_LIBDIR ${MUMPS_LIBRARY} PATH)
    # mumps_common library
    find_library(MUMPS_LIBRARY_COMMON NAMES mumps_common_${_ARCH} PATHS ${_MUMPS_LIBDIR} NO_DEFAULT_PATH)
    # esmumps/dmumps/tesmumps
    find_library(ESMUMPS_LIBRARY NAMES esmumps_TYRONE64 PATHS ${_MUMPS_LIBDIR} NO_DEFAULT_PATH)
    find_library(DMUMPS_LIBRARY NAMES dmumps_${_ARCH} PATHS ${_MUMPS_LIBDIR} NO_DEFAULT_PATH)
    find_library(TESMUMPS_LIBRARY NAMES ptesmumps_${_ARCH} PATHS ${_MUMPS_LIBDIR} NO_DEFAULT_PATH)
    # scotch and pord libraries
    find_library(PORD_LIBRARY NAMES pord_${_ARCH} PATHS ${_MUMPS_LIBDIR} NO_DEFAULT_PATH)
    find_library(SCOTCH_LIBRARY NAMES scotch_${_ARCH} PATHS ${_MUMPS_LIBDIR} NO_DEFAULT_PATH)
    find_library(SCOTCHERR_LIBRARY NAMES scotcherr_TYRONE64 PATHS ${_MUMPS_LIBDIR} NO_DEFAULT_PATH)
    
    set(MUMPS_LIBRARY ${MUMPS_LIBRARY} ${ESMUMPS_LIBRARY} ${DMUMPS_LIBRARY} ${TESMUMPS_LIBRARY}
                         ${MUMPS_LIBRARY_COMMON} ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY} ${PORD_LIBRARY}
                         ${ESMUMPS_LIBRARY} ${DMUMPS_LIBRARY} ${TESMUMPS_LIBRARY})
     
  # Combine MUMPS and its dependencies
  set(MUMPS_LIBRARY ${MUMPS_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
       ${PARMETIS_LIBRARIES} ${METIS_LIBRARIES})
  
  # Handling of standard arguments.
  include(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set MUMPS_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(MUMPS  DEFAULT_MSG
                                    MUMPS_LIBRARY MUMPS_INCLUDE_DIR)

  mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_LIBRARY)
      
  # Set two non-cache variables to be use when including and linking against MUMPS.          
  set(MUMPS_LIBRARIES ${MUMPS_LIBRARY})
  set(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})

endif (NOT MUMPS_FOUND)

