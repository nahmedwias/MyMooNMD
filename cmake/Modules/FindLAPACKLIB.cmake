# ===================================================================
# This is FindLAPACK.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 17 June 2015
# searching for a LAPACK lib in the system 
# if found, this will define
#  LAPACK_FOUND - System has LAPACK
#  LAPACK_INCLUDE_DIRS - The LAPACK include directories
#  LAPACK_LIBRARIES - The libraries needed to use LAPACK
# History:
# 2015/08/18 Clemens Bartsch: Added TODO. Added NO_DEFAULT_PATH
#	keyword. Slightly changed the output.
# 2015/08/20 Clemens Bartsch: Placed standard handling of arguments
#         outside the if(..._library) block - that was a bug.
#
#	TODO The module needs a rework.	Which libraries/paths MUST be
#	found, which are optional? What should a default search find and set?
#	Why is the include path set to equal the library path? 
# ===================================================================
if(LAPACK_INCLUDES AND LAPACK_LIBRARIES)
  set(LAPACK_FIND_QUIETLY TRUE)
endif(LAPACK_INCLUDES AND LAPACK_LIBRARIES)

# TODO CB 2015/08/18 Add code to search at user-defined places (see e.g. FindACMLBLAS.cmake)
if(FIND_USER_LAPACK)
  message("Find user Lapack is not yet implemented.")
endif(FIND_USER_LAPACK)

# Search for the needed libraries exclusively in the ParMooN EXT_LIB directory.
if(NOT LAPACK_FOUND)
  message("Searching in ParMooN EXT_LIB path.")
  find_library(SCALAPACK_LIBRARY NAMES scalapack_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS NO_DEFAULT_PATH) 
  find_library(BLACS_LIBRARY NAMES blacs_MPI_${ARCH}-0 PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS NO_DEFAULT_PATH)     
  find_library(BLACS_F77LIBRARY NAMES blacsF77init_MPI_${ARCH}-0 PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS NO_DEFAULT_PATH) 
  
  # Combine libraries into LAPACK_LIBRARY. TODO Here the combination is unclear to me.
  if(SCALAPACK_LIBRARY)
    set(LAPACK_LIBRARY  ${LAPACK_LIBRARY}  ${SCALAPACK_LIBRARY} ) 
  else(SCALAPACK_LIBRARY)   
    set(LAPACK_LIBRARY FALSE)
    message(FATAL_ERROR "SCALAPACK_LIBRARY: ${SCALAPACK_LIBRARY} not found !") 
  endif(SCALAPACK_LIBRARY)
    
  if(BLACS_LIBRARY)
    set(LAPACK_LIBRARY  ${LAPACK_LIBRARY}  ${BLACS_LIBRARY} ) 
  else(BLACS_LIBRARY)   
    set(LAPACK_LIBRARY FALSE)
  endif(BLACS_LIBRARY)    
    
  if(BLACS_F77LIBRARY)
    set(LAPACK_LIBRARY  ${LAPACK_LIBRARY}  ${BLACS_F77LIBRARY} ) 
  else(BLACS_F77LIBRARY)   
    set(LAPACK_LIBRARY FALSE)
  endif(BLACS_F77LIBRARY)   
  
  # This is a hack - include directory is set to the library path.
  # To use LAPACK we don't need its include files (#extern blocks
  # in ParMooN code...), but we want it to be set nevertheless,
  # for standard handling.
  #set(LAPACK_INCLUDE_DIR ${LAPACK_LIBRARY}) 
  
  # Handling of standard arguments.
  include(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set LAPACK_FOUND to TRUE
  # if all listed variables are TRUE
  #find_package_handle_standard_args(LAPACK  DEFAULT_MSG
  #                                  LAPACK_LIBRARY LAPACK_INCLUDE_DIR)
  find_package_handle_standard_args(LAPACK  DEFAULT_MSG LAPACK_LIBRARY)

  #mark_as_advanced(LAPACK_INCLUDE_DIR LAPACK_LIBRARY)
  mark_as_advanced(LAPACK_LIBRARY)
  
  # Set two non-cache variables to be use when including and linking against LAPACK.
  set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})
  #set(LAPACK_INCLUDE_DIRS ${LAPACK_INCLUDE_DIR})
  
endif(NOT LAPACK_FOUND)



