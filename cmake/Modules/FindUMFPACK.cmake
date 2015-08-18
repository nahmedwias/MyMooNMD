# ===================================================================
# This is FindUMFPACK.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a UMFPACK lib in the system 
# if found, this will define
#  UMFPACK_FOUND - System has UMFPACK
#  UMFPACK_INCLUDE_DIRS - The UMFPACK include directories
#  UMFPACK_LIBRARIES - The libraries needed to use UMFPACK
# History:
# 2015/08/18 Clemens Bartsch: Replaced variable USE_SYSTEM_UMFPACK by
#	FIND_USER_UMFPACK. Added NO_DEFAULT_PATH and REQUIRED option.
#   Slightly altered output.
# ===================================================================
if(UMFPACK_INCLUDES AND UMFPACK_LIBRARIES)
  set(UMFPACK_FIND_QUIETLY TRUE)
endif(UMFPACK_INCLUDES AND UMFPACK_LIBRARIES)

if(NOT UMFPACK_FOUND)

  # Search for the Library in standard non-ParMooN paths and in those specified after PATHS.
  if(FIND_USER_UMFPACK)
    message("Searching in default and user paths.")
    find_path(UMFPACK_INCLUDE_DIR   umfpack.h PATHS $ENV{UMFPACKDIR}/include ${CMAKE_INCLUDE_PATH})
    find_library(UMFPACK_LIBRARY NAMES umfpack PATHS $ENV{UMFPACKDIR}/lib ${CMAKE_LIBRARY_PATH})
    get_filename_component(UMFPACK_LIBDIR ${UMFPACK_LIBRARY} PATH)
    find_library(UMFPACK_LIBRARY_SUITESE NAMES suitesparseconfig PATHS ${UMFPACK_LIBDIR})
    find_library(UMFPACK_LIBRARY_AMD NAMES amd PATHS ${UMFPACK_LIBDIR})
  endif(FIND_USER_UMFPACK)
  
  # Search for the library exclusively in the ParMooN EXT_LIB path.     
  if(NOT UMFPACK_LIBRARY)
    message("Searching in ParMooN EXT_LIB path.")
    find_path(UMFPACK_INCLUDE_DIR  umfpack.h PATHS ${PARMOON_EXTLIB_PATH}/UMFPACK/Include NO_DEFAULT_PATH)
    find_library(UMFPACK_LIBRARY NAMES umfpack_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/UMFPACK/Lib NO_DEFAULT_PATH) 
    get_filename_component(UMFPACK_LIBDIR ${UMFPACK_LIBRARY} PATH)
    find_library(UMFPACK_LIBRARY_SUITESE NAMES suitesparseconfig_${ARCH} PATHS ${UMFPACK_LIBDIR} NO_DEFAULT_PATH)
    find_library(UMFPACK_LIBRARY_AMD NAMES amd_${ARCH} PATHS ${UMFPACK_LIBDIR} NO_DEFAULT_PATH)
  endif(NOT UMFPACK_LIBRARY)
  
  if(UMFPACK_LIBRARY)  
    # Find libraries mumps depends on. Start with default search.
    find_library(SCALAPACK_LIB NAMES scalapack PATHS ${UMFPACK_LIBDIR} ${CMAKE_LIBRARY_PATH})
    find_library(PARMETIS_LIB NAMES parmetis PATHS ${UMFPACK_LIBDIR} ${CMAKE_LIBRARY_PATH} ${PARMOON_EXTLIB_PATH}/Metis)
    find_library(METIS_LIB NAMES metis PATHS ${UMFPACK_LIBDIR} ${CMAKE_LIBRARY_PATH}) 
    # Search for the dependencies in ParMooN EXT_LIB exclusively.   
    if(NOT SCALAPACK_LIB)
      find_library(SCALAPACK_LIB NAMES scalapack_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS REQUIRED NO_DEFAULT_PATH)
    endif()
    if(NOT PARMETIS_LIB)
      find_library(PARMETIS_LIB NAMES parmetis_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Metis REQUIRED NO_DEFAULT_PATH)
    endif()
    if(NOT METIS_LIB)
      find_library(METIS_LIB NAMES metis_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Metis REQUIRED NO_DEFAULT_PATH)   
    endif()      
    
    # combine mumps and its deps    
    if(UMFPACK_LIBRARY_SUITESE AND UMFPACK_LIBRARY_AMD)
      set(UMFPACK_LIBRARY  ${UMFPACK_LIBRARY} ${UMFPACK_LIBRARY_AMD} ${UMFPACK_LIBRARY_SUITESE} ) 
    else(UMFPACK_LIBRARY_SUITESE AND UMFPACK_LIBRARY_AMD)   
      set(UMFPACK_LIBRARY FALSE)
    endif(UMFPACK_LIBRARY_SUITESE AND UMFPACK_LIBRARY_AMD)
    
    # set UMFPACK
    if(UMFPACK_LIBRARY)
      include(FindPackageHandleStandardArgs)
    
      set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARY})
      set(UMFPACK_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIR})

      # handle the QUIETLY and REQUIRED arguments and set UMFPACK_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(UMFPACK  DEFAULT_MSG
                                        UMFPACK_LIBRARY UMFPACK_INCLUDE_DIR)

      mark_as_advanced(UMFPACK_INCLUDE_DIR UMFPACK_LIBRARY)
    endif(UMFPACK_LIBRARY)  
  endif()

endif()


