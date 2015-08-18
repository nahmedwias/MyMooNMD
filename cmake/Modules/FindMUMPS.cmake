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
# ===================================================================
if(MUMPS_INCLUDES AND MUMPS_LIBRARIES)
  set(MUMPS_FIND_QUIETLY TRUE)
endif(MUMPS_INCLUDES AND MUMPS_LIBRARIES)

if(NOT MUMPS_FOUND)
	
  # Search for the Library in standard non-ParMooN paths and in those specified after PATHS.
  if(FIND_USER_MUMPS) 
    message("Searching in default and user paths.")
    find_path(MUMPS_INCLUDE_DIR  mumps_compat.h PATHS $ENV{MUMPSDIR}/include ${CMAKE_INCLUDE_PATH})
    find_library(MUMPS_LIBRARY NAMES dmumps PATHS $ENV{MUMPSDIR}/lib ${CMAKE_LIBRARY_PATH})
    get_filename_component(MUMPS_LIBDIR ${MUMPS_LIBRARY} PATH)
    find_library(MUMPS_LIBRARY_COMMON NAMES mumps_common PATHS ${MUMPS_LIBDIR})
   endif(FIND_USER_MUMPS)
  
  # Search for the library exclusively in the ParMooN EXT_LIB path.
  if(NOT MUMPS_LIBRARY)
    message("Searching in ParMooN EXT_LIB path. Selected architecture ARCH=${ARCH}")
    find_path(MUMPS_INCLUDE_DIR  mumps_compat.h PATHS ${PARMOON_EXTLIB_PATH}/MUMPS/include NO_DEFAULT_PATH)
    find_library(MUMPS_LIBRARY NAMES dmumps_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/MUMPS/lib NO_DEFAULT_PATH) 
    get_filename_component(MUMPS_LIBDIR ${MUMPS_LIBRARY} PATH)
    find_library(MUMPS_LIBRARY_COMMON NAMES mumps_common_${ARCH} PATHS ${MUMPS_LIBDIR} NO_DEFAULT_PATH)
  endif(NOT MUMPS_LIBRARY)
  
  if(MUMPS_LIBRARY)  

    # Find dependencies for MUMPS.
    
    # First search in default and user paths.
    find_library(SCALAPACK_LIB NAMES scalapack PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})
    find_library(PARMETIS_LIB NAMES parmetis PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})
    find_library(METIS_LIB NAMES metis PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})
    
    # If not found there, take the libraries supplied by ParMooN.
    if(NOT SCALAPACK_LIB)
      find_library(SCALAPACK_LIB NAMES scalapack_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS NO_DEFAULT_PATH)
    endif()
    if(NOT PARMETIS_LIB)
      find_library(PARMETIS_LIB NAMES parmetis_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Metis NO_DEFAULT_PATH)
    endif()
    if(NOT METIS_LIB)
      find_library(METIS_LIB NAMES metis_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Metis NO_DEFAULT_PATH)   
    endif()      
    
    # Combine MUMPS and its dependencies.  
    if(MUMPS_LIBRARY_COMMON AND SCALAPACK_LIB AND METIS_LIB)
      set(MUMPS_LIBRARY  ${MUMPS_LIBRARY} ${MUMPS_LIBRARY_COMMON} ${SCALAPACK_LIB} ${PARMETIS_LIB} ${METIS_LIB})
    else(MUMPS_LIBRARY_COMMON AND SCALAPACK_LIB AND METIS_LIB)
      #This case should not happen, because the libraries ship with ParMooN.  
      if(NOT SCALAPACK_LIB)
        message("MUMPS found in system but scalapack not found in the search path")
      endif()
      if(NOT PARMETIS_LIB)
        message("MUMPS found in system but parmetis not found in the search path")
      endif()
      if(NOT METIS_LIB)
        message("MUMPS found in system but metis not found in the search path")
      endif()      
      set(MUMPS_LIBRARY FALSE)
    endif(MUMPS_LIBRARY_COMMON AND SCALAPACK_LIB AND METIS_LIB)
    
    # set MUMPS
    if(MUMPS_LIBRARY)
      include(FindPackageHandleStandardArgs)
    
      set(MUMPS_LIBRARIES ${MUMPS_LIBRARY})
      set(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})

      # handle the QUIETLY and REQUIRED arguments and set MUMPS_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(MUMPS  DEFAULT_MSG
                                        MUMPS_LIBRARY MUMPS_INCLUDE_DIR)

      mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_LIBRARY)
    endif(MUMPS_LIBRARY)  
  endif(MUMPS_LIBRARY)

endif(NOT MUMPS_FOUND)


