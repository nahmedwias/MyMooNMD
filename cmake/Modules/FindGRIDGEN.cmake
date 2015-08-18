# ===================================================================
# This is FindGRIDGEN.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a GRIDGEN lib in the system 
# if found, this will define
#  GRIDGEN_FOUND - System has GRIDGEN
#  GRIDGEN_INCLUDE_DIRS - The GRIDGEN include directories
#  GRIDGEN_LIBRARIES - The libraries needed to use GRIDGEN
# 2015/08/18 Clemens Bartsch: Introduced variable FIND_USER_GRIDGEN
#	for activating non-ParMooN search. Added NO_DEFAULT_PATH
#	keyword. Slightly changed the output. Added comments.
# ===================================================================
if(GRIDGEN_INCLUDES AND GRIDGEN_LIBRARIES)
  set(GRIDGEN_FIND_QUIETLY TRUE)
endif(GRIDGEN_INCLUDES AND GRIDGEN_LIBRARIES)

if(NOT GRIDGEN_FOUND)

  # Search for the library in standard non-ParMooN paths and in those specified after PATHS.
  if(FIND_USER_GRIDGEN)
    message("Searching in default and user paths.")
    find_path(GRIDGEN_INCLUDE_DIR   gridgen.h PATHS $ENV{GRIDGENDIR}/include ${CMAKE_INCLUDE_PATH})
    find_library(GRIDGEN_LIBRARY NAMES gridgen PATHS $ENV{GRIDGENDIR}/lib ${CMAKE_LIBRARY_PATH})
    get_filename_component(GRIDGEN_LIBDIR ${GRIDGEN_LIBRARY} PATH)
  endif(FIND_USER_GRIDGEN)
  
  # Search for the library exclusively in the ParMooN EXT_LIB path.
  if(NOT GRIDGEN_LIBRARY)
    message("Searching in ParMooN EXT_LIB path. Selected architecture ARCH=${ARCH}")
    find_path(GRIDGEN_INCLUDE_DIR  gridgen.h PATHS ${PARMOON_EXTLIB_PATH}/GridGen NO_DEFAULT_PATH)
    find_library(GRIDGEN_LIBRARY NAMES gridgen_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/GridGen NO_DEFAULT_PATH)
  endif(NOT GRIDGEN_LIBRARY)
  
  if(GRIDGEN_LIBRARY)      
    # Set GRIDGEN library and directory variables
    include(FindPackageHandleStandardArgs)
    
    set(GRIDGEN_LIBRARIES ${GRIDGEN_LIBRARY})
    set(GRIDGEN_INCLUDE_DIRS ${GRIDGEN_INCLUDE_DIR})

    # handle the QUIETLY and REQUIRED arguments and set GRIDGEN_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args(GRIDGEN  DEFAULT_MSG
                                      GRIDGEN_LIBRARY GRIDGEN_INCLUDE_DIR)

    mark_as_advanced(GRIDGEN_INCLUDE_DIR GRIDGEN_LIBRARY)
  endif(GRIDGEN_LIBRARY)

endif(NOT GRIDGEN_FOUND)


