# ===================================================================
# This is FindTETGEN.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a TETGEN lib in the system 
# if found, this will define
#  TETGEN_FOUND - System has TETGEN
#  TETGEN_INCLUDE_DIRS - The TETGEN include directories
#  TETGEN_LIBRARIES - The libraries needed to use TETGEN
# 2015/08/18 Clemens Bartsch: Introduced variable FIND_USER_TETGEN
#	for activating non-ParMooN search. Added NO_DEFAULT_PATH
#	keyword. Slightly changed the output. Added comments.
# ===================================================================
if(TETGEN_INCLUDES AND TETGEN_LIBRARIES)
  set(TETGEN_FIND_QUIETLY TRUE)
endif(TETGEN_INCLUDES AND TETGEN_LIBRARIES)

if(NOT TETGEN_FOUND)
  # Search for the library in standard non-ParMooN paths and in those specified after PATHS.
  if(FIND_USER_TETGEN)
    message("Searching in default and user paths.")
    find_path(TETGEN_INCLUDE_DIR   tetgen.h PATHS $ENV{TETGENDIR}/include ${CMAKE_INCLUDE_PATH})
    find_library(TETGEN_LIBRARY NAMES tet PATHS $ENV{TETGENDIR}/lib ${CMAKE_LIBRARY_PATH})
    get_filename_component(TETGEN_LIBDIR ${TETGEN_LIBRARY} PATH)
   endif(FIND_USER_TETGEN)   
 
   # Search for the library exclusively in the ParMooN EXT_LIB path.
  if(NOT TETGEN_LIBRARY)
    message("Searching in ParMooN EXT_LIB path. Selected architecture ARCH=${ARCH}")
    find_path(TETGEN_INCLUDE_DIR  tetgen.h PATHS ${PARMOON_EXTLIB_PATH}/tetgen)
    find_library(TETGEN_LIBRARY NAMES tet_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/tetgen)
  endif(NOT TETGEN_LIBRARY)
  
  if(TETGEN_LIBRARY)      
    # Set TETGEN library and directory variables.
    include(FindPackageHandleStandardArgs)
    
    set(TETGEN_LIBRARIES ${TETGEN_LIBRARY})
    set(TETGEN_INCLUDE_DIRS ${TETGEN_INCLUDE_DIR})

    # handle the QUIETLY and REQUIRED arguments and set TETGEN_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args(TETGEN  DEFAULT_MSG
                                      TETGEN_LIBRARY TETGEN_INCLUDE_DIR)

    mark_as_advanced(TETGEN_INCLUDE_DIR TETGEN_LIBRARY)
  endif(TETGEN_LIBRARY)

endif()


