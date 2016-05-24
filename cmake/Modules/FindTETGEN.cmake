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
# 2015/08/20 Clemens Bartsch: Placed standard handling of arguments
#         outside the if(..._library) block - that was a bug.
# ===================================================================
if(TETGEN_INCLUDES AND TETGEN_LIBRARIES)
  set(TETGEN_FIND_QUIETLY TRUE)
endif(TETGEN_INCLUDES AND TETGEN_LIBRARIES)

if(NOT TETGEN_FOUND)
  # Search for the library in standard non-ParMooN paths and in those specified after PATHS.
  if(FIND_USER_TETGEN)
    message("Searching in default and user paths.")
    find_path(TETGEN_INCLUDE_DIR   tetgen.h PATHS $ENV{TETGENDIR}/include)
    find_library(TETGEN_LIBRARY NAMES tet PATHS $ENV{TETGENDIR}/lib)
    get_filename_component(_TETGEN_LIBDIR ${TETGEN_LIBRARY} PATH)
   endif(FIND_USER_TETGEN)   
 
   # Search for the library exclusively in the ParMooN EXT_LIB path.
  if(NOT TETGEN_LIBRARY)
    message("Searching in ParMooN EXT_LIB path. Selected architecture _ARCH=${_ARCH}")
    find_path(TETGEN_INCLUDE_DIR  tetgen.h PATHS ${PARMOON_EXTLIB_PATH}/tetgen)
    find_library(TETGEN_LIBRARY NAMES tet_${_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/tetgen)
  endif(NOT TETGEN_LIBRARY)
   
  # Handling of standard arguments.
  include(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set TETGEN_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(TETGEN  DEFAULT_MSG
                                    TETGEN_LIBRARY TETGEN_INCLUDE_DIR)
  mark_as_advanced(TETGEN_INCLUDE_DIR TETGEN_LIBRARY)
    
  # Set two non-cache variables to be use when including and linking against TETGEN.
  set(TETGEN_LIBRARIES ${TETGEN_LIBRARY})
  message("TETGEN_LIBRARIES" ${TETGEN_LIBRARIES})
  set(TETGEN_INCLUDE_DIRS ${TETGEN_INCLUDE_DIR})

endif()


