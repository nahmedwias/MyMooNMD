# ===================================================================
# This is FindTRIANGLE.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a TRIANGLE lib in the system 
# if found, this will define
#  TRIANGLE_FOUND - System has TRIANGLE
#  TRIANGLE_INCLUDE_DIRS - The TRIANGLE include directories
#  TRIANGLE_LIBRARIES - The libraries needed to use TRIANGLE
# 2015/08/18 Clemens Bartsch: Introduced variable FIND_USER_TRIANGLE
#	for activating non-ParMooN search. Added NO_DEFAULT_PATH
#	keyword. Slightly changed the output. Added comments.
# 2015/08/20 Clemens Bartsch: Placed standard handling of arguments
#         outside the if(..._library) block - that was a bug.
# ===================================================================
if(TRIANGLE_INCLUDES AND TRIANGLE_LIBRARIES)
  set(TRIANGLE_FIND_QUIETLY TRUE)
endif(TRIANGLE_INCLUDES AND TRIANGLE_LIBRARIES)

if(NOT TRIANGLE_FOUND)

  # Search for the library in standard non-ParMooN paths and in those specified after PATHS.
  if(FIND_USER_TRIANGLE)
    message("Searching in default and user paths.")
    find_path(TRIANGLE_INCLUDE_DIR   triangle.h PATHS $ENV{TRIANGLEDIR}/include)
    find_library(TRIANGLE_LIBRARY NAMES triangle_${_ARCH} PATHS $ENV{TRIANGLEDIR}/lib)
    get_filename_component(_TRIANGLE_LIBDIR ${TRIANGLE_LIBRARY} PATH)
  endif(FIND_USER_TRIANGLE)
  
  # Search for the library exclusively in the ParMooN EXT_LIB path.
  if(NOT TRIANGLE_LIBRARY)
    message("Searching in ParMooN EXT_LIB path. Selected architecture _ARCH=${_ARCH}")
    find_path(TRIANGLE_INCLUDE_DIR  triangle.h PATHS ${PARMOON_EXTLIB_PATH}/Triangle NO_DEFAULT_PATH)
    find_library(TRIANGLE_LIBRARY NAMES triangle_${_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Triangle NO_DEFAULT_PATH)
    message("Found library: " ${TRIANGLE_LIBRARY})
  endif(NOT TRIANGLE_LIBRARY)

  # Handling of standard arguments.
  include(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set TRIANGLE_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(TRIANGLE  DEFAULT_MSG
                                      TRIANGLE_LIBRARY TRIANGLE_INCLUDE_DIR)
  mark_as_advanced(TRIANGLE_INCLUDE_DIR TRIANGLE_LIBRARY)

  # Set two non-cache variables to be use when including and linking against BLAS.
  set(TRIANGLE_LIBRARIES ${TRIANGLE_LIBRARY})
  set(TRIANGLE_INCLUDE_DIRS ${TRIANGLE_INCLUDE_DIR}) 

endif(NOT TRIANGLE_FOUND)


