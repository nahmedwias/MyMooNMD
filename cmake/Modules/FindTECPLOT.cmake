# ===================================================================
# This is FindTECPLOT.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a TECPLOT lib in the system 
# if found, this will define
#  TECPLOT_FOUND - System has TECPLOT
#  TECPLOT_INCLUDE_DIRS - The TECPLOT include directories
#  TECPLOT_LIBRARIES - The libraries needed to use TECPLOT
# 2015/08/18 Clemens Bartsch: Introduced variable FIND_USER_TECPLOT
#	for activating non-ParMooN search. Added NO_DEFAULT_PATH
#	keyword. Slightly changed the output. Added comments.
# ===================================================================
if(TECPLOT_INCLUDES AND TECPLOT_LIBRARIES)
  set(TECPLOT_FIND_QUIETLY TRUE)
endif(TECPLOT_INCLUDES AND TECPLOT_LIBRARIES)

if(NOT TECPLOT_FOUND)

  # Search for the library in standard non-ParMooN paths and in those specified after PATHS.
  if(FIND_USER_TECPLOT)
    message("Searching in default and user paths.")
    find_path(TECPLOT_INCLUDE_DIR   TECIO.h PATHS $ENV{TECPLOTDIR}/include)
    find_library(TECPLOT_LIBRARY NAMES tecio PATHS $ENV{TECPLOTDIR}/lib)
    get_filename_component(_TECPLOT_LIBDIR ${TECPLOT_LIBRARY} PATH)
  endif(FIND_USER_TECPLOT)
     
  # Search for the library exclusively in the ParMooN EXT_LIB path.   
  if(NOT TECPLOT_LIBRARY)
    message("Searching in ParMooN EXT_LIB path. Selected architecture _ARCH=${_ARCH}")
    find_path(TECPLOT_INCLUDE_DIR  TECIO.h PATHS ${PARMOON_EXTLIB_PATH}/tecplot/include NO_DEFAULT_PATH)
    find_library(TECPLOT_LIBRARY NAMES tecio_${_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/tecplot/lib NO_DEFAULT_PATH)
  endif(NOT TECPLOT_LIBRARY)
   
  # Handling of standard arguments.
  include(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set TECPLOT_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(TECPLOT  DEFAULT_MSG
                                    TECPLOT_LIBRARY TECPLOT_INCLUDE_DIR)
  mark_as_advanced(TECPLOT_INCLUDE_DIR TECPLOT_LIBRARY)

  # Set two non-cache variables to be use when including and linking against TECPLOT.
  set(TECPLOT_LIBRARIES ${TECPLOT_LIBRARY})
  set(TECPLOT_INCLUDE_DIRS ${TECPLOT_INCLUDE_DIR})

endif(NOT TECPLOT_FOUND)


