# ===================================================================
# This is FindTRIANGLE.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a TRIANGLE lib in the system 
# if found, this will define
#  TRIANGLE_FOUND - System has TRIANGLE
#  TRIANGLE_INCLUDE_DIRS - The TRIANGLE include directories
#  TRIANGLE_LIBRARIES - The libraries needed to use TRIANGLE
# ===================================================================
if(TRIANGLE_INCLUDES AND TRIANGLE_LIBRARIES)
  set(TRIANGLE_FIND_QUIETLY TRUE)
endif(TRIANGLE_INCLUDES AND TRIANGLE_LIBRARIES)

if(NOT TRIANGLE_FOUND)
  find_path(TRIANGLE_INCLUDE_DIR   triangle.h PATHS $ENV{TRIANGLEDIR}/include ${CMAKE_INCLUDE_PATH})
  find_library(TRIANGLE_LIBRARY NAMES triangle PATHS $ENV{TRIANGLEDIR}/lib ${CMAKE_LIBRARY_PATH})
  get_filename_component(TRIANGLE_LIBDIR ${TRIANGLE_LIBRARY} PATH)
     
  if(NOT TRIANGLE_LIBRARY)
    message("TRIANGLE not found in the system, so checking the availability in ParMooN for the selected ARCH=${ARCH}")
    find_path(TRIANGLE_INCLUDE_DIR  triangle.h PATHS ${PARMOON_EXTLIB_PATH}/GridGen)
    find_library(TRIANGLE_LIBRARY NAMES gridgen_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/GridGen)
  endif(NOT TRIANGLE_LIBRARY)
  
  if(TRIANGLE_LIBRARY)      
    # set TRIANGLE
    if(TRIANGLE_LIBRARY)
      include(FindPackageHandleStandardArgs)
    
      set(TRIANGLE_LIBRARIES ${TRIANGLE_LIBRARY})
      set(TRIANGLE_INCLUDE_DIRS ${TRIANGLE_INCLUDE_DIR})

      # handle the QUIETLY and REQUIRED arguments and set TRIANGLE_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(TRIANGLE  DEFAULT_MSG
                                        TRIANGLE_LIBRARY TRIANGLE_INCLUDE_DIR)

      mark_as_advanced(TRIANGLE_INCLUDE_DIR TRIANGLE_LIBRARY)
    endif(TRIANGLE_LIBRARY)  
  endif()

endif()


