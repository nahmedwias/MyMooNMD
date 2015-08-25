# ===================================================================
# This is FindBLAS.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a BLAS lib in the system 
# if found, this will define
#  BLAS_FOUND - System has BLAS
#  BLAS_INCLUDE_DIRS - The BLAS include directories
#  BLAS_LIBRARIES - The libraries needed to use BLAS
#
# History:
# 2015/08/18 Clemens Bartsch: Introduced variable FIND_USER_BLAS 
#	for activating non-ParMooN search. Added NO_DEFAULT_PATH
#	keyword. Slightly changed the output.
# 2015/08/20 Clemens Bartsch: Placed standard handling of arguments
#         outside the if(..._library) block - that was a bug.
# ===================================================================
if(BLAS_INCLUDES AND BLAS_LIBRARIES)
  set(BLAS_FIND_QUIETLY TRUE)
endif(BLAS_INCLUDES AND BLAS_LIBRARIES)


if(NOT BLAS_FOUND)
 
  # Search for the Library in standard non-ParMooN paths and in those specified after PATHS.
  if(FIND_USER_BLAS)
    message("Searching in default and user paths.")
    find_path(BLAS_INCLUDE_DIR   acml.h PATHS $ENV{ACMBLASDIR}/include )
    find_library(BLAS_LIBRARY NAMES acml PATHS $ENV{BLASDIR}/lib)
    get_filename_component(_BLAS_LIBDIR ${BLAS_LIBRARY} PATH)
    find_library(BLAS_LIBRARY_MP NAMES acml_mv PATHS ${_BLAS_LIBDIR})
    if(NOT BLAS_LIBRARY_MP)
      find_library(BLAS_LIBRARY_MP NAMES acml_mp PATHS ${_BLAS_LIBDIR})
    endif(NOT BLAS_LIBRARY_MP)  
  endif(FIND_USER_BLAS)
  
  # Search for the library exclusively in the ParMooN EXT_LIB path.
  if(NOT BLAS_LIBRARY)
    message("Searching in ParMooN EXT_LIB path.")
    find_path(BLAS_INCLUDE_DIR  acml.h PATHS ${PARMOON_EXTLIB_PATH}/ACML/gfortran64/include NO_DEFAULT_PATH)
    find_library(BLAS_LIBRARY NAMES acml PATHS ${PARMOON_EXTLIB_PATH}/ACML/gfortran64/lib NO_DEFAULT_PATH) 
    get_filename_component(_BLAS_LIBDIR ${BLAS_LIBRARY} PATH)
    find_library(BLAS_LIBRARY_MP NAMES acml_mv PATHS ${_BLAS_LIBDIR} NO_DEFAULT_PATH)
    if(NOT BLAS_LIBRARY_MP)
      find_library(BLAS_LIBRARY_MP NAMES acml_mp PATHS ${_BLAS_LIBDIR} NO_DEFAULT_PATH)  
    endif(NOT BLAS_LIBRARY_MP)       
  endif(NOT BLAS_LIBRARY)
  
  if(BLAS_LIBRARY)  
    # combine blas and its deps    
    if(BLAS_LIBRARY_MP)
      set(BLAS_LIBRARY  ${BLAS_LIBRARY}  ${BLAS_LIBRARY_MP} ) 
    else(BLAS_LIBRARY_MP)   
      set(BLAS_LIBRARY FALSE)
    endif(BLAS_LIBRARY_MP)
  endif(BLAS_LIBRARY) 
  
  # Handling of standard arguments.
  include(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set BLAS_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(BLAS  DEFAULT_MSG
                                    BLAS_LIBRARY BLAS_INCLUDE_DIR)
  mark_as_advanced(BLAS_INCLUDE_DIR BLAS_LIBRARY)

  # Set two non-cache variables to be use when including and linking against BLAS.
  set(BLAS_LIBRARIES ${BLAS_LIBRARY})
  set(BLAS_INCLUDE_DIRS ${BLAS_INCLUDE_DIR})


endif()


