# ===================================================================
# This is FindACMLBLAS.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a ACMLBLAS lib in the system 
# if found, this will define
#  ACMLBLAS_FOUND - System has ACMLBLAS
#  ACMLBLAS_INCLUDE_DIRS - The ACMLBLAS include directories
#  ACMLBLAS_LIBRARIES - The libraries needed to use ACMLBLAS
# ===================================================================
if(ACMLBLAS_INCLUDES AND ACMLBLAS_LIBRARIES)
  set(ACMLBLAS_FIND_QUIETLY TRUE)
endif(ACMLBLAS_INCLUDES AND ACMLBLAS_LIBRARIES)

if(NOT ACMLBLAS_FOUND)
 
  find_path(ACMLBLAS_INCLUDE_DIR   acml.h PATHS $ENV{ACMLBLASDIR}/include ${CMAKE_INCLUDE_PATH})
  find_library(ACMLBLAS_LIBRARY NAMES acml PATHS $ENV{ACMLBLASDIR}/lib ${CMAKE_LIBRARY_PATH})
  get_filename_component(ACMLBLAS_LIBDIR ${ACMLBLAS_LIBRARY} PATH)
  find_library(ACMLBLAS_LIBRARY_MP NAMES acml_mv PATHS ${ACMLBLAS_LIBDIR})
   if(NOT ACMLBLAS_LIBRARY_MP)
      find_library(ACMLBLAS_LIBRARY_MP NAMES acml_mp PATHS ${ACMLBLAS_LIBDIR})
    endif(NOT ACMLBLAS_LIBRARY_MP)  
      
  if(NOT ACMLBLAS_LIBRARY)
    message("ACMLBLAS not found in the system, so checking the availability in ParMooN for the selected ARCH=${ARCH}")
    find_path(ACMLBLAS_INCLUDE_DIR  acml.h PATHS ${PARMOON_EXTLIB_PATH}/ACML/gfortran64/include)
    find_library(ACMLBLAS_LIBRARY NAMES acml PATHS ${PARMOON_EXTLIB_PATH}/ACML/gfortran64/lib) 
    get_filename_component(ACMLBLAS_LIBDIR ${ACMLBLAS_LIBRARY} PATH)
    find_library(ACMLBLAS_LIBRARY_MP NAMES acml_mv PATHS ${ACMLBLAS_LIBDIR})
    if(NOT ACMLBLAS_LIBRARY_MP)
      find_library(ACMLBLAS_LIBRARY_MP NAMES acml_mp PATHS ${ACMLBLAS_LIBDIR})  
    endif(NOT ACMLBLAS_LIBRARY_MP)       
  endif(NOT ACMLBLAS_LIBRARY)
  
  if(ACMLBLAS_LIBRARY)  
    # combine mumps and its deps    
    if(ACMLBLAS_LIBRARY_MP)
      set(ACMLBLAS_LIBRARY  ${ACMLBLAS_LIBRARY}  ${ACMLBLAS_LIBRARY_MP} ) 
    else(ACMLBLAS_LIBRARY_MP)   
      set(ACMLBLAS_LIBRARY FALSE)
    endif(ACMLBLAS_LIBRARY_MP)
    
    # set ACMLBLAS
    if(ACMLBLAS_LIBRARY)
      include(FindPackageHandleStandardArgs)
    
      set(ACMLBLAS_LIBRARIES ${ACMLBLAS_LIBRARY})
      set(ACMLBLAS_INCLUDE_DIRS ${ACMLBLAS_INCLUDE_DIR})

      # handle the QUIETLY and REQUIRED arguments and set ACMLBLAS_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(ACMLBLAS  DEFAULT_MSG
                                        ACMLBLAS_LIBRARY ACMLBLAS_INCLUDE_DIR)

      mark_as_advanced(ACMLBLAS_INCLUDE_DIR ACMLBLAS_LIBRARY)
    endif(ACMLBLAS_LIBRARY)  
  endif()

endif()


