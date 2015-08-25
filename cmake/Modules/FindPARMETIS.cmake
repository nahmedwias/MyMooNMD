# - Try to find PARMETIS parallel graph partitioning package.
# written by Clemens Bartsch
# date: 24 August 2015
# This is a rudimentary find module, which does only look for the PARMETIS library shipped with ParMooN.
# TODO Extend to finding system PARMETIS.
#
# Once done this will define
#  PARMETIS_FOUND - System has PARMETIS
#  PARMETIS_INCLUDE_DIRS - The PARMETIS include directories
#  PARMETIS_LIBRARIES - The libraries needed to use PARMETIS
#  TODO(?) PARMETIS_DEFINITIONS - Compiler switches required for using PARMETIS

if (PARMETIS_LIBRARIES AND PARMETIS_INCLUDE_DIRS)
    message("PARMETIS package already found. Not searching again.")
else(PARMETIS_LIBRARIES AND PARMETIS_INCLUDE_DIRS)
    
    #Check if an MPI implementation is there. If not so, find one with the default module. 
    find_package(MPI REQUIRED)
    
    if(FIND_USER_PARMETIS)
        message("Find user PARMETIS is not yet implemented.")
    endif(FIND_USER_PARMETIS)
      
    # Search for the library exclusively in the ParMooN EXT_LIB path.
    if(NOT PARMETIS_LIBRARY OR NOT PARMETIS_INCLUDE_DIR)
        message("Searching in ParMooN EXT_LIB path. Selected architecture ARCH=${ARCH}")
        find_path(PARMETIS_INCLUDE_DIR  parmetis.h PATHS ${PARMOON_EXTLIB_PATH}/Metis/Lib NO_DEFAULT_PATH)
        find_library(PARMETIS_LIBRARY NAMES parmetis_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Metis NO_DEFAULT_PATH)
        message("Found library: " ${PARMETIS_LIBRARY})
    endif(NOT PARMETIS_LIBRARY OR NOT PARMETIS_INCLUDE_DIR)
    
    set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${MPI_LIBRARIES} )
    set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR} ${MPI_INCLUDE_DIRS})
    
    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args(PARMETIS  DEFAULT_MSG
                                  PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR)

    mark_as_advanced(PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY )

endif(PARMETIS_LIBRARIES AND PARMETIS_INCLUDE_DIRS)