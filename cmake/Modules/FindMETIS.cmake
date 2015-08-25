# - Try to find METIS graph partitioning package.
# written by Clemens Bartsch
# date: 24 August 2015
# This is a rudimentary find module, which does only look for the METIS library shipped with ParMooN.
# TODO Extend to finding system METIS.
#
# Once done this will define
#  METIS_FOUND - System has METIS
#  METIS_INCLUDE_DIRS - The METIS include directories
#  METIS_LIBRARIES - The libraries needed to use METIS
#  TODO(?) METIS_DEFINITIONS - Compiler switches required for using METIS

if (METIS_LIBRARIES AND METIS_INCLUDE_DIRS)
    message("METIS package already found. Not searching again.")
else(METIS_LIBRARIES AND METIS_INCLUDE_DIRS)

    if(FIND_USER_METIS)
      message("Find user METIS is not yet implemented.")
    endif(FIND_USER_METIS)
    
      # Search for the library exclusively in the ParMooN EXT_LIB path.
    if(NOT METIS_LIBRARY OR NOT METIS_INCLUDE_DIR)
        message("Searching in ParMooN EXT_LIB path. Selected architecture ARCH=${ARCH}")
        find_path(METIS_INCLUDE_DIR  metis.h PATHS ${PARMOON_EXTLIB_PATH}/Metis/Lib NO_DEFAULT_PATH)
        find_library(METIS_LIBRARY NAMES metis_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Metis NO_DEFAULT_PATH)
        message("Found library: " ${METIS_LIBRARY})
    endif(NOT METIS_LIBRARY OR NOT METIS_INCLUDE_DIR)
    
    set(METIS_LIBRARIES ${METIS_LIBRARY} )
    set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR} )

    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args(METIS  DEFAULT_MSG
                                  METIS_LIBRARY METIS_INCLUDE_DIR)

    mark_as_advanced(METIS_INCLUDE_DIR METIS_LIBRARY )

endif(METIS_LIBRARIES AND METIS_INCLUDE_DIRS)