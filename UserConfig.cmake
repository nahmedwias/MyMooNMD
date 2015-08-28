# ===================================================================
# This is a user configuration file for ParMooN, Version 1.1,
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 05 June 2015
#
# Change History 
# 2015/08/20 Clemens Bartsch Reworked the file with the directive:
# the user should set all cached variables only once and never edit
# anything in CMakeCache.txt. The only variables that may be edited
# here are the local variables _OUTPUT_DIR_PATH and _PARALLEL_TYPE.
# ===================================================================

################################################################################
# Setting User Cache variables.
#
# These variables should be set only once, before running cmake the first time.
# Changing them here later has no effect, and changing them in the cache might
# lead to unwanted behavior. So if you really need to change any of them later,
# delete CMakeCache.txt and directory CMakeFiles and rerun cmake.
################################################################################

# Enable verbose output from MakeFiles. Would be set by "project"
# command, defaults to FALSE.
set(CMAKE_VERBOSE_MAKEFILE TRUE CACHE STRING "Do not edit in CMakeCache.txt!")

###### Variables dealing with the finding of external libraries.

# Additional non-standard search path for extern include folders.
# Note: This path is among those searched by find_path by default!
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "my/special/path" CACHE PATH "Do not edit in CMakeCache.txt!")

# Additional non-standard search path for extern libraries.
# Note: This path is among those searched by find_library by default!
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "my/special/path" CACHE PATH "Do not edit in CMakeCache.txt!")

# If true: use cmake default blas search module before using ParMooN's own. 
#set(DEFAULT_MODULE_BLAS FALSE CACHE STRING "Do not edit in CMakeCache.txt!")

# If true: use cmake default lapack search module before using ParMooN's own. 
#set(DEFAULT_MODULE_LAPACK FALSE CACHE STRING "Do not edit in CMakeCache.txt!")

#The following variables take effect only when using ParMooN's search modules.
# Try to find user blas library before using ParMooN's own.
set(FIND_USER_BLAS FALSE CACHE STRING "Do not edit in CMakeCache.txt!")

# Try to find user lapack library before using ParMooN's own.
set(FIND_USER_LAPACK FALSE CACHE STRING "Do not edit in CMakeCache.txt!")

# Try to find user umfpack library before using ParMooN's own.
set(FIND_USER_UMFPACK FALSE CACHE STRING "Do not edit in CMakeCache.txt!")

# Try to find user mumps library before using ParMooN's own.
set(FIND_USER_MUMPS FALSE CACHE STRING "Do not edit in CMakeCache.txt!")

# Try to find user tetgen library before using ParMooN's own.
set(FIND_USER_TETGEN FALSE CACHE STRING "Do not edit in CMakeCache.txt!")

# Try to find user tecplot library before using ParMooN's own.
set(FIND_USER_TECPLOT FALSE CACHE STRING "Do not edit in CMakeCache.txt!")

# Try to find user triangle library before using ParMooN's own.
set(FIND_USER_TRIANGLE FALSE CACHE STRING "Do not edit in CMakeCache.txt!")

################################################################################
# Setting User non-Cache variables.
#
# These are non-cache variables which can be changed here deliberately. They
# will take effect on every new run of cmake, without deleting anything before.
################################################################################
# Selection of the program type. Chose between SEQUENTIAL MPI OMPONLY HYBRID,
# entering anything else will lead to undefined behaviour. TODO Catch that case!
set(_PARALLEL_TYPE "MPI")

# Select the output directory for the binaries. Make sure the directory exists in your file system!
set(_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/")