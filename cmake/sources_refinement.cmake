# CMakeLists.txt for subdirectory Refinement of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/Refinement")

# TODO Replace this by a call-every-file approach, which is the recommended approach.
# (background: with globbing CMake will not notice when there are any files added
#  after its first run)
file(GLOB_RECURSE REF_SOURCES "${PROJECT_SOURCE_DIR}/src/Refinement/*.C")


list(APPEND PARMOON_SOURCES_2D ${REF_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${REF_SOURCES})