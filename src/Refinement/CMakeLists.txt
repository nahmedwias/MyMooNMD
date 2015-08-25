# CMakeLists.txt for subdirectory Refinement of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# TODO This include might be superfluous - it is done in the main CMakeLists.
include_directories("${CMAKE_SOURCE_DIR}/include/Refinement")

# TODO Replace this by a call-every-file approach, which is the recommended approach.
# (background: with globbing CMake will not notice when there are any files added
#  after its first run)
file(GLOB_RECURSE REF_SOURCES "${PROJECT_SOURCE_DIR}/src/Refinement/*.C")


# Define two static libraries. Each has its own precompiler flag (2D/3D)
add_library(ref_2d STATIC ${REF_SOURCES})
target_compile_definitions(ref_2d PUBLIC -D__2D__)

add_library(ref_3d STATIC ${REF_SOURCES})
target_compile_definitions(ref_3d PUBLIC -D__3D__)

