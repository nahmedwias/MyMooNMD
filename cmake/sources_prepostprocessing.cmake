# CMakeLists.txt for subdirectory System of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
#

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/PrePostProcessing")

list(APPEND PrePostProcess "${PROJECT_SOURCE_DIR}/src/PrePostProcessing/PrePost_Cylinder_Square.C")
# Source files used in 2D and 3D.

list(APPEND PARMOON_SOURCES_3D ${PrePostProcess} )
