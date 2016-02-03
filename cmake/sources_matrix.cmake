# CMakeLists.txt for subdirectory Matrix of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/25 Ulrich Wilbrandt: Adjust to Clemens changes to the cmake system
#

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/Matrix")

# Source files to be added to the 2D and 3D library.
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/BlockMatrix.C")
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/BlockFEMatrix.C")
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/BlockVector.C")


list(APPEND PARMOON_SOURCES_2D ${MATRIX_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${MATRIX_SOURCES})
