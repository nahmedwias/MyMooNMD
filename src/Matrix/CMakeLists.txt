# CMakeLists.txt for subdirectory Matrix of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/25 Ulrich Wilbrandt: Adjust to Clemens changes to the cmake system
#

include_directories("${CMAKE_SOURCE_DIR}/include/Matrix")

# Source files to be added to the 2D library.
list(APPEND MATRIX_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockMatrix.C")
list(APPEND MATRIX_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockMatrixCD2D.C")
list(APPEND MATRIX_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockMatrixDarcy2D.C")
list(APPEND MATRIX_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockMatrixNSE2D.C")
list(APPEND MATRIX_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockPattern.C")
list(APPEND MATRIX_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockVector.C")
#list(APPEND MATRIX_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockVector2D.C")

list(APPEND MATRIX_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockMatrix.C")
list(APPEND MATRIX_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockMatrix3D.C")
list(APPEND MATRIX_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockMatrixCD3D.C")
list(APPEND MATRIX_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockPattern.C")
list(APPEND MATRIX_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockVector.C")
#list(APPEND MATRIX_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Matrix/BlockVector3D.C")


# Define two static libraries. Each has its own precompiler flag (2D/3D). 
add_library(matrix_2d STATIC ${MATRIX_SOURCES_2D})
target_compile_definitions(matrix_2d PUBLIC -D__2D__)

add_library(matrix_3d STATIC ${MATRIX_SOURCES_3D})
target_compile_definitions(matrix_3d PUBLIC -D__3D__)
