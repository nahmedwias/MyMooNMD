# CMakeLists.txt for subdirectory QuadFormulas of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# TODO This include might be superfluous - it is done in the main CMakeLists.
include_directories("${CMAKE_SOURCE_DIR}/include/QuadFormulas")

# Source files used in 2D and 3D (CB that correct? Could be seperated, probably).
list(APPEND QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/QuadFormula.C")
list(APPEND QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/QuadFormula1D.C")
list(APPEND QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/QuadFormula2D.C")
list(APPEND QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/QuadFormula3D.C")
list(APPEND QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/QuadFormulaHexa.C")
list(APPEND QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/QuadFormulaQuad.C")
list(APPEND QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/QuadFormulaTetra.C")
list(APPEND QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/QuadFormulaTria.C")

# Define two static libraries. Each has its own precompiler flag (2D/3D)
add_library(quad_2d STATIC ${QUAD_SOURCES})
target_compile_definitions(quad_2d PUBLIC -D__2D__)
 
add_library(quad_3d STATIC ${QUAD_SOURCES})
target_compile_definitions(quad_3d PUBLIC -D__3D__)

