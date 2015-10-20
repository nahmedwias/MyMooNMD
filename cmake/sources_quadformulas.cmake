# CMakeLists.txt for subdirectory QuadFormulas of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# Include header files. 
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


list(APPEND PARMOON_SOURCES_2D ${QUAD_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${QUAD_SOURCES})