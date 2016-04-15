# CMakeLists.txt for subdirectory Solver of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/Solver")

# Source files to be added to the 2D and 3D library.
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Solver.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/DirectSolver.C")

list(APPEND PARMOON_SOURCES_2D ${SOLVER_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${SOLVER_SOURCES})
