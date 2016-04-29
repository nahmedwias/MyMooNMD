# CMakeLists.txt for subdirectory Solver of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/Solver")

# Source files to be added to the 2D and 3D library.
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/DirectSolver.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_bicgstab.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_cg.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_cgs.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_gmres.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_jacobi.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_richardson.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Saddle_point_preconditioner.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Solver.C")

list(APPEND PARMOON_SOURCES_2D ${SOLVER_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${SOLVER_SOURCES})
