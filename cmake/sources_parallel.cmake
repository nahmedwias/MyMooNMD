# CMakeLists.txt for subdirectory Parallel of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# TODO This include might be superfluous - it is done in the main CMakeLists.
include_directories("${CMAKE_SOURCE_DIR}/include/Parallel")

list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/MeshPartition.C") 
list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/MumpsSolver.C")
list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParDirectSolver.C")
list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParDiso.C")
list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParFEMapper3D.C")
list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParFECommunicator3D.C")

# Define two static libraries. Each has its own precompiler flag (2D/3D)
add_library(par_2d STATIC ${PAR_SOURCES})
target_compile_definitions(par_2d PUBLIC -D__2D__)

add_library(par_3d STATIC ${PAR_SOURCES})
target_compile_definitions(par_3d PUBLIC -D__3D__)
