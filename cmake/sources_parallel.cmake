# CMakeLists.txt for subdirectory Parallel of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/Parallel")

list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/MeshPartition.C") 
list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/MumpsSolver.C")
list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParDirectSolver.C")
list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParDiso.C")
list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParFEMapper3D.C")
list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParFECommunicator3D.C")


list(APPEND PARMOON_SOURCES_2D ${PAR_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${PAR_SOURCES})