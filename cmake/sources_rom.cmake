# CMakeLists.txt for subdirectory Solver of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/ReducedOrderModels")

# Source files to be added to the 2D and 3D library.
list(APPEND ROM_SOURCES "${PROJECT_SOURCE_DIR}/src/ReducedOrderModels/SnapshotsCollector.C")
list(APPEND ROM_SOURCES "${PROJECT_SOURCE_DIR}/src/ReducedOrderModels/POD.C")
list(APPEND ROM_SOURCES "${PROJECT_SOURCE_DIR}/src/ReducedOrderModels/TimeConvectionDiffusionPOD.C")
list(APPEND ROM_SOURCES "${PROJECT_SOURCE_DIR}/src/ReducedOrderModels/ROM.C")
list(APPEND ROM_SOURCES "${PROJECT_SOURCE_DIR}/src/ReducedOrderModels/TimeConvectionDiffusionROM.C")

list(APPEND PARMOON_SOURCES_2D ${ROM_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${ROM_SOURCES})
