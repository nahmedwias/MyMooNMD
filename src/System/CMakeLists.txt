# CMakeLists.txt for subdirectory System of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# TODO This include might be superfluous - it is done in the main CMakeLists.
include_directories("${CMAKE_SOURCE_DIR}/include/System")

# Source files to be added to both 2D and 3D library.
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/Assemble.C") 

# Source files to be added to the 2D library.
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/CD2D.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/NSE2D.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/SystemMat2D.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/SystemMatDarcy2D.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/SystemMatNSE2D.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/SystemMatScalar2D.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/SystemMatTimeScalar2D.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/SystemMatTNSE2D.C")
# TODO CB Will the following files be replaced/renamed?
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/SystemNSE2D.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/SystemTNSE2D.C")

# Source files to be added to the 3D library.
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/AssembleMat3D.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/SystemMatScalar3D.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/SystemCD3D.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/SystemNSE3D.C")
# TODO CB Will the following file be replaced/renamed?
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/SystemTCD3D.C")
#list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/SystemMatTimeScalar3D.C") #CB Why not?
#list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/SystemMatNSE3D.C") #CB Why not?

# Define two static libraries. Each has its own precompiler flag (2D/3D). 
add_library(sys_2d STATIC ${SYST_SOURCES} ${SYST_SOURCES_2D})
target_compile_definitions(sys_2d PUBLIC -D__2D__)

add_library(sys_3d STATIC ${SYST_SOURCES} ${SYST_SOURCES_3D})
target_compile_definitions(sys_3d PUBLIC -D__3D__)
