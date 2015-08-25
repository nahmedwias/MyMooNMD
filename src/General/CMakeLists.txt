# CMakeLists.txt for subdirectory General of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# TODO This include might be superfluous - it is done in the main CMakeLists.
include_directories("${CMAKE_SOURCE_DIR}/include/General")

# Source files used in 2D and 3D.
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/Database.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/MainUtilities.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/MooNMD_Io.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/ReadParam.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/Utilities.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/Vector.C") 

# Define two static libraries. Each has its own precompiler flag (2D/3D)
add_library(gen_2d STATIC ${GENERAL_SOURCES})
target_compile_definitions(gen_2d PUBLIC -D__2D__)

add_library(gen_3d STATIC ${GENERAL_SOURCES})
target_compile_definitions(gen_3d PUBLIC -D__3D__)