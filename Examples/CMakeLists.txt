# CMakeLists.txt for subdirectory Examples of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#


# Include header files from subdirectories.
include_directories("${CMAKE_SOURCE_DIR}/Examples"
                    "${CMAKE_SOURCE_DIR}/Examples/CD_2D"
                    "${CMAKE_SOURCE_DIR}/Examples/CD_3D"
                    "${CMAKE_SOURCE_DIR}/Examples/Darcy_2D"
                    "${CMAKE_SOURCE_DIR}/Examples/NSE_2D"
                    "${CMAKE_SOURCE_DIR}/Examples/NSE_3D"
                    "${CMAKE_SOURCE_DIR}/Examples/TCD_2D"
                    "${CMAKE_SOURCE_DIR}/Examples/TCD_3D"
                    "${CMAKE_SOURCE_DIR}/Examples/TNSE_2D")

# Source files to be added to the 2D library.
list(APPEND EXAMPLE_SOURCES_2D "${CMAKE_CURRENT_SOURCE_DIR}/Example2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${CMAKE_CURRENT_SOURCE_DIR}/Example_CD2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${CMAKE_CURRENT_SOURCE_DIR}/Example_Darcy2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${CMAKE_CURRENT_SOURCE_DIR}/Example_NSE2D.C")

# Source files to be added to the 3D library.
list(APPEND EXAMPLE_SOURCES_3D "${CMAKE_CURRENT_SOURCE_DIR}/Example3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${CMAKE_CURRENT_SOURCE_DIR}/Example_CD3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${CMAKE_CURRENT_SOURCE_DIR}/Example_NSE3D.C")

# Define two static libraries. Each has its own precompiler flag (2D/3D). 
add_library(examples_2d STATIC ${EXAMPLE_SOURCES_2D})
target_compile_definitions(examples_2d PUBLIC -D__2D__)

add_library(examples_3d STATIC ${EXAMPLE_SOURCES_3D})
target_compile_definitions(examples_3d PUBLIC -D__3D__)