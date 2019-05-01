# CMakeLists.txt for subdirectory AssembleRoutines of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2017/03/22 Najib Alia: Rework to separate FE and AssembleRoutines files.
#

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/PrePostProcessing")

list(APPEND ASSEMBLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/PrePostProcessing/ChannelFlowRoutines.C")
list(APPEND PARMOON_SOURCES_3D ${ASSEMBLE_SOURCES} ${ASSEMBLE_SOURCES_3D})
