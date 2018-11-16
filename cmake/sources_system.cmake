# CMakeLists.txt for subdirectory System of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/System")

# Source files to be added to the 2D library.
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/ConvectionDiffusion.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/NSE2D.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/Darcy.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/Time_CD2D.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/Time_NSE2D.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/Brinkman2D.C")

# Source files to be added to the 3D library.
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/ConvectionDiffusion.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/NSE3D.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/Darcy.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/Time_CD3D.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/Time_NSE3D.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/Brinkman3D.C")

# Source files to be added to both 2D and 3D library
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/Residuals.C")


list(APPEND PARMOON_SOURCES_2D ${SYST_SOURCES} ${SYST_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${SYST_SOURCES} ${SYST_SOURCES_3D})
