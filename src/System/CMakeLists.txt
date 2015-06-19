
include_directories("${CMAKE_SOURCE_DIR}/include/System")

list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/Assemble.C") 

if("${GEO}" STREQUAL "2D")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatNSE2D.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatScalar2D.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatTimeScalar2D.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatTimeScalar2D_ALE.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatTNSE2D.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatTNSE2D_ALE.C")
elseif("${GEO}" STREQUAL "3D")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/AssembleMat3D.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatScalar3D.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatTimeScalar3D.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatNSE3D.C")
endif("${GEO}" STREQUAL "2D")

add_library(sys STATIC ${SYST_SOURCES})
