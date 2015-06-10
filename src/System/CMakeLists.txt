
include_directories("${CMAKE_SOURCE_DIR}/include/System")

list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/Assemble.C") 
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/AssembleMat3D.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatScalar3D.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatTimeScalar3D.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/SystemMatNSE3D.C")

add_library(sys STATIC ${SYST_SOURCES})
