# CMakeLists.txt for subdirectory FE of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# Include header files. 
include_directories("${CMAKE_SOURCE_DIR}/include/ROM")

# We Need Boost For Brush
find_package(Boost 1.47.0 REQUIRED COMPONENTS regex filesystem system program_options)
IF(Boost_FOUND)
    MESSAGE("Boost includes located in: ${Boost_INCLUDE_DIRS}")
    MESSAGE("Boost libraries located in: ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIRS})
    list(APPEND _EXTERN_LIBRARIES ${Boost_LIBRARIES})
ENDIF(Boost_FOUND)

#set( ENV{BLA_VENDOR} "ATLAS" )
find_package( BLAS REQUIRED )

# Source files used in 2D and 3D.
list(APPEND _ROM_SOURCES "${PROJECT_SOURCE_DIR}/src/ROM/POD.C")
list(APPEND _ROM_SOURCES "${PROJECT_SOURCE_DIR}/src/ROM/SNAPS.C")
list(APPEND _ROM_SOURCES "${PROJECT_SOURCE_DIR}/src/ROM/POD_TCDR2D.C")
list(APPEND _ROM_SOURCES "${PROJECT_SOURCE_DIR}/src/ROM/ROM_TCDR2D.C")
list(APPEND _ROM_SOURCES "${PROJECT_SOURCE_DIR}/src/ROM/ROM.C")

list(APPEND PARMOON_SOURCES_2D ${_ROM_SOURCES})
