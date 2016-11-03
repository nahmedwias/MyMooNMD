# CMakeLists.txt for subdirectory Examples of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#


# Include header files from subdirectories.
include_directories("${CMAKE_SOURCE_DIR}/include/Examples"
                    "${CMAKE_SOURCE_DIR}/include/Examples/CD_2D"
                    "${CMAKE_SOURCE_DIR}/include/Examples/CD_3D"
                    "${CMAKE_SOURCE_DIR}/include/Examples/Darcy_2D"
                    "${CMAKE_SOURCE_DIR}/include/Examples/NSE_2D"
                    "${CMAKE_SOURCE_DIR}/include/Examples/NSE_3D"
                    "${CMAKE_SOURCE_DIR}/include/Examples/TCD_2D"
                    "${CMAKE_SOURCE_DIR}/include/Examples/TCD_3D"
                    "${CMAKE_SOURCE_DIR}/include/Examples/TNSE_2D"
		    "${CMAKE_SOURCE_DIR}/include/Examples/Brinkman_2D"
		    "${CMAKE_SOURCE_DIR}/include/Examples/Brinkman_3D")

# Source files to be added to the 2D library.
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_NonStationary2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_CD2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_TimeCD2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_Darcy2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_NSE2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_TimeNSE2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_Brinkman2D.C")

# Source files to be added to the 3D library.
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example_NonStationary3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example_CD3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example_TimeCD3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example_NSE3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example_TimeNSE3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example_Brinkman3D.C")

list(APPEND PARMOON_SOURCES_2D ${EXAMPLE_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${EXAMPLE_SOURCES_3D})
