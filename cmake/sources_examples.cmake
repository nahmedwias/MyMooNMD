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
		    "${CMAKE_SOURCE_DIR}/user_projects/include/Examples"
		    "${CMAKE_SOURCE_DIR}/user_projects/include/Examples/CD2D"
		    "${CMAKE_SOURCE_DIR}/user_projects/include/Examples/NSE2D"
		    "${CMAKE_SOURCE_DIR}/user_projects/include/Examples/Time_NSE2D"
		    "${CMAKE_SOURCE_DIR}/user_projects/include/Examples/Time_CD2D")

list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/BoundaryCondition.C")
list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/BoundaryData.C")
list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/PDECoefficients.C")
list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/AnalyticalFunction.C")
list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/Example.C")
list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/Example_Brinkman.C")
list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/Example_ConvDiff.C")
list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/Example_Darcy.C")
list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/Example_NavierStokes.C")
list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/Example_NonStationary.C")
list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/Example_Time_ConvDiff.C")
list(APPEND EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/src/Examples/Example_Time_NavierStokes.C")

# Source files to be added to the 2D library.
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_NonStationary2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_CD2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_TimeCD2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_Darcy2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_NSE2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_TimeNSE2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Examples/Example_Brinkman2D.C")
list(APPEND EXAMPLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/Examples/Example_TimeLinElastic2D.C")

# Source files to be added to the 3D library.
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example_NonStationary3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example_CD3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example_TimeCD3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example_NSE3D.C")
list(APPEND EXAMPLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Examples/Example_TimeNSE3D.C")


list(APPEND PARMOON_SOURCES_2D ${EXAMPLE_SOURCES_2D} ${EXAMPLE_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${EXAMPLE_SOURCES_3D} ${EXAMPLE_SOURCES})
