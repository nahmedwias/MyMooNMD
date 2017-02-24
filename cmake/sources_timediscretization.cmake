include_directories("${CMAKE_SOURCE_DIR}/include/TimeDiscretization")

# Source files used in 2D and 3D.
list(APPEND TD_SOURCES "${PROJECT_SOURCE_DIR}/src/TimeDiscretization/Time_NSE2D_Merged.C")
list(APPEND TD_SOURCES "${PROJECT_SOURCE_DIR}/src/TimeDiscretization/RungeKuttaTable.C")
list(APPEND TD_SOURCES "${PROJECT_SOURCE_DIR}/src/TimeDiscretization/TimeDiscretization.C")
list(APPEND TD_SOURCES "${PROJECT_SOURCE_DIR}/src/TimeDiscretization/assemble_routine_tnse2D_smagorinsky.C")


list(APPEND PARMOON_SOURCES_2D ${TD_SOURCES} )
list(APPEND PARMOON_SOURCES_3D ${TD_SOURCES} )