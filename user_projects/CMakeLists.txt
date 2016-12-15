# this is the place to put the code which cmake needs to compile your project.
# Feel free to create a subdirectory for it, if you like.
#
# If you have header files in this directory, you can do
#    include_directories("${PROJECT_SOURCE_DIR}/user_projects")
# Than the compile will be able to find these header files. If you have 
# additional source files, i.e. source files without the main method, you can
# create a library for your project 
#    list(APPEND PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/project_file1.C")
#    list(APPEND PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/project_file2.C")
#    add_library(PROJECT_LIB_2D STATIC ${PROJECT_SOURCES_2D})
#    target_compile_definitions(PROJECT_LIB_2D PUBLIC -D__2D__)
# Then you can create an executable as follows
#    add_executable(project ${PROJECT_SOURCE_DIR}/user_projects/project_file_with_main.C)
#    set_target_properties(project PROPERTIES OUTPUT_NAME parMooN_project_${PARMOON_PARALLEL_TYPE})
#    target_compile_definitions(project PUBLIC -D__2D__)
#    target_link_libraries(project PROJECT_LIB_2D parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
# 
# Then cmake will create a makefile with a target 'project', so typing 
#    make project
# will compile your project. The executable will be called parMooN_project_X
# where 'X' is the build type, e.g. SEQUENTIAL.
####################################################################################################

# Here begins the code which enables usage of our project specific Coupled CDR stuff.
# Clemens Bartsch, 2015/10/27

# Include header files from inc/ and examples/ directory. 
include_directories("${PROJECT_SOURCE_DIR}/user_projects/inc")
include_directories("${PROJECT_SOURCE_DIR}/user_projects/examples")

# We Need Boost For Brush
find_package(Boost 1.47.0 REQUIRED COMPONENTS regex filesystem system program_options)
IF(Boost_FOUND)
    MESSAGE("Boost includes located in: ${Boost_INCLUDE_DIRS}")
    MESSAGE("Boost libraries located in: ${Boost_LIBRARIES}")
ENDIF(Boost_FOUND)

# Shuffle everything into a static project library.
list(APPEND _PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/BrushWrapper.C")
list(APPEND _PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/CoupledCDR_2D.C")
list(APPEND _PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/Coupled_Time_CDR_2D.C")
list(APPEND _PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/Example_CoupledCDR2D.C")
list(APPEND _PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/Example_TimeCoupledCDR2D.C")
list(APPEND _PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/ReactionCoupling.C")

add_library(PROJECT_LIB_2D STATIC ${_PROJECT_SOURCES_2D})
target_compile_definitions(PROJECT_LIB_2D PUBLIC -D__2D__)

# Brush will need OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Include the Brush solver libraries..
list(APPEND BRUSH_LIBS ${BRUSH_LIB_DIR}/libbrush.a)
list(APPEND BRUSH_LIBS ${BRUSH_LIB_DIR}/libcamxml.a)
list(APPEND BRUSH_LIBS ${BRUSH_LIB_DIR}/libchemkinReader.a)
list(APPEND BRUSH_LIBS ${BRUSH_LIB_DIR}/libcomostrings.a)
list(APPEND BRUSH_LIBS ${BRUSH_LIB_DIR}/libgeometry.a)
list(APPEND BRUSH_LIBS ${BRUSH_LIB_DIR}/libsprog.a)
list(APPEND BRUSH_LIBS ${BRUSH_LIB_DIR}/libsweep.a)
list(APPEND BRUSH_LIBS ${MOON_GEO_LIB_DIR}/libgeomoon.a)
list(APPEND BRUSH_LIBS ${Boost_LIBRARIES})
list(APPEND BRUSH_LIBS ${BRUSH_LIBS}) #double it, to be on the safe side
#...and the corresponding source code.
list(APPEND BRUSH_INCLUDES "${BRUSH_SOURCE_DIR}/brush/include/")
include_directories(${BRUSH_INCLUDES})

# Add an executable target from CoupledCDR.C
add_executable(ccdr2d ${PROJECT_SOURCE_DIR}/user_projects/CoupledCDR2D_ParMooN.C)
set_target_properties(ccdr2d PROPERTIES OUTPUT_NAME parMooN_ccdr2d_${PARMOON_PARALLEL_TYPE})
target_compile_definitions(ccdr2d PUBLIC -D__2D__)
target_link_libraries(ccdr2d PROJECT_LIB_2D parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

# Add an executable target from Coupled_Time_CDR2D_ParMooN.C
add_executable(time_ccdr2d ${PROJECT_SOURCE_DIR}/user_projects/Coupled_Time_CDR2D_ParMooN.C)
set_target_properties(time_ccdr2d PROPERTIES OUTPUT_NAME parMooN_time_ccdr2d_${PARMOON_PARALLEL_TYPE})
target_compile_definitions(time_ccdr2d PUBLIC -D__2D__)
target_link_libraries(time_ccdr2d PROJECT_LIB_2D parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

# Add an executable target from Flow_And_Concentrations_ParMooN.c++
add_executable(flow_and_conc_2d ${PROJECT_SOURCE_DIR}/user_projects/Flow_Concentration_Particles_ParMooN.c++)
set_target_properties(flow_and_conc_2d PROPERTIES OUTPUT_NAME parMooN_flow_and_concentrations_2d_${PARMOON_PARALLEL_TYPE})
target_compile_definitions(flow_and_conc_2d PUBLIC -D__2D__)
target_link_libraries(flow_and_conc_2d PROJECT_LIB_2D parmoon_2d_${PARMOON_PARALLEL_TYPE} ${BRUSH_LIBS} ${_EXTERN_LIBRARIES})

