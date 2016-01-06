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

# Shuffle everything into a static project library.
list(APPEND _PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/CoupledCDR_2D.C")
list(APPEND _PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/Coupled_Time_CDR_2D.C")
list(APPEND _PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/ReactionCoupling.C")
list(APPEND _PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/Example_CoupledCDR2D.C")
list(APPEND _PROJECT_SOURCES_2D "${PROJECT_SOURCE_DIR}/user_projects/src/FEM_FCT_Copy.C")
add_library(PROJECT_LIB_2D STATIC ${_PROJECT_SOURCES_2D})
target_compile_definitions(PROJECT_LIB_2D PUBLIC -D__2D__)

# Add an executable target from CoupledCDR.C
add_executable(ccdr2d ${PROJECT_SOURCE_DIR}/user_projects/CoupledCDR2D_ParMooN.C)
set_target_properties(ccdr2d PROPERTIES OUTPUT_NAME parMooN_ccdr2d_${PARMOON_PARALLEL_TYPE})
target_compile_definitions(ccdr2d PUBLIC -D__2D__)
target_link_libraries(ccdr2d PROJECT_LIB_2D parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})



