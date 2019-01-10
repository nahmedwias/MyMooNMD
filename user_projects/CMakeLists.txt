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

set(non_newt "${PROJECT_SOURCE_DIR}/user_projects")

include_directories(${non_newt})
#source files
list(APPEND NON_NEWT2D "${non_newt}/CoupledNavierStokesStress.cpp")
list(APPEND NON_NEWT2D "${non_newt}/Coupled_NS_Stress_Local_routines.cpp")

add_library(NON_NEWT_LIB_2D STATIC ${NON_NEWT2D})
target_compile_definitions(NON_NEWT_LIB_2D PUBLIC -D__2D__)

### TNSE2D ###
# Standard executable target.
add_executable(tnse2d_coupled ${non_newt}/CoupledNSE2D_Stress_ParMooN.cpp)
set_target_properties(tnse2d_coupled PROPERTIES OUTPUT_NAME parMooN_tnse2d_coupled_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tnse2d_coupled PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(tnse2d_coupled NON_NEWT_LIB_2D parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})