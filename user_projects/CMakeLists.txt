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

# Here begins the code which enables usage of our project specific POD/ROM CDR stuff.
# Swetlana code copied and modified by Naveed 

set(pod_rom "${PROJECT_SOURCE_DIR}/user_projects")

include_directories(${pod_rom})
include_directories("${pod_rom}/inc")
include_directories("${pod_rom}/src")

# We Need Boost For POD
find_package(Boost 1.47.0 REQUIRED COMPONENTS regex filesystem system program_options)
IF(Boost_FOUND)
    MESSAGE("Boost includes located in: ${Boost_INCLUDE_DIRS}")
    MESSAGE("Boost libraries located in: ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIRS})
    list(APPEND _EXTERN_LIBRARIES ${Boost_LIBRARIES})
ENDIF(Boost_FOUND)

#set( ENV{BLA_VENDOR} "ATLAS" )
find_package( BLAS REQUIRED )

#source files 
list(APPEND ROM_SOURCES_2D "${pod_rom}/src/SNAPS.C")
list(APPEND ROM_SOURCES_2D "${pod_rom}/src/POD.C")
list(APPEND ROM_SOURCES_2D "${pod_rom}/src/Temprary.cpp")
list(APPEND ROM_SOURCES_2D "${pod_rom}/src/TimeConvectionDiffusionPOD.C")
list(APPEND ROM_SOURCES_2D "${pod_rom}/src/TimeConvectionDiffusionROM.C")


add_library(POD_LIB_2D STATIC ${ROM_SOURCES_2D})
target_compile_definitions(POD_LIB_2D PUBLIC -D__2D__)

# Add an executable target from TCDR2D_ParMooN.C
add_executable(tcdr2d ${pod_rom}/TCDR2D_ParMooN.C)
set_target_properties(tcdr2d PROPERTIES OUTPUT_NAME 
parMooN_tcdr2d_${PARMOON_PARALLEL_TYPE})
target_compile_definitions(tcdr2d PUBLIC -D__2D__)
target_link_libraries(tcdr2d POD_LIB_2D parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

# Add an executable target from TCD
add_executable(tcdr2d_pod ${pod_rom}/TCDR2D_POD_ParMooN.C)
set_target_properties(tcdr2d_pod PROPERTIES OUTPUT_NAME 
parMooN_tcdr2d_pod_${PARMOON_PARALLEL_TYPE})
target_compile_definitions(tcdr2d_pod PUBLIC -D__2D__)
target_link_libraries(tcdr2d_pod POD_LIB_2D parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
