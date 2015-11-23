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
#    target_link_libraries(project PROJECT_LIB_2D ${_PARMOON_LIBRARIES_2D} ${_EXTERN_LIBRARIES})
# 
# Then cmake will create a makefile with a target 'project', so typing 
#    make project
# will compile your project. The executable will be called parMooN_project_X
# where 'X' is the build type, e.g. SEQUENTIAL.
include_directories("${PROJECT_SOURCE_DIR}/user_projects/include")

function(add_memcheck_test name binary)
    set(memcheck_command "valgrind --leak-check=yes --num-callers=4 --show-possibly-lost=no --undef-value-errors=yes --track-origins=yes")
    separate_arguments(memcheck_command)
    # add_test(${name} ${binary} ${ARGN})
    add_test(memcheck_${name} ${memcheck_command} ./${binary} ${ARGN})
endfunction(add_memcheck_test)

function(set_memcheck_test_properties name)
    set_tests_properties(${name} ${ARGN})
    set_tests_properties(memcheck_${name} ${ARGN})
endfunction(set_memcheck_test_properties)

add_executable(cd2d_moritz ${PROJECT_SOURCE_DIR}/user_projects/CD2D_ParMooN.C)
set_target_properties(cd2d_moritz PROPERTIES OUTPUT_NAME cd2d_moritz)
target_compile_definitions(cd2d_moritz PUBLIC -D__2D__)
set_target_properties(cd2d_moritz PROPERTIES LINK_FLAGS "-v")
target_link_libraries(cd2d_moritz parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})


add_executable(nse2d_moritz ${PROJECT_SOURCE_DIR}/user_projects/NSE2D_ParMooN.C)
set_target_properties(nse2d_moritz PROPERTIES OUTPUT_NAME nse2d_moritz)
target_compile_definitions(nse2d_moritz PUBLIC -D__2D__)
set_target_properties(nse2d_moritz PROPERTIES LINK_FLAGS "-v")
target_link_libraries(nse2d_moritz parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

add_executable(localized_refinement localized/LocalizedDomainRefinement.C localized/LocalizedErrorIndicator.cpp)
set_target_properties(localized_refinement PROPERTIES OUTPUT_NAME localized_refinement)
target_compile_definitions(localized_refinement PUBLIC -D__2D__)
target_link_libraries(localized_refinement parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})

add_memcheck_test("valgrind-nse2d" nse2d_moritz)
# add_custom_target("valgrind-nse2d_moritz" COMMAND valgrind -q --leak-check=no --num-callers=4 --show-possibly-lost=no --undef-value-errors=yes --track-origins=yes ${CMAKE_BINARY_DIR}/user_projects/nse2d_moritz)