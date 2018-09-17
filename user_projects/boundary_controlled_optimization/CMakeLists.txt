
set(bco_dir "${PROJECT_SOURCE_DIR}/user_projects/boundary_controlled_optimization")
include_directories(${bco_dir})

# list of sources for boundary controlled optimization
list(APPEND BCOSRC "${bco_dir}/BoundaryControlledOptimization.cpp" )
list(APPEND BCOSRC "${bco_dir}/NSE2D_Adjoint.cpp" )
list(APPEND BCOSRC "${bco_dir}/Time_BoundaryControlledOptimization.cpp" )
list(APPEND BCOSRC "${bco_dir}/Time_NSE2D_Adjoint.cpp" )

add_library(BCO_LIB_2D STATIC ${BCOSRC})
target_compile_definitions(BCO_LIB_2D PUBLIC -D__2D__)
set_target_properties(BCO_LIB_2D PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

set(all_bco_libs BCO_LIB_2D parmoon_2d_${PARMOON_PARALLEL_TYPE} 
    ${_EXTERN_LIBRARIES} ${NLOPT_LIBRARIES})

add_executable(bco2d ${bco_dir}/main.cpp)
set_target_properties(bco2d PROPERTIES OUTPUT_NAME 
parMooN_bco2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(bco2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(bco2d ${all_bco_libs})

add_executable(tbco2d ${bco_dir}/TBCO2D_ParMooN.cpp)
set_target_properties(tbco2d PROPERTIES OUTPUT_NAME 
parMooN_tbco2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(tbco2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(tbco2d ${all_bco_libs})


#add_executable(ParMooN_test_bco2d ${bco_dir}/regression_test.C)
#target_link_libraries(ParMooN_test_bco2d ${all_bco_libs})
#add_test(stokes_darcy ParMooN_test_bco2d)
#list(APPEND parmoon_tests ParMooN_test_bco2d)
#set(parmoon_tests ${parmoon_tests} PARENT_SCOPE)


