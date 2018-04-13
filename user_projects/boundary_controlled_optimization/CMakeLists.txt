
set(bco_dir "${PROJECT_SOURCE_DIR}/user_projects/boundary_controlled_optimization")

# include 'nlopt' as an external library, 
# see https://nlopt.readthedocs.io/en/latest/
set(NLOPT_DIR CACHE STRING "The path to the directory to find nlopt in.")
message(STATUS "NLOPT_DIR ${NLOPT_DIR}")
if(NLOPT_DIR STREQUAL "")
  message(WARNING "Please provide a path to where 'nlopt' is installed")
  message(FATAL_ERROR "Missing NLOPT_DIR")
endif(NLOPT_DIR STREQUAL "")
include("${NLOPT_DIR}/lib/cmake/nlopt/NLoptConfig.cmake")
include_directories(${NLOPT_INCLUDE_DIRS})


include_directories(${bco_dir})

# list of sources for boundary controlled optimization
list(APPEND BCOSRC "${bco_dir}/BoundaryControlledOptimization.cpp" )
list(APPEND BCOSRC "${bco_dir}/NSE2D_Adjoint.cpp" )

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


#add_executable(ParMooN_test_bco2d ${bco_dir}/regression_test.C)
#target_link_libraries(ParMooN_test_bco2d ${all_bco_libs})
#add_test(stokes_darcy ParMooN_test_bco2d)
#list(APPEND parmoon_tests ParMooN_test_bco2d)
#set(parmoon_tests ${parmoon_tests} PARENT_SCOPE)


