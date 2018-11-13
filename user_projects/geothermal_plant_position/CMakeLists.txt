set(gppo_dir "${PROJECT_SOURCE_DIR}/user_projects/geothermal_plant_position")
include_directories(${gppo_dir})

# list of sources for boundary controlled optimization
list(APPEND GTPSRC "${gppo_dir}/GeothermalPlantsPositionOptimization.cpp" )
list(APPEND GTPSRC "${gppo_dir}/CD2D_Temperature.cpp" )

add_library(GTP_LIB_2D STATIC ${GTPSRC})
target_compile_definitions(GTP_LIB_2D PUBLIC -D__2D__)
set_target_properties(GTP_LIB_2D PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

set(all_gppo_libs GTP_LIB_2D parmoon_2d_${PARMOON_PARALLEL_TYPE} 
    ${_EXTERN_LIBRARIES} ${NLOPT_LIBRARIES})

if(${PARMOON_PARALLEL_TYPE} STREQUAL SEQUENTIAL OR ${PARMOON_PARALLEL_TYPE} STREQUAL OMPONLY)

  add_executable(gppo2d ${gppo_dir}/GTP2D_ParMooN.cpp)
  set_target_properties(gppo2d PROPERTIES OUTPUT_NAME 
  parMooN_gppo2d_${PARMOON_PARALLEL_TYPE})
  # Add the correct geometry flag to the target.
  target_compile_definitions(gppo2d PUBLIC -D__2D__)
  # Link in the required libraries.
  target_link_libraries(gppo2d ${all_gppo_libs})
  set(parmoon_main_programs ${parmoon_main_programs} PARENT_SCOPE)
  
endif(${PARMOON_PARALLEL_TYPE} STREQUAL SEQUENTIAL OR ${PARMOON_PARALLEL_TYPE} STREQUAL OMPONLY)

