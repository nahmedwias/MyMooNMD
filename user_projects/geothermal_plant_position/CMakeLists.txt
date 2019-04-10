set(gppo_dir "${PROJECT_SOURCE_DIR}/user_projects/geothermal_plant_position")
include_directories(${gppo_dir})

# list of sources for geothermal plant position optimization
list(APPEND GTPSRC "${gppo_dir}/GeothermalPlantsPositionOptimization.cpp" )
#list(APPEND GTPSRC "${gppo_dir}/TCD2D_Temperature.cpp" )
list(APPEND GTPSRC "${gppo_dir}/TCD_Temperature.cpp" )
#list(APPEND GTPSRC "${gppo_dir}/NSE2D_GPPO.cpp" )
list(APPEND GTPSRC "${gppo_dir}/NSE_GPPO.cpp" )

#list(APPEND GPPOSRC_3D "${gppo_dir}/GeothermalPlantsPositionOptimization3D.cpp" )
list(APPEND GPPOSRC_3D "${gppo_dir}/GeothermalPlantsPositionOptimization.cpp" )
#list(APPEND GPPOSRC_3D "${gppo_dir}/CD3D_Temperature.cpp" )
list(APPEND GPPOSRC_3D "${gppo_dir}/TCD_Temperature.cpp" )
#list(APPEND GPPOSRC_3D "${gppo_dir}/NSE3D_GPPO.cpp" )
list(APPEND GPPOSRC_3D "${gppo_dir}/NSE_GPPO.cpp" )

# ONLY IN SEQ CASE
if(${PARMOON_PARALLEL_TYPE} STREQUAL SEQUENTIAL OR ${PARMOON_PARALLEL_TYPE} STREQUAL OMPONLY)

add_library(GTP_LIB_2D STATIC ${GTPSRC})
target_compile_definitions(GTP_LIB_2D PUBLIC -D__2D__)
set_target_properties(GTP_LIB_2D PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

add_library(GPPO_LIB_3D STATIC ${GPPOSRC_3D})
target_compile_definitions(GPPO_LIB_3D PUBLIC -D__3D__)
set_target_properties(GPPO_LIB_3D PROPERTIES POSITION_INDEPENDENT_CODE TRUE)


set(all_gppo_libs GTP_LIB_2D parmoon_2d_${PARMOON_PARALLEL_TYPE} 
    ${PARMOON_EXTERNAL_LIBRARIES} ${NLOPT_LIBRARIES})

set(all_gppo_libs_3D GPPO_LIB_3D parmoon_3d_${PARMOON_PARALLEL_TYPE} 
    ${PARMOON_EXTERNAL_LIBRARIES} ${NLOPT_LIBRARIES})

add_executable(gppo2d ${gppo_dir}/GPPO2D_ParMooN.cpp)
set_target_properties(gppo2d PROPERTIES OUTPUT_NAME 
parMooN_gppo2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(gppo2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(gppo2d ${all_gppo_libs})
list(APPEND parmoon_main_programs gppo2d)

add_executable(gppo3d ${gppo_dir}/GPPO3D_ParMooN.cpp)
set_target_properties(gppo3d PROPERTIES OUTPUT_NAME 
parMooN_gppo3d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(gppo3d PUBLIC -D__3D__)
# Link in the required libraries.
target_link_libraries(gppo3d ${all_gppo_libs_3D})
list(APPEND parmoon_main_programs gppo3d)

set(parmoon_main_programs ${parmoon_main_programs} PARENT_SCOPE)

endif(${PARMOON_PARALLEL_TYPE} STREQUAL SEQUENTIAL OR ${PARMOON_PARALLEL_TYPE} STREQUAL OMPONLY)

