
set(stoda_dir "${PROJECT_SOURCE_DIR}/Main_Users/StokesDarcy2D")

include_directories(${stoda_dir})
include_directories("${stoda_dir}/Examples" )
include_directories("${stoda_dir}/solvers" )

# list of sources for Stokes-Darcy (StoDa) 2D
list(APPEND STODA_SOURCES_2D "${stoda_dir}/Examples/Example_StokesDarcy2D.C" )
list(APPEND STODA_SOURCES_2D "${stoda_dir}/auxiliaryFunctions.C" )
list(APPEND STODA_SOURCES_2D "${stoda_dir}/DarcyProblem.C" )
list(APPEND STODA_SOURCES_2D "${stoda_dir}/InterfaceFunction.C" )
list(APPEND STODA_SOURCES_2D "${stoda_dir}/StokesDarcy2D.C" )
list(APPEND STODA_SOURCES_2D "${stoda_dir}/StokesProblem.C" )
list(APPEND STODA_SOURCES_2D "${stoda_dir}/solvers/UmfpackSolverWrapper.C" )


add_library(STODA_LIB_2D STATIC ${STODA_SOURCES_2D})
target_compile_definitions(STODA_LIB_2D PUBLIC -D__2D__)

add_executable(sd2d ${stoda_dir}/main.C)
set_target_properties(sd2d PROPERTIES OUTPUT_NAME 
parMooN_sd2d_${PARMOON_PARALLEL_TYPE})
# Add the correct geometry flag to the target.
target_compile_definitions(sd2d PUBLIC -D__2D__)
# Link in the required libraries.
target_link_libraries(sd2d STODA_LIB_2D ${_PARMOON_LIBRARIES_2D} 
${_EXTERN_LIBRARIES})
