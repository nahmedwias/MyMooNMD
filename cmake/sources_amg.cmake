

include_directories("${CMAKE_SOURCE_DIR}/include/AMG")

file(GLOB_RECURSE AMG_SOURCES "${PROJECT_SOURCE_DIR}/src/AMG/*.c")

list(APPEND PARMOON_SOURCES_2D ${AMG_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${AMG_SOURCES})

