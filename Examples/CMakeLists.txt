

include_directories("${CMAKE_SOURCE_DIR}/Examples")

# does this work?
#file(GLOB_RECURSE EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/Examples/*${GEO}.C")

if("${GEO}" STREQUAL "2D")
  file(GLOB_RECURSE EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/Examples/*2D.C")
elseif("${GEO}" STREQUAL "3D")
  file(GLOB_RECURSE EXAMPLE_SOURCES "${PROJECT_SOURCE_DIR}/Examples/*3D.C")
endif("${GEO}" STREQUAL "2D")

add_library(examples STATIC ${EXAMPLE_SOURCES})

