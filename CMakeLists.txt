
cmake_minimum_required (VERSION 2.8)

project (ParMooN)

# The version number.
set (ParMooN_VERSION_MAJOR 1)
set (ParMooN_VERSION_MINOR 0)


include(${CMAKE_SOURCE_DIR}/config.cmake)


if (NOT CONFIGURED_ONCE)
    set(CMAKE_CXX_FLAGS "${warnings}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_C_FLAGS   "${warnings}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif()

include_directories("${CMAKE_SOURCE_DIR}/include/AMG")
include_directories("${CMAKE_SOURCE_DIR}/include/Refinement")
include_directories("${CMAKE_SOURCE_DIR}/include/QuadFormulas")
include_directories("${CMAKE_SOURCE_DIR}/include/Geometry")
include_directories("${CMAKE_SOURCE_DIR}/include/General")
include_directories("${CMAKE_SOURCE_DIR}/include/FE")
include_directories("${CMAKE_SOURCE_DIR}/include/System")
include_directories("${CMAKE_SOURCE_DIR}/include/PBE")
include_directories("${CMAKE_SOURCE_DIR}/Examples")


add_subdirectory(src/General)  
add_subdirectory(src/Geometry)          

# link_directories(${CMAKE_BINARY_DIR}/src/General/obj1D)

# MESSAGE("===================================="  )        
# MESSAGE("ParMooN using the compiler: '${CMAKE_CXX_COMPILER_ID}'"  )  
# MESSAGE("C compiler  : '${CMAKE_C_COMPILER}'"  )  
# MESSAGE("C++ compiler: '${CMAKE_CXX_COMPILER}'"  )  
# MESSAGE(" C++ flags  : '${PARALLEL_TYPE}'" )
# MESSAGE("===================================="  ) 



message("UMFPACK_USED ='${UMFPACK_USED}'")
message("MUMPS_USED ='${MUMPS_USED}'")

# FOREACH(var ${CMAKE_SYSTEM_PREFIX_PATH})
#     MESSAGE(STATUS "CMAKE_SYSTEM_PREFIX_PATH = " ${var})
# ENDFOREACH(var)
# MESSAGE(STATUS)

#         

# if (USE_MUMPS)
#   include_directories ("${PROJECT_SOURCE_DIR}/EXT_LIB/MUMPS/")
#   add_subdirectory (EXT_LIB/MUMPS/)
#   set (EXTRA_LIBS ${EXTRA_LIBS} EXT_LIB/MUMPS/)
# endif (USE_MUMPS)        
        
add_executable(parmoon 2DPrograms/CD2D_ParMooN.C)

target_link_libraries(parmoon PRIVATE "gen" "geo" )

#  add_library(mylib src1.c src2.c)
#  add_executable(myexe main1.c)
#   set_property(
#    DIRECTORY
#    PROPERTY COMPILE_DEFINITIONS A AV=1
#    )
#  set_property(
#    TARGET mylib
#    PROPERTY COMPILE_DEFINITIONS B BV=2
#    )
#  set_property(
#    SOURCE src1.c
#    PROPERTY COMPILE_DEFINITIONS C CV=3
#    )
#  set_property(
#    SOURCE src1.c
#    APPEND PROPERTY COMPILE_DEFINITIONS D DV=4
#    )
# will build the source files with these definitions:
#  src1.c:   -DA -DAV=1 -DB -DBV=2 -DC -DCV=3
#  src2.c:   -DA -DAV=1 -DB -DBV=2
#  main2.c:  -DA -DAV=1

