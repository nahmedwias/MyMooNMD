###############################################################################
# getting information from the repository and write into a file which can be 
# included during compilation.

set(PARMOON_HG_REVISION "unknown")
set(PARMOON_HG_BRANCH "unknown")
set(PARMOON_LOCAL_CHANGES "unknown")

find_package(Hg)
if(HG_FOUND)
  # find the revision number and find out if there are local changes
  execute_process(COMMAND ${HG_EXECUTABLE} id -i
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                  OUTPUT_VARIABLE _info
                  RESULT_VARIABLE _result
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(${_result} EQUAL 0)
    # there is a '+' at the end of the revision number in case there are local
    # changes
    # remove the + if it exists to get the revision number
    string(REPLACE "+" "" PARMOON_HG_REVISION ${_info})
    # if there are local changes, set the variable correspondingly
    if(${_info} MATCHES "\\+$")
      set(PARMOON_LOCAL_CHANGES "true")
    else()
      set(PARMOON_LOCAL_CHANGES "false")
    endif()
  endif()
  
  # find the branch name
  execute_process(COMMAND ${HG_EXECUTABLE} branch
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                  OUTPUT_VARIABLE _info
                  RESULT_VARIABLE _result
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(${_result} EQUAL 0)
    set(PARMOON_HG_BRANCH ${_info})
  endif()
endif()

configure_file(${CMAKE_SOURCE_DIR}/cmake/ParMooN_repository_info.h.in
               ${CMAKE_SOURCE_DIR}/include/General/ParMooN_repository_info.h
               @ONLY)
