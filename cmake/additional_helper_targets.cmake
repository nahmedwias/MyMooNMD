###############################################################################
# Additional target to perform clang-format run
# Check if clang-format is available
find_program(clang_format
             NAMES clang-format clang-format-4.0 clang-format-3.6)
if(clang_format-NOTFOUND)
  add_custom_target(clang-format
                  COMMENT "It seems cmake could not find a clang-format
                  executable. Therefore no formatting is done!" VERBATIM)
else()
  # Get all project files
  file(GLOB_RECURSE ALL_SOURCE_FILES RELATIVE 
     ${CMAKE_CURRENT_SOURCE_DIR}/documentation *.cpp *.hpp)
  add_custom_target(
    clang-format
    COMMAND ${clang_format} -style=file -i ${ALL_SOURCE_FILES}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/documentation)
endif(clang_format-NOTFOUND)

###############################################################################
# Additional target to perform cppcheck (Note: many IDEs integrate this nicely)
# Check if cpp_check is available
find_program(cpp_check NAMES cppcheck)
if(NOT cpp_check)
  add_custom_target(cppcheck
                    COMMENT "It seems cmake could not find a cpp_check 
                    executable. Therefore no check was done!" VERBATIM)
else()
  set (cppcheck-outfile cppcheck.txt)
  # the following createes a file 'compile_commands.json' which makes using
  # cppcheck especially simple
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  add_custom_target(
    cppcheck
    COMMAND ${cpp_check} --enable=all --project=compile_commands.json 
                         --output-file=${cppcheck-outfile} --std=c++11
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "A file ${cppcheck-outfile} will be created in the build directory!"
    VERBATIM)
  add_custom_command(
    TARGET cppcheck POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan
    "The results of cppcheck are in the file: "
    "${CMAKE_CURRENT_BINARY_DIR}/${cppcheck-outfile}")
 endif(NOT cpp_check)
 
