# unfortunately we can not use 'make test' as we would like to as it does not 
# set the dependencies correctly. That means entering 'make test' will not 
# compile the tests, even if they have changed or never been compiled before. 
# 
# Here we create a custom target 'check' which first compiles all tests and its
# dependencies if necessary and then runs the tests. All you need to do enter 
# is then 'make check'.
#
# In order for this to work we collect each executable which is a test in a 
# list of tests called 'parmoon_tests'. The check target is then a custom 
# target which calls ctest and depends on all the executables.
# See also https://cmake.org/Wiki/CMakeEmulateMakeCheck
#
# Currently "make check" runs a different set of tests, depending
# on the current parallel type. Make sure all do run.

if(${PARMOON_PARALLEL_TYPE} STREQUAL SEQUENTIAL OR ${PARMOON_PARALLEL_TYPE} STREQUAL OMPONLY)
  
  ###################################################
  ### A function with which to register tests in  ###   
  ### the "check" target for sequential build     ###
  ###################################################
  function(register_executable exec_name file_name)
    
    add_executable(${exec_name} ${file_name})
    target_link_libraries(${exec_name} parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
    set(parmoon_tests ${parmoon_tests} ${exec_name} PARENT_SCOPE)
    
  endfunction(register_executable)

  ###################################################
  ####### Tests for parallel type: SEQUENTIAL #######
  ###################################################
  
  #parameter database test
  add_executable(ParMooN_test_parameter_database parameter_database_test.c++)
  target_link_libraries(ParMooN_test_parameter_database parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(parameter_database ParMooN_test_parameter_database)
  list(APPEND parmoon_tests ParMooN_test_parameter_database)
  
  #matrix test
  add_executable(ParMooN_test_matrix matrix_test.c++)
  target_link_libraries(ParMooN_test_matrix parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(matrix ParMooN_test_matrix)
  list(APPEND parmoon_tests ParMooN_test_matrix)
  
  #block matrix test
  add_executable(ParMooN_test_block_matrix block_matrix_test.c++)
  target_link_libraries(ParMooN_test_block_matrix parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(block_matrix ParMooN_test_block_matrix)
  list(APPEND parmoon_tests ParMooN_test_block_matrix)
  
  #block fe matrix test
  add_executable(ParMooN_test_block_fe_matrix block_fe_matrix_test.c++)
  target_link_libraries(ParMooN_test_block_fe_matrix parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(block_fe_matrix ParMooN_test_block_fe_matrix)
  list(APPEND parmoon_tests ParMooN_test_block_fe_matrix)
 
  #cd2d test
  add_executable(ParMooN_test_cd2d cd2d_test.c++)
  target_link_libraries(ParMooN_test_cd2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(cd2d ParMooN_test_cd2d)
  list(APPEND parmoon_tests ParMooN_test_cd2d)
	
  #tcd2d test
  add_executable(ParMooN_test_tcd2d tcd2d_test.c++)
  target_link_libraries(ParMooN_test_tcd2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(tcd2d ParMooN_test_tcd2d)
  list(APPEND parmoon_tests ParMooN_test_tcd2d)
  
  #tcd2d with fem fct test
  add_executable(ParMooN_test_tcd2d_femfct tcd2d_femfct_test.c++)
  target_link_libraries(ParMooN_test_tcd2d_femfct parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(tcd2d_femfct ParMooN_test_tcd2d_femfct)
  list(APPEND parmoon_tests ParMooN_test_tcd2d_femfct)
	
  #darcy2d test
  add_executable(ParMooN_test_darcy2d darcy2d_test.c++)
  target_link_libraries(ParMooN_test_darcy2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(darcy2d ParMooN_test_darcy2d)
  list(APPEND parmoon_tests ParMooN_test_darcy2d)
  
  #nse2d test
  add_executable(ParMooN_test_nse2d nse2d_test.c++)
  target_link_libraries(ParMooN_test_nse2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(nse2d ParMooN_test_nse2d)
  list(APPEND parmoon_tests ParMooN_test_nse2d)
  
  #nse2d multigrid test
  add_executable(ParMooN_test_nse2d_mg nse2d_multigrid_test.c++)
  target_link_libraries(ParMooN_test_nse2d_mg parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(nse2d_mg ParMooN_test_nse2d_mg)
  list(APPEND parmoon_tests ParMooN_test_nse2d_mg)
  
  #tnse2d test (direct solver)
  add_executable(ParMooN_test_tnse2d tnse2d_test.c++)
  target_link_libraries(ParMooN_test_tnse2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(tnse2d ParMooN_test_tnse2d)
  list(APPEND parmoon_tests ParMooN_test_tnse2d)
	
  # test for hdiv elements
  add_executable(ParMooN_test_hdiv Hdiv_element_test.c++)
  target_link_libraries(ParMooN_test_hdiv parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(hdiv ParMooN_test_hdiv)
  list(APPEND parmoon_tests ParMooN_test_hdiv)
    
  #cd3d test sequential
  add_executable(ParMooN_test_cd3d_seq cd3d_test.c++)
  target_link_libraries(ParMooN_test_cd3d_seq parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(cd3d_seq ParMooN_test_cd3d_seq)
  list(APPEND parmoon_tests ParMooN_test_cd3d_seq)
  
  #cd3d solver tests
  set(_EXEC_NAME ParMooN_test_cd3d_solver_seq)
  set(_FILE_NAME cd3d_solvers_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  #umfpack
  add_test(NAME cd3d_seq_umfpack COMMAND ParMooN_test_cd3d_solver_seq umfpack)
  #richardson iteration with jacobi
  add_test(NAME cd3d_seq_jacobi COMMAND ParMooN_test_cd3d_solver_seq jacobi)
  #richardson iteration with multigrid
  add_test(NAME cd3d_seq_mg COMMAND ParMooN_test_cd3d_solver_seq multigrid)

  #nse3d solver tests
  set(_EXEC_NAME ParMooN_test_nse3d_solver_seq)
  set(_FILE_NAME nse3d_solvers_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  #umfpack
  add_test(NAME nse3d_seq_umfpack COMMAND ParMooN_test_nse3d_solver_seq umfpack)
  #multigrid
  add_test(NAME nse3d_seq_mg COMMAND ParMooN_test_nse3d_solver_seq multigrid)
  #lsc
  add_test(NAME nse3d_seq_lsc COMMAND ParMooN_test_nse3d_solver_seq lsc)

  
  # direct solver test
  add_executable(ParMooN_test_direct_solver direct_solver_test.c++)
  target_link_libraries(ParMooN_test_direct_solver parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(direct_solver ParMooN_test_direct_solver)
  list(APPEND parmoon_tests ParMooN_test_direct_solver)
  
  # iterative solver test
  add_executable(ParMooN_test_iterative_solver iterative_solver_test.c++)
  target_link_libraries(ParMooN_test_iterative_solver parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(iterative_solver ParMooN_test_iterative_solver)
  list(APPEND parmoon_tests ParMooN_test_iterative_solver)
	
elseif(${PARMOON_PARALLEL_TYPE} STREQUAL MPI)	
  
  #####################################################
  ### A function with which to register executables ###   
  ### in the "check" target for mpi build              ###
  #####################################################
  function(register_executable exec_name file_name)
    
    add_executable(${exec_name} ${file_name})
    target_link_libraries(${exec_name} parmoon_3d_MPI ${_EXTERN_LIBRARIES})
    set(parmoon_tests ${parmoon_tests} ${exec_name} PARENT_SCOPE)
    
  endfunction(register_executable)
  
  
  ###################################################
  ####### Tests for parallel type: MPI #######
  ###################################################
  
  # cd3d test mpi (run with 4 processes)
  add_executable(ParMooN_test_cd3d_mpi cd3d_test.c++)
  target_link_libraries(ParMooN_test_cd3d_mpi parmoon_3d_MPI ${_EXTERN_LIBRARIES})
  add_test(NAME cd3d_mpi COMMAND mpirun -n 4 ParMooN_test_cd3d_mpi )
  list(APPEND parmoon_tests ParMooN_test_cd3d_mpi)
  
  # 1st mumps solver test mpi (run with 2 processes)
  add_executable(ParMooN_test_mumps_1_mpi mumps_test_1.c++)
  target_link_libraries(ParMooN_test_mumps_1_mpi ${_EXTERN_LIBRARIES})
  add_test(NAME mumps_1_mpi COMMAND mpirun -n 2 ParMooN_test_mumps_1_mpi )
  list(APPEND parmoon_tests ParMooN_test_mumps_1_mpi)
  
  # 2nd mumps solver test mpi (run with 4 processes)
  add_executable(ParMooN_test_mumps_2_mpi mumps_test_2.c++)
  target_link_libraries(ParMooN_test_mumps_2_mpi ${_EXTERN_LIBRARIES})
  add_test(NAME mumps_2_mpi COMMAND mpirun -n 4 ParMooN_test_mumps_2_mpi )
  list(APPEND parmoon_tests ParMooN_test_mumps_2_mpi)
  
  # mumps wrapper test mpi (run with 2 processes)
  add_executable(ParMooN_test_mumps_wrapper mumps_wrapper_test.c++)
  target_link_libraries(ParMooN_test_mumps_wrapper parmoon_3d_MPI ${_EXTERN_LIBRARIES})
  add_test(NAME mumps_wrapper COMMAND mpirun -n 2 ParMooN_test_mumps_wrapper )
  list(APPEND parmoon_tests ParMooN_test_mumps_wrapper)
  
  #cd3d solver tests
  set(_EXEC_NAME ParMooN_test_cd3d_solver_mpi)
  set(_FILE_NAME cd3d_solvers_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  # mumps
  add_test(NAME cd3d_mpi_mumps COMMAND mpirun -n 3 ParMooN_test_cd3d_solver_mpi mumps)
  #jacobi - CB: THIS DOES STILL FAIL - FIX THE ITERATIVE MPI!
  add_test(NAME cd3d_mpi_jacobi COMMAND mpirun -n 3 ParMooN_test_cd3d_solver_mpi jacobi)
  #multigrid - CB: THIS DOES STILL FAIL - FIX THE ITERATIVE MPI!
  add_test(NAME cd3d_mpi_mg COMMAND mpirun -n 3 ParMooN_test_cd3d_solver_mpi multigrid)
  
  #nse3d solver tests
  set(_EXEC_NAME ParMooN_test_nse3d_solver_mpi)
  set(_FILE_NAME nse3d_solvers_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  #mumps
  add_test(NAME nse3d_mpi_mumps COMMAND mpirun -n 3 ParMooN_test_nse3d_solver_mpi mumps)


endif(${PARMOON_PARALLEL_TYPE} STREQUAL SEQUENTIAL OR ${PARMOON_PARALLEL_TYPE} STREQUAL OMPONLY)
	
####################################################
#                THE CHECK TARGET                  #
#   making this target results in building and     #
#   running all tests corresponding to the current #
#   parallel type. Use it excessively!             #
####################################################
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure DEPENDS ${parmoon_tests})


