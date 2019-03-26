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
  ### the "check" target for 3D sequential build  ###
  ###################################################
  function(register_executable exec_name file_name)
    
    add_executable(${exec_name} ${file_name})
    target_link_libraries(${exec_name} parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
    set(parmoon_tests ${parmoon_tests} ${exec_name} PARENT_SCOPE)
    
  endfunction(register_executable)
  
  ###################################################
  ### A function with which to register tests in  ###   
  ### the "check" target for 2D sequential build  ###
  ###################################################
  function(register_2D_executable exec_name file_name)
    
    add_executable(${exec_name} ${file_name})
    target_link_libraries(${exec_name} parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
    set(parmoon_tests ${parmoon_tests} ${exec_name} PARENT_SCOPE)
    
  endfunction(register_2D_executable)

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
  
  #dense matrix test
  add_executable(ParMooN_test_dense_matrix dense_matrix_test.c++)
  target_link_libraries(ParMooN_test_dense_matrix parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(dense_matrix ParMooN_test_dense_matrix)
  list(APPEND parmoon_tests ParMooN_test_dense_matrix)
  
  # analytical function test
  add_executable(ParMooN_test_analytical_function analytical_function_test.c++)
  target_link_libraries(ParMooN_test_analytical_function 
                        parmoon_2d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(analytical_function ParMooN_test_analytical_function)
  list(APPEND parmoon_tests ParMooN_test_analytical_function)
  
  # boundary condition test (2d and 3d)
  add_executable(ParMooN_test_boundary_condition2d boundary_condition_test.c++)
  target_link_libraries(ParMooN_test_boundary_condition2d
                        parmoon_2d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(boundary_condition2d ParMooN_test_boundary_condition2d)
  list(APPEND parmoon_tests ParMooN_test_boundary_condition2d)
  add_executable(ParMooN_test_boundary_condition3d boundary_condition_test.c++)
  target_link_libraries(ParMooN_test_boundary_condition3d
                        parmoon_3d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(boundary_condition3d ParMooN_test_boundary_condition3d)
  list(APPEND parmoon_tests ParMooN_test_boundary_condition3d)
  
  # boundary data test (2d and 3d)
  add_executable(ParMooN_test_boundary_data2d boundary_data_test.c++)
  target_link_libraries(ParMooN_test_boundary_data2d
                        parmoon_2d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(boundary_data2d ParMooN_test_boundary_data2d)
  list(APPEND parmoon_tests ParMooN_test_boundary_data2d)
  add_executable(ParMooN_test_boundary_data3d boundary_data_test.c++)
  target_link_libraries(ParMooN_test_boundary_data3d
                        parmoon_3d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(boundary_data3d ParMooN_test_boundary_data3d)
  list(APPEND parmoon_tests ParMooN_test_boundary_data3d)
  
  add_executable(ParMooN_test_scott_vogelius_2d scott_vogelius_test.c++)
  target_link_libraries(ParMooN_test_scott_vogelius_2d
                        parmoon_2d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(scott_vogelius_2d ParMooN_test_scott_vogelius_2d)
  list(APPEND parmoon_tests ParMooN_test_scott_vogelius_2d)
  add_executable(ParMooN_test_scott_vogelius_3d scott_vogelius_test.c++)
  target_link_libraries(ParMooN_test_scott_vogelius_3d
                        parmoon_3d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(scott_vogelius_3d ParMooN_test_scott_vogelius_3d)
  list(APPEND parmoon_tests ParMooN_test_scott_vogelius_3d)
  
  # function evaluation test
  add_executable(ParMooN_test_function_evaluation function_evaluation_test.c++)
  target_link_libraries(ParMooN_test_function_evaluation
                        parmoon_2d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(function_evaluation ParMooN_test_function_evaluation)
  list(APPEND parmoon_tests ParMooN_test_function_evaluation)
  
  # read/write fe function
  add_executable(ParMooN_test_read_write_fe_function read_write_fe_function_test.c++)
  target_link_libraries(ParMooN_test_read_write_fe_function
                        parmoon_3d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(read_write_fe_function ParMooN_test_read_write_fe_function)
  list(APPEND parmoon_tests ParMooN_test_read_write_fe_function)
  
  #cd2d test
  add_executable(ParMooN_test_cd2d cd2d_test.c++)
  target_link_libraries(ParMooN_test_cd2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(cd2d ParMooN_test_cd2d)
  list(APPEND parmoon_tests ParMooN_test_cd2d)
  
  #cd2d algebraic flux correction(afc) test
  add_executable(ParMooN_test_cd2d_afc cd2d_afc_test.c++)
  target_link_libraries(ParMooN_test_cd2d_afc parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(cd2d_afc ParMooN_test_cd2d_afc)
  list(APPEND parmoon_tests ParMooN_test_cd2d_afc)
	
  #cd2d test
  add_executable(ParMooN_test_cd2d_adaptive cd2d_adaptive_test.c++)
  target_link_libraries(ParMooN_test_cd2d_adaptive parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(cd2d_adaptive ParMooN_test_cd2d_adaptive)
  list(APPEND parmoon_tests ParMooN_test_cd2d_adaptive)

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
  add_executable(ParMooN_test_darcy2d darcy_test.c++)
  target_link_libraries(ParMooN_test_darcy2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(darcy2d ParMooN_test_darcy2d)
  list(APPEND parmoon_tests ParMooN_test_darcy2d)
  add_executable(ParMooN_test_darcy3d darcy_test.c++)
  target_link_libraries(ParMooN_test_darcy3d parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(darcy3d ParMooN_test_darcy3d)
  list(APPEND parmoon_tests ParMooN_test_darcy3d)

  #brinkman2d test
#   add_executable(ParMooN_test_brinkman2d brinkman2d_test.c++)
#   target_link_libraries(ParMooN_test_brinkman2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
#   add_test(brinkman2d ParMooN_test_brinkman2d)
#   list(APPEND parmoon_tests ParMooN_test_brinkman2d)

  #nse2d test
  add_executable(ParMooN_test_nse2d nse2d_test.c++)
  target_link_libraries(ParMooN_test_nse2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(nse2d ParMooN_test_nse2d)
  add_test(NAME nse2d_all COMMAND ParMooN_test_nse2d testall)
  list(APPEND parmoon_tests ParMooN_test_nse2d)
  
  #stokes equal order test
  add_executable(ParMooN_test_stokes_equal_order2d stokes_equal_order_test.c++)
  target_link_libraries(ParMooN_test_stokes_equal_order2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(stokes_equal_order2d ParMooN_test_stokes_equal_order2d)
  list(APPEND parmoon_tests ParMooN_test_stokes_equal_order2d)
  add_executable(ParMooN_test_stokes_equal_order3d stokes_equal_order_test.c++)
  target_link_libraries(ParMooN_test_stokes_equal_order3d parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(stokes_equal_order3d ParMooN_test_stokes_equal_order3d)
  list(APPEND parmoon_tests ParMooN_test_stokes_equal_order3d)
  
  #lsc / (bc-lsc!) test
  add_executable(ParMooN_test_lsc lsc_test.c++)
  target_link_libraries(ParMooN_test_lsc parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(lsc2d ParMooN_test_lsc)
  list(APPEND parmoon_tests ParMooN_test_lsc)
   
  #tnse2d test (umfpack and lsc, multigrid's parameters are not set up in the test yet)
  add_executable(ParMooN_test_tnse2d tnse2d_test.c++)
  target_link_libraries(ParMooN_test_tnse2d parmoon_2d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})  
  #umfpack
  add_test(NAME tnse2d_umfpack COMMAND ParMooN_test_tnse2d umfpack)
  #lsc 
  add_test(NAME tnse2d_lsc COMMAND ParMooN_test_tnse2d lsc)
  list(APPEND parmoon_tests ParMooN_test_tnse2d)
	
  # test for hdiv elements
  add_executable(ParMooN_test_hdiv_2d Hdiv_element_test.c++)
  target_link_libraries(ParMooN_test_hdiv_2d parmoon_2d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(hdiv_2d ParMooN_test_hdiv_2d)
  list(APPEND parmoon_tests ParMooN_test_hdiv_2d)
  add_executable(ParMooN_test_hdiv_3d Hdiv_element_test.c++)
  target_link_libraries(ParMooN_test_hdiv_3d parmoon_3d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(hdiv_3d ParMooN_test_hdiv_3d)
  list(APPEND parmoon_tests ParMooN_test_hdiv_3d)
  
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
  # petsc
  add_test(NAME cd3d_seq_petsc COMMAND ParMooN_test_cd3d_solver_seq petsc)
  
  #cd3d algebraic flux correction(afc) test
  add_executable(ParMooN_test_cd3d_afc cd3d_afc_test.c++)
  target_link_libraries(ParMooN_test_cd3d_afc parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(cd3d_afc ParMooN_test_cd3d_afc)
  list(APPEND parmoon_tests ParMooN_test_cd3d_afc)
  
  #tcd3d solver tests
  set(_EXEC_NAME ParMooN_test_tcd3d_solver_seq)
  set(_FILE_NAME tcd3d_solvers_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  #umfpack
  add_test(NAME tcd3d_seq_umfpack COMMAND ParMooN_test_tcd3d_solver_seq umfpack)
  #richardson iteration with jacobi
  add_test(NAME tcd3d_seq_jacobi COMMAND ParMooN_test_tcd3d_solver_seq jacobi)
  #richardson iteration with multigrid
  add_test(NAME tcd3d_seq_mg COMMAND ParMooN_test_tcd3d_solver_seq multigrid)
  # petsc
  add_test(NAME tcd3d_seq_petsc COMMAND ParMooN_test_tcd3d_solver_seq petsc)
  
  #tcd3d algebraic flux correction tests
  set(_EXEC_NAME ParMooN_test_tcd3d_femfct_seq)
  set(_FILE_NAME tcd3d_femfct_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  add_test(NAME tcd3d_femfct_seq COMMAND ParMooN_test_tcd3d_femfct_seq)

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
  #multigrid with cell_vanka_jacobi smoother (works in SEQ and MPI)
  add_test(NAME nse3d_seq_cell_vanka_seq COMMAND ParMooN_test_nse3d_solver_seq cell_vanka_jacobi)
  # petsc
  add_test(NAME nse3d_seq_petsc COMMAND ParMooN_test_nse3d_solver_seq petsc)
  
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
  
  #tnse3d tests
  set(_EXEC_NAME ParMooN_test_tnse3d_seq)
  set(_FILE_NAME tnse3d_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  #umfpack
  add_test(NAME tnse3d_seq_umfpack COMMAND ParMooN_test_tnse3d_seq umfpack)
  add_test(NAME tnse3d_seq_umfpack_all COMMAND ParMooN_test_tnse3d_seq umfpack testall)
  #multigrid
  add_test(NAME tnse3d_seq_multigrid COMMAND ParMooN_test_tnse3d_seq multigrid)
  add_test(NAME tnse3d_seq_multigrid_all COMMAND ParMooN_test_tnse3d_seq multigrid testall)
  #lsc 
  add_test(NAME tnse3d_seq_lsc COMMAND ParMooN_test_tnse3d_seq lsc)
  add_test(NAME tnse3d_seq_lsc_all COMMAND ParMooN_test_tnse3d_seq lsc testall)
  
    
  #multigrid tests
  set(_EXEC_NAME ParMooN_test_multigrid_2d)
  set(_FILE_NAME multigrid_test.c++)
  register_2D_executable(${_EXEC_NAME} ${_FILE_NAME})
  set(_EXEC_NAME ParMooN_test_multigrid_3d)
  set(_FILE_NAME multigrid_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  #test 1
  add_test(NAME multigrid_2d COMMAND ParMooN_test_multigrid_2d)
  add_test(NAME multigrid_3d COMMAND ParMooN_test_multigrid_3d)

 #mesh3d tests
  set(_EXEC_NAME ParMooN_test_mesh3d_seq)
  set(_FILE_NAME mesh3d_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  #umfpack
  add_test(NAME mesh3d COMMAND ParMooN_test_mesh3d_seq umfpack)
    
  #sandwich test (sandwish grid)
  add_executable(ParMooN_test_sandwich sandwich_test.c++)
  target_link_libraries(ParMooN_test_sandwich parmoon_3d_${PARMOON_PARALLEL_TYPE} ${_EXTERN_LIBRARIES})
  add_test(sandwich ParMooN_test_sandwich)
  list(APPEND parmoon_tests ParMooN_test_sandwich)

  #lineseval tests 2d and 3d
  add_executable(ParMooN_test_linesEval2d linesEval_test.c++)
  target_link_libraries(ParMooN_test_linesEval2d
                        parmoon_2d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(lineseval2d ParMooN_test_linesEval2d)
  list(APPEND parmoon_tests ParMooN_test_linesEval2d)
  add_executable(ParMooN_test_linesEval3d linesEval_test.c++)
  target_link_libraries(ParMooN_test_linesEval3d
                        parmoon_3d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(lineseval3d ParMooN_test_linesEval3d)
  list(APPEND parmoon_tests ParMooN_test_linesEval3d)
  
  # collection tests 2d and 3d
  add_executable(ParMooN_test_collection_2d collection_test.c++)
  target_link_libraries(ParMooN_test_collection_2d
                        parmoon_2d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(collection2d ParMooN_test_collection_2d)
  list(APPEND parmoon_tests ParMooN_test_collection_2d)
  add_executable(ParMooN_test_collection_3d collection_test.c++)
  target_link_libraries(ParMooN_test_collection_3d
                        parmoon_3d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(collection3d ParMooN_test_collection_3d)
  list(APPEND parmoon_tests ParMooN_test_collection_3d)

  #snapshots, POD and ROM tests 2d
  add_executable(ParMooN_test_snaps_pod_rom2d snaps_pod_rom_test.c++)
  target_link_libraries(ParMooN_test_snaps_pod_rom2d
                        parmoon_2d_${PARMOON_PARALLEL_TYPE}
                        ${_EXTERN_LIBRARIES})
  add_test(snaps_pod_rom2d ParMooN_test_snaps_pod_rom2d)
  list(APPEND parmoon_tests ParMooN_test_snaps_pod_rom2d)

elseif(${PARMOON_PARALLEL_TYPE} STREQUAL MPI)	
  
  #####################################################
  ### A function with which to register executables ###   
  ### in the "check" target for mpi build           ###
  #####################################################
  function(register_executable exec_name file_name)
    
    add_executable(${exec_name} ${file_name})
    target_link_libraries(${exec_name} parmoon_3d_MPI ${_EXTERN_LIBRARIES})
    set(parmoon_tests ${parmoon_tests} ${exec_name} PARENT_SCOPE)
    
  endfunction(register_executable)
  
  
  ###################################################
  ####### Tests for parallel type: MPI #######
  ###################################################
  



  # 1st mumps solver test mpi (run with 2 processes)
  add_executable(ParMooN_test_mumps_1_mpi mumps_test_1.c++)
  target_link_libraries(ParMooN_test_mumps_1_mpi ${_EXTERN_LIBRARIES})
  add_test(NAME mumps_1_mpi COMMAND ${PETSC_MPIEXEC} -n 2 ParMooN_test_mumps_1_mpi )
  list(APPEND parmoon_tests ParMooN_test_mumps_1_mpi)
  
  # 2nd mumps solver test mpi (run with 4 processes)
  add_executable(ParMooN_test_mumps_2_mpi mumps_test_2.c++)
  target_link_libraries(ParMooN_test_mumps_2_mpi ${_EXTERN_LIBRARIES})
  add_test(NAME mumps_2_mpi COMMAND ${PETSC_MPIEXEC} -n 4 ParMooN_test_mumps_2_mpi )
  list(APPEND parmoon_tests ParMooN_test_mumps_2_mpi)
  
  # mumps wrapper test mpi (run with 2 processes)
  add_executable(ParMooN_test_mumps_wrapper mumps_wrapper_test.c++)
  target_link_libraries(ParMooN_test_mumps_wrapper parmoon_3d_MPI ${_EXTERN_LIBRARIES})
  add_test(NAME mumps_wrapper COMMAND ${PETSC_MPIEXEC} -n 2 ParMooN_test_mumps_wrapper )
  list(APPEND parmoon_tests ParMooN_test_mumps_wrapper)
  
  # domain partitioning read-in and write-out test
  register_executable(ParMooN_test_domain_partitioning_read_write domain_partitioning_read_write_test.c++)
  add_test(NAME partition_read_write COMMAND ${PETSC_MPIEXEC} -n 4 ParMooN_test_domain_partitioning_read_write)
  
  #cd3d solver tests
  set(_EXEC_NAME ParMooN_test_cd3d_solver_mpi)
  set(_FILE_NAME cd3d_solvers_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  # mumps
  add_test(NAME cd3d_mpi_mumps COMMAND ${PETSC_MPIEXEC} -n 3 ParMooN_test_cd3d_solver_mpi mumps)
  #jacobi
  add_test(NAME cd3d_mpi_jacobi COMMAND ${PETSC_MPIEXEC} -n 3 ParMooN_test_cd3d_solver_mpi jacobi)
  # using defualt petsc solver
  add_test(NAME cd3d_mpi_petsc COMMAND ${PETSC_MPIEXEC} -n 3 ParMooN_test_cd3d_solver_mpi 
           petsc)
  # mg with iterative solver fgmres and smoother SOR
  add_test(NAME cd3d_mpi_mg COMMAND ${PETSC_MPIEXEC} -n 3 ParMooN_test_cd3d_solver_mpi 
           multigrid)
  
  #tcd3d solver tests 
  set(_EXEC_NAME ParMooN_test_tcd3d_solver_mpi)
  set(_FILE_NAME tcd3d_solvers_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  # mumps
  add_test(NAME tcd3d_mpi_mumps COMMAND ${PETSC_MPIEXEC} -n 3 
           ParMooN_test_tcd3d_solver_mpi mumps)
  add_test(NAME tcd3d_mpi_jacobi COMMAND ${PETSC_MPIEXEC} -n 3 
           ParMooN_test_tcd3d_solver_mpi jacobi)
  # mg with iterative solver fgmres and smoother SOR
  add_test(NAME tcd3d_mpi_mg COMMAND ${PETSC_MPIEXEC} -n 3 ParMooN_test_tcd3d_solver_mpi 
           multigrid)
  
  #tcd3d algebraic flux correction tests with MPI
  set(_EXEC_NAME ParMooN_test_tcd3d_femfct_mpi)
  set(_FILE_NAME tcd3d_femfct_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  # only tested with Mumps so far
  add_test(NAME tcd3d_femfct_mpi COMMAND ${PETSC_MPIEXEC} -n 2 ParMooN_test_tcd3d_femfct_mpi)
  
  #nse3d solver tests
  set(_EXEC_NAME ParMooN_test_nse3d_solver_mpi)
  set(_FILE_NAME nse3d_solvers_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  #mumps
  add_test(NAME nse3d_mpi_mumps COMMAND ${PETSC_MPIEXEC} -n 3 ParMooN_test_nse3d_solver_mpi mumps)
  #multigrid with cell_vanka_jacobi smoother (works in SEQ and MPI), with iterative solver fgmres
  add_test(NAME nse3d_mpi_cell_vanka_jac COMMAND ${PETSC_MPIEXEC} -n 3 ParMooN_test_nse3d_solver_mpi cell_vanka_jacobi)
  #parallelized lsc
  add_test(NAME nse3d_mpi_lsc COMMAND ${PETSC_MPIEXEC} -n 3 ParMooN_test_nse3d_solver_mpi lsc)
  
  #tnse3d tests
  set(_EXEC_NAME ParMooN_test_tnse3d_mpi)
  set(_FILE_NAME tnse3d_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  #mumps
  add_test(NAME tnse3d_mpi_mumps COMMAND ${PETSC_MPIEXEC} -n 2 ParMooN_test_tnse3d_mpi mumps)
  add_test(NAME tnse3d_mpi_mumps_all COMMAND ${PETSC_MPIEXEC} -n 2 ParMooN_test_tnse3d_mpi mumps testall)
  #lsc
  add_test(NAME tnse3d_mpi_lsc COMMAND ${PETSC_MPIEXEC} -n 3 ParMooN_test_tnse3d_mpi lsc)
  add_test(NAME tnse3d_mpi_lsc_all COMMAND ${PETSC_MPIEXEC} -n 3 ParMooN_test_tnse3d_mpi lsc testall)
  
  
  #read_write_fe_function_mpi test
  set(_EXEC_NAME ParMooN_test_read_write_fe_function_mpi)
  set(_FILE_NAME read_write_fe_function_test.c++)
  register_executable(${_EXEC_NAME} ${_FILE_NAME})
  add_test(NAME read_write_fe_function_mpi COMMAND ${PETSC_MPIEXEC} -n 2 ParMooN_test_read_write_fe_function_mpi)
  
  #lineseval (2d and 3d) mpi tests (run with 4 processes)
  add_executable(ParMooN_test_linesEval3d linesEval_test.c++)
  target_link_libraries(ParMooN_test_linesEval3d
                        parmoon_3d_MPI ${_EXTERN_LIBRARIES})
  add_test(NAME lineseval3d_mpi COMMAND ${PETSC_MPIEXEC} -n 4
           ParMooN_test_linesEval3d)
  list(APPEND parmoon_tests ParMooN_test_linesEval3d)

endif(${PARMOON_PARALLEL_TYPE} STREQUAL SEQUENTIAL OR ${PARMOON_PARALLEL_TYPE} STREQUAL OMPONLY)

set(parmoon_tests ${parmoon_tests} PARENT_SCOPE)
