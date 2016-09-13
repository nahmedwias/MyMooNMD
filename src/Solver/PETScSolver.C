#include <PETScSolver.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>

#ifdef _MPI
#include <mpi.h>
#include <ParFECommunicator3D.h>
#endif // _MPI
#include <algorithm>

/**
 * Adding PETSc parameters from the database to the PETSc arguments string.
 *
 * Supported parameters are
 *  - maximum number of iterations
 *  - absolute residual tolerance
 *  - residual reduction
 *  - gmres restart
 *  - damping factor for Richardson iteration
 *
 * If a parameters is set in this function it will take precedence over
 * the parameter in the db["petsc_arguments"] string, e.g. if
 * your .dat file reads
 *
 * max_n_iterations: 1000
 *
 * [...]
 *
 * petsc_arguments: -pc_type lu -ksp_max_it 5
 *
 * The 1000 will be set in this function and your code will run with
 * 1000 iterations maximum instead of 5.
 *
 * @param db[in] database from the solver
 * @return PETSc arguments conforming string
 */
std::string addParameters2str(const ParameterDatabase& db)
{
  std::string petsc_args = "";
  if(db.contains("petsc_arguments"))
  {
    petsc_args = db["petsc_arguments"].value_as_string();
  }
  // adding maximum number of iterations
  petsc_args.append(" -ksp_max_it "+db["max_n_iterations"].value_as_string());
  // adding absolute tolerance
  petsc_args.append(" -ksp_atol "+db["residual_tolerance"].value_as_string());
  // adding reduction tolerance
  petsc_args.append(" -ksp_rtol "+db["residual_reduction"].value_as_string());
  // adding gmres restart
  petsc_args.append(" -ksp_gmres_restart "+db["gmres_restart"].value_as_string());
  // adding damping for Richardson iteration
  petsc_args.append(" -ksp_richardson_scale "+db["damping_factor"].value_as_string());
  return petsc_args;
}

/**
 * @brief Converts a std::string to a vector of char* by splitting words at 
 * white spaces.
 * 
 * The idea is that the returned vector mimics the standard arguments of a C/C++
 * main function: int main(int argc, char** argv). These arguments are needed to
 * properly call PetscInitialize.
 * 
 * @param[in] s input string
 * @return vector of char* containing all words + nullptr
 */
std::vector<char*> str_to_vector_char_p(const std::string &s)
{
  // number of arguments, the first is the program name by default
  size_t n_arg = 1;
  // a position in the string s
  auto current_position = s.find(" ");
  // count the white spaces and therefore the number of arguments
  while(current_position != std::string::npos)
  {
    n_arg++;
    current_position = s.find(" ", current_position + 1);
  }
  // the number of white spaces is smaller than the number arguments by 1
  n_arg++;
  
  // store where words (separated by words) are
  std::vector<size_t> word_positions;
  word_positions.reserve(n_arg+1);
  // zero is always the first position
  word_positions.push_back((size_t)0);

  // fill the vector word_positions
  current_position = s.find(" ");
  while(current_position != std::string::npos)
  {
    word_positions.push_back(current_position + 1);
    current_position = s.find(" ", current_position + 1);
  }

  // push the length of the string as last "position"
  // now word_positions[i+1]-word_positions[i] is valid for i = 0, ... n_arg-1
  word_positions.push_back(s.size()+1);

  // plus 1 here for the terminating nullptr, this will be returned
  std::vector<char*> argv(n_arg+1);

  // PETSc expects the program name at first.
  char* progname = new char[8];
  strcpy(progname, "./dummy");
  argv[0] = &progname[0];

  for(size_t i = 1; i < n_arg; ++i)
  {
    size_t wordLength = word_positions[i] - word_positions[i-1] - 1;
    argv[i] = new char[wordLength+1];
    strncpy( argv[i], &s.c_str()[word_positions[i-1]], wordLength );
    // add terminating \0
    argv[i][wordLength] = '\0';
  }

  // set the nullptr at the end
  // demanded by the C++ Standard for main function argument argv
  argv[n_arg] = nullptr;
  return argv;
}

PETScSolver::PETScSolver(const BlockFEMatrix& matrix,
                         const ParameterDatabase& db)
 : petsc_mat(nullptr)
#ifdef _MPI
  , comms_(matrix.get_communicators())
#endif // _MPI
{
  // we will create petsc_mat (type: Mat). All entries of 'matrix' will be
  // copied to it. That copying is really slow if petsc_mat has not been created
  // with a known sparsity structure. This would cause many memory allocations 
  // and copies. Instead all that PETSc needs to know at first is the number of 
  // non-zero entries in each row. This is therefore computed first here.
#if (defined _MPI && defined __2D__)
  static_assert(false, "you can not use MPI in 2D (yet)");
#endif
  // sanity check:
  if(matrix.get_n_total_rows() != matrix.get_n_total_columns())
    ErrThrow("PETScSolver: non-square matrix ", matrix.get_n_total_rows(),
             ",", matrix.get_n_total_columns());
  size_t n_block_rows = matrix.get_n_cell_rows();
  size_t n_block_cols = matrix.get_n_cell_columns();
#ifdef _MPI
  if(n_block_rows * n_block_cols != 1)
    ErrThrow("in MPI the PETScSolver only works for scalar problems currently");
#endif

  // indicate if matrix is stored as a transposed, this is set when calling
  // matrix.get_block(...)
  bool transposed;
  // indicate if we have a saddle point problem
  /// @todo save this information in the BlockMatrix
  bool saddlepoint_problem = false;
  // this is considered a saddle point problem if the last block is all zero
  auto block = matrix.get_block(n_block_rows-1, n_block_cols-1, transposed);
  if(block->GetNorm() == 0.0)
  {
    saddlepoint_problem = true;
  }

  // call PetscInitialize
  {
    // the following is necessary to properly call PetscInitialize
    std::string petsc_args = addParameters2str(db);
    if(saddlepoint_problem)
    {
      // PETSc can figure out which block is zero. Also this is the only way to
      // use specialized saddle point solvers
      petsc_args.append(" -pc_fieldsplit_detect_saddle_point");
    }
    auto char_vector = str_to_vector_char_p(petsc_args);
    char ** params_pointer = char_vector.data();
    int argc = char_vector.size();

    PetscInitialize(&argc, &params_pointer, (char*)0, nullptr);

    // properly delete all the char* which where created
    // in str_to_vector_char_p(...)
    for(auto cp : char_vector)
      delete [] cp;

    // Is pc_type LU, ILU or CHOLESKY? I.e. do we want to use a direct solver
    bool use_direct_petsc = false;
    auto npos = std::string::npos;
    if ( petsc_args.find("-pc_type lu", 0) != npos
      || petsc_args.find("-pc_type ilu", 0) != npos
      || petsc_args.find("-pc_type cholesky", 0) != npos)
    {
      use_direct_petsc = true;
    }

    // for block matrices the direct solver does not seem to work in PETSc
    if((saddlepoint_problem || n_block_rows*n_block_cols != 1)
        && use_direct_petsc)
    {
      ErrThrow("solver_type is petsc, with a direct solver for a saddle point problem. "
               "Please set the solver_type to direct when attempting to solve "
               "a saddle point problem with a direct solver.");
    }
  }

  // reserve space for enough sub matrices for PETSc
  sub_petsc_mats.resize(n_block_rows * n_block_cols);

  // copy the entries from the blocks in 'matrix' to petsc_mat:
  // loop over all block rows
  for(size_t block_row = 0; block_row < n_block_rows;
      ++block_row)
  {
    // loop over all block columns
    for(size_t block_col = 0; block_col < n_block_cols; ++block_col)
    {
      // the current block (sparse matrix)
      auto block = matrix.get_block(block_row, block_col, transposed);
      if (transposed)
      {
        /// @todo enable transposed blocks to be used in the PETScSolver class
        ErrThrow("It is not yet possible to use the PETScSolver with a "
                 "BlockFEMatrix which stores a block as transposed.");
      }
      
      size_t block_index = block_col + block_row * n_block_cols;
      bool isOnDiag = (block_row == block_col);
      FEBlock2PETScBlock(block, isOnDiag, sub_petsc_mats[block_index]);
    }
  }

  // if scalar problem
  if (n_block_cols == 1 && n_block_rows == 1)
  {
    // the vector sub_petsc_mats is not needed in this case at all
    petsc_mat = sub_petsc_mats[0];
  }
  else
  {
    // This is the PETSc way of having a BlockMatrix, direct solvers wont work
    MatCreateNest(PETSC_COMM_WORLD, n_block_rows, NULL, n_block_cols, NULL,
                  &sub_petsc_mats[0], &petsc_mat);
  }
  
  // print out the petsc matrix to console
  //PetscViewer viewer;
  //PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  //PetscViewerSetType(viewer, PETSCVIEWERASCII);
  //PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_COMMON);
  //MatView(petsc_mat, viewer);
  //PetscViewerDestroy(&viewer);
  //PetscReal norm;
  //MatNorm(petsc_mat, NORM_FROBENIUS, &norm);
  //Output::print("norm mat ", norm, "  ", 
  //              matrix.get_block(0, 0, transposed)->GetNorm(-2));
  
  // create solver and preconditioner objects
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, petsc_mat, petsc_mat);
  
  // PETSc preconditioner object
  PC  pc;
  // set preconditioner context (PC)
  KSPGetPC(ksp, &pc);
  
  // set the options from the string passed to PetscInitialize which includes
  // the values in the parameter "petsc_arguments".
  KSPSetFromOptions(ksp);
  PCSetFromOptions(pc);
  
  // some PETSc information about the solver
  /// @todo nicer output for the PETSc solver
  KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
}

/* ************************************************************************** */
PETScSolver::PETScSolver(const BlockMatrix& matrix, const ParameterDatabase& db)
{
  ErrThrow("can not construct a PETScSolver using a BlockMatrix. I need a "
           "BlockFEMatrix instead. It would be possible to implement this "
           "here, in case it is really needed.");
}


/* ************************************************************************** */
PETScSolver::PETScSolver(PETScSolver && other)
{
  // the type of this->petsc_mat is 'Mat' which is in fact a pointer
  // destroy the old petsc_mat
  MatDestroy(&petsc_mat);
  // copy the pointer
  this->petsc_mat = other.petsc_mat;
  // reset to some zero matrix in the other object
  int n_rows = 0;
  int n_cols = 0;
  int n_nonzero_per_row = 0;
  MatCreateSeqAIJ(PETSC_COMM_WORLD, n_rows, n_cols, n_nonzero_per_row, nullptr,
                  &other.petsc_mat);
}

/* ************************************************************************** */
void PETScSolver::FEBlock2PETScBlock(std::shared_ptr<const FEMatrix> feblock,
                                     bool isOnDiag, Mat &sub_mat)
{
  size_t n_rows = feblock->GetN_Rows();
  size_t n_active_rows = feblock->GetActiveBound();

  // unfortunately the sequential and parallel code differs a lot here.
  
#ifndef _MPI 
  // sequential case
  // 1) create matrix
  // number of non-zeros for each row
  std::vector<int> nnz(n_rows, 0);
  // fill vector nnz
  for(size_t row = 0; row < n_rows; ++row)
  {
    nnz[row] += feblock->get_n_entries_in_row(row);
  }
  size_t n_cols = feblock->GetN_Columns();
  // create a matrix with known non-zero distribution among the rows
  MatCreateSeqAIJ(PETSC_COMM_WORLD, n_rows, n_cols, 0, &nnz[0], &sub_mat);
  
  // 2) fill matrix
  // get raw data of the matrix
  const int * row_ptr = feblock->GetRowPtr();
  const int * col_ptr = feblock->GetKCol();
  const double * entries = feblock->GetEntries();
  for(size_t row = 0; row < n_rows; ++row)
  {
    size_t n_entries_in_row = nnz[row];
    int row_index = row;
    // only add this row if it is active or in the diagonal block
    if(row < n_active_rows || isOnDiag)
      MatSetValues(sub_mat, 1, &row_index, n_entries_in_row,
          &col_ptr[row_ptr[row]], entries+row_ptr[row],
          INSERT_VALUES);
  }
#else
  // MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
  // 1) create matrix
  
  // communicator for the row space (test space) of this matrix
  const TParFECommunicator3D& comm_row = 
    feblock->GetTestSpace3D()->get_communicator();
  // for PETSc only the master rows belong to a (process local) matrix
  int n_master_rows = comm_row.GetN_Master();
  size_t n_global_dof = comm_row.get_n_global_dof();
  std::vector<int> nnz_diagonal(n_master_rows, 0);
  std::vector<int> nnz_off_diagonal(n_master_rows, 0);
  const int * masters = comm_row.GetMaster();
  // fill vectors nnz_diagonal and nnz_off_diagonal
  // loop over all rows
  for(size_t row = 0, master_index = 0; row < n_rows; ++row)
  {
    // skip slave rows
    if(masters[row] != my_rank)
      continue;
    
    // number of nonzero entries in a master column in this row
    size_t n_masters_in_row = 0;
    int begin_row = feblock->GetRowPtr()[row];
    int end_row = feblock->GetRowPtr()[row+1];
    // loop over all entries in this row
    for(int index = begin_row; index < end_row; index++)
    {
      int column = feblock->GetKCol()[index]; // column index
      if(masters[column] == my_rank)
        n_masters_in_row++;
    }
    nnz_diagonal.at(master_index) = n_masters_in_row;
    nnz_off_diagonal.at(master_index) = feblock->get_n_entries_in_row(row)
                                        - n_masters_in_row;
    master_index++;
  }
  //Output::print("n_master_rows: ", n_master_rows, "   n_rows ", n_rows, 
  //              "  n_cols ", feblock->GetN_Columns(), "  n_global_dof ", 
  //              n_global_dof);
  MatCreateAIJ(PETSC_COMM_WORLD, n_master_rows, n_master_rows, 
               (int)n_global_dof, (int)n_global_dof, 0, &nnz_diagonal[0], 0, 
               &nnz_off_diagonal[0], &sub_mat);
  
  // 2) fill matrix
  // get raw data of the matrix
  const std::vector<int>& row_ptr = feblock->get_row_array();
  const int * col_ptr = feblock->GetKCol();
  const double * entries = feblock->GetEntries();
  const int * local2global = comm_row.Get_Local2Global();
  // loop over all rows, fill each row at once
  for(size_t row_local = 0; row_local < n_rows; ++row_local)
  {
    // skip slave rows
    if(masters[row_local] != my_rank)
      continue;
    
    int n_entries_in_row = feblock->get_n_entries_in_row(row_local);
    // number of nonzero entries in a master column in this row
    int begin_row = row_ptr[row_local];
    int end_row = row_ptr.at(row_local+1);
    std::vector<int> columns_global(n_entries_in_row, 0);
    // we need to pass the global column indices to PETSc, we can not use the 
    // local column indices from the feblock.
    // loop over all entries in this row
    for(int index = begin_row; index < end_row; index++)
    {
      // local (to this process) column index
      int column_local = col_ptr[index]; // column index
      // write global column index
      columns_global[index-begin_row] = local2global[column_local];
      //Output::print("entry (", local2global[row_local], ",", 
      //              columns_global[index-begin_row], ")_global = (", 
      //              row_local, ",", column_local, ") = ", entries[index]);
    }
    int row_global = local2global[row_local];
    if(row_local < n_active_rows || isOnDiag)
      MatSetValues(sub_mat, 1, &row_global, n_entries_in_row,
          &columns_global[0], entries+begin_row, INSERT_VALUES);
  }
#endif // _MPI
  MatAssemblyBegin(sub_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(sub_mat, MAT_FINAL_ASSEMBLY);
}

/* ************************************************************************** */
PETScSolver& PETScSolver::operator=(PETScSolver && other)
{
  // the type of this->petsc_mat is 'Mat' which is in fact a pointer
  // destroy the old petsc_mat
  MatDestroy(&petsc_mat);
  // copy the pointer
  this->petsc_mat = other.petsc_mat;
  // reset to some zero matrix in the other object
  int n_rows = 0;
  int n_cols = 0;
  int n_nonzero_per_row = 0;
  MatCreateSeqAIJ(PETSC_COMM_WORLD, n_rows, n_cols, n_nonzero_per_row, nullptr,
                  &other.petsc_mat);
  return *this;
}

/* ************************************************************************** */
PETScSolver::~PETScSolver()
{
  KSPDestroy(&ksp);
  // if we have just one sub matrix
  // it is the petsc_mat
  if (sub_petsc_mats.size() > 1)
  {
    for (size_t i = 0; i < sub_petsc_mats.size(); ++i)
    {
      MatDestroy(&sub_petsc_mats[i]);
    }
  }
  MatDestroy(&petsc_mat);
  // remember you have to call 'PetscFinalize();' at the very end of the program
}

/* ************************************************************************** */
void PETScSolver::solve(const BlockVector& rhs, BlockVector& solution)
{
  /// @todo check if matrix sizes are ok with solution and rhs sizes
  /// MatGetSize() does not work for MatNest mat type, which is used for
  /// saddle point problems.

  size_t n_local = solution.length();

  if(n_local != rhs.length())
  {
    ErrThrow("PETScSolver::solve: size of the rhs and the solution are not equal, ",
                 rhs.length(), " != ", n_local);
  }
  
  // length of vector
  int n_global = n_local;
  int n_local_masters = n_local;
#ifdef _MPI
  n_local_masters = comms_.at(0)->GetN_Master();
  n_global = comms_[0]->get_n_global_dof();
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  const int * masters = comms_[0]->GetMaster();
#endif
  
  // create PETSc vectors
  Vec x, b;
  VecCreate(PETSC_COMM_WORLD, &x);
  PetscObjectSetName((PetscObject) x, "Solution");
  VecSetSizes(x, n_local_masters, n_global);
  VecSetFromOptions(x);
  VecDuplicate(x, &b);
  PetscObjectSetName((PetscObject) b, "Rhs");
  
  // copy into PETSc vectors
  for(size_t i = 0; i < n_local; ++i)
  {
    int index = i; // index in vector (conversion to int)
#ifdef _MPI
    // skip slave dofs
    if(masters[i] != my_rank)
      continue;
    index = comms_[0]->Get_Local2Global()[i];
#endif // _MPI
    VecSetValue(b, index, rhs[i], INSERT_VALUES);
    VecSetValue(x, index, solution[i], INSERT_VALUES);
  }
  VecAssemblyBegin(x);
  VecAssemblyBegin(b);
  VecAssemblyEnd(x);
  VecAssemblyEnd(b);
  
  // solve the linear system using PETSc
  KSPSolve(ksp, b, x);
  
  PetscInt its;
  KSPGetIterationNumber(ksp, &its);
#ifdef _MPI
  if(my_rank == 0)
#endif
  Output::print("some PETSc solver: number of iterations: ", its);
  
  // copy back to solution:
  for(size_t i = 0; i < n_local; ++i)
  {
    int index = i; // index in vector (conversion to int)
#ifdef _MPI
    // skip slave dofs
    if(masters[i] != my_rank)
      continue;
    index = comms_[0]->Get_Local2Global()[i];
#endif // _MPI
    VecGetValues(x, 1, &index, &solution[i]);
  }
  
  VecDestroy(&b);
  VecDestroy(&x);
}

