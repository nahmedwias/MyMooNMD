#include <PETScSolver.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>

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
{
#ifdef MPI
  ErrThrow("PETScSolver not yet implemented for MPI ");
#endif
  
  // we will create petsc_mat (type: Mat). All entries of 'matrix' will be
  // copied to it. That copying is really slow if petsc_mat has not been created
  // with a known sparsity structure. This would cause many memory allocations 
  // and copies. Instead all that PETSc needs to know at first is the number of 
  // non-zero entries in each row. This is computed first here.

  size_t n_rows = matrix.get_n_total_rows();
  size_t n_cols = matrix.get_n_total_columns();
  if(n_rows != n_cols)
    ErrThrow("PETScSolver: non-square matrix");
  size_t n_block_rows = matrix.get_n_cell_rows();
  size_t n_block_cols = matrix.get_n_cell_columns();

  // indicate if matrix is stored as a transposed
  bool transposed;
  
  // indicate if we have a saddle point problem
  /// @todo save this information in the BlockMatrix
  bool saddlepoint_problem = false;

  // check for saddle point problems
  auto block = matrix.get_block(n_block_rows-1, n_block_cols-1, transposed);
  if(block->GetNorm() == 0.0)
  {
    // this seems to be a saddle point problem.
    saddlepoint_problem = true;
  }

  // call PetscInitialize
  {
    // the following is necessary to properly call PetscInitialize
    std::string petsc_args = addParameters2str(db);

    if(saddlepoint_problem)
    {
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

    /*
     * Is pc_type LU, ILU or CHOLESKY?
     * I.e. do we want to use a direct solver
     * for a saddle point problem?
     */
    bool use_direct_petsc = false;

    auto npos = std::string::npos;

    if ( petsc_args.find("-pc_type lu", 0) != npos
      || petsc_args.find("-pc_type ilu", 0) != npos
      || petsc_args.find("-pc_type cholesky", 0) != npos)
    {
      use_direct_petsc = true;
    }

    if(saddlepoint_problem && use_direct_petsc)
    {
      ErrThrow("solver_type is petsc, with a direct solver for a saddle point problem. "
               "Please set the solver_type to direct when attempting to solve "
               "a saddle point problem with a direct solver.");
    }
  }

  // now we know how much blocks we have in the BlockFEMatrix
  // and can reserve space for enough sub matrices for PETSc
  sub_petsc_mats.resize(n_block_rows * n_block_cols);

  // now we need to copy the entries from the blocks in 'matrix' to petsc_mat:
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
        ErrThrow("It is not yet possible to use the PETScSolver with a "
                 "BlockFEMatrix which stores a block as transposed.");
      }
      
      size_t block_index = block_col + block_row * n_block_cols;

      bool isOnDiag = block_row==block_col;

      FEBlock2PETScBlock(block, isOnDiag,
                         sub_petsc_mats[block_index]);
    }
  }

  // if scalar problem
  if (n_block_cols == 1 && n_block_rows == 1)
  {
    petsc_mat = sub_petsc_mats[0];
  }
  else
  {
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
  size_t n_cols = feblock->GetN_Columns();
  size_t n_entries = feblock->GetN_Entries();

  // number of non-zeros for each row
  std::vector<int> nnz(n_rows, 0);
  // fill vector nnz
  for(size_t row = 0; row < n_rows; ++row)
  {
    nnz[row] += feblock->get_n_entries_in_row(row);
  }

  // create a matrix with known non-zero distribution among the rows
  MatCreateSeqAIJ(PETSC_COMM_WORLD, n_rows, n_cols, n_entries, &nnz[0],
                  &sub_mat);

  // get raw data of the matrix
  const int * row_ptr = feblock->GetRowPtr();
  const int * col_ptr = feblock->GetKCol();
  const double * entries = feblock->GetEntries();

  for(size_t row = 0; row < n_rows; ++row)
  {
  	size_t n_entries_in_row = row_ptr[row+1] - row_ptr[row];
  	int row_index = row;
  	// only add this row if it is active or in the diagonal block
  	if(row < n_active_rows || isOnDiag)
  		MatSetValues(sub_mat, 1, &row_index, n_entries_in_row,
  		    &col_ptr[row_ptr[row]], entries+row_ptr[row],
  				INSERT_VALUES);
  }

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

  PetscInt m = rhs.length();
  PetscInt n = solution.length();

  if(n != m)
  {
    ErrThrow("PETScSolver::solve: size of the rhs and the solution are not equal, ",
                 m, " != ", n);
  }
  
  // create PETSc vectors
  Vec x, b;
  VecCreate(PETSC_COMM_WORLD, &x);
  PetscObjectSetName((PetscObject) x, "Solution");
  VecSetSizes(x, PETSC_DECIDE, solution.length());
  VecSetFromOptions(x);
  VecDuplicate(x, &b);
  PetscObjectSetName((PetscObject) b, "Rhs");
  
  // copy into PETSc vectors
  for(int i = 0; i < n; ++i)
  {
    VecSetValue(b, i, rhs[i], INSERT_VALUES);
    VecSetValue(x, i, solution[i], INSERT_VALUES);
  }
  
  // create solver and preconditioner objects
  KSP ksp;         /* linear solver context */
  PC  pc;          /* preconditioner context */
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, petsc_mat, petsc_mat);

  // set preconditioner context (PC)
  KSPGetPC(ksp, &pc);
  
  // set the options from the string passed to PetscInitialize which includes
  // the values in the parameter "petsc_arguments".
  KSPSetFromOptions(ksp);
  PCSetFromOptions(pc);

  KSPSolve(ksp, b, x);

  // some PETSc information about the solver
  /// @todo nicer output for the PETSc solver
  // KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
  PetscInt its;
  KSPGetIterationNumber(ksp, &its);
  Output::print("some PETSc solver: number of iterations: ", its);
  
  // copy back to solution:
  for(int i = 0; i < n; ++i)
  {
    VecGetValues(x, 1, (int*)&i, &solution[i]);
  }
  
  VecDestroy(&b);
  VecDestroy(&x);
  KSPDestroy(&ksp);
}

