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
 * Converts a std::string to a char** by splitting words at white spaces.
 *
 * @param[in] s input string
 * @param[out] argc number of words in the input string
 * @param[out] argv char** containing all words + nullptr
 */
void str2charpp(const std::string &s, int &argc, char*** argv)
{
  // zero is always the first position
  std::vector<size_t> wordPositions;
  wordPositions.push_back((size_t)0);

  auto current_position = s.find(" ");

  // argc = 1 for the "progname"
  argc = 1;

  // count the white spaces and save the positions next to them
  while(current_position != std::string::npos)
  {
	argc++;
	wordPositions.push_back(current_position + 1);
	current_position = s.find(" ", current_position + 1);
  }

  // push the length of the string as last "position"
  // now position[i+1]-position[i] is valid for i = 0, ... argc
  argc++;
  wordPositions.push_back(s.size()+1);

  // one extra entry for the file name. PETSc expects the program name at first
  char* progname = new char[8];
  progname[0] = '.'; progname[1] = '/'; progname[2] = 'd';
  progname[3] = 'u'; progname[4] = 'm'; progname[5] = 'm';
  progname[6] = 'y'; progname[7] = '\0';

  // plus 1 here for the terminating NULL
  *argv = (char**)malloc(sizeof(char*)*(argc+1));

  // filename first
  (*argv)[0] = &progname[0];

  for (auto i = 1; i < argc; ++i)
  {
	size_t wordLength = wordPositions[i]-wordPositions[i-1] - 1;

	(*argv)[i] = (char*)malloc(sizeof(char)*(wordLength+1));
	strncpy( (*argv)[i], &s.c_str()[wordPositions[i-1]], wordLength );

	// add terminating \0
	(*argv)[i][wordLength] = '\0';
  }

  // set the nullptr at the end
  // demanded by the C++ Standard for main function argument argv
  (*argv)[argc] = nullptr;

}


PETScSolver::PETScSolver(const BlockFEMatrix& matrix, const ParameterDatabase& db)
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

  std::string petsc_args = addParameters2str(db);

  Output::print("Out new argument string\n", petsc_args, "\n");

  char ** params_pointer;
  int argc;

  str2charpp(petsc_args, argc, &params_pointer);

  PetscInitialize(&argc, &params_pointer, (char*)0, nullptr);
  
  size_t n_rows = matrix.get_n_total_rows();
  size_t n_cols = matrix.get_n_total_columns();
  size_t n_entries = matrix.get_n_total_entries();
  if(n_rows != n_cols)
    ErrThrow("PETScSolver: non-square matrix");
  size_t n_block_rows = matrix.get_n_cell_rows();
  size_t n_block_cols = matrix.get_n_cell_columns();
  
  // number of non-zeros for each row
  std::vector<int> nnz(n_rows, 0);
  // indicate if matrix is stored as a transposed
  bool transposed;
  // fill vector nnz
  for(size_t block_row = 0, row_offset = 0; block_row < n_block_rows; 
      ++block_row)
  {
    for(size_t block_col = 0; block_col < n_block_cols; ++block_col)
    {
      auto block = matrix.get_block(block_row, block_col, transposed);
      if(transposed)
      {
        ErrThrow("It is not yet possible to use the PETScSolver with a "
                 "BlockFEMatrix which stores a block as transposed.");
      }
      size_t n_rows = block->GetN_Rows();
      for(size_t row = 0; row < n_rows; ++row)
      {
        nnz[row+row_offset] += block->get_n_entries_in_row(row);
      }
      if(block_col == n_block_cols-1)
      {
        row_offset += n_rows;
      }
    }
  }
  
  // create a matrix with known non-zero distribution among the rows
  MatCreateSeqAIJ(PETSC_COMM_WORLD, n_rows, n_cols, n_entries, &nnz[0],
                  &petsc_mat);
  // MatSetUp must be called before MatSetValues (says PETSc documentation)
  MatSetUp(petsc_mat);
  
  // now we need to copy the entries from the blocks in 'matrix' to petsc_mat:
  // loop over all blocks rows
  for(size_t block_row = 0, row_offset = 0; block_row < n_block_rows;
      ++block_row)
  {
    size_t column_offset = 0;
    // loop over all block columns
    for(size_t block_col = 0; block_col < n_block_cols; ++block_col)
    {
      // the current block (sparse matrix)
      auto block = matrix.get_block(block_row, block_col, transposed);
      // transposed is false, we already checked this earlier
      
      // insert row-wise
      size_t n_rows = block->GetN_Rows();
      size_t n_active_rows = block->GetActiveBound();
      // get raw data of the matrix
      const int * row_ptr = block->GetRowPtr();
      const int * col_ptr = block->GetKCol();
      const double * entries = block->GetEntries();
      // copy columns into a new vector, needed to change the global column 
      // index. If this matrix only has one block this is not necessary.
      std::vector<int> colums(col_ptr, col_ptr + block->GetN_Entries());
      std::for_each(colums.begin(), colums.end(),
                    [column_offset](int & c){c += column_offset;});
      
      for(size_t row = 0; row < n_rows; ++row)
      {
        size_t n_entries_in_row = row_ptr[row+1] - row_ptr[row];
        int row_index = row + row_offset;
        // only add this row if it is active or in the diagonal block
        if(row < n_active_rows || block_row == block_col)
          MatSetValues(petsc_mat, 1, &row_index, n_entries_in_row, 
                       &colums[row_ptr[row]], entries+row_ptr[row], 
                       INSERT_VALUES);
      }
      
      column_offset += block->GetN_Columns();
      if(block_col == n_block_cols-1)
      {
        row_offset += n_rows;
      }
    }
  }
  
  // the following two functions must be called after MatSetValues (says PETSc 
  // documentation)
  MatAssemblyBegin(petsc_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(petsc_mat, MAT_FINAL_ASSEMBLY);
  
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
  MatDestroy(&petsc_mat);
  // remember you have to call 'PetscFinalize();' at the very end of the program
}

/* ************************************************************************** */
void PETScSolver::solve(const BlockVector& rhs, BlockVector& solution)
{
  PetscInt n, m;
  MatGetSize(petsc_mat, &m, &n);
  if((int)rhs.length() != m)
  {
    ErrThrow("PETScSolver::solver: size of the right hand side is not "
             "suitable ", rhs.length(), " != ", m);
  }
  if((int)solution.length() != n)
  {
    ErrThrow("PETScSolver::solver: size of the solution is not suitable ",
             solution.length(), " != ", n);
  }
  
  // create PETSc vectors
  Vec x, b;
  VecCreate(PETSC_COMM_WORLD, &x);
  PetscObjectSetName((PetscObject) x, "Solution");
  VecSetSizes(x, PETSC_DECIDE, solution.length());
  VecSetFromOptions(x);
  VecDuplicate(x, &b);
  
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
  KSPGetPC(ksp, &pc);
  
  // http://scicomp.stackexchange.com/questions/513/why-is-my-iterative-linear-
  // solver-not-converging
  // http://scicomp.stackexchange.com/questions/7288/which-preconditioners-and-
  // solver-in-petsc-for-indefinite-symmetric-systems-sho

  PCSetFromOptions(pc);

  // obviously this needs to be controlled from outside via some database 
  // entries
  double relative_tolerance = 1.e-14;
  double absolute_tolerance = 1.e-20;
  int max_iterations = 1000;
  KSPSetTolerances(ksp, relative_tolerance, absolute_tolerance, PETSC_DEFAULT,
                   max_iterations);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp,b,x);

  // some PETSc information about the solver
  KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
  PetscInt its;
  KSPGetIterationNumber(ksp,&its);
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

