// =======================================================================
// @(#)DirectSolver.h
//
// Purpose:     solve equation system by direct solver
//
// Author:      Gunar Matthies (06.09.05)
//
// History:     start of implementation 06.09.05 (Gunar Matthies)
//
// =======================================================================

#include <MainUtilities.h>
#include <DirectSolver.h>
#include <Database.h>
#include "umfpack.h"
#include <BlockMatrix.h>

void handle_error_umfpack(int ierror)
{
  if (ierror == UMFPACK_OK)
    return;

  switch(ierror)
  {
    //WARNINGS
    case UMFPACK_WARNING_singular_matrix:
      Output::print<2>("umfpack Warning: Matrix is singular!");
      break;
    case UMFPACK_WARNING_determinant_underflow:
      Output::print<2>("umfpack Warning: Determinant smaller than eps");
      break;
    case UMFPACK_WARNING_determinant_overflow:
      Output::print<2>("umfpack Warning: Determinant is larger than IEEE Inf");
      break;
    //ERRORS
    case UMFPACK_ERROR_out_of_memory:
      ErrThrow("umfpack: Out of Memory");
      break;
    case UMFPACK_ERROR_invalid_Numeric_object:
      ErrThrow("umfpack: Invalid numeric factorization object");
      break;
    case UMFPACK_ERROR_invalid_Symbolic_object:
      ErrThrow("umfpack: Invalid symbolic factorization object");
      break;
    case UMFPACK_ERROR_argument_missing:
      ErrThrow("umfpack: Argument Missing.");
      break;
    case UMFPACK_ERROR_n_nonpositive:
      ErrThrow("umfpack: Matrix dimensions not positive.");
      break;
    case UMFPACK_ERROR_invalid_matrix:
      ErrThrow("umfpack: Invalid Matrix Structure.");
      break;
    case UMFPACK_ERROR_different_pattern:
      ErrThrow("umfpack: Different sparse pattern.");
      break;
    case UMFPACK_ERROR_invalid_system:
      ErrThrow("umfpack: Invalid system provided with sys.");
      break;
    case UMFPACK_ERROR_invalid_permutation:
      ErrThrow("umfpack: Invalid permutation vector.");
      break;
    case UMFPACK_ERROR_file_IO:
      ErrThrow("umfpack: Fille IO error.");
      break;
    case UMFPACK_ERROR_internal_error:
      ErrThrow("umfpack: Internal error.");
      break;
    
    default:
      ErrThrow("umfpack: unkown error. Error number ", ierror); break;
    break;
  }
}

/** ************************************************************************ */
DirectSolver::DirectSolver(std::shared_ptr<TMatrix> matrix, 
                           DirectSolver::DirectSolverTypes type)
 : type(type), matrix(matrix), symbolic(nullptr), numeric(nullptr), cols(), 
   rows()
{
  Output::print<3>("constructing a DirectSolver object");
  if(!matrix->is_square())
  {
    ErrThrow("unable to factorize a non-square matrix ", matrix->GetN_Rows(),
             "  ", matrix->GetN_Columns());
  }
  
  if(type == DirectSolverTypes::pardiso)
  {
    ErrThrow("Pardiso does not yet work");
  }
  
  // the threshold is rather small here, it should furthermore depend on the 
  // dimension (2 or 3) and the polynomial degree (and possibly more).
  if(this->matrix->GetN_Rows() > 2e5)
  {
    this->cols.resize(this->matrix->GetN_Entries(), 0);
    this->rows.resize(this->matrix->GetN_Rows()+1, 0);
    for(size_t i = 0; i < this->matrix->GetN_Entries(); ++i)
      this->cols[i] = this->matrix->GetKCol()[i];
    for(size_t i = 0; i < this->matrix->GetN_Rows()+1; ++i)
      this->rows[i] = this->matrix->GetRowPtr()[i];
  }
  this->symetric_factorize();
  this->numeric_factorize();
}

/** ************************************************************************ */
DirectSolver::DirectSolver(const BlockMatrix& matrix, 
                           DirectSolver::DirectSolverTypes type)
 : DirectSolver(matrix.get_combined_matrix(), type)
{
}

/** ************************************************************************ */
DirectSolver::DirectSolver(DirectSolver&& other)
 : type(other.type), matrix(other.matrix), symbolic(other.symbolic), 
   numeric(other.numeric), cols(std::move(other.cols)), 
   rows(std::move(other.rows))
{
  other.symbolic = nullptr;
  other.numeric = nullptr;
  Output::print<4>("DirectSolver::DirectSolver(DirectSolver&&)");
}

/** ************************************************************************ */
class DirectSolver& DirectSolver::operator=(DirectSolver&& other)
{
  this->type = other.type;
  this->matrix = other.matrix;
  this->symbolic = other.symbolic;
  this->numeric = other.numeric;
  this->cols = std::move(other.cols);
  this->rows = std::move(other.rows);
  other.symbolic = nullptr;
  other.numeric = nullptr;
  Output::print<4>("DirectSolver::operator=(DirectSolver&&)");
  return *this;
}

/** ************************************************************************ */
DirectSolver::~DirectSolver()
{
  switch(type)
  {
    case DirectSolver::DirectSolverTypes::umfpack:
      if(this->cols.size() == 0)
      {
        // using int for indices
        umfpack_di_free_symbolic(&symbolic);
        umfpack_di_free_numeric(&numeric);
      }
      else
      {
        // using long for indices
        umfpack_dl_free_symbolic(&symbolic);
        umfpack_dl_free_numeric(&numeric);
      }
      break;
    case DirectSolver::DirectSolverTypes::pardiso:
      ErrThrow("Pardiso does not yet work");
      break;
      default:
      ErrThrow("unknown DirectSolverTypes ",
               static_cast<typename 
                 std::underlying_type<DirectSolverTypes>::type> (type));
  }
  Output::print<3>("destructed a DirectSolver object");
}

/** ************************************************************************ */
void DirectSolver::symetric_factorize()
{
  int n_eq = matrix->GetN_Rows();
  switch(type)
  {
    case DirectSolverTypes::umfpack:
    {
      // symbolic factorization
      if(this->cols.size() == 0)
      {
        // using int for indices
        int error = umfpack_di_symbolic(n_eq, n_eq, matrix->GetRowPtr(),
                                        matrix->GetKCol(), matrix->GetEntries(),
                                        &symbolic, nullptr, nullptr);
        handle_error_umfpack(error);
      }
      else
      {
        // using long for indices
        int error = umfpack_dl_symbolic(n_eq, n_eq, &this->rows[0], 
                                        &this->cols[0], matrix->GetEntries(), 
                                        &symbolic, nullptr, nullptr);
        handle_error_umfpack(error);
      }
      break;
    }
    case DirectSolverTypes::pardiso:
    {
      ErrThrow("Pardiso does not yet work");
      break;
    }
    default:
      ErrThrow("unknown DirectSolverTypes ",
               static_cast<typename 
                 std::underlying_type<DirectSolverTypes>::type> (type));
      break;
  }
}

/** ************************************************************************ */
void DirectSolver::numeric_factorize()
{
  switch(type)
  {
    case DirectSolverTypes::umfpack:
    {
      double Info[UMFPACK_INFO];
      double Control[UMFPACK_CONTROL];
      umfpack_di_defaults(Control);
      
      if(this->cols.size() == 0)
      {
        // using int for indices
        int error = umfpack_di_numeric(matrix->GetRowPtr(), matrix->GetKCol(),
                                       matrix->GetEntries(), symbolic, &numeric,
                                       Control, Info);
        handle_error_umfpack(error);
      }
      else
      {
        // using long for indices
        int error = umfpack_dl_numeric(&this->rows[0], &this->cols[0],
                                       matrix->GetEntries(), symbolic, &numeric,
                                       Control, Info);
        handle_error_umfpack(error);
      }
      break;
    }
    case DirectSolverTypes::pardiso:
    {
      ErrThrow("Pardiso does not yet work");
      break;
    }
    default:
      ErrThrow("unknown DirectSolverTypes ",
               static_cast<typename 
                 std::underlying_type<DirectSolverTypes>::type> (type));
      break;
  }
}

/** ************************************************************************ */
void DirectSolver::solve(const double* rhs, double* solution)
{
  Output::print<3>("solving using a direct solver");
  switch(type)
  {
    case DirectSolverTypes::umfpack:
    {
      // symbolic factorization
      if(this->cols.size() == 0)
      {
        // using int for indices
        int error = umfpack_di_solve(UMFPACK_At, matrix->GetRowPtr(), 
                                     matrix->GetKCol(), matrix->GetEntries(),
                                     solution, rhs, numeric, nullptr, nullptr);
        handle_error_umfpack(error);
      }
      else
      {
        // using long for indices
        int error = umfpack_dl_solve(UMFPACK_At, &this->rows[0], 
                                     &this->cols[0], matrix->GetEntries(),  
                                     solution, rhs, numeric, nullptr, nullptr);
        handle_error_umfpack(error);
      }
      break;
    }
    case DirectSolverTypes::pardiso:
    {
      ErrThrow("Pardiso does not yet work");
      break;
    }
    default:
      ErrThrow("unknown DirectSolverTypes ",
               static_cast<typename 
                 std::underlying_type<DirectSolverTypes>::type> (type));
      break;
  }
}

/** ************************************************************************ */
void DirectSolver::solve(const BlockVector& rhs, BlockVector& solution)
{
  if(  rhs.length() != this->matrix->GetN_Rows() 
    || solution.length() != this->matrix->GetN_Columns())
    ErrThrow("solution or right hand side vector has wrong size. ",
             "Size of the matrix: ", this->matrix->GetN_Rows(), " x ", 
             this->matrix->GetN_Columns(),"\t rhs size: ", rhs.length(), 
             "\tsolution size: ", solution.length());

  solve(rhs.get_entries(), solution.get_entries());
}

/** ************************************************************************ */
/// @note everything below this line is to be deleted

void DirectSolver_old(TSquareMatrix *matrix, double *rhs, double *sol)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  int N_Eqn;
  const int *Row, *KCol;
  double *Values;
  void *Symbolic, *Numeric;

  N_Eqn = matrix->GetN_Columns();
  Row = matrix->GetRowPtr();
  KCol = matrix->GetKCol();
  Values = matrix->GetEntries();

  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
  {
    // sort matrix
    Output::print("umfpack: reordering of the columns will be performed", matrix->GetColOrder());
    Output::print("umfpack: no back ordering implemented !!!");
    matrix->reorderMatrix();
  }
 
  t1 = GetTime();
  ret = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, KCol, Values,
    &Symbolic, NULL, NULL);
  t2 = GetTime();
  // error occured
  if (ret!=0)
  {
    ErrThrow("error in umfpack_di_symbolic ", ret);
  }

  ret = umfpack_di_numeric(Row, KCol, Values, Symbolic,
    &Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic);
  t3 = GetTime();
  // error occured
  if (ret!=0)
  {
    ErrThrow("error in umfpack_di_numeric ", ret);
  }

  ret = umfpack_di_solve(UMFPACK_At, Row, KCol, Values,
    sol, rhs, Numeric, NULL, NULL);
  umfpack_di_free_numeric(&Numeric);
  t4 = GetTime();
  if (ret!=0)
  {
    ErrThrow("error in umfpack_di_solve ", ret);
  }
//  OutPut("umfpack: " << ret << " " << t4-t1 << " sec." << endl);
}

// rb_flag = 0 ==> allocieren und LU-Zerlegung
// rb_flag = 1 ==> nur vorw./rueckw.
// rb_flag = 2 ==> speicher wieder freigeben
// rb_flag = 3 ==> allocieren, LU-Zerl, Freigabe
void DirectSolver_old(TSquareMatrix2D *sqmatrixA11,
                      TSquareMatrix2D *sqmatrixA12,
                      TSquareMatrix2D *sqmatrixA21,
                      TSquareMatrix2D *sqmatrixA22,
                      TMatrix2D *matrixB1T, TMatrix2D *matrixB2T,
                      TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                      double *rhs, double *sol, int rb_flag)
{
  const int *KColA, *RowPtrA;
  const int *KColB, *RowPtrB;
  const int *KColBT, *RowPtrBT;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;
  double *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  int N_, N_U, N_P, N_Entries;
  static double *Entries;
  static int *KCol, *RowPtr;
  double *null = (double *) NULL;
  static void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;
  int verbose = TDatabase::ParamDB->SC_VERBOSE;
  double sum = 0;

  if (rb_flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
    return;
  }

  Output::print<3>("rb_flag: ", rb_flag);
  
  if (rb_flag==0 || rb_flag==3)
  {
    t1 = GetTime();
    N_U = sqmatrixA11->GetN_Rows();
    N_P = matrixB1->GetN_Rows();
    N_ = 2*N_U + N_P;
    N_Active = sqmatrixA11->GetActiveBound();

    KColA = sqmatrixA11->GetKCol();
    RowPtrA = sqmatrixA11->GetRowPtr();

    KColB = matrixB1->GetKCol();
    RowPtrB = matrixB1->GetRowPtr();

    KColBT = matrixB1T->GetKCol();
    RowPtrBT = matrixB1T->GetRowPtr();

    EntriesA11 = sqmatrixA11->GetEntries();
    EntriesA12 = sqmatrixA12->GetEntries();
    EntriesA21 = sqmatrixA21->GetEntries();
    EntriesA22 = sqmatrixA22->GetEntries();

    EntriesB1 = matrixB1->GetEntries();
    EntriesB2 = matrixB2->GetEntries();
    EntriesB1T = matrixB1T->GetEntries();
    EntriesB2T = matrixB2T->GetEntries();

    N_Entries = 4*RowPtrA[N_U] + 2*RowPtrB[N_P] + 2*RowPtrBT[N_U];
    Entries = new double[N_Entries];
    KCol = new int[N_Entries];

    RowPtr = new int[N_+1];
    RowPtr[0] = 0;

    pos = 0;

    for(i=0;i<N_U;i++)
    {
      begin = RowPtrA[i];
      end = RowPtrA[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesA11[j];
        KCol[pos] = KColA[j];
        pos++;

        Entries[pos] = (i<N_Active)?EntriesA12[j]:0;
        KCol[pos] = KColA[j]+N_U;
        pos++;
      }

      if(i<N_Active)
      {
        begin = RowPtrBT[i];
        end = RowPtrBT[i+1];
        for(j=begin;j<end;j++)
        {
          Entries[pos] = EntriesB1T[j];
          KCol[pos] = KColBT[j]+2*N_U;
          pos++;
        }
      }
      RowPtr[i+1] = pos;
    }

    for(i=0;i<N_U;i++)
    {
      begin = RowPtrA[i];
      end = RowPtrA[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = (i<N_Active)?EntriesA21[j]:0;
        KCol[pos] = KColA[j];
        pos++;
        Entries[pos] = EntriesA22[j];
        KCol[pos] = KColA[j]+N_U;
        pos++;
      }

      if(i<N_Active)
      {
        begin = RowPtrBT[i];
        end = RowPtrBT[i+1];
        for(j=begin;j<end;j++)
        {
          Entries[pos] = EntriesB2T[j];
          KCol[pos] = KColBT[j]+2*N_U;
          pos++;
        }
      }
      RowPtr[N_U+i+1] = pos;
    }

    for(i=0;i<N_P;i++)
    {
      begin = RowPtrB[i];
      end = RowPtrB[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB1[j];
        KCol[pos] = KColB[j];
        pos++;
      }
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB2[j];
        KCol[pos] = KColB[j]+N_U;
        pos++;
      }
      RowPtr[2*N_U+i+1] = pos;
    }

    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    {
      // pressure constant
      begin = RowPtr[2*N_U];
      end = RowPtr[2*N_U+1];
      for(j=begin+1;j<end;j++)
        Entries[j] = 0;
      Entries[begin] = 1;
      KCol[begin] = 2*N_U;
      rhs[2*N_U] = 0;
    }

    // sort matrix
    for(i=0;i<N_;i++)
    {
      begin=RowPtr[i];
      end=RowPtr[i+1];

      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];      value = Entries[j];
            KCol[j] = KCol[k]; Entries[j] = Entries[k];
            KCol[k] = l;       Entries[k] = value;
          }                      // endif
        }                        // endfor k
      }                          // endfor j
    }                            // endfor i

    /*
    for(i=0;i<N_;i++)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
        cout << i << " " << KCol[j] << " " << Entries[j] << endl;
    }
    */

    t2 = GetTime();
    ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
    if (ret!=0)
    {
       ErrThrow("WARNING: symbolic: ", ret);
    }
    t3 = GetTime();

    ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
    if (ret!=0)
    {
      ErrThrow("WARNING: numeric: ", ret);
    }
    t4 = GetTime();
    umfpack_di_free_symbolic(&Symbolic);
  }

  t4 = GetTime();
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if (ret!=0)
  {
    ErrThrow("WARNING: solve: ", ret);
  }

  t5 = GetTime();

  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
  }

  Output::print<3>("UMFPACK Time:");
  Output::print<3>("  data prep: ", t2-t1, "s ");
  Output::print<3>("  symbolic: ", t3-t2, "s ");
  Output::print<3>("  numeric: ", t4-t3, "s ");
  Output::print<3>("  solve: ", t5-t4, "s ");
  Output::print<3>("UMFPACK total time: ", t5-t1, "s ");
  

  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}

void DirectSolver_old(TSquareMatrix2D *sqmatrixA, TMatrix2D *matrixB1,
                      TMatrix2D *matrixB2, double *rhs, double *sol,
                      int rb_flag)
{
  const int *KColA, *RowPtrA;
  const int *KColB, *RowPtrB;
  double *EntriesA, *EntriesB1, *EntriesB2;
  int N_, N_U, N_P, N_B, N_Entries;
  static double *Entries;
  static int *KCol, *RowPtr;
  double *null = (double *) NULL;
  static void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5, sum;
  
  if (rb_flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
      Output::print<3>("rb_flag: ", rb_flag);
    return;
  }

  Output::print<3>("rb_flag: ", rb_flag);
  
  if (rb_flag==0 || rb_flag==3)
  {
    t1 = GetTime();
    // get information from the matrices
    // size
    N_U = sqmatrixA->GetN_Rows();
    N_P = matrixB1->GetN_Rows();
    N_ = 2*N_U + N_P;
    N_Active = sqmatrixA->GetActiveBound();
    // pointer to the index arrays
    KColA = sqmatrixA->GetKCol();
    RowPtrA = sqmatrixA->GetRowPtr();

    KColB = matrixB1->GetKCol();
    RowPtrB = matrixB1->GetRowPtr();

    // entries
    EntriesA = sqmatrixA->GetEntries();

    EntriesB1 = matrixB1->GetEntries();
    EntriesB2 = matrixB2->GetEntries();

    // allocate arrays for structure of combined matrix
    // total number of entries
    N_Entries = 2*RowPtrA[N_U] + 4*RowPtrB[N_P];
    Entries = new double[N_Entries];
    KCol = new int[N_Entries];
    RowPtr = new int[N_+1];
    RowPtr[0] = 0;
    N_B = RowPtrB[N_P];

    pos = 0;
    // fill combined matrix
    for(i=0;i<N_U;i++)
    {
      // first velocity component
      begin = RowPtrA[i];
      end = RowPtrA[i+1];
      for(j=begin;j<end;j++)
      {
        // A11
        Entries[pos] = EntriesA[j];
        KCol[pos] = KColA[j];
        pos++;
      }
      // B1T
      if(i<N_Active)
      {
        // this is quit inefficient, think about more efficient solutions
        // later
        // loop over column indices of matrix B1
        for (k=0;k< N_P; k++)
        {
          begin = RowPtrB[k];
          end = RowPtrB[k+1];

          // if column index equal to i
          for(l=begin;l<end;l++)
          {
            if (KColB[l] == i)
            {
              Entries[pos] = EntriesB1[l];
              KCol[pos] = k+2*N_U;
              pos++;
            }
          }
        }
      }
      RowPtr[i+1] = pos;
    }

    // second velocity component
    for(i=0;i<N_U;i++)
    {
      begin = RowPtrA[i];
      end = RowPtrA[i+1];
      for(j=begin;j<end;j++)
      {
        // A22
        Entries[pos] = EntriesA[j];
        KCol[pos] = KColA[j]+N_U;
        pos++;
      }
      // B2T
      if(i<N_Active)
      {
        // this is quit inefficient, think about more efficient solutions
        // later
        // loop over column indices of matrix B1
        for (k=0;k< N_P; k++)
        {
          begin = RowPtrB[k];
          end = RowPtrB[k+1];

          // if column index equal to i
          for(l=begin;l<end;l++)
          {
            if (KColB[l] == i)
            {
              Entries[pos] = EntriesB2[l];
              KCol[pos] = k+2*N_U;
              pos++;
            }
          }
        }
      }
      RowPtr[N_U+i+1] = pos;
    }

    // pressure
    for(i=0;i<N_P;i++)
    {
      // B1
      begin = RowPtrB[i];
      end = RowPtrB[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB1[j];
        KCol[pos] = KColB[j];
        pos++;
      }
      // B2
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB2[j];
        KCol[pos] = KColB[j]+N_U;
        pos++;
      }
      RowPtr[2*N_U+i+1] = pos;
    }

    /*  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    {
      // pressure constant
      begin = RowPtr[2*N_U];
      end = RowPtr[2*N_U+1];
      for(j=begin+1;j<end;j++)
        Entries[j] = 0;
      Entries[begin] = 1;
      KCol[begin] = 2*N_U;
      rhs[2*N_U] = 0;
    }
    */

    // sort matrix
    for(i=0;i<N_;i++)
    {
      begin=RowPtr[i];
      end=RowPtr[i+1];

      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];      value = Entries[j];
            KCol[j] = KCol[k]; Entries[j] = Entries[k];
            KCol[k] = l;       Entries[k] = value;
          }                        // endif
        }                          // endfor k
      }                            // endfor j
    }                              // endfor i

    /*
    for(i=0;i<N_;i++)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
        cout << i << " " << KCol[j] << " " << Entries[j] << endl;
    }
    */

//     t2 = GetTime();

    ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
//     OutPut("symbolic: " << ret << endl);
//     t3 = GetTime();
    ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
//     OutPut("numeric: " << ret << endl);
//     t4 = GetTime();
    umfpack_di_free_symbolic(&Symbolic);
  }

  t4 = GetTime();
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  Output::print<2>("solve: ", ret);
  t5 = GetTime();

  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
  }
  
//   cout << "UMFPACK:";
//   cout << "  data prep: " << t2-t1 << " ";
//   cout << "  symbolic: " << t3-t2 << " ";
//   cout << "  numeric: " << t4-t3 << " ";
//   cout << "  solve: " << t5-t4 << endl;
//   cout << "UMFPACK total time: " << t5-t1 << endl;

  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}

#ifdef __3D__
//****************************************************************************/
//
// for NSTYPE == 4
// flag = 0 ==> allocate and LU factorization, solve
// flag = 1 ==> only forward and backward solve
// flag = 2 ==> free memory
// flag = 3 ==> allocate, LU factorization, solve  and free memory
// flag = 4 ==> free memory
//
//****************************************************************************/
void DirectSolver_old(TSquareMatrix3D *sqmatrixA11,
                      TSquareMatrix3D *sqmatrixA12,
                      TSquareMatrix3D *sqmatrixA13,
                      TSquareMatrix3D *sqmatrixA21,
                      TSquareMatrix3D *sqmatrixA22,
                      TSquareMatrix3D *sqmatrixA23,
                      TSquareMatrix3D *sqmatrixA31,
                      TSquareMatrix3D *sqmatrixA32,
                      TSquareMatrix3D *sqmatrixA33,
                      TMatrix3D *matrixB1T, TMatrix3D *matrixB2T,
                      TMatrix3D *matrixB3T, TMatrix3D *matrixB1,
                      TMatrix3D *matrixB2, TMatrix3D *matrixB3,
                      double *rhs, double *sol, int flag)
{
  Output::print<4>("umf3d");
  int *KColA, *RowPtrA;
  int *KColB, *RowPtrB;
  int *KColBT, *RowPtrBT;
  double *EntriesA11, *EntriesA12, *EntriesA13, *EntriesA21;
  double *EntriesA22, *EntriesA23, *EntriesA31, *EntriesA32, *EntriesA33;
  double *EntriesB1, *EntriesB2,  *EntriesB3, *EntriesB1T, *EntriesB2T, *EntriesB3T;
  int N_, N_U, N_P, N_Entries;
  static double *Entries;
  static int *KCol, *RowPtr;
  double *null = (double *) NULL;
  static void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;
  int verbose = TDatabase::ParamDB->SC_VERBOSE;

  if (flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
    return;
  }

  Output::print<4>("flag: ", flag);
  
  t1 = GetTime();
  if (flag==0 || flag==3)
  {
  N_U = sqmatrixA11->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 3*N_U + N_P;
  N_Active = sqmatrixA11->GetActiveBound();

  KColA = sqmatrixA11->GetKCol();
  RowPtrA = sqmatrixA11->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  KColBT = matrixB1T->GetKCol();
  RowPtrBT = matrixB1T->GetRowPtr();

  EntriesA11 = sqmatrixA11->GetEntries();
  EntriesA12 = sqmatrixA12->GetEntries();
  EntriesA13 = sqmatrixA13->GetEntries();
  EntriesA21 = sqmatrixA21->GetEntries();
  EntriesA22 = sqmatrixA22->GetEntries();
  EntriesA23 = sqmatrixA23->GetEntries();
  EntriesA31 = sqmatrixA31->GetEntries();
  EntriesA32 = sqmatrixA32->GetEntries();
  EntriesA33 = sqmatrixA33->GetEntries();

  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();
  EntriesB3 = matrixB3->GetEntries();
  EntriesB1T = matrixB1T->GetEntries();
  EntriesB2T = matrixB2T->GetEntries();
  EntriesB3T = matrixB3T->GetEntries();

  N_Entries = 9*RowPtrA[N_U] + 3*RowPtrB[N_P] + 3*RowPtrBT[N_U];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];

  RowPtr = new int[N_+1];
  RowPtr[0] = 0;

  pos = 0;

  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      // A11
      Entries[pos] = EntriesA11[j];
      KCol[pos] = KColA[j];
      pos++;
      // A12
      Entries[pos] = (i<N_Active)?EntriesA12[j]:0;
      KCol[pos] = KColA[j]+N_U;
      pos++;
      // A13
      Entries[pos] = (i<N_Active)?EntriesA13[j]:0;
      KCol[pos] = KColA[j]+2*N_U;
      pos++;
    }

    if(i<N_Active)
    {
      // B1T
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB1T[j];
        KCol[pos] = KColBT[j]+3*N_U;
        pos++;
      }
    }
    RowPtr[i+1] = pos;
  }

  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      // A21
      Entries[pos] = (i<N_Active)?EntriesA21[j]:0;
      KCol[pos] = KColA[j];
      pos++;
      // A22
      Entries[pos] = EntriesA22[j];
      KCol[pos] = KColA[j]+N_U;
      pos++;
      // A23
      Entries[pos] = (i<N_Active)?EntriesA23[j]:0;
      KCol[pos] = KColA[j]+2*N_U;
      pos++;
    }

    if(i<N_Active)
    {
      // B2T
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB2T[j];
        KCol[pos] = KColBT[j]+3*N_U;
        pos++;
      }
    }
    RowPtr[N_U+i+1] = pos;
  }

  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      // A31
      Entries[pos] = (i<N_Active)?EntriesA31[j]:0;
      KCol[pos] = KColA[j];
      pos++;
      // A32
      Entries[pos] = (i<N_Active)?EntriesA32[j]:0;
      KCol[pos] = KColA[j]+N_U;
      pos++;
      // A33
      Entries[pos] = EntriesA33[j];
      KCol[pos] = KColA[j]+2*N_U;
      pos++;
    }

    if(i<N_Active)
    {
      // B3T
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB3T[j];
        KCol[pos] = KColBT[j]+3*N_U;
        pos++;
      }
    }
    RowPtr[2*N_U+i+1] = pos;
  }

  for(i=0;i<N_P;i++)
  {
    begin = RowPtrB[i];
    end = RowPtrB[i+1];
    for(j=begin;j<end;j++)
    {
      // B1
      Entries[pos] = EntriesB1[j];
      KCol[pos] = KColB[j];
      pos++;
    }
    for(j=begin;j<end;j++)
    {
      // B2
      Entries[pos] = EntriesB2[j];
      KCol[pos] = KColB[j]+N_U;
      pos++;
    }
    for(j=begin;j<end;j++)
    {
      // B3
      Entries[pos] = EntriesB3[j];
      KCol[pos] = KColB[j]+2*N_U;
      pos++;
    }
    RowPtr[3*N_U+i+1] = pos;
  }

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    // pressure constant
    begin = RowPtr[3*N_U];
    end = RowPtr[3*N_U+1];
    for(j=begin+1;j<end;j++)
      Entries[j] = 0;
    Entries[begin] = 1;
    KCol[begin] = 3*N_U;
    rhs[3*N_U] = 0;
  }
  
  // sort matrix
  for(i=0;i<N_;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_;i++)
  {
    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << i << " " << KCol[j] << " " << Entries[j] << endl;
  }
  */

  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
  if (ret!=0)
  {
    Output::print("WARNING: symbolic: ", ret);
  }
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  if (ret!=0)
  {
    Output::print("WARNING: numeric: ", ret);
  }
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  }
  t4 = GetTime();
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if (ret!=0)
  {
    Output::print("WARNING: solve: ", ret);
  }
  t5 = GetTime();

  if (flag==2 || flag==3)
  {
      umfpack_di_free_numeric(&Numeric);
      delete [] Entries;
      delete [] KCol;
      delete [] RowPtr;
  }

  Output::print<3>("UMFPACK Time:");
  if (flag==0 || flag==3)
  {
    Output::print<3>("  data prep: ", t2-t1, "s ");
    Output::print<3>("  symbolic: ", t3-t2, "s ");
    Output::print<3>("  numeric: ", t4-t3, "s ");
  }
  Output::print<3>("  solve: ", t5-t4, "s ");
  Output::print<3>("UMFPACK total time: ", t5-t1, "s ");
  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}

#endif // __3D__
