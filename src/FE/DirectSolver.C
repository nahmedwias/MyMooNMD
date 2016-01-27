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
#include <ColoredBlockMatrix.h>

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
 : type(type), matrix(matrix), symbolic(nullptr), numeric(nullptr), 
   isFortranShifted(false)
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
  
  this->symetric_factorize();
  this->numeric_factorize();
}

/** ************************************************************************ */
DirectSolver::DirectSolver(const ColoredBlockMatrix& matrix, 
DirectSolver::DirectSolverTypes type)
 : DirectSolver(matrix.get_combined_matrix(), type)
{
}

/** ************************************************************************ */
DirectSolver::~DirectSolver()
{
  switch(type)
  {
    case DirectSolver::DirectSolverTypes::umfpack:
      umfpack_di_free_symbolic(&symbolic);
      umfpack_di_free_numeric(&numeric);
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
      int error = umfpack_di_symbolic(n_eq, n_eq, matrix->GetRowPtr(), 
                                      matrix->GetKCol(), matrix->GetEntries(), 
                                      &symbolic, nullptr, nullptr);
      handle_error_umfpack(error);
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
  
      int error = umfpack_di_numeric(matrix->GetRowPtr(), matrix->GetKCol(),
                                     matrix->GetEntries(), symbolic, &numeric,
                                     Control, Info);
      handle_error_umfpack(error);
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
void DirectSolver::fortranShift()
{
  unsigned int n_eq = this->matrix->GetN_Rows();
  unsigned int n_en = this->matrix->GetN_Entries();
  if (!isFortranShifted)
  {
    for(unsigned int i = 0; i < n_eq+1; i++)
    {
      this->matrix->GetRowPtr()[i] += 1;
    }
    for (unsigned int i = 0; i < n_en; i++)
    {
      this->matrix->GetKCol()[i] += 1;
    }
    isFortranShifted = true;
  }
  else
  {
    for(unsigned int i = 0; i < n_eq+1; i++)
    {
      this->matrix->GetRowPtr()[i] -= 1;
    }
    for (unsigned int i = 0; i < n_en; i++)
    {
      this->matrix->GetKCol()[i] -= 1;
    }
    isFortranShifted = false;
  }
}

/** ************************************************************************ */
void DirectSolver::solve(double* rhs, double* solution)
{
  Output::print<3>("solving using a direct solver");
  switch(type)
  {
    case DirectSolverTypes::umfpack:
    {
      // symbolic factorization
      int error = umfpack_di_solve(UMFPACK_At, matrix->GetRowPtr(), 
                                   matrix->GetKCol(), matrix->GetEntries(),  
                                   solution, rhs, numeric, nullptr, nullptr);
      handle_error_umfpack(error);
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
void DirectSolver::solve(BlockVector& rhs, BlockVector& solution)
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

void UMFPACK_return(int ret)
{

  if (TDatabase::ParamDB->SC_VERBOSE>=2)
  {
      OutPut("solved " << ret << endl);
  }
  if (ret == 0) 
  {
      OutPut("solved sucessfully (" << ret << ")" << endl);
  }
  else if (ret ==   1) 
  {
      OutPut("solved failed (" << ret << "): Matrix singular" << endl);
  }
  else if (ret ==   2)
  {
      OutPut("solved failed (" << ret << "): det(A) != 0 but < eps" << endl);
  }
  else if (ret ==   3)
  {
      OutPut("solved failed (" << ret << "): det(A) != 0 but > inf" << endl);
  }
  else if (ret ==  -1)
  {
      OutPut("solved failed (" << ret << "): Not enough memory!" << endl);
  }
  else if (ret ==  -3)
  {
      OutPut("solved failed (" << ret << "): Used Numeric object is invalided!" << endl);
  }
  else if (ret ==  -4)
  {
      OutPut("solved failed (" << ret << "): Used Symbolic object is invalided!" << endl);
  }
  else if (ret ==  -5)
  {
      OutPut("solved failed (" << ret << "): Argument missing!" << endl);
  }
  else if (ret ==  -6)
  {
      OutPut("solved failed (" << ret << "): Number of rows and columns must be greater 0!" << endl);
  }
  else if (ret ==  -8)
  {
      OutPut("solved failed (" << ret << "): Invalidid matrix!" << endl);
  }
  else if (ret == -11)
  {
      OutPut("solved failed (" << ret << "): Different pattern" << endl);
  }
  else if (ret == -13)
  {
      OutPut("solved failed (" << ret << "): Invalidid system!" << endl);
  }
  else if (ret == -15)
  {
      OutPut("solved failed (" << ret << "): Invalidid permutation!" << endl);
  }
  else if (ret == -17)
  {
      OutPut("solved failed (" << ret << "): Error due to I/O!!!" << endl);
  }
  else if (ret == -911)
  {
      OutPut("solved failed (" << ret << "): Internal error!!!" << endl);
  }

  if (ret != 0) exit(4711);
}

/*******************************************************************/
/*        SCALAR PROBLEMS                                          */
/*******************************************************************/
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol)
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

/*******************************************************************/
/*        SCALAR PROBLEMS                                          */
/*******************************************************************/
void DirectSolver(TSquareMatrix2D *matrix, double *rhs, double *sol)
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
    Output::print("umfpack: reordering of the columns will be performed");
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
  Output::print("umfpack: ", ret, " ", t4-t1, " sec.");
}


// rb_flag = 0 ==> allocation and LU-decomposition forward/backward.
// rb_flag = 1 ==> only forward/backward.
// rb_flag = 2 ==> forward/backward and free up memory
// rb_flag = 3 ==> allocation, LU-decomposition, forward/backward, free up memory
// rb_flag = 4 ==> only free up memory
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, double *&Values,
                   int *&KCol, int *&Row, void *&Symbolic, void *&Numeric, int rb_flag)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  int N_Eqn, N_Entries;
  const int *Row_orig, *KCol_orig;
  double *Values_orig;
//   void *Symbolic, *Numeric;

//   static double *Values;
//   static int *KCol, *Row;
//   static void *Symbolic, *Numeric;  
//     
  double *null = (double *) NULL;

  
  if (rb_flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Values;
    delete [] KCol;
    delete [] Row;
    return;
  }  
  
  Output::print<3>("rb_flag: ", rb_flag);
  
  
 if (rb_flag==0 || rb_flag==3)
 {  
  N_Eqn = matrix->GetN_Columns();
  Row_orig = matrix->GetRowPtr();
  KCol_orig = matrix->GetKCol();
  Values_orig = matrix->GetEntries();
  
  N_Entries = Row_orig[N_Eqn];     
  KCol = new int[N_Entries];
  Row = new int[N_Eqn+1];
  Values = new double[N_Entries];
  
  
  memcpy(Values, Values_orig, N_Entries*SizeOfDouble);
  memcpy(KCol, KCol_orig, N_Entries*SizeOfInt);
  memcpy(Row, Row_orig, (N_Eqn+1)*SizeOfInt);
   
   
  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
  {
    // sort matrix
    Output::print("umfpack: reordering of the columns will be performed");
    Output::print("umfpack: no back ordering implemented !!!");

    for(i=0;i<N_Eqn;i++)
    {
      begin=Row[i];
      end=Row[i+1];
      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];     value = Values[j];
            KCol[j] = KCol[k]; Values[j] = Values[k];
            KCol[k] = l;       Values[k] = value;
          }                      // endif
        }                        // endfor k
      }                          // endfor j
    }                            // endfor i
  }
 
  t1 = GetTime();
  ret = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, KCol, Values, &Symbolic, NULL, NULL);
  t2 = GetTime();
  // error occured
  if (ret!=0)
  {
    ErrThrow("error in umfpack_di_symbolic ", ret);
  }

  ret = umfpack_di_numeric(Row, KCol, Values, Symbolic, &Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic);
  t3 = GetTime();
  // error occured
  if (ret!=0)
  {
    ErrThrow("error in umfpack_di_numeric ", ret);
  }
 } // if (rb_flag==0 || rb_flag==3)
 
  
  ret = umfpack_di_solve(UMFPACK_At, Row, KCol, Values, sol, rhs, Numeric, NULL, NULL);
//   umfpack_di_free_numeric(&Numeric);
  t4 = GetTime();
  if (ret!=0)
  {
    ErrThrow("error in umfpack_di_solve ", ret);
  }
  
  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Values;
    delete [] KCol;
    delete [] Row;
  }
  
//  Output::print<1>("umfpack: ", ret, " ", t4-t1, " sec.");
}




/*******************************************************************/
/*        SCALAR PROBLEMS with multiple rhs                        */
/*******************************************************************/
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs, int N_Rhs_Disp)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  int N_Eqn;
  const int *Row, *KCol;
  double *Values, *Sol, *Rhs;
  void *Symbolic, *Numeric;

  N_Eqn = matrix->GetN_Columns();
  Row = matrix->GetRowPtr();
  KCol = matrix->GetKCol();
  Values = matrix->GetEntries();

  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
  {
    // sort matrix
    Output::print("umfpack: reordering of the columns will be performed");
    Output::print("umfpack: no back ordering implemented !!!");
    matrix->reorderMatrix();
  }

  //t1 = GetTime();
  ret = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, KCol, Values, &Symbolic, NULL, NULL);
  //t2 = GetTime();
  // error occured
  if (ret!=0)
  {
    ErrThrow("error in umfpack_di_symbolic ", ret);
  }

  ret = umfpack_di_numeric(Row, KCol, Values, Symbolic, &Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic);
  //t3 = GetTime();
  // error occured
  if (ret!=0)
  {
    ErrThrow("error in umfpack_di_numeric ", ret);
  }

 for(i=N_Rhs_Disp; i<N_Rhs; i++)
  {
   //Sol = sol[i];
   Sol = sol+i*N_Eqn;
   Rhs = rhs+i*N_Eqn;

   ret = umfpack_di_solve(UMFPACK_At, Row, KCol, Values, Sol, Rhs, Numeric, NULL, NULL);

   //t4 = GetTime();
   if (ret!=0)
   {
     ErrThrow("error in umfpack_di_solve ", ret);
   }
  } // for(i=0; i<N_Rhs; i++)

  umfpack_di_free_numeric(&Numeric);

//  Output::print<1>("umfpack: ", ret, " ", t4-t1, " sec.");
}



// rb_flag = 0 ==> allocation and LU-decomposition forward/backward.
// rb_flag = 1 ==> only forward/backward.
// rb_flag = 2 ==> forward/backward and free up memory
// rb_flag = 3 ==> allocation, LU-decomposition, forward/backward, free up memory
// rb_flag = 4 ==> only free up memory
/*******************************************************************/
/*        SCALAR PROBLEMS with multiple rhs                        */
/*******************************************************************/
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs, int N_Rhs_Disp, double *&Values,
                   int *&KCol, int *&Row, void *&Symbolic, void *&Numeric, int rb_flag)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  int N_Eqn, N_Entries;
  const int *Row_orig, *KCol_orig;
  double *Values_orig, *Sol, *Rhs;
//   void *Symbolic, *Numeric;
 
  if (rb_flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Values;
    delete [] KCol;
    delete [] Row;
    return;
  }    
 
 
  Output::print<3>("rb_flag: ", rb_flag);
  
  N_Eqn = matrix->GetN_Columns(); 
 if (rb_flag==0 || rb_flag==3)
 {   
  Row_orig = matrix->GetRowPtr();
  KCol_orig = matrix->GetKCol();
  Values_orig = matrix->GetEntries();

  N_Entries = Row_orig[N_Eqn];     
  KCol = new int[N_Entries];
  Row = new int[N_Eqn+1];
  Values = new double[N_Entries];
    
  memcpy(Values, Values_orig, N_Entries*SizeOfDouble);
  memcpy(KCol, KCol_orig, N_Entries*SizeOfInt);
  memcpy(Row, Row_orig, (N_Eqn+1)*SizeOfInt);
  
  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
  {
    // sort matrix
    Output::print("umfpack: reordering of the columns will be performed");
    Output::print("umfpack: no back ordering implemented !!!");

    for(i=0;i<N_Eqn;i++)
    {
      begin=Row[i];
      end=Row[i+1];
      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];     value = Values[j];
            KCol[j] = KCol[k]; Values[j] = Values[k];
            KCol[k] = l;       Values[k] = value;
          }                      // endif
        }                        // endfor k
      }                          // endfor j
    }                            // endfor i
  }

  //t1 = GetTime();
  ret = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, KCol, Values, &Symbolic, NULL, NULL);
  //t2 = GetTime();
  // error occured
  if (ret!=0)
  {
    ErrThrow("error in umfpack_di_symbolic ", ret);
  }

  ret = umfpack_di_numeric(Row, KCol, Values, Symbolic, &Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic);
  //t3 = GetTime();
  // error occured
  if (ret!=0)
  {
    ErrThrow("error in umfpack_di_numeric ", ret);
  }
 }//if (rb_flag==0 || rb_flag==3)
 
 
 for(i=N_Rhs_Disp; i<N_Rhs; i++)
  {
   //Sol = sol[i];
   Sol = sol+i*N_Eqn;
   Rhs = rhs+i*N_Eqn;

   ret = umfpack_di_solve(UMFPACK_At, Row, KCol, Values, Sol, Rhs, Numeric, NULL, NULL);

   //t4 = GetTime();
   if (ret!=0)
   {
     ErrThrow("error in umfpack_di_solve ", ret);
   }
  } // for(i=0; i<N_Rhs; i++)

//   umfpack_di_free_numeric(&Numeric);

  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Values;
    delete [] KCol;
    delete [] Row;
  }

//  Output::print<1>("umfpack: ", ret, " ", t4-t1, " sec.");
}



void DirectSolverLong(TSquareMatrix *matrix, double *rhs, double *sol)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  long N_Eqn;
  long *Row, *KCol;
  const int *row, *kcol;
  double *Values;
  void *Symbolic, *Numeric;

  N_Eqn = matrix->GetN_Columns();
  row = matrix->GetRowPtr();
  kcol = matrix->GetKCol();
  Values = matrix->GetEntries();

  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
  {
      // sort matrix
      Output::print("umfpack: reordering of the columns will be performed");
      Output::print("umfpack: no back ordering implemented !!!");

      matrix->reorderMatrix();
  }

  Row = new long[N_Eqn+1];
  for(i=0;i<=N_Eqn;i++)
    Row[i] = row[i];

  end = Row[N_Eqn];
  KCol = new long[end];
  for(i=0;i<end;i++)
    KCol[i] = kcol[i];

  t1 = GetTime();
  ret = umfpack_dl_symbolic(N_Eqn, N_Eqn, Row, KCol, Values,
                          &Symbolic, NULL, NULL);
  t2 = GetTime();
 
  ret = umfpack_dl_numeric(Row, KCol, Values, Symbolic,
                          &Numeric, NULL, NULL);
  umfpack_dl_free_symbolic(&Symbolic);
  t3 = GetTime();
 
  ret = umfpack_dl_solve(UMFPACK_At, Row, KCol, Values,
                       sol, rhs, Numeric, NULL, NULL);
  umfpack_dl_free_numeric(&Numeric);
  delete [] KCol;
  delete [] Row;
  t4 = GetTime();
//   OutPut("solve: " << ret << " " << t4-t3 << endl);
//   OutPut("long umfpack: " << ret << " " << t4-t1 << " sec." << endl);
}


// Solver for PDAE2D
//
// rb_flag = 0 ==> allocieren und LU-Zerlegung
// rb_flag = 1 ==> only forward/backward.
// rb_flag = 2 ==> forward/backward and free up memory
// rb_flag = 3 ==> allocieren, LU-Zerl, Freigabe
void DirectSolver(const TSquareMatrix2D *sqmatrixA11, const TSquareMatrix2D *sqmatrixA12,
                  const TSquareMatrix2D *sqmatrixA21, const TSquareMatrix2D *sqmatrixA22,
double *rhs1, double *rhs2, double *sol1, double *sol2, int rb_flag)
{
  const int *KColA, *RowPtrA;
  const double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;
  double *sol, *rhs;
  int N_, N_U, N_Entries;
  static double *Entries;
  static int *KCol, *RowPtr;
  double *null = (double *) NULL;
  static void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  N_U = sqmatrixA11->GetN_Rows();
  N_ = 2*N_U;

  // copy sol and rhs piece by piece
  sol = new double[N_];
  rhs = new double[N_];
  for (i=0;i<N_U;i++)
  {
    sol[i] = sol1[i];
    rhs[i] = rhs1[i];
  }
  for (i=N_U;i<N_;i++)
  {
    sol[i] = sol2[i-N_U];
    rhs[i] = rhs2[i-N_U];
  }

  Output::print<3>("rb_flag: ", rb_flag);
  if (rb_flag==0 || rb_flag==3)
  {
    t1 = GetTime();
    N_Active = sqmatrixA11->GetActiveBound();

    KColA = sqmatrixA11->GetKCol();
    RowPtrA = sqmatrixA11->GetRowPtr();

    EntriesA11 = sqmatrixA11->GetEntries();
    EntriesA12 = sqmatrixA12->GetEntries();
    EntriesA21 = sqmatrixA21->GetEntries();
    EntriesA22 = sqmatrixA22->GetEntries();

    N_Entries = 4*RowPtrA[N_U];
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

      RowPtr[N_U+i+1] = pos;
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
//   OutPut("solve: " << ret << endl);
  t5 = GetTime();

  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
  }

  // copy sol and rhs piece by piece
  for (i=0;i<N_U;i++)
    sol1[i] = sol[i];
  for (i=N_U;i<N_;i++)
    sol2[i-N_U] = sol[i];

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


// rb_flag = 0 ==> allocieren und LU-Zerlegung
// rb_flag = 1 ==> nur vorw./rueckw.
// rb_flag = 2 ==> speicher wieder freigeben
// rb_flag = 3 ==> allocieren, LU-Zerl, Freigabe
void DirectSolver(const TSquareMatrix2D *sqmatrixA11, const TSquareMatrix2D *sqmatrixA12,
                  const TSquareMatrix2D *sqmatrixA21, const TSquareMatrix2D *sqmatrixA22,
                  const TMatrix2D *matrixB1T, const TMatrix2D *matrixB2T,
                  const TMatrix2D *matrixB1, const  TMatrix2D *matrixB2,
double *rhs, double *sol, int rb_flag)
{
  const int *KColA, *RowPtrA;
  const int *KColB, *RowPtrB;
  const int *KColBT, *RowPtrBT;
  const double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;
  const double *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
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

void DirectSolver(const TSquareMatrix2D *sqmatrixA,
                  const TMatrix2D *matrixB1T, const TMatrix2D *matrixB2T,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
                  const TMatrix2D *matrixC,
double *rhs, double *sol)
{
  const int *KColA, *RowPtrA;
  const int *KColB, *RowPtrB;
  const int *KColBT, *RowPtrBT;
  const int *KColC, *RowPtrC;
  const double *EntriesA, *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T, *EntriesC;
  int N_, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  N_U = sqmatrixA->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 2*N_U + N_P;
  N_Active = sqmatrixA->GetActiveBound();

  KColA = sqmatrixA->GetKCol();
  RowPtrA = sqmatrixA->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  KColBT = matrixB1T->GetKCol();
  RowPtrBT = matrixB1T->GetRowPtr();

  KColC = matrixC->GetKCol();
  RowPtrC = matrixC->GetRowPtr();

  EntriesA = sqmatrixA->GetEntries();
  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();
  EntriesB1T = matrixB1T->GetEntries();
  EntriesB2T = matrixB2T->GetEntries();
  EntriesC = matrixC->GetEntries();

  N_Entries = 2*RowPtrA[N_U] + 2*RowPtrB[N_P] + 2*RowPtrBT[N_U] + RowPtrC[N_P];
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
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j];
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
      Entries[pos] = EntriesA[j];
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

    begin = RowPtrC[i];
    end = RowPtrC[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesC[j];
      KCol[pos] = KColC[j] + 2*N_U;
      pos++;
    }

    RowPtr[2*N_U+i+1] = pos;
  }

  /*
    // check residual
    {
      double *res, val, sum;
      res = new double[N_];
      sum = 0;
      for(i=0;i<N_;i++)
      {
        val = rhs[i];
        begin = RowPtr[i];
        end = RowPtr[i+1];
        for(j=begin;j<end;j++)
          val -= Entries[j]*sol[KCol[j]];
        res[i] = val;
        sum += val*val;
        if(fabs(val)>1e-12)
        {
          cout << setw(5) << i << setw(5) << i-2*N_U << setw(25) << res[i] << setw(25) << rhs[i] << endl;
        }
      }
      cout << "defect: " << sqrt(sum) << endl;
    }
  */

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    // pressure constant
    begin = RowPtr[2*N_U];
    end = RowPtr[2*N_U+1];
    k = begin;
    for(j=begin+1;j<end;j++)
    {
      Entries[j] = 0;
      // check whether col (2*N_U) is already in this row
      if(KCol[j] == 2*N_U) k = j;
    }
    Entries[k] = 1;
    KCol[k] = 2*N_U;
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
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_;i++)
  {
    cout << "row: " << i << " " << RowPtr[i] << ".." << RowPtr[i+1] << endl;
    for(j=RowPtr[i];j<RowPtr[i+1]-1;j++)
      if(KCol[j]>=KCol[j+1]) cout << "Error: " << j << endl;

    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << i << " " << KCol[j] << " " << Entries[j] << endl;
  }
  */

  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
//   OutPut("symbolic: " << ret << endl);
//   t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
//   OutPut("numeric: " << ret << endl);
//   t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries, sol, rhs, Numeric, null, null);
//   OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

//   cout << "UMFPACK:";
//   cout << "  data prep: " << t2-t1 << " ";
//   cout << "  symbolic: " << t3-t2 << " ";
//   cout << "  numeric: " << t4-t3 << " ";
//   cout << "  solve: " << t5-t4 << endl;
//   cout << "UMFPACK total time: " << t5-t1 << endl;
}


//****************************************************************************/
//
// for NSTYPE == 1
//
//****************************************************************************/

void DirectSolver(const TSquareMatrix2D *sqmatrixA,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
                  double *rhs, double *sol)
{
  const int *KColA, *RowPtrA;
  const int *KColB, *RowPtrB;
  const double *EntriesA, *EntriesB1, *EntriesB2;
  int N_, N_U, N_P, N_B, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5, sum;

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

  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
//   OutPut("symbolic: " << ret << endl);
//   t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
//   OutPut("numeric: " << ret << endl);
//   t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries, sol, rhs, Numeric, null, null);
//   OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

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

void DirectSolver(const TSquareMatrix2D *sqmatrixA,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
                  double *rhs, double *sol, int rb_flag)
{
  const int *KColA, *RowPtrA;
  const int *KColB, *RowPtrB;
  const double *EntriesA, *EntriesB1, *EntriesB2;
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

//****************************************************************************/
//
// for NSTYPE == 2
//
//****************************************************************************/

void DirectSolver(const TSquareMatrix2D *sqmatrixA,
                  const TMatrix2D *matrixB1T, const  TMatrix2D *matrixB2T,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
double *rhs, double *sol)
{
  const int *KColA, *RowPtrA;
  const int *KColB, *RowPtrB;
  const int *KColBT, *RowPtrBT;
  const double *EntriesA, *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  int N_, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  N_U = sqmatrixA->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 2*N_U + N_P;
  N_Active = sqmatrixA->GetActiveBound();

  KColA = sqmatrixA->GetKCol();
  RowPtrA = sqmatrixA->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  KColBT = matrixB1T->GetKCol();
  RowPtrBT = matrixB1T->GetRowPtr();

  EntriesA = sqmatrixA->GetEntries();
  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();
  EntriesB1T = matrixB1T->GetEntries();
  EntriesB2T = matrixB2T->GetEntries();

  N_Entries = 2*RowPtrA[N_U] + 2*RowPtrB[N_P] + 2*RowPtrBT[N_U];
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
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j];
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
      Entries[pos] = EntriesA[j];
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

//   t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
//   OutPut("symbolic: " << ret << endl);
//   t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
//   OutPut("numeric: " << ret << endl);
//   t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
//   OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
//   t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

//   cout << "UMFPACK:";
//   cout << "  data prep: " << t2-t1 << " ";
//   cout << "  symbolic: " << t3-t2 << " ";
//   cout << "  numeric: " << t4-t3 << " ";
//   cout << "  solve: " << t5-t4 << endl;
//   cout << "UMFPACK total time: " << t5-t1 << endl;
}

void DirectSolver(const TSquareMatrix2D *sqmatrixA,
                  const TMatrix2D *matrixB1T, const TMatrix2D *matrixB2T,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
double *rhs, double *sol, int rb_flag)
{
  const int *KColA, *RowPtrA;
  const int *KColB, *RowPtrB;
  const int *KColBT, *RowPtrBT;
  const double *EntriesA, *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  int N_, N_U, N_P, N_Entries;
  static double *Entries;
  static int *KCol, *RowPtr;
  double *null = (double *) NULL;
  static void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  Output::print<3>("rb_flag: ", rb_flag);
  if (rb_flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;

    return;
  }

  if (rb_flag==0 || rb_flag==3)
  {
    t1 = GetTime();
    N_U = sqmatrixA->GetN_Rows();
    N_P = matrixB1->GetN_Rows();
    N_ = 2*N_U + N_P;
    N_Active = sqmatrixA->GetActiveBound();

    KColA = sqmatrixA->GetKCol();
    RowPtrA = sqmatrixA->GetRowPtr();

    KColB = matrixB1->GetKCol();
    RowPtrB = matrixB1->GetRowPtr();

    KColBT = matrixB1T->GetKCol();
    RowPtrBT = matrixB1T->GetRowPtr();

    EntriesA = sqmatrixA->GetEntries();
    EntriesB1 = matrixB1->GetEntries();
    EntriesB2 = matrixB2->GetEntries();
    EntriesB1T = matrixB1T->GetEntries();
    EntriesB2T = matrixB2T->GetEntries();

    N_Entries = 2*RowPtrA[N_U] + 2*RowPtrB[N_P] + 2*RowPtrBT[N_U];
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
        Entries[pos] = EntriesA[j];
        KCol[pos] = KColA[j];
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
        Entries[pos] = EntriesA[j];
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

    /*
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

    t2 = GetTime();

    ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
    UMFPACK_return(ret);

    t3 = GetTime();
    ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
    UMFPACK_return(ret);

    t4 = GetTime();
    umfpack_di_free_symbolic(&Symbolic);
  }

  t4 = GetTime();
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  UMFPACK_return(ret);

  t5 = GetTime();

  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
  }
  
  Output::print<2>("UMFPACK time:");
  Output::print<2>("  data prep : ", t2-t1);
  Output::print<2>("  symbolic: ", t3-t2);
  Output::print<2>("  numeric: ", t4-t3);
  Output::print<2>("  solve: ", t5-t4);
  Output::print<2>("UMFPACK total time: ", t5-t1);
}

//****************************************************************************/
//
// for NSTYPE == 4
//
//****************************************************************************/

void DirectSolver(const TSquareMatrix2D *sqmatrixA11, const TSquareMatrix2D *sqmatrixA12,
                  const TSquareMatrix2D *sqmatrixA21, const TSquareMatrix2D *sqmatrixA22,
                  const TMatrix2D *matrixB1T, const TMatrix2D *matrixB2T,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
double *rhs, double *sol)
{
  const int *KColA, *RowPtrA;
  const int *KColB, *RowPtrB;
  const int *KColBT, *RowPtrBT;
  const double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;
  const double *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  int N_, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

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
  if(ret!=0)
    Output::print("symbolic: ", ret);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  if(ret!=0)  
    Output::print("numeric: ", ret);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if(ret!=0)    
    Output::print("solve: ", ret);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();
// exit(0);
  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

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


//****************************************************************************/
//
// for NSTYPE == 14
//
//****************************************************************************/

void DirectSolver(const TSquareMatrix2D *sqmatrixA11, const TSquareMatrix2D *sqmatrixA12,
                  const TSquareMatrix2D *sqmatrixA21, const TSquareMatrix2D *sqmatrixA22,
                  const TSquareMatrix2D *sqmatrixC,
                  const TMatrix2D *matrixB1T, const TMatrix2D *matrixB2T,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
                  double *rhs, double *sol)
{
  const int *KColA, *RowPtrA, *KColC, *RowPtrC;
  const int *KColB, *RowPtrB;
  const int *KColBT, *RowPtrBT;
  const double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22, *EntriesC;
  const double *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  int N_, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  // get information from the matrices
  // size
  N_U = sqmatrixA11->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 2*N_U + N_P;
  N_Active = sqmatrixA11->GetActiveBound();
  // pointer to the index arrays
  KColA = sqmatrixA11->GetKCol();
  RowPtrA = sqmatrixA11->GetRowPtr();

  KColC = sqmatrixC->GetKCol();
  RowPtrC = sqmatrixC->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  KColBT = matrixB1T->GetKCol();
  RowPtrBT = matrixB1T->GetRowPtr();
  // entries
  EntriesA11 = sqmatrixA11->GetEntries();
  EntriesA12 = sqmatrixA12->GetEntries();
  EntriesA21 = sqmatrixA21->GetEntries();
  EntriesA22 = sqmatrixA22->GetEntries();
  EntriesC = sqmatrixC->GetEntries();

  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();
  EntriesB1T = matrixB1T->GetEntries();
  EntriesB2T = matrixB2T->GetEntries();

  // allocate arrays for structure of combined matrix
  // total number of entries
  N_Entries = 4*RowPtrA[N_U] + RowPtrC[N_P] + 2*RowPtrB[N_P] + 2*RowPtrBT[N_U];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];
  RowPtr = new int[N_+1];
  RowPtr[0] = 0;

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
      Entries[pos] = EntriesA11[j];
      KCol[pos] = KColA[j];
      pos++;
      // A12
      Entries[pos] = (i<N_Active)?EntriesA12[j]:0;
      KCol[pos] = KColA[j]+N_U;
      pos++;
    }
    // B1T
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
  // second velocity component
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
    }
    // B2T
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
    // C
    begin = RowPtrC[i];
    end = RowPtrC[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesC[j];
      KCol[pos] = KColC[j]+2*N_U;
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



  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
  Output::print<1>("symbolic: ", ret);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  Output::print<1>("numeric: ", ret);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  Output::print<1>("solve: ", ret);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete  []RowPtr;

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

void DirectSolver(const TSquareMatrix2D *sqmatrixA, const TSquareMatrix2D *sqmatrixC,
                  const TMatrix2D *matrixBT, const TMatrix2D *matrixB,
                  double *rhs, double *sol)
{
  const int *KColA, *RowPtrA, *KColC, *RowPtrC;
  const int *KColB, *RowPtrB;
  const int *KColBT, *RowPtrBT;
  const double *EntriesA, *EntriesC, *EntriesB, *EntriesBT;
  int N_DOF, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  // get information from the matrices
  // size
  N_U = sqmatrixA->GetN_Rows();
  N_P = matrixB->GetN_Rows();
  N_DOF = N_U + N_P;
  N_Active = sqmatrixA->GetActiveBound();
  
  // pointer to the index arrays
  KColA = sqmatrixA->GetKCol();
  RowPtrA = sqmatrixA->GetRowPtr();

  KColC = sqmatrixC->GetKCol();
  RowPtrC = sqmatrixC->GetRowPtr();

  KColB = matrixB->GetKCol();
  RowPtrB = matrixB->GetRowPtr();

  KColBT = matrixBT->GetKCol();
  RowPtrBT = matrixBT->GetRowPtr();
  // entries
  EntriesA = sqmatrixA->GetEntries();
  EntriesC = sqmatrixC->GetEntries();

  EntriesB = matrixB->GetEntries();
  EntriesBT = matrixBT->GetEntries();
  
  // allocate arrays for structure of combined matrix
  // total number of entries
  N_Entries = RowPtrA[N_U] + RowPtrC[N_P] + RowPtrB[N_P] + RowPtrBT[N_U];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];
  RowPtr = new int[N_DOF+1];
  RowPtr[0] = 0;

  pos = 0;
  // fill combined matrix
  for(i=0;i<N_U;i++)
  {
    // first velocity component
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      // A
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j];
      pos++;
    }
    // BT
    //if(i<N_Active)
    {
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesBT[j];
        KCol[pos] = KColBT[j]+N_U;
        pos++;
      }
    }
    RowPtr[i+1] = pos;
  }
  // pressure
  for(i=0;i<N_P;i++)
  {
    // B
    begin = RowPtrB[i];
    end = RowPtrB[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB[j];
      KCol[pos] = KColB[j];
      pos++;
    }
    // C
    begin = RowPtrC[i];
    end = RowPtrC[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesC[j];
      KCol[pos] = KColC[j]+N_U;
      pos++;
    }
    RowPtr[N_U+i+1] = pos;
  }

  // sort matrix
  for(i=0;i<N_DOF;i++)
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
  for(i=0;i<N_DOF;i++)
  {
    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << "A(" << i+1 << "," << KCol[j]+1 << ") = " << Entries[j] << ";\n";
  }
  for(i=0;i<N_DOF;i++)
  {
    OutPut("rhs(" << i+1 << ") = " << rhs[i] << endl);
  }
  */
  
  t2 = GetTime();

  ret = umfpack_di_symbolic(N_DOF, N_DOF, RowPtr, KCol, Entries, &Symbolic, 
                            null, null);
  if(ret!=0)
    Output::print("symbolic: ", ret);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, 
                           null);
  if(ret!=0)
    Output::print("numeric: ", ret);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if(ret!=0)
    Output::print("solve: ", ret);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

  Output::print<3>("UMFPACK Time:");
  Output::print<3>("  data prep: ", t2-t1, "s ");
  Output::print<3>("  symbolic: ", t3-t2, "s ");
  Output::print<3>("  numeric: ", t4-t3, "s ");
  Output::print<3>("  solve: ", t5-t4, "s ");
  Output::print<3>("UMFPACK total time: ", t5-t1, "s ");
}


#ifdef __3D__
//****************************************************************************/
//
// for NSTYPE == 2
//
//****************************************************************************/

void DirectSolver(TSquareMatrix3D *sqmatrixA,
TMatrix3D *matrixB1T, TMatrix3D *matrixB2T,
TMatrix3D *matrixB3T,
TMatrix3D *matrixB1,  TMatrix3D *matrixB2,
TMatrix3D *matrixB3,
double *rhs, double *sol)
{
  int *KColA, *RowPtrA;
  int *KColB, *RowPtrB;
  int *KColBT, *RowPtrBT;
  double *EntriesA;
  double *EntriesB1, *EntriesB2, *EntriesB3;
  double *EntriesB1T, *EntriesB2T, *EntriesB3T;
  int N_, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  N_U = sqmatrixA->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 3*N_U + N_P;
  N_Active = sqmatrixA->GetActiveBound();

  KColA = sqmatrixA->GetKCol();
  RowPtrA = sqmatrixA->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  KColBT = matrixB1T->GetKCol();
  RowPtrBT = matrixB1T->GetRowPtr();

  EntriesA = sqmatrixA->GetEntries();

  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();
  EntriesB3 = matrixB3->GetEntries();

  EntriesB1T = matrixB1T->GetEntries();
  EntriesB2T = matrixB2T->GetEntries();
  EntriesB3T = matrixB3T->GetEntries();

  N_Entries = 3*RowPtrA[N_U] + 3*RowPtrB[N_P] + 3*RowPtrBT[N_U];
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
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j];
      pos++;
    }

    if(i<N_Active)
    {
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
      Entries[pos] = EntriesA[j];
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
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j]+2*N_U;
      pos++;
    }

    if(i<N_Active)
    {
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
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB3[j];
      KCol[pos] = KColB[j]+2*N_U;
      pos++;
    }
    RowPtr[3*N_U+i+1] = pos;
  }

double len = 0;
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    // pressure constant
    begin = RowPtr[3*N_U];
    end = RowPtr[3*N_U+1];
    
    len=end- begin;
//     cout << "No. Internal Pressure row entries " <<len-- <<endl;
    
    for(j=begin+1;j<end;j++)
      Entries[j] = 0;
    Entries[begin] = 1;
    KCol[begin] = 3*N_U;
    rhs[3*N_U] = 0;
  }

//   cout << "Total entries : " <<  RowPtr[N_]-len << endl;
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


//   for(i=0;i<N_;i++)
//   {
//     cout << "=====" << RowPtr[i] << " " << RowPtr[i+1] << endl;
//     for(j=RowPtr[i];j<RowPtr[i+1];j++)
//       cout << i << " " << KCol[j] << " " << Entries[j] << endl;
//   }

//   for(i=0;i<N_U;i++)
//   cout << "Rhs " << i << " " << rhs[i] <<" " << rhs[i+1] <<  " " << rhs[i+2] << endl;

//   for(i=0;i<N_P;i++)
//   cout << "Rhs " << i << " " << rhs[3*N_U+i] << endl;
// 

// exit(0);
  

  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
  if(ret!=0)
    Output::print("symbolic: ", ret);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  if(ret!=0)
    Output::print("numeric: ", ret);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if(ret!=0)
    Output::print("solve: ", ret);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

//   cout << "UMFPACK:";
//   cout << "  data prep: " << t2-t1 << " ";
//   cout << "  symbolic: " << t3-t2 << " ";
//   cout << "  numeric: " << t4-t3 << " ";
//   cout << "  solve: " << t5-t4 << endl;
//   cout << "UMFPACK total time: " << t5-t1 << endl;
}


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

void DirectSolver(TSquareMatrix3D *sqmatrixA11, TSquareMatrix3D *sqmatrixA12,
TSquareMatrix3D *sqmatrixA13,
TSquareMatrix3D *sqmatrixA21, TSquareMatrix3D *sqmatrixA22,
TSquareMatrix3D *sqmatrixA23,
TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
TSquareMatrix3D *sqmatrixA33,
TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T,
TMatrix3D *matrixB1,  TMatrix3D *matrixB2, TMatrix3D *matrixB3,
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

void DirectSolver(TSquareMatrix3D **sqmatrices, int n_row, int n_column,
                  double *sol, double *rhs)
{
  int N_Active;
  int N_Rows;
  int N_Entries, N_Row, begin, end, pos, l;
  int *RowPtr, *rowptr;
  int *KCol, *kcol;
  double *Entries, *entries, value;
  
  N_Active = sqmatrices[0]->GetActiveBound();
  N_Rows   = sqmatrices[0]->GetN_Rows();
  
  N_Entries = n_row*n_column*sqmatrices[0]->GetN_Entries();
  
  Entries = new double [N_Entries];
  RowPtr  = new int [N_Rows*n_row+1];
  KCol    = new int [N_Entries];
  
  N_Row = N_Rows*n_row;
  
  pos = 0;
  RowPtr[0] = 0;
  
  for (int i=0;i<n_row;++i)
  {
    for (int row=0;row<N_Rows;++row)
    {
      for (int j=0;j<n_column;++j)
      {
        if ( i != j && row >= N_Active ) continue;
 
        rowptr = sqmatrices[i*n_column+j]->GetRowPtr();
        kcol   = sqmatrices[i*n_column+j]->GetKCol();
        entries = sqmatrices[i*n_column+j]->GetEntries();
 
        begin = rowptr[row];
        end   = rowptr[row+1];

        for (int loc_pos=begin;loc_pos<end;++loc_pos)
         {
          Entries[pos] = entries[loc_pos];
          KCol[pos] = kcol[loc_pos] + j*N_Rows;
          ++pos;
         }
      }
      RowPtr[i*N_Rows+row+1] = pos;
    }
  }
  
   // sort matrix
  for(int i=0;i<N_Row;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(int j=begin;j<end;j++)
    {
      for(int k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }       
  
  void *Symbolic, *Numeric;
  int ret;
  
  ret = umfpack_di_symbolic(N_Row, N_Row, RowPtr, KCol, Entries, &Symbolic, NULL, NULL);
  if ( ret != 0) 
  {
    ErrThrow("symbolic: ", ret);
  }
//   t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, NULL, NULL);
  if ( ret != 0 )
  {
    ErrThrow("numeric: ", ret);
  }
//   t4 = GetTime();
  
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, NULL, NULL);
  if ( ret != 0) 
  {
    ErrThrow("solve: ", ret);
  }
  umfpack_di_free_symbolic(&Symbolic);
  umfpack_di_free_numeric(&Numeric);
  
  delete [] Entries;
  delete [] RowPtr;
  delete [] KCol;
}


#endif
