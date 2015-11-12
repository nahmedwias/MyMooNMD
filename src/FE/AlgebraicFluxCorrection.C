#include <Database.h>
#include <LinAlg.h>
#include <Solver.h>
//#ifdef __2D__
#include <SquareMatrix2D.h>
#include <DiscreteForm2D.h>
#include <FEFunction2D.h>
#include <FEDatabase2D.h>
#include <FE2D.h>
//#endif
//#ifdef __3D__
//#include <SquareMatrix3D.h>
//#include <DiscreteForm3D.h>
//#include <FEFunction3D.h>
//#include <FEDatabase3D.h>
//#include <FE3D.h>
//#endif
#include <AlgebraicFluxCorrection.h>
#include <IsoBoundEdge.h>
#include <BoundComp.h>

#include <stdlib.h>
#include <string.h>
#include <MooNMD_Io.h>

//**************************************************************
//compute the lumped matrix
//output is a vector
//**************************************************************
//#ifdef __2D__
void AlgebraicFluxCorrection::LumpMassMatrixToVector(TSquareMatrix2D *M, double *lump_mass)
//#endif
//#ifdef __3D__
//void LumpMassMatrixToVector(TSquareMatrix3D *M, double *lump_mass)
//#endif
{
  double *Entries;
  int *RowPtr, i, j, rows, j0, j1;

  RowPtr        = M->GetRowPtr();
  Entries       = M->GetEntries();
  rows          = M->GetN_Rows();

  memset(lump_mass, 0, rows*SizeOfDouble);
  for (i=0; i<rows; i++)
  {
    lump_mass[i]=0.0;
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];

    for (j=j0;j<j1;j++)
      lump_mass[i] += Entries[j];

    if(lump_mass[i]==0)
    {
      OutPut("zero entry in lumped matrix "<< i << " "  << lump_mass[i] << endl);
      exit(4711);
    }
  }
}


/*******************************************************************************/
//
// FEM_TVD_ForConvDiff for steady-state cdr equations
// following D. Kuzmin (2007)
//
// inputs : *sqmatrix         - system matrix
//          N_U               - total number of dof
//          N_Active          - active dof, should be the same as N_U
//          *matrix_D_Entries - entries of matrix D
//          *sol              - current solution
//          N_neum_to_diri    - number of Dirichlet dof which are treated as
//                              Neumann dof, MUST BE ORDERED !!!
//          *neum_to_diri     - array which contains the indices of Dirichlet
//                              dof which are treated as Neumann dof
//          compute_matrix_D  - flag which says if matrix_D_Entries should be
//                              computed
//
// outputs: *rhs              - right hand side of FEM-TVD problem
//          *matrix_D_Entries - entries of matrix D if compute_matrix_D == 1
//
/*******************************************************************************/
//#ifdef __2D__
void AlgebraicFluxCorrection::FEM_TVD_ForConvDiff(TSquareMatrix2D *sqmatrix, int N_U, int N_Active,
double *matrix_D_Entries, double *sol, double *rhs,
int N_neum_to_diri, int *neum_to_diri, int compute_matrix_D)
//#endif
//#ifdef __3D__
//void FEM_TVD_ForConvDiff(TSquareMatrix3D *sqmatrix, int N_U, int N_Active,
//double *matrix_D_Entries, double *sol, double *rhs,
//int N_neum_to_diri, int *neum_to_diri, int compute_matrix_D)
//#endif
{
  int  *RowPtr, *ColInd, N_Entries;
  int i,j,j0,j1,j2,j3,jj,index;
  double nenner, zaehler;
  double *Entries, *F;
  double *P_plus, *P_minus, *Q_plus, *Q_minus, *R_plus, *R_minus;

  if (N_Active < N_U)
  {
    OutPut("N_Active < N_U ("<< N_Active<< "<"<< N_U <<
      ") !!! FOR APPLYING ALGEBRAIC FLUX CORRECTION, THE BOUNDARY CONDITIONS SHOULD BE NEUMANN !!!" << endl);
    exit(4711);
  }
  // get pointers to columns, rows and entries of matrix A
  ColInd = sqmatrix->GetKCol();
  RowPtr = sqmatrix->GetRowPtr();
  Entries = sqmatrix->GetEntries();
  N_Entries = sqmatrix->GetN_Entries();

  // allocate memory for array F
  F = new double[N_Entries+6*N_U];
  memset(F, 0, (N_Entries+6*N_U)*SizeOfDouble);
  P_plus = F + N_Entries;
  P_minus = P_plus + N_U;
  Q_plus = P_minus + N_U;
  Q_minus = Q_plus + N_U;
  R_plus = Q_minus + N_U;
  R_minus = R_plus + N_U;


	if(compute_matrix_D){
		computeArtificialDiffusionMatrix(*sqmatrix, matrix_D_Entries,N_U);
	}


  // add this matrix to A giving \tilde A (Entries)
  // this is the matrix with the properties of an M matrix
  Daxpy(N_Entries, 1.0, matrix_D_Entries, Entries);

  // compute matrix F
  // loop over all rows
  for(i=0;i<N_Active;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];

    for(j=j0;j<j1;j++)
    {
      // column
      index = ColInd[j];
      // d_ij (u_i - u_j)
      F[j] = matrix_D_Entries[j] * (sol[index]-sol[i]);
    }
  }
  // matrix F is computed

  // compute flux limiters
  // loop over all rows
  for(i=0;i<N_Active;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    for(j=j0;j<j1;j++)
    {
      if (Entries[j] > 0)
        continue;
      // column
      index = ColInd[j];
      // check transposed entry -> jj
      // diagonal
      if (index==i)
        continue;
      j2 = RowPtr[index];
      j3 = RowPtr[index+1];
      for (jj=j2;jj<j3;jj++)
      {
        if (ColInd[jj]==i)
        {
          //OutPut(Entries[j] << " " << Entries[jj] <<endl);
          break;
        }
      }
      // check upwind condition
      // this ensures that the 'link' between i and index is treated only once
      if (Entries[jj] > Entries[j])
        continue;
      // only the active part of the matrix
      if (F[j] > 0)
      {
        P_plus[i] += F[j];
        if (index<N_U)
          Q_plus[index] += F[j];
        Q_minus[i] -= F[j];
      }
      if (F[j] < 0)
      {
        P_minus[i] += F[j];
        Q_plus[i] -= F[j];
        if (index<N_U)
          Q_minus[index] +=  F[j];
      }
    }                                             // end loop j
  }

  // apply the nodal correction factor evaluated at the upwind node i
  // loop over all nodes
  if (TDatabase::ParamDB->P9!=4711) //CB FIXME Replace that by some okay parameter!
  { // original but discontinuous proposal
    for(i=0;i<N_U;i++)
    {
      if (fabs(P_plus[i])>0)
      {
        R_plus[i] = Q_plus[i]/P_plus[i];
        if (R_plus[i] >1)
          R_plus[i] = 1;
      }
      if (fabs(P_minus[i])>0)
      {
        R_minus[i] = Q_minus[i]/P_minus[i];
        if (R_minus[i] >1)
          R_minus[i] = 1;
      }
      //OutPut(" P " << P_plus[i] << " " <<  P_minus[i] << " ");
      //OutPut(R_plus[i] << " " <<  R_minus[i] << endl);
    }
  }
  else
  {
    // continuous proposal

    for(i=0;i<N_U;i++)
    {
      zaehler =  Q_plus[i];
      if (-Q_minus[i] < zaehler)
        zaehler = -Q_minus[i];
      nenner = 1e-32;
      if (P_plus[i] > nenner)
        nenner = P_plus[i];
      if (-P_minus[i] > nenner)
        nenner = -P_minus[i];
      //OutPut(zaehler << " " << nenner );
      R_plus[i] = zaehler/nenner;
      if (R_plus[i] > 1)
        R_plus[i] = 1;
      R_minus[i] = R_plus[i];
      //OutPut("new " << i << " P " << P_plus[i] << " "  << P_minus[i] <<" Q " << Q_plus[i] << " "  << Q_minus[i] <<
      //" R " << R_plus[i] << " " << R_minus[i] << endl);
    }
  }


  // treat Dirichlet nodes
  for (j=0;j<N_neum_to_diri;j++)
  {
    i=neum_to_diri[j];
    R_plus[i] = 1;
    R_minus[i] = 1;
  }

  // apply the flux limiters
  // loop over all rows
  for(i=0;i<N_Active;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    for(j=j0;j<j1;j++)
    {
      // column
      index = ColInd[j];
      if (index==i)
        continue;
       // this should not happen
      if (Entries[j] > 0)
      {
        OutPut("positive entry in FEMTVD " << i << " " << j << " " << Entries[j] << endl);
        exit(4711);
        continue;
      }

      // check transposed entry
      j2 = RowPtr[index];
      j3 = RowPtr[index+1];
      for (jj=j2;jj<j3;jj++)
      {
        if (ColInd[jj]==i)
        {
          //OutPut(Entries[j] << " " << Entries[jj] <<endl);
          break;
        }
      }

      if (TDatabase::ParamDB->P8!=4711) //CB FIXME Replace that by some okay parameter!
      {
        // original, symmetric application
        // check upwind condition
        // this ensures that the 'link' between i and index is treated only once
        if (Entries[jj] > Entries[j])
          continue;
        //OutPut(R_plus[i] << " " << R_minus[i] << " : " << R_plus[index] << " " << R_minus[index] << "::");
        // compute contribution to rhs
        if (F[j] > 0)
        {
          F[j] = R_plus[i]*F[j];
        }
        else
          F[j] = R_minus[i]*F[j];
        // update rhs of current row
        rhs[i] += F[j];
        // update rhs wrt to current column
        // note that F[j] = -F[jj] and alpha_j = alpha_jj (symmetry of alpha matrix)
        if (index<N_Active)
          rhs[index] -= F[j];
      }
      else
      { // nonsymmetric application
        // compute contribution to rhs
        if (F[j] > 0)
        {
          F[j] = R_plus[i]*F[j];
        }
        else
          F[j] = R_minus[i]*F[j];
        // update rhs of current row
        rhs[i] += F[j];
      }
    }
  }
  delete [] F;
}

//**************************************************************
//compute system matrix for FEM-FCT
//output is a vector
//**************************************************************
//#ifdef __2D__
void AlgebraicFluxCorrection::FEM_FCT_SystemMatrix(TSquareMatrix2D *M_C, TSquareMatrix2D *A,
double *lump_mass,int N_U)
//#endif
//#ifdef __3D__
//void FEM_FCT_SystemMatrix(TSquareMatrix3D *M_C, TSquareMatrix3D *A,
//double *lump_mass,int N_U)
//#endif
{
  //int *ColInd, *RowPtr, N_Entries, *ColInd_M, *RowPtr_M, N_Entries_M;
  int i,j,j0,j1,index;
  //double *Entries,*Entries_M, cfl;
  double delta_t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1; //appears in the system matrix
  double theta2 = TDatabase::TimeDB->THETA2; //needed only for checking cfl (or checking, whether cfl has to be checked...)

//  // get pointers to columns, rows and entries of matrix M_C
  int* ColInd_M = M_C->GetKCol();
  int* RowPtr_M = M_C->GetRowPtr();
  double*  Entries_M = M_C->GetEntries();
//  N_Entries_M = M_C->GetN_Entries();

  // get pointers to columns, rows and entries of matrix A
  //ColInd = A->GetKCol();
  //RowPtr = A->GetRowPtr();
  double* Entries = A->GetEntries();
  //N_Entries = A->GetN_Entries();

  for(i=0;i<N_U;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr_M[i];
    j1 = RowPtr_M[i+1];

    for(j=j0;j<j1;j++)
    {
      index = ColInd_M[j];
      if(i == index)
      {
        //diagonale
        Entries_M[j] = lump_mass[i]+delta_t*theta1*Entries[j];
        double cfl = lump_mass[i]/Entries[j];
        if (theta2 >0)
        {
          cfl /= theta2;
          if (delta_t > cfl)
          {
            OutPut(lump_mass[i] << " " << "cfl violated: cfl " << cfl << " delta t:" << delta_t <<endl);
          }
        }
      }
      else
      {
        Entries_M[j]= delta_t*theta1*Entries[j];
        if ( Entries_M[j]>1e-10)
          OutPut( Entries_M[j] << " ");

      }
    }
  }
}


//**************************************************************
// MINMOD prelimiter
//**************************************************************
double AlgebraicFluxCorrection::MinMod(double a, double b)
{
  if (a*b < 0)
  {
    return 0.0;
  }
  //CB 2015/07/27 should not this return a negative value, when a<0 and b<0?
  // int sign = a < 0 ? -1 : 1;

  if (fabs(a) <  fabs(b))
  {
    return a; //return sign * a;
  }
  else
  {
    return b; //return sign * b;
  }
}


/*******************************************************************************/
//
// FCT-FEM algorithm
// following D. Kuzmin, M. M"oller (2005) (nonlinear scheme)
//           D. Kuzmin (2008) (linear scheme)
//
// inputs : M_C              - consistent mass matrix
//          A                - stiffness matrix
//          N_U              - number of all unknowns
//          N_Active         - number of active unknows (should be the same
//                             as N_U !!!)
//          lump_mass        - lumped mass matrix, stored as array
//          matrix_D_Entries - set only once in each discrete time, in this
//                             routine
//          sol              - array for solution, contains the current
//                             approximation on the solution
//          oldsol           - array with solution from previous discrete time
//          rhs              - array with right hand side f from current
//                             discrete time
//          rhs_old          - rhs f from previous discrete time
//          tilde_u          - low order solution at (Delta t)/2
//          N_neum_to_diri   - number of dofs which should become finally
//                             Dirichlet nodes, MUST BE ORDERED !!!
//          neum_to_diri     - array containing the indices of the dofs which
//                             should become Dirichlet nodes
//          neum_to_diri_bdry- array containing the number of the boundary part
//          neum_to_diri_param - array containing the parameter of the boundary part
//          compute_matrix_D - flag which says if to compute matrix_D_entries
//          BoundaryValue    - pointer to function containing the boundary values
//          BoundaryValues   - contains the boundary values if no pointer is specified
//
// output : matrix_D_Entries - if compute_matrix_D is true
//          B                - array for right hand side for the solver
//          tilde_u          - if compute_matrix_D is true (CB: fuer nichtlinearen Fall nur dann, sonst immer)
//
/*******************************************************************************/

/*! New documentation for the input and output.
 *
 * Apply different variants of the FEM-FCT algebraic flux correction algorithm.
 *
 * @param[out] tilde_u Allocated memory for the low-order auxiliary solution.
 * As far as I see this is only used internally TODO Put construction and destruction into the algorithm,
 * remove as parameter.
 *
 * @param[in] compute_matrix_D - Whether to compute the artificial diffusion matrix D (non-zero)
 * or whether it already is computed and handed over as parameter.
 *
 *
 *
 */
//#ifdef __2D__
void AlgebraicFluxCorrection::FEM_FCT_ForConvDiff(TSquareMatrix2D *M_C, TSquareMatrix2D *A,
			 int N_U, int N_Active, //da diese ohnehin gleich sein muessen und die matrix speziell assembliert werden muss reicht es, hier nur N_U zu uebergeben.
			 double *lump_mass, double *matrix_D_Entries,
			 double *sol, double *oldsol, double *B, double *rhs,
			 double *rhs_old, double *tilde_u,
			 int N_neum_to_diri, int *neum_to_diri, int *neum_to_diri_bdry, double *neum_to_diri_param,
			 int compute_matrix_D, BoundValueFunct2D *BoundaryValue, double *BoundaryValues) //compute_matrix_D kann man sich klemmen, wenn man checkt ob matrix_D_Entries NULL ist
//#endif
//#ifdef __3D__
//void FEM_FCT_ForConvDiff(TSquareMatrix3D *M_C,TSquareMatrix3D *A,
//			 int N_U, int N_Active,
//			 double *lump_mass, double *matrix_D_Entries,
//			 double *sol, double *oldsol,
//			 double *B, double *rhs, double *rhs_old,
//			 double *tilde_u,
//			 int N_neum_to_diri, int *neum_to_diri,
//			 double *neum_to_diri_x, double *neum_to_diri_y, double *neum_to_diri_z,
//			 int compute_matrix_D,
//			 BoundValueFunct3D *BoundaryValue,
//			 double *BoundaryValues)
//#endif
{
  int i,j,j0,j1,index;
  int solver_param[3];
  double eps = 1e-10, help;
  double val, val1;
  double solver_param_d[1];
  double delta_t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;

  // get pointers to columns, rows and entries of matrix M_C
  int* ColInd_M = M_C->GetKCol();
  int* RowPtr_M = M_C->GetRowPtr();
  double* Entries_M = M_C->GetEntries();
  //int N_Entries_M = M_C->GetN_Entries();

  // get pointers to columns, rows and entries of matrix A
  int* ColInd = A->GetKCol();
  int* RowPtr = A->GetRowPtr();
  double* Entries = A->GetEntries();
  int N_Entries = A->GetN_Entries();

  //Allocate space for auxiliary quantities.
  double* res = new double[N_Entries];
  double* alpha = new double[N_Entries];
  double* P_plus = new double[N_U];
  double* P_minus = new double[N_U];
  double* Q_plus = new double[N_U];
  double* Q_minus = new double[N_U];
  double* R_plus = new double[N_U];
  double* R_minus = new double[N_U];
  double* aux_v = new double[N_U];

  //STEP 0: Compute artificial diffusion matrix D and add to A if need be.
  if(compute_matrix_D){
	  computeArtificialDiffusionMatrix(*A, matrix_D_Entries,N_U);
  }
  // Add the artificial diffusion matrix to A giving A+D=L.
  // TODO Should this also only be done once? Differs in the implementation of FEM_TVD and FEM_FCT!!!
  Daxpy(N_Entries, 1.0, matrix_D_Entries, Entries);


  //STEP 1: Computation of the auxiliary solution depending on the chosen
  // TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME.

  // 0 - nonlinear scheme, with minmod prelimiting
  //
  // "linear schemes only for Crank-Nicolson" Worauf genau bezieht sich das?
  //
  // 1 -simplest linear scheme, Kuzmin 2009.
  // 2,3,4 - other linear schemes, but which exactly??

  switch (TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME)
  {
    // simplest linear scheme, from the preprint of Kuzmin (2008)
    // linear schemes only for Crank-Nicolson
    case 1:
      // L u_{k-1}
      MatVect(A,oldsol,aux_v);

      // M_L^(-1)(f_{k-1} - L u_{k-1})
      for(i=0;i<N_U;i++)
      {
    	  aux_v[i] = (rhs_old[i] - aux_v[i])/lump_mass[i];
      }
      //OutPut("auxv0 " << Ddot(N_U,aux_v,aux_v) <<
      // " " << Ddot(N_U,lump_mass,lump_mass) <<
      // " " <<  Ddot(N_U,rhs_old,rhs_old) <<
      // " " <<  Ddot(N_U,oldsol,oldsol) << endl);
      val = delta_t/2.0;
      for(i=0;i<N_U;i++)
      {
        tilde_u[i] = oldsol[i] + val * aux_v[i];
      }

      TDatabase::TimeDB->CURRENTTIME -= val;
      // set correct boundary conditions for intermediate solution
      for (j=0 ;j<N_neum_to_diri;j++)
      {
        i=neum_to_diri[j];
//#ifdef __2D__
	if (BoundaryValue!=NULL)
	{
	    BoundaryValue(neum_to_diri_bdry[j], neum_to_diri_param[j],
			  tilde_u[i]);
	}
	else
	{
	    tilde_u[i] =   BoundaryValues[j];
	}
//#endif
//#ifdef __3D__
//	if (BoundaryValue!=NULL)
//	{
//	    BoundaryValue(neum_to_diri_x[j], neum_to_diri_y[j], neum_to_diri_z[j],
//			  tilde_u[i]);
//	}
//	else
//	{
//	    tilde_u[i] =   BoundaryValues[j];
//	}
//#endif
        aux_v[i] = 2*(tilde_u[i] - oldsol[i])/delta_t;
      }
      TDatabase::TimeDB->CURRENTTIME += val;
      // compute matrix res
      // loop over all rows
      for(i=0;i<N_U;i++)
      {
        // i-th row of sqmatrix
        j0 = RowPtr_M[i];
        j1 = RowPtr_M[i+1];

        for(j=j0;j<j1;j++)
        {
          // column
          index = ColInd_M[j];
          if (index==i)
            continue;
          val = - matrix_D_Entries[j] * (tilde_u[i]-tilde_u[index]);
          res[j]= Entries_M[j] * (aux_v[i]-aux_v[index]) + val;
          // prelimiting with MINMOD
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==1)||
            (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
            res[j] = delta_t*MinMod(res[j]/delta_t, val);
          // more prelimiting
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==2)||
            (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
          {
            if (res[j]*(tilde_u[index]-tilde_u[i])>0)
              res[j] = 0;
          }
          res[j] *= delta_t;
        }
      }



      break; //here the "simplest linear scheme" ends with a casebreak.
    case 2:
    case 3:
    case 4:
      if (TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME > 2)
      {
        theta1 = 0.5;
        theta2 = 1-theta1;
        // compute auxiliary solution, eq. (55) from [Kuz09]
        // rhs comes to aux_v
        // L u_{k-1}, (note: A = -L)
        MatVect(A,oldsol,aux_v);
        // delta_t Lu_{k-1}/2
        val = -theta2 * delta_t;
        Dscal(N_U, val, aux_v);
      }
      val = theta1 * delta_t;
      for(i=0;i<N_U;i++)
      {
        // (M  + delta_t L/2) u_{k-1}
        aux_v[i] += oldsol[i]*lump_mass[i];
        // (M  + delta_t L/2) u_{k-1} + delta_t/2 rhs_old
        aux_v[i] += val * rhs_old[i];
        // (M  + delta_t L/2) u_{k-1} + delta_t/2 rhs_old + delta_t/2 rhs
        aux_v[i] += val * rhs[i];
      }

      // matrix comes to A
      for (i=0;i<N_U;i++)
      {
        // i-th row of sqmatrix
        j0 = RowPtr[i];
        j1 = RowPtr[i+1];
        // j-th column
        for(j=j0;j<j1;j++)
        {
          // column
          index = ColInd[j];
          if (index==i)
          {
            Entries[j] = lump_mass[i] + val * Entries[j];
          }
          else
          {
            Entries[j] *= val;
          }
        }
      }

      // prepare solver
      // use bicgstab with ssor preconditioner
      solver_param[0] = TDatabase::ParamDB->SC_SOLVER_SCALAR;
      TDatabase::ParamDB->SC_SOLVER_SCALAR = 13;
      solver_param[1] = TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR;
      TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = 3;
      solver_param[2] = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
      TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = 100;
      solver_param_d[0] = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
      TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = 1e-15;
      Solver(A, aux_v, tilde_u);

      // reset solver
      TDatabase::ParamDB->SC_SOLVER_SCALAR = solver_param[0];
      TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = solver_param[1];
      TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = solver_param[2];
      TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = solver_param_d[0];

      // reset matrix
      for (i=0;i<N_U;i++)
      {
        // i-th row of sqmatrix
        j0 = RowPtr[i];
        j1 = RowPtr[i+1];
        // j-th column
        for(j=j0;j<j1;j++)
        {
          // column
          index = ColInd[j];
          if (index==i)
          {
            Entries[j] = (Entries[j] - lump_mass[i])/val;
          }
          else
          {
            Entries[j] /= val;
          }
        }
      }

      val = delta_t;
      TDatabase::TimeDB->CURRENTTIME -= val;
      // set correct boundary conditions for auxiliary solution
      for (j=0 ;j<N_neum_to_diri;j++)
      {
        i=neum_to_diri[j];
//#ifdef __2D__
	if (BoundaryValue!=NULL)
	{
	    BoundaryValue(neum_to_diri_bdry[j], neum_to_diri_param[j],
			  tilde_u[i]);
	}
	else
	{
	    tilde_u[i] =   BoundaryValues[j];
	}
//#endif
//#ifdef __3D__
//	if (BoundaryValue!=NULL)
//	{
//	    BoundaryValue(neum_to_diri_x[j], neum_to_diri_y[j], neum_to_diri_z[j],
//			  tilde_u[i]);
//	}
//	else
//	{
//	    tilde_u[i] =   BoundaryValues[j];
//	}
//#endif
      }

      TDatabase::TimeDB->CURRENTTIME += val;

      // intermediate solution should be nonnegative
      for(i=0;i<N_U;i++)
      {
       if (tilde_u[i] < 0)
        OutPut("int sol neg " << tilde_u[i] << endl);
      }

      switch (TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME)
      {
        case 2:
          for(i=0;i<N_U;i++)
          {
            aux_v[i] = tilde_u[i] - oldsol[i];
          }
          break;
        case 3:
          // eq. (60) from [Kuz09]
          // L tilde_u
          MatVect(A,tilde_u,aux_v);
          // approximation of derivative of tilde_u
          for(i=0;i<N_U;i++)
          {
            aux_v[i] = -aux_v[i]/lump_mass[i];
            //if (delta_t * aux_v[i] + tilde_u[i]  < -1e-13)
            //OutPut(delta_t * aux_v[i] + tilde_u[i] << " " << endl);
          }
          break;
        case 4:
          // restore matrix A
          Daxpy(N_Entries, -1.0, matrix_D_Entries, Entries);
          // A tilde_u
          MatVect(A, tilde_u, R_minus);
          // restore matrix D
          Daxpy(N_Entries, 1.0, matrix_D_Entries, Entries);
          // prepare solver
          // use bicgstab with ssor preconditioner
          solver_param[0] = TDatabase::ParamDB->SC_SOLVER_SCALAR;
          TDatabase::ParamDB->SC_SOLVER_SCALAR = 13;
          solver_param[1] = TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR;
          TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = 3;
          solver_param[2] = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
          TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = 100;
          solver_param_d[0] = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
          TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = 1e-11;

          Solver(M_C, R_minus, aux_v);

          // reset solver
          TDatabase::ParamDB->SC_SOLVER_SCALAR = solver_param[0];
          TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = solver_param[1];
          TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = solver_param[2];
          TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = solver_param_d[0];
          memset(R_minus, 0 , N_U*SizeOfDouble);
          break;
      }

      // HAS TO BE CHECKED
      /*for (j=0 ;j<N_neum_to_diri;j++)
      {
          i=neum_to_diri[j];
          aux_v[i] = 0;
	  }*/
      // compute matrix res, eq. (57) in [Kuz09]
      // loop over all rows
      for(i=0;i<N_U;i++)
      {
        // i-th row of sqmatrix
        j0 = RowPtr_M[i];
        j1 = RowPtr_M[i+1];

        for(j=j0;j<j1;j++)
        {
          // column
          index = ColInd_M[j];
          if (index==i)
            continue;
          switch (TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME)
          {
            case 2:
              val = - matrix_D_Entries[j] * (tilde_u[i]-tilde_u[index]);
	      val *= delta_t;
              res[j]= 2 * Entries_M[j] * (aux_v[i]-aux_v[index]) + val;
	      break;
            case 3:
            case 4:
              theta1 = TDatabase::TimeDB->THETA1;
              theta2 = TDatabase::TimeDB->THETA2;
              val = -theta1*tilde_u[i] - theta2*oldsol[i];
              val += theta1*tilde_u[index] + theta2*oldsol[index];
              val *= matrix_D_Entries[j];
              res[j] = Entries_M[j] * (aux_v[i]-aux_v[index]) + val;
              res[j] *= delta_t;
              /*val = -theta1*tilde_u[i] - theta2*oldsol[i];
              val += theta1*tilde_u[index] + theta2*oldsol[index];
              val *= delta_t * matrix_D_Entries[j];
              res[j] = Entries_M[j] *(tilde_u[i] - oldsol[i]);
              res[j] += Entries_M[j] *(tilde_u[index] - oldsol[index]);
              res[j] += val;*/
              break;
          }
          //prelimiting with MINMOD
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==1)||
              (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
          {
             val =  -matrix_D_Entries[j]*(tilde_u[i]-tilde_u[index]);
              res[j] = delta_t*MinMod(res[j]/delta_t, val);
          }
          // more prelimiting
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==2)||
              (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
          {
	    if (res[j]*(tilde_u[index]-tilde_u[i])>0)
	      res[j] = 0;
	  }
        }
      }
      break;
    case 0:
      //=======================================================
      //nonlinear scheme
      //=======================================================
      if (compute_matrix_D)
      {
	  // reset parameter if necessary
	  if (TDatabase::ParamDB->FEM_FCT_PRELIMITING==0)
	  {
	      TDatabase::ParamDB->FEM_FCT_PRELIMITING = 2;
	      OutPut("nonlinear FEM-FCT works badly without prelimiting, "
		     << " set ParamDB->FEM_FCT_PRELIMITING to "
		     << TDatabase::ParamDB->FEM_FCT_PRELIMITING << endl);
	  }
        // L u_{k-1}
        MatVect(A,oldsol,aux_v);

        // M_L^(-1)(f_{k-1} - L u_{k-1})
        val = delta_t/2.0;
        for(i=0;i<N_U;i++)
        {
          val1 = (rhs_old[i] - aux_v[i])/lump_mass[i];
          tilde_u[i] = oldsol[i] + val * val1;
        }
        //auxiliary variable is computed
        TDatabase::TimeDB->CURRENTTIME -= val;
        // set correct boundary conditions for intermediate solution
        for (j=0;j<N_neum_to_diri;j++)
        {
          i=neum_to_diri[j];
//#ifdef __2D__
	if (BoundaryValue!=NULL)
	{
          BoundaryValue(neum_to_diri_bdry[j], neum_to_diri_param[j],
            tilde_u[i]);
	}
	else
	{
	    tilde_u[i] =   BoundaryValues[j];
	}

          //aux_v[i] = 2*(tilde_u[i] - oldsol[i])/delta_t;
//#endif
//#ifdef __3D__
//          // WRONG
//	if (BoundaryValue!=NULL)
//	{
//          BoundaryValue(neum_to_diri_x[j], neum_to_diri_y[j], neum_to_diri_z[j],
//            tilde_u[i]);
//	}
//	else
//	{
//	    tilde_u[i] =   BoundaryValues[j];
//	}
//#endif
        }
        TDatabase::TimeDB->CURRENTTIME += val;
      }
      // compute matrix res
      // loop over all rows
      val = theta1 * delta_t;
      val1 = theta2 * delta_t;
      for(i=0;i<N_U;i++)
      {
        // i-th row of sqmatrix
        j0 = RowPtr_M[i];
        j1 = RowPtr_M[i+1];

        for(j=j0;j<j1;j++)
        {
          // column
          index = ColInd_M[j];

          res[j]= -Entries_M[j] * ((sol[index]-sol[i]) - (oldsol[index] - oldsol[i]))
            + val * matrix_D_Entries[j] * (sol[index] - sol[i])
            + val1 * matrix_D_Entries[j] * (oldsol[index]- oldsol[i]);
          // prelimiting with MINMOD, THIS HAS TO BE CHECKED
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==1)||
            (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
            res[j] = delta_t*MinMod(res[j]/delta_t, -matrix_D_Entries[j] * (tilde_u[i]-tilde_u[index]));
          // more prelimiting
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==2)||
            (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
          {
            if (res[j]*(tilde_u[index]-tilde_u[i])>0)
              res[j] = 0;
          }
        }
      }
      break;
    default:
      OutPut("FEM_FCT_LINEAR_TYPE " <<
        TDatabase::ParamDB->FEM_FCT_LINEAR_TYPE << " DOES NOT EXIST !!!"<<endl);
      exit(4711);
  }

  //=======================================================
  //Zalesaks limiter
  //=======================================================

  for(i=0;i<N_U;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr_M[i];
    j1 = RowPtr_M[i+1];
    for(j=j0;j<j1;j++)
    {
      index = ColInd_M[j];
      // only the active part of the matrix
      if(i != index)
      {
        // if (i==1111)
        //  OutPut(i << " " << index << " " << res[j] << " : ");
        if (res[j] > 0.0)
          P_plus[i] += res[j];
        if (res[j] <= 0.0)
          P_minus[i] += res[j];

        help = tilde_u[index]-tilde_u[i];

        // if (i==1111)
        //    OutPut(tilde_u[index] << " " << tilde_u[i] << " " << help << endl);
        if (help > Q_plus[i])
          Q_plus[i]= help;

        if (help < Q_minus[i])
          Q_minus[i]= help;
      }
    }                                             // end loop j
  }
  //for(i=0;i<N_U;i++)
  //    OutPut(i << " "<< P_plus[i] << " " <<  P_minus[i] << " " << Q_plus[i] << " " <<  Q_minus[i] << endl);
  //exit(1);
  /*
    for(i=0;i<N_U;i++)
    {
        if (tilde_u[i] + Q_minus[i] < t_min)
      OutPut("tilde1 " << i << " "  << tilde_u[i] + Q_minus[i] <<
       " " << tilde_u[i]  << " "  <<   Q_minus[i] <<
       " " << t_min << endl);
    }
  */

  for(i=0;i<N_U;i++)
  {
    if(fabs(P_plus[i]) == 0.0)
      R_plus[i] = 1.0;
    else
    {
      //OutPut(Q_plus[i] << " "  << P_plus[i] << " " << lump_mass[i] << " ");
      help = (lump_mass[i] * Q_plus[i])/P_plus[i];
      if(help < 1.0)
        R_plus[i] = help;
      else
        R_plus[i] = 1.0;
    }
    if(fabs(P_minus[i]) == 0.0)
      R_minus[i] = 1.0;
    else
    {
      //OutPut(Q_minus[i] << " "  << P_minus[i] << " " << lump_mass[i] << " ");
      help = (lump_mass[i] * Q_minus[i])/P_minus[i];
      if(help < 1.0)
        R_minus[i] = help;
      else
        R_minus[i] = 1.0;
    }
  }

  // treat Dirichlet nodes
  for (j=0;j<N_neum_to_diri;j++)
  {
    i = neum_to_diri[j];
    R_plus[i] = 1;
    R_minus[i] = 1;
  }

  for(i=0;i<N_U;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];

    for(j=j0;j<j1;j++)
    {
      index = ColInd[j];
      if(res[j] > 0.0)
      {
        //Initialisation
        alpha[j] = R_plus[i];
        if(alpha[j] > R_minus[index])
          alpha[j] = R_minus[index];
      }
      if(res[j]<=0.0)
      {
        //initialisation
        alpha[j] = R_minus[i];
        if(alpha[j] > R_plus[index])
          alpha[j] = R_plus[index];
      }
      // clipping, see Kuzmin (2009), end of Section 5
      //if ((fabs(Q_plus[i])< eps)&&(fabs(Q_minus[i])< eps))
      //  alpha[j] = 0;
    }                                             //end loop j
  }

  //=======================================================
  //correct right hand side
  //=======================================================

  // Crank--Nicolson
  if ((fabs(theta2-0.5)<eps) && (fabs(theta3-0.5)<eps)&&
    (TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME<2)) //linear scheme 1, nonlinear scheme 0, only for Crank-Nicolson
  {
    for(i=0;i<N_U;i++)
    {
      val=0.0;
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      for(j=j0;j<j1;j++)
      {
        index = ColInd[j];
        if(i != index)
        {
          val +=  alpha[j] * res[j];
        }
      }
      B[i] = tilde_u[i] * lump_mass[i] + theta4*delta_t*rhs[i] + val;
    }
  }
  else //linear schemes 2,3,4; other time steps?
  {
    // L u_{k-1}
    MatVect(A,oldsol,aux_v);
    val1 = theta2*delta_t;
    for(i=0;i<N_U;i++)
    {
      val =0.0;
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      for(j=j0;j<j1;j++)
      {
        index = ColInd[j];
        if(i != index)
        {
          val +=  alpha[j] * res[j];
        }
      }
      B[i] = lump_mass[i] * oldsol[i] - val1 * aux_v[i] + delta_t * theta3 * rhs_old[i]
        + delta_t * theta4*rhs[i] + val;
      /* if (B[i] < -1e-10)
      {
          OutPut(i << " " << B[i] <<  " " << lump_mass[i] * oldsol[i]  << " " <<
           -val1 * aux_v[i] << " " <<
           delta_t * theta3 * rhs_old[i] << " " <<
           delta_t * theta4*rhs[i] << " " <<  val << endl);
           }*/
    }
  }

  //OutPut(" vecB " << Ddot(N_U,B,B) << endl);
  //Clean up.
  delete[] res;
  delete[] alpha;
  delete[] P_plus;
  delete[] P_minus;
  delete[] Q_plus;
  delete[] Q_minus;
  delete[] R_plus;
  delete[] R_minus;
  delete[] aux_v;
}

/*!
 * Compute the artificial diffusion matrix. Put here to make FEM-FCT algo shorter.
 * @param[in] A The stiffness matrix which need additional diffusion.
 * @param[out] matrix_D_Entries Memory space to write the entries of D, the artificial diffusion matrix.
 * @param[in] N_U The number of degrees of freedom, i.e. number of rows of A.
 */
void AlgebraicFluxCorrection::computeArtificialDiffusionMatrix(const TSquareMatrix2D& A, double* matrix_D_Entries, int N_U){

	  // get pointers to columns, rows and entries of matrix A
	  int* ColInd = A.GetKCol();
	  int* RowPtr = A.GetRowPtr();
	  double* Entries = A.GetEntries();
	  int N_Entries = A.GetN_Entries();

		memset(matrix_D_Entries , 0, N_Entries*SizeOfDouble);
		// loop over all rows
		for(int i=0;i<N_U;i++)
		{
			// i-th row of sqmatrix
			int j0 = RowPtr[i];
			int j1 = RowPtr[i+1];
			// compute first the matrix D
			for(int j=j0;j<j1;j++)
			{
				// column
				int index = ColInd[j];
				// only the active part of the matrix
				//if (index>=N_U)
				//	  continue;
				// only off-diagonals
				if (index!=i)
				{
					if (Entries[j] > 0)
						matrix_D_Entries[j] = -Entries[j];
					// now check the transposed entry
					int j2 = RowPtr[index];
					int j3 = RowPtr[index+1];
					for (int jj=j2;jj<j3;jj++)
					{
						if (ColInd[jj]==i)
						{
							if (-Entries[jj]<matrix_D_Entries[j])
								matrix_D_Entries[j] = -Entries[jj];
							break;
						}
					}
				}
			}
		}

		// compute diagonal entry of D
		// loop over all rows
		for(int i=0;i<N_U;i++)
		{
			// i-th row of sqmatrix
			int j0 = RowPtr[i];
			int j1 = RowPtr[i+1];
			double val = 0.0;
			// add all entries of i-th row
			int jj = -1;
			for(int j=j0;j<j1;j++)
			{
				val +=  matrix_D_Entries[j];
				int index = ColInd[j];
				if (index==i)
					jj = j; //Hold the place of the diagonal entry in the entries array. Is this reached exactly once ?
			}
			matrix_D_Entries[jj] = -val;
		}
}

/*! New documentation for the input and output.
 *
 * Apply different variants of the FEM-FCT algebraic flux correction algorithm.
 *
 * @param[out] tilde_u Allocated memory for the low-order auxiliary solution.
 * As far as I see this is only used internally TODO Put construction and destruction into the algorithm,
 * remove as parameter.
 *
 * @param[in] compute_matrix_D - Whether to compute the artificial diffusion matrix D (non-zero)
 * or whether it already is computed and handed over as parameter.
 *
 *
 *
 */
void AlgebraicFluxCorrection::FEM_FCT_SimpleLinear(TSquareMatrix2D *M_C,TSquareMatrix2D *A,
			 int N_U, int N_Active,
			 double *lump_mass, double *matrix_D_Entries, int compute_matrix_D,
			 double *sol, double *oldsol,
			 double *B, double *rhs, double *rhs_old,
			 int N_neum_to_diri, int *neum_to_diri,
			 double *BoundaryValues){

	double eps = 1e-10;
	double delta_t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
	double theta1 = TDatabase::TimeDB->THETA1;
	double theta2 = TDatabase::TimeDB->THETA2;
	double theta3 = TDatabase::TimeDB->THETA3;
	double theta4 = TDatabase::TimeDB->THETA4;

	// get pointers to columns, rows and entries of matrix M_C
	int* ColInd_M = M_C->GetKCol();
	int* RowPtr_M = M_C->GetRowPtr();
	double* Entries_M = M_C->GetEntries();

	// get pointers to columns, rows and entries of matrix A (later L); plus number of entries.
	int* ColInd = A->GetKCol();
	int* RowPtr = A->GetRowPtr();
	double* Entries = A->GetEntries();
	int N_Entries = A->GetN_Entries();

	//Allocate space for auxiliary quantities and initialize with 0.
	double* res = new double[N_Entries]();
	double* alpha = new double[N_Entries]();
	double* P_plus = new double[N_U]();
	double* P_minus = new double[N_U]();
	double* Q_plus = new double[N_U]();
	double* Q_minus = new double[N_U]();
	double* R_plus = new double[N_U]();
	double* R_minus = new double[N_U]();
	double* aux_v = new double[N_U]();

	//STEP 0: Compute artificial diffusion matrix D and add to A if need be.
	if(compute_matrix_D){
		computeArtificialDiffusionMatrix(*A, matrix_D_Entries,N_U);
		// Add the artificial diffusion matrix to A giving A+D=L.
		// This is only done once, too - from here on the matrix A is equipped with lots of additional diffusion.
		Daxpy(N_Entries, 1.0, matrix_D_Entries, Entries);
	}

	//STEP 1: Computation of the auxiliary solution. Here I stripped the code down
	// to the "simplest linear scheme", case 1 in the original algorithm.

	double* tilde_u = new double[N_U]; //memory for the auxiliary solution

	// L u_{k-1}
	MatVect(A,oldsol,aux_v); //multiply A*oldsol and store in aux_v

	// M_L^(-1)(f_{k-1} - L u_{k-1})
	for(int i=0;i<N_U;i++)
	{
		aux_v[i] = (rhs_old[i] - aux_v[i])/lump_mass[i];
	}

	double halftimestep = delta_t/2.0;
	for(int i=0;i<N_U;i++)
	{
		tilde_u[i] = oldsol[i] + halftimestep * aux_v[i];
	}

	// set correct boundary conditions for intermediate solution
	for (int j=0 ;j<N_neum_to_diri;j++)
	{
		int i=neum_to_diri[j];
		//sicherstellen, dass die richtigen bdry values uebergeben werden, siehe Handmitschrift 27.07.2015.
		tilde_u[i] =   BoundaryValues[j];
		//aux_v wird auf einen Zwischenstand zurueck gesetzt
		aux_v[i] = 2*(tilde_u[i] - oldsol[i])/delta_t;
	}
	//Computation of auxiliary solution is complete.

	//STEP 2: compute matrix res (?) - raw antidiffusive fluxes
	// loop over all rows of mass matrix
	for(int i=0;i<N_U;i++)
	{
		// i-th row of mass matrix
		int j0 = RowPtr_M[i];
		int j1 = RowPtr_M[i+1];

		for(int j=j0;j<j1;j++)
		{
			// column
			int index = ColInd_M[j];

			if (index==i) //nothing to do if this is a diagonal entry (?)
				continue;

			double val = - matrix_D_Entries[j] * (tilde_u[i]-tilde_u[index]);
			res[j]= Entries_M[j] * (aux_v[i]-aux_v[index]) + val;

			// prelimiting with MINMOD
			if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==1)||
					(TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
				res[j] = delta_t*MinMod(res[j]/delta_t, val);
			// prelimiting by setting diffusive fluxes in gradient direction to 0
			if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==2)||
					(TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
			{
				if (res[j]*(tilde_u[index]-tilde_u[i])>0)
					res[j] = 0;
			}
			// It seems there are 3 prelimiting strategies implemented. 1 is only minmod, 2
			// is setting diffusive fluxes in gradient direction to 0 and 3 is combining both.
			// If FEM_FCT_PRELIMITING is neither 1,2 nor 3 no prelimiting as applied.

			res[j] *= delta_t;
		}
	}

	//STEP 3: Apply Zalesak's limiter, input quantities are "res" and "tilde_u".

	// 3.1 	Compute Ps and Qs (negative/positive diffusive/antidiffusive fluxes),
	//		where the Qs are "distances to local extrema of the auxiliary solution"
	for(int i=0;i<N_U;i++)
	{
		// i-th row of mass matrix
		int j0 = RowPtr_M[i];
		int j1 = RowPtr_M[i+1];
		for(int j=j0;j<j1;j++)
		{
			int index = ColInd_M[j];

			if(i != index) //diagonal entries are skipped
			{
				if (res[j] > 0.0)
					P_plus[i] += res[j]; //P_plus was initialized with 0
				if (res[j] <= 0.0)
					P_minus[i] += res[j]; //P_minus was initialized with 0

				double help = tilde_u[index]-tilde_u[i];

				if (help > Q_plus[i]) //Q_plus was initialized with 0
					Q_plus[i]= help;

				if (help < Q_minus[i]) //Q_minus was initialized with 0
					Q_minus[i]= help;
			}
		}                                             // end loop j
	}

	// 3.2	Compute Rs (nodal correction factors)
	for(int i=0;i<N_U;i++)
	{
		if(fabs(P_plus[i]) == 0.0)
			R_plus[i] = 1.0;
		else
		{
			//OutPut(Q_plus[i] << " "  << P_plus[i] << " " << lump_mass[i] << " ");
			double help = (lump_mass[i] * Q_plus[i])/P_plus[i];
			if(help < 1.0)
				R_plus[i] = help;
			else
				R_plus[i] = 1.0;
		}
		if(fabs(P_minus[i]) == 0.0)
			R_minus[i] = 1.0;
		else
		{
			//OutPut(Q_minus[i] << " "  << P_minus[i] << " " << lump_mass[i] << " ");
			double help = (lump_mass[i] * Q_minus[i])/P_minus[i];
			if(help < 1.0)
				R_minus[i] = help;
			else
				R_minus[i] = 1.0;
		}
	}

	// treat Dirichlet nodes (e.g. Kuzmin & Moeller 2005 S. 27) TODO "the same should be performed at outflow boundaries"
	for (int j=0;j<N_neum_to_diri;j++)
	{
		int i = neum_to_diri[j];
		R_plus[i] = 1;
		R_minus[i] = 1;
	}

	// 3.3 Compute alphas (final correction factors)
	for(int i=0;i<N_U;i++)
	{
		// i-th row of sqmatrix
		int j0 = RowPtr[i];
		int j1 = RowPtr[i+1];

		for(int j=j0;j<j1;j++)
		{
			int index = ColInd[j];
			if(res[j] > 0.0)
			{
				//Initialisation
				alpha[j] = R_plus[i];
				if(alpha[j] > R_minus[index])
					alpha[j] = R_minus[index];
			}
			if(res[j]<=0.0)
			{
				//initialisation
				alpha[j] = R_minus[i];
				if(alpha[j] > R_plus[index])
					alpha[j] = R_plus[index];
			}
			// clipping, see Kuzmin (2009), end of Section 5
			//if ((fabs(Q_plus[i])< eps)&&(fabs(Q_minus[i])< eps))
			//  alpha[j] = 0;
		}                                             //end loop j
	}

	// STEP 4: Calculate the right hand side.
	// This is, apart from the addition of a lot of diffusion to A, the only place where smething is actually changed.
	// but strange enough, the new rhs with the additional limited antidiffusive fluxes is stored in "B",
	// not in rhs. How is one supposed to use this class from the outside???

	//Die rechte Seite-Berechnung fuer simple linear (1) und nonlinear (0) ist nur fuer C-N implementiert (???)
	// Crank--Nicolson bzw. theta2 und theta3 sollen 0.5 sein - fuer C-N muss das fuer alle 4 thetas gelten....
	// Ich wuerde diese ganze Struktur gerne aendern. Warum soll ueberhaupt "B" vollgeschrieben werden,
	// statt die ganze rechte Seite zu aendern? Und warum wird die Addition von theta4*delta_t*rhs[i]
	// hier vorgenommen? Dafuer sorgt doch im besten Fall (intermediate Version) schon die Assemblierungsmethode...
	// TODO Change the whole set up.
	if ((fabs(theta2-0.5)<eps) && (fabs(theta3-0.5)<eps)) //Why ask for theta 2 and theta3 ??
	{
		for(int i=0;i<N_U;i++)
		{
			double val=0.0;
			int j0 = RowPtr[i];
			int j1 = RowPtr[i+1];
			for(int j=j0;j<j1;j++)
			{
				int index = ColInd[j];
				if(i != index)
				{
					val +=  alpha[j] * res[j];
				}
			}
			B[i] = tilde_u[i] * lump_mass[i] + theta4*delta_t*rhs[i] + val;
		}
	}

	//Clean up.
	delete[] res; delete[] alpha; delete[] P_plus; delete[] P_minus;
	delete[] Q_plus; delete[] Q_minus; delete[] R_plus; delete[] R_minus;
	delete[] aux_v;
	delete[] tilde_u;
}

void AlgebraicFluxCorrection::correctDirichletRows(TSquareMatrix2D& MatrixA)
{
	//hold pointers to row, kcol, entries array
	int* RowPtr_A      = MatrixA.GetRowPtr();
	int* KCol_A        = MatrixA.GetKCol();
	double* Entries_A  = MatrixA.GetEntries();

	//determine first and one-after-last dirichlet rows
	int diriHighBound = MatrixA.GetFESpace()->GetDirichletBound();
	int diriLowBound = diriHighBound - MatrixA.GetFESpace()->GetN_Dirichlet();

	// loop over rows and set them to unity-vectors
	for (size_t rowIndex = diriLowBound;
			rowIndex < diriHighBound ;++rowIndex)
	{
		int l0 = RowPtr_A[rowIndex];
		int l1 = RowPtr_A[rowIndex+1];
		for (int l=l0;l<l1;l++)
		{
			// diagonal entry
			if (KCol_A[l]==rowIndex)
				Entries_A[l] = 1;
			else
				Entries_A[l] = 0;
		}
	}
}
