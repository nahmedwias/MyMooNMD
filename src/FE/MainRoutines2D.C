// =======================================================================
// @(#)MainRoutines.C
//
// Purpose: contains routines which are called from the main program 
//
// Author: Volker John 
//
// History: start of implementation 22.09.2009
//
// =======================================================================

#include <Constants.h>
#include <Database.h>
#include <AuxParam2D.h>
#include <DirectSolver.h>
#include <Solver.h>
#include <ItMethod.h>
#include <FEDatabase2D.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <JacobiIte.h>
#include <MultiGridScaIte.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <Output2D.h>
#include <ConvDiff2D_Routines.h>
#include <Assemble2D.h>
#include <MainUtilities.h>
#include <MainRoutines2D.h>
#include <Upwind.h>

#include <stdlib.h>
// #include <malloc.h>
#include <string.h>
#include <sstream>

/******************************************************************************/
// SetParametersCDAdapt2D()
// sets parameters of the data base for the main program CDAdapt2D.C
/******************************************************************************/

void SetParametersCDAdapt2D()
{
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== KLR02_3)
    TDatabase::ParamDB->SOLD_S = 0;
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== LP96)
  {
    OutPut("SOLD_PARAMETER_TYPE == LP96 should be used with higher quadrature rule,"<<endl);
    OutPut("since right hand side is in general not linear !!!"<<endl);
  }

  if(TDatabase::ParamDB->DISCTYPE != LOCAL_PROJECTION)
  {
    // switch off all local projection terms
    TDatabase::ParamDB->LP_FULL_GRADIENT = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;

    TDatabase::ParamDB->LP_STREAMLINE = 0;
    TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
    TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
  }

  if(TDatabase::ParamDB->LP_FULL_GRADIENT)
  {
    if(TDatabase::ParamDB->LP_STREAMLINE)
    {
      TDatabase::ParamDB->LP_STREAMLINE = 0;
      TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
      TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
      OutPut("local projection stabilisation in streamline direction ");
      OutPut("is switched off due to stabilisation of full gradient." << endl);
    }
  }

  if(TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  if(TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;
  

  if ((TDatabase::ParamDB->SDFEM_TYPE == 100)&&(!TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&
      (TDatabase::ParamDB->DISCTYPE==SDFEM))
  {
      TDatabase::ParamDB->SDFEM_TYPE = 2;
      OutPut("Changed Database::ParamDB->SDFEM_TYPE to " << TDatabase::ParamDB->SDFEM_TYPE
	     << " since no adjoint problem is solved !!! "<<endl);
  } 
  if ((TDatabase::ParamDB->SDFEM_TYPE != 100)&&(TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&
      (TDatabase::ParamDB->DISCTYPE==SDFEM))
  {
      TDatabase::ParamDB->SDFEM_TYPE = 100;
      OutPut("Changed Database::ParamDB->SDFEM_TYPE to " << TDatabase::ParamDB->SDFEM_TYPE
	     << " since adjoint problem is solved !!! "<<endl);
  } 
  if ((TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM==4))
    TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS = 1;
  // SUPG 
  if ((TDatabase::ParamDB->DISCTYPE==SDFEM)&&(TDatabase::ParamDB->SOLD_TYPE==0))
  {
      // this excludes some not wished side effects
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
  } 
  if  (!(TDatabase::ParamDB->DISCTYPE==SDFEM))
    {
      TDatabase::ParamDB->SOLD_TYPE = 0;
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE =0;
    }
  if  ((TDatabase::ParamDB->DISCTYPE==SDFEM)&&(TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD))
    {
      TDatabase::ParamDB->SDFEM_TYPE = 0;
      TDatabase::ParamDB->DELTA0 =  TDatabase::ParamDB->DELTA1 = 0;
      OutPut("FEM-TVD: switched stabilization off!" << endl);
    }

  TDatabase::ParamDB->NSTYPE = 0;

  if (TDatabase::ParamDB->DISCTYPE==CIP)
    {
      TDatabase::ParamDB->DISCTYPE=GALERKIN;
      TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS = 1;
    }
  if (TDatabase::ParamDB->DISCTYPE==DG)
    {
      TDatabase::ParamDB->DISCTYPE=GALERKIN;
      TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS = 2;
      if ( TDatabase::ParamDB->ANSATZ_ORDER < 10)
	TDatabase::ParamDB->ANSATZ_ORDER = -TDatabase::ParamDB->ANSATZ_ORDER-10;
      else 
	// P elements on quads
	TDatabase::ParamDB->ANSATZ_ORDER = -10*TDatabase::ParamDB->ANSATZ_ORDER;
      if (TDatabase::ParamDB->ESTIMATE_ERRORS)
	{
	  TDatabase::ParamDB->ESTIMATE_ERRORS = 0;
	  OutPut("Error estimation does not work for DG !!!"<< endl);
	}
    }
}

/******************************************************************************/
// Solver
// solves linear system
// output in sol
// scalar problems: ns_type == 0
/******************************************************************************/
void Solver(TSquareMatrix **sqmatrices, TMatrix **matrices,
	    double *rhs, double *sol, 
	    MatVecProc *MatVect, DefectProc *Defect,
	    TMultiGrid2D *MG, 
	    int N_Unknowns, int ns_type)
{
  int solver_type, prec_type;
  double *itmethod_sol, *itmethod_rhs, t1, t2;
  TItMethod *itmethod, *prec;
  TSquareMatrix2D *A;

  t1 = GetTime();
  
  if (!ns_type)
  {
      prec_type = TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR;
      solver_type = TDatabase::ParamDB->SC_SOLVER_SCALAR;
  }
  else
  {
      prec_type = TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE;
      solver_type = TDatabase::ParamDB->SC_SOLVER_SADDLE;
  }

  switch(ns_type)
  {
      case 0:
	  A = (TSquareMatrix2D*) sqmatrices[0];
	  break;
  }
  /* if (TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE == 0)
    {
      TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE =  TDatabase::ParamDB->SOLVER_TYPE;
      TDatabase::ParamDB->SOLVER_TYPE = 2;
      }*/
  switch(TDatabase::ParamDB->SOLVER_TYPE)
  {
      case DIRECT:
	  switch(ns_type)
	  {
	      case 0:
		  DirectSolver(A, rhs, sol);
		  break;
	  }
	  //TDatabase::ParamDB->SOLVER_TYPE = 1;
	  break;
	  
      case AMG_SOLVE:
	  switch(ns_type)
	  {
	      case 0:
		  Solver(A, rhs, sol);
		  break;
	  }
	  break;
	  
      case GMG:
	  switch (prec_type)
	  {
	      case 1:
		  prec = new TJacobiIte(MatVect, Defect, NULL,
					0, N_Unknowns, 1);
		  break;
	      case 5:
		  prec = new TMultiGridScaIte(MatVect, Defect, NULL,
					      0, N_Unknowns, MG, 0);
		  break;
	      default:
		  OutPut("Unknown preconditioner !!!" << endl);
		  exit(4711);
	  }
	  switch (solver_type)
	  {
	      case 11:
		  itmethod = new TFixedPointIte(MatVect, Defect, prec,
						0, N_Unknowns, 1);
		  if (prec_type == 5)
		  {
		      itmethod_sol = new double[N_Unknowns];
		      itmethod_rhs = new double[N_Unknowns];
		      memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
		      memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
		  }
		  else
		  {
		      itmethod_sol = sol;
		      itmethod_rhs = rhs;
		  }
		  break;
	      case 16:
		  itmethod = new TFgmresIte(MatVect, Defect, prec,
					    0, N_Unknowns, 1);
		  if (prec_type == 5)
		  {
		      itmethod_sol = new double[N_Unknowns];
		      itmethod_rhs = new double[N_Unknowns];
		      memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
		      memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
		  }
		  else
		  {
		      itmethod_sol = sol;
		      itmethod_rhs = rhs;
		  }
		  break;
	      default:
		  OutPut("Unknown solver !!!" << endl);
		  exit(4711);
	  }
	  // solve linear system
	  itmethod->Iterate(sqmatrices,NULL,itmethod_sol,itmethod_rhs);
	  
	  delete prec;
	  delete itmethod;
	  
	  switch (solver_type)
	  {
	      case 11:
		  if (prec_type == 5)
		  {
		      memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
		      memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
		      delete itmethod_sol;
		      delete itmethod_rhs;
		  }
		  break;
	      case 16:
		  if (prec_type == 5)
		  {
		      memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
		      memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
		      delete itmethod_sol;
		      delete itmethod_rhs;
		  }
		  break;
	  }
	  break;
  }
  t2 = GetTime();
  if (TDatabase::ParamDB->SC_VERBOSE>1)
    OutPut("time for solving: " << t2-t1 << endl);
}

/******************************************************************************/
// ComputeErrorEstimate
// computes residual based error estimates
// output in eta_K and estimated_global_error
/******************************************************************************/

void ComputeErrorEstimate(TCollection *coll, TFEFunction2D *u,
			  CoeffFct2D *Coeffs, BoundCondFunct2D **BoundaryConditions,
			  BoundValueFunct2D **BoundaryValues,
			  double* &eta_K,
			  double *maximal_local_error,
			  double *estimated_global_error,
			  double l2, double h1, int N_Unknowns)
{
    int N_Cells, bdry_tvd = 0, bdry_ad = 0;
    double t1, t2;
    TAuxParam2D *aux;
    TCD2DErrorEstimator *CDErrorEstimator;
    TFESpace2D *fesp[2];
    MultiIndex2D Derivatives_All[5] = { D10, D01, D00, D20, D02 };
 
    // set correct boundary conditions
    if (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
    {
	bdry_ad = TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM;
	TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM = 0;
    }
    if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD);
    {
	bdry_tvd = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
	TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
    }

    N_Cells = coll->GetN_Cells();
    // allocate arrays for local estimate
    eta_K = new double[N_Cells];

    CDErrorEstimator = new TCD2DErrorEstimator(TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION,
					     u,
					     TDatabase::ParamDB->ERROR_CONTROL);
    t1 = GetTime();
    aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    fesp[0] = u->GetFESpace2D();
    CDErrorEstimator->GetErrorEstimate(5, Derivatives_All,
				       Coeffs, BoundaryConditions,
				       BoundaryValues, aux, 1, fesp,
				       eta_K, &maximal_local_error[0], &estimated_global_error[0]);
    delete aux;
    delete CDErrorEstimator;
    t2 = GetTime();

    if (TDatabase::ParamDB->SC_VERBOSE>1)
      {
	OutPut("time for error estimation: " << t2-t1 << endl);
	
	OutPut("gradient estimate " << estimated_global_error[0] << endl);
	OutPut("L2 estimate " << estimated_global_error[2]);
	if ((l2>0)&&(TDatabase::ParamDB->MEASURE_ERRORS))
	  {
	    OutPut("   effectivity index " << estimated_global_error[2]/l2 << endl);
	  }
	else
	  {
	    OutPut(endl);
	  }
	OutPut("H1 estimate " << estimated_global_error[1]);
	if ((h1>0)&&(TDatabase::ParamDB->MEASURE_ERRORS))
	  {
	    OutPut("   effectivity index " << estimated_global_error[1]/h1 << endl);
	  }
	else
	  {
	    OutPut(endl);
	  }
	OutPut("energy estimate " << estimated_global_error[3]);
	if ((l2*l2+h1*h1>0)   &&(TDatabase::ParamDB->MEASURE_ERRORS))
	  {
	    OutPut("   effectivity index " << estimated_global_error[3]/
		   sqrt(l2*l2+h1*h1/TDatabase::ParamDB->RE_NR));// << endl);
	  }
	else
	  {
	      ;//OutPut(endl);
	  }
	OutPut(" (without edge res.) " << estimated_global_error[4]);
	OutPut(endl);
	
	if (TDatabase::ParamDB->MEASURE_ERRORS)
	  OutPut("eff_ind " << N_Unknowns << " " << estimated_global_error[2]/l2
		 << " " << estimated_global_error[1]/h1 <<
		 " " << estimated_global_error[3]/
		 sqrt(l2*l2+h1*h1/TDatabase::ParamDB->RE_NR) << endl);
      }
    // reset parameters
      if (bdry_ad)
	TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM = bdry_ad;

      if (bdry_tvd)
	  TDatabase::ParamDB->SOLD_PARAMETER_TYPE = bdry_tvd;
}

/******************************************************************************/
// ComputeAdjointMatrix
// computes the adjoint matrix of a matrix A
// it is assumed that both matrices have the same structure
/******************************************************************************/

void ComputeAdjointMatrix(TSquareMatrix2D *A,TSquareMatrix2D *AT)
{
    int i, j, k, n, m, ii, k1, k2, index, *row, *col;
    double *a, *at, val;

    m = A->GetN_Rows();
    n = A->GetN_Columns();
    row = A->GetRowPtr();
    col = A->GetKCol();
    a = A->GetEntries();
    at = AT->GetEntries();
  
    // loop over A
    j = row[0];
    for (i=0;i<m;i++)
    {
        ii = row[i+1];
	for (;j<ii;j++)
	{
	    index = col[j];
	    val = a[j];
	    // row in AT
	    k1 = row[index];
	    k2 = row[index+1];
	    for (k=k1;k<k2;k++)
	    {
		// compare column of AT with row of A
		if (col[k] == i)
		{
		    at[k] = val;
		    break;
		}
	    }
	}
    }
}

/******************************************************************************/
// OutputData2D
// writes the outputs
/******************************************************************************/
void OutputData2D(std::ostringstream& os, TOutput2D *Output,
		  int counter)
{
    if(TDatabase::ParamDB->WRITE_GRAPE)
    {
      os.seekp(std::ios::beg);
      os << TDatabase::ParamDB->GRAPEBASENAME << counter << ".dat" << ends;
      Output->WriteGrape(os.str().c_str());
    }
    if(TDatabase::ParamDB->WRITE_GMV)
    {
      os.seekp(std::ios::beg);
      os << TDatabase::ParamDB->GMVBASENAME  << counter << ".gmv" << ends;
      Output->WriteGMV(os.str().c_str());
    }

    if(TDatabase::ParamDB->WRITE_GNU)
    {
      os.seekp(std::ios::beg);
      os <<  TDatabase::ParamDB->GNUBASENAME  << counter << ".gnu" << ends;
      Output->WriteGnuplot(os.str().c_str());
    }

    if(TDatabase::ParamDB->WRITE_VTK)
    {
      os.seekp(std::ios::beg);
      os <<  TDatabase::ParamDB->VTKBASENAME << counter << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
    }
}

/******************************************************************************/
// PrepareAdjointProblem2D
// preparation of assembling and solution of adjoint problem
// the array for the piecewise constant stabilization parameters is filled
//    with standard values
/******************************************************************************/
void PrepareAdjointProblem2D(TCollection *coll,
			     TFEFunction2D* &u_adjoint,
			     TSquareMatrix2D* &sqmatrixAadjoint,
			     TFESpace2D* &pw_const_param_space,
			     TFEFunction2D* &pw_const_param_fe,
			     TFEFunction2D* &pw_const_param_deriv_fe,
			     TFESpace2D *velocity_space,
			     TSquareStructure2D *sqstructureA,
			     BoundCondFunct2D *BoundCondition,
			     CoeffFct2D *Coeffs,
			     double* &sol_adjoint,
			     double* &rhs_edge,
			     double* &pw_const_param,
			     double* &pw_const_param_deriv,
			     double* &pw_const_param_old,
			     int N_U, int N_Cells)
{
  int i, j, N_Edges, type;
  double *coeff, x, y, xs, ys;
  char UadjointString[] = "u_adjoint";
  char PwString[] = "pw_const_supg_param";
  char PwdString[] = "pw_const_supg_param_deriv";
  char Name[] = "name";
  char Description[] = "space for pw constant parameters";
  TBaseCell *cell;

  // coefficients of the equation
  coeff = new double[13];
  // allocate objects and arrays for solution and right hand side of adjoint problem
  sol_adjoint = new double[N_U];
  memset(sol_adjoint, 0, N_U*SizeOfDouble);
  u_adjoint = new TFEFunction2D(velocity_space, UadjointString, UadjointString, sol_adjoint, N_U);
  rhs_edge = new double[N_U];
  // allocate matrix for adjoint problem
  sqmatrixAadjoint = new TSquareMatrix2D(sqstructureA);
  // allocate objects and array that contains the piecewise constant stabilization parameters
  pw_const_param = new double[N_Cells];
  memset(pw_const_param, 0, N_Cells*SizeOfDouble);
  TDatabase::ParamDB->INTERNAL_P1_Array = pw_const_param;
  // piecewise constant finite element space
  pw_const_param_space = new TFESpace2D(coll, Name, Description, BoundCondition, 0, NULL);
  pw_const_param_fe = new TFEFunction2D(pw_const_param_space, PwString, PwString, pw_const_param, N_Cells);
  // allocate objects and array that contains the piecewise constant derivatives of the stabilization parameters
  pw_const_param_deriv = new double[N_Cells];
  pw_const_param_deriv_fe = new TFEFunction2D(pw_const_param_space, 
					      PwdString, PwdString, 
					      pw_const_param_deriv, N_Cells);
  // array for the previous stabilization parameters
  pw_const_param_old = new double[N_Cells];
  
  // initialization: fill pw_const_param_fe with the standard stabilization parameters
  type = TDatabase::ParamDB->SDFEM_TYPE;
  TDatabase::ParamDB->SDFEM_TYPE = 2;
  // loop over all cells
  for (i=0;i<N_Cells;i++)
  {
      cell = coll->GetCell(i);
      N_Edges=cell->GetN_Edges();
      // prepare data for computing mesh size in convection direction
      for (j=0;j<N_Edges;j++)
      {
	  TDatabase::ParamDB->INTERNAL_VERTEX_X[j] = cell->GetVertex(j)->GetX();
	  TDatabase::ParamDB->INTERNAL_VERTEX_Y[j] = cell->GetVertex(j)->GetY();
      }
      if (N_Edges==3)
	  TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
      xs = ys = 0;
      for (j=0;j<N_Edges;j++)
      {
	  // compute coordinates
	  x = cell->GetVertex(j)->GetX();
	  y = cell->GetVertex(j)->GetY();
	  xs += x;
	  ys += y;
      }
      xs /= N_Edges;
      ys /= N_Edges;
      Coeffs(1, &xs, &ys, NULL, &coeff);
      pw_const_param[i] = 
	  Compute_SDFEM_delta(1, coeff[0], coeff[1], coeff[2], coeff[3], 0);
  }
  TDatabase::ParamDB->SDFEM_TYPE = type;
  delete coeff;
}

/******************************************************************************/
// SolveAdjointProblem2D
//  - assembles right hand side for adjoint problem
//  - solves adjoint problem
//  - computes derivatives of error estimator wrt stabilization parameters
/******************************************************************************/
void SolveAdjointProblem2D(TCollection *coll,
			   TDiscreteForm2D *DiscreteForm,
			   TFESpace2D *velocity_space,
			   TFEFunction2D *velo,
			   TFEFunction2D *u_adjoint,
			   TSquareMatrix2D* &sqmatrixAadjoint,
			   CoeffFct2D *Coeffs,
			   BoundCondFunct2D **BoundaryConditions,
			   BoundValueFunct2D **BoundaryValues,
			   double *rhs, double *rhs_edge,
			   double *sol_adjoint,
			   double *pw_const_param_deriv,
			   int N_U, int N_Active, int N_neum_to_diri,
			   int *neum_to_diri, int *neum_to_diri_bdry,
			   double *neum_to_diri_param)
{
    int solver_type, i;
    TAuxParam2D *aux;
    TSquareMatrix2D *SQMATRICES[1];
    TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
    MultiIndex2D FEMultiIndex_All_Deriv[5] = { D00, D10, D01, D20, D02 };
    MultiIndex2D FEMultiIndex_Sol[1] = { D00 };

    // assemble right hand side for adjoint problem
    // coercivity constant for rhs of adjoint problem
    TDatabase::ParamDB->INTERNAL_COERCIVITY = EstimateCoercivityConstant(coll, Coeffs);
    memset(rhs, 0, N_U*SizeOfDouble);
    // set aux object
    switch (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
    {
      case 100: 
	aux =  new TAuxParam2D(N_FESpaces_Sol, N_Fct_Sol, N_ParamFct_Sol,
			       N_FEValues_Sol, &velocity_space, &velo,
			       Fct_Sol,
			       FEFctIndex_Sol,  FEMultiIndex_Sol,
			       N_Params_Sol, BeginParam_Sol);
	break;
    default: 
      aux =  new TAuxParam2D(N_FESpaces_All_Deriv, N_Fct_All_Deriv, N_ParamFct_All_Deriv,
			     N_FEValues_All_Deriv, &velocity_space, &velo,
			     Fct_All_Deriv,
			     FEFctIndex_All_Deriv, FEMultiIndex_All_Deriv,
			     N_Params_All_Deriv, BeginParam_All_Deriv);
      break;
    }

    // this is just a dummy for the cell measure
    solver_type = TDatabase::ParamDB->CELL_MEASURE;
    // diameter
    TDatabase::ParamDB->CELL_MEASURE = 0;
    // assemble with homogeneous boundary conditions
    // IMPORTANT: bc. of primal problem to rule out Dirichlet boundaries, if set
    Assemble2D(1, &velocity_space,
	       0, NULL,
	       0, NULL,
	       1, &rhs, &velocity_space,
	       DiscreteForm,
	       BoundaryConditions,
	       BoundaryValues+2,
	       aux);
    // reset parameter
    TDatabase::ParamDB->CELL_MEASURE = solver_type;    
    // compute the contributions from the jumps across the edges
    // >= 100 -- adjoint problem with known error
    switch (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
      {
      case 1:
	memset(rhs_edge, 0, N_U*SizeOfDouble);
	JumpTermsForAdjointProblem(velocity_space, velo, Coeffs, BoundaryConditions[0], rhs_edge);
	// sum up both contributions
	Daxpy(N_U, 1, rhs, rhs_edge);
	break;
      default:
	memcpy(rhs_edge, rhs, N_U*SizeOfDouble);
	break;
      }

    // set Dirichlet boundary conditions
    SQMATRICES[0] = sqmatrixAadjoint;
    SetDirichletNodesFromNeumannNodes(SQMATRICES, rhs_edge, sol_adjoint,
				      N_neum_to_diri, neum_to_diri,
				      neum_to_diri_bdry, neum_to_diri_param,
				      BoundaryValues[2]);
    // this is just a dummy
    solver_type = TDatabase::ParamDB->SOLVER_TYPE;
    // change solver to direct solver
    TDatabase::ParamDB->SOLVER_TYPE = 2;
    Solver(sqmatrices, NULL,
	   rhs_edge, sol_adjoint,
	   NULL, NULL,
	   NULL, N_U, 0);
    // reset parameter
    TDatabase::ParamDB->SOLVER_TYPE = solver_type;

    // compute derivatives of error estimator wrt stabilization parameters
    ComputeDerivativeOfEstimator2D(coll, Coeffs, velo, u_adjoint, pw_const_param_deriv);
}

/******************************************************************************/
// ComputeDerivativeOfEstimator2D
/******************************************************************************/
void ComputeDerivativeOfEstimator2D(TCollection *coll,
				    CoeffFct2D *Coeffs,
				    TFEFunction2D *velo,
				    TFEFunction2D *u_adjoint,
				    double *pw_const_param_deriv)
{
  int i,j, N_Cells, N_Edges;
  double area, integral, xs, ys, x[4], y[4], val_velo[4], val_u_ad[4], *coeff;
  TBaseCell *cell;

  // coefficients of the equation
  coeff = new double[13];
  // number of mesh cells
  N_Cells = coll->GetN_Cells();
  // loop over all mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get cell no. i
    cell = coll->GetCell(i);
    // get number of edges
    N_Edges=cell->GetN_Edges();
    // get area of mesh cell
    area = cell->GetMeasure();
    // compute coordinates of vertices
    for (j=0;j<N_Edges;j++)
    {
      x[j] = cell->GetVertex(j)->GetX();
      y[j] = cell->GetVertex(j)->GetY();
    }
    // initialize integral
    integral = 0;
    // loop over the edges for edge mid point rule
    for (j=0;j<N_Edges;j++)
    {
	// barycenter of the edge
	xs = (x[j]+x[(j+1)%N_Edges])/2.0;
	ys = (y[j]+y[(j+1)%N_Edges])/2.0;
	// get value and derivatives of current solution of primal problem
	velo->FindGradientLocal(cell,i,xs,ys,val_velo);
	// get value and derivatives of solution of adjoint problem
	u_adjoint->FindGradientLocal(cell,i,xs,ys,val_u_ad);
	// get coefficients of the equation
	Coeffs(1, &xs, &ys, NULL, &coeff);
	// linearization of error estimator
	integral += (coeff[1]*val_velo[1] + coeff[2]*val_velo[2] + 
		     coeff[3]*val_velo[0] - coeff[4])*
	    (coeff[1]*val_u_ad[1] + coeff[2]*val_u_ad[2]);
    }
    // scaling for edge midpoint rule
    pw_const_param_deriv[i]  = -integral*area/3;
  }
  delete coeff;
}


void Assemble_CD_2D(TCollection *coll,
		    TFESpace2D **USpaces, TFEFunction2D **UArray,
		    double **RhsArray, TSquareMatrix2D **MatricesA,
		    CoeffFct2D *Coeffs,
		    BoundCondFunct2D **BoundaryConditions,
		    BoundValueFunct2D **BoundaryValues,
		    TDiscreteForm2D **DiscreteForms,
		    TSquareMatrix2D *sqmatrixAadjoint,
		    CheckWrongNeumannNodesFunct2D **CheckWrongNeumannNodesFct,
		    int *N_Uarray,
		    int low, int mg_level, int mg_type, int i,
		    int &N_neum_to_diri, int* &neum_to_diri,
		    int* &neum_to_diri_bdry, 
		    double* &neum_to_diri_param)
{
    int k, N_U, N_Active, N_NonActive, type;
    double *rhs, *RHSs[1];
    TFESpace2D *fesp[1], *ferhs[1];
    TSquareMatrix2D *SQMATRICES[1];
    TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
    TAuxParam2D *aux;
    TDiscreteForm2D *DiscreteForm;
    BoundCondFunct2D *BdCond;
    BoundValueFunct2D *BdVal;

    for(k=low;k<=mg_level;k++)
    {
      rhs = RhsArray[k];
      N_U = N_Uarray[k];
      N_Active = USpaces[k]->GetActiveBound();
      N_NonActive = N_U - N_Active;

      RHSs[0] = rhs;
      memset(rhs, 0, N_U*SizeOfDouble);

      // find discrete form
      if ((mg_type==1) && (k<i+1))
      {
        DiscreteForm = DiscreteForms[2];
	OutPut("Upwind ");
      }
      else
        switch(TDatabase::ParamDB->DISCTYPE)
        {
          case GALERKIN:
          case LOCAL_PROJECTION:
	      DiscreteForm = DiscreteForms[0];
	      OutPut("Galerkin ");
	      break;

          case SDFEM:
            DiscreteForm = DiscreteForms[1];
	    type = TDatabase::ParamDB->SDFEM_TYPE;
	    if ((TDatabase::ParamDB->SDFEM_TYPE == 100)&&(k<mg_level))
	    {
		TDatabase::ParamDB->SDFEM_TYPE = 2;
	    }
            break;

          case UPWIND:
	      DiscreteForm = DiscreteForms[2];
	      OutPut("Upwind ");
	      break;

        default:
          OutPut("Unknown DISCTYPE" << endl);
          exit(4711);;
      }

      fesp[0] = USpaces[k];
      ferhs[0] = USpaces[k];
      // initialize matrix
      SQMATRICES[0] = MatricesA[k];
      SQMATRICES[0]->Reset();

      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      SetPolynomialDegree();
      BdCond = BoundaryConditions[0];
      BdVal = BoundaryValues[0];
      // set homogeneous Neumann boundary conditions
      if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)
	  || (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM))
      {
	  BdCond = BoundaryConditions[1];
	  BdVal = BoundaryValues[1];
      }
      // assemble
      Assemble2D(1, fesp,
		 1, SQMATRICES,
		 0, NULL,
		 1, RHSs, ferhs,
		 DiscreteForm,
		 &BdCond,
		 &BdVal,
		 aux);
      if (TDatabase::ParamDB->SC_VERBOSE>1)
	OutPut("Assembling done on mg-level: "<< mg_level<<endl);

      // DG and CIP
      switch (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)
	{
	case 1:  Assemble2D_CIP(Coeffs,1, fesp,
				1, SQMATRICES,
				0, NULL,
				1, RHSs, ferhs,
				BoundaryConditions,
				BoundaryValues,
				aux);
	  break;
	case 2:  Assemble2D_DG(Coeffs,1, fesp,
			       1, SQMATRICES,
			       0, NULL,
			       1, RHSs, ferhs,
			       BoundaryConditions,
			       BoundaryValues,
			       aux);
	  break;
	}

      // reset SUPG parameter
      if (TDatabase::ParamDB->DISCTYPE == SDFEM)
      {
	  TDatabase::ParamDB->SDFEM_TYPE = type;
      }	  
      // LPS ???
      if(TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
      {
	  if(TDatabase::ParamDB->LP_FULL_GRADIENT)
	      ;//UltraLocalProjection(MatricesA[k], FALSE, Coefficients[0]);
	  
	  if(TDatabase::ParamDB->LP_STREAMLINE)
	  {
	      OutPut("local projection stabilisation in streamline direction ");
	      OutPut("is currently not available." << endl);
	      exit(4711);
	  }
      }
      // apply upwind
      if (DiscreteForm == DiscreteForms[2])
      {
	  UpwindForConvDiff(Coeffs, SQMATRICES[0],RHSs[0],
			    fesp[0],DiscreteForms[2],NULL,NULL,0);
	  cout << "UPWINDING DONE : level " << k << endl;
      }                                           // endif
      delete aux;

      // compute adjoint matrix, only on the fines level
      if ((TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&(k==mg_level))
      {
	  ComputeAdjointMatrix(MatricesA[k],sqmatrixAadjoint);
      }
      if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)
	  || (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM))
      {
	  CheckWrongNeumannNodesFct[0](coll, USpaces[k], N_neum_to_diri, neum_to_diri,
				       neum_to_diri_bdry, neum_to_diri_param);
	  if (N_neum_to_diri)
	      SetDirichletNodesFromNeumannNodes(SQMATRICES, RHSs[0], UArray[k]->GetValues(),
						N_neum_to_diri, neum_to_diri,
						neum_to_diri_bdry, neum_to_diri_param,
						BoundaryValues[0]);
      }
    }	  
}

/******************************************************************************/
// restrict parameters to admitted values 
/******************************************************************************/
void  RestrictParametersToAdmissibleValues(TCollection *coll,
					   TFESpace2D *pw_const_param_space,
                                           CoeffFct2D *Coeffs,
					   double* pw_const_param)
{
    int i, j, N_Cells, N_Edges, index, order;
    int *GlobalNumbers, *BeginIndex, *DOF;
  double omega = TDatabase::ParamDB->INTERNAL_COERCIVITY, val1, val2, x, y;
  double *coeff, c_inv2 = 1, hK, lower_bound;
  TBaseCell *cell;
  FE2D CurrentElement;

  N_Cells = coll->GetN_Cells();
  GlobalNumbers = pw_const_param_space->GetGlobalNumbers();
  BeginIndex = pw_const_param_space->GetBeginIndex();
  order = TDatabase::ParamDB->ANSATZ_ORDER;

  coeff = new double[13];
  // loop over the cells
  for (i=0;i<N_Cells;i++)
  {
      cell = coll->GetCell(i);
      hK = cell->GetDiameter();
      lower_bound = 1e-12 * hK * hK;
      // dof of current mesh cell
      DOF = GlobalNumbers + BeginIndex[i];
      index = DOF[0];
      // apply lower bound
      if (pw_const_param[index]<=lower_bound)
      {
	  pw_const_param[index] = lower_bound;
	  continue;
      }
      // compute maximal value
      // compute first barycenter
      N_Edges = cell->GetN_Edges();
      x = y = 0;
      for (j=0; j<N_Edges; j++)
      {
	  x += cell->GetVertex(j)->GetX();
	  y += cell->GetVertex(j)->GetY();
      }
      x /= N_Edges;
      y /= N_Edges;
      // get coefficients
      Coeffs(1, &x, &y, NULL, &coeff);
      // finite element on the mesh cell
      CurrentElement = pw_const_param_space->GetFE2D(i, cell);
      switch(CurrentElement)
      {
	  // P_1, Q_1
	  case C_P1_2D_T_A:
	  case C_Q1_2D_Q_A:
	  case C_Q1_2D_Q_M:
	      c_inv2 = 1.0;
	      break;
	case C_P2_2D_T_A:
	    c_inv2 = 48.0;
	    break;
	case C_Q2_2D_Q_A:
	case C_Q2_2D_Q_M:
	    c_inv2 = 24.0;
	    break;
	case C_P3_2D_T_A:
	    c_inv2 = (435+sqrt(26025.0))/4.0;
	    break;
	case C_Q3_2D_Q_A:
	case C_Q3_2D_Q_M:
	    c_inv2 = (244+sqrt(9136.0))/3.0;
	    break;
	default:
	    c_inv2 = (435+sqrt(26025.0))/4.0;
	    break;
      }
	    
      // first parameter
      val1 = hK * hK;
      val1 /= (coeff[0]*c_inv2);
      // this bound not necessary for linear and bilinear fe (on parallelograms)
      if ((CurrentElement==C_P1_2D_T_A) || (CurrentElement==C_Q1_2D_Q_A) || (CurrentElement==C_Q1_2D_Q_M))
	  val1 = 1e36;
      // there is a second bound for omega > 0
      if (omega > 0)
      {
	  // second parameter
	  if (coeff[3]!=0)
	      val2 = omega/(coeff[3]*coeff[3]);
	  else
	      val2 = 1e36;
	  // put minimum on val1
	  if (val2 < val1)
	      val1 = val2/2.0;
	  else
	      val1 /= 2.0;
      }
      // apply upper bound
      if (pw_const_param[index] > val1)
      {
	  pw_const_param[index] = val1;
      }
  }
  delete coeff;
}

/******************************************************************************/
// computes the wrong Neumann nodes for problems on the unit square with
// Dirichlet boundary conditions
/******************************************************************************/

void CheckWrongNeumannNodesUnitSquareDiri(TCollection *Coll, TFESpace2D *fespace,
			    int &N_neum_to_diri, int* &neum_to_diri,
			    int* &neum_to_diri_bdry, 
			    double* &neum_to_diri_param)
{
    const int max_entries = 4096;  
    int i, j, N_, min_val, type;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[4], tmp_diri[max_entries], tmp_bdry[max_entries];
  double x[4], y[4], eps = 1e-6, tmp_param[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE2D CurrentElement;

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();

  diri_counter = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    found = 0;
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      boundary_vertices[j] = 0;
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j]);
      // vertex on the upper lid
      if ((fabs(x[j])<eps)||(fabs(y[j])<eps)||(fabs(1-x[j])<eps)||(fabs(1-y[j])<eps))
      {
	   boundary_vertices[j] = 1;
	   found++;
      }
    }
    // no cell with edge with vertex on the boundary
    if (found<2) 
	continue;
    // finite element on the mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = TFEDatabase2D::GetN_BaseFunctFromFE2D(CurrentElement);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = global_numbers+begin_index[i];
    switch(CurrentElement)
    {
	// P_1, Q_1
	case C_P1_2D_T_A:
	case C_Q1_2D_Q_A:
	case C_Q1_2D_Q_M:
	    for (j=0;j<N_V;j++)
	    {
		// vertex on the boundary
		if (boundary_vertices[j])
		{
		    if (CurrentElement==C_P1_2D_T_A)
			tmp_diri[diri_counter] = dof[j];
		    else
		    {
			if (j<2)
			    tmp_diri[diri_counter] = dof[j];
			else
			{
			    if (j==2)
				tmp_diri[diri_counter] = dof[3];
			    else
				tmp_diri[diri_counter] = dof[2];
			}
		    }
		    if (diri_counter > max_entries)
		    {
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		    }
		    if (fabs(y[j])<eps)
		    {
			tmp_bdry[diri_counter] = 0;
			tmp_param[diri_counter] = x[j];
		    }
		    if (fabs(1-y[j])<eps)
		    {
			tmp_bdry[diri_counter] = 2;
			tmp_param[diri_counter] = 1-x[j];
		    }
		    if (fabs(x[j])<eps)
		    {
			tmp_bdry[diri_counter] = 3;
			tmp_param[diri_counter] = 1-y[j];
		    }
		    if (fabs(1-x[j])<eps)
		    {
			tmp_bdry[diri_counter] = 1;
			tmp_param[diri_counter] = y[j];
		    }
		    diri_counter++;
		}
	    }
	    break;
	// P_2, Q_2
	case C_P2_2D_T_A:
	case C_Q2_2D_Q_A:
	case C_Q2_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[4];
                       tmp_diri[diri_counter+2] = dof[5];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[5];
                       tmp_diri[diri_counter+2] = dof[8];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[5];
                       tmp_diri[diri_counter+1] = dof[3];
                       tmp_diri[diri_counter+2] = dof[0];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[8];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[6];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[6];
                     tmp_diri[diri_counter+1] = dof[3];
                     tmp_diri[diri_counter+2] = dof[0];
                   break;

                }
              
		if (diri_counter+2 > max_entries)
		{
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		}
		// boundary at y=0
		if ((fabs(y[j])<eps)&&(fabs(y[(j+1)%N_V])<eps))
		  {
		    tmp_bdry[diri_counter] = 0;
		    tmp_bdry[diri_counter+1] = 0;
		    tmp_bdry[diri_counter+2] = 0;
		    tmp_param[diri_counter] = x[j];
		    tmp_param[diri_counter+1] = (x[j] + x[(j+1)%N_V])/2.0;
		    tmp_param[diri_counter+2] = x[(j+1)%N_V];
		  }
		// boundary at y = 1
		if ((fabs(1-y[j])<eps)&&(fabs(1-y[(j+1)%N_V])<eps))
		  {
		    tmp_bdry[diri_counter] = 2;
		    tmp_bdry[diri_counter+1] = 2;
		    tmp_bdry[diri_counter+2] = 2;
		    tmp_param[diri_counter] = 1-x[j];
		    tmp_param[diri_counter+1] = (1-x[j] + 1-x[(j+1)%N_V])/2.0;
		    tmp_param[diri_counter+2] = 1-x[(j+1)%N_V];
		  }
		// boundary at x = 0
		if ((fabs(x[j])<eps)&&(fabs(x[(j+1)%N_V])<eps))
		  {
		    tmp_bdry[diri_counter] = 3;
		    tmp_bdry[diri_counter+1] = 3;
		    tmp_bdry[diri_counter+2] = 3;
		    tmp_param[diri_counter] = 1-y[j];
		    tmp_param[diri_counter+1] = (1-y[j] + 1-y[(j+1)%N_V])/2.0;
		    tmp_param[diri_counter+2] = 1-y[(j+1)%N_V];
		  }
		// boundary at x = 1
		if ((fabs(1-x[j])<eps)&&(fabs(1-x[(j+1)%N_V])<eps))
		  {
		    tmp_bdry[diri_counter] = 1;
		    tmp_bdry[diri_counter+1] = 1;
		    tmp_bdry[diri_counter+2] = 1;
		    tmp_param[diri_counter] = y[j];
		    tmp_param[diri_counter+1] = (y[j] + y[(j+1)%N_V])/2.0;
		    tmp_param[diri_counter+2] = y[(j+1)%N_V];
		  }
		diri_counter +=3;
	      }
	    }
	    break;
	// P_3, Q_3
	case C_P3_2D_T_A:
	case C_Q3_2D_Q_A:
	case C_Q3_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;

               // P3: local dof 0, 1, 2, 3 are on the boundary
               // Q3: local dof 0, 1, 2, 3 are on the boundary
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
		     tmp_diri[diri_counter+3] = dof[3];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[6];
                       tmp_diri[diri_counter+2] = dof[8];
		       tmp_diri[diri_counter+3] = dof[9];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[11];
		       tmp_diri[diri_counter+3] = dof[15];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[9];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[4];
                       tmp_diri[diri_counter+3] = dof[0];
		     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[15];
                       tmp_diri[diri_counter+1] = dof[14];
                       tmp_diri[diri_counter+2] = dof[13];
			tmp_diri[diri_counter+3] = dof[12];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[12];
                     tmp_diri[diri_counter+1] = dof[8];
                     tmp_diri[diri_counter+2] = dof[4];
		     tmp_diri[diri_counter+3] = dof[0];
                   break;
                }
              
		if (diri_counter+3 > max_entries)
		{
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		}

		if ((fabs(y[j])<eps)&&(fabs(y[(j+1)%N_V])<eps))
		{
		    tmp_bdry[diri_counter] = 0;
		    tmp_bdry[diri_counter+1] = 0;
		    tmp_bdry[diri_counter+2] = 0;
		    tmp_bdry[diri_counter+3] = 0;
		    tmp_param[diri_counter] = x[j];
		    tmp_param[diri_counter+1] = 2*x[j]/3.0 + x[(j+1)%N_V]/3.0;
		    tmp_param[diri_counter+2] = x[j]/3.0 + 2*x[(j+1)%N_V]/3.0; 
		    tmp_param[diri_counter+3]= x[(j+1)%N_V];
		}
		if ((fabs(1-y[j])<eps)&&(fabs(1-y[(j+1)%N_V])<eps))
		{
		    tmp_bdry[diri_counter] = 2;
		    tmp_bdry[diri_counter+1] = 2;
		    tmp_bdry[diri_counter+2] = 2;
		    tmp_bdry[diri_counter+3] = 2;
		    tmp_param[diri_counter] = 1-x[j];
		    tmp_param[diri_counter+1] = 2*(1-x[j])/3.0 + (1-x[(j+1)%N_V])/3.0;
		    tmp_param[diri_counter+2] = (1-x[j])/3.0 + 2*(1-x[(j+1)%N_V])/3.0;
		    tmp_param[diri_counter+3] = 1-x[(j+1)%N_V];
		}
		    if ((fabs(x[j])<eps)&&(fabs(x[(j+1)%N_V])<eps))
		    {
			tmp_bdry[diri_counter] = 3;
			tmp_bdry[diri_counter+1] = 3;
			tmp_bdry[diri_counter+2] = 3;
			tmp_bdry[diri_counter+3] = 3;
			tmp_param[diri_counter] = 1-y[j];
			tmp_param[diri_counter+1] = 2*(1-y[j])/3.0 + (1-y[(j+1)%N_V])/3.0;
			tmp_param[diri_counter+2] = (1-y[j])/3.0 + (1-y[(j+1)%N_V])/3.0;
			tmp_param[diri_counter+3] = 1-y[(j+1)%N_V];
		    }
		    if ((fabs(1-x[j])<eps)&&(fabs(1-x[(j+1)%N_V])<eps))
		    {
			tmp_bdry[diri_counter] = 1;
			tmp_bdry[diri_counter+1] = 1;
			tmp_bdry[diri_counter+2] = 1;
			tmp_bdry[diri_counter+3] = 1;
			tmp_param[diri_counter] = y[j];
			tmp_param[diri_counter+1] = 2*y[j]/3.0 + y[(j+1)%N_V]/3.0;
			tmp_param[diri_counter+2] = y[j]/3.0 + 2*y[(j+1)%N_V]/3.0;
			tmp_param[diri_counter+3] = y[(j+1)%N_V];
		    }
		    diri_counter +=4;
		}
	    }
	    break;
	default:
	    OutPut("CheckNeumannNodesForVelocity not implemented for element "
		   << CurrentElement << endl);
	    OutPut("code can be run without CheckNeumannNodesForVelocity, just delete the exit" << endl);
	    exit(4711);
    }	    
  }
 
  // condense
  for (i=0;i<diri_counter;i++)
  {
      if (tmp_diri[i] == -1)
	  continue;
      diri_counter_1++;
      for (j=i+1;j<diri_counter;j++)
      {
	  if (tmp_diri[i] == tmp_diri[j])
	  {
	      tmp_diri[j] = -1;
	  }
      }
  }

  //OutPut("CheckNeumannNodesForVelocity: N_neum_to_diri " << diri_counter_1 << endl);
  N_neum_to_diri = diri_counter_1;
  // allocate array for the indices
  neum_to_diri = new int[diri_counter_1];
  // allocate array for the corresponding boundary numbers
  neum_to_diri_bdry = new int[diri_counter_1];
  // allocate array for the corresponding boundary parameters
  neum_to_diri_param = new double[diri_counter_1];
  // fill array and sort
  for (i=0;i<diri_counter_1;i++)
  {
      min_val = tmp_diri[0];
      found = 0;
      for (j=1;j<diri_counter;j++)
      {
	  if ((tmp_diri[j]>0) && ((tmp_diri[j] < min_val) || 
				  (min_val == -1)))
	  {
	       min_val =  tmp_diri[j];
	       found = j;
	  }
      }
      neum_to_diri[i] = tmp_diri[found];
      neum_to_diri_bdry[i] = tmp_bdry[found];
      neum_to_diri_param[i] = tmp_param[found];
      tmp_diri[found] = -1;
  }
}
