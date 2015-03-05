// =======================================================================
// @(#)MainRoutines.h        
//
// Purpose: contains routines which are called from the main program 
//
// Author: Volker John 
//
// History: start of implementation 22.09.2009
//
// =======================================================================

#ifndef __MAINROUTINES2D__
#define __MAINROUTINES2D__

void SetParametersCDAdapt2D();

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
		    double* &neum_to_diri_param);

void Solver(TSquareMatrix **sqmatrices, TMatrix **matrices,
	    double *rhs, double *sol, 
	    MatVecProc *MatVect, DefectProc *Defect,
	    TMultiGrid2D *MG, 
	    int N_Unknowns, int ns_type);

void OutputData2D(std::ostringstream& os, TOutput2D *Output, int counter);

void ComputeErrorEstimate(TCollection *coll, TFEFunction2D *u,
			  CoeffFct2D *Coeffs, BoundCondFunct2D **BoundaryConditions,
			  BoundValueFunct2D **BoundaryValues,
			  double* &eta_K,
			  double *maximal_local_error,
			  double *estimated_global_error,
			  double l2, double h1, int N_Unknowns);

void ComputeAdjointMatrix(TSquareMatrix2D *A,TSquareMatrix2D *AT);


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
			     int N_U, int N_Cells);

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
			   double *neum_to_diri_param);

void ComputeDerivativeOfEstimator2D(TCollection *coll,
				    CoeffFct2D *Coeffs,
				    TFEFunction2D *velo,
				    TFEFunction2D *u_adjoint,
				    double *pw_const_param_deriv);

void  RestrictParametersToAdmissibleValues(TCollection *coll,
					   TFESpace2D *pw_const_param_space,
                                           CoeffFct2D *Coeffs,
					   double* pw_const_param);

void CheckWrongNeumannNodesUnitSquareDiri(TCollection *Coll, TFESpace2D *fespace,
					  int &N_neum_to_diri, int* &neum_to_diri,
					  int* &neum_to_diri_bdry, 
					  double* &neum_to_diri_param);

					  /** all derivatives of current solution as params */
void Params_All_Deriv(double *in, double *out);

int N_FESpaces_All_Deriv = 1;
int N_Fct_All_Deriv = 1;
int N_ParamFct_All_Deriv = 1;
int N_FEValues_All_Deriv = 5;
int N_Params_All_Deriv = 5;
int FEFctIndex_All_Deriv[5] = { 0, 0, 0, 0, 0 };
ParamFct *Fct_All_Deriv[1] = { Params_All_Deriv };
int BeginParam_All_Deriv[1] = { 0 };

/** all derivatives of current solution as params */
void Params_Sol(double *in, double *out);

int N_FESpaces_Sol = 1;
int N_Fct_Sol = 1;
int N_ParamFct_Sol = 1;
int N_FEValues_Sol = 1;
int N_Params_Sol = 1;
int FEFctIndex_Sol[1] = { 0 };
ParamFct *Fct_Sol[1] = { Params_Sol };
int BeginParam_Sol[1] = { 0 };

#endif
