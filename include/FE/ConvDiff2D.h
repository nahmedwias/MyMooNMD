// ======================================================================
// @(#)ConvDiff2D.h        1.13 10/19/99
//
// common declaration for all convection diffusion problems
// ======================================================================

#ifndef __CONVDIFF2D__
#define __CONVDIFF2D__


// part for standard Galerkin
int N_Terms = 3;
MultiIndex2D Derivatives[3] = { D10, D01, D00 };
int SpacesNumbers[3] = { 0, 0, 0 };

// part for SDFEM
int N_Terms_SD = 5;
MultiIndex2D Derivatives_SD[5] = { D10, D01, D00, D20, D02 };
int SpacesNumbers_SD[5] = { 0, 0, 0, 0, 0 };

// part for UPWIND with lumping of reaction term and rhs
int N_Terms_UPW1 = 2;
MultiIndex2D Derivatives_UPW1[2] = { D10, D01 };
int SpacesNumbers_UPW1[2] = { 0, 0 };

// part for UPWIND without lumping of reaction term and rhs
int N_Terms_UPW2 = 3;
MultiIndex2D Derivatives_UPW2[3] = { D10, D01, D00 };
int SpacesNumbers_UPW2[3] = { 0, 0, 0 };

// part for rhs
int N_Terms_rhs = 1;
MultiIndex2D Derivatives_rhs[1] = { D00 };
int SpacesNumbers_rhs[1] = { 0 };

// part for all
int CD_N_Matrices = 1;
int CD_RowSpace[1] = { 0 };
int CD_ColumnSpace[1] = { 0 };
int CD_N_Rhs = 1;
int CD_RhsSpace[1] = { 0 };

MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
MultiIndex2D ZeroDerivative[1] = { D00};

void BilinearAssemble(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void BilinearAssemble_Axial3D(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);


void BilinearAssemble_SD(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void BilinearAssemble_GLS(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void BilinearAssemble_UPW1(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void BilinearAssemble_UPW2(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SOLD(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SOLD_Orthogonal(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void RhsAssemble_LP96(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
		      double ***LocMatrices, double **LocRhs);

void BilinearAssemble_MH_Kno06(double Mult, double *coeff, double *param,
			       double hK,
			       double **OrigValues, int *N_BaseFuncts,
			       double ***LocMatrices, double **LocRhs);

void RhsAssemble_RhsAdjointEnergyEstimate(double Mult, double *coeff, double *param,
					  double hK,
					  double **OrigValues, int *N_BaseFuncts,
					  double ***LocMatrices, double **LocRhs);

void RhsAssemble_RhsAdjointL2Error(double Mult, double *coeff, double *param,
				   double hK,
				   double **OrigValues, int *N_BaseFuncts,
				   double ***LocMatrices, double **LocRhs);

void RhsAssemble_RhsAdjointH1Error(double Mult, double *coeff, double *param,
				   double hK,
				   double **OrigValues, int *N_BaseFuncts,
				   double ***LocMatrices, double **LocRhs);

void EdgeStabilization(TFESpace2D *fespace, 
		       TFEFunction2D *u, 
		       CoeffFct2D *Coeffs,
		       double *rhs,
		       int time_dependent,
		       double *time_step,
		       TFEFunction2D *old_u);

// parameters for DC/CD shock capturing scheme
void DC_CD_Params(double *in, double *out);

int DC_CD_N_FESpaces = 1;
int DC_CD_N_Fct = 2;
int DC_CD_N_ParamFct = 1;
int DC_CD_N_FEValues = 2;
int DC_CD_N_Params = 2;
int DC_CD_FEFctIndex[2] = { 0, 1 };
MultiIndex2D DC_CD_FEMultiIndex[2] = { D00, D00 };
ParamFct *DC_CD_Fct[1] = { DC_CD_Params };
int DC_CD_BeginParam[1] = { 0 };

// parameters for SC_2 shock capturing scheme
void SC_2_Params(double *in, double *out);

int SC_2_N_FESpaces = 2;
int SC_2_N_Fct = 3;
int SC_2_N_ParamFct = 1;
int SC_2_N_FEValues = 4;
int SC_2_N_Params = 4;
int SC_2_FEFctIndex[4] = { 0, 1, 2, 2 };
MultiIndex2D SC_2_FEMultiIndex[4] = { D00, D00, D10, D01 };
ParamFct *SC_2_Fct[1] = { SC_2_Params };
int SC_2_BeginParam[1] = { 0 };

// parameters for SOLD schemes
void SOLD_Params(double *in, double *out);

int SOLD_N_FESpaces = 2;
int SOLD_N_Fct = 3;
int SOLD_N_ParamFct = 1;
int SOLD_N_FEValues = 7;
int SOLD_N_Params = 7;
int SOLD_FEFctIndex[7] = { 0, 0, 0, 0, 0, 1, 2 };
MultiIndex2D SOLD_FEMultiIndex[7] = { D00, D10, D01, D20, D02, D00, D00 };
ParamFct *SOLD_Fct[1] = { SOLD_Params };
int SOLD_BeginParam[1] = { 0 };

// ========================================================================
// write parameters for latex table
// ========================================================================

void ParamsForLatexTable();
void SetSoldParameters(int i);
void EvaluateResults(int *nonlinite, double *scale, double *globmin,
double *globmax, double *bdrmin, double *bdrmax,
double *bdrh11, double *maxprinciple_cells, double *maxprinciple_max,
int ij_end);
void EvaluateResults(int *nonlinite, double *scale, double *globmin,
double *globmax, double *intlayer_width, double *bdrlay_osc,
double *bdrlay_smear, int ij_end);

#endif
