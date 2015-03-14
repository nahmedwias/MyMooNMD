#ifndef __NSE2D_PARAMROUT__
#define __NSE2D_PARAMROUT__
#include <NSE2D_FixPo.h>

// ======================================================================
// setting for error calculation for all types
// ======================================================================
MultiIndex2D NSAllDerivatives[3] = { D00, D10, D01 };
MultiIndex2D NSErrorEstiamtorU_Derivatives[8] = { D10, D01, D00, D20, D02, 
                                         D10, D01, D00  };
MultiIndex2D NSErrorEstiamtorU_DerivativesEstimator[5] = { D10, D01, D00, D20, D02};
MultiIndex2D NSErrorEstiamtorP_Derivatives[3] = { D10, D01, D00 };
MultiIndex2D NSErrorEstiamtorP_DerivativesEstimator[2] = { D10, D01 };

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
int NSN_FESpacesVelo = 1;
int NSN_FctVelo = 2;
int NSN_ParamFctVelo = 1;
int NSN_FEValuesVelo = 2;
int NSN_ParamsVelo = 2;
int NSFEFctIndexVelo[2] = { 0, 1 };
MultiIndex2D NSFEMultiIndexVelo[2] = { D00, D00 };
ParamFct *NSFctVelo[1] = { NSParamsVelo };
int NSBeginParamVelo[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void NSParamsVelo_GradVelo(double *in, double *out);

int NSN_FESpacesVelo_GradVelo = 1;
int NSN_FctVelo_GradVelo = 2;
int NSN_ParamFctVelo_GradVelo = 1;
int NSN_FEValuesVelo_GradVelo = 6;
int NSN_ParamsVelo_GradVelo = 6;
int NSFEFctIndexVelo_GradVelo[6] = { 0, 1, 0, 1, 0, 1 };
MultiIndex2D NSFEMultiIndexVelo_GradVelo[6] = { D00, D00, D10, D10, D01, D01 };
ParamFct *NSFctVelo_GradVelo[1] = { NSParamsVelo_GradVelo };
int NSBeginParamVelo_GradVelo[1] = { 0 };

// ======================================================================
// auxiliary problem
// ======================================================================
void NSAuxProblem(double Mult, double *coeff, 
                   double *param, double hK, 
                   double **OrigValues, int *N_BaseFuncts,
                   double ***LocMatrices, double **LocRhs);


int NSAuxProblemN_Terms = 3;
MultiIndex2D NSAuxProblemDerivatives[3] = { D10, D01, D00};
int NSAuxProblemSpaceNumbers[3] = { 0, 0, 0};
int NSAuxProblemN_Matrices = 1;
int NSAuxProblemRowSpace[1] = { 0 };
int NSAuxProblemColumnSpace[1] = { 0 };
int NSAuxProblemN_Rhs = 2;
int NSAuxProblemRhsSpace[2] = { 0, 0 };

void NSParamsVeloExact(double *in, double *out);
ParamFct *NSFctVeloExact[1] = { NSParamsVeloExact };

// ======================================================================
// auxiliary problem for differential filter
// ======================================================================
void Filter_Galerkin(double Mult, double *coeff, 
                     double *param, double hK, 
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);
// ======================================================================
//  declarations for auxiliary problem for differential filter
//      one matrix 
//      two rhs
// ======================================================================

int Filter_N_Terms = 3;
MultiIndex2D Filter_Derivatives[3] = { D10, D01, D00};
int Filter_SpaceNumbers[3] = { 0, 0, 0};
int Filter_N_Matrices = 1;
int Filter_RowSpace[1] = { 0 };
int Filter_ColumnSpace[1] = { 0 };
int Filter_N_Rhs = 2;
int Filter_RhsSpace[2] = { 0, 0 };

// ========================================================================
// parameters: separated pressure
// ========================================================================
void NSParamsPressSep(double *in, double *out);

int NSN_FESpacesPressSep = 1;
int NSN_FctPressSep = 1;
int NSN_ParamFctPressSep = 1;
int NSN_FEValuesPressSep = 2;
int NSN_ParamsPressSep = 2;
int NSFEFctIndexPressSep[2] = { 0, 0 };
MultiIndex2D NSFEMultiIndexPressSep[2] = { D10, D01 };
ParamFct *NSFctPressSep[1] = { NSParamsPressSep };
int NSBeginParamPressSep[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
void NSParamsVeloAxialSymm3D(double *in, double *out);

int NSN_FESpacesVeloAxialSymm3D = 1;
int NSN_FctVeloAxialSymm3D = 2;
int NSN_ParamFctVeloAxialSymm3D = 1;
int NSN_FEValuesVeloAxialSymm3D = 2;
int NSN_ParamsVeloAxialSymm3D = 3;
int NSFEFctIndexVeloAxialSymm3D[2] = { 0, 1 };
MultiIndex2D NSFEMultiIndexVeloAxialSymm3D[2] = { D00, D00 };
ParamFct *NSFctVeloAxialSymm3D[1] = { NSParamsVeloAxialSymm3D };
int NSBeginParamVeloAxialSymm3D[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void NSParamsVelo_GradVeloAxialSymm3D(double *in, double *out);

int NSN_FESpacesVelo_GradVeloAxialSymm3D = 1;
int NSN_FctVelo_GradVeloAxialSymm3D = 2;
int NSN_ParamFctVelo_GradVeloAxialSymm3D = 1;
int NSN_FEValuesVelo_GradVeloAxialSymm3D = 6;
int NSN_ParamsVelo_GradVeloAxialSymm3D = 7;
int NSFEFctIndexVelo_GradVeloAxialSymm3D[6] = { 0, 1, 0, 1, 0, 1 };
MultiIndex2D NSFEMultiIndexVelo_GradVeloAxialSymm3D[6] = { D00, D00, D10, D10, D01, D01 };
ParamFct *NSFctVelo_GradVeloAxialSymm3D[1] = { NSParamsVelo_GradVeloAxialSymm3D };
int NSBeginParamVelo_GradVeloAxialSymm3D[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void NSParamsVelo_Press(double *in, double *out);

int NSN_FESpacesVelo_Press = 2;
int NSN_FctVelo_Press = 3;
int NSN_ParamFctVelo_Press = 1;
int NSN_FEValuesVelo_Press = 4;
int NSN_ParamsVelo_Press = 4;
int NSFEFctIndexVelo_Press[4] = { 0, 1, 2, 2 };
MultiIndex2D NSFEMultiIndexVelo_Press[6] = { D00, D00, D10, D10};
ParamFct *NSFctVelo_Press[1] = { NSParamsVelo_Press};
int NSBeginParamVelo_Press[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, temperature
// ========================================================================
/*
void NSParamsVeloTemp(double *in, double *out);

int NSN_FESpacesVeloTemp = 3;
int NSN_FctVeloTemp = 3;
int NSN_ParamFctVeloTemp = 1;
int NSN_FEValuesVeloTemp = 3;
int NSN_ParamsVeloTemp = 3;
int NSFEFctIndexVeloTemp[3] = { 0, 1, 2 };
MultiIndex2D NSFEMultiIndexVeloTemp[3] = { D00, D00, D00 };
ParamFct *NSFctVeloTemp[1] = { NSParamsVeloTemp };
int NSBeginParamVeloTemp[1] = { 0 };
*/

#endif // __NSE2D_PARAMROUT__
