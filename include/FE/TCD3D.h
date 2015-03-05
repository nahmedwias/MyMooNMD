// ======================================================================
// TCD3D.h
//
// common declaration for time dependent convection diffusion problems
// ======================================================================

#ifndef __TIMECONVDIFF3D__
#define __TIMECONVDIFF3D__

// MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };

// ======================================================================
// definitions for assembling the mass matrix and rhs
// ======================================================================
int N_Terms_MatrixMRhs = 1;
MultiIndex3D Derivatives_MatrixMRhs[1] = { D000 };
int SpacesNumbers_MatrixMRhs[1] = { 0 };
int N_Matrices_MatrixMRhs = 1;
int RowSpace_MatrixMRhs[1] = { 0 };
int ColumnSpace_MatrixMRhs[1] = { 0 };
int N_Rhs_MatrixMRhs = 1;
int RhsSpace_MatrixMRhs[1] = { 0 };

void MatrixMRhsAssemble(double Mult, double *coeff, double *param,
                           double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);

int N_Terms_MatrixMRhs_SUPG = 4;
MultiIndex3D Derivatives_MatrixMRhs_SUPG[4] = { D100, D010, D001, D000 };
int SpacesNumbers_MatrixMRhs_SUPG[4] = { 0, 0, 0, 0 };
int N_Matrices_MatrixMRhs_SUPG = 1;
int RowSpace_MatrixMRhs_SUPG[1] = { 0 };
int ColumnSpace_MatrixMRhs_SUPG[1] = { 0 };
int N_Rhs_MatrixMRhs_SUPG = 1;
int RhsSpace_MatrixMRhs_SUPG[1] = { 0 };

int N_Matrices_MatricesAKRhs_SOLD = 3;
int RowSpace_MatricesAKRhs_SOLD[3] = { 0, 0, 0 };
int ColumnSpace_MatricesAKRhs_SOLD[3] = { 0, 0, 0 };

void MatrixMRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the matrices A, K and rhs
// ======================================================================
int N_Terms_MatricesAKRhs_SUPG = 4;
MultiIndex3D Derivatives_MatricesAKRhs_SUPG[4] = { D100, D010, D001, D000 };
int SpacesNumbers_MatricesAKRhs_SUPG[4] = { 0, 0, 0, 0 };
int N_Matrices_MatricesAKRhs_SUPG = 2;
int RowSpace_MatricesAKRhs_SUPG[2] = { 0, 0 };
int ColumnSpace_MatricesAKRhs_SUPG[2] = { 0, 0 };
int N_Rhs_MatricesAKRhs_SUPG = 1;
int RhsSpace_MatricesAKRhs_SUPG[1] = { 0 };

void MatricesAKRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);

void RhsAssemble_SUPG(double Mult, double *coeff, double *param,
                      double hK, 
                      double **OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the matrix A and rhs
// ======================================================================
int N_Terms_MatrixARhs = 4;
MultiIndex3D Derivatives_MatrixARhs[4] = { D100, D010, D001, D000 };
int SpacesNumbers_MatrixARhs[4] = { 0, 0, 0, 0 };
int N_Matrices_MatrixARhs = 1;
int RowSpace_MatrixARhs[1] = { 0 };
int ColumnSpace_MatrixARhs[1] = { 0 };
int N_Rhs_MatrixARhs = 1;
int RhsSpace_MatrixARhs[1] = { 0 };

void MatrixARhsAssemble(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);

void MatrixAUpwindRhsAssemble(double Mult, double *coeff, double *param,
                              double hK, 
                              double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the rhs
// ======================================================================
int N_Terms_Rhs = 1;
MultiIndex3D Derivatives_Rhs[1] = { D000 };
int SpacesNumbers_Rhs[1] = { 0 };
int N_Matrices_Rhs = 0;
int *RowSpace_Rhs = NULL;
int *ColumnSpace_Rhs = NULL;
int N_Rhs_Rhs = 1;
int RhsSpace_Rhs[1] = { 0 };

void RhsAssemble(double Mult, double *coeff, double *param,
                 double hK, 
                 double **OrigValues, int *N_BaseFuncts,
                 double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the rhs
// ======================================================================
int N_Terms_Rhs_SUPG = 4;
MultiIndex3D Derivatives_Rhs_SUPG[4] = { D100, D010, D001, D000 };
int SpacesNumbers_Rhs_SUPG[4] = { 0, 0, 0, 0 };
int N_Matrices_Rhs_SUPG = 0;
int *RowSpace_Rhs_SUPG = NULL;
int *ColumnSpace_Rhs_SUPG = NULL;
int N_Rhs_Rhs_SUPG = 1;
int RhsSpace_Rhs_SUPG[1] = { 0 };

void RhsAssemble_SUPG(double Mult, double *coeff, double *param,
                      double hK, 
                      double **OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs);
// ======================================================================
// definitions for assembling the mass matrix for bulk problem
// ======================================================================
int N_Terms_MatrixM_Bulk = 1;
MultiIndex3D Derivatives_MatrixM_Bulk[1] = { D000 };
int SpacesNumbers_MatrixM_Bulk[1] = { 0 };
int N_Matrices_MatrixM_Bulk = 1;
int RowSpace_MatrixM_Bulk[1] = { 0 };
int ColumnSpace_MatrixM_Bulk[1] = { 0 };
int N_Rhs_MatrixM_Bulk = 0;
int *RhsSpace_MatrixM_Bulk = NULL;

void MatrixMAssemble_Bulk3D(double Mult, double *coeff, double *param,
			  double hK, 
			  double **OrigValues, int *N_BaseFuncts,
			  double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the matrices A for bulk problem
// ======================================================================

int N_Terms_MatricesA_SUPG_Bulk = 4;
MultiIndex3D Derivatives_MatricesA_SUPG_Bulk[4] = { D100, D010, D001, D000 };
int SpacesNumbers_MatricesA_SUPG_Bulk[4] = { 0, 0, 0, 0  };
int N_Matrices_MatricesA_SUPG_Bulk = 2;
int RowSpace_MatricesA_SUPG_Bulk[2] = { 0, 0 };
int ColumnSpace_MatricesA_SUPG_Bulk[2] = { 0, 0 };
int N_Rhs_MatricesA_SUPG_Bulk = 0;
int *RhsSpace_MatricesA_SUPG_Bulk = NULL;

void MatricesA_Assemble_SUPG_Bulk3D(double Mult, double *coeff, double *param,
                                    double hK,
                                    double **OrigValues, int *N_BaseFuncts,
                                    double ***LocMatrices, double **LocRhs);

int N_Matrices_MatricesA_Galerkin_Bulk = 1;
int RowSpace_MatricesA_Galerkin_Bulk[1] = { 0 };
int ColumnSpace_MatricesA_Galerkin_Bulk[1] = { 0 };

void MatricesA_Assemble_Bulk3D(double Mult, double *coeff, double *param,
                                  double hK,
                                  double **OrigValues, int *N_BaseFuncts,
                                  double ***LocMatrices, double **LocRhs);

void MatricesA_Assemble_Galerkin_Bulk3D(double Mult, double *coeff, double *param,
                                  double hK,
                                  double **OrigValues, int *N_BaseFuncts,
                                  double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the rhs for bulk problem
// ======================================================================

int N_Terms_Rhs_SUPG_Bulk = 4;
MultiIndex3D Derivatives_Rhs_SUPG_Bulk[4] = { D100, D010, D001, D000 };
int SpacesNumbers_Rhs_SUPG_Bulk[4] = { 0, 0, 0, 0  };
int N_Matrices_Rhs_SUPG_Bulk = 0;
int *RowSpace_Rhs_SUPG_Bulk = NULL;
int *ColumnSpace_Rhs_SUPG_Bulk = NULL;
int N_Rhs_Rhs_SUPG_Bulk = 1;
int RhsSpace_Rhs_SUPG_Bulk[1] = { 0 };

void Rhs_Assemble_SUPG_Bulk3D(double Mult, double *coeff, double *param,
                              double hK,
                              double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);

int N_Terms_Rhs_Galerkin_Bulk = 1;
MultiIndex3D Derivatives_Rhs_Galerkin_Bulk[1] = { D000 };
int SpacesNumbers_Rhs_Galerkin_Bulk[1] = { 0 };

void Rhs_Assemble_Bulk3D(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs);

// ======================================================================
// parameter routine settings
// ======================================================================
void TimeCDParamsVeloField(double *in, double *out);

int TimeCDParamsVeloFieldN_FESpaces = 1;
int TimeCDParamsVeloFieldN_Fct = 3;
int TimeCDParamsVeloFieldN_ParamFct = 1;
int TimeCDParamsVeloFieldN_FEValues = 3;
int TimeCDParamsVeloFieldN_Params = 3;
int TimeCDParamsVeloFieldFEFctIndex[3] = { 0, 1, 2 };
MultiIndex3D TimeCDParamsVeloFieldFEMultiIndex[3] = { D000, D000, D000};
ParamFct *TimeCDParamsVeloFieldFct[1] = { TimeCDParamsVeloField };
int TimeCDParamsVeloFieldBeginParam[1] = { 0 };

// ======================================================================
// parameter routine settings
// SOLD methods
// ======================================================================
void TimeCDParamsSOLD(double *in, double *out);

int TimeCDParamsSOLDN_FESpaces = 1;
int TimeCDParamsSOLDN_Fct = 2;
int TimeCDParamsSOLDN_ParamFct = 1;
int TimeCDParamsSOLDN_FEValues = 8;
int TimeCDParamsSOLDN_Params = 8;
int TimeCDParamsSOLDFEFctIndex[8] = { 0, 0, 0, 0, 1, 1, 1, 1};
MultiIndex3D TimeCDParamsSOLDFEMultiIndex[8] = { D000, D100, D010, D001,
						 D000, D100, D010, D001};

ParamFct *TimeCDParamsSOLDFct[1] = { TimeCDParamsSOLD };
int TimeCDParamsSOLDBeginParam[1] = { 0 };

// ======================================================================
// parameters for bulk problem
// ======================================================================

void TimeCDParamsBulk(double *in, double *out);

int TimeCDParamsBulkN_FESpaces = 2;
int TimeCDParamsBulkN_Fct = 4;
int TimeCDParamsBulkN_ParamFct = 1;
int TimeCDParamsBulkN_FEValues = 4;
int TimeCDParamsBulkN_Params = 4;
int TimeCDParamsBulkFEFctIndex[4] = { 0, 1, 2, 3 };
MultiIndex3D TimeCDParamsBulkFEMultiIndex[4] = { D000, D000, D000, D000 };
ParamFct *TimeCDParamsBulkFct[1] = { TimeCDParamsBulk };
int TimeCDParamsBulkBeginParam[1] = { 0 };

void TimeCDParamsBulk_Cc(double *in, double *out);

int TimeCDParamsBulk_CcN_FESpaces = 5; // velocity, C_a (C_b), integral_c_C
int TimeCDParamsBulk_CcN_Fct = 7; // u_1, u_2, u_3, C_a, C_b, C_c, integral_c_C
int TimeCDParamsBulk_CcN_ParamFct = 1; // number of ParamRout
int TimeCDParamsBulk_CcN_FEValues = 7; // u_1, u_2, u_3, C_a, C_b, C_c_old, integral_c_C
int TimeCDParamsBulk_CcN_Params = 7;
int TimeCDParamsBulk_CcFEFctIndex[7] = { 0, 1, 2, 3, 4, 5, 6 };
MultiIndex3D TimeCDParamsBulk_CcFEMultiIndex[7] = { D000, D000, D000, D000, D000, D000, D000 };
ParamFct *TimeCDParamsBulk_CcFct[1] = { TimeCDParamsBulk_Cc };
int TimeCDParamsBulk_CcBeginParam[1] = { 0 };




// ======================================================================
//
// definitions for assembling the mass matrix for urea synthesis
//
// ======================================================================

int N_Terms_MatrixM_Urea = 1;
MultiIndex3D Derivatives_MatrixM_Urea[1] = { D000 };
int SpacesNumbers_MatrixM_Urea[1] = { 0 };
int N_Matrices_MatrixM_Urea = 1;
int RowSpace_MatrixM_Urea[1] = { 0 };
int ColumnSpace_MatrixM_Urea[1] = { 0 };
int N_Rhs_MatrixM_Urea = 0;
int *RhsSpace_MatrixM_Urea = NULL;

// ======================================================================
// definitions for assembling the matrices A for urea problem
// ======================================================================

int N_Terms_MatricesA_SUPG_Urea = 4;
MultiIndex3D Derivatives_MatricesA_SUPG_Urea[4] = { D100, D010, D001, D000 };
int SpacesNumbers_MatricesA_SUPG_Urea[4] = { 0, 0, 0, 0  };
int N_Matrices_MatricesA_SUPG_Urea = 2;
int RowSpace_MatricesA_SUPG_Urea[2] = { 0, 0 };
int ColumnSpace_MatricesA_SUPG_Urea[2] = { 0, 0 };
int N_Rhs_MatricesA_SUPG_Urea = 0;
int *RhsSpace_MatricesA_SUPG_Urea = NULL;

int N_Matrices_MatricesA_Galerkin_Urea = 1;
int RowSpace_MatricesA_Galerkin_Urea[1] = { 0 };
int ColumnSpace_MatricesA_Galerkin_Urea[1] = { 0 };

// ======================================================================
// definitions for assembling the rhs for bulk problem
// ======================================================================

int N_Terms_Rhs_SUPG_Urea = 4;
MultiIndex3D Derivatives_Rhs_SUPG_Urea[4] = { D100, D010, D001, D000 };
int SpacesNumbers_Rhs_SUPG_Urea[4] = { 0, 0, 0, 0  };
int N_Matrices_Rhs_SUPG_Urea = 0;
int *RowSpace_Rhs_SUPG_Urea = NULL;
int *ColumnSpace_Rhs_SUPG_Urea = NULL;
int N_Rhs_Rhs_SUPG_Urea = 1;
int RhsSpace_Rhs_SUPG_Urea[1] = { 0 };
// assembling routine same as in BULK

int N_Terms_Rhs_Galerkin_Urea = 1;
MultiIndex3D Derivatives_Rhs_Galerkin_Urea[1] = { D000 };
int SpacesNumbers_Rhs_Galerkin_Urea[1] = { 0 };
// assembling routine same as in BULK


void TimeCDParamsUrea(double *in, double *out);

int TimeCDParamsUreaN_FESpaces = 1;
int TimeCDParamsUreaN_Fct = 3;
int TimeCDParamsUreaN_ParamFct = 1;
int TimeCDParamsUreaN_FEValues = 3;
int TimeCDParamsUreaN_Params = 3;
int TimeCDParamsUreaFEFctIndex[3] = { 0, 1, 2};
MultiIndex3D TimeCDParamsUreaFEMultiIndex[3] = { D000, D000, D000 };
ParamFct *TimeCDParamsUreaFct[1] = { TimeCDParamsUrea };
int TimeCDParamsUreaBeginParam[1] = { 0 };

void TimeCDParamsUrea_conc(double *in, double *out);

int TimeCDParamsUrea_concN_FESpaces = 4; // conc, velocity, temp, integral_conce
int TimeCDParamsUrea_concN_Fct = 6; // u_1, u_2, u_3, conc, temp, integral_conce
int TimeCDParamsUrea_concN_ParamFct = 1; // number of ParamRout
int TimeCDParamsUrea_concN_FEValues = 6; // u_1, u_2, u_3, conc, temp integral_c_C
int TimeCDParamsUrea_concN_Params = 6;
int TimeCDParamsUrea_concFEFctIndex[6] = { 0, 1, 2, 3, 4, 5};
MultiIndex3D TimeCDParamsUrea_concFEMultiIndex[6] = { D000, D000, D000, D000, D000, D000};
ParamFct *TimeCDParamsUrea_concFct[1] = { TimeCDParamsUrea_conc };
int TimeCDParamsUrea_concBeginParam[1] = { 0 };

void TimeCDParamsUrea_conc_mat(double *in, double *out);

int TimeCDParamsUrea_conc_matN_FESpaces = 1; // velocity
int TimeCDParamsUrea_conc_matN_Fct = 3; // u_1, u_2, u_3
int TimeCDParamsUrea_conc_matN_ParamFct = 1; // number of ParamRout
int TimeCDParamsUrea_conc_matN_FEValues = 3; // u_1, u_2, u_3
int TimeCDParamsUrea_conc_matN_Params = 3;
int TimeCDParamsUrea_conc_matFEFctIndex[3] = { 0, 1, 2};
MultiIndex3D TimeCDParamsUrea_conc_matFEMultiIndex[3] = { D000, D000, D000};
ParamFct *TimeCDParamsUrea_conc_matFct[1] = { TimeCDParamsUrea_conc_mat };
int TimeCDParamsUrea_conc_matBeginParam[1] = { 0 };

#endif
