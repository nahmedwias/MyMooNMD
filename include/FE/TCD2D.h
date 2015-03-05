// ======================================================================
// TCD2D.h
//
// common declaration for time dependent convection diffusion problems
// ======================================================================

#ifndef __TIMECONVDIFF2D__
#define __TIMECONVDIFF2D__

// ======================================================================
// definitions for assembling the mass matrix and rhs
// ======================================================================
int N_Terms_MatrixMRhs = 1;
MultiIndex2D Derivatives_MatrixMRhs[1] = { D00 };
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

void MatrixMRhsAssemble_Axial3D(double Mult, double *coeff, double *param,
                           double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);


int N_Terms_MatrixMRhs_SUPG = 3;
MultiIndex2D Derivatives_MatrixMRhs_SUPG[3] = { D10, D01, D00 };
int SpacesNumbers_MatrixMRhs_SUPG[3] = { 0, 0, 0 };
int N_Matrices_MatrixMRhs_SUPG = 1;
int RowSpace_MatrixMRhs_SUPG[1] = { 0 };
int ColumnSpace_MatrixMRhs_SUPG[1] = { 0 };
int N_Rhs_MatrixMRhs_SUPG = 1;
int RhsSpace_MatrixMRhs_SUPG[1] = { 0 };

void MatrixMRhsAssemble_SUPG_Axial3D(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);

void MatrixMRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the matrices A, K and rhs
// ======================================================================
int N_Terms_MatricesAKRhs_SUPG = 5;
MultiIndex2D Derivatives_MatricesAKRhs_SUPG[5] = { D10, D01, D00, D20, D02 };
int SpacesNumbers_MatricesAKRhs_SUPG[5] = { 0, 0, 0, 0, 0  };
int N_Matrices_MatricesAKRhs_SUPG = 2;
int RowSpace_MatricesAKRhs_SUPG[2] = { 0, 0 };
int ColumnSpace_MatricesAKRhs_SUPG[2] = { 0, 0 };
int N_Rhs_MatricesAKRhs_SUPG = 1;
int RhsSpace_MatricesAKRhs_SUPG[1] = { 0 };

int N_Matrices_MatricesAKRhs_SOLD = 3;
int RowSpace_MatricesAKRhs_SOLD[3] = { 0, 0, 0 };
int ColumnSpace_MatricesAKRhs_SOLD[3] = { 0, 0, 0 };

void MatricesAKRhsAssemble_SUPG_Axial3D(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);

void MatricesAKRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the matrix A and rhs
// ======================================================================
int N_Terms_MatrixARhs = 3;
MultiIndex2D Derivatives_MatrixARhs[3] = { D10, D01, D00};
int SpacesNumbers_MatrixARhs[4] = { 0, 0, 0  };
int N_Matrices_MatrixARhs = 1;
int RowSpace_MatrixARhs[1] = { 0 };
int ColumnSpace_MatrixARhs[1] = { 0 };
int N_Rhs_MatrixARhs = 1;
int RhsSpace_MatrixARhs[1] = { 0 };

void MatrixARhsAssemble(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);

void MatrixARhsAssemble_Axial3D(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the matrices for VMM
// ======================================================================

// number of fe functions 
int N_Terms_Matrices_VMM = 4;

// derivatives needed for fe functions 
// D10 = derivative w.r.t. x
// D01 = derivative w.r.t. y
// D00 = function itself 
MultiIndex2D Derivatives_Matrices_VMM[4] = { D10, D01, D10, D01};

// fe spaces where the functions belong to
// 0 - fine space
// 1 - coarse space
int SpacesNumbers_Matrices_VMM[4] = { 0, 0, 1, 1};

// number of matrices to assemble (M, B, C)
int N_Matrices_Matrices_VMM = 3;

// fe element space which determines the number of rows in 
// the matrices 
// M - coarse space
// B - fine space
// C - coarse space
int RowSpace_Matrices_VMM[3] = { 1, 0, 1 };


// fe element space which determines the number of columns in 
// the matrices 
// M - coarse space
// B - coarse space
// C - fine space
int ColumnSpace_Matrices_VMM[3] = { 1, 1, 0 };

// number of right hand sides
int N_Rhs_Matrices_VMM = 0;

// fe space which determines the length of the rhs vector
int *RhsSpace_Matrices_VMM = NULL;

// routine to assemble the matrices and the rhs
void MatricesAssemble_VMM(double Mult, double *coeff, double *param,
                          double hK, 
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the matrices for VMM on different grids
// as in the paper by Kaya and Layton (2002)
// ======================================================================

// number of fe functions 
int N_Terms_Matrices_VMM_KL02 = 2;

// derivatives needed for fe functions 
// D10 = derivative w.r.t. x
// D01 = derivative w.r.t. y
// D00 = function itself 
MultiIndex2D Derivatives_Matrices_VMM_KL02[2] = { D10, D01};

// fe spaces where the functions belong to
// 0 - fine space
// 1 - coarse space
int SpacesNumbers_Matrices_VMM_KL02[2] = { 0, 0};

// number of matrices to assemble (M, B, C)
int N_Matrices_Matrices_VMM_KL02 = 5;

// fe element space which determines the number of rows in 
// the matrices 
// M - coarse space
// B1 - fine space
// B2 - fine space
// C1 - coarse space
// C2 - coarse space
int RowSpace_Matrices_VMM_KL02[5] = { 1, 0, 0, 1, 1};


// fe element space which determines the number of columns in 
// the matrices 
// M - coarse space
// B1 - coarse space
// B2 - coarse space
// C1 - fine space
// C2 - fine space
int ColumnSpace_Matrices_VMM_KL02[5] = { 1, 1, 1, 0, 0 };

// number of right hand sides
int N_Rhs_Matrices_VMM_KL02 = 0;

// fe space which determines the length of the rhs vector
int *RhsSpace_Matrices_VMM_KL02 = NULL;

// routine to assemble the matrices and the rhs
void MatricesAssemble_VMM_KL02(double Mult, double *coeff, double *param,
                          double hK, 
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs);


// ======================================================================
// definitions for assembling the mass matrix for bulk problem
// ======================================================================
int N_Terms_MatrixM_Bulk = 1;
MultiIndex2D Derivatives_MatrixM_Bulk[1] = { D00 };
int SpacesNumbers_MatrixM_Bulk[1] = { 0 };
int N_Matrices_MatrixM_Bulk = 1;
int RowSpace_MatrixM_Bulk[1] = { 0 };
int ColumnSpace_MatrixM_Bulk[1] = { 0 };
int N_Rhs_MatrixM_Bulk = 0;
int *RhsSpace_MatrixM_Bulk = NULL;

void MatrixMAssemble_Bulk(double Mult, double *coeff, double *param,
			  double hK, 
			  double **OrigValues, int *N_BaseFuncts,
			  double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the matrices A for bulk problem
// ======================================================================

int N_Terms_MatricesA_SUPG_Bulk = 3;
MultiIndex2D Derivatives_MatricesA_SUPG_Bulk[3] = { D10, D01, D00 };
int SpacesNumbers_MatricesA_SUPG_Bulk[3] = { 0, 0, 0  };
int N_Matrices_MatricesA_SUPG_Bulk = 2;
int RowSpace_MatricesA_SUPG_Bulk[2] = { 0, 0 };
int ColumnSpace_MatricesA_SUPG_Bulk[2] = { 0, 0 };
int N_Rhs_MatricesA_SUPG_Bulk = 0;
int *RhsSpace_MatricesA_SUPG_Bulk = NULL;

void MatricesA_Assemble_SUPG_Bulk(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs);

int N_Matrices_MatricesA_Galerkin_Bulk = 1;
int RowSpace_MatricesA_Galerkin_Bulk[1] = { 0 };
int ColumnSpace_MatricesA_Galerkin_Bulk[1] = { 0 };

void MatricesA_Assemble_Bulk(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs);

void MatricesA_Assemble_Galerkin_Bulk(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the rhs for bulk problem
// ======================================================================

int N_Terms_Rhs_SUPG_Bulk = 3;
MultiIndex2D Derivatives_Rhs_SUPG_Bulk[3] = { D10, D01, D00 };
int SpacesNumbers_Rhs_SUPG_Bulk[3] = { 0, 0, 0  };
int N_Matrices_Rhs_SUPG_Bulk = 0;
int *RowSpace_Rhs_SUPG_Bulk = NULL;
int *ColumnSpace_Rhs_SUPG_Bulk = NULL;
int N_Rhs_Rhs_SUPG_Bulk = 1;
int RhsSpace_Rhs_SUPG_Bulk[1] = { 0 };

void Rhs_Assemble_SUPG_Bulk(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs);

int N_Terms_Rhs_Galerkin_Bulk = 1;
MultiIndex2D Derivatives_Rhs_Galerkin_Bulk[1] = { D00 };
int SpacesNumbers_Rhs_Galerkin_Bulk[1] = { 0  };

void Rhs_Assemble_Bulk(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs);
// ======================================================================
// definitions for assembling the matrix M, A and rhs for moving mesh
// ======================================================================
int N_Terms_MatrixMARhs = 3;
MultiIndex2D Derivatives_MatrixMARhs[3] = { D10, D01, D00};
int SpacesNumbers_MatrixMARhs[3] = { 0, 0, 0  };
int N_Matrices_MatrixMARhs = 2;
int RowSpace_MatrixMARhs[2] = { 0, 0 };
int ColumnSpace_MatrixMARhs[2] = { 0, 0};
int N_Rhs_MatrixMARhs = 1;
int RhsSpace_MatrixMARhs[1] = { 0 };

void MatrixMARhsAssemble_Axial3D(double Mult, double *coeff, double *param,
                                 double hK, 
                                 double **OrigValues, int *N_BaseFuncts,
                                 double ***LocMatrices, double **LocRhs);

void MatrixMARhsAssemble(double Mult, double *coeff, double *param,
                         double hK, 
                         double **OrigValues, int *N_BaseFuncts,
                         double ***LocMatrices, double **LocRhs);

int N_Terms_MatrixMARhs_SUPG = 5;
MultiIndex2D Derivatives_MatrixMARhs_SUPG[5] = { D10, D01, D00, D20, D02 };
int SpacesNumbers_MatrixMARhs_SUPG[5] = { 0, 0, 0, 0, 0  };
int N_Matrices_MatrixMARhs_SUPG = 3;
int RowSpace_MatrixMARhs_SUPG[3] = { 0, 0, 0 };
int ColumnSpace_MatrixMARhs_SUPG[3] = { 0, 0, 0};
int N_Rhs_MatrixMARhs_SUPG = 1;
int RhsSpace_MatrixMARhs_SUPG[1] = { 0 };


void MatrixMARhsAssemble_SUPG(double Mult, double *coeff, double *param,
                              double hK, 
                              double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);


// int N_Terms_MatricesAKRhs_SUPG = 5;
// MultiIndex2D Derivatives_MatricesAKRhs_SUPG[5] = { D10, D01, D00, D20, D02 };
// int SpacesNumbers_MatricesAKRhs_SUPG[5] = { 0, 0, 0, 0, 0  };
// int N_Matrices_MatricesAKRhs_SUPG = 2;
// int RowSpace_MatricesAKRhs_SUPG[2] = { 0, 0 };
// int ColumnSpace_MatricesAKRhs_SUPG[2] = { 0, 0 };
// int N_Rhs_MatricesAKRhs_SUPG = 1;
// int RhsSpace_MatricesAKRhs_SUPG[1] = { 0 };
// 
// int N_Matrices_MatricesAKRhs_SOLD = 3;
// int RowSpace_MatricesAKRhs_SOLD[3] = { 0, 0, 0 };
// int ColumnSpace_MatricesAKRhs_SOLD[3] = { 0, 0, 0 };
// 
// void MatricesAKRhsAssemble_SUPG(double Mult, double *coeff, double *param,
//                             double hK, 
//                             double **OrigValues, int *N_BaseFuncts,
//                             double ***LocMatrices, double **LocRhs);



// ======================================================================
// parameter routine settings
// SOLD methods
// ======================================================================
void TimeCDParamsSOLD(double *in, double *out);

int TimeCDParamsSOLDN_FESpaces = 1;
int TimeCDParamsSOLDN_Fct = 2;
int TimeCDParamsSOLDN_ParamFct = 1;
int TimeCDParamsSOLDN_FEValues = 10;
int TimeCDParamsSOLDN_Params = 10;
int TimeCDParamsSOLDFEFctIndex[10] = { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
MultiIndex2D TimeCDParamsSOLDFEMultiIndex[10] = { D00, D10, D01, D20, D02, 
						 D00, D10, D01, D20, D02};

ParamFct *TimeCDParamsSOLDFct[1] = { TimeCDParamsSOLD };
int TimeCDParamsSOLDBeginParam[1] = { 0 };

// ======================================================================
// parameters for bulk problem
// ======================================================================

void TimeCDParamsBulk(double *in, double *out);

int TimeCDParamsBulkN_FESpaces = 2;
int TimeCDParamsBulkN_Fct = 3;
int TimeCDParamsBulkN_ParamFct = 1;
int TimeCDParamsBulkN_FEValues = 3;
int TimeCDParamsBulkN_Params = 3;
int TimeCDParamsBulkFEFctIndex[3] = { 0, 1, 2 };
MultiIndex2D TimeCDParamsBulkFEMultiIndex[3] = { D00, D00, D00 };
ParamFct *TimeCDParamsBulkFct[1] = { TimeCDParamsBulk };
int TimeCDParamsBulkBeginParam[1] = { 0 };

void TimeCDParamsBulk_SOLD(double *in, double *out);

int TimeCDParamsBulk_SOLDN_FESpaces = 3;
int TimeCDParamsBulk_SOLDN_Fct = 5;
int TimeCDParamsBulk_SOLDN_ParamFct = 1;
int TimeCDParamsBulk_SOLDN_FEValues = 10;
int TimeCDParamsBulk_SOLDN_Params = 13;
int TimeCDParamsBulk_SOLDFEFctIndex[10] = { 0, 1, 2, 3, 4, 3, 3, 5, 5, 5 };
MultiIndex2D TimeCDParamsBulk_SOLDFEMultiIndex[10] = { D00, D00, D00, D00, D00, D10, D01, D00,
						       D10, D01};
ParamFct *TimeCDParamsBulk_SOLDFct[1] = { TimeCDParamsBulk_SOLD };
int TimeCDParamsBulk_SOLDBeginParam[1] = { 0 };

void TimeCDParamsBulk_Cc(double *in, double *out);

int TimeCDParamsBulk_CcN_FESpaces = 5; // C_a, C_b, velo, C_c, integral_c_C
int TimeCDParamsBulk_CcN_Fct = 6; // C_a, C_b, u_1, u_2, C_c, integral_c_C
int TimeCDParamsBulk_CcN_ParamFct = 1; // number of ParamRout
int TimeCDParamsBulk_CcN_FEValues = 6; // C_a, C_b, u_1, u_2, C_c, integral_c_C
int TimeCDParamsBulk_CcN_Params = 6;
int TimeCDParamsBulk_CcFEFctIndex[6] = { 0, 1, 2, 3, 4, 5 };
MultiIndex2D TimeCDParamsBulk_CcFEMultiIndex[6] = { D00, D00, D00, D00, D00, D00 };
ParamFct *TimeCDParamsBulk_CcFct[1] = { TimeCDParamsBulk_Cc };
int TimeCDParamsBulk_CcBeginParam[1] = { 0 };

void TimeCDParamsBulk_SOLD_Cc(double *in, double *out);

int TimeCDParamsBulk_SOLD_CcN_FESpaces = 5; // C_a, C_b, velo, C_c, integral_c_C
int TimeCDParamsBulk_SOLD_CcN_Fct = 9; // C_a, C_b, u_1, u_2, C_c, integral_c_C
                                       // C_c_old, C_a_old, C_b_old
int TimeCDParamsBulk_SOLD_CcN_ParamFct = 1; // number of ParamRout
int TimeCDParamsBulk_SOLD_CcN_FEValues = 13; // C_a, C_b, u_1, u_2, C_c, integral_c_C
int TimeCDParamsBulk_SOLD_CcN_Params = 13;
int TimeCDParamsBulk_SOLD_CcFEFctIndex[13] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 6, 6, 4, 4  };
MultiIndex2D TimeCDParamsBulk_SOLD_CcFEMultiIndex[13] = { D00, D00, D00, D00, D00, D00,
							 D00, D00, D00, D10, D01, D10, D01};
ParamFct *TimeCDParamsBulk_SOLD_CcFct[1] = { TimeCDParamsBulk_SOLD_Cc };
int TimeCDParamsBulk_SOLD_CcBeginParam[1] = { 0 };

void TimeCDParamsBulk_mom(double *in, double *out);

int TimeCDParamsBulk_momN_FESpaces = 3; // vel, c_C, mom_{k-1} 
int TimeCDParamsBulk_momN_Fct = 4; // u_1, u_2,  c_C, mom_{k-1} 
int TimeCDParamsBulk_momN_ParamFct = 1; // number of ParamRout
int TimeCDParamsBulk_momN_FEValues = 4; // u_1, u_2,  c_C, mom_{k-1} 
int TimeCDParamsBulk_momN_Params = 4;
int TimeCDParamsBulk_momFEFctIndex[4] = { 0, 1, 2, 3};
MultiIndex2D TimeCDParamsBulk_momFEMultiIndex[4] = { D00, D00, D00, D00};
ParamFct *TimeCDParamsBulk_momFct[1] = { TimeCDParamsBulk_mom };
int TimeCDParamsBulk_momBeginParam[1] = { 0 };

void TimeCDParamsBulk_SOLD_mom(double *in, double *out);

int TimeCDParamsBulk_SOLD_momN_FESpaces = 3; // velo, C_c
int TimeCDParamsBulk_SOLD_momN_Fct = 4; // u_1, u_2, C_c, C_c_old
int TimeCDParamsBulk_SOLD_momN_ParamFct = 1; // number of ParamRout
int TimeCDParamsBulk_SOLD_momN_FEValues = 4; //  u_1, u_2, C_c, C_c_old
int TimeCDParamsBulk_SOLD_momN_Params = 10;
int TimeCDParamsBulk_SOLD_momFEFctIndex[4] = { 0, 1, 2, 3 };
MultiIndex2D TimeCDParamsBulk_SOLD_momFEMultiIndex[4] = { D00, D00, D00, D00};
ParamFct *TimeCDParamsBulk_SOLD_momFct[1] = { TimeCDParamsBulk_SOLD_mom };
int TimeCDParamsBulk_SOLD_momBeginParam[1] = { 0 };

#
#endif
