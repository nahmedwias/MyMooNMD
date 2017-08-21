// ======================================================================
// TCD3D.h
//
// common declaration for time dependent convection diffusion problems
// ======================================================================

#ifndef __TIMECONVDIFF3D__
#define __TIMECONVDIFF3D__

// MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };

// ======================================================================
// definitions for assembling the mass matrix, stiffness matrix and rhs 
// Galerkin 
// ======================================================================
void MatrixMARhsAssemble(double Mult, double *coeff, double *param,
                           double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);

void MatrixARhsAssemble(double Mult, double *coeff, double *param,
                           double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);
// ======================================================================
// definitions for assembling the mass matrix, stiffness matrix and rhs 
// SUPG: NOTE: Mass matrix also include the term from time 
// discretization
// ======================================================================
void MatricesMARhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);

// ======================================================================
// ======================================================================
void MatricesMAKRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);
// ======================================================================
// definitions for assembling the matrices A, K and rhs
// ======================================================================
// TODO: the local routine is still here which is for SOLD method which
// also includes the part from the SUPG method, One have to addopt the routine
// to be fitted in the class Time_CD3D and LocalAssembling3D
void MatricesAKSRhsAssemble_SUPG(double Mult, double *coeff, double *param,
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
void MatrixAUpwindRhsAssemble(double Mult, double *coeff, double *param,
                              double hK, 
                              double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the rhs
// ======================================================================
void RhsAssemble(double Mult, double *coeff, double *param,
                 double hK, 
                 double **OrigValues, int *N_BaseFuncts,
                 double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the rhs
// ======================================================================
void RhsAssemble_SUPG(double Mult, double *coeff, double *param,
                      double hK, 
                      double **OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs);

// ======================================================================
//  definitions for assembling the matrix M, A group fem-fct
// ======================================================================
void MatrixMAGroupFEMAssemble(double Mult, double *coeff, double *param,
                              double hK, 
                              double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);
        
// ======================================================================
// definitions for assembling the matrix C1, C2, C3, R for group fem-fct
// ======================================================================
void MatrixGroupFEMAssemble(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);      

// ======================================================================
// definitions for assembling the mass matrix for bulk problem
// ======================================================================
void MatrixMAssemble_Bulk3D(double Mult, double *coeff, double *param,
			  double hK, 
			  double **OrigValues, int *N_BaseFuncts,
			  double ***LocMatrices, double **LocRhs);

// ======================================================================
// definitions for assembling the matrices A for bulk problem
// ======================================================================
void MatricesA_Assemble_SUPG_Bulk3D(double Mult, double *coeff, double *param,
                                    double hK,
                                    double **OrigValues, int *N_BaseFuncts,
                                    double ***LocMatrices, double **LocRhs);

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
void Rhs_Assemble_SUPG_Bulk3D(double Mult, double *coeff, double *param,
                              double hK,
                              double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);

void Rhs_Assemble_Bulk3D(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs);

// ======================================================================
// parameter routine settings
// ======================================================================
void TimeCDParamsVeloField(double *in, double *out);

// ======================================================================
// parameter routine settings
// SOLD methods
// ======================================================================
void TimeCDParamsSOLD(double *in, double *out);

// ======================================================================
// parameter routine settings
// for explicit reaction
// ======================================================================
void TimeCDParamsSolution(double *in, double *out);

// ======================================================================
// parameters for bulk problem
// ======================================================================

void TimeCDParamsBulk(double *in, double *out);


void TimeCDParamsBulk_Cc(double *in, double *out);

// ======================================================================
//
// definitions for assembling the mass matrix for urea synthesis
//
// ======================================================================

// ======================================================================
// definitions for assembling the matrices A for urea problem
// ======================================================================
// ======================================================================
// definitions for assembling the rhs for bulk problem
// ======================================================================
void TimeCDParamsUrea(double *in, double *out);

void TimeCDParamsUrea_conc2(double *in, double *out);

void TimeCDParamsUrea_temp2(double *in, double *out);

void TimeCDParamsUrea_conc_mat(double *in, double *out);

#endif // __TIMECONVDIFF3D__
