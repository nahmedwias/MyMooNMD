/*
 * FEM_FCT_Copy.h
 *
 * Fem-tvd and fem-fct methods for scalar cdr equations. Copied
 * from Volkers's implementation and supposed to be adapted adapted.
 *
 * Outcommented 3D parts. Want to take advantage of Eclipse Indexing.
 *
 *  Created on: Jul 27, 2015
 *      Author: bartsch
 */

#ifndef FEM_FCT_COPY_H_
#define FEM_FCT_COPY_H_

namespace AlgebraicFluxCorrection {

//**************************************************************
//compute the lumped matrix
//output is a vector
//**************************************************************
//#ifdef __2D__
void LumpMassMatrixToVector(TSquareMatrix2D *M, double *lump_mass);
//#endif
//#ifdef __3D__
//void LumpMassMatrixToVector(TSquareMatrix3D *M, double *lump_mass);
//#endif

/*******************************************************************************/
//
// FEM_TVD_ForConvDiff for steady-state cdr equations
// following D. Kuzmin (2007)
//
/*******************************************************************************/
//#ifdef __2D__
void FEM_TVD_ForConvDiff(TSquareMatrix2D *sqmatrix, int N_U, int N_Active,
double *matrix_D_Entries, double *sol, double *rhs,
			 int N_neum_to_diri, int *neum_to_diri, int compute_matrix_D);
//#endif
//#ifdef __3D__
//void FEM_TVD_ForConvDiff(TSquareMatrix3D *sqmatrix, int N_U, int N_Active,
//double *matrix_D_Entries, double *sol, double *rhs,
//int N_neum_to_diri, int *neum_to_diri, int compute_matrix_D);
//#endif

//**************************************************************
//compute system matrix for FEM-FCT
//output is a vector
//**************************************************************
//#ifdef __2D__
void FEM_FCT_SystemMatrix(TSquareMatrix2D *M_C, TSquareMatrix2D *A,
			  double *lump_mass,int N_U);
//#endif
//#ifdef __3D__
//void FEM_FCT_SystemMatrix(TSquareMatrix3D *M_C, TSquareMatrix3D *A,
//			  double *lump_mass,int N_U);
//#endif

//**************************************************************
// MINMOD prelimiter
//**************************************************************
double MinMod(double a, double b);

/*******************************************************************************/
//
// FCT-FEM algorithm
// following D. Kuzmin, M. M"oller (2005) (nonlinear scheme)
//           D. Kuzmin (2008) (linear scheme)
//
/*******************************************************************************/
//#ifdef __2D__
void FEM_FCT_ForConvDiff(TSquareMatrix2D *M_C,TSquareMatrix2D *A,
			 int N_U, int N_Active,
			 double *lump_mass, double *matrix_D_Entries,
			 double *sol, double *oldsol,
			 double *B, double *rhs, double *rhs_old,
			 double *tilde_u,
			 int N_neum_to_diri, int *neum_to_diri,
			 int *neum_to_diri_bdry,
			 double *neum_to_diri_param,
			 int compute_matrix_D, BoundValueFunct2D *BoundaryValue,
			 double *BoundaryValues);
//#endif
//#ifdef __3D__
//void FEM_FCT_ForConvDiff(TSquareMatrix3D *M_C,TSquareMatrix3D *A,
//    int N_U, int N_Active,
//    double *lump_mass, double *matrix_D_Entries,
//    double *sol, double *oldsol,
//    double *B, double *rhs, double *rhs_old,
//    double *tilde_u,
//    int N_neum_to_diri, int *neum_to_diri,
//    double *neum_to_diri_x, double *neum_to_diri_y, double *neum_to_diri_z,
//    int compute_matrix_D, BoundValueFunct3D *BoundaryValue,
//    double *BoundaryValues);
//#endif

// This is an auxiliary function for tests etc., Containing only the simplest linear scheme
// for calculation of the auxiliary solution.
void FEM_FCT_SimpleLinear(TSquareMatrix2D *M_C,TSquareMatrix2D *A,
			 int N_U, int N_Active,
			 double *lump_mass, double *matrix_D_Entries, int compute_matrix_D,
			 double *sol, double *oldsol,
			 double *B, double *rhs, double *rhs_old,
			 int N_neum_to_diri, int *neum_to_diri,
			 double *BoundaryValues);

void computeArtificialDiffusionMatrix(const TSquareMatrix2D& A, double* matrix_D_Entries, int N_U);

}




#endif /* FEM_FCT_COPY_H_ */
