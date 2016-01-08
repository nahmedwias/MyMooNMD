/*
 * AlgebraicFluxCorrection.h
 *
 * Implementation of FEM-TVD and FEM-FCT methods, which are algebraic
 * flux correction algorithms for scalar cdr equations. The methods are
 * copied from the old implementation in FEM_TVD_FCT.C and currently
 * undergo refactoring.
 *
 * @date: 2015/11/11
 * @author: Volker John, Ellen Schmeyer, Clemens Bartsch
 */

#ifndef ALGEBRAIC_FLUX_CORRECTION_H_
#define ALGEBRAIC_FLUX_CORRECTION_H_

namespace AlgebraicFluxCorrection {


/**
 * FEM-TVD ('Total Variation Diminishing') algorithm for steady-state
 * convection-diffusion-reaction problems, following
 *
 * Kuzmin, D.(2007): Algebraic Flux Correction for Finite Element
 * Discretizations of Coupled Systems.
 *
 * Note that for the algorithm to operate properly, the system matrix
 * must have been constructed as if it only had active degrees of freedom.
 * This is not checked, because when reaching here the matrix is a purely
 * algebraic object without knowledge of FE spaces.
 * To regain the correct Dirichlet rows, one has to call
 * correct_dirichlet_rows upon the system matrix afterwards.
 *
 * @param[in,out] sqmatrix The matrix to undergo algebraic flux correction.
 *                         Must be square.
 * @param[in, out] matrix_D The artificial diffusion matrix, as a vector of its
 * entries. Gets completely recomputed if must_compute_D.
 * @param[in] must_compute_D flag which says if artificial diffusion matrix ("D")
 * is given or has to be (re-)computed.
 * @param[in] sol The current solution vector.
 * @param[out] rhs The current right hand side vector. Gets modified accordingly.
 * @param[in] neum_to_diri array which contains the indices of actual
 *  Dirichlet dof which were but treated as Neumann dof
 * @param[in] continuous_proposal If true, a slightly different version
 * of the algorithm is used (continuous proposal in appllication of nodal
 * correction factors).
 * @param[in] nonsymmetric_application If true, a slightly different version
 * of the algorithm is used (nonsymmetric application of flux correction)
 *
 * @note The parameters continuous_proposal and nonsymmetric_application
 * were formerly controlled by global Database entries "P8" and "P9".
 * The current setup is to be maintaiend only until we have a new
 * parameter passing framework which avoids global scope.
 */
void fem_tvd_algorithm(
    TMatrix& sqmatrix,
    std::vector<double>& matrix_D,
    bool must_compute_D,
    const std::vector<double>& sol,
    std::vector<double>& rhs,
    const std::vector<int>& neum_to_diri,
    bool continuous_proposal = false,
    bool nonsymmetric_application = false);

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

/**
 * Sets all Dirichlet rows of the input Matrix to 0 and 1 on diagonal.
 * In the assembling of a matrix which is supposed to be treated with
 * AFC all dofs are treated as active. So after the AFC completed, the
 * actual dirichlet rows have to be reset. This is achieved here.
 *
 * @param[in,out] MatrixA The matrix to be corrected. Must be square.
 */
void correct_dirichlet_rows(FEMatrix& MatrixA);

/*!
 * Compute the artificial diffusion matrix. Imagine this to be private.
 * Put here to make FEM-FCT algo shorter.
 *
 * @param[in,out] A The stiffness matrix which needs additional
 * diffusion. Must be square.
 * @param[out] matrix_D A vector to write the entries
 * of D, the artificial diffusion matrix to.
 */
void compute_artificial_diffusion_matrix(
    const TMatrix& A, std::vector<double>& matrix_D);

/**
 * @brief Lump a mass matrix row wise.
 *
 * @param[in] M The mass matrix to be lumped. Must be square.
 *
 * @param[out] lump_mass A vector containing the diagonal entries
 *  of the lumped mass matrix. Must be of according length.
 *
 * @note This method is untested in ParMooN!
 */
void lump_mass_matrix_to_vector(const TMatrix& M, std::vector<double>& lump_mass);

/**
 * @brief The MinMod Flux Prelimiter.
 * @param a first value
 * @param b second value
 * @return The result of the MinMod Prelimiter.
 */
double MinMod(double a, double b);




}




#endif /* ALGEBRAIC_FLUX_CORRECTION_H_ */
