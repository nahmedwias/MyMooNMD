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
 * Control parameters which prelimiter to use (flux limiter before application
 * of Zalesak's).
 * NONE: do not use a prelimiter
 * MIN_MOD: Min_Mod prelimiter.
 * GRAD_DIRECTION: setting anit-diffusive fluxes in gradient direction to 0
 * BOTH: apply both - first MIN_MOD, then GRAD_DIRECTION.
 */
enum class Prelimiter
{
    NONE, MIN_MOD, GRAD_DIRECTION, BOTH
};

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
 * @param[in,out] system_matrix The matrix to undergo algebraic flux correction.
 *                         Must be square.
 * @param[in] sol The current solution vector.
 * @param[out] rhs The current right hand side vector. Gets modified.
 * @param[in] neum_to_diri array which contains the indices of actual
 *  Dirichlet dof which were but treated as Neumann dof. CB: That paramerter
 *  is unused at the moment!
 * @param[in] continuous_proposal If true, a slightly different version
 * of the algorithm is used (continuous proposal in appllication of nodal
 * correction factors).
 * @param[in] nonsymmetric_application If true, a slightly different version
 * of the algorithm is used (nonsymmetric application of flux correction)
 */
void fem_tvd_algorithm(
    TMatrix& system_matrix,
    const std::vector<double>& sol,
    std::vector<double>& rhs,
    const std::vector<int>& neum_to_diri,
    bool continuous_proposal = false,
    bool nonsymmetric_application = false);

/**
 * @brief Apply the linear Crank-Nicolson FEM-FCT algorithm for
 * time-dependent convection-diffusion-reaction problems to a CDR system.
 *
 * The algorithm is the FEM-FCT variant which is described in Chapter 7 of
 * D.Kuzmin(2009): Explicit and implicit FEM-FCT algorithms with flux linearization.
 * J Comput. Phys. 228(7) pp. 2517--2534.
 *
 *
 * @param[in] mass Consistent mass matrix.
 * @param[in, out] stiff High order transport matrix (
 * before modifications due to time stepping!). Gets modified to be the
 * system matrix of the system which has to be solved finally (so will include
 * modifications due to flux correction AND to timestepping). Must have been
 * assembled with the INTERNAL_FULL_MATRIX_STRUCTURE switch on (but that will
 * not be checked by the algorithm).
 * @param oldsol[in] The solution vector of the last time step.
 * @param rhs[in] The current right hand side (before modifications due to
 * timestepping!). Will be modified as to be the the right hand side of the system
 * which has to be solved finally (so it will include
 * modifications due to flux correction AND to timestepping).
 * @param[in, out] rhs_old The right hand side of the last time step
 * (without modification). As output, stores the right hand side of the current
 * time step without modifications.
 * @param[in] delta_t The current time step length.
 * @param[in] neum_to_diri array which contains the indices of actual
 * Dirichlet dof which were but treated as Neumann dof due to
 * INTERNAL_FULL_MATRIX_STRUCTURE.
 * @param[in] prelim Which prelimiter to use, if at all.
 */
void crank_nicolson_fct(
    const TMatrix& mass, TMatrix& stiff,
    const std::vector<double>& oldsol,
    std::vector<double>& rhs, std::vector<double>& rhs_old,
    double delta_t,
    const std::vector<int>& neum_to_diri,
    AlgebraicFluxCorrection::Prelimiter prelim = AlgebraicFluxCorrection::Prelimiter::NONE
    );

/**
 * Sets all Dirichlet rows of the input Matrix to 0 and 1 on diagonal.
 * In the assembling of a matrix which is supposed to be treated with
 * AFC all dofs are treated as active. So after the AFC completed, the
 * actual dirichlet rows have to be reset. This is achieved here.
 *
 * @param[in,out] MatrixA The matrix to be corrected. Must be square.
 */
void correct_dirichlet_rows(FEMatrix& MatrixA);

//////////////// Consider all following methods 'private' ///////////////////////

/**
 * Compute the entire system matrix of a FEM-FCT corrected system, after
 * its mass matrix got lumped and its stiffness matrix flux corrected.
 *
 * Also performs checking of the cfl condition in case an explicit
 * method has been used as predictor ("intermediate solution").
 *
 * @note Makes sense only for (at least partly) implicit schemes.
 * @param[in, out] system_matrix The systems stiffness matrix A, with additional
 * diffusion already added to it. Must be square. Will be turned to
 * S = M_(\text{lumped}) + \Delta_t * \theta_1 * S.
 * @param[in] lumped_mass_matrix The row wise lumped mass matrix of the system,
 * as a vector of its diagonal entries. Must have length equal to the stiffness
 * matrix' dimension.
 * @param[in] tau The current time step length.
 * @param theta_1 The impliciteness parameter. Must be 0.5 so far
 * (and defaults to that value).
 *  * @param theta_2 The expliciteness parameter. Must be 0.5 so far
 * (and defaults to that value).
 */
void fem_fct_compute_system_matrix(
    TMatrix& system_matrix,
    const std::vector<double>& lumped_mass_matrix,
    double tau, double theta_1 = 0.5, double theta_2 = 0.5);

/*!
 * Compute the artificial diffusion matrix to a given stiffness matrix.
 *
 * @param[in,out] A The stiffness matrix which needs additional
 * diffusion. Must be square.
 * @param[out] matrix_D A vector to write the entries
 * of D, the artificial diffusion matrix to. (Consider TMatrix
 * with same TStructure as A!)
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
 */
void lump_mass_matrix_to_vector(const TMatrix& M, std::vector<double>& lump_mass);

/** Check whether the setup fulfils the cfl-like condition. Write a warning if not so.
 * tau is timestep, theta2 is expliciteness parameter (Kuzmin 2009, eq (11) )
 */
void check_cfl_condition(double mass_diag_entry,
                         double system_diag_entry,
                         double tau, double theta2);

/**
 * @brief The MinMod Flux Prelimiter.
 * @param a first value
 * @param b second value
 * @return The result of the MinMod Prelimiter.
 */
double MinMod(double a, double b);

/**
 * Apply Zalesaks flux limiter as described in Kuzmin 2009 pp.6-7
 * to compute flux correction factors. Used with FEM-FCT for time-
 * dependent CDR problems.
 *
 * @param[out] alphas The flux correction factors, row by row (sorted like
 * the entries array of the mass_matrix!).
 * @param[in] mass_matrix The complete mass matrix.
 * @param[in] lumped_mass The row-wise lumped mass matrix as vector
 * of its diagonal values.
 * @param[in] raw_fluxes The raw fluxes, sorted like alphas. Must have been
 * multiplied with current timesteplength!
 * @param[in] sol_approx The approximate solution used to calculate
 * the flux correction.
 * @param[in] neum_to_diri See FEM-FCT, same story.
 * @note I'm pretty sure we can get rid of the whole "neum_to_diri"
 * thing when passing an FEMatrix here which got assembled with Neumann
 * treatment of dirichlet rows.
 */
void ZalesaksFluxLimiter(
    std::vector<double>& alphas,
    const TMatrix& mass_matrix,
    const std::vector<double>& lumped_mass,
    const std::vector<double>& raw_fluxes,
    const std::vector<double>& sol_approx,
    const std::vector<int>& neum_to_diri
    );

}




#endif /* ALGEBRAIC_FLUX_CORRECTION_H_ */
