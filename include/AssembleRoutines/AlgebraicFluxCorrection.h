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

#include <ParameterDatabase.h>

#include <vector>

class TMatrix;
class FEMatrix;

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
 * Sets up and returns a default algebraic flux correction database,
 * which contains all control parameters necessary for an algebraic
 * flux correction. All of them are initialized with default values.
 */
ParameterDatabase default_afc_database();

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
 * D.Kuzmin(2009): Explicit and implicit FEM-FCT algorithms with flux
 * linearization.
 * J Comput. Phys. 228(7) pp. 2517--2534.
 *
 * TODO CB: In the algorithm, the explicit part of the C-N scheme is computed
 * using the transport operator K_{k} of the current time step and not K_{k-1}
 * of the previous timestep.
 * Could one run into trouble when the transport is time dependent?
 *
 * @param[in] mass Consistent mass matrix.
 * @param[in, out] stiff High order transport matrix (
 * before modifications due to time stepping!). Gets modified to be the
 * system matrix of the system which has to be solved finally (so will include
 * modifications due to flux correction AND to timestepping). Must have been
 * assembled with the INTERNAL_FULL_MATRIX_STRUCTURE switch on (but that will
 * not be checked by the algorithm).
 * @param[in] oldsol The solution vector of the last time step.
 * @param[in] rhs The current right hand side (before modifications due to
 * timestepping!). Will be modified as to be the the right hand side of the =
 * system which has to be solved finally (so it will include
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
    AlgebraicFluxCorrection::Prelimiter prelim =
        AlgebraicFluxCorrection::Prelimiter::NONE
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

}




#endif /* ALGEBRAIC_FLUX_CORRECTION_H_ */
