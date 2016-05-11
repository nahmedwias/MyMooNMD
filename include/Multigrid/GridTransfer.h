/**
 * @file GridTransfer.h
 *
 * Gathers grid transfer operations needed for multgrid, formerly these were
 * declared in LinAlg.h and defined in MultigridComponents2D and
 * MultigridComponents3D.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_GRIDTRANSFER_H_
#define INCLUDE_MULTIGRID_GRIDTRANSFER_H_

#include <vector>
#include <cstddef>

//forward declarations
class TFESpace2D;
class TFESpace3D;

namespace GridTransfer
{

#ifdef __2D__

void Prolongate(
    const TFESpace2D& CoarseSpace, const TFESpace2D& FineSpace,
    const double* CoarseFunction, size_t n_coarse_dofs,
    double* FineFunction, size_t n_fine_dofs);

void DefectRestriction(
    const TFESpace2D& CoarseSpace, const TFESpace2D& FineSpace,
    double* CoarseFunction, size_t n_coarse_dofs,
    const double* FineFunction, size_t n_fine_dofs);

/** function restriction from level+1 to level */
void RestrictFunction(
    const TFESpace2D& CoarseSpace, const TFESpace2D& FineSpace,
    double* CoarseFunction, size_t n_coarse_dofs,
    const double* FineFunction, size_t n_fine_dofs);

#endif
#ifdef __3D__

/** prolongate */
void Prolongate(
    const TFESpace3D& CoarseSpace, const TFESpace3D& FineSpace,
    const double* CoarseFunction, size_t n_coarse_dofs,
    double* FineFunction, size_t n_fine_dofs);

/** defect restriction from level+1 to level */
void DefectRestriction(
    const TFESpace3D& CoarseSpace, const TFESpace3D& FineSpace,
    double* CoarseFunction, size_t n_coarse_dofs,
    const double* FineFunction, size_t n_fine_dofs);

/** function restriction from level+1 to level */
void RestrictFunction(
    const TFESpace3D& CoarseSpace, const TFESpace3D& FineSpace,
    double* CoarseFunction, size_t n_coarse_dofs,
    const double* FineFunction, size_t n_fine_dofs);

#endif
}




#endif /* INCLUDE_MULTIGRID_GRIDTRANSFER_H_ */
