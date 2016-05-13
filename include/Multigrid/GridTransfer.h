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

/**
 * Restrict defect from fine function to coarse function -
 * maps into dual space. For an in-depth explanation see lecture notes
 * of Volker John.
 */
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

/**
 * Restricts the function on the finest grid to all other grids
 * successively.
 *
 * @note No checks whether the spaces form an actual hierarchy are
 * performed - we rely on "garbage in, garbage out" here.
 *
 * @param[in] space_hierarchy An ordered hierarchy of FE Spaces, finest first.
 * @param[in,out] function_entries The functions to be filled.
 * @param[in] The lengths of the functions. Must match the FESpaces number of dofs.
 *
 * TODO First function should be const.
 */
void RestrictFunctionRepeatedly(
#ifdef __2D__
  std::vector<const TFESpace2D*> space_hierarchy,
#elif __3D__
  std::vector<const TFESpace3D*> space_hierarchy,
#endif
  std::vector<double*> function_entries,
  std::vector<size_t> function_n_dofs);
}

#endif /* INCLUDE_MULTIGRID_GRIDTRANSFER_H_ */
