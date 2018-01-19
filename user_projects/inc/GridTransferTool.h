/*
 * GridTransferTool.h
 *
 *  Created on: May 19, 2017
 *      Author: bartsch
 */

#ifndef USER_PROJECTS_INC_GRIDTRANSFERTOOL_H_
#define USER_PROJECTS_INC_GRIDTRANSFERTOOL_H_

#ifdef __XD__
#define TFEFunctionXD TFEFunctionXD
#define TFESpaceXD TFESpaceXD
#define GetFESpaceXD() GetFESpaceXD()
#elif defined(__3D__)
#define TFEFunctionXD TFEFunction3D
#define TFESpaceXD TFESpace3D
#define GetFESpaceXD() GetFESpace3D()
#endif

#include <memory>
#include <valarray>

class FEFunctionInterpolator;
class TCollection;
class TFEFunctionXD;
class TFESpaceXD;

enum class GridTransferType{ MultiGrid, Interpolation};

/**
 * A tool which wraps up several different approaches to transfer fe functions
 * between different grids. This calls for a derivation hierarchy, but this is
 * too much work at this early stage. Instead, two different types of transfer
 * are crammed into one class.
 *
 * If the type MultiGrid is used, the used grids must be in parent-child relation,
 * with a maximum distance of 1 refinement level.
 * No check of this is performed and it holds: garbage in, garbage out.
 * Sorry for that.
 */
class GridTransferTool
{
  public:
    GridTransferTool(
        const TFESpaceXD* to_space,
        GridTransferType type);

    void transfer(const TFEFunctionXD& input_fct, TFEFunctionXD& output_fct,
                  std::vector<double>& output_fct_values) const;

  private:

    /// The grid to which this tool maps. Should,
    /// of course, be a shared or weak pointer instead.
    const TFESpaceXD* to_space_;

    /// The transfer type.
    GridTransferType type_;

    // A vector of function interpolators.
    mutable std::vector<std::shared_ptr<FEFunctionInterpolator>> interpolators_;


};



#endif /* USER_PROJECTS_INC_GRIDTRANSFERTOOL_H_ */
