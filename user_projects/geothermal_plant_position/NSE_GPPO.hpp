#ifndef NSE_GPPO_H
#define NSE_GPPO_H

#include "NavierStokes.h"

#ifdef __2D__
#include "FEVectFunct2D.h"
#include "Example_NSE2D.h"
#include "FEFunction2D.h"
#else
#include "FEVectFunct3D.h"
#include "Example_NSE3D.h"
#include "FEFunction3D.h"
#endif

#include "Domain.h"
#include "ParameterDatabase.h"
#include "Collection.h"
#include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!
#include "LPS_scott_zhang.h"

#include "templateNames.h"

/** ************************************************************************ */
template<int d>
class NSE_GPPO : public NavierStokes<d>
{
protected:

#ifdef __2D__
  std::shared_ptr<const TFESpace2D> coefficient_function_FEspace;
  BlockVector coefficient_function_vector;
  TFEFunction2D coefficient_function;
#endif


public:

  using FEFunction = typename Template_names<d>::FEFunction;
  using FEVectFunct = typename Template_names<d>::FEVectFunct;
  using DoubleFunction = typename Template_names<d>::DoubleFunction;
  using BoundaryValuesFunction = typename Template_names<d>::BoundaryValuesFunction;
  using FESpace = typename Template_names<d>::FESpace;
  using Example_NSE = typename Template_names<d>::Example_NSE;

  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryConditionFunction  = typename Template_names<d>::BoundaryConditionFunction;

 //using LocalAssembling = LocalAssembling2D
 // using Assemble = Assemble2D




#ifdef __2D__
  NSE_GPPO(const TDomain &domain, const ParameterDatabase& param_db, const Example_NSE& example);

  void assemble_with_coefficient_fct(TFEFunction2D* coefficient_function = nullptr);

  FEFunction & get_coefficient_function()
  {
    return this->coefficient_function;
  }

#else
  NSE_GPPO(TDomain &domain, const ParameterDatabase& param_db, Example_NSE example);

#endif



};

#endif // NSE_GPPO_H
