#ifndef NSE2D_GPPO_H
#define NSE2D_GPPO_H


#include "NSE2D.h"
#include "FEVectFunct2D.h"
#include "Domain.h"
#include "ParameterDatabase.h"
#include "Collection.h"
#include "Example_NSE2D.h"
#include "FEFunction2D.h"
#include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!
#include "LPS_scott_zhang.h"



class NSE2D_GPPO : public NSE2D
{

protected:

  std::shared_ptr<const TFESpace2D> coefficient_function_FEspace;
  BlockVector coefficient_function_vector;
  TFEFunction2D coefficient_function;

public:

  NSE2D_GPPO(const TDomain &domain, const ParameterDatabase& param_db, const Example_NSE2D& example);

  void assemble_with_coefficient_fct(TFEFunction2D* coefficient_function = nullptr);

  TFEFunction2D & get_coefficient_function()
  {
    return this->coefficient_function;
  }

};


#endif // NSE2D_GPPO_H
