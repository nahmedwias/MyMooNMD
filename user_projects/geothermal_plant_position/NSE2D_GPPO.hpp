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
	public: 


NSE2D_GPPO(const TDomain &domain, const ParameterDatabase& param_db, const Example_NSE2D& example);


void read_coefficient_function(const ParameterDatabase& param_db, TCollection *coll, TFEFunction2D *coefficient_function_ptr);

void assemble(TFEFunction2D* coefficient_function);

};


#endif // NSE2D_GPPO_H
