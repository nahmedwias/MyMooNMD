#ifndef NSE3D_GPPO_H
#define NSE3D_GPPO_H
#include "NSE3D.h"
#include "FEVectFunct3D.h"
#include "Domain.h"
#include "FEFunction3D.h"
#include "Collection.h"
#include "Example_NSE3D.h"
#include "ParameterDatabase.h"
#include "Collection.h"
#include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!
#include "LPS_scott_zhang.h"


class NSE3D_GPPO : public NSE3D
{
	public: 

 NSE3D_GPPO(TDomain &domain, const ParameterDatabase& param_db, Example_NSE3D example);


//void read_coefficient_function(std::string read_coefficient_function_directory, TCollection *coll, TFEFunction3D *coefficient_function_ptr);


//void assemble(TFEFunction3D* coefficient_function = nullptr);


};












#endif // NSE3D_GPPO_H
