#ifndef TCD_TEMPERATURE_H
#define TCD_TEMPERATURE_H

#include "TimeConvectionDiffusion.h"



#ifdef __2D__
#include "Example_TimeCD2D.h"
#include "FEFunction2D.h"
#else
#include "Example_TimeCD3D.h"
#include "FEFunction3D.h"
#endif

#include "Domain.h"
//class ParameterDatabase;

/** ************************************************************************ */
template<int d>
class TCD_Temperature : public TimeConvectionDiffusion<d>
{
public:

  using FEFunction = typename Template_names<d>::FEFunction;
  using FEVectFunct = typename Template_names<d>::FEVectFunct;
  using FESpace = typename Template_names<d>::FESpace;
  using Example_TimeCD = typename Template_names<d>::Example_TimeCD;
  using BoundaryValuesFunction = typename Template_names<d>::BoundaryValuesFunction;

  //constexpr static char required_database_name[] = "TCD parameter database";



#ifdef _MPI
  TCD_Temperature(TDomain& domain,
          const ParameterDatabase& param_db,
          Example_TimeCD example,
          int maxSubDomainPerDof);
#else
  TCD_Temperature(const TDomain& domain,
          const ParameterDatabase& param_db,
          Example_TimeCD example);
#endif

  void assemble(const FEVectFunct& convection, const double * x, double nu);

  void reset_for_output();
};


#endif // TCD_TEMPERATURE_H
