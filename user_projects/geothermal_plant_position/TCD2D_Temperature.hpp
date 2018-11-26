#ifndef TCD2D_TEMPERATURE_H
#define TCD2D_TEMPERATURE_H

#include "Time_CD2D.h"
#include "FEVectFunct2D.h"

//class ParameterDatabase;

/** ************************************************************************ */
class TCD2D_Temperature : public Time_CD2D
{
public:

 // constexpr static char required_database_name[] = "TCD2D parameter database";

  TCD2D_Temperature(const TDomain& domain, const ParameterDatabase& param_db,
          Example_TimeCD2D example);

  void assemble(const TFEVectFunct2D& convection, const double * x, double nu);

  void reset_for_output();
};


#endif // TCD2D_TEMPERATURE_H
