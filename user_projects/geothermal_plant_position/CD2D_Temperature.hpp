#ifndef CD2D_TEMPERATURE_H
#define CD2D_TEMPERATURE_H
#include "Time_CD2D.h"
#include "FEVectFunct2D.h"

class CD2D_Temperature : public Time_CD2D
{
  public: 
  CD2D_Temperature(const TDomain& domain, const ParameterDatabase& param_db, 
                   const Example_TimeCD2D& example);
  
  void assemble(const TFEVectFunct2D& convection, const double * x, double nu);
  
  void reset_for_output();
};


#endif // CD2D_TEMPERATURE_H
