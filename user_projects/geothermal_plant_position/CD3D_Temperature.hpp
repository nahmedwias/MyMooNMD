#ifndef CD3D_TEMPERATURE_H
#define CD3D_TEMPERATURE_H
#include "Time_CD3D.h"
#include "FEVectFunct3D.h"

class CD3D_Temperature : public Time_CD3D
{
  public: 
  CD3D_Temperature(std::list<TCollection* > collections, //const TDomain& domain,
 const ParameterDatabase& param_db, 
                   const Example_TimeCD3D& example);
  
  void assemble(const TFEVectFunct3D& convection, const double * x, double nu);
  
  void reset_for_output();
};


#endif // CD3D_TEMPERATURE_H
