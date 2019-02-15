#ifndef LOCALASSEMBLECOUPLEDSTRESSNS_H
#define LOCALASSEMBLECOUPLEDSTRESSNS_H

#include "LocalAssembling.h"

// temprary for 2D, TODO: generalize the class
class LocalAssembleCoupledStressNS : public LocalAssembling<2>
{
public:
  // constructor 
  LocalAssembleCoupledStressNS(ParameterDatabase db, CoeffFct2D coeffs);
};

#endif // LOCALASSEMBLECOUPLEDSTRESSNS_H