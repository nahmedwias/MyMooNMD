#ifndef VOLUMEOFFLUID_H
#define VOLUMEOFFLUID_H

#include "vector"
#include "../include/General/templateNames.h"

#ifdef __2D__
#include "../include/FE/FEFunction2D.h"
#include "../include/FE/FEVectFunct2D.h"
#else
#include "../include/FE/FEFunction3D.h"
#include "../include/FE/FEVectFunct3D.h"
#endif
/**
 * @todo write docs
 */

template <int d>
class VolumeOfFluid
{
public:
  using FEFunction = typename Template_names<d>::FEFunction;
  using FEVectFunct = typename Template_names<d>::FEVectFunct;
  using FESpace = typename Template_names<d>::FESpace;
  using Example_TimeNSE = typename Template_names<d>::Example_TimeNSE;
  
protected:
  VolumeOfFluid();
};

#endif // VOLUMEOFFLUID_H
