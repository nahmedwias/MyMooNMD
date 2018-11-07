#include "RefTrans2D.h"
#include "MooNMD_Io.h"

void TRefTrans2D::PiolaMapOrigFromRef(int N_Functs, double *refD00,
                                      double *origD00)
{
  ErrThrow("Piola Map not defined for this element");
};

void TRefTrans2D::PiolaMapOrigFromRef(int N_Functs, double *refD10,
                                      double *refD01, double *origD10,
                                      double *origD01)
{
  ErrThrow("Piola Map not defined for this element");
};
