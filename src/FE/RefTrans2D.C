#include "RefTrans2D.h"
#include "MooNMD_Io.h"

void TRefTrans2D::PiolaMapOrigFromRef(double, double, int, const double *,
                                      double *)
{
  ErrThrow("Piola Map not defined for this element");
}

void TRefTrans2D::PiolaMapOrigFromRef(double, double, int, const double *,
                                      const double *, const double *,
                                      double *, double *)
{
  ErrThrow("Piola Map not defined for this element");
}
