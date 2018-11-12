#include "RefTrans2D.h"
#include "MooNMD_Io.h"

void TRefTrans2D::PiolaMapOrigFromRef(double ci, double eta, int N_Functs,
                                      const double *refD00, double *origD00)
{
  ErrThrow("Piola Map not defined for this element");
};

void TRefTrans2D::PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                                      const double *refD00,
                                      const double *refD10,
                                      const double *refD01, double *origD10,
                                      double *origD01)
{
  ErrThrow("Piola Map not defined for this element");
};
