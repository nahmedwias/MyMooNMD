#include <Residuals.h>
#include <cmath>
#include <iomanip>


Residuals::Residuals()
  : impulsResidual(1e10), massResidual(1e10), fullResidual(1e10)
{}

Residuals::Residuals(double imR, double maR)
 : impulsResidual(std::sqrt(imR)), massResidual(std::sqrt(maR)),
   fullResidual(std::sqrt(imR + maR))
{
}

std::ostream& operator<<(std::ostream& s, const Residuals& n)
{
  s << std::setw(14) << n.impulsResidual << "\t" << std::setw(14)
        << n.massResidual << "\t" << std::setw(14) << n.fullResidual;
  return s;
}
