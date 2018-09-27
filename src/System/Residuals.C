/**
 * @file Implementation of struct Residuals, declared in src/system/Residuals.C
 */
#include <Residuals.h>
#include <cmath>
#include <iomanip>
#ifdef _MPI
#include <mpi.h>
#endif

///@brief standard constructor, initialize with large numbers (1e10).
Residuals::Residuals()
: impulsResidual(1e10), massResidual(1e10), fullResidual(1e10)
  {}

/// @brief constructor given the \e square of the impuls and mass
/// residuals - in MPI case, both must be stored consistently
Residuals::Residuals(double imR, double maR)
: impulsResidual(0), massResidual(0),
  fullResidual(0)
{
  impulsResidual = sqrt(imR);
  massResidual = sqrt(maR);
  fullResidual = sqrt(imR + maR);
}

/// @brief Write out the three numbers to a stream.
std::ostream& operator<<(std::ostream& s, const Residuals& n)
{
  s << std::setw(14) << n.impulsResidual << "\t" << std::setw(14)
        << n.massResidual << "\t" << std::setw(14) << n.fullResidual;
  return s;
}
