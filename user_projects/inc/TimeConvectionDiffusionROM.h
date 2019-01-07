#ifndef TIMECONVECTIONDIFFUSIONROM_H
#define TIMECONVECTIONDIFFUSIONROM_H


#include "ROM.h"
#include <BlockFEMatrix.h>
#include <BlockVector.h>

#ifdef __2D__
#include "Example_TimeCD2D.h"
#include "FEFunction2D.h"
#else
#include "Example_TimeCD3D.h"
#include "FEFunction3D.h"
#endif

#include "DataWriter.h"
#include "templateNames.h"
#include "TimeDiscretizations.h"
#include "Solver.h"

template <int d>
class TimeConvectionDiffusionROM : public ROM
{
public:
  using FEFunction = typename Template_names<d>::FEFunction;
  using FESpace = typename Template_names<d>::FESpace;
  using Example_TimeCD = typename Template_names<d>::Example_TimeCD;
  
  TimeConvectionDiffusionROM();
  
private:
  /** @brief Finite Element space */
  FESpace space;
  /** @brief Gramian matrix (needed for reduction of solution) */
  BlockFEMatrix gramian_matrix;
  /** @brief right hand side vector */
  BlockVector rhs;
  /** @brief full-order solution */
  BlockVector solution;
  /** @brief Finite element function */
  FEFunction fe_function;
  
  /** @brief a local parameter database which controls this class
   * 
   * The database given to the constructor will be merged into this one. Only
   * parameters which are of interest to this class are stored (and the
   * default ParMooN parameters). Note that this usually does not include
   * other parameters such as solver parameters. Those are only in the
   * Solver object.
   */
  ParameterDatabase db;
  
  /** @brief a solver object which will solve the linear system
   * 
   * Storing it means that for a direct solver we also store the factorization
   * which is usually not necessary.
   */
  Solver<BlockFEMatrix, BlockVector> solver;
  
  /** @brief Definition of the used example */
  const Example_TimeCD2D example;
  /** @brief Reduced system matrix */
  ublas::matrix<double> sys_mat_;
  /** @brief Reduced system rhs */
  ublas::vector<double> sys_rhs_;
  /** @brief ROM solution */
  ublas::vector<double> red_sol_;
  
  /** @brief Reduced mass matrix */
  ublas::matrix<double> mass_mat_;
  /** @brief Reduced convection-diffusion-reaction matrix */
  ublas::matrix<double> cdr_mat_;
  /** @brief Reduced convection-diffusion-reaction matrix times snaps mean */
  ublas::vector<double> cdr_mat_mean_;
  /** @brief Reduced source term (of its variational formulation)*/
  ublas::vector<double> reduced_rhs_;
  
  /// @brief time stepping scheme object to access everything
  TimeDiscretization time_stepping_scheme;
  
  /** @brief store the errors to compute accumulated error norms */
  std::vector<double> errors;
  
  /// @brief class for handling (time dependent) output 
  DataWriter<d> outputWriter;
  
};

#endif // TIMECONVECTIONDIFFUSIONROM_H