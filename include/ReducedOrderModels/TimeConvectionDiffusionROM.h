#ifndef TIMECONVECTIONDIFFUSIONROM_H
#define TIMECONVECTIONDIFFUSIONROM_H

// boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "ROM.h"
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include "deque"
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
#include "LocalAssembling.h"

template<int d>
class TimeConvectionDiffusionROM : public ROM
{
public:
  using FEFunction = typename Template_names<d>::FEFunction;
  using FESpace = typename Template_names<d>::FESpace;
  using Example_TimeCD = typename Template_names<d>::Example_TimeCD;
  /** @brief The standard constructor, can be used for multigrid and 
   * non-multigrid.
   *
   * @param[in] Domain collection of cells
   * @param[in] param_db A parameter database with parameters concerning this
   *                     class or any of its members (fe space, solver,
   *                     assemble,...)
   * @param[in] example The example which is to be calculated.
   */
  TimeConvectionDiffusionROM(const ParameterDatabase& param_db,
			     const Example_TimeCD& ex);
  
  /** @brief return a database with all parameters necessary for 
   * time-dependent convection-diffusion (tcd) probems
   */
  static ParameterDatabase default_tcd_rom_database();
  
  /** @brief Compute ROM initial condition
   * 
   * With parameter 'rom_init_regularized' set to false, ROM initial condition
   * is computed in the standard way (L2 projection into POD space). Otherwise,
   * the regularized ROM initial condition (see PhD Thesis of S.Giere, Sec. 3.2.2.) is
   * computed. The filter width in the Helmholtz equation is controlled by database
   * parameter 'differential_filter_width'.
   */
  void compute_initial_solution();
  
  /** @brief assemble matrices
   * stiffness matrix and mass matrix for the CDR problem
   * also
   */
  void assemble_matrices_rhs(bool mat_depend_time = false);
  
  /** @brief Compute system matrix
   * 
   * System matrix sys_mat is computed and its inverse is pre-computed for a
   * more efficient solution of the linear system (only valid for
   * time-independent problem coefficients).
   */
  void set_system_matrix();
  
  /** @brief Compute system rhs
   * 
   * @note Note that assembling in every time step in very inefficient in the
   * ROM context. Try to avoid it, e.g., by storing all FE source terms at all
   * times and reduce them offline and online just use the reduced vectors.
   * Alternatively, for source terms given in the separated time-space form,
   * one could pre-assemble the space part and reduce it offline and within the
   * time loop just multiply it with the corresponding temporal coefficient
   * (e.g., see PhD Thesis of S.Giere, Sec. 4.2.).
   *
   * @param reassemble Flag indicating if the source term has to be reassembled
   */
  void set_system_rhs(bool reassemble=true);
  /** @brief solve the system */
  void solve();
  
  /** @brief measure errors and write solution
   * 
   */
  void output(int m);
  
  
  // getters and setters
  TimeDiscretization& get_time_stepping_scheme()
  { return time_stepping_scheme; }
  const TimeDiscretization& get_time_stepping_scheme() const
  { return time_stepping_scheme; }
  
private:
  /** @brief store a complete system on a particular grid.
   * 
   * This combines a matrix, rhs, solution, spaces and functions 
   * needed to describe a Time CDR problem in 2D
   */
   struct System_per_grid
   {
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
     /** @brief constructor*/
     System_per_grid(const Example_TimeCD& example, TCollection& coll, int order);
     
     /**
      * Special member functions mostly deleted,
      * for struct takes ownership of the bad
      * classes TFEFunction2D and TFESpace2D.
      */
     //! Delete copy constructor.
     System_per_grid(const System_per_grid&) = delete;
     
     //! Delete move constructor.
     System_per_grid(System_per_grid&&) = delete;
     
     //! Delete copy assignment operator.
     System_per_grid& operator=(const System_per_grid&) = delete;
     
     //! Delete move assignment operator.
     System_per_grid& operator=(System_per_grid&&) = delete;
     
     //! Default destructor. Most likely causes memory leaks.
     ~System_per_grid() = default;
   };
    
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
  
  /** @brief a complete system on each grid 
    *
    * Note that the size of this deque is at least one and larger than that
    * only in case of multigrid (when it holds as many systems as there are
    * multigrid levels).
    */
  std::deque<System_per_grid> systems;
  
  /** @brief Definition of the used example */
  const Example_TimeCD example;
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
  
  /** @brief set parameters in database
   * 
   * This functions checks if the parameters in the database are meaningful 
   * and resets them otherwise. The hope is that after calling this function
   * this class is fully functional. 
   * 
   * If some parameters are set to unsupported values, an error occurs and 
   * throws an exception.
   */
  void set_parameters();
  /**
   * This wraps up the tedious call to Assemble2D and is only used
   * to avoid code duping. Is really not written very sophisticated,
   * use it with care.
   * @param block_mat should be one system's stiffness or mass matrix.
   * @param la A fittingly constructed LocalAssemble2D object which
   * is responsible for the assembling of the stiffness matrix and right hand
   * side.
   * @param assemble_both If true, both stiffness (+rhs) and mass matrix are
   * assembled, if false only stiffness matrix and rhs.
   */
  void assemble_and_reduce(LocalAssembling_type type);
  
  // temporry matrix
  BlockFEMatrix temp_matrix;
  
};

#endif // TIMECONVECTIONDIFFUSIONROM_H
