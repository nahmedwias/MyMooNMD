#ifndef COUPLEDNAVIERSTOKESSTRESS_H
#define COUPLEDNAVIERSTOKESSTRESS_H

#include "BlockFEMatrix.h"
#ifdef __2D__
#include "Example_NSE2D.h"
#include "FEFunction2D.h"
#include "FEVectFunct2D.h"
#else
#include "Example_NSE3D.h"
#include "FEFunction3D.h"
#include "FEVectFunct3D.h"
#endif
#include "BlockVector.h"
#include "ParameterDatabase.h"
#include "Solver.h"
#include "templateNames.h"
#include "LocalAssembling.h"
#include "DataWriter.h"
#include "Residuals.h"
#include "MainUtilities.h" // FixedSizeQueue
#include "TimeDiscretizations.h"
#include "Example_CoupledNS_Stress.h"


template<int d>
class CoupledNavierStokesStress
{
public:
  using FEFunction = typename Template_names<d>::FEFunction;
  using FEVectFunct = typename Template_names<d>::FEVectFunct;
  using FESpace = typename Template_names<d>::FESpace;
  using Example_NSE = typename Template_names<d>::Example_NSE;
  
  /** @brief see the other constructor */
  CoupledNavierStokesStress(const TDomain& domain, const ParameterDatabase& param_db);
  
  /** @brief constructor
   * 
   * The domain must have been refined a couple of times already. On the
   * finest level the finite element spaces and functions as well as
   * matrices, solution and right hand side vectors are initialized.
   *
   * @param domain the computational domain to get the grid(s)
   * @param param_db parameters controlling this class
   * @param example The example to use
   */
  CoupledNavierStokesStress(const TDomain& domain, const ParameterDatabase& param_db,
                     const Example_CoupledNS_Stress& ex);
  
  static ParameterDatabase default_coupled_database();
  
  /** @brief Assemble those parts which do not contain nonlinearities
   * like a Stokes or Brinkman problem. When solving a Navier-Stokes problem,
   * then this must be called once before entering the nonlinear loop.
   */
  void assemble_linear_terms();
  
  // getters and setters
  const BlockFEMatrix & get_matrix() const
  { return this->systems.front().matrix; }
  BlockFEMatrix & get_matrix()
  { return this->systems.front().matrix; }
  
  const BlockVector & get_rhs() const
  { return this->systems.front().rhs; }
  BlockVector & get_rhs()
  { return this->systems.front().rhs; }
  
  const FEFunction & get_stress_xx() const
  {return this->systems.front().stress1;}
  FEFunction & get_stress_xx() 
  {return this->systems.front().stress1;}
  
  const FEFunction & get_stress_xy() const
  {return this->systems.front().stress2;}
  FEFunction & get_stress_xy() 
  {return this->systems.front().stress2;}
  
  const FEFunction & get_stress_yy() const
  {return this->systems.front().stress3;}
  FEFunction & get_stress_yy() 
  {return this->systems.front().stress3;}
  
  const FEVectFunct & get_velocity() const
  { return this->systems.front().u; }
  FEVectFunct & get_velocity()
  { return this->systems.front().u; }
  
  // try not to use this as it is not const
  FEFunction *get_velocity_component(int i);
  
  const FEFunction & get_pressure() const
  { return this->systems.front().p; }
  FEFunction & get_pressure()
  { return this->systems.front().p; }
  
  const FESpace & get_velocity_space() const
  { return *this->systems.front().velocity_space.get(); }
  const FESpace & get_pressure_space() const
  { return *this->systems.front().pressure_space.get(); }
  
  const BlockVector & get_solution() const
  { return this->systems.front().solution; }
  BlockVector & get_solution()
  { return this->systems.front().solution; }
  
  const LoopInfo& get_it_solver_info()
  {return solver.get_solver_loop_info();}
  
  unsigned int get_size() const
  { return this->systems.front().solution.length(); }
  
  const Example_CoupledNS_Stress & get_example() const
  { return example; }
protected:
  /// @brief default copy constructor (useful in derived classes)
  CoupledNavierStokesStress(const CoupledNavierStokesStress<d> &) = default;
  
  /** @brief store a complete system on a paticular grid.
   * 
   * This combines a matrix, rhs, solution, spaces and functions
   * needed to describe a Time Navier-Stokes problem in 2D
   */
  struct System_per_grid
  {
    /** @brief Finite element space for stress */
    std::shared_ptr<FESpace> stress_space;
    /** @brief Finite element space for velocity */
    std::shared_ptr<FESpace> velocity_space;
    /** @brief Finite element space for pressure */
    std::shared_ptr<FESpace> pressure_space;
    
    /** @brief system matrix
     *  T11    0     0    S11' S12' 0                    
     *  0    T22     0    S21' S22' 0
     *  0      0    T33   S31' S32' 0
     *  S11  S12    S13   A11  A12  B1T
     *  S21  S22    S23   A21  A22  B2T
     *  0      0     0    B1   B2   C              
     */
    BlockFEMatrix matrix;
    
    /** @brief mass matrix: this will be the standard mass matrix
     * for the standard Galerkin scheme. However, for the SUPG/RBVMS
     * schems, this includes all the terms which are related to the
     * time derivatives in the fully discrete scheme Eq: (45)
     * Ahmed, Rebollo, John and Rubino (2015)
     *  [ M11  M12  0 ]
     *  [ M21  M22  0 ]
     *  [ 0     0   0 ]
     * This additionaly created due to the struture of the
     * residual based Variational Multiscale Method:
     */
    // BlockFEMatrix mass_matrix;
    
    /** @brief right hand side vector*/
    BlockVector rhs;
    /** @brief solution vector */
    BlockVector solution;
    /** @brief Finite element function for stress */
    FEFunction stress1;
    FEFunction stress2;
    FEFunction stress3;
#ifdef __3D__
    FEFunction stress4;
    FEFunction stress5;
    FEFunction stress6;
#endif
    /** @brief Finite element function for velocity*/
    FEVectFunct u;
    /** @brief Finite element function for pressure*/
    FEFunction p;
    
    
    
    /** @brief constructor*/
    System_per_grid(const Example_CoupledNS_Stress& example, TCollection& coll,
                    std::tuple<int,int,int> order);
    
    /** @brief copy constructor*/
    System_per_grid(const System_per_grid&);
    
    //! Delete move constructor. No moves allowed.
    System_per_grid( System_per_grid&& ) = delete;
    
    //! Delete copy assignment operator. No copies allowed.
    System_per_grid& operator=( const System_per_grid& ) = delete;
    
    //! Default move assignment operator. No moves allowed.
    System_per_grid& operator=( System_per_grid&& ) = delete;
    
    //! Default destructor. Does most likely cause memory leaks.
    ~System_per_grid() = default;
  };
  
  /** @brief a local parameter database which controls this class
   * 
   *
   * The database given to the constructor will be merged into this one. Only
   * parameters which are of interest to this class are stored (and the
   * default ParMooN parameters). Note that this usually does not include
   * other parameters such as solver parameters. Those are only in the
   * Solver object.
   */
  ParameterDatabase db;
   
  /** @brief a complete system on each grid
  *
  * Note that the size of this deque is at least one and larger only in case
  * of multigrid.
  */
  std::deque<System_per_grid> systems;
  
  /** @brief Definition of the used example */
  Example_CoupledNS_Stress example;
  
  /** @brief a solver object which will solve the linear system
   * 
   * Storing it means that for a direct solver we also store the factorization
   * which is usually not necessary.
   */
  Solver<BlockFEMatrix, BlockVector> solver;
  
  /** @brief output object */
  DataWriter<d> outputWriter;
    
  /// this is for setting the discretization type globally, since the DISCTYPE 
  /// is removed fromt the global database but it's not simple/possible to use
  /// different space discretization in the LocalAssembling2D
  int space_disc_global;
  
  /** @brief set parameters in database
   * 
   * This functions checks if the parameters in the database are meaningful
   * and resets them otherwise. The hope is that after calling this function
   * this class is fully functional.
   *
   * If some parameters are set to unsupported values, an error occurs and
   * throws an exception.
   */
  void check_and_set_parameters();

  /** @brief write some information (number of cells, dofs, ...) */
  void output_problem_size_info() const;
    
  /** @brief get velocity and pressure space*/
  void get_velocity_pressure_orders(
    std::tuple<int, int, int> &velocity_pressure_orders);
  
  /// @brief set the spaces depending on disc types
  void set_arrays_stress(CoupledNavierStokesStress<d>::System_per_grid& s,
                  std::vector<const FESpace*> &spaces,
                  std::vector<const FESpace*>& spaces_rhs,
                  std::vector<FEFunction*>& functions);
  /// @brief set the spaces depending on disc types
  void set_arrays_ns(CoupledNavierStokesStress<d>::System_per_grid& s,
                  std::vector<const FESpace*> &spaces,
                  std::vector<const FESpace*>& spaces_rhs,
                  std::vector<FEFunction*>& functions);
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD; 
  /// @brief set the matrices and right hand side depending on the
  /// assemling routines, for stress part
  void set_matrices_rhs_stress(CoupledNavierStokesStress<d>::System_per_grid& s,
                        LocalAssembling_type type,
                        std::vector<SquareMatrixD*> &sqMat,
                        std::vector<MatrixD*> &reMat,
                        std::vector<double*> &rhs_array);
  /// @brief set the matrices and right hand side depending on the
  /// assemling routines, for navier stokes part
  void set_matrices_rhs_ns(CoupledNavierStokesStress<d>::System_per_grid& s,
                        LocalAssembling_type type,
                        std::vector<SquareMatrixD*> &sqMat,
                        std::vector<MatrixD*> &reMat,
                        std::vector<double*> &rhs_array);
  
  /// @brief this routines wraps up the call to Assemble2D
  void call_assembling_routine_stress(CoupledNavierStokesStress<d>::System_per_grid& s);
  
  /// @brief this routines wraps up the call to Assemble2D
  void call_assembling_routine_ns(CoupledNavierStokesStress<d>::System_per_grid& s, 
                                  LocalAssembling_type type);
};

#endif // COUPLEDNAVIERSTOKESSTRESS_H