/** ************************************************************************
 * 
 * @name         Time_LinElastic2D
 * @brief        store everything needed to solve a 2D linear elastic
 *               problem
 *               
 * @author       Najib Alia
 * @History      16.03.2017
**************************************************************************/

#ifndef __Time_LinElastic2D__
#define __Time_LinElastic2D__

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <FESpace2D.h>
#include <FEVectFunct2D.h>
#include <Solver.h>
#include <PostProcessing2D.h>

#include <Example_TimeLinElastic2D.h>

#include <MainUtilities.h>
#include <ParameterDatabase.h>
#include <vector>
#include <deque>
#include <utility>
#include <array>

using namespace std;

class Time_LinElastic2D
{
  protected:
    /** @brief store a complete system on a particular grid.
     *
     * This combines a matrix, rhs, solution, spaces and functions
     * needed to describe a Time Linear Elasticity problem in 2D
     */
    struct System_per_grid
    {
      /** @brief Finite element space for displacement */
      TFESpace2D fe_space_;
      /** @brief stiffness matrix K
       *  [ K11  K12 ]
       *  [ K21  K22 ]
       */
      BlockFEMatrix stiffness_matrix_;
      /** @brief mass matrix M
       *  [ M11  0   ]
       *  [ 0    M22 ]
       */
      BlockFEMatrix mass_matrix_;
      /** @brief right hand side vector*/
      BlockVector rhs_;
      /** @brief solution vector = displacement */
      BlockVector solution_;
      /** @brief Finite element vector for displacement*/
      TFEVectFunct2D u_;
      /** @brief Finite element function for Lame coefficient Lambda */
      BlockVector vector_lambda_;
      TFEFunction2D lambda_;
      /** @brief Finite element function for Lame coefficient Lambda */
      BlockVector vector_mu_;
      TFEFunction2D mu_;

      /** @brief constructor */
      System_per_grid(const Example_TimeLinElastic2D& example,
                      TCollection& coll);

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
    ParameterDatabase db_;

    /** @brief Definition of the used example */
    const Example_TimeLinElastic2D example_;

    /** @brief a solver object which will solve the linear system
     *
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver_;

    /** @brief class for output handling */
    PostProcessing2D outputwriter_;

    /** @brief a complete system on each grid
     *
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<System_per_grid> systems_;

    /** @brief store errors  */
//    std::vector<double> errors_;

    /** @brief write some information (number of cells, dofs, ...) */
    void output_problem_size_info() const;


  public:
    /** @brief constructor*/
    Time_LinElastic2D(const TDomain& domain,
                      const ParameterDatabase& param_db,
                      int reference_id = -4711);

    /** @brief constructor*/
    Time_LinElastic2D(const TDomain& domain,
                      const ParameterDatabase& param_db,
                      const Example_TimeLinElastic2D& ex,
                      int reference_id = -4711);

    /** @brief Assemble all the matrices and rhs before the time iterations
     *
     * This includes the assembling of: Stiff_matrix, Mass_Matrix,
     * rhs
     */
    void assemble_initial_time();



    /*************************************************************/
   /**
    * Special member functions mostly deleted
    * ...needs to be optimized
    */
   //! Delete copy constructor.
    Time_LinElastic2D(const Time_LinElastic2D&) = delete;

   //! Delete move constructor.
    Time_LinElastic2D(Time_LinElastic2D&&) = delete;

   //! Delete copy assignment operator.
    Time_LinElastic2D& operator=(const Time_LinElastic2D&) = delete;

   //! Delete move assignment operator.
    Time_LinElastic2D& operator=(Time_LinElastic2D&&) = delete;

   //! Default destructor. Most likely causes memory leaks.
   ~Time_LinElastic2D() = default;

   /*************************************************************/
   /* Getters and setters */
   const TFESpace2D & get_fe_space() const
   { return this->systems_.front().fe_space_; }
   const BlockVector& get_solution() const
   { return this->systems_.front().solution_; }
   const BlockVector& get_rhs() const
   { return this->systems_.front().rhs_; }
   const BlockVector& get_lamecoefficients_lambda() const
   { return this->systems_.front().vector_lambda_; }
   const BlockVector& get_lamecoefficients_mu() const
   { return this->systems_.front().vector_mu_; }
   const TFEVectFunct2D & get_displacement() const
   { return this->systems_.front().u_; }
   // try not to use this as it is not const
   TFEFunction2D *get_displacement_component(int i)
   { return (i==0) ? this->systems_.front().u_.GetComponent(0)
                   : this->systems_.front().u_.GetComponent(1); }
   unsigned int get_size() const
   { return this->systems_.front().solution_.length(); }
   const Example_TimeLinElastic2D & get_example() const
   { return example_; }
   const ParameterDatabase & get_db() const
   { return db_; }
   const BlockFEMatrix & get_stiff_matrix() const
   { return this->systems_.front().stiffness_matrix_; }
   const BlockFEMatrix & get_mass_matrix() const
   { return this->systems_.front().mass_matrix_; }

};

#endif // __Time_LinElastic2D__
